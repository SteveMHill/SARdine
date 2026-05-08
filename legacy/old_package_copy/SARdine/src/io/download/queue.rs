#![allow(dead_code, unused_variables)]
//! Download queue management system

use crate::types::{SarError, SarResult};
use log::debug;
use std::collections::{HashMap, VecDeque};
use std::sync::{Arc, Mutex};
use std::time::{Duration, Instant};

/// Helper macro to recover from poisoned mutexes
macro_rules! lock_or_recover {
    ($mutex:expr) => {
        $mutex.lock().unwrap_or_else(|poisoned| {
            log::warn!("Mutex poisoned, recovering: {}", poisoned);
            poisoned.into_inner()
        })
    };
}

// Constants
const DEFAULT_TASK_CLEANUP_AGE_SECONDS: u64 = 3600; // 1 hour
const MIN_PRODUCT_ID_LENGTH: usize = 50;

/// Download task priority
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum Priority {
    Low = 0,
    Normal = 1,
    High = 2,
    Urgent = 3,
}

/// Download task status
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum TaskStatus {
    Pending,
    Downloading,
    Completed,
    Failed,
    Paused,
    Cancelled,
}

/// Download task
#[derive(Debug, Clone)]
pub struct DownloadTask {
    pub id: String,
    pub url: String,
    pub output_path: String,
    pub priority: Priority,
    pub status: TaskStatus,
    pub created_at: Instant,
    pub started_at: Option<Instant>,
    pub completed_at: Option<Instant>,
    pub downloaded_bytes: u64,
    pub total_bytes: Option<u64>,
    pub error: Option<String>,
    pub retry_count: u32,
    pub max_retries: u32,
}

impl DownloadTask {
    pub fn new(id: String, url: String, output_path: String, priority: Priority) -> Self {
        Self {
            id,
            url,
            output_path,
            priority,
            status: TaskStatus::Pending,
            created_at: Instant::now(),
            started_at: None,
            completed_at: None,
            downloaded_bytes: 0,
            total_bytes: None,
            error: None,
            retry_count: 0,
            max_retries: 3,
        }
    }

    pub fn progress_percentage(&self) -> f64 {
        if let Some(total) = self.total_bytes {
            if total > 0 {
                return (self.downloaded_bytes as f64 / total as f64) * 100.0;
            }
        }
        0.0
    }
}

/// Download queue manager
pub struct DownloadQueue {
    tasks: Arc<Mutex<HashMap<String, DownloadTask>>>,
    queue: Arc<Mutex<VecDeque<String>>>, // Task IDs in priority order
    active_downloads: Arc<Mutex<usize>>,
    max_concurrent: usize,
}

impl DownloadQueue {
    pub fn new(max_concurrent: usize) -> Self {
        Self {
            tasks: Arc::new(Mutex::new(HashMap::new())),
            queue: Arc::new(Mutex::new(VecDeque::new())),
            active_downloads: Arc::new(Mutex::new(0)),
            max_concurrent,
        }
    }

    /// Add a task to the queue
    pub fn add_task(&self, task: DownloadTask) -> SarResult<()> {
        let id = task.id.clone();
        let priority = task.priority;

        {
            let mut tasks = lock_or_recover!(self.tasks);
            tasks.insert(id.clone(), task);
        }

        // Insert into queue based on priority
        // FIXED: Avoid nested locks by collecting task IDs first
        let task_ids: Vec<String> = {
            let queue = lock_or_recover!(self.queue);
            queue.iter().cloned().collect()
        };

        // Find insertion position without holding queue lock
        let insert_pos = {
            let tasks = lock_or_recover!(self.tasks);
            task_ids
                .iter()
                .position(|task_id| {
                    tasks
                        .get(task_id)
                        .map(|t| t.priority < priority)
                        .unwrap_or(false)
                })
                .unwrap_or(task_ids.len())
        };

        // Insert into queue
        let mut queue = lock_or_recover!(self.queue);
        queue.insert(insert_pos, id);
        Ok(())
    }

    /// Get next task to process
    pub fn next_task(&self) -> Option<DownloadTask> {
        let active = *lock_or_recover!(self.active_downloads);

        if active >= self.max_concurrent {
            return None;
        }

        // FIXED: Avoid nested locks by collecting task IDs first
        let task_ids: Vec<(usize, String)> = {
            let queue = lock_or_recover!(self.queue);
            queue
                .iter()
                .enumerate()
                .map(|(i, id)| (i, id.clone()))
                .collect()
        };

        // Find next pending task
        for (i, task_id) in task_ids {
            let mut tasks = lock_or_recover!(self.tasks);
            if let Some(task) = tasks.get_mut(&task_id) {
                if task.status == TaskStatus::Pending {
                    task.status = TaskStatus::Downloading;
                    task.started_at = Some(Instant::now());
                    let task_clone = task.clone();
                    drop(tasks);

                    // Remove from queue and update active count
                    let mut queue = lock_or_recover!(self.queue);
                    queue.remove(i);
                    *lock_or_recover!(self.active_downloads) += 1;
                    return Some(task_clone);
                }
            }
        }

        None
    }

    /// Update task progress
    pub fn update_progress(&self, task_id: &str, downloaded: u64, total: Option<u64>) {
        let mut tasks = lock_or_recover!(self.tasks);
        if let Some(task) = tasks.get_mut(task_id) {
            task.downloaded_bytes = downloaded;
            if let Some(total_bytes) = total {
                task.total_bytes = Some(total_bytes);
            }
        }
    }

    /// Mark task as completed
    pub fn complete_task(&self, task_id: &str) {
        {
            let mut tasks = lock_or_recover!(self.tasks);
            if let Some(task) = tasks.get_mut(task_id) {
                task.status = TaskStatus::Completed;
                task.completed_at = Some(Instant::now());
            }
        }
        let mut active = lock_or_recover!(self.active_downloads);
        *active = active.saturating_sub(1);
    }

    /// Mark task as failed
    pub fn fail_task(&self, task_id: &str, error: String) {
        let should_retry = {
            let mut tasks = lock_or_recover!(self.tasks);
            if let Some(task) = tasks.get_mut(task_id) {
                task.status = TaskStatus::Failed;
                task.error = Some(error);
                task.retry_count += 1;

                // Retry if under max retries
                if task.retry_count < task.max_retries {
                    task.status = TaskStatus::Pending;
                    task.error = None;
                    true
                } else {
                    false
                }
            } else {
                false
            }
        };

        if should_retry {
            // Re-add to queue
            let mut queue = lock_or_recover!(self.queue);
            queue.push_back(task_id.to_string());
        }

        let mut active = lock_or_recover!(self.active_downloads);
        *active = active.saturating_sub(1);
    }

    /// Pause a task
    pub fn pause_task(&self, task_id: &str) -> SarResult<()> {
        let was_downloading = {
            let mut tasks = lock_or_recover!(self.tasks);
            if let Some(task) = tasks.get_mut(task_id) {
                let was_downloading = task.status == TaskStatus::Downloading;
                // Can pause both Downloading and Pending tasks
                if was_downloading || task.status == TaskStatus::Pending {
                    task.status = TaskStatus::Paused;
                }
                Ok(was_downloading)
            } else {
                Err(SarError::Processing(format!("Task not found: {}", task_id)))
            }
        }?;

        if was_downloading {
            let mut active = lock_or_recover!(self.active_downloads);
            *active = active.saturating_sub(1);
        }
        Ok(())
    }

    /// Resume a paused task
    pub fn resume_task(&self, task_id: &str) -> SarResult<()> {
        {
            let mut tasks = lock_or_recover!(self.tasks);
            if let Some(task) = tasks.get_mut(task_id) {
                if task.status == TaskStatus::Paused {
                    task.status = TaskStatus::Pending;
                } else {
                    return Err(SarError::Processing(format!(
                        "Task is not paused: {}",
                        task_id
                    )));
                }
            } else {
                return Err(SarError::Processing(format!("Task not found: {}", task_id)));
            }
        }

        // Re-add to queue
        let mut queue = lock_or_recover!(self.queue);
        queue.push_back(task_id.to_string());
        Ok(())
    }

    /// Cancel a task
    pub fn cancel_task(&self, task_id: &str) -> SarResult<()> {
        let was_downloading = {
            let mut tasks = lock_or_recover!(self.tasks);
            if let Some(task) = tasks.get_mut(task_id) {
                // FIXED: Check status before changing it
                let was_downloading = task.status == TaskStatus::Downloading;
                task.status = TaskStatus::Cancelled;
                Ok(was_downloading)
            } else {
                Err(SarError::Processing(format!("Task not found: {}", task_id)))
            }
        }?;

        if was_downloading {
            let mut active = lock_or_recover!(self.active_downloads);
            *active = active.saturating_sub(1);
        }
        Ok(())
    }

    /// Get task status
    pub fn get_task(&self, task_id: &str) -> Option<DownloadTask> {
        lock_or_recover!(self.tasks).get(task_id).cloned()
    }

    /// Get all tasks
    pub fn get_all_tasks(&self) -> Vec<DownloadTask> {
        lock_or_recover!(self.tasks).values().cloned().collect()
    }

    /// Get queue statistics
    pub fn get_stats(&self) -> (usize, usize, usize, usize) {
        let tasks = lock_or_recover!(self.tasks);
        let pending = tasks
            .values()
            .filter(|t| t.status == TaskStatus::Pending)
            .count();
        let downloading = tasks
            .values()
            .filter(|t| t.status == TaskStatus::Downloading)
            .count();
        let completed = tasks
            .values()
            .filter(|t| t.status == TaskStatus::Completed)
            .count();
        let failed = tasks
            .values()
            .filter(|t| t.status == TaskStatus::Failed)
            .count();
        (pending, downloading, completed, failed)
    }

    /// Remove old completed/failed tasks to prevent memory growth
    pub fn cleanup_old_tasks(&self, max_age_seconds: Option<u64>) {
        let max_age = max_age_seconds.unwrap_or(DEFAULT_TASK_CLEANUP_AGE_SECONDS);
        let cutoff = Instant::now() - Duration::from_secs(max_age);
        let mut tasks = lock_or_recover!(self.tasks);

        let initial_count = tasks.len();
        tasks.retain(|_, task| {
            match task.status {
                TaskStatus::Completed | TaskStatus::Failed => task
                    .completed_at
                    .or(task.started_at)
                    .map(|t| t > cutoff)
                    .unwrap_or(true),
                _ => true, // Keep pending/downloading/paused/cancelled tasks
            }
        });

        let removed = initial_count - tasks.len();
        if removed > 0 {
            debug!("Cleaned up {} old tasks from queue", removed);
        }
    }
}
