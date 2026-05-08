//! Download progress reporting

/// Progress callback function type
/// Arguments: (downloaded_bytes, total_bytes, bytes_per_second)
pub type ProgressCallback = Box<dyn Fn(u64, u64, f64) + Send + Sync>;

/// Progress state for tracking download progress
#[derive(Debug, Clone)]
pub struct ProgressState {
    pub downloaded: u64,
    pub total: u64,
    pub bytes_per_second: f64,
    pub percentage: f64,
    pub elapsed_seconds: f64,
    pub estimated_remaining_seconds: f64,
}

impl ProgressState {
    pub fn new(total: u64) -> Self {
        Self {
            downloaded: 0,
            total,
            bytes_per_second: 0.0,
            percentage: 0.0,
            elapsed_seconds: 0.0,
            estimated_remaining_seconds: 0.0,
        }
    }

    pub fn update(&mut self, downloaded: u64, elapsed_seconds: f64) {
        self.downloaded = downloaded;
        self.elapsed_seconds = elapsed_seconds;

        if elapsed_seconds > 0.0 {
            self.bytes_per_second = downloaded as f64 / elapsed_seconds;
        }

        if self.total > 0 {
            self.percentage = (downloaded as f64 / self.total as f64) * 100.0;
            if self.bytes_per_second > 0.0 {
                let remaining = self.total.saturating_sub(downloaded);
                self.estimated_remaining_seconds = remaining as f64 / self.bytes_per_second;
            }
        }
    }
}
