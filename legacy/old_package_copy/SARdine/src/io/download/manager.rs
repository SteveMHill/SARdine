#![allow(dead_code, unused_variables)]
//! Download manager - unified coordinator for all download types

use crate::types::{SarError, SarResult};
use chrono::{DateTime, Utc};
use std::collections::HashMap;
use std::path::PathBuf;

use super::cache::CacheManager;
use super::dem::DEMDownloader;
use super::orbit::OrbitDownloader;
use super::products::ProductDownloader;
use super::progress::ProgressCallback;
use super::providers::ProviderRegistry;
use super::queue::{DownloadQueue, Priority};
use std::sync::Arc;

/// Download configuration
#[derive(Debug, Clone)]
pub struct DownloadConfig {
    pub cache_dir: PathBuf,
    pub enabled_sources: Vec<String>, // ["esa", "asf", "aws", "cds"]
    pub max_parallel_downloads: usize,
    pub retry_attempts: u32,
    pub timeout_seconds: u64,
    pub verify_checksums: bool,
    pub auto_download: bool, // Auto-download missing data in pipeline
    pub provider_configs: HashMap<String, ProviderConfig>,
}

/// Provider-specific configuration
#[derive(Debug, Clone)]
pub struct ProviderConfig {
    pub username: Option<String>,
    pub password: Option<String>,
    pub token: Option<String>,
    pub base_url: Option<String>,
}

impl Default for DownloadConfig {
    fn default() -> Self {
        let cache_dir = std::env::var("SARDINE_DOWNLOAD_CACHE")
            .map(PathBuf::from)
            .unwrap_or_else(|_| {
                dirs::home_dir()
                    .unwrap_or_else(|| PathBuf::from("."))
                    .join(".sardine")
                    .join("cache")
            });

        let enabled_sources = std::env::var("SARDINE_DOWNLOAD_SOURCES")
            .map(|s| s.split(',').map(|x| x.trim().to_string()).collect())
            .unwrap_or_else(|_| vec!["aws".to_string(), "esa".to_string()]);

        Self {
            cache_dir,
            enabled_sources,
            max_parallel_downloads: 3,
            retry_attempts: 3,
            timeout_seconds: 3600, // 1 hour for large files
            verify_checksums: true,
            auto_download: false,
            provider_configs: HashMap::new(),
        }
    }
}

/// Download result
#[derive(Debug)]
pub struct DownloadResult {
    pub path: PathBuf,
    pub size: u64,
    pub source: String,
}

/// Unified download manager
pub struct DownloadManager {
    config: DownloadConfig,
    cache: Arc<CacheManager>,
    providers: Arc<ProviderRegistry>,
    product_downloader: ProductDownloader,
    orbit_downloader: OrbitDownloader,
    dem_downloader: DEMDownloader,
    queue: Arc<DownloadQueue>,
}

impl DownloadManager {
    /// Create a new download manager
    pub fn new(config: DownloadConfig) -> SarResult<Self> {
        let cache = Arc::new(CacheManager::new(&config.cache_dir)?);
        let providers = Arc::new(ProviderRegistry::new(&config)?);

        let product_downloader = ProductDownloader::new(&cache, &providers, &config)?;
        let orbit_downloader = OrbitDownloader::new(&cache, &providers, &config)?;
        let dem_downloader = DEMDownloader::new(&cache, &providers, &config)?;
        let queue = Arc::new(DownloadQueue::new(config.max_parallel_downloads));

        Ok(Self {
            config,
            cache,
            providers,
            product_downloader,
            orbit_downloader,
            dem_downloader,
            queue,
        })
    }

    /// Download a Sentinel-1 product
    pub fn download_product(
        &mut self,
        product_id: &str,
        progress: Option<ProgressCallback>,
    ) -> SarResult<DownloadResult> {
        // Check cache first
        if let Some(entry) = self.cache.get_product(product_id) {
            if CacheManager::validate_cached_file(&entry.path)? {
                log::info!("Using cached product: {}", entry.path.display());
                return Ok(DownloadResult {
                    path: entry.path.clone(),
                    size: entry.size,
                    source: "cache".to_string(),
                });
            }
        }

        // Download product
        let result = self.product_downloader.download(product_id, progress)?;

        // Add to cache
        let cache_entry = CacheManager::get_cache_entry(&result.path)?;
        self.cache
            .add_product(product_id.to_string(), cache_entry.clone());

        Ok(result)
    }

    /// Download orbit file for a product
    pub fn download_orbit(
        &mut self,
        product_id: &str,
        start_time: DateTime<Utc>,
        progress: Option<ProgressCallback>,
    ) -> SarResult<DownloadResult> {
        // Check cache first
        let orbit_key = format!("{}_{}", product_id, start_time.format("%Y%m%d"));
        if let Some(entry) = self.cache.get_orbit(&orbit_key) {
            if CacheManager::validate_cached_file(&entry.path)? {
                log::info!("Using cached orbit file: {}", entry.path.display());
                return Ok(DownloadResult {
                    path: entry.path.clone(),
                    size: entry.size,
                    source: "cache".to_string(),
                });
            }
        }

        // Download orbit
        let result = self
            .orbit_downloader
            .download(product_id, start_time, progress)?;

        // Add to cache
        let cache_entry = CacheManager::get_cache_entry(&result.path)?;
        self.cache.add_orbit(orbit_key, cache_entry.clone());

        Ok(result)
    }

    /// Download DEM tiles for a bounding box
    pub fn download_dem(
        &mut self,
        bbox: &crate::types::BoundingBox,
        progress: Option<ProgressCallback>,
    ) -> SarResult<Vec<DownloadResult>> {
        let results = self.dem_downloader.download_tiles(bbox, progress)?;

        // Add to cache
        for result in &results {
            if let Some(file_name) = result.path.file_stem().and_then(|s| s.to_str()) {
                if let Ok(cache_entry) = CacheManager::get_cache_entry(&result.path) {
                    self.cache.add_dem_tile(file_name.to_string(), cache_entry);
                }
            }
        }

        Ok(results)
    }

    /// Search for products (when provider supports it)
    pub fn search_products(
        &self,
        params: &super::products::SearchParams,
    ) -> SarResult<Vec<super::products::ProductMetadata>> {
        self.product_downloader.search(params)
    }

    /// Download multiple products in parallel
    pub fn download_products_batch(
        &mut self,
        product_ids: &[String],
    ) -> SarResult<Vec<(String, DownloadResult)>> {
        use rayon::prelude::*;
        use std::sync::Mutex;

        // Use the product downloader directly (it's thread-safe)
        let product_downloader = &self.product_downloader;
        let cache = self.cache.clone();

        let results: Mutex<Vec<(String, SarResult<DownloadResult>)>> = Mutex::new(Vec::new());

        product_ids.par_iter().for_each(|product_id| {
            // Check cache first
            let result = if let Some(entry) = cache.get_product(product_id) {
                if CacheManager::validate_cached_file(&entry.path).unwrap_or(false) {
                    Ok(DownloadResult {
                        path: entry.path.clone(),
                        size: entry.size,
                        source: "cache".to_string(),
                    })
                } else {
                    product_downloader.download(product_id, None)
                }
            } else {
                product_downloader.download(product_id, None)
            };

            // Add to cache if successful
            if let Ok(ref download_result) = result {
                if let Ok(cache_entry) = CacheManager::get_cache_entry(&download_result.path) {
                    cache.add_product(product_id.clone(), cache_entry);
                }
            }

            let mut results_guard = results.lock().unwrap_or_else(|poisoned| {
                log::warn!("Mutex poisoned in batch download, recovering");
                poisoned.into_inner()
            });
            results_guard.push((product_id.clone(), result));
        });

        let results = results.into_inner().unwrap_or_else(|poisoned| {
            log::warn!("Mutex poisoned when extracting batch results, recovering");
            poisoned.into_inner()
        });
        let mut successes = Vec::new();
        let mut errors = Vec::new();

        for (product_id, result) in results {
            match result {
                Ok(download_result) => successes.push((product_id, download_result)),
                Err(e) => errors.push((product_id, e)),
            }
        }

        if !errors.is_empty() {
            log::warn!("Some batch downloads failed: {} errors", errors.len());
            for (id, err) in &errors {
                log::warn!("  {}: {}", id, err);
            }
        }

        if successes.is_empty() {
            return Err(SarError::Processing(format!(
                "All {} downloads failed",
                product_ids.len()
            )));
        }

        Ok(successes)
    }

    /// Get cache manager reference
    pub fn cache(&self) -> &CacheManager {
        &self.cache
    }

    /// Get cache manager reference (Arc is already interior mutable)
    pub fn cache_ref(&self) -> &CacheManager {
        &self.cache
    }

    /// Check if auto-download is enabled
    pub fn auto_download(&self) -> bool {
        self.config.auto_download
    }

    /// Get download queue reference
    pub fn queue(&self) -> &DownloadQueue {
        &self.queue
    }

    /// Add product download to queue
    pub fn queue_product_download(
        &self,
        product_id: String,
        priority: Priority,
    ) -> SarResult<String> {
        use super::queue::DownloadTask;
        let task_id = format!("product_{}", product_id);
        let output_path = self
            .cache
            .products_dir()
            .join(format!("{}.zip", product_id))
            .to_string_lossy()
            .to_string();

        // Resolve URL first (simplified - would need actual URL resolution)
        let url = format!("placeholder://{}", product_id); // Would be resolved later

        let task = DownloadTask::new(task_id.clone(), url, output_path, priority);
        self.queue.add_task(task)?;
        Ok(task_id)
    }
}
