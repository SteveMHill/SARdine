//! Unified caching system for all download types

use crate::types::{SarError, SarResult};
use std::collections::HashMap;
use std::fs;
use std::path::{Path, PathBuf};
use std::sync::{Arc, Mutex};
use std::time::SystemTime;

/// Cache entry metadata
#[derive(Debug, Clone)]
pub struct CacheEntry {
    pub path: PathBuf,
    pub size: u64,
    pub modified: SystemTime,
    pub checksum: Option<String>,
}

/// Unified cache manager for products, orbits, and DEMs
#[derive(Clone)]
pub struct CacheManager {
    base_dir: PathBuf,
    product_cache: Arc<Mutex<HashMap<String, CacheEntry>>>,
    orbit_cache: Arc<Mutex<HashMap<String, CacheEntry>>>,
    dem_cache: Arc<Mutex<HashMap<String, CacheEntry>>>,
}

impl CacheManager {
    /// Create a new cache manager
    pub fn new(base_dir: &Path) -> SarResult<Self> {
        let base_dir = base_dir.to_path_buf();
        fs::create_dir_all(&base_dir).map_err(SarError::Io)?;

        // Create subdirectories
        fs::create_dir_all(base_dir.join("products")).map_err(SarError::Io)?;
        fs::create_dir_all(base_dir.join("orbits")).map_err(SarError::Io)?;
        fs::create_dir_all(base_dir.join("dem")).map_err(SarError::Io)?;

        Ok(Self {
            base_dir,
            product_cache: Arc::new(Mutex::new(HashMap::new())),
            orbit_cache: Arc::new(Mutex::new(HashMap::new())),
            dem_cache: Arc::new(Mutex::new(HashMap::new())),
        })
    }

    /// Get cache directory for products
    pub fn products_dir(&self) -> PathBuf {
        self.base_dir.join("products")
    }

    /// Get cache directory for orbits
    pub fn orbits_dir(&self) -> PathBuf {
        self.base_dir.join("orbits")
    }

    /// Get cache directory for DEMs
    pub fn dem_dir(&self) -> PathBuf {
        self.base_dir.join("dem")
    }

    /// Check if a product is cached
    pub fn get_product(&self, product_id: &str) -> Option<CacheEntry> {
        self.product_cache
            .lock()
            .unwrap_or_else(|poisoned| poisoned.into_inner())
            .get(product_id)
            .cloned()
    }

    /// Check if an orbit file is cached
    pub fn get_orbit(&self, orbit_key: &str) -> Option<CacheEntry> {
        self.orbit_cache
            .lock()
            .unwrap_or_else(|poisoned| poisoned.into_inner())
            .get(orbit_key)
            .cloned()
    }

    /// Check if a DEM tile is cached
    pub fn get_dem_tile(&self, tile_name: &str) -> Option<CacheEntry> {
        self.dem_cache
            .lock()
            .unwrap_or_else(|poisoned| poisoned.into_inner())
            .get(tile_name)
            .cloned()
    }

    /// Add a product to cache
    pub fn add_product(&self, product_id: String, entry: CacheEntry) {
        if let Ok(mut cache) = self.product_cache.lock() {
            cache.insert(product_id, entry);
        }
    }

    /// Add an orbit file to cache
    pub fn add_orbit(&self, orbit_key: String, entry: CacheEntry) {
        if let Ok(mut cache) = self.orbit_cache.lock() {
            cache.insert(orbit_key, entry);
        }
    }

    /// Add a DEM tile to cache
    pub fn add_dem_tile(&self, tile_name: String, entry: CacheEntry) {
        if let Ok(mut cache) = self.dem_cache.lock() {
            cache.insert(tile_name, entry);
        }
    }

    /// Validate a cached file exists and is readable
    pub fn validate_cached_file(path: &Path) -> SarResult<bool> {
        if !path.exists() {
            return Ok(false);
        }

        let metadata = fs::metadata(path).map_err(SarError::Io)?;
        if metadata.len() == 0 {
            return Ok(false);
        }

        Ok(true)
    }

    /// Get cache entry for a file
    pub fn get_cache_entry(path: &Path) -> SarResult<CacheEntry> {
        let metadata = fs::metadata(path).map_err(SarError::Io)?;
        Ok(CacheEntry {
            path: path.to_path_buf(),
            size: metadata.len(),
            modified: metadata.modified().map_err(SarError::Io)?,
            checksum: None,
        })
    }

    /// Scan cache directory and populate cache entries
    pub fn scan_cache(&mut self) -> SarResult<()> {
        // Scan products
        let products_dir = self.products_dir();
        if products_dir.exists() {
            let mut cache = self
                .product_cache
                .lock()
                .unwrap_or_else(|poisoned| poisoned.into_inner());
            for entry in fs::read_dir(&products_dir).map_err(SarError::Io)? {
                let entry = entry.map_err(SarError::Io)?;
                let path = entry.path();
                if path.is_file() {
                    if let Some(file_name) = path.file_stem().and_then(|s| s.to_str()) {
                        if let Ok(cache_entry) = Self::get_cache_entry(&path) {
                            cache.insert(file_name.to_string(), cache_entry);
                        }
                    }
                }
            }
        }

        // Scan orbits
        let orbits_dir = self.orbits_dir();
        if orbits_dir.exists() {
            let mut cache = self
                .orbit_cache
                .lock()
                .unwrap_or_else(|poisoned| poisoned.into_inner());
            for entry in fs::read_dir(&orbits_dir).map_err(SarError::Io)? {
                let entry = entry.map_err(SarError::Io)?;
                let path = entry.path();
                if path.is_file() {
                    if let Some(file_name) = path.file_stem().and_then(|s| s.to_str()) {
                        if let Ok(cache_entry) = Self::get_cache_entry(&path) {
                            cache.insert(file_name.to_string(), cache_entry);
                        }
                    }
                }
            }
        }

        // Scan DEM tiles
        let dem_dir = self.dem_dir();
        if dem_dir.exists() {
            let mut cache = self
                .dem_cache
                .lock()
                .unwrap_or_else(|poisoned| poisoned.into_inner());
            for entry in fs::read_dir(&dem_dir).map_err(SarError::Io)? {
                let entry = entry.map_err(SarError::Io)?;
                let path = entry.path();
                if path.is_file() {
                    if let Some(file_name) = path.file_stem().and_then(|s| s.to_str()) {
                        if let Ok(cache_entry) = Self::get_cache_entry(&path) {
                            cache.insert(file_name.to_string(), cache_entry);
                        }
                    }
                }
            }
        }

        Ok(())
    }

    /// Clear cache (optional: specify type)
    pub fn clear_cache(&self, cache_type: Option<&str>) -> SarResult<()> {
        match cache_type {
            Some("products") => {
                let products_dir = self.products_dir();
                if products_dir.exists() {
                    fs::remove_dir_all(&products_dir).map_err(SarError::Io)?;
                    fs::create_dir_all(&products_dir).map_err(SarError::Io)?;
                }
                if let Ok(mut cache) = self.product_cache.lock() {
                    cache.clear();
                }
            }
            Some("orbits") => {
                let orbits_dir = self.orbits_dir();
                if orbits_dir.exists() {
                    fs::remove_dir_all(&orbits_dir).map_err(SarError::Io)?;
                    fs::create_dir_all(&orbits_dir).map_err(SarError::Io)?;
                }
                if let Ok(mut cache) = self.orbit_cache.lock() {
                    cache.clear();
                }
            }
            Some("dem") => {
                let dem_dir = self.dem_dir();
                if dem_dir.exists() {
                    fs::remove_dir_all(&dem_dir).map_err(SarError::Io)?;
                    fs::create_dir_all(&dem_dir).map_err(SarError::Io)?;
                }
                if let Ok(mut cache) = self.dem_cache.lock() {
                    cache.clear();
                }
            }
            None => {
                // Clear all
                self.clear_cache(Some("products"))?;
                self.clear_cache(Some("orbits"))?;
                self.clear_cache(Some("dem"))?;
            }
            _ => {
                return Err(SarError::Processing(format!(
                    "Unknown cache type: {}",
                    cache_type.unwrap()
                )));
            }
        }
        Ok(())
    }
}
