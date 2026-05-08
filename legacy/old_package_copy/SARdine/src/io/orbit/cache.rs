use std::path::{Path, PathBuf};

use chrono::{DateTime, Utc};

use crate::types::{OrbitData, SarError, SarResult};

use super::interp::validate_orbit_coverage;
use super::{OrbitReader, OrbitType};

pub struct OrbitCache {
    cache_dir: PathBuf,
}

impl OrbitCache {
    pub fn new(cache_dir: PathBuf) -> Self {
        Self { cache_dir }
    }

    pub fn get_cached_orbit(&self, product_id: &str, orbit_type: OrbitType) -> Option<PathBuf> {
        let filename = format!("{}_{}.EOF", product_id, orbit_type);
        let path = self.cache_dir.join(filename);
        if path.exists() {
            Some(path)
        } else {
            None
        }
    }

    pub fn cache_orbit(
        &self,
        product_id: &str,
        orbit_type: OrbitType,
        content: &str,
    ) -> SarResult<PathBuf> {
        std::fs::create_dir_all(&self.cache_dir)?;
        let filename = format!("{}_{}.EOF", product_id, orbit_type);
        let path = self.cache_dir.join(filename);
        std::fs::write(&path, content)?;
        Ok(path)
    }
}

pub struct OrbitManager {
    cache: OrbitCache,
}

impl OrbitManager {
    pub fn new(cache_dir: PathBuf) -> Self {
        Self {
            cache: OrbitCache::new(cache_dir),
        }
    }

    pub fn get_orbit_data(
        &self,
        product_id: &str,
        start_time: DateTime<Utc>,
    ) -> SarResult<OrbitData> {
        let orbit_type = OrbitReader::determine_orbit_type(start_time);

        // FIXED: Check cache first (in main cache directory)
        if let Some(cached_path) = self.cache.get_cached_orbit(product_id, orbit_type) {
            log::info!("Using cached orbit file: {}", cached_path.display());
            match OrbitReader::read_orbit_file(&cached_path) {
                Ok(orbit) => {
                    // CRITICAL: Validate temporal coverage to prevent wrong orbit usage
                    if let Err(e) = validate_orbit_coverage(&orbit, start_time) {
                        log::error!("❌ SAFETY: Cached orbit file failed temporal validation: {}. This could corrupt results!", e);
                        return Err(SarError::Processing(format!(
                            "Cached orbit file {} failed temporal validation: {}. Remove invalid cache file and download correct orbit.",
                            cached_path.display(), e
                        )));
                    } else {
                        log::info!(
                            "✅ Cached orbit file validated: covers target time {}",
                            start_time
                        );
                        return Ok(orbit);
                    }
                }
                Err(e) => {
                    log::error!("❌ Failed to read cached orbit file: {}", e);
                    return Err(e);
                }
            }
        }

        // FIXED: Also check orbits/ subdirectory (where DownloadManager stores files)
        let orbits_subdir = self.cache.cache_dir.join("orbits");
        let filename = format!("{}_{}.EOF", product_id, orbit_type);
        let orbits_path = orbits_subdir.join(&filename);
        if orbits_path.exists() {
            log::info!(
                "Using orbit file from orbits/ subdirectory: {}",
                orbits_path.display()
            );
            match OrbitReader::read_orbit_file(&orbits_path) {
                Ok(orbit) => {
                    // CRITICAL: Validate temporal coverage to prevent wrong orbit usage
                    if let Err(e) = validate_orbit_coverage(&orbit, start_time) {
                        log::error!("❌ SAFETY: Orbits subdirectory file failed temporal validation: {}. This could corrupt results!", e);
                        return Err(SarError::Processing(format!(
                            "Orbit file in orbits/ subdirectory {} failed temporal validation: {}. Remove invalid file and download correct orbit.",
                            orbits_path.display(), e
                        )));
                    } else {
                        log::info!(
                            "✅ Orbits subdirectory file validated: covers target time {}",
                            start_time
                        );
                        return Ok(orbit);
                    }
                }
                Err(e) => {
                    log::error!(
                        "❌ Failed to read orbit file from orbits/ subdirectory: {}",
                        e
                    );
                    return Err(e);
                }
            }
        }

        // FIXED: Check for RESORB fallback if we were looking for POEORB
        if matches!(orbit_type, OrbitType::POEORB) {
            let resorb_filename = format!("{}_{}.EOF", product_id, OrbitType::RESORB);
            let resorb_path = self.cache.cache_dir.join(&resorb_filename);
            if resorb_path.exists() {
                log::info!(
                    "Using cached RESORB file as fallback: {}",
                    resorb_path.display()
                );
                match OrbitReader::read_orbit_file(&resorb_path) {
                    Ok(orbit) => {
                        // CRITICAL: Validate temporal coverage to prevent wrong orbit usage
                        if let Err(e) = validate_orbit_coverage(&orbit, start_time) {
                            log::error!("❌ SAFETY: Cached RESORB fallback failed temporal validation: {}. This could corrupt results!", e);
                            return Err(SarError::Processing(format!(
                                "Cached RESORB fallback {} failed temporal validation: {}. Remove invalid cache file and download correct orbit.",
                                resorb_path.display(), e
                            )));
                        } else {
                            log::info!(
                                "✅ Cached RESORB fallback validated: covers target time {}",
                                start_time
                            );
                            return Ok(orbit);
                        }
                    }
                    Err(e) => {
                        log::error!("❌ Failed to read cached RESORB fallback: {}", e);
                        // Continue to try orbits/ subdirectory
                    }
                }
            }
            // Also check orbits/ subdirectory for RESORB
            let resorb_orbits_path = orbits_subdir.join(&resorb_filename);
            if resorb_orbits_path.exists() {
                log::info!(
                    "Using RESORB file from orbits/ subdirectory: {}",
                    resorb_orbits_path.display()
                );
                match OrbitReader::read_orbit_file(&resorb_orbits_path) {
                    Ok(orbit) => {
                        // CRITICAL: Validate temporal coverage to prevent wrong orbit usage
                        if let Err(e) = validate_orbit_coverage(&orbit, start_time) {
                            log::error!("❌ SAFETY: RESORB from orbits/ subdirectory failed temporal validation: {}. This could corrupt results!", e);
                            return Err(SarError::Processing(format!(
                                "RESORB file in orbits/ subdirectory {} failed temporal validation: {}. Remove invalid file and download correct orbit.",
                                resorb_orbits_path.display(), e
                            )));
                        } else {
                            log::info!("✅ RESORB from orbits/ subdirectory validated: covers target time {}", start_time);
                            return Ok(orbit);
                        }
                    }
                    Err(e) => {
                        log::error!("❌ Failed to read RESORB from orbits/ subdirectory: {}", e);
                        // Continue to final error
                    }
                }
            }
        }

        // FIXED: Don't try to download using old method - return error instructing to use DownloadManager
        return Err(crate::types::SarError::Processing(format!(
            "Orbit file not found in cache: {}. Please use DownloadManager to download the orbit file first. Expected locations:\n  - {}\n  - {}",
            filename,
            self.cache.cache_dir.join(&filename).display(),
            orbits_path.display()
        )));
    }

    pub fn has_orbit_cached(&self, product_id: &str, orbit_type: OrbitType) -> bool {
        self.cache
            .get_cached_orbit(product_id, orbit_type)
            .is_some()
    }

    pub fn get_cache_dir(&self) -> &Path {
        &self.cache.cache_dir
    }

    pub fn download_and_cache_orbit_public(
        &self,
        product_id: &str,
        start_time: DateTime<Utc>,
        orbit_type: OrbitType,
    ) -> SarResult<PathBuf> {
        let orbit_data =
            OrbitReader::download_orbit_file(product_id, start_time, orbit_type, None)?;

        let mut content = String::new();
        for sv in &orbit_data.state_vectors {
            content.push_str(&format!(
                "UTC={} X={} Y={} Z={} VX={} VY={} VZ={}\n",
                sv.time.format("%Y-%m-%dT%H:%M:%S%.6f"),
                sv.position[0],
                sv.position[1],
                sv.position[2],
                sv.velocity[0],
                sv.velocity[1],
                sv.velocity[2]
            ));
        }

        self.cache.cache_orbit(product_id, orbit_type, &content)
    }
}
