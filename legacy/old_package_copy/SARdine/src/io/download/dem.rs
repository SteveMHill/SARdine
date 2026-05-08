#![allow(dead_code, unused_variables)]
//! DEM tile downloader

use crate::types::{BoundingBox, SarError, SarResult};
use std::path::PathBuf;
use std::sync::Arc;

use super::cache::CacheManager;
use super::manager::{DownloadConfig, DownloadResult};
use super::progress::ProgressCallback;
use super::providers::ProviderRegistry;
use super::utils::{download_from_url, extract_from_zip, is_zip_content};

/// DEM downloader
pub struct DEMDownloader {
    cache: Arc<CacheManager>,
    providers: Arc<ProviderRegistry>,
    config: DownloadConfig,
}

impl DEMDownloader {
    pub fn new(
        cache: &Arc<CacheManager>,
        providers: &Arc<ProviderRegistry>,
        config: &DownloadConfig,
    ) -> SarResult<Self> {
        Ok(Self {
            cache: cache.clone(),
            providers: providers.clone(),
            config: config.clone(),
        })
    }

    /// Download DEM tiles for a bounding box
    pub fn download_tiles(
        &self,
        bbox: &BoundingBox,
        _progress: Option<ProgressCallback>,
    ) -> SarResult<Vec<DownloadResult>> {
        // Calculate required SRTM tiles
        let tiles = self.calculate_srtm_tiles(bbox);
        let mut downloaded_files = Vec::new();

        let dem_dir = self.cache.dem_dir();
        std::fs::create_dir_all(&dem_dir).map_err(SarError::Io)?;

        for tile in tiles {
            // Check if file already exists
            let mut existing_path = None;
            for ext in &["hgt", "tif", "tiff"] {
                let test_path = dem_dir.join(format!("{}.{}", tile, ext));
                if test_path.exists() {
                    log::info!(
                        "DEM tile {} already exists, skipping download",
                        test_path.display()
                    );
                    existing_path = Some(test_path);
                    break;
                }
            }

            if let Some(path) = existing_path {
                if let Ok(entry) = CacheManager::get_cache_entry(&path) {
                    downloaded_files.push(DownloadResult {
                        path,
                        size: entry.size,
                        source: "cache".to_string(),
                    });
                }
                continue;
            }

            // Try downloading from multiple sources
            if let Some(path) = self.try_download_tile(&tile, &dem_dir)? {
                downloaded_files.push(DownloadResult {
                    path,
                    size: 0, // Will be set properly
                    source: "download".to_string(),
                });
            }
        }

        if downloaded_files.is_empty() {
            return Err(SarError::Processing(
                "Failed to download any DEM tiles from available sources.".to_string(),
            ));
        }

        Ok(downloaded_files)
    }

    fn try_download_tile(&self, tile: &str, output_dir: &PathBuf) -> SarResult<Option<PathBuf>> {
        // Multiple DEM data sources (in order of preference)
        let mut sources = Vec::new();

        // AWS USGS SRTM 1 Arcsecond Tiles (most reliable, no auth required)
        if let Some(dir) = self.get_skadi_directory(tile) {
            sources.push(format!(
                "https://s3.amazonaws.com/elevation-tiles-prod/skadi/{}/{}.hgt.gz",
                dir, tile
            ));
        }

        // Add backup sources
        sources.extend(vec![format!(
            "https://e4ftl01.cr.usgs.gov/MEASURES/SRTMGL1.003/2000.02.11/{}.SRTMGL1.hgt.zip",
            tile
        )]);

        for url in &sources {
            log::info!("Attempting to download DEM tile from: {}", url);
            match self.download_and_extract_tile(url, output_dir, tile) {
                Ok(path) => {
                    log::info!("Successfully downloaded DEM tile: {}", path.display());
                    return Ok(Some(path));
                }
                Err(e) => {
                    log::debug!("Download failed: {}, trying next source", e);
                }
            }
        }

        Ok(None)
    }

    fn download_and_extract_tile(
        &self,
        url: &str,
        output_dir: &PathBuf,
        tile: &str,
    ) -> SarResult<PathBuf> {
        let bytes = download_from_url(url, None, self.config.timeout_seconds)?;

        let output_path = if url.ends_with(".gz") {
            // Handle gzipped files
            use flate2::read::GzDecoder;
            use std::io::Read;
            let mut decoder = GzDecoder::new(&bytes[..]);
            let mut decompressed = Vec::new();
            decoder
                .read_to_end(&mut decompressed)
                .map_err(SarError::Io)?;
            let path = output_dir.join(format!("{}.hgt", tile));
            std::fs::write(&path, &decompressed).map_err(SarError::Io)?;
            path
        } else if url.ends_with(".zip") || is_zip_content(&bytes) {
            // Handle ZIP files
            let hgt_bytes = extract_from_zip(&bytes, ".hgt")?;
            let path = output_dir.join(format!("{}.hgt", tile));
            std::fs::write(&path, &hgt_bytes).map_err(SarError::Io)?;
            path
        } else {
            // Assume raw HGT file
            let path = output_dir.join(format!("{}.hgt", tile));
            std::fs::write(&path, &bytes).map_err(SarError::Io)?;
            path
        };

        Ok(output_path)
    }

    fn calculate_srtm_tiles(&self, bbox: &BoundingBox) -> Vec<String> {
        let mut tiles = Vec::new();

        // SRTM tiles are 1x1 degree
        let min_lat = bbox.min_lat.floor() as i32;
        let max_lat = bbox.max_lat.floor() as i32;
        let min_lon = bbox.min_lon.floor() as i32;
        let max_lon = bbox.max_lon.floor() as i32;

        for lat in min_lat..=max_lat {
            for lon in min_lon..=max_lon {
                let lat_prefix = if lat >= 0 { "N" } else { "S" };
                let lon_prefix = if lon >= 0 { "E" } else { "W" };

                let tile_name = format!(
                    "{}{:02}{}{:03}",
                    lat_prefix,
                    lat.abs(),
                    lon_prefix,
                    lon.abs()
                );

                tiles.push(tile_name);
            }
        }

        tiles
    }

    fn get_skadi_directory(&self, tile: &str) -> Option<String> {
        if tile.len() >= 3 {
            let lat_part = &tile[0..3];
            Some(lat_part.to_string())
        } else {
            None
        }
    }
}
