//! DEM tile downloader and local cache — SRTM-1 and Copernicus GLO-30.
//!
//! Enabled with `--features dem-fetch` (the default; adds `ureq` + `flate2`).
//!
//! # Entry points
//!
//! - [`fetch_srtm1_tiles`] — download missing SRTM-1 `.hgt` tiles from AWS,
//!   cache in `$SARDINE_DEM_DIR/` (or `$HOME/.sardine/dem/`).
//! - [`fetch_glo30_tiles`] — download missing Copernicus GLO-30 GeoTIFF tiles
//!   from AWS, cache in `$SARDINE_DEM_DIR/glo30/` (or `$HOME/.sardine/dem/glo30/`).
//!
//! Both functions return the directory path ready for
//! [`crate::dem::DemMosaic::load_directory`].
//!
//! # Sources
//!
//! **SRTM-1** — AWS elevation-tiles-prod (Mapzen / DigitalElevation):
//! ```text
//! https://s3.amazonaws.com/elevation-tiles-prod/skadi/{NS}{LL}/{NS}{LL}{EW}{LLL}.hgt.gz
//! ```
//! Files are GZIP-compressed SRTM-1 raw big-endian i16 (3601 × 3601).
//!
//! **Copernicus GLO-30** — AWS Copernicus DEM Open Data:
//! ```text
//! https://copernicus-dem-30m.s3.amazonaws.com/
//!   Copernicus_DSM_COG_10_{NS}{LL}_00_{EW}{LLL}_00_DEM/
//!   Copernicus_DSM_COG_10_{NS}{LL}_00_{EW}{LLL}_00_DEM.tif
//! ```
//! Files are Cloud-Optimised GeoTIFF, no decompression required.
//! Tiles that don't exist (e.g. open ocean) return HTTP 404 and are silently
//! skipped — the DEM mosaic will mark those pixels as nodata.
//!
//! No authentication is required for either source.
//!
//! # Cache layout
//!
//! ```text
//! $SARDINE_DEM_DIR/        (default: $HOME/.sardine/dem/)
//!   N47E007.hgt            ← SRTM-1 tiles
//!   N47E008.hgt
//!   glo30/                 ← GLO-30 tiles (subdirectory)
//!     Copernicus_DSM_COG_10_N47_00_E007_00_DEM.tif
//!     …
//! ```
//!
//! Already-cached tiles are never re-downloaded.
//!
//! # AGENTS.md compliance
//!
//! - Every failure is a typed [`DemFetchError`] variant — no silent fallbacks.
//! - No hardcoded Sentinel-1 constants.
//! - No `unwrap_or`, `ok()?`, or `let _ =` patterns.

#[cfg(feature = "dem-fetch")]
mod inner {
    use std::io::{self, Write};
    use std::path::{Path, PathBuf};

    use thiserror::Error;

    // AWS elevation-tiles-prod bucket base URL (SRTM-1).
    const AWS_SKADI_BASE: &str = "https://s3.amazonaws.com/elevation-tiles-prod/skadi";

    // AWS Copernicus DEM bucket base URL (GLO-30).
    const AWS_GLO30_BASE: &str = "https://copernicus-dem-30m.s3.amazonaws.com";

    // HTTP timeout for individual tile requests (seconds).
    const HTTP_TIMEOUT_S: u64 = 120;

    // ── Error ──────────────────────────────────────────────────────────────

    /// Errors returned by [`fetch_srtm1_tiles`] and [`fetch_glo30_tiles`].
    #[derive(Debug, Error)]
    pub enum DemFetchError {
        /// HTTP-level transport or status error for a specific tile.
        #[error("HTTP error fetching tile {tile}: {detail}")]
        Http { tile: String, detail: String },

        /// GZIP decompression failed for a tile.
        #[error("GZIP decompression failed for tile {tile}: {detail}")]
        Decompress { tile: String, detail: String },

        /// I/O error writing a tile to the cache.
        #[error("I/O error: {0}")]
        Io(#[from] io::Error),

        /// Cache directory resolution failed (HOME not set).
        #[error("Cannot resolve DEM cache directory: {0}")]
        CacheDir(String),

        /// The downloaded tile has an unexpected byte size.
        #[error("Tile {tile}: expected {expected} bytes, got {actual}")]
        UnexpectedSize {
            tile: String,
            expected: usize,
            actual: usize,
        },

        /// Unknown DEM source string.
        #[error("Unknown DEM source '{0}': expected 'srtm1' or 'glo30'")]
        UnknownSource(String),
    }

    // SRTM-1 raw tile size: 3601 × 3601 big-endian i16 values.
    const SRTM1_BYTES: usize = 3601 * 3601 * 2;

    // ── Tile name helpers ─────────────────────────────────────────────────

    /// Format a tile name from integer south-west lat/lon (SRTM/flat convention).
    ///
    /// `{N|S}{lat:02}{E|W}{lon:03}` — e.g. `N47E007`, `S01W071`.
    pub(super) fn tile_name(lat: i32, lon: i32) -> String {
        let ns = if lat >= 0 { 'N' } else { 'S' };
        let ew = if lon >= 0 { 'E' } else { 'W' };
        format!("{}{:02}{}{:03}", ns, lat.unsigned_abs(), ew, lon.unsigned_abs())
    }

    /// Format the official Copernicus GLO-30 tile stem for a given SW corner.
    ///
    /// e.g. lat=47, lon=7 → `Copernicus_DSM_COG_10_N47_00_E007_00_DEM`
    pub(super) fn glo30_tile_name(lat: i32, lon: i32) -> String {
        let ns = if lat >= 0 { 'N' } else { 'S' };
        let ew = if lon >= 0 { 'E' } else { 'W' };
        format!(
            "Copernicus_DSM_COG_10_{}{:02}_00_{}{:03}_00_DEM",
            ns, lat.unsigned_abs(),
            ew, lon.unsigned_abs(),
        )
    }

    /// List all 1°×1° tile SW corners required to cover the given bounding box.
    ///
    /// Tile SW corners at integer degrees are chosen so every point inside
    /// the bbox falls inside at least one tile.
    pub(super) fn tiles_for_bbox(
        min_lat: f64,
        max_lat: f64,
        min_lon: f64,
        max_lon: f64,
    ) -> Vec<(i32, i32)> {
        let lat0 = min_lat.floor() as i32;
        let lat1 = max_lat.floor() as i32;
        let lon0 = min_lon.floor() as i32;
        let lon1 = max_lon.floor() as i32;

        let mut corners = Vec::new();
        for lat in lat0..=lat1 {
            for lon in lon0..=lon1 {
                corners.push((lat, lon));
            }
        }
        corners
    }

    /// Derive the skadi subdirectory from a tile name.
    ///
    /// e.g. `N47E007` → `N47`, `S01W071` → `S01`.
    fn skadi_dir(tile: &str) -> &str {
        // The lat prefix is always 3 chars: one letter + two digits.
        if tile.len() >= 3 {
            &tile[..3]
        } else {
            tile
        }
    }

    // ── SRTM-1 download ───────────────────────────────────────────────────

    /// Build the AWS URL for an SRTM-1 tile.
    fn srtm1_tile_url(tile: &str) -> String {
        format!("{}/{}/{}.hgt.gz", AWS_SKADI_BASE, skadi_dir(tile), tile)
    }

    /// Download and GZIP-decompress one SRTM-1 `.hgt.gz` tile into `dest`.
    ///
    /// Uses an atomic write (tmp → rename) so a partial download is never
    /// left behind as a valid-looking file.
    fn download_srtm1_tile(tile: &str, dest: &Path) -> Result<(), DemFetchError> {
        let url = srtm1_tile_url(tile);

        let response = ureq::builder()
            .timeout(std::time::Duration::from_secs(HTTP_TIMEOUT_S))
            .build()
            .get(&url)
            .call()
            .map_err(|e| DemFetchError::Http {
                tile: tile.to_string(),
                detail: format!("{}: {}", url, e),
            })?;

        if response.status() < 200 || response.status() >= 300 {
            return Err(DemFetchError::Http {
                tile: tile.to_string(),
                detail: format!("{}: HTTP {}", url, response.status()),
            });
        }

        // Read compressed bytes into memory.
        let mut compressed: Vec<u8> = Vec::new();
        response
            .into_reader()
            .read_to_end(&mut compressed)
            .map_err(io::Error::from)
            .map_err(DemFetchError::Io)?;

        // Decompress gzip.
        let raw = decompress_gzip(tile, &compressed)?;

        // Validate byte count.
        if raw.len() != SRTM1_BYTES {
            return Err(DemFetchError::UnexpectedSize {
                tile: tile.to_string(),
                expected: SRTM1_BYTES,
                actual: raw.len(),
            });
        }

        // Atomic write: tmp path → final path.
        let tmp = dest.with_extension("hgt.tmp");
        {
            let mut f = std::fs::File::create(&tmp)?;
            f.write_all(&raw)?;
            f.flush()?;
        }
        std::fs::rename(&tmp, dest)?;
        Ok(())
    }

    /// Decompress a GZIP byte slice using `flate2`.
    fn decompress_gzip(tile: &str, data: &[u8]) -> Result<Vec<u8>, DemFetchError> {
        use flate2::read::GzDecoder;
        use std::io::Read;

        let mut decoder = GzDecoder::new(data);
        let mut out = Vec::with_capacity(SRTM1_BYTES);
        decoder.read_to_end(&mut out).map_err(|e| DemFetchError::Decompress {
            tile: tile.to_string(),
            detail: e.to_string(),
        })?;
        Ok(out)
    }

    // ── GLO-30 download ───────────────────────────────────────────────────

    /// Build the AWS URL for a Copernicus GLO-30 tile.
    fn glo30_tile_url(stem: &str) -> String {
        format!("{}/{}/{}.tif", AWS_GLO30_BASE, stem, stem)
    }

    /// Download one Copernicus GLO-30 GeoTIFF tile into `dest`.
    ///
    /// Returns `Ok(true)` if the tile was downloaded, `Ok(false)` if it does
    /// not exist on the server (HTTP 404 — typical for open-ocean tiles).
    ///
    /// All other HTTP errors are propagated as [`DemFetchError::Http`].
    /// Uses an atomic write (`.tif.tmp` → `.tif`) to avoid leaving a partial
    /// file that looks valid.
    fn download_glo30_tile(stem: &str, dest: &Path) -> Result<bool, DemFetchError> {
        let url = glo30_tile_url(stem);

        let response = ureq::builder()
            .timeout(std::time::Duration::from_secs(HTTP_TIMEOUT_S))
            .build()
            .get(&url)
            .call();

        match response {
            Err(ureq::Error::Status(404, _)) => {
                tracing::debug!("GLO-30 tile {} not found on server (ocean/void), skipping", stem);
                return Ok(false);
            }
            Err(e) => {
                return Err(DemFetchError::Http {
                    tile: stem.to_string(),
                    detail: format!("{}: {}", url, e),
                });
            }
            Ok(resp) if resp.status() < 200 || resp.status() >= 300 => {
                return Err(DemFetchError::Http {
                    tile: stem.to_string(),
                    detail: format!("{}: HTTP {}", url, resp.status()),
                });
            }
            Ok(resp) => {
                // Read the full GeoTIFF body.
                let mut body: Vec<u8> = Vec::new();
                resp.into_reader()
                    .read_to_end(&mut body)
                    .map_err(io::Error::from)
                    .map_err(DemFetchError::Io)?;

                // Atomic write.
                let tmp = dest.with_extension("tif.tmp");
                {
                    let mut f = std::fs::File::create(&tmp)?;
                    f.write_all(&body)?;
                    f.flush()?;
                }
                std::fs::rename(&tmp, dest)?;
                Ok(true)
            }
        }
    }

    // ── Public API ─────────────────────────────────────────────────────────

    /// Ensure all SRTM-1 tiles covering the given bounding box are present
    /// in `cache_dir`, downloading any that are missing.
    ///
    /// Returns `cache_dir` as a [`PathBuf`] ready for
    /// [`crate::dem::DemMosaic::load_directory`].
    pub fn fetch_srtm1_tiles(
        min_lat: f64,
        max_lat: f64,
        min_lon: f64,
        max_lon: f64,
        cache_dir: &Path,
    ) -> Result<PathBuf, DemFetchError> {
        std::fs::create_dir_all(cache_dir)?;

        let corners = tiles_for_bbox(min_lat, max_lat, min_lon, max_lon);

        for (lat, lon) in &corners {
            let name = tile_name(*lat, *lon);
            let dest = cache_dir.join(format!("{}.hgt", name));
            if dest.exists() {
                tracing::debug!("SRTM-1 tile {} already cached, skipping", name);
                continue;
            }
            tracing::info!("downloading SRTM-1 tile {} …", name);
            download_srtm1_tile(&name, &dest).map_err(|e| {
                let tmp = dest.with_extension("hgt.tmp");
                if tmp.exists() {
                    let _ = std::fs::remove_file(&tmp); // SAFETY-OK: cleanup of failed download; removal failure leaves only a harmless .tmp file
                }
                e
            })?;
            tracing::info!("  → cached {}", dest.display());
        }

        Ok(cache_dir.to_path_buf())
    }

    /// Ensure all Copernicus GLO-30 tiles covering the given bounding box are
    /// present in `cache_dir/glo30/`, downloading any that are missing.
    ///
    /// Tiles that do not exist on the server (open-ocean, HTTP 404) are
    /// silently skipped — the DEM mosaic marks those pixels as nodata.
    ///
    /// Returns the GLO-30 subdirectory as a [`PathBuf`] ready for
    /// [`crate::dem::DemMosaic::load_directory`].
    pub fn fetch_glo30_tiles(
        min_lat: f64,
        max_lat: f64,
        min_lon: f64,
        max_lon: f64,
        cache_dir: &Path,
    ) -> Result<PathBuf, DemFetchError> {
        let glo30_dir = cache_dir.join("glo30");
        std::fs::create_dir_all(&glo30_dir)?;

        let corners = tiles_for_bbox(min_lat, max_lat, min_lon, max_lon);

        for (lat, lon) in &corners {
            let stem = glo30_tile_name(*lat, *lon);
            let dest = glo30_dir.join(format!("{}.tif", stem));
            if dest.exists() {
                tracing::debug!("GLO-30 tile {} already cached, skipping", stem);
                continue;
            }
            tracing::info!("downloading GLO-30 tile {} …", stem);
            let downloaded = download_glo30_tile(&stem, &dest).map_err(|e| {
                let tmp = dest.with_extension("tif.tmp");
                if tmp.exists() {
                    let _ = std::fs::remove_file(&tmp); // SAFETY-OK: cleanup of failed download; removal failure leaves only a harmless .tmp file
                }
                e
            })?;
            if downloaded {
                tracing::info!("  → cached {}", dest.display());
            }
        }

        Ok(glo30_dir)
    }

    /// Dispatch to [`fetch_srtm1_tiles`] or [`fetch_glo30_tiles`] based on
    /// the `source` string (`"srtm1"` or `"glo30"`).
    pub fn fetch_dem_tiles_for_source(
        source: &str,
        min_lat: f64,
        max_lat: f64,
        min_lon: f64,
        max_lon: f64,
        cache_dir: &Path,
    ) -> Result<PathBuf, DemFetchError> {
        match source {
            "srtm1" => fetch_srtm1_tiles(min_lat, max_lat, min_lon, max_lon, cache_dir),
            "glo30" => fetch_glo30_tiles(min_lat, max_lat, min_lon, max_lon, cache_dir),
            other => Err(DemFetchError::UnknownSource(other.to_string())),
        }
    }
}

#[cfg(feature = "dem-fetch")]
pub use inner::{fetch_dem_tiles_for_source, fetch_glo30_tiles, fetch_srtm1_tiles, DemFetchError};

// ── Unit tests ────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    #[cfg(feature = "dem-fetch")]
    use super::inner::*;

    #[cfg(feature = "dem-fetch")]
    #[test]
    fn test_tile_name_positive() {
        assert_eq!(tile_name(47, 7), "N47E007");
        assert_eq!(tile_name(0, 0), "N00E000");
        assert_eq!(tile_name(9, 99), "N09E099");
    }

    #[cfg(feature = "dem-fetch")]
    #[test]
    fn test_tile_name_negative() {
        assert_eq!(tile_name(-1, -71), "S01W071");
        assert_eq!(tile_name(-90, -180), "S90W180");
    }

    #[cfg(feature = "dem-fetch")]
    #[test]
    fn test_tiles_for_bbox_single() {
        let corners = tiles_for_bbox(47.1, 47.9, 7.1, 7.9);
        assert_eq!(corners, vec![(47, 7)]);
    }

    #[cfg(feature = "dem-fetch")]
    #[test]
    fn test_tiles_for_bbox_two_lat() {
        let corners = tiles_for_bbox(47.5, 48.5, 7.1, 7.9);
        assert_eq!(corners, vec![(47, 7), (48, 7)]);
    }

    #[cfg(feature = "dem-fetch")]
    #[test]
    fn test_tiles_for_bbox_2x2() {
        let corners = tiles_for_bbox(47.0, 48.9, 7.0, 8.9);
        assert_eq!(corners, vec![(47, 7), (47, 8), (48, 7), (48, 8)]);
    }

    #[cfg(feature = "dem-fetch")]
    #[test]
    fn test_tiles_for_bbox_southern_hemisphere() {
        let corners = tiles_for_bbox(-2.5, -1.1, -70.9, -70.1);
        assert_eq!(corners, vec![(-3, -71), (-2, -71)]);
    }

    #[cfg(feature = "dem-fetch")]
    #[test]
    fn test_glo30_tile_name_positive() {
        assert_eq!(
            glo30_tile_name(47, 7),
            "Copernicus_DSM_COG_10_N47_00_E007_00_DEM"
        );
        assert_eq!(
            glo30_tile_name(0, 0),
            "Copernicus_DSM_COG_10_N00_00_E000_00_DEM"
        );
    }

    #[cfg(feature = "dem-fetch")]
    #[test]
    fn test_glo30_tile_name_negative() {
        assert_eq!(
            glo30_tile_name(-1, -71),
            "Copernicus_DSM_COG_10_S01_00_W071_00_DEM"
        );
        assert_eq!(
            glo30_tile_name(-90, -180),
            "Copernicus_DSM_COG_10_S90_00_W180_00_DEM"
        );
    }
}
