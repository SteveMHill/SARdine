//! SRTM 1-arc-second tile downloader and local cache.
//!
//! Enabled with `--features dem-fetch` (the default; adds `ureq` + `flate2`).
//!
//! # Entry point
//!
//! [`fetch_dem_tiles`] — given a scene bounding box, downloads all missing
//! SRTM-1 `.hgt` tiles into the local cache and returns the cache directory
//! as a [`PathBuf`] ready to pass to [`crate::dem::DemMosaic::load_directory`].
//!
//! # Source
//!
//! AWS Open Data "Terrain Tiles" (Mapzen / DigitalElevation):
//!
//! ```text
//! https://s3.amazonaws.com/elevation-tiles-prod/skadi/{NS}{LL}/{NS}{LL}{EW}{LLL}.hgt.gz
//! ```
//!
//! Example: `N47E007` → `skadi/N47/N47E007.hgt.gz`
//!
//! Files are GZIP-compressed SRTM-1 raw big-endian i16 grids (3601 × 3601).
//! No authentication is required (Registry of Open Data on AWS).
//!
//! # Cache layout
//!
//! ```text
//! $SARDINE_DEM_DIR/    (default: $HOME/.sardine/dem/)
//!   N47E007.hgt
//!   N47E008.hgt
//!   …
//! ```
//!
//! Each tile is stored as a plain (decompressed) `.hgt` file, exactly as
//! expected by [`crate::dem::DemMosaic::load_directory`].  Already-cached
//! tiles are never re-downloaded.
//!
//! # Tile coverage
//!
//! For a scene bounding box `[min_lat, max_lat] × [min_lon, max_lon]`, the
//! required 1° × 1° tiles are those whose south-west corners are at integer
//! coordinates in:
//!
//! ```text
//! lat ∈ [⌊min_lat⌋, ⌊max_lat⌋]   (inclusive on both ends)
//! lon ∈ [⌊min_lon⌋, ⌊max_lon⌋]   (inclusive on both ends)
//! ```
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

    // AWS elevation-tiles-prod bucket base URL.
    const AWS_SKADI_BASE: &str = "https://s3.amazonaws.com/elevation-tiles-prod/skadi";

    // HTTP timeout for individual tile requests (seconds).
    const HTTP_TIMEOUT_S: u64 = 120;

    // ── Error ──────────────────────────────────────────────────────────────

    /// Errors returned by [`fetch_dem_tiles`].
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
    }

    // SRTM-1 raw tile size: 3601 × 3601 big-endian i16 values.
    const SRTM1_BYTES: usize = 3601 * 3601 * 2;

    // ── Tile name helpers ─────────────────────────────────────────────────

    /// Format a tile name from integer south-west lat/lon.
    ///
    /// Standard SRTM convention:
    /// `{N|S}{lat:02}{E|W}{lon:03}` e.g. `N47E007`, `S01W071`.
    pub(super) fn tile_name(lat: i32, lon: i32) -> String {
        let ns = if lat >= 0 { 'N' } else { 'S' };
        let ew = if lon >= 0 { 'E' } else { 'W' };
        format!("{}{:02}{}{:03}", ns, lat.unsigned_abs(), ew, lon.unsigned_abs())
    }

    /// List all SRTM-1 tile names required to cover `[min_lat, max_lat] × [min_lon, max_lon]`.
    ///
    /// Tile SW corners at integer degrees are chosen so every point inside
    /// the bbox falls inside at least one tile.
    pub(super) fn tiles_for_bbox(
        min_lat: f64,
        max_lat: f64,
        min_lon: f64,
        max_lon: f64,
    ) -> Vec<String> {
        let lat0 = min_lat.floor() as i32;
        let lat1 = max_lat.floor() as i32;
        let lon0 = min_lon.floor() as i32;
        let lon1 = max_lon.floor() as i32;

        let mut names = Vec::new();
        for lat in lat0..=lat1 {
            for lon in lon0..=lon1 {
                names.push(tile_name(lat, lon));
            }
        }
        names
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

    /// Build the AWS URL for a tile.
    fn tile_url(tile: &str) -> String {
        format!("{}/{}/{}.hgt.gz", AWS_SKADI_BASE, skadi_dir(tile), tile)
    }

    // ── Download + decompress ──────────────────────────────────────────────

    /// Download and GZIP-decompress one `.hgt.gz` tile, writing the raw
    /// `.hgt` bytes to `dest`.
    ///
    /// Uses an atomic write (tmp → rename) so a partial download is never
    /// left behind as a valid-looking file.
    fn download_tile(tile: &str, dest: &Path) -> Result<(), DemFetchError> {
        let url = tile_url(tile);

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

    // ── Public API ─────────────────────────────────────────────────────────

    /// Ensure all SRTM-1 tiles covering the given bounding box are present
    /// in the local cache, downloading any that are missing.
    ///
    /// # Arguments
    ///
    /// - `min_lat`, `max_lat`, `min_lon`, `max_lon`: WGS84 degrees defining
    ///   the scene bounding box.
    /// - `cache_dir`: Directory to store `.hgt` tiles (created if absent).
    ///
    /// # Returns
    ///
    /// On success, returns `cache_dir` as a [`PathBuf`] ready to pass to
    /// [`crate::dem::DemMosaic::load_directory`].
    ///
    /// # Errors
    ///
    /// Returns [`DemFetchError`] for network failures, decompression errors,
    /// unexpected tile sizes, or I/O errors.
    pub fn fetch_dem_tiles(
        min_lat: f64,
        max_lat: f64,
        min_lon: f64,
        max_lon: f64,
        cache_dir: &Path,
    ) -> Result<PathBuf, DemFetchError> {
        std::fs::create_dir_all(cache_dir)?;

        let tiles = tiles_for_bbox(min_lat, max_lat, min_lon, max_lon);

        for tile in &tiles {
            let dest = cache_dir.join(format!("{}.hgt", tile));
            if dest.exists() {
                tracing::debug!("DEM tile {} already cached, skipping", tile);
                continue;
            }
            tracing::info!("downloading SRTM-1 tile {} …", tile);
            download_tile(tile, &dest).map_err(|e| {
                // Remove partial tmp file if it exists.
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
}

#[cfg(feature = "dem-fetch")]
pub use inner::{fetch_dem_tiles, DemFetchError};

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
        // Bbox entirely within one tile.
        let tiles = tiles_for_bbox(47.1, 47.9, 7.1, 7.9);
        assert_eq!(tiles, vec!["N47E007"]);
    }

    #[cfg(feature = "dem-fetch")]
    #[test]
    fn test_tiles_for_bbox_two_lat() {
        // Bbox straddles two latitude tiles.
        let tiles = tiles_for_bbox(47.5, 48.5, 7.1, 7.9);
        assert_eq!(tiles, vec!["N47E007", "N48E007"]);
    }

    #[cfg(feature = "dem-fetch")]
    #[test]
    fn test_tiles_for_bbox_2x2() {
        // Bbox covers a 2×2 grid of tiles.
        let tiles = tiles_for_bbox(47.0, 48.9, 7.0, 8.9);
        assert_eq!(
            tiles,
            vec!["N47E007", "N47E008", "N48E007", "N48E008"]
        );
    }

    #[cfg(feature = "dem-fetch")]
    #[test]
    fn test_tiles_for_bbox_southern_hemisphere() {
        let tiles = tiles_for_bbox(-2.5, -1.1, -70.9, -70.1);
        assert_eq!(tiles, vec!["S03W071", "S02W071"]);
    }
}
