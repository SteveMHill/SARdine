//! EGM96 and EGM2008 geoid grid downloader and local cache.
//!
//! Enabled with `--features geoid-fetch` (adds the `ureq` HTTP dependency).
//!
//! # Entry points
//!
//! * [`fetch_egm96`] — ensures the EGM96 geoid grid is present in the local
//!   cache directory, downloading it if necessary, then loads and returns an
//!   [`Egm96Grid`].
//! * [`fetch_egm2008`] — same, for [`Egm2008Grid`] (EGM2008 at 2.5°; use
//!   with Copernicus GLO-30 DEMs).
//!
//! # Source
//!
//! The NGA EGM96 15-arc-minute geoid grid is distributed by the PROJ project
//! as a Cloud Optimised GeoTIFF on the PROJ datum-grid CDN:
//!
//! ```text
//! https://cdn.proj.org/us_nga_egm96_15.tif
//! ```
//!
//! The GeoTIFF is georeferenced via `ModelPixelScaleTag` (33550) and
//! `ModelTiepointTag` (33922), with a single f32 band of undulation values
//! in metres.  It is parsed with the `tiff` crate (same approach as the
//! GLO-30 DEM tiles in `dem.rs`) and downsampled to the 73×144 (2.5°) grid
//! stored by [`Egm96Grid`].
//!
//! # Cache layout
//!
//! ```text
//! $SARDINE_GEOID_DIR/   (default: $HOME/.sardine/geoid/)
//!   us_nga_egm96_15.tif  raw PROJ CDN download (kept for re-use)
//!   egm96_2p5deg.bin     converted 73×144 f32 LE compact binary (fast reload)
//! ```
//!
//! If `egm96_2p5deg.bin` already exists it is loaded directly without
//! re-downloading.
//!
//! # AGENTS.md compliance
//!
//! - Every failure is a typed `GeoidFetchError` variant — no silent fallbacks.
//! - No hardcoded Sentinel-1 constants.
//! - No `unwrap_or`, `ok()?`, or `let _ =` patterns.

#[cfg(feature = "geoid-fetch")]
mod inner {
    use std::io::{self, Write};
    use std::path::{Path, PathBuf};

    use thiserror::Error;

    use crate::geoid::{Egm96Grid, Egm2008Grid, GeoidError};

    // ── Constants ──────────────────────────────────────────────────────────

    // PROJ CDN EGM96 15-arc-minute Cloud Optimised GeoTIFF.
    // Single f32 band of undulation values in metres, georeferenced via
    // ModelPixelScaleTag (34264) and ModelTiepointTag (33922).
    const COG_URL: &str = "https://cdn.proj.org/us_nga_egm96_15.tif";

    const COG_FILENAME: &str = "us_nga_egm96_15.tif";
    const BIN_FILENAME: &str = "egm96_2p5deg.bin";

    // Minimum source-grid dimensions for a 15-arc-minute (0.25°) global grid.
    // Actual PROJ CDN file is 720 rows × 1440 cols.
    const COG_MIN_NLAT: usize = 720;
    const COG_MIN_NLON: usize = 1440;

    // ── EGM2008 constants ───────────────────────────────────────────────────

    // PROJ CDN EGM2008 2.5-arc-minute Cloud Optimised GeoTIFF.
    const COG_URL_EGM2008: &str = "https://cdn.proj.org/us_nga_egm08_25.tif";
    const COG_FILENAME_EGM2008: &str = "us_nga_egm08_25.tif";
    const BIN_FILENAME_EGM2008: &str = "egm2008_2p5deg.bin";

    // Minimum source-grid dimensions for a 2.5-arc-minute (0.04167°) global grid.
    // Actual PROJ CDN file is 4320 rows × 8640 cols.
    const COG_MIN_NLAT_EGM2008: usize = 4320;
    const COG_MIN_NLON_EGM2008: usize = 8640;

    // 2.5° target grid used by Egm96Grid.
    const EGM_NLAT: usize = 73; // (90 - -90) / 2.5 + 1
    const EGM_NLON: usize = 144; // 360 / 2.5

    // HTTP timeout.
    const HTTP_TIMEOUT_S: u64 = 120;

    // ── Error ──────────────────────────────────────────────────────────────

    /// Errors returned by [`fetch_egm96`] and [`fetch_egm2008`].
    #[derive(Debug, Error)]
    pub enum GeoidFetchError {
        /// HTTP-level transport or status error.
        #[error("HTTP error fetching EGM96 grid: {0}")]
        Http(String),

        /// I/O error writing the cache file.
        #[error("I/O error: {0}")]
        Io(#[from] io::Error),

        /// The downloaded GTX file does not have the expected structure.
        #[error("EGM96 GTX file is malformed: {detail}")]
        MalformedGtx { detail: String },

        /// Grid loading failed after successful download/cache hit.
        #[error("EGM96 grid load error: {0}")]
        GridLoad(#[from] GeoidError),
    }

    // ── Public entry point ─────────────────────────────────────────────────

    /// Ensure the EGM96 2.5° grid is cached locally and load it.
    ///
    /// # Cache directory
    ///
    /// 1. Checks `$SARDINE_GEOID_DIR` environment variable.
    /// 2. Falls back to `$HOME/.sardine/geoid/`.
    ///
    /// If neither resolves, returns `GeoidFetchError::Io`.
    ///
    /// # Caching behaviour
    ///
    /// - If `egm96_2p5deg.bin` exists in the cache directory, load it
    ///   directly (fast path, ~42 KiB).
    /// - If not, download `egm96_15.gtx` (~47 MiB), downsample to 73×144,
    ///   write `egm96_2p5deg.bin`, then load that.
    ///
    /// The GTX file is kept in the cache so subsequent calls with a missing
    /// `.bin` (e.g. after a version change) can convert without re-downloading.
    pub fn fetch_egm96() -> Result<Egm96Grid, GeoidFetchError> {
        let cache_dir = resolve_cache_dir()?;
        std::fs::create_dir_all(&cache_dir)?;

        let bin_path = cache_dir.join(BIN_FILENAME);

        if bin_path.exists() {
            return Ok(Egm96Grid::load_binary(&bin_path)?);
        }

        // Need to download and convert.
        let cog_path = cache_dir.join(COG_FILENAME);
        if !cog_path.exists() {
            tracing::info!("Downloading EGM96 15-arc-minute grid from PROJ CDN …");
            download_to_file(COG_URL, &cog_path)?;
            tracing::info!("Download complete: {}", cog_path.display());
        } else {
            tracing::info!("Using cached GeoTIFF: {}", cog_path.display());
        }

        tracing::info!("Converting GeoTIFF → 2.5° compact binary …");
        let grid = cog_to_egm96_grid(&cog_path)?;
        save_binary(&grid, &bin_path)?;
        tracing::info!("Saved: {}", bin_path.display());

        Ok(grid)
    }

    /// Ensure the EGM2008 2.5° grid is cached locally and load it.
    ///
    /// Uses the same cache directory as [`fetch_egm96`]
    /// (`$SARDINE_GEOID_DIR` or `$HOME/.sardine/geoid/`).
    ///
    /// The PROJ CDN source file (`us_nga_egm08_25.tif`, ~120 MiB) is
    /// downloaded once, downsampled to a 73×144 compact binary
    /// (`egm2008_2p5deg.bin`, ~41 KiB), and cached for fast subsequent loads.
    ///
    /// Use with Copernicus GLO-30 DEM tiles, whose vertical datum is EGM2008.
    pub fn fetch_egm2008() -> Result<Egm2008Grid, GeoidFetchError> {
        let cache_dir = resolve_cache_dir()?;
        std::fs::create_dir_all(&cache_dir)?;

        let bin_path = cache_dir.join(BIN_FILENAME_EGM2008);

        if bin_path.exists() {
            return Ok(Egm2008Grid::load_binary(&bin_path)?);
        }

        let cog_path = cache_dir.join(COG_FILENAME_EGM2008);
        if !cog_path.exists() {
            tracing::info!("Downloading EGM2008 2.5-arc-minute grid from PROJ CDN …");
            download_to_file(COG_URL_EGM2008, &cog_path)?;
            tracing::info!("Download complete: {}", cog_path.display());
        } else {
            tracing::info!("Using cached GeoTIFF: {}", cog_path.display());
        }

        tracing::info!("Converting EGM2008 GeoTIFF → 2.5° compact binary …");
        let grid = cog_to_egm2008_grid(&cog_path)?;
        save_binary_egm2008(&grid, &bin_path)?;
        tracing::info!("Saved: {}", bin_path.display());

        Ok(grid)
    }

    // ── Cache directory resolution ─────────────────────────────────────────

    fn resolve_cache_dir() -> Result<PathBuf, GeoidFetchError> {
        if let Ok(dir) = std::env::var("SARDINE_GEOID_DIR") {
            return Ok(PathBuf::from(dir));
        }
        let home = std::env::var("HOME").map_err(|_| {
            GeoidFetchError::Io(io::Error::new(
                io::ErrorKind::NotFound,
                "HOME environment variable not set and SARDINE_GEOID_DIR not set",
            ))
        })?;
        Ok(PathBuf::from(home).join(".sardine").join("geoid"))
    }

    // ── HTTP download ──────────────────────────────────────────────────────

    fn download_to_file(url: &str, dest: &Path) -> Result<(), GeoidFetchError> {
        let response = ureq::builder()
            .timeout(std::time::Duration::from_secs(HTTP_TIMEOUT_S))
            .build()
            .get(url)
            .call()
            .map_err(|e| GeoidFetchError::Http(format!("{url}: {e}")))?;

        let status = response.status();
        if status < 200 || status >= 300 {
            return Err(GeoidFetchError::Http(format!("{url}: HTTP {status}")));
        }

        // Write via temp file then atomic rename.
        let tmp = dest.with_extension("tmp");
        {
            let mut file = std::fs::File::create(&tmp)?;
            let mut reader = response.into_reader();
            io::copy(&mut reader, &mut file)?;
            file.flush()?;
        }
        std::fs::rename(&tmp, dest)?;
        Ok(())
    }

    // ── GeoTIFF → Egm96Grid conversion ────────────────────────────────────

    /// Parse the PROJ CDN Cloud Optimised GeoTIFF and downsample to the
    /// 73×144 (2.5°) grid.
    ///
    /// Reads georeferencing from `ModelPixelScaleTag` (33550) and
    /// `ModelTiepointTag` (33922), then nearest-neighbour resamples to the
    /// 2.5° output grid.  Uses the same `tiff`-crate approach as the GLO-30
    /// DEM loader in `dem.rs`.
    fn cog_to_egm96_grid(path: &Path) -> Result<Egm96Grid, GeoidFetchError> {
        use std::io::BufReader;
        use tiff::decoder::{Decoder, DecodingResult};
        use tiff::tags::Tag;

        const TAG_MODEL_PIXEL_SCALE: u16 = 33550;
        const TAG_MODEL_TIEPOINT: u16 = 33922;

        let file = std::fs::File::open(path)?;
        let mut dec = Decoder::new(BufReader::new(file)).map_err(|e| GeoidFetchError::MalformedGtx {
            detail: e.to_string(),
        })?;

        let (width_u32, height_u32) = dec.dimensions().map_err(|e| GeoidFetchError::MalformedGtx {
            detail: e.to_string(),
        })?;
        let n_lon = width_u32 as usize;
        let n_lat = height_u32 as usize;

        // Require minimum dimensions for a 15-arc-minute global grid.
        if n_lat < COG_MIN_NLAT || n_lon < COG_MIN_NLON {
            return Err(GeoidFetchError::MalformedGtx {
                detail: format!(
                    "EGM96 grid {n_lat}×{n_lon} is too small; expected ≥ {COG_MIN_NLAT}×{COG_MIN_NLON}"
                ),
            });
        }

        // ModelPixelScaleTag: [scale_lon, scale_lat, 0]  (all positive)
        let scale = dec
            .get_tag_f64_vec(Tag::Unknown(TAG_MODEL_PIXEL_SCALE))
            .map_err(|e| GeoidFetchError::MalformedGtx {
                detail: format!("ModelPixelScaleTag: {e}"),
            })?;
        if scale.len() < 2 {
            return Err(GeoidFetchError::MalformedGtx {
                detail: format!(
                    "ModelPixelScaleTag has {} elements, expected ≥ 2",
                    scale.len()
                ),
            });
        }
        let d_lon = scale[0]; // degrees per column
        let d_lat = scale[1]; // degrees per row (positive; rows go N→S)

        // ModelTiepointTag: [raster_i, raster_j, 0, world_lon, world_lat, 0]
        let tp = dec
            .get_tag_f64_vec(Tag::Unknown(TAG_MODEL_TIEPOINT))
            .map_err(|e| GeoidFetchError::MalformedGtx {
                detail: format!("ModelTiepointTag: {e}"),
            })?;
        if tp.len() < 6 {
            return Err(GeoidFetchError::MalformedGtx {
                detail: format!("ModelTiepointTag has {} elements, expected ≥ 6", tp.len()),
            });
        }
        // Centre of pixel (0, 0): add 0.5 pixel from the corner-referenced tiepoint.
        let lon0 = tp[3] + (0.5 - tp[0]) * d_lon;
        let lat0 = tp[4] - (0.5 - tp[1]) * d_lat;

        // Read pixels.
        let image = dec.read_image().map_err(|e| GeoidFetchError::MalformedGtx {
            detail: e.to_string(),
        })?;
        let pixels: Vec<f32> = match image {
            DecodingResult::F32(v) => v,
            DecodingResult::I16(v) => v.iter().map(|&x| x as f32).collect(),
            DecodingResult::U16(v) => v.iter().map(|&x| x as f32).collect(),
            other => {
                return Err(GeoidFetchError::MalformedGtx {
                    detail: format!("unsupported pixel type: {other:?}"),
                });
            }
        };

        // Nearest-neighbour downsample to the 73×144 (2.5°) target grid.
        //
        // Egm96Grid row r → lat = 90° − r × 2.5°
        // Egm96Grid col c → lon = c × 2.5°  (0° to 357.5°)
        //
        // Source raster: row 0 is northernmost, rows go south.
        //   row_f = (lat0 − lat) / d_lat
        //   col_f = (lon  − lon0) / d_lon
        let mut arr = [0f32; EGM_NLAT * EGM_NLON];
        for r in 0..EGM_NLAT {
            let lat = 90.0 - r as f64 * 2.5;
            let row_f = (lat0 - lat) / d_lat;
            let row_src = (row_f.round() as isize).clamp(0, (n_lat - 1) as isize) as usize;
            for c in 0..EGM_NLON {
                let lon = c as f64 * 2.5;
                let col_f = (lon - lon0) / d_lon;
                let col_src = (col_f.round() as isize).clamp(0, (n_lon - 1) as isize) as usize;
                arr[r * EGM_NLON + c] = pixels[row_src * n_lon + col_src];
            }
        }

        // Construct Egm96Grid via binary round-trip (internal field is not pub).
        let mut buf: Vec<u8> = Vec::with_capacity(EGM_NLAT * EGM_NLON * 4);
        for &v in &arr {
            buf.extend_from_slice(&v.to_le_bytes());
        }
        let tmp_path = path.with_file_name("egm96_tmp_convert.bin");
        std::fs::write(&tmp_path, &buf)?;
        let grid = Egm96Grid::load_binary(&tmp_path).map_err(GeoidFetchError::GridLoad)?;
        std::fs::remove_file(&tmp_path).map_err(GeoidFetchError::Io)?;

        Ok(grid)
    }

    // ── Save compact binary (EGM96) ────────────────────────────────────────

    fn save_binary(grid: &Egm96Grid, path: &Path) -> Result<(), GeoidFetchError> {
        // Re-export via sample queries at every grid node.
        // EGM_NLAT × EGM_NLON = 10 512 values.
        let mut buf: Vec<u8> = Vec::with_capacity(EGM_NLAT * EGM_NLON * 4);
        for r in 0..EGM_NLAT {
            let lat = 90.0 - r as f64 * 2.5;
            for c in 0..EGM_NLON {
                let lon = c as f64 * 2.5;
                // Query the grid at exact nodes — these will be identical to the
                // stored values (bilinear at integer coords returns exact values).
                let v = grid.undulation_m(lat, lon) as f32;
                buf.extend_from_slice(&v.to_le_bytes());
            }
        }
        let tmp = path.with_extension("tmp");
        std::fs::write(&tmp, &buf)?;
        std::fs::rename(&tmp, path)?;
        Ok(())
    }

    // ── EGM2008 COG parsing ────────────────────────────────────────────────

    /// Parse the PROJ CDN EGM2008 2.5-arc-minute COG and downsample to the
    /// 73×144 (2.5°) target grid used by [`Egm2008Grid`].
    ///
    /// The source file (`us_nga_egm08_25.tif`) has the same GeoTIFF structure
    /// as the EGM96 file but is much larger (4320×8640 pixels).
    fn cog_to_egm2008_grid(path: &Path) -> Result<Egm2008Grid, GeoidFetchError> {
        use std::io::BufReader;
        use tiff::decoder::{Decoder, DecodingResult};
        use tiff::tags::Tag;

        const TAG_MODEL_PIXEL_SCALE: u16 = 33550;
        const TAG_MODEL_TIEPOINT: u16 = 33922;

        let file = std::fs::File::open(path)?;
        let mut dec = Decoder::new(BufReader::new(file)).map_err(|e| GeoidFetchError::MalformedGtx {
            detail: e.to_string(),
        })?;

        let (width_u32, height_u32) = dec.dimensions().map_err(|e| GeoidFetchError::MalformedGtx {
            detail: e.to_string(),
        })?;
        let n_lon = width_u32 as usize;
        let n_lat = height_u32 as usize;

        if n_lat < COG_MIN_NLAT_EGM2008 || n_lon < COG_MIN_NLON_EGM2008 {
            return Err(GeoidFetchError::MalformedGtx {
                detail: format!(
                    "EGM2008 grid {n_lat}×{n_lon} is too small; \
                     expected ≥ {COG_MIN_NLAT_EGM2008}×{COG_MIN_NLON_EGM2008}"
                ),
            });
        }

        let scale = dec
            .get_tag_f64_vec(Tag::Unknown(TAG_MODEL_PIXEL_SCALE))
            .map_err(|e| GeoidFetchError::MalformedGtx {
                detail: format!("ModelPixelScaleTag: {e}"),
            })?;
        if scale.len() < 2 {
            return Err(GeoidFetchError::MalformedGtx {
                detail: format!(
                    "ModelPixelScaleTag has {} elements, expected ≥ 2",
                    scale.len()
                ),
            });
        }
        let d_lon = scale[0];
        let d_lat = scale[1];

        let tp = dec
            .get_tag_f64_vec(Tag::Unknown(TAG_MODEL_TIEPOINT))
            .map_err(|e| GeoidFetchError::MalformedGtx {
                detail: format!("ModelTiepointTag: {e}"),
            })?;
        if tp.len() < 6 {
            return Err(GeoidFetchError::MalformedGtx {
                detail: format!("ModelTiepointTag has {} elements, expected ≥ 6", tp.len()),
            });
        }
        let lon0 = tp[3] + (0.5 - tp[0]) * d_lon;
        let lat0 = tp[4] - (0.5 - tp[1]) * d_lat;

        // The EGM2008 COG is large (~120 MiB). Use unlimited decoder limits.
        use tiff::decoder::Limits;
        let file2 = std::fs::File::open(path)?;
        let mut dec2 = Decoder::new(BufReader::new(file2))
            .map_err(|e| GeoidFetchError::MalformedGtx { detail: e.to_string() })?
            .with_limits(Limits::unlimited());
        let image = dec2.read_image().map_err(|e| GeoidFetchError::MalformedGtx {
            detail: e.to_string(),
        })?;
        let pixels: Vec<f32> = match image {
            DecodingResult::F32(v) => v,
            DecodingResult::I16(v) => v.iter().map(|&x| x as f32).collect(),
            DecodingResult::U16(v) => v.iter().map(|&x| x as f32).collect(),
            other => {
                return Err(GeoidFetchError::MalformedGtx {
                    detail: format!("unsupported pixel type: {other:?}"),
                });
            }
        };

        // Nearest-neighbour downsample to the 73×144 (2.5°) target grid.
        let mut arr = [0f32; EGM_NLAT * EGM_NLON];
        for r in 0..EGM_NLAT {
            let lat = 90.0 - r as f64 * 2.5;
            let row_f = (lat0 - lat) / d_lat;
            let row_src = (row_f.round() as isize).clamp(0, (n_lat - 1) as isize) as usize;
            for c in 0..EGM_NLON {
                let lon = c as f64 * 2.5;
                let col_f = (lon - lon0) / d_lon;
                let col_src = (col_f.round() as isize).clamp(0, (n_lon - 1) as isize) as usize;
                arr[r * EGM_NLON + c] = pixels[row_src * n_lon + col_src];
            }
        }

        // Construct Egm2008Grid via binary round-trip (undulations_m is private).
        let mut buf: Vec<u8> = Vec::with_capacity(EGM_NLAT * EGM_NLON * 4);
        for &v in &arr {
            buf.extend_from_slice(&v.to_le_bytes());
        }
        let tmp_path = path.with_file_name("egm2008_tmp_convert.bin");
        std::fs::write(&tmp_path, &buf)?;
        let grid = Egm2008Grid::load_binary(&tmp_path).map_err(GeoidFetchError::GridLoad)?;
        std::fs::remove_file(&tmp_path).map_err(GeoidFetchError::Io)?;

        Ok(grid)
    }

    fn save_binary_egm2008(grid: &Egm2008Grid, path: &Path) -> Result<(), GeoidFetchError> {
        let mut buf: Vec<u8> = Vec::with_capacity(EGM_NLAT * EGM_NLON * 4);
        for r in 0..EGM_NLAT {
            let lat = 90.0 - r as f64 * 2.5;
            for c in 0..EGM_NLON {
                let lon = c as f64 * 2.5;
                let v = grid.undulation_m(lat, lon) as f32;
                buf.extend_from_slice(&v.to_le_bytes());
            }
        }
        let tmp = path.with_extension("tmp");
        std::fs::write(&tmp, &buf)?;
        std::fs::rename(&tmp, path)?;
        Ok(())
    }

} // end mod inner

// ── Feature-gated public surface ──────────────────────────────────

#[cfg(feature = "geoid-fetch")]
pub use inner::{fetch_egm96, fetch_egm2008, GeoidFetchError};
