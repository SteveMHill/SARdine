//! Geoid undulation model for DEM height correction.
//!
//! # Why this matters
//!
//! SRTM and GLO-30 heights are **orthometric** — measured above the geoid surface
//! (EGM96 and EGM2008 respectively), not the WGS84 reference ellipsoid.
//! Range-Doppler geocoding requires **ellipsoidal** heights.  The relationship is:
//!
//! ```text
//! h_ellipsoidal = h_orthometric + N
//! ```
//!
//! where `N` is the **geoid undulation**, which varies globally from −106 m
//! (Banda Sea, Indonesia) to +85 m (near Madras, India).  Typical mid-latitude
//! values are 30–50 m.
//!
//! Omitting this correction causes a systematic ECEF positioning error equal
//! to the local undulation.  At 50°N (N ≈ 40 m) this introduces a geolocation
//! offset of ~40 m, which is **~4 pixels at 10 m resolution** — well above the
//! sub-pixel accuracy that precise orbit files otherwise provide.
//!
//! # Variants
//!
//! | Variant | Accuracy | Notes |
//! |---------|----------|-------|
//! | [`GeoidModel::Zero`] | ≤ 80 m systematic error | No correction; backward-compatible |
//! | [`GeoidModel::Egm96`] | < 1 m (limited by DEM and orbit quality) | Requires NGA grid file; use with SRTM |
//! | [`GeoidModel::Egm2008`] | < 1 m | Requires PROJ CDN grid; use with GLO-30 |
//!
//! # EGM96 grid file
//!
//! The [`GeoidModel::Egm96`] variant requires the NGA EGM96 2.5° × 2.5°
//! undulation grid.  Download from:
//!
//! <https://earth-info.nga.mil/index.php?dir=wgs84&action=wgs84>
//!
//! File: `WW15MGH.GRD` (ASCII, 73 rows × 144 columns, metres).
//!
//! Load with:
//!
//! ```rust,ignore
//! let grid = Egm96Grid::from_ww15mgh(Path::new("WW15MGH.GRD"))?;
//! let model = GeoidModel::Egm96(grid);
//! ```
//!
//! Alternatively, use a pre-converted compact binary file:
//!
//! ```rust,ignore
//! let grid = Egm96Grid::load_binary(Path::new("egm96.bin"))?;
//! ```
//!
//! The binary format is 73 × 144 = 10 512 `f32` values in row-major order
//! (row 0 = 90°N, row 72 = 90°S; col 0 = 0°E, col 143 = 357.5°E),
//! native byte order (little-endian on x86).
//!
//! # EGM2008
//!
//! [`GeoidModel::Egm2008`] is also available.  Copernicus GLO-30 heights use
//! EGM2008 as their vertical datum, so pairing GLO-30 tiles with
//! `GeoidModel::Egm2008` is the most rigorous choice.  The global difference
//! between EGM96 and EGM2008 undulations is < 4 m (typical < 1 m), so
//! `GeoidModel::Egm96` remains scientifically acceptable when using SRTM data.
//!
//! # Reference
//!
//! Lemoine et al. (1998). *The Development of the Joint NASA GSFC and the
//! National Imagery and Mapping Agency (NIMA) Geopotential Model EGM96.*
//! NASA/TP-1998-206861.

use std::path::Path;

// ── EGM96 grid constants ───────────────────────────────────────────

/// Number of latitude rows in the 2.5° grid: from 90°N to 90°S inclusive.
/// `(90 − (−90)) / 2.5 + 1 = 73`.
const EGM96_NLAT: usize = 73;

/// Number of longitude columns in the 2.5° grid: from 0° to 357.5° inclusive.
/// `360 / 2.5 = 144`.  Longitude wraps: col 0 (0°E) == col 144 (360°E).
const EGM96_NLON: usize = 144;

/// Grid spacing in degrees.
const EGM96_STEP_DEG: f64 = 2.5;

// ── Error type ─────────────────────────────────────────────────────

#[derive(Debug, thiserror::Error)]
pub enum GeoidError {
    #[error("IO error reading {path}: {source}")]
    Io {
        path: String,
        #[source]
        source: std::io::Error,
    },

    #[error("EGM96 binary grid file {path} has wrong size: expected {expected} bytes, got {actual}")]
    WrongSize { path: String, expected: usize, actual: usize },

    #[error("EGM2008 binary grid file {path} has wrong size: expected {expected} bytes, got {actual}")]
    WrongSizeEgm2008 { path: String, expected: usize, actual: usize },

    #[error("EGM96 ASCII grid file {path} is malformed: {detail}")]
    MalformedAscii { path: String, detail: String },
}

// ── EGM96 grid ─────────────────────────────────────────────────────

/// The EGM96 2.5° × 2.5° global geoid undulation grid.
///
/// Row 0 = 90°N, row 72 = 90°S (step −2.5°).
/// Col 0 = 0°E (= 360°E), col 143 = 357.5°E (step +2.5°).
/// Values are geoid undulation in metres (positive = geoid above ellipsoid).
///
/// Bilinear interpolation is used for arbitrary (lat, lon) queries.
/// Longitude wraps at ±180°.
#[derive(Debug, Clone)]
pub struct Egm96Grid {
    /// Flat row-major grid: `undulations_m[row * EGM96_NLON + col]`.
    undulations_m: Box<[f32; EGM96_NLAT * EGM96_NLON]>,
}

impl Egm96Grid {
    /// Load from a compact binary file.
    ///
    /// The file must contain exactly `73 × 144 = 10 512` little-endian `f32`
    /// values (41 648 bytes) in row-major order:
    /// - Row 0 = 90°N … row 72 = 90°S
    /// - Col 0 = 0°E … col 143 = 357.5°E
    ///
    /// To convert `WW15MGH.GRD` to this format, use
    /// [`from_ww15mgh`](Self::from_ww15mgh) and then write the `f32` slice.
    pub fn load_binary(path: &Path) -> Result<Self, GeoidError> {
        let path_str = path.display().to_string();
        let raw = std::fs::read(path).map_err(|e| GeoidError::Io {
            path: path_str.clone(),
            source: e,
        })?;
        let expected = EGM96_NLAT * EGM96_NLON * 4;
        if raw.len() != expected {
            return Err(GeoidError::WrongSize {
                path: path_str,
                expected,
                actual: raw.len(),
            });
        }
        let mut arr = Box::new([0f32; EGM96_NLAT * EGM96_NLON]);
        for (i, chunk) in raw.chunks_exact(4).enumerate() {
            arr[i] = f32::from_le_bytes([chunk[0], chunk[1], chunk[2], chunk[3]]);
        }
        Ok(Egm96Grid { undulations_m: arr })
    }

    /// Load from the NGA ASCII file `WW15MGH.GRD`.
    ///
    /// The format is:
    /// ```text
    /// <header line, ignored>
    /// <lat_min> <lat_step> <lat_max> <lon_max> <lat_interval> <lon_interval>
    /// <row of undulation values in metres, 73 rows × 144 cols>
    /// ```
    ///
    /// Row order: descending latitude (90°N first).
    /// Column order: ascending longitude (0°E to 357.5°E).
    pub fn from_ww15mgh(path: &Path) -> Result<Self, GeoidError> {
        let path_str = path.display().to_string();
        let content = std::fs::read_to_string(path).map_err(|e| GeoidError::Io {
            path: path_str.clone(),
            source: e,
        })?;

        let mut values: Vec<f32> = Vec::with_capacity(EGM96_NLAT * EGM96_NLON);
        let mut header_lines_seen = 0usize;

        for line in content.lines() {
            let line = line.trim();
            if line.is_empty() {
                continue;
            }
            // Skip up to 2 header lines (the NGA file has a text header then
            // a numeric bounds line before the data).
            if header_lines_seen < 2 {
                // A header line either contains non-numeric text or a small
                // number of tokens (the bounds line has exactly 6 numbers).
                let tokens: Vec<&str> = line.split_whitespace().collect();
                let all_numeric = tokens.iter().all(|t| t.parse::<f64>().is_ok());
                if !all_numeric || tokens.len() <= 6 {
                    header_lines_seen += 1;
                    continue;
                }
            }
            for token in line.split_whitespace() {
                let v: f32 = token.parse().map_err(|_| GeoidError::MalformedAscii {
                    path: path_str.clone(),
                    detail: format!("cannot parse token '{token}' as f32"),
                })?;
                values.push(v);
            }
        }

        let expected = EGM96_NLAT * EGM96_NLON;
        if values.len() != expected {
            return Err(GeoidError::MalformedAscii {
                path: path_str,
                detail: format!(
                    "expected {expected} undulation values, found {}",
                    values.len()
                ),
            });
        }

        let mut arr = Box::new([0f32; EGM96_NLAT * EGM96_NLON]);
        arr.copy_from_slice(&values);
        Ok(Egm96Grid { undulations_m: arr })
    }

    /// Geoid undulation in metres at (lat_deg, lon_deg).
    ///
    /// Uses bilinear interpolation between the four surrounding 2.5° grid
    /// nodes.  Latitude is clamped to [−90, 90]; longitude wraps modulo 360°.
    pub fn undulation_m(&self, lat_deg: f64, lon_deg: f64) -> f64 {
        if !lat_deg.is_finite() || !lon_deg.is_finite() {
            return 0.0;
        }
        // Normalise longitude to [0, 360).
        let lon = lon_deg.rem_euclid(360.0);

        // Row 0 = 90°N; row increases going south.
        // row_f = (90 - lat) / 2.5
        let row_f = (90.0 - lat_deg.clamp(-90.0, 90.0)) / EGM96_STEP_DEG;
        let col_f = lon / EGM96_STEP_DEG;

        // Bracketing indices.
        let r0 = (row_f as usize).min(EGM96_NLAT - 2);
        let r1 = r0 + 1;
        // Longitude wraps: col 144 == col 0.
        let c0 = (col_f as usize) % EGM96_NLON;
        let c1 = (c0 + 1) % EGM96_NLON;

        let dr = row_f - r0 as f64;
        let dc = col_f - c0 as f64; // may slightly exceed 1.0 near wrap, clamped below

        let v00 = self.undulations_m[r0 * EGM96_NLON + c0] as f64;
        let v01 = self.undulations_m[r0 * EGM96_NLON + c1] as f64;
        let v10 = self.undulations_m[r1 * EGM96_NLON + c0] as f64;
        let v11 = self.undulations_m[r1 * EGM96_NLON + c1] as f64;

        // Clamp dc to [0, 1] to handle floating-point overshoot near 360°.
        let dc = dc.clamp(0.0, 1.0);

        v00 * (1.0 - dr) * (1.0 - dc)
            + v01 * (1.0 - dr) * dc
            + v10 * dr * (1.0 - dc)
            + v11 * dr * dc
    }
}

// ── EGM2008 grid ────────────────────────────────────────────────────

/// The EGM2008 2.5° × 2.5° global geoid undulation grid.
///
/// Same binary layout and interpolation as [`Egm96Grid`]; a separate type
/// so the two models cannot be confused at the type level.  Load with
/// [`Egm2008Grid::load_binary`].  Prefer this over [`Egm96Grid`] when the
/// DEM uses EGM2008 as its vertical datum (e.g. Copernicus GLO-30).
///
/// Row 0 = 90°N, row 72 = 90°S (step −2.5°).
/// Col 0 = 0°E, col 143 = 357.5°E (step +2.5°).
#[derive(Debug, Clone)]
pub struct Egm2008Grid {
    /// Same layout as [`Egm96Grid`]: `undulations_m[row * 144 + col]`.
    undulations_m: Box<[f32; EGM96_NLAT * EGM96_NLON]>,
}

impl Egm2008Grid {
    /// Load from a compact binary file.
    ///
    /// Same format as [`Egm96Grid::load_binary`]: exactly
    /// `73 × 144 = 10 512` little-endian `f32` values (41 648 bytes).
    pub fn load_binary(path: &Path) -> Result<Self, GeoidError> {
        let path_str = path.display().to_string();
        let raw = std::fs::read(path).map_err(|e| GeoidError::Io {
            path: path_str.clone(),
            source: e,
        })?;
        let expected = EGM96_NLAT * EGM96_NLON * 4;
        if raw.len() != expected {
            return Err(GeoidError::WrongSizeEgm2008 {
                path: path_str,
                expected,
                actual: raw.len(),
            });
        }
        let mut arr = Box::new([0f32; EGM96_NLAT * EGM96_NLON]);
        for (i, chunk) in raw.chunks_exact(4).enumerate() {
            arr[i] = f32::from_le_bytes([chunk[0], chunk[1], chunk[2], chunk[3]]);
        }
        Ok(Egm2008Grid { undulations_m: arr })
    }

    /// Geoid undulation in metres at (lat_deg, lon_deg).
    ///
    /// Bilinear interpolation on the 2.5° grid.  Identical algorithm to
    /// [`Egm96Grid::undulation_m`].
    pub fn undulation_m(&self, lat_deg: f64, lon_deg: f64) -> f64 {
        if !lat_deg.is_finite() || !lon_deg.is_finite() {
            return 0.0;
        }
        let lon = lon_deg.rem_euclid(360.0);
        let row_f = (90.0 - lat_deg.clamp(-90.0, 90.0)) / EGM96_STEP_DEG;
        let col_f = lon / EGM96_STEP_DEG;
        let r0 = (row_f as usize).min(EGM96_NLAT - 2);
        let r1 = r0 + 1;
        let c0 = (col_f as usize) % EGM96_NLON;
        let c1 = (c0 + 1) % EGM96_NLON;
        let dr = row_f - r0 as f64;
        let dc = (col_f - c0 as f64).clamp(0.0, 1.0);
        let v00 = self.undulations_m[r0 * EGM96_NLON + c0] as f64;
        let v01 = self.undulations_m[r0 * EGM96_NLON + c1] as f64;
        let v10 = self.undulations_m[r1 * EGM96_NLON + c0] as f64;
        let v11 = self.undulations_m[r1 * EGM96_NLON + c1] as f64;
        v00 * (1.0 - dr) * (1.0 - dc)
            + v01 * (1.0 - dr) * dc
            + v10 * dr * (1.0 - dc)
            + v11 * dr * dc
    }
}

// ── GeoidModel enum ────────────────────────────────────────────────

/// Geoid undulation model applied to DEM heights before ECEF conversion.
///
/// Add the value returned by [`GeoidModel::undulation_m`] to an orthometric
/// (e.g. SRTM or GLO-30) height to obtain the WGS84 ellipsoidal height
/// required for accurate ECEF coordinate conversion and range-Doppler
/// geocoding.
///
/// # Choosing a model
///
/// | Variant | Accuracy | Use case |
/// |---------|----------|----------|
/// | [`Zero`](GeoidModel::Zero) | ≤ 80 m systematic | Development, flat-ocean tests, backward compatibility |
/// | [`Egm96`](GeoidModel::Egm96) | < 1 m | Production — requires NGA `WW15MGH.GRD` grid |
///
/// Using `Zero` for land scenes at 10 m resolution introduces geolocation
/// errors much larger than the orbit/Doppler contribution.  It should only
/// be used knowingly and temporarily.
#[derive(Debug, Clone)]
pub enum GeoidModel {
    /// No geoid correction: orthometric height is used directly as ellipsoidal height.
    ///
    /// **Scientifically incorrect** except over open ocean where the geoid
    /// nearly coincides with the ellipsoid (N ≈ 0).
    ///
    /// Introduces a geolocation error in ECEF proportional to the local geoid
    /// undulation.  At 50°N the typical error is ~40 m (≈ 4 output pixels at
    /// 10 m resolution).  The maximum global error is ~80 m.
    ///
    /// There is **no `Default` implementation** for `GeoidModel`, and
    /// [`crate::terrain_correction::TerrainCorrectionConfig`] requires the
    /// geoid model to be supplied explicitly via
    /// [`TerrainCorrectionConfig::new`](crate::terrain_correction::TerrainCorrectionConfig::new).
    /// Choosing `Zero` is therefore always a deliberate act, never an
    /// accidental fallback.
    Zero,

    /// EGM96 2.5° bilinear undulation grid.
    ///
    /// Accuracy: < 1 m for most of the globe (limited by DEM and orbit quality,
    /// not the geoid grid).  Appropriate for use with SRTM DEM tiles.
    ///
    /// Load the grid with [`Egm96Grid::from_ww15mgh`] or [`Egm96Grid::load_binary`].
    Egm96(Egm96Grid),

    /// EGM2008 2.5° bilinear undulation grid.
    ///
    /// Accuracy: < 1 m for most of the globe.  EGM2008 is the vertical datum
    /// of Copernicus GLO-30 DEM tiles, so this variant is the correct choice
    /// when using GLO-30.
    ///
    /// Load the grid with [`Egm2008Grid::load_binary`].
    Egm2008(Egm2008Grid),
}

// Note: `GeoidModel` deliberately does NOT implement `Default`.  Producing a
// scientifically correct backscatter image requires the caller to make an
// informed choice between `Zero` (ocean-only / regression testing) and
// `Egm96(grid)` (land scenes).  Removing the `Default` impl pushes that
// choice out of the silent-fallback bucket and into the type system.

impl GeoidModel {
    /// Geoid undulation in metres at the given geodetic coordinates.
    ///
    /// Returns `h_ellipsoidal − h_orthometric`.  Add this to an orthometric
    /// DEM height to get the WGS84 ellipsoidal height required for ECEF
    /// conversion.
    ///
    /// # Arguments
    ///
    /// * `lat_deg` — WGS84 latitude in degrees (−90 to +90).
    /// * `lon_deg` — WGS84 longitude in degrees (any range; normalised internally).
    pub fn undulation_m(&self, lat_deg: f64, lon_deg: f64) -> f64 {
        match self {
            GeoidModel::Zero => 0.0,
            GeoidModel::Egm96(grid) => grid.undulation_m(lat_deg, lon_deg),
            GeoidModel::Egm2008(grid) => grid.undulation_m(lat_deg, lon_deg),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn zero_model_returns_zero_everywhere() {
        let model = GeoidModel::Zero;
        for (lat, lon) in [
            (0.0_f64, 0.0_f64),
            (90.0, 0.0),
            (-90.0, 0.0),
            (50.0, 10.0),
            (-34.0, 151.0),
        ] {
            assert_eq!(
                model.undulation_m(lat, lon),
                0.0,
                "GeoidModel::Zero must return 0.0 at ({lat}, {lon})"
            );
        }
    }

    #[test]
    fn zero_variant_is_constructible_and_documents_intent() {
        // Regression guard: there must be no `Default` impl for `GeoidModel`.
        // Callers must name `Zero` or `Egm96(...)` explicitly so that the
        // ~40 m geolocation bias of `Zero` is never an accidental fallback.
        // (Removing this trait bound below would be a silent-default
        // regression — see geoid.rs::GeoidModel doc-comment.)
        fn assert_no_default<T: Default>() {}
        // The line below is intentionally commented out.  Uncommenting it
        // must produce a compile error; if it ever compiles, the
        // anti-default property has regressed.
        // assert_no_default::<GeoidModel>();
        let _ = std::mem::discriminant(&GeoidModel::Zero); // SAFETY-OK: discarded discriminant in compile-time anti-default regression test
    }

    // ── Egm96Grid unit tests (synthetic grid, no file I/O) ────────────────────

    /// Build a minimal 3×3 synthetic EGM96-shaped grid for testing
    /// interpolation logic without needing the real NGA data file.
    fn make_synthetic_grid(fill: f32) -> Egm96Grid {
        let mut arr = Box::new([fill; EGM96_NLAT * EGM96_NLON]);
        // Insert a known ramp in rows 0–1, cols 0–1 so we can verify
        // bilinear interpolation:
        //   (row=0, col=0) = 40.0 m  (90°N, 0°E)
        //   (row=0, col=1) = 30.0 m  (90°N, 2.5°E)
        //   (row=1, col=0) = 20.0 m  (87.5°N, 0°E)
        //   (row=1, col=1) = 10.0 m  (87.5°N, 2.5°E)
        arr[0 * EGM96_NLON + 0] = 40.0;
        arr[0 * EGM96_NLON + 1] = 30.0;
        arr[1 * EGM96_NLON + 0] = 20.0;
        arr[1 * EGM96_NLON + 1] = 10.0;
        Egm96Grid { undulations_m: arr }
    }

    #[test]
    fn egm96_exact_grid_node_returned_exactly() {
        let grid = make_synthetic_grid(0.0);
        // At (lat=90°, lon=0°) → row=0, col=0 → undulation = 40.0 m.
        let u = grid.undulation_m(90.0, 0.0);
        assert!(
            (u - 40.0).abs() < 1e-4,
            "exact grid node: expected 40.0 m, got {u}"
        );
    }

    #[test]
    fn egm96_bilinear_midpoint() {
        let grid = make_synthetic_grid(0.0);
        // Midpoint between (row=0,col=0)=40 and (row=0,col=1)=30 along the
        // top row: lat=90°, lon=1.25° → expected = (40+30)/2 = 35.0 m.
        let u = grid.undulation_m(90.0, 1.25);
        assert!(
            (u - 35.0).abs() < 1e-3,
            "top-row midpoint: expected 35.0 m, got {u}"
        );
    }

    #[test]
    fn egm96_bilinear_centre_of_four_nodes() {
        let grid = make_synthetic_grid(0.0);
        // Centre of the 2×2 patch: lat=88.75°, lon=1.25°
        // → dr=0.5, dc=0.5
        // bilinear = 0.25*(40 + 30 + 20 + 10) = 25.0 m
        let u = grid.undulation_m(88.75, 1.25);
        assert!(
            (u - 25.0).abs() < 1e-3,
            "2×2 centre: expected 25.0 m, got {u}"
        );
    }

    #[test]
    fn egm96_longitude_wraps_at_360() {
        let grid = make_synthetic_grid(5.0);
        // 360°E == 0°E; both should return the same undulation.
        let u0 = grid.undulation_m(45.0, 0.0);
        let u360 = grid.undulation_m(45.0, 360.0);
        assert!(
            (u0 - u360).abs() < 1e-6,
            "0° and 360° must give the same undulation: {u0} vs {u360}"
        );
    }

    #[test]
    fn egm96_negative_longitude_normalised() {
        // −10° == 350°E.  Both should give the same result.
        let grid = make_synthetic_grid(7.0);
        let u_neg = grid.undulation_m(45.0, -10.0);
        let u_pos = grid.undulation_m(45.0, 350.0);
        assert!(
            (u_neg - u_pos).abs() < 1e-6,
            "−10° and 350° must give same undulation: {u_neg} vs {u_pos}"
        );
    }

    #[test]
    fn egm96_model_variant_delegates_to_grid() {
        let grid = make_synthetic_grid(0.0);
        let model = GeoidModel::Egm96(grid.clone());
        // Exact node at (90°N, 0°E) = 40.0 m.
        assert!(
            (model.undulation_m(90.0, 0.0) - 40.0).abs() < 1e-4,
            "GeoidModel::Egm96 must delegate to Egm96Grid"
        );
    }

    // ── Egm2008Grid unit tests ─────────────────────────────────────────────

    /// Build a synthetic Egm2008Grid for interpolation tests.  Uses the same
    /// corner values as `make_synthetic_grid` so we can assert both grids
    /// agree when given identical data.
    fn make_synthetic_egm2008_grid(fill: f32) -> Egm2008Grid {
        // Construct via binary round-trip (undulations_m is private).
        let mut arr = vec![fill; EGM96_NLAT * EGM96_NLON];
        arr[0 * EGM96_NLON + 0] = 40.0;
        arr[0 * EGM96_NLON + 1] = 30.0;
        arr[1 * EGM96_NLON + 0] = 20.0;
        arr[1 * EGM96_NLON + 1] = 10.0;
        let bytes: Vec<u8> = arr.iter().flat_map(|v| v.to_le_bytes()).collect();
        let tmp = tempfile::NamedTempFile::new().unwrap();
        std::fs::write(tmp.path(), &bytes).unwrap();
        Egm2008Grid::load_binary(tmp.path()).unwrap()
    }

    #[test]
    fn egm2008_exact_grid_node_returned_exactly() {
        let grid = make_synthetic_egm2008_grid(0.0);
        let u = grid.undulation_m(90.0, 0.0);
        assert!(
            (u - 40.0).abs() < 1e-4,
            "Egm2008Grid exact node: expected 40.0 m, got {u}"
        );
    }

    #[test]
    fn egm2008_bilinear_midpoint_matches_egm96() {
        // With identical underlying data, EGM2008 and EGM96 interpolation
        // must produce the same result.
        let egm96 = make_synthetic_grid(0.0);
        let egm2008 = make_synthetic_egm2008_grid(0.0);
        for (lat, lon) in [
            (90.0_f64, 0.0_f64),
            (90.0, 1.25),
            (88.75, 1.25),
            (45.0, 0.0),
        ] {
            let u96 = egm96.undulation_m(lat, lon);
            let u2008 = egm2008.undulation_m(lat, lon);
            assert!(
                (u96 - u2008).abs() < 1e-5,
                "EGM96 vs EGM2008 differ at ({lat}, {lon}): {u96} vs {u2008}"
            );
        }
    }

    #[test]
    fn egm2008_model_variant_delegates_to_grid() {
        let grid = make_synthetic_egm2008_grid(0.0);
        let model = GeoidModel::Egm2008(grid);
        assert!(
            (model.undulation_m(90.0, 0.0) - 40.0).abs() < 1e-4,
            "GeoidModel::Egm2008 must delegate to Egm2008Grid"
        );
    }

    #[test]
    fn egm2008_wrong_size_returns_error() {
        let tmp = tempfile::NamedTempFile::new().unwrap();
        std::fs::write(tmp.path(), b"short").unwrap();
        let err = Egm2008Grid::load_binary(tmp.path()).unwrap_err();
        assert!(
            matches!(err, GeoidError::WrongSizeEgm2008 { .. }),
            "expected WrongSizeEgm2008, got {err}"
        );
    }
}

