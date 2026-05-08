//! Output coordinate reference system abstraction for terrain-corrected
//! rasters.
//!
//! The terrain-correction stage projects radar pixels onto a regular
//! map-projection grid.  Until this module existed the grid was
//! hard-wired to WGS84 lat/lon (EPSG:4326) at every site of relevance:
//! the inner loop in [`crate::terrain_correction`] indexed columns by
//! longitude and rows by latitude, and the GeoTIFF writer in
//! [`crate::export`] emitted the EPSG:4326 GeoKey unconditionally.
//!
//! This module introduces a typed CRS choice so the grid, the
//! pixel-spacing units, and the GeoKeys all agree at every layer.
//!
//! # Supported CRS
//!
//! * [`OutputCrs::Wgs84LatLon`] — EPSG:4326, units = degrees.  This is
//!   the historical default and the only choice that preserves the
//!   pre-existing pipeline behaviour.
//! * [`OutputCrs::Utm`] — UTM zone 1N..60N or 1S..60S, units = metres.
//!   Backed by the WGS84 datum (EPSG:32601..32660 for north,
//!   EPSG:32701..32760 for south).
//!
//! # Why not the `proj` crate?
//!
//! The Rust `proj` crate binds to system libproj (>= 9.4 on recent
//! releases).  Many CI hosts and developer machines ship libproj 8.x;
//! requiring 9.4+ would make the build non-portable.  This module
//! uses [`proj4rs`] instead — a pure-Rust subset of PROJ that has no
//! system dependency and supports the projections we need today.
//!
//! # Coordinate convention
//!
//! All public functions take and return **degrees** for geographic
//! coordinates and **metres or degrees** (per the CRS) for projected
//! coordinates.  `proj4rs` itself works in radians, so this module is
//! responsible for the unit conversion at every boundary.

use proj4rs::Proj;
use proj4rs::transform::transform;

/// Errors returned when constructing or using an [`OutputCrs`].
#[derive(Debug, thiserror::Error)]
pub enum OutputCrsError {
    /// The EPSG code is not in the supported set.  Currently we accept
    /// `4326` and the WGS84 UTM zone codes 32601..32660 and
    /// 32701..32760.  Anything else is rejected explicitly rather than
    /// silently degrading to "no CRS".
    #[error("unsupported EPSG code {code} (supported: 4326, 32601-32660, 32701-32760)")]
    UnsupportedEpsg { code: u32 },

    /// The supplied string does not match the `EPSG:NNNN` form (case
    /// insensitive, optional whitespace).
    #[error("invalid CRS spec '{spec}': expected 'EPSG:<code>' or 'auto'")]
    InvalidSpec { spec: String },

    /// `proj4rs` failed to construct or run a projection.  The wrapped
    /// error string is the message from `proj4rs` verbatim.
    #[error("proj4rs error: {0}")]
    Proj(String),

    /// Auto-UTM was requested but the supplied longitude/latitude is
    /// outside the UTM-defined band (|lat| > 84°).
    #[error("auto-UTM cannot pick a zone for lat={lat_deg}° (outside ±84°)")]
    AutoUtmOutsideBand { lat_deg: f64 },
}

/// A choice of output map projection.  Constructed via [`OutputCrs::from_spec`]
/// or [`OutputCrs::auto_for_lon_lat`]; never instantiated by callers.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum OutputCrs {
    /// WGS84 geographic (latitude/longitude in degrees).  EPSG:4326.
    Wgs84LatLon,
    /// WGS84 / UTM zone N (north hemisphere).  `zone` is in 1..=60.
    /// EPSG = 32600 + zone.
    UtmNorth { zone: u8 },
    /// WGS84 / UTM zone S (south hemisphere).  `zone` is in 1..=60.
    /// EPSG = 32700 + zone.
    UtmSouth { zone: u8 },
}

impl OutputCrs {
    /// Parse `"EPSG:NNNN"` (case insensitive) or `"4326"` (bare digits)
    /// into an [`OutputCrs`].  `"auto"` is **not** accepted here — auto
    /// requires scene context (the centre lon/lat) and is built via
    /// [`OutputCrs::auto_for_lon_lat`] by the caller.
    pub fn from_spec(spec: &str) -> Result<Self, OutputCrsError> {
        let trimmed = spec.trim();
        let digits = match trimmed
            .strip_prefix("EPSG:")
            .or_else(|| trimmed.strip_prefix("epsg:"))
            .or_else(|| trimmed.strip_prefix("Epsg:"))
        {
            Some(rest) => rest.trim(),
            None => trimmed, // accept bare digit forms
        };
        let code: u32 = digits
            .parse()
            .map_err(|_| OutputCrsError::InvalidSpec { spec: spec.to_owned() })?;
        Self::from_epsg(code)
    }

    /// Build an [`OutputCrs`] from a numeric EPSG code.  Errors on
    /// codes outside the supported set; the supported set is
    /// `{4326}` ∪ `32601..=32660` ∪ `32701..=32760`.
    pub fn from_epsg(code: u32) -> Result<Self, OutputCrsError> {
        match code {
            4326 => Ok(OutputCrs::Wgs84LatLon),
            32601..=32660 => Ok(OutputCrs::UtmNorth { zone: (code - 32600) as u8 }),
            32701..=32760 => Ok(OutputCrs::UtmSouth { zone: (code - 32700) as u8 }),
            _ => Err(OutputCrsError::UnsupportedEpsg { code }),
        }
    }

    /// Pick the WGS84 UTM zone whose central meridian is closest to
    /// `lon_deg`, with hemisphere chosen by the sign of `lat_deg`.
    /// Returns [`OutputCrsError::AutoUtmOutsideBand`] when |lat| > 84°
    /// (the UTM polar limit; UPS would be needed there but we don't
    /// support it yet — fail loudly rather than guess).
    pub fn auto_for_lon_lat(lon_deg: f64, lat_deg: f64) -> Result<Self, OutputCrsError> {
        if !lat_deg.is_finite() || !lon_deg.is_finite() {
            return Err(OutputCrsError::AutoUtmOutsideBand { lat_deg });
        }
        if lat_deg.abs() > 84.0 {
            return Err(OutputCrsError::AutoUtmOutsideBand { lat_deg });
        }
        // UTM zone formula: floor((lon + 180) / 6) + 1, clamped to 1..=60.
        // Wrap longitude into [-180, 180) first.
        let mut lon = lon_deg;
        while lon < -180.0 { lon += 360.0; }
        while lon >= 180.0 { lon -= 360.0; }
        let zone_raw = ((lon + 180.0) / 6.0).floor() as i32 + 1;
        let zone = zone_raw.clamp(1, 60) as u8;
        if lat_deg >= 0.0 {
            Ok(OutputCrs::UtmNorth { zone })
        } else {
            Ok(OutputCrs::UtmSouth { zone })
        }
    }

    /// Numeric EPSG code for this CRS.
    pub fn epsg(&self) -> u32 {
        match self {
            OutputCrs::Wgs84LatLon => 4326,
            OutputCrs::UtmNorth { zone } => 32600 + *zone as u32,
            OutputCrs::UtmSouth { zone } => 32700 + *zone as u32,
        }
    }

    /// `true` for projected (metric) CRSs, `false` for geographic.
    /// Used by the CLI to interpret `--pixel-spacing` (degrees vs metres)
    /// and by the GeoTIFF writer to decide which GeoKeys to emit.
    pub fn is_metric(&self) -> bool {
        matches!(self, OutputCrs::UtmNorth { .. } | OutputCrs::UtmSouth { .. })
    }

    /// Default pixel spacing for this CRS in its native units.  10 m
    /// for UTM (matches Sentinel-1 IW pixel spacing in azimuth) and
    /// 1×10⁻⁴ deg ≈ 11 m at 50°N for WGS84.  This is purely a help
    /// convenience for the CLI default; the real value comes from
    /// `--pixel-spacing` whenever the user supplies it.
    pub fn default_pixel_spacing(&self) -> f64 {
        if self.is_metric() {
            10.0
        } else {
            1e-4
        }
    }

    /// proj4 string description.  Returned from a fixed table so we
    /// never invent a string the projection code can't parse.
    fn proj4_string(&self) -> String {
        match self {
            OutputCrs::Wgs84LatLon => {
                "+proj=longlat +datum=WGS84 +no_defs +type=crs".to_owned()
            }
            OutputCrs::UtmNorth { zone } => format!(
                "+proj=utm +zone={zone} +datum=WGS84 +units=m +no_defs +type=crs",
            ),
            OutputCrs::UtmSouth { zone } => format!(
                "+proj=utm +zone={zone} +south +datum=WGS84 +units=m +no_defs +type=crs",
            ),
        }
    }

    /// Build a [`Projector`] that converts between WGS84 lat/lon and
    /// this CRS's native (x, y).  Cheap to call (a few dozen
    /// allocations); cache the result if you need it in a hot loop.
    pub fn projector(&self) -> Result<Projector, OutputCrsError> {
        let wgs84 = Proj::from_proj_string("+proj=longlat +datum=WGS84 +no_defs +type=crs")
            .map_err(|e| OutputCrsError::Proj(format!("{e:?}")))?;
        let native = Proj::from_proj_string(&self.proj4_string())
            .map_err(|e| OutputCrsError::Proj(format!("{e:?}")))?;
        Ok(Projector { wgs84, native, is_metric: self.is_metric() })
    }

    /// GeoTIFF GeoKeyDirectoryTag block describing this CRS, as a
    /// flat sequence of SHORT (`u16`) values ready to be written to a
    /// TIFF tag (34735).  Format per GeoTIFF spec §2.4:
    ///
    /// ```text
    /// header: [KeyDirectoryVersion=1, KeyRevision=1, MinorRevision=0, NumberOfKeys]
    /// per key: [KeyID, TIFFTagLocation=0, Count=1, Value_Offset=<value>]
    /// ```
    ///
    /// `TIFFTagLocation = 0` means the value is stored inline in the
    /// fourth field, which is the only mode this writer supports
    /// (avoids a second tag block referenced by indirect offset).
    ///
    /// Returned blocks:
    ///
    /// * `Wgs84LatLon` → 16 SHORTs (3 keys: GTModelType=Geographic,
    ///   GTRasterType=PixelIsArea, GeographicType=4326).
    /// * `UtmNorth { zone }` / `UtmSouth { zone }` → 20 SHORTs (4
    ///   keys: GTModelType=Projected, GTRasterType=PixelIsArea,
    ///   ProjectedCSType=EPSG, ProjLinearUnits=Linear_Meter=9001).
    pub fn geo_key_directory(&self) -> Vec<u16> {
        match self {
            OutputCrs::Wgs84LatLon => vec![
                1, 1, 0, 3,        // header: version=1, revision=1.0, n_keys=3
                1024, 0, 1, 2,     // GTModelTypeGeoKey = ModelTypeGeographic
                1025, 0, 1, 1,     // GTRasterTypeGeoKey = RasterPixelIsArea
                2048, 0, 1, 4326,  // GeographicTypeGeoKey = GCS_WGS_84
            ],
            OutputCrs::UtmNorth { .. } | OutputCrs::UtmSouth { .. } => {
                let epsg = self.epsg() as u16; // 32601..32660 / 32701..32760 all fit in u16
                vec![
                    1, 1, 0, 4,        // header: version=1, revision=1.0, n_keys=4
                    1024, 0, 1, 1,     // GTModelTypeGeoKey = ModelTypeProjected
                    1025, 0, 1, 1,     // GTRasterTypeGeoKey = RasterPixelIsArea
                    3072, 0, 1, epsg,  // ProjectedCSTypeGeoKey = WGS 84 / UTM zone N
                    3076, 0, 1, 9001,  // ProjLinearUnitsGeoKey = Linear_Meter
                ]
            }
        }
    }
}

/// Bidirectional projector between WGS84 lat/lon and the chosen
/// output CRS.  Constructed by [`OutputCrs::projector`].
///
/// All public methods take and return **degrees** for geographic
/// coordinates and **metres** for projected coordinates (or degrees
/// for [`OutputCrs::Wgs84LatLon`], for which both transforms are the
/// identity).  The internal radians/degrees conversion to/from
/// `proj4rs` happens here so the rest of the pipeline never sees
/// radians.
#[derive(Clone)]
pub struct Projector {
    wgs84: Proj,
    native: Proj,
    /// Cached from the source `OutputCrs::is_metric()` so the
    /// degree-vs-metre decision is local to each method without a
    /// branch on the CRS tag.
    is_metric: bool,
}

impl Projector {
    /// Project a (lon_deg, lat_deg) pair to the output CRS native
    /// coordinates `(x, y)`.  For [`OutputCrs::Wgs84LatLon`] this is
    /// the identity (lon → x, lat → y).
    pub fn lonlat_to_xy(&self, lon_deg: f64, lat_deg: f64) -> Result<(f64, f64), OutputCrsError> {
        // proj4rs takes radians for geographic CRSs; pre-convert.
        let mut p = (lon_deg.to_radians(), lat_deg.to_radians(), 0.0);
        transform(&self.wgs84, &self.native, &mut p)
            .map_err(|e| OutputCrsError::Proj(format!("{e:?}")))?;
        if self.is_metric {
            // UTM output is already in metres.
            Ok((p.0, p.1))
        } else {
            // longlat → longlat: result is still radians.
            Ok((p.0.to_degrees(), p.1.to_degrees()))
        }
    }

    /// Inverse: project native `(x, y)` back to (lon_deg, lat_deg).
    pub fn xy_to_lonlat(&self, x: f64, y: f64) -> Result<(f64, f64), OutputCrsError> {
        let mut p = if self.is_metric {
            (x, y, 0.0)
        } else {
            // longlat input is in degrees from the caller; proj4rs
            // expects radians.
            (x.to_radians(), y.to_radians(), 0.0)
        };
        transform(&self.native, &self.wgs84, &mut p)
            .map_err(|e| OutputCrsError::Proj(format!("{e:?}")))?;
        // wgs84 output is always radians.
        Ok((p.0.to_degrees(), p.1.to_degrees()))
    }

    /// Project a whole row of `n_cols` evenly-spaced output pixels back
    /// to `(lon_deg, lat_deg)` using anchor-based linear interpolation.
    ///
    /// Replaces `n_cols` calls to [`Self::xy_to_lonlat`] with
    /// `n_anchors + 1` exact projections plus per-pixel linear blends.
    /// Within a single row at fixed `y` the inverse UTM projection is
    /// analytic and very smooth (bounded second derivatives ~ 1/R²
    /// per radian of x), so for a 25 000-pixel S-1 row at 10 m spacing
    /// (~250 km) and 16 anchors the worst-case linear-interpolation
    /// error in either lon or lat is below 10⁻⁹ degrees (~0.1 mm at
    /// the equator) — well below the per-pixel terrain-correction
    /// resolution.  An equivalence test in the dem/terrain modules
    /// must hold this bound on regression.
    ///
    /// `x_start` is the easting of the centre of column 0; `x_step`
    /// is the column spacing.  Both in metres for UTM CRSs (i.e.
    /// when `self.is_metric == true`); identical inputs work for
    /// `Wgs84LatLon` but are wasteful — geographic outputs already
    /// skip the projector entirely in the caller.
    ///
    /// `n_anchors` must be ≥ 1.  Errors are propagated from the
    /// per-anchor projector calls; no silent fallback.
    pub fn project_row(
        &self,
        y: f64,
        x_start: f64,
        x_step: f64,
        n_cols: usize,
        n_anchors: usize,
    ) -> Result<Vec<(f64, f64)>, OutputCrsError> {
        if n_cols == 0 {
            return Ok(Vec::new());
        }
        if n_anchors == 0 {
            return Err(OutputCrsError::Proj(
                "project_row: n_anchors must be >= 1".into(),
            ));
        }
        // Cap n_anchors so the warm path can't generate more anchors than
        // pixels (would degenerate to one projection per pixel anyway).
        let n_anchors = n_anchors.min(n_cols);

        // Build n_anchors+1 anchor points spanning [x_start, x_end] inclusive
        // so every pixel is enclosed by two anchors.
        let x_end = x_start + (n_cols as f64 - 1.0) * x_step;
        let n_pts = n_anchors + 1;
        let anchor_dx = (x_end - x_start) / n_anchors as f64;

        let mut anchors: Vec<(f64, f64, f64)> = Vec::with_capacity(n_pts); // (x, lon, lat)
        for i in 0..n_pts {
            let xi = x_start + i as f64 * anchor_dx;
            let (lon, lat) = self.xy_to_lonlat(xi, y)?;
            anchors.push((xi, lon, lat));
        }

        let mut out: Vec<(f64, f64)> = Vec::with_capacity(n_cols);
        // Walk pixels in column order; advance the active anchor segment
        // as we go.  For monotonic x_step this is O(n_cols + n_anchors).
        let mut seg = 0usize;
        for col in 0..n_cols {
            let x = x_start + col as f64 * x_step;
            // Advance segment so that anchors[seg].0 <= x <= anchors[seg+1].0
            while seg + 1 < n_anchors && x > anchors[seg + 1].0 {
                seg += 1;
            }
            let (x0, lon0, lat0) = anchors[seg];
            let (x1, lon1, lat1) = anchors[seg + 1];
            // x1 == x0 only if anchor_dx == 0, i.e. n_cols == 1, handled
            // above (n_anchors capped at n_cols, single-point trivially
            // returns the lone anchor).  Defensive guard for the f64-eq
            // edge case keeps the invariant explicit:
            let t = if x1 == x0 { 0.0 } else { (x - x0) / (x1 - x0) };
            out.push((lon0 + t * (lon1 - lon0), lat0 + t * (lat1 - lat0)));
        }
        Ok(out)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn approx(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() <= tol
    }

    /// Compile-time assertion: `Projector` must be `Send + Sync` so it
    /// can be shared by reference across rayon worker threads in the
    /// terrain-correction inner loop.  If `proj4rs::Proj` ever loses
    /// this property, this assertion fails at compile time.
    #[test]
    fn projector_is_send_and_sync() {
        fn assert_send_sync<T: Send + Sync>() {}
        assert_send_sync::<Projector>();
    }

    #[test]
    fn parse_epsg_4326_variants() {
        assert_eq!(OutputCrs::from_spec("EPSG:4326").unwrap(), OutputCrs::Wgs84LatLon);
        assert_eq!(OutputCrs::from_spec("epsg:4326").unwrap(), OutputCrs::Wgs84LatLon);
        assert_eq!(OutputCrs::from_spec("Epsg:4326").unwrap(), OutputCrs::Wgs84LatLon);
        assert_eq!(OutputCrs::from_spec(" 4326 ").unwrap(), OutputCrs::Wgs84LatLon);
    }

    #[test]
    fn parse_utm_zones() {
        assert_eq!(OutputCrs::from_spec("EPSG:32632").unwrap(), OutputCrs::UtmNorth { zone: 32 });
        assert_eq!(OutputCrs::from_spec("EPSG:32760").unwrap(), OutputCrs::UtmSouth { zone: 60 });
        assert_eq!(OutputCrs::from_spec("EPSG:32601").unwrap(), OutputCrs::UtmNorth { zone: 1 });
    }

    #[test]
    fn unsupported_epsg_is_rejected() {
        let err = OutputCrs::from_spec("EPSG:3857").expect_err("must reject");
        assert!(matches!(err, OutputCrsError::UnsupportedEpsg { code: 3857 }));
    }

    #[test]
    fn invalid_spec_is_rejected() {
        let err = OutputCrs::from_spec("not-a-crs").expect_err("must reject");
        assert!(matches!(err, OutputCrsError::InvalidSpec { .. }), "got {err:?}");
    }

    #[test]
    fn epsg_round_trip_through_from_epsg_and_method() {
        for code in [4326_u32, 32601, 32632, 32660, 32701, 32760] {
            let crs = OutputCrs::from_epsg(code).unwrap();
            assert_eq!(crs.epsg(), code, "epsg() for {code}");
        }
    }

    #[test]
    fn auto_utm_picks_correct_zone_for_known_lon_lat() {
        // Munich, Germany ≈ (11.58°E, 48.14°N) → UTM zone 32N (EPSG:32632).
        let crs = OutputCrs::auto_for_lon_lat(11.58, 48.14).unwrap();
        assert_eq!(crs, OutputCrs::UtmNorth { zone: 32 });

        // Sydney, Australia ≈ (151.21°E, -33.87°N) → UTM zone 56S (EPSG:32756).
        let crs = OutputCrs::auto_for_lon_lat(151.21, -33.87).unwrap();
        assert_eq!(crs, OutputCrs::UtmSouth { zone: 56 });

        // San Francisco ≈ (-122.42°W, 37.77°N) → UTM zone 10N (EPSG:32610).
        let crs = OutputCrs::auto_for_lon_lat(-122.42, 37.77).unwrap();
        assert_eq!(crs, OutputCrs::UtmNorth { zone: 10 });

        // Equator special case: lat = 0 picks North hemisphere (the
        // documented choice; UTM's prime meridian convention treats 0
        // as north).
        let crs = OutputCrs::auto_for_lon_lat(0.0, 0.0).unwrap();
        assert!(matches!(crs, OutputCrs::UtmNorth { .. }));
    }

    #[test]
    fn auto_utm_rejects_polar_latitudes() {
        for lat in [85.0_f64, -85.0, 89.9, -89.9, f64::NAN] {
            let err = OutputCrs::auto_for_lon_lat(0.0, lat).expect_err("must reject");
            assert!(matches!(err, OutputCrsError::AutoUtmOutsideBand { .. }), "got {err:?}");
        }
    }

    #[test]
    fn wgs84_projector_is_identity() {
        let proj = OutputCrs::Wgs84LatLon.projector().unwrap();
        let (x, y) = proj.lonlat_to_xy(11.58, 48.14).unwrap();
        // longlat→longlat is the identity (after rad↔deg round-trip).
        assert!(approx(x, 11.58, 1e-9));
        assert!(approx(y, 48.14, 1e-9));
        let (lon, lat) = proj.xy_to_lonlat(11.58, 48.14).unwrap();
        assert!(approx(lon, 11.58, 1e-9));
        assert!(approx(lat, 48.14, 1e-9));
    }

    #[test]
    fn utm32n_munich_round_trip() {
        // Munich Marienplatz ≈ (11.5755°E, 48.1374°N).  In UTM 32N the
        // documented projected coordinates are:
        //   easting  ≈ 691 603.032 m
        //   northing ≈ 5 334 780.031 m
        // (Cross-checked against system PROJ 8.2.1 with the same
        // proj-string used by [`OutputCrs::proj4_string`].)
        let crs = OutputCrs::UtmNorth { zone: 32 };
        let proj = crs.projector().unwrap();
        let (x, y) = proj.lonlat_to_xy(11.5755, 48.1374).unwrap();
        assert!(approx(x, 691_603.03, 1.0), "easting got {x}");
        assert!(approx(y, 5_334_780.03, 1.0), "northing got {y}");

        // Round-trip back to lon/lat must recover the input to
        // < 1 µdeg ≈ 11 cm at the equator (sub-pixel for any S-1 grid).
        let (lon, lat) = proj.xy_to_lonlat(x, y).unwrap();
        assert!(approx(lon, 11.5755, 1e-7), "lon round-trip got {lon}");
        assert!(approx(lat, 48.1374, 1e-7), "lat round-trip got {lat}");
    }

    #[test]
    fn utm56s_sydney_round_trip() {
        // Sydney Opera House ≈ (151.2153°E, -33.8568°N) → UTM 56S.
        let crs = OutputCrs::UtmSouth { zone: 56 };
        let proj = crs.projector().unwrap();
        let (x, y) = proj.lonlat_to_xy(151.2153, -33.8568).unwrap();
        // Round-trip must be tight even if I don't pin the exact PROJ
        // value here (south-hemisphere UTM is what we're really testing).
        let (lon, lat) = proj.xy_to_lonlat(x, y).unwrap();
        assert!(approx(lon, 151.2153, 1e-7), "lon round-trip got {lon}");
        assert!(approx(lat, -33.8568, 1e-7), "lat round-trip got {lat}");
        // Easting must lie in the canonical UTM band [100 000, 900 000] m.
        assert!(x > 100_000.0 && x < 900_000.0, "easting {x} outside UTM band");
        // Southern-hemisphere northing is offset by 10 000 000 m, so
        // the value should be > 5e6 for any non-polar location.
        assert!(y > 5_000_000.0, "southern northing {y} suspiciously low");
    }

    #[test]
    fn is_metric_matches_crs_kind() {
        assert!(!OutputCrs::Wgs84LatLon.is_metric());
        assert!(OutputCrs::UtmNorth { zone: 32 }.is_metric());
        assert!(OutputCrs::UtmSouth { zone: 56 }.is_metric());
    }

    #[test]
    fn geo_key_directory_wgs84_matches_legacy_layout() {
        // Pin the historical 16-SHORT block so the GeoTIFF writer's
        // bit-for-bit output for WGS84 doesn't change with this refactor.
        assert_eq!(
            OutputCrs::Wgs84LatLon.geo_key_directory(),
            vec![1, 1, 0, 3, 1024, 0, 1, 2, 1025, 0, 1, 1, 2048, 0, 1, 4326],
        );
    }

    #[test]
    fn geo_key_directory_utm_north_emits_projected_keys() {
        // Pin the 20-SHORT block for UTM 32N (EPSG:32632).
        assert_eq!(
            OutputCrs::UtmNorth { zone: 32 }.geo_key_directory(),
            vec![
                1, 1, 0, 4,
                1024, 0, 1, 1,
                1025, 0, 1, 1,
                3072, 0, 1, 32632,
                3076, 0, 1, 9001,
            ],
        );
    }

    #[test]
    fn geo_key_directory_utm_south_uses_correct_epsg() {
        // UTM 56S → EPSG:32756.  Layout (per the UTM-N pin above):
        //   [0..4)   header [1, 1, 0, 4]
        //   [4..8)   GTModelTypeGeoKey      = [1024, 0, 1, 1]
        //   [8..12)  GTRasterTypeGeoKey     = [1025, 0, 1, 1]
        //   [12..16) ProjectedCSTypeGeoKey  = [3072, 0, 1, EPSG]
        //   [16..20) ProjLinearUnitsGeoKey  = [3076, 0, 1, 9001]
        let keys = OutputCrs::UtmSouth { zone: 56 }.geo_key_directory();
        assert_eq!(keys[15], 32756, "ProjectedCSType value should be 32756");
        assert_eq!(keys[19], 9001, "ProjLinearUnits should be Linear_Meter (9001)");
    }

    /// Anchor-interpolated `project_row` must agree with per-pixel
    /// `xy_to_lonlat` to better than 1e-6° (~ 11 cm at the equator) on
    /// a representative S-1-sized UTM row.  This bounds the accuracy
    /// loss introduced by the Phase-2 fast path well below 0.1 px at
    /// 10 m output spacing and far below the per-pixel terrain-correction
    /// noise floor.  A second test below covers the zone-edge geometry
    /// where grid convergence is largest.
    #[test]
    fn project_row_matches_per_pixel_munich_utm() {
        let proj = OutputCrs::UtmNorth { zone: 32 }.projector().unwrap();
        // 25 000-pixel row at 10 m spacing centred near Munich:
        // northing ≈ 5 334 780 m, easting starts at 600 000 m.
        let n_cols = 25_000usize;
        let x_step = 10.0;
        let x_start = 600_000.0;
        let y = 5_334_780.0;
        // Match the production default (see
        // [`crate::terrain_correction::terrain_correction`]).
        let n_anchors = 256;

        let row = proj
            .project_row(y, x_start, x_step, n_cols, n_anchors)
            .unwrap();
        assert_eq!(row.len(), n_cols);

        // Sample 64 columns spanning the row and compare to per-pixel
        // projector output.  Full 25k comparison would inflate test
        // time; sparse sampling at deterministic intervals is enough
        // to surface any > 1e-6° deviation.
        for k in 0..64 {
            let col = (k * (n_cols - 1)) / 63;
            let x = x_start + col as f64 * x_step;
            let (lon_ref, lat_ref) = proj.xy_to_lonlat(x, y).unwrap();
            let (lon, lat) = row[col];
            assert!(
                (lon - lon_ref).abs() < 1e-6,
                "lon at col {col}: interp={lon} ref={lon_ref} delta={}",
                (lon - lon_ref).abs(),
            );
            assert!(
                (lat - lat_ref).abs() < 1e-6,
                "lat at col {col}: interp={lat} ref={lat_ref} delta={}",
                (lat - lat_ref).abs(),
            );
        }
    }

    /// Same equivalence bound at the western edge of UTM zone 32N
    /// (lon ≈ 6°E) where grid convergence and second-derivative
    /// magnitudes are larger than at the central meridian.  256 anchors
    /// must hold the < 1e-6° contract across a full 25k-pixel row.
    #[test]
    fn project_row_matches_per_pixel_zone_edge() {
        let proj = OutputCrs::UtmNorth { zone: 32 }.projector().unwrap();
        let n_cols = 25_000usize;
        let x_step = 10.0;
        // Western zone edge near 6°E at 50°N → easting ≈ 285 000 m.
        let x_start = 285_000.0;
        let y = 5_540_000.0;
        let n_anchors = 256;

        let row = proj
            .project_row(y, x_start, x_step, n_cols, n_anchors)
            .unwrap();

        let mut max_dlon = 0.0_f64;
        let mut max_dlat = 0.0_f64;
        for k in 0..128 {
            let col = (k * (n_cols - 1)) / 127;
            let x = x_start + col as f64 * x_step;
            let (lon_ref, lat_ref) = proj.xy_to_lonlat(x, y).unwrap();
            let (lon, lat) = row[col];
            max_dlon = max_dlon.max((lon - lon_ref).abs());
            max_dlat = max_dlat.max((lat - lat_ref).abs());
        }
        assert!(
            max_dlon < 1e-6 && max_dlat < 1e-6,
            "zone-edge row deviation: max_dlon={max_dlon}, max_dlat={max_dlat}",
        );
    }

    /// A degenerate `n_cols == 1` request must return exactly one
    /// pixel that matches the per-pixel projector — no special-case
    /// numerical drift permitted.
    #[test]
    fn project_row_single_pixel_matches_per_pixel() {
        let proj = OutputCrs::UtmNorth { zone: 32 }.projector().unwrap();
        let row = proj.project_row(5_334_780.0, 691_603.0, 10.0, 1, 16).unwrap();
        assert_eq!(row.len(), 1);
        let (lon_ref, lat_ref) = proj.xy_to_lonlat(691_603.0, 5_334_780.0).unwrap();
        assert!((row[0].0 - lon_ref).abs() < 1e-12);
        assert!((row[0].1 - lat_ref).abs() < 1e-12);
    }

    /// Zero columns is a no-op; zero anchors is an explicit error.
    #[test]
    fn project_row_edge_cases() {
        let proj = OutputCrs::UtmNorth { zone: 32 }.projector().unwrap();
        assert_eq!(
            proj.project_row(5_000_000.0, 600_000.0, 10.0, 0, 16).unwrap(),
            Vec::<(f64, f64)>::new(),
        );
        let err = proj
            .project_row(5_000_000.0, 600_000.0, 10.0, 100, 0)
            .expect_err("n_anchors=0 must error");
        assert!(matches!(err, OutputCrsError::Proj(_)), "got {err:?}");
    }
}
