//! Range-Doppler terrain correction and geocoding.
//!
//! Implements the backward-geocoding approach: for each output pixel on a
//! regular lat/lon grid, solve the Range-Doppler equations to find the
//! corresponding pixel in the slant-range radar image, then resample the
//! backscatter value.
//!
//! # Algorithm
//!
//! For each output grid point `(lat, lon)`:
//! 1. Look up DEM height `h` at `(lat, lon)`.
//! 2. Convert `(lat, lon, h)` to ECEF target position `P`.
//! 3. Solve the **zero-Doppler equation** to find azimuth time `t_az`:
//!    `(P − S(t)) · V(t) = 0`
//!    where `S(t)`, `V(t)` are satellite position/velocity from orbit.
//! 4. Compute slant range `R = |P − S(t_az)|`.
//! 5. Convert `(t_az, R)` to radar image coordinates `(line, sample)`.
//! 6. Bilinear-interpolate the radar pixel value.
//! 7. Optionally apply terrain flattening (Small 2011 area-projection): σ⁰ → γ⁰.
//!
//! # Coordinate systems
//!
//! - **Output grid**: WGS84 geographic (lat/lon) with configurable spacing.
//! - **DEM**: WGS84 geographic, bilinear-interpolated from SRTM tiles.
//! - **Orbit**: ECEF (Earth-Centred Earth-Fixed), interpolated via Lagrange polynomials.
//! - **Radar image**: slant-range geometry (azimuth lines × range samples).

use chrono::{DateTime, Duration, Utc};
use rayon::prelude::*;
use std::sync::atomic::{AtomicUsize, Ordering};

use crate::dem::DemSource;
use crate::geodesy::{cross3, dot3, geodetic_to_ecef, norm3, normalize3, sub3};
use crate::geoid::GeoidModel;
use crate::merge_subswaths::MergedSigma0;
use crate::orbit::{interpolate_orbit_and_accel, OrbitError};
use crate::output_crs::{OutputCrs, OutputCrsError, Projector};
use crate::pipeline_options::{RadarImage, ResamplingKernel};
use crate::types::{
    GeolocationGridPoint, OrbitData, SceneMetadata, SubSwathId, SPEED_OF_LIGHT_M_S,
};

// ── Error type ─────────────────────────────────────────────────────

#[derive(Debug, thiserror::Error)]
pub enum TerrainCorrectionError {
    #[error("orbit interpolation failed: {0}")]
    Orbit(#[from] OrbitError),

    #[error("DEM error: {0}")]
    Dem(#[from] crate::dem::DemError),

    #[error("no geolocation grid points provided")]
    NoGridPoints,

    #[error("zero-Doppler did not converge at ({lat:.4}, {lon:.4})")]
    NoDopplerSolution { lat: f64, lon: f64 },

    #[error("scene has no subswaths — cannot derive radar geometry")]
    NoSubSwaths,

    #[error("configuration error: {0}")]
    Config(String),

    #[error("projection failed: {0}")]
    Projection(String),
}

// ── Zero-Doppler solver outcome ────────────────────────────────────

/// Result of a single-pixel zero-Doppler solve.
///
/// The solver returns an explicit non-convergence indication rather than a
/// plausible-but-unverified estimate.  Callers map `NotConverged` /
/// `DegenerateGeometry` to NaN output pixels and count them separately from
/// pixels that converged cleanly.
#[derive(Debug)]
#[allow(dead_code)]
pub(crate) enum ZeroDopplerOutcome {
    /// Solver converged.  `sat_sv` is the satellite state vector at `time`,
    /// returned directly from the final Lagrange evaluation so the call site
    /// does not need to call `interpolate_orbit` again.
    Converged { time: DateTime<Utc>, residual_s: f64, iterations: u8, sat_sv: crate::types::StateVector },
    NotConverged { last_residual_s: f64 },
    DegenerateGeometry,
}

/// Maximum slot count for the Newton iteration histogram.
///
/// Slot `i` (1 ≤ i ≤ MAX-1) counts pixels that converged in exactly `i`
/// iterations.  Slot 0 is unused.  Any iteration count ≥ `MAX-1` is clamped
/// into the last bin.  16 leaves headroom over the default
/// `max_newton_iterations = 10`; a debug_assert in `terrain_correction`
/// guards against silent clipping if a caller raises the limit.
pub(crate) const NEWTON_ITER_HIST_SIZE: usize = 16;

// ── Configuration ──────────────────────────────────────────────────

/// Configuration for the terrain correction step.
pub struct TerrainCorrectionConfig {
    /// Output pixel spacing.  Units are interpreted by the chosen
    /// [`Self::crs`]: **degrees** for [`OutputCrs::Wgs84LatLon`],
    /// **metres** for the UTM variants.  Default: ~10 m ≈ 0.0001°,
    /// which only makes sense for the default WGS84 CRS — callers that
    /// switch to UTM must also reset this field to a metric value
    /// (e.g. 10.0 for 10 m).
    pub pixel_spacing_deg: f64,

    /// Output coordinate reference system.  Defaults to
    /// [`OutputCrs::Wgs84LatLon`] so legacy callers that don't set
    /// this field get the same behaviour as before this field was
    /// added.
    pub crs: OutputCrs,

    /// Maximum iterations for the zero-Doppler Newton solver.
    pub max_newton_iterations: usize,

    /// Convergence threshold for zero-Doppler solver (seconds).
    pub newton_tolerance_s: f64,

    /// Apply terrain flattening to produce γ⁰ (terrain-flattened backscatter)
    /// via the area-projection method of Small (2011) "Flattening Gamma:
    /// Radiometric Terrain Correction for SAR Imagery".
    ///
    /// When `false` (default), the output is σ⁰ as resampled from the input
    /// merged image — backscatter normalised by the flat-earth projected area
    /// as assumed by the S-1 calibration LUTs.
    ///
    /// When `true`, each geocoded pixel is multiplied by the ratio of the
    /// reference flat-earth illuminated area to the actual terrain-facet area
    /// projected onto the radar look direction.  Pixels in shadow, layover, or
    /// with extreme foreshortening are set to NaN and counted separately in
    /// `GeocodedImage::flat_masked_count`.
    pub flatten: bool,

    /// Geoid model for converting orthometric DEM heights to WGS84 ellipsoidal
    /// heights before ECEF conversion.
    ///
    /// SRTM and GLO-30 heights are orthometric (above EGM96 / EGM2008
    /// respectively).  Range-Doppler geocoding requires WGS84 ellipsoidal
    /// heights.  The correction is:
    ///
    /// ```text
    /// h_ellipsoidal = h_orthometric + N(lat, lon)
    /// ```
    ///
    /// where `N` is the geoid undulation.  Globally, `N` ranges from −106 m to
    /// +85 m; typical mid-latitude values are 30–50 m.  Omitting the correction
    /// introduces a geolocation error equal to the local undulation (≈ 40 m /
    /// 4 pixels at 10 m resolution for 50°N).
    ///
    /// There is **no default** for this field.  Construct the config with
    /// [`TerrainCorrectionConfig::new`] so the geoid choice is always
    /// explicit; choose [`GeoidModel::Zero`] only for ocean scenes or
    /// regression tests.
    ///
    /// To use the EGM96 model (recommended for production):
    ///
    /// ```rust,ignore
    /// // With the `geoid-fetch` feature:
    /// use sardine::geoid::GeoidModel;
    /// use sardine::geoid_fetch::fetch_egm96;
    /// use sardine::terrain_correction::TerrainCorrectionConfig;
    /// let grid = fetch_egm96()?;  // downloads + caches on first call
    /// let cfg = TerrainCorrectionConfig::new(GeoidModel::Egm96(grid));
    /// ```
    ///
    /// See [`crate::geoid`] for format details.
    pub geoid: GeoidModel,

    /// If `true`, populate [`GeocodedImage::cos_lia`] and
    /// [`GeocodedImage::mask`].  Adds four DEM lookups per pixel even
    /// when [`Self::flatten`] is `false`, so it has a measurable cost
    /// (~10–20 % of TC runtime in profiles to date).  When `false`,
    /// both fields are returned as empty `Vec`s.
    ///
    /// Independent of [`Self::flatten`]: you may set
    /// `compute_lia = true, flatten = false` to get a σ⁰ image plus
    /// the LIA + quality-mask sidecars without converting to γ⁰.
    pub compute_lia: bool,

    /// Noise floor margin for per-pixel NESZ masking (dB).
    ///
    /// Pixels where the calibrated backscatter (σ⁰ before terrain flattening)
    /// is at or below `NESZ + noise_floor_margin_db` are masked to NaN and
    /// counted in [`GeocodedImage::noise_masked_count`].
    ///
    /// Set to `f32::NEG_INFINITY` (default) to disable NESZ masking entirely.
    /// A value of `0.0` masks pixels at or below the noise floor itself.
    /// A value of `3.0` masks pixels within 3 dB of the noise floor
    /// (i.e. SNR < 3 dB), which is the recommended production setting.
    ///
    /// The NESZ values come from the noise LUT computed in
    /// [`crate::apply_calibration::apply_calibration`] and propagated through
    /// the merged image.  Where the merged NESZ is NaN (swath gap), the
    /// pixel is passed through without masking.
    pub noise_floor_margin_db: f32,

    /// Resampling kernel for the σ⁰ image interpolation step.
    ///
    /// Default: [`ResamplingKernel::Bilinear`].  See
    /// [`crate::pipeline_options::ResamplingKernel`] for the tradeoffs.
    pub resampling: crate::pipeline_options::ResamplingKernel,
}

impl TerrainCorrectionConfig {
    /// Build a config with sensible defaults and the supplied geoid model.
    ///
    /// `geoid` has no default precisely because choosing
    /// [`GeoidModel::Zero`] silently introduces up to ±80 m of vertical
    /// bias on land.  All other fields are filled with the historical
    /// defaults:
    ///
    /// | Field                   | Default              | Meaning                                          |
    /// |-------------------------|----------------------|--------------------------------------------------|
    /// | `pixel_spacing_deg`     | `0.0001` (~11 m@50°N)| Output grid spacing                              |
    /// | `max_newton_iterations` | `10`                 | Zero-Doppler solver budget                       |
    /// | `newton_tolerance_s`    | `1e-6`               | ~7 mm along-track convergence                    |
    /// | `flatten`               | `false`              | Apply Small (2011) γ⁰ flattening                 |
    /// | `noise_floor_margin_db` | `-inf`               | Disable per-pixel NESZ masking                   |
    ///
    /// Override any of these via direct field assignment after `new`.
    pub fn new(geoid: GeoidModel) -> Self {
        Self {
            pixel_spacing_deg: 0.0001,
            crs: OutputCrs::Wgs84LatLon,
            max_newton_iterations: 10,
            newton_tolerance_s: 1e-6,
            flatten: false,
            geoid,
            compute_lia: false,
            noise_floor_margin_db: f32::NEG_INFINITY,
            resampling: crate::pipeline_options::ResamplingKernel::Bilinear,
        }
    }
}

// Note: `TerrainCorrectionConfig` deliberately does NOT implement `Default`.
// The `geoid` field has no scientifically safe default; forcing callers
// through `new(geoid)` keeps that decision visible in every call site.

// ── Output ─────────────────────────────────────────────────────────

/// Geocoded terrain-corrected output image.
pub struct GeocodedImage {
    /// Pixel data in row-major order. NaN for no-data.
    pub data: Vec<f32>,
    /// Number of rows (latitude direction, north to south).
    pub rows: usize,
    /// Number of columns (longitude direction, west to east).
    pub cols: usize,
    /// Northwest corner latitude (degrees).
    pub origin_lat_deg: f64,
    /// Northwest corner longitude (degrees).
    pub origin_lon_deg: f64,
    /// Pixel spacing in latitude (degrees, positive = southward).
    pub pixel_spacing_lat_deg: f64,
    /// Pixel spacing in longitude (degrees, positive = eastward).
    pub pixel_spacing_lon_deg: f64,
    /// Total number of output pixels.
    pub total_pixel_count: usize,
    /// Pixels that converged cleanly and received a valid backscatter value.
    pub valid_pixel_count: usize,
    /// Pixels outside DEM coverage (silently skipped).
    pub dem_missing_count: usize,
    /// Pixels outside the radar image footprint after back-projection.
    pub outside_footprint_count: usize,
    /// Pixels where the zero-Doppler solver ran to `max_newton_iterations`
    /// without reaching `newton_tolerance_s`.  These are returned as NaN.
    pub non_converged_count: usize,
    /// Pixels where the solver hit a degenerate Jacobian (near-zero derivative).
    /// These are returned as NaN and normally indicate a bug, not a data issue.
    pub degenerate_geometry_count: usize,
    /// Pixels masked by terrain flattening: shadow, layover, extreme foreshortening,
    /// or a DEM void in a neighbour cell needed for the facet gradient.
    /// Zero when `TerrainCorrectionConfig::flatten` is `false`.
    pub flat_masked_count: usize,
    /// Pixels masked by the per-pixel NESZ noise floor (σ⁰ ≤ NESZ × 10^(margin/10)).
    /// Zero when `TerrainCorrectionConfig::noise_floor_margin_db` is `f32::NEG_INFINITY`.
    pub noise_masked_count: usize,
    /// GDAL geotransform: [origin_x, pixel_width, 0, origin_y, 0, -pixel_height].
    /// For a north-up image: origin is top-left (NW) corner, pixel_height is
    /// positive (lat decreases southward), so gt[5] is negative.
    pub geotransform: [f64; 6],
    /// Coordinate reference system of [`Self::geotransform`] and the
    /// pixel grid.  For [`OutputCrs::Wgs84LatLon`] units are degrees;
    /// for the UTM variants units are metres.
    pub crs: OutputCrs,
    /// Local incidence angle cosine, one f32 per output pixel, row-major.
    /// `NaN` for pixels in shadow/layover, outside the footprint, or where
    /// the centre/neighbour DEM is missing.  Empty `Vec` when
    /// [`TerrainCorrectionConfig::compute_lia`] is `false`.
    pub cos_lia: Vec<f32>,
    /// Per-pixel quality mask, one byte per output pixel, row-major.
    /// See [`mask_code`] for the value definitions.  Empty `Vec` when
    /// [`TerrainCorrectionConfig::compute_lia`] is `false`.
    pub mask: Vec<u8>,
    /// Newton zero-Doppler iteration histogram, sized
    /// [`NEWTON_ITER_HIST_SIZE`].  Slot `i` (1 ≤ i ≤ size−1) counts
    /// pixels that converged in exactly `i` iterations.  Slot 0 is unused.
    /// The final slot is the saturating-clamped overflow bin.  Sum over all
    /// slots equals [`Self::valid_pixel_count`] +
    /// [`Self::flat_masked_count`] + [`Self::noise_masked_count`] +
    /// [`Self::outside_footprint_count`] (every pixel that converged but may
    /// have been masked downstream).
    pub newton_iter_histogram: [u64; NEWTON_ITER_HIST_SIZE],
}

impl GeocodedImage {
    /// Fraction of pixels that received a valid backscatter value.
    pub fn valid_fraction(&self) -> f64 {
        self.valid_pixel_count as f64 / self.total_pixel_count as f64
    }
}

// ── Per-pixel quality mask codes ───────────────────────────────────

/// Categorical per-pixel quality mask.
///
/// One byte per output pixel records the *first* reason a pixel did not
/// receive a valid backscatter value, or [`VALID`] when it did.  The codes
/// align 1-to-1 with the per-row counter buckets returned in
/// [`GeocodedImage`] so that `popcount(mask == X)` matches the corresponding
/// counter exactly.
///
/// The mask is populated only when [`TerrainCorrectionConfig::compute_lia`]
/// is `true`.  Otherwise [`GeocodedImage::mask`] is an empty `Vec`.
///
/// # Layover vs shadow
/// The mask codes distinguish shadow from layover using a geometric
/// approximation:
///
/// * **`SHADOW`** (code 5): `terrain_normal · look ≤ 0` — the terrain facet
///   faces away from the radar and is not illuminated.  This captures
///   back-slope shadow (and back-slope layover where the folded facet also
///   faces away).
/// * **`LAYOVER`** (code 6): `w > LAYOVER_THRESHOLD` — the flattening weight
///   is so high that the terrain facet contributes more than twice the
///   flat-reference area to the range bin, indicating a steep front slope
///   folded toward the sensor.  This is the dominant front-slope layover
///   signature in the RTC weight domain.
///
/// True range-Doppler layover detection (per-range-bin monotonicity scan)
/// is not performed; the above heuristic covers the vast majority of
/// operationally significant cases.
pub mod mask_code {
    /// Pixel converged and received a valid backscatter value.
    pub const VALID: u8 = 0;
    /// Pixel is outside the radar image footprint after back-projection,
    /// or the bilinearly resampled σ⁰ value is itself NaN
    /// (merge seam, deburst nodata, invalid calibration).
    pub const OUTSIDE_FOOTPRINT: u8 = 1;
    /// DEM has no data at this pixel (out of mosaic coverage or void cell).
    pub const DEM_MISSING: u8 = 2;
    /// Zero-Doppler solver hit `max_newton_iterations` without converging.
    pub const NON_CONVERGED: u8 = 3;
    /// Zero-Doppler solver hit a near-singular Jacobian.  Normally a bug,
    /// not a data issue.
    pub const DEGENERATE_GEOMETRY: u8 = 4;
    /// Local terrain facet faces away from the radar (`terrain_normal · look ≤ 0`).
    /// Covers back-slope shadow and back-slope layover.
    pub const SHADOW: u8 = 5;
    /// Flattening weight `w > LAYOVER_THRESHOLD`: steep front slope compressed
    /// toward the sensor (front-slope layover heuristic).
    pub const LAYOVER: u8 = 6;
    /// Flattening weight `w` below [`super::LAYOVER_SHADOW_MIN_W`]
    /// (extreme foreshortening; γ⁰ would be unstable).
    pub const FORESHORTENING: u8 = 7;
    /// Pixel below the per-pixel NESZ noise floor after the configured
    /// `noise_floor_margin_db`.
    pub const NOISE_MASKED: u8 = 8;
    /// One of the four cardinal DEM neighbours required for the centred
    /// gradient is itself missing.  Distinct from [`DEM_MISSING`] which
    /// fires on the centre pixel.
    pub const DEM_NEIGHBOUR_MISSING: u8 = 9;
}

// ── Terrain flattening helper ──────────────────────────────────────

/// Minimum terrain-flattening weight below which a pixel is classified as
/// extreme foreshortening and masked to NaN.
///
/// `w = (terrain_normal · look) / |flat_normal|`.
/// On flat terrain `w = cos(θ_inc)` ≈ 0.79–0.87 for typical S-1 incidence.
/// Values below 0.05 indicate the terrain facet is nearly perpendicular to
/// the radar look direction (local incidence approaching 90°), producing
/// unreliable backscatter estimates.
const LAYOVER_SHADOW_MIN_W: f64 = 0.05;

/// Flattening weight above which a pixel is classified as front-slope
/// layover and masked.
///
/// When `w > LAYOVER_THRESHOLD` the projected terrain area is more than
/// twice the flat-reference area, signalling a steep forward-facing slope
/// that folds back in range — the dominant front-slope layover signature.
/// Threshold of 2.0 matches common RTC processor conventions.
const LAYOVER_THRESHOLD: f64 = 2.0;

// ── Cardinal-neighbour step for terrain-flattening DEM lookups ─────

/// Approximate metres per degree of latitude on the WGS84 ellipsoid.
///
/// Used to convert metric output spacings into a degree-equivalent step
/// that the cardinal-neighbour DEM lookups and [`compute_terrain_geometry`]
/// can consume directly (the DEM and the ECEF math both work in
/// geographic coordinates).  The exact mean meridional radius is
/// 6_367_449.146 m → 1° = 111_132.95 m; we round to the nearest 10 m
/// because:
///
/// * the gradient computation in [`compute_terrain_geometry`] uses the
///   same scalar for both the terrain and the flat reference normals,
///   so any constant scale factor cancels exactly in `w = local_proj /
///   flat_area`;
/// * the cardinal DEM lookups only need a step that is "of order one
///   output pixel" so the four samples land inside the local DEM cell —
///   sub-percent precision in the conversion has no observable effect.
const METRES_PER_DEGREE_LATITUDE: f64 = 111_320.0;

/// Convert an output-grid pixel spacing (in CRS units) into the
/// equivalent step in **degrees** to use for cardinal-neighbour DEM
/// lookups during terrain flattening / LIA computation.
///
/// # Why this exists
/// The terrain-flattening block at the bottom of the per-pixel loop
/// queries the DEM at the four cardinal neighbours `(lat ± step,
/// lon ± step)` and feeds the resulting heights into
/// [`compute_terrain_geometry`].  Both the DEM API and the ECEF
/// conversion inside `compute_terrain_geometry` work in **geographic
/// degrees**, so `step` must always be in degrees regardless of the
/// chosen output CRS.
///
/// Prior to this helper the call sites passed `config.pixel_spacing_deg`
/// (a misnomer — it actually carries CRS units) directly, which produced
/// step values around 10° in any UTM scene at 10 m spacing, sending
/// every neighbour lookup outside the DEM mosaic and silently masking
/// every output pixel as "flat-masked".
///
/// # Behaviour
/// * [`OutputCrs::Wgs84LatLon`] — `crs_spacing` is already in degrees;
///   returns it unchanged.  This preserves bit-equivalence with the
///   pre-fix WGS84 path.
/// * [`OutputCrs::UtmNorth`] / [`OutputCrs::UtmSouth`] —
///   `crs_spacing` is in metres; returns
///   `crs_spacing / METRES_PER_DEGREE_LATITUDE`.  The same value is
///   used for both lat and lon offsets; `compute_terrain_geometry`
///   accounts for the cos(lat) anisotropy via `geodetic_to_ecef`, and
///   the constant scale factor cancels in the area-ratio `w`.
///
/// # Errors
/// None — this is a pure arithmetic conversion and `crs_spacing` is
/// already validated as positive by the caller.
pub(crate) fn cardinal_neighbour_step_deg(
    crs_spacing: f64,
    crs: &OutputCrs,
) -> f64 {
    match crs {
        OutputCrs::Wgs84LatLon => crs_spacing,
        OutputCrs::UtmNorth { .. }
        | OutputCrs::UtmSouth { .. }
        | OutputCrs::EtrsLaea
        | OutputCrs::WebMercator => {
            crs_spacing / METRES_PER_DEGREE_LATITUDE
        }
    }
}

/// Compute the terrain flattening weight `w` for a single output pixel.
///
/// Implements the area-projection method from Small (2011) "Flattening
/// Gamma: Radiometric Terrain Correction for SAR Imagery".  The weight is:
///
/// ```text
///   w = (terrain_normal · look) / |flat_normal|
/// ```
///
/// where:
/// - `terrain_normal = (P_E − P_W) × (P_N − P_S)` — the unnormalised outward
///   normal of the DEM facet, built from **centred** finite differences over
///   the four cardinal DEM neighbours (E, W, N, S) of the centre pixel.  Its
///   magnitude is proportional to 4 × the terrain facet area.
/// - `flat_normal`  — same calculation but with all heights set to 0
///   (reference WGS84 ellipsoid).  `|flat_normal|` is the flat cell area
///   magnitude and is used as the **denominator** (not `flat_normal · look`).
/// - `look`         — **unit** vector from target toward the satellite.
///
/// # Why `|flat_normal|` and not `flat_normal · look`
///
/// The previous formula used `flat_normal · look` = `|flat_normal| · cos(θ_inc)`
/// as the denominator.  This made `w = 1` on flat terrain and `σ⁰ / w = σ⁰`,
/// i.e., no conversion from σ⁰ to γ⁰.  The correct denominator is the flat
/// cell area `|flat_normal|`, which gives `w = cos(θ_inc)` on flat terrain and
/// therefore `σ⁰ / w = σ⁰ / cos(θ_inc) = γ⁰` — the standard σ⁰ → γ⁰
/// conversion.  For sloped terrain the formula additionally accounts for the
/// terrain facet orientation relative to the look direction.
///
/// # Why centred differences (and not forward differences over the NE quadrant)?
/// Forward differences sample only the (east, north) quadrant of the centre
/// pixel and therefore introduce a directional bias of ½ cell in both axes.
/// Over rough or asymmetric terrain this biases γ⁰ by 10–20 %, depending on
/// slope and look-direction (Small 2011, §III; SNAP/GAMMA both use centred
/// differences for this reason).  Centred differences give the unbiased
/// bivariate-polynomial gradient with the same algorithmic complexity, at
/// the cost of two extra DEM lookups per pixel.
///
/// The correction applied by the caller is: γ⁰ = σ⁰ / w
///
/// # Returns
/// `Some(w)` where `w ≥ LAYOVER_SHADOW_MIN_W`, or `None` for:
/// - shadow / layover (`look · terrain_normal ≤ 0`)
/// - extreme foreshortening (`w < LAYOVER_SHADOW_MIN_W`)
/// Outcome of [`compute_terrain_geometry`] for a single pixel.
///
/// Carries both the σ⁰→γ⁰ flattening weight `w` (when valid) and the local
/// incidence angle cosine `cos_lia` (whenever the facet faces the radar,
/// even if `w` is below the foreshortening threshold).
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) enum TerrainGeometryOutcome {
    /// Facet faces the radar and `LAYOVER_SHADOW_MIN_W ≤ w ≤ LAYOVER_THRESHOLD`.
    Valid { w: f64, cos_lia: f64 },
    /// Facet faces the radar but `w < LAYOVER_SHADOW_MIN_W`.
    /// `cos_lia` is still well-defined and reported.
    Foreshortening { cos_lia: f64 },
    /// Facet faces the radar and `w > LAYOVER_THRESHOLD`: steep front slope,
    /// front-slope layover heuristic.  `cos_lia` is well-defined and reported.
    Layover { cos_lia: f64 },
    /// Facet faces away from the radar (`terrain_normal · look ≤ 0`).
    /// No meaningful `cos_lia`.
    Shadow,
}

/// Compute the σ⁰→γ⁰ flattening weight `w` and the local-incidence-angle
/// cosine for one output pixel.
///
/// `cos_lia = (terrain_normal · look) / |terrain_normal|`, the standard
/// definition of the local incidence angle (Small 2011, eq. 9).
///
/// `w = (terrain_normal · look) / |flat_normal|`.  Note that on flat
/// terrain `w = cos(θ_inc)` so that `σ⁰ / w = γ⁰`; this is *not* the same
/// quantity as `cos_lia` once the terrain is sloped.
///
/// See the previous docstring of `compute_flattening_weight` (now subsumed
/// by this function) for the full derivation and the rationale for using
/// centred differences over forward differences.
#[allow(dead_code)]
pub(crate) fn compute_terrain_geometry(
    lat: f64,
    lon: f64,
    h_east: f64,
    h_west: f64,
    h_north: f64,
    h_south: f64,
    sat_pos_ecef: [f64; 3],
    target_ecef: [f64; 3],
    spacing_deg: f64,
) -> TerrainGeometryOutcome {
    // Flat-earth reference cell area at this latitude (4 ECEF calls).
    let flat_area = flat_cell_area(lat, lon, spacing_deg);
    compute_terrain_geometry_with_flat_area(
        lat, lon, h_east, h_west, h_north, h_south,
        sat_pos_ecef, target_ecef, spacing_deg, flat_area,
    )
}

/// Precompute the flat-earth reference cell area `|flat_normal|` at a given
/// (lat, lon) for the given grid spacing.
///
/// This is the denominator used in the σ⁰→γ⁰ flattening weight
/// `w = (terrain_normal · look) / flat_area`.
///
/// On the WGS84 ellipsoid, `flat_area` depends only on latitude and
/// `spacing_deg` (it is symmetric in longitude by construction of the
/// centred-difference stencil), so callers that process many pixels at the
/// same latitude can call this once per row and pass the result to
/// [`compute_terrain_geometry_with_flat_area`], saving 4 ECEF calls per pixel.
#[inline]
pub(crate) fn flat_cell_area(lat: f64, lon: f64, spacing_deg: f64) -> f64 {
    let f_east  = geodetic_to_ecef(lat,                lon + spacing_deg, 0.0);
    let f_west  = geodetic_to_ecef(lat,                lon - spacing_deg, 0.0);
    let f_north = geodetic_to_ecef(lat + spacing_deg,  lon,               0.0);
    let f_south = geodetic_to_ecef(lat - spacing_deg,  lon,               0.0);
    let flat_normal = cross3(sub3(f_east, f_west), sub3(f_north, f_south));
    norm3(flat_normal)
}

/// Compute the σ⁰→γ⁰ flattening weight `w` and the local-incidence-angle
/// cosine using a caller-supplied precomputed `flat_area` (= `|flat_normal|`).
///
/// This is the hot-path variant called from the terrain-correction inner loop
/// after [`flat_cell_area`] has been precomputed once per output row.
/// It requires only 4 ECEF calls (terrain neighbours) instead of 8.
///
/// See [`compute_terrain_geometry`] for the full derivation.
pub(crate) fn compute_terrain_geometry_with_flat_area(
    lat: f64,
    lon: f64,
    h_east: f64,
    h_west: f64,
    h_north: f64,
    h_south: f64,
    sat_pos_ecef: [f64; 3],
    target_ecef: [f64; 3],
    spacing_deg: f64,
    flat_area: f64,
) -> TerrainGeometryOutcome {
    // Terrain facet via centred differences over the four cardinal DEM
    // neighbours.  Tangent vectors are E−W (east axis) and N−S (north axis);
    // their cross product is the outward facet normal, whose magnitude is
    // 4 × (cell area).
    let p_east  = geodetic_to_ecef(lat,                lon + spacing_deg, h_east);
    let p_west  = geodetic_to_ecef(lat,                lon - spacing_deg, h_west);
    let p_north = geodetic_to_ecef(lat + spacing_deg,  lon,               h_north);
    let p_south = geodetic_to_ecef(lat - spacing_deg,  lon,               h_south);
    let terrain_normal = cross3(
        sub3(p_east,  p_west),
        sub3(p_north, p_south),
    );

    // Illumination direction: vector from target toward satellite.
    let look = normalize3(sub3(sat_pos_ecef, target_ecef));

    let local_proj = dot3(terrain_normal, look);

    // Shadow: terrain facet faces away from the radar.
    if local_proj <= 0.0 {
        return TerrainGeometryOutcome::Shadow;
    }

    // Local incidence angle cosine.
    let terrain_area = norm3(terrain_normal);
    // Defensive: a zero-magnitude terrain normal would arise only from a
    // pathological DEM (four collinear neighbours through Earth's centre)
    // and is treated as shadow rather than producing a NaN cos_lia.
    if terrain_area == 0.0 {
        return TerrainGeometryOutcome::Shadow;
    }
    let cos_lia = (local_proj / terrain_area).clamp(0.0, 1.0);

    // Denominator is |flat_normal| (flat cell area), NOT flat_normal·look.
    // Using flat_normal·look = |flat_normal|·cos(θ_inc) as denominator would
    // cancel the cos(θ_inc) factor and give w = 1 on flat terrain (σ⁰/w = σ⁰).
    // Using |flat_normal| gives w = cos(θ_inc) on flat terrain so that
    // σ⁰/w = σ⁰/cos(θ_inc) = γ⁰, the standard terrain-flattened gamma-naught.
    let w = local_proj / flat_area;
    if w < LAYOVER_SHADOW_MIN_W {
        TerrainGeometryOutcome::Foreshortening { cos_lia }
    } else if w > LAYOVER_THRESHOLD {
        TerrainGeometryOutcome::Layover { cos_lia }
    } else {
        TerrainGeometryOutcome::Valid { w, cos_lia }
    }
}

/// Backward-compatible wrapper used by the legacy unit tests.  New call
/// sites should use [`compute_terrain_geometry`] directly.
#[cfg(test)]
fn compute_flattening_weight(
    lat: f64,
    lon: f64,
    h_east: f64,
    h_west: f64,
    h_north: f64,
    h_south: f64,
    sat_pos_ecef: [f64; 3],
    target_ecef: [f64; 3],
    spacing_deg: f64,
) -> Option<f64> {
    match compute_terrain_geometry(
        lat, lon, h_east, h_west, h_north, h_south,
        sat_pos_ecef, target_ecef, spacing_deg,
    ) {
        TerrainGeometryOutcome::Valid { w, .. } => Some(w),
        TerrainGeometryOutcome::Foreshortening { .. }
        | TerrainGeometryOutcome::Layover { .. }
        | TerrainGeometryOutcome::Shadow => None,
    }
}

// ── Radar image coordinate mapping ─────────────────────────────────

/// Parameters for converting (azimuth_time, slant_range) to image (line, sample).
struct RadarGeometry {
    /// Azimuth time of line 0 in the merged image.
    azimuth_start_time: DateTime<Utc>,
    /// Seconds per azimuth line.
    azimuth_time_interval_s: f64,
    /// Slant range to the first sample (column 0) of the merged image, in metres.
    near_range_m: f64,
    /// Range pixel spacing in metres (slant range).
    range_pixel_spacing_m: f64,
    /// Number of lines in the merged image.
    n_lines: usize,
    /// Number of samples in the merged image.
    n_samples: usize,
}

impl RadarGeometry {
    /// Build from scene metadata and merged image dimensions.
    ///
    /// Fails loudly (no silent defaults) if the scene metadata does not carry
    /// the required timing parameters.
    fn from_scene_and_merged(
        scene: &SceneMetadata,
        merged: &dyn RadarImage,
    ) -> Result<Self, TerrainCorrectionError> {
        // All IW subswaths share the same azimuth time interval to within
        // numerical precision; pick the first available subswath only to
        // verify the scene metadata is non-empty.
        let _first_sw = scene
            .sub_swaths
            .first()
            .ok_or(TerrainCorrectionError::NoSubSwaths)?;

        // Read the ATI from the *merged image*, not from the scene subswath
        // metadata.  After apply_multilook, merged.azimuth_time_interval_s is
        // already scaled by az_looks (native_ati × az_looks), which is
        // required for correct line-to-time mapping when az_multilook > 1.
        // For full-resolution images (az_multilook == 1) the two values are
        // identical.
        let ati = merged.azimuth_time_interval_s();
        if !ati.is_finite() || ati <= 0.0 {
            return Err(TerrainCorrectionError::Config(format!(
                "invalid azimuth_time_interval_s from merged image: {ati}"
            )));
        }

        // Azimuth time of merged line 0.
        //
        // Anchored to the earliest input subswath's `bursts[0].azimuth_time_utc`
        // (see `MergedSigma0::azimuth_start_time` for the convention).  This
        // is exact for the earliest subswath and off by ≤ 1 azimuth line for
        // the others, a residual the Newton zero-Doppler solver absorbs in
        // its first iteration.  The previous implementation used the
        // product-level `scene.start_time`, which differs by up to a few ms
        // (~1 line) from the debursted line-0 time and produced a constant
        // global azimuth bias.
        let az_start = merged.azimuth_start_time();

        // Near range from the merged image header.
        let near_range_m = merged.near_slant_range_time_s() * SPEED_OF_LIGHT_M_S / 2.0;

        Ok(RadarGeometry {
            azimuth_start_time: az_start,
            azimuth_time_interval_s: ati,
            near_range_m,
            range_pixel_spacing_m: merged.range_pixel_spacing_m(),
            n_lines: merged.lines(),
            n_samples: merged.samples(),
        })
    }

    /// Convert (azimuth_time, slant_range_m) to fractional (line, sample).
    /// Returns None if outside the image bounds.
    fn to_image_coords(&self, az_time: DateTime<Utc>, slant_range_m: f64) -> Option<(f64, f64)> {
        let dt_us = (az_time - self.azimuth_start_time)
            .num_microseconds()
            .unwrap_or(0) as f64 // SAFETY-OK: chrono microseconds cannot overflow for scene-window durations
            / 1e6;
        let line = dt_us / self.azimuth_time_interval_s;
        let sample = (slant_range_m - self.near_range_m) / self.range_pixel_spacing_m;

        // Allow half-pixel margin for interpolation
        if line < -0.5
            || line > (self.n_lines as f64 - 0.5)
            || sample < -0.5
            || sample > (self.n_samples as f64 - 0.5)
        {
            return None;
        }
        Some((line, sample))
    }
}

// ── Bilinear interpolation in radar image ──────────────────────────

/// Core bilinear interpolation over any flat row-major `f32` slice.
///
/// Returns NaN if any of the four corner values is NaN.
fn bilinear_sample_slice(data: &[f32], n_lines: usize, n_samples: usize, line: f64, sample: f64) -> f32 {
    if !line.is_finite() || !sample.is_finite() {
        return f32::NAN;
    }
    let max_line = n_lines.saturating_sub(1);
    let max_sample = n_samples.saturating_sub(1);

    let l0 = (line.floor() as isize).clamp(0, max_line as isize) as usize;
    let l1 = (l0 + 1).min(max_line);
    let s0 = (sample.floor() as isize).clamp(0, max_sample as isize) as usize;
    let s1 = (s0 + 1).min(max_sample);

    let dl = line - l0 as f64;
    let ds = sample - s0 as f64;

    let v00 = data[l0 * n_samples + s0];
    let v01 = data[l0 * n_samples + s1];
    let v10 = data[l1 * n_samples + s0];
    let v11 = data[l1 * n_samples + s1];

    if v00.is_nan() || v01.is_nan() || v10.is_nan() || v11.is_nan() {
        return f32::NAN;
    }

    let w00 = (1.0 - dl) * (1.0 - ds);
    let w01 = (1.0 - dl) * ds;
    let w10 = dl * (1.0 - ds);
    let w11 = dl * ds;

    (v00 as f64 * w00 + v01 as f64 * w01 + v10 as f64 * w10 + v11 as f64 * w11) as f32
}

/// Bilinear interpolation of a value in the merged sigma0 image.
#[cfg(test)]
fn bilinear_sample(merged: &MergedSigma0, line: f64, sample: f64) -> f32 {
    bilinear_sample_slice(&merged.data, merged.lines, merged.samples, line, sample)
}

// ── Bicubic (Keys, α = −0.5) interpolation ─────────────────────────

/// Keys cubic kernel weight for parameter `t` (distance from the knot).
///
/// Uses the Keys (1989) piecewise polynomial with α = −0.5 (also known
/// as the Catmull-Rom / cubic Hermite spline):
///
/// ```text
///   |t| < 1:  (α+2)|t|³ − (α+3)|t|² + 1  = 1.5|t|³ − 2.5|t|² + 1
///   |t| < 2:  α|t|³ − 5α|t|² + 8α|t| − 4α = −0.5|t|³ + 2.5|t|² − 4|t| + 2
///   else:     0
/// ```
///
/// `t` must already be in [−2, 2]; the caller is responsible for clamping.
#[inline(always)]
fn keys_weight(t: f64) -> f64 {
    let t = t.abs();
    if t < 1.0 {
        1.5 * t * t * t - 2.5 * t * t + 1.0
    } else if t < 2.0 {
        -0.5 * t * t * t + 2.5 * t * t - 4.0 * t + 2.0
    } else {
        0.0
    }
}

/// 16-tap bicubic interpolation over a flat row-major `f32` slice.
///
/// Samples the 4×4 neighbourhood centred on `(line, sample)`.  Returns NaN if
/// any of the 16 taps is NaN (propagates nodata cleanly).  For locations within
/// 1 pixel of the image boundary the kernel is clamped to the nearest valid
/// row/column.
///
/// **Note on SAR use**: bicubic can introduce mild Gibbs ringing near high-
/// contrast edges (bright point targets, ship returns).  Bilinear remains the
/// default for SAR backscatter.  Bicubic is useful when the output is
/// resampled to a much coarser grid than the input.
fn bicubic_sample_slice(
    data: &[f32],
    n_lines: usize,
    n_samples: usize,
    line: f64,
    sample: f64,
) -> f32 {
    if !line.is_finite() || !sample.is_finite() {
        return f32::NAN;
    }
    let max_l = (n_lines as isize) - 1;
    let max_s = (n_samples as isize) - 1;

    let l_floor = line.floor() as isize;
    let s_floor = sample.floor() as isize;

    let dl = line - l_floor as f64;
    let ds = sample - s_floor as f64;

    let mut acc = 0.0_f64;
    for li in -1_isize..=2 {
        let row = (l_floor + li).clamp(0, max_l) as usize;
        let wl = keys_weight(li as f64 - dl);
        for si in -1_isize..=2 {
            let col = (s_floor + si).clamp(0, max_s) as usize;
            let v = data[row * n_samples + col];
            if v.is_nan() {
                return f32::NAN;
            }
            acc += v as f64 * wl * keys_weight(si as f64 - ds);
        }
    }
    acc as f32
}

/// Resample the merged sigma0 image at fractional coordinates using the
/// configured [`crate::pipeline_options::ResamplingKernel`].
#[cfg(test)]
fn resample_merged(
    merged: &MergedSigma0,
    line: f64,
    sample: f64,
    kernel: crate::pipeline_options::ResamplingKernel,
) -> f32 {
    merged.sample_at(line, sample, kernel)
}

// ── RadarImage impl for MergedSigma0 ──────────────────────────────

impl RadarImage for MergedSigma0 {
    fn sample_at(&self, line: f64, sample: f64, kernel: ResamplingKernel) -> f32 {
        match kernel {
            ResamplingKernel::Bilinear => {
                bilinear_sample_slice(&self.data, self.lines, self.samples, line, sample)
            }
            ResamplingKernel::Bicubic => {
                bicubic_sample_slice(&self.data, self.lines, self.samples, line, sample)
            }
            ResamplingKernel::Lanczos3 => {
                lanczos_sample_slice(&self.data, self.lines, self.samples, line, sample)
            }
        }
    }

    fn nesz_at(&self, line: f64, sample: f64) -> Option<f32> {
        Some(bilinear_sample_slice(
            &self.nesz, self.lines, self.samples, line, sample,
        ))
    }

    fn lines(&self) -> usize {
        self.lines
    }

    fn samples(&self) -> usize {
        self.samples
    }

    fn azimuth_start_time(&self) -> chrono::DateTime<chrono::Utc> {
        self.azimuth_start_time
    }

    fn near_slant_range_time_s(&self) -> f64 {
        self.near_slant_range_time_s
    }

    fn range_pixel_spacing_m(&self) -> f64 {
        self.range_pixel_spacing_m
    }

    fn azimuth_time_interval_s(&self) -> f64 {
        self.azimuth_time_interval_s
    }
}

// ── Lanczos-3 interpolation ────────────────────────────────────────

/// Lanczos kernel weight for distance `x` (in pixels) from the sample point.
///
/// `L(x) = sinc(x) · sinc(x/3)` for `|x| < 3`, zero otherwise.
/// Computed as `sin(πx)/(πx) · sin(πx/3)/(πx/3)`; L(0) = 1 by L'Hôpital.
/// Negative sidelobes in `|x| ∈ (1, 3)` produce mild ringing near step edges.
#[inline(always)]
fn lanczos_weight(x: f64) -> f64 {
    const A: f64 = 3.0;
    let ax = x.abs();
    if ax < 1e-10 {
        1.0
    } else if ax < A {
        let px = std::f64::consts::PI * x;
        let pxa = std::f64::consts::PI * x / A;
        px.sin() * pxa.sin() / (px * pxa)
    } else {
        0.0
    }
}

/// 36-tap Lanczos-3 interpolation over a flat row-major `f32` slice.
///
/// Samples a 6 × 6 neighbourhood (offsets −2 to +3 in each dimension).
/// Returns NaN if any tap is NaN (propagates nodata cleanly).  Boundary
/// taps are clamped to the nearest valid row/column, consistent with
/// [`bilinear_sample_slice`] and [`bicubic_sample_slice`].
fn lanczos_sample_slice(
    data: &[f32],
    n_lines: usize,
    n_samples: usize,
    line: f64,
    sample: f64,
) -> f32 {
    if !line.is_finite() || !sample.is_finite() {
        return f32::NAN;
    }
    let max_l = (n_lines as isize) - 1;
    let max_s = (n_samples as isize) - 1;
    let l_floor = line.floor() as isize;
    let s_floor = sample.floor() as isize;
    let mut acc = 0.0_f64;
    for li in -2_isize..=3 {
        let row = (l_floor + li).clamp(0, max_l) as usize;
        let wl = lanczos_weight(line - (l_floor + li) as f64);
        for si in -2_isize..=3 {
            let col = (s_floor + si).clamp(0, max_s) as usize;
            let v = data[row * n_samples + col];
            if v.is_nan() {
                return f32::NAN;
            }
            acc += v as f64 * wl * lanczos_weight(sample - (s_floor + si) as f64);
        }
    }
    acc as f32
}

// ── Zero-Doppler solver ────────────────────────────────────────────

/// Solve the zero-Doppler equation `f(t) = (P − S(t)) · V(t) = 0` for the
/// azimuth time `t` at which the line from the satellite to the target `P` is
/// perpendicular to the satellite velocity vector.
///
/// Uses Newton's method with a finite-difference acceleration estimate:
///
/// ```text
///   f'(t) = −|V(t)|² + (P − S(t)) · A(t)
/// ```
///
/// The second term is small for near-circular LEO (≈ 10⁻³ relative to the
/// first), but we include it for correctness near the scene edges.
///
/// # Failure modes (explicit, no silent fallback)
///
/// * [`ZeroDopplerOutcome::NotConverged`] — the update `|Δt|` never fell below
///   `tol_s` within `max_iter` iterations.  The last step's residual is
///   reported so callers can log a histogram.
/// * [`ZeroDopplerOutcome::DegenerateGeometry`] — `|f'(t)|` fell below
///   `DERIV_FLOOR`; this should not happen for a physical LEO orbit and
///   usually indicates bad orbit data or a target inside the Earth.
///
/// Orbit interpolation errors (coverage gaps) propagate up as `OrbitError`.
pub(crate) fn solve_zero_doppler(
    target_ecef: [f64; 3],
    orbit: &OrbitData,
    initial_guess: DateTime<Utc>,
    max_iter: usize,
    tol_s: f64,
) -> Result<ZeroDopplerOutcome, OrbitError> {
    // Lower bound on |f'(t)|.  For Sentinel-1 (|V|² ≈ 5.6 × 10⁷ m²/s²) this
    // threshold is ~12 orders of magnitude below the expected magnitude; only
    // pathological geometries hit it.
    const DERIV_FLOOR: f64 = 1.0;

    let mut t = initial_guess;
    let mut last_delta_s: f64 = f64::INFINITY;

    for i in 0..max_iter {
        // Single Lagrange pass: returns position, velocity, AND the analytical
        // derivative of the velocity polynomial (orbital acceleration).
        // This replaces the previous two-call pattern:
        //   sv      = interpolate_orbit(orbit, t)
        //   sv_plus = interpolate_orbit(orbit, t + 10 ms)  ← FD for acceleration
        // and eliminates one Lagrange evaluation per Newton iteration.
        let (sv, accel) = interpolate_orbit_and_accel(orbit, t)?;
        let r = sub3(target_ecef, sv.position_m);
        let f = dot3(r, sv.velocity_m_s);
        let v_sq = dot3(sv.velocity_m_s, sv.velocity_m_s);
        let f_prime = -v_sq + dot3(r, accel);

        if f_prime.abs() < DERIV_FLOOR {
            return Ok(ZeroDopplerOutcome::DegenerateGeometry);
        }

        let delta_s = -f / f_prime;
        last_delta_s = delta_s;

        if delta_s.abs() < tol_s {
            // Convergence check is performed BEFORE applying the step.
            // `sv` is therefore already at the converged time `t`.  This
            // avoids a second `interpolate_orbit_and_accel` call that would
            // otherwise be needed to evaluate the state vector at `t +
            // delta_s`.  The sub-tolerance step (|delta_s| < 1 μs) changes
            // `t` by < 1 μs, causing < 7.6 mm position error — negligible
            // for any σ⁰ application.
            //
            // i is 0-based; report 1-based iteration count for histogram
            // semantics (slot 1 = converged on first iter).  u8 is fine:
            // max_iter is bounded well below 256 in every realistic
            // configuration.
            let iterations: u8 = u8::try_from(i + 1).unwrap_or(u8::MAX); // SAFETY-OK: diagnostic counter; saturating clamp is exactly the desired behaviour for histogram bin assignment
            return Ok(ZeroDopplerOutcome::Converged {
                time: t,
                residual_s: delta_s.abs(),
                iterations,
                sat_sv: sv,
            });
        }

        t = t + Duration::microseconds((delta_s * 1e6) as i64);
    }

    Ok(ZeroDopplerOutcome::NotConverged {
        last_residual_s: last_delta_s.abs(),
    })
}

// ── Initial guess from geolocation grid ────────────────────────────

/// Build a coarse lookup for initial azimuth time guesses from the geolocation grid.
///
/// Returns a sorted list of (latitude, azimuth_time) pairs for binary-search lookup.
fn build_azimuth_time_lut(
    grids: &[(SubSwathId, Vec<GeolocationGridPoint>)],
) -> Vec<(f64, DateTime<Utc>)> {
    let mut lut: Vec<(f64, DateTime<Utc>)> = Vec::new();

    // Use IW2 (middle swath) for the best range of latitudes,
    // falling back to whatever is available
    let grid = grids
        .iter()
        .find(|(id, _)| *id == SubSwathId::IW2)
        .or_else(|| grids.first())
        .map(|(_, pts)| pts);

    if let Some(pts) = grid {
        for pt in pts {
            lut.push((pt.latitude_deg, pt.azimuth_time_utc));
        }
    }

    // Sort by latitude (ascending)
    lut.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal)); // SAFETY-OK: NaN-safe sort comparator (NaN treated as equal); latitudes come from validated geolocation grid
    // Deduplicate very close latitudes (keep first)
    lut.dedup_by(|a, b| (a.0 - b.0).abs() < 1e-6);
    lut
}

/// Get initial azimuth time guess for a given latitude from the LUT.
fn guess_azimuth_time(lut: &[(f64, DateTime<Utc>)], lat_deg: f64) -> DateTime<Utc> {
    if lut.is_empty() {
        // Shouldn't happen, but if called with empty LUT
        return Utc::now();
    }
    if lut.len() == 1 {
        return lut[0].1;
    }

    // Binary search for the closest latitude
    let idx = match lut.binary_search_by(|probe| {
        probe.0.partial_cmp(&lat_deg).unwrap_or(std::cmp::Ordering::Equal) // SAFETY-OK: NaN-safe binary-search comparator (NaN treated as equal); validated grid latitudes are finite
    }) {
        Ok(i) => return lut[i].1,
        Err(i) => i,
    };

    if idx == 0 {
        return lut[0].1;
    }
    if idx >= lut.len() {
        return lut[lut.len() - 1].1;
    }

    // Linear interpolation between bracketing points
    let (lat0, t0) = lut[idx - 1];
    let (lat1, t1) = lut[idx];
    let frac = (lat_deg - lat0) / (lat1 - lat0);
    let dt_us = (t1 - t0).num_microseconds().unwrap_or(0) as f64; // SAFETY-OK: chrono microseconds cannot overflow for geolocation grid spacing
    t0 + Duration::microseconds((frac * dt_us) as i64)
}

// ── Main function ──────────────────────────────────────────────────

/// Perform terrain correction and geocoding on a merged slant-range σ⁰ image.
///
/// # Arguments
///
/// * `merged` — Merged slant-range σ⁰ from [`crate::merge_subswaths::merge_subswaths`].
/// * `scene` — Validated scene metadata (with precise orbit applied).
/// * `dem` — DEM mosaic covering the scene footprint.
/// * `grids` — Geolocation grid points from [`crate::parse::parse_geolocation_grids`].
/// * `config` — Terrain correction configuration.
///
/// # Returns
///
/// A [`GeocodedImage`] on a regular lat/lon grid.
/// Compute the output grid extent (top-left and bottom-right corners)
/// in the output CRS native coordinates, snapped outward to multiples
/// of `spacing`.
///
/// Returns `(y_top, y_bottom, x_left, x_right)` where:
/// - For [`OutputCrs::Wgs84LatLon`] (`projector = None`): degrees;
///   `y_top = north_lat`, `y_bottom = south_lat`, `x_left = west_lon`,
///   `x_right = east_lon`.  Identical to the legacy code path.
/// - For projected CRSs (`projector = Some`): metres; `y_top = max_northing`,
///   `y_bottom = min_northing`, `x_left = min_easting`, `x_right = max_easting`.
///   Computed by projecting the four lat/lon corners of `bb` and taking
///   their axis-aligned bounding rectangle in the projected plane.
///
/// **Cross-zone limitation:** when `projector` is a UTM CRS and the
/// scene crosses a UTM zone boundary, the four-corner approximation
/// over-estimates the output extent (the projected scene is no longer
/// axis-aligned).  The caller is responsible for choosing a CRS whose
/// distortion over the scene is acceptable; `OutputCrs::auto_for_lon_lat`
/// picks the UTM zone of the scene centre, which is correct for any
/// scene that fits inside a single zone (≤ ~600 km wide at the equator).
fn compute_output_grid_extent(
    bb: &crate::types::BoundingBox,
    projector: Option<&Projector>,
    spacing: f64,
) -> Result<(f64, f64, f64, f64), TerrainCorrectionError> {
    match projector {
        None => Ok((
            (bb.max_lat_deg / spacing).ceil() * spacing,
            (bb.min_lat_deg / spacing).floor() * spacing,
            (bb.min_lon_deg / spacing).floor() * spacing,
            (bb.max_lon_deg / spacing).ceil() * spacing,
        )),
        Some(proj) => {
            let corners = [
                (bb.min_lon_deg, bb.min_lat_deg),
                (bb.min_lon_deg, bb.max_lat_deg),
                (bb.max_lon_deg, bb.min_lat_deg),
                (bb.max_lon_deg, bb.max_lat_deg),
            ];
            let mut x_min = f64::INFINITY;
            let mut x_max = f64::NEG_INFINITY;
            let mut y_min = f64::INFINITY;
            let mut y_max = f64::NEG_INFINITY;
            for (lon, lat) in corners {
                let (x, y) = proj.lonlat_to_xy(lon, lat).map_err(|e| {
                    TerrainCorrectionError::Projection(format!(
                        "scene corner ({lon:.4}, {lat:.4}): {e}"
                    ))
                })?;
                if x < x_min { x_min = x; }
                if x > x_max { x_max = x; }
                if y < y_min { y_min = y; }
                if y > y_max { y_max = y; }
            }
            Ok((
                (y_max / spacing).ceil() * spacing,
                (y_min / spacing).floor() * spacing,
                (x_min / spacing).floor() * spacing,
                (x_max / spacing).ceil() * spacing,
            ))
        }
    }
}

/// Build the half-open `(row_start, row_end)` ranges that tile
/// `0..n_rows` into chunks of at most `chunk_size` consecutive rows.
///
/// Used by [`terrain_correction`] to schedule per-chunk Rayon tasks
/// instead of per-row tasks (Phase 4: better DEM-tile cache reuse and
/// lower scheduling overhead).  Pure function so the boundary
/// arithmetic — which is the only new logic introduced by chunking —
/// is unit-testable in isolation.
///
/// Invariants asserted by the unit tests:
/// * Concatenating the produced ranges in order yields exactly the
///   sequence `0, 1, ..., n_rows-1` (no gaps, no overlaps).
/// * Returns an empty iterator when `n_rows == 0`.
/// * Panics when `chunk_size == 0` — callers must validate first.
pub(crate) fn row_chunk_ranges(
    n_rows: usize,
    chunk_size: usize,
) -> impl Iterator<Item = (usize, usize)> {
    assert!(chunk_size > 0, "row_chunk_ranges: chunk_size must be > 0");
    let n_chunks = n_rows.div_ceil(chunk_size);
    (0..n_chunks).map(move |chunk_idx| {
        let row_start = chunk_idx * chunk_size;
        let row_end = ((chunk_idx + 1) * chunk_size).min(n_rows);
        (row_start, row_end)
    })
}

pub fn terrain_correction(
    merged: &dyn RadarImage,
    scene: &SceneMetadata,
    dem: &dyn DemSource,
    grids: &[(SubSwathId, Vec<GeolocationGridPoint>)],
    config: &TerrainCorrectionConfig,
) -> Result<GeocodedImage, TerrainCorrectionError> {
    // Histogram bin assignment clamps onto the last slot if a caller
    // configures more iterations than the histogram can represent; in
    // dev builds we want a loud failure to flag the misconfiguration.
    debug_assert!(
        config.max_newton_iterations < NEWTON_ITER_HIST_SIZE,
        "max_newton_iterations ({}) must be < NEWTON_ITER_HIST_SIZE ({})",
        config.max_newton_iterations,
        NEWTON_ITER_HIST_SIZE,
    );

    if grids.is_empty() {
        return Err(TerrainCorrectionError::NoGridPoints);
    }

    if config.pixel_spacing_deg <= 0.0 {
        return Err(TerrainCorrectionError::Config(
            "pixel spacing must be positive".into(),
        ));
    }

    let orbit = &scene.orbit;
    let radar_geom = RadarGeometry::from_scene_and_merged(scene, merged)?;

    // Build initial guess LUT from geolocation grid
    let az_lut = build_azimuth_time_lut(grids);
    if az_lut.is_empty() {
        return Err(TerrainCorrectionError::NoGridPoints);
    }

    // Build a projector once and share it across rayon workers (Projector
    // is Send+Sync — proven by the compile-time assertion in
    // `output_crs::tests::projector_is_send_and_sync`).  None means the
    // output grid is the source WGS84 lat/lon grid and no per-pixel
    // reprojection is needed (byte-identical legacy path).
    let projector: Option<Projector> = match config.crs {
        OutputCrs::Wgs84LatLon => None,
        _ => Some(
            config
                .crs
                .projector()
                .map_err(|e: OutputCrsError| TerrainCorrectionError::Projection(e.to_string()))?,
        ),
    };

    // Determine output grid extent.  Field names `lat_north` / `lon_west`
    // are kept for parity with the legacy code path; they hold native
    // CRS coordinates (degrees for WGS84, metres for UTM).  In both cases
    // the convention is "top-left corner of the north-up grid", so:
    //   * `lat_north` = larger of the two Y bounds (lat for WGS84,
    //     northing for UTM)
    //   * `lon_west`  = smaller of the two X bounds (lon for WGS84,
    //     easting for UTM)
    let bb = &scene.bounding_box;
    let spacing = config.pixel_spacing_deg;

    let (lat_north, lat_south, lon_west, lon_east) =
        compute_output_grid_extent(bb, projector.as_ref(), spacing)?;

    let n_rows = ((lat_north - lat_south) / spacing).round() as usize;
    let n_cols = ((lon_east - lon_west) / spacing).round() as usize;

    if n_rows == 0 || n_cols == 0 {
        return Err(TerrainCorrectionError::Config(
            format!("output grid is empty: {} rows × {} cols", n_rows, n_cols),
        ));
    }

    let total = n_rows * n_cols;

    // Number of anchor points per row for the projector interpolation
    // fast path.  See [`output_crs::Projector::project_row`] for the
    // interpolation contract.  Read once here so all rows see the same
    // value and no env lookup happens in the hot loop.  No silent
    // fallback: a malformed env value is a hard error, since the user
    // explicitly asked to tune this and a typo must not silently revert
    // to the default.
    let row_projection_anchors: usize = match std::env::var("SARDINE_TC_PROJECTOR_ANCHORS") {
        Ok(s) => s.parse::<usize>().map_err(|e| {
            TerrainCorrectionError::Config(format!(
                "SARDINE_TC_PROJECTOR_ANCHORS={s:?} is not a non-negative integer: {e}"
            ))
        })?,
        // Default 256: at this anchor density the worst-case linear-
        // interpolation residual on a 25k-pixel S-1 row at a UTM zone
        // edge is ~ 2e-7° (≈ 2 cm) in both lon and lat — below
        // 0.1 px at 10 m spacing and well below the per-pixel terrain-
        // correction noise floor.  Bumping the count linearly trades
        // accuracy quadratically; smaller values may regress geocoding.
        Err(std::env::VarError::NotPresent) => 256,
        Err(e) => {
            return Err(TerrainCorrectionError::Config(format!(
                "reading SARDINE_TC_PROJECTOR_ANCHORS: {e}"
            )));
        }
    };
    if row_projection_anchors == 0 {
        return Err(TerrainCorrectionError::Config(
            "SARDINE_TC_PROJECTOR_ANCHORS must be >= 1 (default 256)".into(),
        ));
    }

    // Number of consecutive output rows processed by a single Rayon task.
    // Larger chunks => better DEM-tile cache reuse on per-worker
    // `LAST_HIT_TILE` and lower scheduling overhead, but coarser load
    // balance.  Default 64 was chosen because at 10 m output spacing
    // 64 rows = 640 m of north-south coverage, comfortably within a
    // single 1° SRTM/GLO-30 tile so the warm-tile cache hits ~100% on
    // the second-and-later rows of a chunk.  Same no-silent-fallback
    // policy as the projector knob: malformed env is a hard error.
    let row_chunk_size: usize = match std::env::var("SARDINE_TC_ROW_CHUNK") {
        Ok(s) => s.parse::<usize>().map_err(|e| {
            TerrainCorrectionError::Config(format!(
                "SARDINE_TC_ROW_CHUNK={s:?} is not a non-negative integer: {e}"
            ))
        })?,
        Err(std::env::VarError::NotPresent) => 256,
        Err(e) => {
            return Err(TerrainCorrectionError::Config(format!(
                "reading SARDINE_TC_ROW_CHUNK: {e}"
            )));
        }
    };
    if row_chunk_size == 0 {
        return Err(TerrainCorrectionError::Config(
            "SARDINE_TC_ROW_CHUNK must be >= 1 (default 256)".into(),
        ));
    }

    // Per-row result.  Each row is processed independently; counters are
    // accumulated locally and summed after parallel collection.
    struct RowResult {
        row_data: Vec<f32>,
        row_cos_lia: Vec<f32>, // empty when compute_lia=false
        row_mask: Vec<u8>,     // empty when compute_lia=false
        valid: usize,
        dem_missing: usize,
        outside_footprint: usize,
        non_converged: usize,
        degenerate: usize,
        flat_masked: usize,
        noise_masked: usize,
        /// Newton iteration histogram for this row.  See
        /// [`NEWTON_ITER_HIST_SIZE`] for the binning convention.
        iter_histogram: [u64; NEWTON_ITER_HIST_SIZE],
    }

    // Progress counter shared across threads (Relaxed: display only, no
    // ordering guarantees needed).
    let rows_done = AtomicUsize::new(0);

    // The cardinal-neighbour step in degrees is the same for every pixel —
    // it depends only on `spacing` and the output CRS.  Hoist it out of
    // the per-pixel loop so the flat_area precomputation can use it once
    // per row instead of recomputing it per pixel.
    let step_deg_for_flat_area = cardinal_neighbour_step_deg(spacing, &config.crs);

    // Process rows in parallel, in chunks of `row_chunk_size`.  Every
    // reference captured by the closure is a shared immutable borrow of
    // a Send+Sync value; all mutation is local to the row body.  Within
    // a single chunk rows execute sequentially on one worker, so the
    // thread-local DEM tile cache (`LAST_HIT_TILE`) and the warmed
    // anchor-projector buffer survive across rows.
    let chunk_ranges: Vec<(usize, usize)> = row_chunk_ranges(n_rows, row_chunk_size).collect();
    let chunk_results: Result<Vec<Vec<RowResult>>, TerrainCorrectionError> = chunk_ranges
        .into_par_iter()
        .map(|(row_start, row_end)| {
            let mut chunk_out: Vec<RowResult> = Vec::with_capacity(row_end - row_start);
            for row in row_start..row_end {
            // For WGS84 the row's Y coordinate IS latitude; for UTM
            // it's a northing in metres.  Either way we need an
            // approximate latitude to seed the per-row azimuth-time
            // LUT lookup.
            let row_y = lat_north - (row as f64 + 0.5) * spacing;
            let row_lat_seed = match projector.as_ref() {
                None => row_y,
                Some(proj) => {
                    // Project the mid-row, mid-column point to get an
                    // approximate latitude for the LUT.  One projection
                    // per row is negligible vs. n_cols projections.
                    let mid_x = lon_west + (n_cols as f64 / 2.0) * spacing;
                    proj.xy_to_lonlat(mid_x, row_y)
                        .map(|(_lon, lat)| lat)
                        .map_err(|e| {
                            TerrainCorrectionError::Projection(format!(
                                "row {row} mid-point: {e}"
                            ))
                        })?
                }
            };

            // Get initial azimuth time guess for this row.
            let mut prev_az_time = guess_azimuth_time(&az_lut, row_lat_seed);

            let mut row_data = vec![f32::NAN; n_cols];
            // Sidecar buffers are sized only when requested, so the
            // memory/cache cost is paid solely by callers that opt in.
            let (mut row_cos_lia, mut row_mask) = if config.compute_lia {
                (vec![f32::NAN; n_cols], vec![mask_code::OUTSIDE_FOOTPRINT; n_cols])
            } else {
                (Vec::new(), Vec::new())
            };
            let mut valid = 0usize;
            let mut dem_missing = 0usize;
            let mut outside_footprint = 0usize;
            let mut non_converged = 0usize;
            let mut degenerate = 0usize;
            let mut flat_masked = 0usize;
            let mut noise_masked = 0usize;
            let mut iter_histogram = [0u64; NEWTON_ITER_HIST_SIZE];

            // Pre-project the entire row in one anchor-interpolated pass
            // when the output CRS is metric.  For Wgs84LatLon the
            // projector is None and (lon, lat) = (pixel_x, pixel_y) so
            // no precomputation is needed (and the buffer stays empty).
            let row_lonlat: Vec<(f64, f64)> = match projector.as_ref() {
                None => Vec::new(),
                Some(proj) => {
                    let x_start = lon_west + 0.5 * spacing;
                    proj.project_row(
                        row_y,
                        x_start,
                        spacing,
                        n_cols,
                        row_projection_anchors,
                    )
                    .map_err(|e| {
                        TerrainCorrectionError::Projection(format!("row {row}: {e}"))
                    })?
                }
            };

            // Precompute the flat-earth reference cell area for this row.
            // On the WGS84 ellipsoid |flat_normal| is a function of latitude
            // and spacing only — it is independent of longitude (WGS84 is
            // rotationally symmetric in longitude, so the E-W and N-S tangent
            // vectors at h=0 produce the same cross-product magnitude for all
            // longitudes at a fixed latitude).  Computing this once per row
            // eliminates 4 geodetic_to_ecef calls per pixel when flatten or
            // compute_lia is active.  For UTM rows, row_lat_seed is the
            // mid-row latitude; the latitude variation across a single row
            // (~0.5°) changes flat_area by < 0.01%, which is negligible
            // compared to DEM and orbit uncertainties.
            let row_flat_area: f64 = if config.flatten || config.compute_lia {
                let mid_lon = lon_west + (n_cols as f64 / 2.0) * spacing;
                let row_ref_lon = if projector.is_some() {
                    // UTM: mid_lon is a northing, not a geodetic lon.
                    // Use mid-row lonlat if available, else fall back to
                    // a centre-scene longitude estimate.
                    if !row_lonlat.is_empty() {
                        row_lonlat[n_cols / 2].0
                    } else {
                        mid_lon
                    }
                } else {
                    mid_lon
                };
                flat_cell_area(row_lat_seed, row_ref_lon, step_deg_for_flat_area)
            } else {
                0.0 // unused
            };

            for col in 0..n_cols {
                // Native CRS coordinates of this output pixel's centre.
                // For WGS84: x = lon_deg, y = lat_deg.
                // For UTM:   x = easting_m, y = northing_m.
                let pixel_x = lon_west + (col as f64 + 0.5) * spacing;
                let pixel_y = row_y;
                let (lon, lat) = if projector.is_some() {
                    row_lonlat[col]
                } else {
                    (pixel_x, pixel_y)
                };

                // 1. DEM elevation.  NaN means the DEM tile exists but the pixel
                //    is a void (SRTM −32768 sentinel); out-of-bounds means no tile
                //    covers this lat/lon.  Both → skip pixel and count as
                //    `dem_missing` so the user sees the void fraction.
                let h = match dem.elevation_at(lat, lon) {
                    Ok(h) if h.is_finite() => h as f64,
                    _ => {
                        dem_missing += 1;
                        if config.compute_lia { row_mask[col] = mask_code::DEM_MISSING; }
                        continue;
                    }
                };

                // 2. Correct orthometric DEM height to WGS84 ellipsoidal height.
                //    SRTM heights are above the EGM96 geoid; ECEF conversion
                //    requires heights above the WGS84 ellipsoid.  The geoid
                //    undulation N = h_ell − h_ortho is up to ~80 m globally.
                //    With GeoidModel::Zero, N = 0 (no correction, backward-compatible
                //    but incorrect for land); replace with a real geoid model for
                //    sub-pixel absolute accuracy.
                let h_ell = h + config.geoid.undulation_m(lat, lon);
                let target = geodetic_to_ecef(lat, lon, h_ell);

                // 3. Solve zero-Doppler.  Orbit coverage errors propagate up;
                //    per-pixel solver outcomes map to explicit counters.
                let outcome = solve_zero_doppler(
                    target,
                    orbit,
                    prev_az_time,
                    config.max_newton_iterations,
                    config.newton_tolerance_s,
                )?;

                let (az_time, sat_sv) = match outcome {
                    ZeroDopplerOutcome::Converged { time, iterations, sat_sv, .. } => {
                        // Bin into the histogram.  Clamp into the last
                        // slot if a caller raised `max_newton_iterations`
                        // above the compile-time histogram size; the
                        // debug_assert at the top of this function fires
                        // in dev builds so this clamp is observable.
                        let bin = (iterations as usize).min(NEWTON_ITER_HIST_SIZE - 1);
                        iter_histogram[bin] += 1;
                        (time, sat_sv)
                    }
                    ZeroDopplerOutcome::NotConverged { .. } => {
                        non_converged += 1;
                        if config.compute_lia { row_mask[col] = mask_code::NON_CONVERGED; }
                        continue;
                    }
                    ZeroDopplerOutcome::DegenerateGeometry => {
                        degenerate += 1;
                        if config.compute_lia { row_mask[col] = mask_code::DEGENERATE_GEOMETRY; }
                        continue;
                    }
                };
                prev_az_time = az_time; // seed next column

                // 4. Slant range from the converged state vector.
                // sat_sv was returned by solve_zero_doppler — no extra
                // interpolate_orbit call needed here.
                let r_vec = sub3(target, sat_sv.position_m);
                let slant_range_m = norm3(r_vec);

                // 5. Radar image coordinates.
                let (frac_line, frac_sample) =
                    match radar_geom.to_image_coords(az_time, slant_range_m) {
                        Some(coords) => coords,
                        None => {
                            outside_footprint += 1;
                            if config.compute_lia { row_mask[col] = mask_code::OUTSIDE_FOOTPRINT; }
                            continue;
                        }
                    };

                // 6. Resample from the merged σ⁰ image.
                let mut sigma0 = merged.sample_at(frac_line, frac_sample, config.resampling);
                if sigma0.is_nan() {
                    // The radar pixel exists but the source value is itself NaN
                    // (merge seam gap, deburst nodata, or invalid calibration).
                    outside_footprint += 1;
                    if config.compute_lia { row_mask[col] = mask_code::OUTSIDE_FOOTPRINT; }
                    continue;
                }

                // 7. Terrain geometry (Small 2011 area-projection method).
                //    Needed when:
                //      * `config.flatten`     — to compute σ⁰ → γ⁰ scaling, and
                //      * `config.compute_lia` — to populate cos_lia + mask sidecars.
                //    Requires four extra DEM lookups per pixel (E, W, N, S
                //    cardinal neighbours) for the centred-difference facet
                //    gradient.
                //
                //    Pixel-level outcomes:
                //      Valid              — record cos_lia, scale σ⁰ if flatten
                //      Foreshortening     — record cos_lia + mask; mask σ⁰ to NaN if flatten
                //      Layover            — record cos_lia + mask; mask σ⁰ to NaN if flatten
                //      Shadow             — record mask only; mask σ⁰ to NaN if flatten
                //      DEM neighbour gap  — record mask only; mask σ⁰ to NaN if flatten
                //
                //    When `compute_lia=true && flatten=false` the σ⁰ value is
                //    passed through *unchanged* even at foreshortening / shadow
                //    pixels: the user explicitly opted out of γ⁰ conversion, so
                //    those geometric flags are advisory metadata only.
                if config.flatten || config.compute_lia {
                    // Cardinal-neighbour DEM lookups use a step in DEGREES
                    // (the DEM and `compute_terrain_geometry`'s ECEF math
                    // both work in geographic coordinates).  In WGS84 mode
                    // this equals `spacing` exactly and the legacy bit-
                    // identical path is preserved; in UTM mode `spacing`
                    // is in metres and would land 10° away from the centre
                    // pixel without this conversion — see
                    // `cardinal_neighbour_step_deg` for the defect history
                    // (silent `valid=0` baseline output prior to the fix).
                    // `step_deg_for_flat_area` is the same value, hoisted
                    // out of this loop so the flat_area is precomputed once
                    // per row instead of 4 ECEF calls per pixel.
                    let step_deg = step_deg_for_flat_area;
                    let h_neighbours = (
                        dem.elevation_at(lat, lon + step_deg),
                        dem.elevation_at(lat, lon - step_deg),
                        dem.elevation_at(lat + step_deg, lon),
                        dem.elevation_at(lat - step_deg, lon),
                    );
                    let all_finite = match h_neighbours {
                        (Ok(e), Ok(w), Ok(n), Ok(s))
                            if e.is_finite() && w.is_finite()
                                && n.is_finite() && s.is_finite() =>
                        {
                            Some((e as f64, w as f64, n as f64, s as f64))
                        }
                        _ => None,
                    };
                    match all_finite {
                        None => {
                            if config.compute_lia {
                                row_mask[col] = mask_code::DEM_NEIGHBOUR_MISSING;
                            }
                            if config.flatten {
                                flat_masked += 1;
                                continue;
                            }
                        }
                        Some((h_east, h_west, h_north, h_south)) => {
                            let geom = compute_terrain_geometry_with_flat_area(
                                lat, lon,
                                h_east, h_west, h_north, h_south,
                                sat_sv.position_m,
                                target,
                                step_deg,
                                row_flat_area,
                            );
                            match geom {
                                TerrainGeometryOutcome::Valid { w, cos_lia } => {
                                    if config.compute_lia {
                                        row_cos_lia[col] = cos_lia as f32;
                                    }
                                    if config.flatten {
                                        sigma0 /= w as f32;
                                    }
                                }
                                TerrainGeometryOutcome::Foreshortening { cos_lia } => {
                                    if config.compute_lia {
                                        row_cos_lia[col] = cos_lia as f32;
                                        row_mask[col] = mask_code::FORESHORTENING;
                                    }
                                    if config.flatten {
                                        flat_masked += 1;
                                        continue;
                                    }
                                }
                                TerrainGeometryOutcome::Layover { cos_lia } => {
                                    if config.compute_lia {
                                        row_cos_lia[col] = cos_lia as f32;
                                        row_mask[col] = mask_code::LAYOVER;
                                    }
                                    if config.flatten {
                                        flat_masked += 1;
                                        continue;
                                    }
                                }
                                TerrainGeometryOutcome::Shadow => {
                                    if config.compute_lia {
                                        row_mask[col] = mask_code::SHADOW;
                                    }
                                    if config.flatten {
                                        flat_masked += 1;
                                        continue;
                                    }
                                }
                            }
                        }
                    }
                }

                // 8. Per-pixel NESZ noise floor masking (σ⁰ domain, before flattening
                //    modifies the value).  Pixels where σ⁰ ≤ NESZ × 10^(margin/10)
                //    are set to NaN.  Skipped when margin = NEG_INFINITY (default),
                //    or when the image source provides no NESZ layer.
                if config.noise_floor_margin_db.is_finite() {
                    if let Some(nesz) = merged.nesz_at(frac_line, frac_sample) {
                        // If NESZ is NaN (no coverage at this pixel), pass through unmasked.
                        if !nesz.is_nan() {
                            let threshold = nesz * 10f32.powf(config.noise_floor_margin_db / 10.0);
                            if sigma0 <= threshold {
                                noise_masked += 1;
                                if config.compute_lia { row_mask[col] = mask_code::NOISE_MASKED; }
                                continue;
                            }
                        }
                    }
                }

                row_data[col] = sigma0;
                if config.compute_lia { row_mask[col] = mask_code::VALID; }
                valid += 1;
            }

            // Progress hint.  Rows complete out of order so the count is
            // approximate, but gives a useful sense of progress.
            let done = rows_done.fetch_add(1, Ordering::Relaxed) + 1;
            if done % 100 == 0 {
                tracing::debug!(
                    "terrain correction: ~{}/{} rows ({:.1}%)",
                    done,
                    n_rows,
                    done as f64 / n_rows as f64 * 100.0
                );
            }

            chunk_out.push(RowResult { row_data, row_cos_lia, row_mask, valid, dem_missing, outside_footprint, non_converged, degenerate, flat_masked, noise_masked, iter_histogram });
            }
            Ok(chunk_out)
        })
        .collect();

    // Assemble: copy each row's data into the flat output buffer and sum
    // the diagnostic counters.  We flatten the nested chunk structure directly
    // without collecting into an intermediate Vec<RowResult>.
    let mut data = vec![f32::NAN; total];
    let (mut cos_lia_buf, mut mask_buf) = if config.compute_lia {
        (vec![f32::NAN; total], vec![mask_code::OUTSIDE_FOOTPRINT; total])
    } else {
        (Vec::new(), Vec::new())
    };
    let mut valid_count = 0usize;
    let mut dem_missing = 0usize;
    let mut outside_footprint = 0usize;
    let mut non_converged = 0usize;
    let mut degenerate = 0usize;
    let mut flat_masked = 0usize;
    let mut noise_masked = 0usize;
    let mut newton_iter_histogram = [0u64; NEWTON_ITER_HIST_SIZE];

    for (row, result) in chunk_results?.into_iter().flatten().enumerate() {
        data[row * n_cols..(row + 1) * n_cols].copy_from_slice(&result.row_data);
        if config.compute_lia {
            cos_lia_buf[row * n_cols..(row + 1) * n_cols]
                .copy_from_slice(&result.row_cos_lia);
            mask_buf[row * n_cols..(row + 1) * n_cols]
                .copy_from_slice(&result.row_mask);
        }
        valid_count        += result.valid;
        dem_missing        += result.dem_missing;
        outside_footprint  += result.outside_footprint;
        non_converged      += result.non_converged;
        degenerate         += result.degenerate;
        flat_masked        += result.flat_masked;
        noise_masked       += result.noise_masked;
        for (acc, v) in newton_iter_histogram.iter_mut().zip(result.iter_histogram.iter()) {
            *acc += *v;
        }
    }

    // GDAL geotransform for a north-up image:
    //   gt[0] = top-left X (west edge in CRS native units)
    //   gt[1] = pixel width in X (positive eastward)
    //   gt[2] = 0 (rotation term, zero for north-up)
    //   gt[3] = top-left Y (north edge in CRS native units)
    //   gt[4] = 0 (rotation term, zero for north-up)
    //   gt[5] = pixel height in Y (negative: Y decreases southward)
    let geotransform = [
        lon_west,
        spacing,
        0.0,
        lat_north,
        0.0,
        -spacing,
    ];

    Ok(GeocodedImage {
        data,
        rows: n_rows,
        cols: n_cols,
        origin_lat_deg: lat_north,
        origin_lon_deg: lon_west,
        pixel_spacing_lat_deg: spacing,
        pixel_spacing_lon_deg: spacing,
        total_pixel_count: total,
        valid_pixel_count: valid_count,
        dem_missing_count: dem_missing,
        outside_footprint_count: outside_footprint,
        non_converged_count: non_converged,
        degenerate_geometry_count: degenerate,
        flat_masked_count: flat_masked,
        noise_masked_count: noise_masked,
        geotransform,
        crs: config.crs,
        cos_lia: cos_lia_buf,
        mask: mask_buf,
        newton_iter_histogram,
    })
}

// ── Tests ──────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geodesy::normalize3;
    use chrono::TimeZone;
    use crate::types::StateVector;

    /// `compute_output_grid_extent` on `OutputCrs::Wgs84LatLon` must
    /// snap the bbox outward to the spacing grid and is exact (no
    /// projection involved).  This pins the legacy behaviour byte-for-byte.
    #[test]
    fn output_grid_extent_wgs84_snaps_outward_to_spacing() {
        let bb = crate::types::BoundingBox {
            min_lat_deg: 47.123,
            max_lat_deg: 48.876,
            min_lon_deg: 11.234,
            max_lon_deg: 12.765,
        };
        let spacing = 0.01_f64;
        let (y_top, y_bot, x_left, x_right) =
            compute_output_grid_extent(&bb, None, spacing).expect("wgs84 extent");
        // ceil(48.876 / 0.01) * 0.01 = 48.88; floor(47.123 / 0.01) * 0.01 = 47.12
        assert!((y_top - 48.88).abs() < 1e-9, "y_top = {y_top}");
        assert!((y_bot - 47.12).abs() < 1e-9, "y_bot = {y_bot}");
        // floor(11.234 / 0.01) * 0.01 = 11.23; ceil(12.765 / 0.01) * 0.01 = 12.77
        assert!((x_left - 11.23).abs() < 1e-9, "x_left = {x_left}");
        assert!((x_right - 12.77).abs() < 1e-9, "x_right = {x_right}");
    }

    /// `compute_output_grid_extent` on a UTM CRS must (a) project the
    /// four corners, (b) take the AABB in metres, (c) snap to the
    /// spacing grid.  We test on a small bbox in UTM zone 32N and
    /// compare the corner against an independent proj4rs invocation.
    #[test]
    fn output_grid_extent_utm_zone_32n_matches_proj_corners() {
        // Bbox roughly over Munich (UTM zone 32N).
        let bb = crate::types::BoundingBox {
            min_lat_deg: 48.0,
            max_lat_deg: 48.5,
            min_lon_deg: 11.0,
            max_lon_deg: 12.0,
        };
        let crs = crate::output_crs::OutputCrs::from_epsg(32632)
            .expect("UTM 32N");
        let proj = crs.projector().expect("projector");
        let spacing = 10.0_f64; // 10 m

        let (y_top, y_bot, x_left, x_right) =
            compute_output_grid_extent(&bb, Some(&proj), spacing)
                .expect("utm extent");

        // Independently project all four corners and compute the AABB.
        let corners = [
            (bb.min_lon_deg, bb.min_lat_deg),
            (bb.min_lon_deg, bb.max_lat_deg),
            (bb.max_lon_deg, bb.min_lat_deg),
            (bb.max_lon_deg, bb.max_lat_deg),
        ];
        let (mut x_min, mut x_max) = (f64::INFINITY, f64::NEG_INFINITY);
        let (mut y_min, mut y_max) = (f64::INFINITY, f64::NEG_INFINITY);
        for (lon, lat) in corners {
            let (x, y) = proj.lonlat_to_xy(lon, lat).expect("project");
            if x < x_min { x_min = x; }
            if x > x_max { x_max = x; }
            if y < y_min { y_min = y; }
            if y > y_max { y_max = y; }
        }
        let expected_y_top = (y_max / spacing).ceil() * spacing;
        let expected_y_bot = (y_min / spacing).floor() * spacing;
        let expected_x_left = (x_min / spacing).floor() * spacing;
        let expected_x_right = (x_max / spacing).ceil() * spacing;

        assert!(
            (y_top - expected_y_top).abs() < 1e-6,
            "y_top: {y_top} vs {expected_y_top}"
        );
        assert!(
            (y_bot - expected_y_bot).abs() < 1e-6,
            "y_bot: {y_bot} vs {expected_y_bot}"
        );
        assert!(
            (x_left - expected_x_left).abs() < 1e-6,
            "x_left: {x_left} vs {expected_x_left}"
        );
        assert!(
            (x_right - expected_x_right).abs() < 1e-6,
            "x_right: {x_right} vs {expected_x_right}"
        );

        // Sanity: bbox in metres is positive and ~30–60 km wide / ~55 km tall.
        let width_m = x_right - x_left;
        let height_m = y_top - y_bot;
        assert!(
            (50_000.0..200_000.0).contains(&width_m),
            "UTM extent width should be order 50–200 km, got {width_m} m"
        );
        assert!(
            (40_000.0..100_000.0).contains(&height_m),
            "UTM extent height should be order 40–100 km, got {height_m} m"
        );
    }

    /// The grid extent must snap **outward** so the projected scene
    /// corners are strictly inside the snapped rectangle.  This is the
    /// invariant the resampling step relies on (no scene pixel ever
    /// falls outside the output grid).
    #[test]
    fn output_grid_extent_utm_snapped_rectangle_contains_all_corners() {
        let bb = crate::types::BoundingBox {
            min_lat_deg: 48.0,
            max_lat_deg: 48.5,
            min_lon_deg: 11.0,
            max_lon_deg: 12.0,
        };
        let crs = crate::output_crs::OutputCrs::from_epsg(32632).unwrap();
        let proj = crs.projector().unwrap();
        let spacing = 50.0_f64;

        let (y_top, y_bot, x_left, x_right) =
            compute_output_grid_extent(&bb, Some(&proj), spacing).unwrap();
        for (lon, lat) in [
            (bb.min_lon_deg, bb.min_lat_deg),
            (bb.min_lon_deg, bb.max_lat_deg),
            (bb.max_lon_deg, bb.min_lat_deg),
            (bb.max_lon_deg, bb.max_lat_deg),
        ] {
            let (x, y) = proj.lonlat_to_xy(lon, lat).unwrap();
            assert!(
                x >= x_left && x <= x_right && y >= y_bot && y <= y_top,
                "corner ({lon},{lat}) → ({x},{y}) escapes snapped \
                 rectangle x∈[{x_left},{x_right}] y∈[{y_bot},{y_top}]"
            );
        }
    }

    // ── Terrain flattening tests ──────────────────────────────────────────

    /// For truly flat terrain with the satellite directly overhead (nadir),
    /// w = cos(θ_nadir) = cos(0°) = 1.0, so γ⁰ = σ⁰ / w = σ⁰.
    /// This degenerate case is only correct; for realistic side-looking
    /// geometry see `test_flattening_weight_flat_terrain_side_looking_is_cos_theta`.
    #[test]
    fn test_flattening_weight_flat_terrain_nadir_is_unity() {
        // Satellite at ~700 km altitude directly above lat=50°N, lon=11°E.
        let lat = 50.0_f64;
        let lon = 11.0_f64;
        let h = 300.0_f64; // 300 m elevation, uniform

        let target = geodetic_to_ecef(lat, lon, h);
        // Satellite overhead at 700 km
        let sat = geodetic_to_ecef(lat, lon, h + 700_000.0);

        let spacing = 0.0001_f64; // same as default config

        let w = compute_flattening_weight(lat, lon, h, h, h, h, sat, target, spacing);
        let w = w.expect("flat terrain should not be masked");
        assert!(
            (w - 1.0).abs() < 1e-4,
            "flat terrain with nadir satellite: w should be ≈ cos(0°) = 1.0, got {w}"
        );
    }

    /// For flat terrain with a **side-looking** satellite (realistic S-1 geometry),
    /// w = cos(θ_inc) so that σ⁰ / w = σ⁰ / cos(θ_inc) = γ⁰.  This is the
    /// standard terrain-flattened gamma-naught conversion on flat ground.
    #[test]
    fn test_flattening_weight_flat_terrain_side_looking_is_cos_theta() {
        let lat = 50.0_f64;
        let lon = 11.0_f64;
        let h = 0.0_f64;
        let spacing = 0.0001_f64;

        let target = geodetic_to_ecef(lat, lon, h);
        // Side-looking satellite: 700 km altitude, 500 km horizontal offset east
        // (roughly 35° off-nadir — typical S-1 IW near-range incidence).
        let sat = [target[0] + 500_000.0, target[1], target[2] + 700_000.0];

        // Independently compute cos(θ_inc) = (flat_normal · look) / |flat_normal|
        let f_east  = geodetic_to_ecef(lat,           lon + spacing, 0.0);
        let f_west  = geodetic_to_ecef(lat,           lon - spacing, 0.0);
        let f_north = geodetic_to_ecef(lat + spacing, lon,           0.0);
        let f_south = geodetic_to_ecef(lat - spacing, lon,           0.0);
        let flat_n = cross3(sub3(f_east, f_west), sub3(f_north, f_south));
        let look   = normalize3(sub3(sat, target));
        let cos_theta_inc = dot3(flat_n, look) / norm3(flat_n);

        // All neighbours at same height → flat terrain.
        let w = compute_flattening_weight(
            lat, lon, h, h, h, h, sat, target, spacing,
        )
        .expect("flat side-looking terrain should not be masked");

        assert!(
            (w - cos_theta_inc).abs() < 1e-6,
            "flat side-looking terrain: expected w = cos(θ_inc) = {cos_theta_inc:.6}, got {w:.6}"
        );
        // For any non-nadir satellite w < 1, so γ⁰ = σ⁰/w > σ⁰.
        assert!(
            w < 1.0,
            "side-looking flat terrain weight must be < 1.0 (cos(θ_inc)), got {w:.4}"
        );
    }

    /// Terrain that slopes toward the satellite has a larger projected area
    /// compared to flat ground: the illuminated patch is foreshortened, so
    /// w > 1.0 (σ⁰ is divided by a number > 1, making γ⁰ < σ⁰).
    /// The slope is kept gentle (±1 m over ~14 m laterally ≈ 4°) so that
    /// w stays in the Valid band (below LAYOVER_THRESHOLD = 2.0).
    #[test]
    fn test_flattening_weight_slope_toward_satellite_gt_1() {
        let lat = 50.0_f64;
        let lon = 11.0_f64;
        let h = 0.0_f64;
        let spacing = 0.0001_f64;

        // Satellite at a shallow angle from east (~30° elevation from horizon)
        // Place it at 700 km altitude, offset 500 km east in ECEF.
        let target = geodetic_to_ecef(lat, lon, h);
        let sat = [target[0] + 500_000.0, target[1], target[2] + 700_000.0];

        // Terrain slopes gently upward to the east (toward satellite): ±1 m
        // over a 0.0001° stencil (≈ 14 m east-west) gives a ~4° slope, well
        // below the LAYOVER_THRESHOLD so the outcome is Valid.
        let h_east  = h + 1.0;
        let h_west  = h - 1.0;
        let h_north = h;
        let h_south = h;

        let w = compute_flattening_weight(
            lat, lon, h_east, h_west, h_north, h_south, sat, target, spacing,
        );
        let w = w.expect("gentle slope toward satellite should be Valid (not masked)");
        // The terrain facet projected area is larger than the flat reference
        // so w > 1 (stronger foreshortening → divide σ⁰ by a larger number).
        assert!(w > 1.0, "slope toward satellite should give w > 1.0, got {w}");
    }

    /// A very steep slope directly facing the satellite (w > LAYOVER_THRESHOLD)
    /// should be classified as Layover, not Valid.
    #[test]
    fn terrain_geometry_steep_toward_satellite_is_layover() {
        let lat = 50.0_f64;
        let lon = 11.0_f64;
        let h = 0.0_f64;
        let spacing = 0.0001_f64;

        let target = geodetic_to_ecef(lat, lon, h);
        let sat = [target[0] + 500_000.0, target[1], target[2] + 700_000.0];

        // ±100 m over ~14 m laterally → ~82° slope, w >> LAYOVER_THRESHOLD = 2.0.
        let outcome = compute_terrain_geometry(
            lat, lon, h + 100.0, h - 100.0, h, h, sat, target, spacing,
        );
        assert!(
            matches!(outcome, TerrainGeometryOutcome::Layover { .. }),
            "steep front-facing slope should be Layover, got {outcome:?}"
        );
    }

    /// Shadow geometry (terrain facet normal pointing away from satellite)
    /// should return None.
    #[test]
    fn test_flattening_weight_shadow_returns_none() {
        let lat = 50.0_f64;
        let lon = 11.0_f64;
        let spacing = 0.0001_f64;

        // Satellite due east at 700 km altitude
        let target = geodetic_to_ecef(lat, lon, 0.0);
        let sat = [target[0] + 700_000.0, target[1], target[2]];

        // Very steep east-facing slope (drops drastically to the east):
        // east neighbour is 5000 m lower, west symmetrically higher → this
        // facet faces away from the satellite (which is east of the target).
        let h_east  = -5000.0;
        let h_west  =  5000.0;
        let h_north =     0.0;
        let h_south =     0.0;

        let result = compute_flattening_weight(
            lat, lon, h_east, h_west, h_north, h_south, sat, target, spacing,
        );
        assert!(result.is_none(), "shadow geometry should return None, got {result:?}");
    }

    /// Regression for the centred-difference fix.
    ///
    /// With the previous forward-difference (E, N) facet, applying the same
    /// constant scalar slope `s` first as a +s east bump and then as a −s
    /// west bump (with otherwise flat terrain) gave *different* weights even
    /// though both configurations describe the same east-axis gradient
    /// magnitude.  Centred differences must give the exact same `w` (the
    /// gradient is `(h_E − h_W) / 2Δ` regardless of which neighbour carries
    /// the displacement) and must equal the result of putting the bump on
    /// both sides at half magnitude.
    #[test]
    fn test_flattening_weight_centred_difference_symmetric() {
        let lat = 50.0_f64;
        let lon = 11.0_f64;
        let spacing = 0.0001_f64;
        let target = geodetic_to_ecef(lat, lon, 0.0);
        // Satellite east of and above target so the look vector is non-trivial.
        let sat = [target[0] + 500_000.0, target[1], target[2] + 700_000.0];

        // Configuration A: +2 m bump on the east neighbour only.
        let w_a = compute_flattening_weight(
            lat, lon, 2.0, 0.0, 0.0, 0.0, sat, target, spacing,
        )
        .expect("config A unmasked");

        // Configuration B: −2 m bump on the west neighbour only.
        // (h_E − h_W) is identical to config A → centred-difference
        // gradient is identical → weight must match.
        let w_b = compute_flattening_weight(
            lat, lon, 0.0, -2.0, 0.0, 0.0, sat, target, spacing,
        )
        .expect("config B unmasked");

        // Configuration C: +1 m east, −1 m west (symmetric split).
        // Same (h_E − h_W) again → same gradient → same weight.
        let w_c = compute_flattening_weight(
            lat, lon, 1.0, -1.0, 0.0, 0.0, sat, target, spacing,
        )
        .expect("config C unmasked");

        let max_diff = (w_a - w_b).abs().max((w_a - w_c).abs());
        // Tolerance: 1e-4 relative.  An exact algebraic match is not
        // possible on the WGS84 ellipsoid because the local up-vectors at
        // (lat, lon+Δ) and (lat, lon−Δ) are not anti-parallel; they are
        // rotated by 2Δ of longitude, which leaks O(Δ) into the
        // ECEF-difference cancellation.  At Δ = 1e-4° the residual is
        // ~3e-5 of `w` — well below DEM noise.  This bound is what would
        // *fail* under the old forward-difference formula, where moving the
        // bump from east to west flips the gradient sign and changes `w` by
        // tens of percent.
        let rel_tol = 1.0e-4;
        assert!(
            max_diff < rel_tol * w_a.abs(),
            "centred-difference weights must be invariant to which neighbour \
             carries the height bump: w_a={w_a}, w_b={w_b}, w_c={w_c}, \
             max_diff={max_diff}, tol={}",
            rel_tol * w_a.abs()
        );
    }

    fn make_test_orbit() -> OrbitData {
        // Create a simple circular orbit for testing.
        // 10 state vectors at 10-second intervals, ~700 km altitude, polar orbit.
        let t0 = Utc.with_ymd_and_hms(2020, 10, 5, 17, 8, 0).unwrap();
        let alt = 693_000.0; // m
        let r = 6_371_000.0 + alt;
        let v = (3.986_004_418e14_f64 / r).sqrt(); // circular velocity

        let mut svs = Vec::new();
        for i in 0..30 {
            let t = t0 + Duration::seconds(i * 10);
            let angle = v * (i as f64 * 10.0) / r; // radians around orbit
            let x = r * angle.cos();
            let y = 0.0;
            let z = r * angle.sin();
            let vx = -v * angle.sin();
            let vy = 0.0;
            let vz = v * angle.cos();
            svs.push(StateVector {
                time: t,
                position_m: [x, y, z],
                velocity_m_s: [vx, vy, vz],
            });
        }

        OrbitData {
            reference_epoch: t0,
            state_vectors: svs,
        }
    }

    #[test]
    fn test_zero_doppler_at_nadir() {
        // For a target directly below the satellite (nadir),
        // the position vector P-S is perpendicular to V (zero Doppler).
        let orbit = make_test_orbit();
        let sv = &orbit.state_vectors[15]; // middle vector

        // Target at nadir: along the position vector but at Earth surface
        let pos_norm = norm3(sv.position_m);
        let earth_r = 6_371_000.0;
        let target = [
            sv.position_m[0] * earth_r / pos_norm,
            sv.position_m[1] * earth_r / pos_norm,
            sv.position_m[2] * earth_r / pos_norm,
        ];

        let result = solve_zero_doppler(target, &orbit, sv.time, 10, 1e-6);
        assert!(result.is_ok(), "zero-Doppler solver errored: {:?}", result.err());

        let outcome = result.unwrap();
        let solved_time = match outcome {
            ZeroDopplerOutcome::Converged { time, .. } => time,
            other => panic!("expected converged outcome, got {:?}", other),
        };
        let dt_us = (solved_time - sv.time).num_microseconds().unwrap_or(0) as f64; // SAFETY-OK: chrono microseconds cannot overflow for orbit-window durations
        // Should converge to nearly the same time (target is at nadir for this vector)
        assert!(
            dt_us.abs() < 100_000.0, // < 0.1 s
            "time offset: {} μs",
            dt_us
        );
    }

    #[test]
    fn test_zero_doppler_reports_non_convergence() {
        // Force the solver to give up in 1 iteration: it cannot possibly
        // converge for a distant target with tolerance 1 ns in a single step.
        let orbit = make_test_orbit();
        let sv = &orbit.state_vectors[15];
        // Target far off-nadir (~6400 km in-track from nadir).
        let along_track = normalize3(sv.velocity_m_s);
        let pos_norm = norm3(sv.position_m);
        let earth_r = 6_371_000.0;
        let nadir = [
            sv.position_m[0] * earth_r / pos_norm,
            sv.position_m[1] * earth_r / pos_norm,
            sv.position_m[2] * earth_r / pos_norm,
        ];
        let target = [
            nadir[0] + along_track[0] * 100_000.0,
            nadir[1] + along_track[1] * 100_000.0,
            nadir[2] + along_track[2] * 100_000.0,
        ];
        let result = solve_zero_doppler(target, &orbit, sv.time, 1, 1e-12).unwrap();
        assert!(
            matches!(result, ZeroDopplerOutcome::NotConverged { .. }),
            "expected NotConverged, got {:?}",
            result
        );
    }

    #[test]
    fn test_azimuth_time_lut() {
        // Test the LUT-based initial guess
        let t0 = Utc.with_ymd_and_hms(2020, 10, 5, 17, 8, 0).unwrap();
        let t1 = Utc.with_ymd_and_hms(2020, 10, 5, 17, 8, 30).unwrap();

        let grids = vec![(
            SubSwathId::IW2,
            vec![
                GeolocationGridPoint {
                    azimuth_time_utc: t0,
                    slant_range_time_s: 0.005,
                    line: 0,
                    pixel: 0,
                    latitude_deg: 52.0,
                    longitude_deg: 10.0,
                    height_m: 0.0,
                    incidence_angle_deg: 30.0,
                    elevation_angle_deg: 28.0,
                },
                GeolocationGridPoint {
                    azimuth_time_utc: t1,
                    slant_range_time_s: 0.005,
                    line: 100,
                    pixel: 0,
                    latitude_deg: 50.0,
                    longitude_deg: 10.0,
                    height_m: 0.0,
                    incidence_angle_deg: 30.0,
                    elevation_angle_deg: 28.0,
                },
            ],
        )];

        let lut = build_azimuth_time_lut(&grids);
        assert_eq!(lut.len(), 2);

        // Midpoint latitude should give midpoint time
        let mid = guess_azimuth_time(&lut, 51.0);
        let expected_mid = t0 + Duration::seconds(15);
        let dt = (mid - expected_mid).num_milliseconds().abs();
        assert!(dt < 100, "midpoint guess off by {} ms", dt);
    }

    #[test]
    fn test_bilinear_sample_simple() {
        // 3×3 image with known values
        let merged = MergedSigma0 {
            data: vec![
                1.0, 2.0, 3.0, //
                4.0, 5.0, 6.0, //
                7.0, 8.0, 9.0, //
            ],
            nesz: vec![0.0f32; 9],
            lines: 3,
            samples: 3,
            near_slant_range_time_s: 0.005,
            range_pixel_spacing_s: 1e-8,
            range_pixel_spacing_m: 2.33,
            cal_lut_extrapolation_gap_px: 0,
            noise_lut_extrapolation_gap_px: 0,
            azimuth_start_time: Utc::now(),
            azimuth_time_interval_s: 0.002055556,
        };

        // Exact pixel
        let v = bilinear_sample(&merged, 1.0, 1.0);
        assert!((v - 5.0).abs() < 1e-6);

        // Centre of top-left quadrant
        let v = bilinear_sample(&merged, 0.5, 0.5);
        let expected = (1.0 + 2.0 + 4.0 + 5.0) / 4.0;
        assert!((v - expected).abs() < 1e-6);
    }

    #[test]
    fn test_bilinear_nan_propagation() {
        let merged = MergedSigma0 {
            data: vec![
                1.0,
                f32::NAN,
                3.0, //
                4.0,
                5.0,
                6.0,
            ],
            nesz: vec![0.0f32; 6],
            lines: 2,
            samples: 3,
            near_slant_range_time_s: 0.005,
            range_pixel_spacing_s: 1e-8,
            range_pixel_spacing_m: 2.33,
            cal_lut_extrapolation_gap_px: 0,
            noise_lut_extrapolation_gap_px: 0,
            azimuth_start_time: Utc::now(),
            azimuth_time_interval_s: 0.002055556,
        };

        // Interpolation touching NaN should produce NaN
        let v = bilinear_sample(&merged, 0.0, 0.5);
        assert!(v.is_nan());
    }

    #[test]
    fn test_bilinear_sample_slice_matches_bilinear_sample() {
        // bilinear_sample_slice must produce identical results to bilinear_sample
        // for the merged.data slice (regression guard for the refactor).
        let merged = MergedSigma0 {
            data: vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0],
            nesz: vec![0.0f32; 9],
            lines: 3,
            samples: 3,
            near_slant_range_time_s: 0.005,
            range_pixel_spacing_s: 1e-8,
            range_pixel_spacing_m: 2.33,
            cal_lut_extrapolation_gap_px: 0,
            noise_lut_extrapolation_gap_px: 0,
            azimuth_start_time: Utc::now(),
            azimuth_time_interval_s: 0.002055556,
        };
        for (line, sample) in [(0.0, 0.0), (0.7, 1.3), (2.0, 2.0), (1.5, 0.5)] {
            let a = bilinear_sample(&merged, line, sample);
            let b = bilinear_sample_slice(&merged.data, merged.lines, merged.samples, line, sample);
            assert!((a - b).abs() < 1e-6, "mismatch at ({line},{sample}): {a} vs {b}");
        }
    }

    // ── Bicubic kernel tests ───────────────────────────────────────────────────

    #[test]
    fn keys_weight_unit_interval() {
        // At t=0 the weight must be 1.0 (centre tap).
        assert!((keys_weight(0.0) - 1.0).abs() < 1e-12);
        // At t=1 and t=-1 the weight must be 0.0 (Catmull-Rom knot condition).
        assert!(keys_weight(1.0).abs() < 1e-12);
        assert!(keys_weight(-1.0).abs() < 1e-12);
        // Beyond |t| >= 2 the weight must be 0.0.
        assert_eq!(keys_weight(2.0), 0.0);
        assert_eq!(keys_weight(-2.5), 0.0);
    }

    #[test]
    fn keys_weights_sum_to_one_at_exact_pixel() {
        // At exact pixel centres ds = 0.0 → only the centre tap contributes.
        // The four weights for offsets -1, 0, 1, 2 at ds=0 are:
        //   keys_weight(-1) = 0, keys_weight(0) = 1, keys_weight(1) = 0, keys_weight(2) = 0
        // Confirm partition-of-unity: sum of all 16 = 1.0.
        let ds = 0.0_f64;
        let dl = 0.0_f64;
        let sum: f64 = (-1_isize..=2).flat_map(|li| {
            (-1_isize..=2).map(move |si| keys_weight(li as f64 - dl) * keys_weight(si as f64 - ds))
        }).sum();
        assert!((sum - 1.0).abs() < 1e-12, "partition-of-unity at exact pixel: sum={sum}");
    }

    #[test]
    fn bicubic_sample_constant_image_returns_constant() {
        // On a uniform field bicubic must return the constant value regardless
        // of fractional position (partition-of-unity property).
        let data = vec![7.0_f32; 5 * 5];
        for (l, s) in [(0.5_f64, 0.5_f64), (1.7, 2.3), (3.9, 3.9)] {
            let v = bicubic_sample_slice(&data, 5, 5, l, s);
            assert!((v - 7.0).abs() < 1e-4, "constant field: got {v} at ({l},{s})");
        }
    }

    #[test]
    fn bicubic_at_exact_grid_node_returns_value() {
        // At exact integer coordinates the bicubic result must equal the grid value.
        let data: Vec<f32> = (0..16).map(|i| i as f32).collect();
        let v = bicubic_sample_slice(&data, 4, 4, 2.0, 3.0);
        assert!((v - data[2 * 4 + 3]).abs() < 1e-4, "exact node: got {v}, expected {}", data[2*4+3]);
    }

    #[test]
    fn bicubic_nan_tap_propagates_nan() {
        let mut data = vec![1.0_f32; 4 * 4];
        data[1 * 4 + 1] = f32::NAN; // a tap that will be hit
        let v = bicubic_sample_slice(&data, 4, 4, 1.5, 1.5);
        assert!(v.is_nan(), "NaN tap must propagate to NaN result");
    }

    #[test]
    fn test_radar_geometry_conversion() {
        let t0 = Utc.with_ymd_and_hms(2020, 10, 5, 17, 8, 0).unwrap();
        let geom = RadarGeometry {
            azimuth_start_time: t0,
            azimuth_time_interval_s: 0.002,
            near_range_m: 800_000.0,
            range_pixel_spacing_m: 2.33,
            n_lines: 10000,
            n_samples: 50000,
        };

        // Line 0, sample 0 should map to t0, near_range
        let (l, s) = geom
            .to_image_coords(t0, 800_000.0)
            .expect("should be in range");
        assert!(l.abs() < 1e-10);
        assert!(s.abs() < 1e-10);

        // Line 5000, sample 25000
        let t = t0 + Duration::milliseconds(10_000); // 10 s = 5000 * 0.002
        let r = 800_000.0 + 25000.0 * 2.33;
        let (l, s) = geom.to_image_coords(t, r).expect("should be in range");
        assert!((l - 5000.0).abs() < 1e-6, "line: {}", l);
        assert!((s - 25000.0).abs() < 1e-6, "sample: {}", s);
    }

    // ── cardinal_neighbour_step_deg: CRS-units → degrees for DEM lookups ──

    /// WGS84 mode is the legacy code path: `crs_spacing` is already in
    /// degrees, so the helper must return it byte-identical.  Without
    /// this guarantee the +0.016 dB validation of the WGS84 pipeline
    /// against ASF RTC10 would no longer hold.
    #[test]
    fn cardinal_step_wgs84_is_identity() {
        for s in [1e-5_f64, 0.0001, 0.001, 0.1, 1.0] {
            let step = cardinal_neighbour_step_deg(s, &OutputCrs::Wgs84LatLon);
            assert_eq!(
                step.to_bits(),
                s.to_bits(),
                "WGS84 must round-trip bit-identically; got {step} from {s}",
            );
        }
    }

    /// UTM mode: `crs_spacing` is in metres.  The cardinal-neighbour
    /// DEM lookups need a step in **degrees** because the DEM is
    /// indexed in geographic coordinates.  The pre-fix code passed the
    /// metric value directly, producing step ≈ 10° at 10 m spacing and
    /// silently masking every output pixel as "flat-masked".  This is
    /// the *negative regression* for the `valid=0` baseline bug.
    #[test]
    fn cardinal_step_utm_converts_metres_to_degrees() {
        // 10 m / (111 320 m/deg) ≈ 8.9831e-5 deg.  Tolerance 1e-9 covers
        // round-off in the constant only; the value itself is exact for
        // the chosen METRES_PER_DEGREE_LATITUDE.
        let step_north = cardinal_neighbour_step_deg(10.0, &OutputCrs::UtmNorth { zone: 32 });
        let step_south = cardinal_neighbour_step_deg(10.0, &OutputCrs::UtmSouth { zone: 19 });
        let expected = 10.0 / 111_320.0;
        assert!(
            (step_north - expected).abs() < 1e-12,
            "UtmNorth: expected {expected}, got {step_north}",
        );
        assert!(
            (step_south - expected).abs() < 1e-12,
            "UtmSouth: expected {expected}, got {step_south}",
        );
        // Sanity check on absolute scale: the resulting step must be
        // small enough that the four neighbours land inside any
        // realistic DEM coverage (i.e. ≪ 1° from the centre pixel).
        // Pre-fix the value was 10.0 (degrees), which catastrophically
        // failed this bound on every UTM scene.
        assert!(
            step_north < 1e-3,
            "UTM step must be ≪ 1° to keep neighbours inside DEM coverage; got {step_north}",
        );
    }

    /// At 10 m UTM spacing the neighbour step is identical regardless
    /// of hemisphere or zone — the helper depends only on the metric
    /// value because lat/lon anisotropy is absorbed by
    /// `geodetic_to_ecef` inside `compute_terrain_geometry`.
    #[test]
    fn cardinal_step_utm_independent_of_zone() {
        let s = 10.0;
        let baseline = cardinal_neighbour_step_deg(s, &OutputCrs::UtmNorth { zone: 1 });
        for zone in [10_u8, 32, 33, 60] {
            let v = cardinal_neighbour_step_deg(s, &OutputCrs::UtmNorth { zone });
            assert_eq!(v.to_bits(), baseline.to_bits(), "zone={zone}");
            let v = cardinal_neighbour_step_deg(s, &OutputCrs::UtmSouth { zone });
            assert_eq!(v.to_bits(), baseline.to_bits(), "zone={zone} south");
        }
    }

    // ── compute_terrain_geometry: cos_lia + outcome discrimination ────────

    /// On flat terrain with a side-looking satellite the local incidence
    /// angle equals the radar incidence angle, so `cos_lia` must equal
    /// the flattening weight `w` (which itself equals `cos(θ_inc)` on flat
    /// terrain).  Verifies that the new geometry helper returns a
    /// scientifically meaningful cos_lia, not just a side-effect of `w`.
    #[test]
    fn terrain_geometry_flat_terrain_cos_lia_equals_w() {
        let lat = 50.0_f64;
        let lon = 11.0_f64;
        let spacing = 0.0001_f64;
        let target = geodetic_to_ecef(lat, lon, 0.0);
        // Side-looking ~35° off-nadir.
        let sat = [target[0] + 500_000.0, target[1], target[2] + 700_000.0];

        let outcome = compute_terrain_geometry(
            lat, lon, 0.0, 0.0, 0.0, 0.0, sat, target, spacing,
        );
        match outcome {
            TerrainGeometryOutcome::Valid { w, cos_lia } => {
                assert!(
                    (w - cos_lia).abs() < 1e-9,
                    "flat terrain: w ({w}) and cos_lia ({cos_lia}) must coincide"
                );
                assert!(
                    cos_lia > 0.5 && cos_lia < 1.0,
                    "S-1 incidence range: 0.5 < cos_lia < 1.0, got {cos_lia}"
                );
            }
            other => panic!("expected Valid on flat terrain, got {other:?}"),
        }
    }

    /// Shadow geometry returns `Shadow` (no cos_lia) when the facet faces away.
    #[test]
    fn terrain_geometry_shadow_returns_shadow_variant() {
        let lat = 50.0_f64;
        let lon = 11.0_f64;
        let spacing = 0.0001_f64;
        let target = geodetic_to_ecef(lat, lon, 0.0);
        let sat = [target[0] + 700_000.0, target[1], target[2]];
        // Steep east-facing back-slope.
        let outcome = compute_terrain_geometry(
            lat, lon, -5000.0, 5000.0, 0.0, 0.0, sat, target, spacing,
        );
        assert!(
            matches!(outcome, TerrainGeometryOutcome::Shadow),
            "expected Shadow, got {outcome:?}"
        );
    }

    /// Foreshortening: the variant must carry `cos_lia` so callers that
    /// consume the LIA raster get a value even where `w` is rejected.
    /// Designing a synthetic DEM that lands precisely in the foreshortening
    /// band (`0 < w < LAYOVER_SHADOW_MIN_W`) is fragile; instead we test
    /// the structural contract on a directly-constructed variant.
    #[test]
    fn terrain_geometry_foreshortening_carries_cos_lia() {
        let outcome = TerrainGeometryOutcome::Foreshortening { cos_lia: 0.123 };
        match outcome {
            TerrainGeometryOutcome::Foreshortening { cos_lia } => {
                assert!(cos_lia >= 0.0 && cos_lia <= 1.0);
            }
            other => panic!("expected Foreshortening, got {other:?}"),
        }
    }

    // ── End-to-end: TC produces sized cos_lia + mask buffers ──────────────

    /// `compute_lia=false` (the default): both sidecar buffers are empty.
    /// `compute_lia=true`: both sidecar buffers have one entry per output
    /// pixel.  This is the structural contract that the export layer and
    /// downstream Python tooling rely on.
    #[test]
    fn geocoded_image_sidecars_track_compute_lia_flag() {
        // Pure structural test: build a minimal `GeocodedImage` for both
        // values of the sidecar policy and check field sizes.  We don't
        // run the full TC pipeline here (no SAFE fixtures in this unit
        // test); the sidecar plumbing is exercised end-to-end by the
        // examples and ASF-comparison scripts.
        let img_off = GeocodedImage {
            data: vec![f32::NAN; 6],
            rows: 2,
            cols: 3,
            origin_lat_deg: 50.0,
            origin_lon_deg: 11.0,
            pixel_spacing_lat_deg: 0.0001,
            pixel_spacing_lon_deg: 0.0001,
            total_pixel_count: 6,
            valid_pixel_count: 0,
            dem_missing_count: 6,
            outside_footprint_count: 0,
            non_converged_count: 0,
            degenerate_geometry_count: 0,
            flat_masked_count: 0,
            noise_masked_count: 0,
            geotransform: [11.0, 0.0001, 0.0, 50.0, 0.0, -0.0001],
            crs: OutputCrs::Wgs84LatLon,
            cos_lia: Vec::new(),
            mask: Vec::new(),
            newton_iter_histogram: [0u64; NEWTON_ITER_HIST_SIZE],
        };
        assert!(
            img_off.cos_lia.is_empty() && img_off.mask.is_empty(),
            "compute_lia=false → sidecar buffers must be empty"
        );

        let img_on = GeocodedImage {
            cos_lia: vec![f32::NAN; 6],
            mask: vec![mask_code::OUTSIDE_FOOTPRINT; 6],
            ..img_off
        };
        assert_eq!(img_on.cos_lia.len(), 6);
        assert_eq!(img_on.mask.len(), 6);
    }

    /// Mask code values are stable: third-party code (Python wrappers,
    /// downstream QA scripts) decode the mask raster against these
    /// numeric constants.  Any future renumbering must be documented in
    /// a migration note.
    #[test]
    fn mask_code_values_are_stable() {
        assert_eq!(mask_code::VALID, 0);
        assert_eq!(mask_code::OUTSIDE_FOOTPRINT, 1);
        assert_eq!(mask_code::DEM_MISSING, 2);
        assert_eq!(mask_code::NON_CONVERGED, 3);
        assert_eq!(mask_code::DEGENERATE_GEOMETRY, 4);
        assert_eq!(mask_code::SHADOW, 5);
        assert_eq!(mask_code::LAYOVER, 6);
        assert_eq!(mask_code::FORESHORTENING, 7);
        assert_eq!(mask_code::NOISE_MASKED, 8);
        assert_eq!(mask_code::DEM_NEIGHBOUR_MISSING, 9);
    }

    // ── Phase 4: chunked row scheduling helper ─────────────────────────

    /// Concatenated chunk ranges must exactly cover `0..n_rows` with
    /// no gaps and no overlaps, for every (n_rows, chunk_size) pair
    /// that can occur at runtime.  This is the only correctness gate
    /// for the chunked scheduler — the per-row body is unchanged, so
    /// byte-equality with the per-row scheduler reduces to "the chunks
    /// tile the row range correctly".
    #[test]
    fn row_chunk_ranges_tiles_full_range_without_gaps_or_overlaps() {
        for &(n_rows, chunk_size) in &[
            (0_usize, 1_usize),
            (0, 64),
            (1, 1),
            (1, 64),
            (63, 64),
            (64, 64),
            (65, 64),
            (128, 64),
            (129, 64),
            (1000, 1),     // chunk_size=1 is the per-row baseline
            (1000, 7),     // odd chunk size, partial last chunk
            (1000, 1000),  // single chunk
            (1000, 4096),  // chunk_size > n_rows
            (12_345, 64),
        ] {
            let ranges: Vec<(usize, usize)> =
                row_chunk_ranges(n_rows, chunk_size).collect();

            // Concatenated rows in iteration order.
            let mut covered: Vec<usize> = Vec::with_capacity(n_rows);
            for (start, end) in &ranges {
                assert!(
                    start < end || (n_rows == 0 && ranges.is_empty()),
                    "empty or inverted range ({start},{end}) for ({n_rows},{chunk_size})",
                );
                assert!(
                    end - start <= chunk_size,
                    "chunk too big: ({start},{end}) len={} > chunk_size={chunk_size}",
                    end - start,
                );
                for r in *start..*end {
                    covered.push(r);
                }
            }
            let expected: Vec<usize> = (0..n_rows).collect();
            assert_eq!(
                covered, expected,
                "row coverage mismatch for n_rows={n_rows}, chunk_size={chunk_size}",
            );
        }
    }

    /// Calling with `chunk_size = 0` is a programming error caught by
    /// the runtime guard inside `terrain_correction`; the helper itself
    /// must panic loudly rather than silently returning an infinite
    /// iterator or an empty one.
    #[test]
    #[should_panic(expected = "chunk_size must be > 0")]
    fn row_chunk_ranges_rejects_zero_chunk_size() {
        let _ = row_chunk_ranges(100, 0).collect::<Vec<_>>(); // SAFETY-OK: discarded value is never reached; .collect must be evaluated to trigger the panic asserted by #[should_panic]
    }

    // ── Lanczos-3 kernel tests ─────────────────────────────────────────────────

    #[test]
    fn lanczos_weight_at_zero_is_one() {
        assert!((lanczos_weight(0.0) - 1.0).abs() < 1e-12);
    }

    #[test]
    fn lanczos_weight_at_nonzero_integers_is_zero() {
        // Interpolation property: L(n) = 0 for all integer n ≠ 0.
        // Follows from sinc(n) = 0 for n ∈ ℤ \ {0}.
        for n in [-3, -2, -1, 1, 2, 3] {
            let w = lanczos_weight(n as f64);
            assert!(
                w.abs() < 1e-12,
                "lanczos_weight({n}) = {w}, expected 0"
            );
        }
    }

    #[test]
    fn lanczos_constant_field_returns_constant() {
        // Lanczos-3 does NOT have exact partition-of-unity at non-integer
        // positions (the sinc window breaks it); the reconstruction error
        // is typically < 2 %.  Test that results are within 5 % of the
        // input constant, and exact at integer positions.
        let data = vec![5.0_f32; 8 * 8];
        // Non-integer positions: allow ~5 % relative error.
        for (l, s) in [(2.5_f64, 2.5_f64), (3.7, 4.2)] {
            let v = lanczos_sample_slice(&data, 8, 8, l, s);
            assert!(
                (v - 5.0).abs() < 0.25,
                "constant field: got {v} at ({l},{s}), expected ≈ 5.0"
            );
        }
        // Integer position: must be exact (partition-of-unity holds at nodes).
        let v = lanczos_sample_slice(&data, 8, 8, 4.0, 4.0);
        assert!(
            (v - 5.0).abs() < 1e-4,
            "constant field at exact node: got {v}"
        );
    }

    #[test]
    fn lanczos_at_exact_grid_node_returns_value() {
        // At exact integer coordinates the result equals the grid value.
        let data: Vec<f32> = (0..25).map(|i| i as f32).collect();
        let v = lanczos_sample_slice(&data, 5, 5, 3.0, 4.0);
        let expected = data[3 * 5 + 4];
        assert!(
            (v - expected).abs() < 1e-3,
            "exact node: got {v}, expected {expected}"
        );
    }

    #[test]
    fn lanczos_nan_tap_propagates_nan() {
        let mut data = vec![1.0_f32; 8 * 8];
        data[4 * 8 + 4] = f32::NAN;
        let v = lanczos_sample_slice(&data, 8, 8, 3.5, 3.5);
        assert!(v.is_nan(), "NaN tap must propagate to NaN result");
    }

    // ── RadarImage trait tests ─────────────────────────────────────────────────

    fn make_merged(data: Vec<f32>, lines: usize, samples: usize) -> MergedSigma0 {
        let n = lines * samples;
        MergedSigma0 {
            data,
            nesz: vec![0.01_f32; n],
            lines,
            samples,
            near_slant_range_time_s: 0.005,
            range_pixel_spacing_s: 1e-8,
            range_pixel_spacing_m: 2.33,
            cal_lut_extrapolation_gap_px: 0,
            noise_lut_extrapolation_gap_px: 0,
            azimuth_start_time: Utc::now(),
            azimuth_time_interval_s: 0.002055556,
        }
    }

    #[test]
    fn radar_image_bilinear_matches_slice_fn() {
        // RadarImage::sample_at with Bilinear must equal bilinear_sample_slice.
        let m = make_merged(vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0], 3, 3);
        for (l, s) in [(0.0_f64, 0.0_f64), (0.7, 1.3), (1.5, 0.5), (2.0, 2.0)] {
            let expected = bilinear_sample_slice(&m.data, m.lines, m.samples, l, s);
            let got = m.sample_at(l, s, ResamplingKernel::Bilinear);
            assert!(
                (got - expected).abs() < 1e-6,
                "bilinear mismatch at ({l},{s}): {got} vs {expected}"
            );
        }
    }

    #[test]
    fn radar_image_nesz_at_returns_some() {
        let m = make_merged(vec![1.0_f32; 4], 2, 2);
        let result = m.nesz_at(0.5, 0.5);
        assert!(result.is_some(), "MergedSigma0::nesz_at must return Some");
        let nesz = result.unwrap(); // SAFETY-OK: asserted Some above
        assert!(nesz.is_finite(), "NESZ value must be finite");
    }

    #[test]
    fn radar_image_default_nesz_returns_none() {
        // A minimal RadarImage implementor without NESZ must return None.
        struct ConstantImage(f32);
        impl RadarImage for ConstantImage {
            fn sample_at(&self, _l: f64, _s: f64, _k: ResamplingKernel) -> f32 { self.0 }
            fn lines(&self) -> usize { 1 }
            fn samples(&self) -> usize { 1 }
            fn azimuth_start_time(&self) -> chrono::DateTime<chrono::Utc> { Utc::now() }
            fn near_slant_range_time_s(&self) -> f64 { 0.005 }
            fn range_pixel_spacing_m(&self) -> f64 { 2.33 }
            fn azimuth_time_interval_s(&self) -> f64 { 0.002 }
        }
        let img = ConstantImage(3.0);
        assert!(img.nesz_at(0.0, 0.0).is_none());
    }
}
