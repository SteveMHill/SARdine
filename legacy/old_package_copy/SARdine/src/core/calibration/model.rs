#![allow(dead_code, unused_variables)]
use ndarray::Array2;

use crate::core::calibration::CalibrationType;
use crate::core::global_clamp_debug::ClampDebug;
use crate::core::metadata::{
    CalibrationScale as ParserCalibrationScale, CalibrationUnits as ParserCalibrationUnits,
    CompactCalibrationData as ParserCompactCalibrationData,
    CompactGeometricData as ParserCompactGeometricData, CompactNoiseData as ParserCompactNoiseData,
    CompactTimingData as ParserCompactTimingData,
};
use crate::types::SarResult;

/// Minimum valid denoised power (before calibration). If denoised power falls below this
/// (noise > signal), the pixel is invalid.
pub const MIN_VALID_POWER: f32 = 1e-10;

/// Guard-rail threshold used when reporting extreme magnitudes.
pub const POWER_ALERT_THRESHOLD: f32 = 4.0e9;

/// Valid sample range per azimuth line (for burst edge handling). Used to clip LUT access
/// to avoid garbage coefficients at invalid samples.
///
/// # Bound Semantics
/// **INCLUSIVE bounds**: `(first_valid, last_valid)` where both are valid indices.
/// This differs from most other bounds in this codebase which use exclusive end.
///
/// | Element | Convention | Description |
/// |---------|------------|-------------|
/// | `first_valid` | inclusive | First valid sample index |
/// | `last_valid` | **inclusive** | Last valid sample index (NOT one-past-end) |
///
/// # Example
/// ```ignore
/// // Line with valid samples 100..=24999 (24900 valid samples):
/// ranges[line] = (100, 24999);  // 24999 IS a valid index
///
/// // Iteration pattern:
/// for sample in ranges[line].0..=ranges[line].1 {
///     // process valid sample
/// }
/// ```
///
/// # Why Inclusive?
/// Calibration LUT access patterns check `sample <= last_valid`, making inclusive
/// bounds more natural for boundary checks.
#[derive(Debug, Clone)]
pub struct ValidSampleRanges {
    /// Per-line valid sample range: `(first_valid, last_valid)` **both inclusive**.
    pub ranges: Vec<(usize, usize)>,
}

impl ValidSampleRanges {
    /// Create from first/last valid sample arrays.
    ///
    /// # Note on Bound Semantics
    /// This struct uses **INCLUSIVE** bounds where both `first` and `last` are valid indices.
    /// Iteration should use `first..=last`, NOT `first..last`.
    pub fn from_arrays(first_valid: &[i32], last_valid: &[i32]) -> Self {
        let ranges = first_valid
            .iter()
            .zip(last_valid.iter())
            .map(|(first, last)| {
                let start = (*first).max(0) as usize;
                let end = (*last).max(0) as usize;

                // Validate that bounds make sense after clamping.
                // Both being negative (e.g., -1/-1) means the entire line is outside the
                // valid burst window; treat as zero-width by using a reversed pair (1, 0)
                // so that every `col >= start && col <= end` check returns false.
                if *first < 0 && *last < 0 {
                    return (1, 0); // Sentinel: no valid pixels on this line
                }
                // In all other cases end should be >= start. Log a warning in release
                // builds so annotation corruption surfaces without panicking.
                if end < start {
                    log::warn!(
                        "ValidSampleRanges: end ({}) < start ({}) after clamping \
                         (annotation first={}, last={}); treating line as all-valid",
                        end, start, first, last
                    );
                    return (0, end);
                }

                (start, end)
            })
            .collect();
        Self { ranges }
    }

    /// Create "all valid" range (no clipping) for non-burst data.
    pub fn all_valid(num_lines: usize, num_samples: usize) -> Self {
        let ranges = vec![(0, num_samples - 1); num_lines];
        Self { ranges }
    }
}

/// Clamp calibration gains to sane bounds to prevent astronomical outputs.
#[inline]
pub fn sane_gain(g: f32) -> f32 {
    if !g.is_finite() {
        log::error!("⚠️  Non-finite calibration coefficient encountered: {}", g);
        // Scientific requirement (Jan 2026): never silently convert invalid
        // calibration gains into zeros. Downstream code treats NaN as
        // "no data" and diagnostics will surface the fraction of
        // non-finite outputs.
        return f32::NAN;
    }

    if g <= 0.0 {
        log::warn!(
            "⚠️  Non-positive calibration coefficient encountered ({}); treating as invalid gain",
            g
        );
        return f32::NAN;
    }

    g
}

/// Calibration vector from Sentinel-1 XML.
#[derive(Debug, Clone)]
pub struct CalibrationVector {
    pub azimuth_time: String,
    pub line: i32, // Can be negative for pre-burst calibration data
    pub pixels: Vec<usize>,
    pub sigma_nought: Vec<f32>,
    pub beta_nought: Vec<f32>,
    pub gamma: Vec<f32>,
    pub dn: Vec<f32>,
    // Flatness detection flags for robust calibration
    pub beta_flat: bool,
    pub sigma_flat: bool,
    pub gamma_flat: bool,
}

impl CalibrationVector {
    /// Recover corrupted calibration coefficients using physical relationships.
    /// Returns true if any corruption was detected and recovered.
    pub fn recover_corrupted_coefficients(
        &mut self,
        min_incidence_deg: f32,
        max_incidence_deg: f32,
        swath_width: usize,
    ) -> bool {
        use super::parsing::{
            derive_beta_from_sigma, derive_gamma_from_sigma, derive_sigma_from_beta,
        };

        let mut recovered = false;

        // Priority 1: Recover Beta0 from Sigma0 (most common corruption in archived data)
        if self.beta_flat && !self.sigma_flat && !self.sigma_nought.is_empty() {
            log::warn!(
                "🔧 Recovering corrupted Beta0 from Sigma0 at line {} ({} pixels)",
                self.line,
                self.pixels.len()
            );
            self.beta_nought = derive_beta_from_sigma(
                &self.sigma_nought,
                &self.pixels,
                min_incidence_deg,
                max_incidence_deg,
                swath_width,
            );
            self.beta_flat = false;
            recovered = true;
        }

        // Priority 2: Recover Gamma0 from Sigma0
        if self.gamma_flat && !self.sigma_flat && !self.sigma_nought.is_empty() {
            log::warn!(
                "🔧 Recovering corrupted Gamma0 from Sigma0 at line {} ({} pixels)",
                self.line,
                self.pixels.len()
            );
            self.gamma = derive_gamma_from_sigma(
                &self.sigma_nought,
                &self.pixels,
                min_incidence_deg,
                max_incidence_deg,
                swath_width,
            );
            self.gamma_flat = false;
            recovered = true;
        }

        // Priority 3: Recover Sigma0 from Beta0 (rare, but possible)
        if self.sigma_flat && !self.beta_flat && !self.beta_nought.is_empty() {
            log::warn!(
                "🔧 Recovering corrupted Sigma0 from Beta0 at line {} ({} pixels)",
                self.line,
                self.pixels.len()
            );
            self.sigma_nought = derive_sigma_from_beta(
                &self.beta_nought,
                &self.pixels,
                min_incidence_deg,
                max_incidence_deg,
                swath_width,
            );
            self.sigma_flat = false;
            recovered = true;
        }

        if recovered {
            log::info!(
                "✅ Calibration recovery complete at line {} (beta_flat={}, sigma_flat={}, gamma_flat={})",
                self.line,
                self.beta_flat,
                self.sigma_flat,
                self.gamma_flat
            );
        }

        recovered
    }
}

/// Thermal noise vector from Sentinel-1 noise XML (range component).
#[derive(Debug, Clone)]
pub struct NoiseVector {
    pub azimuth_time: String,
    pub azimuth_time_seconds: f64, // Parsed azimuth time in seconds since reference
    pub line: f64,
    pub range_pixels: Vec<f64>,
    pub noise_range_lut: Vec<f32>,
}

/// Azimuth noise vector from Sentinel-1 noise XML (IPF 3.x+).
/// Represents along-track noise variation which is typically smaller than range noise
/// but can be significant for specific applications.
#[derive(Debug, Clone)]
pub struct NoiseAzimuthVector {
    pub swath: String,
    pub first_azimuth_line: usize,
    pub first_range_sample: usize,
    pub last_azimuth_line: usize,
    pub last_range_sample: usize,
    pub lines: Vec<usize>,
    pub noise_azimuth_lut: Vec<f32>,
}

/// Antenna pattern vector from annotation XML.
#[derive(Debug, Clone)]
pub struct AntennaPatternVector {
    pub azimuth_time: String,
    pub line: i32,
    pub pixels: Vec<usize>,
    pub values: Vec<f64>,
}

/// Incidence angle model for deriving β⁰/γ⁰ from σ⁰.
pub trait IncidenceAngleModel: Send + Sync + std::fmt::Debug {
    fn alpha_ellipsoid(&self, line_slc: i32, col_slc: usize) -> SarResult<f32>;
    fn clone_model(&self) -> Box<dyn IncidenceAngleModel>;
}

/// Simple ellipsoid incidence angle model derived from annotation metadata.
#[derive(Debug, Clone)]
pub struct EllipsoidIncidenceModel {
    pub min_incidence_rad: f32,
    pub max_incidence_rad: f32,
    pub swath_width: usize,
}

impl EllipsoidIncidenceModel {
    pub fn new(min_incidence_deg: f32, max_incidence_deg: f32, swath_width: usize) -> Self {
        Self {
            min_incidence_rad: min_incidence_deg.to_radians(),
            max_incidence_rad: max_incidence_deg.to_radians(),
            swath_width,
        }
    }
}

impl IncidenceAngleModel for EllipsoidIncidenceModel {
    fn alpha_ellipsoid(&self, _line_slc: i32, col_slc: usize) -> SarResult<f32> {
        let range_fraction = (col_slc as f32) / (self.swath_width as f32).max(1.0);
        let range_fraction = range_fraction.dbg_clamp(0.0, 1.0, "calibrate_range_fraction");

        let incidence_angle = self.min_incidence_rad
            + range_fraction * (self.max_incidence_rad - self.min_incidence_rad);

        Ok(incidence_angle)
    }

    fn clone_model(&self) -> Box<dyn IncidenceAngleModel> {
        Box::new(self.clone())
    }
}

/// Pre-computed calibration lookup table for fast access.
#[derive(Debug, Clone)]
pub struct CalibrationLUT {
    pub sigma_values: Array2<f32>,
    pub beta_values: Array2<f32>,
    pub gamma_values: Array2<f32>,
    pub dn_values: Array2<f32>,
    pub is_precomputed: bool,
}

impl CalibrationLUT {
    pub fn new(dims: (usize, usize)) -> Self {
        let (height, width) = dims;
        Self {
            sigma_values: Array2::zeros((height, width)),
            beta_values: Array2::zeros((height, width)),
            gamma_values: Array2::zeros((height, width)),
            dn_values: Array2::ones((height, width)),
            is_precomputed: false,
        }
    }
}

/// Separable calibration model: K(line, pixel) ≈ A(line) × R(pixel).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum NoiseLutMode {
    Full2D,
    RangeOnly,
    AzimuthInterpolated,
    Disabled,
}

/// Pre-computed thermal noise lookup table for denoising. Some modes keep the
/// dense grid; lightweight modes store compressed representations instead.
#[derive(Debug, Clone)]
pub struct NoiseLUT {
    pub noise_values: Array2<f32>,
    pub(crate) range_profile: Option<Vec<f32>>, // Range-only model
    pub(crate) precomputed_rows: Option<Vec<Vec<f32>>>, // Per-vector range rows
    pub(crate) azimuth_brackets: Option<Vec<AzimuthBracketCache>>, // Cached azimuth weights
    pub azimuth_axis: Vec<f64>,
    pub range_axis: Vec<f64>,
    pub height: usize,
    pub width: usize,
    pub mode: NoiseLutMode,
    pub is_precomputed: bool,
}

impl NoiseLUT {
    pub fn new(dims: (usize, usize)) -> Self {
        let (height, width) = dims;
        Self {
            noise_values: Array2::zeros((height, width)),
            range_profile: None,
            precomputed_rows: None,
            azimuth_brackets: None,
            azimuth_axis: Vec::new(),
            range_axis: Vec::new(),
            height,
            width,
            mode: NoiseLutMode::Full2D,
            is_precomputed: false,
        }
    }

    pub fn disabled(dims: (usize, usize)) -> Self {
        let (height, width) = dims;
        Self {
            noise_values: Array2::zeros((0, 0)),
            range_profile: None,
            precomputed_rows: None,
            azimuth_brackets: None,
            azimuth_axis: vec![],
            range_axis: vec![],
            height,
            width,
            mode: NoiseLutMode::Disabled,
            is_precomputed: true,
        }
    }
}

/// Antenna pattern correction LUT for removing azimuth scalloping.
#[derive(Debug, Clone)]
pub struct AntennaPatternLUT {
    pub pattern_values: Array2<f32>,
    pub is_precomputed: bool,
}

impl AntennaPatternLUT {
    pub fn new(dims: (usize, usize)) -> Self {
        let (height, width) = dims;
        Self {
            pattern_values: Array2::ones((height, width)),
            is_precomputed: false,
        }
    }

    pub fn unity(dims: (usize, usize)) -> Self {
        let (height, width) = dims;
        Self {
            pattern_values: Array2::ones((height, width)),
            is_precomputed: true,
        }
    }
}

/// Calibration coefficients for radiometric correction.
#[derive(Debug)]
pub struct CalibrationCoefficients {
    pub vectors: Vec<CalibrationVector>,
    pub swath: String,
    pub polarization: String,
    pub product_first_line_utc_time: String,
    pub product_last_line_utc_time: String,
    pub lut: Option<CalibrationLUT>,
    pub coordinate_mapper: Option<CalibrationCoordinateMapper>,
    pub incidence_model: Option<Box<dyn IncidenceAngleModel>>,
    pub incidence_lut_max_grid: Option<usize>,
    pub abs_const: f64,
    pub antenna_pattern_lut: Option<AntennaPatternLUT>,
    pub antenna_pattern_vectors: Vec<AntennaPatternVector>,
    pub valid_sample_ranges: Option<ValidSampleRanges>,
}

impl CalibrationCoefficients {
    pub fn new() -> Self {
        Self {
            vectors: Vec::new(),
            swath: String::new(),
            polarization: String::new(),
            product_first_line_utc_time: String::new(),
            product_last_line_utc_time: String::new(),
            lut: None,
            coordinate_mapper: None,
            incidence_model: None,
            incidence_lut_max_grid: None,
            abs_const: 1.0,
            antenna_pattern_lut: None,
            antenna_pattern_vectors: Vec::new(),
            valid_sample_ranges: None,
        }
    }

    pub fn precompute_lut(&mut self, image_dims: (usize, usize)) -> SarResult<()> {
        crate::core::calibration::precompute_calibration_lut(self, image_dims)
    }

    pub fn get_calibration_value(
        &self,
        line: i32,
        pixel: usize,
        cal_type: CalibrationType,
    ) -> SarResult<f32> {
        crate::core::calibration::get_calibration_value(self, line, pixel, cal_type)
    }

    pub fn precompute_antenna_pattern_lut(&mut self, image_dims: (usize, usize)) -> SarResult<()> {
        crate::core::calibration::precompute_antenna_pattern_lut(self, image_dims)
    }

    pub fn parse_and_set_antenna_pattern(
        &mut self,
        xml: &str,
        range_sampling_rate_hz: Option<f64>,
    ) -> SarResult<()> {
        let patterns =
            crate::core::calibration::parse_antenna_pattern_from_xml(xml, range_sampling_rate_hz)?;
        self.antenna_pattern_vectors = patterns;
        Ok(())
    }

    pub fn set_incidence_model(&mut self, model: Box<dyn IncidenceAngleModel>) -> SarResult<()> {
        self.incidence_model = Some(model);
        Ok(())
    }

    pub fn set_incidence_lut_max_grid(&mut self, max_grid: usize) {
        self.incidence_lut_max_grid = Some(max_grid.max(1));
    }

    pub fn set_coordinate_mapper(&mut self, mapper: CalibrationCoordinateMapper) -> SarResult<()> {
        self.coordinate_mapper = Some(mapper);
        Ok(())
    }

    pub fn create_auto_coordinate_mapper(
        &self,
        burst_start_line: i32,
        image_start_line: i32,
        image_width: usize,
    ) -> SarResult<CalibrationCoordinateMapper> {
        let line_offset = image_start_line - burst_start_line;
        Ok(CalibrationCoordinateMapper {
            burst_start_line,
            image_start_line,
            line_offset,
            slc_range_offset: 0.0,
            slc_range_scale: 1.0,
        })
    }
}

impl Clone for CalibrationCoefficients {
    fn clone(&self) -> Self {
        Self {
            vectors: self.vectors.clone(),
            swath: self.swath.clone(),
            polarization: self.polarization.clone(),
            product_first_line_utc_time: self.product_first_line_utc_time.clone(),
            product_last_line_utc_time: self.product_last_line_utc_time.clone(),
            lut: self.lut.clone(),
            coordinate_mapper: self.coordinate_mapper.clone(),
            incidence_model: self.incidence_model.as_ref().map(|m| m.clone_model()),
            incidence_lut_max_grid: self.incidence_lut_max_grid,
            abs_const: self.abs_const,
            antenna_pattern_lut: self.antenna_pattern_lut.clone(),
            antenna_pattern_vectors: self.antenna_pattern_vectors.clone(),
            valid_sample_ranges: self.valid_sample_ranges.clone(),
        }
    }
}

impl Default for CalibrationCoefficients {
    fn default() -> Self {
        Self::new()
    }
}

/// Coordinate mapper for transforming post-processing coordinates back to original SLC grid.
#[derive(Debug, Clone)]
pub struct CalibrationCoordinateMapper {
    pub burst_start_line: i32,
    pub image_start_line: i32,
    pub line_offset: i32,
    pub slc_range_offset: f64,
    pub slc_range_scale: f64,
}

impl CalibrationCoordinateMapper {
    pub fn map_to_slc_coordinates(&self, line: i32, pixel: usize) -> (i32, usize) {
        let slc_line = line + self.line_offset;
        let slc_pixel = self.slc_range_offset + (pixel as f64) * self.slc_range_scale;
        (slc_line, slc_pixel.round() as usize)
    }
}

/// Coordinate mapper for transforming (row, col) → (azimuth_coord, range_coord) for noise LUTs.
#[derive(Debug, Clone)]
pub struct NoiseCoordinateMapper {
    pub burst_start_line: f64,
    pub burst_start_time: f64,
    pub use_time_axis: bool,
}

impl NoiseCoordinateMapper {
    pub fn new(burst_start_line: f64, burst_start_time: f64, use_time_axis: bool) -> Self {
        Self {
            burst_start_line,
            burst_start_time,
            use_time_axis,
        }
    }

    pub fn map_coordinates(&self, full_azimuth: f64, full_range: f64) -> (f64, f64) {
        // When use_time_axis is true, full_azimuth is expected to be in time (seconds)
        // and we subtract burst_start_time. Otherwise, it's a line index and we
        // subtract burst_start_line.
        let azimuth_coord = if self.use_time_axis {
            full_azimuth - self.burst_start_time
        } else {
            full_azimuth - self.burst_start_line
        };
        (azimuth_coord, full_range)
    }
}

/// Thermal noise coefficients for noise removal.
#[derive(Debug, Clone)]
pub struct NoiseCoefficients {
    pub vectors: Vec<NoiseVector>,
    /// Azimuth noise vectors (IPF 3.x+ products only, may be empty for older products)
    pub azimuth_vectors: Vec<NoiseAzimuthVector>,
    pub swath: String,
    pub polarization: String,
    pub product_first_line_utc_time: String,
    pub product_last_line_utc_time: String,
    pub lut: Option<NoiseLUT>,
    pub coordinate_mapper: Option<NoiseCoordinateMapper>,
}

impl NoiseCoefficients {
    pub fn new() -> Self {
        Self {
            vectors: Vec::new(),
            azimuth_vectors: Vec::new(),
            swath: String::new(),
            polarization: String::new(),
            product_first_line_utc_time: String::new(),
            product_last_line_utc_time: String::new(),
            lut: None,
            coordinate_mapper: None,
        }
    }

    pub fn precompute_lut(&mut self, image_dims: (usize, usize)) -> SarResult<()> {
        crate::core::calibration::precompute_noise_lut(self, image_dims)
    }

    pub fn validate_vectors(&self) -> SarResult<()> {
        Ok(())
    }
}

#[derive(Clone, Copy, Debug)]
pub(crate) struct AzimuthBracketCache {
    pub(crate) lower: usize,
    pub(crate) upper: usize,
    pub(crate) weight: f32,
}

/// Pre-computed interpolation cache for reuse across multiple polarizations.
#[derive(Debug, Clone)]
pub struct SharedCalibrationCache {
    pub image_dims: (usize, usize),
    pub subswath_id: String,
    pub(crate) azimuth_brackets: Vec<AzimuthBracketCache>,
    pub slc_lines: Vec<i32>,
    pub slc_pixels: Vec<usize>,
    pub created_at: std::time::Instant,
}

impl SharedCalibrationCache {
    pub fn new(
        coefficients: &CalibrationCoefficients,
        image_dims: (usize, usize),
        subswath_id: String,
    ) -> SarResult<Self> {
        let (height, width) = image_dims;

        let slc_lines: Vec<i32> = (0..height)
            .map(|row| {
                if let Some(ref mapper) = coefficients.coordinate_mapper {
                    mapper.map_to_slc_coordinates(row as i32, 0).0
                } else {
                    row as i32
                }
            })
            .collect();

        let slc_pixels: Vec<usize> = (0..width)
            .map(|col| {
                if let Some(ref mapper) = coefficients.coordinate_mapper {
                    mapper.map_to_slc_coordinates(0, col).1
                } else {
                    col
                }
            })
            .collect();

        let azimuth_brackets = slc_lines
            .iter()
            .map(|&slc_line| {
                let mut lower = 0usize;
                let mut upper = coefficients.vectors.len() - 1;

                for (idx, vector) in coefficients.vectors.iter().enumerate() {
                    if vector.line <= slc_line {
                        lower = idx;
                    }
                    if vector.line >= slc_line {
                        upper = idx;
                        break;
                    }
                }

                if lower == upper {
                    return AzimuthBracketCache {
                        lower,
                        upper,
                        weight: 0.0,
                    };
                }

                let line1 = coefficients.vectors[lower].line as f32;
                let line2 = coefficients.vectors[upper].line as f32;
                let weight = if (line2 - line1).abs() > f32::EPSILON {
                    (slc_line as f32 - line1) / (line2 - line1)
                } else {
                    0.0
                };

                AzimuthBracketCache {
                    lower,
                    upper,
                    weight: weight.dbg_clamp(0.0, 1.0, "azimuth_bracket_weight"),
                }
            })
            .collect::<Vec<_>>();

        Ok(Self {
            image_dims,
            subswath_id,
            azimuth_brackets,
            slc_lines,
            slc_pixels,
            created_at: std::time::Instant::now(),
        })
    }
}

/// Type aliases to bridge existing parsed data structures into the new calibration module layout.
pub type CompactCalibrationData = ParserCompactCalibrationData;
pub type CompactNoiseData = ParserCompactNoiseData;
pub type CompactGeometricData = ParserCompactGeometricData;
pub type CompactTimingData = ParserCompactTimingData;
pub type CalibrationUnits = ParserCalibrationUnits;
pub type CalibrationScale = ParserCalibrationScale;
