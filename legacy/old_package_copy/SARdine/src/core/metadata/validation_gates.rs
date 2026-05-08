#![allow(dead_code, unused_variables)]
//! # SARdine Validation Gates
//!
//! Comprehensive validation framework enforcing critical processing quality checks:
//! - Power preservation in deburst/merge operations
//! - Uncovered pixel rate validation
//! - NESZ (Noise Equivalent Sigma Zero) compliance
//! - Point Target Analysis (PTA) validation
//! - Geometry accuracy validation
//! - LUT domain boundary checks

use crate::types::{SarError, SarResult};
use log::{error, info, warn};
use ndarray::Array2;

/// Validation gate results with detailed diagnostics
#[derive(Debug, Clone)]
pub struct ValidationResult {
    pub passed: bool,
    pub test_name: String,
    pub measured_value: f64,
    pub threshold: f64,
    pub tolerance: f64,
    pub diagnostic_message: String,
    pub recommendations: Vec<String>,
}

impl ValidationResult {
    pub fn pass(test_name: &str, measured: f64, threshold: f64, message: &str) -> Self {
        Self {
            passed: true,
            test_name: test_name.to_string(),
            measured_value: measured,
            threshold,
            tolerance: 0.0,
            diagnostic_message: message.to_string(),
            recommendations: Vec::new(),
        }
    }

    pub fn fail(
        test_name: &str,
        measured: f64,
        threshold: f64,
        tolerance: f64,
        message: &str,
    ) -> Self {
        Self {
            passed: false,
            test_name: test_name.to_string(),
            measured_value: measured,
            threshold,
            tolerance,
            diagnostic_message: message.to_string(),
            recommendations: Vec::new(),
        }
    }

    pub fn with_recommendations(mut self, recommendations: Vec<String>) -> Self {
        self.recommendations = recommendations;
        self
    }
}

/// Comprehensive validation gates framework
pub struct ValidationGates {
    /// Enable/disable specific validation checks
    pub enable_power_preservation: bool,
    pub enable_uncovered_pixel_check: bool,
    pub enable_nesz_check: bool,
    pub enable_pta_check: bool,
    pub enable_geometry_check: bool,
    pub enable_lut_domain_check: bool,

    /// Validation thresholds
    pub power_preservation_tolerance: f64, // ±1% default
    pub uncovered_pixel_threshold: f64, // 0.2% default
    pub nesz_tolerance_db: f64,         // ±1-2 dB default
    pub geometry_ce90_threshold: f64,   // 10m default

    /// Point Target Analysis specs for Sentinel-1
    pub pslr_threshold_db: f64, // Peak Sidelobe Ratio
    pub islr_threshold_db: f64,          // Integrated Sidelobe Ratio
    pub irw_threshold_m: f64,            // Impulse Response Width
    pub radiometric_gain_tolerance: f64, // Radiometric gain accuracy
}

impl Default for ValidationGates {
    fn default() -> Self {
        Self {
            enable_power_preservation: true,
            enable_uncovered_pixel_check: true,
            enable_nesz_check: true,
            enable_pta_check: true,
            enable_geometry_check: true,
            enable_lut_domain_check: true,

            power_preservation_tolerance: 0.01, // ±1%
            uncovered_pixel_threshold: 0.002,   // 0.2%
            nesz_tolerance_db: 2.0,             // ±2 dB
            geometry_ce90_threshold: 10.0,      // 10m CE90

            // Sentinel-1 Point Target Analysis specifications
            pslr_threshold_db: -13.0,        // PSLR better than -13 dB
            islr_threshold_db: -10.0,        // ISLR better than -10 dB
            irw_threshold_m: 20.0,           // IRW within 20m
            radiometric_gain_tolerance: 0.5, // ±0.5 dB radiometric accuracy
        }
    }
}

impl ValidationGates {
    /// Create validation gates with custom thresholds
    pub fn with_thresholds(
        power_tolerance: f64,
        uncovered_threshold: f64,
        nesz_tolerance: f64,
        geometry_threshold: f64,
    ) -> Self {
        Self {
            power_preservation_tolerance: power_tolerance,
            uncovered_pixel_threshold: uncovered_threshold,
            nesz_tolerance_db: nesz_tolerance,
            geometry_ce90_threshold: geometry_threshold,
            ..Default::default()
        }
    }

    /// 1. Power preservation test after deburst/merge
    /// Validates that mean(|S|²) in overlap windows pre/post match within ±1%
    pub fn validate_power_preservation(
        &self,
        pre_deburst_power: &Array2<f32>,
        post_deburst_power: &Array2<f32>,
        overlap_regions: &[(usize, usize, usize, usize)], // (row_start, row_end, col_start, col_end)
    ) -> SarResult<ValidationResult> {
        if !self.enable_power_preservation {
            return Ok(ValidationResult::pass(
                "Power Preservation",
                0.0,
                0.0,
                "Validation disabled",
            ));
        }

        let mut pre_power_sum = 0.0;
        let mut post_power_sum = 0.0;
        let mut total_pixels = 0;

        for &(row_start, row_end, col_start, col_end) in overlap_regions {
            for row in row_start..row_end {
                for col in col_start..col_end {
                    if row < pre_deburst_power.nrows()
                        && col < pre_deburst_power.ncols()
                        && row < post_deburst_power.nrows()
                        && col < post_deburst_power.ncols()
                    {
                        pre_power_sum += pre_deburst_power[[row, col]] as f64;
                        post_power_sum += post_deburst_power[[row, col]] as f64;
                        total_pixels += 1;
                    }
                }
            }
        }

        if total_pixels == 0 {
            return Err(SarError::InvalidParameter(
                "No valid overlap pixels found for power preservation test".to_string(),
            ));
        }

        let pre_mean_power = pre_power_sum / total_pixels as f64;
        let post_mean_power = post_power_sum / total_pixels as f64;

        // Add epsilon to protect against division by zero / denormals
        let eps = 1e-12;
        let power_ratio = post_mean_power / (pre_mean_power + eps);
        let power_error = (power_ratio - 1.0).abs();

        info!(
            "🔋 Power Preservation: pre={:.6e}, post={:.6e}, ratio={:.6}, error={:.4}%",
            pre_mean_power,
            post_mean_power,
            power_ratio,
            power_error * 100.0
        );

        if power_error <= self.power_preservation_tolerance {
            Ok(ValidationResult::pass(
                "Power Preservation",
                power_error,
                self.power_preservation_tolerance,
                &format!(
                    "Power preserved within ±{:.1}% tolerance (measured: {:.4}%)",
                    self.power_preservation_tolerance * 100.0,
                    power_error * 100.0
                ),
            ))
        } else {
            Ok(ValidationResult::fail(
                "Power Preservation",
                power_error,
                self.power_preservation_tolerance,
                self.power_preservation_tolerance,
                &format!(
                    "Power loss detected: {:.4}% exceeds ±{:.1}% tolerance",
                    power_error * 100.0,
                    self.power_preservation_tolerance * 100.0
                ),
            )
            .with_recommendations(vec![
                "Check deburst weighting functions for energy normalization".to_string(),
                "Verify overlap region blending algorithms".to_string(),
                "Inspect burst boundary handling for power leakage".to_string(),
            ]))
        }
    }

    /// 2. Uncovered pixel rate after deburst
    /// Validates that <0.2% of pixels remain uncovered after deburst
    pub fn validate_uncovered_pixels(
        &self,
        debursted_image: &Array2<f32>,
        valid_mask: Option<&Array2<bool>>,
    ) -> SarResult<ValidationResult> {
        if !self.enable_uncovered_pixel_check {
            return Ok(ValidationResult::pass(
                "Uncovered Pixels",
                0.0,
                0.0,
                "Validation disabled",
            ));
        }

        let total_pixels = debursted_image.len();
        let mut uncovered_pixels = 0;

        // Count uncovered pixels (NaN, zero, or invalid)
        for (i, &pixel_value) in debursted_image.iter().enumerate() {
            let is_uncovered =
                pixel_value.is_nan() || pixel_value == 0.0 || pixel_value.is_infinite();

            // Check against valid mask if provided
            let mask_invalid = if let Some(mask) = valid_mask {
                let row = i / debursted_image.ncols();
                let col = i % debursted_image.ncols();
                !mask[[row, col]]
            } else {
                false
            };

            if is_uncovered || mask_invalid {
                uncovered_pixels += 1;
            }
        }

        let uncovered_rate = uncovered_pixels as f64 / total_pixels as f64;

        info!(
            "🕳️  Uncovered Pixels: {}/{} ({:.4}% rate)",
            uncovered_pixels,
            total_pixels,
            uncovered_rate * 100.0
        );

        if uncovered_rate <= self.uncovered_pixel_threshold {
            Ok(ValidationResult::pass(
                "Uncovered Pixels",
                uncovered_rate,
                self.uncovered_pixel_threshold,
                &format!(
                    "Uncovered pixel rate {:.4}% within {:.1}% threshold",
                    uncovered_rate * 100.0,
                    self.uncovered_pixel_threshold * 100.0
                ),
            ))
        } else {
            Ok(ValidationResult::fail(
                "Uncovered Pixels",
                uncovered_rate,
                self.uncovered_pixel_threshold,
                self.uncovered_pixel_threshold,
                &format!(
                    "High uncovered pixel rate: {:.4}% exceeds {:.1}% threshold",
                    uncovered_rate * 100.0,
                    self.uncovered_pixel_threshold * 100.0
                ),
            )
            .with_recommendations(vec![
                "Check burst timing parameters for gap alignment".to_string(),
                "Verify azimuth antenna pattern overlap regions".to_string(),
                "Inspect TOPS acquisition parameters for coverage".to_string(),
                "Review deburst algorithm for edge handling".to_string(),
            ]))
        }
    }

    /// 3. NESZ (Noise Equivalent Sigma Zero) check
    /// Validates far-range ocean σ⁰ close to mission NESZ ±1–2 dB
    pub fn validate_nesz_compliance(
        &self,
        sigma0_image: &Array2<f32>,
        range_pixel_spacing: f64,
        nesz_lut: &[f32], // Mission NESZ lookup table
        ocean_mask: Option<&Array2<bool>>,
    ) -> SarResult<ValidationResult> {
        if !self.enable_nesz_check {
            return Ok(ValidationResult::pass(
                "NESZ Compliance",
                0.0,
                0.0,
                "Validation disabled",
            ));
        }
        // Heuristic dB domain assertions (sigma0 and NESZ LUT should be dB values here)
        if !sigma0_image.is_empty() {
            let suspicious_sigma = sigma0_image
                .iter()
                .take(32)
                .filter(|v| v.is_finite())
                .any(|v| *v > 100.0);
            if suspicious_sigma {
                warn!("⚠️ NESZ validation: sigma0_image appears to contain linear values (>100). Expected dB domain.");
            }
        }
        if !nesz_lut.is_empty() {
            let suspicious_nesz = nesz_lut.iter().take(32).any(|v| *v > 100.0);
            if suspicious_nesz {
                warn!("⚠️ NESZ validation: nesz_lut appears to contain linear values (>100). Expected dB domain.");
            }
        }

        let ncols = sigma0_image.ncols();
        let far_range_start = (ncols as f64 * 0.8) as usize; // Far range (last 20%)

        let mut ocean_sigma0_sum = 0.0;
        let mut ocean_pixel_count = 0;
        let mut nesz_sum = 0.0;

        // Sample far-range ocean pixels
        for row in 0..sigma0_image.nrows() {
            for col in far_range_start..ncols {
                let is_ocean = if let Some(mask) = ocean_mask {
                    mask[[row, col]]
                } else {
                    // Simple heuristic: low backscatter likely ocean
                    sigma0_image[[row, col]] < -15.0 // dB threshold for ocean
                };

                if is_ocean && !sigma0_image[[row, col]].is_nan() {
                    ocean_sigma0_sum += sigma0_image[[row, col]] as f64;
                    ocean_pixel_count += 1;

                    // Get NESZ value for this range
                    if col < nesz_lut.len() {
                        nesz_sum += nesz_lut[col] as f64;
                    }
                }
            }
        }

        if ocean_pixel_count == 0 {
            return Ok(ValidationResult::fail(
                "NESZ Compliance",
                0.0,
                self.nesz_tolerance_db,
                self.nesz_tolerance_db,
                "No ocean pixels found in far range for NESZ validation",
            )
            .with_recommendations(vec![
                "Verify ocean mask or use different validation area".to_string(),
                "Check sigma0 calibration for reasonable ocean values".to_string(),
            ]));
        }

        let mean_ocean_sigma0 = ocean_sigma0_sum / ocean_pixel_count as f64;
        let mean_nesz = nesz_sum / ocean_pixel_count as f64;
        let nesz_difference = (mean_ocean_sigma0 - mean_nesz).abs();

        info!(
            "🌊 NESZ Check: ocean_σ⁰={:.2} dB, mission_NESZ={:.2} dB, diff={:.2} dB",
            mean_ocean_sigma0, mean_nesz, nesz_difference
        );

        if nesz_difference <= self.nesz_tolerance_db {
            Ok(ValidationResult::pass(
                "NESZ Compliance",
                nesz_difference,
                self.nesz_tolerance_db,
                &format!(
                    "Ocean σ⁰ within ±{:.1} dB of mission NESZ (diff: {:.2} dB)",
                    self.nesz_tolerance_db, nesz_difference
                ),
            ))
        } else {
            Ok(ValidationResult::fail(
                "NESZ Compliance",
                nesz_difference,
                self.nesz_tolerance_db,
                self.nesz_tolerance_db,
                &format!(
                    "Ocean σ⁰ deviation {:.2} dB exceeds ±{:.1} dB NESZ tolerance",
                    nesz_difference, self.nesz_tolerance_db
                ),
            )
            .with_recommendations(vec![
                "Check radiometric calibration accuracy".to_string(),
                "Verify noise subtraction in calibration".to_string(),
                "Inspect antenna pattern correction".to_string(),
                "Review thermal noise removal algorithm".to_string(),
            ]))
        }
    }

    /// 4. Point Target Analysis (PTA) validation
    /// Validates PSLR, ISLR, IRW, and radiometric gain within S-1 specs
    pub fn validate_point_target_analysis(
        &self,
        pslr_db: f64,
        islr_db: f64,
        irw_meters: f64,
        radiometric_gain_error_db: f64,
    ) -> SarResult<Vec<ValidationResult>> {
        let mut results = Vec::new();

        if !self.enable_pta_check {
            results.push(ValidationResult::pass(
                "Point Target Analysis",
                0.0,
                0.0,
                "PTA validation disabled",
            ));
            return Ok(results);
        }

        // PSLR (Peak Sidelobe Ratio) check
        if pslr_db <= self.pslr_threshold_db {
            results.push(ValidationResult::pass(
                "PSLR",
                pslr_db,
                self.pslr_threshold_db,
                &format!(
                    "PSLR {:.2} dB meets requirement (<{:.1} dB)",
                    pslr_db, self.pslr_threshold_db
                ),
            ));
        } else {
            results.push(
                ValidationResult::fail(
                    "PSLR",
                    pslr_db,
                    self.pslr_threshold_db,
                    1.0,
                    &format!(
                        "PSLR {:.2} dB exceeds {:.1} dB threshold",
                        pslr_db, self.pslr_threshold_db
                    ),
                )
                .with_recommendations(vec![
                    "Check focusing algorithm parameters".to_string(),
                    "Verify azimuth processing quality".to_string(),
                    "Inspect range compression settings".to_string(),
                ]),
            );
        }

        // ISLR (Integrated Sidelobe Ratio) check
        if islr_db <= self.islr_threshold_db {
            results.push(ValidationResult::pass(
                "ISLR",
                islr_db,
                self.islr_threshold_db,
                &format!(
                    "ISLR {:.2} dB meets requirement (<{:.1} dB)",
                    islr_db, self.islr_threshold_db
                ),
            ));
        } else {
            results.push(
                ValidationResult::fail(
                    "ISLR",
                    islr_db,
                    self.islr_threshold_db,
                    1.0,
                    &format!(
                        "ISLR {:.2} dB exceeds {:.1} dB threshold",
                        islr_db, self.islr_threshold_db
                    ),
                )
                .with_recommendations(vec![
                    "Review windowing function selection".to_string(),
                    "Check sidelobe suppression algorithms".to_string(),
                    "Verify processing gain settings".to_string(),
                ]),
            );
        }

        // IRW (Impulse Response Width) check
        if irw_meters <= self.irw_threshold_m {
            results.push(ValidationResult::pass(
                "IRW",
                irw_meters,
                self.irw_threshold_m,
                &format!(
                    "IRW {:.1}m meets requirement (<{:.1}m)",
                    irw_meters, self.irw_threshold_m
                ),
            ));
        } else {
            results.push(
                ValidationResult::fail(
                    "IRW",
                    irw_meters,
                    self.irw_threshold_m,
                    2.0,
                    &format!(
                        "IRW {:.1}m exceeds {:.1}m threshold",
                        irw_meters, self.irw_threshold_m
                    ),
                )
                .with_recommendations(vec![
                    "Check spatial resolution processing".to_string(),
                    "Verify bandwidth utilization".to_string(),
                    "Inspect multilook parameters".to_string(),
                ]),
            );
        }

        // Radiometric gain accuracy check
        let radiometric_error = radiometric_gain_error_db.abs();
        if radiometric_error <= self.radiometric_gain_tolerance {
            results.push(ValidationResult::pass(
                "Radiometric Gain",
                radiometric_error,
                self.radiometric_gain_tolerance,
                &format!(
                    "Radiometric error {:.2} dB within ±{:.1} dB tolerance",
                    radiometric_error, self.radiometric_gain_tolerance
                ),
            ));
        } else {
            results.push(
                ValidationResult::fail(
                    "Radiometric Gain",
                    radiometric_error,
                    self.radiometric_gain_tolerance,
                    self.radiometric_gain_tolerance,
                    &format!(
                        "Radiometric error {:.2} dB exceeds ±{:.1} dB tolerance",
                        radiometric_error, self.radiometric_gain_tolerance
                    ),
                )
                .with_recommendations(vec![
                    "Verify calibration constant accuracy".to_string(),
                    "Check antenna pattern correction".to_string(),
                    "Review radiometric calibration algorithm".to_string(),
                ]),
            );
        }

        info!(
            "🎯 PTA Results: PSLR={:.2}dB, ISLR={:.2}dB, IRW={:.1}m, Radio_err={:.2}dB",
            pslr_db, islr_db, irw_meters, radiometric_error
        );

        Ok(results)
    }

    /// 5. Geometry accuracy validation
    /// Validates known corner reflectors or coastline edges within <10m CE90 after terrain correction
    pub fn validate_geometry_accuracy(
        &self,
        measured_positions: &[(f64, f64)],  // (lat, lon) pairs
        reference_positions: &[(f64, f64)], // Known accurate positions
        coordinate_system: &str,
    ) -> SarResult<ValidationResult> {
        if !self.enable_geometry_check {
            return Ok(ValidationResult::pass(
                "Geometry Accuracy",
                0.0,
                0.0,
                "Validation disabled",
            ));
        }

        if measured_positions.len() != reference_positions.len() {
            return Err(SarError::InvalidParameter(
                "Measured and reference position arrays must have same length".to_string(),
            ));
        }

        if measured_positions.is_empty() {
            return Ok(ValidationResult::fail(
                "Geometry Accuracy",
                0.0,
                self.geometry_ce90_threshold,
                self.geometry_ce90_threshold,
                "No reference points available for geometry validation",
            ));
        }

        let mut azimuth_errors = Vec::new();
        let mut range_errors = Vec::new();
        let mut total_errors = Vec::new();

        for (measured, reference) in measured_positions.iter().zip(reference_positions.iter()) {
            // Convert lat/lon differences to meters (approximate)
            let lat_diff_m = (measured.0 - reference.0) * 111320.0; // ~111.32 km per degree latitude
                                                                    // Use mean latitude for longitudinal meter conversion (reduces bias)
            let mean_lat = 0.5 * (measured.0 + reference.0);
            let lon_diff_m = (measured.1 - reference.1) * 111320.0 * mean_lat.to_radians().cos();

            let total_error_m = (lat_diff_m.powi(2) + lon_diff_m.powi(2)).sqrt();

            azimuth_errors.push(lat_diff_m);
            range_errors.push(lon_diff_m);
            total_errors.push(total_error_m);
        }

        // Calculate CE90 (90% circular error)
        total_errors.sort_by(|a, b| a.partial_cmp(b).unwrap());
        let ce90_index = ((total_errors.len() as f64) * 0.9) as usize;
        let ce90_error = total_errors[ce90_index.min(total_errors.len() - 1)];

        // Calculate mean errors for diagnostic
        let mean_azimuth_error: f64 =
            azimuth_errors.iter().sum::<f64>() / azimuth_errors.len() as f64;
        let mean_range_error: f64 = range_errors.iter().sum::<f64>() / range_errors.len() as f64;

        info!(
            "📍 Geometry Check: CE90={:.2}m, mean_az_err={:.2}m, mean_rg_err={:.2}m ({} points)",
            ce90_error,
            mean_azimuth_error,
            mean_range_error,
            measured_positions.len()
        );

        if ce90_error <= self.geometry_ce90_threshold {
            Ok(ValidationResult::pass(
                "Geometry Accuracy",
                ce90_error,
                self.geometry_ce90_threshold,
                &format!(
                    "CE90 {:.2}m within {:.1}m threshold ({} reference points)",
                    ce90_error,
                    self.geometry_ce90_threshold,
                    measured_positions.len()
                ),
            ))
        } else {
            let mut recommendations = vec![
                "Check orbit accuracy and timing parameters".to_string(),
                "Verify DEM quality and resolution".to_string(),
                "Review terrain correction algorithm".to_string(),
            ];

            // Specific diagnostics based on error patterns
            if mean_azimuth_error.abs() > mean_range_error.abs() {
                recommendations
                    .push("Azimuth error dominant - check timing/PRF parameters".to_string());
            } else {
                recommendations
                    .push("Range error dominant - check range sampling/delay".to_string());
            }

            Ok(ValidationResult::fail(
                "Geometry Accuracy",
                ce90_error,
                self.geometry_ce90_threshold,
                self.geometry_ce90_threshold,
                &format!(
                    "CE90 {:.2}m exceeds {:.1}m threshold (az_err={:.2}m, rg_err={:.2}m)",
                    ce90_error, self.geometry_ce90_threshold, mean_azimuth_error, mean_range_error
                ),
            )
            .with_recommendations(recommendations))
        }
    }

    /// 6. LUT domain audit
    /// Validates all lookups are inside min/max line/pixel bounds
    pub fn validate_lut_domain_bounds(
        &self,
        line_indices: &[usize],
        pixel_indices: &[usize],
        max_lines: usize,
        max_pixels: usize,
        operation_name: &str,
    ) -> SarResult<ValidationResult> {
        if !self.enable_lut_domain_check {
            return Ok(ValidationResult::pass(
                "LUT Domain Check",
                0.0,
                0.0,
                "Validation disabled",
            ));
        }

        let mut out_of_bounds_count = 0;
        let mut max_line_violation = 0;
        let mut max_pixel_violation = 0;

        // Check line indices
        for &line_idx in line_indices {
            if line_idx >= max_lines {
                out_of_bounds_count += 1;
                max_line_violation = max_line_violation.max(line_idx);
            }
        }

        // Check pixel indices
        for &pixel_idx in pixel_indices {
            if pixel_idx >= max_pixels {
                out_of_bounds_count += 1;
                max_pixel_violation = max_pixel_violation.max(pixel_idx);
            }
        }

        let total_indices = line_indices.len() + pixel_indices.len();
        let violation_rate = out_of_bounds_count as f64 / total_indices as f64;

        if out_of_bounds_count == 0 {
            info!(
                "✅ LUT Domain ({}): All {} indices within bounds [0..{}, 0..{}]",
                operation_name, total_indices, max_lines, max_pixels
            );

            Ok(ValidationResult::pass(
                "LUT Domain Check",
                0.0,
                0.0,
                &format!(
                    "All {} indices within bounds for {}",
                    total_indices, operation_name
                ),
            ))
        } else {
            error!(
                "❌ LUT Domain ({}): {} out of {} indices out of bounds ({:.2}%)",
                operation_name,
                out_of_bounds_count,
                total_indices,
                violation_rate * 100.0
            );

            if max_line_violation > 0 {
                error!(
                    "   Max line violation: {} >= {}",
                    max_line_violation, max_lines
                );
            }
            if max_pixel_violation > 0 {
                error!(
                    "   Max pixel violation: {} >= {}",
                    max_pixel_violation, max_pixels
                );
            }

            Err(SarError::InvalidParameter(format!(
                "LUT domain violation in {}: {} indices out of bounds. Max violations: line={} (max={}), pixel={} (max={}). This indicates extrapolation which can cause processing artifacts.",
                operation_name, out_of_bounds_count, max_line_violation, max_lines, max_pixel_violation, max_pixels
            )))
        }
    }

    /// Run comprehensive validation suite
    pub fn run_validation_suite(
        &self,
        pre_deburst_power: Option<&Array2<f32>>,
        post_deburst_power: Option<&Array2<f32>>,
        overlap_regions: Option<&[(usize, usize, usize, usize)]>,
        debursted_image: Option<&Array2<f32>>,
        sigma0_image: Option<&Array2<f32>>,
        nesz_lut: Option<&[f32]>,
        pta_metrics: Option<(f64, f64, f64, f64)>, // (PSLR, ISLR, IRW, radiometric_error)
        geometry_points: Option<(&[(f64, f64)], &[(f64, f64)])>, // (measured, reference)
    ) -> SarResult<Vec<ValidationResult>> {
        let mut all_results = Vec::new();

        info!("🛡️  Starting SARdine Validation Gate Suite");

        // 1. Power preservation
        if let (Some(pre), Some(post), Some(regions)) =
            (pre_deburst_power, post_deburst_power, overlap_regions)
        {
            match self.validate_power_preservation(pre, post, regions) {
                Ok(result) => all_results.push(result),
                Err(e) => {
                    error!("Power preservation validation failed: {}", e);
                    all_results.push(ValidationResult::fail(
                        "Power Preservation",
                        0.0,
                        0.0,
                        0.0,
                        &format!("Validation error: {}", e),
                    ));
                }
            }
        }

        // 2. Uncovered pixels
        if let Some(image) = debursted_image {
            match self.validate_uncovered_pixels(image, None) {
                Ok(result) => all_results.push(result),
                Err(e) => {
                    error!("Uncovered pixel validation failed: {}", e);
                    all_results.push(ValidationResult::fail(
                        "Uncovered Pixels",
                        0.0,
                        0.0,
                        0.0,
                        &format!("Validation error: {}", e),
                    ));
                }
            }
        }

        // 3. NESZ compliance
        if let (Some(sigma0), Some(nesz)) = (sigma0_image, nesz_lut) {
            match self.validate_nesz_compliance(sigma0, 10.0, nesz, None) {
                Ok(result) => all_results.push(result),
                Err(e) => {
                    error!("NESZ validation failed: {}", e);
                    all_results.push(ValidationResult::fail(
                        "NESZ Compliance",
                        0.0,
                        0.0,
                        0.0,
                        &format!("Validation error: {}", e),
                    ));
                }
            }
        }

        // 4. Point Target Analysis
        if let Some((pslr, islr, irw, radio_err)) = pta_metrics {
            match self.validate_point_target_analysis(pslr, islr, irw, radio_err) {
                Ok(mut results) => all_results.append(&mut results),
                Err(e) => {
                    error!("PTA validation failed: {}", e);
                    all_results.push(ValidationResult::fail(
                        "Point Target Analysis",
                        0.0,
                        0.0,
                        0.0,
                        &format!("Validation error: {}", e),
                    ));
                }
            }
        }

        // 5. Geometry accuracy
        if let Some((measured, reference)) = geometry_points {
            match self.validate_geometry_accuracy(measured, reference, "WGS84") {
                Ok(result) => all_results.push(result),
                Err(e) => {
                    error!("Geometry validation failed: {}", e);
                    all_results.push(ValidationResult::fail(
                        "Geometry Accuracy",
                        0.0,
                        0.0,
                        0.0,
                        &format!("Validation error: {}", e),
                    ));
                }
            }
        }

        // Summary
        let passed_count = all_results.iter().filter(|r| r.passed).count();
        let total_count = all_results.len();

        if passed_count == total_count {
            info!("✅ Validation Suite: ALL {} checks PASSED", total_count);
        } else {
            warn!(
                "⚠️  Validation Suite: {}/{} checks passed ({} failed)",
                passed_count,
                total_count,
                total_count - passed_count
            );
        }

        Ok(all_results)
    }
}

/// Validation gate result summary
#[derive(Debug)]
pub struct ValidationSummary {
    pub total_tests: usize,
    pub passed_tests: usize,
    pub failed_tests: usize,
    pub critical_failures: Vec<String>,
    pub warnings: Vec<String>,
    pub overall_success: bool,
}

impl ValidationSummary {
    pub fn from_results(results: &[ValidationResult]) -> Self {
        let total_tests = results.len();
        let passed_tests = results.iter().filter(|r| r.passed).count();
        let failed_tests = total_tests - passed_tests;

        let critical_failures: Vec<String> = results
            .iter()
            .filter(|r| !r.passed && r.test_name.contains("LUT Domain"))
            .map(|r| r.diagnostic_message.clone())
            .collect();

        let warnings: Vec<String> = results
            .iter()
            .filter(|r| !r.passed && !r.test_name.contains("LUT Domain"))
            .map(|r| format!("{}: {}", r.test_name, r.diagnostic_message))
            .collect();

        Self {
            total_tests,
            passed_tests,
            failed_tests,
            critical_failures: critical_failures.clone(),
            warnings,
            overall_success: critical_failures.is_empty() && failed_tests == 0,
        }
    }

    pub fn log_summary(&self) {
        info!(
            "📊 Validation Summary: {}/{} tests passed",
            self.passed_tests, self.total_tests
        );

        if !self.critical_failures.is_empty() {
            for failure in &self.critical_failures {
                error!("💥 CRITICAL: {}", failure);
            }
        }

        if !self.warnings.is_empty() {
            for warning in &self.warnings {
                warn!("⚠️  WARNING: {}", warning);
            }
        }

        if self.overall_success {
            info!("🎉 All validation gates PASSED - processing quality verified!");
        } else {
            error!("❌ Validation gates FAILED - review diagnostics before proceeding");
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use ndarray::Array2;

    #[test]
    fn test_power_preservation_validation() {
        let gates = ValidationGates::default();

        // Create test data with ~1% power difference
        let pre_power = Array2::from_elem((100, 100), 1.0f32);
        let mut post_power = Array2::from_elem((100, 100), 1.01f32); // 1% increase

        let overlap_regions = vec![(10, 90, 10, 90)]; // Central region

        let result = gates
            .validate_power_preservation(&pre_power, &post_power, &overlap_regions)
            .unwrap();
        assert!(result.passed, "1% power difference should pass");

        // Test failure case
        post_power.fill(1.02f32); // 2% increase - should fail
        let result = gates
            .validate_power_preservation(&pre_power, &post_power, &overlap_regions)
            .unwrap();
        assert!(!result.passed, "2% power difference should fail");
    }

    #[test]
    fn test_uncovered_pixel_validation() {
        let gates = ValidationGates::default();

        // Test with mostly valid pixels
        let mut image = Array2::from_elem((1000, 1000), 1.0f32);
        image[(0, 0)] = 0.0; // 1 uncovered pixel out of 1M = 0.0001%

        let result = gates.validate_uncovered_pixels(&image, None).unwrap();
        assert!(result.passed, "Very low uncovered rate should pass");

        // Test with too many uncovered pixels
        for i in 0..3000 {
            // 3000 uncovered pixels = 0.3%
            let row = i / 1000;
            let col = i % 1000;
            image[(row, col)] = 0.0;
        }

        let result = gates.validate_uncovered_pixels(&image, None).unwrap();
        assert!(
            !result.passed,
            "0.3% uncovered rate should fail (>0.2% threshold)"
        );
    }

    #[test]
    fn test_lut_domain_validation() {
        let gates = ValidationGates::default();

        // Test valid indices
        let line_indices = vec![0, 10, 99];
        let pixel_indices = vec![0, 50, 199];

        let result = gates
            .validate_lut_domain_bounds(&line_indices, &pixel_indices, 100, 200, "test_operation")
            .unwrap();
        assert!(result.passed, "Valid indices should pass");

        // Test out-of-bounds indices
        let bad_line_indices = vec![0, 10, 100]; // 100 >= max_lines(100)
        let result = gates.validate_lut_domain_bounds(
            &bad_line_indices,
            &pixel_indices,
            100,
            200,
            "test_operation",
        );
        assert!(result.is_err(), "Out-of-bounds indices should return error");
    }

    #[test]
    fn test_validation_summary() {
        let results = vec![
            ValidationResult::pass("Test1", 0.5, 1.0, "Passed"),
            ValidationResult::fail("Test2", 1.5, 1.0, 0.1, "Failed"),
            ValidationResult::pass("Test3", 0.8, 1.0, "Passed"),
        ];

        let summary = ValidationSummary::from_results(&results);
        assert_eq!(summary.total_tests, 3);
        assert_eq!(summary.passed_tests, 2);
        assert_eq!(summary.failed_tests, 1);
        assert!(!summary.overall_success);
    }
}
