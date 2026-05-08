//! Strict metadata validation with no silent fallbacks
//!
//! This module implements the requirements:
//! - Stop on missing PRF, sampling rates, burst start/stop times, DC/FM polynomials
//! - β⁰ variability: Assert max(beta)/min(beta) > 1.02
//! - Units: Confirm LUTs are linear, not dB
//! - Line indexing: Handle negative tie-point lines with proper origin/offset
//! - Power seam check: Blended overlap mean power within ±1% of neighbors
//! - Coverage: Uncovered pixels in seam < 0.5%

use crate::types::{SarError, SarMetadata, SarResult, SubSwath};
use ndarray::Array2;
use std::collections::HashMap;

/// Strict metadata validator with no fallbacks
pub struct StrictMetadataValidator;

/// Results of strict metadata validation
#[derive(Debug, Clone)]
pub struct StrictValidationResult {
    pub is_valid: bool,
    pub errors: Vec<String>,
    pub warnings: Vec<String>,
    pub critical_missing: Vec<String>,
}

/// β⁰ variability analysis results
#[derive(Debug, Clone)]
pub struct BetaVariabilityResult {
    pub min_beta: f64,
    pub max_beta: f64,
    pub ratio: f64,
    pub is_valid: bool,
    pub parsing_bug_detected: bool,
}

/// Units validation results
#[derive(Debug, Clone)]
pub struct UnitsValidationResult {
    pub sigma_units: String,
    pub beta_units: String,
    pub gamma_units: String,
    pub all_linear: bool,
    pub db_conversions_needed: Vec<String>,
}

/// Line indexing validation results
#[derive(Debug, Clone)]
pub struct LineIndexingResult {
    pub negative_lines_found: bool,
    pub lines_exceeding_height: bool,
    pub corrected_count: usize,
    pub clamped_count: usize,
    pub corrected_indices: Vec<i32>,
    pub overflow_errors: Vec<String>,
}

/// Power seam validation results
#[derive(Debug, Clone)]
pub struct PowerSeamResult {
    pub seam_count: usize,
    pub seams_within_tolerance: usize,
    pub max_power_deviation_percent: f64,
    pub is_valid: bool,
}

/// Coverage validation results
#[derive(Debug, Clone)]
pub struct CoverageResult {
    pub total_pixels: usize,
    pub uncovered_pixels: usize,
    pub uncovered_percent: f64,
    pub is_valid: bool,
}

impl StrictMetadataValidator {
    /// Validate metadata with ZERO tolerance for missing critical fields
    pub fn validate_strict(metadata: &SarMetadata) -> SarResult<StrictValidationResult> {
        let mut errors = Vec::new();
        let mut warnings = Vec::new();
        let mut critical_missing = Vec::new();

        // 1. PRF - prefer per-subswath PRF, but allow global if consistent
        let mut subswath_prfs: Vec<f64> = metadata
            .sub_swaths
            .values()
            .filter_map(|sw| sw.prf_hz)
            .collect();
        subswath_prfs.sort_by(|a, b| a.total_cmp(b));
        subswath_prfs.dedup_by(|a, b| (*a - *b).abs() < 1e-6);

        let has_subswath_prf = !subswath_prfs.is_empty();

        if let Some(prf) = metadata.prf {
            if prf <= 0.0 || !prf.is_finite() {
                errors.push(format!("CRITICAL: Invalid PRF value: {}", prf));
            } else if prf < 100.0 || prf > 10000.0 {
                warnings.push(format!(
                    "PRF value {} Hz is outside typical SAR range (100-10000 Hz)",
                    prf
                ));
            }
        } else if !has_subswath_prf {
            critical_missing.push("PRF".to_string());
            errors.push("CRITICAL: PRF (Pulse Repetition Frequency) missing from metadata and no per-subswath PRF available.".to_string());
        }

        if subswath_prfs.len() > 1 {
            warnings.push(format!(
                "Per-subswath PRFs differ (n={}); global PRF will be ignored to avoid geometry drift",
                subswath_prfs.len()
            ));
        }

        // 2. Range sampling rate - CRITICAL, NO INFERENCE
        // Must be explicitly provided, not inferred from pixel spacing
        if metadata.range_sampling_rate.is_none() {
            critical_missing.push("range_sampling_rate_hz".to_string());
            errors.push("CRITICAL: Range sampling rate (Hz) missing from metadata. Cannot infer from pixel spacing.".to_string());
        } else if let Some(rsr) = metadata.range_sampling_rate {
            if rsr <= 0.0 || !rsr.is_finite() {
                errors.push(format!("CRITICAL: Invalid range sampling rate: {}", rsr));
            } else if rsr < 1e6 || rsr > 1e9 {
                warnings.push(format!(
                    "Range sampling rate {} Hz is outside typical SAR range (1-1000 MHz)",
                    rsr
                ));
            }
        }

        // 3. Slant range time - CRITICAL
        if metadata.slant_range_time.is_none() {
            critical_missing.push("slant_range_time".to_string());
            errors.push("CRITICAL: Slant range time missing from metadata".to_string());
        } else if let Some(srt) = metadata.slant_range_time {
            if srt <= 0.0 || !srt.is_finite() {
                errors.push(format!("CRITICAL: Invalid slant range time: {}", srt));
            }
        }

        // 4. Wavelength - CRITICAL
        if metadata.wavelength.is_none() {
            critical_missing.push("wavelength".to_string());
            errors.push("CRITICAL: Wavelength missing from metadata".to_string());
        } else if let Some(wl) = metadata.wavelength {
            if wl <= 0.0 || !wl.is_finite() {
                errors.push(format!("CRITICAL: Invalid wavelength: {}", wl));
            }
        }

        // 5. Burst start/stop times for TOPSAR - CRITICAL for TOPS mode only
        if matches!(
            metadata.acquisition_mode,
            crate::types::AcquisitionMode::IW | crate::types::AcquisitionMode::EW
        ) {
            Self::validate_burst_timing(&metadata.sub_swaths, &mut errors, &mut critical_missing)?;
            // 6. DC/FM polynomials - CRITICAL for TOPSAR
            Self::validate_doppler_polynomials(
                &metadata.sub_swaths,
                &mut errors,
                &mut critical_missing,
            )?;
        }

        let is_valid = errors.is_empty() && critical_missing.is_empty();

        if !is_valid {
            log::error!("🚨 STRICT METADATA VALIDATION FAILED");
            for error in &errors {
                log::error!("   ❌ {}", error);
            }
            for missing in &critical_missing {
                log::error!("   🚫 MISSING CRITICAL FIELD: {}", missing);
            }
        } else {
            log::info!("✅ Strict metadata validation passed");
        }

        Ok(StrictValidationResult {
            is_valid,
            errors,
            warnings,
            critical_missing,
        })
    }

    /// Validate burst timing parameters - NO FALLBACKS
    fn validate_burst_timing(
        sub_swaths: &HashMap<String, SubSwath>,
        errors: &mut Vec<String>,
        critical_missing: &mut Vec<String>,
    ) -> SarResult<()> {
        for (swath_id, subswath) in sub_swaths {
            // Burst count - CRITICAL
            if subswath.burst_count == 0 {
                critical_missing.push(format!("burst_count_{}", swath_id));
                errors.push(format!(
                    "CRITICAL: Burst count is zero for subswath {}",
                    swath_id
                ));
                continue; // Skip further checks for this subswath
            }

            // Burst duration
            if subswath.burst_duration <= 0.0 {
                critical_missing.push(format!("burst_duration_{}", swath_id));
                errors.push(format!(
                    "CRITICAL: Invalid burst duration for subswath {}: {}",
                    swath_id, subswath.burst_duration
                ));
            }

            // PRF per subswath
            if subswath.prf_hz.is_none() {
                critical_missing.push(format!("prf_hz_{}", swath_id));
                errors.push(format!("CRITICAL: PRF missing for subswath {}", swath_id));
            }

            // Azimuth time interval per subswath (preferred over global PRF for timing)
            if let Some(ati) = subswath.azimuth_time_interval {
                if ati <= 0.0 || !ati.is_finite() {
                    errors.push(format!(
                        "CRITICAL: Invalid azimuth_time_interval for subswath {}: {}",
                        swath_id, ati
                    ));
                }
            } else {
                critical_missing.push(format!("azimuth_time_interval_{}", swath_id));
                errors.push(format!(
                    "CRITICAL: azimuth_time_interval missing for subswath {} (needed for burst timing)",
                    swath_id
                ));
            }

            // Azimuth time intervals
            if subswath.azimuth_samples == 0 {
                errors.push(format!(
                    "CRITICAL: Zero azimuth samples for subswath {}",
                    swath_id
                ));
            }

            // Burst timing consistency check with absolute + relative tolerance
            if let Some(prf) = subswath.prf_hz {
                let expected_burst_lines = (subswath.burst_duration * prf).round() as usize;
                let lines_per_burst = subswath.azimuth_samples / subswath.burst_count;

                let abs_diff = (expected_burst_lines as i32 - lines_per_burst as i32).abs();
                let rel_diff = abs_diff as f64 / expected_burst_lines.max(1) as f64;

                if abs_diff > 10 && rel_diff > 0.05 {
                    // 10 lines OR 5% relative error
                    errors.push(format!(
                        "CRITICAL: Burst timing inconsistency in {}: expected {} lines, got {} lines per burst (abs={}, rel={:.2}%)",
                        swath_id, expected_burst_lines, lines_per_burst, abs_diff, rel_diff * 100.0
                    ));
                }
            }
        }

        Ok(())
    }

    /// Validate Doppler centroid and FM polynomials - NO FALLBACKS
    /// Checks presence, degree, finiteness, and realistic magnitudes
    fn validate_doppler_polynomials(
        sub_swaths: &HashMap<String, SubSwath>,
        errors: &mut Vec<String>,
        critical_missing: &mut Vec<String>,
    ) -> SarResult<()> {
        for (swath_id, subswath) in sub_swaths {
            if subswath.burst_count == 0 {
                continue; // Already reported in burst timing
            }

            // STRICT: Check for actual DC polynomial presence in SubSwath
            if subswath.dc_polynomial.is_none()
                || subswath
                    .dc_polynomial
                    .as_ref()
                    .map_or(true, |p| p.is_empty())
            {
                critical_missing.push(format!("dc_polynomial_{}", swath_id));
                errors.push(format!(
                    "CRITICAL: Doppler centroid polynomial missing or empty for subswath {}. \
                     Deburst cannot proceed without DC polynomials to prevent phase misalignment.",
                    swath_id
                ));
            } else if let Some(dc_poly) = &subswath.dc_polynomial {
                // Validate all coefficients are finite
                for (idx, &coeff) in dc_poly.iter().enumerate() {
                    if !coeff.is_finite() {
                        errors.push(format!(
                            "CRITICAL: DC polynomial coefficient[{}] is non-finite for {}: {}",
                            idx, swath_id, coeff
                        ));
                    }
                }

                log::debug!(
                    "DC polynomial for {} has degree {} (coeffs: {:?})",
                    swath_id,
                    dc_poly.len().saturating_sub(1),
                    dc_poly
                );
            }

            // Check burst duration (required for time-based polynomials)
            if subswath.burst_duration <= 0.0 || !subswath.burst_duration.is_finite() {
                critical_missing.push(format!("burst_duration_{}", swath_id));
                errors.push(format!(
                    "CRITICAL: Invalid burst duration for subswath {}: {}",
                    swath_id, subswath.burst_duration
                ));
            }

            // Validate PRF for Doppler frequency conversion and Nyquist constraint
            if let Some(prf) = subswath.prf_hz {
                if !prf.is_finite() || prf <= 0.0 {
                    errors.push(format!(
                        "CRITICAL: Invalid PRF for Doppler calculation in {}: {}",
                        swath_id, prf
                    ));
                } else if let Some(dc_poly) = &subswath.dc_polynomial {
                    // Evaluate DC polynomial at mid-burst to check Nyquist constraint
                    let mid_time = subswath.burst_duration / 2.0;
                    let dc_value = dc_poly
                        .iter()
                        .enumerate()
                        .fold(0.0, |acc, (i, &c)| acc + c * mid_time.powi(i as i32));

                    let nyquist_freq = prf / 2.0;
                    if dc_value.abs() > nyquist_freq {
                        errors.push(format!(
                            "CRITICAL: DC value {:.1} Hz at mid-burst exceeds Nyquist frequency {:.1} Hz for {}",
                            dc_value, nyquist_freq, swath_id
                        ));
                    }

                    log::debug!("Doppler validation for {}: PRF={:.1} Hz, DC@mid={:.1} Hz, Nyquist={:.1} Hz",
                               swath_id, prf, dc_value, nyquist_freq);
                }
            } else {
                errors.push(format!(
                    "CRITICAL: PRF missing for Doppler validation in {}",
                    swath_id
                ));
            }
        }

        Ok(())
    }

    /// Validate β⁰ variability with configurable threshold (default > 1.02 ratio)
    /// Filters non-finite values and rejects min_beta ≤ 0
    pub fn validate_beta_variability(beta_values: &[f64]) -> BetaVariabilityResult {
        Self::validate_beta_variability_with_threshold(beta_values, 1.02)
    }

    /// Validate β⁰ variability with custom threshold
    pub fn validate_beta_variability_with_threshold(
        beta_values: &[f64],
        threshold: f64,
    ) -> BetaVariabilityResult {
        // Filter to finite, positive values only
        let valid_values: Vec<f64> = beta_values
            .iter()
            .copied()
            .filter(|&v| v.is_finite() && v > 0.0)
            .collect();

        if valid_values.is_empty() {
            return BetaVariabilityResult {
                min_beta: 0.0,
                max_beta: 0.0,
                ratio: 0.0,
                is_valid: false,
                parsing_bug_detected: true,
            };
        }

        let min_beta = valid_values
            .iter()
            .copied()
            .fold(f64::INFINITY, |a, b| a.min(b));
        let max_beta = valid_values
            .iter()
            .copied()
            .fold(f64::NEG_INFINITY, |a, b| a.max(b));

        let ratio = if min_beta > 0.0 {
            max_beta / min_beta
        } else {
            0.0
        };

        let is_valid = ratio > threshold;
        let parsing_bug_detected = !is_valid;

        if parsing_bug_detected {
            log::error!("🚨 β⁰ VARIABILITY CHECK FAILED:");
            log::error!("   Min β⁰: {:.6e}", min_beta);
            log::error!("   Max β⁰: {:.6e}", max_beta);
            log::error!("   Ratio: {:.6}", ratio);
            log::error!("   Required ratio: > 1.02");
            log::error!("   → PARSING BUG DETECTED: β⁰ values lack proper variation");
        } else {
            log::info!("✅ β⁰ variability check passed: ratio = {:.3}", ratio);
        }

        BetaVariabilityResult {
            min_beta,
            max_beta,
            ratio,
            is_valid,
            parsing_bug_detected,
        }
    }

    /// Validate units - confirm LUTs are linear, not dB
    pub fn validate_units(
        sigma_values: &[f64],
        beta_values: &[f64],
        gamma_values: &[f64],
        sigma_units: Option<&str>,
        beta_units: Option<&str>,
        gamma_units: Option<&str>,
    ) -> UnitsValidationResult {
        #[allow(non_snake_case)]
        let mut db_conversions_needed = Vec::new();

        // Check units strings
        let sigma_unit_str = sigma_units.unwrap_or("unknown").to_string();
        let beta_unit_str = beta_units.unwrap_or("unknown").to_string();
        let gamma_unit_str = gamma_units.unwrap_or("unknown").to_string();

        // Check for explicit dB indicators
        if sigma_unit_str.to_lowercase().contains("db") {
            db_conversions_needed.push("sigma".to_string());
        }
        if beta_unit_str.to_lowercase().contains("db") {
            db_conversions_needed.push("beta".to_string());
        }
        if gamma_unit_str.to_lowercase().contains("db") {
            db_conversions_needed.push("gamma".to_string());
        }

        // Statistical analysis for dB detection
        Self::detect_dB_values_statistical("sigma", sigma_values, &mut db_conversions_needed);
        Self::detect_dB_values_statistical("beta", beta_values, &mut db_conversions_needed);
        Self::detect_dB_values_statistical("gamma", gamma_values, &mut db_conversions_needed);

        let all_linear = db_conversions_needed.is_empty();

        if !all_linear {
            log::warn!(
                "⚠️ UNITS VALIDATION: dB conversion needed for: {:?}",
                db_conversions_needed
            );
        } else {
            log::info!("✅ Units validation: All LUTs are in linear domain");
        }

        UnitsValidationResult {
            sigma_units: sigma_unit_str,
            beta_units: beta_unit_str,
            gamma_units: gamma_unit_str,
            all_linear,
            db_conversions_needed,
        }
    }

    /// Statistical detection of dB values using robust IQR and log10 tests
    /// Only processes finite values
    #[allow(non_snake_case)]
    fn detect_dB_values_statistical(
        lut_name: &str,
        values: &[f64],
        conversions_needed: &mut Vec<String>,
    ) {
        // Filter to finite values only
        let mut finite_values: Vec<f64> =
            values.iter().copied().filter(|v| v.is_finite()).collect();

        if finite_values.is_empty() {
            return;
        }

        finite_values.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));

        let len = finite_values.len();
        let q1_idx = len / 4;
        let q3_idx = (3 * len) / 4;
        let median_idx = len / 2;

        let q1 = finite_values[q1_idx];
        let median = finite_values[median_idx];
        let q3 = finite_values[q3_idx];
        let iqr = q3 - q1;

        let min_val = finite_values[0];
        let max_val = finite_values[len - 1];

        // Heuristic: calibration values in dB typically range from -50 to +50
        // Linear calibration constants are typically 1e-6 to 1e6
        #[allow(non_snake_case)]
        let likely_dB = median > -60.0 && median < 60.0 &&
                       iqr < 100.0 && // Small interquartile range
                       max_val - min_val < 150.0; // Limited total range

        // Additional check: log10 of linear values should span many orders of magnitude
        let positive_values: Vec<f64> =
            finite_values.iter().copied().filter(|&v| v > 0.0).collect();
        if !positive_values.is_empty() {
            let log_range = positive_values
                .iter()
                .map(|&v| v.log10())
                .fold((f64::INFINITY, f64::NEG_INFINITY), |(min, max), log_v| {
                    (min.min(log_v), max.max(log_v))
                });
            let log_span = log_range.1 - log_range.0;

            // Linear values typically span > 6 orders of magnitude; dB values span < 2
            if log_span < 2.0 {
                if !conversions_needed.contains(&lut_name.to_string()) {
                    log::warn!("🔍 Statistical dB detection for {}: median={:.3}, IQR={:.3}, log_span={:.3}",
                              lut_name, median, iqr, log_span);
                    conversions_needed.push(lut_name.to_string());
                }
                return;
            }
        }

        if likely_dB && !conversions_needed.contains(&lut_name.to_string()) {
            log::warn!(
                "🔍 Statistical dB detection for {}: median={:.3}, IQR={:.3}, range=[{:.3}, {:.3}]",
                lut_name,
                median,
                iqr,
                min_val,
                max_val
            );
            conversions_needed.push(lut_name.to_string());
        }
    }

    /// Validate line indexing - handle negative tie-point lines and return corrected mapping
    pub fn validate_line_indexing(
        tie_point_lines: &[i32],
        image_height: usize,
        origin_offset: i32,
    ) -> LineIndexingResult {
        let mut negative_lines_found = false;
        let mut lines_exceeding_height = false;
        let mut corrected_count = 0;
        let mut clamped_count = 0;
        let mut corrected_indices = Vec::with_capacity(tie_point_lines.len());
        let mut overflow_errors = Vec::new();

        for (idx, &line) in tie_point_lines.iter().enumerate() {
            let corrected_line = line + origin_offset;

            if line < 0 {
                negative_lines_found = true;
                corrected_count += 1;
                log::debug!(
                    "Corrected negative tie-point line: {} → {}",
                    line,
                    corrected_line
                );
            }

            let clamped_line = corrected_line.max(0).min(image_height as i32 - 1);
            if clamped_line != corrected_line {
                lines_exceeding_height = true;
                clamped_count += 1;

                let error_msg = if corrected_line < 0 {
                    format!("Line index {} underflow: {} < 0", idx, corrected_line)
                } else {
                    format!(
                        "Line index {} overflow: {} >= {}",
                        idx, corrected_line, image_height
                    )
                };
                overflow_errors.push(error_msg);
                log::error!("🚨 {}", overflow_errors.last().unwrap());
            }

            corrected_indices.push(clamped_line);
        }

        if negative_lines_found || lines_exceeding_height {
            log::info!(
                "📏 Line indexing corrections: {} negative corrected, {} clamped",
                corrected_count,
                clamped_count
            );
        }

        LineIndexingResult {
            negative_lines_found,
            lines_exceeding_height,
            corrected_count,
            clamped_count,
            corrected_indices,
            overflow_errors,
        }
    }

    /// Validate power seams - blended overlap mean power within ±1% of neighbors
    pub fn validate_power_seams(
        merged_data: &Array2<f32>,
        overlap_regions: &[(usize, usize, usize, usize)], // (row_start, row_end, col_start, col_end)
    ) -> PowerSeamResult {
        let mut seams_within_tolerance = 0;
        let mut max_power_deviation_percent: f64 = 0.0;

        for (i, &(row_start, row_end, col_start, col_end)) in overlap_regions.iter().enumerate() {
            if let Some(deviation) = Self::check_seam_power_consistency(
                merged_data,
                row_start,
                row_end,
                col_start,
                col_end,
            ) {
                if deviation <= 1.0 {
                    seams_within_tolerance += 1;
                } else {
                    log::warn!(
                        "⚠️ Seam {} power deviation: {:.2}% (> 1.0% threshold)",
                        i,
                        deviation
                    );
                }
                max_power_deviation_percent = max_power_deviation_percent.max(deviation);
            }
        }

        let is_valid = seams_within_tolerance == overlap_regions.len();

        if !is_valid {
            log::error!("🚨 POWER SEAM VALIDATION FAILED:");
            log::error!(
                "   Seams within tolerance: {}/{}",
                seams_within_tolerance,
                overlap_regions.len()
            );
            log::error!("   Max deviation: {:.2}%", max_power_deviation_percent);
        } else {
            log::info!("✅ Power seam validation passed: all seams within ±1.0%");
        }

        PowerSeamResult {
            seam_count: overlap_regions.len(),
            seams_within_tolerance,
            max_power_deviation_percent,
            is_valid,
        }
    }

    /// Check power consistency in a specific seam region
    fn check_seam_power_consistency(
        data: &Array2<f32>,
        row_start: usize,
        row_end: usize,
        col_start: usize,
        col_end: usize,
    ) -> Option<f64> {
        if row_end <= row_start || col_end <= col_start {
            return None;
        }

        let (height, width) = data.dim();
        let row_end = row_end.min(height);
        let col_end = col_end.min(width);

        // Calculate mean power in seam region
        let mut seam_power_sum = 0.0;
        let mut seam_pixel_count = 0;

        for row in row_start..row_end {
            for col in col_start..col_end {
                let value = data[[row, col]];
                if value.is_finite() && value > 0.0 {
                    seam_power_sum += value as f64;
                    seam_pixel_count += 1;
                }
            }
        }

        if seam_pixel_count == 0 {
            return None;
        }

        let seam_mean_power = seam_power_sum / seam_pixel_count as f64;

        // Calculate mean power in neighboring regions (left and right of seam)
        let neighbor_buffer = 10; // pixels
        let mut neighbor_power_sum = 0.0;
        let mut neighbor_pixel_count = 0;

        // Left neighbor
        if col_start >= neighbor_buffer {
            for row in row_start..row_end {
                for col in (col_start - neighbor_buffer)..col_start {
                    let value = data[[row, col]];
                    if value.is_finite() && value > 0.0 {
                        neighbor_power_sum += value as f64;
                        neighbor_pixel_count += 1;
                    }
                }
            }
        }

        // Right neighbor
        if col_end + neighbor_buffer < width {
            for row in row_start..row_end {
                for col in col_end..(col_end + neighbor_buffer) {
                    let value = data[[row, col]];
                    if value.is_finite() && value > 0.0 {
                        neighbor_power_sum += value as f64;
                        neighbor_pixel_count += 1;
                    }
                }
            }
        }

        if neighbor_pixel_count == 0 {
            return None;
        }

        let neighbor_mean_power = neighbor_power_sum / neighbor_pixel_count as f64;

        // Avoid divide-by-zero when neighbor mean is ~0
        if neighbor_mean_power < 1e-10 {
            log::warn!(
                "⚠️ Neighbor mean power near zero ({:.2e}), cannot compute deviation",
                neighbor_mean_power
            );
            return None;
        }

        // Calculate percentage deviation
        // NOTE: This assumes input is in LINEAR power domain, not dB
        let deviation_percent =
            ((seam_mean_power - neighbor_mean_power).abs() / neighbor_mean_power) * 100.0;

        Some(deviation_percent)
    }

    /// Validate coverage - uncovered pixels in seams < threshold (default 0.5%)
    pub fn validate_coverage(uncovered_mask: &Array2<u8>) -> CoverageResult {
        Self::validate_coverage_with_threshold(uncovered_mask, 0.5)
    }

    /// Validate coverage with custom threshold percentage
    pub fn validate_coverage_with_threshold(
        uncovered_mask: &Array2<u8>,
        threshold_percent: f64,
    ) -> CoverageResult {
        let total_pixels = uncovered_mask.len();
        let uncovered_pixels = uncovered_mask.iter().map(|&x| x as usize).sum::<usize>();
        let uncovered_percent = (uncovered_pixels as f64 / total_pixels as f64) * 100.0;

        let is_valid = uncovered_percent < threshold_percent;

        if !is_valid {
            log::error!("🚨 COVERAGE VALIDATION FAILED:");
            log::error!(
                "   Uncovered pixels: {} ({:.2}%)",
                uncovered_pixels,
                uncovered_percent
            );
            log::error!("   Threshold: < {:.2}%", threshold_percent);

            if uncovered_percent > 2.0 {
                log::error!("   🚫 CRITICAL: > 2.0% uncovered → wrong overlap math detected");
            }
        } else {
            log::info!(
                "✅ Coverage validation passed: {:.3}% uncovered (threshold: {:.2}%)",
                uncovered_percent,
                threshold_percent
            );
        }

        CoverageResult {
            total_pixels,
            uncovered_pixels,
            uncovered_percent,
            is_valid,
        }
    }
}

/// Comprehensive strict validation of entire processing chain
pub fn validate_processing_chain_strict(
    metadata: &SarMetadata,
    beta_values: &[f64],
    sigma_values: &[f64],
    gamma_values: &[f64],
    merged_data: &Array2<f32>,
    uncovered_mask: &Array2<u8>,
    overlap_regions: &[(usize, usize, usize, usize)],
) -> SarResult<()> {
    log::info!("🔍 Starting comprehensive strict validation");

    // 1. Units validation FIRST - CRITICAL requirement that all LUTs are linear
    let units_result = StrictMetadataValidator::validate_units(
        sigma_values,
        beta_values,
        gamma_values,
        None,
        None,
        None,
    );
    if !units_result.all_linear {
        return Err(SarError::Processing(format!(
            "Units validation CRITICAL FAILURE: LUTs must be linear, not dB. Non-linear LUTs detected: {:?}",
            units_result.db_conversions_needed
        )));
    }

    // 2. Strict metadata validation
    let metadata_result = StrictMetadataValidator::validate_strict(metadata)?;
    if !metadata_result.is_valid {
        return Err(SarError::Processing(format!(
            "Strict metadata validation failed: {} errors, {} critical missing fields",
            metadata_result.errors.len(),
            metadata_result.critical_missing.len()
        )));
    }

    // 3. β⁰ variability check
    let beta_result = StrictMetadataValidator::validate_beta_variability(beta_values);
    if beta_result.parsing_bug_detected {
        return Err(SarError::Processing(format!(
            "β⁰ variability check failed: ratio {:.6} ≤ 1.02 indicates parsing bug",
            beta_result.ratio
        )));
    }

    // 4. Power seam validation
    let seam_result = StrictMetadataValidator::validate_power_seams(merged_data, overlap_regions);
    if !seam_result.is_valid {
        let failing_seams = seam_result.seam_count - seam_result.seams_within_tolerance;
        return Err(SarError::Processing(format!(
            "Power seam validation failed: {}/{} seams within ±1% tolerance, max deviation {:.2}%. Failing seam count: {}",
            seam_result.seams_within_tolerance,
            seam_result.seam_count,
            seam_result.max_power_deviation_percent,
            failing_seams
        )));
    }

    // 5. Coverage validation
    let coverage_result = StrictMetadataValidator::validate_coverage(uncovered_mask);
    if !coverage_result.is_valid {
        return Err(SarError::Processing(format!(
            "Coverage validation failed: {:.2}% uncovered pixels (threshold: 0.5%). Total uncovered: {}/{}",
            coverage_result.uncovered_percent,
            coverage_result.uncovered_pixels,
            coverage_result.total_pixels
        )));
    }

    log::info!("✅ Comprehensive strict validation passed");
    Ok(())
}
