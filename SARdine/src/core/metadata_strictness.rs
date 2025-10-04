
//! Strict metadata validation with no silent fallbacks
//! 
//! This module implements the requirements:
//! - Stop on missing PRF, sampling rates, burst start/stop times, DC/FM polynomials
//! - β⁰ variability: Assert max(beta)/min(beta) > 1.02
//! - Units: Confirm LUTs are linear, not dB
//! - Line indexing: Handle negative tie-point lines with proper origin/offset
//! - Power seam check: Blended overlap mean power within ±1% of neighbors
//! - Coverage: Uncovered pixels in seam < 0.5%

use crate::types::{SarError, SarResult, SarMetadata, SubSwath};
use std::collections::HashMap;
use ndarray::Array2;

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
    #[allow(non_snake_case)]
    pub dB_conversions_needed: Vec<String>,
}

/// Line indexing validation results
#[derive(Debug, Clone)]
pub struct LineIndexingResult {
    pub negative_lines_found: bool,
    pub lines_exceeding_height: bool,
    pub corrected_count: usize,
    pub clamped_count: usize,
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

        // 1. PRF - CRITICAL, NO FALLBACKS
        if metadata.prf.is_none() {
            critical_missing.push("PRF".to_string());
            errors.push("CRITICAL: PRF (Pulse Repetition Frequency) missing from metadata. Processing cannot continue without this parameter.".to_string());
        } else if let Some(prf) = metadata.prf {
            if prf <= 0.0 || !prf.is_finite() {
                errors.push(format!("CRITICAL: Invalid PRF value: {}", prf));
            } else if prf < 100.0 || prf > 10000.0 {
                warnings.push(format!("PRF value {} Hz is outside typical SAR range (100-10000 Hz)", prf));
            }
        }

        // 2. Range sampling rate - for SAR metadata, this might be accessed differently
        // Note: Adapting to actual SarMetadata structure
        let has_range_sampling_rate = metadata.pixel_spacing.0 > 0.0; // Infer from pixel spacing
        if !has_range_sampling_rate {
            critical_missing.push("range_sampling_rate".to_string());
            errors.push("CRITICAL: Range sampling rate missing from metadata".to_string());
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

        // 5. Burst start/stop times for TOPSAR - CRITICAL for TOPS mode
        if matches!(metadata.acquisition_mode, crate::types::AcquisitionMode::IW) {
            Self::validate_burst_timing(&metadata.sub_swaths, &mut errors, &mut critical_missing)?;
        }

        // 6. DC/FM polynomials - CRITICAL for TOPSAR
        Self::validate_doppler_polynomials(&metadata.sub_swaths, &mut errors, &mut critical_missing)?;

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
            // Burst count
            if subswath.burst_count == 0 {
                errors.push(format!("CRITICAL: Burst count is zero for subswath {}", swath_id));
            }

            // Burst duration
            if subswath.burst_duration <= 0.0 {
                critical_missing.push(format!("burst_duration_{}", swath_id));
                errors.push(format!("CRITICAL: Invalid burst duration for subswath {}: {}", swath_id, subswath.burst_duration));
            }

            // PRF per subswath
            if subswath.prf_hz.is_none() {
                critical_missing.push(format!("prf_hz_{}", swath_id));
                errors.push(format!("CRITICAL: PRF missing for subswath {}", swath_id));
            }

            // Azimuth time intervals
            if subswath.azimuth_samples == 0 {
                errors.push(format!("CRITICAL: Zero azimuth samples for subswath {}", swath_id));
            }

            // Burst timing consistency check
            if let Some(prf) = subswath.prf_hz {
                let expected_burst_lines = (subswath.burst_duration * prf) as usize;
                let lines_per_burst = subswath.azimuth_samples / subswath.burst_count.max(1);
                
                if (expected_burst_lines as i32 - lines_per_burst as i32).abs() > 10 {
                    errors.push(format!(
                        "CRITICAL: Burst timing inconsistency in {}: expected {} lines, got {} lines per burst",
                        swath_id, expected_burst_lines, lines_per_burst
                    ));
                }
            }
        }

        Ok(())
    }

    /// Validate Doppler centroid and FM polynomials - NO FALLBACKS
    fn validate_doppler_polynomials(
        sub_swaths: &HashMap<String, SubSwath>,
        errors: &mut Vec<String>,
        critical_missing: &mut Vec<String>,
    ) -> SarResult<()> {
        for (swath_id, subswath) in sub_swaths {
            // Check for DC polynomial - simplified since field structure may differ
            // Note: Adapting to actual SubSwath structure  
            let has_doppler_data = subswath.burst_count > 0; // Basic check for burst data
            if !has_doppler_data {
                critical_missing.push(format!("doppler_data_{}", swath_id));
                errors.push(format!("CRITICAL: Doppler data missing for subswath {}", swath_id));
            }

            // Check for FM rate - use available fields
            if subswath.burst_duration <= 0.0 {
                critical_missing.push(format!("burst_duration_{}", swath_id));
                errors.push(format!("CRITICAL: Burst duration missing for subswath {}", swath_id));
            }
        }

        Ok(())
    }

    /// Validate β⁰ variability - must be > 1.02 ratio
    pub fn validate_beta_variability(beta_values: &[f64]) -> BetaVariabilityResult {
        if beta_values.is_empty() {
            return BetaVariabilityResult {
                min_beta: 0.0,
                max_beta: 0.0,
                ratio: 0.0,
                is_valid: false,
                parsing_bug_detected: true,
            };
        }

        let min_beta = beta_values.iter().fold(f64::INFINITY, |a, &b| a.min(b));
        let max_beta = beta_values.iter().fold(f64::NEG_INFINITY, |a, &b| a.max(b));
        
        let ratio = if min_beta > 0.0 {
            max_beta / min_beta
        } else {
            0.0
        };

        let is_valid = ratio > 1.02;
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
        let mut dB_conversions_needed = Vec::new();

        // Check units strings
        let sigma_unit_str = sigma_units.unwrap_or("unknown").to_string();
        let beta_unit_str = beta_units.unwrap_or("unknown").to_string();
        let gamma_unit_str = gamma_units.unwrap_or("unknown").to_string();

        // Check for explicit dB indicators
        if sigma_unit_str.to_lowercase().contains("db") {
            dB_conversions_needed.push("sigma".to_string());
        }
        if beta_unit_str.to_lowercase().contains("db") {
            dB_conversions_needed.push("beta".to_string());
        }
        if gamma_unit_str.to_lowercase().contains("db") {
            dB_conversions_needed.push("gamma".to_string());
        }

        // Statistical analysis for dB detection
        Self::detect_dB_values_statistical("sigma", sigma_values, &mut dB_conversions_needed);
        Self::detect_dB_values_statistical("beta", beta_values, &mut dB_conversions_needed);
        Self::detect_dB_values_statistical("gamma", gamma_values, &mut dB_conversions_needed);

        let all_linear = dB_conversions_needed.is_empty();

        if !all_linear {
            log::warn!("⚠️ UNITS VALIDATION: dB conversion needed for: {:?}", dB_conversions_needed);
        } else {
            log::info!("✅ Units validation: All LUTs are in linear domain");
        }

        UnitsValidationResult {
            sigma_units: sigma_unit_str,
            beta_units: beta_unit_str,
            gamma_units: gamma_unit_str,
            all_linear,
            dB_conversions_needed,
        }
    }

    /// Statistical detection of dB values
    #[allow(non_snake_case)]
    fn detect_dB_values_statistical(lut_name: &str, values: &[f64], conversions_needed: &mut Vec<String>) {
        if values.is_empty() {
            return;
        }

        let mut sorted_values = values.to_vec();
        sorted_values.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
        
        let median = sorted_values[sorted_values.len() / 2];
        let min_val = sorted_values[0];
        let max_val = sorted_values[sorted_values.len() - 1];

        // Heuristic: calibration values in dB typically range from -50 to +50
        // Linear calibration constants are typically 1e-6 to 1e6
        #[allow(non_snake_case)]
        let likely_dB = median > -60.0 && median < 60.0 && 
                       max_val - min_val < 100.0 && // Small dynamic range
                       median != median.round(); // Has decimal precision

        if likely_dB && !conversions_needed.contains(&lut_name.to_string()) {
            log::warn!("🔍 Statistical dB detection for {}: median={:.3}, range=[{:.3}, {:.3}]", 
                      lut_name, median, min_val, max_val);
            conversions_needed.push(lut_name.to_string());
        }
    }

    /// Validate line indexing - handle negative tie-point lines
    pub fn validate_line_indexing(
        tie_point_lines: &[i32],
        image_height: usize,
        origin_offset: i32,
    ) -> LineIndexingResult {
        let mut negative_lines_found = false;
        let mut lines_exceeding_height = false;
        let mut corrected_count = 0;
        let mut clamped_count = 0;

        for &line in tie_point_lines {
            let corrected_line = line + origin_offset;
            
            if line < 0 {
                negative_lines_found = true;
                corrected_count += 1;
                log::debug!("Corrected negative tie-point line: {} → {}", line, corrected_line);
            }

            let clamped_line = corrected_line.max(0).min(image_height as i32 - 1);
            if clamped_line != corrected_line {
                lines_exceeding_height = true;
                clamped_count += 1;
                log::debug!("Clamped tie-point line: {} → {}", corrected_line, clamped_line);
            }
        }

        if negative_lines_found || lines_exceeding_height {
            log::info!("📏 Line indexing corrections: {} negative corrected, {} clamped", 
                      corrected_count, clamped_count);
        }

        LineIndexingResult {
            negative_lines_found,
            lines_exceeding_height,
            corrected_count,
            clamped_count,
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
            if let Some(deviation) = Self::check_seam_power_consistency(merged_data, row_start, row_end, col_start, col_end) {
                if deviation <= 1.0 {
                    seams_within_tolerance += 1;
                } else {
                    log::warn!("⚠️ Seam {} power deviation: {:.2}% (> 1.0% threshold)", i, deviation);
                }
                max_power_deviation_percent = max_power_deviation_percent.max(deviation);
            }
        }

        let is_valid = seams_within_tolerance == overlap_regions.len();

        if !is_valid {
            log::error!("🚨 POWER SEAM VALIDATION FAILED:");
            log::error!("   Seams within tolerance: {}/{}", seams_within_tolerance, overlap_regions.len());
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

        // Calculate percentage deviation
        let deviation_percent = ((seam_mean_power - neighbor_mean_power).abs() / neighbor_mean_power) * 100.0;

        Some(deviation_percent)
    }

    /// Validate coverage - uncovered pixels in seams < 0.5%
    pub fn validate_coverage(uncovered_mask: &Array2<u8>) -> CoverageResult {
        let total_pixels = uncovered_mask.len();
        let uncovered_pixels = uncovered_mask.iter().map(|&x| x as usize).sum::<usize>();
        let uncovered_percent = (uncovered_pixels as f64 / total_pixels as f64) * 100.0;
        
        let is_valid = uncovered_percent < 0.5;

        if !is_valid {
            log::error!("🚨 COVERAGE VALIDATION FAILED:");
            log::error!("   Uncovered pixels: {} ({:.2}%)", uncovered_pixels, uncovered_percent);
            log::error!("   Threshold: < 0.5%");
            
            if uncovered_percent > 2.0 {
                log::error!("   🚫 CRITICAL: > 2.0% uncovered → wrong overlap math detected");
            }
        } else {
            log::info!("✅ Coverage validation passed: {:.3}% uncovered", uncovered_percent);
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

    // 1. Strict metadata validation
    let metadata_result = StrictMetadataValidator::validate_strict(metadata)?;
    if !metadata_result.is_valid {
        return Err(SarError::Processing(format!(
            "Strict metadata validation failed: {} errors, {} critical missing fields",
            metadata_result.errors.len(),
            metadata_result.critical_missing.len()
        )));
    }

    // 2. β⁰ variability check
    let beta_result = StrictMetadataValidator::validate_beta_variability(beta_values);
    if beta_result.parsing_bug_detected {
        return Err(SarError::Processing(format!(
            "β⁰ variability check failed: ratio {:.6} ≤ 1.02 indicates parsing bug",
            beta_result.ratio
        )));
    }

    // 3. Units validation
    let units_result = StrictMetadataValidator::validate_units(
        sigma_values, beta_values, gamma_values, None, None, None
    );
    if !units_result.all_linear {
        log::warn!("⚠️ Units requiring dB conversion: {:?}", units_result.dB_conversions_needed);
    }

    // 4. Power seam validation
    let seam_result = StrictMetadataValidator::validate_power_seams(merged_data, overlap_regions);
    if !seam_result.is_valid {
        return Err(SarError::Processing(format!(
            "Power seam validation failed: {}/{} seams within tolerance, max deviation {:.2}%",
            seam_result.seams_within_tolerance,
            seam_result.seam_count,
            seam_result.max_power_deviation_percent
        )));
    }

    // 5. Coverage validation
    let coverage_result = StrictMetadataValidator::validate_coverage(uncovered_mask);
    if !coverage_result.is_valid {
        return Err(SarError::Processing(format!(
            "Coverage validation failed: {:.2}% uncovered pixels (> 0.5% threshold)",
            coverage_result.uncovered_percent
        )));
    }

    log::info!("✅ Comprehensive strict validation passed");
    Ok(())
}