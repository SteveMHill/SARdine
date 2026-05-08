/*!
 * STEP-2 Diagnostics and Invariants Module
 *
 * This module implements scientific invariants and diagnostic metrics for the
 * TOPSAR IW deburst pipeline. It provides gates to prevent regressions and
 * tools to diagnose remaining sources of burst-to-burst power variation.
 *
 * Scope:
 * - Coverage metrics (pixel accounting)
 * - Contribution histograms (burst overlap validation)
 * - Overlap radiometry agreement (pre-blend power comparison)
 * - Burst power CV with proper masking
 * - DC/FM interpolation bracketing validation
 * - Deramp effectiveness measurement
 * - Calibration LUT integrity checks
 *
 * Usage:
 * - Enable via SARDINE_STEP2_DIAGNOSTICS=1
 * - Configure thresholds via env or DiagnosticsConfig
 * - Outputs JSON summary: STEP2_DIAGNOSTICS_SUMMARY.json
 */

use crate::types::{SarError, SarResult};
use serde::{Deserialize, Serialize};
use std::path::PathBuf;

// ============================================================================
// CONFIGURATION
// ============================================================================

/// STEP-2 diagnostics configuration
#[derive(Debug, Clone)]
pub struct DiagnosticsConfig {
    /// Enable diagnostics (default: from env SARDINE_STEP2_DIAGNOSTICS)
    pub enabled: bool,
    
    /// Minimum coverage fraction (default: 0.99)
    pub min_coverage_fraction: f64,
    
    /// Maximum overlap delta in dB (default: 1.0)
    pub max_overlap_delta_db: f64,
    
    /// Maximum burst power CV for homogeneous ROI (default: 0.08)
    pub max_burst_cv: f64,
    
    /// Warning threshold for CV (default: 0.15)
    pub warn_burst_cv: f64,
    
    /// Maximum DC value in Hz (sanity check, default: 200.0)
    pub max_dc_hz: f64,
    
    /// Minimum deramp reduction factor (default: 5.0)
    pub min_deramp_reduction: f64,
    
    /// Range edge margin for masking (fraction, default: 0.05 = 5%)
    pub range_edge_margin: f64,
    
    /// Azimuth edge margin for masking (fraction, default: 0.01 = 1%)
    pub azimuth_edge_margin: f64,
    
    /// Strict mode: panic on invariant failures (default: false)
    pub strict_mode: bool,
    
    /// Output directory for JSON summary (default: current dir)
    pub output_dir: PathBuf,
}

impl Default for DiagnosticsConfig {
    fn default() -> Self {
        let enabled = std::env::var("SARDINE_STEP2_DIAGNOSTICS")
            .unwrap_or_else(|_| "0".to_string())
            == "1";
        
        let strict_mode = std::env::var("SARDINE_DIAGNOSTICS_STRICT")
            .unwrap_or_else(|_| "0".to_string())
            == "1";
        
        Self {
            enabled,
            min_coverage_fraction: 0.99,
            max_overlap_delta_db: 1.0,
            max_burst_cv: 0.08,
            warn_burst_cv: 0.15,
            max_dc_hz: 200.0,
            min_deramp_reduction: 5.0,
            range_edge_margin: 0.05,
            azimuth_edge_margin: 0.01,
            strict_mode,
            output_dir: PathBuf::from("."),
        }
    }
}

// ============================================================================
// DIAGNOSTIC TYPES
// ============================================================================

/// Coverage metrics for a debursted array
/// 
/// IMPORTANT: For Sentinel-1 TOPSAR IW SLC, each burst has invalid range margins
/// (first_valid_sample/last_valid_sample) and azimuth edges. Therefore, full-canvas
/// fill is expected to be ~92-95% and must NOT be used as a correctness gate.
/// 
/// We compute TWO metrics:
/// - full_canvas: Informational only (shows actual canvas utilization)
/// - valid_window: Gated metric computed within the canonical valid range window
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CoverageMetrics {
    // Full canvas metrics (informational, non-gating)
    pub total_pixels_full: usize,
    pub covered_pixels_full: usize,
    pub fraction_full: f64,
    
    // Valid window metrics (gating)
    pub total_pixels_valid: Option<usize>,
    pub covered_pixels_valid: Option<usize>,
    pub fraction_valid: Option<f64>,
    
    // Expected fill range from metadata
    pub expected_fill_range: Option<f64>,
    
    // Diagnostics
    pub has_interior_holes: bool,
    pub valid_window_available: bool,
    
    // Pass/fail (based on valid window if available, else null)
    pub pass: bool,
    pub fail_reason: Option<String>,
}

/// Contribution count statistics
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ContributionHistogram {
    pub count_0: usize, // No burst contributed
    pub count_1: usize, // Single burst
    pub count_2: usize, // Two bursts (overlap)
    pub count_invalid: usize, // More than 2 (ERROR)
    pub pass: bool,
}

/// Overlap zone radiometry for a burst pair
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct OverlapMetrics {
    pub burst_a_idx: usize,
    pub burst_b_idx: usize,
    pub overlap_start_line: usize,
    pub overlap_end_line: usize,
    pub mean_power_a: f64,
    pub mean_power_b: f64,
    pub ratio_b_over_a: f64,
    pub delta_db: f64,
    pub pass: bool,
}

/// Per-burst power statistics with masking
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BurstPowerStats {
    pub burst_idx: usize,
    pub mean_power: f64,
    pub std_power: f64,
    pub median_power: f64,
    pub valid_pixel_count: usize,
    pub deviation_from_global_pct: f64,
}

/// Burst power CV summary
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BurstPowerCV {
    pub burst_stats: Vec<BurstPowerStats>,
    pub global_mean: f64,
    pub global_std: f64,
    pub cv: f64, // coefficient of variation (std / mean)
    pub pass: bool,
}

/// DC/FM interpolation diagnostics for a single time point
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DCFMBracketDiag {
    pub target_time: f64,
    pub lower_time: f64,
    pub upper_time: f64,
    pub weight: f64,
    pub is_extrapolated: bool,
    pub dc_value_hz: f64,
    pub fm_value_hz_per_s: f64,
    pub sanity_pass: bool,
}

/// Deramp effectiveness metrics
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DerampEffectiveness {
    pub sample_count: usize,
    pub median_slope_before: f64,
    pub median_slope_after: f64,
    pub reduction_factor: f64,
    pub pass: bool,
}

/// Calibration LUT integrity check
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CalibrationLUTDiag {
    pub calibration_type: String, // "sigma0", "beta0", "gamma0"
    pub lut_array_name: String,   // Which array was actually used
    pub lut_min: f64,
    pub lut_median: f64,
    pub lut_max: f64,
    pub has_invalid: bool,
    pub has_zeros: bool,
    pub pass: bool,
}

/// Per-subswath diagnostics summary
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SubswathDiagnostics {
    pub subswath_name: String,
    pub polarization: String,
    pub burst_count: usize,
    
    // Coverage
    pub coverage: CoverageMetrics,
    pub contribution_histogram: ContributionHistogram,
    
    // Overlap radiometry
    pub overlap_metrics: Vec<OverlapMetrics>,
    
    // Burst power
    pub burst_power_cv: BurstPowerCV,
    
    // DC/FM bracketing (sampled time points)
    pub dc_fm_diagnostics: Vec<DCFMBracketDiag>,
    
    // Deramp
    pub deramp_effectiveness: Option<DerampEffectiveness>,
    
    // Calibration
    pub calibration_lut: Option<CalibrationLUTDiag>,
}

/// Top-level STEP-2 diagnostics summary
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Step2DiagnosticsSummary {
    pub product_id: String,
    pub timestamp: String,
    pub config: DiagnosticsConfigSummary,
    pub subswaths: Vec<SubswathDiagnostics>,
    pub overall_pass: bool,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DiagnosticsConfigSummary {
    pub min_coverage_fraction: f64,
    pub max_overlap_delta_db: f64,
    pub max_burst_cv: f64,
    pub strict_mode: bool,
}

// ============================================================================
// DIAGNOSTICS COLLECTOR
// ============================================================================

/// Main diagnostics collector for STEP-2
pub struct Step2Diagnostics {
    pub config: DiagnosticsConfig,
    subswath_diags: Vec<SubswathDiagnostics>,
    product_id: String,
}

impl Step2Diagnostics {
    /// Create a new diagnostics collector
    pub fn new(product_id: String, config: DiagnosticsConfig) -> Self {
        if config.enabled {
            log::info!("STEP2|INIT: Diagnostics ENABLED (strict={})", config.strict_mode);
        }
        
        Self {
            config,
            subswath_diags: Vec::new(),
            product_id,
        }
    }
    
    /// Check if diagnostics are enabled
    pub fn is_enabled(&self) -> bool {
        self.config.enabled
    }
    
    /// Add subswath diagnostics
    pub fn add_subswath(&mut self, diag: SubswathDiagnostics) {
        self.subswath_diags.push(diag);
    }
    
    /// Generate final summary
    pub fn summarize(&self) -> Step2DiagnosticsSummary {
        let overall_pass = self.subswath_diags.iter().all(|sw| {
            sw.coverage.pass
                && sw.contribution_histogram.pass
                && sw.overlap_metrics.iter().all(|o| o.pass)
                && sw.burst_power_cv.pass
                && sw.dc_fm_diagnostics.iter().all(|d| d.sanity_pass)
                && sw.deramp_effectiveness.as_ref().map_or(true, |d| d.pass)
                && sw.calibration_lut.as_ref().map_or(true, |c| c.pass)
        });
        
        Step2DiagnosticsSummary {
            product_id: self.product_id.clone(),
            timestamp: chrono::Utc::now().to_rfc3339(),
            config: DiagnosticsConfigSummary {
                min_coverage_fraction: self.config.min_coverage_fraction,
                max_overlap_delta_db: self.config.max_overlap_delta_db,
                max_burst_cv: self.config.max_burst_cv,
                strict_mode: self.config.strict_mode,
            },
            subswaths: self.subswath_diags.clone(),
            overall_pass,
        }
    }
    
    /// Write summary to JSON file
    pub fn write_json(&self, filename: &str) -> SarResult<()> {
        let summary = self.summarize();
        let path = self.config.output_dir.join(filename);
        
        let json = serde_json::to_string_pretty(&summary)
            .map_err(|e| SarError::Processing(format!("Failed to serialize diagnostics: {}", e)))?;
        
        std::fs::write(&path, json)?;
        
        log::info!("STEP2|OUTPUT: Wrote diagnostics to {}", path.display());
        
        if !summary.overall_pass && self.config.strict_mode {
            return Err(SarError::Validation(
                "STEP-2 diagnostics FAILED in strict mode".to_string(),
            ));
        }
        
        Ok(())
    }
}

// ============================================================================
// INVARIANT COMPUTATION FUNCTIONS
// ============================================================================

/// Compute coverage metrics from contribution count array
///
/// Computes TWO coverage metrics:
/// 1. Full canvas: Total coverage (informational, never gates)
/// 2. Valid window: Coverage within canonical valid range (gating metric)
///
/// # Arguments
/// * `contribution_count` - Array where each pixel = number of bursts that contributed
/// * `valid_range_window` - Optional (first_valid_col, last_valid_col) from burst metadata
/// * `expected_fill` - Optional expected fill fraction from metadata analysis
/// * `config` - Diagnostics configuration
pub fn compute_coverage_metrics(
    contribution_count: &ndarray::Array2<u8>,
    valid_range_window: Option<(usize, usize)>,
    expected_fill: Option<f64>,
    config: &DiagnosticsConfig,
) -> CoverageMetrics {
    // Full canvas metrics (informational)
    let total_pixels_full = contribution_count.len();
    let covered_pixels_full = contribution_count.iter().filter(|&&c| c > 0).count();
    let fraction_full = covered_pixels_full as f64 / total_pixels_full as f64;
    
    log::info!(
        "STEP2|COVERAGE|FULL: {}/{} pixels ({:.2}%) - informational only",
        covered_pixels_full,
        total_pixels_full,
        fraction_full * 100.0
    );
    
    // Valid window metrics (gating)
    let (total_pixels_valid, covered_pixels_valid, fraction_valid, valid_window_available) =
        if let Some((first_col, last_col)) = valid_range_window {
            if last_col > first_col && last_col <= contribution_count.shape()[1] {
                // Crop to valid window
                let valid_window = contribution_count.slice(ndarray::s![.., first_col..last_col]);
                let total_valid = valid_window.len();
                let covered_valid = valid_window.iter().filter(|&&c| c > 0).count();
                let frac_valid = covered_valid as f64 / total_valid as f64;
                
                log::info!(
                    "STEP2|COVERAGE|VALID: {}/{} pixels ({:.2}%) in window [{}:{}] - GATING",
                    covered_valid,
                    total_valid,
                    frac_valid * 100.0,
                    first_col,
                    last_col
                );
                
                (Some(total_valid), Some(covered_valid), Some(frac_valid), true)
            } else {
                log::warn!("STEP2|COVERAGE: Invalid window range [{}, {}], using full canvas", first_col, last_col);
                (None, None, None, false)
            }
        } else {
            log::warn!("STEP2|COVERAGE: No valid range window metadata, cannot compute gating metric");
            (None, None, None, false)
        };
    
    // Check for interior holes
    let has_interior_holes = detect_interior_holes(contribution_count);
    
    // Pass/fail logic: only gate on valid window if available
    let (pass, fail_reason) = if let Some(frac_valid) = fraction_valid {
        if frac_valid >= config.min_coverage_fraction && !has_interior_holes {
            (true, None)
        } else if has_interior_holes {
            (false, Some("Interior holes detected".to_string()))
        } else {
            (false, Some(format!("Valid window coverage {:.2}% < {:.2}%", 
                frac_valid * 100.0, config.min_coverage_fraction * 100.0)))
        }
    } else {
        // No valid window - cannot gate, mark as failed with reason
        (false, Some("Missing valid range metadata".to_string()))
    };
    
    if !pass {
        log::warn!(
            "STEP2|COVERAGE: FAIL - {}",
            fail_reason.as_ref().unwrap_or(&"Unknown".to_string())
        );
        if config.strict_mode {
            log::error!("STEP2|COVERAGE: Strict mode enabled - this would fail the pipeline");
        }
    }
    
    CoverageMetrics {
        total_pixels_full,
        covered_pixels_full,
        fraction_full,
        total_pixels_valid,
        covered_pixels_valid,
        fraction_valid,
        expected_fill_range: expected_fill,
        has_interior_holes,
        valid_window_available,
        pass,
        fail_reason,
    }
}

/// Compute contribution histogram
pub fn compute_contribution_histogram(
    contribution_count: &ndarray::Array2<u8>,
    config: &DiagnosticsConfig,
) -> ContributionHistogram {
    let mut count_0 = 0;
    let mut count_1 = 0;
    let mut count_2 = 0;
    let mut count_invalid = 0;
    
    for &c in contribution_count.iter() {
        match c {
            0 => count_0 += 1,
            1 => count_1 += 1,
            2 => count_2 += 1,
            _ => count_invalid += 1,
        }
    }
    
    let pass = count_invalid == 0;
    
    log::info!(
        "STEP2|CONTRIB: Histogram=[0:{}, 1:{}, 2:{}, invalid:{}], PASS={}",
        count_0,
        count_1,
        count_2,
        count_invalid,
        pass
    );
    
    if !pass {
        log::error!("STEP2|CONTRIB: Found {} pixels with >2 burst contributions!", count_invalid);
        if config.strict_mode {
            panic!("STEP2|CONTRIB: Invariant violation - contribution count > 2");
        }
    }
    
    ContributionHistogram {
        count_0,
        count_1,
        count_2,
        count_invalid,
        pass,
    }
}

/// Detect interior holes (uncovered pixels not on boundary)
fn detect_interior_holes(contribution_count: &ndarray::Array2<u8>) -> bool {
    let shape = contribution_count.shape();
    let nrows = shape[0];
    let ncols = shape[1];
    
    if nrows < 3 || ncols < 3 {
        return false;
    }
    
    // Check interior region (exclude first/last row/col)
    for i in 1..(nrows - 1) {
        for j in 1..(ncols - 1) {
            if contribution_count[[i, j]] == 0 {
                // Check if surrounded by covered pixels
                let neighbors = [
                    contribution_count[[i - 1, j]],
                    contribution_count[[i + 1, j]],
                    contribution_count[[i, j - 1]],
                    contribution_count[[i, j + 1]],
                ];
                
                if neighbors.iter().all(|&n| n > 0) {
                    log::warn!("STEP2|COVERAGE: Interior hole detected at ({}, {})", i, j);
                    return true;
                }
            }
        }
    }
    
    false
}

/// Compute overlap zone radiometry before blending
///
/// # Arguments
/// * `burst_a_data` - Burst A intensity data in overlap region
/// * `burst_b_data` - Burst B intensity data in overlap region
/// * `burst_a_idx`, `burst_b_idx` - Burst indices
/// * `overlap_lines` - (start_line, end_line) in output coordinates
/// * `config` - Diagnostics configuration
pub fn compute_overlap_metrics(
    burst_a_data: &ndarray::Array2<f32>,
    burst_b_data: &ndarray::Array2<f32>,
    burst_a_idx: usize,
    burst_b_idx: usize,
    overlap_lines: (usize, usize),
    config: &DiagnosticsConfig,
) -> OverlapMetrics {
    // Compute mean power (exclude NaN/inf and zeros)
    let mean_a = compute_valid_mean(burst_a_data);
    let mean_b = compute_valid_mean(burst_b_data);
    
    let ratio = if mean_a > 0.0 { mean_b / mean_a } else { 0.0 };
    let delta_db = if mean_a > 0.0 && mean_b > 0.0 {
        10.0 * (mean_b / mean_a).log10()
    } else {
        f64::NAN
    };
    
    let pass = delta_db.is_finite() && delta_db.abs() <= config.max_overlap_delta_db;
    
    log::info!(
        "STEP2|OVERLAP: Burst {}-{} overlap lines [{}-{}]: mean_A={:.2}, mean_B={:.2}, Δ={:.2}dB, PASS={}",
        burst_a_idx,
        burst_b_idx,
        overlap_lines.0,
        overlap_lines.1,
        mean_a,
        mean_b,
        delta_db,
        pass
    );
    
    if !pass && delta_db.is_finite() {
        log::warn!(
            "STEP2|OVERLAP: Overlap delta {:.2}dB exceeds threshold {:.2}dB",
            delta_db,
            config.max_overlap_delta_db
        );
    }
    
    OverlapMetrics {
        burst_a_idx,
        burst_b_idx,
        overlap_start_line: overlap_lines.0,
        overlap_end_line: overlap_lines.1,
        mean_power_a: mean_a,
        mean_power_b: mean_b,
        ratio_b_over_a: ratio,
        delta_db,
        pass,
    }
}

/// Compute valid mean (exclude NaN, inf, and zeros)
fn compute_valid_mean(data: &ndarray::Array2<f32>) -> f64 {
    let valid_values: Vec<f64> = data
        .iter()
        .filter(|&&v| v.is_finite() && v > 0.0)
        .map(|&v| v as f64)
        .collect();
    
    if valid_values.is_empty() {
        return 0.0;
    }
    
    valid_values.iter().sum::<f64>() / valid_values.len() as f64
}

/// Compute per-burst power statistics with edge masking
///
/// # Arguments
/// * `burst_data` - Single burst intensity data
/// * `burst_idx` - Burst index
/// * `first_valid_sample` - First valid range sample
/// * `last_valid_sample` - Last valid range sample
/// * `config` - Diagnostics configuration
pub fn compute_burst_power_stats(
    burst_data: &ndarray::Array2<f32>,
    burst_idx: usize,
    first_valid_sample: usize,
    last_valid_sample: usize,
    config: &DiagnosticsConfig,
) -> BurstPowerStats {
    let shape = burst_data.shape();
    let nlines = shape[0];
    let nsamples = shape[1];
    
    // Apply margins
    let range_margin_pixels = ((nsamples as f64) * config.range_edge_margin) as usize;
    let azimuth_margin_lines = ((nlines as f64) * config.azimuth_edge_margin) as usize;
    
    let range_start = first_valid_sample.max(range_margin_pixels);
    let range_end = last_valid_sample.min(nsamples - range_margin_pixels);
    let az_start = azimuth_margin_lines;
    let az_end = nlines.saturating_sub(azimuth_margin_lines);
    
    // Collect valid values in masked ROI
    let mut values = Vec::new();
    for i in az_start..az_end {
        for j in range_start..range_end {
            let val = burst_data[[i, j]];
            if val.is_finite() && val > 0.0 {
                values.push(val as f64);
            }
        }
    }
    
    let valid_pixel_count = values.len();
    let mean_power = if !values.is_empty() {
        values.iter().sum::<f64>() / values.len() as f64
    } else {
        0.0
    };
    
    let std_power = if values.len() > 1 {
        let variance = values.iter()
            .map(|&v| (v - mean_power).powi(2))
            .sum::<f64>() / (values.len() - 1) as f64;
        variance.sqrt()
    } else {
        0.0
    };
    
    // Median
    values.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let median_power = if !values.is_empty() {
        values[values.len() / 2]
    } else {
        0.0
    };
    
    BurstPowerStats {
        burst_idx,
        mean_power,
        std_power,
        median_power,
        valid_pixel_count,
        deviation_from_global_pct: 0.0, // Will be filled by CV computation
    }
}

/// Compute burst power CV from per-burst stats
pub fn compute_burst_power_cv(
    mut burst_stats: Vec<BurstPowerStats>,
    config: &DiagnosticsConfig,
) -> BurstPowerCV {
    if burst_stats.is_empty() {
        return BurstPowerCV {
            burst_stats: vec![],
            global_mean: 0.0,
            global_std: 0.0,
            cv: 0.0,
            pass: false,
        };
    }
    
    // Global mean
    let global_mean = burst_stats.iter()
        .map(|s| s.mean_power)
        .sum::<f64>() / burst_stats.len() as f64;
    
    // Global std (population std of burst means)
    let global_std = if burst_stats.len() > 1 {
        let variance = burst_stats.iter()
            .map(|s| (s.mean_power - global_mean).powi(2))
            .sum::<f64>() / burst_stats.len() as f64;
        variance.sqrt()
    } else {
        0.0
    };
    
    let cv = if global_mean > 0.0 {
        global_std / global_mean
    } else {
        0.0
    };
    
    // Fill in deviation percentages
    for stat in &mut burst_stats {
        stat.deviation_from_global_pct = if global_mean > 0.0 {
            ((stat.mean_power - global_mean) / global_mean) * 100.0
        } else {
            0.0
        };
    }
    
    let pass = cv <= config.max_burst_cv;
    
    log::info!(
        "STEP2|BURST_CV: {} bursts, global_mean={:.2}, CV={:.4} ({:.2}%), PASS={}",
        burst_stats.len(),
        global_mean,
        cv,
        cv * 100.0,
        pass
    );
    
    if cv > config.warn_burst_cv {
        log::warn!(
            "STEP2|BURST_CV: High CV={:.2}% exceeds warning threshold {:.2}%",
            cv * 100.0,
            config.warn_burst_cv * 100.0
        );
        
        // Log per-burst deviations
        for stat in &burst_stats {
            log::warn!(
                "  Burst {}: mean={:.2}, dev={:+.1}%",
                stat.burst_idx,
                stat.mean_power,
                stat.deviation_from_global_pct
            );
        }
    }
    
    BurstPowerCV {
        burst_stats,
        global_mean,
        global_std,
        cv,
        pass,
    }
}

/// Validate DC/FM interpolation bracketing
///
/// # Arguments
/// * `target_time` - Target azimuth time
/// * `lower_time`, `upper_time` - Bracketing estimate times
/// * `weight` - Interpolation weight
/// * `dc_value` - Interpolated DC value in Hz
/// * `fm_value` - Interpolated FM value in Hz/s
/// * `config` - Diagnostics configuration
pub fn validate_dc_fm_bracket(
    target_time: f64,
    lower_time: f64,
    upper_time: f64,
    weight: f64,
    dc_value: f64,
    fm_value: f64,
    config: &DiagnosticsConfig,
) -> DCFMBracketDiag {
    // Check if extrapolated
    let is_extrapolated = target_time < lower_time || target_time > upper_time;
    
    // Sanity check DC magnitude
    let dc_sanity = dc_value.abs() < config.max_dc_hz;
    
    let sanity_pass = !is_extrapolated && dc_sanity;
    
    if is_extrapolated {
        log::warn!(
            "STEP2|DC_BRACKET: EXTRAPOLATION detected at t={:.3}: [{:.3}, {:.3}]",
            target_time,
            lower_time,
            upper_time
        );
    }
    
    if !dc_sanity {
        log::error!(
            "STEP2|DC_BRACKET: DC value {:.1} Hz exceeds sanity limit {} Hz",
            dc_value,
            config.max_dc_hz
        );
    }
    
    log::debug!(
        "STEP2|DC_BRACKET: t={:.3}, bracket=[{:.3}, {:.3}], w={:.3}, DC={:.1}Hz, FM={:.1}Hz/s, PASS={}",
        target_time,
        lower_time,
        upper_time,
        weight,
        dc_value,
        fm_value,
        sanity_pass
    );
    
    DCFMBracketDiag {
        target_time,
        lower_time,
        upper_time,
        weight,
        is_extrapolated,
        dc_value_hz: dc_value,
        fm_value_hz_per_s: fm_value,
        sanity_pass,
    }
}

/// Measure deramp effectiveness (phase slope reduction)
///
/// # Arguments
/// * `phase_before` - Phase array before deramp (sampled)
/// * `phase_after` - Phase array after deramp (sampled)
/// * `config` - Diagnostics configuration
pub fn measure_deramp_effectiveness(
    phase_before: &ndarray::Array1<f32>,
    phase_after: &ndarray::Array1<f32>,
    config: &DiagnosticsConfig,
) -> DerampEffectiveness {
    let sample_count = phase_before.len();
    
    // Compute median absolute phase slope
    let slope_before = compute_median_phase_slope(phase_before);
    let slope_after = compute_median_phase_slope(phase_after);
    
    let reduction_factor = if slope_after > 0.0 {
        slope_before / slope_after
    } else {
        f64::INFINITY
    };
    
    let pass = reduction_factor >= config.min_deramp_reduction;
    
    log::info!(
        "STEP2|DERAMP: {} samples, slope before={:.3}, after={:.3}, reduction={:.1}x, PASS={}",
        sample_count,
        slope_before,
        slope_after,
        reduction_factor,
        pass
    );
    
    if !pass {
        log::warn!(
            "STEP2|DERAMP: Low reduction factor {:.1}x < threshold {:.1}x",
            reduction_factor,
            config.min_deramp_reduction
        );
    }
    
    DerampEffectiveness {
        sample_count,
        median_slope_before: slope_before,
        median_slope_after: slope_after,
        reduction_factor,
        pass,
    }
}

/// Compute median absolute phase slope
fn compute_median_phase_slope(phase: &ndarray::Array1<f32>) -> f64 {
    if phase.len() < 2 {
        return 0.0;
    }
    
    let mut slopes = Vec::new();
    for i in 1..phase.len() {
        let slope = (phase[i] - phase[i - 1]).abs();
        if slope.is_finite() {
            slopes.push(slope as f64);
        }
    }
    
    if slopes.is_empty() {
        return 0.0;
    }
    
    slopes.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    slopes[slopes.len() / 2]
}

/// Validate calibration LUT integrity
///
/// # Arguments
/// * `calibration_type` - "sigma0", "beta0", or "gamma0"
/// * `lut_array` - The actual LUT array used
/// * `lut_array_name` - Name of the array (for logging)
/// * `config` - Diagnostics configuration
pub fn validate_calibration_lut(
    calibration_type: &str,
    lut_array: &ndarray::Array2<f32>,
    lut_array_name: &str,
    config: &DiagnosticsConfig,
) -> CalibrationLUTDiag {
    // Expected LUT array names
    let expected_array = match calibration_type {
        "sigma0" => "sigma",
        "beta0" => "beta",
        "gamma0" => "gamma",
        _ => "unknown",
    };
    
    // Check for name mismatch
    let name_match = lut_array_name.to_lowercase().contains(expected_array);
    if !name_match {
        log::error!(
            "STEP2|CAL_LUT: Calibration type '{}' expects '{}' array, got '{}'",
            calibration_type,
            expected_array,
            lut_array_name
        );
    }
    
    // Compute stats on valid pixels
    let mut valid_values: Vec<f64> = lut_array
        .iter()
        .filter(|&&v| v.is_finite() && v > 0.0)
        .map(|&v| v as f64)
        .collect();
    
    let has_invalid = lut_array.iter().any(|&v| !v.is_finite());
    let has_zeros = lut_array.iter().any(|&v| v == 0.0);
    
    valid_values.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    
    let lut_min = valid_values.first().copied().unwrap_or(0.0);
    let lut_median = if !valid_values.is_empty() {
        valid_values[valid_values.len() / 2]
    } else {
        0.0
    };
    let lut_max = valid_values.last().copied().unwrap_or(0.0);
    
    let pass = name_match && !has_invalid && lut_median > 0.0;
    
    log::info!(
        "STEP2|CAL_LUT: type='{}', array='{}', range=[{:.2e}, {:.2e}, {:.2e}], invalid={}, zeros={}, PASS={}",
        calibration_type,
        lut_array_name,
        lut_min,
        lut_median,
        lut_max,
        has_invalid,
        has_zeros,
        pass
    );
    
    if !pass && config.strict_mode {
        panic!(
            "STEP2|CAL_LUT: Calibration LUT validation failed for type '{}' with array '{}'",
            calibration_type, lut_array_name
        );
    }
    
    CalibrationLUTDiag {
        calibration_type: calibration_type.to_string(),
        lut_array_name: lut_array_name.to_string(),
        lut_min,
        lut_median,
        lut_max,
        has_invalid,
        has_zeros,
        pass,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use ndarray::Array2;
    
    #[test]
    fn test_contribution_histogram() {
        let config = DiagnosticsConfig::default();
        
        // Valid histogram: only 0, 1, 2
        let mut data = Array2::zeros((10, 10));
        data[[0, 0]] = 1;
        data[[1, 1]] = 2;
        
        let hist = compute_contribution_histogram(&data, &config);
        assert!(hist.pass);
        assert_eq!(hist.count_invalid, 0);
        
        // Invalid: contribution count > 2
        data[[2, 2]] = 3;
        let hist = compute_contribution_histogram(&data, &config);
        assert!(!hist.pass);
        assert_eq!(hist.count_invalid, 1);
    }
    
    #[test]
    fn test_overlap_metrics() {
        let config = DiagnosticsConfig::default();
        
        // Similar power
        let burst_a = Array2::from_elem((10, 10), 100.0);
        let burst_b = Array2::from_elem((10, 10), 105.0);
        
        let metrics = compute_overlap_metrics(
            &burst_a,
            &burst_b,
            0,
            1,
            (0, 10),
            &config,
        );
        
        assert!(metrics.pass);
        assert!(metrics.delta_db.abs() < 1.0);
        
        // Different power (should fail)
        let burst_c = Array2::from_elem((10, 10), 200.0);
        let metrics = compute_overlap_metrics(
            &burst_a,
            &burst_c,
            0,
            2,
            (0, 10),
            &config,
        );
        
        assert!(!metrics.pass);
        assert!(metrics.delta_db.abs() > 1.0);
    }
    
    #[test]
    fn test_dc_fm_bracket_validation() {
        let config = DiagnosticsConfig::default();
        
        // Good bracketing
        let diag = validate_dc_fm_bracket(
            50.0,
            45.0,
            55.0,
            0.5,
            30.0,
            -2100.0,
            &config,
        );
        
        assert!(!diag.is_extrapolated);
        assert!(diag.sanity_pass);
        
        // Extrapolation
        let diag = validate_dc_fm_bracket(
            60.0,
            45.0,
            55.0,
            1.5,
            30.0,
            -2100.0,
            &config,
        );
        
        assert!(diag.is_extrapolated);
        assert!(!diag.sanity_pass);
        
        // DC out of range
        let diag = validate_dc_fm_bracket(
            50.0,
            45.0,
            55.0,
            0.5,
            250.0, // > 200 Hz
            -2100.0,
            &config,
        );
        
        assert!(!diag.sanity_pass);
    }
    
    #[test]
    fn test_calibration_lut_validation() {
        let config = DiagnosticsConfig::default();
        
        // Correct usage: sigma0 with sigma array
        let sigma_lut = Array2::from_elem((10, 10), 0.5);
        let diag = validate_calibration_lut(
            "sigma0",
            &sigma_lut,
            "sigma",
            &config,
        );
        assert!(diag.pass);
        
        // Wrong usage: sigma0 with beta array
        let beta_lut = Array2::from_elem((10, 10), 0.3);
        let diag = validate_calibration_lut(
            "sigma0",
            &beta_lut,
            "beta",
            &config,
        );
        assert!(!diag.pass);
        
        // Invalid values
        let mut invalid_lut = Array2::from_elem((10, 10), 0.5);
        invalid_lut[[0, 0]] = f32::NAN;
        let diag = validate_calibration_lut(
            "beta0",
            &invalid_lut,
            "beta",
            &config,
        );
        assert!(!diag.pass);
    }
    
    #[test]
    fn test_coverage_metrics_two_metric_approach() {
        let config = DiagnosticsConfig::default();
        
        // Simulate TOPSAR scenario: full canvas has ~92% fill, valid window has 99%
        // Create 100x1000 array with contribution counts
        let mut data = Array2::zeros((100, 1000));
        
        // Fill valid window [50:950] with high coverage
        for i in 0..100 {
            for j in 50..950 {
                data[[i, j]] = 1; // All covered
            }
        }
        
        // Leave edges [0:50] and [950:1000] mostly empty (typical TOPSAR invalidity)
        // This gives us ~90% full canvas, but 100% valid window
        
        let valid_window = Some((50, 950));
        let expected_fill = Some(0.90); // 900/1000 columns
        
        let metrics = compute_coverage_metrics(&data, valid_window, expected_fill, &config);
        
        // Full canvas should be ~90% (informational)
        assert!((metrics.fraction_full - 0.90).abs() < 0.01);
        assert_eq!(metrics.total_pixels_full, 100 * 1000);
        
        // Valid window should be 100% (gating)
        assert!(metrics.valid_window_available);
        assert_eq!(metrics.fraction_valid, Some(1.0));
        assert_eq!(metrics.total_pixels_valid, Some(100 * 900));
        
        // Should PASS based on valid window, not full canvas
        assert!(metrics.pass);
        assert!(metrics.fail_reason.is_none());
        assert_eq!(metrics.expected_fill_range, Some(0.90));
        
        // Now test failure case: low coverage in valid window
        let mut sparse_data = Array2::zeros((100, 1000));
        for i in 0..100 {
            for j in 50..150 { // Only 100 columns in valid window
                sparse_data[[i, j]] = 1;
            }
        }
        
        let sparse_metrics = compute_coverage_metrics(&sparse_data, valid_window, expected_fill, &config);
        
        // Should FAIL because valid window coverage < 99%
        assert!(!sparse_metrics.pass);
        assert!(sparse_metrics.fail_reason.is_some());
        
        // Without valid window metadata, should fail with appropriate reason
        let no_window_metrics = compute_coverage_metrics(&data, None, None, &config);
        assert!(!no_window_metrics.pass);
        assert_eq!(no_window_metrics.fail_reason, Some("Missing valid range metadata".to_string()));
    }
}
