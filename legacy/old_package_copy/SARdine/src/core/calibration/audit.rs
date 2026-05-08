use std::env;
use std::fs;
use std::path::{Path, PathBuf};

use ndarray::Array2;
use serde::Serialize;

use crate::core::calibration::{
    CalibrationCoefficients, CalibrationType, ValidSampleRanges,
};
use crate::core::calibration::model::MIN_VALID_POWER;
use crate::types::{SarError, SarResult};

#[derive(Serialize)]
struct ScalarStats {
    min: f64,
    max: f64,
    mean: f64,
    finite_fraction: f64,
}

#[derive(Serialize)]
struct NonFiniteCounts {
    total: u64,
    finite: u64,
    non_finite: u64,
    nan: u64,
    pos_inf: u64,
    neg_inf: u64,
}

#[derive(Serialize)]
struct HistogramBucket {
    label: String,
    count: u64,
}

#[derive(Serialize)]
struct InvariantResults {
    sig1_non_finite_only_where_input_invalid: bool,
    sig2_calibrated_fraction_not_much_worse_than_power: bool,
    sig3_gains_finite_and_positive: bool,
    notes: Vec<String>,
}

#[derive(Serialize)]
struct BurstSegmentSummary {
    start_line: usize,
    end_line: usize,
    mean_sigma0_db: f64,
    mean_ratio_db: f64,
    zero_fraction: f64,
    nan_fraction: f64,
    finite_fraction: f64,
}

#[derive(Serialize)]
struct BurstBoundaryBandSummary {
    boundary_line: usize,
    pre_mean_sigma0_db: f64,
    post_mean_sigma0_db: f64,
    delta_sigma0_db: f64,
    pre_zero_fraction: f64,
    post_zero_fraction: f64,
    pre_nan_fraction: f64,
    post_nan_fraction: f64,
    pre_mean_ratio_db: f64,
    post_mean_ratio_db: f64,
}

#[derive(Serialize)]
struct BurstDiagnosticsJson {
    central_start_sample: usize,
    central_end_sample: usize,
    detection_threshold_db: f64,
    band_half_height_lines: usize,
    segments: Vec<BurstSegmentSummary>,
    boundaries: Vec<BurstBoundaryBandSummary>,
}

#[derive(Serialize)]
struct SigmaAuditJson {
    product_id: String,
    subswath: String,
    polarization: String,
    calibration_type: String,
    start_time_utc: Option<String>,
    noise_removal_applied: bool,

    // Grid / mask
    lines: usize,
    samples: usize,
    mask_valid_pixels: u64,
    mask_invalid_pixels: u64,

    // Statistics restricted to mask
    power_stats: ScalarStats,
    gain_stats: ScalarStats,
    calibrated_stats: ScalarStats,

    // Non-finite summaries inside mask
    power_nonfinite: NonFiniteCounts,
    gain_nonfinite: NonFiniteCounts,
    calibrated_nonfinite: NonFiniteCounts,

    // Simple histograms (log-scaled buckets)
    power_hist: Vec<HistogramBucket>,
    gain_hist: Vec<HistogramBucket>,
    calibrated_hist: Vec<HistogramBucket>,

    invariants: InvariantResults,
    root_cause_guess: String,
    burst_diagnostics: Option<BurstDiagnosticsJson>,
}

fn parse_bool_env(var: &str) -> bool {
    match env::var(var) {
        Ok(v) => {
            let v = v.to_ascii_lowercase();
            v == "1" || v == "true" || v == "yes" || v == "on"
        }
        Err(_) => false,
    }
}

fn histogram_buckets() -> Vec<(f32, f32, &'static str)> {
    vec![
        (f32::MIN_POSITIVE, 1.0e-10, "(0,1e-10)"),
        (1.0e-10, 1.0e-8, "[1e-10,1e-8)"),
        (1.0e-8, 1.0e-6, "[1e-8,1e-6)"),
        (1.0e-6, 1.0e-4, "[1e-6,1e-4)"),
        (1.0e-4, 1.0e-2, "[1e-4,1e-2)"),
        (1.0e-2, 1.0, "[1e-2,1)"),
        (1.0, 1.0e2, "[1,1e2)"),
        (1.0e2, f32::INFINITY, "[1e2,inf)"),
    ]
}

fn ensure_dir(path: &Path) -> SarResult<()> {
    if let Err(e) = fs::create_dir_all(path) {
        return Err(SarError::Io(e));
    }
    Ok(())
}

/// Run Sigma0 calibration audit and, if configured, write a JSON summary.
///
/// This is intentionally conservative: if `SARDINE_CALIB_AUDIT_DIR` is not set,
/// the function becomes a no-op. Strict failure behaviour is controlled via
/// `SARDINE_CALIB_AUDIT_STRICT`.
pub fn sigma_audit_and_maybe_write_json(
    product_id: &str,
    subswath: &str,
    polarization: &str,
    calibration_type: CalibrationType,
    start_time_utc: Option<&str>,
    noise_removal_applied: bool,
    power_data: &Array2<f32>,
    calibrated_data: &Array2<f32>,
    calibration: &CalibrationCoefficients,
) -> SarResult<()> {
    // Only active for Sigma0
    if calibration_type != CalibrationType::Sigma0 {
        return Ok(());
    }

    let audit_dir = match env::var("SARDINE_CALIB_AUDIT_DIR") {
        Ok(v) if !v.trim().is_empty() => PathBuf::from(v),
        _ => {
            // Audit disabled if no directory configured
            return Ok(());
        }
    };

    let strict = parse_bool_env("SARDINE_CALIB_AUDIT_STRICT");

    let (lines, samples) = power_data.dim();

    let lut = calibration.lut.as_ref().ok_or_else(|| {
        SarError::Processing("Sigma0 calibration audit requires precomputed LUT".to_string())
    })?;

    let effective_ranges: Option<&ValidSampleRanges> = calibration.valid_sample_ranges.as_ref();

    let mut mask_valid_pixels: u64 = 0;

    // Stats initialisation
    let mut power_min = f64::INFINITY;
    let mut power_max = f64::NEG_INFINITY;
    let mut power_sum = 0.0f64;
    let mut power_finite = 0u64;

    let mut gain_min = f64::INFINITY;
    let mut gain_max = f64::NEG_INFINITY;
    let mut gain_sum = 0.0f64;
    let mut gain_finite = 0u64;

    let mut cal_min = f64::INFINITY;
    let mut cal_max = f64::NEG_INFINITY;
    let mut cal_sum = 0.0f64;
    let mut cal_finite = 0u64;

    // Non-finite counters (inside mask)
    let mut power_nan = 0u64;
    let mut power_pos_inf = 0u64;
    let mut power_neg_inf = 0u64;

    let mut gain_nan = 0u64;
    let mut gain_pos_inf = 0u64;
    let mut gain_neg_inf = 0u64;

    let mut cal_nan = 0u64;
    let mut cal_pos_inf = 0u64;
    let mut cal_neg_inf = 0u64;

    // Invariant-related counters
    let mut nonfinite_calib_with_finite_inputs = 0u64;

    let buckets = histogram_buckets();
    let mut power_hist_counts = vec![0u64; buckets.len()];
    let mut gain_hist_counts = vec![0u64; buckets.len()];
    let mut cal_hist_counts = vec![0u64; buckets.len()];

    for row in 0..lines {
        let in_range = effective_ranges.and_then(|ranges| ranges.ranges.get(row));

        for col in 0..samples {
            let in_mask = match in_range {
                Some((start, end_valid)) => col >= *start && col <= *end_valid,
                None => effective_ranges.is_none(),
            };

            if !in_mask {
                continue;
            }
            mask_valid_pixels += 1;

            let p = power_data[[row, col]];
            let g = lut.sigma_values[[row, col]];
            let c = calibrated_data[[row, col]];

            // Power stats
            if p.is_finite() {
                power_finite += 1;
                let v = p as f64;
                power_sum += v;
                if v < power_min {
                    power_min = v;
                }
                if v > power_max {
                    power_max = v;
                }
            } else {
                if p.is_nan() {
                    power_nan += 1;
                } else if p.is_infinite() {
                    if p.is_sign_positive() {
                        power_pos_inf += 1;
                    } else {
                        power_neg_inf += 1;
                    }
                }
            }

            // Gain stats (after LUT & clamping)
            if g.is_finite() {
                let v = g as f64;
                gain_finite += 1;
                gain_sum += v;
                if v < gain_min {
                    gain_min = v;
                }
                if v > gain_max {
                    gain_max = v;
                }
            } else {
                if g.is_nan() {
                    gain_nan += 1;
                } else if g.is_infinite() {
                    if g.is_sign_positive() {
                        gain_pos_inf += 1;
                    } else {
                        gain_neg_inf += 1;
                    }
                }
            }

            // Calibrated stats
            if c.is_finite() {
                cal_finite += 1;
                let v = c as f64;
                cal_sum += v;
                if v < cal_min {
                    cal_min = v;
                }
                if v > cal_max {
                    cal_max = v;
                }
            } else {
                if c.is_nan() {
                    cal_nan += 1;
                } else if c.is_infinite() {
                    if c.is_sign_positive() {
                        cal_pos_inf += 1;
                    } else {
                        cal_neg_inf += 1;
                    }
                }

                // SIG-1 style: calibrated non-finite while both inputs finite
                if p.is_finite() && g.is_finite() {
                    nonfinite_calib_with_finite_inputs += 1;
                }
            }

            // Histograms (positive finite only)
            if p.is_finite() && p > 0.0 {
                for (idx, (lo, hi, _)) in buckets.iter().enumerate() {
                    if p >= *lo && p < *hi {
                        power_hist_counts[idx] += 1;
                        break;
                    }
                }
            }

            if g.is_finite() && g > 0.0 {
                for (idx, (lo, hi, _)) in buckets.iter().enumerate() {
                    if g >= *lo && g < *hi {
                        gain_hist_counts[idx] += 1;
                        break;
                    }
                }
            }

            if c.is_finite() && c > 0.0 {
                for (idx, (lo, hi, _)) in buckets.iter().enumerate() {
                    if c >= *lo && c < *hi {
                        cal_hist_counts[idx] += 1;
                        break;
                    }
                }
            }
        }
    }

    let total_mask_pixels = mask_valid_pixels.max(1); // prevent div-by-zero

    let power_stats = ScalarStats {
        min: if power_min.is_finite() { power_min } else { 0.0 },
        max: if power_max.is_finite() { power_max } else { 0.0 },
        mean: if power_finite > 0 {
            power_sum / (power_finite as f64)
        } else {
            0.0
        },
        finite_fraction: (power_finite as f64) / (total_mask_pixels as f64),
    };

    let gain_stats = ScalarStats {
        min: if gain_min.is_finite() { gain_min } else { 0.0 },
        max: if gain_max.is_finite() { gain_max } else { 0.0 },
        mean: if gain_finite > 0 {
            gain_sum / (gain_finite as f64)
        } else {
            0.0
        },
        finite_fraction: (gain_finite as f64) / (total_mask_pixels as f64),
    };

    let calibrated_stats = ScalarStats {
        min: if cal_min.is_finite() { cal_min } else { 0.0 },
        max: if cal_max.is_finite() { cal_max } else { 0.0 },
        mean: if cal_finite > 0 {
            cal_sum / (cal_finite as f64)
        } else {
            0.0
        },
        finite_fraction: (cal_finite as f64) / (total_mask_pixels as f64),
    };

    let power_nonfinite = NonFiniteCounts {
        total: total_mask_pixels,
        finite: power_finite,
        non_finite: total_mask_pixels - power_finite,
        nan: power_nan,
        pos_inf: power_pos_inf,
        neg_inf: power_neg_inf,
    };

    let gain_nonfinite = NonFiniteCounts {
        total: total_mask_pixels,
        finite: gain_finite,
        non_finite: total_mask_pixels - gain_finite,
        nan: gain_nan,
        pos_inf: gain_pos_inf,
        neg_inf: gain_neg_inf,
    };

    let calibrated_nonfinite = NonFiniteCounts {
        total: total_mask_pixels,
        finite: cal_finite,
        non_finite: total_mask_pixels - cal_finite,
        nan: cal_nan,
        pos_inf: cal_pos_inf,
        neg_inf: cal_neg_inf,
    };

    let power_hist: Vec<HistogramBucket> = buckets
        .iter()
        .zip(power_hist_counts.iter())
        .map(|((_, _, label), count)| HistogramBucket {
            label: (*label).to_string(),
            count: *count,
        })
        .collect();

    let gain_hist: Vec<HistogramBucket> = buckets
        .iter()
        .zip(gain_hist_counts.iter())
        .map(|((_, _, label), count)| HistogramBucket {
            label: (*label).to_string(),
            count: *count,
        })
        .collect();

    let calibrated_hist: Vec<HistogramBucket> = buckets
        .iter()
        .zip(cal_hist_counts.iter())
        .map(|((_, _, label), count)| HistogramBucket {
            label: (*label).to_string(),
            count: *count,
        })
        .collect();

    // Invariants (heuristic, tuned for debugging rather than hard science here)
    let mut notes = Vec::new();

    let sig1_ok = nonfinite_calib_with_finite_inputs == 0;
    if !sig1_ok {
        notes.push(format!(
            "SIG-1: {} calibrated pixels non-finite despite finite inputs",
            nonfinite_calib_with_finite_inputs
        ));
    }

    // Allow up to 10 percentage points drop from power to calibrated finite fraction
    let max_allowed_drop = 0.10;
    let finite_power_frac = power_stats.finite_fraction;
    let finite_calib_frac = calibrated_stats.finite_fraction;
    let drop = (finite_power_frac - finite_calib_frac).max(0.0);
    let sig2_ok = drop <= max_allowed_drop;
    if !sig2_ok {
        notes.push(format!(
            "SIG-2: finite fraction drop {:.3} exceeds allowed {:.3}",
            drop, max_allowed_drop
        ));
    }

    // Require essentially all gains finite & positive inside mask
    let sig3_ok = gain_nonfinite.non_finite == 0;
    if !sig3_ok {
        notes.push(format!(
            "SIG-3: {} non-finite gain entries inside mask",
            gain_nonfinite.non_finite
        ));
    }

    let invariants = InvariantResults {
        sig1_non_finite_only_where_input_invalid: sig1_ok,
        sig2_calibrated_fraction_not_much_worse_than_power: sig2_ok,
        sig3_gains_finite_and_positive: sig3_ok,
        notes,
    };

    let root_cause_guess = if !sig1_ok {
        "calibration_lut_or_math".to_string()
    } else if !sig2_ok {
        "mask_mismatch_or_excessive_clamping".to_string()
    } else if !sig3_ok {
        "lut_coverage_or_parsing_issue".to_string()
    } else {
        "ok".to_string()
    };

    let burst_diagnostics =
        compute_burst_diagnostics(power_data, calibrated_data, effective_ranges, lines, samples);

    let audit_json = SigmaAuditJson {
        product_id: product_id.to_string(),
        subswath: subswath.to_string(),
        polarization: polarization.to_string(),
        calibration_type: match calibration_type {
            CalibrationType::Sigma0 => "sigma0".to_string(),
            CalibrationType::Beta0 => "beta0".to_string(),
            CalibrationType::Gamma0 => "gamma0".to_string(),
            CalibrationType::Dn => "dn".to_string(),
        },
        start_time_utc: start_time_utc.map(|s| s.to_string()),
        noise_removal_applied,
        lines,
        samples,
        mask_valid_pixels,
        mask_invalid_pixels: (lines as u64 * samples as u64).saturating_sub(mask_valid_pixels),
        power_stats,
        gain_stats,
        calibrated_stats,
        power_nonfinite,
        gain_nonfinite,
        calibrated_nonfinite,
        power_hist,
        gain_hist,
        calibrated_hist,
        invariants,
        root_cause_guess: root_cause_guess.clone(),
        burst_diagnostics,
    };

    ensure_dir(&audit_dir)?;

    let fname = format!(
        "sigma_audit_{}_{}_{}.json",
        subswath.to_lowercase(),
        polarization.to_lowercase(),
        product_id
    );
    let out_path = audit_dir.join(fname);

    let json_str = serde_json::to_string_pretty(&audit_json).map_err(|e| {
        SarError::Processing(format!("Failed to serialize sigma audit JSON: {}", e))
    })?;

    if let Err(e) = fs::write(&out_path, json_str) {
        return Err(SarError::Io(e));
    }

    log::info!(
        "[SIG-AUDIT] {} {}: mask_valid={} finite_power_frac={:.3} finite_calib_frac={:.3} root_cause={}",
        subswath,
        polarization,
        mask_valid_pixels,
        finite_power_frac,
        finite_calib_frac,
        root_cause_guess,
    );

    // If strict mode is enabled and any invariant fails, raise a hard error.
    if strict {
        if !audit_json
            .invariants
            .sig1_non_finite_only_where_input_invalid
            || !audit_json
                .invariants
                .sig2_calibrated_fraction_not_much_worse_than_power
            || !audit_json.invariants.sig3_gains_finite_and_positive
        {
            return Err(SarError::Processing(
                "Sigma0 calibration audit invariants failed (strict mode)".to_string(),
            ));
        }
    }

    Ok(())
}

/// Apply a conservative, multiplicative equalization across burst seams for
/// Sigma0, driven by the same logic used for burst diagnostics. This mutates
/// the calibrated Sigma0 array in-place while preserving NaN/zero semantics.
///
/// The algorithm is intentionally cautious:
/// - Only consider seams where both sides have high finite fraction and low
///   zero_fraction (to avoid coverage ramps and noise-only regions).
/// - Require |Δσ0| >= MIN_JUMP_DB before touching any seam.
/// - Enforce that σ0/|DN|^2 is essentially unchanged across the seam.
/// - Cluster nearby boundaries into a single physical seam and cap the
///   correction magnitude per seam.
pub(crate) fn equalize_sigma0_burst_seams_inplace(
    power_data: &Array2<f32>,
    calibrated_data: &mut Array2<f32>,
    effective_ranges: Option<&ValidSampleRanges>,
) {
    let (lines, samples) = calibrated_data.dim();
    if lines < 4 || samples < 4 {
        return;
    }

    let diagnostics = match compute_burst_diagnostics(
        power_data,
        calibrated_data,
        effective_ranges,
        lines,
        samples,
    ) {
        Some(d) => d,
        None => return,
    };

    if diagnostics.boundaries.is_empty() {
        return;
    }

    // Heuristic thresholds – tuned to only touch clear ≈1 dB stripes while
    // leaving smaller, acceptable seams untouched.
    const MIN_JUMP_DB: f64 = 0.7; // require at least this jump to touch seam
    const MAX_CORRECTION_DB: f64 = 0.8; // clamp per-seam correction magnitude
    const MIN_FINITE_FRACTION: f64 = 0.8; // both sides must be well-populated
    const MAX_ZERO_FRACTION: f64 = 0.2; // exclude noise-floor/coverage ramps
    const MAX_RATIO_DELTA_DB: f64 = 0.02; // σ0/|DN|^2 must be essentially flat

    struct SeamCluster<'a> {
        start_line: usize,
        boundaries: Vec<&'a BurstBoundaryBandSummary>,
    }

    let mut clusters: Vec<SeamCluster<'_>> = Vec::new();

    for b in &diagnostics.boundaries {
        if !b.delta_sigma0_db.is_finite() {
            continue;
        }

        // Approximate finite fraction from zeros + NaNs.
        let pre_finite = 1.0 - b.pre_zero_fraction - b.pre_nan_fraction;
        let post_finite = 1.0 - b.post_zero_fraction - b.post_nan_fraction;

        if pre_finite < MIN_FINITE_FRACTION || post_finite < MIN_FINITE_FRACTION {
            continue;
        }
        if b.pre_zero_fraction > MAX_ZERO_FRACTION || b.post_zero_fraction > MAX_ZERO_FRACTION {
            continue;
        }

        let delta = b.delta_sigma0_db;
        if delta.abs() < MIN_JUMP_DB {
            continue;
        }

        let ratio_delta = (b.pre_mean_ratio_db - b.post_mean_ratio_db).abs();
        if ratio_delta > MAX_RATIO_DELTA_DB {
            continue;
        }

        if let Some(last) = clusters.last_mut() {
            if let Some(prev) = last.boundaries.last() {
                if b.boundary_line
                    <= prev.boundary_line.saturating_add(diagnostics.band_half_height_lines)
                {
                    last.boundaries.push(b);
                    continue;
                }
            }
        }

        clusters.push(SeamCluster {
            start_line: b.boundary_line,
            boundaries: vec![b],
        });
    }

    if clusters.is_empty() {
        return;
    }

    let n_clusters = clusters.len();
    for idx in 0..n_clusters {
        let start_row = clusters[idx].start_line.min(lines.saturating_sub(1));
        let end_row = if idx + 1 < n_clusters {
            let next_start = clusters[idx + 1].start_line;
            next_start.saturating_sub(1).min(lines.saturating_sub(1))
        } else {
            lines.saturating_sub(1)
        };

        if end_row <= start_row {
            continue;
        }

        let mut deltas: Vec<f64> = clusters[idx]
            .boundaries
            .iter()
            .map(|b| b.delta_sigma0_db)
            .filter(|d| d.is_finite())
            .collect();
        if deltas.is_empty() {
            continue;
        }
        deltas.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
        let mid = deltas.len() / 2;
        let mut delta_db = if deltas.len() % 2 == 1 {
            deltas[mid]
        } else {
            0.5 * (deltas[mid - 1] + deltas[mid])
        };

        if !delta_db.is_finite() || delta_db.abs() < MIN_JUMP_DB {
            continue;
        }

        if delta_db.abs() > MAX_CORRECTION_DB {
            delta_db = if delta_db.is_sign_positive() {
                MAX_CORRECTION_DB
            } else {
                -MAX_CORRECTION_DB
            };
        }

        let factor: f32 = 10f32.powf((-delta_db as f32) / 10.0);

        log::info!(
            "[SIG-AUDIT|BURST-EQ] Applying {:.3} dB correction (factor {:.6}) to rows {}..={} (Sigma0)",
            delta_db,
            factor,
            start_row,
            end_row
        );

        for row in start_row..=end_row {
            if row >= lines {
                break;
            }
            let mut row_view = calibrated_data.row_mut(row);
            for v in row_view.iter_mut() {
                if v.is_finite() && *v > 0.0 {
                    *v *= factor;
                }
            }
        }
    }
}

fn compute_burst_diagnostics(
    power_data: &Array2<f32>,
    calibrated_data: &Array2<f32>,
    effective_ranges: Option<&ValidSampleRanges>,
    lines: usize,
    samples: usize,
) -> Option<BurstDiagnosticsJson> {
    if lines < 4 || samples < 4 {
        return None;
    }

    // Determine a global central-range window in range (columns)
    let (global_start, global_end) = if let Some(ranges) = effective_ranges {
        let mut min_start = samples;
        let mut max_end = 0usize;
        for (row, (start, end)) in ranges.ranges.iter().enumerate() {
            if row >= lines {
                break;
            }
            min_start = min_start.min(*start);
            max_end = max_end.max(*end);
        }
        if min_start >= samples || max_end <= min_start {
            (0usize, samples.saturating_sub(1))
        } else {
            (min_start, max_end)
        }
    } else {
        (0usize, samples.saturating_sub(1))
    };

    if global_end <= global_start {
        return None;
    }

    let width = global_end - global_start + 1;
    let central_start = global_start + width / 5; // ~20% in from left
    let central_end = global_start + (4 * width) / 5; // ~80% in from left

    if central_end <= central_start {
        return None;
    }

    let detection_threshold_db = 0.2_f64;
    let band_half_height_lines: usize = 16;

    // Per-row aggregates over the central window
    let mut row_total: Vec<u64> = vec![0; lines];
    let mut row_zero: Vec<u64> = vec![0; lines];
    let mut row_nan: Vec<u64> = vec![0; lines];
    let mut row_finite: Vec<u64> = vec![0; lines];
    let mut row_sigma_sum_db: Vec<f64> = vec![0.0; lines];
    let mut row_ratio_sum_db: Vec<f64> = vec![0.0; lines];

    for row in 0..lines {
        let range = effective_ranges.and_then(|ranges| ranges.ranges.get(row));

        for col in central_start..=central_end {
            if let Some((start, end_valid)) = range {
                if col < *start || col > *end_valid {
                    continue;
                }
            }

            let p = power_data[[row, col]];
            let c = calibrated_data[[row, col]];
            row_total[row] += 1;

            if !p.is_finite() || p <= MIN_VALID_POWER {
                row_zero[row] += 1;
                continue;
            }

            if !c.is_finite() || c <= 0.0 {
                if c.is_nan() {
                    row_nan[row] += 1;
                }
                continue;
            }

            row_finite[row] += 1;
            let sigma_db = 10.0 * (c as f64).log10();
            let ratio = (c / p).max(f32::MIN_POSITIVE);
            let ratio_db = 10.0 * (ratio as f64).log10();
            row_sigma_sum_db[row] += sigma_db;
            row_ratio_sum_db[row] += ratio_db;
        }
    }

    let mut row_mean_sigma_db: Vec<f64> = vec![f64::NAN; lines];
    let mut row_mean_ratio_db: Vec<f64> = vec![f64::NAN; lines];
    for row in 0..lines {
        if row_finite[row] > 0 {
            let denom = row_finite[row] as f64;
            row_mean_sigma_db[row] = row_sigma_sum_db[row] / denom;
            row_mean_ratio_db[row] = row_ratio_sum_db[row] / denom;
        }
    }

    // Detect candidate burst boundaries where mean sigma0 dB jumps
    let min_valid_per_row: u64 = 16;
    let mut boundary_rows: Vec<usize> = Vec::new();
    for row in 1..lines {
        if row_finite[row] < min_valid_per_row || row_finite[row - 1] < min_valid_per_row {
            continue;
        }
        let prev = row_mean_sigma_db[row - 1];
        let curr = row_mean_sigma_db[row];
        if !prev.is_finite() || !curr.is_finite() {
            continue;
        }
        let diff = curr - prev;
        if diff.abs() >= detection_threshold_db {
            // Avoid duplicating very close boundaries
            if boundary_rows.last().copied().map_or(true, |last| row > last + 1) {
                boundary_rows.push(row);
            }
        }
    }

    // Build segment summaries between boundaries
    let mut segments: Vec<BurstSegmentSummary> = Vec::new();
    let mut seg_start: usize = 0;
    let add_segment = |start: usize, end: usize,
                           segments: &mut Vec<BurstSegmentSummary>| {
        if end < start || end >= lines {
            return;
        }
        let mut total: u64 = 0;
        let mut zeros: u64 = 0;
        let mut nans: u64 = 0;
        let mut finite: u64 = 0;
        let mut sigma_sum: f64 = 0.0;
        let mut ratio_sum: f64 = 0.0;
        for row in start..=end {
            total += row_total[row];
            zeros += row_zero[row];
            nans += row_nan[row];
            finite += row_finite[row];
            sigma_sum += row_sigma_sum_db[row];
            ratio_sum += row_ratio_sum_db[row];
        }
        if total == 0 {
            return;
        }
        let finite_f = finite as f64;
        let mean_sigma0_db = if finite > 0 {
            sigma_sum / finite_f
        } else {
            f64::NAN
        };
        let mean_ratio_db = if finite > 0 {
            ratio_sum / finite_f
        } else {
            f64::NAN
        };
        let total_f = total as f64;
        let zero_fraction = zeros as f64 / total_f;
        let nan_fraction = nans as f64 / total_f;
        let finite_fraction = finite as f64 / total_f;
        segments.push(BurstSegmentSummary {
            start_line: start,
            end_line: end,
            mean_sigma0_db,
            mean_ratio_db,
            zero_fraction,
            nan_fraction,
            finite_fraction,
        });
    };

    for &b in &boundary_rows {
        let end = b.saturating_sub(1);
        add_segment(seg_start, end, &mut segments);
        seg_start = b;
    }
    if seg_start < lines {
        add_segment(seg_start, lines - 1, &mut segments);
    }

    // Build boundary-band diagnostics around each detected boundary
    let mut boundaries: Vec<BurstBoundaryBandSummary> = Vec::new();
    for &b in &boundary_rows {
        if b == 0 || b >= lines {
            continue;
        }
        let pre_start = b.saturating_sub(band_half_height_lines).min(lines.saturating_sub(1));
        let pre_end = b.saturating_sub(1);
        let post_start = b;
        let post_end = (b + band_half_height_lines - 1).min(lines.saturating_sub(1));

        let summarize_band = |start: usize, end: usize| -> (f64, f64, f64, f64, f64) {
            if end < start {
                return (f64::NAN, 0.0, 0.0, f64::NAN, f64::NAN);
            }
            let mut total: u64 = 0;
            let mut zeros: u64 = 0;
            let mut nans: u64 = 0;
            let mut finite: u64 = 0;
            let mut sigma_sum: f64 = 0.0;
            let mut ratio_sum: f64 = 0.0;
            for row in start..=end {
                total += row_total[row];
                zeros += row_zero[row];
                nans += row_nan[row];
                finite += row_finite[row];
                sigma_sum += row_sigma_sum_db[row];
                ratio_sum += row_ratio_sum_db[row];
            }
            if total == 0 {
                return (f64::NAN, 0.0, 0.0, f64::NAN, f64::NAN);
            }
            let finite_f = finite as f64;
            let mean_sigma0_db = if finite > 0 {
                sigma_sum / finite_f
            } else {
                f64::NAN
            };
            let mean_ratio_db = if finite > 0 {
                ratio_sum / finite_f
            } else {
                f64::NAN
            };
            let total_f = total as f64;
            let zero_fraction = zeros as f64 / total_f;
            let nan_fraction = nans as f64 / total_f;
            (mean_sigma0_db, mean_ratio_db, zero_fraction, nan_fraction, finite_f / total_f)
        };

        let (pre_mean_sigma0_db, pre_mean_ratio_db, pre_zero_fraction, pre_nan_fraction, _pre_finite_frac) =
            summarize_band(pre_start, pre_end);
        let (post_mean_sigma0_db, post_mean_ratio_db, post_zero_fraction, post_nan_fraction, _post_finite_frac) =
            summarize_band(post_start, post_end);

        let delta_sigma0_db = if pre_mean_sigma0_db.is_finite() && post_mean_sigma0_db.is_finite()
        {
            post_mean_sigma0_db - pre_mean_sigma0_db
        } else {
            f64::NAN
        };

        boundaries.push(BurstBoundaryBandSummary {
            boundary_line: b,
            pre_mean_sigma0_db,
            post_mean_sigma0_db,
            delta_sigma0_db,
            pre_zero_fraction,
            post_zero_fraction,
            pre_nan_fraction,
            post_nan_fraction,
            pre_mean_ratio_db,
            post_mean_ratio_db,
        });
    }

    Some(BurstDiagnosticsJson {
        central_start_sample: central_start,
        central_end_sample: central_end,
        detection_threshold_db,
        band_half_height_lines,
        segments,
        boundaries,
    })
}

