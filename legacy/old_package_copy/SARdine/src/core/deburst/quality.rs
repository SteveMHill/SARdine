#![allow(dead_code)]
use crate::core::deburst::burst_ops::valid_window;
use crate::core::deburst::geometry::BurstInfo;
use crate::types::SarComplex;
use ndarray::Array2;

/// Compute total power of a complex array.
pub(crate) fn calculate_total_power(data: &Array2<SarComplex>) -> f64 {
    data.iter()
        .map(|&sample| (sample.re as f64).powi(2) + (sample.im as f64).powi(2))
        .sum()
}

/// Per-burst power diagnostics for phase tracking.
pub(crate) fn calculate_burst_power_diagnostics(
    burst_info: &[BurstInfo],
    slc_data: &Array2<SarComplex>,
) -> Vec<(usize, f64, f64)> {
    let mut diagnostics = Vec::new();

    for (burst_idx, burst) in burst_info.iter().enumerate() {
        let burst_lines = burst.lines();
        let burst_samples = burst.end_sample.saturating_sub(burst.start_sample) + 1;

        let mut burst_power = 0.0_f64;
        let mut valid_pixels = 0_usize;

        for line_in_burst in 0..burst_lines {
            let src_line = burst.start_line + line_in_burst;
            if src_line >= slc_data.nrows() {
                continue;
            }

            let (valid_start, valid_end) = crate::core::deburst::burst_ops::valid_window(
                line_in_burst,
                &burst.first_valid_sample,
                &burst.last_valid_sample,
                burst_samples,
            );

            for col in valid_start..valid_end {
                let src_col = burst.start_sample + col;
                if src_col >= slc_data.ncols() {
                    continue;
                }

                let sample = slc_data[[src_line, src_col]];
                burst_power += (sample.re as f64).powi(2) + (sample.im as f64).powi(2);
                valid_pixels += 1;
            }
        }

        let mean_power = if valid_pixels > 0 {
            burst_power / valid_pixels as f64
        } else {
            0.0
        };

        diagnostics.push((burst_idx, burst_power, mean_power));
    }

    diagnostics
}

/// Log diagnostics for debugging deramp issues.
pub(crate) fn log_power_diagnostics(diagnostics: &[(usize, f64, f64)]) {
    if diagnostics.is_empty() {
        return;
    }

    let total_power: f64 = diagnostics.iter().map(|(_, p, _)| p).sum();
    let mean_powers: Vec<f64> = diagnostics.iter().map(|(_, _, mp)| *mp).collect();
    let mean_of_means = mean_powers.iter().sum::<f64>() / mean_powers.len() as f64;

    let variance: f64 = mean_powers
        .iter()
        .map(|mp| (mp - mean_of_means).powi(2))
        .sum::<f64>()
        / mean_powers.len() as f64;
    let std_dev = variance.sqrt();
    let coeff_variation = if mean_of_means > 0.0 {
        std_dev / mean_of_means
    } else {
        0.0
    };

    log::info!("📊 Enhancement #4: Per-burst power diagnostics:");
    log::info!("   Total power across all bursts: {:.3e}", total_power);
    log::info!("   Mean power per pixel (average): {:.3e}", mean_of_means);
    log::info!(
        "   Power variation across bursts (CV): {:.3}%",
        coeff_variation * 100.0
    );

    for (burst_idx, burst_power, mean_power) in diagnostics {
        let deviation = ((mean_power - mean_of_means) / mean_of_means * 100.0).abs();
        let status = if deviation < 5.0 {
            "✅"
        } else if deviation < 15.0 {
            "⚠️ "
        } else {
            "❌"
        };

        log::info!(
            "   {} Burst {}: power={:.3e}, mean={:.3e}, deviation={:.1}%",
            status,
            burst_idx + 1,
            burst_power,
            mean_power,
            deviation
        );
    }

    log::info!(
        "   Power CV summary: {:.2}% (threshold 15%)",
        coeff_variation * 100.0
    );

    if coeff_variation > 0.15 {
        log::warn!("⚠️  High power variation across bursts (CV > 15%)");
        log::warn!("   Possible causes:");
        log::warn!("   - Incorrect DC/FM polynomial reference time");
        log::warn!("   - Missing range-dependent deramp");
        log::warn!("   - Invalid deramp parameters from annotation");
    }
}

/// Evaluate polynomial using Horner's rule for efficiency
/// OPTIMIZATION #32: Avoids repeated t.powi() calls - O(n) multiplications instead of O(n²)
fn eval_poly(poly: &[f64], t: f64) -> f64 {
    // Horner's method: c0 + c1*t + c2*t^2 = c0 + t*(c1 + t*c2)
    poly.iter().rev().fold(0.0, |acc, &c| acc * t + c)
}

/// Check Doppler centroid and FM consistency at burst overlaps.
/// Flags bursts whose DC/FM coefficients differ beyond small thresholds.
///
/// NOTE (Jan 2026): DC/FM polynomials are RANGE polynomials, not azimuth polynomials!
/// The independent variable is slant range time (dt = range_time - t0), NOT azimuth time.
/// We compare the polynomial coefficients directly since they should be similar for
/// temporally adjacent bursts if per-burst selection is working correctly.
pub(crate) fn log_doppler_fm_consistency(burst_info: &[BurstInfo]) {
    if burst_info.len() < 2 {
        return;
    }

    // Compare DC c0 coefficients directly (constant term in Hz)
    // For properly selected per-burst polynomials, adjacent bursts should have similar c0
    const DC_C0_THRESH_HZ: f64 = 50.0; // Allow some variation for scene Doppler changes
    const FM_C0_THRESH_HZ_S: f64 = 100.0;

    for w in burst_info.windows(2) {
        let a = &w[0];
        let b = &w[1];

        // Compare c0 coefficients directly (no polynomial evaluation needed)
        let dc_c0_a = a.dc_polynomial.first().copied().unwrap_or(0.0);
        let dc_c0_b = b.dc_polynomial.first().copied().unwrap_or(0.0);
        let fm_c0_a = a.fm_polynomial.first().copied().unwrap_or(0.0);
        let fm_c0_b = b.fm_polynomial.first().copied().unwrap_or(0.0);

        let dc_diff = dc_c0_a - dc_c0_b;
        let fm_diff = fm_c0_a - fm_c0_b;

        let status_dc = if dc_diff.abs() > DC_C0_THRESH_HZ {
            "⚠️"
        } else {
            "✅"
        };
        let status_fm = if fm_diff.abs() > FM_C0_THRESH_HZ_S {
            "⚠️"
        } else {
            "✅"
        };

        // Check if bursts have identical coefficients (indicates subswath-wide fallback)
        let coeffs_identical = a.dc_polynomial == b.dc_polynomial;
        let identity_warning = if coeffs_identical {
            " [IDENTICAL - check per-burst selection]"
        } else {
            ""
        };

        log::info!(
            "{} Burst {}→{} DC c0 diff: {:.3} Hz (a={:.3}, b={:.3}){}",
            status_dc,
            a.burst_id,
            b.burst_id,
            dc_diff,
            dc_c0_a,
            dc_c0_b,
            identity_warning
        );
        log::info!(
            "{} Burst {}→{} FM c0 diff: {:.3} Hz/s",
            status_fm,
            a.burst_id,
            b.burst_id,
            fm_diff
        );

        if dc_diff.abs() > DC_C0_THRESH_HZ || fm_diff.abs() > FM_C0_THRESH_HZ_S {
            log::warn!(
                "⚠️  Doppler/FM coefficient discontinuity ({}→{}): Δdc_c0={:.2} Hz, Δfm_c0={:.2} Hz/s",
                a.burst_id,
                b.burst_id,
                dc_diff,
                fm_diff
            );
        }
    }
}

fn line_mean_intensity(
    burst: &BurstInfo,
    line_in_burst: usize,
    slc_data: &Array2<SarComplex>,
) -> Option<f64> {
    let burst_width = burst.end_sample.saturating_sub(burst.start_sample) + 1;
    let (valid_start, valid_end) = valid_window(
        line_in_burst,
        &burst.first_valid_sample,
        &burst.last_valid_sample,
        burst_width,
    );
    if valid_end <= valid_start {
        return None;
    }

    let src_line = burst.start_line + line_in_burst;
    if src_line >= slc_data.nrows() {
        return None;
    }

    let mut sum = 0.0;
    let mut count = 0usize;
    for col in valid_start..valid_end {
        let src_col = burst.start_sample + col;
        if src_col >= slc_data.ncols() {
            continue;
        }
        let s = slc_data[[src_line, src_col]];
        sum += (s.re as f64).powi(2) + (s.im as f64).powi(2);
        count += 1;
    }

    if count == 0 {
        None
    } else {
        Some(sum / count as f64)
    }
}

fn best_shift(a: &[f64], b: &[f64], max_shift: i32) -> Option<f32> {
    if a.is_empty() || b.is_empty() {
        return None;
    }

    let mut best_corr = f64::NEG_INFINITY;
    let mut best_shift = 0i32;

    for s in -max_shift..=max_shift {
        let mut num = 0.0;
        let mut na = 0;
        for (i, &va) in a.iter().enumerate() {
            let j = i as i32 + s;
            if j < 0 || j as usize >= b.len() {
                continue;
            }
            num += va * b[j as usize];
            na += 1;
        }
        if na == 0 {
            continue;
        }
        let corr = num / na as f64;
        if corr > best_corr {
            best_corr = corr;
            best_shift = s;
        }
    }

    // Optional sub-pixel refinement via quadratic fit around peak
    if best_shift > -max_shift && best_shift < max_shift {
        let s0 = best_shift;
        let f_m1 = correlation_at_shift(a, b, s0 - 1);
        let f_0 = correlation_at_shift(a, b, s0);
        let f_p1 = correlation_at_shift(a, b, s0 + 1);
        let denom = f_m1 - 2.0 * f_0 + f_p1;
        if denom.abs() > 1e-6 {
            let delta = 0.5 * (f_m1 - f_p1) / denom;
            let refined = s0 as f64 + delta;
            return Some(refined as f32);
        }
    }

    Some(best_shift as f32)
}

fn correlation_at_shift(a: &[f64], b: &[f64], shift: i32) -> f64 {
    let mut num = 0.0;
    let mut na = 0;
    for (i, &va) in a.iter().enumerate() {
        let j = i as i32 + shift;
        if j < 0 || j as usize >= b.len() {
            continue;
        }
        num += va * b[j as usize];
        na += 1;
    }
    if na == 0 {
        0.0
    } else {
        num / na as f64
    }
}

/// Estimate small azimuth shifts between adjacent bursts using overlap means (ESD-lite).
/// Returns per-burst shift in lines (burst 0 fixed at 0). Positive shift delays burst k+1.
pub(crate) fn estimate_overlap_shifts(
    burst_info: &[BurstInfo],
    slc_data: &Array2<SarComplex>,
    blend_lines: usize,
) -> Vec<f32> {
    let mut shifts = vec![0.0f32; burst_info.len()];
    if burst_info.len() < 2 || blend_lines == 0 {
        return shifts;
    }

    let max_shift = 2; // search ±2 lines

    for i in 0..burst_info.len().saturating_sub(1) {
        let a = &burst_info[i];
        let b = &burst_info[i + 1];
        let overlap = blend_lines.min(a.lines()).min(b.lines());
        if overlap < 8 {
            continue;
        }

        let mut seq_a = Vec::with_capacity(overlap);
        let mut seq_b = Vec::with_capacity(overlap);
        for k in 0..overlap {
            let idx_a = a.lines().saturating_sub(overlap) + k;
            let idx_b = k;
            if let (Some(ma), Some(mb)) = (
                line_mean_intensity(a, idx_a, slc_data),
                line_mean_intensity(b, idx_b, slc_data),
            ) {
                seq_a.push(ma);
                seq_b.push(mb);
            }
        }

        if seq_a.len() < 4 || seq_b.len() < 4 {
            continue;
        }

        if let Some(shift) = best_shift(&seq_a, &seq_b, max_shift) {
            shifts[i + 1] = shift;
            log::info!(
                "ESD-lite overlap shift Burst {}→{}: {:.3} lines (window {})",
                a.burst_id,
                b.burst_id,
                shift,
                overlap
            );
        }
    }

    shifts
}
