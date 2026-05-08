//! STEP-2 diagnostics for IW TOPSAR deburst.
//!
//! This module computes coverage, contribution histograms, overlap radiometry,
//! DC/FM bracketing sanity, deramp effectiveness, and calibration LUT integrity
//! metrics. It is gated by `SARDINE_STEP2_DIAGNOSTICS` to avoid runtime cost
//! unless explicitly enabled.

use crate::core::deburst::burst_ops::valid_window;
use crate::core::deburst::iw_deburst::types::BurstInfo;
use crate::core::deburst::iw_deburst::types::Step2Diagnostics;
use crate::core::deburst::iw_deburst::DerampPhaseStats;
use crate::types::SarComplex;
use ndarray::{s, Array2};
use serde::Serialize;

#[derive(Debug, Clone, Serialize)]
pub(crate) struct OverlapMetric {
    pub burst_a: usize,
    pub burst_b: usize,
    pub mean_power_a: f64,
    pub mean_power_b: f64,
    pub ratio: f64,
    pub delta_db: f64,
}

#[derive(Debug, Clone, Serialize)]
pub(crate) struct DcBracketLog {
    pub burst_id: usize,
    pub target_time: f64,
    pub lower_time: f64,
    pub upper_time: f64,
    pub weight: f64,
    pub extrapolated: bool,
    pub dc_hz: f64,
}

/// Collect coverage and contribution histogram from hit-count mask.
pub(crate) fn coverage_stats(hit_count: &Array2<u16>) -> (usize, usize, f64, [usize; 3]) {
    let total = hit_count.len();
    let mut covered = 0usize;
    let mut hist = [0usize; 3]; // 0,1,2+

    for &v in hit_count.iter() {
        match v {
            0 => hist[0] += 1,
            1 => {
                hist[1] += 1;
                covered += 1;
            }
            _ => {
                hist[2] += 1;
                covered += 1;
            }
        }
    }

    let uncovered = total.saturating_sub(covered);
    let frac = if total == 0 { 0.0 } else { covered as f64 / total as f64 };
    (covered, uncovered, frac, hist)
}

/// Compute overlap radiometry using pre-blend SLC data.
pub(crate) fn overlap_metrics(
    burst_info: &[BurstInfo],
    slc: &Array2<SarComplex>,
    blend_lines: usize,
    sample_margin_frac: f32,
) -> Vec<OverlapMetric> {
    let mut out = Vec::new();
    if burst_info.len() < 2 || blend_lines == 0 {
        return out;
    }

    for w in burst_info.windows(2) {
        let a = &w[0];
        let b = &w[1];

        let overlap = blend_lines.min(a.lines()).min(b.lines());
        if overlap == 0 {
            continue;
        }

        // Restrict to valid window and crop margins to avoid edge noise.
        let burst_width = |burst: &BurstInfo| burst.end_sample.saturating_sub(burst.start_sample) + 1;
        let width_a = burst_width(a);
        let width_b = burst_width(b);

        let margin_a = (width_a as f32 * sample_margin_frac).round() as usize;
        let margin_b = (width_b as f32 * sample_margin_frac).round() as usize;

        let mut sum_a = 0.0;
        let mut sum_b = 0.0;
        let mut n_a = 0usize;
        let mut n_b = 0usize;

        for k in 0..overlap {
            let line_a = a.lines().saturating_sub(overlap) + k;
            let line_b = k;

            let (va0, va1) = valid_window(line_a, &a.first_valid_sample, &a.last_valid_sample, width_a);
            let (vb0, vb1) = valid_window(line_b, &b.first_valid_sample, &b.last_valid_sample, width_b);

            let va0 = va0.saturating_add(margin_a).min(width_a);
            let va1 = va1.saturating_sub(margin_a).max(va0);
            let vb0 = vb0.saturating_add(margin_b).min(width_b);
            let vb1 = vb1.saturating_sub(margin_b).max(vb0);

            if va1 > va0 {
                let src_line = a.start_line + line_a;
                if src_line < slc.nrows() {
                    let row = slc.slice(s![src_line, va0 + a.start_sample..va1 + a.start_sample]);
                    for s in row.iter() {
                        sum_a += (s.re as f64).powi(2) + (s.im as f64).powi(2);
                    }
                    n_a += row.len();
                }
            }

            if vb1 > vb0 {
                let src_line = b.start_line + line_b;
                if src_line < slc.nrows() {
                    let row = slc.slice(s![src_line, vb0 + b.start_sample..vb1 + b.start_sample]);
                    for s in row.iter() {
                        sum_b += (s.re as f64).powi(2) + (s.im as f64).powi(2);
                    }
                    n_b += row.len();
                }
            }
        }

        let mean_a = if n_a > 0 { sum_a / n_a as f64 } else { 0.0 };
        let mean_b = if n_b > 0 { sum_b / n_b as f64 } else { 0.0 };
        let ratio = if mean_b > 0.0 { mean_a / mean_b } else { 0.0 };
        let delta_db = if mean_a > 0.0 && mean_b > 0.0 {
            10.0 * (mean_a / mean_b).log10()
        } else {
            0.0
        };

        out.push(OverlapMetric {
            burst_a: a.burst_id,
            burst_b: b.burst_id,
            mean_power_a: mean_a,
            mean_power_b: mean_b,
            ratio,
            delta_db,
        });
    }

    out
}

/// Compute masked per-burst mean power for CV estimation.
pub(crate) fn per_burst_means(
    burst_info: &[BurstInfo],
    slc: &Array2<SarComplex>,
    sample_margin_frac: f32,
    line_margin_frac: f32,
) -> Vec<f64> {
    let mut means = Vec::new();

    for burst in burst_info {
        let lines = burst.lines();
        let width = burst.end_sample.saturating_sub(burst.start_sample) + 1;
        if lines == 0 || width == 0 {
            means.push(0.0);
            continue;
        }

        let line_margin = ((lines as f32) * line_margin_frac).round() as usize;
        let sample_margin = ((width as f32) * sample_margin_frac).round() as usize;

        let mut sum = 0.0;
        let mut n = 0usize;

        for line_in_burst in line_margin..lines.saturating_sub(line_margin) {
            let (v0, v1) = valid_window(line_in_burst, &burst.first_valid_sample, &burst.last_valid_sample, width);
            let v0 = v0.saturating_add(sample_margin).min(width);
            let v1 = v1.saturating_sub(sample_margin).max(v0);
            if v1 <= v0 {
                continue;
            }
            let src_line = burst.start_line + line_in_burst;
            if src_line >= slc.nrows() {
                continue;
            }
            let row = slc.slice(s![src_line, v0 + burst.start_sample..v1 + burst.start_sample]);
            for s in row.iter() {
                sum += (s.re as f64).powi(2) + (s.im as f64).powi(2);
            }
            n += row.len();
        }

        let mean = if n > 0 { sum / n as f64 } else { 0.0 };
        means.push(mean);
    }

    means
}

pub(crate) fn coeff_of_variation(means: &[f64]) -> f64 {
    if means.is_empty() {
        return 0.0;
    }
    let m = means.iter().sum::<f64>() / means.len() as f64;
    if m == 0.0 {
        return 0.0;
    }
    let var = means.iter().map(|x| (x - m).powi(2)).sum::<f64>() / means.len() as f64;
    (var.sqrt() / m).abs()
}

/// Bundle diagnostics and optionally serialize.
pub(crate) fn build_step2_diagnostics(
    burst_info: &[BurstInfo],
    slc: &Array2<SarComplex>,
    hit_count: &Array2<u16>,
    deramp_stats: Option<&DerampPhaseStats>,
    blend_lines: usize,
    sample_margin_frac: f32,
    line_margin_frac: f32,
    coverage_thresh: f64,
    overlap_db_thresh: f64,
    cv_warn_thresh: f64,
) -> Step2Diagnostics {
    let (covered, uncovered, frac, hist) = coverage_stats(hit_count);
    let overlap = overlap_metrics(burst_info, slc, blend_lines, sample_margin_frac);
    let means = per_burst_means(burst_info, slc, sample_margin_frac, line_margin_frac);
    let cv = coeff_of_variation(&means);

    for m in &overlap {
        let status = if m.delta_db.abs() <= overlap_db_thresh { "✅" } else { "⚠️" };
        log::info!(
            "STEP2|OVERLAP {} bursts {}→{} delta_db={:.3}dB ratio={:.4}",
            status, m.burst_a, m.burst_b, m.delta_db, m.ratio
        );
    }

    let coverage_ok = frac >= coverage_thresh;
    log::info!(
        "STEP2|COVERAGE {} covered={} uncovered={} frac={:.5} hist0/1/2+={}/{}/{}",
        if coverage_ok { "✅" } else { "❌" },
        covered,
        uncovered,
        frac,
        hist[0],
        hist[1],
        hist[2]
    );

    let cv_pct = cv * 100.0;
    let cv_status = if cv <= cv_warn_thresh { "✅" } else { "⚠️" };
    log::info!("STEP2|POWER_CV {} cv={:.2}% (threshold {:.2}%)", cv_status, cv_pct, cv_warn_thresh * 100.0);

    if let Some(d) = deramp_stats {
        log::info!(
            "STEP2|DERAMP reduction_factor={:.3} pre_slope={:.6} post_slope={:.6} samples={} lines={}",
            d.reduction_factor, d.pre_slope, d.post_slope, d.samples_used, d.lines_used
        );
    }

    Step2Diagnostics {
        coverage_fraction: frac,
        covered_pixels: covered,
        uncovered_pixels: uncovered,
        contrib_histogram: hist,
        overlap_metrics: overlap,
        burst_means: means,
        cv,
        deramp_stats: deramp_stats.cloned(),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::deburst::iw_deburst::types::BurstInfo;
    use crate::types::SarComplex;

    fn make_burst(burst_id: usize, start_line: usize, end_line: usize, start_sample: usize, end_sample: usize) -> BurstInfo {
        let lines = end_line.saturating_sub(start_line) + 1;
        BurstInfo {
            burst_id,
            start_line,
            end_line,
            start_sample,
            end_sample,
            azimuth_time: String::new(),
            sensing_time: String::new(),
            first_valid_sample: vec![start_sample as i32; lines],
            last_valid_sample: vec![end_sample as i32; lines],
            byte_offset: 0,
            // Extended fields
            azimuth_fm_rate: 0.0,
            azimuth_steering_rate: 0.0,
            slant_range_time: 0.0,
            doppler_centroid: 0.0,
            azimuth_bandwidth: 0.0,
            range_sampling_rate: 0.0,
            range_pixel_spacing: 1.0,
            azimuth_pixel_spacing: 1.0,
            azimuth_time_interval: 0.0,
            dc_polynomial: vec![0.0],
            fm_polynomial: vec![0.0],
            dc_range_poly: None,
            fm_range_poly: None,
            dc_polynomial_t0: None,
            fm_polynomial_t0: None,
            burst_reference_time_seconds: None,
            burst_azimuth_time_seconds: None,
            burst_start_time_utc: None,
            next_burst_start_time_utc: None,
            dc_selection_lower_idx: None,
            dc_selection_upper_idx: None,
            dc_selection_weight: None,
            fm_selection_lower_idx: None,
            fm_selection_upper_idx: None,
            fm_selection_weight: None,
        }
    }

    #[test]
    fn test_coverage_histogram() {
        let hit = Array2::from_shape_vec((2, 3), vec![0u16, 1, 2, 0, 1, 1]).unwrap();
        let (covered, uncovered, frac, hist) = coverage_stats(&hit);
        assert_eq!(covered, 4);
        assert_eq!(uncovered, 2);
        assert!((frac - 4.0 / 6.0).abs() < 1e-6);
        assert_eq!(hist, [2, 3, 1]);
    }

    #[test]
    fn test_overlap_metric_basic() {
        // Two bursts with simple overlaps; use unity power so ratio == 1.
        let burst_a = make_burst(0, 0, 3, 0, 2);
        let burst_b = make_burst(1, 3, 6, 0, 2);
        let burst_info = vec![burst_a, burst_b];

        let mut data = Array2::<SarComplex>::zeros((6, 3));
        for v in data.iter_mut() {
            *v = SarComplex { re: 1.0, im: 0.0 };
        }

        let m = overlap_metrics(&burst_info, &data, 2, 0.0);
        assert_eq!(m.len(), 1);
        assert!((m[0].delta_db).abs() < 1e-6);
        assert!((m[0].ratio - 1.0).abs() < 1e-6);
    }
}
