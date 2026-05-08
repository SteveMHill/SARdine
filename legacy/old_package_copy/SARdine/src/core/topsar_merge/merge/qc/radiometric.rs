#![allow(dead_code)]
//! Radiometric consistency validation and gain computation.
//!
//! Provides overlap gain calculation for radiometric equalization
//! between adjacent subswaths.

use ndarray::Array2;
use std::cmp::Ordering;
use std::collections::HashMap;
use std::env;

use crate::types::{SarRealImage, SarResult};

use super::super::overlap::OverlapRegion;

/// Sampling stride for radiometric consistency checks
const OVERLAP_SAMPLE_STEP_RADIO: usize = 10;

/// Radiometric consistency threshold in dB
const RADIOMETRIC_CONSISTENCY_DB_THRESHOLD: f32 = 0.2;

// Default gain clamp bounds for overlap equalization. Can be overridden via env.
// Tighter bounds (≈ -1 dB to +1 dB) keep overlap normalization modest and
// avoid introducing bright bands when cross-swath offsets are large.
const GAIN_CLAMP_MIN_DEFAULT: f32 = 0.8;
const GAIN_CLAMP_MAX_DEFAULT: f32 = 1.25;

fn env_usize(name: &str) -> Option<usize> {
    env::var(name)
        .ok()
        .and_then(|v| v.parse::<usize>().ok())
        .filter(|&v| v > 0)
}

fn env_f32(name: &str) -> Option<f32> {
    env::var(name)
        .ok()
        .and_then(|v| v.parse::<f32>().ok())
        .filter(|&v| v.is_finite())
}

fn gain_clamp_bounds() -> (f32, f32) {
    let min = env_f32("SARDINE_OVERLAP_GAIN_CLAMP_MIN").unwrap_or(GAIN_CLAMP_MIN_DEFAULT);
    let max = env_f32("SARDINE_OVERLAP_GAIN_CLAMP_MAX").unwrap_or(GAIN_CLAMP_MAX_DEFAULT);
    if max > min {
        (min, max)
    } else {
        (GAIN_CLAMP_MIN_DEFAULT, GAIN_CLAMP_MAX_DEFAULT)
    }
}

fn overlap_sample_step_radio() -> usize {
    env_usize("SARDINE_OVERLAP_SAMPLE_STEP_RADIO").unwrap_or(OVERLAP_SAMPLE_STEP_RADIO)
}

fn overlap_sample_step_fine() -> usize {
    // Used for denser sampling in overlap gain computation
    env_usize("SARDINE_OVERLAP_SAMPLE_STEP_FINE").unwrap_or(4)
}

/// Proper Pool-Adjacent-Violators (PAV) isotonic regression.
///
/// Returns a monotonically non-decreasing sequence of the same length as input.
/// This is a correct implementation that preserves output length.
pub fn isotonic_non_decreasing(values: &[f64]) -> Vec<f64> {
    if values.is_empty() {
        return Vec::new();
    }

    let n = values.len();

    // Each block stores (sum, count, start_idx, end_idx_exclusive)
    // We'll merge blocks when violating monotonicity
    struct Block {
        sum: f64,
        count: usize,
        start: usize,
        end: usize, // exclusive
    }

    impl Block {
        fn mean(&self) -> f64 {
            self.sum / self.count as f64
        }
    }

    let mut blocks: Vec<Block> = values
        .iter()
        .enumerate()
        .map(|(i, &v)| Block {
            sum: v,
            count: 1,
            start: i,
            end: i + 1,
        })
        .collect();

    // Pool adjacent violators
    let mut i = 0;
    while i + 1 < blocks.len() {
        if blocks[i].mean() > blocks[i + 1].mean() {
            // Merge blocks i and i+1
            blocks[i].sum += blocks[i + 1].sum;
            blocks[i].count += blocks[i + 1].count;
            blocks[i].end = blocks[i + 1].end;
            blocks.remove(i + 1);

            // Check backward for new violations
            if i > 0 {
                i -= 1;
            }
        } else {
            i += 1;
        }
    }

    // Expand blocks back to original indices
    let mut result = vec![0.0; n];
    for block in &blocks {
        let mean = block.mean();
        for idx in block.start..block.end {
            result[idx] = mean;
        }
    }

    result
}

/// Smooth gain curve using median filter and isotonic regression.
pub fn smooth_gain_curve(raw: &[f32]) -> Vec<f32> {
    if raw.is_empty() {
        return Vec::new();
    }
    let n = raw.len();
    let (clamp_min, clamp_max) = gain_clamp_bounds();

    // Convert to log domain for robust statistics
    let log_vals: Vec<f64> = raw.iter().map(|v| (*v as f64).max(1e-4).ln()).collect();

    // Median filter with outlier-robust averaging
    let window = n.min(21).max(5);
    let mut smoothed = vec![0.0f64; n];

    for i in 0..n {
        let start = i.saturating_sub(window / 2);
        let end = (i + window / 2 + 1).min(n);
        let slice = &log_vals[start..end];

        let mut med_vec = slice.to_vec();
        med_vec.sort_by(|a, b| a.partial_cmp(b).unwrap_or(Ordering::Equal));
        let median = med_vec[med_vec.len() / 2];

        // Outlier-robust weighted mean
        let mad = med_vec
            .iter()
            .map(|v| (v - median).abs())
            .fold(0.0, f64::max);
        let delta = mad.max(0.05);

        let mut num = 0.0;
        let mut den = 0.0;
        for &v in slice {
            let r = v - median;
            let w = if r.abs() <= delta {
                1.0
            } else {
                delta / r.abs()
            };
            num += w * v;
            den += w;
        }
        smoothed[i] = if den > 0.0 { num / den } else { median };
    }

    // Apply isotonic regression for monotonicity
    let mono = isotonic_non_decreasing(&smoothed);

    // Convert back from log domain with clamping
    mono.iter()
        .map(|v| v.exp().clamp(clamp_min as f64, clamp_max as f64) as f32)
        .collect()
}

/// Compute per-overlap radiometric gain to equalize swath2 to swath1.
pub fn compute_overlap_gains(
    overlap_regions: &[OverlapRegion],
    subswath_data: &HashMap<String, SarRealImage>,
    mask_data: Option<&HashMap<String, Array2<u8>>>,
) -> SarResult<Vec<Vec<f32>>> {
    let mut gains = Vec::with_capacity(overlap_regions.len());

    for overlap in overlap_regions {
        log::info!(
            "🧭 Overlap {}-{}: azimuth {}..{}, sw1 range {}..{}, sw2 range {}..{}",
            overlap.swath1_id,
            overlap.swath2_id,
            overlap.azimuth_start,
            overlap.azimuth_end,
            overlap.swath1_range_start,
            overlap.swath1_range_end,
            overlap.swath2_range_start,
            overlap.swath2_range_end
        );

        let sw1 = match subswath_data.get(&overlap.swath1_id) {
            Some(s) => s,
            None => {
                log::warn!(
                    "Missing subswath {} for overlap gain compute",
                    overlap.swath1_id
                );
                gains.push(vec![1.0]);
                continue;
            }
        };
        let sw2 = match subswath_data.get(&overlap.swath2_id) {
            Some(s) => s,
            None => {
                log::warn!(
                    "Missing subswath {} for overlap gain compute",
                    overlap.swath2_id
                );
                gains.push(vec![1.0]);
                continue;
            }
        };

        let az_start = overlap.azimuth_start.min(sw1.nrows().saturating_sub(1));
        let az_end = overlap
            .azimuth_end
            .min(sw1.nrows().saturating_sub(1))
            .min(sw2.nrows().saturating_sub(1));
        if az_end <= az_start {
            gains.push(vec![1.0]);
            continue;
        }

        let rg1_start = overlap
            .swath1_range_start
            .min(sw1.ncols().saturating_sub(1));
        let rg1_end = overlap.swath1_range_end.min(sw1.ncols().saturating_sub(1));
        let rg2_start = overlap
            .swath2_range_start
            .min(sw2.ncols().saturating_sub(1));
        let rg2_end = overlap.swath2_range_end.min(sw2.ncols().saturating_sub(1));

        if rg1_end <= rg1_start || rg2_end <= rg2_start {
            gains.push(vec![1.0]);
            continue;
        }

        // Compute power floor to avoid low-SNR edges biasing gains
        let (thresh1, thresh2) = compute_power_thresholds(
            sw1, sw2, az_start, az_end, rg1_start, rg1_end, rg2_start, rg2_end,
        );

        // STATE-OF-THE-ART: Dynamic least-squares method for radiometric normalization
        // Based on research: "Dynamic least-squares method for additive noise removal in Sentinel-1 TOPSAR data"
        // This method adjusts scaling parameters for each subswath individually to minimize radiometric discrepancies

        // Step 1: Compute block-level gains for better spatial consistency
        let block_size = 50; // Process in blocks for spatial consistency
        let num_blocks = ((az_end - az_start + 1) + block_size - 1) / block_size;
        let mut block_gains: Vec<f32> = Vec::new();

        for block_idx in 0..num_blocks {
            let block_az_start = az_start + block_idx * block_size;
            let block_az_end = (block_az_start + block_size).min(az_end + 1);

            // Collect all valid ratios in this block
            let mut block_ratios: Vec<f32> = Vec::new();
            let sample_step = overlap_sample_step_fine(); // Denser sampling for better statistics

            for az in block_az_start..block_az_end {
                let row1 = sw1.row(az);
                let row2 = sw2.row(az);
                let mask_row1 = mask_data
                    .and_then(|m| m.get(&overlap.swath1_id))
                    .map(|a| a.row(az));
                let mask_row2 = mask_data
                    .and_then(|m| m.get(&overlap.swath2_id))
                    .map(|a| a.row(az));

                let len = (rg1_end - rg1_start)
                    .min(rg2_end - rg2_start)
                    .min(row1.len());
                if len == 0 {
                    continue;
                }

                let mut idx = 0usize;
                while rg1_start + idx < rg1_start + len && rg2_start + idx < rg2_start + len {
                    let v1 = row1[rg1_start + idx] as f64;
                    let v2 = row2[rg2_start + idx] as f64;

                    let mask_ok1 = mask_row1
                        .as_ref()
                        .map(|m| m[rg1_start + idx] > 0)
                        .unwrap_or(true);
                    let mask_ok2 = mask_row2
                        .as_ref()
                        .map(|m| m[rg2_start + idx] > 0)
                        .unwrap_or(true);

                    if v1.is_finite()
                        && v2.is_finite()
                        && v1 as f32 > thresh1
                        && v2 as f32 > thresh2
                        && mask_ok1
                        && mask_ok2
                    {
                        let ratio = (v1 / v2) as f32;
                        if ratio.is_finite() && ratio > 0.0 {
                            block_ratios.push(ratio);
                        }
                    }
                    idx += sample_step;
                }
            }

            // Dynamic least-squares: Use robust median for block gain
            let block_gain = if block_ratios.is_empty() {
                1.0
            } else {
                compute_trimmed_median(&mut block_ratios).clamp(0.5, 2.0)
            };
            block_gains.push(block_gain);
        }

        // Step 2: Compute per-row gains with block-level guidance
        let sample_step = overlap_sample_step_fine(); // Denser sampling
        let mut row_gains: Vec<f32> = Vec::new();

        for az in az_start..=az_end {
            let row1 = sw1.row(az);
            let row2 = sw2.row(az);
            let mask_row1 = mask_data
                .and_then(|m| m.get(&overlap.swath1_id))
                .map(|a| a.row(az));
            let mask_row2 = mask_data
                .and_then(|m| m.get(&overlap.swath2_id))
                .map(|a| a.row(az));

            let len = (rg1_end - rg1_start)
                .min(rg2_end - rg2_start)
                .min(row1.len());
            if len == 0 {
                // Use block-level gain if available
                let block_idx =
                    ((az - az_start) / block_size).min(block_gains.len().saturating_sub(1));
                row_gains.push(block_gains.get(block_idx).copied().unwrap_or(1.0));
                continue;
            }

            let mut ratios: Vec<f32> = Vec::new();
            let mut idx = 0usize;
            while rg1_start + idx < rg1_start + len && rg2_start + idx < rg2_start + len {
                let v1 = row1[rg1_start + idx] as f64;
                let v2 = row2[rg2_start + idx] as f64;

                let mask_ok1 = mask_row1
                    .as_ref()
                    .map(|m| m[rg1_start + idx] > 0)
                    .unwrap_or(true);
                let mask_ok2 = mask_row2
                    .as_ref()
                    .map(|m| m[rg2_start + idx] > 0)
                    .unwrap_or(true);

                if v1.is_finite()
                    && v2.is_finite()
                    && v1 as f32 > thresh1
                    && v2 as f32 > thresh2
                    && mask_ok1
                    && mask_ok2
                {
                    let ratio = (v1 / v2) as f32;
                    if ratio.is_finite() && ratio > 0.0 {
                        ratios.push(ratio);
                    }
                }
                idx += sample_step;
            }

            let gain = if ratios.is_empty() {
                // Fallback to block-level gain
                let block_idx =
                    ((az - az_start) / block_size).min(block_gains.len().saturating_sub(1));
                block_gains.get(block_idx).copied().unwrap_or(1.0)
            } else {
                // STATE-OF-THE-ART: Dynamic least-squares approach
                // Use trimmed median for robustness, but blend with block-level gain for spatial consistency
                let trimmed_median = compute_trimmed_median(&mut ratios);
                let block_idx =
                    ((az - az_start) / block_size).min(block_gains.len().saturating_sub(1));
                let block_gain = block_gains.get(block_idx).copied().unwrap_or(1.0);

                // Weighted combination: 70% row-level, 30% block-level for spatial smoothness
                let blended_gain = 0.7 * trimmed_median + 0.3 * block_gain;

                // Validate and clamp
                let (clamp_min, clamp_max) = gain_clamp_bounds();
                if blended_gain < clamp_min || blended_gain > clamp_max {
                    log::warn!(
                        "⚠️  Extreme gain ratio {:.3} for row {} in overlap {}-{} (row={:.3}, block={:.3}, ratios: {} samples)",
                        blended_gain, az, overlap.swath1_id, overlap.swath2_id, trimmed_median, block_gain, ratios.len()
                    );
                    blended_gain.clamp(clamp_min, clamp_max)
                } else {
                    blended_gain
                }
            };

            row_gains.push(gain);
        }

        if !row_gains.is_empty() {
            let row_gains = smooth_gain_curve(&row_gains);
            log_gain_stats(&overlap.swath1_id, &overlap.swath2_id, &row_gains);
            gains.push(row_gains);
        } else {
            gains.push(vec![1.0]);
        }
    }

    Ok(gains)
}

fn compute_power_thresholds(
    sw1: &SarRealImage,
    sw2: &SarRealImage,
    az_start: usize,
    az_end: usize,
    rg1_start: usize,
    rg1_end: usize,
    rg2_start: usize,
    rg2_end: usize,
) -> (f32, f32) {
    let mut samples1 = Vec::new();
    let mut samples2 = Vec::new();

    for az in (az_start..=az_end).step_by(32) {
        let row1 = sw1.row(az);
        let row2 = sw2.row(az);
        let len = (rg1_end - rg1_start)
            .min(rg2_end - rg2_start)
            .min(row1.len());

        let mut idx = 0usize;
        while rg1_start + idx < rg1_start + len && rg2_start + idx < rg2_start + len {
            let v1 = row1[rg1_start + idx];
            let v2 = row2[rg2_start + idx];
            if v1.is_finite() && v1 > 0.0 {
                samples1.push(v1);
            }
            if v2.is_finite() && v2 > 0.0 {
                samples2.push(v2);
            }
            idx += 64;
        }
    }

    let percentile_10 = |mut v: Vec<f32>| -> f32 {
        if v.is_empty() {
            return 0.0;
        }
        v.sort_by(|a, b| a.partial_cmp(b).unwrap_or(Ordering::Equal));
        let pos = ((v.len() as f32) * 0.1).floor() as usize;
        v[pos.min(v.len().saturating_sub(1))]
    };

    (percentile_10(samples1), percentile_10(samples2))
}

fn compute_trimmed_median(ratios: &mut Vec<f32>) -> f32 {
    let (clamp_min, clamp_max) = gain_clamp_bounds();
    ratios.sort_by(|a, b| a.partial_cmp(b).unwrap_or(Ordering::Equal));

    let len = ratios.len();
    let lo = (len as f32 * 0.10).floor() as usize;
    let hi = (len as f32 * 0.90).ceil() as usize;
    let lo = lo.min(len.saturating_sub(1));
    let hi = hi.max(lo + 1).min(len);
    let trimmed = &ratios[lo..hi];

    let mid = trimmed.len() / 2;
    let median = if trimmed.len() % 2 == 0 && trimmed.len() >= 2 {
        0.5 * (trimmed[mid - 1] + trimmed[mid])
    } else if !trimmed.is_empty() {
        trimmed[mid]
    } else {
        1.0
    };

    median.clamp(clamp_min, clamp_max)
}

fn log_gain_stats(swath1_id: &str, swath2_id: &str, gains: &[f32]) {
    let mut sorted = gains.to_vec();
    sorted.sort_by(|a, b| a.partial_cmp(b).unwrap_or(Ordering::Equal));
    let median = sorted[sorted.len() / 2];
    let mean = sorted.iter().copied().sum::<f32>() / sorted.len() as f32;
    let min = *sorted.first().unwrap_or(&1.0);
    let max = *sorted.last().unwrap_or(&1.0);
    log::info!(
        "   Gain stats {}-{}: mean={:.3}, median={:.3}, min={:.3}, max={:.3}, rows={}",
        swath1_id,
        swath2_id,
        mean,
        median,
        min,
        max,
        gains.len()
    );
}

/// Validate radiometric consistency in overlap regions.
pub fn validate_radiometric_consistency(
    overlaps: &[OverlapRegion],
    subswath_data: &HashMap<String, SarRealImage>,
    subswaths: &[crate::types::SubSwath],
) -> SarResult<()> {
    let subswath_lookup: HashMap<&str, &crate::types::SubSwath> =
        subswaths.iter().map(|s| (s.id.as_str(), s)).collect();

    for overlap in overlaps {
        if let (Some(swath1_data), Some(swath2_data)) = (
            subswath_data.get(&overlap.swath1_id),
            subswath_data.get(&overlap.swath2_id),
        ) {
            let Some(swath1_meta) = subswath_lookup.get(overlap.swath1_id.as_str()) else {
                continue;
            };
            let Some(swath2_meta) = subswath_lookup.get(overlap.swath2_id.as_str()) else {
                continue;
            };

            let mut differences = Vec::new();

            // FIX for bug G: Simplified validation without redundant overlap_range calls
            let az_start = overlap.azimuth_start;
            let az_end = overlap.azimuth_end;

            if az_end <= az_start {
                log::warn!(
                    "⚠️ Invalid azimuth range for overlap {}-{}: {}..{}",
                    overlap.swath1_id,
                    overlap.swath2_id,
                    az_start,
                    az_end
                );
                continue;
            }

            let sw1_start = overlap.swath1_range_start;
            let sw1_end = overlap.swath1_range_end;
            let sw2_start = overlap.swath2_range_start;
            let sw2_end = overlap.swath2_range_end;

            if sw1_end <= sw1_start || sw2_end <= sw2_start {
                continue;
            }

            let step = overlap_sample_step_radio();
            for row in (az_start..az_end).step_by(step) {
                for col1 in (sw1_start..sw1_end).step_by(step) {
                    let Some(col1_offset) = col1.checked_sub(sw1_start) else {
                        continue;
                    };
                    let Some(col2) = sw2_start.checked_add(col1_offset) else {
                        continue;
                    };
                    if col2 >= sw2_end {
                        continue;
                    }

                    let Some(row_local_1) = row.checked_sub(swath1_meta.first_line_global) else {
                        continue;
                    };
                    let Some(row_local_2) = row.checked_sub(swath2_meta.first_line_global) else {
                        continue;
                    };
                    let Some(col_local_1) = col1.checked_sub(swath1_meta.first_sample_global)
                    else {
                        continue;
                    };
                    let Some(col_local_2) = col2.checked_sub(swath2_meta.first_sample_global)
                    else {
                        continue;
                    };

                    if row_local_1 < swath1_data.nrows()
                        && col_local_1 < swath1_data.ncols()
                        && row_local_2 < swath2_data.nrows()
                        && col_local_2 < swath2_data.ncols()
                    {
                        let val1 = swath1_data[[row_local_1, col_local_1]];
                        let val2 = swath2_data[[row_local_2, col_local_2]];

                        if val1 > 0.0 && val2 > 0.0 {
                            let db_diff = 10.0 * (val1 / val2).log10();
                            differences.push(db_diff);
                        }
                    }
                }
            }

            if !differences.is_empty() {
                differences.sort_by(|a, b| a.total_cmp(b));
                let median_diff = differences[differences.len() / 2].abs();

                log::info!(
                    "📊 Radiometric consistency {}-{}: {:.3} dB median difference",
                    overlap.swath1_id,
                    overlap.swath2_id,
                    median_diff
                );

                if median_diff > RADIOMETRIC_CONSISTENCY_DB_THRESHOLD {
                    log::warn!(
                        "⚠️  Radiometric inconsistency detected: {:.3} dB > {:.1} dB threshold",
                        median_diff,
                        RADIOMETRIC_CONSISTENCY_DB_THRESHOLD
                    );
                }
            }
        }
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn isotonic_preserves_length() {
        let input = vec![5.0, 3.0, 4.0, 2.0, 6.0];
        let output = isotonic_non_decreasing(&input);

        assert_eq!(output.len(), input.len(), "Output length must match input");
    }

    #[test]
    fn isotonic_is_monotonic() {
        let input = vec![5.0, 3.0, 4.0, 2.0, 6.0, 1.0, 7.0];
        let output = isotonic_non_decreasing(&input);

        for i in 1..output.len() {
            assert!(
                output[i] >= output[i - 1],
                "Output should be non-decreasing: output[{}]={} < output[{}]={}",
                i - 1,
                output[i - 1],
                i,
                output[i]
            );
        }
    }

    #[test]
    fn isotonic_already_monotonic() {
        let input = vec![1.0, 2.0, 3.0, 4.0, 5.0];
        let output = isotonic_non_decreasing(&input);

        assert_eq!(output, input, "Already monotonic input should be unchanged");
    }

    #[test]
    fn isotonic_single_element() {
        let input = vec![3.14];
        let output = isotonic_non_decreasing(&input);

        assert_eq!(output, vec![3.14]);
    }

    #[test]
    fn isotonic_empty() {
        let input: Vec<f64> = vec![];
        let output = isotonic_non_decreasing(&input);

        assert!(output.is_empty());
    }

    #[test]
    fn isotonic_all_decreasing() {
        let input = vec![5.0, 4.0, 3.0, 2.0, 1.0];
        let output = isotonic_non_decreasing(&input);

        // All should be the mean: (5+4+3+2+1)/5 = 3.0
        let expected_mean = 3.0;
        for &v in &output {
            assert!((v - expected_mean).abs() < 1e-10);
        }
        assert_eq!(output.len(), 5);
    }

    #[test]
    fn smooth_gain_preserves_length() {
        let input = vec![1.1, 0.9, 1.0, 1.2, 0.8, 1.1];
        let output = smooth_gain_curve(&input);

        assert_eq!(output.len(), input.len());
    }

    #[test]
    fn smooth_gain_is_bounded() {
        let input = vec![0.1, 5.0, 1.0, 0.01, 10.0];
        let output = smooth_gain_curve(&input);

        for &v in &output {
            assert!(v >= 0.5 && v <= 2.0, "Gain {} should be in [0.5, 2.0]", v);
        }
    }
}
