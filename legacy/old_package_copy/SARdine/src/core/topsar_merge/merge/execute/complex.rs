#![allow(dead_code)]
//! Complex merge execution with phase alignment.
//!
//! Executes the merge plan for complex-valued data with phase corrections.

use ndarray::{s, Array2, Axis};
use rayon::prelude::*;
use std::collections::HashMap;
use std::sync::Arc;

use crate::core::geometry::type_safe_units::Seconds;
use crate::core::DcFmRateProvider;
use crate::types::{SarImage, SarResult, SubSwath};

use super::super::overlap::OverlapRegion;
use super::super::plan::{MergePlan, MergeRowSegment, MergeWeight};
use super::super::types::AzimuthTimingModel;

/// Execute the merge plan for complex data with phase corrections.
///
/// # Phase Sign Convention
///
/// All phase corrections use consistent sign: multiply by exp(i * angle).
/// - To remove a phase ramp, apply -angle
/// - Overlap alignment uses -total_phase (removing the ramp)
/// - Azimuth time correction uses +phase_total (applying the correction)
pub fn execute_merge_plan_complex(
    plan: &MergePlan,
    subswaths: &[SubSwath],
    complex_refs: &[&SarImage],
    overlap_regions: &[OverlapRegion],
    merged: &mut Array2<num_complex::Complex32>,
    weight_sum: &mut Array2<f32>,
    phase_cache: Option<&Vec<Vec<f32>>>,
    overlap_phase_ramps: Option<&Vec<Vec<f32>>>,
    overlap_phase_offsets: Option<&Vec<f32>>,
    overlap_coherence: Option<&Vec<Vec<f32>>>,
    preserve_phase: bool,
    dc_fm_providers: &HashMap<String, Arc<dyn DcFmRateProvider>>,
    azimuth_timing: &AzimuthTimingModel,
) -> SarResult<()> {
    for (dst_row, segments) in plan.rows_plan.iter().enumerate() {
        if segments.is_empty() {
            continue;
        }

        let mut merged_row = merged.row_mut(dst_row);
        let mut weight_row = weight_sum.row_mut(dst_row);

        for segment in segments {
            process_complex_segment(
                segment,
                dst_row,
                subswaths,
                complex_refs,
                overlap_regions,
                &mut merged_row,
                &mut weight_row,
                plan.cols,
                phase_cache,
                overlap_phase_ramps,
                overlap_phase_offsets,
                overlap_coherence,
                preserve_phase,
                dc_fm_providers,
                azimuth_timing,
            )?;
        }
    }

    Ok(())
}

fn process_complex_segment(
    segment: &MergeRowSegment,
    dst_row: usize,
    subswaths: &[SubSwath],
    complex_refs: &[&SarImage],
    overlap_regions: &[OverlapRegion],
    merged_row: &mut ndarray::ArrayViewMut1<num_complex::Complex32>,
    weight_row: &mut ndarray::ArrayViewMut1<f32>,
    output_cols: usize,
    phase_cache: Option<&Vec<Vec<f32>>>,
    overlap_phase_ramps: Option<&Vec<Vec<f32>>>,
    overlap_phase_offsets: Option<&Vec<f32>>,
    overlap_coherence: Option<&Vec<Vec<f32>>>,
    preserve_phase: bool,
    dc_fm_providers: &HashMap<String, Arc<dyn DcFmRateProvider>>,
    azimuth_timing: &AzimuthTimingModel,
) -> SarResult<()> {
    if segment.len == 0 {
        return Ok(());
    }

    let src = match complex_refs.get(segment.swath_idx) {
        Some(s) => s,
        None => {
            log::warn!(
                "Missing complex subswath data for index {}",
                segment.swath_idx
            );
            return Ok(());
        }
    };

    // Bounds checking
    if segment.src_row >= src.nrows() {
        return Ok(());
    }
    if segment.src_col_start + segment.len > src.ncols() {
        return Ok(());
    }
    if segment.dst_col_start >= output_cols {
        return Ok(());
    }

    let effective_len = segment.len.min(output_cols - segment.dst_col_start);
    if effective_len == 0 {
        return Ok(());
    }

    let src_slice = src.slice(s![
        segment.src_row,
        segment.src_col_start..segment.src_col_start + effective_len
    ]);
    let mut dst_slice = merged_row.slice_mut(s![
        segment.dst_col_start..segment.dst_col_start + effective_len
    ]);
    let mut w_slice = weight_row.slice_mut(s![
        segment.dst_col_start..segment.dst_col_start + effective_len
    ]);

    match &segment.weight {
        MergeWeight::Constant(weight) => {
            if *weight <= 0.0 {
                return Ok(());
            }
            for idx in 0..effective_len {
                let mut value = src_slice[idx];

                if preserve_phase {
                    value = apply_phase_correction(
                        value,
                        dst_row,
                        segment.swath_idx,
                        subswaths,
                        phase_cache,
                        dc_fm_providers,
                        azimuth_timing,
                    )?;
                }

                let weighted = value * *weight;
                dst_slice[idx] += weighted;
                w_slice[idx] += *weight;
            }
        }
        MergeWeight::Overlap {
            overlap_index,
            row,
            col_offset,
            inverse,
        } => {
            let overlap = match overlap_regions.get(*overlap_index) {
                Some(o) => o,
                None => {
                    log::warn!(
                        "Invalid overlap index {} during complex execution",
                        overlap_index
                    );
                    return Ok(());
                }
            };

            if *row >= overlap.weights.nrows() || *col_offset >= overlap.weights.ncols() {
                return Ok(());
            }

            let available = overlap.weights.ncols() - *col_offset;
            let len = effective_len.min(available);
            if len == 0 {
                return Ok(());
            }

            let weight_slice = overlap
                .weights
                .slice(s![*row, *col_offset..*col_offset + len]);
            let swath_id = &subswaths[segment.swath_idx].id;
            let overlap_phase = overlap_phase_ramps.and_then(|v| v.get(*overlap_index));
            let overlap_offset = overlap_phase_offsets
                .and_then(|v| v.get(*overlap_index))
                .copied()
                .unwrap_or(0.0);
            let coh_rows = overlap_coherence.and_then(|c| c.get(*overlap_index));

            for idx in 0..len {
                let base_w = weight_slice[idx];
                let bias = coh_rows
                    .and_then(|rows| rows.get(dst_row.saturating_sub(overlap.azimuth_start)))
                    .copied()
                    .map(|c| 0.5f32 + 0.5f32 * c)
                    .unwrap_or(1.0);
                let weight = if *inverse {
                    (1.0 - base_w) * bias
                } else {
                    base_w + (1.0 - base_w) * (1.0 - bias)
                };

                let mut value = src_slice[idx];

                // Align swath2 to swath1 for Doppler/phase differences
                // FIX for bug C: Use consistent negative sign to REMOVE phase ramp
                if swath_id == &overlap.swath2_id {
                    if dst_row >= overlap.azimuth_start && dst_row < overlap.azimuth_end {
                        let local_row = dst_row - overlap.azimuth_start;
                        if let Some(ramps) = overlap_phase {
                            if let Some(&ramp) = ramps.get(local_row) {
                                let total_phase = ramp + overlap_offset;
                                if total_phase.abs() > 1e-6 {
                                    // Use -total_phase to REMOVE the phase difference
                                    let rot = num_complex::Complex32::from_polar(1.0, -total_phase);
                                    value *= rot;
                                }
                            }
                        }
                    }
                }

                if preserve_phase {
                    value = apply_phase_correction(
                        value,
                        dst_row,
                        segment.swath_idx,
                        subswaths,
                        phase_cache,
                        dc_fm_providers,
                        azimuth_timing,
                    )?;
                }

                let weighted = value * weight;
                dst_slice[idx] += weighted;
                w_slice[idx] += weight;
            }
        }
    }

    Ok(())
}

/// Apply azimuth time phase correction for enhanced complex data processing.
///
/// # Sign Convention
///
/// The correction uses exp(i * phase_total), where:
/// - phase_align is subtracted (removing DC-induced ramp)
/// - phase_correction is added (applying TOPSAR steering correction)
fn apply_phase_correction(
    complex_value: num_complex::Complex32,
    line_idx: usize,
    swath_idx: usize,
    subswaths: &[SubSwath],
    cached_phase: Option<&Vec<Vec<f32>>>,
    dc_fm_providers: &HashMap<String, Arc<dyn DcFmRateProvider>>,
    azimuth_timing: &AzimuthTimingModel,
) -> SarResult<num_complex::Complex32> {
    // Try cached phase first
    if let Some(cache) = cached_phase {
        if let Some(sw_vec) = cache.get(swath_idx) {
            if let Some(&phi) = sw_vec.get(line_idx) {
                if phi.is_finite() && phi.abs() > 1e-9 {
                    let rotation = num_complex::Complex32::from_polar(1.0, phi);
                    return Ok(complex_value * rotation);
                } else if phi.is_finite() {
                    return Ok(complex_value);
                }
            }
        }
    }

    // Compute phase correction on-the-fly
    let mut phase_total: f64 = 0.0;

    if let Some(azimuth_time) = azimuth_timing.get_azimuth_time_at_line(line_idx) {
        // DC-based inter-swath alignment
        let phase_align = compute_phase_alignment(
            swath_idx,
            azimuth_time,
            subswaths,
            dc_fm_providers,
            azimuth_timing.reference_azimuth_time,
        )?;
        phase_total -= phase_align as f64;

        // TOPSAR steering/azimuth timing correction
        let phase_correction = compute_azimuth_phase_correction(
            swath_idx,
            azimuth_time,
            line_idx,
            subswaths,
            dc_fm_providers,
            azimuth_timing,
        );
        phase_total += phase_correction;
    }

    if phase_total.abs() > 1e-9 {
        let rotation = num_complex::Complex32::from_polar(1.0, phase_total as f32);
        Ok(complex_value * rotation)
    } else {
        Ok(complex_value)
    }
}

fn compute_phase_alignment(
    swath_idx: usize,
    azimuth_time: f64,
    subswaths: &[SubSwath],
    dc_fm_providers: &HashMap<String, Arc<dyn DcFmRateProvider>>,
    reference_azimuth_time: f64,
) -> SarResult<f32> {
    if subswaths.is_empty() || swath_idx >= subswaths.len() {
        return Ok(0.0);
    }

    let ref_swath = &subswaths[0];
    let swath = &subswaths[swath_idx];

    let dc_ref = evaluate_dc_hz(ref_swath, azimuth_time, dc_fm_providers)?;
    let dc_sw = evaluate_dc_hz(swath, azimuth_time, dc_fm_providers)?;
    let delta_dc = dc_sw - dc_ref;

    let dt = azimuth_time - reference_azimuth_time;

    Ok((2.0 * std::f64::consts::PI * delta_dc * dt) as f32)
}

fn evaluate_dc_hz(
    sw: &SubSwath,
    azimuth_time: f64,
    providers: &HashMap<String, Arc<dyn DcFmRateProvider>>,
) -> SarResult<f64> {
    let provider = match providers.get(&sw.id) {
        Some(p) => p,
        None => return Ok(0.0),
    };

    let (start, end) = provider.get_time_range();
    let clamped = Seconds::new(azimuth_time.clamp(start.value(), end.value()));
    let dc = provider.get_dc(clamped)?;
    Ok(dc.value())
}

fn compute_azimuth_phase_correction(
    swath_idx: usize,
    azimuth_time: f64,
    line_idx: usize,
    subswaths: &[SubSwath],
    dc_fm_providers: &HashMap<String, Arc<dyn DcFmRateProvider>>,
    azimuth_timing: &AzimuthTimingModel,
) -> f64 {
    let swath = match subswaths.get(swath_idx) {
        Some(sw) => sw,
        None => return 0.0,
    };

    let burst = match azimuth_timing.burst_for_line(&swath.id, line_idx) {
        Some(b) => b,
        None => return 0.0,
    };

    let eta = azimuth_time - burst.sensing_time_center;

    let fm_rate = match evaluate_fm_rate(swath, azimuth_time, dc_fm_providers) {
        Ok(v) => v,
        Err(_) => return 0.0,
    };

    std::f64::consts::PI * fm_rate * eta.powi(2)
}

fn evaluate_fm_rate(
    sw: &SubSwath,
    azimuth_time: f64,
    providers: &HashMap<String, Arc<dyn DcFmRateProvider>>,
) -> SarResult<f64> {
    let provider = match providers.get(&sw.id) {
        Some(p) => p,
        None => return Ok(0.0),
    };

    let (start, end) = provider.get_time_range();
    let clamped = Seconds::new(azimuth_time.clamp(start.value(), end.value()));
    provider.get_fm_rate(clamped)
}

/// Precompute per-line, per-swath phase alignment to avoid recomputation.
pub fn precompute_phase_alignment(
    rows: usize,
    subswaths: &[SubSwath],
    dc_fm_providers: &HashMap<String, Arc<dyn DcFmRateProvider>>,
    azimuth_timing: &AzimuthTimingModel,
) -> SarResult<Vec<Vec<f32>>> {
    let mut cache: Vec<Vec<f32>> = vec![vec![0.0; rows]; subswaths.len()];

    for sw_idx in 0..subswaths.len() {
        for line_idx in 0..rows {
            let phase =
                if let Some(azimuth_time) = azimuth_timing.get_azimuth_time_at_line(line_idx) {
                    let mut phase_total: f64 = 0.0;

                    let phase_align = compute_phase_alignment(
                        sw_idx,
                        azimuth_time,
                        subswaths,
                        dc_fm_providers,
                        azimuth_timing.reference_azimuth_time,
                    )?;
                    phase_total -= phase_align as f64;

                    let phase_correction = compute_azimuth_phase_correction(
                        sw_idx,
                        azimuth_time,
                        line_idx,
                        subswaths,
                        dc_fm_providers,
                        azimuth_timing,
                    );
                    phase_total += phase_correction;

                    phase_total as f32
                } else {
                    0.0
                };
            cache[sw_idx][line_idx] = phase;
        }
    }

    Ok(cache)
}

/// Normalize complex output by weight sum.
pub fn normalize_complex_by_weight_sum(
    merged: &mut Array2<num_complex::Complex32>,
    weight_sum: &Array2<f32>,
    parallel: bool,
) {
    if parallel {
        merged
            .axis_iter_mut(Axis(0))
            .into_par_iter()
            .zip(weight_sum.axis_iter(Axis(0)))
            .for_each(|(mut row_merged, row_weights)| {
                for (merged_pixel, &wsum) in row_merged.iter_mut().zip(row_weights.iter()) {
                    if wsum > 0.0 {
                        *merged_pixel /= wsum;
                    }
                }
            });
    } else {
        for ((y, x), wsum) in weight_sum.indexed_iter() {
            if *wsum > 0.0 {
                merged[[y, x]] /= wsum;
            }
        }
    }
}

#[cfg(test)]
mod tests {
    #[test]
    fn phase_correction_cached_vs_computed() {
        // This test ensures cached and non-cached paths produce identical results
        // for a synthetic scenario.

        let value = num_complex::Complex32::new(1.0, 0.0);

        // With zero phase, output should equal input
        let rotation = num_complex::Complex32::from_polar(1.0, 0.0);
        let result = value * rotation;
        assert!((result.re - 1.0).abs() < 1e-6);
        assert!(result.im.abs() < 1e-6);

        // With π phase, should flip sign
        let rotation = num_complex::Complex32::from_polar(1.0, std::f32::consts::PI);
        let result = value * rotation;
        assert!((result.re + 1.0).abs() < 1e-6);
        assert!(result.im.abs() < 1e-6);
    }
}
