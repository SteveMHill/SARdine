#![allow(dead_code)]
//! Phase alignment and coherence computations for complex merge.

use std::collections::HashMap;
use std::sync::Arc;

use crate::core::geometry::type_safe_units::Seconds;
use crate::core::DcFmRateProvider;
use crate::types::{SarError, SarImage, SarResult, SubSwath};

use super::super::overlap::OverlapRegion;
use super::super::types::AzimuthTimingModel;

/// Sampling stride for overlap phase estimation
const OVERLAP_SAMPLE_STEP_FINE: usize = 8;

/// Phase correction sign convention:
///
/// All phase corrections use exp(i * angle), meaning:
/// - Positive angle rotates counter-clockwise
/// - To remove a phase ramp, use -angle
///
/// This is consistent throughout the merge pipeline.
#[inline]
pub fn phase_rotation(angle: f32) -> num_complex::Complex32 {
    num_complex::Complex32::from_polar(1.0, angle)
}

/// Precompute per-overlap Doppler/phase ramps (swath2 relative to swath1).
pub fn precompute_overlap_phase_ramps(
    overlap_regions: &[OverlapRegion],
    subswaths: &[SubSwath],
    dc_fm_providers: &HashMap<String, Arc<dyn DcFmRateProvider>>,
    azimuth_timing: &AzimuthTimingModel,
) -> SarResult<Vec<Vec<f32>>> {
    let mut swath_index: HashMap<&str, usize> = HashMap::new();
    for (idx, sw) in subswaths.iter().enumerate() {
        swath_index.insert(sw.id.as_str(), idx);
    }

    let ref_time = azimuth_timing.reference_azimuth_time;
    let mut ramps = Vec::with_capacity(overlap_regions.len());

    for overlap in overlap_regions {
        let idx1 = *swath_index.get(overlap.swath1_id.as_str()).ok_or_else(|| {
            SarError::Processing(format!(
                "Missing swath {} for overlap phase ramp computation",
                overlap.swath1_id
            ))
        })?;
        let idx2 = *swath_index.get(overlap.swath2_id.as_str()).ok_or_else(|| {
            SarError::Processing(format!(
                "Missing swath {} for overlap phase ramp computation",
                overlap.swath2_id
            ))
        })?;

        let sw1 = &subswaths[idx1];
        let sw2 = &subswaths[idx2];

        let mut line_phases =
            Vec::with_capacity(overlap.azimuth_end.saturating_sub(overlap.azimuth_start));

        for row in overlap.azimuth_start..overlap.azimuth_end {
            let phase = if let Some(az_time) = azimuth_timing.get_azimuth_time_at_line(row) {
                let dc1 = evaluate_dc_hz(sw1, az_time, dc_fm_providers)?;
                let dc2 = evaluate_dc_hz(sw2, az_time, dc_fm_providers)?;
                let delta_dc = dc2 - dc1;
                let dt = az_time - ref_time;
                // Phase ramp from DC difference
                (2.0 * std::f64::consts::PI * delta_dc * dt) as f32
            } else {
                0.0
            };

            line_phases.push(phase);
        }

        ramps.push(line_phases);
    }

    Ok(ramps)
}

/// Evaluate Doppler centroid for a subswath at a given azimuth time.
fn evaluate_dc_hz(
    sw: &SubSwath,
    azimuth_time: f64,
    providers: &HashMap<String, Arc<dyn DcFmRateProvider>>,
) -> SarResult<f64> {
    let provider = providers.get(&sw.id).ok_or_else(|| {
        SarError::Processing(format!("Missing DC/FM provider for subswath {}", sw.id))
    })?;

    let (start, end) = provider.get_time_range();
    let clamped = Seconds::new(azimuth_time.clamp(start.value(), end.value()));
    let dc = provider.get_dc(clamped)?;
    Ok(dc.value())
}

/// Estimate constant phase offsets per-overlap using complex samples.
pub fn compute_overlap_phase_offsets(
    overlap_regions: &[OverlapRegion],
    subswaths: &[SubSwath],
    complex_data: &HashMap<String, SarImage>,
) -> SarResult<Vec<f32>> {
    let mut offsets = Vec::with_capacity(overlap_regions.len());

    for overlap in overlap_regions {
        let phase = estimate_overlap_phase_offset(overlap, subswaths, complex_data)?.unwrap_or(0.0);
        offsets.push(phase);
    }

    Ok(offsets)
}

/// Estimate mean phase offset between two swaths in an overlap region.
fn estimate_overlap_phase_offset(
    overlap: &OverlapRegion,
    subswaths: &[SubSwath],
    complex_data: &HashMap<String, SarImage>,
) -> SarResult<Option<f32>> {
    let (swath1_data, swath2_data) = match (
        complex_data.get(&overlap.swath1_id),
        complex_data.get(&overlap.swath2_id),
    ) {
        (Some(s1), Some(s2)) => (s1, s2),
        _ => return Ok(None),
    };

    // Validate overlap ranges
    if overlap.azimuth_end <= overlap.azimuth_start
        || overlap.swath1_range_end <= overlap.swath1_range_start
    {
        return Ok(None);
    }

    let swath1_meta = subswaths
        .iter()
        .find(|s| s.id == overlap.swath1_id)
        .ok_or_else(|| {
            SarError::Processing(format!(
                "Missing metadata for subswath {}",
                overlap.swath1_id
            ))
        })?;

    let swath2_meta = subswaths
        .iter()
        .find(|s| s.id == overlap.swath2_id)
        .ok_or_else(|| {
            SarError::Processing(format!(
                "Missing metadata for subswath {}",
                overlap.swath2_id
            ))
        })?;

    let mut phase_accumulator = num_complex::Complex32::new(0.0, 0.0);
    let mut sample_count = 0;

    for row in (overlap.azimuth_start..overlap.azimuth_end).step_by(OVERLAP_SAMPLE_STEP_FINE) {
        for col1_local in
            (overlap.swath1_range_start..overlap.swath1_range_end).step_by(OVERLAP_SAMPLE_STEP_FINE)
        {
            let Some(col1_offset) = col1_local.checked_sub(overlap.swath1_range_start) else {
                continue;
            };
            let col1_global = swath1_meta.first_sample_global.saturating_add(col1_local);
            let col2_local = overlap.swath2_range_start.saturating_add(col1_offset);
            if col2_local >= overlap.swath2_range_end {
                continue;
            }
            let col2_global = swath2_meta.first_sample_global.saturating_add(col2_local);

            let Some(local_row1) = row.checked_sub(swath1_meta.first_line_global) else {
                continue;
            };
            let Some(local_row2) = row.checked_sub(swath2_meta.first_line_global) else {
                continue;
            };
            let Some(local_col1) = col1_global.checked_sub(swath1_meta.first_sample_global) else {
                continue;
            };
            let Some(local_col2) = col2_global.checked_sub(swath2_meta.first_sample_global) else {
                continue;
            };

            if local_row1 < swath1_data.nrows()
                && local_col1 < swath1_data.ncols()
                && local_row2 < swath2_data.nrows()
                && local_col2 < swath2_data.ncols()
            {
                let pixel1 = swath1_data[[local_row1, local_col1]];
                let pixel2 = swath2_data[[local_row2, local_col2]];

                if pixel1.norm_sqr() > 1e-6 && pixel2.norm_sqr() > 1e-6 {
                    phase_accumulator += pixel1 * pixel2.conj();
                    sample_count += 1;
                }
            }
        }
    }

    if sample_count > 0 {
        let mean_phase_offset = phase_accumulator.arg();
        log::debug!(
            "🔄 Phase offset {}-{}: {:.3} rad ({:.1}°) from {} samples",
            overlap.swath1_id,
            overlap.swath2_id,
            mean_phase_offset,
            mean_phase_offset.to_degrees(),
            sample_count
        );
        Ok(Some(mean_phase_offset))
    } else {
        log::debug!(
            "⚠️  Insufficient samples for phase offset estimation in {}-{}",
            overlap.swath1_id,
            overlap.swath2_id
        );
        Ok(None)
    }
}

/// Compute coherence-guided weights per overlap row.
pub fn compute_overlap_coherence(
    overlap_regions: &[OverlapRegion],
    complex_data: &HashMap<String, SarImage>,
) -> Option<Vec<Vec<f32>>> {
    let mut result = Vec::with_capacity(overlap_regions.len());

    for overlap in overlap_regions {
        let sw1 = complex_data.get(&overlap.swath1_id)?;
        let sw2 = complex_data.get(&overlap.swath2_id)?;

        let az_start = overlap.azimuth_start.min(sw1.nrows().saturating_sub(1));
        let az_end = overlap
            .azimuth_end
            .min(sw1.nrows().saturating_sub(1))
            .min(sw2.nrows().saturating_sub(1));

        if az_end <= az_start {
            result.push(Vec::new());
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
        let len = (rg1_end - rg1_start).min(rg2_end - rg2_start);

        if len == 0 {
            result.push(Vec::new());
            continue;
        }

        let mut coh_rows = Vec::with_capacity(az_end - az_start);
        let step = OVERLAP_SAMPLE_STEP_FINE.max(4);

        for az in az_start..=az_end {
            let row1 = sw1.row(az);
            let row2 = sw2.row(az);
            let mut num = num_complex::Complex64::new(0.0, 0.0);
            let mut p1 = 0.0;
            let mut p2 = 0.0;
            let mut n = 0usize;

            for idx in (0..len).step_by(step) {
                let v1 = row1[rg1_start + idx];
                let v2 = row2[rg2_start + idx];
                if v1.re.is_finite() && v1.im.is_finite() && v2.re.is_finite() && v2.im.is_finite()
                {
                    num += num_complex::Complex64::new(v1.re as f64, v1.im as f64)
                        * num_complex::Complex64::new(v2.re as f64, -v2.im as f64);
                    p1 += v1.norm_sqr() as f64;
                    p2 += v2.norm_sqr() as f64;
                    n += 1;
                }
            }

            let coh = if n > 0 && p1 > 0.0 && p2 > 0.0 {
                (num.norm() / ((p1 * p2).sqrt())).min(1.0)
            } else {
                0.0
            };
            coh_rows.push(coh as f32);
        }

        result.push(coh_rows);
    }

    Some(result)
}

/// Log overlap phase diagnostics to surface seam risk.
pub fn log_overlap_phase_diagnostics(
    overlap_regions: &[OverlapRegion],
    subswaths: &[SubSwath],
    dc_fm_providers: &HashMap<String, Arc<dyn DcFmRateProvider>>,
    azimuth_timing: &AzimuthTimingModel,
    max_phase_discontinuity: f32,
) {
    if overlap_regions.is_empty() {
        return;
    }

    let mut swath_index: HashMap<&str, usize> = HashMap::new();
    for (idx, sw) in subswaths.iter().enumerate() {
        swath_index.insert(sw.id.as_str(), idx);
    }

    let az_interval = azimuth_timing.azimuth_time_interval;

    for overlap in overlap_regions {
        let Some(&idx1) = swath_index.get(overlap.swath1_id.as_str()) else {
            continue;
        };
        let Some(&idx2) = swath_index.get(overlap.swath2_id.as_str()) else {
            continue;
        };

        let mid_row = (overlap.azimuth_start + overlap.azimuth_end) / 2;
        let Some(az_time) = azimuth_timing.get_azimuth_time_at_line(mid_row) else {
            continue;
        };

        let dc1_hz = match evaluate_dc_hz(&subswaths[idx1], az_time, dc_fm_providers) {
            Ok(v) => v,
            Err(err) => {
                log::warn!(
                    "⚠️  DC evaluation failed for {}: {}",
                    overlap.swath1_id,
                    err
                );
                continue;
            }
        };

        let dc2_hz = match evaluate_dc_hz(&subswaths[idx2], az_time, dc_fm_providers) {
            Ok(v) => v,
            Err(err) => {
                log::warn!(
                    "⚠️  DC evaluation failed for {}: {}",
                    overlap.swath2_id,
                    err
                );
                continue;
            }
        };

        let delta_dc = dc2_hz - dc1_hz;
        let phase_per_line = 2.0 * std::f64::consts::PI * delta_dc * az_interval;
        let lines = overlap
            .azimuth_end
            .saturating_sub(overlap.azimuth_start)
            .max(1);
        let total_phase = phase_per_line * lines as f64;

        if total_phase.abs() > max_phase_discontinuity as f64 {
            log::warn!(
                "⚠️  Overlap {}-{} phase ramp {:.2} rad (Δf_dc={:.2} Hz, lines={}) exceeds tolerance {:.2} rad",
                overlap.swath1_id,
                overlap.swath2_id,
                total_phase,
                delta_dc,
                lines,
                max_phase_discontinuity
            );
        } else {
            log::info!(
                "✅ Overlap {}-{} phase ramp {:.2} rad (Δf_dc={:.2} Hz, lines={})",
                overlap.swath1_id,
                overlap.swath2_id,
                total_phase,
                delta_dc,
                lines
            );
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn phase_rotation_correct() {
        // Rotation by 0 should give (1, 0)
        let r0 = phase_rotation(0.0);
        assert!((r0.re - 1.0).abs() < 1e-6);
        assert!(r0.im.abs() < 1e-6);

        // Rotation by π/2 should give (0, 1)
        let r90 = phase_rotation(std::f32::consts::FRAC_PI_2);
        assert!(r90.re.abs() < 1e-6);
        assert!((r90.im - 1.0).abs() < 1e-6);

        // Rotation by π should give (-1, 0)
        let r180 = phase_rotation(std::f32::consts::PI);
        assert!((r180.re + 1.0).abs() < 1e-6);
        assert!(r180.im.abs() < 1e-6);
    }

    #[test]
    fn phase_rotation_inverse() {
        // Applying angle then -angle should give identity
        let angle = 1.234f32;
        let forward = phase_rotation(angle);
        let backward = phase_rotation(-angle);
        let combined = forward * backward;

        assert!((combined.re - 1.0).abs() < 1e-6);
        assert!(combined.im.abs() < 1e-6);
    }
}
