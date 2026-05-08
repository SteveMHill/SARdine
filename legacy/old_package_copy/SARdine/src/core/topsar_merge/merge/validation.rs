//! Validation helpers for TOPSAR merge processing
//!
//! This module contains pre-merge validation functions for:
//! - DC polynomial prerequisites
//! - DC/FM provider coverage
//! - Grid alignment between subswaths
//! - Radiometric consistency in overlap regions

use crate::core::geometry::dc_fm_provider::DcFmRateProvider;
use crate::types::{SarError, SarRealImage, SarResult, SubSwath};
use std::collections::HashMap;
use std::sync::Arc;

use super::overlap::overlap_range;
use super::{
    AzimuthTimingModel, TopsarMerge, OVERLAP_SAMPLE_STEP_RADIO,
    RADIOMETRIC_CONSISTENCY_DB_THRESHOLD,
};

/// Validate that all subswaths have the DC polynomial prerequisites for merge.
///
/// Checks:
/// - Each subswath has a non-empty DC polynomial
/// - Each subswath has dc_polynomial_t0 reference time
/// - Optionally warns if azimuth_time_interval is missing
pub(crate) fn validate_dc_prerequisites(subswaths: &[SubSwath]) -> SarResult<()> {
    for sw in subswaths {
        match &sw.dc_polynomial {
            None => {
                return Err(SarError::Processing(format!(
                    "❌ CRITICAL: Subswath {} lacks DC polynomial!\n\
                     This indicates deburst did NOT use DC-aware processing.\n\
                     Inter-subswath alignment will be INCORRECT.\n\
                     Solution: Re-run deburst with DC extraction enabled.",
                    sw.id
                )));
            }
            Some(coeffs) if coeffs.is_empty() => {
                return Err(SarError::Processing(format!(
                    "❌ CRITICAL: Subswath {} has empty DC polynomial!\n\
                     This indicates DC extraction failed during deburst.\n\
                     Inter-subswath alignment will be INCORRECT.",
                    sw.id
                )));
            }
            Some(coeffs) => {
                log::info!(
                    "✅ Subswath {} has DC polynomial: {} coefficients",
                    sw.id,
                    coeffs.len()
                );
            }
        }

        if sw.dc_polynomial_t0.is_none() {
            return Err(SarError::Processing(format!(
                "❌ CRITICAL: Subswath {} missing dc_polynomial_t0 reference time; DC evaluation would silently fall back to zero.",
                sw.id
            )));
        }

        // Also check azimuth time interval (needed for timing-based alignment)
        if sw.azimuth_time_interval.is_none() {
            log::warn!(
                "⚠️ Subswath {} lacks azimuth time interval - will use PRF fallback",
                sw.id
            );
        }
    }

    log::info!("✅ DC prerequisites validated - proceeding with DC-aware merge");
    Ok(())
}

/// Validate that DC/FM providers cover all burst timing ranges.
///
/// Ensures that the polynomial providers can evaluate DC/FM rates for
/// all azimuth times encountered in the burst timing model.
pub(crate) fn validate_dc_fm_providers(
    providers: &HashMap<String, Arc<dyn DcFmRateProvider>>,
    azimuth_timing: &AzimuthTimingModel,
    subswaths: &[SubSwath],
) -> SarResult<()> {
    // Mirror build_dc_fm_provider_for_swath tolerance: allow ~2 ms drift
    const AZ_TIME_EPS: f64 = 2.0e-3;

    for sw in subswaths {
        let provider = providers.get(&sw.id).ok_or_else(|| {
            SarError::Processing(format!(
                "Missing DC/FM provider for subswath {} during coverage validation",
                sw.id
            ))
        })?;

        let (start, end) = provider.get_time_range();
        for burst in azimuth_timing
            .burst_timing
            .iter()
            .filter(|b| b.subswath_id == sw.id)
        {
            let line_count = burst
                .last_line_merged
                .saturating_sub(burst.first_line_merged)
                .saturating_add(1);
            let expected_end = burst.azimuth_time_start
                + (line_count.saturating_sub(1)) as f64 * burst.azimuth_time_interval;

            if start.value() > burst.azimuth_time_start + AZ_TIME_EPS
                || end.value() + AZ_TIME_EPS < expected_end
            {
                return Err(SarError::Processing(format!(
                    "DC/FM provider for {} does not cover burst [{:.6}, {:.6}] s (expected_end {:.6} from lines={}, dt={:.6}); provider [{:.6}, {:.6}]",
                    sw.id,
                    burst.azimuth_time_start,
                    burst.azimuth_time_end,
                    expected_end,
                    line_count,
                    burst.azimuth_time_interval,
                    start.value(),
                    end.value()
                )));
            }
        }
    }

    Ok(())
}

/// Validate that IW subswaths are on a unified range/azimuth grid.
///
/// Reports misalignment diagnostics to detect grid inconsistencies.
/// CRITICAL: Helps identify when ESA annotation coordinates are incorrect.
pub(crate) fn validate_grid_alignment(subswaths: &[SubSwath]) -> SarResult<()> {
    if subswaths.len() < 2 {
        return Ok(()); // Single swath, no alignment needed
    }

    log::info!("🔍 GRID ALIGNMENT VALIDATION");
    log::info!("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━");

    const C: f64 = crate::constants::physical::SPEED_OF_LIGHT_M_S;

    // Use first swath as reference
    let ref_swath = &subswaths[0];
    let ref_tau0 = ref_swath.slant_range_time;
    let ref_range_spacing = ref_swath.range_pixel_spacing;

    let mut max_range_misalignment: f64 = 0.0;
    let mut max_azimuth_misalignment: f64 = 0.0;

    for (i, swath) in subswaths.iter().enumerate().skip(1) {
        // Calculate expected range offset in pixels
        // Δn_r = round((τ₀,ᵢ - τ₀,ref) × c/2 / pixel_spacing)
        let delta_tau = swath.slant_range_time - ref_tau0;
        let delta_range_m = delta_tau * (C / 2.0);
        let expected_delta_n_r = (delta_range_m / ref_range_spacing).round();

        // Calculate actual offset from annotation
        let actual_delta_n_r =
            swath.first_sample_global as f64 - ref_swath.first_sample_global as f64;

        let range_misalignment = (actual_delta_n_r - expected_delta_n_r).abs();
        max_range_misalignment = max_range_misalignment.max(range_misalignment);

        // Calculate azimuth misalignment (if timing available)
        if let (Some(_ref_ati), Some(_swath_ati)) =
            (ref_swath.azimuth_time_interval, swath.azimuth_time_interval)
        {
            // In TOPS IW mode, all subswaths should share same azimuth grid
            // So expected offset should be close to zero (modulo burst alignment)
            let expected_delta_n_a = 0.0_f64; // Simplified - proper calculation needs burst timing
            let actual_delta_n_a =
                swath.first_line_global as f64 - ref_swath.first_line_global as f64;

            let azimuth_misalignment = (actual_delta_n_a - expected_delta_n_a).abs();
            max_azimuth_misalignment = f64::max(max_azimuth_misalignment, azimuth_misalignment);
        }

        log::info!("📊 {} → {}:", ref_swath.id, swath.id);
        log::info!(
            "   Slant range time: {:.9} s → {:.9} s (Δτ={:.6} ms)",
            ref_tau0,
            swath.slant_range_time,
            delta_tau * 1000.0
        );
        log::info!(
            "   Range offset: expected={:.1}, actual={:.1}, Δ={:.2} samples",
            expected_delta_n_r,
            actual_delta_n_r,
            range_misalignment
        );

        // Quality assessment
        if range_misalignment > 2.0 {
            log::warn!("   ⚠️  RANGE MISALIGNMENT > 2 samples - may cause visible seams!");
        } else if range_misalignment > 0.5 {
            log::warn!(
                "   ⚠️  Range misalignment {:.2} samples detected",
                range_misalignment
            );
        } else {
            log::info!("   ✅ Range grid aligned (Δ < 0.5 samples)");
        }

        // Check for critical misalignment in first_sample_global
        if swath.first_sample_global == 0 && ref_swath.first_sample_global == 0 {
            log::warn!("   ⚠️  Both subswaths have first_sample_global=0 - annotation may be missing grid info");
        }

        // Suppress unused variable warning
        let _ = i;
    }

    log::info!("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━");
    log::info!("📈 Grid Alignment Summary:");
    log::info!(
        "   Max range misalignment: {:.2} samples",
        max_range_misalignment
    );
    log::info!(
        "   Max azimuth misalignment: {:.2} lines",
        max_azimuth_misalignment
    );

    // Fail if critical misalignment detected
    if max_range_misalignment > 5.0 {
        return Err(SarError::Processing(
            format!("CRITICAL: Range grid misalignment {:.1} samples exceeds threshold (5.0) - IW merge will produce artifacts", 
                max_range_misalignment)
        ));
    }

    if max_range_misalignment > 1.0 {
        log::warn!(
            "⚠️  Range misalignment {:.2} samples may cause minor radiometric discontinuities",
            max_range_misalignment
        );
    } else {
        log::info!("✅ Grid alignment validated - proceeding with merge");
    }

    Ok(())
}

impl TopsarMerge {
    /// Validate radiometric consistency in overlap regions between subswaths.
    ///
    /// Samples pixels in overlap regions and compares dB differences to detect
    /// calibration or processing inconsistencies that could cause visible seams.
    pub(crate) fn validate_radiometric_consistency(
        &self,
        subswath_data: &HashMap<String, SarRealImage>,
    ) -> SarResult<()> {
        let subswath_lookup: HashMap<&str, &SubSwath> =
            self.subswaths.iter().map(|s| (s.id.as_str(), s)).collect();

        for overlap in &self.overlap_regions {
            if let (Some(swath1_data), Some(swath2_data)) = (
                subswath_data.get(&overlap.swath1_id),
                subswath_data.get(&overlap.swath2_id),
            ) {
                let Some(swath1_meta) = subswath_lookup.get(overlap.swath1_id.as_str()) else {
                    log::warn!(
                        "⚠️  Overlap references unknown subswath {}",
                        overlap.swath1_id
                    );
                    continue;
                };
                let Some(swath2_meta) = subswath_lookup.get(overlap.swath2_id.as_str()) else {
                    log::warn!(
                        "⚠️  Overlap references unknown subswath {}",
                        overlap.swath2_id
                    );
                    continue;
                };

                let mut differences = Vec::new();

                log::debug!(
                    "Overlap {}-{}: azimuth {}..{}, swath1 {}..{}, swath2 {}..{}",
                    overlap.swath1_id,
                    overlap.swath2_id,
                    overlap.azimuth_start,
                    overlap.azimuth_end,
                    overlap.swath1_range_start,
                    overlap.swath1_range_end,
                    overlap.swath2_range_start,
                    overlap.swath2_range_end
                );

                // Use safe overlap range validation
                let azimuth_overlap = overlap_range(
                    (overlap.azimuth_start, overlap.azimuth_end),
                    (overlap.azimuth_start, overlap.azimuth_end),
                );
                let swath1_range_overlap = overlap_range(
                    (overlap.swath1_range_start, overlap.swath1_range_end),
                    (overlap.swath1_range_start, overlap.swath1_range_end),
                );
                let swath2_range_overlap = overlap_range(
                    (overlap.swath2_range_start, overlap.swath2_range_end),
                    (overlap.swath2_range_start, overlap.swath2_range_end),
                );

                if azimuth_overlap.is_none()
                    || swath1_range_overlap.is_none()
                    || swath2_range_overlap.is_none()
                {
                    log::warn!(
                        "⚠️ Overlap region {}-{} has invalid bounds (azimuth {}..{}, swath1 {}..{}, swath2 {}..{}); skipping radiometric consistency check",
                        overlap.swath1_id,
                        overlap.swath2_id,
                        overlap.azimuth_start,
                        overlap.azimuth_end,
                        overlap.swath1_range_start,
                        overlap.swath1_range_end,
                        overlap.swath2_range_start,
                        overlap.swath2_range_end
                    );
                    continue;
                }

                // Sample pixels in overlap region for comparison with safe iteration
                let sample_step = OVERLAP_SAMPLE_STEP_RADIO; // Sample every Nth pixel for efficiency

                // Use safe range bounds that were already validated
                let (az_start, az_end) = azimuth_overlap.unwrap();
                let (sw1_start, sw1_end) = swath1_range_overlap.unwrap();
                let (sw2_start, sw2_end) = swath2_range_overlap.unwrap();

                for row in (az_start..az_end).step_by(sample_step) {
                    for col1 in (sw1_start..sw1_end).step_by(sample_step) {
                        let Some(col1_offset) = col1.checked_sub(sw1_start) else {
                            continue;
                        };
                        let Some(col2) = sw2_start.checked_add(col1_offset) else {
                            continue;
                        };

                        // Ensure col2 is within valid range for swath2
                        if col2 >= sw2_end {
                            continue;
                        }

                        // EXPERT FIX: Convert global coordinates to local coordinates for each swath
                        let Some(row_local_1) = row.checked_sub(swath1_meta.first_line_global)
                        else {
                            continue;
                        };
                        let Some(row_local_2) = row.checked_sub(swath2_meta.first_line_global)
                        else {
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

                    // Compute calibration ratio (linear scale) for diagnostic
                    let mean_diff_db = differences.iter().sum::<f32>() / differences.len() as f32;
                    let calib_ratio = 10.0_f32.powf(mean_diff_db / 10.0);

                    log::info!(
                        "📊 Radiometric consistency {}-{}: {:.3} dB median difference, calibration ratio={:.4}",
                        overlap.swath1_id,
                        overlap.swath2_id,
                        median_diff,
                        calib_ratio
                    );

                    // Detailed diagnostic for SNAP vs SARdine comparison
                    if (calib_ratio - 1.0).abs() > 0.02 {
                        log::warn!(
                            "⚠️  IW calibration mismatch detected: {}/{} power ratio = {:.4} (ideal: 1.0 ± 0.02)\n\
                             With cosine feathering, this creates gradient stripes across the overlap.\n\
                             Consider SARDINE_MERGE_MODE=sharp to use SNAP-style midpoint cutoff.",
                            overlap.swath1_id,
                            overlap.swath2_id,
                            calib_ratio
                        );
                    }

                    if median_diff > RADIOMETRIC_CONSISTENCY_DB_THRESHOLD {
                        log::warn!(
                            "⚠️  Radiometric inconsistency detected: {:.3} dB > {:.1} dB threshold",
                            median_diff,
                            RADIOMETRIC_CONSISTENCY_DB_THRESHOLD
                        );
                        // Could apply gain correction here if needed
                    }
                }
            }
        }

        Ok(())
    }
}
