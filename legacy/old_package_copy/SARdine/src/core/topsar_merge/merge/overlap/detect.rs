#![allow(dead_code, unused_variables)]
#![allow(dead_code)]
//! Overlap detection and region calculation.
//!
//! Detects overlapping regions between adjacent IW subswaths using
//! both global grid coordinates and slant range fallback.

use crate::types::{SarError, SarResult, SubSwath};

use super::super::types::{BlendingMethod, OverlapQuality, OverlapRegion};
use super::weights::{create_complementary_cosine_weights, create_sharp_cutoff_weights};
use super::super::assert_subswath_geometry;

/// Safe overlap range calculation - returns None if no overlap exists.
///
/// Given two ranges `a = (start, end)` and `b = (start, end)`, computes
/// the intersection. Returns `None` if ranges do not overlap.
#[inline]
pub fn overlap_range(a: (usize, usize), b: (usize, usize)) -> Option<(usize, usize)> {
    let start = a.0.max(b.0);
    let end = a.1.min(b.1);
    if end > start {
        Some((start, end))
    } else {
        None
    }
}

/// Detect overlap regions between adjacent IW subswaths.
///
/// Processes IW1-IW2 and IW2-IW3 overlaps using global grid coordinates
/// with slant-range fallback for robustness.
pub fn detect_subswath_overlaps(
    subswaths: &[SubSwath],
    feather_width: usize,
) -> SarResult<Vec<OverlapRegion>> {
    log::info!("🎯 Detecting subswath overlaps using global grid analysis");
    let mut overlaps = Vec::new();

    #[cfg(debug_assertions)]
    for sw in subswaths {
        assert_subswath_geometry(sw);
    }

    // Find IW subswaths and sort by ID for proper ordering
    let mut iw_swaths: Vec<_> = subswaths
        .iter()
        .filter(|sw| sw.id.starts_with("IW"))
        .collect();
    iw_swaths.sort_by(|a, b| a.id.cmp(&b.id));

    // Detect IW1-IW2 overlap
    if let (Some(iw1), Some(iw2)) = (
        iw_swaths.iter().find(|sw| sw.id == "IW1"),
        iw_swaths.iter().find(|sw| sw.id == "IW2"),
    ) {
        if let Some(overlap) = calculate_overlap_region(iw1, iw2, feather_width)? {
            log::info!(
                "✓ IW1-IW2 overlap: range [{}-{}] × [{}-{}], azimuth [{}-{}]",
                overlap.swath1_range_start,
                overlap.swath1_range_end,
                overlap.swath2_range_start,
                overlap.swath2_range_end,
                overlap.azimuth_start,
                overlap.azimuth_end
            );
            overlaps.push(overlap);
        }
    }

    // Detect IW2-IW3 overlap
    if let (Some(iw2), Some(iw3)) = (
        iw_swaths.iter().find(|sw| sw.id == "IW2"),
        iw_swaths.iter().find(|sw| sw.id == "IW3"),
    ) {
        if let Some(overlap) = calculate_overlap_region(iw2, iw3, feather_width)? {
            log::info!(
                "✓ IW2-IW3 overlap: range [{}-{}] × [{}-{}], azimuth [{}-{}]",
                overlap.swath1_range_start,
                overlap.swath1_range_end,
                overlap.swath2_range_start,
                overlap.swath2_range_end,
                overlap.azimuth_start,
                overlap.azimuth_end
            );
            overlaps.push(overlap);
        }
    }

    log::info!("📊 Detected {} overlap regions", overlaps.len());
    Ok(overlaps)
}

/// Detect overlap regions with explicit blending method selection.
///
/// This version allows choosing between cosine feathering (default) and
/// SNAP-compatible sharp cutoff (midpoint selection).
///
/// # Arguments
///
/// * `subswaths` - List of subswaths to analyze
/// * `feather_width` - Width of feather zone (ignored for SnapCompatible)
/// * `blending_method` - Blending method to use for overlap weights
pub fn detect_subswath_overlaps_with_blending(
    subswaths: &[SubSwath],
    feather_width: usize,
    blending_method: &BlendingMethod,
) -> SarResult<Vec<OverlapRegion>> {
    let use_sharp_cutoff = matches!(blending_method, BlendingMethod::SnapCompatible);

    if use_sharp_cutoff {
        log::info!(
            "🎯 Detecting subswath overlaps using SNAP-compatible sharp cutoff (no blending)"
        );
    } else {
        log::info!("🎯 Detecting subswath overlaps using cosine feathering");
    }

    let mut overlaps = Vec::new();

    #[cfg(debug_assertions)]
    for sw in subswaths {
        assert_subswath_geometry(sw);
    }

    // Find IW subswaths and sort by ID for proper ordering
    let mut iw_swaths: Vec<_> = subswaths
        .iter()
        .filter(|sw| sw.id.starts_with("IW"))
        .collect();
    iw_swaths.sort_by(|a, b| a.id.cmp(&b.id));

    // Detect IW1-IW2 overlap
    if let (Some(iw1), Some(iw2)) = (
        iw_swaths.iter().find(|sw| sw.id == "IW1"),
        iw_swaths.iter().find(|sw| sw.id == "IW2"),
    ) {
        if let Some(overlap) =
            calculate_overlap_region_with_blending(iw1, iw2, feather_width, use_sharp_cutoff)?
        {
            log::info!(
                "✓ IW1-IW2 overlap: range [{}-{}] × [{}-{}], azimuth [{}-{}], mode={}",
                overlap.swath1_range_start,
                overlap.swath1_range_end,
                overlap.swath2_range_start,
                overlap.swath2_range_end,
                overlap.azimuth_start,
                overlap.azimuth_end,
                if use_sharp_cutoff {
                    "sharp_cutoff"
                } else {
                    "cosine"
                }
            );
            overlaps.push(overlap);
        }
    }

    // Detect IW2-IW3 overlap
    if let (Some(iw2), Some(iw3)) = (
        iw_swaths.iter().find(|sw| sw.id == "IW2"),
        iw_swaths.iter().find(|sw| sw.id == "IW3"),
    ) {
        if let Some(overlap) =
            calculate_overlap_region_with_blending(iw2, iw3, feather_width, use_sharp_cutoff)?
        {
            log::info!(
                "✓ IW2-IW3 overlap: range [{}-{}] × [{}-{}], azimuth [{}-{}], mode={}",
                overlap.swath1_range_start,
                overlap.swath1_range_end,
                overlap.swath2_range_start,
                overlap.swath2_range_end,
                overlap.azimuth_start,
                overlap.azimuth_end,
                if use_sharp_cutoff {
                    "sharp_cutoff"
                } else {
                    "cosine"
                }
            );
            overlaps.push(overlap);
        }
    }

    log::info!("📊 Detected {} overlap regions", overlaps.len());
    Ok(overlaps)
}

/// Calculate overlap region between two adjacent subswaths.
///
/// Uses global grid first, with slant-range fallback for robustness.
pub fn calculate_overlap_region(
    swath1: &SubSwath,
    swath2: &SubSwath,
    feather_width: usize,
) -> SarResult<Option<OverlapRegion>> {
    #[cfg(debug_assertions)]
    {
        assert_subswath_geometry(swath1);
        assert_subswath_geometry(swath2);
    }

    // 1) Try fast path using already-aligned global sample coordinates.
    let global_start = swath1.first_sample_global.max(swath2.first_sample_global);
    let global_end = swath1.last_sample_global.min(swath2.last_sample_global);

    // Check for no overlap (when global_end <= global_start)
    if global_end <= global_start {
        log::debug!(
            "No overlap between {} and {} (no global range intersection)",
            swath1.id,
            swath2.id
        );
        return Ok(None);
    }

    let (swath1_final_start, swath1_final_end, swath2_final_start, swath2_final_end) = (
        global_start.saturating_sub(swath1.first_sample_global),
        global_end.saturating_sub(swath1.first_sample_global),
        global_start.saturating_sub(swath2.first_sample_global),
        global_end.saturating_sub(swath2.first_sample_global),
    );

    // Calculate common width
    let width1 = swath1_final_end.saturating_sub(swath1_final_start);
    let width2 = swath2_final_end.saturating_sub(swath2_final_start);
    let common_width = width1.min(width2);

    if common_width == 0 {
        log::debug!(
            "No overlap between {} and {} (zero width)",
            swath1.id,
            swath2.id
        );
        return Ok(None);
    }

    // Adjust end indices to use common width
    let swath1_final_end = swath1_final_start
        .checked_add(common_width)
        .ok_or_else(|| SarError::Processing("Swath1 overlap range overflow".to_string()))?;
    let swath2_final_end = swath2_final_start
        .checked_add(common_width)
        .ok_or_else(|| SarError::Processing("Swath2 overlap range overflow".to_string()))?;

    // Azimuth overlap using GLOBAL coordinates
    // FIX: Use UNION (min/max) instead of INTERSECTION (max/min) to ensure ALL rows
    // where both swaths have range overlap get proper blending weights.
    // Previously, rows outside the strict azimuth intersection got Constant(1.0) weight
    // from BOTH swaths, causing doubled intensity (bright stripes).
    let azimuth_start = swath1.first_line_global.min(swath2.first_line_global);
    let azimuth_end = swath1.last_line_global.max(swath2.last_line_global);

    if azimuth_end <= azimuth_start {
        log::debug!(
            "No azimuth overlap between {} and {} ({}..{})",
            swath1.id,
            swath2.id,
            azimuth_start,
            azimuth_end
        );
        return Ok(None);
    }

    let azimuth_height = azimuth_end.saturating_sub(azimuth_start);

    // Create weights with correct monotonic direction
    // FIX: Always use full overlap width for feathering to eliminate hard cuts
    // The feather_width parameter is kept for API compatibility but full-width feathering
    // is used internally for seamless blending (matches SNAP behavior)
    let weights = create_complementary_cosine_weights(common_width, azimuth_height, common_width);

    let quality_metrics = OverlapQuality::default();

    Ok(Some(OverlapRegion {
        swath1_id: swath1.id.clone(),
        swath2_id: swath2.id.clone(),
        swath1_range_start: swath1_final_start,
        swath1_range_end: swath1_final_end,
        swath2_range_start: swath2_final_start,
        swath2_range_end: swath2_final_end,
        azimuth_start,
        azimuth_end,
        weights,
        quality_metrics,
    }))
}

/// Calculate overlap region with selectable blending method.
///
/// # Arguments
///
/// * `swath1` - First (left) subswath
/// * `swath2` - Second (right) subswath  
/// * `feather_width` - Width of feather zone (ignored for sharp cutoff)
/// * `use_sharp_cutoff` - If true, use SNAP-style midpoint selection instead of cosine
pub fn calculate_overlap_region_with_blending(
    swath1: &SubSwath,
    swath2: &SubSwath,
    feather_width: usize,
    use_sharp_cutoff: bool,
) -> SarResult<Option<OverlapRegion>> {
    #[cfg(debug_assertions)]
    {
        assert_subswath_geometry(swath1);
        assert_subswath_geometry(swath2);
    }

    // 1) Try fast path using already-aligned global sample coordinates.
    let global_start = swath1.first_sample_global.max(swath2.first_sample_global);
    let global_end = swath1.last_sample_global.min(swath2.last_sample_global);

    // Check for no overlap (when global_end <= global_start)
    if global_end <= global_start {
        log::debug!(
            "No overlap between {} and {} (no global range intersection)",
            swath1.id,
            swath2.id
        );
        return Ok(None);
    }

    let (swath1_final_start, swath1_final_end, swath2_final_start, swath2_final_end) = (
        global_start.saturating_sub(swath1.first_sample_global),
        global_end.saturating_sub(swath1.first_sample_global),
        global_start.saturating_sub(swath2.first_sample_global),
        global_end.saturating_sub(swath2.first_sample_global),
    );

    // Calculate common width
    let width1 = swath1_final_end.saturating_sub(swath1_final_start);
    let width2 = swath2_final_end.saturating_sub(swath2_final_start);
    let common_width = width1.min(width2);

    if common_width == 0 {
        log::debug!(
            "No overlap between {} and {} (zero width)",
            swath1.id,
            swath2.id
        );
        return Ok(None);
    }

    // Adjust end indices to use common width
    let swath1_final_end = swath1_final_start
        .checked_add(common_width)
        .ok_or_else(|| SarError::Processing("Swath1 overlap range overflow".to_string()))?;
    let swath2_final_end = swath2_final_start
        .checked_add(common_width)
        .ok_or_else(|| SarError::Processing("Swath2 overlap range overflow".to_string()))?;

    // Azimuth overlap using GLOBAL coordinates
    // FIX: Use UNION (min/max) instead of INTERSECTION (max/min) to ensure ALL rows
    // where both swaths have range overlap get proper blending weights.
    let azimuth_start = swath1.first_line_global.min(swath2.first_line_global);
    let azimuth_end = swath1.last_line_global.max(swath2.last_line_global);

    if azimuth_end <= azimuth_start {
        log::debug!(
            "No azimuth overlap between {} and {} ({}..{})",
            swath1.id,
            swath2.id,
            azimuth_start,
            azimuth_end
        );
        return Ok(None);
    }

    let azimuth_height = azimuth_end.saturating_sub(azimuth_start);

    // Create weights based on blending method selection
    let weights = if use_sharp_cutoff {
        // SNAP-style: sharp midpoint cutoff - no gradient, no calibration mismatch amplification
        create_sharp_cutoff_weights(common_width, azimuth_height)
    } else {
        // Default: full-width cosine feathering
        create_complementary_cosine_weights(common_width, azimuth_height, common_width)
    };

    let quality_metrics = OverlapQuality::default();

    Ok(Some(OverlapRegion {
        swath1_id: swath1.id.clone(),
        swath2_id: swath2.id.clone(),
        swath1_range_start: swath1_final_start,
        swath1_range_end: swath1_final_end,
        swath2_range_start: swath2_final_start,
        swath2_range_end: swath2_final_end,
        azimuth_start,
        azimuth_end,
        weights,
        quality_metrics,
    }))
}

/// Audit overlap geometry to surface zero-width or very thin regions early.
pub fn audit_overlap_regions(overlaps: &[OverlapRegion]) -> SarResult<()> {
    if overlaps.is_empty() {
        log::warn!("⚠️  No overlap regions detected; merged output will have hard seams");
        return Ok(());
    }

    for ov in overlaps {
        if ov.swath1_id == ov.swath2_id {
            return Err(SarError::Processing(format!(
                "CRITICAL: overlap lists identical swaths {} vs {}",
                ov.swath1_id, ov.swath2_id
            )));
        }

        let w1 = ov.swath1_range_end.saturating_sub(ov.swath1_range_start);
        let w2 = ov.swath2_range_end.saturating_sub(ov.swath2_range_start);
        let az = ov.azimuth_end.saturating_sub(ov.azimuth_start);

        if w1 == 0 || w2 == 0 || az == 0 {
            return Err(SarError::Processing(format!(
                "CRITICAL: overlap {}-{} has zero extent (w1={}, w2={}, az={})",
                ov.swath1_id, ov.swath2_id, w1, w2, az
            )));
        }

        if w1.min(w2) < 8 {
            log::warn!(
                "⚠️  Thin overlap {}-{}: width1={} width2={} az_height={}; seams likely",
                ov.swath1_id,
                ov.swath2_id,
                w1,
                w2,
                az
            );
        }
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use ndarray::Array2;

    fn make_swath(id: &str, first_sample: usize, last_sample: usize) -> SubSwath {
        SubSwath {
            id: id.to_string(),
            burst_count: 1,
            lines_per_burst: 100,
            range_samples: last_sample - first_sample,
            azimuth_samples: 100,
            first_line_global: 0,
            last_line_global: 100,
            first_sample_global: first_sample,
            last_sample_global: last_sample,
            full_range_samples: last_sample - first_sample,
            valid_first_line: Some(0),
            valid_last_line: Some(100),
            valid_first_sample: Some(first_sample),
            valid_last_sample: Some(last_sample),
            range_pixel_spacing: 2.3,
            azimuth_pixel_spacing: 14.0,
            slant_range_time: 0.005,
            burst_duration: 2.7,
            near_range_m: 800000.0,
            prf_hz: Some(1680.0),
            dc_polynomial: Some(vec![0.0, 0.0, 0.0]),
            azimuth_time_interval: Some(0.000595),
            dc_polynomial_t0: Some(0.0),
            fm_rate_estimates: None,
        }
    }

    #[test]
    fn overlap_detection_simple() {
        // IW1: samples 0-100, IW2: samples 80-180 → overlap at 80-100 (20 samples)
        let iw1 = make_swath("IW1", 0, 100);
        let iw2 = make_swath("IW2", 80, 180);

        let overlap = calculate_overlap_region(&iw1, &iw2, 20)
            .unwrap()
            .expect("should detect overlap");

        assert_eq!(overlap.swath1_range_start, 80); // local to IW1: 80
        assert_eq!(overlap.swath1_range_end, 100); // local to IW1: 100
        assert_eq!(overlap.swath2_range_start, 0); // local to IW2: 0
        assert_eq!(overlap.swath2_range_end, 20); // local to IW2: 20

        // Weights should be correct dimension
        assert_eq!(overlap.weights.nrows(), 100); // azimuth extent
        assert_eq!(overlap.weights.ncols(), 20); // range extent
    }

    #[test]
    fn no_overlap_when_separated() {
        let iw1 = make_swath("IW1", 0, 100);
        let iw2 = make_swath("IW2", 200, 300);

        let overlap = calculate_overlap_region(&iw1, &iw2, 20).unwrap();
        assert!(overlap.is_none());
    }

    #[test]
    fn audit_catches_zero_extent() {
        let bad_overlap = OverlapRegion {
            swath1_id: "IW1".to_string(),
            swath2_id: "IW2".to_string(),
            swath1_range_start: 50,
            swath1_range_end: 50, // zero width
            swath2_range_start: 0,
            swath2_range_end: 10,
            azimuth_start: 0,
            azimuth_end: 100,
            weights: Array2::zeros((100, 0)),
            quality_metrics: OverlapQuality::default(),
        };

        let result = audit_overlap_regions(&[bad_overlap]);
        assert!(result.is_err());
    }
}
