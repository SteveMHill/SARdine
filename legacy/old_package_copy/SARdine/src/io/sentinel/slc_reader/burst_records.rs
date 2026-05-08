//! Burst record construction with correct TOPS geometry
//!
//! Fixes the following issues:
//! - A2: TOPS burst range geometry - bursts tile in azimuth, NOT range
//! - A3: Valid sample bounds use exclusive upper bound
//! - A1: Timing uses orbit-reference-relative seconds

use super::time_utils::{parse_iso8601_to_datetime_utc, seconds_since_epoch};
use crate::types::{BurstRecord, SubSwath};
use chrono::{DateTime, Utc};

/// Indicates which time domain was used for burst timing
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BurstTimingDomain {
    /// Time computed relative to orbit reference epoch (preferred)
    OrbitRelative,
    /// Time from annotation azimuth_anx_time (fallback when orbit epoch unavailable)
    AzimuthAnxTime,
    /// No timing information available
    Unavailable,
}

/// Build burst records for a subswath with CORRECT TOPS geometry
///
/// **Critical fixes applied:**
/// 1. **Range window is CONSTANT per subswath** - bursts tile in azimuth only
/// 2. **Valid sample bounds are EXCLUSIVE** - last_valid_sample_exclusive = lastValidSample + 1
/// 3. **Timing uses orbit-reference-relative seconds** when epoch is available
///
/// # Arguments
/// * `annotation` - The parsed annotation ProductRoot
/// * `subswath_id` - The subswath identifier (e.g., "IW1", "IW2", "IW3")
/// * `subswath` - The SubSwath geometry information
/// * `orbit_ref_epoch` - Optional orbit reference epoch for timing conversion
///
/// # Returns
/// Vector of BurstRecord with corrected geometry
pub fn build_burst_records_for_subswath_fixed(
    annotation: &crate::io::annotation::ProductRoot,
    subswath_id: &str,
    subswath: &SubSwath,
    orbit_ref_epoch: Option<DateTime<Utc>>,
) -> Vec<BurstRecord> {
    let swath_timing = match &annotation.swath_timing {
        Some(timing) => timing,
        None => return Vec::new(),
    };

    let bursts = match swath_timing
        .burst_list
        .as_ref()
        .and_then(|list| list.bursts.as_ref())
    {
        Some(list) if !list.is_empty() => list,
        _ => return Vec::new(),
    };

    // Calculate total lines in subswath
    let total_lines = subswath
        .last_line_global
        .saturating_sub(subswath.first_line_global)
        .saturating_add(1);

    if total_lines == 0 {
        return Vec::new();
    }

    let total_bursts = bursts.len().max(1);

    // Lines per burst from annotation or computed
    let lines_per_burst = swath_timing
        .lines_per_burst
        .map(|v| v as usize)
        .filter(|&v| v > 0)
        .unwrap_or_else(|| std::cmp::max(1, total_lines / total_bursts));

    // =========================================================================
    // FIX A2: RANGE WINDOW IS CONSTANT FOR ALL BURSTS IN A SUBSWATH
    // In TOPS/IW SLC, bursts tile in AZIMUTH, not range.
    // All bursts share the same range extent.
    // =========================================================================
    let constant_first_sample = subswath.first_sample_global;
    let constant_last_sample = subswath.last_sample_global;

    log::debug!(
        "📐 {} burst geometry: {} bursts, {} lines/burst, range=[{}, {}] (constant)",
        subswath_id,
        total_bursts,
        lines_per_burst,
        constant_first_sample,
        constant_last_sample
    );

    let strict = std::env::var("SARDINE_STRICT")
        .map(|v| v == "1" || v.to_lowercase() == "true")
        .unwrap_or(false);

    if strict && orbit_ref_epoch.is_none() {
        panic!(
            "SARDINE_STRICT=1 requires orbit_ref_epoch for burst timing in build_burst_records_for_subswath_fixed (subswath={})",
            subswath_id
        );
    }

    let mut current_line = subswath.first_line_global;
    let mut records = Vec::with_capacity(bursts.len());

    for (idx, burst) in bursts.iter().enumerate() {
        // Azimuth (line) window advances per burst
        let first_line = current_line;
        let mut last_line = first_line.saturating_add(lines_per_burst).saturating_sub(1);

        // Clamp to subswath bounds for last burst
        if idx == total_bursts - 1 || last_line > subswath.last_line_global {
            last_line = subswath.last_line_global;
        }
        current_line = last_line.saturating_add(1);

        // =====================================================================
        // FIX A3: VALID SAMPLE BOUNDS USE EXCLUSIVE UPPER BOUND
        // first_valid_sample: inclusive (first valid sample index)
        // last_valid_sample_exclusive: exclusive (one past last valid sample)
        // =====================================================================
        let (first_valid, last_valid_exclusive) = compute_valid_sample_bounds(
            &burst.first_valid_sample,
            &burst.last_valid_sample,
            constant_first_sample,
            constant_last_sample,
        );

        // =====================================================================
        // FIX A1: TIMING USES ORBIT-REFERENCE-RELATIVE SECONDS
        // =====================================================================
        let (azimuth_time_rel_orbit, azimuth_time_absolute, timing_domain) =
            compute_burst_timing(burst, orbit_ref_epoch);

        if strict && timing_domain != BurstTimingDomain::OrbitRelative {
            panic!(
                "SARDINE_STRICT=1 requires OrbitRelative burst timing; got {:?} for subswath={} burst_index={}",
                timing_domain,
                subswath_id,
                idx
            );
        }

        if idx == 0 {
            log::info!(
                "🕐 {} burst timing domain: {:?}",
                subswath_id,
                timing_domain
            );
        }

        records.push(BurstRecord {
            subswath_id: subswath_id.to_string(),
            burst_index: idx,
            // Azimuth bounds (advance per burst)
            first_line_global: first_line,
            last_line_global: last_line,
            // Range bounds (CONSTANT for all bursts in subswath)
            start_sample_global: constant_first_sample,
            end_sample_global: constant_last_sample,
            // Valid bounds with exclusive upper bound
            first_valid_sample: Some(first_valid),
            last_valid_sample: Some(last_valid_exclusive), // NOTE: This is now EXCLUSIVE
            first_valid_line: Some(first_line),
            last_valid_line: Some(last_line),
            // Timing
            azimuth_time_rel_orbit,
            azimuth_time_absolute,
        });
    }

    if strict && !records.is_empty() {
        let first_record = &records[0];
        let last_record = &records[records.len() - 1];

        // Azimuth tiling: bursts must cover the subswath exactly with no gaps/overlaps
        assert_eq!(
            first_record.first_line_global,
            subswath.first_line_global,
            "SARDINE_STRICT=1: first burst line ({}) must equal subswath.first_line_global ({}) for subswath {}",
            first_record.first_line_global,
            subswath.first_line_global,
            subswath_id
        );
        assert_eq!(
            last_record.last_line_global,
            subswath.last_line_global,
            "SARDINE_STRICT=1: last burst line ({}) must equal subswath.last_line_global ({}) for subswath {}",
            last_record.last_line_global,
            subswath.last_line_global,
            subswath_id
        );

        let expected_first_sample = subswath.first_sample_global;
        let expected_last_sample = subswath.last_sample_global;

        let mut prev_last_line = first_record.last_line_global;
        for (i, record) in records.iter().enumerate() {
            // Constant range window for all bursts in subswath
            assert_eq!(
                record.start_sample_global,
                expected_first_sample,
                "SARDINE_STRICT=1: burst {} in subswath {} has start_sample_global={} but expected {}",
                i,
                subswath_id,
                record.start_sample_global,
                expected_first_sample
            );
            assert_eq!(
                record.end_sample_global,
                expected_last_sample,
                "SARDINE_STRICT=1: burst {} in subswath {} has end_sample_global={} but expected {}",
                i,
                subswath_id,
                record.end_sample_global,
                expected_last_sample
            );

            // For all but the first burst, enforce exact tiling in azimuth
            if i > 0 {
                let expected_first_line = prev_last_line
                    .checked_add(1)
                    .expect("azimuth line index overflow in burst tiling check");
                assert_eq!(
                    record.first_line_global,
                    expected_first_line,
                    "SARDINE_STRICT=1: bursts must tile contiguously in azimuth; burst {} in subswath {} starts at line {} but expected {} (previous last line was {})",
                    i,
                    subswath_id,
                    record.first_line_global,
                    expected_first_line,
                    prev_last_line
                );
                prev_last_line = record.last_line_global;
            }
        }
    }

    log::info!(
        "✅ {} built {} burst records: lines={}..{}, samples={}..{} (constant)",
        subswath_id,
        records.len(),
        subswath.first_line_global,
        subswath.last_line_global,
        constant_first_sample,
        constant_last_sample
    );

    records
}

/// Compute valid sample bounds with exclusive upper bound
///
/// Returns (first_valid_sample, last_valid_sample_exclusive)
fn compute_valid_sample_bounds(
    first_valid_samples: &[i32],
    last_valid_samples: &[i32],
    subswath_first_sample: usize,
    subswath_last_sample: usize,
) -> (usize, usize) {
    // Find minimum first valid sample across all lines
    let min_first = first_valid_samples
        .iter()
        .filter(|&&v| v >= 0)
        .map(|&v| subswath_first_sample.saturating_add(v as usize))
        .min()
        .unwrap_or(subswath_first_sample);

    // Find maximum last valid sample across all lines
    let max_last_inclusive = last_valid_samples
        .iter()
        .filter(|&&v| v >= 0)
        .map(|&v| subswath_first_sample.saturating_add(v as usize))
        .max()
        .unwrap_or(subswath_last_sample);

    // Convert to exclusive upper bound
    let max_last_exclusive = max_last_inclusive.saturating_add(1);

    // Clamp to subswath range window
    let first_valid = min_first.max(subswath_first_sample);
    let last_valid_exclusive = max_last_exclusive.min(subswath_last_sample.saturating_add(1));

    (first_valid, last_valid_exclusive)
}

/// Compute burst timing with proper domain handling
///
/// Returns (azimuth_time_rel_orbit, azimuth_time_absolute, timing_domain)
fn compute_burst_timing(
    burst: &crate::io::annotation::Burst,
    orbit_ref_epoch: Option<DateTime<Utc>>,
) -> (Option<f64>, Option<f64>, BurstTimingDomain) {
    // Try to parse absolute burst time from annotation
    let burst_datetime = burst
        .azimuth_time
        .as_ref()
        .or(burst.sensing_time.as_ref())
        .and_then(|t| parse_iso8601_to_datetime_utc(t));

    // Compute absolute Unix seconds (for legacy interchange)
    let azimuth_time_absolute = burst_datetime
        .map(|dt| dt.timestamp() as f64 + (dt.timestamp_subsec_nanos() as f64) * 1e-9);

    // Compute orbit-relative time if we have both epoch and burst datetime
    if let (Some(epoch), Some(burst_dt)) = (orbit_ref_epoch, burst_datetime) {
        let rel_s = seconds_since_epoch(epoch, burst_dt);
        if rel_s.is_finite() {
            return (
                Some(rel_s),
                azimuth_time_absolute,
                BurstTimingDomain::OrbitRelative,
            );
        } else {
            log::warn!(
                "⚠️  Non-finite orbit-relative time computed; falling back to azimuth_anx_time"
            );
        }
    }

    // Fallback: use azimuth_anx_time from annotation
    if burst.azimuth_anx_time.is_finite() && burst.azimuth_anx_time > 0.0 {
        if orbit_ref_epoch.is_some() {
            log::warn!(
                "⚠️  Could not parse burst datetime; using annotation azimuth_anx_time={:.6}s",
                burst.azimuth_anx_time
            );
        }
        return (
            Some(burst.azimuth_anx_time),
            azimuth_time_absolute,
            BurstTimingDomain::AzimuthAnxTime,
        );
    }

    // No timing available
    (None, azimuth_time_absolute, BurstTimingDomain::Unavailable)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_valid_sample_bounds_exclusive() {
        // If annotation says lastValidSample = 99, we should get exclusive = 100
        let first_valid = vec![10, 12, 11];
        let last_valid = vec![95, 99, 97];

        let (first, last_excl) = compute_valid_sample_bounds(
            &first_valid,
            &last_valid,
            0,   // subswath_first_sample
            100, // subswath_last_sample
        );

        assert_eq!(first, 10, "Should use minimum first_valid");
        assert_eq!(
            last_excl, 100,
            "Should be max(last_valid) + 1 = 99 + 1 = 100"
        );
    }

    #[test]
    fn test_valid_sample_bounds_with_invalid_markers() {
        // -1 values should be filtered out
        let first_valid = vec![-1, 10, -1];
        let last_valid = vec![-1, 50, -1];

        let (first, last_excl) = compute_valid_sample_bounds(&first_valid, &last_valid, 0, 100);

        assert_eq!(first, 10);
        assert_eq!(last_excl, 51); // 50 + 1
    }

    #[test]
    fn test_valid_sample_bounds_clamped_to_subswath() {
        // Valid samples should be clamped to subswath range
        let first_valid = vec![5];
        let last_valid = vec![150]; // Beyond subswath

        let (first, last_excl) = compute_valid_sample_bounds(
            &first_valid,
            &last_valid,
            10,  // subswath_first_sample
            100, // subswath_last_sample
        );

        assert_eq!(
            first, 15,
            "Should be offset by subswath_first but clamped to >= 10"
        );
        // 5 + 10 = 15, which is >= 10, so no clamping needed for first
        // Actually first = max(5+10, 10) = max(15, 10) = 15 ✓

        assert_eq!(last_excl, 101, "Should be clamped to subswath_last + 1");
        // 150 + 10 = 160 (global), then +1 = 161, but clamped to 100+1 = 101 ✓
    }
}
