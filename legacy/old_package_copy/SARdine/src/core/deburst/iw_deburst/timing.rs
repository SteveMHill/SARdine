//! Timing utilities for IW TOPSAR deburst processing
//!
//! This module provides per-line azimuth timing calculations used during
//! Doppler centroid and FM rate polynomial evaluation.

/// Holds per-line timing derived from annotation (preferred over spacing/velocity)
#[derive(Clone, Copy)]
pub(crate) struct LineTiming {
    /// seconds relative to polynomial reference time (annotation-aware)
    pub(crate) t_az: f64,
}

/// Build per-line timings with correct annotation reference time
///
/// CRITICAL: Doppler/FM polynomials in annotation are referenced to a specific time origin
/// (usually orbit epoch or product reference time). Evaluating with incorrect t causes phase errors.
pub(crate) fn build_line_timing_with_offset(
    lines: usize,
    az_time_interval_s: f64,
    time_offset_s: f64,
) -> Vec<LineTiming> {
    let center = (lines as f64 - 1.0) * 0.5;
    (0..lines)
        .map(|l| {
            let t_burst_relative = (l as f64 - center) * az_time_interval_s;
            LineTiming {
                t_az: t_burst_relative + time_offset_s,
            }
        })
        .collect()
}

// REMOVED: build_line_timing (legacy, use build_line_timing_with_offset for proper reference time handling)
