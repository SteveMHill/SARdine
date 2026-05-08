#![allow(dead_code)]
//! Weighting and blending functions for TOPSAR burst mosaicking
//!
//! This module provides the windowing and weight computation functions used
//! for seamless burst blending during TOPSAR deburst processing. The core
//! algorithm uses complementary cosine-squared (cos²) weights to ensure
//! perfect energy preservation across burst overlap regions.
//!
//! # References
//! - ESA S1-TN-MDA-52-7445: "TOPSAR Debursting" Section 3.2 (complementary blending)

/// Respect valid samples per line (first/last) with floor-safe indexing
///
/// Note: Handles edge case where last_valid[line] = -1 by clamping to 0.
/// This relies on the (a <= b) check to return empty window (0,0).
/// If you see 1-pixel slivers, check raw annotation values here.
pub(crate) fn valid_window(
    line_in_burst: usize,
    first_valid: &[i32],
    last_valid: &[i32],
    width: usize,
) -> (usize, usize) {
    if line_in_burst >= first_valid.len() || line_in_burst >= last_valid.len() {
        return (0, 0);
    }

    let first_raw = first_valid[line_in_burst];
    let last_raw = last_valid[line_in_burst];

    // CRITICAL FIX: -1 means NO valid samples (not "up to pixel 0")
    // Return empty window immediately to prevent spurious 1-pixel slivers
    // Check BOTH first_raw and last_raw for -1 (invalid marker)
    if first_raw < 0 || last_raw < 0 {
        return (0, 0);
    }

    // Now clamp positive values
    let mut a = first_raw.max(0) as usize;
    let mut b = last_raw as usize; // last_raw is >= 0 here

    a = a.min(width.saturating_sub(1));
    b = b.min(width.saturating_sub(1));

    if a <= b {
        (a, b + 1)
    } else {
        // Log suspicious cases for debugging annotation issues
        if first_raw >= 0 && last_raw >= 0 && first_raw <= last_raw {
            log::debug!(
                "Valid window inverted after clamping: line={}, first_raw={}, last_raw={}, width={}, a={}, b={}",
                line_in_burst, first_raw, last_raw, width, a, b
            );
        }
        (0, 0)
    }
}

/// Complementary overlap blending (pairwise cos²) with hit-count mask
#[inline]
pub(crate) fn w_cos2(u01: f32) -> f32 {
    let x = std::f32::consts::FRAC_PI_2 * u01.clamp(0.0, 1.0);
    let c = x.cos();
    c * c
}

/// Compute pairwise complementary weights for exact overlap region
/// Ensures perfect energy preservation: w₁ + w₂ = 1.0 at every overlap pixel
///
/// This enforces strict complementarity for seamless blending without seams.
///
/// # Arguments
/// * `overlap_len` - Total number of lines in overlap region
/// * `line_in_overlap` - Current line index within overlap (0 = start of overlap)
///
/// # Returns
/// (w_current, w_next) where w_current + w_next = 1.0 (guaranteed)
///
/// # References
/// - ESA S1-TN-MDA-52-7445: "TOPSAR Debursting" Section 3.2 (complementary blending)
///
/// **Enhancement #2 (2025-10-04):** Pairwise weight enforcement for perfect energy conservation
pub(crate) fn compute_pairwise_weights(overlap_len: usize, line_in_overlap: usize) -> (f32, f32) {
    if overlap_len < 2 {
        return (1.0, 0.0); // No overlap
    }

    // Normalized position in overlap: 0.0 at start, 1.0 at end
    let u = if overlap_len > 1 {
        line_in_overlap as f32 / (overlap_len - 1) as f32
    } else {
        0.5
    };

    let w_current = w_cos2(u); // Current burst weight (fades out)
    let w_next = 1.0 - w_current; // Next burst weight (fades in)

    // Safety check: ensure perfect complementarity
    debug_assert!(
        (w_current + w_next - 1.0).abs() < 1e-6,
        "Weights not complementary: {} + {} = {} (line {} of {})",
        w_current,
        w_next,
        w_current + w_next,
        line_in_overlap,
        overlap_len
    );

    (w_current, w_next)
}

/// Returns line weight for burst k, and (optionally) 1-w for k+1
pub(crate) fn overlap_weight(
    line_in_burst: usize,
    lines: usize,
    blend_lines: usize,
    is_prev: bool,
) -> f32 {
    if blend_lines == 0 {
        return 1.0;
    }
    let head = line_in_burst;
    let tail = lines.saturating_sub(1).saturating_sub(line_in_burst);
    if head < blend_lines {
        let u = head as f32 / blend_lines as f32;
        if is_prev {
            w_cos2(u)
        } else {
            1.0 - w_cos2(u)
        }
    } else if tail < blend_lines {
        let u = tail as f32 / blend_lines as f32;
        if is_prev {
            1.0 - w_cos2(u)
        } else {
            w_cos2(u)
        }
    } else {
        1.0
    }
}
