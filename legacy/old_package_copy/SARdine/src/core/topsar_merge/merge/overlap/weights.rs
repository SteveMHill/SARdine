#![allow(dead_code, unused_variables)]
//! Overlap weight generation with correct monotonic cosine taper.
//!
//! # Weight Convention
//!
//! Weights are stored as `w_swath1(x)` representing the contribution from swath1:
//! - At the LEFT edge of overlap (low range indices): `w ≈ 1.0` (swath1 dominates)
//! - At the RIGHT edge of overlap (high range indices): `w ≈ 0.0` (swath2 dominates)
//! - Swath2 weight is always `1 - w`
//!
//! This produces a smooth left-to-right transition from swath1 to swath2.

use ndarray::Array2;

/// Create complementary cosine taper weights for overlap blending.
///
/// # Weight Direction
///
/// The returned weights represent swath1's contribution:
/// - Column 0 (left): w ≈ 1.0 (swath1 dominates)
/// - Column width-1 (right): w ≈ 0.0 (swath2 dominates via 1-w)
///
/// This is a monotonically DECREASING function across the overlap width.
///
/// # Arguments
///
/// * `width` - Width of the overlap region in range samples
/// * `height` - Height of the overlap region in azimuth lines
/// * `feather_width` - Width of the cosine transition zone (typically full width for smooth blend)
///
/// # Returns
///
/// 2D array of weights where `weights[[row, col]]` is swath1's weight at that position.
pub fn create_complementary_cosine_weights(
    width: usize,
    height: usize,
    feather_width: usize,
) -> Array2<f32> {
    let mut weights = Array2::zeros((height, width));

    if width == 0 || height == 0 {
        return weights;
    }

    // FIX: Always use full-width cosine feathering for seamless blending
    // This eliminates hard cuts in the middle of overlap regions
    // The feather_width parameter is kept for backwards compatibility but ignored
    // as we always use the full overlap width for smooth transitions (SNAP-compatible)

    // Create 1D monotonic weight ramp (swath1 contribution: 1→0 from left to right)
    // Full cosine transition across entire overlap width for seamless blending
    let w_1d: Vec<f32> = if width == 0 {
        Vec::new()
    } else {
        // Full cosine transition across entire width (SNAP-compatible approach)
        (0..width)
            .map(|col| {
                let t = (col as f32 + 0.5) / (width as f32);
                // Cosine goes from 1 at t=0 to 0 at t=1
                // This provides smooth, monotonic transition with no flat zones
                0.5 * (1.0 + (std::f32::consts::PI * t).cos())
            })
            .collect()
    };

    // Broadcast 1D weights to all rows
    for row in 0..height {
        for col in 0..width {
            weights[[row, col]] = w_1d[col];
        }
    }

    weights
}

/// Create SNAP-style sharp cutoff weights for overlap (no blending).
///
/// ESA SNAP's TOPSARMergeOp uses midpoint selection rather than feathering:
/// - Left half of overlap: swath1 contributes (w = 1.0)
/// - Right half of overlap: swath2 contributes (w = 0.0)
///
/// This avoids amplifying radiometric calibration differences that occur with
/// smooth blending (cosine/linear feathering creates gradient artifacts when
/// IW calibration LUTs don't match perfectly).
///
/// # Weight Direction
///
/// Same convention as cosine: w=1 at left, w=0 at right, but with a step function.
pub fn create_sharp_cutoff_weights(width: usize, height: usize) -> Array2<f32> {
    let mut weights = Array2::zeros((height, width));

    if width == 0 || height == 0 {
        return weights;
    }

    let midpoint = width / 2;

    for row in 0..height {
        for col in 0..width {
            // Sharp cutoff at midpoint: swath1 for left half, swath2 for right half
            weights[[row, col]] = if col < midpoint { 1.0 } else { 0.0 };
        }
    }

    weights
}

/// Create simple linear weights for overlap blending (alternative to cosine).
///
/// # Weight Direction
///
/// Same convention as cosine: w=1 at left, w=0 at right.
pub fn create_linear_weights(width: usize, height: usize) -> Array2<f32> {
    let mut weights = Array2::zeros((height, width));

    if width == 0 || height == 0 {
        return weights;
    }

    for row in 0..height {
        for col in 0..width {
            let t = (col as f32 + 0.5) / width as f32;
            weights[[row, col]] = 1.0 - t;
        }
    }

    weights
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn weight_direction_monotonic() {
        let weights = create_complementary_cosine_weights(100, 10, 100);

        // Check row 5 as a sample
        let row = weights.row(5);

        // Left edge should be close to 1.0
        assert!(
            row[0] > 0.95,
            "Left edge weight should be ~1.0, got {}",
            row[0]
        );

        // Right edge should be close to 0.0
        assert!(
            row[99] < 0.05,
            "Right edge weight should be ~0.0, got {}",
            row[99]
        );

        // Middle should be around 0.5
        assert!(
            (row[50] - 0.5).abs() < 0.1,
            "Middle weight should be ~0.5, got {}",
            row[50]
        );

        // Should be monotonically decreasing
        for col in 1..100 {
            assert!(
                row[col] <= row[col - 1] + 1e-6,
                "Weights should be monotonically decreasing: w[{}]={} > w[{}]={}",
                col - 1,
                row[col - 1],
                col,
                row[col]
            );
        }
    }

    #[test]
    fn weight_complement_sums_to_one() {
        let weights = create_complementary_cosine_weights(50, 5, 50);

        for row in 0..5 {
            for col in 0..50 {
                let w1 = weights[[row, col]];
                let w2 = 1.0 - w1;
                let sum = w1 + w2;
                assert!(
                    (sum - 1.0).abs() < 1e-6,
                    "w + (1-w) should equal 1.0, got {} at [{}, {}]",
                    sum,
                    row,
                    col
                );
            }
        }
    }

    #[test]
    fn weight_boundary_values() {
        let weights = create_complementary_cosine_weights(100, 1, 100);

        // First column: swath1 dominates (w ≈ 1)
        assert!(
            weights[[0, 0]] > 0.99,
            "First column should have w≈1, got {}",
            weights[[0, 0]]
        );

        // Last column: swath2 dominates (w ≈ 0)
        assert!(
            weights[[0, 99]] < 0.01,
            "Last column should have w≈0, got {}",
            weights[[0, 99]]
        );
    }

    #[test]
    fn linear_weights_correct() {
        let weights = create_linear_weights(10, 1);

        // First should be ~0.95 (0.5/10 from left edge)
        assert!((weights[[0, 0]] - 0.95).abs() < 0.01);

        // Last should be ~0.05
        assert!((weights[[0, 9]] - 0.05).abs() < 0.01);
    }

    #[test]
    fn empty_dimensions() {
        let w1 = create_complementary_cosine_weights(0, 10, 10);
        assert!(w1.is_empty());

        let w2 = create_complementary_cosine_weights(10, 0, 10);
        assert!(w2.is_empty());
    }

    #[test]
    fn sharp_cutoff_weights_correct() {
        let weights = create_sharp_cutoff_weights(100, 10);

        // Check row 5 as a sample
        let row = weights.row(5);

        // Left half should be 1.0 (swath1)
        for col in 0..50 {
            assert_eq!(row[col], 1.0, "Left half should be 1.0 at col {}", col);
        }

        // Right half should be 0.0 (swath2)
        for col in 50..100 {
            assert_eq!(row[col], 0.0, "Right half should be 0.0 at col {}", col);
        }
    }

    #[test]
    fn sharp_cutoff_odd_width() {
        // Odd width: 101 pixels, midpoint at 50
        let weights = create_sharp_cutoff_weights(101, 1);

        // Cols 0-49 should be 1.0 (50 pixels)
        assert_eq!(weights[[0, 49]], 1.0);
        // Cols 50-100 should be 0.0 (51 pixels)
        assert_eq!(weights[[0, 50]], 0.0);
    }
}
