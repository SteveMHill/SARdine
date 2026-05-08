//! Centralized DEM Gradient Operators
//!
//! This module provides unified gradient computation for terrain processing,
//! implementing the DEM gradient centralization (DUP-4) from the structural audit.
//!
//! # Supported Operators
//!
//! - **Horn 3×3**: Industry standard (most noise-robust) - recommended for SAR
//! - **Sobel 3×3**: Equivalent to Horn mathematically, explicit name for familiarity
//! - **Central Difference**: Simple 4-point stencil (faster, less noise-robust)
//!
//! # Scientific Reference
//!
//! Horn, B.K.P. (1981) "Hill shading and the reflectance map"
//! Proceedings of the IEEE, 69(1), 14-47
//!
//! # Usage
//!
//! ```rust,ignore
//! use sardine::core::terrain_correction::gradient::{compute_gradients, GradientOperator};
//!
//! let (dz_dx, dz_dy) = compute_gradients(&dem, GradientOperator::Horn3x3, pixel_spacing)?;
//! ```

use ndarray::Array2;
use rayon::prelude::*;

use crate::types::{SarError, SarResult};

/// Gradient operator selection
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum GradientOperator {
    /// Horn 3×3 operator (industry standard, most noise-robust)
    /// Weights: center=2, corners=1, total weight=8 per direction
    /// Reference: Horn 1981
    #[default]
    Horn3x3,

    /// Sobel 3×3 operator (mathematically equivalent to Horn)
    /// Provided as explicit name for users familiar with Sobel terminology
    Sobel3x3,

    /// Central difference (simple 4-point stencil)
    /// Faster but less noise-robust
    CentralDifference,
}

impl GradientOperator {
    /// Get operator label for logging/metadata
    pub fn label(&self) -> &'static str {
        match self {
            GradientOperator::Horn3x3 => "Horn 3×3",
            GradientOperator::Sobel3x3 => "Sobel 3×3",
            GradientOperator::CentralDifference => "Central Difference",
        }
    }
}

/// Compute DEM gradients using the specified operator.
///
/// # Arguments
///
/// * `dem` - 2D DEM array (rows × cols)
/// * `operator` - Gradient operator to use
/// * `dx` - Pixel spacing in X direction (meters)
/// * `dy` - Pixel spacing in Y direction (meters)
///
/// # Returns
///
/// * Tuple of (∂z/∂x, ∂z/∂y) arrays
///
/// # Errors
///
/// Returns error if DEM is too small (< 3×3)
pub fn compute_gradients(
    dem: &Array2<f32>,
    operator: GradientOperator,
    dx: f64,
    dy: f64,
) -> SarResult<(Array2<f32>, Array2<f32>)> {
    let (rows, cols) = dem.dim();

    if rows < 3 || cols < 3 {
        return Err(SarError::Processing(
            "DEM too small for gradient computation (minimum 3×3)".to_string(),
        ));
    }

    match operator {
        GradientOperator::Horn3x3 | GradientOperator::Sobel3x3 => {
            compute_horn_gradients(dem, dx, dy)
        }
        GradientOperator::CentralDifference => compute_central_difference_gradients(dem, dx, dy),
    }
}

/// Compute gradients using Horn 3×3 operator (parallel).
///
/// Formula for ∂z/∂x:
/// ```text
/// dz_dx = (z[i-1,j+1] + 2*z[i,j+1] + z[i+1,j+1] - z[i-1,j-1] - 2*z[i,j-1] - z[i+1,j-1]) / (8 * dx)
/// ```
fn compute_horn_gradients(
    dem: &Array2<f32>,
    dx: f64,
    dy: f64,
) -> SarResult<(Array2<f32>, Array2<f32>)> {
    let (rows, cols) = dem.dim();
    let mut dz_dx = Array2::<f32>::zeros((rows, cols));
    let mut dz_dy = Array2::<f32>::zeros((rows, cols));

    let dx_f = dx as f32;
    let dy_f = dy as f32;
    let inv_8dx = 1.0 / (8.0 * dx_f);
    let inv_8dy = 1.0 / (8.0 * dy_f);

    // Parallel interior computation
    let row_indices: Vec<usize> = (1..rows - 1).collect();
    let row_results: Vec<Vec<(usize, f32, f32)>> = row_indices
        .par_iter()
        .map(|&i| {
            let mut local_results = Vec::with_capacity(cols - 2);
            for j in 1..cols - 1 {
                // Horn 3×3 gradient for X (right column - left column)
                let px = (dem[[i - 1, j + 1]] + 2.0 * dem[[i, j + 1]] + dem[[i + 1, j + 1]]
                    - dem[[i - 1, j - 1]]
                    - 2.0 * dem[[i, j - 1]]
                    - dem[[i + 1, j - 1]])
                    * inv_8dx;

                // Horn 3×3 gradient for Y (bottom row - top row)
                let py = (dem[[i + 1, j - 1]] + 2.0 * dem[[i + 1, j]] + dem[[i + 1, j + 1]]
                    - dem[[i - 1, j - 1]]
                    - 2.0 * dem[[i - 1, j]]
                    - dem[[i - 1, j + 1]])
                    * inv_8dy;

                local_results.push((j, px, py));
            }
            local_results
        })
        .collect();

    // Apply results
    for (i, results) in row_indices.iter().zip(row_results.iter()) {
        for &(j, px, py) in results {
            dz_dx[[*i, j]] = px;
            dz_dy[[*i, j]] = py;
        }
    }

    // Fill edges using nearest interior values
    fill_gradient_edges(&mut dz_dx, &mut dz_dy);

    Ok((dz_dx, dz_dy))
}

/// Compute gradients using central differences (parallel).
///
/// Formula for ∂z/∂x:
/// ```text
/// dz_dx = (z[i, j+1] - z[i, j-1]) / (2 * dx)
/// ```
fn compute_central_difference_gradients(
    dem: &Array2<f32>,
    dx: f64,
    dy: f64,
) -> SarResult<(Array2<f32>, Array2<f32>)> {
    let (rows, cols) = dem.dim();
    let mut dz_dx = Array2::<f32>::zeros((rows, cols));
    let mut dz_dy = Array2::<f32>::zeros((rows, cols));

    let dx_f = dx as f32;
    let dy_f = dy as f32;
    let inv_2dx = 1.0 / (2.0 * dx_f);
    let inv_2dy = 1.0 / (2.0 * dy_f);

    // Parallel interior computation
    let row_indices: Vec<usize> = (1..rows - 1).collect();
    let row_results: Vec<Vec<(usize, f32, f32)>> = row_indices
        .par_iter()
        .map(|&i| {
            let mut local_results = Vec::with_capacity(cols - 2);
            for j in 1..cols - 1 {
                // Central difference for X
                let px = (dem[[i, j + 1]] - dem[[i, j - 1]]) * inv_2dx;

                // Central difference for Y
                let py = (dem[[i + 1, j]] - dem[[i - 1, j]]) * inv_2dy;

                local_results.push((j, px, py));
            }
            local_results
        })
        .collect();

    // Apply results
    for (i, results) in row_indices.iter().zip(row_results.iter()) {
        for &(j, px, py) in results {
            dz_dx[[*i, j]] = px;
            dz_dy[[*i, j]] = py;
        }
    }

    // Fill edges using nearest interior values
    fill_gradient_edges(&mut dz_dx, &mut dz_dy);

    Ok((dz_dx, dz_dy))
}

/// Fill edge pixels with nearest interior values.
fn fill_gradient_edges(dz_dx: &mut Array2<f32>, dz_dy: &mut Array2<f32>) {
    let (rows, cols) = dz_dx.dim();

    // Top and bottom rows
    for j in 0..cols {
        dz_dx[[0, j]] = dz_dx[[1, j.min(cols - 2).max(1)]];
        dz_dy[[0, j]] = dz_dy[[1, j.min(cols - 2).max(1)]];
        dz_dx[[rows - 1, j]] = dz_dx[[rows - 2, j.min(cols - 2).max(1)]];
        dz_dy[[rows - 1, j]] = dz_dy[[rows - 2, j.min(cols - 2).max(1)]];
    }

    // Left and right columns (excluding corners already filled)
    for i in 1..rows - 1 {
        dz_dx[[i, 0]] = dz_dx[[i, 1]];
        dz_dy[[i, 0]] = dz_dy[[i, 1]];
        dz_dx[[i, cols - 1]] = dz_dx[[i, cols - 2]];
        dz_dy[[i, cols - 1]] = dz_dy[[i, cols - 2]];
    }
}

/// Compute surface normal from gradients.
///
/// Surface normal: n = (-∂z/∂x, -∂z/∂y, 1), normalized to unit length.
///
/// # Arguments
///
/// * `dz_dx` - Gradient in X direction
/// * `dz_dy` - Gradient in Y direction
/// * `i`, `j` - Pixel indices
///
/// # Returns
///
/// * Unit normal vector [nx, ny, nz]
pub fn gradient_to_normal(dz_dx: f32, dz_dy: f32) -> [f32; 3] {
    let nx = -dz_dx;
    let ny = -dz_dy;
    let nz = 1.0_f32;

    let norm = (nx * nx + ny * ny + nz * nz).sqrt();
    if norm > 1e-6 {
        [nx / norm, ny / norm, nz / norm]
    } else {
        [0.0, 0.0, 1.0] // Flat surface fallback
    }
}

/// Compute slope magnitude from gradients.
///
/// Slope = √(1 + (∂z/∂x)² + (∂z/∂y)²)
///
/// This is used in area-projection RTC normalization.
pub fn slope_magnitude(dz_dx: f32, dz_dy: f32) -> f32 {
    (1.0 + dz_dx * dz_dx + dz_dy * dz_dy).sqrt()
}

#[cfg(test)]
mod tests {
    use super::*;
    use ndarray::array;

    #[test]
    fn test_horn_gradient_flat_dem() {
        // Flat DEM at elevation 100m
        let dem = Array2::<f32>::from_elem((5, 5), 100.0);
        let (dz_dx, dz_dy) =
            compute_gradients(&dem, GradientOperator::Horn3x3, 10.0, 10.0).unwrap();

        // Gradients should be zero for flat surface
        for i in 1..4 {
            for j in 1..4 {
                assert!(
                    dz_dx[[i, j]].abs() < 1e-6,
                    "Expected zero gradient, got {}",
                    dz_dx[[i, j]]
                );
                assert!(
                    dz_dy[[i, j]].abs() < 1e-6,
                    "Expected zero gradient, got {}",
                    dz_dy[[i, j]]
                );
            }
        }
    }

    #[test]
    fn test_horn_gradient_east_slope() {
        // DEM with linear slope in X direction: z = x * 10 (10m rise per 10m run = 45° slope)
        let mut dem = Array2::<f32>::zeros((5, 5));
        for i in 0..5 {
            for j in 0..5 {
                dem[[i, j]] = j as f32 * 10.0;
            }
        }

        let (dz_dx, dz_dy) =
            compute_gradients(&dem, GradientOperator::Horn3x3, 10.0, 10.0).unwrap();

        // Gradient in X should be 1.0 (10m / 10m spacing)
        // Gradient in Y should be 0.0
        for i in 1..4 {
            for j in 1..4 {
                assert!(
                    (dz_dx[[i, j]] - 1.0).abs() < 0.01,
                    "Expected dz_dx=1.0, got {}",
                    dz_dx[[i, j]]
                );
                assert!(
                    dz_dy[[i, j]].abs() < 0.01,
                    "Expected dz_dy=0.0, got {}",
                    dz_dy[[i, j]]
                );
            }
        }
    }

    #[test]
    fn test_central_diff_matches_horn_for_linear_slope() {
        // Linear slope should give same result for both operators
        let mut dem = Array2::<f32>::zeros((5, 5));
        for i in 0..5 {
            for j in 0..5 {
                dem[[i, j]] = j as f32 * 10.0;
            }
        }

        let (horn_dx, _) = compute_gradients(&dem, GradientOperator::Horn3x3, 10.0, 10.0).unwrap();
        let (cd_dx, _) =
            compute_gradients(&dem, GradientOperator::CentralDifference, 10.0, 10.0).unwrap();

        // For linear slope, both should give same result
        for i in 1..4 {
            for j in 1..4 {
                assert!(
                    (horn_dx[[i, j]] - cd_dx[[i, j]]).abs() < 0.01,
                    "Horn: {}, CentralDiff: {}",
                    horn_dx[[i, j]],
                    cd_dx[[i, j]]
                );
            }
        }
    }

    #[test]
    fn test_slope_magnitude() {
        // 45° slope in X direction
        let slope = slope_magnitude(1.0, 0.0);
        assert!((slope - 2.0_f32.sqrt()).abs() < 0.001);

        // Flat surface
        let flat = slope_magnitude(0.0, 0.0);
        assert!((flat - 1.0).abs() < 0.001);
    }

    #[test]
    fn test_gradient_to_normal() {
        // Flat surface should have normal pointing up
        let normal = gradient_to_normal(0.0, 0.0);
        assert!((normal[2] - 1.0).abs() < 0.001);

        // 45° slope in X should have nx = -1/sqrt(2), nz = 1/sqrt(2)
        let normal = gradient_to_normal(1.0, 0.0);
        let expected_component = 1.0 / 2.0_f32.sqrt();
        assert!(
            (normal[0] + expected_component).abs() < 0.01,
            "Expected nx=-0.707, got {}",
            normal[0]
        );
        assert!(
            (normal[2] - expected_component).abs() < 0.01,
            "Expected nz=0.707, got {}",
            normal[2]
        );
    }

    #[test]
    fn test_dem_too_small() {
        let small = Array2::<f32>::zeros((2, 2));
        let result = compute_gradients(&small, GradientOperator::Horn3x3, 10.0, 10.0);
        assert!(result.is_err());
    }
}
