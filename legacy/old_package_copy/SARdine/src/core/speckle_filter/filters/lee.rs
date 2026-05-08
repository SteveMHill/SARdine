#![allow(dead_code, unused_variables)]
//! Lee filter family implementations
//!
//! Includes the original Lee filter, Enhanced Lee, Refined Lee, and Lee Sigma filters.
//! These adaptive filters use local statistics to balance noise reduction with edge preservation.

use crate::core::speckle_filter::stats::IntegralImage;
use crate::types::SarResult;
use ndarray::{s, Array2};
use rayon::prelude::*;

/// Lee filter implementations for SpeckleFilter
impl super::SpeckleFilter {
    /// Apply Lee filter (adaptive)
    /// Implements the adaptive Lee filter based on local statistics
    /// Reference: Lee, J.S. (1980). Digital Image Enhancement and Noise Filtering by Use of Local Statistics
    /// CORRECTED: Uses proper literature-based weight formula for scientific accuracy
    pub(super) fn apply_lee_filter(&self, image: &Array2<f32>) -> SarResult<Array2<f32>> {
        let (height, width) = image.dim();
        let mut filtered = Array2::zeros((height, width));
        let half_window = self.params.window_size / 2;

        log::debug!(
            "Applying Lee filter with window size {}",
            self.params.window_size
        );

        // Build integral image for large windows
        let integral = if self.params.window_size >= 9 {
            Some(IntegralImage::new(image))
        } else {
            None
        };

        // OPTIMIZATION: Parallelize row processing with rayon (Phase 3)
        // Each row's computation is independent once integral image is built
        let num_looks = self.params.num_looks;
        let window_size = self.params.window_size;
        let cu = 1.0 / num_looks.sqrt(); // Theoretical CV for fully developed speckle
        const EPSILON: f32 = 1e-6;

        // Process interior rows in parallel
        let row_results: Vec<(usize, Vec<f32>)> = (half_window..(height - half_window))
            .into_par_iter()
            .map(|i| {
                let mut row_data = Vec::with_capacity(width - 2 * half_window);
                for j in half_window..(width - half_window) {
                    let center_value = image[[i, j]];

                    let result = if !center_value.is_finite() || center_value <= 0.0 {
                        center_value
                    } else {
                        // Calculate local statistics using optimized method
                        let (local_mean, local_variance) = if window_size >= 9 {
                            if let Some(ref integral_img) = integral {
                                let i_start = i.saturating_sub(half_window);
                                let i_end = (i + half_window).min(height - 1);
                                let j_start = j.saturating_sub(half_window);
                                let j_end = (j + half_window).min(width - 1);
                                integral_img.window_stats(i_start, j_start, i_end, j_end)
                            } else {
                                calculate_local_statistics_inline(
                                    image,
                                    i,
                                    j,
                                    half_window,
                                    height,
                                    width,
                                )
                            }
                        } else {
                            calculate_local_statistics_inline(
                                image,
                                i,
                                j,
                                half_window,
                                height,
                                width,
                            )
                        };

                        if local_mean <= 0.0 {
                            center_value
                        } else {
                            // Lee filter coefficient of variation
                            let cv = local_variance.sqrt() / local_mean;

                            // Apply CORRECTED Lee filter formula
                            if cv <= cu + EPSILON {
                                // Homogeneous area - use local mean
                                local_mean
                            } else {
                                // Heterogeneous area - weighted combination
                                // CORRECTED FORMULA: weight = (cv² - cu²) / (cv² + cu²)
                                let weight = (cv * cv - cu * cu) / (cv * cv + cu * cu);
                                (local_mean + weight * (center_value - local_mean)).max(0.0)
                            }
                        }
                    };
                    row_data.push(result);
                }
                (i, row_data)
            })
            .collect();

        // Write results back to output array
        for (i, row_data) in row_results {
            for (j_offset, &value) in row_data.iter().enumerate() {
                filtered[[i, half_window + j_offset]] = value;
            }
        }

        // Handle borders by copying original values
        self.handle_borders(image, &mut filtered, half_window);

        Ok(filtered)
    }

    /// Apply Enhanced Lee filter with edge detection
    /// More sophisticated version that preserves edges better
    /// CORRECTED: Implements robust numerical handling for edge cases
    pub(super) fn apply_enhanced_lee_filter(&self, image: &Array2<f32>) -> SarResult<Array2<f32>> {
        let (height, width) = image.dim();
        let mut filtered = Array2::zeros((height, width));
        let half_window = self.params.window_size / 2;

        log::debug!(
            "Applying Enhanced Lee filter with window size {}",
            self.params.window_size
        );

        // Build integral image for large windows
        let integral = if self.params.window_size >= 9 {
            Some(IntegralImage::new(image))
        } else {
            None
        };

        for i in half_window..(height - half_window) {
            for j in half_window..(width - half_window) {
                let center_value = image[[i, j]];

                if !center_value.is_finite() || center_value <= 0.0 {
                    filtered[[i, j]] = center_value;
                    continue;
                }

                // Calculate local statistics using optimized method
                let (local_mean, local_variance) = self.calculate_local_statistics_optimized(
                    image,
                    i,
                    j,
                    half_window,
                    integral.as_ref(),
                );

                if local_mean <= 0.0 {
                    filtered[[i, j]] = center_value;
                    continue;
                }

                // Enhanced Lee with robust coefficient of variation limits and numerical stability
                let cv = if local_variance > 0.0 && local_mean > 0.0 {
                    local_variance.sqrt() / local_mean
                } else {
                    0.0
                };

                let cu = 1.0 / self.params.num_looks.sqrt();
                let cmax = 1.73; // Maximum CV for speckle noise

                // Epsilon for numerical stability near critical points
                const EPSILON: f32 = 1e-6;
                const MIN_VARIANCE: f32 = 1e-8;

                let result = if local_variance < MIN_VARIANCE {
                    // Very low variance - homogeneous area
                    local_mean
                } else if cv <= cu + EPSILON {
                    // Homogeneous area - use local mean
                    local_mean
                } else if cv < cmax - EPSILON {
                    // Heterogeneous area - apply Enhanced Lee with bounded weights
                    let cv_sq = cv * cv;
                    let cu_sq = cu * cu;

                    // Prevent division by very small numbers
                    let denominator = cv_sq * (1.0 + cu_sq);
                    if denominator > EPSILON {
                        let weight = (cv_sq - cu_sq) / denominator;
                        // Clamp weight to reasonable bounds
                        let bounded_weight = weight.max(0.0).min(1.0);
                        local_mean + bounded_weight * (center_value - local_mean)
                    } else {
                        local_mean
                    }
                } else {
                    // Strong edge or point target - preserve original
                    center_value
                };

                filtered[[i, j]] = result.max(0.0);
            }
        }

        self.handle_borders(image, &mut filtered, half_window);
        Ok(filtered)
    }

    /// Apply Lee Sigma filter (edge-preserving)
    /// Uses sigma-based edge detection to preserve linear features
    pub(super) fn apply_lee_sigma_filter(&self, image: &Array2<f32>) -> SarResult<Array2<f32>> {
        let (height, width) = image.dim();
        let mut filtered = Array2::zeros((height, width));
        let half_window = self.params.window_size / 2;

        log::debug!(
            "Applying Lee Sigma filter with window size {}",
            self.params.window_size
        );

        for i in half_window..(height - half_window) {
            for j in half_window..(width - half_window) {
                let center_value = image[[i, j]];

                if !center_value.is_finite() || center_value <= 0.0 {
                    filtered[[i, j]] = center_value;
                    continue;
                }

                // Extract local window
                let window = image.slice(s![
                    i - half_window..=i + half_window,
                    j - half_window..=j + half_window
                ]);

                let local_mean = self.calculate_window_mean(&window);
                let local_std = self.calculate_window_variance(&window, local_mean).sqrt();

                if local_mean <= 0.0 {
                    filtered[[i, j]] = center_value;
                    continue;
                }

                // Sigma range for homogeneous pixels
                let sigma_range = self.params.edge_threshold * local_std;
                let lower_bound = local_mean - sigma_range;
                let upper_bound = local_mean + sigma_range;

                // Count pixels within sigma range
                let mut sum = 0.0;
                let mut count = 0;

                for value in window.iter() {
                    if *value >= lower_bound
                        && *value <= upper_bound
                        && value.is_finite()
                        && *value > 0.0
                    {
                        sum += *value;
                        count += 1;
                    }
                }

                let result = if count > 0 {
                    sum / count as f32
                } else {
                    center_value
                };

                filtered[[i, j]] = result.max(0.0);
            }
        }

        self.handle_borders(image, &mut filtered, half_window);
        Ok(filtered)
    }

    /// Apply Refined Lee filter with improved edge detection
    /// CORRECTED: Includes robust numerical handling for edge cases
    pub(super) fn apply_refined_lee_filter(&self, image: &Array2<f32>) -> SarResult<Array2<f32>> {
        let (height, width) = image.dim();
        let mut filtered = Array2::zeros((height, width));
        let half_window = self.params.window_size / 2;

        log::debug!(
            "Applying Refined Lee filter with window size {}",
            self.params.window_size
        );

        // Build integral image for large windows
        let integral = if self.params.window_size >= 9 {
            Some(IntegralImage::new(image))
        } else {
            None
        };

        for i in half_window..(height - half_window) {
            for j in half_window..(width - half_window) {
                let center_value = image[[i, j]];

                if !center_value.is_finite() || center_value <= 0.0 {
                    filtered[[i, j]] = center_value;
                    continue;
                }

                // Calculate directional gradients for edge detection
                let edge_strength = self.calculate_edge_strength(image, i, j);

                if edge_strength > self.params.edge_threshold {
                    // Strong edge - preserve original value
                    filtered[[i, j]] = center_value;
                    continue;
                }

                // Calculate local statistics using optimized method
                let (local_mean, local_variance) = self.calculate_local_statistics_optimized(
                    image,
                    i,
                    j,
                    half_window,
                    integral.as_ref(),
                );

                if local_mean <= 0.0 {
                    filtered[[i, j]] = center_value;
                    continue;
                }

                // Robust Lee filter processing with numerical stability
                let cv = if local_variance > 0.0 && local_mean > 0.0 {
                    local_variance.sqrt() / local_mean
                } else {
                    0.0
                };

                let cu = 1.0 / self.params.num_looks.sqrt();

                // Numerical stability constants
                const EPSILON: f32 = 1e-6;
                const MIN_VARIANCE: f32 = 1e-8;

                let result = if local_variance < MIN_VARIANCE {
                    // Very low variance - use local mean
                    local_mean
                } else if cv <= cu + EPSILON {
                    // Homogeneous area
                    local_mean
                } else {
                    // Heterogeneous area - apply corrected Lee filter
                    let cv_sq = cv * cv;
                    let cu_sq = cu * cu;

                    // Robust weight calculation with bounds checking
                    let denominator = cv_sq + cu_sq;
                    if denominator > EPSILON {
                        let weight = (cv_sq - cu_sq) / denominator;
                        let bounded_weight = weight.max(0.0).min(1.0);
                        local_mean + bounded_weight * (center_value - local_mean)
                    } else {
                        local_mean
                    }
                };

                filtered[[i, j]] = result.max(0.0);
            }
        }

        self.handle_borders(image, &mut filtered, half_window);
        Ok(filtered)
    }

    /// Apply Lee filter to a tile with cache optimization
    pub(super) fn apply_lee_filter_to_tile(
        &self,
        tile: &Array2<f32>,
        filtered: &mut Array2<f32>,
        integral: Option<&IntegralImage>,
    ) -> SarResult<()> {
        let (height, width) = tile.dim();
        let half_window = self.params.window_size / 2;
        let cu = 1.0 / self.params.num_looks.sqrt();
        const EPSILON: f32 = 1e-6;

        for i in 0..height {
            for j in 0..width {
                let center_value = tile[[i, j]];

                if !center_value.is_finite() || center_value <= 0.0 {
                    filtered[[i, j]] = center_value;
                    continue;
                }

                let (local_mean, local_variance) = if let Some(integral_img) = integral {
                    // Use integral image if available (for large windows)
                    let i_start = i.saturating_sub(half_window);
                    let i_end = (i + half_window).min(height - 1);
                    let j_start = j.saturating_sub(half_window);
                    let j_end = (j + half_window).min(width - 1);
                    integral_img.window_stats(i_start, j_start, i_end, j_end)
                } else {
                    // Use direct calculation for smaller windows
                    self.calculate_local_statistics_fast(tile, i, j, half_window)
                };

                if local_mean <= 0.0 {
                    filtered[[i, j]] = center_value;
                    continue;
                }

                let cv = local_variance.sqrt() / local_mean;

                let result = if cv <= cu + EPSILON {
                    local_mean
                } else {
                    let weight = (cv * cv - cu * cu) / (cv * cv + cu * cu);
                    local_mean + weight * (center_value - local_mean)
                };

                filtered[[i, j]] = result.max(0.0);
            }
        }

        Ok(())
    }

    /// Apply Enhanced Lee filter to a tile
    pub(super) fn apply_enhanced_lee_filter_to_tile(
        &self,
        tile: &Array2<f32>,
        filtered: &mut Array2<f32>,
        integral: Option<&IntegralImage>,
    ) -> SarResult<()> {
        let (height, width) = tile.dim();
        let half_window = self.params.window_size / 2;
        let cu = 1.0 / self.params.num_looks.sqrt();
        let cmax = 1.73;
        const EPSILON: f32 = 1e-6;

        for i in 0..height {
            for j in 0..width {
                let center_value = tile[[i, j]];

                if !center_value.is_finite() || center_value <= 0.0 {
                    filtered[[i, j]] = center_value;
                    continue;
                }

                let (local_mean, local_variance) = if let Some(integral_img) = integral {
                    let i_start = i.saturating_sub(half_window);
                    let i_end = (i + half_window).min(height - 1);
                    let j_start = j.saturating_sub(half_window);
                    let j_end = (j + half_window).min(width - 1);
                    integral_img.window_stats(i_start, j_start, i_end, j_end)
                } else {
                    self.calculate_local_statistics_fast(tile, i, j, half_window)
                };

                if local_mean <= 0.0 {
                    filtered[[i, j]] = center_value;
                    continue;
                }

                let cv = local_variance.sqrt() / local_mean;

                let result = if cv <= cu + EPSILON {
                    local_mean
                } else if cv < cmax {
                    let weight = (cv * cv - cu * cu) / (cv * cv * (1.0 + cu * cu));
                    local_mean + weight * (center_value - local_mean)
                } else {
                    center_value
                };

                filtered[[i, j]] = result;
            }
        }

        Ok(())
    }

    /// Enhanced Lee filter with chunked processing for better cache performance
    pub(super) fn apply_enhanced_lee_filter_chunked(
        &self,
        image: &Array2<f32>,
    ) -> SarResult<Array2<f32>> {
        log::debug!("Applying Enhanced Lee filter with chunked processing");

        let (height, width) = image.dim();
        let mut filtered = Array2::zeros((height, width));
        let half_window = self.params.window_size / 2;

        let cu = 1.0 / self.params.num_looks.sqrt();
        let cmax = 1.73;

        // Build integral image for large windows
        let integral = if self.params.window_size >= 9 {
            Some(IntegralImage::new(image))
        } else {
            None
        };

        // Process in chunks for better cache performance
        let chunk_size = 512;

        for i_chunk in (0..height).step_by(chunk_size) {
            let i_end = (i_chunk + chunk_size).min(height);

            for j_chunk in (0..width).step_by(chunk_size) {
                let j_end = (j_chunk + chunk_size).min(width);

                // Process chunk
                for i in i_chunk..i_end {
                    for j in j_chunk..j_end {
                        let center_value = image[[i, j]];

                        if !center_value.is_finite() || center_value <= 0.0 {
                            filtered[[i, j]] = center_value;
                            continue;
                        }

                        let (local_mean, local_variance) = self
                            .calculate_local_statistics_optimized(
                                image,
                                i,
                                j,
                                half_window,
                                integral.as_ref(),
                            );

                        if local_mean <= 0.0 {
                            filtered[[i, j]] = center_value;
                            continue;
                        }

                        let cv = local_variance.sqrt() / local_mean;

                        let result = if cv <= cu {
                            local_mean
                        } else if cv < cmax {
                            let weight = (cv * cv - cu * cu) / (cv * cv * (1.0 + cu * cu));
                            local_mean + weight * (center_value - local_mean)
                        } else {
                            center_value
                        };

                        filtered[[i, j]] = result;
                    }
                }
            }
        }

        Ok(filtered)
    }

    /// Parallel Enhanced Lee filter using Rayon (if available)
    #[cfg(feature = "parallel")]
    pub(super) fn apply_enhanced_lee_filter_parallel(
        &self,
        image: &Array2<f32>,
    ) -> SarResult<Array2<f32>> {
        use rayon::prelude::*;

        log::debug!("Applying Enhanced Lee filter with optimized parallel processing");

        let (height, width) = image.dim();
        let mut filtered = Array2::zeros((height, width));
        let half_window = self.params.window_size / 2;

        let cu = 1.0 / self.params.num_looks.sqrt();
        let cmax = 1.73;

        // Build integral image for large windows
        let integral = if self.params.window_size >= 9 {
            Some(IntegralImage::new(image))
        } else {
            None
        };

        // Determine optimal chunk size based on image dimensions and system capabilities
        let num_threads = rayon::current_num_threads();
        let total_pixels = height * width;
        let target_chunk_size = (total_pixels / (num_threads * 4)).max(1024);

        // Create row-based chunks for better cache locality
        let rows_per_chunk = (target_chunk_size / width).max(1).min(height);

        log::debug!(
            "Parallel processing: {} threads, {} rows per chunk",
            num_threads,
            rows_per_chunk
        );

        // Process in parallel chunks
        let chunk_results: Vec<_> = (0..height)
            .step_by(rows_per_chunk)
            .collect::<Vec<_>>()
            .into_par_iter()
            .map(|start_row| {
                let end_row = (start_row + rows_per_chunk).min(height);
                let mut chunk_result = Vec::with_capacity((end_row - start_row) * width);

                for i in start_row..end_row {
                    for j in 0..width {
                        let center_value = image[[i, j]];

                        let result = if !center_value.is_finite() || center_value <= 0.0 {
                            center_value
                        } else {
                            let (local_mean, local_variance) =
                                if let Some(ref integral_img) = integral {
                                    let i_start = i.saturating_sub(half_window);
                                    let i_end = (i + half_window).min(height - 1);
                                    let j_start = j.saturating_sub(half_window);
                                    let j_end = (j + half_window).min(width - 1);
                                    integral_img.window_stats(i_start, j_start, i_end, j_end)
                                } else {
                                    self.calculate_local_statistics_fast(image, i, j, half_window)
                                };

                            if local_mean <= 0.0 {
                                center_value
                            } else {
                                let cv = local_variance.sqrt() / local_mean;

                                if cv <= cu {
                                    local_mean
                                } else if cv < cmax {
                                    let weight = (cv * cv - cu * cu) / (cv * cv * (1.0 + cu * cu));
                                    local_mean + weight * (center_value - local_mean)
                                } else {
                                    center_value
                                }
                            }
                        };

                        chunk_result.push((i, j, result));
                    }
                }

                (start_row, chunk_result)
            })
            .collect();

        // Assign results back to filtered array
        for (start_row, chunk_data) in chunk_results {
            for (i, j, value) in chunk_data {
                filtered[[i, j]] = value;
            }
        }

        Ok(filtered)
    }

    #[cfg(not(feature = "parallel"))]
    pub(super) fn apply_enhanced_lee_filter_parallel(
        &self,
        image: &Array2<f32>,
    ) -> SarResult<Array2<f32>> {
        // Fallback to chunked processing if parallel feature is not available
        self.apply_enhanced_lee_filter_chunked(image)
    }
}

/// Inline helper for local statistics calculation (used in parallel context)
/// Avoids borrowing self in parallel closures
#[inline]
fn calculate_local_statistics_inline(
    image: &Array2<f32>,
    center_i: usize,
    center_j: usize,
    half_window: usize,
    height: usize,
    width: usize,
) -> (f32, f32) {
    let mut sum = 0.0f32;
    let mut sum_sq = 0.0f32;
    let mut count = 0u32;

    let i_start = center_i.saturating_sub(half_window);
    let i_end = (center_i + half_window + 1).min(height);
    let j_start = center_j.saturating_sub(half_window);
    let j_end = (center_j + half_window + 1).min(width);

    for i in i_start..i_end {
        for j in j_start..j_end {
            let pixel_val = image[[i, j]];
            if pixel_val.is_finite() && pixel_val >= 0.0 {
                sum += pixel_val;
                sum_sq += pixel_val * pixel_val;
                count += 1;
            }
        }
    }

    if count > 1 {
        let mean = sum / count as f32;
        let variance = (sum_sq / count as f32) - (mean * mean);
        (mean, variance.max(0.0))
    } else {
        (0.0, 0.0)
    }
}
