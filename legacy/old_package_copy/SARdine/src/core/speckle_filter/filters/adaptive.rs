//! Adaptive filter implementations
//!
//! Includes Frost and Gamma MAP filters that use adaptive weighting based on local statistics.

use crate::core::speckle_filter::stats::IntegralImage;
use crate::types::SarResult;
use ndarray::{s, Array2};
use rayon::prelude::*;

/// Adaptive filter implementations for SpeckleFilter
impl super::SpeckleFilter {
    /// Apply Gamma MAP filter (Maximum A Posteriori)
    /// Optimal filter based on Gamma distribution assumption for SAR speckle
    /// Reference: Lopes et al. (1993). Adaptive speckle filters and scene heterogeneity
    /// CORRECTED: Uses bounded weight formulation to prevent numerical instabilities
    pub(super) fn apply_gamma_map_filter(&self, image: &Array2<f32>) -> SarResult<Array2<f32>> {
        let (height, width) = image.dim();
        let mut filtered = Array2::zeros((height, width));
        let half_window = self.params.window_size / 2;

        log::debug!(
            "Applying Gamma MAP filter with window size {}",
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
        let cu = 1.0 / num_looks.sqrt();
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
                                calculate_local_statistics_inline_adaptive(
                                    image,
                                    i,
                                    j,
                                    half_window,
                                    height,
                                    width,
                                )
                            }
                        } else {
                            calculate_local_statistics_inline_adaptive(
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
                            // Gamma distribution parameters
                            let cv = local_variance.sqrt() / local_mean;

                            // CORRECTED: Bounded Gamma-MAP estimation
                            if cv <= cu + EPSILON {
                                // Homogeneous area - use local mean
                                local_mean
                            } else {
                                // Heterogeneous area - bounded MAP estimate
                                // CORRECTED FORMULA: weight = (cv² - cu²) / (cv² * (1 + cu²))
                                let weight = (cv * cv - cu * cu) / (cv * cv * (1.0 + cu * cu));
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

        self.handle_borders(image, &mut filtered, half_window);
        Ok(filtered)
    }

    /// Apply Frost filter (exponential weighting)
    /// Applies exponential weighting based on distance and local statistics
    /// Reference: Frost et al. (1982). A model for radar images and its application to adaptive digital filtering
    pub(super) fn apply_frost_filter(&self, image: &Array2<f32>) -> SarResult<Array2<f32>> {
        let (height, width) = image.dim();
        let mut filtered = Array2::zeros((height, width));
        let half_window = self.params.window_size / 2;

        log::debug!(
            "Applying Frost filter with window size {}",
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
                let local_variance = self.calculate_window_variance(&window, local_mean);

                if local_mean <= 0.0 {
                    filtered[[i, j]] = center_value;
                    continue;
                }

                // Frost filter parameters
                let cv = local_variance.sqrt() / local_mean;
                let k = self.params.damping_factor * cv * cv;

                let mut weighted_sum = 0.0;
                let mut weight_sum = 0.0;

                // Apply exponential weighting
                for wi in 0..self.params.window_size {
                    for wj in 0..self.params.window_size {
                        let value = window[[wi, wj]];
                        if value.is_finite() && value >= 0.0 {
                            let di = (wi as i32 - half_window as i32).abs() as f32;
                            let dj = (wj as i32 - half_window as i32).abs() as f32;
                            let distance = (di * di + dj * dj).sqrt();

                            let weight = (-k * distance).exp();
                            weighted_sum += weight * value;
                            weight_sum += weight;
                        }
                    }
                }

                let result = if weight_sum > 0.0 {
                    weighted_sum / weight_sum
                } else {
                    center_value
                };

                filtered[[i, j]] = result.max(0.0);
            }
        }

        self.handle_borders(image, &mut filtered, half_window);
        Ok(filtered)
    }

    /// Apply Gamma-MAP filter to a tile
    pub(super) fn apply_gamma_map_filter_to_tile(
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
                } else {
                    let weight = (cv * cv - cu * cu) / (cv * cv * (1.0 + cu * cu));
                    local_mean + weight * (center_value - local_mean)
                };

                filtered[[i, j]] = result.max(0.0);
            }
        }

        Ok(())
    }
}

/// Inline helper for local statistics calculation (used in parallel context)
/// Avoids borrowing self in parallel closures
#[inline]
fn calculate_local_statistics_inline_adaptive(
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
