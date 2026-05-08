//! Utility functions for speckle filtering
//!
//! Includes local statistics calculation, border handling, edge detection,
//! and other helper functions used by the filter implementations.

use crate::core::speckle_filter::stats::IntegralImage;
use ndarray::Array2;

/// Utility implementations for SpeckleFilter
impl super::SpeckleFilter {
    /// Calculate local statistics for a window (legacy, slower)
    #[allow(dead_code)]
    pub(super) fn calculate_local_statistics(
        &self,
        image: &Array2<f32>,
        center_i: usize,
        center_j: usize,
        half_window: usize,
    ) -> (f32, f32) {
        let (height, width) = image.dim();
        let mut values = Vec::new();

        for wi in 0..self.params.window_size {
            for wj in 0..self.params.window_size {
                let ii = center_i as i32 + wi as i32 - half_window as i32;
                let jj = center_j as i32 + wj as i32 - half_window as i32;

                if ii >= 0 && ii < height as i32 && jj >= 0 && jj < width as i32 {
                    let pixel_val = image[[ii as usize, jj as usize]];
                    if pixel_val.is_finite() && pixel_val >= 0.0 {
                        values.push(pixel_val);
                    }
                }
            }
        }

        if values.is_empty() {
            return (0.0, 0.0);
        }

        // Calculate mean
        let mean = values.iter().sum::<f32>() / values.len() as f32;

        // Calculate variance
        let variance = if values.len() > 1 {
            values.iter().map(|v| (v - mean) * (v - mean)).sum::<f32>() / (values.len() - 1) as f32
        } else {
            0.0
        };

        (mean, variance)
    }

    /// Fast local statistics calculation with optimized loop
    pub(super) fn calculate_local_statistics_fast(
        &self,
        image: &Array2<f32>,
        center_i: usize,
        center_j: usize,
        half_window: usize,
    ) -> (f32, f32) {
        let (height, width) = image.dim();
        let mut sum = 0.0;
        let mut sum_sq = 0.0;
        let mut count = 0;

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

    /// Ultra-fast local statistics using integral images for large windows (O(1) per pixel)
    /// Automatically switches to integral image method for window sizes >= 9x9
    pub(super) fn calculate_local_statistics_optimized(
        &self,
        image: &Array2<f32>,
        center_i: usize,
        center_j: usize,
        half_window: usize,
        integral: Option<&IntegralImage>,
    ) -> (f32, f32) {
        // Use integral images for large windows (>= 9x9)
        if self.params.window_size >= 9 {
            if let Some(integral_img) = integral {
                let (height, width) = image.dim();

                let i_start = center_i.saturating_sub(half_window);
                let i_end = (center_i + half_window).min(height - 1);
                let j_start = center_j.saturating_sub(half_window);
                let j_end = (center_j + half_window).min(width - 1);

                return integral_img.window_stats(i_start, j_start, i_end, j_end);
            }
        }

        // Fallback to fast direct calculation for smaller windows
        self.calculate_local_statistics_fast(image, center_i, center_j, half_window)
    }

    /// Calculate window mean for a 2D array slice
    pub(super) fn calculate_window_mean(&self, window: &ndarray::ArrayView2<f32>) -> f32 {
        let mut sum = 0.0;
        let mut count = 0;

        for value in window.iter() {
            if value.is_finite() && *value >= 0.0 {
                sum += *value;
                count += 1;
            }
        }

        if count > 0 {
            sum / count as f32
        } else {
            0.0
        }
    }

    /// Calculate window variance for a 2D array slice
    pub(super) fn calculate_window_variance(
        &self,
        window: &ndarray::ArrayView2<f32>,
        mean: f32,
    ) -> f32 {
        let mut sum_sq_diff = 0.0;
        let mut count = 0;

        for value in window.iter() {
            if value.is_finite() && *value >= 0.0 {
                let diff = *value - mean;
                sum_sq_diff += diff * diff;
                count += 1;
            }
        }

        if count > 1 {
            sum_sq_diff / (count - 1) as f32
        } else {
            0.0
        }
    }

    /// Handle border pixels by copying original values
    pub(super) fn handle_borders(
        &self,
        original: &Array2<f32>,
        filtered: &mut Array2<f32>,
        border_size: usize,
    ) {
        let (height, width) = original.dim();

        // Top and bottom borders
        for i in 0..border_size {
            for j in 0..width {
                filtered[[i, j]] = original[[i, j]];
                if height > border_size {
                    filtered[[height - 1 - i, j]] = original[[height - 1 - i, j]];
                }
            }
        }

        // Left and right borders
        for i in border_size..(height - border_size) {
            for j in 0..border_size {
                filtered[[i, j]] = original[[i, j]];
                if width > border_size {
                    filtered[[i, width - 1 - j]] = original[[i, width - 1 - j]];
                }
            }
        }
    }

    /// Calculate edge strength using gradient magnitude
    pub(super) fn calculate_edge_strength(&self, image: &Array2<f32>, i: usize, j: usize) -> f32 {
        let (height, width) = image.dim();

        if i == 0 || i >= height - 1 || j == 0 || j >= width - 1 {
            return 0.0;
        }

        // Simple gradient calculation
        let gx = image[[i, j + 1]] - image[[i, j - 1]];
        let gy = image[[i + 1, j]] - image[[i - 1, j]];

        (gx * gx + gy * gy).sqrt()
    }

    /// Calculate adaptive statistics for Enhanced Lee filter
    #[allow(dead_code)]
    pub(super) fn calculate_adaptive_statistics(
        &self,
        window: &ndarray::ArrayView2<f32>,
        _center_i: usize,
        _center_j: usize,
    ) -> (f32, f32) {
        // For now, use standard statistics
        // In a full implementation, this would include directional analysis
        let mean = self.calculate_window_mean(window);
        let variance = self.calculate_window_variance(window, mean);
        (mean, variance)
    }
}
