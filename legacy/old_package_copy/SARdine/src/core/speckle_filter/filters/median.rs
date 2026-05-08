#![allow(dead_code, unused_variables)]
//! Median filter implementations with optimizations
//!
//! Includes standard median filter, quickselect optimization, and histogram-based
//! median for 8-bit intensity images.

use crate::types::SarResult;
use ndarray::Array2;
use rayon::prelude::*;

/// Median filter implementations for SpeckleFilter
impl super::SpeckleFilter {
    /// Apply median filter
    /// OPTIMIZATION #82: Parallelized row processing using rayon
    pub(super) fn apply_median_filter(&self, image: &Array2<f32>) -> SarResult<Array2<f32>> {
        log::debug!("Applying parallel median filter");

        let (height, width) = image.dim();
        let half_window = self.params.window_size / 2;
        let window_size = self.params.window_size;
        let window_capacity = window_size * window_size;

        // Process rows in parallel - each row is independent
        let rows: Vec<Vec<f32>> = (0..height)
            .into_par_iter()
            .map(|i| {
                let mut row_result = vec![0.0f32; width];
                for j in 0..width {
                    let mut window_values = Vec::with_capacity(window_capacity);

                    // Collect window values
                    for wi in 0..window_size {
                        for wj in 0..window_size {
                            let ii = i as i32 + wi as i32 - half_window as i32;
                            let jj = j as i32 + wj as i32 - half_window as i32;

                            if ii >= 0 && ii < height as i32 && jj >= 0 && jj < width as i32 {
                                let pixel_val = image[[ii as usize, jj as usize]];
                                if pixel_val.is_finite() && pixel_val >= 0.0 {
                                    window_values.push(pixel_val);
                                }
                            }
                        }
                    }

                    if !window_values.is_empty() {
                        window_values.sort_by(|a, b| a.total_cmp(b));
                        row_result[j] = window_values[window_values.len() / 2];
                    } else {
                        row_result[j] = image[[i, j]];
                    }
                }
                row_result
            })
            .collect();

        // Assemble result from parallel row results
        let mut filtered = Array2::zeros((height, width));
        for (i, row) in rows.into_iter().enumerate() {
            for (j, val) in row.into_iter().enumerate() {
                filtered[[i, j]] = val;
            }
        }

        Ok(filtered)
    }

    /// Optimized median filter with efficient sorting
    pub(super) fn apply_median_filter_optimized(
        &self,
        image: &Array2<f32>,
    ) -> SarResult<Array2<f32>> {
        log::debug!("Applying optimized median filter");

        let (height, width) = image.dim();
        let mut filtered = Array2::zeros((height, width));
        let half_window = self.params.window_size / 2;

        // Check if image appears to be 8-bit intensity (values in [0, 255])
        let is_8bit = Self::is_8bit_intensity(image);

        for i in 0..height {
            for j in 0..width {
                if is_8bit {
                    // Use histogram-based median for 8-bit images
                    filtered[[i, j]] =
                        Self::histogram_median_8bit(image, i, j, self.params.window_size);
                } else {
                    // Use quickselect median for floating-point images
                    let mut window_values =
                        Vec::with_capacity(self.params.window_size * self.params.window_size);

                    let i_start = i.saturating_sub(half_window);
                    let i_end = (i + half_window + 1).min(height);
                    let j_start = j.saturating_sub(half_window);
                    let j_end = (j + half_window + 1).min(width);

                    for wi in i_start..i_end {
                        for wj in j_start..j_end {
                            let pixel_val = image[[wi, wj]];
                            if pixel_val.is_finite() && pixel_val >= 0.0 {
                                window_values.push(pixel_val);
                            }
                        }
                    }

                    if !window_values.is_empty() {
                        filtered[[i, j]] = Self::quickselect_median(&mut window_values);
                    } else {
                        filtered[[i, j]] = image[[i, j]];
                    }
                }
            }
        }

        Ok(filtered)
    }

    /// Apply median filter to a tile
    pub(super) fn apply_median_filter_to_tile(
        &self,
        tile: &Array2<f32>,
        filtered: &mut Array2<f32>,
    ) -> SarResult<()> {
        let (height, width) = tile.dim();
        let half_window = self.params.window_size / 2;
        let is_8bit = Self::is_8bit_intensity(tile);

        for i in 0..height {
            for j in 0..width {
                if is_8bit {
                    filtered[[i, j]] =
                        Self::histogram_median_8bit(tile, i, j, self.params.window_size);
                } else {
                    let mut window_values = Vec::new();

                    let i_start = i.saturating_sub(half_window);
                    let i_end = (i + half_window + 1).min(height);
                    let j_start = j.saturating_sub(half_window);
                    let j_end = (j + half_window + 1).min(width);

                    for wi in i_start..i_end {
                        for wj in j_start..j_end {
                            let pixel_val = tile[[wi, wj]];
                            if pixel_val.is_finite() && pixel_val >= 0.0 {
                                window_values.push(pixel_val);
                            }
                        }
                    }

                    filtered[[i, j]] = if !window_values.is_empty() {
                        Self::quickselect_median(&mut window_values)
                    } else {
                        tile[[i, j]]
                    };
                }
            }
        }

        Ok(())
    }

    /// Fast median computation using quickselect algorithm for floating-point data
    /// O(n) average case vs O(n log n) for full sorting
    pub(super) fn quickselect_median(data: &mut [f32]) -> f32 {
        if data.is_empty() {
            return 0.0;
        }
        let n = data.len();
        let k = n / 2;
        Self::quickselect(data, k)
    }

    /// Quickselect implementation for finding k-th element
    fn quickselect(data: &mut [f32], k: usize) -> f32 {
        if data.len() == 1 {
            return data[0];
        }

        let pivot_idx = Self::partition(data);

        if k == pivot_idx {
            data[k]
        } else if k < pivot_idx {
            Self::quickselect(&mut data[..pivot_idx], k)
        } else {
            Self::quickselect(&mut data[pivot_idx + 1..], k - pivot_idx - 1)
        }
    }

    /// Partition function for quickselect
    fn partition(data: &mut [f32]) -> usize {
        let pivot = data[data.len() - 1];
        let mut i = 0;

        for j in 0..data.len() - 1 {
            if data[j] <= pivot {
                data.swap(i, j);
                i += 1;
            }
        }
        data.swap(i, data.len() - 1);
        i
    }

    /// Histogram-based median for 8-bit intensity images (O(1) after preprocessing)
    /// This method is optimal for images with limited dynamic range
    pub(super) fn histogram_median_8bit(
        image: &Array2<f32>,
        center_i: usize,
        center_j: usize,
        window_size: usize,
    ) -> f32 {
        let half_window = window_size / 2;
        let (height, width) = image.dim();
        let mut histogram = [0u32; 256];
        let mut total_pixels = 0u32;

        // Build histogram for the window
        let i_start = center_i.saturating_sub(half_window);
        let i_end = (center_i + half_window + 1).min(height);
        let j_start = center_j.saturating_sub(half_window);
        let j_end = (center_j + half_window + 1).min(width);

        for i in i_start..i_end {
            for j in j_start..j_end {
                let pixel_val = image[[i, j]];
                if pixel_val.is_finite() && pixel_val >= 0.0 && pixel_val <= 255.0 {
                    // Round and clamp to [0,255] for robust binning
                    let bin = (pixel_val.round() as i32).clamp(0, 255) as usize;
                    histogram[bin] += 1;
                    total_pixels += 1;
                }
            }
        }

        if total_pixels == 0 {
            return image[[center_i, center_j]];
        }

        // Find median using histogram
        let median_pos = total_pixels / 2;
        let mut cumulative = 0u32;

        for (bin, &count) in histogram.iter().enumerate() {
            cumulative += count;
            if cumulative > median_pos {
                return bin as f32;
            }
        }

        // Fallback
        image[[center_i, center_j]]
    }

    /// Check if image appears to be 8-bit intensity data
    pub(super) fn is_8bit_intensity(image: &Array2<f32>) -> bool {
        // Sample a subset of pixels to check range
        let (height, width) = image.dim();
        let sample_size = 1000.min(height * width);
        let mut sample_count = 0;
        let mut in_range_count = 0;

        for i in (0..height).step_by((height / 100).max(1)) {
            for j in (0..width).step_by((width / 100).max(1)) {
                let val = image[[i, j]];
                if val.is_finite() {
                    sample_count += 1;
                    if val >= 0.0 && val <= 255.0 && val.fract() == 0.0 {
                        in_range_count += 1;
                    }
                    if sample_count >= sample_size {
                        break;
                    }
                }
            }
            if sample_count >= sample_size {
                break;
            }
        }

        if sample_count > 0 {
            (in_range_count as f32 / sample_count as f32) > 0.9
        } else {
            false
        }
    }
}
