//! Speckle filter implementations
//!
//! This module provides various speckle filtering algorithms for SAR imagery:
//! - Mean filter (simple averaging)
//! - Median filter (rank filter with quickselect optimization)
//! - Lee filter (adaptive)
//! - Enhanced Lee filter (with edge preservation)
//! - Lee Sigma filter (edge-preserving)
//! - Frost filter (exponential weighting)
//! - Gamma MAP filter (Maximum A Posteriori)
//! - Refined Lee filter (with edge detection)

mod adaptive;
mod lee;
mod median;
mod parallel;
mod utils;

use crate::core::speckle_filter::config::SpeckleFilterParams;
use crate::core::speckle_filter::stats::IntegralImage;
use crate::types::{SarError, SarResult};
use ndarray::Array2;

/// Available speckle filter types
#[derive(Debug, Clone, Copy)]
pub enum SpeckleFilterType {
    /// Mean filter (simple averaging)
    Mean,
    /// Median filter (rank filter)
    Median,
    /// Lee filter (adaptive)
    Lee,
    /// Enhanced Lee filter
    EnhancedLee,
    /// Lee Sigma filter (edge-preserving)
    LeeSigma,
    /// Frost filter (exponential weighting)
    Frost,
    /// Gamma MAP filter (Maximum A Posteriori)
    GammaMAP,
    /// Refined Lee filter
    RefinedLee,
}

/// Speckle filter processor
pub struct SpeckleFilter {
    pub(super) params: SpeckleFilterParams,
}

impl SpeckleFilter {
    /// Create a new speckle filter with default parameters
    pub fn new() -> Self {
        Self {
            params: SpeckleFilterParams::default(),
        }
    }

    /// Create a speckle filter with custom parameters
    pub fn with_params(params: SpeckleFilterParams) -> Self {
        Self { params }
    }

    /// Apply speckle filtering to SAR image
    pub fn apply_filter(
        &self,
        image: &Array2<f32>,
        filter_type: SpeckleFilterType,
    ) -> SarResult<Array2<f32>> {
        self.params.validate_window()?;

        log::info!("Applying {:?} speckle filter", filter_type);
        log::debug!("Filter parameters: {:?}", self.params);

        // Validate input
        let (height, width) = image.dim();
        if height < self.params.window_size || width < self.params.window_size {
            return Err(SarError::Processing(format!(
                "Image size {}x{} is too small for window size {}",
                height, width, self.params.window_size
            )));
        }

        let filtered = match filter_type {
            SpeckleFilterType::Mean => self.apply_mean_filter(image)?,
            SpeckleFilterType::Median => self.apply_median_filter(image)?,
            SpeckleFilterType::Lee => self.apply_lee_filter(image)?,
            SpeckleFilterType::EnhancedLee => self.apply_enhanced_lee_filter(image)?,
            SpeckleFilterType::LeeSigma => self.apply_lee_sigma_filter(image)?,
            SpeckleFilterType::Frost => self.apply_frost_filter(image)?,
            SpeckleFilterType::GammaMAP => self.apply_gamma_map_filter(image)?,
            SpeckleFilterType::RefinedLee => self.apply_refined_lee_filter(image)?,
        };

        log::info!("Speckle filtering completed successfully");
        Ok(filtered)
    }

    /// Apply mean filter (simple averaging)
    fn apply_mean_filter(&self, image: &Array2<f32>) -> SarResult<Array2<f32>> {
        log::debug!("Applying mean filter");

        let (height, width) = image.dim();
        let mut filtered = Array2::zeros((height, width));
        let half_window = self.params.window_size / 2;

        for i in 0..height {
            for j in 0..width {
                let mut sum = 0.0;
                let mut count = 0;

                // Process window
                for wi in 0..self.params.window_size {
                    for wj in 0..self.params.window_size {
                        let ii = i as i32 + wi as i32 - half_window as i32;
                        let jj = j as i32 + wj as i32 - half_window as i32;

                        if ii >= 0 && ii < height as i32 && jj >= 0 && jj < width as i32 {
                            let pixel_val = image[[ii as usize, jj as usize]];
                            if pixel_val.is_finite() && pixel_val >= 0.0 {
                                sum += pixel_val;
                                count += 1;
                            }
                        }
                    }
                }

                filtered[[i, j]] = if count > 0 {
                    sum / count as f32
                } else {
                    image[[i, j]]
                };
            }
        }

        Ok(filtered)
    }

    /// Optimized mean filter with better memory access patterns
    #[allow(dead_code)]
    fn apply_mean_filter_optimized(&self, image: &Array2<f32>) -> SarResult<Array2<f32>> {
        log::debug!("Applying optimized mean filter");

        let (height, width) = image.dim();
        let mut filtered = Array2::zeros((height, width));
        let half_window = self.params.window_size / 2;

        // Process in cache-friendly chunks
        let chunk_size = 256;

        for i_chunk in (0..height).step_by(chunk_size) {
            let i_end = (i_chunk + chunk_size).min(height);

            for j_chunk in (0..width).step_by(chunk_size) {
                let j_end = (j_chunk + chunk_size).min(width);

                for i in i_chunk..i_end {
                    for j in j_chunk..j_end {
                        let (mean, _) =
                            self.calculate_local_statistics_fast(image, i, j, half_window);
                        filtered[[i, j]] = if mean > 0.0 { mean } else { image[[i, j]] };
                    }
                }
            }
        }

        Ok(filtered)
    }

    /// Apply mean filter to a tile
    pub(super) fn apply_mean_filter_to_tile(
        &self,
        tile: &Array2<f32>,
        filtered: &mut Array2<f32>,
        integral: Option<&IntegralImage>,
    ) -> SarResult<()> {
        let (height, width) = tile.dim();
        let half_window = self.params.window_size / 2;

        for i in 0..height {
            for j in 0..width {
                let (mean, _) = if let Some(integral_img) = integral {
                    let i_start = i.saturating_sub(half_window);
                    let i_end = (i + half_window).min(height - 1);
                    let j_start = j.saturating_sub(half_window);
                    let j_end = (j + half_window).min(width - 1);
                    integral_img.window_stats(i_start, j_start, i_end, j_end)
                } else {
                    self.calculate_local_statistics_fast(tile, i, j, half_window)
                };

                filtered[[i, j]] = if mean > 0.0 { mean } else { tile[[i, j]] };
            }
        }

        Ok(())
    }

    /// Multi-scale speckle filtering
    pub fn apply_multiscale_filter(
        &self,
        image: &Array2<f32>,
        filter_type: SpeckleFilterType,
        scales: &[usize],
    ) -> SarResult<Array2<f32>> {
        log::info!(
            "Applying multi-scale speckle filtering with {} scales",
            scales.len()
        );

        let mut result = image.clone();

        for (i, &scale) in scales.iter().enumerate() {
            log::debug!("Processing scale {} with window size {}", i + 1, scale);

            // Create filter with current scale
            let mut scale_params = self.params.clone();
            scale_params.window_size = scale;
            let scale_filter = SpeckleFilter::with_params(scale_params);

            // Apply filter at current scale
            result = scale_filter.apply_filter(&result, filter_type)?;
        }

        log::info!("Multi-scale filtering completed");
        Ok(result)
    }

    /// Estimate number of looks from image statistics
    /// OPTIMIZATION #104: Single-pass mean and variance using Welford's algorithm
    /// This is numerically stable and eliminates the need to store all values
    pub fn estimate_number_of_looks(image: &Array2<f32>) -> SarResult<f32> {
        log::debug!("Estimating number of looks from image statistics (single-pass)");

        // Welford's online algorithm for single-pass mean and variance
        let mut count = 0u64;
        let mut mean = 0.0f64;
        let mut m2 = 0.0f64; // Sum of squared deviations

        for val in image.iter() {
            if val.is_finite() && *val >= 0.0 {
                count += 1;
                let delta = *val as f64 - mean;
                mean += delta / count as f64;
                let delta2 = *val as f64 - mean;
                m2 += delta * delta2;
            }
        }

        if count < 2 {
            return Err(SarError::Processing(
                "Insufficient valid pixels for ENL estimation".to_string(),
            ));
        }

        // Sample variance = m2 / (n - 1)
        let variance = m2 / (count - 1) as f64;

        // Number of looks = mean² / variance
        let num_looks = if variance > 0.0 {
            (mean * mean / variance) as f32
        } else {
            f32::MAX // Near-constant image has effectively infinite looks
        };

        log::info!(
            "Estimated number of looks: {:.2} (from {} valid pixels)",
            num_looks,
            count
        );
        Ok(num_looks.max(1.0))
    }

    /// Performance benchmark for speckle filters
    pub fn benchmark_speckle_filters(&self, image: &Array2<f32>) -> SarResult<()> {
        use std::time::Instant;

        log::info!("=== Speckle Filter Performance Benchmark ===");
        let (height, width) = image.dim();
        let total_pixels = height * width;
        log::info!("Image size: {}x{} = {} pixels", height, width, total_pixels);

        // Test Enhanced Lee filter
        log::info!("Testing Enhanced Lee filter variants...");

        // Standard implementation
        let start = Instant::now();
        let _result1 = self.apply_enhanced_lee_filter(image)?;
        let standard_time = start.elapsed();
        log::info!(
            "Standard Enhanced Lee: {:.2} seconds ({:.0} pixels/sec)",
            standard_time.as_secs_f64(),
            total_pixels as f64 / standard_time.as_secs_f64()
        );

        // Chunked implementation
        let start = Instant::now();
        let _result2 = self.apply_enhanced_lee_filter_chunked(image)?;
        let chunked_time = start.elapsed();
        log::info!(
            "Chunked Enhanced Lee: {:.2} seconds ({:.0} pixels/sec)",
            chunked_time.as_secs_f64(),
            total_pixels as f64 / chunked_time.as_secs_f64()
        );

        // Parallel implementation (if available)
        #[cfg(feature = "parallel")]
        {
            let start = Instant::now();
            let _result3 = self.apply_enhanced_lee_filter_parallel(image)?;
            let parallel_time = start.elapsed();
            log::info!(
                "Parallel Enhanced Lee: {:.2} seconds ({:.0} pixels/sec)",
                parallel_time.as_secs_f64(),
                total_pixels as f64 / parallel_time.as_secs_f64()
            );

            let parallel_speedup = standard_time.as_secs_f64() / parallel_time.as_secs_f64();
            log::info!("Parallel speedup: {:.2}x", parallel_speedup);
        }

        let chunked_speedup = standard_time.as_secs_f64() / chunked_time.as_secs_f64();
        log::info!("Chunked speedup: {:.2}x", chunked_speedup);
        log::info!("=== Benchmark Complete ===");

        Ok(())
    }
}

impl Default for SpeckleFilter {
    fn default() -> Self {
        Self::new()
    }
}
