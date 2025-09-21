use crate::types::{SarError, SarResult};
use ndarray::{Array2, s};

/// Integral image structure for O(1) window statistics
/// Precomputed sum and sum-of-squares for efficient window calculations
#[derive(Debug, Clone)]
struct IntegralImage {
    /// Integral of the image (prefix sums)
    sum_table: Array2<f64>,
    /// Integral of squared values for variance calculation
    sum_sq_table: Array2<f64>,
    /// Count table for valid pixels (non-zero, finite)
    count_table: Array2<u32>,
}

impl IntegralImage {
    /// Create integral image from input data
    fn new(image: &Array2<f32>) -> Self {
        let (height, width) = image.dim();
        let mut sum_table = Array2::zeros((height + 1, width + 1));
        let mut sum_sq_table = Array2::zeros((height + 1, width + 1));
        let mut count_table = Array2::zeros((height + 1, width + 1));

        // Build integral tables using dynamic programming
        for i in 1..=height {
            for j in 1..=width {
                let pixel = image[[i - 1, j - 1]];
                let (val, val_sq, cnt) = if pixel.is_finite() && pixel > 0.0 {
                    (pixel as f64, (pixel * pixel) as f64, 1u32)
                } else {
                    (0.0, 0.0, 0u32)
                };

                sum_table[[i, j]] = val 
                    + sum_table[[i - 1, j]] 
                    + sum_table[[i, j - 1]] 
                    - sum_table[[i - 1, j - 1]];

                sum_sq_table[[i, j]] = val_sq 
                    + sum_sq_table[[i - 1, j]] 
                    + sum_sq_table[[i, j - 1]] 
                    - sum_sq_table[[i - 1, j - 1]];

                count_table[[i, j]] = cnt 
                    + count_table[[i - 1, j]] 
                    + count_table[[i, j - 1]] 
                    - count_table[[i - 1, j - 1]];
            }
        }

        Self {
            sum_table,
            sum_sq_table,
            count_table,
        }
    }

    /// Calculate window mean and variance in O(1) time
    /// Returns (mean, variance) for the specified window
    fn window_stats(&self, row1: usize, col1: usize, row2: usize, col2: usize) -> (f32, f32) {
        // Ensure bounds are valid
        if row2 < row1 || col2 < col1 {
            return (0.0, 0.0);
        }
        
        // Ensure indices are within bounds
        let (height, width) = (self.sum_table.nrows() - 1, self.sum_table.ncols() - 1);
        if row2 >= height || col2 >= width {
            return (0.0, 0.0);
        }

        let sum = self.sum_table[[row2 + 1, col2 + 1]]
            - self.sum_table[[row1, col2 + 1]]
            - self.sum_table[[row2 + 1, col1]]
            + self.sum_table[[row1, col1]];

        let sum_sq = self.sum_sq_table[[row2 + 1, col2 + 1]]
            - self.sum_sq_table[[row1, col2 + 1]]
            - self.sum_sq_table[[row2 + 1, col1]]
            + self.sum_sq_table[[row1, col1]];

        // Use i64 arithmetic to avoid overflow/underflow
        let count = (self.count_table[[row2 + 1, col2 + 1]] as i64)
            - (self.count_table[[row1, col2 + 1]] as i64)
            - (self.count_table[[row2 + 1, col1]] as i64)
            + (self.count_table[[row1, col1]] as i64);
        
        let count = count.max(0) as u32;

        if count == 0 {
            return (0.0, 0.0);
        }

        let n = count as f64;
        let mean = sum / n;
        let variance = if count > 1 {
            (sum_sq / n - mean * mean).max(0.0)
        } else {
            0.0
        };

        (mean as f32, variance as f32)
    }
}

/// Speckle filtering parameters
#[derive(Debug, Clone)]
pub struct SpeckleFilterParams {
    /// Filter window size (must be odd)
    pub window_size: usize,
    /// Number of looks (for adaptive filters)
    pub num_looks: f32,
    /// Threshold for edge detection (Lee sigma filter)
    pub edge_threshold: f32,
    /// Damping factor (for Gamma MAP filter)
    pub damping_factor: f32,
    /// Coefficient of variation threshold
    pub cv_threshold: f32,
    /// Tile size for cache-optimized processing (0 = auto)
    pub tile_size: usize,
}

impl Default for SpeckleFilterParams {
    fn default() -> Self {
        Self {
            window_size: 7,          // 7x7 window
            num_looks: 1.0,          // Single look
            edge_threshold: 0.5,     // Edge detection threshold
            damping_factor: 1.0,     // No damping
            cv_threshold: 0.5,       // Standard CV threshold
            tile_size: 0,            // Auto tile size selection
        }
    }
}

/// Tile size configuration for cache optimization
#[derive(Debug, Clone, Copy)]
pub enum TileSize {
    /// Small tiles for better cache locality
    Small = 64,
    /// Medium tiles for balanced performance
    Medium = 128,
    /// Large tiles for reduced overhead
    Large = 256,
    /// Auto-select based on image size
    Auto = 0,
}

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
    params: SpeckleFilterParams,
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

        // Validate window size is odd
        if self.params.window_size % 2 == 0 {
            return Err(SarError::Processing(
                "Window size must be odd".to_string()
            ));
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

    /// Apply speckle filtering with tiled processing for optimal cache performance
    pub fn apply_filter_tiled(
        &self,
        image: &Array2<f32>,
        filter_type: SpeckleFilterType,
        tile_size: Option<usize>,
    ) -> SarResult<Array2<f32>> {
        log::info!("Applying {:?} speckle filter with tiled processing", filter_type);
        
        let (height, width) = image.dim();
        let effective_tile_size = self.determine_optimal_tile_size(height, width, tile_size);
        
        log::debug!("Using tile size: {}x{}", effective_tile_size, effective_tile_size);

        // For small images, use standard processing
        if height <= effective_tile_size && width <= effective_tile_size {
            return self.apply_filter(image, filter_type);
        }

        let mut result = Array2::zeros((height, width));
        let half_window = self.params.window_size / 2;

        // Build integral image if needed for large windows
        let integral = if self.params.window_size >= 9 {
            Some(IntegralImage::new(image))
        } else {
            None
        };

        // Process image in tiles with overlap for boundary handling
        for tile_row in (0..height).step_by(effective_tile_size) {
            for tile_col in (0..width).step_by(effective_tile_size) {
                let tile_end_row = (tile_row + effective_tile_size).min(height);
                let tile_end_col = (tile_col + effective_tile_size).min(width);

                // Extract tile with padding for window operations
                let padded_start_row = tile_row.saturating_sub(half_window);
                let padded_end_row = (tile_end_row + half_window).min(height);
                let padded_start_col = tile_col.saturating_sub(half_window);
                let padded_end_col = (tile_end_col + half_window).min(width);

                let tile = image.slice(s![
                    padded_start_row..padded_end_row,
                    padded_start_col..padded_end_col
                ]).to_owned();

                // Process tile
                let filtered_tile = self.apply_filter_to_tile(&tile, filter_type, integral.as_ref())?;

                // Calculate offsets for copying results back
                let copy_start_row = half_window.min(tile_row - padded_start_row);
                let copy_start_col = half_window.min(tile_col - padded_start_col);
                let copy_height = tile_end_row - tile_row;
                let copy_width = tile_end_col - tile_col;

                // Copy results back to main result array
                result.slice_mut(s![
                    tile_row..tile_end_row,
                    tile_col..tile_end_col
                ]).assign(&filtered_tile.slice(s![
                    copy_start_row..copy_start_row + copy_height,
                    copy_start_col..copy_start_col + copy_width
                ]));
            }
        }

        log::info!("Tiled speckle filtering completed successfully");
        Ok(result)
    }

    /// Determine optimal tile size based on image dimensions and cache considerations
    fn determine_optimal_tile_size(&self, height: usize, width: usize, user_tile_size: Option<usize>) -> usize {
        if let Some(size) = user_tile_size {
            return size;
        }

        if self.params.tile_size > 0 {
            return self.params.tile_size;
        }

        // Auto-select tile size based on image size and memory considerations
        let total_pixels = height * width;
        
        if total_pixels < 1_000_000 {
            // Small images: use smaller tiles for better cache locality
            TileSize::Small as usize
        } else if total_pixels < 10_000_000 {
            // Medium images: balanced approach
            TileSize::Medium as usize
        } else {
            // Large images: use larger tiles to reduce overhead
            TileSize::Large as usize
        }
    }

    /// Apply filter to a single tile with optimized processing
    fn apply_filter_to_tile(
        &self,
        tile: &Array2<f32>,
        filter_type: SpeckleFilterType,
        integral: Option<&IntegralImage>,
    ) -> SarResult<Array2<f32>> {
        let (tile_height, tile_width) = tile.dim();
        let mut filtered = Array2::zeros((tile_height, tile_width));
        let half_window = self.params.window_size / 2;

        match filter_type {
            SpeckleFilterType::Lee => {
                self.apply_lee_filter_to_tile(tile, &mut filtered, integral)
            },
            SpeckleFilterType::EnhancedLee => {
                self.apply_enhanced_lee_filter_to_tile(tile, &mut filtered, integral)
            },
            SpeckleFilterType::GammaMAP => {
                self.apply_gamma_map_filter_to_tile(tile, &mut filtered, integral)
            },
            SpeckleFilterType::Mean => {
                self.apply_mean_filter_to_tile(tile, &mut filtered, integral)
            },
            SpeckleFilterType::Median => {
                self.apply_median_filter_to_tile(tile, &mut filtered)
            },
            _ => {
                // For other filters, use standard implementation
                return self.apply_filter(tile, filter_type);
            }
        }?;

        Ok(filtered)
    }

    /// Apply Lee filter to a tile with cache optimization
    fn apply_lee_filter_to_tile(
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
    fn apply_enhanced_lee_filter_to_tile(
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

    /// Apply Gamma-MAP filter to a tile
    fn apply_gamma_map_filter_to_tile(
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

    /// Apply mean filter to a tile
    fn apply_mean_filter_to_tile(
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

    /// Apply median filter to a tile
    fn apply_median_filter_to_tile(
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
                    filtered[[i, j]] = Self::histogram_median_8bit(tile, i, j, self.params.window_size);
                } else {
                    let mut window_values = Vec::new();

                    let i_start = i.saturating_sub(half_window);
                    let i_end = (i + half_window + 1).min(height);
                    let j_start = j.saturating_sub(half_window);
                    let j_end = (j + half_window + 1).min(width);

                    for wi in i_start..i_end {
                        for wj in j_start..j_end {
                            let pixel_val = tile[[wi, wj]];
                            if pixel_val.is_finite() && pixel_val > 0.0 {
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
                            if pixel_val.is_finite() && pixel_val > 0.0 {
                                sum += pixel_val;
                                count += 1;
                            }
                        }
                    }
                }

                filtered[[i, j]] = if count > 0 { sum / count as f32 } else { image[[i, j]] };
            }
        }

        Ok(filtered)
    }

    /// Enhanced Lee filter with chunked processing for better cache performance
    fn apply_enhanced_lee_filter_chunked(&self, image: &Array2<f32>) -> SarResult<Array2<f32>> {
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

                        let (local_mean, local_variance) = self.calculate_local_statistics_optimized(image, i, j, half_window, integral.as_ref());
                        
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
    fn apply_enhanced_lee_filter_parallel(&self, image: &Array2<f32>) -> SarResult<Array2<f32>> {
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
        let target_chunk_size = (total_pixels / (num_threads * 4)).max(1024); // At least 1KB per chunk
        
        // Create row-based chunks for better cache locality
        let rows_per_chunk = (target_chunk_size / width).max(1).min(height);
        
        log::debug!("Parallel processing: {} threads, {} rows per chunk", num_threads, rows_per_chunk);

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
                            let (local_mean, local_variance) = if let Some(ref integral_img) = integral {
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
    fn apply_enhanced_lee_filter_parallel(&self, image: &Array2<f32>) -> SarResult<Array2<f32>> {
        // Fallback to chunked processing if parallel feature is not available
        self.apply_enhanced_lee_filter_chunked(image)
    }

    /// NUMA-aware parallel processing for very large images
    #[cfg(feature = "parallel")]
    pub fn apply_filter_numa_optimized(
        &self,
        image: &Array2<f32>,
        filter_type: SpeckleFilterType,
    ) -> SarResult<Array2<f32>> {
        use rayon::prelude::*;
        
        log::info!("Applying {:?} speckle filter with NUMA-optimized parallel processing", filter_type);
        
        let (height, width) = image.dim();
        let total_pixels = height * width;
        
        // Use NUMA optimization only for very large images
        if total_pixels < 50_000_000 {
            return self.apply_filter_tiled(image, filter_type, None);
        }
        
        let num_threads = rayon::current_num_threads();
        let numa_regions = num_threads.min(8); // Limit NUMA regions
        let rows_per_region = height / numa_regions;
        
        log::debug!("NUMA processing: {} regions, {} rows per region", numa_regions, rows_per_region);
        
        // Process in NUMA-aware regions
        let region_results: Vec<_> = (0..numa_regions)
            .into_par_iter()
            .map(|region_id| {
                let start_row = region_id * rows_per_region;
                let end_row = if region_id == numa_regions - 1 {
                    height
                } else {
                    (region_id + 1) * rows_per_region
                };
                
                // Extract region with padding for window operations
                let half_window = self.params.window_size / 2;
                let padded_start = start_row.saturating_sub(half_window);
                let padded_end = (end_row + half_window).min(height);
                
                let region_slice = image.slice(s![padded_start..padded_end, ..]).to_owned();
                
                // Apply filter to region
                let filtered_region = self.apply_filter_tiled(&region_slice, filter_type, Some(256))?;
                
                // Calculate copy parameters
                let copy_start = half_window.min(start_row - padded_start);
                let copy_height = end_row - start_row;
                
                Ok::<_, SarError>((start_row, copy_start, copy_height, filtered_region))
            })
            .collect::<Result<Vec<_>, _>>()?;
        
        // Combine results
        let mut result = Array2::zeros((height, width));
        for (start_row, copy_start, copy_height, region_data) in region_results {
            result.slice_mut(s![start_row..start_row + copy_height, ..])
                .assign(&region_data.slice(s![copy_start..copy_start + copy_height, ..]));
        }
        
        log::info!("NUMA-optimized parallel processing completed successfully");
        Ok(result)
    }

    #[cfg(not(feature = "parallel"))]
    pub fn apply_filter_numa_optimized(
        &self,
        image: &Array2<f32>,
        filter_type: SpeckleFilterType,
    ) -> SarResult<Array2<f32>> {
        // Fallback to tiled processing if parallel feature is not available
        self.apply_filter_tiled(image, filter_type, None)
    }

    /// Apply median filter
    fn apply_median_filter(&self, image: &Array2<f32>) -> SarResult<Array2<f32>> {
        log::debug!("Applying median filter");
        
        let (height, width) = image.dim();
        let mut filtered = Array2::zeros((height, width));
        let half_window = self.params.window_size / 2;

        for i in 0..height {
            for j in 0..width {
                let mut window_values = Vec::new();

                // Collect window values
                for wi in 0..self.params.window_size {
                    for wj in 0..self.params.window_size {
                        let ii = i as i32 + wi as i32 - half_window as i32;
                        let jj = j as i32 + wj as i32 - half_window as i32;

                        if ii >= 0 && ii < height as i32 && jj >= 0 && jj < width as i32 {
                            let pixel_val = image[[ii as usize, jj as usize]];
                            if pixel_val.is_finite() && pixel_val > 0.0 {
                                window_values.push(pixel_val);
                            }
                        }
                    }
                }

                if !window_values.is_empty() {
                    window_values.sort_by(|a, b| a.partial_cmp(b).unwrap());
                    filtered[[i, j]] = window_values[window_values.len() / 2];
                } else {
                    filtered[[i, j]] = image[[i, j]];
                }
            }
        }

        Ok(filtered)
    }

    /// Fast median computation using quickselect algorithm for floating-point data
    /// O(n) average case vs O(n log n) for full sorting
    fn quickselect_median(data: &mut [f32]) -> f32 {
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
    fn histogram_median_8bit(image: &Array2<f32>, center_i: usize, center_j: usize, window_size: usize) -> f32 {
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
                    let bin = pixel_val as usize;
                    if bin < 256 {
                        histogram[bin] += 1;
                        total_pixels += 1;
                    }
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

    /// Optimized mean filter with better memory access patterns
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
                        let (mean, _) = self.calculate_local_statistics_fast(image, i, j, half_window);
                        filtered[[i, j]] = if mean > 0.0 { mean } else { image[[i, j]] };
                    }
                }
            }
        }

        Ok(filtered)
    }

    /// Optimized median filter with efficient sorting
    fn apply_median_filter_optimized(&self, image: &Array2<f32>) -> SarResult<Array2<f32>> {
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
                    filtered[[i, j]] = Self::histogram_median_8bit(image, i, j, self.params.window_size);
                } else {
                    // Use quickselect median for floating-point images
                    let mut window_values = Vec::with_capacity(self.params.window_size * self.params.window_size);

                    let i_start = i.saturating_sub(half_window);
                    let i_end = (i + half_window + 1).min(height);
                    let j_start = j.saturating_sub(half_window);
                    let j_end = (j + half_window + 1).min(width);

                    for wi in i_start..i_end {
                        for wj in j_start..j_end {
                            let pixel_val = image[[wi, wj]];
                            if pixel_val.is_finite() && pixel_val > 0.0 {
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

    /// Check if image appears to be 8-bit intensity data
    fn is_8bit_intensity(image: &Array2<f32>) -> bool {
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

    /// Apply Lee filter (adaptive)
    /// Implements the adaptive Lee filter based on local statistics
    /// Reference: Lee, J.S. (1980). Digital Image Enhancement and Noise Filtering by Use of Local Statistics
    /// CORRECTED: Uses proper literature-based weight formula for scientific accuracy
    fn apply_lee_filter(&self, image: &Array2<f32>) -> SarResult<Array2<f32>> {
        let (height, width) = image.dim();
        let mut filtered = Array2::zeros((height, width));
        let half_window = self.params.window_size / 2;

        log::debug!("Applying Lee filter with window size {}", self.params.window_size);

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
                let (local_mean, local_variance) = self.calculate_local_statistics_optimized(image, i, j, half_window, integral.as_ref());

                if local_mean <= 0.0 {
                    filtered[[i, j]] = center_value;
                    continue;
                }

                // Lee filter coefficient of variation
                let cv = local_variance.sqrt() / local_mean;
                let cu = 1.0 / self.params.num_looks.sqrt(); // Theoretical CV for fully developed speckle

                // Epsilon for numerical stability near cu
                const EPSILON: f32 = 1e-6;

                // Apply CORRECTED Lee filter formula
                let result = if cv <= cu + EPSILON {
                    // Homogeneous area - use local mean
                    local_mean
                } else {
                    // Heterogeneous area - weighted combination
                    // CORRECTED FORMULA: weight = (cv² - cu²) / (cv² + cu²)
                    // Previous incorrect: weight = 1.0 - cu²/cv²
                    let weight = (cv * cv - cu * cu) / (cv * cv + cu * cu);
                    local_mean + weight * (center_value - local_mean)
                };

                filtered[[i, j]] = result.max(0.0); // Ensure non-negative
            }
        }

        // Handle borders by copying original values
        self.handle_borders(image, &mut filtered, half_window);

        Ok(filtered)
    }

    /// Apply Enhanced Lee filter with edge detection
    /// More sophisticated version that preserves edges better
    /// CORRECTED: Implements robust numerical handling for edge cases
    fn apply_enhanced_lee_filter(&self, image: &Array2<f32>) -> SarResult<Array2<f32>> {
        let (height, width) = image.dim();
        let mut filtered = Array2::zeros((height, width));
        let half_window = self.params.window_size / 2;

        log::debug!("Applying Enhanced Lee filter with window size {}", self.params.window_size);

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
                let (local_mean, local_variance) = self.calculate_local_statistics_optimized(image, i, j, half_window, integral.as_ref());

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

    /// Apply Gamma MAP filter (Maximum A Posteriori)
    /// Optimal filter based on Gamma distribution assumption for SAR speckle
    /// Reference: Lopes et al. (1993). Adaptive speckle filters and scene heterogeneity
    /// CORRECTED: Uses bounded weight formulation to prevent numerical instabilities
    fn apply_gamma_map_filter(&self, image: &Array2<f32>) -> SarResult<Array2<f32>> {
        let (height, width) = image.dim();
        let mut filtered = Array2::zeros((height, width));
        let half_window = self.params.window_size / 2;

        log::debug!("Applying Gamma MAP filter with window size {}", self.params.window_size);

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
                let (local_mean, local_variance) = self.calculate_local_statistics_optimized(image, i, j, half_window, integral.as_ref());

                if local_mean <= 0.0 {
                    filtered[[i, j]] = center_value;
                    continue;
                }

                // Gamma distribution parameters
                let cv = local_variance.sqrt() / local_mean;
                let cu = 1.0 / self.params.num_looks.sqrt();

                // Epsilon for numerical stability
                const EPSILON: f32 = 1e-6;

                // CORRECTED: Bounded Gamma-MAP estimation
                let result = if cv <= cu + EPSILON {
                    // Homogeneous area - use local mean
                    local_mean
                } else {
                    // Heterogeneous area - bounded MAP estimate
                    // CORRECTED FORMULA: weight = (cv² - cu²) / (cv² * (1 + cu²))
                    // This prevents unbounded growth and numerical instabilities
                    let weight = (cv * cv - cu * cu) / (cv * cv * (1.0 + cu * cu));
                    local_mean + weight * (center_value - local_mean)
                };

                filtered[[i, j]] = result.max(0.0);
            }
        }

        self.handle_borders(image, &mut filtered, half_window);
        Ok(filtered)
    }

    /// Apply Lee Sigma filter (edge-preserving)
    /// Uses sigma-based edge detection to preserve linear features
    fn apply_lee_sigma_filter(&self, image: &Array2<f32>) -> SarResult<Array2<f32>> {
        let (height, width) = image.dim();
        let mut filtered = Array2::zeros((height, width));
        let half_window = self.params.window_size / 2;

        log::debug!("Applying Lee Sigma filter with window size {}", self.params.window_size);

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
                    if *value >= lower_bound && *value <= upper_bound && value.is_finite() && *value > 0.0 {
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

    /// Apply Frost filter (exponential weighting)
    /// Applies exponential weighting based on distance and local statistics
    /// Reference: Frost et al. (1982). A model for radar images and its application to adaptive digital filtering
    fn apply_frost_filter(&self, image: &Array2<f32>) -> SarResult<Array2<f32>> {
        let (height, width) = image.dim();
        let mut filtered = Array2::zeros((height, width));
        let half_window = self.params.window_size / 2;

        log::debug!("Applying Frost filter with window size {}", self.params.window_size);

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
                        if value.is_finite() && value > 0.0 {
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

    /// Apply Refined Lee filter with improved edge detection
    /// CORRECTED: Includes robust numerical handling for edge cases
    fn apply_refined_lee_filter(&self, image: &Array2<f32>) -> SarResult<Array2<f32>> {
        let (height, width) = image.dim();
        let mut filtered = Array2::zeros((height, width));
        let half_window = self.params.window_size / 2;

        log::debug!("Applying Refined Lee filter with window size {}", self.params.window_size);

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
                let (local_mean, local_variance) = self.calculate_local_statistics_optimized(image, i, j, half_window, integral.as_ref());

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

    /// Calculate local statistics for a window
    #[allow(dead_code)]
    fn calculate_local_statistics(&self, image: &Array2<f32>, center_i: usize, center_j: usize, half_window: usize) -> (f32, f32) {
        let (height, width) = image.dim();
        let mut values = Vec::new();

        for wi in 0..self.params.window_size {
            for wj in 0..self.params.window_size {
                let ii = center_i as i32 + wi as i32 - half_window as i32;
                let jj = center_j as i32 + wj as i32 - half_window as i32;

                if ii >= 0 && ii < height as i32 && jj >= 0 && jj < width as i32 {
                    let pixel_val = image[[ii as usize, jj as usize]];
                    if pixel_val.is_finite() && pixel_val > 0.0 {
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
            values.iter()
                .map(|v| (v - mean) * (v - mean))
                .sum::<f32>() / (values.len() - 1) as f32
        } else {
            0.0
        };

        (mean, variance)
    }

    /// Fast local statistics calculation with optimized loop
    fn calculate_local_statistics_fast(&self, image: &Array2<f32>, center_i: usize, center_j: usize, half_window: usize) -> (f32, f32) {
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
                if pixel_val.is_finite() && pixel_val > 0.0 {
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
    fn calculate_local_statistics_optimized(&self, image: &Array2<f32>, center_i: usize, center_j: usize, half_window: usize, integral: Option<&IntegralImage>) -> (f32, f32) {
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

    /// Multi-scale speckle filtering
    pub fn apply_multiscale_filter(
        &self,
        image: &Array2<f32>,
        filter_type: SpeckleFilterType,
        scales: &[usize],
    ) -> SarResult<Array2<f32>> {
        log::info!("Applying multi-scale speckle filtering with {} scales", scales.len());
        
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
    pub fn estimate_number_of_looks(image: &Array2<f32>) -> SarResult<f32> {
        log::debug!("Estimating number of looks from image statistics");
        
        let (height, width) = image.dim();
        let mut values = Vec::new();
        
        // Collect valid pixel values
        for i in 0..height {
            for j in 0..width {
                let val = image[[i, j]];
                if val.is_finite() && val > 0.0 {
                    values.push(val);
                }
            }
        }
        
        if values.is_empty() {
            return Err(SarError::Processing("No valid pixels found".to_string()));
        }
        
        // Calculate mean and variance
        let mean = values.iter().sum::<f32>() / values.len() as f32;
        let variance = values.iter()
            .map(|v| (v - mean) * (v - mean))
            .sum::<f32>() / (values.len() - 1) as f32;
        
        // Number of looks = mean² / variance
        let num_looks = (mean * mean) / variance;
        
        log::info!("Estimated number of looks: {:.2}", num_looks);
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
        log::info!("Standard Enhanced Lee: {:.2} seconds ({:.0} pixels/sec)", 
                  standard_time.as_secs_f64(), 
                  total_pixels as f64 / standard_time.as_secs_f64());
        
        // Chunked implementation
        let start = Instant::now();
        let _result2 = self.apply_enhanced_lee_filter_chunked(image)?;
        let chunked_time = start.elapsed();
        log::info!("Chunked Enhanced Lee: {:.2} seconds ({:.0} pixels/sec)", 
                  chunked_time.as_secs_f64(),
                  total_pixels as f64 / chunked_time.as_secs_f64());
        
        // Parallel implementation (if available)
        #[cfg(feature = "parallel")]
        {
            let start = Instant::now();
            let _result3 = self.apply_enhanced_lee_filter_parallel(image)?;
            let parallel_time = start.elapsed();
            log::info!("Parallel Enhanced Lee: {:.2} seconds ({:.0} pixels/sec)", 
                      parallel_time.as_secs_f64(),
                      total_pixels as f64 / parallel_time.as_secs_f64());
            
            let parallel_speedup = standard_time.as_secs_f64() / parallel_time.as_secs_f64();
            log::info!("Parallel speedup: {:.2}x", parallel_speedup);
        }
        
        let chunked_speedup = standard_time.as_secs_f64() / chunked_time.as_secs_f64();
        log::info!("Chunked speedup: {:.2}x", chunked_speedup);
        log::info!("=== Benchmark Complete ===");
        
        Ok(())
    }

    /// Calculate window mean for a 2D array slice
    fn calculate_window_mean(&self, window: &ndarray::ArrayView2<f32>) -> f32 {
        let mut sum = 0.0;
        let mut count = 0;
        
        for value in window.iter() {
            if value.is_finite() && *value > 0.0 {
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
    fn calculate_window_variance(&self, window: &ndarray::ArrayView2<f32>, mean: f32) -> f32 {
        let mut sum_sq_diff = 0.0;
        let mut count = 0;
        
        for value in window.iter() {
            if value.is_finite() && *value > 0.0 {
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
    fn handle_borders(&self, original: &Array2<f32>, filtered: &mut Array2<f32>, border_size: usize) {
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
    fn calculate_edge_strength(&self, image: &Array2<f32>, i: usize, j: usize) -> f32 {
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
    fn calculate_adaptive_statistics(&self, window: &ndarray::ArrayView2<f32>, _center_i: usize, _center_j: usize) -> (f32, f32) {
        // For now, use standard statistics
        // In a full implementation, this would include directional analysis
        let mean = self.calculate_window_mean(window);
        let variance = self.calculate_window_variance(window, mean);
        (mean, variance)
    }
}

impl Default for SpeckleFilter {
    fn default() -> Self {
        Self::new()
    }
}
