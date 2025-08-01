use crate::types::{SarError, SarResult};
use ndarray::{Array2, s};

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
}

impl Default for SpeckleFilterParams {
    fn default() -> Self {
        Self {
            window_size: 7,          // 7x7 window
            num_looks: 1.0,          // Single look
            edge_threshold: 0.5,     // Edge detection threshold
            damping_factor: 1.0,     // No damping
            cv_threshold: 0.5,       // Standard CV threshold
        }
    }
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

    /// Apply speckle filtering with optimized algorithm selection
    pub fn apply_filter_optimized(
        &self,
        image: &Array2<f32>,
        filter_type: SpeckleFilterType,
    ) -> SarResult<Array2<f32>> {
        log::info!("Applying OPTIMIZED {:?} speckle filter", filter_type);
        log::debug!("Filter parameters: {:?}", self.params);

        let (height, width) = image.dim();
        let total_pixels = height * width;

        // Choose optimization strategy based on image size
        let filtered = match filter_type {
            SpeckleFilterType::EnhancedLee => {
                if total_pixels > 10_000_000 {  // > 10M pixels
                    self.apply_enhanced_lee_filter_parallel(image)?
                } else if total_pixels > 1_000_000 {  // 1-10M pixels
                    self.apply_enhanced_lee_filter_chunked(image)?
                } else {
                    self.apply_enhanced_lee_filter(image)?  // Use standard for small images
                }
            },
            SpeckleFilterType::Mean => {
                if total_pixels > 5_000_000 {
                    self.apply_mean_filter_optimized(image)?
                } else {
                    self.apply_mean_filter(image)?
                }
            },
            SpeckleFilterType::Median => {
                if total_pixels > 5_000_000 {
                    self.apply_median_filter_optimized(image)?
                } else {
                    self.apply_median_filter(image)?
                }
            },
            _ => {
                // For other filters, use standard implementation
                self.apply_filter(image, filter_type)?
            }
        };

        log::info!("Optimized speckle filtering completed successfully");
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

                        let (local_mean, local_variance) = self.calculate_local_statistics_fast(image, i, j, half_window);
                        
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
        
        log::debug!("Applying Enhanced Lee filter with parallel processing");
        
        let (height, width) = image.dim();
        let mut filtered = Array2::zeros((height, width));
        let half_window = self.params.window_size / 2;

        let cu = 1.0 / self.params.num_looks.sqrt();
        let cmax = 1.73;

        // Create index pairs for parallel processing
        let indices: Vec<(usize, usize)> = (0..height)
            .flat_map(|i| (0..width).map(move |j| (i, j)))
            .collect();

        let results: Vec<(usize, usize, f32)> = indices
            .into_par_iter()
            .map(|(i, j)| {
                let center_value = image[[i, j]];
                
                let result = if !center_value.is_finite() || center_value <= 0.0 {
                    center_value
                } else {
                    let (local_mean, local_variance) = self.calculate_local_statistics_fast(image, i, j, half_window);
                    
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
                
                (i, j, result)
            })
            .collect();

        // Assign results back to filtered array
        for (i, j, value) in results {
            filtered[[i, j]] = value;
        }

        Ok(filtered)
    }

    #[cfg(not(feature = "parallel"))]
    fn apply_enhanced_lee_filter_parallel(&self, image: &Array2<f32>) -> SarResult<Array2<f32>> {
        // Fallback to chunked processing if parallel feature is not available
        self.apply_enhanced_lee_filter_chunked(image)
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

        for i in 0..height {
            for j in 0..width {
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
                    // Sort and find median
                    window_values.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
                    let mid = window_values.len() / 2;
                    filtered[[i, j]] = window_values[mid];
                } else {
                    filtered[[i, j]] = image[[i, j]];
                }
            }
        }

        Ok(filtered)
    }

    /// Apply Lee filter (adaptive)
    /// Implements the adaptive Lee filter based on local statistics
    /// Reference: Lee, J.S. (1980). Digital Image Enhancement and Noise Filtering by Use of Local Statistics
    fn apply_lee_filter(&self, image: &Array2<f32>) -> SarResult<Array2<f32>> {
        let (height, width) = image.dim();
        let mut filtered = Array2::zeros((height, width));
        let half_window = self.params.window_size / 2;

        log::debug!("Applying Lee filter with window size {}", self.params.window_size);

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

                // Calculate local statistics
                let local_mean = self.calculate_window_mean(&window);
                let local_variance = self.calculate_window_variance(&window, local_mean);

                if local_mean <= 0.0 {
                    filtered[[i, j]] = center_value;
                    continue;
                }

                // Lee filter coefficient of variation
                let cv = local_variance.sqrt() / local_mean;
                let cu = 1.0 / self.params.num_looks.sqrt(); // Theoretical CV for fully developed speckle

                // Apply Lee filter formula
                let result = if cv <= cu {
                    // Homogeneous area - use local mean
                    local_mean
                } else {
                    // Heterogeneous area - weighted combination
                    let weight = 1.0 - cu * cu / (cv * cv);
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
    fn apply_enhanced_lee_filter(&self, image: &Array2<f32>) -> SarResult<Array2<f32>> {
        let (height, width) = image.dim();
        let mut filtered = Array2::zeros((height, width));
        let half_window = self.params.window_size / 2;

        log::debug!("Applying Enhanced Lee filter with window size {}", self.params.window_size);

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

                // Calculate directional statistics for edge detection
                let (local_mean, local_variance) = self.calculate_adaptive_statistics(&window, i, j);

                if local_mean <= 0.0 {
                    filtered[[i, j]] = center_value;
                    continue;
                }

                // Enhanced Lee with coefficient of variation limits
                let cv = local_variance.sqrt() / local_mean;
                let cu = 1.0 / self.params.num_looks.sqrt();
                let cmax = 1.73; // Maximum CV for speckle noise

                let result = if cv <= cu {
                    // Homogeneous area
                    local_mean
                } else if cv < cmax {
                    // Heterogeneous area - apply Enhanced Lee
                    let weight = (cv * cv - cu * cu) / (cv * cv * (1.0 + cu * cu));
                    local_mean + weight * (center_value - local_mean)
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
    fn apply_gamma_map_filter(&self, image: &Array2<f32>) -> SarResult<Array2<f32>> {
        let (height, width) = image.dim();
        let mut filtered = Array2::zeros((height, width));
        let half_window = self.params.window_size / 2;

        log::debug!("Applying Gamma MAP filter with window size {}", self.params.window_size);

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

                // Gamma distribution parameters
                let cv = local_variance.sqrt() / local_mean;
                let cu = 1.0 / self.params.num_looks.sqrt();

                // MAP estimation for Gamma distribution
                let alpha = (1.0 + cu * cu) / (cv * cv - cu * cu).max(1e-10);
                let k = alpha / local_mean;

                let result = if cv <= cu {
                    // Homogeneous area
                    local_mean
                } else {
                    // MAP estimate
                    (alpha - 1.0) / k.max(1e-10)
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
    fn apply_refined_lee_filter(&self, image: &Array2<f32>) -> SarResult<Array2<f32>> {
        let (height, width) = image.dim();
        let mut filtered = Array2::zeros((height, width));
        let half_window = self.params.window_size / 2;

        log::debug!("Applying Refined Lee filter with window size {}", self.params.window_size);

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

                // Apply standard Lee filter for non-edge areas
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

                let cv = local_variance.sqrt() / local_mean;
                let cu = 1.0 / self.params.num_looks.sqrt();

                let result = if cv <= cu {
                    local_mean
                } else {
                    let weight = 1.0 - cu * cu / (cv * cv);
                    local_mean + weight * (center_value - local_mean)
                };

                filtered[[i, j]] = result.max(0.0);
            }
        }

        self.handle_borders(image, &mut filtered, half_window);
        Ok(filtered)
    }

    /// Calculate local statistics for a window
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
