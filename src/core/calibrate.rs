use crate::types::{SarError, SarImage, SarRealImage, SarResult, Polarization};
use ndarray::{Array1, Array2, Zip, Axis};
use regex::Regex;
use rayon::prelude::*;

/// Calibration vector from Sentinel-1 XML
#[derive(Debug, Clone)]
pub struct CalibrationVector {
    pub azimuth_time: String,
    pub line: usize,
    pub pixels: Vec<usize>,
    pub sigma_nought: Vec<f32>,
    pub beta_nought: Vec<f32>,
    pub gamma: Vec<f32>,
    pub dn: Vec<f32>,
}

/// Pre-computed calibration lookup table for fast access
#[derive(Debug, Clone)]
pub struct CalibrationLUT {
    pub sigma_values: Array2<f32>,
    pub beta_values: Array2<f32>, 
    pub gamma_values: Array2<f32>,
    pub dn_values: Array2<f32>,
    pub is_precomputed: bool,
}

impl CalibrationLUT {
    pub fn new(dims: (usize, usize)) -> Self {
        let (height, width) = dims;
        Self {
            sigma_values: Array2::zeros((height, width)),
            beta_values: Array2::zeros((height, width)),
            gamma_values: Array2::zeros((height, width)), 
            dn_values: Array2::ones((height, width)),
            is_precomputed: false,
        }
    }
}

/// Calibration coefficients for radiometric correction
#[derive(Debug, Clone)]
pub struct CalibrationCoefficients {
    pub vectors: Vec<CalibrationVector>,
    pub swath: String,
    pub polarization: String,
    pub product_first_line_utc_time: String,
    pub product_last_line_utc_time: String,
    pub lut: Option<CalibrationLUT>,
}

impl CalibrationCoefficients {
    pub fn new() -> Self {
        Self {
            vectors: Vec::new(),
            swath: String::new(),
            polarization: String::new(),
            product_first_line_utc_time: String::new(),
            product_last_line_utc_time: String::new(),
            lut: None,
        }
    }

    /// Pre-compute calibration lookup table for entire image (MAJOR OPTIMIZATION)
    pub fn precompute_lut(&mut self, image_dims: (usize, usize)) -> SarResult<()> {
        log::info!("Pre-computing calibration LUT for {}x{} image", image_dims.0, image_dims.1);
        let start_time = std::time::Instant::now();
        
        let (height, width) = image_dims;
        let mut lut = CalibrationLUT::new(image_dims);
        
        // Sort vectors by line for faster access
        self.vectors.sort_by_key(|v| v.line);
        
        // More efficient LUT computation: only compute the values we need
        // Use chunked parallel processing to reduce memory pressure
        let chunk_size = std::cmp::min(1024, height / 4); // Process in 1024-line chunks max
        
        log::debug!("Using chunk size of {} lines for LUT computation", chunk_size);
        
        for chunk_start in (0..height).step_by(chunk_size) {
            let chunk_end = std::cmp::min(chunk_start + chunk_size, height);
            
            // Process this chunk in parallel
            let mut sigma_chunk = Array2::zeros((chunk_end - chunk_start, width));
            let mut beta_chunk = Array2::zeros((chunk_end - chunk_start, width));
            let mut gamma_chunk = Array2::zeros((chunk_end - chunk_start, width));
            
            sigma_chunk.axis_iter_mut(Axis(0))
                .into_par_iter()
                .enumerate()
                .for_each(|(local_i, mut row)| {
                    let global_i = chunk_start + local_i;
                    for j in 0..width {
                        if let Ok(val) = self.get_calibration_value_fast(global_i, j, CalibrationType::Sigma0) {
                            row[j] = val;
                        }
                    }
                });
            
            beta_chunk.axis_iter_mut(Axis(0))
                .into_par_iter()
                .enumerate()
                .for_each(|(local_i, mut row)| {
                    let global_i = chunk_start + local_i;
                    for j in 0..width {
                        if let Ok(val) = self.get_calibration_value_fast(global_i, j, CalibrationType::Beta0) {
                            row[j] = val;
                        }
                    }
                });
                
            gamma_chunk.axis_iter_mut(Axis(0))
                .into_par_iter()
                .enumerate()
                .for_each(|(local_i, mut row)| {
                    let global_i = chunk_start + local_i;
                    for j in 0..width {
                        if let Ok(val) = self.get_calibration_value_fast(global_i, j, CalibrationType::Gamma0) {
                            row[j] = val;
                        }
                    }
                });
            
            // Copy chunks to main LUT
            for local_i in 0..(chunk_end - chunk_start) {
                let global_i = chunk_start + local_i;
                for j in 0..width {
                    lut.sigma_values[[global_i, j]] = sigma_chunk[[local_i, j]];
                    lut.beta_values[[global_i, j]] = beta_chunk[[local_i, j]];
                    lut.gamma_values[[global_i, j]] = gamma_chunk[[local_i, j]];
                }
            }
            
            log::debug!("Completed LUT chunk {}-{}/{}", chunk_start, chunk_end, height);
        }
        
        lut.is_precomputed = true;
        self.lut = Some(lut);
        
        let duration = start_time.elapsed();
        log::info!("Calibration LUT pre-computation completed in {:.3}s", duration.as_secs_f64());
        Ok(())
    }

    /// Get calibration value from pre-computed LUT (ULTRA FAST)
    pub fn get_calibration_value_from_lut(
        &self, 
        line: usize,
        pixel: usize,
        cal_type: CalibrationType
    ) -> SarResult<f32> {
        if let Some(ref lut) = self.lut {
            if !lut.is_precomputed {
                return Err(SarError::Processing("LUT not pre-computed".to_string()));
            }
            
            let values = match cal_type {
                CalibrationType::Sigma0 => &lut.sigma_values,
                CalibrationType::Beta0 => &lut.beta_values,
                CalibrationType::Gamma0 => &lut.gamma_values,
                CalibrationType::Dn => &lut.dn_values,
            };
            
            if line < values.nrows() && pixel < values.ncols() {
                Ok(values[[line, pixel]])
            } else {
                Err(SarError::Processing("Pixel coordinates out of bounds".to_string()))
            }
        } else {
            // Fallback to interpolation
            self.get_calibration_value_fast(line, pixel, cal_type)
        }
    }

    /// Get calibration values for a specific pixel using bilinear interpolation
    pub fn get_calibration_value(
        &self,
        line: usize,
        pixel: usize,
        cal_type: CalibrationType,
    ) -> SarResult<f32> {
        if self.vectors.is_empty() {
            return Err(SarError::Processing(
                "No calibration vectors available".to_string(),
            ));
        }

        // Find the two vectors surrounding this line
        let mut before_idx = 0;
        let mut after_idx = self.vectors.len() - 1;
        
        for (i, vector) in self.vectors.iter().enumerate() {
            if vector.line <= line {
                before_idx = i;
            }
            if vector.line >= line && after_idx == self.vectors.len() - 1 {
                after_idx = i;
                break;
            }
        }

        if before_idx == after_idx {
            // Exact line match or single vector
            return self.interpolate_pixel_value(&self.vectors[before_idx], pixel, cal_type);
        }

        // Bilinear interpolation between two vectors
        let before_vector = &self.vectors[before_idx];
        let after_vector = &self.vectors[after_idx];
        
        let before_value = self.interpolate_pixel_value(before_vector, pixel, cal_type)?;
        let after_value = self.interpolate_pixel_value(after_vector, pixel, cal_type)?;
        
        // Linear interpolation between lines
        if after_vector.line == before_vector.line {
            Ok(before_value)
        } else {
            let weight = (line - before_vector.line) as f32 / 
                        (after_vector.line - before_vector.line) as f32;
            Ok(before_value * (1.0 - weight) + after_value * weight)
        }
    }
    
    /// Interpolate calibration value for a specific pixel within a vector
    fn interpolate_pixel_value(
        &self,
        vector: &CalibrationVector,
        pixel: usize,
        cal_type: CalibrationType,
    ) -> SarResult<f32> {
        let values = match cal_type {
            CalibrationType::Sigma0 => &vector.sigma_nought,
            CalibrationType::Beta0 => &vector.beta_nought,
            CalibrationType::Gamma0 => &vector.gamma,
            CalibrationType::Dn => &vector.dn,
        };
        
        if values.is_empty() || vector.pixels.is_empty() {
            return Err(SarError::Processing(
                "Empty calibration vector".to_string(),
            ));
        }
        
        // Find pixel positions surrounding the target pixel
        let mut before_idx = 0;
        let mut after_idx = vector.pixels.len() - 1;
        
        for (i, &pix) in vector.pixels.iter().enumerate() {
            if pix <= pixel {
                before_idx = i;
            }
            if pix >= pixel && after_idx == vector.pixels.len() - 1 {
                after_idx = i;
                break;
            }
        }
        
        if before_idx == after_idx {
            // Exact pixel match or single value
            return Ok(values[before_idx]);
        }
        
        // Linear interpolation between pixels
        let before_pixel = vector.pixels[before_idx];
        let after_pixel = vector.pixels[after_idx];
        
        if after_pixel == before_pixel {
            Ok(values[before_idx])
        } else {
            let weight = (pixel - before_pixel) as f32 / 
                        (after_pixel - before_pixel) as f32;
            Ok(values[before_idx] * (1.0 - weight) + values[after_idx] * weight)
        }
    }

    /// Pre-compute calibration coefficients for faster access (OPTIMIZATION)
    pub fn precompute_coefficients(&mut self, image_dims: (usize, usize)) -> SarResult<()> {
        log::info!("Pre-computing calibration coefficients for {}x{} image", image_dims.0, image_dims.1);
        
        // Sort vectors by line for faster binary search
        self.vectors.sort_by_key(|v| v.line);
        
        // Pre-validate all vectors
        for (i, vector) in self.vectors.iter().enumerate() {
            if vector.pixels.is_empty() || vector.sigma_nought.is_empty() {
                return Err(SarError::Processing(
                    format!("Empty calibration vector at index {}", i)
                ));
            }
        }
        
        log::debug!("Calibration coefficients pre-computation completed");
        Ok(())
    }

    /// Get calibration value with optimized lookup (FAST VERSION)
    pub fn get_calibration_value_fast(
        &self,
        line: usize,
        pixel: usize,
        cal_type: CalibrationType,
    ) -> SarResult<f32> {
        if self.vectors.is_empty() {
            return Err(SarError::Processing(
                "No calibration vectors available".to_string(),
            ));
        }

        // Binary search for surrounding vectors (O(log n) instead of O(n))
        let (before_idx, after_idx) = self.find_surrounding_vectors_fast(line)?;

        if before_idx == after_idx {
            // Exact line match or single vector
            return self.interpolate_pixel_value_fast(&self.vectors[before_idx], pixel, cal_type);
        }

        // Fast bilinear interpolation between two vectors
        let before_vector = &self.vectors[before_idx];
        let after_vector = &self.vectors[after_idx];
        
        let before_value = self.interpolate_pixel_value_fast(before_vector, pixel, cal_type)?;
        let after_value = self.interpolate_pixel_value_fast(after_vector, pixel, cal_type)?;
        
        // Linear interpolation between lines
        if after_vector.line == before_vector.line {
            Ok(before_value)
        } else {
            let weight = (line - before_vector.line) as f32 / 
                        (after_vector.line - before_vector.line) as f32;
            Ok(before_value * (1.0 - weight) + after_value * weight)
        }
    }

    /// Binary search for surrounding vectors (much faster than linear search)
    fn find_surrounding_vectors_fast(&self, line: usize) -> SarResult<(usize, usize)> {
        // Binary search for the vector just before or at the line
        let mut left = 0;
        let mut right = self.vectors.len();
        
        while left < right {
            let mid = (left + right) / 2;
            if self.vectors[mid].line <= line {
                left = mid + 1;
            } else {
                right = mid;
            }
        }
        
        let before_idx = if left > 0 { left - 1 } else { 0 };
        let after_idx = if left < self.vectors.len() { left } else { self.vectors.len() - 1 };
        
        Ok((before_idx, after_idx))
    }

    /// Fast pixel interpolation with binary search for pixels
    fn interpolate_pixel_value_fast(
        &self,
        vector: &CalibrationVector,
        pixel: usize,
        cal_type: CalibrationType,
    ) -> SarResult<f32> {
        let values = match cal_type {
            CalibrationType::Sigma0 => &vector.sigma_nought,
            CalibrationType::Beta0 => &vector.beta_nought,
            CalibrationType::Gamma0 => &vector.gamma,
            CalibrationType::Dn => &vector.dn,
        };
        
        if values.is_empty() || vector.pixels.is_empty() {
            return Err(SarError::Processing("Empty calibration vector".to_string()));
        }
        
        // Binary search for surrounding pixels
        let mut left = 0;
        let mut right = vector.pixels.len();
        
        while left < right {
            let mid = (left + right) / 2;
            if vector.pixels[mid] <= pixel {
                left = mid + 1;
            } else {
                right = mid;
            }
        }
        
        let before_idx = if left > 0 { left - 1 } else { 0 };
        let after_idx = if left < vector.pixels.len() { left } else { vector.pixels.len() - 1 };
        
        if before_idx == after_idx {
            return Ok(values[before_idx]);
        }
        
        // Linear interpolation between pixels
        let before_pixel = vector.pixels[before_idx];
        let after_pixel = vector.pixels[after_idx];
        
        if after_pixel == before_pixel {
            Ok(values[before_idx])
        } else {
            let weight = (pixel - before_pixel) as f32 / 
                        (after_pixel - before_pixel) as f32;
            Ok(values[before_idx] * (1.0 - weight) + values[after_idx] * weight)
        }
    }
}

/// Radiometric calibration processor
pub struct CalibrationProcessor {
    coefficients: CalibrationCoefficients,
    calibration_type: CalibrationType,
}

/// Types of radiometric calibration
#[derive(Debug, Clone, Copy)]
pub enum CalibrationType {
    Sigma0,  // Radar cross section per unit area
    Beta0,   // Radar brightness  
    Gamma0,  // Backscatter coefficient
    Dn,      // Digital numbers (uncalibrated)
}

impl CalibrationProcessor {
    /// Create a new calibration processor
    pub fn new(coefficients: CalibrationCoefficients, calibration_type: CalibrationType) -> Self {
        Self {
            coefficients,
            calibration_type,
        }
    }

    /// Apply radiometric calibration to SLC data (ORIGINAL VERSION)
    pub fn calibrate(&self, slc_data: &SarImage) -> SarResult<SarRealImage> {
        log::info!("Applying radiometric calibration: {:?}", self.calibration_type);
        
        let (azimuth_lines, range_samples) = slc_data.dim();
        log::debug!("Input dimensions: {} x {}", azimuth_lines, range_samples);

        // Calculate intensity from complex SLC data (|SLC|^2)
        let mut intensity = Array2::zeros((azimuth_lines, range_samples));
        
        for ((i, j), &slc_pixel) in slc_data.indexed_iter() {
            intensity[[i, j]] = slc_pixel.norm_sqr(); // |SLC|^2
        }

        // Apply calibration lookup table
        let calibrated = match self.calibration_type {
            CalibrationType::Dn => intensity, // No calibration for DN
            _ => {
                // Apply proper calibration using the lookup table
                let mut calibrated = Array2::zeros((azimuth_lines, range_samples));
                for ((i, j), &intensity_val) in intensity.indexed_iter() {
                    let cal_value = self.coefficients.get_calibration_value(i, j, self.calibration_type)?;
                    calibrated[[i, j]] = intensity_val * cal_value;
                }
                calibrated
            }
        };

        log::info!("Calibration completed. Output range: {:.2e} to {:.2e}",
                  calibrated.iter().cloned().fold(f32::INFINITY, f32::min),
                  calibrated.iter().cloned().fold(f32::NEG_INFINITY, f32::max));

        Ok(calibrated)
    }

    /// Apply radiometric calibration with minimal, targeted optimization
    pub fn calibrate_optimized(&mut self, slc_data: &SarImage) -> SarResult<SarRealImage> {
        log::info!("Applying MINIMAL calibration optimization: {:?}", self.calibration_type);
        let start_time = std::time::Instant::now();
        
        let (azimuth_lines, range_samples) = slc_data.dim();
        log::debug!("Input dimensions: {} x {}", azimuth_lines, range_samples);

        // The ONLY optimization: avoid repeated vector sorting and use faster access patterns
        // This is much simpler and avoids all the expensive overhead
        
        // Sort vectors once for binary search (this is fast)
        self.coefficients.vectors.sort_by_key(|v| v.line);

        // Pre-allocate output array (this is necessary anyway)
        let mut calibrated = Array2::zeros((azimuth_lines, range_samples));

        match self.calibration_type {
            CalibrationType::Dn => {
                // DN: Just calculate intensity - no calibration needed
                for ((i, j), &slc_pixel) in slc_data.indexed_iter() {
                    calibrated[[i, j]] = slc_pixel.norm_sqr();
                }
            },
            _ => {
                // For other calibration types: use the existing fast interpolation but avoid repeated work
                // The key insight: batch process by line to reduce function call overhead
                for i in 0..azimuth_lines {
                    // Find surrounding vectors once per line (much more efficient)
                    let (before_idx, after_idx) = self.coefficients.find_surrounding_vectors_fast(i)?;
                    
                    for j in 0..range_samples {
                        let slc_pixel = slc_data[[i, j]];
                        let intensity = slc_pixel.norm_sqr();
                        
                        // Use existing fast interpolation with cached vector indices
                        let cal_value = if before_idx == after_idx {
                            // Same vector - direct interpolation
                            self.coefficients.interpolate_pixel_value_fast(
                                &self.coefficients.vectors[before_idx], j, self.calibration_type
                            )?
                        } else {
                            // Bilinear interpolation between vectors
                            let before_vector = &self.coefficients.vectors[before_idx];
                            let after_vector = &self.coefficients.vectors[after_idx];
                            
                            let before_value = self.coefficients.interpolate_pixel_value_fast(
                                before_vector, j, self.calibration_type
                            )?;
                            let after_value = self.coefficients.interpolate_pixel_value_fast(
                                after_vector, j, self.calibration_type
                            )?;
                            
                            // Linear interpolation between lines
                            if after_vector.line == before_vector.line {
                                before_value
                            } else {
                                let weight = (i - before_vector.line) as f32 / 
                                           (after_vector.line - before_vector.line) as f32;
                                before_value * (1.0 - weight) + after_value * weight
                            }
                        };
                        
                        calibrated[[i, j]] = intensity * cal_value;
                    }
                }
            }
        }

        let total_time = start_time.elapsed();
        log::info!("MINIMAL calibration optimization completed in {:.3}s. Output range: {:.2e} to {:.2e}",
                  total_time.as_secs_f64(),
                  calibrated.iter().cloned().fold(f32::INFINITY, f32::min),
                  calibrated.iter().cloned().fold(f32::NEG_INFINITY, f32::max));

        Ok(calibrated)
    }
    
    /// Calibration using pre-computed LUT for very large images
    fn calibrate_with_lut(&mut self, slc_data: &SarImage) -> SarResult<SarRealImage> {
        let (azimuth_lines, range_samples) = slc_data.dim();
        
        // Pre-compute calibration LUT if not already done
        if self.coefficients.lut.is_none() {
            log::debug!("Pre-computing LUT for large image");
            self.coefficients.precompute_lut((azimuth_lines, range_samples))?;
        }

        // Fast parallel intensity calculation using SIMD
        let mut intensity = Array2::zeros((azimuth_lines, range_samples));
        Zip::from(&mut intensity)
            .and(slc_data)
            .par_for_each(|intensity_val, &slc_pixel| {
                *intensity_val = slc_pixel.norm_sqr();
            });

        // Apply calibration using pre-computed LUT
        let calibrated = match self.calibration_type {
            CalibrationType::Dn => intensity, // No calibration for DN
            _ => {
                if let Some(ref lut) = self.coefficients.lut {
                    let cal_values = match self.calibration_type {
                        CalibrationType::Sigma0 => &lut.sigma_values,
                        CalibrationType::Beta0 => &lut.beta_values,
                        CalibrationType::Gamma0 => &lut.gamma_values,
                        CalibrationType::Dn => &lut.dn_values,
                    };
                    
                    // SIMD element-wise multiplication
                    let mut calibrated = Array2::zeros((azimuth_lines, range_samples));
                    Zip::from(&mut calibrated)
                        .and(&intensity)
                        .and(cal_values)
                        .par_for_each(|out, &intensity_val, &cal_val| {
                            *out = intensity_val * cal_val;
                        });
                    
                    calibrated
                } else {
                    return Err(SarError::Processing("Calibration LUT not available".to_string()));
                }
            }
        };

        Ok(calibrated)
    }

    /// Apply calibration using chunked processing for memory efficiency
    pub fn calibrate_chunked(&mut self, slc_data: &SarImage, chunk_size: usize) -> SarResult<SarRealImage> {
        log::info!("Applying CHUNKED calibration: {:?} (chunk_size: {})", self.calibration_type, chunk_size);
        let start_time = std::time::Instant::now();
        
        let (azimuth_lines, range_samples) = slc_data.dim();
        let mut calibrated = Array2::zeros((azimuth_lines, range_samples));
        
        // Pre-compute calibration LUT if needed
        if self.coefficients.lut.is_none() {
            self.coefficients.precompute_lut((azimuth_lines, range_samples))?;
        }
        
        // Process in chunks to reduce memory usage
        let num_chunks = (azimuth_lines + chunk_size - 1) / chunk_size;
        
        for chunk_idx in 0..num_chunks {
            let start_row = chunk_idx * chunk_size;
            let end_row = std::cmp::min(start_row + chunk_size, azimuth_lines);
            
            log::debug!("Processing chunk {}/{}: rows {}-{}", chunk_idx + 1, num_chunks, start_row, end_row);
            
            // Process this chunk
            for i in start_row..end_row {
                for j in 0..range_samples {
                    let slc_pixel = slc_data[[i, j]];
                    let intensity = slc_pixel.norm_sqr();
                    
                    let calibrated_val = match self.calibration_type {
                        CalibrationType::Dn => intensity,
                        _ => {
                            let cal_value = self.coefficients.get_calibration_value_from_lut(i, j, self.calibration_type)?;
                            intensity * cal_value
                        }
                    };
                    
                    calibrated[[i, j]] = calibrated_val;
                }
            }
        }
        
        let total_time = start_time.elapsed();
        log::info!("CHUNKED calibration completed in {:.3}s", total_time.as_secs_f64());
        
        Ok(calibrated)
    }
}
