use crate::types::{SarError, SarImage, SarRealImage, SarResult, Polarization};
use ndarray::{Array1, Array2, Zip};
use regex::Regex;

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

/// Calibration coefficients for radiometric correction
#[derive(Debug, Clone)]
pub struct CalibrationCoefficients {
    pub vectors: Vec<CalibrationVector>,
    pub swath: String,
    pub polarization: String,
    pub product_first_line_utc_time: String,
    pub product_last_line_utc_time: String,
}

impl CalibrationCoefficients {
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

    /// Apply radiometric calibration to SLC data
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
            CalibrationType::Sigma0 => self.apply_sigma0_calibration_optimized(&intensity)?,
            CalibrationType::Beta0 => self.apply_beta0_calibration_optimized(&intensity)?,
            CalibrationType::Gamma0 => self.apply_gamma0_calibration_optimized(&intensity)?,
            CalibrationType::Dn => intensity, // No calibration for DN
        };

        log::info!("Calibration completed. Output range: {:.2e} to {:.2e}",
                  calibrated.iter().cloned().fold(f32::INFINITY, f32::min),
                  calibrated.iter().cloned().fold(f32::NEG_INFINITY, f32::max));

        Ok(calibrated)
    }

    /// Apply radiometric calibration to SLC data (OPTIMIZED VERSION)
    pub fn calibrate_optimized(&self, slc_data: &SarImage) -> SarResult<SarRealImage> {
        log::info!("Applying optimized radiometric calibration: {:?}", self.calibration_type);
        
        let (azimuth_lines, range_samples) = slc_data.dim();
        log::debug!("Input dimensions: {} x {}", azimuth_lines, range_samples);

        // Pre-calculate intensity from complex SLC data (|SLC|^2)
        let intensity = slc_data.mapv(|slc_pixel| slc_pixel.norm_sqr());

        // Apply optimized calibration based on type
        let calibrated = match self.calibration_type {
            CalibrationType::Sigma0 => self.apply_sigma0_calibration_optimized(&intensity)?,
            CalibrationType::Beta0 => self.apply_beta0_calibration_optimized(&intensity)?,
            CalibrationType::Gamma0 => self.apply_gamma0_calibration_optimized(&intensity)?,
            CalibrationType::Dn => intensity, // No calibration for DN
        };

        log::info!("Optimized calibration completed. Output range: {:.2e} to {:.2e}",
                  calibrated.iter().cloned().fold(f32::INFINITY, f32::min),
                  calibrated.iter().cloned().fold(f32::NEG_INFINITY, f32::max));

        Ok(calibrated)
    }

    /// Apply Sigma-0 calibration with pre-computed lookup table (FAST)
    fn apply_sigma0_calibration_optimized(&self, intensity: &SarRealImage) -> SarResult<SarRealImage> {
        log::debug!("Applying optimized Sigma-0 calibration");
        
        let (azimuth_lines, range_samples) = intensity.dim();
        
        // Pre-compute calibration lookup table for entire image
        log::debug!("Pre-computing calibration lookup table...");
        let cal_lut = self.build_calibration_lut(azimuth_lines, range_samples, CalibrationType::Sigma0)?;
        
        // Vectorized calibration: intensity / cal_coeff^2
        log::debug!("Applying vectorized calibration...");
        let calibrated = Zip::from(intensity).and(&cal_lut).map_collect(|intensity_val, &cal_coeff| {
            if cal_coeff > 0.0 {
                intensity_val / (cal_coeff * cal_coeff)
            } else {
                0.0
            }
        });

        Ok(calibrated)
    }

    /// Apply Beta-0 calibration (optimized)
    fn apply_beta0_calibration_optimized(&self, intensity: &SarRealImage) -> SarResult<SarRealImage> {
        log::debug!("Applying optimized Beta-0 calibration");
        
        let (azimuth_lines, range_samples) = intensity.dim();
        let cal_lut = self.build_calibration_lut(azimuth_lines, range_samples, CalibrationType::Beta0)?;
        
        let calibrated = Zip::from(intensity).and(&cal_lut).map_collect(|intensity_val, &cal_coeff| {
            if cal_coeff > 0.0 {
                intensity_val / (cal_coeff * cal_coeff)
            } else {
                0.0
            }
        });

        Ok(calibrated)
    }

    /// Apply Gamma-0 calibration (optimized)
    fn apply_gamma0_calibration_optimized(&self, intensity: &SarRealImage) -> SarResult<SarRealImage> {
        log::debug!("Applying optimized Gamma-0 calibration");
        
        let (azimuth_lines, range_samples) = intensity.dim();
        let cal_lut = self.build_calibration_lut(azimuth_lines, range_samples, CalibrationType::Gamma0)?;
        
        let calibrated = Zip::from(intensity).and(&cal_lut).map_collect(|intensity_val, &cal_coeff| {
            if cal_coeff > 0.0 {
                intensity_val / (cal_coeff * cal_coeff)
            } else {
                0.0
            }
        });

        Ok(calibrated)
    }

    /// Build calibration lookup table for entire image (KEY OPTIMIZATION)
    fn build_calibration_lut(
        &self,
        azimuth_lines: usize,
        range_samples: usize,
        cal_type: CalibrationType,
    ) -> SarResult<Array2<f32>> {
        log::debug!("Building {}x{} calibration LUT", azimuth_lines, range_samples);
        
        // Strategy 1: Sparse interpolation (for large images)
        if azimuth_lines * range_samples > 50_000_000 {  // > 50M pixels
            return self.build_sparse_calibration_lut(azimuth_lines, range_samples, cal_type);
        }
        
        // Strategy 2: Full LUT for smaller images
        let mut cal_lut = Array2::zeros((azimuth_lines, range_samples));
        
        // Process in chunks for better cache performance
        let chunk_size = 1000;
        
        for azimuth_chunk in (0..azimuth_lines).step_by(chunk_size) {
            let azimuth_end = (azimuth_chunk + chunk_size).min(azimuth_lines);
            
            for range_chunk in (0..range_samples).step_by(chunk_size) {
                let range_end = (range_chunk + chunk_size).min(range_samples);
                
                // Process chunk
                for i in azimuth_chunk..azimuth_end {
                    for j in range_chunk..range_end {
                        cal_lut[[i, j]] = self.coefficients.get_calibration_value(i, j, cal_type)?;
                    }
                }
            }
        }
        
        Ok(cal_lut)
    }

    /// Build sparse calibration lookup table for very large images
    fn build_sparse_calibration_lut(
        &self,
        azimuth_lines: usize,
        range_samples: usize,
        cal_type: CalibrationType,
    ) -> SarResult<Array2<f32>> {
        log::debug!("Building sparse calibration LUT (large image optimization)");
        
        // Use lower resolution grid and interpolate
        let sparse_factor = 10; // Sample every 10th pixel
        let sparse_az = (azimuth_lines + sparse_factor - 1) / sparse_factor;
        let sparse_rg = (range_samples + sparse_factor - 1) / sparse_factor;
        
        // Build sparse LUT
        let mut sparse_lut = Array2::zeros((sparse_az, sparse_rg));
        for i in 0..sparse_az {
            for j in 0..sparse_rg {
                let orig_i = (i * sparse_factor).min(azimuth_lines - 1);
                let orig_j = (j * sparse_factor).min(range_samples - 1);
                sparse_lut[[i, j]] = self.coefficients.get_calibration_value(orig_i, orig_j, cal_type)?;
            }
        }
        
        // Interpolate to full resolution
        self.interpolate_sparse_lut(&sparse_lut, azimuth_lines, range_samples, sparse_factor)
    }

    /// Interpolate sparse LUT to full resolution
    fn interpolate_sparse_lut(
        &self,
        sparse_lut: &Array2<f32>,
        target_azimuth: usize,
        target_range: usize,
        sparse_factor: usize,
    ) -> SarResult<Array2<f32>> {
        log::debug!("Interpolating sparse LUT to full resolution");
        
        let mut full_lut = Array2::zeros((target_azimuth, target_range));
        
        for i in 0..target_azimuth {
            for j in 0..target_range {
                // Map to sparse coordinates
                let sparse_i = (i as f32) / (sparse_factor as f32);
                let sparse_j = (j as f32) / (sparse_factor as f32);
                
                // Bilinear interpolation
                let i0 = sparse_i.floor() as usize;
                let i1 = (i0 + 1).min(sparse_lut.nrows() - 1);
                let j0 = sparse_j.floor() as usize;
                let j1 = (j0 + 1).min(sparse_lut.ncols() - 1);
                
                let fi = sparse_i - i0 as f32;
                let fj = sparse_j - j0 as f32;
                
                // Bilinear interpolation formula
                let val = sparse_lut[[i0, j0]] * (1.0 - fi) * (1.0 - fj) +
                         sparse_lut[[i1, j0]] * fi * (1.0 - fj) +
                         sparse_lut[[i0, j1]] * (1.0 - fi) * fj +
                         sparse_lut[[i1, j1]] * fi * fj;
                
                full_lut[[i, j]] = val;
            }
        }
        
        Ok(full_lut)
    }

    /// Parallel calibration using Rayon (if available)
    #[cfg(feature = "parallel")]
    pub fn calibrate_parallel(&self, slc_data: &SarImage) -> SarResult<SarRealImage> {
        use rayon::prelude::*;
        
        log::info!("Applying parallel radiometric calibration: {:?}", self.calibration_type);
        
        let (azimuth_lines, range_samples) = slc_data.dim();
        let intensity = slc_data.mapv(|slc_pixel| slc_pixel.norm_sqr());
        
        // Build calibration LUT in parallel
        let cal_lut = self.build_calibration_lut_parallel(azimuth_lines, range_samples, self.calibration_type)?;
        
        // Apply calibration in parallel
        let calibrated = Zip::from(&intensity).and(&cal_lut).par_map_collect(|intensity_val, cal_coeff| {
            if *cal_coeff > 0.0 {
                intensity_val / (cal_coeff * cal_coeff)
            } else {
                0.0
            }
        });
        
        Ok(calibrated)
    }

    #[cfg(feature = "parallel")]
    fn build_calibration_lut_parallel(
        &self,
        azimuth_lines: usize,
        range_samples: usize,
        cal_type: CalibrationType,
    ) -> SarResult<Array2<f32>> {
        use rayon::prelude::*;
        
        log::debug!("Building calibration LUT in parallel");
        
        let coords: Vec<(usize, usize)> = (0..azimuth_lines)
            .flat_map(|i| (0..range_samples).map(move |j| (i, j)))
            .collect();
        
        let values: Result<Vec<f32>, SarError> = coords
            .into_par_iter()
            .map(|(i, j)| self.coefficients.get_calibration_value(i, j, cal_type))
            .collect();
        
        let cal_values = values?;
        
        Array2::from_shape_vec((azimuth_lines, range_samples), cal_values)
            .map_err(|e| SarError::Processing(format!("Shape error: {}", e)))
    }

    /// Convert calibrated data to dB scale
    pub fn to_db(linear_data: &SarRealImage) -> SarRealImage {
        log::debug!("Converting to dB scale");
        
        linear_data.mapv(|x| {
            if x > 0.0 {
                10.0 * x.log10()
            } else {
                -50.0 // Minimum dB value for zero/negative values
            }
        })
    }

    /// Apply thermal noise correction
    pub fn apply_thermal_noise_correction(
        &self,
        intensity: &SarRealImage,
        noise_lut: &Array2<f32>,
    ) -> SarResult<SarRealImage> {
        log::debug!("Applying thermal noise correction");
        
        let (azimuth_lines, range_samples) = intensity.dim();
        let mut corrected = intensity.clone();

        for i in 0..azimuth_lines {
            for j in 0..range_samples {
                let noise_power = self.interpolate_calibration_coefficient(
                    noise_lut, i, j
                )?;
                
                // Subtract noise power
                corrected[[i, j]] = (corrected[[i, j]] - noise_power).max(0.0);
            }
        }

        Ok(corrected)
    }

    /// Interpolate calibration coefficient for a given pixel
    fn interpolate_calibration_coefficient(
        &self,
        lut: &Array2<f32>,
        azimuth_idx: usize,
        range_idx: usize,
    ) -> SarResult<f32> {
        let (lut_az, lut_rg) = lut.dim();
        
        if lut_az == 0 || lut_rg == 0 {
            return Err(SarError::Processing(
                "Empty calibration lookup table".to_string(),
            ));
        }

        // Simple nearest neighbor interpolation for now
        // In practice, you'd use bilinear interpolation
        let az_idx = (azimuth_idx * lut_az / azimuth_idx.max(1)).min(lut_az - 1);
        let rg_idx = (range_idx * lut_rg / range_idx.max(1)).min(lut_rg - 1);
        
        Ok(lut[[az_idx, rg_idx]])
    }

    /// Benchmark calibration performance
    pub fn benchmark_calibration(&self, slc_data: &SarImage) -> SarResult<()> {
        use std::time::Instant;
        
        log::info!("=== Calibration Performance Benchmark ===");
        let (azimuth_lines, range_samples) = slc_data.dim();
        let total_pixels = azimuth_lines * range_samples;
        log::info!("Image size: {}x{} = {} pixels", azimuth_lines, range_samples, total_pixels);
        
        // Test 1: Original method
        log::info!("Testing original calibration method...");
        let start = Instant::now();
        let _result1 = self.calibrate(slc_data)?;
        let original_time = start.elapsed();
        log::info!("Original method: {:.2} seconds ({:.0} pixels/sec)", 
                  original_time.as_secs_f64(), 
                  total_pixels as f64 / original_time.as_secs_f64());
        
        // Test 2: Optimized method
        log::info!("Testing optimized calibration method...");
        let start = Instant::now();
        let _result2 = self.calibrate_optimized(slc_data)?;
        let optimized_time = start.elapsed();
        log::info!("Optimized method: {:.2} seconds ({:.0} pixels/sec)", 
                  optimized_time.as_secs_f64(),
                  total_pixels as f64 / optimized_time.as_secs_f64());
        
        // Test 3: Parallel method (if available)
        #[cfg(feature = "parallel")]
        {
            log::info!("Testing parallel calibration method...");
            let start = Instant::now();
            let _result3 = self.calibrate_parallel(slc_data)?;
            let parallel_time = start.elapsed();
            log::info!("Parallel method: {:.2} seconds ({:.0} pixels/sec)", 
                      parallel_time.as_secs_f64(),
                      total_pixels as f64 / parallel_time.as_secs_f64());
            
            let parallel_speedup = original_time.as_secs_f64() / parallel_time.as_secs_f64();
            log::info!("Parallel speedup: {:.2}x", parallel_speedup);
        }
        
        let speedup = original_time.as_secs_f64() / optimized_time.as_secs_f64();
        log::info!("Optimized speedup: {:.2}x", speedup);
        log::info!("=== Benchmark Complete ===");
        
        Ok(())
    }

    /// Get recommended calibration method based on image size
    pub fn get_recommended_method(&self, image_dims: (usize, usize)) -> &'static str {
        let total_pixels = image_dims.0 * image_dims.1;
        
        match total_pixels {
            0..=1_000_000 => "standard",        // < 1M pixels: standard is fine
            1_000_001..=10_000_000 => "optimized",  // 1-10M pixels: use optimized
            10_000_001..=100_000_000 => "sparse",   // 10-100M pixels: use sparse LUT
            _ => "parallel"                          // > 100M pixels: use parallel
        }
    }
}

/// Parse calibration data from Sentinel-1 calibration XML
pub fn parse_calibration_from_xml(xml_content: &str) -> SarResult<CalibrationCoefficients> {
    log::debug!("Parsing calibration data from XML (length: {})", xml_content.len());
    
    // Debug: show first 1000 chars of XML
    if xml_content.len() > 0 {
        let preview = if xml_content.len() > 1000 { &xml_content[..1000] } else { xml_content };
        log::debug!("XML preview: {}", preview);
    }
    
    // Extract calibration vectors using regex (fallback to XML parsing if needed)
    let vectors = parse_calibration_vectors_regex(xml_content)?;
    
    if vectors.is_empty() {
        log::warn!("No calibration vectors found with regex, trying fallback parsing");
        return Err(SarError::Processing(
            "No calibration vectors found in XML".to_string(),
        ));
    }
    
    // Extract metadata
    let swath = extract_xml_value(xml_content, "swath").unwrap_or_else(|| "IW".to_string());
    let polarization = extract_xml_value(xml_content, "polarisation").unwrap_or_else(|| "VV".to_string());
    let first_time = extract_xml_value(xml_content, "productFirstLineUtcTime").unwrap_or_else(|| "".to_string());
    let last_time = extract_xml_value(xml_content, "productLastLineUtcTime").unwrap_or_else(|| "".to_string());
    
    log::info!("Parsed {} calibration vectors for {}/{}", vectors.len(), swath, polarization);
    
    Ok(CalibrationCoefficients {
        vectors,
        swath,
        polarization,
        product_first_line_utc_time: first_time,
        product_last_line_utc_time: last_time,
    })
}

/// Parse calibration vectors from XML using regex patterns
fn parse_calibration_vectors_regex(xml_content: &str) -> SarResult<Vec<CalibrationVector>> {
    let mut vectors = Vec::new();
    
    log::debug!("XML content length: {} chars", xml_content.len());
    
    // Count calibrationVector tags first
    let tag_count = xml_content.matches("<calibrationVector>").count();
    log::debug!("Found {} <calibrationVector> tags", tag_count);
    
    // First, try using simple string search for calibrationVector tags
    let mut start_pos = 0;
    
    while let Some(start_tag_pos) = xml_content[start_pos..].find("<calibrationVector>") {
        let absolute_start_pos = start_pos + start_tag_pos;
        let search_start = absolute_start_pos + 19; // Length of "<calibrationVector>"
        
        if let Some(end_tag_pos) = xml_content[search_start..].find("</calibrationVector>") {
            let absolute_end_pos = search_start + end_tag_pos;
            // Extract just the content between the tags (without the tags themselves)
            let vector_xml = &xml_content[search_start..absolute_end_pos];
            
            match parse_single_calibration_vector(vector_xml) {
                Ok(vector) => {
                    vectors.push(vector);
                    log::debug!("Successfully parsed calibration vector {}", vectors.len());
                }
                Err(e) => {
                    log::warn!("Failed to parse calibration vector: {}", e);
                }
            }
            
            start_pos = absolute_end_pos + 20; // Move past "</calibrationVector>"
        } else {
            break;
        }
    }
    
    log::debug!("Extracted {} calibration vectors using string search", vectors.len());
    
    // If string search failed, try regex as fallback
    if vectors.is_empty() {
        log::debug!("Trying regex fallback");
        
        let vector_pattern = Regex::new(
            r"(?s)<calibrationVector>(.*?)</calibrationVector>"
        ).map_err(|e| SarError::Processing(format!("Regex error: {}", e)))?;
        
        for captures in vector_pattern.captures_iter(xml_content) {
            if let Some(vector_match) = captures.get(1) {
                let vector_xml = vector_match.as_str();
                log::debug!("Processing vector XML of length: {}", vector_xml.len());
                
                match parse_single_calibration_vector(vector_xml) {
                    Ok(vector) => {
                        vectors.push(vector);
                        log::debug!("Successfully parsed calibration vector {}", vectors.len());
                    }
                    Err(e) => {
                        log::warn!("Failed to parse calibration vector: {}", e);
                    }
                }
            }
        }
        
        log::debug!("Extracted {} calibration vectors using regex fallback", vectors.len());
    }
    
    Ok(vectors)
}

/// Parse a single calibration vector from XML
fn parse_single_calibration_vector(vector_xml: &str) -> SarResult<CalibrationVector> {
    // Extract azimuth time
    let azimuth_time = extract_xml_value(vector_xml, "azimuthTime")
        .ok_or_else(|| SarError::Processing("Missing azimuthTime".to_string()))?;
    
    // Extract line number
    let line_str = extract_xml_value(vector_xml, "line")
        .ok_or_else(|| SarError::Processing("Missing line".to_string()))?;
    let line = line_str.parse::<usize>()
        .map_err(|e| SarError::Processing(format!("Invalid line number: {}", e)))?;
    
    // Extract pixel coordinates
    let pixel_str = extract_xml_value(vector_xml, "pixel")
        .ok_or_else(|| SarError::Processing("Missing pixel".to_string()))?;
    let pixels = parse_space_separated_numbers::<usize>(&pixel_str)?;
    
    // Extract calibration values
    let sigma_str = extract_xml_value(vector_xml, "sigmaNought")
        .ok_or_else(|| SarError::Processing("Missing sigmaNought".to_string()))?;
    let sigma_nought = parse_space_separated_numbers::<f32>(&sigma_str)?;
    
    let beta_str = extract_xml_value(vector_xml, "betaNought")
        .ok_or_else(|| SarError::Processing("Missing betaNought".to_string()))?;
    let beta_nought = parse_space_separated_numbers::<f32>(&beta_str)?;
    
    let gamma_str = extract_xml_value(vector_xml, "gamma")
        .ok_or_else(|| SarError::Processing("Missing gamma".to_string()))?;
    let gamma = parse_space_separated_numbers::<f32>(&gamma_str)?;
    
    let dn_str = extract_xml_value(vector_xml, "dn")
        .ok_or_else(|| SarError::Processing("Missing dn".to_string()))?;
    let dn = parse_space_separated_numbers::<f32>(&dn_str)?;
    
    // Validate array lengths
    if pixels.len() != sigma_nought.len() || 
       pixels.len() != beta_nought.len() || 
       pixels.len() != gamma.len() || 
       pixels.len() != dn.len() {
        return Err(SarError::Processing(
            "Calibration vector arrays have mismatched lengths".to_string(),
        ));
    }
    
    Ok(CalibrationVector {
        azimuth_time,
        line,
        pixels,
        sigma_nought,
        beta_nought,
        gamma,
        dn,
    })
}

/// Parse space-separated numbers from a string
fn parse_space_separated_numbers<T>(input: &str) -> SarResult<Vec<T>>
where
    T: std::str::FromStr,
    T::Err: std::fmt::Display,
{
    input
        .split_whitespace()
        .map(|s| s.parse::<T>().map_err(|e| 
            SarError::Processing(format!("Parse error: {}", e))
        ))
        .collect()
}

/// Extract XML value using regex (fallback for robust parsing)
fn extract_xml_value(xml_content: &str, tag: &str) -> Option<String> {
    // Handle tags with or without attributes
    let pattern = format!(r"<{}\s*[^>]*>\s*([^<]*)\s*</{}>", tag, tag);
    if let Ok(re) = Regex::new(&pattern) {
        if let Some(cap) = re.captures(xml_content) {
            return Some(cap[1].trim().to_string());
        }
    }
    
    // Fallback: simple tag without attributes
    let simple_pattern = format!(r"<{}>\s*([^<]*)\s*</{}>", tag, tag);
    if let Ok(re) = Regex::new(&simple_pattern) {
        if let Some(cap) = re.captures(xml_content) {
            return Some(cap[1].trim().to_string());
        }
    }
    
    None
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::types::SarComplex;
    use ndarray::Array2;

    #[test]
    fn test_calibration_processor() {
        // Create test calibration vectors
        let mut vectors = Vec::new();
        for i in 0..10 {
            vectors.push(CalibrationVector {
                azimuth_time: format!("2020-01-01T00:00:{:02}", i),
                line: i * 100,
                pixels: (0..100).step_by(10).collect(),
                sigma_nought: vec![1000.0; 10],
                beta_nought: vec![800.0; 10],
                gamma: vec![1200.0; 10],
                dn: vec![500.0; 10],
            });
        }
        
        let coefficients = CalibrationCoefficients {
            vectors,
            swath: "IW1".to_string(),
            polarization: "VV".to_string(),
            product_first_line_utc_time: "2020-01-01T00:00:00".to_string(),
            product_last_line_utc_time: "2020-01-01T00:01:00".to_string(),
        };

        let processor = CalibrationProcessor::new(coefficients, CalibrationType::Sigma0);
        let test_slc = Array2::from_elem((100, 100), SarComplex::new(100.0, 0.0));
        
        let result = processor.calibrate(&test_slc);
        assert!(result.is_ok());
        
        let calibrated = result.unwrap();
        assert_eq!(calibrated.dim(), (100, 100));
    }

    #[test]
    fn test_db_conversion() {
        let linear_data = Array2::from_elem((10, 10), 100.0);
        let db_data = CalibrationProcessor::to_db(&linear_data);
        
        // 100.0 in linear scale should be 20 dB
        assert!((db_data[[0, 0]] - 20.0).abs() < 1e-6);
    }

    #[test]
    fn test_calibration_vector_parsing() {
        let test_xml = r#"
        <calibrationVector>
            <azimuthTime>2020-01-03T17:08:15.674828</azimuthTime>
            <line>0</line>
            <pixel>0 40 80 120 160</pixel>
            <sigmaNought>3.339847e+02 3.339192e+02 3.338538e+02 3.337885e+02 3.337232e+02</sigmaNought>
            <betaNought>2.370000e+02 2.370000e+02 2.370000e+02 2.370000e+02 2.370000e+02</betaNought>
            <gamma>3.104379e+02 3.103564e+02 3.102749e+02 3.101935e+02 3.101122e+02</gamma>
            <dn>2.370000e+02 2.370000e+02 2.370000e+02 2.370000e+02 2.370000e+02</dn>
        </calibrationVector>
        "#;
        
        let result = parse_single_calibration_vector(test_xml);
        assert!(result.is_ok());
        
        let vector = result.unwrap();
        assert_eq!(vector.line, 0);
        assert_eq!(vector.pixels.len(), 5);
        assert_eq!(vector.sigma_nought.len(), 5);
        assert!((vector.sigma_nought[0] - 333.9847).abs() < 1e-3);
    }
}
