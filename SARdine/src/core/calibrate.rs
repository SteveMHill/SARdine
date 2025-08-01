use crate::types::{SarError, SarImage, SarRealImage, SarResult};
use ndarray::{Array2, Zip, Axis, s};
use rayon::prelude::*;
use quick_xml::Reader;
use quick_xml::events::Event;

/// Parse calibration data from Sentinel-1 XML file
pub fn parse_calibration_from_xml(xml_content: &str) -> SarResult<CalibrationCoefficients> {
    log::debug!("Parsing calibration XML content ({} bytes)", xml_content.len());
    
    let mut reader = Reader::from_str(xml_content);
    reader.trim_text(true);
    
    let mut calibration = CalibrationCoefficients::new();
    let mut buf = Vec::new();
    
    // State tracking for parsing
    let mut in_ads_header = false;
    let mut in_calibration_vector = false;
    let mut in_calibration_vector_list = false;
    let mut current_vector: Option<CalibrationVector> = None;
    let mut current_tag = String::new();
    let mut text_content = String::new();
    
    loop {
        match reader.read_event_into(&mut buf) {
            Ok(Event::Start(ref e)) => {
                current_tag = String::from_utf8_lossy(e.name().as_ref()).to_string();
                
                match current_tag.as_str() {
                    "adsHeader" => in_ads_header = true,
                    "calibrationVectorList" => in_calibration_vector_list = true,
                    "calibrationVector" => {
                        in_calibration_vector = true;
                        current_vector = Some(CalibrationVector {
                            azimuth_time: String::new(),
                            line: 0,
                            pixels: Vec::new(),
                            sigma_nought: Vec::new(),
                            beta_nought: Vec::new(),
                            gamma: Vec::new(),
                            dn: Vec::new(),
                        });
                    }
                    _ => {}
                }
            }
            Ok(Event::End(ref e)) => {
                let tag_name = String::from_utf8_lossy(e.name().as_ref()).to_string();
                
                match tag_name.as_str() {
                    "adsHeader" => in_ads_header = false,
                    "calibrationVectorList" => in_calibration_vector_list = false,
                    "calibrationVector" => {
                        if let Some(vector) = current_vector.take() {
                            calibration.vectors.push(vector);
                        }
                        in_calibration_vector = false;
                    }
                    _ => {}
                }
                current_tag.clear();
            }
            Ok(Event::Text(ref e)) => {
                text_content = e.unescape().unwrap().to_string();
                
                if in_ads_header {
                    match current_tag.as_str() {
                        "polarisation" => calibration.polarization = text_content.clone(),
                        "swath" => calibration.swath = text_content.clone(),
                        "startTime" => calibration.product_first_line_utc_time = text_content.clone(),
                        "stopTime" => calibration.product_last_line_utc_time = text_content.clone(),
                        _ => {}
                    }
                } else if in_calibration_vector {
                    if let Some(ref mut vector) = current_vector {
                        match current_tag.as_str() {
                            "azimuthTime" => vector.azimuth_time = text_content.clone(),
                            "line" => {
                                vector.line = text_content.parse().unwrap_or(0);
                            }
                            "pixel" => {
                                // Parse space-separated pixel indices
                                vector.pixels = text_content
                                    .split_whitespace()
                                    .filter_map(|s| s.parse().ok())
                                    .collect();
                            }
                            "sigmaNought" => {
                                // Parse space-separated sigma nought values
                                vector.sigma_nought = text_content
                                    .split_whitespace()
                                    .filter_map(|s| s.parse().ok())
                                    .collect();
                            }
                            "betaNought" => {
                                // Parse space-separated beta nought values
                                vector.beta_nought = text_content
                                    .split_whitespace()
                                    .filter_map(|s| s.parse().ok())
                                    .collect();
                            }
                            "gamma" => {
                                // Parse space-separated gamma values
                                vector.gamma = text_content
                                    .split_whitespace()
                                    .filter_map(|s| s.parse().ok())
                                    .collect();
                            }
                            "dn" => {
                                // Parse space-separated DN values
                                vector.dn = text_content
                                    .split_whitespace()
                                    .filter_map(|s| s.parse().ok())
                                    .collect();
                            }
                            _ => {}
                        }
                    }
                }
            }
            Ok(Event::Eof) => break,
            Err(e) => {
                return Err(SarError::Processing(format!("XML parsing error: {}", e)));
            }
            _ => {}
        }
        buf.clear();
    }
    
    log::info!("Parsed {} calibration vectors for swath {} polarization {}", 
               calibration.vectors.len(), calibration.swath, calibration.polarization);
    
    if calibration.vectors.is_empty() {
        return Err(SarError::Processing("No calibration vectors found in XML".to_string()));
    }
    
    Ok(calibration)
}

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

    /// Get reference to calibration data
    pub fn get_calibration_data(&self) -> &CalibrationCoefficients {
        &self.coefficients
    }
}
