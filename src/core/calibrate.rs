use crate::types::{SarComplex, SarError, SarImage, SarReal, SarRealImage, SarResult};
use ndarray::{Array1, Array2};
use std::f32::consts::PI;

/// Calibration coefficients for radiometric correction
#[derive(Debug, Clone)]
pub struct CalibrationCoefficients {
    pub sigma0_lut: Array2<f32>,  // Sigma-0 lookup table
    pub beta0_lut: Array2<f32>,   // Beta-0 lookup table  
    pub gamma0_lut: Array2<f32>,  // Gamma-0 lookup table
    pub dn_lut: Array2<f32>,      // Digital number lookup table
    pub range_times: Array1<f64>, // Slant range times
    pub azimuth_times: Array1<f64>, // Azimuth times
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
            CalibrationType::Sigma0 => self.apply_sigma0_calibration(&intensity)?,
            CalibrationType::Beta0 => self.apply_beta0_calibration(&intensity)?,
            CalibrationType::Gamma0 => self.apply_gamma0_calibration(&intensity)?,
            CalibrationType::Dn => intensity, // No calibration for DN
        };

        log::info!("Calibration completed. Output range: {:.2e} to {:.2e}",
                  calibrated.iter().cloned().fold(f32::INFINITY, f32::min),
                  calibrated.iter().cloned().fold(f32::NEG_INFINITY, f32::max));

        Ok(calibrated)
    }

    /// Apply Sigma-0 calibration
    fn apply_sigma0_calibration(&self, intensity: &SarRealImage) -> SarResult<SarRealImage> {
        log::debug!("Applying Sigma-0 calibration");
        
        let (azimuth_lines, range_samples) = intensity.dim();
        let mut calibrated = Array2::zeros((azimuth_lines, range_samples));

        for i in 0..azimuth_lines {
            for j in 0..range_samples {
                // Interpolate calibration coefficient
                let cal_coeff = self.interpolate_calibration_coefficient(
                    &self.coefficients.sigma0_lut, i, j
                )?;
                
                // Apply calibration: sigma0 = intensity / cal_coeff^2
                if cal_coeff > 0.0 {
                    calibrated[[i, j]] = intensity[[i, j]] / (cal_coeff * cal_coeff);
                } else {
                    calibrated[[i, j]] = 0.0;
                }
            }
        }

        Ok(calibrated)
    }

    /// Apply Beta-0 calibration
    fn apply_beta0_calibration(&self, intensity: &SarRealImage) -> SarResult<SarRealImage> {
        log::debug!("Applying Beta-0 calibration");
        
        let (azimuth_lines, range_samples) = intensity.dim();
        let mut calibrated = Array2::zeros((azimuth_lines, range_samples));

        for i in 0..azimuth_lines {
            for j in 0..range_samples {
                let cal_coeff = self.interpolate_calibration_coefficient(
                    &self.coefficients.beta0_lut, i, j
                )?;
                
                if cal_coeff > 0.0 {
                    calibrated[[i, j]] = intensity[[i, j]] / (cal_coeff * cal_coeff);
                } else {
                    calibrated[[i, j]] = 0.0;
                }
            }
        }

        Ok(calibrated)
    }

    /// Apply Gamma-0 calibration (includes terrain correction factor)
    fn apply_gamma0_calibration(&self, intensity: &SarRealImage) -> SarResult<SarRealImage> {
        log::debug!("Applying Gamma-0 calibration");
        
        let (azimuth_lines, range_samples) = intensity.dim();
        let mut calibrated = Array2::zeros((azimuth_lines, range_samples));

        for i in 0..azimuth_lines {
            for j in 0..range_samples {
                let cal_coeff = self.interpolate_calibration_coefficient(
                    &self.coefficients.gamma0_lut, i, j
                )?;
                
                if cal_coeff > 0.0 {
                    calibrated[[i, j]] = intensity[[i, j]] / (cal_coeff * cal_coeff);
                } else {
                    calibrated[[i, j]] = 0.0;
                }
            }
        }

        Ok(calibrated)
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
}

/// Parse calibration data from Sentinel-1 annotation XML
pub fn parse_calibration_from_xml(xml_content: &str) -> SarResult<CalibrationCoefficients> {
    log::debug!("Parsing calibration data from XML");
    
    // This is a simplified implementation
    // In practice, you'd parse the complete calibration XML structure
    
    // Create dummy calibration coefficients for now
    let dummy_lut = Array2::from_elem((100, 1000), 1000.0);
    let dummy_times = Array1::linspace(0.0, 1.0, 100);
    
    Ok(CalibrationCoefficients {
        sigma0_lut: dummy_lut.clone(),
        beta0_lut: dummy_lut.clone(),
        gamma0_lut: dummy_lut.clone(),
        dn_lut: dummy_lut,
        range_times: dummy_times.clone(),
        azimuth_times: dummy_times,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use ndarray::Array2;

    #[test]
    fn test_calibration_processor() {
        let dummy_lut = Array2::from_elem((10, 10), 1000.0);
        let dummy_times = Array1::linspace(0.0, 1.0, 10);
        
        let coefficients = CalibrationCoefficients {
            sigma0_lut: dummy_lut.clone(),
            beta0_lut: dummy_lut.clone(),
            gamma0_lut: dummy_lut.clone(),
            dn_lut: dummy_lut,
            range_times: dummy_times.clone(),
            azimuth_times: dummy_times,
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
}
