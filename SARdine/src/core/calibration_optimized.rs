
//! Optimized calibration implementation with dense per-column caching
//! 
//! This module implements the calibration speed optimizations:
//! 2.1 Dense per‑column gain cache - interpolate LUTs once into gain[x] vectors  
//! 2.2 One pass over the image - fused complex scale, noise removal, dB conversion
//! 2.3 No binary searches in hot path - use pre-computed dense tables

use crate::types::{SarError, SarResult};
use ndarray::Array2;
use num_complex::Complex32;
use std::sync::atomic::{AtomicUsize, Ordering};

/// Dense per-column calibration cache for ultra-fast processing
#[derive(Debug, Clone)]
pub struct DenseCalibrationCache {
    /// Pre-interpolated β⁰ gains per column [width]
    pub beta_gains: Vec<f32>,
    /// Pre-interpolated σ⁰ gains per column [width]  
    pub sigma_gains: Vec<f32>,
    /// Pre-interpolated γ⁰ gains per column [width]
    pub gamma_gains: Vec<f32>,
    /// Pre-computed sin(θ) per column [width]
    pub sin_theta: Vec<f32>,
    /// Pre-computed cos(θ) per column [width] 
    pub cos_theta: Vec<f32>,
    /// Image width (number of range samples)
    pub width: usize,
    /// Cache validity flag
    pub is_valid: bool,
}

/// Fused calibration operation mode
#[derive(Debug, Clone, Copy)]
pub enum FusedCalibrationMode {
    /// Complex to power + calibration + noise removal
    ComplexToPowerCalibrated,
    /// Power + calibration only (noise already removed)
    PowerCalibrationOnly,
    /// Complex to power + calibration + noise removal + dB conversion
    ComplexToPowerCalibratedDB,
}

/// Optimized calibration processor with dense caching
pub struct OptimizedCalibrationProcessor {
    /// Dense per-column cache
    cache: Option<DenseCalibrationCache>,
    /// Processing statistics
    pixels_processed: AtomicUsize,
    cache_hits: AtomicUsize,
}

/// Calibration LUT interpolation parameters
#[derive(Debug, Clone)]
pub struct CalibrationLUTParams {
    /// Azimuth lines where LUT values are defined
    pub azimuth_lines: Vec<usize>,
    /// Range pixels where LUT values are defined  
    pub range_pixels: Vec<usize>,
    /// σ⁰ calibration values at tie points
    pub sigma_values: Vec<Vec<f64>>,
    /// β⁰ calibration values at tie points
    pub beta_values: Vec<Vec<f64>>,
    /// γ⁰ calibration values at tie points
    pub gamma_values: Vec<Vec<f64>>,
}

impl OptimizedCalibrationProcessor {
    /// Create new optimized processor
    pub fn new() -> Self {
        Self {
            cache: None,
            pixels_processed: AtomicUsize::new(0),
            cache_hits: AtomicUsize::new(0),
        }
    }

    /// 2.1 Build dense per-column gain cache with piecewise linear interpolation
    pub fn build_dense_cache(
        &mut self,
        lut_params: &CalibrationLUTParams,
        image_width: usize,
        image_height: usize,
    ) -> SarResult<()> {
        log::info!("🏗️ Building dense per-column calibration cache for {}x{} image", image_height, image_width);

        if lut_params.azimuth_lines.is_empty() || lut_params.range_pixels.is_empty() {
            return Err(SarError::Processing(
                "Cannot build dense cache: empty LUT parameters".to_string()
            ));
        }

        // Validate LUT dimensions
        let azimuth_count = lut_params.azimuth_lines.len();
        let range_count = lut_params.range_pixels.len();

        if lut_params.sigma_values.len() != azimuth_count ||
           lut_params.beta_values.len() != azimuth_count ||
           lut_params.gamma_values.len() != azimuth_count {
            return Err(SarError::Processing(
                "LUT azimuth dimension mismatch".to_string()
            ));
        }

        for i in 0..azimuth_count {
            if lut_params.sigma_values[i].len() != range_count ||
               lut_params.beta_values[i].len() != range_count ||
               lut_params.gamma_values[i].len() != range_count {
                return Err(SarError::Processing(
                    format!("LUT range dimension mismatch at azimuth index {}", i)
                ));
            }
        }

        // Pre-compute dense interpolation for each column
        let mut beta_gains = vec![0.0f32; image_width];
        let mut sigma_gains = vec![0.0f32; image_width];
        let mut gamma_gains = vec![0.0f32; image_width];
        let mut sin_theta = vec![0.0f32; image_width];
        let mut cos_theta = vec![0.0f32; image_width];

        // For simplicity, use scene-center azimuth line for interpolation
        let center_azimuth_idx = azimuth_count / 2;

        for col in 0..image_width {
            // Interpolate gains for this column using piecewise linear interpolation
            let (sigma_gain, beta_gain, gamma_gain) = self.interpolate_gains_for_column(
                lut_params,
                center_azimuth_idx,
                col as f64,
            );

            sigma_gains[col] = sigma_gain as f32;
            beta_gains[col] = beta_gain as f32;
            gamma_gains[col] = gamma_gain as f32;

            // Pre-compute incidence angle functions (simplified - use range-dependent model)
            let normalized_range = col as f32 / image_width as f32;
            let incidence_angle = 20.0 + 30.0 * normalized_range; // Typical range: 20° to 50°
            let theta_rad = incidence_angle.to_radians();
            
            sin_theta[col] = theta_rad.sin();
            cos_theta[col] = theta_rad.cos();
        }

        // Validate gains for scientific soundness
        Self::validate_dense_gains(&sigma_gains, "σ⁰")?;
        Self::validate_dense_gains(&beta_gains, "β⁰")?;
        Self::validate_dense_gains(&gamma_gains, "γ⁰")?;

        self.cache = Some(DenseCalibrationCache {
            beta_gains,
            sigma_gains,
            gamma_gains,
            sin_theta,
            cos_theta,
            width: image_width,
            is_valid: true,
        });

        log::info!("✅ Dense calibration cache built successfully");
        log::info!("   Cache size: {} columns × 5 arrays = {} elements", 
                  image_width, image_width * 5);

        Ok(())
    }

    /// Interpolate calibration gains for a specific column using piecewise linear
    fn interpolate_gains_for_column(
        &self,
        lut_params: &CalibrationLUTParams,
        azimuth_idx: usize,
        range_position: f64,
    ) -> (f64, f64, f64) {
        let range_pixels = &lut_params.range_pixels;
        let sigma_row = &lut_params.sigma_values[azimuth_idx];
        let beta_row = &lut_params.beta_values[azimuth_idx];
        let gamma_row = &lut_params.gamma_values[azimuth_idx];

        // Find range interpolation bracket
        let mut lower_idx = 0;
        let mut upper_idx = range_pixels.len() - 1;

        for i in 0..range_pixels.len() - 1 {
            if range_position >= range_pixels[i] as f64 && range_position <= range_pixels[i + 1] as f64 {
                lower_idx = i;
                upper_idx = i + 1;
                break;
            }
        }

        // Linear interpolation weight
        let range_lower = range_pixels[lower_idx] as f64;
        let range_upper = range_pixels[upper_idx] as f64;
        let weight = if range_upper > range_lower {
            (range_position - range_lower) / (range_upper - range_lower)
        } else {
            0.0
        };

        // Interpolate each gain type
        let sigma_gain = sigma_row[lower_idx] * (1.0 - weight) + sigma_row[upper_idx] * weight;
        let beta_gain = beta_row[lower_idx] * (1.0 - weight) + beta_row[upper_idx] * weight;
        let gamma_gain = gamma_row[lower_idx] * (1.0 - weight) + gamma_row[upper_idx] * weight;

        (sigma_gain, beta_gain, gamma_gain)
    }

    /// Validate dense gain arrays for scientific soundness
    fn validate_dense_gains(gains: &[f32], name: &str) -> SarResult<()> {
        if gains.is_empty() {
            return Err(SarError::Processing(
                format!("Empty {} gains array", name)
            ));
        }

        let finite_gains: Vec<f32> = gains.iter().cloned().filter(|x| x.is_finite()).collect();
        if finite_gains.is_empty() {
            return Err(SarError::Processing(
                format!("No finite {} gains found", name)
            ));
        }

        let min_gain = finite_gains.iter().fold(f32::INFINITY, |a, &b| a.min(b));
        let max_gain = finite_gains.iter().fold(f32::NEG_INFINITY, |a, &b| a.max(b));
        let mean_gain = finite_gains.iter().sum::<f32>() / finite_gains.len() as f32;

        // Check for flat arrays (parsing bugs)
        let variation_ratio = if mean_gain.abs() > 1e-10 {
            (max_gain - min_gain) / mean_gain.abs()
        } else {
            0.0
        };

        if variation_ratio < 0.01 {
            return Err(SarError::Processing(
                format!("{} gains show insufficient variation (ratio: {:.6}) - parsing bug suspected", 
                        name, variation_ratio)
            ));
        }

        // Check for reasonable value ranges
        if name == "σ⁰" || name == "β⁰" {
            if mean_gain < 1e-10 || mean_gain > 1e10 {
                return Err(SarError::Processing(
                    format!("{} gains have unrealistic mean: {:.3e}", name, mean_gain)
                ));
            }
        }

        log::info!("✅ {} gains validation passed:", name);
        log::info!("   Range: [{:.3e}, {:.3e}]", min_gain, max_gain);
        log::info!("   Mean: {:.3e}", mean_gain);
        log::info!("   Variation ratio: {:.3}", variation_ratio);

        Ok(())
    }

    /// 2.2 One pass over the image - fused operations
    pub fn apply_fused_calibration(
        &self,
        input_data: &Array2<Complex32>,
        noise_lut: Option<&Array2<f32>>,
        calibration_mode: FusedCalibrationMode,
        target_calibration: &str, // "sigma0", "beta0", "gamma0"
    ) -> SarResult<Array2<f32>> {
        let cache = self.cache.as_ref().ok_or_else(|| {
            SarError::Processing("Dense calibration cache not built".to_string())
        })?;

        if !cache.is_valid {
            return Err(SarError::Processing("Dense calibration cache is invalid".to_string()));
        }

        let (height, width) = input_data.dim();
        if width != cache.width {
            return Err(SarError::Processing(
                format!("Image width {} does not match cache width {}", width, cache.width)
            ));
        }

        log::info!("🚀 Applying fused calibration in single pass");
        log::info!("   Mode: {:?}", calibration_mode);
        log::info!("   Target: {}", target_calibration);
        log::info!("   Dimensions: {}x{}", height, width);

        let mut output = Array2::<f32>::zeros((height, width));

        // Select appropriate gain array based on target calibration
        let gains = match target_calibration.to_lowercase().as_str() {
            "sigma0" => &cache.sigma_gains,
            "beta0" => &cache.beta_gains,
            "gamma0" => &cache.gamma_gains,
            _ => return Err(SarError::Processing(
                format!("Unsupported calibration target: {}", target_calibration)
            )),
        };

        // VECTORIZED FUSED KERNEL: Single pass with all operations
        for row in 0..height {
            for col in 0..width {
                let complex_value = input_data[[row, col]];
                
                // Step 1: Complex to power (|I + jQ|²)
                let power_value = complex_value.norm_sqr();
                
                // Step 2: Thermal noise removal (if provided)
                let denoised_power = if let Some(noise) = noise_lut {
                    (power_value - noise[[row, col]]).max(0.0)
                } else {
                    power_value
                };
                
                // Step 3: Radiometric calibration (no table lookup - direct array access)
                let calibrated_linear = denoised_power * gains[col];
                
                // Step 4: Optional σ⁰/γ⁰ scaling using pre-computed angles
                let final_linear = if target_calibration == "gamma0" {
                    // γ⁰ = σ⁰ / cos(θ) scaling
                    calibrated_linear / cache.cos_theta[col]
                } else {
                    calibrated_linear
                };
                
                // Step 5: Optional dB conversion
                let final_value = match calibration_mode {
                    FusedCalibrationMode::ComplexToPowerCalibratedDB => {
                        if final_linear > 0.0 {
                            10.0 * final_linear.log10()
                        } else {
                            -50.0 // Noise floor in dB
                        }
                    },
                    _ => final_linear,
                };
                
                output[[row, col]] = final_value;
            }
        }

        // Update processing statistics
        let pixels_processed = height * width;
        self.pixels_processed.fetch_add(pixels_processed, Ordering::Relaxed);
        self.cache_hits.fetch_add(pixels_processed, Ordering::Relaxed); // All cache hits with dense array

        log::info!("✅ Fused calibration completed in single pass");
        log::info!("   Pixels processed: {}", pixels_processed);
        log::info!("   Cache hit rate: 100% (dense array access)");

        Ok(output)
    }

    /// 2.3 Alternative vectorized implementation for maximum performance
    /// Note: SIMD implementation disabled due to unstable feature requirement
    /*
    #[cfg(feature = "simd")]
    pub fn apply_fused_calibration_simd(
        &self,
        input_data: &Array2<Complex32>,
        noise_lut: Option<&Array2<f32>>,
        target_calibration: &str,
    ) -> SarResult<Array2<f32>> {
        // SIMD implementation would go here when portable_simd is stable
        // For now, fall back to regular implementation
        self.apply_fused_calibration(input_data, noise_lut, 
                                   FusedCalibrationMode::ComplexToPowerCalibrated, 
                                   target_calibration)
    }
    */

    /// Get processing statistics
    pub fn get_statistics(&self) -> (usize, usize, f64) {
        let processed = self.pixels_processed.load(Ordering::Relaxed);
        let hits = self.cache_hits.load(Ordering::Relaxed);
        let hit_rate = if processed > 0 {
            hits as f64 / processed as f64
        } else {
            0.0
        };
        (processed, hits, hit_rate)
    }

    /// Reset processing statistics
    pub fn reset_statistics(&self) {
        self.pixels_processed.store(0, Ordering::Relaxed);
        self.cache_hits.store(0, Ordering::Relaxed);
    }
}

impl Default for OptimizedCalibrationProcessor {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use num_complex::Complex32;

    #[test]
    fn test_dense_cache_validation() {
        let mut processor = OptimizedCalibrationProcessor::new();
        
        // Test with valid LUT parameters
        let lut_params = CalibrationLUTParams {
            azimuth_lines: vec![0, 100, 200],
            range_pixels: vec![0, 50, 100],
            sigma_values: vec![
                vec![1e-3, 2e-3, 3e-3],
                vec![1.1e-3, 2.1e-3, 3.1e-3],
                vec![1.2e-3, 2.2e-3, 3.2e-3],
            ],
            beta_values: vec![
                vec![2e-3, 4e-3, 6e-3],
                vec![2.1e-3, 4.1e-3, 6.1e-3],
                vec![2.2e-3, 4.2e-3, 6.2e-3],
            ],
            gamma_values: vec![
                vec![3e-3, 6e-3, 9e-3],
                vec![3.1e-3, 6.1e-3, 9.1e-3],
                vec![3.2e-3, 6.2e-3, 9.2e-3],
            ],
        };

        let result = processor.build_dense_cache(&lut_params, 100, 200);
        assert!(result.is_ok(), "Dense cache building should succeed with valid parameters");
        
        let cache = processor.cache.as_ref().unwrap();
        assert_eq!(cache.width, 100);
        assert!(cache.is_valid);
        assert_eq!(cache.sigma_gains.len(), 100);
    }

    #[test]
    fn test_fused_calibration() {
        let mut processor = OptimizedCalibrationProcessor::new();
        
        // Create simple test cache manually
        processor.cache = Some(DenseCalibrationCache {
            sigma_gains: vec![1e-3; 10],
            beta_gains: vec![2e-3; 10],
            gamma_gains: vec![3e-3; 10],
            sin_theta: vec![0.5; 10],
            cos_theta: vec![0.866; 10],
            width: 10,
            is_valid: true,
        });

        // Create test complex data
        let mut test_data = Array2::<Complex32>::zeros((5, 10));
        for i in 0..5 {
            for j in 0..10 {
                test_data[[i, j]] = Complex32::new(i as f32 + 1.0, j as f32 + 1.0);
            }
        }

        let result = processor.apply_fused_calibration(
            &test_data,
            None,
            FusedCalibrationMode::ComplexToPowerCalibrated,
            "sigma0",
        );

        assert!(result.is_ok(), "Fused calibration should succeed");
        let output = result.unwrap();
        assert_eq!(output.dim(), (5, 10));
        
        // Check that output values are reasonable (power * gain)
        assert!(output[[0, 0]] > 0.0, "Calibrated values should be positive");
    }
}