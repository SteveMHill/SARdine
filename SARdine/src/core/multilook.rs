use crate::types::{SarError, SarResult};
use ndarray::Array2;

/// Multilooking parameters for speckle reduction
#[derive(Debug, Clone)]
pub struct MultilookParams {
    /// Number of looks in range direction
    pub range_looks: usize,
    /// Number of looks in azimuth direction  
    pub azimuth_looks: usize,
    /// Output pixel spacing in meters (optional, for georeferencing)
    pub output_pixel_spacing: Option<f64>,
}

impl Default for MultilookParams {
    fn default() -> Self {
        Self {
            range_looks: 4,
            azimuth_looks: 1,
            output_pixel_spacing: None,
        }
    }
}

/// Multilook processor for reducing speckle in SAR imagery
pub struct MultilookProcessor {
    params: MultilookParams,
}

impl MultilookProcessor {
    /// Create a new multilook processor
    pub fn new(params: MultilookParams) -> Self {
        Self { params }
    }

    /// Create processor with standard parameters
    pub fn standard() -> Self {
        Self::new(MultilookParams::default())
    }

    /// Apply multilooking to intensity data
    /// 
    /// # Arguments
    /// * `intensity_data` - 2D array of intensity values (sigma0, beta0, etc.)
    /// * `range_spacing` - Original range pixel spacing in meters
    /// * `azimuth_spacing` - Original azimuth pixel spacing in meters
    /// 
    /// # Returns
    /// * Multilooked intensity data and new pixel spacings
    pub fn apply_multilook(
        &self,
        intensity_data: &Array2<f32>,
        range_spacing: f64,
        azimuth_spacing: f64,
    ) -> SarResult<(Array2<f32>, f64, f64)> {
        let (rows, cols) = intensity_data.dim();
        
        log::info!("Applying multilook: {}x{} looks to {}x{} image", 
                   self.params.azimuth_looks, self.params.range_looks, rows, cols);
        
        // Calculate output dimensions
        let out_rows = rows / self.params.azimuth_looks;
        let out_cols = cols / self.params.range_looks;
        
        if out_rows == 0 || out_cols == 0 {
            return Err(SarError::Processing(
                "Multilook parameters too large for input image".to_string(),
            ));
        }
        
        log::info!("Output dimensions: {}x{}", out_rows, out_cols);
        
        // Create output array
        let mut output = Array2::<f32>::zeros((out_rows, out_cols));
        
        // Apply multilooking using averaging
        for out_row in 0..out_rows {
            for out_col in 0..out_cols {
                let mut sum = 0.0f32;
                let mut count = 0;
                
                // Average over the multilook window
                for az_offset in 0..self.params.azimuth_looks {
                    for rg_offset in 0..self.params.range_looks {
                        let in_row = out_row * self.params.azimuth_looks + az_offset;
                        let in_col = out_col * self.params.range_looks + rg_offset;
                        
                        if in_row < rows && in_col < cols {
                            sum += intensity_data[[in_row, in_col]];
                            count += 1;
                        }
                    }
                }
                
                // Store the average
                if count > 0 {
                    output[[out_row, out_col]] = sum / (count as f32);
                }
            }
        }
        
        // Calculate new pixel spacings
        let new_range_spacing = range_spacing * (self.params.range_looks as f64);
        let new_azimuth_spacing = azimuth_spacing * (self.params.azimuth_looks as f64);
        
        log::info!("Multilooking complete: pixel spacing {}m x {}m -> {}m x {}m",
                   range_spacing, azimuth_spacing, new_range_spacing, new_azimuth_spacing);
        
        Ok((output, new_range_spacing, new_azimuth_spacing))
    }

    /// Apply multilooking with spatial filtering (boxcar average)
    /// 
    /// This version applies a more sophisticated multilooking that considers
    /// the exact window for each output pixel
    pub fn apply_multilook_filtered(
        &self,
        intensity_data: &Array2<f32>,
        range_spacing: f64,
        azimuth_spacing: f64,
    ) -> SarResult<(Array2<f32>, f64, f64)> {
        let (rows, cols) = intensity_data.dim();
        
        log::info!("Applying filtered multilook: {}x{} looks to {}x{} image", 
                   self.params.azimuth_looks, self.params.range_looks, rows, cols);
        
        // Calculate output dimensions
        let out_rows = rows / self.params.azimuth_looks;
        let out_cols = cols / self.params.range_looks;
        
        if out_rows == 0 || out_cols == 0 {
            return Err(SarError::Processing(
                "Multilook parameters too large for input image".to_string(),
            ));
        }
        
        // Create output array
        let mut output = Array2::<f32>::zeros((out_rows, out_cols));
        
        // Apply multilooking with proper normalization
        for out_row in 0..out_rows {
            for out_col in 0..out_cols {
                let mut sum = 0.0f64; // Use f64 for better precision in accumulation
                let mut count = 0;
                
                // Define the input window for this output pixel
                let start_row = out_row * self.params.azimuth_looks;
                let end_row = ((out_row + 1) * self.params.azimuth_looks).min(rows);
                let start_col = out_col * self.params.range_looks;
                let end_col = ((out_col + 1) * self.params.range_looks).min(cols);
                
                // Average over the exact window
                for in_row in start_row..end_row {
                    for in_col in start_col..end_col {
                        sum += intensity_data[[in_row, in_col]] as f64;
                        count += 1;
                    }
                }
                
                // Store the average
                if count > 0 {
                    output[[out_row, out_col]] = (sum / (count as f64)) as f32;
                }
            }
        }
        
        // Calculate new pixel spacings
        let new_range_spacing = range_spacing * (self.params.range_looks as f64);
        let new_azimuth_spacing = azimuth_spacing * (self.params.azimuth_looks as f64);
        
        log::info!("Filtered multilooking complete: {}x{} -> {}x{}, spacing: {}m x {}m",
                   rows, cols, out_rows, out_cols, new_range_spacing, new_azimuth_spacing);
        
        Ok((output, new_range_spacing, new_azimuth_spacing))
    }

    /// Calculate equivalent number of looks (ENL) estimate
    /// 
    /// This provides a quality metric for the multilooking
    pub fn estimate_enl(&self, data: &Array2<f32>) -> f32 {
        let mean = data.mean().unwrap_or(0.0);
        let variance = data.mapv(|x| (x - mean).powi(2)).mean().unwrap_or(0.0);
        
        if variance > 1e-10 {  // Use small threshold to avoid division by zero
            mean * mean / variance
        } else {
            f32::MAX  // Return maximum value for uniform data (infinite ENL)
        }
    }

    /// Get the theoretical number of looks
    pub fn theoretical_looks(&self) -> usize {
        self.params.range_looks * self.params.azimuth_looks
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use ndarray::Array;

    #[test]
    fn test_multilook_basic() {
        // Create a 4x4 test image
        let data = Array::from_shape_vec(
            (4, 4),
            vec![
                1.0, 2.0, 3.0, 4.0,
                5.0, 6.0, 7.0, 8.0,
                9.0, 10.0, 11.0, 12.0,
                13.0, 14.0, 15.0, 16.0,
            ],
        ).unwrap();
        
        let params = MultilookParams {
            range_looks: 2,
            azimuth_looks: 2,
            output_pixel_spacing: None,
        };
        
        let processor = MultilookProcessor::new(params);
        let (result, new_range, new_azimuth) = processor
            .apply_multilook(&data, 10.0, 10.0)
            .unwrap();
        
        // Should be 2x2 output
        assert_eq!(result.dim(), (2, 2));
        
        // Check values (should be averages of 2x2 blocks)
        assert_eq!(result[[0, 0]], 3.5); // (1+2+5+6)/4
        assert_eq!(result[[0, 1]], 5.5); // (3+4+7+8)/4
        assert_eq!(result[[1, 0]], 11.5); // (9+10+13+14)/4
        assert_eq!(result[[1, 1]], 13.5); // (11+12+15+16)/4
        
        // Check new spacings
        assert_eq!(new_range, 20.0);
        assert_eq!(new_azimuth, 20.0);
    }

    #[test]
    fn test_multilook_asymmetric() {
        // Test with different looks in each direction
        let data = Array::from_shape_vec(
            (2, 4),
            vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0],
        ).unwrap();
        
        let params = MultilookParams {
            range_looks: 2,
            azimuth_looks: 1,
            output_pixel_spacing: None,
        };
        
        let processor = MultilookProcessor::new(params);
        let (result, _, _) = processor
            .apply_multilook(&data, 10.0, 10.0)
            .unwrap();
        
        // Should be 2x2 output
        assert_eq!(result.dim(), (2, 2));
        
        // Check values
        assert_eq!(result[[0, 0]], 1.5); // (1+2)/2
        assert_eq!(result[[0, 1]], 3.5); // (3+4)/2
        assert_eq!(result[[1, 0]], 5.5); // (5+6)/2
        assert_eq!(result[[1, 1]], 7.5); // (7+8)/2
    }

    #[test]
    fn test_enl_calculation() {
        // Create uniform data (should have high ENL)
        let uniform_data = Array2::<f32>::from_elem((100, 100), 1.0);
        let processor = MultilookProcessor::standard();
        let enl = processor.estimate_enl(&uniform_data);
        
        // Uniform data should have very high ENL (near infinity, but numerically limited)
        assert!(enl > 1000.0);
    }
}
