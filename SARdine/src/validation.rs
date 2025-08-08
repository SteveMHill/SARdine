//! Scientific Parameter Validation Framework
//!
//! Provides comprehensive validation for all SAR processing parameters
//! to ensure scientific accuracy and detect invalid/hardcoded values.

use crate::types::{SarResult, SarError};
use crate::constants::physical::SPEED_OF_LIGHT_M_S;
use std::collections::HashSet;

/// SAR parameter validation ranges based on scientific literature
pub struct ValidationRanges {
    /// Valid wavelength range for SAR systems (meters)
    pub wavelength_range: (f64, f64),
    /// Valid pixel spacing range (meters)  
    pub pixel_spacing_range: (f64, f64),
    /// Valid PRF range (Hz)
    pub prf_range: (f64, f64),
    /// Valid radar frequency range (Hz)
    pub frequency_range: (f64, f64),
}

impl Default for ValidationRanges {
    fn default() -> Self {
        Self {
            // SAR wavelength range: P-band (1m) to Ka-band (8mm)
            wavelength_range: (0.008, 1.0),
            // Pixel spacing: sub-meter to 100m
            pixel_spacing_range: (0.1, 100.0),
            // PRF range: 100 Hz to 10 kHz (typical SAR range)
            prf_range: (100.0, 10_000.0),
            // Radar frequency range: 300 MHz to 40 GHz
            frequency_range: (3e8, 4e10),
        }
    }
}

/// Comprehensive SAR parameter validator
pub struct ParameterValidator {
    ranges: ValidationRanges,
}

impl ParameterValidator {
    pub fn new() -> Self {
        Self {
            ranges: ValidationRanges::default(),
        }
    }
    
    /// Validate wavelength parameter
    pub fn validate_wavelength(&self, wavelength: f64, source: &str) -> SarResult<()> {
        if wavelength < self.ranges.wavelength_range.0 || wavelength > self.ranges.wavelength_range.1 {
            return Err(SarError::InvalidParameter(format!(
                "Wavelength {:.6}m from {} is outside valid SAR range [{:.3}-{:.3}]m", 
                wavelength, source, self.ranges.wavelength_range.0, self.ranges.wavelength_range.1
            )));
        }
        
        // Check for suspicious hardcoded values
        let suspicious_values = [0.055, 0.0555, 0.055465763, 0.05546576, 0.23];
        for &suspicious in &suspicious_values {
            if (wavelength - suspicious).abs() < 1e-8 {
                return Err(SarError::InvalidParameter(format!(
                    "CRITICAL: Wavelength {:.9}m matches known hardcoded value. Must extract from annotation XML radar frequency.", 
                    wavelength
                )));
            }
        }
        
        Ok(())
    }
    
    /// Validate pixel spacing parameters
    pub fn validate_pixel_spacing(&self, range_spacing: f64, azimuth_spacing: f64, source: &str) -> SarResult<()> {
        // Validate range spacing
        if range_spacing < self.ranges.pixel_spacing_range.0 || range_spacing > self.ranges.pixel_spacing_range.1 {
            return Err(SarError::InvalidParameter(format!(
                "Range pixel spacing {:.6}m from {} is outside valid range [{:.1}-{:.1}]m", 
                range_spacing, source, self.ranges.pixel_spacing_range.0, self.ranges.pixel_spacing_range.1
            )));
        }
        
        // Validate azimuth spacing  
        if azimuth_spacing < self.ranges.pixel_spacing_range.0 || azimuth_spacing > self.ranges.pixel_spacing_range.1 {
            return Err(SarError::InvalidParameter(format!(
                "Azimuth pixel spacing {:.6}m from {} is outside valid range [{:.1}-{:.1}]m", 
                azimuth_spacing, source, self.ranges.pixel_spacing_range.0, self.ranges.pixel_spacing_range.1
            )));
        }
        
        // Check for suspicious hardcoded values
        let suspicious_range = [2.3, 2.33, 2.329562, 2.8767];
        let suspicious_azimuth = [14.0, 14.1, 14.059906];
        
        for &suspicious in &suspicious_range {
            if (range_spacing - suspicious).abs() < 1e-6 {
                return Err(SarError::InvalidParameter(format!(
                    "CRITICAL: Range spacing {:.6}m matches hardcoded value {:.6}. Must extract from annotation XML.", 
                    range_spacing, suspicious
                )));
            }
        }
        
        for &suspicious in &suspicious_azimuth {
            if (azimuth_spacing - suspicious).abs() < 1e-6 {
                return Err(SarError::InvalidParameter(format!(
                    "CRITICAL: Azimuth spacing {:.6}m matches hardcoded value {:.6}. Must extract from annotation XML.", 
                    azimuth_spacing, suspicious
                )));
            }
        }
        
        Ok(())
    }
    
    /// Validate radar frequency and derived wavelength consistency
    pub fn validate_frequency_wavelength_consistency(&self, frequency: f64, wavelength: f64) -> SarResult<()> {
        let calculated_wavelength = SPEED_OF_LIGHT_M_S / frequency;
        let difference = (wavelength - calculated_wavelength).abs();
        let relative_error = difference / calculated_wavelength;
        
        if relative_error > 1e-9 {
            return Err(SarError::InvalidParameter(format!(
                "Wavelength {:.9}m inconsistent with frequency {:.3} Hz. Expected wavelength: {:.9}m (error: {:.2e})", 
                wavelength, frequency, calculated_wavelength, relative_error
            )));
        }
        
        Ok(())
    }
    
    /// Comprehensive validation of all SAR parameters
    pub fn validate_all_parameters(
        &self,
        frequency: f64,
        wavelength: f64, 
        range_spacing: f64,
        azimuth_spacing: f64,
        prf: f64,
        source: &str
    ) -> SarResult<()> {
        // Validate individual parameters
        self.validate_wavelength(wavelength, source)?;
        self.validate_pixel_spacing(range_spacing, azimuth_spacing, source)?;
        
        // Validate frequency
        if frequency < self.ranges.frequency_range.0 || frequency > self.ranges.frequency_range.1 {
            return Err(SarError::InvalidParameter(format!(
                "Radar frequency {:.3} Hz from {} is outside valid range [{:.0}-{:.0}] Hz", 
                frequency, source, self.ranges.frequency_range.0, self.ranges.frequency_range.1
            )));
        }
        
        // Validate PRF
        if prf < self.ranges.prf_range.0 || prf > self.ranges.prf_range.1 {
            return Err(SarError::InvalidParameter(format!(
                "PRF {:.3} Hz from {} is outside valid range [{:.0}-{:.0}] Hz", 
                prf, source, self.ranges.prf_range.0, self.ranges.prf_range.1
            )));
        }
        
        // Validate consistency between frequency and wavelength
        self.validate_frequency_wavelength_consistency(frequency, wavelength)?;
        
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_hardcoded_wavelength_detection() {
        let validator = ParameterValidator::new();
        
        // Should detect hardcoded wavelength values
        assert!(validator.validate_wavelength(0.055, "test").is_err());
        assert!(validator.validate_wavelength(0.0555, "test").is_err());
        assert!(validator.validate_wavelength(0.055465763, "test").is_err());
    }
    
    #[test]
    fn test_hardcoded_spacing_detection() {
        let validator = ParameterValidator::new();
        
        // Should detect hardcoded spacing values
        assert!(validator.validate_pixel_spacing(2.3, 10.0, "test").is_err());
        assert!(validator.validate_pixel_spacing(2.329562, 10.0, "test").is_err());
        assert!(validator.validate_pixel_spacing(10.0, 14.0, "test").is_err());
        assert!(validator.validate_pixel_spacing(10.0, 14.059906, "test").is_err());
    }
    
    #[test]
    fn test_valid_parameters() {
        let validator = ParameterValidator::new();
        
        // Valid C-band parameters (frequency: 5.405 GHz)
        let frequency = 5.405e9;
        let wavelength = SPEED_OF_LIGHT_M_S / frequency; // ~0.055465763m
        
        // Should accept valid parameters that don't match hardcoded values
        assert!(validator.validate_all_parameters(
            frequency,
            wavelength,
            2.5, // Slightly different from hardcoded 2.3
            13.8, // Slightly different from hardcoded 14.0
            1000.0,
            "annotation XML"
        ).is_ok());
    }
}

/// Scientific wavelength validation
pub fn validate_wavelength_scientific(wavelength: f64) -> SarResult<()> {
    let ranges = ValidationRanges::default();
    if wavelength >= ranges.wavelength_range.0 && wavelength <= ranges.wavelength_range.1 {
        Ok(())
    } else {
        Err(SarError::InvalidParameter(
            format!("Wavelength {:.6}m outside valid SAR range [{:.3}-{:.1}]m", 
                   wavelength, ranges.wavelength_range.0, ranges.wavelength_range.1)
        ))
    }
}

/// Realistic pixel spacing validation
pub fn validate_pixel_spacing_realistic(range_spacing: f64, azimuth_spacing: f64) -> SarResult<()> {
    let ranges = ValidationRanges::default();
    
    if range_spacing < ranges.pixel_spacing_range.0 || range_spacing > ranges.pixel_spacing_range.1 {
        return Err(SarError::InvalidParameter(
            format!("Range pixel spacing {:.3}m outside valid range [{:.1}-{:.0}]m",
                   range_spacing, ranges.pixel_spacing_range.0, ranges.pixel_spacing_range.1)
        ));
    }
    
    if azimuth_spacing < ranges.pixel_spacing_range.0 || azimuth_spacing > ranges.pixel_spacing_range.1 {
        return Err(SarError::InvalidParameter(
            format!("Azimuth pixel spacing {:.3}m outside valid range [{:.1}-{:.0}]m",
                   azimuth_spacing, ranges.pixel_spacing_range.0, ranges.pixel_spacing_range.1)
        ));
    }
    
    Ok(())
}

/// SAR frequency band validation
pub fn validate_frequency_sar_band(frequency: f64) -> SarResult<()> {
    let ranges = ValidationRanges::default();
    
    if frequency < ranges.frequency_range.0 || frequency > ranges.frequency_range.1 {
        return Err(SarError::InvalidParameter(
            format!("Radar frequency {:.0} Hz outside SAR band [{:.0}-{:.0}] Hz",
                   frequency, ranges.frequency_range.0, ranges.frequency_range.1)
        ));
    }
    
    // Check for common SAR bands
    let c_band = (4.0e9, 8.0e9);  // C-band
    let x_band = (8.0e9, 12.0e9); // X-band
    let l_band = (1.0e9, 2.0e9);  // L-band
    
    if (frequency >= c_band.0 && frequency <= c_band.1) ||
       (frequency >= x_band.0 && frequency <= x_band.1) ||
       (frequency >= l_band.0 && frequency <= l_band.1) {
        Ok(())
    } else {
        Err(SarError::InvalidParameter(
            format!("Frequency {:.3} GHz not in common SAR bands (L/C/X)", frequency / 1e9)
        ))
    }
}

/// Detect hardcoded values that should come from annotations
pub fn validate_no_hardcoded_values(
    wavelength: f64, 
    range_spacing: f64, 
    azimuth_spacing: f64
) -> SarResult<()> {
    // Known hardcoded values that should be eliminated
    let forbidden_wavelengths = [0.055, 0.0555, 0.055465763];
    let forbidden_spacings = [2.3, 14.0, 2.329580, 14.065834];
    
    for &forbidden in &forbidden_wavelengths {
        if (wavelength - forbidden).abs() < 1e-6 {
            return Err(SarError::InvalidParameter(
                format!("Detected hardcoded wavelength {:.6}m - must extract from annotation", wavelength)
            ));
        }
    }
    
    for &forbidden in &forbidden_spacings {
        if (range_spacing - forbidden).abs() < 1e-6 || (azimuth_spacing - forbidden).abs() < 1e-6 {
            return Err(SarError::InvalidParameter(
                format!("Detected hardcoded pixel spacing {:.6}m - must extract from annotation", forbidden)
            ));
        }
    }
    
    Ok(())
}

/// Detect suspicious constants that might be hardcoded
pub fn detect_suspicious_constants(values: &[f64]) -> Vec<String> {
    let mut warnings = Vec::new();
    let mut value_counts = std::collections::HashMap::new();
    
    // Count occurrences of values
    for &value in values {
        *value_counts.entry(format!("{:.6}", value)).or_insert(0) += 1;
    }
    
    // Flag values that appear too frequently (might be hardcoded)
    for (value_str, count) in value_counts {
        if count > 10 {
            warnings.push(format!("Value {} appears {} times - possible hardcoded constant", value_str, count));
        }
    }
    
    warnings
}
