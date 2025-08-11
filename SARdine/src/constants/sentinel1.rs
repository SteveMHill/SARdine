/// Sentinel-1 Satellite Constants
/// 
/// Official satellite parameters for Sentinel-1A and Sentinel-1B
/// Reference: ESA Sentinel-1 Mission Requirements Document (S1-RS-ESA-SY-0007)
/// 
/// IMPORTANT: These are the official design parameters. Actual orbit heights
/// vary slightly and precise values should be extracted from .EOF orbit files
/// when available.

/// Sentinel-1 orbital parameters
pub mod orbital {
    /// Nominal orbit height for Sentinel-1 (meters above Earth surface)
    /// Reference: ESA Sentinel-1 Mission Requirements Document
    /// Note: Actual height varies ±20km during orbit, use .EOF files for precision
    pub const NOMINAL_ORBIT_HEIGHT_M: f64 = 693_000.0;
    
    /// Orbit repeat cycle (days)
    pub const REPEAT_CYCLE_DAYS: u32 = 12;
    
    /// Orbital period (minutes)
    pub const ORBITAL_PERIOD_MIN: f64 = 98.6;
}

/// Sentinel-1 radar system parameters
pub mod radar {
    /// C-band center frequency (Hz)
    /// Reference: ESA Sentinel-1 Product Specification
    pub const CENTER_FREQUENCY_HZ: f64 = 5.405e9;
    
    /// Radar wavelength (meters) - calculated from center frequency
    /// λ = c/f where c = speed of light
    pub const WAVELENGTH_M: f64 = crate::constants::physical::SPEED_OF_LIGHT_M_S / CENTER_FREQUENCY_HZ;
    
    /// Range sampling rate (Hz) - ADC sampling frequency
    /// Reference: Sentinel-1 Product Specification (exact value)
    pub const RANGE_SAMPLING_RATE_HZ: f64 = 64_345_238.0957;
    
    /// Chirp bandwidth (Hz)
    pub const CHIRP_BANDWIDTH_HZ: f64 = 56.5e6;
}

/// Interferometric Wide (IW) mode parameters
pub mod iw_mode {
    /// Number of subswaths in IW mode
    pub const SUBSWATH_COUNT: u8 = 3;
    
    /// Typical incidence angle range (degrees)
    pub const INCIDENCE_ANGLE_MIN_DEG: f64 = 29.1;
    pub const INCIDENCE_ANGLE_MAX_DEG: f64 = 46.0;
    
    /// Swath width (km)
    pub const SWATH_WIDTH_KM: f64 = 250.0;
    
    /// Typical burst duration (seconds)
    pub const BURST_DURATION_S: f64 = 2.758277;
}

/// Data quality parameters
pub mod quality {
    /// Minimum valid backscatter value (dB) - below this is likely noise
    pub const MIN_VALID_BACKSCATTER_DB: f32 = -40.0;
    
    /// Maximum valid backscatter value (dB) - above this is likely clutter
    pub const MAX_VALID_BACKSCATTER_DB: f32 = 10.0;
    
    /// Minimum valid incidence angle (degrees)
    pub const MIN_VALID_INCIDENCE_DEG: f64 = 15.0;
    
    /// Maximum valid incidence angle (degrees)  
    pub const MAX_VALID_INCIDENCE_DEG: f64 = 60.0;
}

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_wavelength_calculation() {
        // Verify C-band wavelength is approximately 5.5 cm
        assert!((radar::WAVELENGTH_M - 0.0555).abs() < 0.001);
    }
    
    #[test]
    fn test_frequency_values() {
        // Verify center frequency is in C-band
        assert!(radar::CENTER_FREQUENCY_HZ > 5.0e9);
        assert!(radar::CENTER_FREQUENCY_HZ < 6.0e9);
    }
    
    #[test]
    fn test_orbital_parameters() {
        // Verify orbit height is in reasonable range for LEO
        assert!(orbital::NOMINAL_ORBIT_HEIGHT_M > 600_000.0);
        assert!(orbital::NOMINAL_ORBIT_HEIGHT_M < 800_000.0);
    }
}
