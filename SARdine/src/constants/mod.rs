//! Scientific Constants Module
//! 
//! Contains only universally accepted physical and geodetic constants
//! with full literature references for scientific reproducibility.
//!
//! IMPORTANT: Sensor-specific parameters (wavelength, pixel spacing, etc.)
//! MUST be extracted from annotation XML files - NOT defined here!

/// Physical Constants
/// Reference: NIST/CODATA 2018 values
pub mod physical {
    /// Speed of light in vacuum (exact, by definition)
    /// Reference: SI base unit definition
    pub const SPEED_OF_LIGHT_M_S: f64 = 299_792_458.0;
}

/// Geodetic Constants  
/// Reference: WGS84 (World Geodetic System 1984)
pub mod geodetic {
    /// WGS84 Earth semi-major axis (meters)
    /// Reference: NIMA TR8350.2, Department of Defense World Geodetic System 1984
    pub const WGS84_SEMI_MAJOR_AXIS_M: f64 = 6_378_137.0;
    
    /// WGS84 Earth flattening factor
    /// Reference: WGS84 definition
    pub const WGS84_FLATTENING: f64 = 1.0 / 298.257223563;
}

/// Mathematical Constants
pub mod math {
    /// Pi (high precision)
    pub const PI: f64 = std::f64::consts::PI;
    
    /// 2 * Pi
    pub const TWO_PI: f64 = 2.0 * std::f64::consts::PI;
}

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_physical_constants() {
        // Verify speed of light matches NIST definition
        assert_eq!(physical::SPEED_OF_LIGHT_M_S, 299_792_458.0);
    }
    
    #[test]
    fn test_geodetic_constants() {
        // Verify WGS84 parameters
        assert_eq!(geodetic::WGS84_SEMI_MAJOR_AXIS_M, 6_378_137.0);
        assert!((geodetic::WGS84_FLATTENING - 1.0/298.257223563).abs() < 1e-12);
    }
}
