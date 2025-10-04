/// Physical Constants
///
/// Fundamental physical constants used in SAR processing
/// All values conform to CODATA 2018 recommendations where applicable
/// Electromagnetic constants
pub mod electromagnetic {
    /// Speed of light in vacuum (m/s)
    /// CODATA 2018 exact value
    pub const SPEED_OF_LIGHT_M_S: f64 = 299_792_458.0;

    /// Permittivity of free space (F/m)
    /// CODATA 2018 exact value
    pub const PERMITTIVITY_VACUUM_F_M: f64 = 8.8541878128e-12;

    /// Permeability of free space (H/m)
    /// CODATA 2018 exact value  
    pub const PERMEABILITY_VACUUM_H_M: f64 = 1.25663706212e-6;
}

/// Earth geodetic constants
pub mod geodetic {
    /// WGS84 Earth semi-major axis (meters)
    /// Reference: WGS 84 Implementation Manual 1991
    pub const WGS84_SEMI_MAJOR_AXIS_M: f64 = 6_378_137.0;

    /// WGS84 Earth semi-minor axis (meters)
    /// Calculated from semi-major axis and flattening
    pub const WGS84_SEMI_MINOR_AXIS_M: f64 = 6_356_752.314245;

    /// WGS84 first eccentricity squared
    pub const WGS84_ECCENTRICITY_SQUARED: f64 = 6.69437999014e-3;

    /// WGS84 flattening factor
    pub const WGS84_FLATTENING: f64 = 1.0 / 298.257223563;

    /// Earth mean radius (meters)
    /// Arithmetic mean of WGS84 semi-axes
    pub const EARTH_MEAN_RADIUS_M: f64 = (WGS84_SEMI_MAJOR_AXIS_M + WGS84_SEMI_MINOR_AXIS_M) / 2.0;

    /// Earth equatorial radius (meters) - same as semi-major axis
    pub const EARTH_EQUATORIAL_RADIUS_M: f64 = WGS84_SEMI_MAJOR_AXIS_M;
}

/// Mathematical constants
pub mod mathematical {
    /// Pi (π)
    pub const PI: f64 = std::f64::consts::PI;

    /// 2 * Pi
    pub const TWO_PI: f64 = 2.0 * PI;

    /// Pi / 2
    pub const PI_2: f64 = PI / 2.0;

    /// Pi / 180 - degrees to radians conversion
    pub const DEG_TO_RAD: f64 = PI / 180.0;

    /// 180 / Pi - radians to degrees conversion
    pub const RAD_TO_DEG: f64 = 180.0 / PI;

    /// Natural logarithm of 10
    pub const LN_10: f64 = std::f64::consts::LN_10;

    /// Base 10 logarithm of e
    pub const LOG10_E: f64 = std::f64::consts::LOG10_E;
}

/// Conversion factors
pub mod conversions {
    use super::mathematical::*;

    /// Convert dB to linear scale: linear = 10^(dB/10)
    pub fn db_to_linear(db: f64) -> f64 {
        10.0_f64.powf(db / 10.0)
    }

    /// Convert linear to dB scale: dB = 10*log10(linear)
    pub fn linear_to_db(linear: f64) -> f64 {
        10.0 * linear.log10()
    }

    /// Convert degrees to radians
    pub fn deg_to_rad(degrees: f64) -> f64 {
        degrees * DEG_TO_RAD
    }

    /// Convert radians to degrees
    pub fn rad_to_deg(radians: f64) -> f64 {
        radians * RAD_TO_DEG
    }
}

// Re-export commonly used constants for convenience
pub use electromagnetic::SPEED_OF_LIGHT_M_S;
pub use geodetic::{EARTH_EQUATORIAL_RADIUS_M, EARTH_MEAN_RADIUS_M, WGS84_SEMI_MAJOR_AXIS_M};
pub use mathematical::{DEG_TO_RAD, PI, RAD_TO_DEG, TWO_PI};

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_speed_of_light() {
        // Verify speed of light is correct
        assert_eq!(electromagnetic::SPEED_OF_LIGHT_M_S, 299_792_458.0);
    }

    // Removed test that asserted on constants to satisfy clippy::assertions-on-constants

    #[test]
    fn test_db_conversions() {
        // Test dB conversion functions
        assert!((conversions::db_to_linear(0.0) - 1.0).abs() < 1e-10);
        assert!((conversions::linear_to_db(1.0) - 0.0).abs() < 1e-10);
        assert!((conversions::db_to_linear(10.0) - 10.0).abs() < 1e-10);
    }

    #[test]
    fn test_angle_conversions() {
        // Test angle conversion functions
        assert!((conversions::deg_to_rad(180.0) - mathematical::PI).abs() < 1e-10);
        assert!((conversions::rad_to_deg(mathematical::PI) - 180.0).abs() < 1e-10);
    }
}
