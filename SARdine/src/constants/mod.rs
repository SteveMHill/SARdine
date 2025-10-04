//! Scientific Constants Module
//!
//! Contains universally accepted physical, geodetic, and mission-specific constants
//! with full literature references for scientific reproducibility.

/// Physical constants (CODATA 2018)
pub mod physical;

/// Sentinel-1 mission-specific constants
pub mod sentinel1;

/// Fast unit conversion functions (emergency bottleneck fix)
pub mod unit_conversion;

// Re-export commonly used constants for convenience
pub use physical::geodetic; // Re-export geodetic module for backward compatibility
pub use physical::{EARTH_EQUATORIAL_RADIUS_M, EARTH_MEAN_RADIUS_M, SPEED_OF_LIGHT_M_S};
pub use sentinel1::{
    orbital::NOMINAL_ORBIT_HEIGHT_M,
    radar::{CENTER_FREQUENCY_HZ, RANGE_SAMPLING_RATE_HZ},
};

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
        assert!((geodetic::WGS84_FLATTENING - 1.0 / 298.257223563).abs() < 1e-12);
    }
}
