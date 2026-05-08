//! Processing Pipeline Constants
//!
//! Contains constants used throughout the SAR processing pipeline.
//! These constants replace magic numbers for better maintainability
//! and scientific traceability.
//!
//! Reference: SARdine Comprehensive Audit Report, Issue #25

/// Bounding box margin constants for terrain correction
pub mod bbox {
    /// Default margin to add around SAR footprint for DEM extraction (degrees)
    /// Reduced from 1km (~0.01°) to 100m (~0.001°) to avoid including excessive
    /// ocean pixels while maintaining numerical stability at boundaries.
    pub const MARGIN_DEG: f64 = 0.001;

    /// Minimum valid bounding box dimension (degrees)
    /// Prevents degenerate bbox with zero or near-zero extent
    pub const MIN_DIMENSION_DEG: f64 = 0.001;
}

/// Coverage validation thresholds
pub mod coverage {
    /// Critical minimum valid pixel percentage - below this, terrain correction fails
    /// SCIENTIFIC FIX: Raised from 1.0% to 75.0% - geocoded SAR should have >75% coverage
    pub const CRITICAL_THRESHOLD_PERCENT: f64 = 75.0;

    /// Warning threshold - below this, emit warning about low coverage  
    /// SCIENTIFIC FIX: Raised from 50.0% to 90.0%
    pub const WARNING_THRESHOLD_PERCENT: f64 = 90.0;

    /// Minimum row coverage fraction for merge grid clipping
    pub const ROW_COVERAGE_MIN_FRAC: f32 = 0.5;
}

/// Multilooking constants
pub mod multilook {
    /// Default multilook factor when metadata is unavailable
    pub const DEFAULT_FACTOR: f64 = 1.0;

    /// Minimum valid multilook factor
    pub const MIN_FACTOR: f64 = 1.0;

    /// Maximum valid multilook factor (prevent excessive downsampling)
    pub const MAX_FACTOR: f64 = 50.0;
}

/// Orbit processing constants
pub mod orbit {
    /// Minimum required orbit state vectors for accurate processing
    pub const MINIMUM_VECTORS: usize = 10;

    /// Time tolerance for orbit interpolation (seconds)
    pub const INTERPOLATION_TOLERANCE_S: f64 = 0.001;

    /// Velocity sanity check range (m/s) - Sentinel-1 should be ~7500-7700 m/s
    pub const MIN_VELOCITY_MPS: f64 = 7400.0;
    pub const MAX_VELOCITY_MPS: f64 = 7800.0;
}

/// Terrain correction constants
pub mod terrain_correction {
    /// Maximum Newton-Raphson iterations for range-Doppler positioning
    pub const MAX_ITERATIONS: u32 = 50;

    /// Convergence tolerance for Newton-Raphson (pixels)
    pub const CONVERGENCE_TOLERANCE: f64 = 0.001;

    /// Maximum DEM elevation for sanity check (meters)
    /// Mt. Everest is ~8849m, add margin for safety
    pub const MAX_VALID_ELEVATION_M: f64 = 9000.0;

    /// Minimum DEM elevation for sanity check (meters)
    /// Dead Sea shore is ~-430m, add margin for safety
    pub const MIN_VALID_ELEVATION_M: f64 = -500.0;
}

/// Speckle filtering constants
pub mod speckle {
    /// Default window size for speckle filter
    pub const DEFAULT_WINDOW_SIZE: usize = 7;

    /// Minimum valid window size
    pub const MIN_WINDOW_SIZE: usize = 3;

    /// Maximum valid window size
    pub const MAX_WINDOW_SIZE: usize = 21;

    /// Equivalent Number of Looks (ENL) for multi-looked data quality assessment
    pub const DEFAULT_ENL: f64 = 4.4;
}

/// Calibration constants
pub mod calibration {
    /// Valid sigma0 range in linear units (not dB)
    /// Corresponds roughly to -40 dB to +10 dB
    pub const MIN_VALID_SIGMA0_LINEAR: f32 = 1e-4;
    pub const MAX_VALID_SIGMA0_LINEAR: f32 = 10.0;

    /// Thermal noise floor estimate (linear units)
    pub const NOISE_FLOOR_LINEAR: f32 = 1e-6;
}

/// dB conversion constants
pub mod db_conversion {
    /// Expected minimum dB value for SAR backscatter
    pub const EXPECTED_MIN_DB: f32 = -50.0;

    /// Expected maximum dB value for SAR backscatter
    pub const EXPECTED_MAX_DB: f32 = 15.0;

    /// Warning threshold - emit warning if values below this
    pub const WARN_MIN_DB: f32 = -40.0;

    /// Warning threshold - emit warning if values above this
    pub const WARN_MAX_DB: f32 = 10.0;

    /// Epsilon value to prevent log(0) in dB conversion
    pub const LOG_EPSILON: f32 = 1e-12;
}

/// Merge/overlap processing constants
pub mod merge {
    /// Maximum seam fill gap (pixels) for thin-gap interpolation
    pub const SEAM_FILL_MAX_GAP: usize = 8;

    /// Minimum fraction of row hits for seam detection
    pub const SEAM_FILL_MIN_ROW_HIT_FRAC: f32 = 0.5;

    /// Minimum fraction of column hits for seam detection
    pub const SEAM_FILL_MIN_COL_HIT_FRAC: f32 = 0.5;

    /// Sampling step for overlap radiometric check (pixels)
    pub const OVERLAP_SAMPLE_STEP_RADIO: usize = 16;

    /// Fine sampling step for overlap phase estimation (pixels)
    pub const OVERLAP_SAMPLE_STEP_FINE: usize = 8;
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_coverage_thresholds_ordered() {
        assert!(coverage::CRITICAL_THRESHOLD_PERCENT < coverage::WARNING_THRESHOLD_PERCENT);
    }

    #[test]
    fn test_multilook_range_valid() {
        assert!(multilook::MIN_FACTOR <= multilook::DEFAULT_FACTOR);
        assert!(multilook::DEFAULT_FACTOR <= multilook::MAX_FACTOR);
    }

    #[test]
    fn test_db_range_valid() {
        assert!(db_conversion::EXPECTED_MIN_DB < db_conversion::EXPECTED_MAX_DB);
        assert!(db_conversion::WARN_MIN_DB < db_conversion::WARN_MAX_DB);
    }

    #[test]
    fn test_elevation_range_valid() {
        assert!(
            terrain_correction::MIN_VALID_ELEVATION_M < terrain_correction::MAX_VALID_ELEVATION_M
        );
    }
}
