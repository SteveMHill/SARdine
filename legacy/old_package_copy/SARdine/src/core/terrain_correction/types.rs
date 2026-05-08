use std::str::FromStr;

use super::config::TerrainCorrectionConfig;

/// Explicit geographic coordinate type to ensure axis order consistency
///
/// Based on expert recommendations to prevent lat/lon confusion throughout the codebase.
/// This enforces the (latitude, longitude) convention consistently.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct LatLon {
    pub lat: f64, // Latitude in degrees [-90, 90]
    pub lon: f64, // Longitude in degrees [-180, 180]
}

impl LatLon {
    /// Create new geographic coordinates with validation
    pub fn new(lat: f64, lon: f64) -> crate::types::SarResult<Self> {
        if lat < -90.0 || lat > 90.0 {
            return Err(crate::types::SarError::Processing(format!(
                "Invalid latitude: {} (must be [-90, 90])",
                lat
            )));
        }
        if lon < -180.0 || lon > 180.0 {
            return Err(crate::types::SarError::Processing(format!(
                "Invalid longitude: {} (must be [-180, 180])",
                lon
            )));
        }
        Ok(LatLon { lat, lon })
    }

    /// Create from tuple ensuring correct axis order
    pub fn from_tuple(coords: (f64, f64)) -> crate::types::SarResult<Self> {
        Self::new(coords.0, coords.1)
    }

    /// Convert to tuple as (lat, lon)
    pub fn to_tuple(self) -> (f64, f64) {
        (self.lat, self.lon)
    }
}

/// 3D position vector
#[derive(Debug, Clone)]
pub struct Position3D {
    pub x: f64,
    pub y: f64,
    pub z: f64,
}

/// 3D velocity vector
#[derive(Debug, Clone)]
pub struct Velocity3D {
    pub x: f64,
    pub y: f64,
    pub z: f64,
}

/// Simple 3D vector for geometric calculations
#[derive(Debug, Clone)]
pub struct Vector3 {
    pub x: f64,
    pub y: f64,
    pub z: f64,
}

/// Ground point with coordinates and elevation
#[derive(Debug, Clone)]
pub struct GroundPoint {
    pub latitude: f64,
    pub longitude: f64,
    pub elevation: f64,
    pub x: f64, // Projected X coordinate
    pub y: f64, // Projected Y coordinate
}

/// Geocoding result for a single pixel
#[derive(Debug, Clone)]
pub struct GeocodedPixel {
    pub sar_range: usize,
    pub sar_azimuth: usize,
    pub ground_point: GroundPoint,
    pub value: f32,
    pub valid: bool,
}

/// Algorithm execution status tracking
#[derive(Debug, Clone)]
pub struct AlgorithmStatus {
    pub algorithm_name: String,
    pub execution_mode: ExecutionMode,
    pub iterations_used: Option<usize>,
    pub convergence_achieved: Option<bool>,
    pub fallback_reason: Option<String>,
    pub processing_time_ms: f64,
}

#[derive(Debug, Clone)]
pub enum ExecutionMode {
    Primary,
    Fallback(String),
    Failed(String),
}

/// Processing metadata for scientific reproducibility
#[derive(Debug, Clone)]
pub struct ProcessingMetadata {
    pub algorithm_statuses: Vec<AlgorithmStatus>,
    pub configuration_used: TerrainCorrectionConfig,
    pub input_validation_results: ValidationResults,
}

#[derive(Debug, Clone)]
pub struct ValidationResults {
    pub bounding_box_valid: bool,
    pub elevation_range_valid: bool,
    pub coordinate_system_valid: bool,
    pub orbit_data_valid: bool,
    pub warnings: Vec<String>,
    pub errors: Vec<String>,
}

/// Interpolation methods for SAR data resampling
#[derive(Debug, Clone, Copy, Default)]
pub enum InterpolationMethod {
    Nearest,
    #[default]
    Bilinear,
    Bicubic,
    /// Sinc interpolation (windowed sinc function) - GAMMA standard
    Sinc,
    /// Lanczos interpolation (Lanczos kernel) - high quality
    Lanczos,
}

impl FromStr for InterpolationMethod {
    type Err = ();
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        Ok(match s.to_lowercase().as_str() {
            "nearest" => InterpolationMethod::Nearest,
            "bilinear" => InterpolationMethod::Bilinear,
            "bicubic" => InterpolationMethod::Bicubic,
            "sinc" => InterpolationMethod::Sinc,
            "lanczos" => InterpolationMethod::Lanczos,
            _ => InterpolationMethod::Bilinear,
        })
    }
}
