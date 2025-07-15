use chrono::{DateTime, Utc};
use ndarray::{Array2, Array3};
use num_complex::Complex;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;

/// Complex-valued SAR data type (I + jQ)
pub type SarComplex = Complex<f32>;

/// Real-valued intensity or amplitude data
pub type SarReal = f32;

/// 2D complex SAR data array (range x azimuth)
pub type SarImage = Array2<SarComplex>;

/// 2D real SAR data array (range x azimuth)
pub type SarRealImage = Array2<SarReal>;

/// 3D SAR data for multi-polarization (pol x range x azimuth)
pub type SarCube = Array3<SarComplex>;

/// Coordinate system enumeration
#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum CoordinateSystem {
    /// Radar coordinates (range, azimuth)
    Radar,
    /// Geographic coordinates (latitude, longitude)
    Geographic,
    /// Projected coordinates (e.g., UTM)
    Projected { epsg: u32 },
}

/// Polarization modes for Sentinel-1
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum Polarization {
    VV,
    VH,
    HV,
    HH,
}

impl std::fmt::Display for Polarization {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Polarization::VV => write!(f, "VV"),
            Polarization::VH => write!(f, "VH"),
            Polarization::HV => write!(f, "HV"),
            Polarization::HH => write!(f, "HH"),
        }
    }
}

/// Sentinel-1 acquisition mode
#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum AcquisitionMode {
    IW, // Interferometric Wide swath
    EW, // Extra Wide swath
    SM, // StripMap
    WV, // Wave
}

/// Sub-swath information for IW mode
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SubSwath {
    pub id: String,                    // IW1, IW2, IW3
    pub burst_count: usize,
    pub range_samples: usize,
    pub azimuth_samples: usize,
    pub range_pixel_spacing: f64,      // meters
    pub azimuth_pixel_spacing: f64,    // meters
    pub slant_range_time: f64,         // seconds
    pub burst_duration: f64,           // seconds
}

/// Orbit state vector
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct StateVector {
    pub time: DateTime<Utc>,
    pub position: [f64; 3],  // [x, y, z] in meters
    pub velocity: [f64; 3],  // [vx, vy, vz] in m/s
}

/// Precise orbit information
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct OrbitData {
    pub state_vectors: Vec<StateVector>,
    pub reference_time: DateTime<Utc>,
}

/// Geospatial bounding box
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BoundingBox {
    pub min_lon: f64,
    pub max_lon: f64,
    pub min_lat: f64,
    pub max_lat: f64,
}

/// Geospatial transformation parameters
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GeoTransform {
    pub top_left_x: f64,
    pub pixel_width: f64,
    pub rotation_x: f64,
    pub top_left_y: f64,
    pub rotation_y: f64,
    pub pixel_height: f64,
}

/// Complete SAR product metadata
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SarMetadata {
    // Product identification
    pub product_id: String,
    pub mission: String,
    pub platform: String,
    pub instrument: String,
    
    // Acquisition parameters
    pub acquisition_mode: AcquisitionMode,
    pub polarizations: Vec<Polarization>,
    pub start_time: DateTime<Utc>,
    pub stop_time: DateTime<Utc>,
    
    // Geometry
    pub bounding_box: BoundingBox,
    pub coordinate_system: CoordinateSystem,
    pub sub_swaths: HashMap<String, SubSwath>,
    
    // Orbit
    pub orbit_data: Option<OrbitData>,
    
    // Processing parameters
    pub range_looks: u32,
    pub azimuth_looks: u32,
    pub pixel_spacing: (f64, f64), // (range, azimuth) in meters
}

/// Processing result with data and metadata
#[derive(Debug, Clone)]
pub struct SarProduct {
    pub metadata: SarMetadata,
    pub data: HashMap<Polarization, SarRealImage>,
    pub coordinate_system: CoordinateSystem,
    pub geo_transform: Option<GeoTransform>,
}

/// Error types for SAR processing
#[derive(Debug, thiserror::Error)]
pub enum SarError {
    #[error("I/O error: {0}")]
    Io(#[from] std::io::Error),
    
    #[error("Invalid data format: {0}")]
    InvalidFormat(String),
    
    #[error("Processing error: {0}")]
    Processing(String),
    
    #[error("Metadata error: {0}")]
    Metadata(String),
    
    #[error("GDAL error: {0}")]
    Gdal(#[from] gdal::errors::GdalError),
    
    #[error("XML parsing error: {0}")]
    XmlParsing(String),
}

/// Result type for SAR operations
pub type SarResult<T> = Result<T, SarError>;
