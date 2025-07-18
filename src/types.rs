use chrono::{DateTime, Utc};
use ndarray::{Array2, Array3};
use num_complex::Complex;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::path::PathBuf;

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
    
    #[error("Invalid input: {0}")]
    InvalidInput(String),
    
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

/// Orbit file availability status
#[derive(Debug, Clone)]
pub struct OrbitStatus {
    /// Product identifier
    pub product_id: String,
    /// Product acquisition start time
    pub start_time: DateTime<Utc>,
    /// Whether orbit data is embedded in SLC
    pub has_embedded: bool,
    /// Primary recommended orbit type
    pub primary_orbit_type: crate::io::orbit::OrbitType,
    /// Whether primary orbit is cached locally
    pub has_primary_cached: bool,
    /// Fallback orbit type
    pub fallback_orbit_type: crate::io::orbit::OrbitType,
    /// Whether fallback orbit is cached locally
    pub has_fallback_cached: bool,
    /// Path to orbit cache directory
    pub cache_dir: PathBuf,
}

impl OrbitStatus {
    /// Check if any orbit data is available (embedded or cached)
    pub fn has_any_orbit(&self) -> bool {
        self.has_embedded || self.has_primary_cached || self.has_fallback_cached
    }
    
    /// Get recommended action for obtaining orbit data
    pub fn recommended_action(&self) -> OrbitAction {
        if self.has_embedded {
            OrbitAction::UseEmbedded
        } else if self.has_primary_cached {
            OrbitAction::UseCached(self.primary_orbit_type)
        } else if self.has_fallback_cached {
            OrbitAction::UseCached(self.fallback_orbit_type)
        } else {
            OrbitAction::Download(self.primary_orbit_type)
        }
    }
}

/// Recommended action for orbit data
#[derive(Debug, Clone, Copy)]
pub enum OrbitAction {
    /// Use orbit data embedded in SLC
    UseEmbedded,
    /// Use cached orbit file of specified type
    UseCached(crate::io::orbit::OrbitType),
    /// Download orbit file of specified type
    Download(crate::io::orbit::OrbitType),
}

/// Burst-specific orbit data for pixel-level satellite position/velocity
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BurstOrbitData {
    /// Satellite positions for each azimuth line [ECEF coordinates in meters]
    pub positions: Vec<[f64; 3]>,
    /// Satellite velocities for each azimuth line [ECEF in m/s]
    pub velocities: Vec<[f64; 3]>,
    /// Azimuth times for each line
    pub azimuth_times: Vec<DateTime<Utc>>,
    /// Burst start time
    pub burst_start_time: DateTime<Utc>,
    /// Time interval between azimuth samples (seconds)
    pub azimuth_time_interval: f64,
}

impl BurstOrbitData {
    /// Get satellite position for a specific azimuth line
    pub fn get_position_at_line(&self, line_idx: usize) -> Option<[f64; 3]> {
        self.positions.get(line_idx).copied()
    }
    
    /// Get satellite velocity for a specific azimuth line
    pub fn get_velocity_at_line(&self, line_idx: usize) -> Option<[f64; 3]> {
        self.velocities.get(line_idx).copied()
    }
    
    /// Get azimuth time for a specific line
    pub fn get_azimuth_time_at_line(&self, line_idx: usize) -> Option<DateTime<Utc>> {
        self.azimuth_times.get(line_idx).copied()
    }
    
    /// Calculate Doppler centroid for a given line and range direction
    pub fn calculate_doppler_at_line(
        &self,
        line_idx: usize,
        look_direction: [f64; 3], // unit vector towards target
        wavelength: f64, // radar wavelength (C-band ≈ 0.055 m)
    ) -> Option<f64> {
        let velocity = self.get_velocity_at_line(line_idx)?;
        
        // Doppler frequency = 2 * (v_sat · look_dir) / λ
        let velocity_dot_look = velocity[0] * look_direction[0] +
                               velocity[1] * look_direction[1] +
                               velocity[2] * look_direction[2];
        
        Some(2.0 * velocity_dot_look / wavelength)
    }
    
    /// Get number of azimuth lines
    pub fn num_lines(&self) -> usize {
        self.positions.len()
    }
}

/// EOF orbit file header information
#[derive(Debug, Clone, Default)]
pub struct EofHeader {
    pub coordinate_system: Option<String>,
    pub time_reference: Option<String>,
    pub file_name: Option<String>,
}

/// Surface normal vector for terrain analysis
#[derive(Debug, Clone, Copy, Default)]
pub struct SurfaceNormal {
    pub x: f64,
    pub y: f64,
    pub z: f64,
}

impl SurfaceNormal {
    pub fn new(x: f64, y: f64, z: f64) -> Self {
        Self { x, y, z }
    }
    
    /// Normalize the vector to unit length
    pub fn normalize(&mut self) {
        let length = (self.x * self.x + self.y * self.y + self.z * self.z).sqrt();
        if length > 0.0 {
            self.x /= length;
            self.y /= length;
            self.z /= length;
        }
    }
    
    /// Calculate dot product with another vector
    pub fn dot(&self, other: &[f64; 3]) -> f64 {
        self.x * other[0] + self.y * other[1] + self.z * other[2]
    }
}

/// Masking workflow configuration
#[derive(Debug, Clone)]
pub struct MaskingWorkflow {
    pub water_mask: bool,
    pub shadow_mask: bool,
    pub layover_mask: bool,
    pub noise_mask: bool,
    pub coherence_threshold: Option<f64>,
    pub intensity_threshold: Option<f64>,
    pub lia_threshold: f64,
    pub dem_threshold: f64,
    pub gamma0_min: f32,
    pub gamma0_max: f32,
}

impl Default for MaskingWorkflow {
    fn default() -> Self {
        Self {
            water_mask: true,
            shadow_mask: true,
            layover_mask: true,
            noise_mask: false,
            coherence_threshold: Some(0.3),
            intensity_threshold: None,
            lia_threshold: 0.1,
            dem_threshold: -100.0,
            gamma0_min: -50.0,
            gamma0_max: 10.0,
        }
    }
}

/// Result of masking operations
#[derive(Debug, Clone)]
pub struct MaskResult {
    pub water_mask: Option<Array2<u8>>,
    pub shadow_mask: Option<Array2<u8>>,
    pub layover_mask: Option<Array2<u8>>,
    pub noise_mask: Option<Array2<u8>>,
    pub combined_mask: Array2<u8>,
    pub lia_cosine: Array2<f32>,
    pub gamma0_mask: Array2<bool>,
    pub dem_mask: Array2<bool>,
    pub lia_mask: Array2<bool>,
    pub valid_pixels: usize,
    pub total_pixels: usize,
    pub coverage_percent: f64,
    pub stats: MaskStats,
}

/// Statistics from masking operations
#[derive(Debug, Clone, Default)]
pub struct MaskStats {
    pub total_pixels: usize,
    pub water_pixels: usize,
    pub shadow_pixels: usize,
    pub layover_pixels: usize,
    pub noise_pixels: usize,
    pub valid_pixels: usize,
}
