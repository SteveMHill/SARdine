use chrono::{DateTime, Utc};
use ndarray::{Array2, Array3};
use num_complex::Complex;
use pyo3::prelude::*;
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

/// Time reference system for SAR processing
///
/// CRITICAL: SAR processing uses multiple time references that must be tracked carefully:
/// - **OrbitEpoch**: Relative to first orbit state vector (for orbit interpolation)
/// - **ProductStart**: Relative to product acquisition start time (for image indexing)
/// - **UtcAbsolute**: Absolute UTC timestamp (for cross-product comparison)
///
/// Mixing these references causes geolocation errors of ~kilometers!
///
/// # References
/// - ESA S1-TN-ESA-GP-0028: "Sentinel-1 Time and Orbit Reference"
/// - Suggestion from external review (2025-10-04): Explicit time reference tracking
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum TimeReference {
    /// Time relative to orbit state vector reference epoch (seconds)
    /// Used for: Orbit interpolation, Newton-Raphson solver
    OrbitEpoch,
    
    /// Time relative to product start time (seconds)
    /// Used for: Image pixel indexing, burst timing
    ProductStart,
    
    /// Absolute UTC timestamp (seconds since Unix epoch)
    /// Used for: Cross-product comparison, annotation parsing
    UtcAbsolute,
}

impl std::fmt::Display for TimeReference {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            TimeReference::OrbitEpoch => write!(f, "OrbitEpoch"),
            TimeReference::ProductStart => write!(f, "ProductStart"),
            TimeReference::UtcAbsolute => write!(f, "UtcAbsolute"),
        }
    }
}

/// Time value with explicit reference tracking
///
/// This structure ensures we never mix time references accidentally.
/// All time conversions must be explicit and traceable.
///
/// # Example
/// ```ignore
/// let orbit_time = TimeValue::new(42.5, TimeReference::OrbitEpoch);
/// let utc_time = orbit_time.to_utc(orbit_reference_epoch);
/// ```
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct TimeValue {
    /// Time value in seconds
    pub value: f64,
    
    /// Reference system for this time value
    pub reference: TimeReference,
}

impl TimeValue {
    /// Create a new time value with explicit reference
    pub fn new(value: f64, reference: TimeReference) -> Self {
        Self { value, reference }
    }
    
    /// Convert to absolute UTC timestamp (seconds since Unix epoch)
    pub fn to_utc_seconds(&self, orbit_epoch: DateTime<Utc>, product_start: DateTime<Utc>) -> f64 {
        match self.reference {
            TimeReference::UtcAbsolute => self.value,
            TimeReference::OrbitEpoch => {
                let epoch_seconds = orbit_epoch.timestamp() as f64
                    + (orbit_epoch.timestamp_subsec_nanos() as f64) * 1e-9;
                epoch_seconds + self.value
            }
            TimeReference::ProductStart => {
                let start_seconds = product_start.timestamp() as f64
                    + (product_start.timestamp_subsec_nanos() as f64) * 1e-9;
                start_seconds + self.value
            }
        }
    }
    
    /// Convert to orbit-relative time (seconds since orbit epoch)
    pub fn to_orbit_relative(&self, orbit_epoch: DateTime<Utc>, product_start: DateTime<Utc>) -> f64 {
        let utc_seconds = self.to_utc_seconds(orbit_epoch, product_start);
        let epoch_seconds = orbit_epoch.timestamp() as f64
            + (orbit_epoch.timestamp_subsec_nanos() as f64) * 1e-9;
        utc_seconds - epoch_seconds
    }
    
    /// Convert to product-relative time (seconds since product start)
    pub fn to_product_relative(&self, orbit_epoch: DateTime<Utc>, product_start: DateTime<Utc>) -> f64 {
        let utc_seconds = self.to_utc_seconds(orbit_epoch, product_start);
        let start_seconds = product_start.timestamp() as f64
            + (product_start.timestamp_subsec_nanos() as f64) * 1e-9;
        utc_seconds - start_seconds
    }
    
    /// Validate time value consistency
    pub fn validate(&self) -> Result<(), String> {
        if !self.value.is_finite() {
            return Err(format!("Time value must be finite: {}", self.value));
        }
        
        // Check reasonable bounds based on reference
        match self.reference {
            TimeReference::UtcAbsolute => {
                // Should be Unix timestamp (> 1e9 for dates after ~2001)
                if self.value < 1e9 {
                    log::warn!(
                        "⚠️  UTC absolute time looks suspicious: {:.3}s (expected > 1e9 for recent dates)",
                        self.value
                    );
                }
            }
            TimeReference::OrbitEpoch | TimeReference::ProductStart => {
                // Relative times should be reasonable for SAR scenes (< 1 hour typical)
                if self.value.abs() > 3600.0 {
                    log::warn!(
                        "⚠️  Relative time looks suspicious: {:.3}s (expected < 3600s for typical scenes)",
                        self.value
                    );
                }
            }
        }
        
        Ok(())
    }
}

// ============================================================================
// TIME CONVERSION UTILITIES
// ============================================================================

/// Convert chrono DateTime to UTC seconds since Unix epoch
pub fn datetime_to_utc_seconds(dt: chrono::DateTime<chrono::Utc>) -> f64 {
    dt.timestamp() as f64 + dt.timestamp_subsec_nanos() as f64 * 1e-9
}

/// Convert UTC seconds since Unix epoch to chrono DateTime
pub fn utc_seconds_to_datetime(seconds: f64) -> chrono::DateTime<chrono::Utc> {
    let secs = seconds.floor() as i64;
    let nanos = ((seconds - secs as f64) * 1e9) as u32;
    chrono::DateTime::from_timestamp(secs, nanos)
        .unwrap_or_else(|| chrono::DateTime::from_timestamp(0, 0).unwrap())
}

/// Timing metadata for SAR products with explicit time references
#[derive(Debug, Clone)]
pub struct ProductTiming {
    /// Product start time (absolute UTC)
    pub product_start_utc: chrono::DateTime<chrono::Utc>,
    
    /// Orbit reference epoch (absolute UTC) - CRITICAL for orbit interpolation
    pub orbit_epoch_utc: chrono::DateTime<chrono::Utc>,
    
    /// Azimuth time interval (PRF period in seconds)
    pub azimuth_time_interval: f64,
    
    /// Product duration (seconds)
    pub product_duration: f64,
}

impl ProductTiming {
    /// Create a TimeValue for a pixel's azimuth time relative to product start
    pub fn pixel_time_from_line(&self, line_index: usize) -> TimeValue {
        let relative_time = line_index as f64 * self.azimuth_time_interval;
        TimeValue {
            value: relative_time,
            reference: TimeReference::ProductStart,
        }
    }
    
    /// Convert product-relative time to orbit-relative time
    pub fn product_to_orbit_time(&self, product_time: f64) -> TimeValue {
        let product_start_orbit_relative = datetime_to_utc_seconds(self.product_start_utc)
            - datetime_to_utc_seconds(self.orbit_epoch_utc);
        
        TimeValue {
            value: product_start_orbit_relative + product_time,
            reference: TimeReference::OrbitEpoch,
        }
    }
    
    /// Convert orbit-relative time to product-relative time
    pub fn orbit_to_product_time(&self, orbit_time: f64) -> TimeValue {
        let product_start_orbit_relative = datetime_to_utc_seconds(self.product_start_utc)
            - datetime_to_utc_seconds(self.orbit_epoch_utc);
        
        TimeValue {
            value: orbit_time - product_start_orbit_relative,
            reference: TimeReference::ProductStart,
        }
    }
}

/// Coordinate system for SAR data
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
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
#[pyclass]
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
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub enum AcquisitionMode {
    IW, // Interferometric Wide swath
    EW, // Extra Wide swath
    SM, // StripMap
    WV, // Wave
}

impl std::fmt::Display for AcquisitionMode {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            AcquisitionMode::IW => write!(f, "IW"),
            AcquisitionMode::EW => write!(f, "EW"),
            AcquisitionMode::SM => write!(f, "SM"),
            AcquisitionMode::WV => write!(f, "WV"),
        }
    }
}

/// Sub-swath information for IW mode
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SubSwath {
    pub id: String, // IW1, IW2, IW3
    pub burst_count: usize,

    // Local sub-swath dimensions
    pub range_samples: usize,
    pub azimuth_samples: usize,

    // Global coordinate reference frame preservation
    // These maintain the full-image coordinates for proper LUT mapping
    pub first_line_global: usize, // First line in full-image coordinates
    pub last_line_global: usize,  // Last line in full-image coordinates
    pub first_sample_global: usize, // First sample in full-image coordinates
    pub last_sample_global: usize, // Last sample in full-image coordinates

    // Physical parameters
    pub range_pixel_spacing: f64,   // meters
    pub azimuth_pixel_spacing: f64, // meters
    pub slant_range_time: f64,      // seconds
    pub burst_duration: f64,        // seconds

    // EXPERT ADDITION: Annotation-derived timing parameters
    pub prf_hz: Option<f64>, // Pulse Repetition Frequency from annotation (Hz)
    
    // DC-AWARE ALIGNMENT: Doppler centroid polynomial for inter-subswath alignment
    // CRITICAL: Must be populated by DC-aware deburst, otherwise TOPSAR merge will fail
    pub dc_polynomial: Option<Vec<f64>>, // Doppler centroid coefficients [c0, c1, c2, ...]
    pub azimuth_time_interval: Option<f64>, // seconds per line (for timing-based alignment)
}

/// Orbit state vector
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct StateVector {
    pub time: DateTime<Utc>,
    pub position: [f64; 3], // [x, y, z] in meters
    pub velocity: [f64; 3], // [vx, vy, vz] in m/s
}

/// Precise orbit information
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct OrbitData {
    pub state_vectors: Vec<StateVector>,
    pub reference_time: DateTime<Utc>,
}

impl OrbitData {
    /// Validate orbit data for time precision and monotonicity
    /// This catches time-related issues that can break interpolation algorithms
    pub fn validate_time_precision(&self) -> Result<(), String> {
        use chrono::{DateTime, Utc};

        // Inline time conversion helpers for validation
        let dt_to_f64 = |dt: DateTime<Utc>| -> f64 {
            dt.timestamp() as f64 + (dt.timestamp_subsec_nanos() as f64) * 1e-9
        };

        if self.state_vectors.is_empty() {
            return Err("No orbit state vectors available".to_string());
        }

        // Check monotonicity
        let mut prev_time = f64::NEG_INFINITY;
        for (i, sv) in self.state_vectors.iter().enumerate() {
            let t = dt_to_f64(sv.time);
            if t <= prev_time {
                return Err(format!(
                    "Orbit times not strictly increasing at index {}: {} <= {}",
                    i, t, prev_time
                ));
            }
            prev_time = t;
        }

        // Calculate precision metrics
        let time_deltas: Vec<f64> = self
            .state_vectors
            .windows(2)
            .map(|w| dt_to_f64(w[1].time) - dt_to_f64(w[0].time))
            .collect();

        if let Some(min_dt) = time_deltas
            .iter()
            .cloned()
            .fold(None, |acc, x| Some(acc.map_or(x, |y| x.min(y))))
        {
            log::info!(
                "✅ Orbit validation: {} state vectors",
                self.state_vectors.len()
            );
            if !self.state_vectors.is_empty() {
                log::info!(
                    "✅ Orbit time span: {:.6} to {:.6} s",
                    dt_to_f64(
                        self.state_vectors
                            .first()
                            .expect("Vector checked for non-empty")
                            .time
                    ),
                    dt_to_f64(
                        self.state_vectors
                            .last()
                            .expect("Vector checked for non-empty")
                            .time
                    )
                );
            }
            log::info!(
                "✅ Orbit min Δt: {:.6} s (expect ~10 s for S1 precise orbits)",
                min_dt
            );

            // Warn if time precision might be insufficient
            if min_dt < 1e-6 {
                log::warn!(
                    "⚠️  Very small orbit time intervals detected: {:.9} s",
                    min_dt
                );
            }
            if min_dt > 20.0 {
                log::warn!("⚠️  Large orbit time intervals detected: {:.1} s", min_dt);
            }
        }

        Ok(())
    }
}

/// Orbit parameters for SAR processing (simplified from OrbitData)
#[derive(Debug, Clone)]
pub struct OrbitParams {
    pub states: Vec<StateVector>,
    pub polynomial_degree: usize,
}

/// Geospatial bounding box
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BoundingBox {
    pub min_lon: f64,
    pub max_lon: f64,
    pub min_lat: f64,
    pub max_lat: f64,
}

impl BoundingBox {
    /// Normalize the bounding box so that min <= max for both lat and lon.
    /// If inversion is detected it will either:
    ///  - Swap and log an error (default behavior)
    ///  - Return an error (if SARDINE_STRICT_BBOX=1)
    pub fn normalize(&mut self) {
        let strict = std::env::var("SARDINE_STRICT_BBOX").ok().as_deref() == Some("1");

        if self.min_lat > self.max_lat {
            if strict {
                panic!(
                    "BoundingBox lat inverted (min_lat > max_lat): {} > {}. Set SARDINE_STRICT_BBOX=0 to auto-correct.",
                    self.min_lat, self.max_lat
                );
            } else {
                log::error!(
                    "⚠️ BoundingBox latitude inverted: min_lat ({}) > max_lat ({}). Swapping to correct.",
                    self.min_lat, self.max_lat
                );
                std::mem::swap(&mut self.min_lat, &mut self.max_lat);
            }
        }
        if self.min_lon > self.max_lon {
            if strict {
                panic!(
                    "BoundingBox lon inverted (min_lon > max_lon): {} > {}. Set SARDINE_STRICT_BBOX=0 to auto-correct.",
                    self.min_lon, self.max_lon
                );
            } else {
                log::error!(
                    "⚠️ BoundingBox longitude inverted: min_lon ({}) > max_lon ({}). Swapping to correct.",
                    self.min_lon, self.max_lon
                );
                std::mem::swap(&mut self.min_lon, &mut self.max_lon);
            }
        }
    }

    /// Return a normalized copy (does not modify self)
    pub fn normalized(mut self) -> Self {
        self.normalize();
        self
    }
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

    // Radar parameters (extracted from annotation XML)
    pub radar_frequency: Option<f64>, // Hz - must be extracted from annotation, never hardcoded
    pub wavelength: Option<f64>,      // meters - must be extracted from annotation
    pub slant_range_time: Option<f64>, // seconds - critical for terrain correction
    pub prf: Option<f64>, // Hz - pulse repetition frequency, critical for terrain correction
    pub range_sampling_rate: Option<f64>, // Hz - range sampling rate, critical for coordinate conversion
    #[serde(default)]
    pub radar_frequency_extracted: bool, // true when value was parsed from annotation XML

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

    #[error("Invalid parameter: {0}")]
    InvalidParameter(String),

    #[error("Parameter error: {0}")]
    ParameterError(String),

    #[error("Processing error: {0}")]
    Processing(String),

    #[error("Data processing error: {0}")]
    DataProcessingError(String),

    #[error("Metadata error: {0}")]
    Metadata(String),

    #[error("Invalid metadata: {0}")]
    InvalidMetadata(String),

    #[error("Missing calibration data: {0}")]
    MissingCalibrationData(String),

    #[error("GDAL error: {0}")]
    Gdal(#[from] gdal::errors::GdalError),

    #[error("XML parsing error: {0}")]
    XmlParsing(String),

    #[error("Missing required parameter: {0}")]
    MissingParameter(String),

    #[error("Missing metadata: {0}")]
    MissingMetadata(String),

    #[error("Numerical computation error: {0}")]
    NumericalError(String),

    #[error("Not implemented: {0}")]
    NotImplemented(String),
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
        wavelength: f64,          // radar wavelength (C-band ≈ 0.055 m)
    ) -> Option<f64> {
        let velocity = self.get_velocity_at_line(line_idx)?;

        // Doppler frequency = 2 * (v_sat · look_dir) / λ
        let velocity_dot_look = velocity[0] * look_direction[0]
            + velocity[1] * look_direction[1]
            + velocity[2] * look_direction[2];

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
    /// Local incidence angle cosine threshold (dimensionless, 0-1)
    pub lia_threshold: f64,
    /// DEM validity threshold in meters
    pub dem_threshold: f64,
    /// Minimum gamma0 threshold in power units (linear scale, not dB)
    pub gamma0_min: f32,
    /// Maximum gamma0 threshold in power units (linear scale, not dB)
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
            lia_threshold: 0.6, // Fixed: Changed from 0.1 (cos 84°) to 0.6 (cos 53°) for reasonable SAR geometry
            dem_threshold: -100.0,
            // POWER DOMAIN THRESHOLDS: Masking now operates in power units for better precision
            // Made more liberal to preserve more data by default
            gamma0_min: 1e-6, // Equivalent to -60 dB: 10^(-60/10) = 1e-6 (very liberal minimum)
            gamma0_max: 100000.0, // Equivalent to 50 dB: 10^(50/10) = 100000 (very liberal maximum)
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
