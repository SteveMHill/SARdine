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
    pub fn to_orbit_relative(
        &self,
        orbit_epoch: DateTime<Utc>,
        product_start: DateTime<Utc>,
    ) -> f64 {
        let utc_seconds = self.to_utc_seconds(orbit_epoch, product_start);
        let epoch_seconds =
            orbit_epoch.timestamp() as f64 + (orbit_epoch.timestamp_subsec_nanos() as f64) * 1e-9;
        utc_seconds - epoch_seconds
    }

    /// Convert to product-relative time (seconds since product start)
    pub fn to_product_relative(
        &self,
        orbit_epoch: DateTime<Utc>,
        product_start: DateTime<Utc>,
    ) -> f64 {
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

/// Azimuth FM-rate polynomial sampled along azimuth time.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct FmRateEstimate {
    /// Absolute azimuth time in UTC for this polynomial (if provided).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub azimuth_time_utc: Option<DateTime<Utc>>,
    /// Azimuth time relative to the orbit reference epoch (seconds).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub azimuth_time_rel_orbit: Option<f64>,
    /// Slant-range reference time τ₀ used by the polynomial (seconds).
    pub slant_range_reference_time: f64,
    /// Polynomial coefficients evaluated as f_m(τ) = Σ c_n (τ - τ₀)^n.
    pub coefficients: Vec<f64>,
}

/// Sub-swath information for IW mode
///
/// # Bound Semantics
/// All bounds use **exclusive end** convention: `[first, last)` means `first..last` in Rust.
///
/// | Field | Convention | Description |
/// |-------|------------|-------------|
/// | `first_line_global` | inclusive | First line in full-image coordinates |
/// | `last_line_global` | **exclusive** | One past the last line |
/// | `first_sample_global` | inclusive | First sample in full-image coordinates |
/// | `last_sample_global` | **exclusive** | One past the last sample |
/// | `valid_*` fields | same pattern | inclusive start, exclusive end |
///
/// # Example
/// ```ignore
/// // A subswath covering lines 0-999 (1000 lines) and samples 100-24999 (24900 samples):
/// SubSwath {
///     first_line_global: 0,
///     last_line_global: 1000,     // exclusive: lines 0..1000
///     first_sample_global: 100,
///     last_sample_global: 25000,  // exclusive: samples 100..25000
///     azimuth_samples: 1000,      // = last_line_global - first_line_global
///     lines_per_burst: 125,       // per-burst azimuth extent (TOPS burst height)
///     range_samples: 24900,       // = last_sample_global - first_sample_global
///     // ...
/// }
/// ```
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SubSwath {
    /// Subswath identifier (IW1, IW2, IW3)
    pub id: String,
    pub burst_count: usize,
    /// Number of azimuth lines per burst (TOPS burst height)
    #[serde(default)]
    pub lines_per_burst: usize,

    // Local sub-swath dimensions
    /// Number of range samples in this subswath
    pub range_samples: usize,
    /// Number of azimuth samples (lines) in this subswath
    pub azimuth_samples: usize,

    // Global coordinate reference frame preservation
    /// First line in full-image coordinates (inclusive)
    pub first_line_global: usize,
    /// Last line in full-image coordinates (exclusive: one past the last line)
    pub last_line_global: usize,
    /// First sample in full-image coordinates (inclusive)
    pub first_sample_global: usize,
    /// Last sample in full-image coordinates (exclusive: one past the last sample)
    pub last_sample_global: usize,
    /// Total native samples including invalid margins
    pub full_range_samples: usize,
    /// First valid line (inclusive)
    pub valid_first_line: Option<usize>,
    /// Last valid line (exclusive: one past the last valid line)
    pub valid_last_line: Option<usize>,
    /// First valid sample (inclusive)
    pub valid_first_sample: Option<usize>,
    /// Last valid sample (exclusive: one past the last valid sample)
    pub valid_last_sample: Option<usize>,

    // Physical parameters
    /// Range pixel spacing in meters
    pub range_pixel_spacing: f64,
    /// Azimuth pixel spacing in meters
    pub azimuth_pixel_spacing: f64,
    /// Slant range time in seconds
    pub slant_range_time: f64,
    /// Burst duration in seconds
    pub burst_duration: f64,
    /// Near-edge slant range in meters
    pub near_range_m: f64,

    /// Instrument Pulse Repetition Frequency from annotation (Hz)
    /// NOTE: For TOPS mode, this is the radar sampling PRF (1451-1717 Hz for S1 IW),
    /// NOT the azimuth line rate of the SLC grid (~486 Hz, stored in azimuth_time_interval)
    pub prf_hz: Option<f64>,

    /// Doppler centroid polynomial coefficients [c0, c1, c2, ...]
    /// CRITICAL: Must be populated by DC-aware deburst for TOPSAR merge
    pub dc_polynomial: Option<Vec<f64>>,
    /// Azimuth time interval (seconds per line) for timing-based alignment
    pub azimuth_time_interval: Option<f64>,
    /// Absolute polynomial reference time (seconds since epoch)
    pub dc_polynomial_t0: Option<f64>,
    /// Azimuth FM-rate estimates sampled along the burst timeline
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub fm_rate_estimates: Option<Vec<FmRateEstimate>>,
}

/// Burst timing record for TOPS processing
///
/// # Bound Semantics
/// All bounds use **exclusive end** convention: `[first, last)` means `first..last` in Rust.
/// - `first_line_global`: inclusive start line
/// - `last_line_global`: **exclusive** end line (one past the last valid line)
/// - `first_valid_sample`: inclusive start sample  
/// - `last_valid_sample`: **exclusive** end sample (one past the last valid sample)
///
/// # Deserialization Compatibility
/// This struct is designed to deserialize from both `BurstTiming` and `BurstRecord` JSON.
/// The `serde(default)` attributes ensure compatibility when deserializing from `BurstRecord`
/// which may have optional or extra fields. The `flatten` annotation captures any additional
/// fields from `BurstRecord` (like `start_sample_global`, `azimuth_time_absolute`) that are
/// not needed for timing calculations but should not cause deserialization failures.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BurstTiming {
    pub subswath_id: String,
    #[serde(default)]
    pub burst_index: usize,
    /// Azimuth time relative to orbit epoch (seconds).
    /// CRITICAL: This field drives all azimuth-to-time conversions in Range-Doppler geocoding.
    /// If null/missing in source JSON, defaults to 0.0 which will cause WRONG geocoding!
    /// The caller MUST validate this is non-zero for IW TOPSAR products.
    #[serde(default, deserialize_with = "deserialize_option_f64_as_f64")]
    pub azimuth_time_rel_orbit: f64,
    /// First line in global coordinates (inclusive)
    #[serde(default)]
    pub first_line_global: usize,
    /// Last line in global coordinates (exclusive: one past the last valid line)
    #[serde(default)]
    pub last_line_global: usize,
    /// First valid sample (inclusive)
    #[serde(default)]
    pub first_valid_sample: Option<usize>,
    /// Last valid sample (exclusive: one past the last valid sample)
    #[serde(default)]
    pub last_valid_sample: Option<usize>,
}

/// Custom deserializer that handles both `f64` and `Option<f64>` (null) JSON values.
/// Returns the value if present, or 0.0 if null/missing.
fn deserialize_option_f64_as_f64<'de, D>(deserializer: D) -> Result<f64, D::Error>
where
    D: serde::Deserializer<'de>,
{
    use serde::Deserialize;
    let opt: Option<f64> = Option::deserialize(deserializer)?;
    Ok(opt.unwrap_or(0.0))
}

/// Burst metadata cached in SarMetadata for downstream consumers
///
/// # Bound Semantics  
/// All bounds use **exclusive end** convention: `[first, last)` means `first..last` in Rust.
/// This matches the Rust range convention and prevents off-by-one errors.
///
/// **IMPORTANT**: When converting from annotation XML (which uses inclusive bounds),
/// add 1 to the last/end values: `last_exclusive = lastValidSample + 1`
///
/// # Example
/// ```ignore
/// // A burst covering lines 0-99 and samples 100-999:
/// BurstRecord {
///     first_line_global: 0,
///     last_line_global: 100,     // exclusive: lines 0..100
///     start_sample_global: 100,
///     end_sample_global: 1000,   // exclusive: samples 100..1000
///     // ...
/// }
/// ```
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BurstRecord {
    pub subswath_id: String,
    pub burst_index: usize,
    /// First line in global coordinates (inclusive)
    pub first_line_global: usize,
    /// Last line in global coordinates (exclusive: one past the last valid line)
    pub last_line_global: usize,
    /// First sample in global coordinates (inclusive)
    pub start_sample_global: usize,
    /// Last sample in global coordinates (exclusive: one past the last valid sample)
    pub end_sample_global: usize,
    /// First valid sample within burst (inclusive)
    pub first_valid_sample: Option<usize>,
    /// Last valid sample within burst (exclusive: one past the last valid sample)
    pub last_valid_sample: Option<usize>,
    /// First valid line within burst (inclusive)
    pub first_valid_line: Option<usize>,
    /// Last valid line within burst (exclusive: one past the last valid line)
    pub last_valid_line: Option<usize>,
    /// Azimuth time relative to orbit epoch (seconds)
    pub azimuth_time_rel_orbit: Option<f64>,
    /// Absolute azimuth time (Unix seconds)
    pub azimuth_time_absolute: Option<f64>,
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
    pub incidence_angle_near: Option<f64>, // degrees - near-range incidence angle from annotation
    pub incidence_angle_far: Option<f64>, // degrees - far-range incidence angle from annotation
    pub incidence_angle_mid_swath: Option<f64>, // degrees - mid-swath incidence angle for RTC normalization
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
    pub burst_records: Vec<BurstRecord>,
    /// Doppler centroid polynomial (t0 + coefficients)
    #[serde(skip_serializing_if = "Option::is_none")]
    pub doppler_centroid: Option<DopplerCentroidPolynomial>,
}

/// Doppler centroid polynomial extracted from annotation metadata
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DopplerCentroidPolynomial {
    pub t0: f64,
    pub coefficients: Vec<f64>,
}

/// Global strict-mode check to forbid non-scientific fallbacks.
/// **Enabled by default** for scientific accuracy. Set SARDINE_STRICT=0 to disable
/// (not recommended for production/research use).
#[inline]
pub fn strict_mode() -> bool {
    std::env::var("SARDINE_STRICT")
        .ok()
        .and_then(|v| v.parse::<u32>().ok())
        .unwrap_or(1)  // Default to strict mode for scientific accuracy
        != 0
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
///
/// These error types provide structured error handling with actionable messages.
/// Each variant includes context to help diagnose and resolve the issue.
#[derive(Debug, thiserror::Error)]
pub enum SarError {
    #[error("I/O error: {0}")]
    Io(#[from] std::io::Error),

    #[error("Invalid data format: {0}. Check that the input file is a valid Sentinel-1 SAFE or ZIP product.")]
    InvalidFormat(String),

    #[error("Invalid input: {0}. Verify the input path exists and is accessible.")]
    InvalidInput(String),

    #[error(
        "Invalid parameter: {0}. Check the parameter against valid ranges in the documentation."
    )]
    InvalidParameter(String),

    #[error("Parameter error: {0}")]
    ParameterError(String),

    #[error("Processing error: {0}. See RUST_LOG=debug for detailed diagnostics.")]
    Processing(String),

    #[error("Data processing error: {0}")]
    DataProcessingError(String),

    #[error("Metadata error: {0}. Ensure the product annotation XML is well-formed.")]
    Metadata(String),

    #[error("Invalid metadata: {0}. This may indicate a non-standard or corrupted SAFE product.")]
    InvalidMetadata(String),

    #[error("Missing calibration data: {0}. Ensure calibration annotation files are present in the SAFE product.")]
    MissingCalibrationData(String),

    #[error("GDAL error: {0}")]
    Gdal(#[from] gdal::errors::GdalError),

    #[error(
        "XML parsing error: {0}. The annotation XML may be malformed or use an unsupported schema."
    )]
    XmlParsing(String),

    #[error("Missing required parameter: {0}. This parameter must be provided or extracted from annotation.")]
    MissingParameter(String),

    #[error("Missing metadata: {0}. Ensure set env vars (SARDINE_SERDE_ONLY, SARDINE_REQUIRE_SUBSWATHS) before import.")]
    MissingMetadata(String),

    #[error(
        "Numerical computation error: {0}. Check input values for NaN/Inf and ensure valid ranges."
    )]
    NumericalError(String),

    #[error("Not implemented: {0}")]
    NotImplemented(String),

    #[error("Coverage error: {0}. Terrain correction produced insufficient valid pixels. Check DEM bbox and timing metadata.")]
    CoverageError(String),

    #[error("Orbit data error: {0}. Ensure precise orbit files are available or embedded orbit data is sufficient.")]
    OrbitError(String),

    #[error("Geocoding error: {0}. Check native_range_pixel_spacing, burst_timings, and orbit vector count.")]
    GeocodingError(String),

    #[error("Validation error: {0}. STEP-2 diagnostic invariant failed.")]
    Validation(String),
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
    /// Optional satellite attitudes as quaternions [w, x, y, z] for each azimuth line
    pub attitudes: Option<Vec<[f64; 4]>>,
    /// Optional geodetic coordinates (lat, lon, h) in radians/meters per azimuth line
    pub geodetic: Option<Vec<[f64; 3]>>,
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

    /// Get quaternion attitude at a specific azimuth line (if available)
    pub fn get_attitude_at_line(&self, line_idx: usize) -> Option<[f64; 4]> {
        self.attitudes
            .as_ref()
            .and_then(|v| v.get(line_idx).copied())
    }

    /// Get geodetic (lat, lon, height) at a specific azimuth line (if available)
    pub fn get_geodetic_at_line(&self, line_idx: usize) -> Option<[f64; 3]> {
        self.geodetic
            .as_ref()
            .and_then(|v| v.get(line_idx).copied())
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
            gamma0_min: 3.16e-4, // Equivalent to -35 dB: 10^(-35/10) = 3.16e-4 (physical minimum)
            gamma0_max: 10.0, // Equivalent to 10 dB: 10^(10/10) = 10 (physical maximum for backscatter)
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
