//! Core type definitions for Sentinel-1 scene metadata.
//!
//! All scientifically required fields are non-optional.
//! Time domains are explicitly tagged. Units are encoded in field names.
//! Bounds use the exclusive-end convention: \[first, last).

use chrono::{DateTime, Utc};
use serde::{Deserialize, Serialize};

/// Speed of light in vacuum (m/s). Used for wavelength and range derivations.
pub const SPEED_OF_LIGHT_M_S: f64 = 299_792_458.0;

/// Minimum number of orbit state vectors required for interpolation.
/// ESA recommendation: 10 vectors provides ~100s of temporal coverage.
pub const MIN_ORBIT_VECTORS: usize = 10;

// ─── Enumerations ────────────────────────────────────────────────────

/// Sentinel-1 satellite mission.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum Mission {
    S1A,
    S1B,
}

impl std::fmt::Display for Mission {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Mission::S1A => write!(f, "S1A"),
            Mission::S1B => write!(f, "S1B"),
        }
    }
}

/// Sentinel-1 acquisition mode.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum AcquisitionMode {
    /// Interferometric Wide swath (3 sub-swaths, TOPS burst mode).
    IW,
    /// Extra Wide swath (5 sub-swaths, TOPS burst mode).
    EW,
    /// StripMap (single continuous swath).
    SM,
    /// Wave mode (small vignettes).
    WV,
}

impl AcquisitionMode {
    /// Whether this mode uses TOPS burst acquisition.
    pub fn is_tops(&self) -> bool {
        matches!(self, AcquisitionMode::IW | AcquisitionMode::EW)
    }
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

/// Radar polarization channel.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
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

/// Sentinel-1 sub-swath identifier.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
pub enum SubSwathId {
    IW1,
    IW2,
    IW3,
    EW1,
    EW2,
    EW3,
    EW4,
    EW5,
}

impl SubSwathId {
    /// Whether this sub-swath belongs to the given acquisition mode.
    pub fn matches_mode(&self, mode: AcquisitionMode) -> bool {
        match mode {
            AcquisitionMode::IW => matches!(self, SubSwathId::IW1 | SubSwathId::IW2 | SubSwathId::IW3),
            AcquisitionMode::EW => matches!(
                self,
                SubSwathId::EW1 | SubSwathId::EW2 | SubSwathId::EW3 | SubSwathId::EW4 | SubSwathId::EW5
            ),
            _ => false,
        }
    }
}

impl std::fmt::Display for SubSwathId {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let s = match self {
            SubSwathId::IW1 => "IW1",
            SubSwathId::IW2 => "IW2",
            SubSwathId::IW3 => "IW3",
            SubSwathId::EW1 => "EW1",
            SubSwathId::EW2 => "EW2",
            SubSwathId::EW3 => "EW3",
            SubSwathId::EW4 => "EW4",
            SubSwathId::EW5 => "EW5",
        };
        write!(f, "{}", s)
    }
}

// ─── Time Domain ─────────────────────────────────────────────────────

/// Seconds since the orbit reference epoch (earliest orbit state vector time).
///
/// This is the canonical time domain for orbit interpolation and zero-Doppler solving.
/// Using absolute UTC or product-relative times instead causes geolocation errors of ~km.
#[derive(Debug, Clone, Copy, PartialEq, PartialOrd, Serialize, Deserialize)]
pub struct OrbitRelativeSeconds(pub f64);

// ─── Orbit State Vectors ─────────────────────────────────────────────

/// Single orbit state vector: satellite position and velocity at a given time.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct StateVector {
    /// Absolute UTC time of this state vector.
    pub time: DateTime<Utc>,
    /// ECEF position \[x, y, z\] in meters.
    pub position_m: [f64; 3],
    /// ECEF velocity \[vx, vy, vz\] in m/s.
    pub velocity_m_s: [f64; 3],
}

/// Precise orbit data with a defined reference epoch.
///
/// The reference epoch is the time of the earliest state vector.
/// Use [`BurstEntry::azimuth_time_rel`] to convert absolute burst times
/// to orbit-relative seconds on demand.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct OrbitData {
    /// Reference epoch: the earliest state vector time.
    pub reference_epoch: DateTime<Utc>,
    /// State vectors sorted by time, strictly monotonically increasing.
    pub state_vectors: Vec<StateVector>,
}

// ─── Geolocation Grid ──────────────────────────────────────────────

/// A single point from the annotation XML `geolocationGrid`.
///
/// The grid is sparse (typically 10 azimuth × 21 range points per subswath)
/// and must be interpolated to per-pixel values when needed.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GeolocationGridPoint {
    /// Azimuth time of this grid point (UTC).
    pub azimuth_time_utc: DateTime<Utc>,
    /// Two-way slant range time in seconds at this grid point.
    pub slant_range_time_s: f64,
    /// Image line index (azimuth).
    pub line: u32,
    /// Image pixel index (range).
    pub pixel: u32,
    /// WGS84 latitude in degrees.
    pub latitude_deg: f64,
    /// WGS84 longitude in degrees.
    pub longitude_deg: f64,
    /// Height above WGS84 ellipsoid in metres.
    pub height_m: f64,
    /// Incidence angle in degrees.
    pub incidence_angle_deg: f64,
    /// Elevation angle in degrees.
    pub elevation_angle_deg: f64,
}

// ─── Bounding Box ────────────────────────────────────────────────────

/// Geographic bounding box in WGS84 degrees.
///
/// **Limitation:** This uses simple min/max and does NOT handle anti-meridian
/// wrapping. Scenes crossing ±180° longitude (e.g. in the Pacific) will
/// produce invalid `min_lon > max_lon` and fail validation.
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct BoundingBox {
    pub min_lat_deg: f64,
    pub max_lat_deg: f64,
    pub min_lon_deg: f64,
    pub max_lon_deg: f64,
}

// ─── Sub-swath Metadata ──────────────────────────────────────────────

/// Per-subswath metadata extracted from annotation XML.
///
/// All bounds use the exclusive-end convention: \[first, last).
/// This matches Rust range semantics and prevents off-by-one errors.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SubSwathMetadata {
    pub id: SubSwathId,

    /// Number of bursts in this subswath (> 0 for TOPS modes).
    pub burst_count: usize,
    /// Azimuth lines per burst (TOPS burst height).
    pub lines_per_burst: usize,
    /// Range samples in this subswath.
    pub range_samples: usize,
    /// Total azimuth samples (lines) in this subswath.
    pub azimuth_samples: usize,

    // ── Full-image coordinate bounds [first, last) ──

    /// First line in full-image coordinates (inclusive).
    pub first_line: usize,
    /// One past the last line (exclusive end).
    pub last_line: usize,
    /// First sample in full-image coordinates (inclusive).
    pub first_sample: usize,
    /// One past the last sample (exclusive end).
    pub last_sample: usize,

    // ── Physical parameters ──

    /// Range pixel spacing in meters (slant range, before multilook).
    pub range_pixel_spacing_m: f64,
    /// Azimuth pixel spacing in meters (before multilook).
    pub azimuth_pixel_spacing_m: f64,
    /// Two-way slant range time to first pixel in seconds.
    pub slant_range_time_s: f64,
    /// Azimuth time interval: seconds per SLC line (1 / azimuth_frequency).
    /// For TOPS IW this is typically ~0.002s, NOT 1/PRF (~0.0006s).
    pub azimuth_time_interval_s: f64,
    /// Instrument pulse repetition frequency in Hz.
    pub prf_hz: f64,
    /// Burst cycle time in seconds: time between consecutive burst starts.
    /// Derived as the median inter-burst azimuth-time interval.
    /// This is NOT the burst illumination duration (which is shorter due to
    /// the inter-burst beam-steering gap).
    pub burst_cycle_time_s: f64,
}

impl SubSwathMetadata {
    /// Near-edge slant range in meters, derived from `slant_range_time_s`.
    pub fn near_range_m(&self) -> f64 {
        self.slant_range_time_s * SPEED_OF_LIGHT_M_S / 2.0
    }
}

// ─── Burst Entry ─────────────────────────────────────────────────────

/// Single burst timing and bounds from annotation XML.
///
/// All bounds use the exclusive-end convention: \[first, last).
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BurstEntry {
    pub subswath_id: SubSwathId,
    pub burst_index: usize,

    /// Absolute UTC azimuth time of this burst's first line.
    /// Stored as absolute time to avoid epoch-coupling bugs when the
    /// orbit reference epoch changes (e.g. after applying precise orbits).
    pub azimuth_time_utc: DateTime<Utc>,

    // ── Line bounds in full-image coordinates [first, last) ──

    /// First line of this burst (inclusive).
    pub first_line: usize,
    /// One past the last line of this burst (exclusive end).
    pub last_line: usize,

    // ── Valid sample extent [first, last) ──

    /// First valid sample within burst (inclusive).
    pub first_valid_sample: usize,
    /// One past the last valid sample (exclusive end).
    pub last_valid_sample: usize,

    /// Index of the source SAFE slice this burst originates from.
    ///
    /// For scenes from a single `.SAFE` file (the default) this is always 0.
    /// For multi-slice assembled scenes (see [`crate::slice_assembly`]) this
    /// identifies which entry in `AssembledScene::safe_paths` owns this burst.
    /// A value of `k` means the burst's SLC data lives in
    /// `safe_paths[k]/measurement/<swath>-<pol>.tiff`.
    #[serde(default)]
    pub slice_index: usize,
}

impl BurstEntry {
    /// Convert absolute azimuth time to seconds relative to an orbit reference epoch.
    ///
    /// Use this at the point of orbit interpolation, NOT at parse time.
    /// This avoids the epoch-coupling bug where changing the orbit
    /// invalidates pre-computed relative times.
    pub fn azimuth_time_rel(&self, reference_epoch: DateTime<Utc>) -> OrbitRelativeSeconds {
        let delta_us = (self.azimuth_time_utc - reference_epoch)
            .num_microseconds()
            .unwrap_or(0); // SAFETY-OK: chrono microseconds cannot overflow for orbit-window durations
        OrbitRelativeSeconds(delta_us as f64 / 1e6)
    }
}

// ─── Scene Metadata (top-level) ──────────────────────────────────────

/// Validated metadata for a Sentinel-1 SLC scene.
///
/// All scientifically required fields are non-optional.
/// Call [`.validated()`](SceneMetadata::validated) to check all invariants.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SceneMetadata {
    /// Product identifier from the SAFE directory name.
    pub product_id: String,
    /// Satellite mission.
    pub mission: Mission,
    /// Acquisition mode.
    pub acquisition_mode: AcquisitionMode,
    /// Polarization channels present in the product. Non-empty.
    pub polarizations: Vec<Polarization>,

    /// Product first line UTC time.
    pub start_time: DateTime<Utc>,
    /// Product last line UTC time.
    pub stop_time: DateTime<Utc>,

    /// Radar carrier frequency in Hz (C-band, ~5.405 GHz for Sentinel-1).
    pub radar_frequency_hz: f64,
    /// Range sampling rate in Hz.
    pub range_sampling_rate_hz: f64,

    /// Geographic bounding box of the scene footprint.
    pub bounding_box: BoundingBox,
    /// Per-subswath metadata. Non-empty for TOPS modes (IW, EW).
    pub sub_swaths: Vec<SubSwathMetadata>,
    /// Per-burst timing and bounds. Non-empty for TOPS modes (IW, EW).
    pub bursts: Vec<BurstEntry>,
    /// Orbit state vectors and reference epoch.
    pub orbit: OrbitData,
}

impl SceneMetadata {
    /// Radar wavelength in meters, derived from `radar_frequency_hz`.
    pub fn wavelength_m(&self) -> f64 {
        SPEED_OF_LIGHT_M_S / self.radar_frequency_hz
    }

    /// Check all invariants. Returns `Ok(self)` if valid, `Err` with all violations otherwise.
    pub fn validated(self) -> Result<Self, crate::validate::ValidationErrors> {
        let errors = crate::validate::check(&self);
        if errors.is_empty() {
            Ok(self)
        } else {
            Err(crate::validate::ValidationErrors(errors))
        }
    }
}
