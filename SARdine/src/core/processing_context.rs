//! Processing Context and Provenance Tracking
//!
//! This module provides a unified context for SAR processing that carries
//! all critical metadata through the processing chain and enables structured
//! provenance tracking for scientific reproducibility.
//!
//! Design follows recommendations from:
//! ai_collaboration_framework/workflows/sentinel1_metadata_annotation_map.md

use chrono::{DateTime, Utc};
use serde::{Deserialize, Serialize};
use std::collections::HashMap;

use crate::types::SarResult;

/// Unified processing context that carries metadata through the entire pipeline
///
/// This struct consolidates all critical parameters that need to be passed between
/// processing stages, preventing parameter loss and enabling full provenance tracking.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ProcessingContext {
    /// Product identification
    pub product_info: ProductInfo,

    /// Time bases and epochs (critical for coordinate conversions)
    pub timing: TimingContext,

    /// Radar instrument parameters
    pub radar: RadarParameters,

    /// Geometric parameters for RTC and geocoding
    pub geometry: GeometryContext,

    /// Burst information (for TOPS modes)
    #[serde(skip_serializing_if = "Option::is_none")]
    pub bursts: Option<BurstContext>,

    /// Calibration metadata
    pub calibration: CalibrationContext,

    /// DEM metadata
    #[serde(skip_serializing_if = "Option::is_none")]
    pub dem: Option<DemContext>,

    /// Processing stages completed
    pub processing_stages: Vec<ProcessingStage>,

    /// Validation results
    pub validation_results: Vec<ProvenanceValidationResult>,

    /// Custom metadata (extensible)
    #[serde(skip_serializing_if = "HashMap::is_empty")]
    pub custom_metadata: HashMap<String, String>,
}

/// Product identification information
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ProductInfo {
    pub platform: String,          // S1A/S1B
    pub acquisition_mode: String,  // IW/SM/EW
    pub product_type: String,      // SLC/GRD
    pub polarisation: String,      // VV/VH/HH/HV
    pub swath: Option<String>,     // IW1/IW2/IW3
    pub product_id: String,
    pub start_time: DateTime<Utc>,
    pub stop_time: DateTime<Utc>,
    pub pass_direction: String, // ASCENDING/DESCENDING
}

/// Timing context - critical for coordinate conversions
///
/// Following Section 2 of metadata annotation map:
/// - t_abs = t_orbit_ref_epoch + azimuth_time (solver)
/// - line = (t_abs - t_product_start) / azimuthTimeInterval
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TimingContext {
    /// Orbit reference epoch (seconds since Unix epoch)
    /// This is the time base for orbit state vectors
    pub orbit_ref_epoch: f64,

    /// Product start time (seconds since Unix epoch)
    /// This is productFirstLineUtcTime from annotation
    pub product_start_time_abs: f64,

    /// Product stop time (seconds since Unix epoch)
    pub product_stop_time_abs: f64,

    /// Product duration (seconds)
    pub product_duration: f64,

    /// Azimuth time interval (seconds per line)
    /// CRITICAL: This is from annotation, NOT 1/PRF for TOPS!
    pub azimuth_time_interval: f64,

    /// Range sampling rate (Hz)
    pub range_sampling_rate: f64,

    /// Slant range time to first pixel (seconds)
    pub slant_range_time: f64,

    /// Total azimuth lines in product
    pub total_azimuth_lines: Option<usize>,

    /// Pulse Repetition Frequency (Hz)
    /// Note: For TOPS, azimuth_time_interval ≠ 1/PRF
    pub prf: f64,
}

/// Radar instrument parameters
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RadarParameters {
    /// Carrier frequency (Hz)
    pub carrier_frequency: f64,

    /// Radar wavelength (meters) = c / carrier_frequency
    pub wavelength: f64,

    /// Range pixel spacing (meters)
    pub range_pixel_spacing: f64,

    /// Azimuth pixel spacing (meters)
    pub azimuth_pixel_spacing: f64,

    /// Speed of light (m/s)
    pub speed_of_light: f64,

    /// Chirp bandwidth (Hz)
    #[serde(skip_serializing_if = "Option::is_none")]
    pub chirp_bandwidth: Option<f64>,

    /// Range bandwidth (Hz)
    #[serde(skip_serializing_if = "Option::is_none")]
    pub range_bandwidth: Option<f64>,
}

/// Geometric context for RTC and terrain correction
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GeometryContext {
    /// Azimuth heading angle (radians)
    /// Derived from orbit velocity or annotation
    pub azimuth_heading_rad: f64,

    /// Ellipsoid incidence angle (radians)
    /// From geolocation grid, used as θᵢ in RTC
    pub ellipsoid_incidence_angle_rad: f64,

    /// Look direction (Right/Left)
    pub look_direction: String,

    /// Scene center latitude (degrees)
    #[serde(skip_serializing_if = "Option::is_none")]
    pub scene_center_lat: Option<f64>,

    /// Scene center longitude (degrees)
    #[serde(skip_serializing_if = "Option::is_none")]
    pub scene_center_lon: Option<f64>,

    /// Doppler centroid model
    #[serde(skip_serializing_if = "Option::is_none")]
    pub doppler_centroid: Option<DopplerCentroidInfo>,
}

/// Doppler centroid information
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DopplerCentroidInfo {
    pub t0: f64,
    pub coeffs: Vec<f64>,
    pub reference_epoch: String, // Time base description
}

/// Burst context for TOPS modes
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BurstContext {
    pub num_bursts: usize,
    pub bursts: Vec<BurstMetadata>,
    pub overlap_windows: Vec<OverlapWindow>,
}

/// Metadata for a single burst
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BurstMetadata {
    pub burst_id: usize,
    pub start_line: usize,
    pub end_line: usize,
    pub sensing_start: DateTime<Utc>,
    pub sensing_stop: DateTime<Utc>,
    pub azimuth_fm_rate: f64,
    pub azimuth_steering_rate: f64,
    pub doppler_centroid: f64,
}

/// Overlap window between bursts
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct OverlapWindow {
    pub burst1_id: usize,
    pub burst2_id: usize,
    pub overlap_lines: usize,
    pub power_preservation_ratio: Option<f64>,
}

/// Calibration context
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CalibrationContext {
    /// Calibration type applied (beta0/sigma0/gamma0)
    pub calibration_type: String,

    /// Calibration constant
    pub calibration_constant: f64,

    /// LUT identifiers
    pub lut_ids: Vec<String>,

    /// Noise removal applied
    pub noise_removal_applied: bool,

    /// NESZ statistics (if computed)
    #[serde(skip_serializing_if = "Option::is_none")]
    pub nesz_stats: Option<NeszStatistics>,

    /// Output radiometry (linear/dB)
    pub output_units: String,
}

/// NESZ (Noise Equivalent Sigma Zero) statistics
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct NeszStatistics {
    pub mean_db: f64,
    pub std_db: f64,
    pub min_db: f64,
    pub max_db: f64,
}

/// DEM context
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DemContext {
    pub dem_source: String,
    pub dem_crs: String,
    pub dem_crs_epsg: u32,
    pub dem_spacing_meters: f32,
    pub vertical_datum: String, // e.g., "EGM96", "WGS84 ellipsoid"
    pub elevation_stats: ElevationStatistics,
    pub dem_bounds: DemBounds,
}

/// DEM elevation statistics
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ElevationStatistics {
    pub min_m: f32,
    pub max_m: f32,
    pub mean_m: f32,
    pub std_m: f32,
}

/// DEM geographic bounds
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DemBounds {
    pub min_lat: f64,
    pub max_lat: f64,
    pub min_lon: f64,
    pub max_lon: f64,
}

/// Processing stage record
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ProcessingStage {
    pub stage_name: String,
    pub timestamp: DateTime<Utc>,
    pub duration_seconds: f64,
    pub parameters: HashMap<String, String>,
    pub status: StageStatus,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub error_message: Option<String>,
}

/// Stage completion status
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub enum StageStatus {
    Success,
    Warning,
    Failed,
    Skipped,
}

/// Validation result record for provenance
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ProvenanceValidationResult {
    pub gate_name: String,
    pub timestamp: DateTime<Utc>,
    pub passed: bool,
    pub metric_value: Option<f64>,
    pub threshold: Option<f64>,
    pub message: String,
}

impl ProcessingContext {
    /// Create a new processing context from annotation and orbit data
    pub fn new(
        product_info: ProductInfo,
        timing: TimingContext,
        radar: RadarParameters,
        geometry: GeometryContext,
    ) -> Self {
        Self {
            product_info,
            timing,
            radar,
            geometry,
            bursts: None,
            calibration: CalibrationContext {
                calibration_type: "none".to_string(),
                calibration_constant: 1.0,
                lut_ids: Vec::new(),
                noise_removal_applied: false,
                nesz_stats: None,
                output_units: "linear".to_string(),
            },
            dem: None,
            processing_stages: Vec::new(),
            validation_results: Vec::new(),
            custom_metadata: HashMap::new(),
        }
    }

    /// Record a processing stage
    pub fn record_stage(
        &mut self,
        stage_name: String,
        duration_seconds: f64,
        parameters: HashMap<String, String>,
        status: StageStatus,
        error_message: Option<String>,
    ) {
        self.processing_stages.push(ProcessingStage {
            stage_name,
            timestamp: Utc::now(),
            duration_seconds,
            parameters,
            status,
            error_message,
        });
    }

    /// Record a validation result
    pub fn record_validation(
        &mut self,
        gate_name: String,
        passed: bool,
        metric_value: Option<f64>,
        threshold: Option<f64>,
        message: String,
    ) {
        self.validation_results.push(ProvenanceValidationResult {
            gate_name,
            timestamp: Utc::now(),
            passed,
            metric_value,
            threshold,
            message,
        });
    }

    /// Check if all validation gates passed
    pub fn all_validations_passed(&self) -> bool {
        self.validation_results.iter().all(|v| v.passed)
    }

    /// Get summary of processing stages
    pub fn stage_summary(&self) -> HashMap<String, StageStatus> {
        self.processing_stages
            .iter()
            .map(|s| (s.stage_name.clone(), s.status.clone()))
            .collect()
    }

    /// Export context to JSON string
    pub fn to_json(&self) -> SarResult<String> {
        serde_json::to_string_pretty(self)
            .map_err(|e| crate::types::SarError::Processing(format!("JSON serialization failed: {}", e)))
    }

    /// Export context to YAML string
    pub fn to_yaml(&self) -> SarResult<String> {
        serde_yaml::to_string(self)
            .map_err(|e| crate::types::SarError::Processing(format!("YAML serialization failed: {}", e)))
    }

    /// Save context to JSON file
    pub fn save_json(&self, path: &std::path::Path) -> SarResult<()> {
        let json = self.to_json()?;
        std::fs::write(path, json)
            .map_err(|e| crate::types::SarError::Processing(format!("Failed to write JSON: {}", e)))?;
        log::info!("✅ Saved processing context to: {}", path.display());
        Ok(())
    }

    /// Save context to YAML file
    pub fn save_yaml(&self, path: &std::path::Path) -> SarResult<()> {
        let yaml = self.to_yaml()?;
        std::fs::write(path, yaml)
            .map_err(|e| crate::types::SarError::Processing(format!("Failed to write YAML: {}", e)))?;
        log::info!("✅ Saved processing context to: {}", path.display());
        Ok(())
    }

    /// Load context from JSON file
    pub fn from_json_file(path: &std::path::Path) -> SarResult<Self> {
        let json = std::fs::read_to_string(path)
            .map_err(|e| crate::types::SarError::Processing(format!("Failed to read JSON: {}", e)))?;
        serde_json::from_str(&json)
            .map_err(|e| crate::types::SarError::Processing(format!("JSON deserialization failed: {}", e)))
    }

    /// Load context from YAML file
    pub fn from_yaml_file(path: &std::path::Path) -> SarResult<Self> {
        let yaml = std::fs::read_to_string(path)
            .map_err(|e| crate::types::SarError::Processing(format!("Failed to read YAML: {}", e)))?;
        serde_yaml::from_str(&yaml)
            .map_err(|e| crate::types::SarError::Processing(format!("YAML deserialization failed: {}", e)))
    }
}

/// Builder for ProcessingContext to simplify construction
pub struct ProcessingContextBuilder {
    product_info: Option<ProductInfo>,
    timing: Option<TimingContext>,
    radar: Option<RadarParameters>,
    geometry: Option<GeometryContext>,
}

impl ProcessingContextBuilder {
    pub fn new() -> Self {
        Self {
            product_info: None,
            timing: None,
            radar: None,
            geometry: None,
        }
    }

    pub fn product_info(mut self, info: ProductInfo) -> Self {
        self.product_info = Some(info);
        self
    }

    pub fn timing(mut self, timing: TimingContext) -> Self {
        self.timing = Some(timing);
        self
    }

    pub fn radar(mut self, radar: RadarParameters) -> Self {
        self.radar = Some(radar);
        self
    }

    pub fn geometry(mut self, geometry: GeometryContext) -> Self {
        self.geometry = Some(geometry);
        self
    }

    pub fn build(self) -> SarResult<ProcessingContext> {
        Ok(ProcessingContext::new(
            self.product_info
                .ok_or_else(|| crate::types::SarError::Processing("Missing product_info".to_string()))?,
            self.timing
                .ok_or_else(|| crate::types::SarError::Processing("Missing timing".to_string()))?,
            self.radar
                .ok_or_else(|| crate::types::SarError::Processing("Missing radar parameters".to_string()))?,
            self.geometry
                .ok_or_else(|| crate::types::SarError::Processing("Missing geometry".to_string()))?,
        ))
    }
}

impl Default for ProcessingContextBuilder {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_processing_context_creation() {
        let ctx = ProcessingContext::new(
            ProductInfo {
                platform: "S1A".to_string(),
                acquisition_mode: "IW".to_string(),
                product_type: "SLC".to_string(),
                polarisation: "VV".to_string(),
                swath: Some("IW1".to_string()),
                product_id: "test".to_string(),
                start_time: Utc::now(),
                stop_time: Utc::now(),
                pass_direction: "ASCENDING".to_string(),
            },
            TimingContext {
                orbit_ref_epoch: 1000000.0,
                product_start_time_abs: 1000010.0,
                product_stop_time_abs: 1000040.0,
                product_duration: 30.0,
                azimuth_time_interval: 0.00334,
                range_sampling_rate: 64e6,
                slant_range_time: 0.005,
                total_azimuth_lines: Some(10000),
                prf: 1700.0,
            },
            RadarParameters {
                carrier_frequency: 5.405e9,
                wavelength: 0.0555,
                range_pixel_spacing: 2.33,
                azimuth_pixel_spacing: 13.9,
                speed_of_light: 299792458.0,
                chirp_bandwidth: None,
                range_bandwidth: None,
            },
            GeometryContext {
                azimuth_heading_rad: 0.0,
                ellipsoid_incidence_angle_rad: 0.61,
                look_direction: "Right".to_string(),
                scene_center_lat: Some(45.0),
                scene_center_lon: Some(10.0),
                doppler_centroid: None,
            },
        );

        assert_eq!(ctx.product_info.platform, "S1A");
        assert_eq!(ctx.timing.prf, 1700.0);
    }

    #[test]
    fn test_stage_recording() {
        let mut ctx = ProcessingContext::new(
            ProductInfo {
                platform: "S1A".to_string(),
                acquisition_mode: "IW".to_string(),
                product_type: "SLC".to_string(),
                polarisation: "VV".to_string(),
                swath: Some("IW1".to_string()),
                product_id: "test".to_string(),
                start_time: Utc::now(),
                stop_time: Utc::now(),
                pass_direction: "ASCENDING".to_string(),
            },
            TimingContext {
                orbit_ref_epoch: 1000000.0,
                product_start_time_abs: 1000010.0,
                product_stop_time_abs: 1000040.0,
                product_duration: 30.0,
                azimuth_time_interval: 0.00334,
                range_sampling_rate: 64e6,
                slant_range_time: 0.005,
                total_azimuth_lines: Some(10000),
                prf: 1700.0,
            },
            RadarParameters {
                carrier_frequency: 5.405e9,
                wavelength: 0.0555,
                range_pixel_spacing: 2.33,
                azimuth_pixel_spacing: 13.9,
                speed_of_light: 299792458.0,
                chirp_bandwidth: None,
                range_bandwidth: None,
            },
            GeometryContext {
                azimuth_heading_rad: 0.0,
                ellipsoid_incidence_angle_rad: 0.61,
                look_direction: "Right".to_string(),
                scene_center_lat: Some(45.0),
                scene_center_lon: Some(10.0),
                doppler_centroid: None,
            },
        );

        ctx.record_stage(
            "calibration".to_string(),
            1.5,
            HashMap::new(),
            StageStatus::Success,
            None,
        );

        assert_eq!(ctx.processing_stages.len(), 1);
        assert_eq!(ctx.processing_stages[0].stage_name, "calibration");
    }

    #[test]
    fn test_json_serialization() {
        let ctx = ProcessingContext::new(
            ProductInfo {
                platform: "S1A".to_string(),
                acquisition_mode: "IW".to_string(),
                product_type: "SLC".to_string(),
                polarisation: "VV".to_string(),
                swath: Some("IW1".to_string()),
                product_id: "test".to_string(),
                start_time: Utc::now(),
                stop_time: Utc::now(),
                pass_direction: "ASCENDING".to_string(),
            },
            TimingContext {
                orbit_ref_epoch: 1000000.0,
                product_start_time_abs: 1000010.0,
                product_stop_time_abs: 1000040.0,
                product_duration: 30.0,
                azimuth_time_interval: 0.00334,
                range_sampling_rate: 64e6,
                slant_range_time: 0.005,
                total_azimuth_lines: Some(10000),
                prf: 1700.0,
            },
            RadarParameters {
                carrier_frequency: 5.405e9,
                wavelength: 0.0555,
                range_pixel_spacing: 2.33,
                azimuth_pixel_spacing: 13.9,
                speed_of_light: 299792458.0,
                chirp_bandwidth: None,
                range_bandwidth: None,
            },
            GeometryContext {
                azimuth_heading_rad: 0.0,
                ellipsoid_incidence_angle_rad: 0.61,
                look_direction: "Right".to_string(),
                scene_center_lat: Some(45.0),
                scene_center_lon: Some(10.0),
                doppler_centroid: None,
            },
        );

        let json = ctx.to_json().unwrap();
        assert!(json.contains("S1A"));
        assert!(json.contains("IW"));
    }
}
