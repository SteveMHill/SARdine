/*!
 * Comprehensive Metadata and Provenance Tracking System
 * 
 * This module implements scientifically robust metadata and provenance tracking for SAR processing
 * including:
 * - Algorithm version tracking and processing parameters
 * - Input data provenance (DEM source, orbit data, calibration parameters)
 * - Uncertainty propagation and quality metrics
 * - Processing history and algorithm execution tracking
 * - Scientific reproducibility metadata
 * 
 * Based on:
 * - ISO 19115: Geographic information metadata standards
 * - CEOS Analysis Ready Data for Land (CARD4L) specifications
 * - Scientific reproducibility best practices
 * - ESA Sentinel-1 product definition standards
 */

use crate::types::*;
use crate::core::quality_assessment::QualityConfig;
use serde::{Serialize, Deserialize};
use chrono::{DateTime, Utc};
use std::collections::HashMap;
// Removed unused imports: std::path::PathBuf

/// Comprehensive processing metadata container
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ProcessingMetadata {
    /// Processing identification
    pub processing_id: String,
    pub processing_timestamp: DateTime<Utc>,
    pub processor_version: String,
    pub algorithm_versions: HashMap<String, String>,
    
    /// Input data provenance
    pub input_provenance: InputProvenance,
    
    /// Processing parameters
    pub processing_parameters: ProcessingParameters,
    
    /// Processing history
    pub processing_history: Vec<ProcessingStep>,
    
    /// Quality assessment results
    pub quality_metrics: QualityMetrics,
    
    /// Uncertainty estimates
    pub uncertainty_estimates: UncertaintyEstimates,
    
    /// Output product metadata
    pub output_metadata: OutputProductMetadata,
}

/// Input data provenance tracking
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct InputProvenance {
    /// SAR data provenance
    pub sar_data: SarDataProvenance,
    
    /// DEM data provenance
    pub dem_data: DemDataProvenance,
    
    /// Orbit data provenance
    pub orbit_data: OrbitDataProvenance,
    
    /// Calibration data provenance
    pub calibration_data: CalibrationProvenance,
}

/// SAR input data provenance
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SarDataProvenance {
    pub product_id: String,
    pub mission: String,
    pub platform: String,
    pub sensor: String,
    pub acquisition_mode: String,
    pub polarizations: Vec<String>,
    pub acquisition_start: DateTime<Utc>,
    pub acquisition_stop: DateTime<Utc>,
    pub orbit_number: u32,
    pub relative_orbit_number: u32,
    pub cycle_number: u32,
    pub slice_number: u32,
    pub product_type: String,
    pub processing_level: String,
    pub file_path: String,
    pub file_size_bytes: u64,
    pub file_checksum: Option<String>,
}

/// DEM data provenance
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DemDataProvenance {
    pub dem_source: String,           // e.g., "SRTM", "Copernicus DEM"
    pub dem_version: String,          // e.g., "SRTM v3.0", "Copernicus DEM 30m"
    pub spatial_resolution: f32,      // meters
    pub vertical_accuracy: f32,       // meters (1-sigma)
    pub horizontal_accuracy: f32,     // meters (1-sigma)
    pub coverage_area: BoundingBox,
    pub download_timestamp: DateTime<Utc>,
    pub tiles_used: Vec<String>,
    pub processing_method: String,    // e.g., "bilinear_interpolation"
    pub void_fill_method: Option<String>,
}

/// Orbit data provenance
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct OrbitDataProvenance {
    pub orbit_file_type: String,     // e.g., "Precise Orbit Ephemerides (POE)"
    pub orbit_file_name: String,
    pub orbit_accuracy: f32,         // meters (3D position accuracy)
    pub validity_start: DateTime<Utc>,
    pub validity_stop: DateTime<Utc>,
    pub creation_date: DateTime<Utc>,
    pub download_source: String,
    pub interpolation_method: String, // e.g., "Lagrange polynomial"
    pub polynomial_degree: u8,
    pub time_sampling: f32,          // seconds
}

/// Calibration data provenance
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CalibrationProvenance {
    pub calibration_source: String,  // e.g., "Sentinel-1 annotation XML"
    pub calibration_type: String,    // e.g., "sigma0", "gamma0", "beta0"
    pub calibration_accuracy: f32,   // dB (1-sigma radiometric accuracy)
    pub antenna_pattern_source: String,
    pub noise_calibration_source: String,
    pub range_spreading_loss_model: String,
    pub calibration_constant_uncertainty: f32, // dB
}

/// Processing parameters configuration
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ProcessingParameters {
    /// Core processing parameters
    pub multilook_factors: (u32, u32), // (range, azimuth)
    pub output_resolution: f32,         // meters
    pub output_projection: String,      // e.g., "WGS84 / UTM zone 33N"
    pub resampling_method: String,      // e.g., "bilinear"
    
    /// Quality masking parameters
    pub quality_config: QualityConfig,
    
    /// Terrain correction parameters
    pub terrain_correction_config: TerrainCorrectionConfig,
    
    /// Speckle filtering parameters
    pub speckle_filter_config: Option<SpeckleFilterConfig>,
    
    /// Custom parameters
    pub custom_parameters: HashMap<String, serde_json::Value>,
}

/// Individual processing step record
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ProcessingStep {
    pub step_number: u32,
    pub step_name: String,
    pub algorithm_name: String,
    pub algorithm_version: String,
    pub start_time: DateTime<Utc>,
    pub end_time: DateTime<Utc>,
    pub processing_duration: f64,    // seconds
    pub input_description: String,
    pub output_description: String,
    pub parameters_used: HashMap<String, serde_json::Value>,
    pub quality_flags: Vec<String>,
    pub warnings: Vec<String>,
    pub algorithm_status: String,    // "success", "fallback", "warning", "error"
}

/// Quality metrics from processing
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct QualityMetrics {
    /// Overall quality assessment
    pub overall_quality_score: f32,  // 0-1
    pub valid_pixel_percentage: f32,
    pub data_coverage_percentage: f32,
    
    /// SNR metrics
    pub mean_snr_db: f32,
    pub snr_histogram: Vec<u32>,     // 10 bins
    pub pixels_below_snr_threshold: u32,
    
    /// Geometric quality
    pub mean_local_incidence_angle: f32,
    pub layover_pixel_count: u32,
    pub shadow_pixel_count: u32,
    pub foreshortening_severity: f32,
    
    /// Radiometric quality
    pub radiometric_accuracy_estimate: f32, // dB
    pub dynamic_range_db: f32,
    pub saturation_pixel_count: u32,
    pub noise_pixel_count: u32,
    
    /// Statistical quality
    pub outlier_pixel_count: u32,
    pub texture_uniformity_index: f32,
    pub statistical_homogeneity: f32,
}

/// Uncertainty propagation estimates
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct UncertaintyEstimates {
    /// Radiometric uncertainties (dB, 1-sigma)
    pub calibration_uncertainty: f32,
    pub speckle_uncertainty: f32,
    pub atmospheric_uncertainty: f32,
    pub terrain_correction_uncertainty: f32,
    pub total_radiometric_uncertainty: f32,
    
    /// Geometric uncertainties (meters, 1-sigma)
    pub dem_uncertainty: f32,
    pub orbit_uncertainty: f32,
    pub coregistration_uncertainty: f32,
    pub total_geometric_uncertainty: f32,
    
    /// Quality score uncertainties
    pub snr_uncertainty: f32,
    pub quality_assessment_confidence: f32,
}

/// Output product metadata
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct OutputProductMetadata {
    pub product_name: String,
    pub product_version: String,
    pub creation_time: DateTime<Utc>,
    pub file_format: String,
    pub coordinate_system: String,
    pub spatial_resolution: f32,
    pub pixel_spacing: (f32, f32),   // (x, y) in meters
    pub image_dimensions: (u32, u32), // (width, height)
    pub data_type: String,
    pub units: String,
    pub nodata_value: Option<f32>,
    pub valid_range: (f32, f32),
    pub file_size_bytes: u64,
}

/// Terrain correction configuration for metadata
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TerrainCorrectionConfig {
    pub correction_type: String,     // "radiometric", "geometric", "both"
    pub reference_incidence_angle: f32,
    pub minimum_incidence_angle: f32,
    pub maximum_incidence_angle: f32,
    pub dem_interpolation_method: String,
    pub geocoding_algorithm: String,  // e.g., "Range-Doppler"
}

/// Speckle filter configuration for metadata
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SpeckleFilterConfig {
    pub filter_type: String,         // e.g., "Lee", "Enhanced Lee", "Gamma MAP"
    pub window_size: u32,
    pub number_of_looks: f32,
    pub preserve_edges: bool,
    pub damping_factor: Option<f32>,
}

/// Metadata manager implementation
pub struct MetadataManager {
    metadata: ProcessingMetadata,
    step_counter: u32,
}

impl MetadataManager {
    /// Create new metadata manager
    pub fn new(processing_id: String) -> Self {
        let metadata = ProcessingMetadata {
            processing_id,
            processing_timestamp: Utc::now(),
            processor_version: env!("CARGO_PKG_VERSION").to_string(),
            algorithm_versions: Self::get_algorithm_versions(),
            input_provenance: InputProvenance::default(),
            processing_parameters: ProcessingParameters::default(),
            processing_history: Vec::new(),
            quality_metrics: QualityMetrics::default(),
            uncertainty_estimates: UncertaintyEstimates::default(),
            output_metadata: OutputProductMetadata::default(),
        };
        
        Self {
            metadata,
            step_counter: 0,
        }
    }
    
    /// Get algorithm versions
    fn get_algorithm_versions() -> HashMap<String, String> {
        let mut versions = HashMap::new();
        versions.insert("sardine".to_string(), env!("CARGO_PKG_VERSION").to_string());
        versions.insert("deburst".to_string(), "2.1.0".to_string());
        versions.insert("calibration".to_string(), "1.5.0".to_string());
        versions.insert("terrain_correction".to_string(), "3.2.0".to_string());
        versions.insert("quality_assessment".to_string(), "1.0.0".to_string());
        versions.insert("speckle_filter".to_string(), "2.0.0".to_string());
        versions
    }
    
    /// Set SAR data provenance
    pub fn set_sar_provenance(&mut self, provenance: SarDataProvenance) {
        self.metadata.input_provenance.sar_data = provenance;
    }
    
    /// Set DEM data provenance
    pub fn set_dem_provenance(&mut self, provenance: DemDataProvenance) {
        self.metadata.input_provenance.dem_data = provenance;
    }
    
    /// Set orbit data provenance
    pub fn set_orbit_provenance(&mut self, provenance: OrbitDataProvenance) {
        self.metadata.input_provenance.orbit_data = provenance;
    }
    
    /// Set processing parameters
    pub fn set_processing_parameters(&mut self, parameters: ProcessingParameters) {
        self.metadata.processing_parameters = parameters;
    }
    
    /// Record processing step
    pub fn record_processing_step(
        &mut self,
        step_name: String,
        algorithm_name: String,
        algorithm_version: String,
        input_description: String,
        output_description: String,
        parameters: HashMap<String, serde_json::Value>,
        quality_flags: Vec<String>,
        warnings: Vec<String>,
        algorithm_status: String,
        processing_duration: f64,
    ) {
        self.step_counter += 1;
        
        let step = ProcessingStep {
            step_number: self.step_counter,
            step_name,
            algorithm_name,
            algorithm_version,
            start_time: Utc::now(),
            end_time: Utc::now(),
            processing_duration,
            input_description,
            output_description,
            parameters_used: parameters,
            quality_flags,
            warnings,
            algorithm_status,
        };
        
        self.metadata.processing_history.push(step.clone());
        log::info!("Recorded processing step {}: {}", self.step_counter, step.step_name);
    }
    
    /// Set quality metrics
    pub fn set_quality_metrics(&mut self, metrics: QualityMetrics) {
        self.metadata.quality_metrics = metrics;
    }
    
    /// Set uncertainty estimates
    pub fn set_uncertainty_estimates(&mut self, uncertainties: UncertaintyEstimates) {
        self.metadata.uncertainty_estimates = uncertainties;
    }
    
    /// Set output product metadata
    pub fn set_output_metadata(&mut self, output_meta: OutputProductMetadata) {
        self.metadata.output_metadata = output_meta;
    }
    
    /// Export metadata as JSON
    pub fn export_json(&self) -> SarResult<String> {
        serde_json::to_string_pretty(&self.metadata)
            .map_err(|e| SarError::Processing(format!("Metadata serialization failed: {}", e)))
    }
    
    /// Export metadata as XML (ISO 19115 compliant)
    pub fn export_xml(&self) -> SarResult<String> {
        // Simplified XML export - could be enhanced with full ISO 19115 compliance
        let mut xml = String::from(r#"<?xml version="1.0" encoding="UTF-8"?>
<metadata xmlns:gmd="http://www.isotc211.org/2005/gmd" xmlns:gco="http://www.isotc211.org/2005/gco">
  <processing_metadata>
"#);
        
        xml.push_str(&format!("    <processing_id>{}</processing_id>\n", self.metadata.processing_id));
        xml.push_str(&format!("    <processing_timestamp>{}</processing_timestamp>\n", 
                              self.metadata.processing_timestamp.format("%Y-%m-%dT%H:%M:%S%.3fZ")));
        xml.push_str(&format!("    <processor_version>{}</processor_version>\n", self.metadata.processor_version));
        
        // Input provenance
        xml.push_str("    <input_provenance>\n");
        xml.push_str(&format!("      <sar_product_id>{}</sar_product_id>\n", self.metadata.input_provenance.sar_data.product_id));
        xml.push_str(&format!("      <dem_source>{}</dem_source>\n", self.metadata.input_provenance.dem_data.dem_source));
        xml.push_str(&format!("      <orbit_file_type>{}</orbit_file_type>\n", self.metadata.input_provenance.orbit_data.orbit_file_type));
        xml.push_str("    </input_provenance>\n");
        
        // Quality metrics
        xml.push_str("    <quality_metrics>\n");
        xml.push_str(&format!("      <overall_quality_score>{:.3}</overall_quality_score>\n", 
                              self.metadata.quality_metrics.overall_quality_score));
        xml.push_str(&format!("      <valid_pixel_percentage>{:.1}</valid_pixel_percentage>\n", 
                              self.metadata.quality_metrics.valid_pixel_percentage));
        xml.push_str(&format!("      <mean_snr_db>{:.1}</mean_snr_db>\n", 
                              self.metadata.quality_metrics.mean_snr_db));
        xml.push_str("    </quality_metrics>\n");
        
        xml.push_str("  </processing_metadata>\n</metadata>");
        
        Ok(xml)
    }
    
    /// Save metadata to file
    pub fn save_to_file(&self, output_path: &str, format: &str) -> SarResult<()> {
        let content = match format.to_lowercase().as_str() {
            "json" => self.export_json()?,
            "xml" => self.export_xml()?,
            _ => return Err(SarError::Processing(format!("Unsupported metadata format: {}", format))),
        };
        
        std::fs::write(output_path, content)
            .map_err(|e| SarError::Io(e))?;
        
        log::info!("Metadata saved to: {}", output_path);
        Ok(())
    }
    
    /// Generate processing summary report
    pub fn generate_processing_summary(&self) -> String {
        let mut summary = String::new();
        
        summary.push_str("=== SAR PROCESSING SUMMARY ===\n\n");
        summary.push_str(&format!("Processing ID: {}\n", self.metadata.processing_id));
        summary.push_str(&format!("Timestamp: {}\n", self.metadata.processing_timestamp.format("%Y-%m-%d %H:%M:%S UTC")));
        summary.push_str(&format!("Processor Version: {}\n\n", self.metadata.processor_version));
        
        summary.push_str("INPUT DATA:\n");
        summary.push_str(&format!("  SAR Product: {}\n", self.metadata.input_provenance.sar_data.product_id));
        summary.push_str(&format!("  Mission: {}\n", self.metadata.input_provenance.sar_data.mission));
        summary.push_str(&format!("  DEM Source: {} ({})\n", 
                                  self.metadata.input_provenance.dem_data.dem_source,
                                  self.metadata.input_provenance.dem_data.dem_version));
        summary.push_str(&format!("  Orbit Data: {}\n\n", self.metadata.input_provenance.orbit_data.orbit_file_type));
        
        summary.push_str("PROCESSING STEPS:\n");
        for step in &self.metadata.processing_history {
            summary.push_str(&format!("  {}. {} ({})\n", step.step_number, step.step_name, step.algorithm_status));
            summary.push_str(&format!("     Duration: {:.2}s\n", step.processing_duration));
            if !step.warnings.is_empty() {
                summary.push_str(&format!("     Warnings: {}\n", step.warnings.join(", ")));
            }
        }
        
        summary.push_str("\nQUALITY METRICS:\n");
        summary.push_str(&format!("  Overall Quality Score: {:.3}\n", self.metadata.quality_metrics.overall_quality_score));
        summary.push_str(&format!("  Valid Pixels: {:.1}%\n", self.metadata.quality_metrics.valid_pixel_percentage));
        summary.push_str(&format!("  Mean SNR: {:.1} dB\n", self.metadata.quality_metrics.mean_snr_db));
        summary.push_str(&format!("  Dynamic Range: {:.1} dB\n", self.metadata.quality_metrics.dynamic_range_db));
        
        summary.push_str("\nUNCERTAINTY ESTIMATES:\n");
        summary.push_str(&format!("  Total Radiometric Uncertainty: ±{:.2} dB\n", 
                                  self.metadata.uncertainty_estimates.total_radiometric_uncertainty));
        summary.push_str(&format!("  Total Geometric Uncertainty: ±{:.1} m\n", 
                                  self.metadata.uncertainty_estimates.total_geometric_uncertainty));
        
        summary
    }
    
    /// Get metadata reference
    pub fn get_metadata(&self) -> &ProcessingMetadata {
        &self.metadata
    }
}

// Default implementations for metadata structures
impl Default for InputProvenance {
    fn default() -> Self {
        Self {
            sar_data: SarDataProvenance::default(),
            dem_data: DemDataProvenance::default(),
            orbit_data: OrbitDataProvenance::default(),
            calibration_data: CalibrationProvenance::default(),
        }
    }
}

impl Default for SarDataProvenance {
    fn default() -> Self {
        Self {
            product_id: "Unknown".to_string(),
            mission: "Unknown".to_string(),
            platform: "Unknown".to_string(),
            sensor: "Unknown".to_string(),
            acquisition_mode: "Unknown".to_string(),
            polarizations: vec![],
            acquisition_start: Utc::now(),
            acquisition_stop: Utc::now(),
            orbit_number: 0,
            relative_orbit_number: 0,
            cycle_number: 0,
            slice_number: 0,
            product_type: "Unknown".to_string(),
            processing_level: "Unknown".to_string(),
            file_path: "Unknown".to_string(),
            file_size_bytes: 0,
            file_checksum: None,
        }
    }
}

impl Default for DemDataProvenance {
    fn default() -> Self {
        Self {
            dem_source: "Unknown".to_string(),
            dem_version: "Unknown".to_string(),
            spatial_resolution: 30.0,
            vertical_accuracy: 10.0,
            horizontal_accuracy: 15.0,
            coverage_area: BoundingBox {
                min_lon: 0.0,
                max_lon: 0.0,
                min_lat: 0.0,
                max_lat: 0.0,
            },
            download_timestamp: Utc::now(),
            tiles_used: vec![],
            processing_method: "Unknown".to_string(),
            void_fill_method: None,
        }
    }
}

impl Default for OrbitDataProvenance {
    fn default() -> Self {
        Self {
            orbit_file_type: "Unknown".to_string(),
            orbit_file_name: "Unknown".to_string(),
            orbit_accuracy: 5.0,
            validity_start: Utc::now(),
            validity_stop: Utc::now(),
            creation_date: Utc::now(),
            download_source: "Unknown".to_string(),
            interpolation_method: "Lagrange polynomial".to_string(),
            polynomial_degree: 5,
            time_sampling: 10.0,
        }
    }
}

impl Default for CalibrationProvenance {
    fn default() -> Self {
        Self {
            calibration_source: "Unknown".to_string(),
            calibration_type: "gamma0".to_string(),
            calibration_accuracy: 0.5,
            antenna_pattern_source: "Unknown".to_string(),
            noise_calibration_source: "Unknown".to_string(),
            range_spreading_loss_model: "Unknown".to_string(),
            calibration_constant_uncertainty: 0.3,
        }
    }
}

impl Default for ProcessingParameters {
    fn default() -> Self {
        Self {
            multilook_factors: (4, 1),
            output_resolution: 10.0,
            output_projection: "WGS84".to_string(),
            resampling_method: "bilinear".to_string(),
            quality_config: QualityConfig::default(),
            terrain_correction_config: TerrainCorrectionConfig::default(),
            speckle_filter_config: None,
            custom_parameters: HashMap::new(),
        }
    }
}

impl Default for QualityMetrics {
    fn default() -> Self {
        Self {
            overall_quality_score: 0.0,
            valid_pixel_percentage: 0.0,
            data_coverage_percentage: 0.0,
            mean_snr_db: 0.0,
            snr_histogram: vec![0; 10],
            pixels_below_snr_threshold: 0,
            mean_local_incidence_angle: 0.0,
            layover_pixel_count: 0,
            shadow_pixel_count: 0,
            foreshortening_severity: 0.0,
            radiometric_accuracy_estimate: 0.0,
            dynamic_range_db: 0.0,
            saturation_pixel_count: 0,
            noise_pixel_count: 0,
            outlier_pixel_count: 0,
            texture_uniformity_index: 0.0,
            statistical_homogeneity: 0.0,
        }
    }
}

impl Default for UncertaintyEstimates {
    fn default() -> Self {
        Self {
            calibration_uncertainty: 0.5,
            speckle_uncertainty: 1.0,
            atmospheric_uncertainty: 0.2,
            terrain_correction_uncertainty: 0.3,
            total_radiometric_uncertainty: 1.2,
            dem_uncertainty: 10.0,
            orbit_uncertainty: 5.0,
            coregistration_uncertainty: 2.0,
            total_geometric_uncertainty: 11.2,
            snr_uncertainty: 1.0,
            quality_assessment_confidence: 0.8,
        }
    }
}

impl Default for OutputProductMetadata {
    fn default() -> Self {
        Self {
            product_name: "Unknown".to_string(),
            product_version: "1.0".to_string(),
            creation_time: Utc::now(),
            file_format: "GeoTIFF".to_string(),
            coordinate_system: "WGS84".to_string(),
            spatial_resolution: 10.0,
            pixel_spacing: (10.0, 10.0),
            image_dimensions: (0, 0),
            data_type: "Float32".to_string(),
            units: "dB".to_string(),
            nodata_value: Some(-9999.0),
            valid_range: (-50.0, 20.0),
            file_size_bytes: 0,
        }
    }
}

impl Default for TerrainCorrectionConfig {
    fn default() -> Self {
        Self {
            correction_type: "both".to_string(),
            reference_incidence_angle: 30.0,
            minimum_incidence_angle: 15.0,
            maximum_incidence_angle: 65.0,
            dem_interpolation_method: "bilinear".to_string(),
            geocoding_algorithm: "Range-Doppler".to_string(),
        }
    }
}
