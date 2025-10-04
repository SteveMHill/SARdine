//! # SARdine Processing Pipeline with Validation Gates
//!
//! Comprehensive processing pipeline demonstrating integration of validation gates
//! at critical processing stages to ensure robust quality control.

use crate::core::validation_gates::{ValidationGates, ValidationResult, ValidationSummary};
use crate::types::{SarError, SarResult};
use ndarray::Array2;
use log::{info, warn};
use std::time::Instant;

/// Processing pipeline with integrated validation gates
pub struct ValidatedProcessingPipeline {
    /// Validation configuration
    pub validation_gates: ValidationGates,
    
    /// Pipeline configuration
    pub fail_fast: bool,           // Stop on first validation failure
    pub log_detailed_results: bool, // Log detailed validation diagnostics
    pub save_validation_report: bool, // Save validation results to file
    
    /// Processing state tracking
    validation_results: Vec<ValidationResult>,
    processing_stages: Vec<String>,
    stage_timing: Vec<f64>,
}

impl Default for ValidatedProcessingPipeline {
    fn default() -> Self {
        Self {
            validation_gates: ValidationGates::default(),
            fail_fast: false,  // Allow processing to continue for diagnostic purposes
            log_detailed_results: true,
            save_validation_report: true,
            validation_results: Vec::new(),
            processing_stages: Vec::new(),
            stage_timing: Vec::new(),
        }
    }
}

impl ValidatedProcessingPipeline {
    /// Create a new validated processing pipeline with custom gates
    pub fn with_validation_gates(gates: ValidationGates) -> Self {
        Self {
            validation_gates: gates,
            ..Default::default()
        }
    }
    
    /// Create strict validation pipeline (fail-fast mode)
    pub fn strict() -> Self {
        Self {
            fail_fast: true,
            validation_gates: ValidationGates {
                power_preservation_tolerance: 0.005,  // Stricter ±0.5%
                uncovered_pixel_threshold: 0.001,     // Stricter 0.1%
                nesz_tolerance_db: 1.0,               // Stricter ±1 dB
                geometry_ce90_threshold: 5.0,         // Stricter 5m CE90
                ..ValidationGates::default()
            },
            ..Default::default()
        }
    }
    
    /// Process SAR data with comprehensive validation at each stage
    pub fn process_sar_with_validation(
        &mut self,
        slc_data: &Array2<num_complex::Complex<f32>>,
        calibration_lut: &Array2<f32>,
        dem_data: &Array2<f32>,
        processing_params: &ProcessingParameters,
    ) -> SarResult<ProcessedSarData> {
        info!("🚀 Starting validated SAR processing pipeline");
        let pipeline_start = Instant::now();
        
        // Stage 1: Pre-processing validation
        self.validate_input_data(slc_data, calibration_lut, dem_data)?;
        
        // Stage 2: Deburst with power preservation validation
        let (debursted_data, deburst_metrics) = self.deburst_with_validation(
            slc_data, &processing_params.deburst_params
        )?;
        
        // Stage 3: Calibration with radiometric validation
        let calibrated_data = self.calibrate_with_validation(
            &debursted_data, calibration_lut, &processing_params.calibration_params
        )?;
        
        // Stage 4: Terrain correction with geometry validation
        let terrain_corrected = self.terrain_correct_with_validation(
            &calibrated_data, dem_data, &processing_params.terrain_params
        )?;
        
        // Stage 5: Final quality assessment
        let final_metrics = self.final_quality_assessment(&terrain_corrected)?;
        
        // Generate comprehensive validation report
        let validation_summary = self.generate_validation_report()?;
        
        let total_time = pipeline_start.elapsed().as_secs_f64();
        info!("✅ Validated processing completed in {:.2}s", total_time);
        
        if !validation_summary.overall_success {
            if self.fail_fast {
                return Err(SarError::Processing(
                    "Pipeline validation failed - check validation report for details".to_string()
                ));
            } else {
                warn!("⚠️  Pipeline completed with validation warnings");
            }
        }
        
        Ok(ProcessedSarData {
            terrain_corrected_data: terrain_corrected,
            deburst_metrics,
            final_metrics,
            validation_summary,
            processing_time_seconds: total_time,
        })
    }
    
    /// Validate input data quality and consistency
    fn validate_input_data(
        &mut self,
        slc_data: &Array2<num_complex::Complex<f32>>,
        calibration_lut: &Array2<f32>,
        dem_data: &Array2<f32>,
    ) -> SarResult<()> {
        let stage_start = Instant::now();
        info!("🔍 Stage 1: Input data validation");
        
        // Check for NaN/infinite values in input data
        let slc_nan_count = slc_data.iter()
            .filter(|c| c.re.is_nan() || c.im.is_nan() || c.re.is_infinite() || c.im.is_infinite())
            .count();
        
        let cal_nan_count = calibration_lut.iter()
            .filter(|&&v| v.is_nan() || v.is_infinite())
            .count();
        
        let dem_nan_count = dem_data.iter()
            .filter(|&&v| v.is_nan() || v.is_infinite())
            .count();
        
        // Validate data consistency
        if slc_data.is_empty() {
            return Err(SarError::InvalidInput("Empty SLC data array".to_string()));
        }
        
        if slc_nan_count > 0 {
            warn!("⚠️  Found {} invalid values in SLC data", slc_nan_count);
        }
        
        if cal_nan_count > 0 {
            warn!("⚠️  Found {} invalid values in calibration LUT", cal_nan_count);
        }
        
        if dem_nan_count > 0 {
            warn!("⚠️  Found {} invalid values in DEM data", dem_nan_count);
        }
        
        // LUT domain validation - check calibration LUT bounds
        let line_indices: Vec<usize> = (0..slc_data.nrows()).collect();
        let pixel_indices: Vec<usize> = (0..slc_data.ncols()).collect();
        
        let lut_result = self.validation_gates.validate_lut_domain_bounds(
            &line_indices,
            &pixel_indices,
            calibration_lut.nrows(),
            calibration_lut.ncols(),
            "calibration_lut"
        )?;
        
        self.validation_results.push(lut_result);
        
        let stage_time = stage_start.elapsed().as_secs_f64();
        self.processing_stages.push("Input Validation".to_string());
        self.stage_timing.push(stage_time);
        
        info!("✅ Input validation completed in {:.2}s", stage_time);
        Ok(())
    }
    
    /// Deburst processing with power preservation validation
    fn deburst_with_validation(
        &mut self,
        slc_data: &Array2<num_complex::Complex<f32>>,
        deburst_params: &DeburstParameters,
    ) -> SarResult<(Array2<num_complex::Complex<f32>>, DeburstMetrics)> {
        let stage_start = Instant::now();
        info!("🔧 Stage 2: Deburst with power preservation validation");
        
        // Calculate pre-deburst power in overlap regions
        let pre_deburst_power = slc_data.map(|c| c.norm_sqr());
        
        // Perform deburst processing (simplified - actual implementation would be more complex)
        let debursted_data = self.simulate_deburst_processing(slc_data, deburst_params)?;
        
        // Calculate post-deburst power
        let post_deburst_power = debursted_data.map(|c| c.norm_sqr());
        
        // Define overlap regions based on burst boundaries
        let overlap_regions = self.identify_overlap_regions(slc_data, deburst_params);
        
        // Validate power preservation
        let power_result = self.validation_gates.validate_power_preservation(
            &pre_deburst_power,
            &post_deburst_power,
            &overlap_regions,
        )?;
        
        self.validation_results.push(power_result.clone());
        
        if !power_result.passed && self.fail_fast {
            return Err(SarError::Processing(
                format!("Deburst power preservation failed: {}", power_result.diagnostic_message)
            ));
        }
        
        // Validate uncovered pixels after deburst
        let uncovered_result = self.validation_gates.validate_uncovered_pixels(
            &post_deburst_power,
            None
        )?;
        
        self.validation_results.push(uncovered_result.clone());
        
        if !uncovered_result.passed && self.fail_fast {
            return Err(SarError::Processing(
                format!("Excessive uncovered pixels: {}", uncovered_result.diagnostic_message)
            ));
        }
        
        let stage_time = stage_start.elapsed().as_secs_f64();
        self.processing_stages.push("Deburst".to_string());
        self.stage_timing.push(stage_time);
        
        let metrics = DeburstMetrics {
            power_preservation_error: power_result.measured_value,
            uncovered_pixel_rate: uncovered_result.measured_value,
            processing_time_seconds: stage_time,
        };
        
        info!("✅ Deburst completed in {:.2}s", stage_time);
        Ok((debursted_data, metrics))
    }
    
    /// Calibration with radiometric validation
    fn calibrate_with_validation(
        &mut self,
        debursted_data: &Array2<num_complex::Complex<f32>>,
        calibration_lut: &Array2<f32>,
        calibration_params: &CalibrationParameters,
    ) -> SarResult<Array2<f32>> {
        let stage_start = Instant::now();
        info!("📏 Stage 3: Calibration with radiometric validation");
        
        // Perform radiometric calibration
        let calibrated_sigma0 = self.simulate_calibration(debursted_data, calibration_lut)?;
        
        // CRITICAL SCIENTIFIC FIX: Fail if real NESZ data not available
        let nesz_lut = self.extract_nesz_from_annotation(calibrated_sigma0.ncols())?;
        
        let nesz_result = self.validation_gates.validate_nesz_compliance(
            &calibrated_sigma0,
            10.0, // 10m range pixel spacing
            &nesz_lut,
            None, // No ocean mask - use heuristic
        )?;
        
        self.validation_results.push(nesz_result.clone());
        
        if !nesz_result.passed && self.fail_fast {
            return Err(SarError::Processing(
                format!("NESZ validation failed: {}", nesz_result.diagnostic_message)
            ));
        }
        
        let stage_time = stage_start.elapsed().as_secs_f64();
        self.processing_stages.push("Calibration".to_string());
        self.stage_timing.push(stage_time);
        
        info!("✅ Calibration completed in {:.2}s", stage_time);
        Ok(calibrated_sigma0)
    }
    
    /// Terrain correction with geometry validation
    fn terrain_correct_with_validation(
        &mut self,
        calibrated_data: &Array2<f32>,
        dem_data: &Array2<f32>,
        terrain_params: &TerrainParameters,
    ) -> SarResult<Array2<f32>> {
        let stage_start = Instant::now();
        info!("🗻 Stage 4: Terrain correction with geometry validation");
        
        // Perform terrain correction
        let terrain_corrected = self.simulate_terrain_correction(calibrated_data, dem_data)?;
        
        // CRITICAL SCIENTIFIC FIX: Fail if real ground control points not available
        let (measured_positions, reference_positions) = self.extract_ground_control_points()?;
        
        let geometry_result = self.validation_gates.validate_geometry_accuracy(
            &measured_positions,
            &reference_positions,
            "WGS84",
        )?;
        
        self.validation_results.push(geometry_result.clone());
        
        if !geometry_result.passed && self.fail_fast {
            return Err(SarError::Processing(
                format!("Geometry validation failed: {}", geometry_result.diagnostic_message)
            ));
        }
        
        let stage_time = stage_start.elapsed().as_secs_f64();
        self.processing_stages.push("Terrain Correction".to_string());
        self.stage_timing.push(stage_time);
        
        info!("✅ Terrain correction completed in {:.2}s", stage_time);
        Ok(terrain_corrected)
    }
    
    /// Final quality assessment with Point Target Analysis
    fn final_quality_assessment(
        &mut self,
        final_data: &Array2<f32>,
    ) -> SarResult<QualityMetrics> {
        let stage_start = Instant::now();
        info!("🎯 Stage 5: Final quality assessment");
        
        // Simulate Point Target Analysis metrics
        let pta_metrics = self.perform_point_target_analysis(final_data)?;
        
        let pta_results = self.validation_gates.validate_point_target_analysis(
            pta_metrics.pslr_db,
            pta_metrics.islr_db,
            pta_metrics.irw_meters,
            pta_metrics.radiometric_error_db,
        )?;
        
        for result in pta_results {
            self.validation_results.push(result.clone());
            if !result.passed && self.fail_fast {
                return Err(SarError::Processing(
                    format!("PTA validation failed: {}", result.diagnostic_message)
                ));
            }
        }
        
        let stage_time = stage_start.elapsed().as_secs_f64();
        self.processing_stages.push("Quality Assessment".to_string());
        self.stage_timing.push(stage_time);
        
        info!("✅ Quality assessment completed in {:.2}s", stage_time);
        Ok(pta_metrics)
    }
    
    /// Generate comprehensive validation report
    fn generate_validation_report(&self) -> SarResult<ValidationSummary> {
        let summary = ValidationSummary::from_results(&self.validation_results);
        
        if self.log_detailed_results {
            summary.log_summary();
            
            // Log detailed stage timing
            info!("⏱️  Processing stage timing:");
            for (stage, time) in self.processing_stages.iter().zip(self.stage_timing.iter()) {
                info!("   {}: {:.2}s", stage, time);
            }
            
            let total_time: f64 = self.stage_timing.iter().sum();
            info!("   Total: {:.2}s", total_time);
        }
        
        Ok(summary)
    }
    
    // Helper simulation methods (would be replaced with actual processing in real implementation)
    
    fn simulate_deburst_processing(
        &self,
        slc_data: &Array2<num_complex::Complex<f32>>,
        _params: &DeburstParameters,
    ) -> SarResult<Array2<num_complex::Complex<f32>>> {
        // Simulate deburst with slight power loss to test validation
        Ok(slc_data.map(|c| c * num_complex::Complex::new(0.995, 0.0))) // 0.5% power loss
    }
    
    fn simulate_calibration(
        &self,
        debursted_data: &Array2<num_complex::Complex<f32>>,
        calibration_lut: &Array2<f32>,
    ) -> SarResult<Array2<f32>> {
        let mut sigma0 = Array2::zeros((debursted_data.nrows(), debursted_data.ncols()));
        
        for ((i, j), value) in debursted_data.indexed_iter() {
            let cal_factor = if i < calibration_lut.nrows() && j < calibration_lut.ncols() {
                calibration_lut[[i, j]]
            } else {
                1.0 // Default calibration
            };
            
            sigma0[[i, j]] = value.norm_sqr() * cal_factor;
        }
        
        // Convert to dB
        Ok(sigma0.map(|&v| 10.0 * v.max(1e-10).log10()))
    }
    
    fn simulate_terrain_correction(
        &self,
        calibrated_data: &Array2<f32>,
        _dem_data: &Array2<f32>,
    ) -> SarResult<Array2<f32>> {
        // Simple simulation - actual terrain correction would be much more complex
        Ok(calibrated_data.clone())
    }
    
    fn identify_overlap_regions(
        &self,
        slc_data: &Array2<num_complex::Complex<f32>>,
        _params: &DeburstParameters,
    ) -> Vec<(usize, usize, usize, usize)> {
        // Simulate overlap regions (in real implementation, would come from burst timing)
        let rows = slc_data.nrows();
        let cols = slc_data.ncols();
        
        vec![
            (rows / 4, rows / 2, 0, cols),       // First overlap region
            (rows / 2, 3 * rows / 4, 0, cols),  // Second overlap region
        ]
    }
    
    // CRITICAL SCIENTIFIC FIX: Remove synthetic data generation
    // NESZ must be derived from actual noise annotation files
    fn extract_nesz_from_annotation(&self, _num_pixels: usize) -> SarResult<Vec<f32>> {
        Err(SarError::Processing(
            "NESZ curve must be extracted from real Sentinel-1 noise annotation XML - no synthetic curves allowed for scientific processing".to_string()
        ))
    }
    
    // CRITICAL SCIENTIFIC FIX: Remove synthetic reference point generation
    // Ground control points must come from actual survey data
    fn extract_ground_control_points(&self) -> SarResult<(Vec<(f64, f64)>, Vec<(f64, f64)>)> {
        Err(SarError::Processing(
            "Ground control points must be loaded from external survey data - no synthetic reference points allowed for scientific processing".to_string()
        ))
    }
    
    // CRITICAL SCIENTIFIC FIX: Remove synthetic point target analysis
    // Point target analysis must be performed on real corner reflectors or transponders
    fn perform_point_target_analysis(&self, _data: &Array2<f32>) -> SarResult<QualityMetrics> {
        Err(SarError::Processing(
            "Point target analysis requires actual corner reflectors or transponders - no synthetic metrics allowed for scientific validation".to_string()
        ))
    }
}

// Supporting data structures

#[derive(Debug, Clone)]
pub struct ProcessingParameters {
    pub deburst_params: DeburstParameters,
    pub calibration_params: CalibrationParameters,
    pub terrain_params: TerrainParameters,
}

#[derive(Debug, Clone)]
pub struct DeburstParameters {
    pub azimuth_window: String,
    pub enable_eap_correction: bool,
}

#[derive(Debug, Clone)]
pub struct CalibrationParameters {
    pub calibration_type: String,
    pub apply_antenna_pattern: bool,
}

#[derive(Debug, Clone)]
pub struct TerrainParameters {
    pub dem_resampling: String,
    pub output_spacing: f64,
}

#[derive(Debug)]
pub struct ProcessedSarData {
    pub terrain_corrected_data: Array2<f32>,
    pub deburst_metrics: DeburstMetrics,
    pub final_metrics: QualityMetrics,
    pub validation_summary: ValidationSummary,
    pub processing_time_seconds: f64,
}

#[derive(Debug, Clone)]
pub struct DeburstMetrics {
    pub power_preservation_error: f64,
    pub uncovered_pixel_rate: f64,
    pub processing_time_seconds: f64,
}

#[derive(Debug, Clone)]
pub struct QualityMetrics {
    pub pslr_db: f64,
    pub islr_db: f64,
    pub irw_meters: f64,
    pub radiometric_error_db: f64,
}

#[cfg(test)]
mod tests {
    use super::*;
    use num_complex::Complex;
    
    #[test]
    fn test_validated_processing_pipeline() {
        let mut pipeline = ValidatedProcessingPipeline::default();
        
        // Create test data
        let slc_data = Array2::from_elem((100, 200), Complex::new(1.0, 0.5));
        let calibration_lut = Array2::from_elem((100, 200), 0.8);
        let dem_data = Array2::from_elem((50, 100), 100.0);
        
        let params = ProcessingParameters {
            deburst_params: DeburstParameters {
                azimuth_window: "Hamming".to_string(),
                enable_eap_correction: true,
            },
            calibration_params: CalibrationParameters {
                calibration_type: "sigma0".to_string(),
                apply_antenna_pattern: true,
            },
            terrain_params: TerrainParameters {
                dem_resampling: "bilinear".to_string(),
                output_spacing: 10.0,
            },
        };
        
        let result = pipeline.process_sar_with_validation(
            &slc_data,
            &calibration_lut,
            &dem_data,
            &params
        );
        
        assert!(result.is_ok(), "Processing pipeline should complete successfully");
        
        let processed_data = result.unwrap();
        assert!(processed_data.validation_summary.total_tests > 0, "Should have validation results");
    }
    
    #[test]
    fn test_strict_validation_pipeline() {
        let mut pipeline = ValidatedProcessingPipeline::strict();
        
        // Test that strict mode has tighter tolerances
        assert!(pipeline.validation_gates.power_preservation_tolerance < 0.01);
        assert!(pipeline.validation_gates.uncovered_pixel_threshold < 0.002);
        assert!(pipeline.validation_gates.geometry_ce90_threshold < 10.0);
        assert!(pipeline.fail_fast);
    }
}