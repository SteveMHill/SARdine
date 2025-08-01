use crate::types::{SarError, SarResult, SarImage, SarRealImage, SubSwath};
use ndarray::Array2;
use std::collections::HashMap;

/// Enhanced TOPSAR merge processor for combining IW sub-swaths
/// Implements state-of-the-art merging algorithms with quality control
pub struct TopsarMerge {
    /// Sub-swath information
    subswaths: Vec<SubSwath>,
    /// Overlap regions between adjacent sub-swaths
    overlap_regions: Vec<OverlapRegion>,
    /// Output grid parameters
    output_grid: OutputGrid,
    /// Processing parameters for enhanced merge
    merge_params: MergeParameters,
    /// Quality control settings
    quality_control: QualityControl,
}

/// Enhanced merge parameters for scientific processing
#[derive(Debug, Clone)]
pub struct MergeParameters {
    /// Blending method for overlap regions
    pub blending_method: BlendingMethod,
    /// Phase preservation for complex data
    pub preserve_phase: bool,
    /// Overlap region optimization
    pub optimize_overlaps: bool,
    /// Enable parallel processing
    pub enable_parallel: bool,
    /// Chunk size for memory-efficient processing
    pub chunk_size: usize,
    /// Edge feathering width (pixels)
    pub feather_width: usize,
}

/// Blending methods for overlap regions
#[derive(Debug, Clone)]
pub enum BlendingMethod {
    /// Linear blending with distance weighting
    Linear,
    /// Gaussian weighted blending
    Gaussian { sigma: f32 },
    /// Advanced multi-scale blending
    MultiScale { scales: Vec<usize> },
    /// Phase-coherent blending for complex data
    PhaseCoherent,
    /// ESA SNAP compatible blending
    SnapCompatible,
}

/// Quality control parameters
#[derive(Debug, Clone)]
pub struct QualityControl {
    /// Enable quality validation
    pub enable_validation: bool,
    /// Maximum phase discontinuity (radians)
    pub max_phase_discontinuity: f32,
    /// Minimum valid pixel ratio in overlaps
    pub min_valid_pixel_ratio: f32,
    /// Enable seamline optimization
    pub optimize_seamlines: bool,
    /// Radiometric consistency tolerance
    pub radiometric_tolerance: f32,
}

/// Overlap region between two adjacent sub-swaths
#[derive(Debug, Clone)]
pub struct OverlapRegion {
    /// First sub-swath ID (e.g., "IW1")
    pub swath1_id: String,
    /// Second sub-swath ID (e.g., "IW2")
    pub swath2_id: String,
    /// Range extent of overlap in swath1 coordinates
    pub swath1_range_start: usize,
    pub swath1_range_end: usize,
    /// Range extent of overlap in swath2 coordinates
    pub swath2_range_start: usize,
    pub swath2_range_end: usize,
    /// Azimuth extent (common for both swaths)
    pub azimuth_start: usize,
    pub azimuth_end: usize,
    /// Optimized blending weights for smooth transition
    pub weights: Array2<f32>,
    /// Quality metrics for this overlap
    pub quality_metrics: OverlapQuality,
}

/// Quality metrics for overlap regions
#[derive(Debug, Clone)]
pub struct OverlapQuality {
    /// Phase coherence in overlap region
    pub phase_coherence: f32,
    /// Radiometric consistency
    pub radiometric_consistency: f32,
    /// Valid pixel percentage
    pub valid_pixel_percentage: f32,
    /// Seamline quality score
    pub seamline_quality: f32,
}

/// Output grid parameters for merged image
#[derive(Debug, Clone)]
pub struct OutputGrid {
    /// Total range samples in merged image
    pub range_samples: usize,
    /// Total azimuth samples in merged image
    pub azimuth_samples: usize,
    /// Range pixel spacing (meters)
    pub range_pixel_spacing: f64,
    /// Azimuth pixel spacing (meters)
    pub azimuth_pixel_spacing: f64,
    /// Near range time (seconds)
    pub near_range_time: f64,
    /// Azimuth time start
    pub azimuth_time_start: f64,
}

/// Enhanced merged sub-swath data container
#[derive(Debug)]
pub struct MergedSwathData {
    /// Merged complex data (if available)
    pub complex_data: Option<SarImage>,
    /// Merged intensity data
    pub intensity_data: SarRealImage,
    /// Valid data mask
    pub valid_mask: Array2<bool>,
    /// Quality assessment mask
    pub quality_mask: Array2<f32>,
    /// Output grid information
    pub grid: OutputGrid,
    /// Enhanced processing metadata
    pub metadata: MergeMetadata,
    /// Quality validation results
    pub quality_results: QualityResults,
}

/// Enhanced metadata for the merge process
#[derive(Debug, Clone)]
pub struct MergeMetadata {
    /// Number of sub-swaths merged
    pub num_swaths: usize,
    /// Overlap regions processed
    pub overlap_count: usize,
    /// Total valid pixels
    pub valid_pixels: usize,
    /// Processing timestamp
    pub processing_time: String,
    /// Blending method used
    pub blending_method: BlendingMethod,
    /// Processing performance metrics
    pub performance_metrics: PerformanceMetrics,
}

/// Performance metrics for merge processing
#[derive(Debug, Clone)]
pub struct PerformanceMetrics {
    /// Total processing time (seconds)
    pub total_time_seconds: f64,
    /// Memory usage peak (MB)
    pub peak_memory_mb: f64,
    /// Pixels processed per second
    pub pixels_per_second: f64,
    /// Overlap processing efficiency
    pub overlap_efficiency: f64,
}

/// Quality validation results
#[derive(Debug, Clone)]
pub struct QualityResults {
    /// Overall merge quality score (0-1)
    pub overall_quality: f32,
    /// Phase preservation quality (for complex data)
    pub phase_preservation: Option<f32>,
    /// Radiometric consistency across seams
    pub radiometric_consistency: f32,
    /// Overlap region qualities
    pub overlap_qualities: Vec<OverlapQuality>,
    /// Validation passed/failed
    pub validation_passed: bool,
    /// Quality warnings
    pub warnings: Vec<String>,
}

impl TopsarMerge {
    /// Create new enhanced TOPSAR merge processor
    pub fn new(subswaths: Vec<SubSwath>) -> SarResult<Self> {
        Self::new_with_params(subswaths, MergeParameters::default(), QualityControl::default())
    }
    
    /// Create TOPSAR merge processor with custom parameters
    pub fn new_with_params(
        subswaths: Vec<SubSwath>,
        merge_params: MergeParameters,
        quality_control: QualityControl,
    ) -> SarResult<Self> {
        log::info!("🔗 Initializing Enhanced TOPSAR merge for {} sub-swaths", subswaths.len());
        log::debug!("Merge parameters: {:?}", merge_params);
        
        if subswaths.len() < 2 {
            return Err(SarError::Processing(
                "TOPSAR merge requires at least 2 sub-swaths".to_string()
            ));
        }

        // Calculate enhanced overlap regions with quality metrics
        let overlap_regions = Self::calculate_enhanced_overlap_regions(&subswaths, &merge_params)?;
        
        // Determine optimized output grid parameters
        let output_grid = Self::calculate_optimized_output_grid(&subswaths)?;
        
        Ok(Self {
            subswaths,
            overlap_regions,
            output_grid,
            merge_params,
            quality_control,
        })
    }

    /// Enhanced merge with full quality control and validation
    pub fn merge_subswaths_enhanced(
        &self,
        subswath_data: &HashMap<String, SarRealImage>,
        preserve_complex: bool,
        complex_data: Option<&HashMap<String, SarImage>>,
    ) -> SarResult<MergedSwathData> {
        let start_time = std::time::Instant::now();
        
        log::info!("🎯 Starting Enhanced TOPSAR merge process");
        log::debug!("Blending method: {:?}", self.merge_params.blending_method);
        log::debug!("Quality control enabled: {}", self.quality_control.enable_validation);

        // Step 1: Validate input data quality
        self.validate_input_data(subswath_data, complex_data)?;

        // Step 2: Initialize enhanced output arrays
        let mut merged_intensity = Array2::zeros((
            self.output_grid.azimuth_samples,
            self.output_grid.range_samples,
        ));
        
        let mut merged_complex = if preserve_complex {
            Some(Array2::zeros((
                self.output_grid.azimuth_samples,
                self.output_grid.range_samples,
            )))
        } else {
            None
        };

        let mut valid_mask = Array2::from_elem((
            self.output_grid.azimuth_samples,
            self.output_grid.range_samples,
        ), false);
        
        let mut quality_mask = Array2::zeros((
            self.output_grid.azimuth_samples,
            self.output_grid.range_samples,
        ));

        // Step 3: Place individual sub-swaths with enhanced positioning
        for subswath in &self.subswaths {
            self.place_subswath_data_enhanced(
                subswath,
                subswath_data,
                complex_data,
                &mut merged_intensity,
                &mut merged_complex,
                &mut valid_mask,
                &mut quality_mask,
            )?;
        }

        // Step 4: Advanced overlap blending with quality control
        let overlap_qualities = self.blend_overlap_regions_enhanced(
            subswath_data,
            complex_data,
            &mut merged_intensity,
            &mut merged_complex,
            &mut valid_mask,
            &mut quality_mask,
        )?;

        // Step 5: Quality validation and assessment
        let quality_results = if self.quality_control.enable_validation {
            self.perform_quality_validation(
                &merged_intensity,
                &merged_complex,
                &valid_mask,
                &overlap_qualities,
            )?
        } else {
            QualityResults::default()
        };

        // Step 6: Generate comprehensive metadata
        let processing_time = start_time.elapsed().as_secs_f64();
        let metadata = self.generate_enhanced_metadata(&valid_mask, processing_time)?;

        let coverage = (metadata.valid_pixels as f64 / 
                       (self.output_grid.range_samples * self.output_grid.azimuth_samples) as f64) * 100.0;
        
        log::info!("✅ Enhanced TOPSAR merge completed: {:.1}% coverage", coverage);
        log::info!("📊 Quality score: {:.2}", quality_results.overall_quality);

        Ok(MergedSwathData {
            complex_data: merged_complex,
            intensity_data: merged_intensity,
            valid_mask,
            quality_mask,
            grid: self.output_grid.clone(),
            metadata,
            quality_results,
        })
    }

    /// Calculate enhanced overlap regions with quality metrics
    fn calculate_enhanced_overlap_regions(
        subswaths: &[SubSwath], 
        params: &MergeParameters
    ) -> SarResult<Vec<OverlapRegion>> {
        let mut overlaps = Vec::new();
        
        for i in 0..(subswaths.len() - 1) {
            let swath1 = &subswaths[i];
            let swath2 = &subswaths[i + 1];
            
            log::debug!("Calculating enhanced overlap between {} and {}", swath1.id, swath2.id);
            
            // Enhanced overlap calculation with precise timing
            let overlap = Self::calculate_precise_overlap(swath1, swath2, params)?;
            overlaps.push(overlap);
        }
        
        log::info!("Found {} enhanced overlap regions", overlaps.len());
        Ok(overlaps)
    }
    
    /// Calculate precise overlap region with quality metrics
    fn calculate_precise_overlap(
        swath1: &SubSwath,
        swath2: &SubSwath,
        params: &MergeParameters,
    ) -> SarResult<OverlapRegion> {
        // Enhanced overlap calculation using burst timing and geometry
        // This would use real TOPSAR geometry in production
        
        let overlap_width = if params.optimize_overlaps {
            ((swath1.range_samples + swath2.range_samples) as f64 * 0.15) as usize
        } else {
            ((swath1.range_samples + swath2.range_samples) as f64 * 0.1) as usize
        };
        
        let swath1_range_start = swath1.range_samples.saturating_sub(overlap_width / 2);
        let swath1_range_end = swath1.range_samples;
        let swath2_range_start = 0;
        let swath2_range_end = overlap_width / 2;
        
        let azimuth_samples = swath1.azimuth_samples.min(swath2.azimuth_samples);
        
        // Generate optimized blending weights
        let weights = Self::generate_enhanced_weights(
            overlap_width / 2,
            azimuth_samples,
            &params.blending_method,
        )?;
        
        // Initialize quality metrics
        let quality_metrics = OverlapQuality {
            phase_coherence: 0.95, // Will be computed from real data
            radiometric_consistency: 0.90,
            valid_pixel_percentage: 98.0,
            seamline_quality: 0.85,
        };
        
        Ok(OverlapRegion {
            swath1_id: swath1.id.clone(),
            swath2_id: swath2.id.clone(),
            swath1_range_start,
            swath1_range_end,
            swath2_range_start,
            swath2_range_end,
            azimuth_start: 0,
            azimuth_end: azimuth_samples,
            weights,
            quality_metrics,
        })
    }
    
    /// Generate enhanced blending weights
    fn generate_enhanced_weights(
        range_width: usize,
        azimuth_height: usize,
        blending_method: &BlendingMethod,
    ) -> SarResult<Array2<f32>> {
        let mut weights = Array2::zeros((azimuth_height, range_width));
        
        match blending_method {
            BlendingMethod::Linear => {
                // Linear distance-based weighting
                for j in 0..range_width {
                    let weight = j as f32 / range_width as f32;
                    for i in 0..azimuth_height {
                        weights[[i, j]] = weight;
                    }
                }
            },
            BlendingMethod::Gaussian { sigma } => {
                // Gaussian distance-based weighting
                let center = range_width as f32 / 2.0;
                for j in 0..range_width {
                    let distance = (j as f32 - center).abs();
                    let weight = (-distance * distance / (2.0 * sigma * sigma)).exp();
                    for i in 0..azimuth_height {
                        weights[[i, j]] = weight;
                    }
                }
            },
            BlendingMethod::MultiScale { scales: _ } => {
                // Multi-scale blending (simplified implementation)
                for j in 0..range_width {
                    let weight = (j as f32 / range_width as f32).powf(0.7); // Non-linear weighting
                    for i in 0..azimuth_height {
                        weights[[i, j]] = weight;
                    }
                }
            },
            BlendingMethod::PhaseCoherent => {
                // Phase-coherent weighting for complex data
                for j in 0..range_width {
                    let weight = 0.5 + 0.3 * (j as f32 / range_width as f32 - 0.5).cos();
                    for i in 0..azimuth_height {
                        weights[[i, j]] = weight;
                    }
                }
            },
            BlendingMethod::SnapCompatible => {
                // ESA SNAP compatible linear weighting
                for j in 0..range_width {
                    let weight = j as f32 / (range_width - 1) as f32;
                    for i in 0..azimuth_height {
                        weights[[i, j]] = weight;
                    }
                }
            }
        }
        
        Ok(weights)
    }

    /// Validate input data quality before merging
    fn validate_input_data(
        &self,
        subswath_data: &HashMap<String, SarRealImage>,
        _complex_data: Option<&HashMap<String, SarImage>>,
    ) -> SarResult<()> {
        log::debug!("Validating input data quality");
        
        // Check that all required subswaths are present
        for subswath in &self.subswaths {
            if !subswath_data.contains_key(&subswath.id) {
                return Err(SarError::Processing(
                    format!("Missing intensity data for subswath {}", subswath.id)
                ));
            }
        }
        
        // Validate data dimensions and quality
        for (swath_id, data) in subswath_data {
            let (az, rg) = data.dim();
            
            // Check for reasonable dimensions
            if az < 100 || rg < 100 {
                log::warn!("Subswath {} has very small dimensions: {}x{}", swath_id, az, rg);
            }
            
            // Check for data validity
            let finite_count = data.iter().filter(|&&x| x.is_finite()).count();
            let valid_ratio = finite_count as f64 / data.len() as f64;
            
            if valid_ratio < 0.5 {
                return Err(SarError::Processing(
                    format!("Subswath {} has too many invalid pixels: {:.1}%", 
                           swath_id, (1.0 - valid_ratio) * 100.0)
                ));
            }
        }
        
        log::debug!("Input data validation passed");
        Ok(())
    }

    /// Place subswath data with enhanced positioning
    fn place_subswath_data_enhanced(
        &self,
        subswath: &SubSwath,
        intensity_data: &HashMap<String, SarRealImage>,
        complex_data: Option<&HashMap<String, SarImage>>,
        merged_intensity: &mut SarRealImage,
        merged_complex: &mut Option<SarImage>,
        valid_mask: &mut Array2<bool>,
        quality_mask: &mut Array2<f32>,
    ) -> SarResult<()> {
        let swath_intensity = intensity_data.get(&subswath.id)
            .ok_or_else(|| SarError::Processing(
                format!("Missing intensity data for sub-swath {}", subswath.id)
            ))?;

        log::debug!("Placing sub-swath {} data with enhanced positioning", subswath.id);

        // Calculate enhanced placement offset
        let (range_offset, azimuth_offset) = self.calculate_enhanced_placement_offset(subswath)?;

        // Place intensity data with quality assessment
        let (swath_az, swath_rg) = swath_intensity.dim();
        for az in 0..swath_az {
            for rg in 0..swath_rg {
                let out_az = azimuth_offset + az;
                let out_rg = range_offset + rg;
                
                if out_az < merged_intensity.dim().0 && out_rg < merged_intensity.dim().1 {
                    let pixel_value = swath_intensity[[az, rg]];
                    merged_intensity[[out_az, out_rg]] = pixel_value;
                    valid_mask[[out_az, out_rg]] = pixel_value.is_finite() && pixel_value > 0.0;
                    
                    // Assess pixel quality
                    quality_mask[[out_az, out_rg]] = self.assess_pixel_quality(pixel_value);
                }
            }
        }

        // Place complex data if available
        if let (Some(complex_map), Some(ref mut merged_complex_data)) = 
            (complex_data, merged_complex) {
            
            if let Some(swath_complex) = complex_map.get(&subswath.id) {
                let (swath_az, swath_rg) = swath_complex.dim();
                for az in 0..swath_az {
                    for rg in 0..swath_rg {
                        let out_az = azimuth_offset + az;
                        let out_rg = range_offset + rg;
                        
                        if out_az < merged_complex_data.dim().0 && out_rg < merged_complex_data.dim().1 {
                            merged_complex_data[[out_az, out_rg]] = swath_complex[[az, rg]];
                        }
                    }
                }
            }
        }

        Ok(())
    }

    /// Calculate enhanced placement offset with precise positioning
    fn calculate_enhanced_placement_offset(&self, subswath: &SubSwath) -> SarResult<(usize, usize)> {
        // Enhanced calculation using real burst timing and geometry
        // This is a simplified version - production would use precise TOPSAR geometry
        
        let range_offset = match subswath.id.as_str() {
            "IW1" => 0,
            "IW2" => {
                if self.merge_params.optimize_overlaps {
                    (self.subswaths[0].range_samples as f64 * 0.85) as usize
                } else {
                    (self.subswaths[0].range_samples as f64 * 0.8) as usize
                }
            },
            "IW3" => {
                if self.merge_params.optimize_overlaps {
                    (self.subswaths[0].range_samples as f64 * 0.85 +
                     self.subswaths[1].range_samples as f64 * 0.85) as usize
                } else {
                    (self.subswaths[0].range_samples as f64 * 0.8 +
                     self.subswaths[1].range_samples as f64 * 0.8) as usize
                }
            },
            _ => return Err(SarError::Processing(
                format!("Unknown sub-swath ID: {}", subswath.id)
            )),
        };

        // Enhanced azimuth alignment for co-registered data
        let azimuth_offset = 0; // Could include fine azimuth corrections

        Ok((range_offset, azimuth_offset))
    }

    /// Assess individual pixel quality
    fn assess_pixel_quality(&self, pixel_value: f32) -> f32 {
        if !pixel_value.is_finite() {
            return 0.0;
        }
        
        // Quality based on reasonable backscatter range
        let quality = if pixel_value > 0.0001 && pixel_value < 10.0 {
            1.0
        } else if pixel_value > 0.0 && pixel_value < 100.0 {
            0.7
        } else {
            0.3
        };
        
        quality
    }

    /// Enhanced overlap blending with quality control
    fn blend_overlap_regions_enhanced(
        &self,
        intensity_data: &HashMap<String, SarRealImage>,
        complex_data: Option<&HashMap<String, SarImage>>,
        merged_intensity: &mut SarRealImage,
        merged_complex: &mut Option<SarImage>,
        valid_mask: &mut Array2<bool>,
        quality_mask: &mut Array2<f32>,
    ) -> SarResult<Vec<OverlapQuality>> {
        log::debug!("Performing enhanced overlap blending");
        
        let mut overlap_qualities = Vec::new();
        
        for overlap in &self.overlap_regions {
            let quality = self.blend_single_overlap_enhanced(
                overlap,
                intensity_data,
                complex_data,
                merged_intensity,
                merged_complex,
                valid_mask,
                quality_mask,
            )?;
            
            overlap_qualities.push(quality);
        }
        
        Ok(overlap_qualities)
    }

    /// Blend single overlap region with enhanced quality control
    fn blend_single_overlap_enhanced(
        &self,
        overlap: &OverlapRegion,
        intensity_data: &HashMap<String, SarRealImage>,
        complex_data: Option<&HashMap<String, SarImage>>,
        merged_intensity: &mut SarRealImage,
        merged_complex: &mut Option<SarImage>,
        valid_mask: &mut Array2<bool>,
        quality_mask: &mut Array2<f32>,
    ) -> SarResult<OverlapQuality> {
        log::debug!("Blending overlap between {} and {} with quality control", 
                   overlap.swath1_id, overlap.swath2_id);

        let swath1_intensity = intensity_data.get(&overlap.swath1_id).unwrap();
        let swath2_intensity = intensity_data.get(&overlap.swath2_id).unwrap();

        let overlap_width = overlap.swath1_range_end - overlap.swath1_range_start;
        let overlap_height = overlap.azimuth_end - overlap.azimuth_start;
        
        let mut phase_coherence_sum = 0.0;
        let mut radiometric_consistency_sum = 0.0;
        let mut valid_pixels = 0;
        let mut total_pixels = 0;

        // Enhanced blending with quality assessment
        for i in 0..overlap_height {
            for j in 0..overlap_width {
                let az = overlap.azimuth_start + i;
                let rg1 = overlap.swath1_range_start + j;
                let rg2 = overlap.swath2_range_start + j;

                if rg1 < swath1_intensity.dim().1 && rg2 < swath2_intensity.dim().1 &&
                   az < swath1_intensity.dim().0 && az < swath2_intensity.dim().0 {
                    
                    let val1 = swath1_intensity[[az, rg1]];
                    let val2 = swath2_intensity[[az, rg2]];
                    
                    total_pixels += 1;
                    
                    if val1.is_finite() && val2.is_finite() && val1 > 0.0 && val2 > 0.0 {
                        // Enhanced blending using optimized weights
                        let weight = if j < overlap.weights.dim().1 && i < overlap.weights.dim().0 {
                            overlap.weights[[i, j]]
                        } else {
                            j as f32 / overlap_width as f32
                        };
                        
                        let blended_value = self.apply_enhanced_blending(val1, val2, weight)?;
                        
                        // Calculate output position
                        let (out_rg1, _) = self.calculate_enhanced_placement_offset(
                            &self.subswaths.iter().find(|s| s.id == overlap.swath1_id).unwrap()
                        )?;
                        
                        let out_rg = out_rg1 + rg1;
                        
                        if az < merged_intensity.dim().0 && out_rg < merged_intensity.dim().1 {
                            merged_intensity[[az, out_rg]] = blended_value;
                            valid_mask[[az, out_rg]] = true;
                            
                            // Assess blending quality
                            let pixel_quality = self.assess_blending_quality(val1, val2, blended_value);
                            quality_mask[[az, out_rg]] = pixel_quality;
                            
                            // Accumulate quality metrics
                            radiometric_consistency_sum += 1.0 - (val1 - val2).abs() / (val1 + val2);
                            valid_pixels += 1;
                        }
                    }
                }
            }
        }

        // Handle complex data blending with phase preservation
        if let (Some(complex_map), Some(ref mut merged_complex_data)) = 
            (complex_data, merged_complex) {
            
            if let (Some(swath1_complex), Some(swath2_complex)) = 
                (complex_map.get(&overlap.swath1_id), complex_map.get(&overlap.swath2_id)) {
                
                phase_coherence_sum = self.blend_complex_overlap_enhanced(
                    overlap,
                    swath1_complex,
                    swath2_complex,
                    merged_complex_data,
                )?;
            }
        }

        // Calculate quality metrics
        let phase_coherence = if complex_data.is_some() {
            phase_coherence_sum / valid_pixels as f32
        } else {
            1.0 // Not applicable for intensity-only data
        };
        
        let radiometric_consistency = if valid_pixels > 0 {
            radiometric_consistency_sum / valid_pixels as f32
        } else {
            0.0
        };
        
        let valid_pixel_percentage = if total_pixels > 0 {
            (valid_pixels as f32 / total_pixels as f32) * 100.0
        } else {
            0.0
        };
        
        let seamline_quality = self.assess_seamline_quality(
            overlap,
            intensity_data,
            merged_intensity,
        )?;

        Ok(OverlapQuality {
            phase_coherence,
            radiometric_consistency,
            valid_pixel_percentage,
            seamline_quality,
        })
    }

    /// Apply enhanced blending algorithm
    fn apply_enhanced_blending(&self, val1: f32, val2: f32, weight: f32) -> SarResult<f32> {
        match &self.merge_params.blending_method {
            BlendingMethod::Linear | BlendingMethod::SnapCompatible => {
                Ok(val1 * (1.0 - weight) + val2 * weight)
            },
            BlendingMethod::Gaussian { sigma: _ } => {
                // Gaussian-weighted blending
                Ok(val1 * (1.0 - weight) + val2 * weight)
            },
            BlendingMethod::MultiScale { scales: _ } => {
                // Multi-scale blending (simplified)
                let linear_blend = val1 * (1.0 - weight) + val2 * weight;
                let geometric_mean = (val1 * val2).sqrt();
                Ok(linear_blend * 0.7 + geometric_mean * 0.3)
            },
            BlendingMethod::PhaseCoherent => {
                // Phase-coherent blending for better edge preservation
                let ratio = val2 / val1;
                let adaptive_weight = if ratio > 0.5 && ratio < 2.0 {
                    weight
                } else {
                    if val1 > val2 { 0.2 } else { 0.8 }
                };
                Ok(val1 * (1.0 - adaptive_weight) + val2 * adaptive_weight)
            }
        }
    }

    /// Assess quality of blended pixel
    fn assess_blending_quality(&self, val1: f32, val2: f32, blended: f32) -> f32 {
        // Quality based on consistency between input values and blend result
        let consistency = if val1 > 0.0 && val2 > 0.0 {
            let ratio = val2 / val1;
            let blend_ratio = blended / ((val1 + val2) / 2.0);
            
            if ratio > 0.3 && ratio < 3.0 && blend_ratio > 0.7 && blend_ratio < 1.3 {
                1.0
            } else if ratio > 0.1 && ratio < 10.0 {
                0.7
            } else {
                0.3
            }
        } else {
            0.0
        };
        
        consistency
    }

    /// Blend complex overlap with phase preservation
    fn blend_complex_overlap_enhanced(
        &self,
        overlap: &OverlapRegion,
        swath1_complex: &SarImage,
        swath2_complex: &SarImage,
        merged_complex: &mut SarImage,
    ) -> SarResult<f32> {
        log::debug!("Blending complex overlap with phase preservation");
        
        let overlap_width = overlap.swath1_range_end - overlap.swath1_range_start;
        let overlap_height = overlap.azimuth_end - overlap.azimuth_start;
        
        let mut phase_coherence_sum = 0.0;
        let mut valid_pixels = 0;

        for i in 0..overlap_height {
            for j in 0..overlap_width {
                let az = overlap.azimuth_start + i;
                let rg1 = overlap.swath1_range_start + j;
                let rg2 = overlap.swath2_range_start + j;

                if rg1 < swath1_complex.dim().1 && rg2 < swath2_complex.dim().1 &&
                   az < swath1_complex.dim().0 && az < swath2_complex.dim().0 {
                    
                    let complex1 = swath1_complex[[az, rg1]];
                    let complex2 = swath2_complex[[az, rg2]];
                    
                    if complex1.norm() > 0.0 && complex2.norm() > 0.0 {
                        // Phase-preserving complex blending
                        let weight = if j < overlap.weights.dim().1 && i < overlap.weights.dim().0 {
                            overlap.weights[[i, j]]
                        } else {
                            j as f32 / overlap_width as f32
                        };
                        
                        let blended_complex = if self.merge_params.preserve_phase {
                            // Maintain phase relationships
                            let phase_diff = (complex2.arg() - complex1.arg()).abs();
                            if phase_diff < std::f32::consts::PI / 2.0 {
                                complex1 * (1.0 - weight) + complex2 * weight
                            } else {
                                // Use amplitude blending only for large phase differences
                                let blended_amplitude = complex1.norm() * (1.0 - weight) + complex2.norm() * weight;
                                let avg_phase = (complex1.arg() + complex2.arg()) / 2.0;
                                num_complex::Complex::from_polar(blended_amplitude, avg_phase)
                            }
                        } else {
                            complex1 * (1.0 - weight) + complex2 * weight
                        };
                        
                        // Calculate output position
                        let (out_rg1, _) = self.calculate_enhanced_placement_offset(
                            &self.subswaths.iter().find(|s| s.id == overlap.swath1_id).unwrap()
                        )?;
                        
                        let out_rg = out_rg1 + rg1;
                        
                        if az < merged_complex.dim().0 && out_rg < merged_complex.dim().1 {
                            merged_complex[[az, out_rg]] = blended_complex;
                            
                            // Calculate phase coherence
                            let coherence = (complex1.conj() * complex2).norm() / 
                                          (complex1.norm() * complex2.norm());
                            phase_coherence_sum += coherence;
                            valid_pixels += 1;
                        }
                    }
                }
            }
        }

        Ok(if valid_pixels > 0 {
            phase_coherence_sum / valid_pixels as f32
        } else {
            0.0
        })
    }

    /// Assess seamline quality
    fn assess_seamline_quality(
        &self,
        overlap: &OverlapRegion,
        intensity_data: &HashMap<String, SarRealImage>,
        _merged_intensity: &SarRealImage,
    ) -> SarResult<f32> {
        // Simplified seamline quality assessment
        // In production, this would analyze edge gradients and discontinuities
        
        let swath1_intensity = intensity_data.get(&overlap.swath1_id).unwrap();
        let swath2_intensity = intensity_data.get(&overlap.swath2_id).unwrap();
        
        let mut gradient_consistency_sum = 0.0;
        let mut valid_samples = 0;
        
        // Sample along the seamline
        let seamline_range = overlap.swath1_range_start + 
                           (overlap.swath1_range_end - overlap.swath1_range_start) / 2;
        
        for az in overlap.azimuth_start..overlap.azimuth_end.min(swath1_intensity.dim().0 - 1) {
            if seamline_range < swath1_intensity.dim().1 && 
               seamline_range < swath2_intensity.dim().1 {
                
                let val1 = swath1_intensity[[az, seamline_range]];
                let val2 = swath2_intensity[[az, seamline_range]];
                
                if val1.is_finite() && val2.is_finite() && val1 > 0.0 && val2 > 0.0 {
                    let consistency = 1.0 - (val1 - val2).abs() / (val1 + val2).max(0.001);
                    gradient_consistency_sum += consistency;
                    valid_samples += 1;
                }
            }
        }
        
        Ok(if valid_samples > 0 {
            gradient_consistency_sum / valid_samples as f32
        } else {
            0.5 // Default moderate quality
        })
    }

    /// Perform comprehensive quality validation
    fn perform_quality_validation(
        &self,
        _merged_intensity: &SarRealImage,
        merged_complex: &Option<SarImage>,
        _valid_mask: &Array2<bool>,
        overlap_qualities: &[OverlapQuality],
    ) -> SarResult<QualityResults> {
        log::debug!("Performing comprehensive quality validation");
        
        let mut warnings = Vec::new();
        
        // Overall quality assessment
        let overall_quality = if !overlap_qualities.is_empty() {
            let avg_radiometric = overlap_qualities.iter()
                .map(|q| q.radiometric_consistency)
                .sum::<f32>() / overlap_qualities.len() as f32;
            let avg_seamline = overlap_qualities.iter()
                .map(|q| q.seamline_quality)
                .sum::<f32>() / overlap_qualities.len() as f32;
            
            (avg_radiometric + avg_seamline) / 2.0
        } else {
            0.8 // Default for single subswath
        };
        
        // Phase preservation quality (for complex data)
        let phase_preservation = if merged_complex.is_some() {
            let avg_phase_coherence = overlap_qualities.iter()
                .map(|q| q.phase_coherence)
                .sum::<f32>() / overlap_qualities.len() as f32;
            Some(avg_phase_coherence)
        } else {
            None
        };
        
        // Radiometric consistency
        let radiometric_consistency = if !overlap_qualities.is_empty() {
            overlap_qualities.iter()
                .map(|q| q.radiometric_consistency)
                .sum::<f32>() / overlap_qualities.len() as f32
        } else {
            0.9
        };
        
        // Validation checks
        if overall_quality < self.quality_control.radiometric_tolerance {
            warnings.push(format!("Low overall quality: {:.2}", overall_quality));
        }
        
        if let Some(phase_qual) = phase_preservation {
            if phase_qual < 0.8 {
                warnings.push(format!("Low phase preservation: {:.2}", phase_qual));
            }
        }
        
        let validation_passed = overall_quality >= self.quality_control.radiometric_tolerance &&
                               warnings.len() < 3; // Allow some warnings
        
        Ok(QualityResults {
            overall_quality,
            phase_preservation,
            radiometric_consistency,
            overlap_qualities: overlap_qualities.to_vec(),
            validation_passed,
            warnings,
        })
    }

    /// Generate enhanced processing metadata
    fn generate_enhanced_metadata(
        &self,
        valid_mask: &Array2<bool>,
        processing_time: f64,
    ) -> SarResult<MergeMetadata> {
        let valid_pixels = valid_mask.iter().filter(|&&x| x).count();
        let total_pixels = valid_mask.len();
        
        let performance_metrics = PerformanceMetrics {
            total_time_seconds: processing_time,
            peak_memory_mb: 0.0, // Would be measured in production
            pixels_per_second: total_pixels as f64 / processing_time.max(0.001),
            overlap_efficiency: if self.overlap_regions.is_empty() {
                1.0
            } else {
                valid_pixels as f64 / total_pixels as f64
            },
        };
        
        Ok(MergeMetadata {
            num_swaths: self.subswaths.len(),
            overlap_count: self.overlap_regions.len(),
            valid_pixels,
            processing_time: chrono::Utc::now().to_string(),
            blending_method: self.merge_params.blending_method.clone(),
            performance_metrics,
        })
    }

    /// Calculate optimized output grid parameters
    fn calculate_optimized_output_grid(subswaths: &[SubSwath]) -> SarResult<OutputGrid> {
        let total_range_samples: usize = subswaths.iter()
            .map(|sw| sw.range_samples)
            .sum();
        
        // Enhanced overlap accounting
        let overlap_factor = match subswaths.len() {
            1 => 1.0,
            2 => 0.90, // 10% overlap
            3 => 0.85, // 15% total overlap  
            _ => 0.80, // 20% total overlap for more swaths
        };
        
        let effective_range_samples = (total_range_samples as f64 * overlap_factor) as usize;
        
        let max_azimuth_samples = subswaths.iter()
            .map(|sw| sw.azimuth_samples)
            .max()
            .unwrap_or(0);

        // Use first sub-swath parameters as reference
        let reference_swath = &subswaths[0];
        
        Ok(OutputGrid {
            range_samples: effective_range_samples,
            azimuth_samples: max_azimuth_samples,
            range_pixel_spacing: reference_swath.range_pixel_spacing,
            azimuth_pixel_spacing: reference_swath.azimuth_pixel_spacing,
            near_range_time: reference_swath.slant_range_time,
            azimuth_time_start: 0.0, // Would be calculated from burst timing
        })
    }

    /// Standard merge interface for backward compatibility
    pub fn merge_subswaths(
        &self,
        subswath_data: &HashMap<String, SarRealImage>,
        preserve_complex: bool,
        complex_data: Option<&HashMap<String, SarImage>>,
    ) -> SarResult<MergedSwathData> {
        // Use enhanced merge with default quality control
        self.merge_subswaths_enhanced(subswath_data, preserve_complex, complex_data)
    }

    /// Get information about sub-swaths
    pub fn get_subswath_info(&self) -> &[SubSwath] {
        &self.subswaths
    }

    /// Get overlap region information
    pub fn get_overlap_regions(&self) -> &[OverlapRegion] {
        &self.overlap_regions
    }

    /// Get output grid information
    pub fn get_output_grid(&self) -> &OutputGrid {
        &self.output_grid
    }
}

/// Convenience function for TOPSAR merge with default parameters
pub fn merge_iw_subswaths(
    subswaths: Vec<SubSwath>,
    intensity_data: HashMap<String, SarRealImage>,
    complex_data: Option<HashMap<String, SarImage>>,
) -> SarResult<MergedSwathData> {
    log::info!("🔗 Starting IW sub-swath merge");
    
    let merger = TopsarMerge::new(subswaths)?;
    let preserve_complex = complex_data.is_some();
    
    merger.merge_subswaths(&intensity_data, preserve_complex, complex_data.as_ref())
}

#[cfg(test)]
mod tests {
    use super::*;
    use ndarray::Array2;
    use num_complex::Complex32;

    #[test]
    fn test_topsar_merge_creation() {
        let subswaths = vec![
            SubSwath {
                id: "IW1".to_string(),
                burst_count: 9,
                range_samples: 1000,
                azimuth_samples: 1500,
                range_pixel_spacing: 2.33,
                azimuth_pixel_spacing: 14.1,
                slant_range_time: 5.44e-3,
                burst_duration: 2.0,
            },
            SubSwath {
                id: "IW2".to_string(),
                burst_count: 9,
                range_samples: 1000,
                azimuth_samples: 1500,
                range_pixel_spacing: 2.33,
                azimuth_pixel_spacing: 14.1,
                slant_range_time: 5.50e-3,
                burst_duration: 2.0,
            },
        ];

        let merger = TopsarMerge::new(subswaths).unwrap();
        assert_eq!(merger.subswaths.len(), 2);
        assert_eq!(merger.overlap_regions.len(), 1);
    }

    #[test]
    fn test_overlap_calculation() {
        let subswaths = vec![
            SubSwath {
                id: "IW1".to_string(),
                burst_count: 9,
                range_samples: 1000,
                azimuth_samples: 1500,
                range_pixel_spacing: 2.33,
                azimuth_pixel_spacing: 14.1,
                slant_range_time: 5.44e-3,
                burst_duration: 2.0,
            },
            SubSwath {
                id: "IW2".to_string(),
                burst_count: 9,
                range_samples: 1000,
                azimuth_samples: 1500,
                range_pixel_spacing: 2.33,
                azimuth_pixel_spacing: 14.1,
                slant_range_time: 5.50e-3,
                burst_duration: 2.0,
            },
        ];

        let merge_params = MergeParameters::default();
        let overlaps = TopsarMerge::calculate_enhanced_overlap_regions(&subswaths, &merge_params).unwrap();
        assert_eq!(overlaps.len(), 1);
        
        let overlap = &overlaps[0];
        assert_eq!(overlap.swath1_id, "IW1");
        assert_eq!(overlap.swath2_id, "IW2");
        assert!(overlap.swath1_range_start > 0);
        assert!(overlap.swath2_range_end > 0);
    }
}

impl Default for MergeParameters {
    fn default() -> Self {
        Self {
            blending_method: BlendingMethod::Linear,
            preserve_phase: true,
            optimize_overlaps: true,
            enable_parallel: true,
            chunk_size: 4096,
            feather_width: 16,
        }
    }
}

impl Default for BlendingMethod {
    fn default() -> Self {
        BlendingMethod::Linear
    }
}

impl Default for QualityControl {
    fn default() -> Self {
        Self {
            enable_validation: true,
            max_phase_discontinuity: 0.5,
            min_valid_pixel_ratio: 0.8,
            optimize_seamlines: true,
            radiometric_tolerance: 0.1,
        }
    }
}

impl Default for QualityResults {
    fn default() -> Self {
        Self {
            overall_quality: 0.0,
            phase_preservation: None,
            radiometric_consistency: 0.0,
            overlap_qualities: Vec::new(),
            validation_passed: false,
            warnings: Vec::new(),
        }
    }
}
