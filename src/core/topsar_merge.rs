use crate::types::{SarError, SarResult, SarImage, SarRealImage, SubSwath, BoundingBox};
use ndarray::{Array2, Array3};
use std::collections::HashMap;

/// TOPSAR merge processor for combining IW sub-swaths
pub struct TopsarMerge {
    /// Sub-swath information
    subswaths: Vec<SubSwath>,
    /// Overlap regions between adjacent sub-swaths
    overlap_regions: Vec<OverlapRegion>,
    /// Output grid parameters
    output_grid: OutputGrid,
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
    /// Blending weights for smooth transition
    pub weights: Array2<f32>,
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

/// Merged sub-swath data container
#[derive(Debug)]
pub struct MergedSwathData {
    /// Merged complex data (if available)
    pub complex_data: Option<SarImage>,
    /// Merged intensity data
    pub intensity_data: SarRealImage,
    /// Valid data mask
    pub valid_mask: Array2<bool>,
    /// Output grid information
    pub grid: OutputGrid,
    /// Processing metadata
    pub metadata: MergeMetadata,
}

/// Metadata for the merge process
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
}

impl TopsarMerge {
    /// Create new TOPSAR merge processor
    pub fn new(subswaths: Vec<SubSwath>) -> SarResult<Self> {
        log::info!("🔗 Initializing TOPSAR merge for {} sub-swaths", subswaths.len());
        
        if subswaths.len() < 2 {
            return Err(SarError::Processing(
                "TOPSAR merge requires at least 2 sub-swaths".to_string()
            ));
        }

        // Calculate overlap regions between adjacent sub-swaths
        let overlap_regions = Self::calculate_overlap_regions(&subswaths)?;
        
        // Determine output grid parameters
        let output_grid = Self::calculate_output_grid(&subswaths)?;
        
        Ok(Self {
            subswaths,
            overlap_regions,
            output_grid,
        })
    }

    /// Merge calibrated sub-swath data into single image
    pub fn merge_subswaths(
        &self,
        subswath_data: &HashMap<String, SarRealImage>,
        preserve_complex: bool,
        complex_data: Option<&HashMap<String, SarImage>>,
    ) -> SarResult<MergedSwathData> {
        log::info!("🎯 Starting TOPSAR merge process");
        log::debug!("Output grid: {}x{} pixels", 
                   self.output_grid.range_samples, 
                   self.output_grid.azimuth_samples);

        // Initialize output arrays
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

        // Step 1: Place individual sub-swaths in output grid
        for subswath in &self.subswaths {
            self.place_subswath_data(
                subswath,
                subswath_data,
                complex_data,
                &mut merged_intensity,
                &mut merged_complex,
                &mut valid_mask,
            )?;
        }

        // Step 2: Blend overlap regions for smooth transitions
        self.blend_overlap_regions(
            subswath_data,
            complex_data,
            &mut merged_intensity,
            &mut merged_complex,
            &mut valid_mask,
        )?;

        // Step 3: Generate metadata
        let metadata = self.generate_metadata(&valid_mask);

        let coverage = (metadata.valid_pixels as f64 / 
                       (self.output_grid.range_samples * self.output_grid.azimuth_samples) as f64) * 100.0;
        
        log::info!("✅ TOPSAR merge completed: {:.1}% coverage", coverage);

        Ok(MergedSwathData {
            complex_data: merged_complex,
            intensity_data: merged_intensity,
            valid_mask,
            grid: self.output_grid.clone(),
            metadata,
        })
    }

    /// Calculate overlap regions between adjacent sub-swaths
    fn calculate_overlap_regions(subswaths: &[SubSwath]) -> SarResult<Vec<OverlapRegion>> {
        let mut overlaps = Vec::new();
        
        for i in 0..(subswaths.len() - 1) {
            let swath1 = &subswaths[i];
            let swath2 = &subswaths[i + 1];
            
            log::debug!("Calculating overlap between {} and {}", swath1.id, swath2.id);
            
            // Simplified overlap calculation - in reality this would use
            // precise burst timing and geometry information
            let overlap = OverlapRegion {
                swath1_id: swath1.id.clone(),
                swath2_id: swath2.id.clone(),
                swath1_range_start: (swath1.range_samples as f64 * 0.8) as usize,
                swath1_range_end: swath1.range_samples,
                swath2_range_start: 0,
                swath2_range_end: (swath2.range_samples as f64 * 0.2) as usize,
                azimuth_start: 0,
                azimuth_end: swath1.azimuth_samples.min(swath2.azimuth_samples),
                weights: Array2::zeros((0, 0)), // Will be calculated during blending
            };
            
            overlaps.push(overlap);
        }
        
        log::info!("Found {} overlap regions", overlaps.len());
        Ok(overlaps)
    }

    /// Calculate output grid parameters
    fn calculate_output_grid(subswaths: &[SubSwath]) -> SarResult<OutputGrid> {
        let total_range_samples: usize = subswaths.iter()
            .map(|sw| sw.range_samples)
            .sum();
        
        // Account for overlaps (simplified)
        let effective_range_samples = (total_range_samples as f64 * 0.85) as usize;
        
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

    /// Place individual sub-swath data in output grid
    fn place_subswath_data(
        &self,
        subswath: &SubSwath,
        intensity_data: &HashMap<String, SarRealImage>,
        complex_data: Option<&HashMap<String, SarImage>>,
        merged_intensity: &mut SarRealImage,
        merged_complex: &mut Option<SarImage>,
        valid_mask: &mut Array2<bool>,
    ) -> SarResult<()> {
        let swath_intensity = intensity_data.get(&subswath.id)
            .ok_or_else(|| SarError::Processing(
                format!("Missing intensity data for sub-swath {}", subswath.id)
            ))?;

        log::debug!("Placing sub-swath {} data", subswath.id);

        // Calculate placement offset in output grid
        let (range_offset, azimuth_offset) = self.calculate_placement_offset(subswath)?;

        // Copy intensity data
        let (swath_az, swath_rg) = swath_intensity.dim();
        for az in 0..swath_az {
            for rg in 0..swath_rg {
                let out_az = azimuth_offset + az;
                let out_rg = range_offset + rg;
                
                if out_az < merged_intensity.dim().0 && out_rg < merged_intensity.dim().1 {
                    merged_intensity[[out_az, out_rg]] = swath_intensity[[az, rg]];
                    valid_mask[[out_az, out_rg]] = true;
                }
            }
        }

        // Copy complex data if available
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

    /// Calculate placement offset for sub-swath in output grid
    fn calculate_placement_offset(&self, subswath: &SubSwath) -> SarResult<(usize, usize)> {
        // Simplified calculation - in reality would use precise timing
        let range_offset = match subswath.id.as_str() {
            "IW1" => 0,
            "IW2" => (self.subswaths[0].range_samples as f64 * 0.8) as usize,
            "IW3" => (self.subswaths[0].range_samples as f64 * 0.8 +
                     self.subswaths[1].range_samples as f64 * 0.8) as usize,
            _ => return Err(SarError::Processing(
                format!("Unknown sub-swath ID: {}", subswath.id)
            )),
        };

        Ok((range_offset, 0)) // Azimuth offset is 0 for co-registered data
    }

    /// Blend overlap regions for smooth transitions
    fn blend_overlap_regions(
        &self,
        intensity_data: &HashMap<String, SarRealImage>,
        complex_data: Option<&HashMap<String, SarImage>>,
        merged_intensity: &mut SarRealImage,
        merged_complex: &mut Option<SarImage>,
        valid_mask: &mut Array2<bool>,
    ) -> SarResult<()> {
        log::debug!("Blending {} overlap regions", self.overlap_regions.len());

        for overlap in &self.overlap_regions {
            self.blend_single_overlap(
                overlap,
                intensity_data,
                complex_data,
                merged_intensity,
                merged_complex,
                valid_mask,
            )?;
        }

        Ok(())
    }

    /// Blend a single overlap region
    fn blend_single_overlap(
        &self,
        overlap: &OverlapRegion,
        intensity_data: &HashMap<String, SarRealImage>,
        complex_data: Option<&HashMap<String, SarImage>>,
        merged_intensity: &mut SarRealImage,
        merged_complex: &mut Option<SarImage>,
        _valid_mask: &mut Array2<bool>,
    ) -> SarResult<()> {
        let swath1_data = intensity_data.get(&overlap.swath1_id)
            .ok_or_else(|| SarError::Processing(
                format!("Missing data for sub-swath {}", overlap.swath1_id)
            ))?;

        let swath2_data = intensity_data.get(&overlap.swath2_id)
            .ok_or_else(|| SarError::Processing(
                format!("Missing data for sub-swath {}", overlap.swath2_id)
            ))?;

        log::debug!("Blending overlap between {} and {}", 
                   overlap.swath1_id, overlap.swath2_id);

        // Calculate blend weights (linear transition)
        let overlap_width = overlap.swath1_range_end - overlap.swath1_range_start;
        
        for az in overlap.azimuth_start..overlap.azimuth_end {
            for i in 0..overlap_width {
                let rg1 = overlap.swath1_range_start + i;
                let rg2 = overlap.swath2_range_start + i;
                
                if az < swath1_data.dim().0 && rg1 < swath1_data.dim().1 &&
                   az < swath2_data.dim().0 && rg2 < swath2_data.dim().1 {
                    
                    // Linear blending weight
                    let weight1 = 1.0 - (i as f32 / overlap_width as f32);
                    let weight2 = i as f32 / overlap_width as f32;
                    
                    let val1 = swath1_data[[az, rg1]];
                    let val2 = swath2_data[[az, rg2]];
                    
                    // Calculate output position
                    let (out_rg1, _) = self.calculate_placement_offset(
                        &self.subswaths.iter().find(|s| s.id == overlap.swath1_id).unwrap()
                    )?;
                    
                    let out_rg = out_rg1 + rg1;
                    
                    if az < merged_intensity.dim().0 && out_rg < merged_intensity.dim().1 {
                        merged_intensity[[az, out_rg]] = val1 * weight1 + val2 * weight2;
                    }
                }
            }
        }

        // Blend complex data if available
        if let (Some(complex_map), Some(ref mut merged_complex_data)) = 
            (complex_data, merged_complex) {
            
            if let (Some(swath1_complex), Some(swath2_complex)) = 
                (complex_map.get(&overlap.swath1_id), complex_map.get(&overlap.swath2_id)) {
                
                for az in overlap.azimuth_start..overlap.azimuth_end {
                    for i in 0..overlap_width {
                        let rg1 = overlap.swath1_range_start + i;
                        let rg2 = overlap.swath2_range_start + i;
                        
                        if az < swath1_complex.dim().0 && rg1 < swath1_complex.dim().1 &&
                           az < swath2_complex.dim().0 && rg2 < swath2_complex.dim().1 {
                            
                            let weight1 = 1.0 - (i as f32 / overlap_width as f32);
                            let weight2 = i as f32 / overlap_width as f32;
                            
                            let val1 = swath1_complex[[az, rg1]];
                            let val2 = swath2_complex[[az, rg2]];
                            
                            let (out_rg1, _) = self.calculate_placement_offset(
                                &self.subswaths.iter().find(|s| s.id == overlap.swath1_id).unwrap()
                            )?;
                            
                            let out_rg = out_rg1 + rg1;
                            
                            if az < merged_complex_data.dim().0 && out_rg < merged_complex_data.dim().1 {
                                merged_complex_data[[az, out_rg]] = 
                                    val1 * weight1 + val2 * weight2;
                            }
                        }
                    }
                }
            }
        }

        Ok(())
    }

    /// Generate processing metadata
    fn generate_metadata(&self, valid_mask: &Array2<bool>) -> MergeMetadata {
        let valid_pixels = valid_mask.iter().filter(|&&v| v).count();
        
        MergeMetadata {
            num_swaths: self.subswaths.len(),
            overlap_count: self.overlap_regions.len(),
            valid_pixels,
            processing_time: chrono::Utc::now().format("%Y-%m-%d %H:%M:%S UTC").to_string(),
        }
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

        let overlaps = TopsarMerge::calculate_overlap_regions(&subswaths).unwrap();
        assert_eq!(overlaps.len(), 1);
        
        let overlap = &overlaps[0];
        assert_eq!(overlap.swath1_id, "IW1");
        assert_eq!(overlap.swath2_id, "IW2");
        assert!(overlap.swath1_range_start > 0);
        assert!(overlap.swath2_range_end > 0);
    }
}
