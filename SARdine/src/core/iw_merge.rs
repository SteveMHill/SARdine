use crate::types::{SarComplex, SarError, SarResult, GeoTransform};
use ndarray::Array2;
use std::collections::HashMap;

/// Sub-swath information for IW merge processing
/// Based on Sentinel-1 IW mode specifications
#[derive(Debug, Clone)]
pub struct SubSwathInfo {
    pub swath_id: String,           // "IW1", "IW2", "IW3"
    pub near_range: f64,            // Near range distance (m)
    pub far_range: f64,             // Far range distance (m)
    pub range_pixel_spacing: f64,   // Range pixel spacing (m)
    pub azimuth_pixel_spacing: f64, // Azimuth pixel spacing (m)
    pub incidence_angle_near: f64,  // Incidence angle at near range (degrees)
    pub incidence_angle_far: f64,   // Incidence angle at far range (degrees)
    pub samples_per_line: usize,    // Number of range samples
    pub lines: usize,               // Number of azimuth lines
    pub range_looks: usize,         // Number of range looks applied
    pub azimuth_looks: usize,       // Number of azimuth looks applied
}

impl SubSwathInfo {
    /// Calculate the slant range for a given sample
    pub fn slant_range_for_sample(&self, sample: usize) -> f64 {
        self.near_range + (sample as f64) * self.range_pixel_spacing
    }
    
    /// Calculate the incidence angle for a given sample (linear interpolation)
    pub fn incidence_angle_for_sample(&self, sample: usize) -> f64 {
        let fraction = (sample as f64) / (self.samples_per_line as f64 - 1.0);
        self.incidence_angle_near + fraction * (self.incidence_angle_far - self.incidence_angle_near)
    }
    
    /// Check if this sub-swath overlaps with another in range
    pub fn overlaps_with(&self, other: &SubSwathInfo) -> bool {
        !(self.far_range < other.near_range || other.far_range < self.near_range)
    }
    
    /// Calculate overlap region with another sub-swath
    pub fn overlap_region(&self, other: &SubSwathInfo) -> Option<(f64, f64)> {
        if !self.overlaps_with(other) {
            return None;
        }
        
        let overlap_near = self.near_range.max(other.near_range);
        let overlap_far = self.far_range.min(other.far_range);
        
        Some((overlap_near, overlap_far))
    }
}

/// Configuration for IW merge processing
#[derive(Debug, Clone)]
pub struct IwMergeConfig {
    pub blend_overlaps: bool,        // Blend overlapping regions between sub-swaths
    pub blend_width_meters: f64,     // Width of blending region in meters
    pub normalize_by_incidence: bool, // Normalize backscatter by incidence angle
    pub preserve_radiometry: bool,   // Preserve absolute radiometric calibration
    pub quality_check: bool,         // Perform quality checks on merged data
    pub output_spacing: f64,         // Output pixel spacing in meters (0 = auto)
}

impl Default for IwMergeConfig {
    fn default() -> Self {
        Self {
            blend_overlaps: true,
            blend_width_meters: 200.0,      // 200m blending width
            normalize_by_incidence: false,  // Preserve original radiometry
            preserve_radiometry: true,
            quality_check: true,
            output_spacing: 0.0,            // Auto-determine from input
        }
    }
}

/// IW Merge processor for combining Sentinel-1 sub-swaths
/// Reference: ESA Sentinel-1 User Handbook, Section 2.3.3
pub struct IwMergeProcessor {
    sub_swaths: HashMap<String, SubSwathInfo>,
    config: IwMergeConfig,
}

impl IwMergeProcessor {
    /// Create a new IW merge processor
    pub fn new(sub_swaths: Vec<SubSwathInfo>, config: IwMergeConfig) -> Self {
        let mut swath_map = HashMap::new();
        for swath in sub_swaths {
            swath_map.insert(swath.swath_id.clone(), swath);
        }
        
        Self {
            sub_swaths: swath_map,
            config,
        }
    }
    
    /// Merge multiple sub-swath images into a single wide-swath image
    /// Implements seamless merging with proper radiometric preservation
    pub fn merge_subswaths(
        &self,
        swath_images: &HashMap<String, Array2<SarComplex>>,
        swath_transforms: &HashMap<String, GeoTransform>,
    ) -> SarResult<(Array2<SarComplex>, GeoTransform)> {
        log::info!("🔗 Starting IW merge processing for {} sub-swaths", swath_images.len());
        
        if swath_images.is_empty() {
            return Err(SarError::Processing("No sub-swath images provided".to_string()));
        }
        
        // Step 1: Determine merge order and overlap regions
        let merge_order = self.determine_merge_order()?;
        log::info!("Merge order: {:?}", merge_order);
        
        // Step 2: Calculate output dimensions and coordinate system
        let (output_lines, output_samples, output_transform) = 
            self.calculate_output_grid(&merge_order, swath_transforms)?;
        
        // Step 3: Initialize output arrays
        let mut merged_image = Array2::zeros((output_lines, output_samples));
        let mut weight_array = Array2::zeros((output_lines, output_samples));
        
        // Step 4: Merge each sub-swath with proper blending
        for swath_id in &merge_order {
            if let (Some(image), Some(transform)) = (
                swath_images.get(swath_id),
                swath_transforms.get(swath_id),
            ) {
                self.merge_single_subswath(
                    image,
                    transform,
                    swath_id,
                    &mut merged_image,
                    &mut weight_array,
                    &output_transform,
                )?;
            }
        }
        
        // Step 5: Normalize overlapping regions
        self.normalize_merged_image(&mut merged_image, &weight_array)?;
        
        // Step 6: Apply final quality checks
        if self.config.quality_check {
            self.apply_quality_checks(&mut merged_image)?;
        }
        
        log::info!("✅ IW merge processing completed successfully");
        Ok((merged_image, output_transform))
    }
    
    /// Determine the optimal order for merging sub-swaths (typically IW1, IW2, IW3)
    fn determine_merge_order(&self) -> SarResult<Vec<String>> {
        let mut swath_ids: Vec<String> = self.sub_swaths.keys().cloned().collect();
        
        // Sort by near range (ascending order)
        swath_ids.sort_by(|a, b| {
            let range_a = self.sub_swaths[a].near_range;
            let range_b = self.sub_swaths[b].near_range;
            range_a.partial_cmp(&range_b).unwrap_or(std::cmp::Ordering::Equal)
        });
        
        log::debug!("Sub-swath merge order determined by range: {:?}", swath_ids);
        Ok(swath_ids)
    }
    
    /// Calculate output grid dimensions and coordinate transform
    fn calculate_output_grid(
        &self,
        merge_order: &[String],
        transforms: &HashMap<String, GeoTransform>,
    ) -> SarResult<(usize, usize, GeoTransform)> {
        if merge_order.is_empty() {
            return Err(SarError::Processing("No sub-swaths to merge".to_string()));
        }
        
        // Find overall bounding box in geographic coordinates
        let mut min_x = f64::INFINITY;
        let mut max_x = f64::NEG_INFINITY;
        let mut min_y = f64::INFINITY;
        let mut max_y = f64::NEG_INFINITY;
        
        for swath_id in merge_order {
            if let (Some(swath_info), Some(transform)) = (
                self.sub_swaths.get(swath_id),
                transforms.get(swath_id),
            ) {
                // Calculate corner coordinates
                let corners = [
                    (0.0, 0.0),
                    (swath_info.samples_per_line as f64, 0.0),
                    (0.0, swath_info.lines as f64),
                    (swath_info.samples_per_line as f64, swath_info.lines as f64),
                ];
                
                for (col, row) in corners.iter() {
                    let x = transform.top_left_x + col * transform.pixel_width;
                    let y = transform.top_left_y + row * transform.pixel_height;
                    
                    min_x = min_x.min(x);
                    max_x = max_x.max(x);
                    min_y = min_y.min(y);
                    max_y = max_y.max(y);
                }
            }
        }
        
        // Determine output pixel spacing
        let pixel_spacing = if self.config.output_spacing > 0.0 {
            self.config.output_spacing
        } else {
            // Use the finest resolution among input sub-swaths
            merge_order
                .iter()
                .filter_map(|id| transforms.get(id))
                .map(|t| t.pixel_width.abs())
                .fold(f64::INFINITY, f64::min)
        };
        
        // Calculate output dimensions
        let output_samples = ((max_x - min_x) / pixel_spacing).ceil() as usize;
        let output_lines = ((max_y - min_y) / pixel_spacing.abs()).ceil() as usize;
        
        // Create output transform
        let output_transform = GeoTransform {
            top_left_x: min_x,
            pixel_width: pixel_spacing,
            rotation_x: 0.0,
            top_left_y: max_y,
            rotation_y: 0.0,
            pixel_height: -pixel_spacing, // Negative for north-up
        };
        
        log::info!("Output grid: {} x {} pixels, spacing: {:.2}m", 
                  output_samples, output_lines, pixel_spacing);
        
        Ok((output_lines, output_samples, output_transform))
    }
    
    /// Merge a single sub-swath into the output grid
    fn merge_single_subswath(
        &self,
        image: &Array2<SarComplex>,
        transform: &GeoTransform,
        swath_id: &str,
        output: &mut Array2<SarComplex>,
        weights: &mut Array2<f32>,
        output_transform: &GeoTransform,
    ) -> SarResult<()> {
        let swath_info = self.sub_swaths.get(swath_id)
            .ok_or_else(|| SarError::Processing(format!("Sub-swath {} info not found", swath_id)))?;
        
        let (input_lines, input_samples) = image.dim();
        let (output_lines, output_samples) = output.dim();
        
        log::info!("Merging sub-swath {} ({} x {} pixels)", swath_id, input_samples, input_lines);
        
        // Process each pixel in the input sub-swath
        for input_row in 0..input_lines {
            for input_col in 0..input_samples {
                // Convert input pixel to geographic coordinates
                let geo_x = transform.top_left_x + (input_col as f64) * transform.pixel_width;
                let geo_y = transform.top_left_y + (input_row as f64) * transform.pixel_height;
                
                // Convert geographic coordinates to output pixel coordinates
                let output_col = ((geo_x - output_transform.top_left_x) / output_transform.pixel_width) as i32;
                let output_row = ((geo_y - output_transform.top_left_y) / output_transform.pixel_height) as i32;
                
                // Check if output pixel is within bounds
                if output_row >= 0 && output_row < output_lines as i32 &&
                   output_col >= 0 && output_col < output_samples as i32 {
                    let out_row = output_row as usize;
                    let out_col = output_col as usize;
                    
                    // Get complex sample value
                    let sample_value = image[[input_row, input_col]];
                    
                    // Calculate blending weight for overlap regions
                    let blend_weight = self.calculate_blend_weight(
                        swath_info,
                        input_col,
                        geo_x,
                    );
                    
                    // Apply radiometric normalization if requested
                    let normalized_value = if self.config.normalize_by_incidence {
                        let inc_angle = swath_info.incidence_angle_for_sample(input_col);
                        sample_value / inc_angle.to_radians().sin() as f32
                    } else {
                        sample_value
                    };
                    
                    // Accumulate weighted sample
                    output[[out_row, out_col]] = output[[out_row, out_col]] + normalized_value * blend_weight;
                    weights[[out_row, out_col]] += blend_weight;
                }
            }
        }
        
        Ok(())
    }
    
    /// Calculate blending weight for overlap regions
    fn calculate_blend_weight(&self, swath_info: &SubSwathInfo, sample: usize, range_coord: f64) -> f32 {
        if !self.config.blend_overlaps {
            return 1.0;
        }
        
        let slant_range = swath_info.slant_range_for_sample(sample);
        let blend_width = self.config.blend_width_meters;
        
        // Calculate distance from swath edges
        let dist_from_near = slant_range - swath_info.near_range;
        let dist_from_far = swath_info.far_range - slant_range;
        
        // Apply cosine-squared tapering near edges
        let edge_distance = dist_from_near.min(dist_from_far);
        
        if edge_distance < blend_width {
            let taper_factor = edge_distance / blend_width;
            (0.5 * (1.0 - (std::f64::consts::PI * taper_factor).cos())) as f32
        } else {
            1.0
        }
    }
    
    /// Normalize the merged image by accumulated weights
    fn normalize_merged_image(
        &self,
        image: &mut Array2<SarComplex>,
        weights: &Array2<f32>,
    ) -> SarResult<()> {
        let (lines, samples) = image.dim();
        let mut normalized_pixels = 0;
        
        for i in 0..lines {
            for j in 0..samples {
                let weight = weights[[i, j]];
                if weight > 0.0 {
                    image[[i, j]] = image[[i, j]] / weight;
                    normalized_pixels += 1;
                }
            }
        }
        
        log::info!("Normalized {} pixels in merged image", normalized_pixels);
        Ok(())
    }
    
    /// Apply quality checks to the merged image
    fn apply_quality_checks(&self, image: &mut Array2<SarComplex>) -> SarResult<()> {
        let (lines, samples) = image.dim();
        let mut invalid_count = 0;
        let total_pixels = lines * samples;
        
        for i in 0..lines {
            for j in 0..samples {
                let sample = image[[i, j]];
                if !sample.re.is_finite() || !sample.im.is_finite() || sample.norm() < 1e-12 {
                    image[[i, j]] = SarComplex::new(1e-8, 0.0);
                    invalid_count += 1;
                }
            }
        }
        
        let invalid_percentage = (invalid_count as f64 / total_pixels as f64) * 100.0;
        log::info!("Quality check: {:.2}% invalid pixels corrected", invalid_percentage);
        
        if invalid_percentage > 50.0 {
            log::warn!("High percentage of invalid pixels ({:.1}%) - check input data quality", invalid_percentage);
        }
        
        Ok(())
    }
}
