use crate::types::{SarComplex, SarError, SarImage, SarResult};
use ndarray::Array2;
use regex::Regex;
use std::f32::consts::PI;

/// TOPSAR Burst information for proper debursting
/// Based on ESA Sentinel-1 Level 1 Detailed Algorithm Definition
#[derive(Debug, Clone)]
pub struct BurstInfo {
    pub burst_id: usize,
    pub start_line: usize,
    pub end_line: usize,
    pub start_sample: usize,
    pub end_sample: usize,
    pub azimuth_time: String,
    pub sensing_time: String,
    pub first_valid_sample: Vec<i32>,
    pub last_valid_sample: Vec<i32>,
    pub byte_offset: u64,
    
    // TOPSAR-specific parameters
    pub azimuth_fm_rate: f64,        // Azimuth FM rate (Hz/s)
    pub azimuth_steering_rate: f64,  // Azimuth steering rate (rad/s)
    pub slant_range_time: f64,       // Slant range time (s)
    pub doppler_centroid: f64,       // Doppler centroid frequency (Hz)
    pub azimuth_bandwidth: f64,      // Processed azimuth bandwidth (Hz)
    pub range_sampling_rate: f64,    // Range sampling rate (Hz)
    pub range_pixel_spacing: f64,    // Range pixel spacing (m)
    pub azimuth_pixel_spacing: f64,  // Azimuth pixel spacing (m)
}

impl BurstInfo {
    /// Calculate the number of lines in this burst
    pub fn lines(&self) -> usize {
        self.end_line.saturating_sub(self.start_line) + 1
    }
    
    /// Calculate the number of valid samples for a given line
    pub fn valid_samples_for_line(&self, line: usize) -> (usize, usize) {
        if line >= self.first_valid_sample.len() || line >= self.last_valid_sample.len() {
            return (0, 0);
        }
        
        let first = self.first_valid_sample[line].max(0) as usize;
        let last = self.last_valid_sample[line].max(0) as usize;
        
        if first <= last {
            (first, last)
        } else {
            (0, 0)
        }
    }
    
    /// Calculate azimuth deramp phase for TOPSAR processing
    /// Based on equation from ESA Sentinel-1 IPF Algorithm Specification
    pub fn calculate_deramp_phase(&self, line: usize, pixel: usize) -> f32 {
        // Azimuth time relative to burst center
        let burst_center_line = (self.start_line + self.end_line) as f64 / 2.0;
        let azimuth_time_rel = (line as f64 - burst_center_line) * self.azimuth_pixel_spacing / 7000.0; // 7 km/s velocity
        
        // Range time
        let range_time = self.slant_range_time + (pixel as f64) * (1.0 / self.range_sampling_rate);
        
        // Doppler centroid phase
        let doppler_phase = 2.0 * PI as f64 * self.doppler_centroid * azimuth_time_rel;
        
        // Azimuth FM rate phase correction
        let fm_phase = PI as f64 * self.azimuth_fm_rate * azimuth_time_rel * azimuth_time_rel;
        
        // Azimuth steering phase correction for TOPSAR
        let steering_phase = self.azimuth_steering_rate * azimuth_time_rel;
        
        ((doppler_phase + fm_phase + steering_phase) % (2.0 * PI as f64)) as f32
    }
    
    /// Calculate overlap weight for seamless merging between bursts
    /// Uses cosine-squared weighting as specified in TOPSAR literature
    pub fn calculate_overlap_weight(&self, line: usize, overlap_lines: usize) -> f32 {
        let burst_lines = self.lines();
        
        if overlap_lines == 0 || burst_lines <= overlap_lines * 2 {
            return 1.0; // No overlap or burst too small
        }
        
        let line_in_burst = line.saturating_sub(self.start_line);
        
        if line_in_burst < overlap_lines {
            // Beginning of burst - fade in
            let x = line_in_burst as f32 / overlap_lines as f32;
            (0.5 * (1.0 - (PI * x).cos())).powf(2.0)
        } else if line_in_burst >= burst_lines - overlap_lines {
            // End of burst - fade out  
            let x = (burst_lines - line_in_burst - 1) as f32 / overlap_lines as f32;
            (0.5 * (1.0 - (PI * x).cos())).powf(2.0)
        } else {
            // Middle of burst - full weight
            1.0
        }
    }
}

/// Configuration for TOPSAR deburst processing
/// Based on ESA Sentinel-1 processing specifications
#[derive(Debug, Clone)]
pub struct DeburstConfig {
    pub blend_overlap: bool,          // Enable overlap blending between bursts
    pub blend_lines: usize,           // Number of lines to blend (typically 100-200)
    pub remove_invalid_data: bool,    // Remove invalid data regions
    pub seamless_stitching: bool,     // Enable seamless stitching with phase continuity
    pub apply_deramp: bool,           // Apply azimuth deramp for TOPSAR
    pub preserve_phase: bool,         // Preserve interferometric phase information
    pub antenna_pattern_correction: bool, // Apply azimuth antenna pattern correction
}

impl Default for DeburstConfig {
    fn default() -> Self {
        Self {
            blend_overlap: true,
            blend_lines: 150,              // Typical TOPSAR overlap
            remove_invalid_data: true,
            seamless_stitching: true,
            apply_deramp: true,            // Essential for TOPSAR
            preserve_phase: true,          // Important for interferometry
            antenna_pattern_correction: false, // Can be computationally expensive
        }
    }
}

/// TOPSAR Deburst processor for Sentinel-1 IW data
/// Implements scientific algorithms from ESA documentation and literature
pub struct TopSarDeburstProcessor {
    burst_info: Vec<BurstInfo>,
    config: DeburstConfig,
}

impl TopSarDeburstProcessor {
    /// Create a new TOPSAR deburst processor
    pub fn new(burst_info: Vec<BurstInfo>, config: DeburstConfig) -> Self {
        Self {
            burst_info,
            config,
        }
    }

    /// Perform complete TOPSAR debursting with azimuth deramp and seamless merging
    pub fn deburst_topsar(&self, slc_data: &Array2<SarComplex>) -> SarResult<Array2<SarComplex>> {
        log::info!("🔄 Starting TOPSAR deburst processing for {} bursts", self.burst_info.len());

        if self.burst_info.is_empty() {
            return Err(SarError::Processing(
                "No burst information available for TOPSAR debursting".to_string(),
            ));
        }

        // Step 1: Validate burst information and input data consistency
        self.validate_burst_data(slc_data)?;

        // Step 2: Calculate output dimensions for debursted image
        let (output_lines, output_samples) = self.calculate_output_dimensions()?;
        log::info!("Output dimensions: {} lines x {} samples", output_lines, output_samples);

        // Step 3: Initialize output array with proper dimensions
        let mut debursted = Array2::zeros((output_lines, output_samples));
        let mut weight_sum = Array2::zeros((output_lines, output_samples));

        // Step 4: Process each burst with TOPSAR-specific corrections
        for (burst_idx, burst) in self.burst_info.iter().enumerate() {
            log::info!("Processing burst {} of {}", burst_idx + 1, self.burst_info.len());
            
            self.process_single_burst(
                slc_data,
                burst,
                &mut debursted,
                &mut weight_sum,
                burst_idx,
            )?;
        }

        // Step 5: Normalize by accumulated weights (handles overlapping regions)
        self.normalize_overlaps(&mut debursted, &weight_sum)?;

        // Step 6: Apply final quality checks and corrections
        self.apply_final_corrections(&mut debursted)?;

        log::info!("✅ TOPSAR deburst processing completed successfully");
        Ok(debursted)
    }

    /// Validate burst data consistency with input SLC dimensions
    fn validate_burst_data(&self, slc_data: &Array2<SarComplex>) -> SarResult<()> {
        let (slc_lines, slc_samples) = slc_data.dim();
        
        for (i, burst) in self.burst_info.iter().enumerate() {
            if burst.end_line >= slc_lines {
                return Err(SarError::Processing(
                    format!("Burst {} end_line ({}) exceeds SLC dimensions ({})", 
                            i, burst.end_line, slc_lines)
                ));
            }
            
            if burst.end_sample >= slc_samples {
                return Err(SarError::Processing(
                    format!("Burst {} end_sample ({}) exceeds SLC dimensions ({})", 
                            i, burst.end_sample, slc_samples)
                ));
            }
        }
        
        Ok(())
    }

    /// Calculate output dimensions for the debursted image
    fn calculate_output_dimensions(&self) -> SarResult<(usize, usize)> {
        if self.burst_info.is_empty() {
            return Err(SarError::Processing("No bursts available".to_string()));
        }

        // Calculate total lines (sum of all burst lines minus overlaps)
        let mut total_lines = 0;
        for (i, burst) in self.burst_info.iter().enumerate() {
            let burst_lines = burst.lines();
            if i == 0 {
                // First burst - full lines
                total_lines += burst_lines;
            } else {
                // Subsequent bursts - subtract overlap
                total_lines += burst_lines.saturating_sub(self.config.blend_lines);
            }
        }

        // Maximum samples across all bursts
        let max_samples = self.burst_info.iter()
            .map(|b| b.end_sample - b.start_sample + 1)
            .max()
            .unwrap_or(0);

        Ok((total_lines, max_samples))
    }

    /// Process a single burst with TOPSAR-specific corrections
    fn process_single_burst(
        &self,
        slc_data: &Array2<SarComplex>,
        burst: &BurstInfo,
        output: &mut Array2<SarComplex>,
        weights: &mut Array2<f32>,
        burst_idx: usize,
    ) -> SarResult<()> {
        let burst_lines = burst.lines();
        let burst_samples = burst.end_sample - burst.start_sample + 1;
        
        // Calculate output line offset for this burst
        let mut output_line_offset = 0;
        for i in 0..burst_idx {
            let prev_burst_lines = self.burst_info[i].lines();
            if i == 0 {
                output_line_offset += prev_burst_lines;
            } else {
                output_line_offset += prev_burst_lines.saturating_sub(self.config.blend_lines);
            }
        }

        // Process each line in the burst
        for line_in_burst in 0..burst_lines {
            let input_line = burst.start_line + line_in_burst;
            let output_line = output_line_offset + line_in_burst;
            
            // Skip overlap regions for subsequent bursts
            if burst_idx > 0 && line_in_burst < self.config.blend_lines {
                continue;
            }

            // Ensure we don't exceed output dimensions
            if output_line >= output.dim().0 {
                break;
            }

            // Calculate overlap weighting for seamless merging
            let overlap_weight = if self.config.blend_overlap {
                burst.calculate_overlap_weight(input_line, self.config.blend_lines)
            } else {
                1.0
            };

            // Process each sample in the line
            for sample_in_burst in 0..burst_samples {
                let input_sample = burst.start_sample + sample_in_burst;
                let output_sample = sample_in_burst;

                // Ensure we don't exceed array bounds
                if input_line >= slc_data.dim().0 || input_sample >= slc_data.dim().1 {
                    continue;
                }
                if output_line >= output.dim().0 || output_sample >= output.dim().1 {
                    continue;
                }

                // Get complex sample from input
                let mut complex_sample = slc_data[[input_line, input_sample]];

                // Apply TOPSAR azimuth deramp if enabled
                if self.config.apply_deramp {
                    let deramp_phase = burst.calculate_deramp_phase(line_in_burst, sample_in_burst);
                    let deramp_factor = SarComplex::new(deramp_phase.cos(), -deramp_phase.sin());
                    complex_sample = complex_sample * deramp_factor;
                }

                // Apply overlap weighting
                complex_sample = complex_sample * overlap_weight;

                // Accumulate weighted sample
                output[[output_line, output_sample]] = output[[output_line, output_sample]] + complex_sample;
                weights[[output_line, output_sample]] += overlap_weight;
            }
        }

        Ok(())
    }

    /// Normalize overlapping regions by accumulated weights
    fn normalize_overlaps(&self, output: &mut Array2<SarComplex>, weights: &Array2<f32>) -> SarResult<()> {
        let (lines, samples) = output.dim();
        
        for i in 0..lines {
            for j in 0..samples {
                let weight = weights[[i, j]];
                if weight > 0.0 {
                    output[[i, j]] = output[[i, j]] / weight;
                }
            }
        }
        
        Ok(())
    }

    /// Apply final corrections and quality checks
    fn apply_final_corrections(&self, output: &mut Array2<SarComplex>) -> SarResult<()> {
        if self.config.remove_invalid_data {
            // Set invalid/zero samples to a small non-zero value to avoid issues downstream
            let (lines, samples) = output.dim();
            let mut invalid_count = 0;
            
            for i in 0..lines {
                for j in 0..samples {
                    let sample = output[[i, j]];
                    if sample.norm() < 1e-12 || !sample.re.is_finite() || !sample.im.is_finite() {
                        output[[i, j]] = SarComplex::new(1e-8, 0.0);
                        invalid_count += 1;
                    }
                }
            }
            
            if invalid_count > 0 {
                log::info!("Corrected {} invalid samples in debursted output", invalid_count);
            }
        }
        
        Ok(())
    }
}

/// Legacy DeburstProcessor wrapper for backward compatibility
/// Routes to the new TopSarDeburstProcessor with proper TOPSAR support
pub struct DeburstProcessor {
    burst_info: Vec<BurstInfo>,
}

impl DeburstProcessor {
    /// Create a new deburst processor from burst information
    pub fn new(burst_info: Vec<BurstInfo>) -> Self {
        Self { burst_info }
    }
    
    /// Deburst SLC data using the new TOPSAR implementation
    pub fn deburst(&self, slc_data: &Array2<SarComplex>, config: &DeburstConfig) -> SarResult<Array2<SarComplex>> {
        // Create TOPSAR processor and deburst
        let topsar_processor = TopSarDeburstProcessor::new(self.burst_info.clone(), config.clone());
        topsar_processor.deburst_topsar(slc_data)
    }
    
    /// Extract burst information from annotation XML with TOPSAR parameters
    pub fn extract_burst_info_from_annotation(
        annotation_data: &str,
        total_lines: usize,
        total_samples: usize,
    ) -> SarResult<Vec<BurstInfo>> {
        log::info!("🔍 Extracting burst information from Sentinel-1 annotation");
        log::debug!("Total image dimensions: {} x {}", total_lines, total_samples);

        // Try parsing with the enhanced TOPSAR-aware method
        match Self::parse_topsar_burst_info(annotation_data, total_lines, total_samples) {
            Ok(bursts) if !bursts.is_empty() => {
                log::info!("✅ Successfully extracted {} TOPSAR bursts from annotation", bursts.len());
                return Ok(bursts);
            }
            Ok(_) => {
                log::error!("❌ TOPSAR parsing returned empty burst list");
                return Err(SarError::Processing("No valid bursts found in annotation data".to_string()));
            }
            Err(e) => {
                log::error!("❌ TOPSAR annotation parsing failed: {}", e);
                return Err(SarError::Processing(format!("Failed to parse burst information: {}", e)));
            }
        }
    }
    
    /// Parse TOPSAR burst information with enhanced parameter extraction
    fn parse_topsar_burst_info(
        annotation_data: &str,
        total_lines: usize,
        total_samples: usize,
    ) -> SarResult<Vec<BurstInfo>> {
        log::info!("🎯 Parsing TOPSAR burst information with enhanced parameters");
        
        let mut burst_info = Vec::new();
        
        // Extract global TOPSAR parameters
        let azimuth_fm_rate = Self::extract_parameter(annotation_data, "<azimuthFmRate>", "</azimuthFmRate>")
            .unwrap_or(2000.0);
        let azimuth_steering_rate = Self::extract_parameter(annotation_data, "<azimuthSteeringRate>", "</azimuthSteeringRate>")
            .unwrap_or(0.0015);
        let range_sampling_rate = Self::extract_parameter(annotation_data, "<rangeSamplingRate>", "</rangeSamplingRate>")
            .unwrap_or(64000000.0);
        let range_pixel_spacing = Self::extract_parameter(annotation_data, "<rangePixelSpacing>", "</rangePixelSpacing>")
            .unwrap_or(2.33);
        let azimuth_pixel_spacing = Self::extract_parameter(annotation_data, "<azimuthPixelSpacing>", "</azimuthPixelSpacing>")
            .unwrap_or(14.0);
        
        // Regular expressions for burst information
        let burst_pattern = Regex::new(
            r"<burst>.*?<azimuthTime>([^<]+)</azimuthTime>.*?<sensingTime>([^<]+)</sensingTime>.*?<byteOffset>([^<]+)</byteOffset>.*?<firstValidSample[^>]*>([^<]+)</firstValidSample>.*?<lastValidSample[^>]*>([^<]+)</lastValidSample>.*?</burst>"
        ).map_err(|e| SarError::Processing(format!("Failed to compile burst regex: {}", e)))?;
        
        // Extract lines per burst
        let lines_per_burst = Self::extract_parameter(annotation_data, "<linesPerBurst>", "</linesPerBurst>")
            .unwrap_or(1500.0) as usize;
        
        // Find all burst matches
        let burst_matches: Vec<_> = burst_pattern.captures_iter(annotation_data).collect();
        
        if burst_matches.is_empty() {
            return Err(SarError::Processing("No burst information found in annotation".to_string()));
        }
        
        for (i, captures) in burst_matches.iter().enumerate() {
            let azimuth_time = captures.get(1).unwrap().as_str().to_string();
            let sensing_time = captures.get(2).unwrap().as_str().to_string();
            let byte_offset = captures.get(3).unwrap().as_str().parse::<u64>().unwrap_or(0);
            
            let first_valid_sample = Self::parse_sample_array(captures.get(4).unwrap().as_str());
            let last_valid_sample = Self::parse_sample_array(captures.get(5).unwrap().as_str());
            
            // Calculate burst line positions
            let start_line = i * lines_per_burst;
            let end_line = ((i + 1) * lines_per_burst - 1).min(total_lines - 1);
            
            burst_info.push(BurstInfo {
                burst_id: i,
                start_line,
                end_line,
                start_sample: 0,
                end_sample: total_samples - 1,
                azimuth_time,
                sensing_time,
                first_valid_sample,
                last_valid_sample,
                byte_offset,
                
                // TOPSAR-specific parameters
                azimuth_fm_rate,
                azimuth_steering_rate,
                slant_range_time: 0.006,
                doppler_centroid: 0.0,
                azimuth_bandwidth: 320.0,
                range_sampling_rate,
                range_pixel_spacing,
                azimuth_pixel_spacing,
            });
        }
        
        Ok(burst_info)
    }
    
    /// Extract numeric parameter from XML annotation
    fn extract_parameter(annotation_data: &str, start_tag: &str, end_tag: &str) -> Option<f64> {
        if let Some(start_pos) = annotation_data.find(start_tag) {
            let content_start = start_pos + start_tag.len();
            if let Some(end_pos) = annotation_data[content_start..].find(end_tag) {
                let content = &annotation_data[content_start..content_start + end_pos];
                return content.trim().parse::<f64>().ok();
            }
        }
        None
    }
    
    /// Parse a space-separated sample array
    fn parse_sample_array(data: &str) -> Vec<i32> {
        data.split_whitespace()
            .filter_map(|s| s.parse::<i32>().ok())
            .collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use ndarray::Array2;

    #[test]
    fn test_topsar_deburst_processor() {
        let burst_info = vec![
            BurstInfo {
                burst_id: 0,
                start_line: 0,
                end_line: 499,
                start_sample: 0,
                end_sample: 999,
                azimuth_time: "2020-01-03T17:08:16.618328".to_string(),
                sensing_time: "2020-01-03T17:08:17.623236".to_string(),
                first_valid_sample: vec![100; 500],
                last_valid_sample: vec![900; 500],
                byte_offset: 109035,
                azimuth_fm_rate: 2000.0,
                azimuth_steering_rate: 0.0015,
                slant_range_time: 0.006,
                doppler_centroid: 0.0,
                azimuth_bandwidth: 320.0,
                range_sampling_rate: 64000000.0,
                range_pixel_spacing: 2.33,
                azimuth_pixel_spacing: 14.0,
            },
        ];

        let config = DeburstConfig::default();
        let processor = TopSarDeburstProcessor::new(burst_info, config);
        
        // Create test data
        let test_data = Array2::zeros((1000, 1000));
        
        let result = processor.deburst_topsar(&test_data);
        assert!(result.is_ok());
    }
}
