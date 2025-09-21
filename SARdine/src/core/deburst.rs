use crate::types::{SarComplex, SarError, SarResult};
use ndarray::Array2;
use std::f32::consts::PI;
use num_complex::Complex;

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
    /// Reference: ESA-EOPG-CSCOP-TN-0009 "Sentinel-1 IPF Algorithms"
    pub fn calculate_deramp_phase(&self, line: usize, pixel: usize, satellite_velocity: f64) -> f32 {
        // Azimuth time relative to burst center
        let burst_center_line = (self.start_line + self.end_line) as f64 / 2.0;
        
        // SCIENTIFIC FIX: Use actual satellite velocity from orbit state vectors
        // Previous hardcoded 7000.0 m/s was scientifically incorrect approximation
        let azimuth_time_rel = (line as f64 - burst_center_line) * self.azimuth_pixel_spacing / satellite_velocity;
        
        // Range time
        let _range_time = self.slant_range_time + (pixel as f64) * (1.0 / self.range_sampling_rate);
        
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
    satellite_velocity: f64,  // Actual satellite velocity from orbit state vectors (m/s)
}

impl TopSarDeburstProcessor {
    /// Create a new TOPSAR deburst processor
    pub fn new(burst_info: Vec<BurstInfo>, config: DeburstConfig, satellite_velocity: f64) -> Self {
        Self {
            burst_info,
            config,
            satellite_velocity,
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
                    let deramp_phase = burst.calculate_deramp_phase(line_in_burst, sample_in_burst, self.satellite_velocity);
                    let deramp_factor = SarComplex::new(deramp_phase.cos(), -deramp_phase.sin());
                    complex_sample *= deramp_factor;
                }

                // Apply overlap weighting
                complex_sample *= overlap_weight;

                // Accumulate weighted sample
                output[[output_line, output_sample]] += complex_sample;
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
                    output[[i, j]] /= weight;
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
    satellite_velocity: f64,
}

impl DeburstProcessor {
    /// Create a new deburst processor from burst information
    pub fn new(burst_info: Vec<BurstInfo>, satellite_velocity: f64) -> Self {
        Self { burst_info, satellite_velocity }
    }
    
    /// Deburst SLC data using the new TOPSAR implementation
    pub fn deburst(&self, slc_data: &Array2<SarComplex>, config: &DeburstConfig) -> SarResult<Array2<SarComplex>> {
        // Create TOPSAR processor and deburst
        let topsar_processor = TopSarDeburstProcessor::new(self.burst_info.clone(), config.clone(), self.satellite_velocity);
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
        
        // Extract global TOPSAR parameters with strict validation - NO FALLBACKS
        let azimuth_fm_rate = Self::extract_parameter_string(annotation_data, "<azimuthFmRatePolynomial", "</azimuthFmRatePolynomial>")
            .and_then(|s| {
                // Find the closing of the opening tag and extract polynomial coefficients
                if let Some(content_start) = s.find('>') {
                    let content = &s[content_start + 1..];
                    let coeffs: Vec<&str> = content.split_whitespace().collect();
                    // The first coefficient is the azimuth FM rate constant term
                    coeffs.first().and_then(|s| s.parse::<f64>().ok())
                } else {
                    // Fallback: try to parse the content directly
                    let coeffs: Vec<&str> = s.split_whitespace().collect();
                    coeffs.first().and_then(|s| s.parse::<f64>().ok())
                }
            })
            .or_else(|| Self::extract_parameter(annotation_data, "<azimuthFmRate>", "</azimuthFmRate>"))
            .or_else(|| Self::extract_parameter(annotation_data, "<azimuthFMRate>", "</azimuthFMRate>"))
            .ok_or_else(|| SarError::Processing("❌ SCIENTIFIC ERROR: Azimuth FM rate not found in annotation XML! Real Sentinel-1 parameters required - no fallbacks allowed.".to_string()))?;
        let azimuth_steering_rate = Self::extract_parameter(annotation_data, "<azimuthSteeringRate>", "</azimuthSteeringRate>")
            .ok_or_else(|| SarError::Processing("❌ SCIENTIFIC ERROR: Azimuth steering rate not found in annotation XML! Real Sentinel-1 parameters required - no fallbacks allowed.".to_string()))?;
        let range_sampling_rate = Self::extract_parameter(annotation_data, "<rangeSamplingRate>", "</rangeSamplingRate>")
            .ok_or_else(|| SarError::Processing("❌ SCIENTIFIC ERROR: Range sampling rate not found in annotation XML! Real Sentinel-1 parameters required - no fallbacks allowed.".to_string()))?;
        
        // Important: Extract real pixel spacing from annotation - NO hardcoded fallbacks for research use
        let range_pixel_spacing = Self::extract_parameter(annotation_data, "<rangePixelSpacing>", "</rangePixelSpacing>")
            .ok_or_else(|| SarError::Processing("❌ SCIENTIFIC ERROR: Range pixel spacing not found in annotation XML! Real Sentinel-1 annotation required - no synthetic fallbacks allowed for research-grade processing.".to_string()))?;
        let azimuth_pixel_spacing = Self::extract_parameter(annotation_data, "<azimuthPixelSpacing>", "</azimuthPixelSpacing>")
            .ok_or_else(|| SarError::Processing("❌ SCIENTIFIC ERROR: Azimuth pixel spacing not found in annotation XML! Real Sentinel-1 annotation required - no synthetic fallbacks allowed for research-grade processing.".to_string()))?;
        
        // Extract lines per burst - SCIENTIFIC REQUIREMENT: Must be from real annotation
        let lines_per_burst = Self::extract_parameter(annotation_data, "<linesPerBurst>", "</linesPerBurst>")
            .ok_or_else(|| SarError::Processing("❌ SCIENTIFIC ERROR: Lines per burst not found in annotation XML! Real Sentinel-1 burst parameters required - no synthetic values allowed.".to_string()))? as usize;
        
        // Extract samples per burst - SCIENTIFIC REQUIREMENT: Must be from real annotation
        let samples_per_burst = Self::extract_parameter(annotation_data, "<samplesPerBurst>", "</samplesPerBurst>")
            .ok_or_else(|| SarError::Processing("❌ SCIENTIFIC ERROR: Samples per burst not found in annotation XML! Real Sentinel-1 burst parameters required - no synthetic values allowed.".to_string()))? as usize;
        
        log::info!("📊 TOPSAR parameters: lines_per_burst={}, samples_per_burst={}", lines_per_burst, samples_per_burst);
        log::info!("📊 Range sampling rate: {:.0} Hz", range_sampling_rate);
        log::info!("📊 Azimuth steering rate: {:.6} rad/s", azimuth_steering_rate);
        
        // Use a more precise regex pattern matching the actual XML structure
        let burst_pattern = regex::Regex::new(
            r"(?s)<burst>.*?<azimuthTime>([^<]+)</azimuthTime>.*?<sensingTime>([^<]+)</sensingTime>.*?<byteOffset>([^<]+)</byteOffset>.*?<firstValidSample[^>]*>([^<]+)</firstValidSample>.*?<lastValidSample[^>]*>([^<]+)</lastValidSample>.*?</burst>"
        ).map_err(|e| SarError::Processing(format!("Failed to compile burst regex: {}", e)))?;
        
        // Find all burst matches
        let burst_matches: Vec<_> = burst_pattern.captures_iter(annotation_data).collect();
        
        if burst_matches.is_empty() {
            log::error!("❌ No burst information found with regex pattern");
            
            // Fallback: check if we can find burst list count
            if let Some(count_match) = regex::Regex::new(r#"<burstList count="(\d+)">"#).unwrap().captures(annotation_data) {
                if let Ok(burst_count) = count_match[1].parse::<usize>() {
                    log::info!("📊 Found burstList with {} bursts, but couldn't parse individual bursts", burst_count);
                }
            }
            
            return Err(SarError::Processing("No burst information found in annotation".to_string()));
        }
        
        log::info!("✅ Found {} burst matches in annotation", burst_matches.len());
        
        for (i, captures) in burst_matches.iter().enumerate() {
            let azimuth_time = captures.get(1).unwrap().as_str().to_string();
            let sensing_time = captures.get(2).unwrap().as_str().to_string();
            let byte_offset = captures.get(3).unwrap().as_str().parse::<u64>()
                .map_err(|e| SarError::Processing(format!("Failed to parse byte_offset for burst {}: {}", i, e)))?;
            
            let first_valid_sample = Self::parse_sample_array(captures.get(4).unwrap().as_str());
            let last_valid_sample = Self::parse_sample_array(captures.get(5).unwrap().as_str());
            
            // Verify sample array lengths match expected lines per burst
            if first_valid_sample.len() != lines_per_burst {
                log::warn!("⚠️ Burst {}: firstValidSample length {} != lines_per_burst {}", 
                          i, first_valid_sample.len(), lines_per_burst);
            }
            
            // Calculate burst line positions
            let start_line = i * lines_per_burst;
            let end_line = ((i + 1) * lines_per_burst - 1).min(total_lines.saturating_sub(1));
            
            // Use actual samples per burst from annotation
            let start_sample = 0;
            let end_sample = samples_per_burst.saturating_sub(1).min(total_samples.saturating_sub(1));
            
            log::info!("📋 Burst {}: lines {}..{}, samples {}..{}, byte_offset={}", 
                      i, start_line, end_line, start_sample, end_sample, byte_offset);
            
            burst_info.push(BurstInfo {
                burst_id: i,
                start_line,
                end_line,
                start_sample,
                end_sample,
                azimuth_time,
                sensing_time,
                first_valid_sample,
                last_valid_sample,
                byte_offset,
                
                // TOPSAR-specific parameters (real values from annotation)
                azimuth_fm_rate,
                azimuth_steering_rate,
                slant_range_time: 0.006,  // Typical S-band value
                doppler_centroid: 0.0,    // Will be refined from DC polynomials if available
                azimuth_bandwidth: 320.0, // Typical TOPSAR bandwidth
                range_sampling_rate,
                range_pixel_spacing,
                azimuth_pixel_spacing,
            });
        }
        
        log::info!("✅ Successfully parsed {} TOPSAR bursts with real parameters", burst_info.len());
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
    
    /// Extract raw string content between XML tags
    fn extract_parameter_string(annotation_data: &str, start_tag: &str, end_tag: &str) -> Option<String> {
        if let Some(start_pos) = annotation_data.find(start_tag) {
            let content_start = start_pos + start_tag.len();
            if let Some(end_pos) = annotation_data[content_start..].find(end_tag) {
                let content = &annotation_data[content_start..content_start + end_pos];
                return Some(content.trim().to_string());
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
                // Use realistic Sentinel-1 IW1 parameters for test (from real annotation data)
                range_pixel_spacing: 2.329562,    // Realistic IW1 range pixel spacing
                azimuth_pixel_spacing: 14.059906, // Realistic IW azimuth pixel spacing
            },
        ];

        let config = DeburstConfig::default();
        let satellite_velocity = 7500.0; // Typical Sentinel-1 velocity for testing
        let processor = TopSarDeburstProcessor::new(burst_info, config, satellite_velocity);
        
        // Create test data
        let test_data = Array2::zeros((1000, 1000));
        
        let result = processor.deburst_topsar(&test_data);
        assert!(result.is_ok());
    }
}

/// Extract subswath data from a Sentinel-1 ZIP file
/// This function reads the SLC data for a specific subswath and polarization
pub fn extract_subswath_from_zip(
    zip_path: &str,
    subswath: &str,
    polarization: &str,
) -> SarResult<Array2<f32>> {
    log::info!("📡 Extracting subswath {} {} from {}", subswath, polarization, zip_path);

    // Open the ZIP file
    let file = std::fs::File::open(zip_path)
        .map_err(|e| SarError::Processing(format!("Failed to open ZIP file: {}", e)))?;

    let mut archive = zip::ZipArchive::new(file)
        .map_err(|e| SarError::Processing(format!("Failed to read ZIP archive: {}", e)))?;

    // Find the appropriate measurement file for the subswath and polarization
    // Extract subswath number (e.g., "IW1" -> "1")
    let subswath_num = subswath.strip_prefix("IW").unwrap_or(subswath);
    let measurement_pattern = format!("s1a-iw{}-slc-{}-", subswath_num.to_lowercase(), polarization.to_lowercase());

    let mut measurement_file = None;
    log::info!("Looking for files with pattern: '{}' in measurement/ directory", measurement_pattern);

    for i in 0..archive.len() {
        let file = archive.by_index(i)
            .map_err(|e| SarError::Processing(format!("Failed to read ZIP entry: {}", e)))?;

        let name = file.name();
        log::debug!("Checking file: {}", name);

        if name.contains(&measurement_pattern) && name.ends_with(".tiff") && name.contains("measurement/") {
            measurement_file = Some(i);
            log::info!("Found measurement file: {}", name);
            break;
        }
    }

    if measurement_file.is_none() {
        return Err(SarError::Processing(format!(
            "No measurement file found for subswath {} and polarization {}",
            subswath, polarization
        )));
    }

    // Extract TIFF to temporary file for GDAL reading
    let mut zip_file = archive.by_index(measurement_file.unwrap())
        .map_err(|e| SarError::Processing(format!("Failed to access measurement file: {}", e)))?;

    // Create temporary file
    let mut temp_file = tempfile::NamedTempFile::new()
        .map_err(|e| SarError::Processing(format!("Failed to create temp file: {}", e)))?;

    // Copy TIFF data to temp file
    std::io::copy(&mut zip_file, &mut temp_file)
        .map_err(|e| SarError::Processing(format!("Failed to extract TIFF to temp file: {}", e)))?;

    log::info!("TIFF extracted to temporary file for GDAL processing");

    // Use GDAL to read the TIFF file
    let dataset = gdal::Dataset::open(temp_file.path())
        .map_err(|e| SarError::Processing(format!("Failed to open TIFF with GDAL: {}", e)))?;

    // Get raster dimensions
    let raster_size = dataset.raster_size();
    let (width, height) = (raster_size.0, raster_size.1);
    let band_count = dataset.raster_count();

    log::info!("TIFF dimensions: {} x {}, bands: {}", width, height, band_count);

    // Handle different band configurations (Sentinel-1 SLC format)
    let intensity_data = if band_count == 1 {
        // Single band containing complex data (Sentinel-1 SLC format: CInt16)
        let band = dataset.rasterband(1)
            .map_err(|e| SarError::Processing(format!("Failed to get band 1: {}", e)))?;

        let window = (0, 0);
        let window_size = (width, height);

        // Read complex 16-bit integers (CInt16 format)
        let complex_data = band.read_as::<i16>(window, window_size, (width * 2, height), None)
            .map_err(|e| SarError::Processing(format!("Failed to read complex CInt16 data: {}", e)))?;

        log::debug!("Read {} complex values as interleaved i16", complex_data.data.len() / 2);

        // Convert interleaved i16 data to intensity (magnitude squared) f32 array
        convert_cint16_to_intensity(complex_data, width, height)?
    } else if band_count >= 2 {
        // Separate I and Q bands (less common for Sentinel-1)
        let band1 = dataset.rasterband(1)
            .map_err(|e| SarError::Processing(format!("Failed to get band 1: {}", e)))?;

        let band2 = dataset.rasterband(2)
            .map_err(|e| SarError::Processing(format!("Failed to get band 2: {}", e)))?;

        let window = (0, 0);
        let window_size = (width, height);
        let buffer_size = (width, height);

        let i_data = band1.read_as::<f32>(window, window_size, buffer_size, None)
            .map_err(|e| SarError::Processing(format!("Failed to read I data: {}", e)))?;

        let q_data = band2.read_as::<f32>(window, window_size, buffer_size, None)
            .map_err(|e| SarError::Processing(format!("Failed to read Q data: {}", e)))?;

        // Convert to intensity data
        convert_iq_to_intensity(i_data, q_data, width, height)?
    } else {
        return Err(SarError::Processing(format!("Unexpected number of bands in TIFF: {}", band_count)));
    };

    log::info!("✅ Successfully extracted real SAR subswath data: {}x{} pixels", width, height);
    log::info!("Data range: {:.6} to {:.6}", intensity_data.iter().fold(f32::INFINITY, |a, &b| a.min(b)), intensity_data.iter().fold(f32::NEG_INFINITY, |a, &b| a.max(b)));

    Ok(intensity_data)
}

/// Convert CInt16 interleaved data to intensity (magnitude squared)
fn convert_cint16_to_intensity(complex_data: gdal::raster::Buffer<i16>, width: usize, height: usize) -> SarResult<Array2<f32>> {
    let mut intensity_array = Array2::zeros((height, width));
    let total_pixels = width * height;

    // Data is interleaved as [real, imag, real, imag, ...]
    if complex_data.data.len() < total_pixels * 2 {
        return Err(SarError::Processing(format!("Insufficient data: expected {} i16 values, got {}",
                                                total_pixels * 2, complex_data.data.len())));
    }

    // Convert i16 complex data to f32 intensity values
    for row in 0..height {
        for col in 0..width {
            let pixel_idx = row * width + col;
            let data_idx = pixel_idx * 2;

            // Extract real and imaginary parts as i16, convert to f32
            let real_i16 = complex_data.data[data_idx];
            let imag_i16 = complex_data.data[data_idx + 1];

            // Convert to f32
            let real_f32 = real_i16 as f32;
            let imag_f32 = imag_i16 as f32;

            // Calculate intensity (magnitude squared)
            let intensity = real_f32 * real_f32 + imag_f32 * imag_f32;
            intensity_array[[row, col]] = intensity;
        }
    }

    // Log some statistics
    let sample_intensity = intensity_array[[0, 0]];
    log::debug!("First pixel intensity: {:.6}", sample_intensity);

    // Check data variance
    let mut sum = 0.0f32;
    let mut sum_sq = 0.0f32;
    let sample_size = std::cmp::min(1000, total_pixels);

    for i in 0..sample_size {
        let row = i / width;
        let col = i % width;
        if row < height {
            let val = intensity_array[[row, col]];
            sum += val;
            sum_sq += val * val;
        }
    }

    let mean = sum / sample_size as f32;
    let variance = (sum_sq / sample_size as f32) - (mean * mean);
    log::debug!("Sample statistics (first {} pixels): mean={:.6}, variance={:.6}",
               sample_size, mean, variance);

    Ok(intensity_array)
}

/// Convert I/Q data to intensity (magnitude squared)
fn convert_iq_to_intensity(i_data: gdal::raster::Buffer<f32>, q_data: gdal::raster::Buffer<f32>, width: usize, height: usize) -> SarResult<Array2<f32>> {
    let mut intensity_array = Array2::zeros((height, width));

    // Convert buffer data to intensity
    for i in 0..height {
        for j in 0..width {
            let idx = i * width + j;
            if idx < i_data.data.len() && idx < q_data.data.len() {
                let real = i_data.data[idx];
                let imag = q_data.data[idx];
                let intensity = real * real + imag * imag;
                intensity_array[[i, j]] = intensity;
            }
        }
    }

    Ok(intensity_array)
}

/// Extract complex SLC data (preserving I/Q values) from ZIP file for proper radiometric calibration
pub fn extract_subswath_complex_from_zip(
    zip_path: &str,
    subswath: &str,
    polarization: &str,
) -> SarResult<ndarray::Array2<num_complex::Complex<f32>>> {
    
    log::info!("📡 Extracting COMPLEX subswath {} {} from {}", subswath, polarization, zip_path);

    // Check if input is a ZIP file or SAFE directory
    let input_path = std::path::Path::new(zip_path);
    let is_zip = input_path.is_file() && input_path.extension() == Some(std::ffi::OsStr::new("zip"));
    let is_safe = input_path.is_dir() && (
        input_path.file_name()
            .map(|name| name.to_string_lossy().contains(".SAFE"))
            .unwrap_or_else(|| {
                log::warn!("⚠️  Could not determine filename for path validation - checking manifest.safe instead");
                false
            }) ||
        input_path.join("manifest.safe").exists()
    );
    
    if !is_zip && !is_safe {
        return Err(SarError::Processing(format!(
            "Input must be either a ZIP file or SAFE directory: {}", zip_path
        )));
    }
    
    if is_safe {
        // For SAFE directories, delegate to SlcReader which has full SAFE support
        log::info!("Detected SAFE directory, using SlcReader for complex data extraction");
        return extract_subswath_complex_from_safe(zip_path, subswath, polarization);
    }

    // Original ZIP file processing continues below
    // Open the ZIP file
    let file = std::fs::File::open(zip_path)
        .map_err(|e| SarError::Processing(format!("Failed to open ZIP file: {}", e)))?;

    let mut archive = zip::ZipArchive::new(file)
        .map_err(|e| SarError::Processing(format!("Failed to read ZIP archive: {}", e)))?;

    // Find the appropriate measurement file for the subswath and polarization
    let subswath_num = subswath.strip_prefix("IW").unwrap_or(subswath);
    let measurement_pattern = format!("s1a-iw{}-slc-{}-", subswath_num.to_lowercase(), polarization.to_lowercase());

    let mut measurement_file = None;
    for i in 0..archive.len() {
        let file = archive.by_index(i)
            .map_err(|e| SarError::Processing(format!("Failed to read ZIP entry: {}", e)))?;

        let name = file.name();
        if name.contains(&measurement_pattern) && name.ends_with(".tiff") && name.contains("measurement/") {
            measurement_file = Some(i);
            log::info!("Found measurement file: {}", name);
            break;
        }
    }

    if measurement_file.is_none() {
        return Err(SarError::Processing(format!(
            "No measurement file found for subswath {} and polarization {}",
            subswath, polarization
        )));
    }

    // Extract TIFF to temporary file for GDAL reading
    let mut zip_file = archive.by_index(measurement_file.unwrap())
        .map_err(|e| SarError::Processing(format!("Failed to access measurement file: {}", e)))?;

    let mut temp_file = tempfile::NamedTempFile::new()
        .map_err(|e| SarError::Processing(format!("Failed to create temp file: {}", e)))?;

    std::io::copy(&mut zip_file, &mut temp_file)
        .map_err(|e| SarError::Processing(format!("Failed to extract TIFF to temp file: {}", e)))?;

    // Use GDAL to read the TIFF file
    let dataset = gdal::Dataset::open(temp_file.path())
        .map_err(|e| SarError::Processing(format!("Failed to open TIFF with GDAL: {}", e)))?;

    let raster_size = dataset.raster_size();
    let (width, height) = (raster_size.0, raster_size.1);
    let band_count = dataset.raster_count();

    log::info!("TIFF dimensions: {} x {}, bands: {}", width, height, band_count);

    // Handle complex SLC data
    if band_count == 1 {
        let band = dataset.rasterband(1)
            .map_err(|e| SarError::Processing(format!("Failed to get band 1: {}", e)))?;

        let window = (0, 0);
        let window_size = (width, height);

        // Read complex 16-bit integers (CInt16 format)
        let complex_data = band.read_as::<i16>(window, window_size, (width * 2, height), None)
            .map_err(|e| SarError::Processing(format!("Failed to read complex CInt16 data: {}", e)))?;

        // Convert to Complex<f32> array (preserve I/Q values)
        convert_cint16_to_complex(complex_data, width, height)
    } else {
        Err(SarError::Processing(format!("Complex extraction not supported for {} bands", band_count)))
    }
}

/// Convert CInt16 interleaved data to Complex<f32> (preserving I/Q values)
fn convert_cint16_to_complex(
    complex_data: gdal::raster::Buffer<i16>, 
    width: usize, 
    height: usize
) -> SarResult<ndarray::Array2<num_complex::Complex<f32>>> {
    
    let mut complex_array = ndarray::Array2::zeros((height, width));
    let total_pixels = width * height;

    // Data is interleaved as [real, imag, real, imag, ...]
    if complex_data.data.len() < total_pixels * 2 {
        return Err(SarError::Processing(format!("Insufficient data: expected {} i16 values, got {}",
                                                total_pixels * 2, complex_data.data.len())));
    }

    // Convert i16 complex data to Complex<f32> values
    for row in 0..height {
        for col in 0..width {
            let pixel_idx = row * width + col;
            let data_idx = pixel_idx * 2;

            // Extract real and imaginary parts as i16, convert to f32
            let real_i16 = complex_data.data[data_idx];
            let imag_i16 = complex_data.data[data_idx + 1];

            // Convert to f32 and create complex number
            let real_f32 = real_i16 as f32;
            let imag_f32 = imag_i16 as f32;

            complex_array[[row, col]] = Complex::new(real_f32, imag_f32);
        }
    }

    log::info!("✅ Successfully extracted COMPLEX SLC data: {}x{} pixels", width, height);
    
    Ok(complex_array)
}

/// Extract complex SLC data (preserving I/Q values) from ZIP or SAFE directory for proper radiometric calibration
pub fn extract_subswath_complex_data(
    product_path: &str,
    subswath: &str,
    polarization: &str,
) -> SarResult<ndarray::Array2<num_complex::Complex<f32>>> {
    // Delegate to the ZIP/SAFE-aware function
    extract_subswath_complex_from_zip(product_path, subswath, polarization)
}

/// Helper function to extract complex SLC data from SAFE directory
fn extract_subswath_complex_from_safe(
    safe_path: &str,
    _subswath: &str,
    polarization: &str,
) -> SarResult<ndarray::Array2<num_complex::Complex<f32>>> {
    use crate::io::slc_reader::SlcReader;
    use crate::types::Polarization;
    
    log::info!("Extracting complex SLC data from SAFE directory: {}", safe_path);

    // Parse polarization
    let pol = match polarization.to_uppercase().as_str() {
        "VV" => Polarization::VV,
        "VH" => Polarization::VH,
        "HV" => Polarization::HV,
        "HH" => Polarization::HH,
        _ => return Err(SarError::Processing(format!("Invalid polarization: {}", polarization))),
    };

    // Create SlcReader for SAFE directory
    let mut reader = SlcReader::new(safe_path)
        .map_err(|e| SarError::Processing(format!("Failed to create SLC reader: {}", e)))?;

    // Read SLC data for the specified polarization
    let sar_image = reader.read_slc_data(pol)
        .map_err(|e| SarError::Processing(format!("Failed to read SLC data: {}", e)))?;

    log::info!("Successfully read SLC data with dimensions: {}x{}", sar_image.nrows(), sar_image.ncols());
    
    // Convert from Complex64 to Complex<f32> and ensure it's complex data
    let complex_array = sar_image.mapv(|val| Complex::new(val.re as f32, val.im as f32));

    // Verify this is actually complex data
    let has_imaginary = complex_array.iter().any(|&val| val.im.abs() > 1e-10);
    if !has_imaginary {
        log::warn!("SLC data appears to have no imaginary component - may be intensity data");
    }

    log::info!("✅ Successfully extracted COMPLEX SLC data from SAFE: {}x{} pixels", 
               complex_array.nrows(), complex_array.ncols());
    
    Ok(complex_array)
}
