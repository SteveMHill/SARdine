/*!
 * Comprehensive SAR Quality Assessment and Masking Module
 * 
 * This module implements scientifically robust quality assessment for SAR data including:
 * - Signal-to-Noise Ratio (SNR) assessment
 * - Geometric distortion detection (layover/shadow/foreshortening)
 * - Per-pixel quality metrics
 * - Multi-criteria quality masking
 * - Statistical quality reporting
 * 
 * Based on:
 * - Freeman, A. (1992). SAR calibration: an overview. IEEE Trans. Geosci. Remote Sens.
 * - Small, D. (2011). Flattening gamma: Radiometric terrain correction for SAR imagery
 * - Kellndorfer, J. et al. (2004). Toward consistent regional-to-global scale vegetation characterization using orbital SAR systems
 */

use crate::types::*;
use ndarray::{Array2, s};
use serde::{Serialize, Deserialize};
// Removed unused imports: Array3, std::collections::HashMap

/// Comprehensive quality assessment configuration
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct QualityConfig {
    /// SNR-based quality assessment
    pub enable_snr_masking: bool,
    pub snr_threshold_db: f32,           // Minimum SNR in dB
    pub noise_equivalent_sigma0_db: f32, // NESZ in dB
    
    /// Geometric quality assessment
    pub enable_geometric_masking: bool,
    pub min_local_incidence_angle: f32,  // degrees
    pub max_local_incidence_angle: f32,  // degrees
    pub foreshortening_threshold: f32,   // ratio threshold
    
    /// Radiometric quality assessment
    pub enable_radiometric_masking: bool,
    pub min_gamma0_db: f32,
    pub max_gamma0_db: f32,
    pub dynamic_range_threshold: f32,    // dB
    
    /// Statistical quality assessment
    pub enable_statistical_masking: bool,
    pub outlier_sigma_threshold: f32,    // standard deviations
    pub texture_variance_threshold: f32,
    
    /// Per-pixel quality scoring
    pub enable_pixel_scoring: bool,
    pub quality_score_threshold: f32,    // 0-1 threshold
    
    /// Speckle coherence assessment
    pub enable_coherence_masking: bool,
    pub coherence_threshold: f32,        // 0-1 threshold
    pub coherence_window_size: usize,
}

impl Default for QualityConfig {
    fn default() -> Self {
        Self {
            // SNR masking - based on Sentinel-1 specifications
            enable_snr_masking: true,
            snr_threshold_db: 10.0,
            noise_equivalent_sigma0_db: -22.0,  // Typical Sentinel-1 NESZ
            
            // Geometric masking - based on SAR imaging geometry literature
            enable_geometric_masking: true,
            min_local_incidence_angle: 15.0,  // Avoid severe layover
            max_local_incidence_angle: 65.0,  // Avoid severe shadow
            foreshortening_threshold: 0.5,    // 50% foreshortening limit
            
            // Radiometric masking - typical SAR backscatter ranges
            enable_radiometric_masking: true,
            min_gamma0_db: -30.0,  // Below this likely noise/water
            max_gamma0_db: 5.0,    // Above this likely corner reflectors
            dynamic_range_threshold: 40.0,  // Minimum dynamic range
            
            // Statistical masking - outlier detection
            enable_statistical_masking: true,
            outlier_sigma_threshold: 3.0,     // 3-sigma outlier detection
            texture_variance_threshold: 10.0,  // High texture variance threshold
            
            // Per-pixel quality scoring
            enable_pixel_scoring: true,
            quality_score_threshold: 0.6,     // 60% quality threshold
            
            // Speckle coherence assessment
            enable_coherence_masking: false,   // Requires complex data
            coherence_threshold: 0.3,
            coherence_window_size: 7,
        }
    }
}

/// Comprehensive quality assessment result
#[derive(Debug, Clone)]
pub struct QualityAssessment {
    /// Individual quality masks
    pub snr_mask: Array2<bool>,
    pub geometric_mask: Array2<bool>,
    pub radiometric_mask: Array2<bool>,
    pub statistical_mask: Array2<bool>,
    pub coherence_mask: Option<Array2<bool>>,
    
    /// Combined quality mask (intersection of all enabled masks)
    pub combined_quality_mask: Array2<bool>,
    
    /// Per-pixel quality scores (0-1)
    pub pixel_quality_scores: Array2<f32>,
    
    /// Quality metrics
    pub snr_map: Array2<f32>,              // SNR in dB
    pub local_incidence_angles: Array2<f32>, // degrees
    pub foreshortening_factors: Array2<f32>, // ratio
    pub texture_measures: Array2<f32>,      // texture variance
    
    /// Quality statistics
    pub statistics: QualityStatistics,
}

/// Quality assessment statistics
#[derive(Debug, Clone)]
pub struct QualityStatistics {
    pub total_pixels: usize,
    pub valid_pixels: usize,
    pub valid_percentage: f32,
    
    // SNR statistics
    pub mean_snr_db: f32,
    pub min_snr_db: f32,
    pub max_snr_db: f32,
    pub snr_below_threshold: usize,
    
    // Geometric statistics
    pub mean_incidence_angle: f32,
    pub layover_pixels: usize,
    pub shadow_pixels: usize,
    pub foreshortened_pixels: usize,
    
    // Radiometric statistics
    pub mean_gamma0_db: f32,
    pub dynamic_range_db: f32,
    pub saturated_pixels: usize,
    pub low_backscatter_pixels: usize,
    
    // Statistical outliers
    pub outlier_pixels: usize,
    pub high_texture_pixels: usize,
    
    // Quality score distribution
    pub mean_quality_score: f32,
    pub quality_score_histogram: Vec<usize>, // 10 bins 0-1
}

/// Main quality assessment implementation
pub struct QualityAssessor;

impl QualityAssessor {
    /// Perform comprehensive quality assessment on SAR data
    pub fn assess_quality(
        gamma0_db: &Array2<f32>,
        dem_data: &Array2<f32>,
        local_incidence_angles: &Array2<f32>,
        orbit_params: &OrbitParams,
        config: &QualityConfig,
    ) -> SarResult<QualityAssessment> {
        let (rows, cols) = gamma0_db.dim();
        
        // Initialize quality assessment
        let mut assessment = QualityAssessment {
            snr_mask: Array2::from_elem((rows, cols), true),
            geometric_mask: Array2::from_elem((rows, cols), true),
            radiometric_mask: Array2::from_elem((rows, cols), true),
            statistical_mask: Array2::from_elem((rows, cols), true),
            coherence_mask: None,
            combined_quality_mask: Array2::from_elem((rows, cols), true),
            pixel_quality_scores: Array2::zeros((rows, cols)),
            snr_map: Array2::zeros((rows, cols)),
            local_incidence_angles: local_incidence_angles.clone(),
            foreshortening_factors: Array2::zeros((rows, cols)),
            texture_measures: Array2::zeros((rows, cols)),
            statistics: QualityStatistics::default(),
        };
        
        // 1. SNR Assessment
        if config.enable_snr_masking {
            Self::assess_snr(gamma0_db, config, &mut assessment)?;
        }
        
        // 2. Geometric Quality Assessment
        if config.enable_geometric_masking {
            Self::assess_geometric_quality(
                local_incidence_angles,
                dem_data,
                config,
                &mut assessment,
            )?;
        }
        
        // 3. Radiometric Quality Assessment
        if config.enable_radiometric_masking {
            Self::assess_radiometric_quality(gamma0_db, config, &mut assessment)?;
        }
        
        // 4. Statistical Quality Assessment
        if config.enable_statistical_masking {
            Self::assess_statistical_quality(gamma0_db, config, &mut assessment)?;
        }
        
        // 5. Texture Analysis
        Self::compute_texture_measures(gamma0_db, &mut assessment)?;
        
        // 6. Per-pixel Quality Scoring
        if config.enable_pixel_scoring {
            Self::compute_pixel_quality_scores(&mut assessment, config)?;
        }
        
        // 7. Combine Quality Masks
        Self::combine_quality_masks(&mut assessment, config)?;
        
        // 8. Compute Quality Statistics
        Self::compute_quality_statistics(&mut assessment, gamma0_db)?;
        
        log::info!("Quality assessment completed:");
        log::info!("  - Valid pixels: {}/{} ({:.1}%)", 
                   assessment.statistics.valid_pixels,
                   assessment.statistics.total_pixels,
                   assessment.statistics.valid_percentage);
        log::info!("  - Mean SNR: {:.1} dB", assessment.statistics.mean_snr_db);
        log::info!("  - Mean quality score: {:.3}", assessment.statistics.mean_quality_score);
        
        Ok(assessment)
    }
    
    /// Assess Signal-to-Noise Ratio (SNR)
    fn assess_snr(
        gamma0_db: &Array2<f32>,
        config: &QualityConfig,
        assessment: &mut QualityAssessment,
    ) -> SarResult<()> {
        let (rows, cols) = gamma0_db.dim();
        
        for i in 0..rows {
            for j in 0..cols {
                let gamma0_val = gamma0_db[[i, j]];
                
                // Calculate SNR relative to noise floor
                // SNR = Signal - NESZ
                let snr_db = gamma0_val - config.noise_equivalent_sigma0_db;
                assessment.snr_map[[i, j]] = snr_db;
                
                // Apply SNR threshold
                assessment.snr_mask[[i, j]] = snr_db >= config.snr_threshold_db 
                    && gamma0_val.is_finite();
            }
        }
        
        Ok(())
    }
    
    /// Assess geometric quality (layover, shadow, foreshortening)
    fn assess_geometric_quality(
        local_incidence_angles: &Array2<f32>,
        dem_data: &Array2<f32>,
        config: &QualityConfig,
        assessment: &mut QualityAssessment,
    ) -> SarResult<()> {
        let (rows, cols) = local_incidence_angles.dim();
        
        for i in 0..rows {
            for j in 0..cols {
                let lia_deg = local_incidence_angles[[i, j]];
                let elevation = dem_data[[i, j]];
                
                // Check local incidence angle bounds
                let lia_valid = lia_deg >= config.min_local_incidence_angle 
                    && lia_deg <= config.max_local_incidence_angle
                    && lia_deg.is_finite()
                    && elevation.is_finite();
                
                // Compute foreshortening factor
                // Foreshortening = cos(local_incidence_angle)
                let foreshortening = if lia_deg.is_finite() {
                    lia_deg.to_radians().cos()
                } else {
                    0.0
                };
                
                assessment.foreshortening_factors[[i, j]] = foreshortening;
                
                // Check foreshortening threshold
                let foreshortening_valid = foreshortening >= config.foreshortening_threshold;
                
                assessment.geometric_mask[[i, j]] = lia_valid && foreshortening_valid;
            }
        }
        
        Ok(())
    }
    
    /// Assess radiometric quality
    fn assess_radiometric_quality(
        gamma0_db: &Array2<f32>,
        config: &QualityConfig,
        assessment: &mut QualityAssessment,
    ) -> SarResult<()> {
        let (rows, cols) = gamma0_db.dim();
        
        // Compute global statistics for dynamic range assessment
        let mut valid_values = Vec::new();
        for &val in gamma0_db.iter() {
            if val.is_finite() {
                valid_values.push(val);
            }
        }
        
        let dynamic_range = if !valid_values.is_empty() {
            let min_val = valid_values.iter().fold(f32::INFINITY, |a, &b| a.min(b));
            let max_val = valid_values.iter().fold(f32::NEG_INFINITY, |a, &b| a.max(b));
            max_val - min_val
        } else {
            0.0
        };
        
        for i in 0..rows {
            for j in 0..cols {
                let gamma0_val = gamma0_db[[i, j]];
                
                // Check radiometric bounds
                let radiometric_valid = gamma0_val >= config.min_gamma0_db
                    && gamma0_val <= config.max_gamma0_db
                    && gamma0_val.is_finite()
                    && dynamic_range >= config.dynamic_range_threshold;
                
                assessment.radiometric_mask[[i, j]] = radiometric_valid;
            }
        }
        
        Ok(())
    }
    
    /// Assess statistical quality (outlier detection)
    fn assess_statistical_quality(
        gamma0_db: &Array2<f32>,
        config: &QualityConfig,
        assessment: &mut QualityAssessment,
    ) -> SarResult<()> {
        let (rows, cols) = gamma0_db.dim();
        
        // Compute global statistics
        let mut valid_values = Vec::new();
        for &val in gamma0_db.iter() {
            if val.is_finite() {
                valid_values.push(val);
            }
        }
        
        if valid_values.is_empty() {
            return Ok(());
        }
        
        let mean = valid_values.iter().sum::<f32>() / valid_values.len() as f32;
        let variance = valid_values.iter()
            .map(|&x| (x - mean).powi(2))
            .sum::<f32>() / valid_values.len() as f32;
        let std_dev = variance.sqrt();
        
        for i in 0..rows {
            for j in 0..cols {
                let gamma0_val = gamma0_db[[i, j]];
                
                // Outlier detection using sigma thresholding
                let z_score = if std_dev > 0.0 && gamma0_val.is_finite() {
                    (gamma0_val - mean).abs() / std_dev
                } else {
                    0.0
                };
                
                let statistical_valid = z_score <= config.outlier_sigma_threshold;
                
                assessment.statistical_mask[[i, j]] = statistical_valid;
            }
        }
        
        Ok(())
    }
    
    /// Compute texture measures using local variance
    fn compute_texture_measures(
        gamma0_db: &Array2<f32>,
        assessment: &mut QualityAssessment,
    ) -> SarResult<()> {
        let (rows, cols) = gamma0_db.dim();
        let window_size = 7; // 7x7 window for texture analysis
        let half_window = window_size / 2;
        
        for i in 0..rows {
            for j in 0..cols {
                let start_i = if i >= half_window { i - half_window } else { 0 };
                let end_i = (i + half_window + 1).min(rows);
                let start_j = if j >= half_window { j - half_window } else { 0 };
                let end_j = (j + half_window + 1).min(cols);
                
                // Extract window
                let window = gamma0_db.slice(s![start_i..end_i, start_j..end_j]);
                
                // Compute local variance as texture measure
                let mut valid_values = Vec::new();
                for &val in window.iter() {
                    if val.is_finite() {
                        valid_values.push(val);
                    }
                }
                
                let texture_variance = if valid_values.len() > 1 {
                    let mean = valid_values.iter().sum::<f32>() / valid_values.len() as f32;
                    valid_values.iter()
                        .map(|&x| (x - mean).powi(2))
                        .sum::<f32>() / (valid_values.len() - 1) as f32
                } else {
                    0.0
                };
                
                assessment.texture_measures[[i, j]] = texture_variance;
            }
        }
        
        Ok(())
    }
    
    /// Compute per-pixel quality scores (0-1)
    fn compute_pixel_quality_scores(
        assessment: &mut QualityAssessment,
        config: &QualityConfig,
    ) -> SarResult<()> {
        let (rows, cols) = assessment.snr_map.dim();
        
        for i in 0..rows {
            for j in 0..cols {
                let mut quality_score = 0.0;
                let mut weight_sum = 0.0;
                
                // SNR component (weight: 0.3)
                if config.enable_snr_masking {
                    let snr_normalized = ((assessment.snr_map[[i, j]] - config.snr_threshold_db) / 20.0)
                        .max(0.0).min(1.0);
                    quality_score += 0.3 * snr_normalized;
                    weight_sum += 0.3;
                }
                
                // Geometric component (weight: 0.25)
                if config.enable_geometric_masking {
                    let lia = assessment.local_incidence_angles[[i, j]];
                    let lia_score = if lia >= config.min_local_incidence_angle && lia <= config.max_local_incidence_angle {
                        // Optimal around 30-40 degrees
                        let optimal_lia = 35.0;
                        1.0 - ((lia - optimal_lia).abs() / optimal_lia).min(1.0)
                    } else {
                        0.0
                    };
                    
                    let foreshortening_score = assessment.foreshortening_factors[[i, j]]
                        .max(0.0).min(1.0);
                    
                    quality_score += 0.25 * (lia_score + foreshortening_score) / 2.0;
                    weight_sum += 0.25;
                }
                
                // Radiometric component (weight: 0.2)
                if config.enable_radiometric_masking && assessment.radiometric_mask[[i, j]] {
                    quality_score += 0.2;
                    weight_sum += 0.2;
                }
                
                // Statistical component (weight: 0.15)
                if config.enable_statistical_masking && assessment.statistical_mask[[i, j]] {
                    quality_score += 0.15;
                    weight_sum += 0.15;
                }
                
                // Texture component (weight: 0.1)
                let texture_score = (1.0 - (assessment.texture_measures[[i, j]] / config.texture_variance_threshold).min(1.0))
                    .max(0.0);
                quality_score += 0.1 * texture_score;
                weight_sum += 0.1;
                
                // Normalize by weight sum
                assessment.pixel_quality_scores[[i, j]] = if weight_sum > 0.0 {
                    quality_score / weight_sum
                } else {
                    0.0
                };
            }
        }
        
        Ok(())
    }
    
    /// Combine all quality masks
    fn combine_quality_masks(
        assessment: &mut QualityAssessment,
        config: &QualityConfig,
    ) -> SarResult<()> {
        let (rows, cols) = assessment.combined_quality_mask.dim();
        
        for i in 0..rows {
            for j in 0..cols {
                let mut is_valid = true;
                
                if config.enable_snr_masking {
                    is_valid &= assessment.snr_mask[[i, j]];
                }
                
                if config.enable_geometric_masking {
                    is_valid &= assessment.geometric_mask[[i, j]];
                }
                
                if config.enable_radiometric_masking {
                    is_valid &= assessment.radiometric_mask[[i, j]];
                }
                
                if config.enable_statistical_masking {
                    is_valid &= assessment.statistical_mask[[i, j]];
                }
                
                if config.enable_pixel_scoring {
                    is_valid &= assessment.pixel_quality_scores[[i, j]] >= config.quality_score_threshold;
                }
                
                if let Some(ref coherence_mask) = assessment.coherence_mask {
                    if config.enable_coherence_masking {
                        is_valid &= coherence_mask[[i, j]];
                    }
                }
                
                assessment.combined_quality_mask[[i, j]] = is_valid;
            }
        }
        
        Ok(())
    }
    
    /// Compute comprehensive quality statistics
    fn compute_quality_statistics(
        assessment: &mut QualityAssessment,
        gamma0_db: &Array2<f32>,
    ) -> SarResult<()> {
        let (rows, cols) = assessment.combined_quality_mask.dim();
        
        let mut stats = QualityStatistics::default();
        stats.total_pixels = rows * cols;
        
        let mut snr_values = Vec::new();
        let mut lia_values = Vec::new();
        let mut gamma0_values = Vec::new();
        let mut quality_scores = Vec::new();
        
        for i in 0..rows {
            for j in 0..cols {
                if assessment.combined_quality_mask[[i, j]] {
                    stats.valid_pixels += 1;
                }
                
                let snr = assessment.snr_map[[i, j]];
                let lia = assessment.local_incidence_angles[[i, j]];
                let gamma0 = gamma0_db[[i, j]];
                let quality_score = assessment.pixel_quality_scores[[i, j]];
                
                if snr.is_finite() {
                    snr_values.push(snr);
                    if snr < 10.0 { // Default SNR threshold
                        stats.snr_below_threshold += 1;
                    }
                }
                
                if lia.is_finite() {
                    lia_values.push(lia);
                    if lia < 15.0 {
                        stats.layover_pixels += 1;
                    } else if lia > 65.0 {
                        stats.shadow_pixels += 1;
                    }
                }
                
                if gamma0.is_finite() {
                    gamma0_values.push(gamma0);
                    if gamma0 < -30.0 {
                        stats.low_backscatter_pixels += 1;
                    } else if gamma0 > 5.0 {
                        stats.saturated_pixels += 1;
                    }
                }
                
                quality_scores.push(quality_score);
                
                if assessment.foreshortening_factors[[i, j]] < 0.5 {
                    stats.foreshortened_pixels += 1;
                }
                
                if !assessment.statistical_mask[[i, j]] {
                    stats.outlier_pixels += 1;
                }
            }
        }
        
        // Compute means
        stats.valid_percentage = (stats.valid_pixels as f32 / stats.total_pixels as f32) * 100.0;
        
        if !snr_values.is_empty() {
            stats.mean_snr_db = snr_values.iter().sum::<f32>() / snr_values.len() as f32;
            stats.min_snr_db = snr_values.iter().fold(f32::INFINITY, |a, &b| a.min(b));
            stats.max_snr_db = snr_values.iter().fold(f32::NEG_INFINITY, |a, &b| a.max(b));
        }
        
        if !lia_values.is_empty() {
            stats.mean_incidence_angle = lia_values.iter().sum::<f32>() / lia_values.len() as f32;
        }
        
        if !gamma0_values.is_empty() {
            stats.mean_gamma0_db = gamma0_values.iter().sum::<f32>() / gamma0_values.len() as f32;
            let min_gamma0 = gamma0_values.iter().fold(f32::INFINITY, |a, &b| a.min(b));
            let max_gamma0 = gamma0_values.iter().fold(f32::NEG_INFINITY, |a, &b| a.max(b));
            stats.dynamic_range_db = max_gamma0 - min_gamma0;
        }
        
        if !quality_scores.is_empty() {
            stats.mean_quality_score = quality_scores.iter().sum::<f32>() / quality_scores.len() as f32;
            
            // Compute quality score histogram (10 bins)
            stats.quality_score_histogram = vec![0; 10];
            for &score in &quality_scores {
                let bin = ((score * 10.0) as usize).min(9);
                stats.quality_score_histogram[bin] += 1;
            }
        }
        
        assessment.statistics = stats;
        Ok(())
    }
}

impl Default for QualityStatistics {
    fn default() -> Self {
        Self {
            total_pixels: 0,
            valid_pixels: 0,
            valid_percentage: 0.0,
            mean_snr_db: 0.0,
            min_snr_db: 0.0,
            max_snr_db: 0.0,
            snr_below_threshold: 0,
            mean_incidence_angle: 0.0,
            layover_pixels: 0,
            shadow_pixels: 0,
            foreshortened_pixels: 0,
            mean_gamma0_db: 0.0,
            dynamic_range_db: 0.0,
            saturated_pixels: 0,
            low_backscatter_pixels: 0,
            outlier_pixels: 0,
            high_texture_pixels: 0,
            mean_quality_score: 0.0,
            quality_score_histogram: vec![0; 10],
        }
    }
}
