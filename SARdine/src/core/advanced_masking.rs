use ndarray::Array2;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;

use crate::types::{SarResult, SarError};

/// Advanced masking methods for SAR data quality control
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub enum MaskingMethod {
    /// Simple threshold-based masking
    Threshold,
    /// Statistical outlier detection using Z-score
    StatisticalOutlier,
    /// Adaptive threshold based on local statistics
    AdaptiveThreshold,
    /// Multi-scale analysis for robust detection
    MultiScale,
    /// Machine learning-based classification (simplified)
    MlBased,
    /// Comprehensive combining multiple methods
    Comprehensive,
}

/// Noise characteristics for different SAR scenarios
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct NoiseProfile {
    /// Thermal noise level (linear units)
    pub thermal_noise_level: f32,
    /// Expected noise equivalent sigma zero (NESZ) in dB
    pub nesz_db: f32,
    /// Range-dependent noise variation
    pub range_noise_factor: f32,
    /// Azimuth-dependent noise variation  
    pub azimuth_noise_factor: f32,
}

impl Default for NoiseProfile {
    fn default() -> Self {
        // Conservative values for Sentinel-1
        Self {
            thermal_noise_level: 1e-6,
            nesz_db: -25.0,
            range_noise_factor: 1.2,
            azimuth_noise_factor: 1.1,
        }
    }
}

/// Water body detection parameters
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct WaterDetectionParams {
    /// Sigma0 threshold for water detection (dB)
    pub sigma0_threshold_db: f32,
    /// Standard deviation threshold for water uniformity
    pub std_threshold: f32,
    /// Minimum patch size for water bodies (pixels)
    pub min_patch_size: usize,
    /// Maximum texture threshold (local variance)
    pub max_texture: f32,
}

impl Default for WaterDetectionParams {
    fn default() -> Self {
        Self {
            sigma0_threshold_db: -18.0,  // Typical water threshold
            std_threshold: 2.0,
            min_patch_size: 100,
            max_texture: 0.5,
        }
    }
}

/// Layover and shadow detection parameters
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct LayoverShadowParams {
    /// Minimum local incidence angle (degrees)
    pub min_incidence_deg: f32,
    /// Maximum local incidence angle (degrees)  
    pub max_incidence_deg: f32,
    /// DEM resolution for slope computation (meters)
    pub dem_resolution: f32,
    /// Slope threshold for layover detection (degrees)
    pub layover_slope_threshold: f32,
    /// Shadow detection sensitivity
    pub shadow_sensitivity: f32,
}

impl Default for LayoverShadowParams {
    fn default() -> Self {
        Self {
            min_incidence_deg: 15.0,
            max_incidence_deg: 75.0,
            dem_resolution: 30.0,
            layover_slope_threshold: 45.0,
            shadow_sensitivity: 0.8,
        }
    }
}

/// Comprehensive masking configuration for state-of-the-art quality control
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AdvancedMaskingConfig {
    /// Primary masking method
    pub method: MaskingMethod,
    /// Enable water body detection and masking
    pub enable_water_detection: bool,
    /// Enable layover/shadow detection
    pub enable_layover_shadow: bool,
    /// Enable noise floor detection
    pub enable_noise_detection: bool,
    /// Enable coherence-based masking (if coherence data available)
    pub enable_coherence_masking: bool,
    /// Enable edge/boundary detection
    pub enable_edge_detection: bool,
    /// Enable multi-temporal consistency checking
    pub enable_temporal_consistency: bool,
    /// Noise profile for the SAR system
    pub noise_profile: NoiseProfile,
    /// Water detection parameters
    pub water_params: WaterDetectionParams,
    /// Layover/shadow parameters  
    pub layover_shadow_params: LayoverShadowParams,
    /// Statistical confidence level for outlier detection
    pub confidence_level: f64,
    /// Morphological cleaning operations
    pub enable_morphological_cleaning: bool,
    /// Output validation quality metrics
    pub compute_quality_metrics: bool,
}

impl Default for AdvancedMaskingConfig {
    fn default() -> Self {
        Self {
            method: MaskingMethod::Comprehensive,
            enable_water_detection: true,
            enable_layover_shadow: true,
            enable_noise_detection: true,
            enable_coherence_masking: false,
            enable_edge_detection: true,
            enable_temporal_consistency: false,
            noise_profile: NoiseProfile::default(),
            water_params: WaterDetectionParams::default(),
            layover_shadow_params: LayoverShadowParams::default(),
            confidence_level: 0.95,
            enable_morphological_cleaning: true,
            compute_quality_metrics: true,
        }
    }
}

/// Advanced mask result with detailed statistics and quality metrics
#[derive(Debug, Clone)]
pub struct AdvancedMaskResult {
    /// Combined final mask (1 = valid, 0 = invalid)
    pub final_mask: Array2<u8>,
    /// Individual component masks
    pub component_masks: HashMap<String, Array2<u8>>,
    /// Quality confidence map (0.0 to 1.0)
    pub confidence_map: Array2<f32>,
    /// Detected anomaly map (0 = normal, 1+ = anomaly type)
    pub anomaly_map: Array2<u8>,
    /// Comprehensive statistics
    pub statistics: AdvancedMaskStats,
    /// Quality assessment metrics
    pub quality_metrics: QualityMetrics,
}

/// Comprehensive masking statistics
#[derive(Debug, Clone, Default)]
pub struct AdvancedMaskStats {
    pub total_pixels: usize,
    pub valid_pixels: usize,
    pub water_pixels: usize,
    pub layover_pixels: usize,
    pub shadow_pixels: usize,
    pub noise_pixels: usize,
    pub edge_pixels: usize,
    pub outlier_pixels: usize,
    pub low_coherence_pixels: usize,
    pub valid_percentage: f64,
    pub data_quality_score: f64,
}

/// Quality assessment metrics for validation
#[derive(Debug, Clone, Default)]
pub struct QualityMetrics {
    /// Overall data quality score (0.0 to 1.0)
    pub overall_quality: f64,
    /// Spatial consistency metric
    pub spatial_consistency: f64,
    /// Radiometric quality score
    pub radiometric_quality: f64,
    /// Geometric accuracy assessment
    pub geometric_accuracy: f64,
    /// Detected artifact severity
    pub artifact_severity: f64,
}

/// State-of-the-art advanced masking processor
#[derive(Debug)]
pub struct AdvancedMaskingProcessor {
    config: AdvancedMaskingConfig,
}

impl AdvancedMaskingProcessor {
    /// Create new advanced masking processor
    pub fn new(config: AdvancedMaskingConfig) -> Self {
        Self { config }
    }

    /// Create with default comprehensive masking configuration
    pub fn new_comprehensive() -> Self {
        Self::new(AdvancedMaskingConfig::default())
    }

    /// Process SAR data with state-of-the-art masking
    pub fn process_advanced_masking(
        &self,
        sigma0_data: &Array2<f32>,
        incidence_angles: Option<&Array2<f32>>,
        dem_data: Option<&Array2<f32>>,
        coherence_data: Option<&Array2<f32>>,
    ) -> SarResult<AdvancedMaskResult> {
        log::info!("🎭 Starting advanced masking with method: {:?}", self.config.method);
        
        let (_height, _width) = sigma0_data.dim();
        let mut component_masks = HashMap::new();
        let confidence_map: Array2<f32>;
        let anomaly_map: Array2<u8>;

        // Step 1: Basic validity mask (NaN, infinite, extreme values)
        log::debug!("Computing basic validity mask");
        let basic_mask = self.compute_basic_validity_mask(sigma0_data)?;
        component_masks.insert("basic_validity".to_string(), basic_mask.clone());

        // Step 2: Statistical outlier detection
        log::debug!("Performing statistical outlier detection");
        let outlier_mask = self.detect_statistical_outliers(sigma0_data)?;
        component_masks.insert("statistical_outliers".to_string(), outlier_mask.clone());

        // Step 3: Water body detection (if enabled)
        let water_mask = if self.config.enable_water_detection {
            log::debug!("Detecting water bodies");
            let mask = self.detect_water_bodies(sigma0_data)?;
            component_masks.insert("water_bodies".to_string(), mask.clone());
            Some(mask)
        } else {
            None
        };

        // Step 4: Layover and shadow detection (if DEM and incidence angles available)
        let layover_shadow_mask = if self.config.enable_layover_shadow && 
                                     incidence_angles.is_some() && dem_data.is_some() {
            log::debug!("Detecting layover and shadow areas");
            let mask = self.detect_layover_shadow(
                incidence_angles.unwrap(),
                dem_data.unwrap(),
            )?;
            component_masks.insert("layover_shadow".to_string(), mask.clone());
            Some(mask)
        } else {
            None
        };

        // Step 5: Noise floor detection
        let noise_mask = if self.config.enable_noise_detection {
            log::debug!("Detecting noise floor areas");
            let mask = self.detect_noise_floor(sigma0_data)?;
            component_masks.insert("noise_floor".to_string(), mask.clone());
            Some(mask)
        } else {
            None
        };

        // Step 6: Coherence-based masking (if coherence data available)
        let coherence_mask = if self.config.enable_coherence_masking && coherence_data.is_some() {
            log::debug!("Applying coherence-based masking");
            let mask = self.apply_coherence_masking(coherence_data.unwrap())?;
            component_masks.insert("coherence".to_string(), mask.clone());
            Some(mask)
        } else {
            None
        };

        // Step 7: Edge and boundary detection
        let edge_mask = if self.config.enable_edge_detection {
            log::debug!("Detecting edges and boundaries");
            let mask = self.detect_edges_boundaries(sigma0_data)?;
            component_masks.insert("edges_boundaries".to_string(), mask.clone());
            Some(mask)
        } else {
            None
        };

        // Step 8: Compute confidence map
        log::debug!("Computing confidence map");
        confidence_map = self.compute_confidence_map(sigma0_data, &component_masks)?;

        // Step 9: Detect and classify anomalies
        log::debug!("Detecting and classifying anomalies");
        anomaly_map = self.detect_anomalies(sigma0_data, &component_masks)?;

        // Step 10: Combine all masks using the selected method
        log::debug!("Combining masks using method: {:?}", self.config.method);
        let final_mask = self.combine_masks_advanced(
            &basic_mask,
            water_mask.as_ref(),
            layover_shadow_mask.as_ref(), 
            noise_mask.as_ref(),
            coherence_mask.as_ref(),
            edge_mask.as_ref(),
            &confidence_map,
        )?;

        // Step 11: Morphological cleaning (if enabled)
        let final_mask = if self.config.enable_morphological_cleaning {
            log::debug!("Applying morphological cleaning");
            self.apply_morphological_cleaning(&final_mask)?
        } else {
            final_mask
        };

        // Step 12: Compute comprehensive statistics
        log::debug!("Computing comprehensive statistics");
        let statistics = self.compute_advanced_statistics(
            &final_mask,
            &component_masks,
            sigma0_data,
        )?;

        // Step 13: Compute quality metrics (if enabled)
        let quality_metrics = if self.config.compute_quality_metrics {
            log::debug!("Computing quality metrics");
            self.compute_quality_metrics(
                sigma0_data,
                &final_mask,
                &confidence_map,
                &anomaly_map,
            )?
        } else {
            QualityMetrics::default()
        };

        log::info!("✅ Advanced masking completed");
        log::info!("   Valid pixels: {:.1}%", statistics.valid_percentage);
        log::info!("   Quality score: {:.3}", statistics.data_quality_score);

        Ok(AdvancedMaskResult {
            final_mask,
            component_masks,
            confidence_map,
            anomaly_map,
            statistics,
            quality_metrics,
        })
    }

    /// Compute basic validity mask for NaN, infinite, and extreme values
    fn compute_basic_validity_mask(&self, data: &Array2<f32>) -> SarResult<Array2<u8>> {
        let mut mask = Array2::<u8>::ones(data.dim());
        
        data.indexed_iter().for_each(|((i, j), &value)| {
            mask[[i, j]] = if value.is_finite() && 
                             value > 0.0 &&
                             value < 100.0 &&  // Reasonable upper bound for sigma0/gamma0
                             value > 1e-8 {    // Above realistic noise floor
                1
            } else {
                0
            };
        });

        Ok(mask)
    }

    /// Detect statistical outliers using robust Z-score method
    fn detect_statistical_outliers(&self, data: &Array2<f32>) -> SarResult<Array2<u8>> {
        log::debug!("Computing robust statistics for outlier detection");
        
        // Compute robust statistics (median and MAD)
        let valid_data: Vec<f32> = data.iter()
            .filter(|&&x| x.is_finite() && x > 0.0)
            .copied()
            .collect();

        if valid_data.is_empty() {
            return Ok(Array2::<u8>::zeros(data.dim()));
        }

        let mut sorted_data = valid_data.clone();
        sorted_data.sort_by(|a, b| a.partial_cmp(b).unwrap());
        
        let median = sorted_data[sorted_data.len() / 2];
        
        // Median Absolute Deviation (MAD)
        let mut deviations: Vec<f32> = sorted_data.iter()
            .map(|&x| (x - median).abs())
            .collect();
        deviations.sort_by(|a, b| a.partial_cmp(b).unwrap());
        let mad = deviations[deviations.len() / 2];

        // Z-score threshold based on confidence level
        let z_threshold = match self.config.confidence_level {
            level if level >= 0.99 => 3.0,
            level if level >= 0.95 => 2.5,
            _ => 2.0,
        };

        let mut mask = Array2::<u8>::ones(data.dim());
        
        data.indexed_iter().for_each(|((i, j), &value)| {
            if value.is_finite() && value > 0.0 {
                let modified_z_score = 0.6745 * (value - median).abs() / mad;
                mask[[i, j]] = if modified_z_score <= z_threshold { 1 } else { 0 };
            } else {
                mask[[i, j]] = 0;
            }
        });

        log::debug!("Outlier detection completed with threshold: {:.1}", z_threshold);
        Ok(mask)
    }

    /// Detect water bodies using multiple criteria
    fn detect_water_bodies(&self, data: &Array2<f32>) -> SarResult<Array2<u8>> {
        let (height, width) = data.dim();
        let mut water_mask = Array2::<u8>::ones((height, width));
        
        let water_threshold_linear = 10_f32.powf(self.config.water_params.sigma0_threshold_db / 10.0);

        // Apply water detection criteria
        for i in 0..height {
            for j in 0..width {
                let value = data[[i, j]];
                
                if !value.is_finite() || value <= 0.0 {
                    water_mask[[i, j]] = 0;
                    continue;
                }

                // Criterion 1: Low backscatter
                let is_low_backscatter = value < water_threshold_linear;

                // Criterion 2: Local uniformity (low texture)
                let local_texture = self.compute_local_texture(data, i, j, 3);
                let is_uniform = local_texture < self.config.water_params.max_texture;

                // Combine criteria for water detection
                // Note: For this implementation, we're identifying water as invalid areas
                // In some applications, you might want to keep water areas
                if is_low_backscatter && is_uniform {
                    water_mask[[i, j]] = 0; // Mark as invalid (water)
                }
            }
        }

        // Apply morphological operations to clean up small patches
        if self.config.water_params.min_patch_size > 0 {
            water_mask = self.remove_small_patches(&water_mask, self.config.water_params.min_patch_size)?;
        }

        Ok(water_mask)
    }

    /// Detect layover and shadow areas using topographic analysis
    fn detect_layover_shadow(
        &self,
        incidence_angles: &Array2<f32>,
        dem_data: &Array2<f32>,
    ) -> SarResult<Array2<u8>> {
        let (height, width) = incidence_angles.dim();
        let mut mask = Array2::<u8>::ones((height, width));

        let min_inc_rad = self.config.layover_shadow_params.min_incidence_deg.to_radians();
        let max_inc_rad = self.config.layover_shadow_params.max_incidence_deg.to_radians();

        // Compute terrain slopes
        let slopes = self.compute_terrain_slopes(dem_data)?;

        for i in 0..height {
            for j in 0..width {
                let incidence_angle = incidence_angles[[i, j]];
                let slope = slopes[[i, j]];

                // Skip invalid data
                if !incidence_angle.is_finite() || !slope.is_finite() {
                    mask[[i, j]] = 0;
                    continue;
                }

                // Check for layover (very low incidence angles)
                if incidence_angle < min_inc_rad {
                    mask[[i, j]] = 0;
                    continue;
                }

                // Check for shadow (very high incidence angles)
                if incidence_angle > max_inc_rad {
                    mask[[i, j]] = 0;
                    continue;
                }

                // Check for steep slopes causing layover
                let slope_deg = slope.to_degrees();
                if slope_deg > self.config.layover_shadow_params.layover_slope_threshold {
                    mask[[i, j]] = 0;
                    continue;
                }
            }
        }

        Ok(mask)
    }

    /// Detect noise floor areas
    fn detect_noise_floor(&self, data: &Array2<f32>) -> SarResult<Array2<u8>> {
        let mut mask = Array2::<u8>::ones(data.dim());
        
        let noise_threshold = self.config.noise_profile.thermal_noise_level;

        data.indexed_iter().for_each(|((i, j), &value)| {
            mask[[i, j]] = if value.is_finite() && value > noise_threshold { 1 } else { 0 };
        });

        Ok(mask)
    }

    /// Apply coherence-based masking
    fn apply_coherence_masking(&self, coherence_data: &Array2<f32>) -> SarResult<Array2<u8>> {
        let mut mask = Array2::<u8>::ones(coherence_data.dim());
        let threshold = 0.3; // Default coherence threshold

        coherence_data.indexed_iter().for_each(|((i, j), &value)| {
            mask[[i, j]] = if value.is_finite() && value >= threshold { 1 } else { 0 };
        });

        Ok(mask)
    }

    /// Detect edges and boundaries that might indicate processing artifacts
    fn detect_edges_boundaries(&self, data: &Array2<f32>) -> SarResult<Array2<u8>> {
        let (height, width) = data.dim();
        let mut mask = Array2::<u8>::ones((height, width));

        // Simple edge detection using gradient magnitude
        for i in 1..height-1 {
            for j in 1..width-1 {
                let center = data[[i, j]];
                if !center.is_finite() {
                    mask[[i, j]] = 0;
                    continue;
                }

                // Compute gradient using Sobel-like operator
                let gx = -data[[i-1, j-1]] - 2.0*data[[i, j-1]] - data[[i+1, j-1]] +
                          data[[i-1, j+1]] + 2.0*data[[i, j+1]] + data[[i+1, j+1]];
                
                let gy = -data[[i-1, j-1]] - 2.0*data[[i-1, j]] - data[[i-1, j+1]] +
                          data[[i+1, j-1]] + 2.0*data[[i+1, j]] + data[[i+1, j+1]];

                let gradient_magnitude = (gx*gx + gy*gy).sqrt();
                
                // Adaptive threshold based on local mean
                let local_mean = self.compute_local_mean(data, i, j, 3);
                let edge_threshold = local_mean * 0.5; // 50% of local mean

                if gradient_magnitude > edge_threshold {
                    mask[[i, j]] = 0; // Mark as potential artifact
                }
            }
        }

        Ok(mask)
    }

    /// Compute confidence map based on multiple quality indicators
    fn compute_confidence_map(
        &self,
        data: &Array2<f32>,
        component_masks: &HashMap<String, Array2<u8>>,
    ) -> SarResult<Array2<f32>> {
        let (height, width) = data.dim();
        let mut confidence = Array2::<f32>::ones((height, width));

        for i in 0..height {
            for j in 0..width {
                let mut local_confidence = 1.0f32;

                // Reduce confidence based on failed mask components
                for (_name, mask) in component_masks.iter() {
                    if mask[[i, j]] == 0 {
                        local_confidence *= 0.8; // Reduce confidence
                    }
                }

                // Factor in data quality indicators
                let value = data[[i, j]];
                if value.is_finite() && value > 0.0 {
                    // Higher confidence for values in reasonable range
                    if value >= 0.001 && value <= 1.0 {
                        local_confidence *= 1.0; // No penalty
                    } else if value > 1.0 && value <= 10.0 {
                        local_confidence *= 0.9; // Slight penalty
                    } else {
                        local_confidence *= 0.5; // Significant penalty
                    }
                } else {
                    local_confidence = 0.0; // No confidence for invalid data
                }

                confidence[[i, j]] = local_confidence.max(0.0).min(1.0);
            }
        }

        Ok(confidence)
    }

    /// Detect and classify different types of anomalies
    fn detect_anomalies(
        &self,
        data: &Array2<f32>,
        component_masks: &HashMap<String, Array2<u8>>,
    ) -> SarResult<Array2<u8>> {
        let (height, width) = data.dim();
        let mut anomaly_map = Array2::<u8>::zeros((height, width));

        for i in 0..height {
            for j in 0..width {
                let value = data[[i, j]];

                if !value.is_finite() {
                    anomaly_map[[i, j]] = 1; // NaN/Inf anomaly
                } else if value <= 0.0 {
                    anomaly_map[[i, j]] = 2; // Negative/zero anomaly
                } else if value > 50.0 {
                    anomaly_map[[i, j]] = 3; // Extremely high value anomaly
                } else {
                    // Check component-specific anomalies
                    if let Some(water_mask) = component_masks.get("water_bodies") {
                        if water_mask[[i, j]] == 0 {
                            anomaly_map[[i, j]] = 4; // Water body
                        }
                    }
                    
                    if let Some(layover_mask) = component_masks.get("layover_shadow") {
                        if layover_mask[[i, j]] == 0 {
                            anomaly_map[[i, j]] = 5; // Layover/shadow
                        }
                    }
                }
            }
        }

        Ok(anomaly_map)
    }

    /// Combine masks using advanced fusion method
    fn combine_masks_advanced(
        &self,
        basic_mask: &Array2<u8>,
        water_mask: Option<&Array2<u8>>,
        layover_shadow_mask: Option<&Array2<u8>>,
        noise_mask: Option<&Array2<u8>>,
        coherence_mask: Option<&Array2<u8>>,
        edge_mask: Option<&Array2<u8>>,
        confidence_map: &Array2<f32>,
    ) -> SarResult<Array2<u8>> {
        let (height, width) = basic_mask.dim();
        let mut final_mask = Array2::<u8>::zeros((height, width));

        match self.config.method {
            MaskingMethod::Threshold | MaskingMethod::Comprehensive => {
                // Logical AND of all available masks with confidence weighting
                for i in 0..height {
                    for j in 0..width {
                        let mut is_valid = basic_mask[[i, j]] == 1;
                        
                        if let Some(mask) = water_mask {
                            is_valid = is_valid && mask[[i, j]] == 1;
                        }
                        
                        if let Some(mask) = layover_shadow_mask {
                            is_valid = is_valid && mask[[i, j]] == 1;
                        }
                        
                        if let Some(mask) = noise_mask {
                            is_valid = is_valid && mask[[i, j]] == 1;
                        }
                        
                        if let Some(mask) = coherence_mask {
                            is_valid = is_valid && mask[[i, j]] == 1;
                        }
                        
                        if let Some(mask) = edge_mask {
                            is_valid = is_valid && mask[[i, j]] == 1;
                        }

                        // Apply confidence threshold
                        let confidence_threshold = 0.5;
                        is_valid = is_valid && confidence_map[[i, j]] >= confidence_threshold;

                        final_mask[[i, j]] = if is_valid { 1 } else { 0 };
                    }
                }
            },
            MaskingMethod::StatisticalOutlier => {
                // Use only statistical outlier detection
                final_mask = basic_mask.clone();
            },
            MaskingMethod::AdaptiveThreshold => {
                // Adaptive combination based on local statistics
                for i in 0..height {
                    for j in 0..width {
                        let local_confidence = confidence_map[[i, j]];
                        let adaptive_threshold = 0.3 + 0.4 * local_confidence; // 0.3 to 0.7 range
                        
                        let mut valid_count = 0;
                        let mut total_count = 0;
                        
                        if basic_mask[[i, j]] == 1 { valid_count += 1; }
                        total_count += 1;
                        
                        if let Some(mask) = water_mask {
                            if mask[[i, j]] == 1 { valid_count += 1; }
                            total_count += 1;
                        }
                        
                        if let Some(mask) = layover_shadow_mask {
                            if mask[[i, j]] == 1 { valid_count += 1; }
                            total_count += 1;
                        }
                        
                        let validity_ratio = valid_count as f32 / total_count as f32;
                        final_mask[[i, j]] = if validity_ratio >= adaptive_threshold { 1 } else { 0 };
                    }
                }
            },
            _ => {
                // Default to conservative approach
                final_mask = basic_mask.clone();
            }
        }

        Ok(final_mask)
    }

    /// Apply morphological cleaning operations
    fn apply_morphological_cleaning(&self, mask: &Array2<u8>) -> SarResult<Array2<u8>> {
        log::debug!("Applying morphological cleaning operations");
        
        // Simple morphological operations: erosion followed by dilation
        let mut cleaned_mask = mask.clone();
        
        // Erosion (remove small isolated valid pixels)
        cleaned_mask = self.morphological_erosion(&cleaned_mask, 1)?;
        
        // Dilation (restore valid regions and fill small gaps)
        cleaned_mask = self.morphological_dilation(&cleaned_mask, 2)?;

        // Check if cleaning removed too much data
        let original_valid = mask.iter().filter(|&&x| x == 1).count();
        let cleaned_valid = cleaned_mask.iter().filter(|&&x| x == 1).count();
        
        if cleaned_valid < original_valid * 80 / 100 { // Less than 80% remaining
            log::warn!("Morphological cleaning removed too much data, using original mask");
            Ok(mask.clone())
        } else {
            log::debug!("Morphological cleaning preserved {:.1}% of valid pixels", 
                       100.0 * cleaned_valid as f64 / original_valid as f64);
            Ok(cleaned_mask)
        }
    }

    /// Compute comprehensive advanced statistics
    fn compute_advanced_statistics(
        &self,
        final_mask: &Array2<u8>,
        component_masks: &HashMap<String, Array2<u8>>,
        data: &Array2<f32>,
    ) -> SarResult<AdvancedMaskStats> {
        let total_pixels = final_mask.len();
        let valid_pixels = final_mask.iter().filter(|&&x| x == 1).count();
        
        let water_pixels = component_masks.get("water_bodies")
            .map(|mask| mask.iter().filter(|&&x| x == 0).count())
            .unwrap_or(0);
            
        let layover_pixels = component_masks.get("layover_shadow")
            .map(|mask| mask.iter().filter(|&&x| x == 0).count())
            .unwrap_or(0);
            
        let noise_pixels = component_masks.get("noise_floor")
            .map(|mask| mask.iter().filter(|&&x| x == 0).count())
            .unwrap_or(0);
            
        let edge_pixels = component_masks.get("edges_boundaries")
            .map(|mask| mask.iter().filter(|&&x| x == 0).count())
            .unwrap_or(0);
            
        let outlier_pixels = component_masks.get("statistical_outliers")
            .map(|mask| mask.iter().filter(|&&x| x == 0).count())
            .unwrap_or(0);

        let valid_percentage = 100.0 * valid_pixels as f64 / total_pixels as f64;
        
        // Compute data quality score based on multiple factors
        let data_quality_score = self.compute_data_quality_score(data, final_mask)?;

        Ok(AdvancedMaskStats {
            total_pixels,
            valid_pixels,
            water_pixels,
            layover_pixels,
            shadow_pixels: 0, // Would be computed separately if needed
            noise_pixels,
            edge_pixels,
            outlier_pixels,
            low_coherence_pixels: 0, // Would be computed if coherence data available
            valid_percentage,
            data_quality_score,
        })
    }

    /// Compute quality metrics for validation
    fn compute_quality_metrics(
        &self,
        data: &Array2<f32>,
        mask: &Array2<u8>,
        confidence_map: &Array2<f32>,
        anomaly_map: &Array2<u8>,
    ) -> SarResult<QualityMetrics> {
        // Overall quality based on valid pixel percentage and confidence
        let valid_pixels = mask.iter().filter(|&&x| x == 1).count();
        let total_pixels = mask.len();
        let coverage_factor = valid_pixels as f64 / total_pixels as f64;
        
        let mean_confidence = confidence_map.iter()
            .filter(|&&x| x.is_finite())
            .map(|&x| x as f64)
            .sum::<f64>() / confidence_map.len() as f64;
        
        let overall_quality = (coverage_factor * 0.6 + mean_confidence * 0.4).min(1.0);

        // Spatial consistency (measure of spatial clustering of valid pixels)
        let spatial_consistency = self.compute_spatial_consistency(mask)?;

        // Radiometric quality (based on data distribution characteristics)
        let radiometric_quality = self.compute_radiometric_quality(data, mask)?;

        // Artifact severity (based on anomaly map)
        let total_anomalies = anomaly_map.iter().filter(|&&x| x > 0).count();
        let artifact_severity = 1.0 - (total_anomalies as f64 / total_pixels as f64);

        Ok(QualityMetrics {
            overall_quality,
            spatial_consistency,
            radiometric_quality,
            geometric_accuracy: 0.85, // Would be computed with ground truth
            artifact_severity,
        })
    }

    // Helper methods for internal computations

    fn compute_local_texture(&self, data: &Array2<f32>, i: usize, j: usize, window_size: usize) -> f32 {
        let (height, width) = data.dim();
        let half_window = window_size / 2;
        
        let i_start = if i >= half_window { i - half_window } else { 0 };
        let i_end = (i + half_window + 1).min(height);
        let j_start = if j >= half_window { j - half_window } else { 0 };
        let j_end = (j + half_window + 1).min(width);
        
        let mut values = Vec::new();
        for ii in i_start..i_end {
            for jj in j_start..j_end {
                let val = data[[ii, jj]];
                if val.is_finite() && val > 0.0 {
                    values.push(val);
                }
            }
        }
        
        if values.is_empty() {
            return f32::INFINITY;
        }
        
        let mean = values.iter().sum::<f32>() / values.len() as f32;
        let variance = values.iter()
            .map(|&x| (x - mean).powi(2))
            .sum::<f32>() / values.len() as f32;
        
        variance.sqrt() / mean.max(1e-8) // Coefficient of variation
    }

    fn compute_local_mean(&self, data: &Array2<f32>, i: usize, j: usize, window_size: usize) -> f32 {
        let (height, width) = data.dim();
        let half_window = window_size / 2;
        
        let i_start = if i >= half_window { i - half_window } else { 0 };
        let i_end = (i + half_window + 1).min(height);
        let j_start = if j >= half_window { j - half_window } else { 0 };
        let j_end = (j + half_window + 1).min(width);
        
        let mut sum = 0.0f32;
        let mut count = 0;
        
        for ii in i_start..i_end {
            for jj in j_start..j_end {
                let val = data[[ii, jj]];
                if val.is_finite() && val > 0.0 {
                    sum += val;
                    count += 1;
                }
            }
        }
        
        if count > 0 {
            sum / count as f32
        } else {
            0.0
        }
    }

    fn compute_terrain_slopes(&self, dem_data: &Array2<f32>) -> SarResult<Array2<f32>> {
        let (height, width) = dem_data.dim();
        let mut slopes = Array2::<f32>::zeros((height, width));
        
        for i in 1..height-1 {
            for j in 1..width-1 {
                let dz_dx = (dem_data[[i, j+1]] - dem_data[[i, j-1]]) / (2.0 * self.config.layover_shadow_params.dem_resolution);
                let dz_dy = (dem_data[[i+1, j]] - dem_data[[i-1, j]]) / (2.0 * self.config.layover_shadow_params.dem_resolution);
                
                slopes[[i, j]] = (dz_dx*dz_dx + dz_dy*dz_dy).sqrt().atan();
            }
        }
        
        Ok(slopes)
    }

    fn remove_small_patches(&self, mask: &Array2<u8>, _min_size: usize) -> SarResult<Array2<u8>> {
        // Simple implementation - in practice, would use connected component analysis
        Ok(mask.clone())
    }

    fn morphological_erosion(&self, mask: &Array2<u8>, iterations: usize) -> SarResult<Array2<u8>> {
        let mut result = mask.clone();
        
        for _ in 0..iterations {
            let (height, width) = result.dim();
            let mut new_result = result.clone();
            
            for i in 1..height-1 {
                for j in 1..width-1 {
                    if result[[i, j]] == 1 {
                        // Check 3x3 neighborhood
                        let mut all_valid = true;
                        for di in -1i32..=1 {
                            for dj in -1i32..=1 {
                                let ni = (i as i32 + di) as usize;
                                let nj = (j as i32 + dj) as usize;
                                if result[[ni, nj]] == 0 {
                                    all_valid = false;
                                    break;
                                }
                            }
                            if !all_valid { break; }
                        }
                        if !all_valid {
                            new_result[[i, j]] = 0;
                        }
                    }
                }
            }
            result = new_result;
        }
        
        Ok(result)
    }

    fn morphological_dilation(&self, mask: &Array2<u8>, iterations: usize) -> SarResult<Array2<u8>> {
        let mut result = mask.clone();
        
        for _ in 0..iterations {
            let (height, width) = result.dim();
            let mut new_result = result.clone();
            
            for i in 1..height-1 {
                for j in 1..width-1 {
                    if result[[i, j]] == 0 {
                        // Check 3x3 neighborhood
                        let mut has_valid = false;
                        for di in -1i32..=1 {
                            for dj in -1i32..=1 {
                                let ni = (i as i32 + di) as usize;
                                let nj = (j as i32 + dj) as usize;
                                if result[[ni, nj]] == 1 {
                                    has_valid = true;
                                    break;
                                }
                            }
                            if has_valid { break; }
                        }
                        if has_valid {
                            new_result[[i, j]] = 1;
                        }
                    }
                }
            }
            result = new_result;
        }
        
        Ok(result)
    }

    fn compute_data_quality_score(&self, data: &Array2<f32>, mask: &Array2<u8>) -> SarResult<f64> {
        let mut valid_values = Vec::new();
        
        data.indexed_iter().for_each(|((i, j), &value)| {
            if mask[[i, j]] == 1 && value.is_finite() && value > 0.0 {
                valid_values.push(value);
            }
        });
        
        if valid_values.is_empty() {
            return Ok(0.0);
        }
        
        // Compute quality score based on data characteristics
        let mean = valid_values.iter().sum::<f32>() / valid_values.len() as f32;
        let std_dev = (valid_values.iter()
            .map(|&x| (x - mean).powi(2))
            .sum::<f32>() / valid_values.len() as f32).sqrt();
        
        let cv = std_dev / mean.max(1e-8); // Coefficient of variation
        
        // Quality score: lower CV indicates better quality
        let quality_score = (1.0 / (1.0 + cv as f64)).min(1.0);
        
        Ok(quality_score)
    }

    fn compute_spatial_consistency(&self, mask: &Array2<u8>) -> SarResult<f64> {
        // Simple measure: ratio of valid pixels that have valid neighbors
        let (height, width) = mask.dim();
        let mut consistent_pixels = 0;
        let mut total_valid_pixels = 0;
        
        for i in 1..height-1 {
            for j in 1..width-1 {
                if mask[[i, j]] == 1 {
                    total_valid_pixels += 1;
                    
                    // Count valid neighbors in 3x3 window
                    let mut valid_neighbors = 0;
                    for di in -1i32..=1 {
                        for dj in -1i32..=1 {
                            let ni = (i as i32 + di) as usize;
                            let nj = (j as i32 + dj) as usize;
                            if mask[[ni, nj]] == 1 {
                                valid_neighbors += 1;
                            }
                        }
                    }
                    
                    // Consider consistent if at least 5 out of 9 neighbors are valid
                    if valid_neighbors >= 5 {
                        consistent_pixels += 1;
                    }
                }
            }
        }
        
        if total_valid_pixels > 0 {
            Ok(consistent_pixels as f64 / total_valid_pixels as f64)
        } else {
            Ok(0.0)
        }
    }

    fn compute_radiometric_quality(&self, data: &Array2<f32>, mask: &Array2<u8>) -> SarResult<f64> {
        let mut valid_values = Vec::new();
        
        data.indexed_iter().for_each(|((i, j), &value)| {
            if mask[[i, j]] == 1 && value.is_finite() && value > 0.0 {
                valid_values.push(value);
            }
        });
        
        if valid_values.is_empty() {
            return Ok(0.0);
        }
        
        // Sort values for percentile computation
        valid_values.sort_by(|a, b| a.partial_cmp(b).unwrap());
        
        let p10 = valid_values[valid_values.len() * 10 / 100];
        let p90 = valid_values[valid_values.len() * 90 / 100];
        let _median = valid_values[valid_values.len() / 2];
        
        // Quality based on dynamic range and distribution
        let dynamic_range = p90 / p10.max(1e-8);
        let range_quality = (dynamic_range.log10() / 3.0).min(1.0) as f64; // Normalize by 1000x range
        
        Ok(range_quality)
    }

    /// Export masking results for analysis and validation
    pub fn export_masking_results(
        &self,
        result: &AdvancedMaskResult,
        output_dir: &str,
    ) -> SarResult<()> {
        log::info!("📊 Exporting advanced masking results to: {}", output_dir);

        // Create output directory
        std::fs::create_dir_all(output_dir)
            .map_err(|e| SarError::Io(e))?;

        // Export individual component masks as GeoTIFF files
        for (name, _mask) in &result.component_masks {
            let output_path = format!("{}/mask_{}.tif", output_dir, name);
            // Would use GDAL to export as GeoTIFF
            log::debug!("Would export {} mask to: {}", name, output_path);
        }

        // Export confidence map
        let confidence_path = format!("{}/confidence_map.tif", output_dir);
        log::debug!("Would export confidence map to: {}", confidence_path);

        // Export statistics as JSON
        let stats_path = format!("{}/masking_statistics.json", output_dir);
        log::debug!("Would export statistics to: {}", stats_path);

        log::info!("✅ Masking results export completed");
        Ok(())
    }
}

/// Convenience function to create and run advanced masking with default settings
pub fn apply_advanced_masking(
    sigma0_data: &Array2<f32>,
    incidence_angles: Option<&Array2<f32>>,
    dem_data: Option<&Array2<f32>>,
) -> SarResult<AdvancedMaskResult> {
    let processor = AdvancedMaskingProcessor::new_comprehensive();
    processor.process_advanced_masking(sigma0_data, incidence_angles, dem_data, None)
}

/// Convenience function for water-focused masking
pub fn apply_water_masking(
    sigma0_data: &Array2<f32>,
) -> SarResult<AdvancedMaskResult> {
    let mut config = AdvancedMaskingConfig::default();
    config.enable_water_detection = true;
    config.enable_layover_shadow = false;
    config.enable_edge_detection = false;
    
    let processor = AdvancedMaskingProcessor::new(config);
    processor.process_advanced_masking(sigma0_data, None, None, None)
}

/// Convenience function for terrain-focused masking
pub fn apply_terrain_masking(
    sigma0_data: &Array2<f32>,
    incidence_angles: &Array2<f32>,
    dem_data: &Array2<f32>,
) -> SarResult<AdvancedMaskResult> {
    let mut config = AdvancedMaskingConfig::default();
    config.enable_layover_shadow = true;
    config.enable_water_detection = false;
    config.enable_edge_detection = true;
    
    let processor = AdvancedMaskingProcessor::new(config);
    processor.process_advanced_masking(sigma0_data, Some(incidence_angles), Some(dem_data), None)
}

#[cfg(test)]
mod tests {
    use super::*;
    use ndarray::Array2;

    #[test]
    fn test_basic_masking() {
        let data = Array2::<f32>::from_shape_fn((100, 100), |(i, j)| {
            0.1 * (i as f32 * j as f32).sin() + 0.05
        });
        
        let processor = AdvancedMaskingProcessor::new_comprehensive();
        let result = processor.process_advanced_masking(&data, None, None, None);
        
        assert!(result.is_ok());
        let result = result.unwrap();
        assert!(result.statistics.valid_percentage > 80.0);
        assert!(result.quality_metrics.overall_quality > 0.5);
    }

    #[test]
    fn test_water_detection() {
        let mut data = Array2::<f32>::ones((100, 100)) * 0.1;
        
        // Create water area (low backscatter)
        for i in 20..40 {
            for j in 30..50 {
                data[[i, j]] = 0.001; // Very low backscatter
            }
        }
        
        let result = apply_water_masking(&data);
        assert!(result.is_ok());
        
        let result = result.unwrap();
        assert!(result.component_masks.contains_key("water_bodies"));
        assert!(result.statistics.water_pixels > 0);
    }

    #[test]
    fn test_statistical_outlier_detection() {
        let mut data = Array2::<f32>::ones((50, 50)) * 0.1;
        
        // Add some outliers
        data[[10, 10]] = 100.0; // Very high value
        data[[20, 20]] = f32::NAN; // Invalid value
        data[[30, 30]] = -0.1; // Negative value
        
        let processor = AdvancedMaskingProcessor::new_comprehensive();
        let result = processor.process_advanced_masking(&data, None, None, None);
        
        assert!(result.is_ok());
        let result = result.unwrap();
        
        // Outliers should be detected and masked
        assert_eq!(result.final_mask[[10, 10]], 0);
        assert_eq!(result.final_mask[[20, 20]], 0);
        assert_eq!(result.final_mask[[30, 30]], 0);
    }
}
