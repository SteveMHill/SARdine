/*!
 * Mask Propagation Module
 * 
 * Implements scientifically-valid mask propagation through the entire SAR processing chain:
 * SLC → Deburst → Calibration → Multilook → Terrain Correction
 * 
 * Key Requirements:
 * - Masks must be carried through all processing stages
 * - Invalid samples (DEM voids, noisy borders, water) propagate forward
 * - Conservative approach: once invalid, stays invalid
 * - Integration with Newton-Raphson to avoid instability
 * 
 * References:
 * - SNAP/ISCE mask handling workflows
 * - ESA S1-TN-ESA-GS-10360: "TOPS Processing"
 */

use crate::types::{SarError, SarResult};
use ndarray::Array2;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;

/// Mask stage identifier for tracking provenance
/// Each variant maps to a bit position in the reason_codes bitfield
#[repr(u32)]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum MaskStage {
    /// Initial SLC invalid samples (firstValidSample/lastValidSample)
    SlcBorders = 0,
    /// DEM void regions (NaN, no-data values)
    DemVoids = 1,
    /// Water/ocean regions from auxiliary data
    WaterMask = 2,
    /// Layover/shadow from terrain analysis  
    TerrainArtifacts = 3,
    /// Noise floor detection
    NoiseFloor = 4,
    /// Invalid after deburst (phase discontinuities)
    Deburst = 5,
    /// Invalid after calibration (LUT edge effects)
    Calibration = 6,
    /// Invalid after multilook (edge effects)
    Multilook = 7,
    /// Invalid after terrain correction (geocoding failures)
    TerrainCorrection = 8,
    /// User-provided external mask
    External = 9,
}

/// Mask with full provenance tracking
#[derive(Debug, Clone)]
pub struct ProvenanceMask {
    /// Current mask state (1 = valid, 0 = invalid)
    pub mask: Array2<u8>,
    /// Per-pixel reason codes (bitfield)
    pub reason_codes: Array2<u32>,
    /// Stage-specific masks for debugging
    pub stage_masks: HashMap<MaskStage, Array2<u8>>,
    /// Processing history
    pub history: Vec<String>,
}

impl ProvenanceMask {
    /// Create initial mask from data dimensions
    pub fn new(height: usize, width: usize) -> Self {
        Self {
            mask: Array2::ones((height, width)),
            reason_codes: Array2::zeros((height, width)),
            stage_masks: HashMap::new(),
            history: vec!["Initialized with all valid".to_string()],
        }
    }

    /// Apply a stage-specific mask (logical AND with existing)
    pub fn apply_stage_mask(
        &mut self,
        stage: MaskStage,
        stage_mask: Array2<u8>,
    ) -> SarResult<()> {
        if self.mask.dim() != stage_mask.dim() {
            return Err(SarError::InvalidInput(format!(
                "Mask dimension mismatch: {:?} vs {:?}",
                self.mask.dim(),
                stage_mask.dim()
            )));
        }
        
        let (height, width) = self.mask.dim();
        if height == 0 || width == 0 {
            return Err(SarError::InvalidInput(
                "Cannot apply mask to zero-sized array".to_string()
            ));
        }

        // Logical AND: invalid in either mask → invalid in combined mask
        // Use bitfield position (1 << stage as u32) not raw enum value
        let stage_bit = 1u32 << (stage as u32);
        
        // Use ndarray::Zip for better performance
        ndarray::Zip::from(&mut self.mask)
            .and(&mut self.reason_codes)
            .and(&stage_mask)
            .for_each(|mask_val, reason, &stage_val| {
                if stage_val == 0 {
                    *mask_val = 0;
                    *reason |= stage_bit;
                }
            });

        // Store stage mask for provenance
        self.stage_masks.insert(stage, stage_mask);

        // Recompute valid percentage from current mask state
        let valid_pixels = self.mask.iter().filter(|&&x| x == 1).count();
        let valid_percentage = 100.0 * valid_pixels as f64 / self.mask.len() as f64;
        let history_entry = format!(
            "Applied {:?} mask: {:.1}% valid remaining",
            stage, valid_percentage
        );
        self.history.push(history_entry.clone());
        log::info!("  📋 {}", history_entry);

        Ok(())
    }

    /// Resample mask to new dimensions (used after multilook)
    /// Uses conservative resampling: any invalid pixel in source block → invalid in output
    pub fn resample_to(&self, new_height: usize, new_width: usize) -> SarResult<Self> {
        let (old_height, old_width) = self.mask.dim();

        let mut new_mask = Array2::zeros((new_height, new_width));
        let mut new_reason_codes = Array2::zeros((new_height, new_width));

        let height_ratio = old_height as f64 / new_height as f64;
        let width_ratio = old_width as f64 / new_width as f64;

        // Conservative resampling: invalid if any contributing pixel is invalid
        for i in 0..new_height {
            for j in 0..new_width {
                // Safe index computation with isize and clamping
                let old_i_start = ((i as f64 * height_ratio).floor() as isize)
                    .max(0)
                    .min(old_height as isize - 1) as usize;
                let old_i_end_raw = ((i + 1) as f64 * height_ratio).ceil() as isize;
                let old_i_end = old_i_end_raw
                    .max(old_i_start as isize + 1)  // Guarantee end > start
                    .min(old_height as isize) as usize;
                
                let old_j_start = ((j as f64 * width_ratio).floor() as isize)
                    .max(0)
                    .min(old_width as isize - 1) as usize;
                let old_j_end_raw = ((j + 1) as f64 * width_ratio).ceil() as isize;
                let old_j_end = old_j_end_raw
                    .max(old_j_start as isize + 1)  // Guarantee end > start
                    .min(old_width as isize) as usize;

                let mut all_valid = true;
                let mut combined_reasons = 0u32;

                // Short-circuit opportunity: if all source pixels are 1, skip inner loop
                for old_i in old_i_start..old_i_end {
                    for old_j in old_j_start..old_j_end {
                        if self.mask[[old_i, old_j]] == 0 {
                            all_valid = false;
                            combined_reasons |= self.reason_codes[[old_i, old_j]];
                            // Continue to accumulate all reasons for provenance
                        }
                    }
                }

                new_mask[[i, j]] = if all_valid { 1 } else { 0 };
                new_reason_codes[[i, j]] = combined_reasons;
            }
        }

        let mut history = self.history.clone();
        history.push(format!(
            "Resampled (conservative: any invalid → invalid) from {}x{} to {}x{}",
            old_height, old_width, new_height, new_width
        ));

        Ok(Self {
            mask: new_mask,
            reason_codes: new_reason_codes,
            stage_masks: HashMap::new(), // Stage masks not resampled
            history,
        })
    }

    /// Get statistics about mask validity
    pub fn statistics(&self) -> MaskStatistics {
        let total_pixels = self.mask.len();
        let valid_pixels = self.mask.iter().filter(|&&x| x == 1).count();
        let invalid_pixels = total_pixels - valid_pixels;

        // Count pixels by reason
        let mut reason_counts: HashMap<MaskStage, usize> = HashMap::new();
        for &code in self.reason_codes.iter() {
            if code != 0 {
                // Check each bit
                for stage_val in 0..32 {
                    if code & (1 << stage_val) != 0 {
                        if let Some(stage) = mask_stage_from_u32(stage_val) {
                            *reason_counts.entry(stage).or_insert(0) += 1;
                        }
                    }
                }
            }
        }

        MaskStatistics {
            total_pixels,
            valid_pixels,
            invalid_pixels,
            valid_percentage: 100.0 * valid_pixels as f64 / total_pixels as f64,
            reason_counts,
        }
    }

    /// Export mask provenance report
    pub fn provenance_report(&self) -> String {
        let stats = self.statistics();
        let mut report = String::new();

        report.push_str("╔══════════════════════════════════════════════════════════════╗\n");
        report.push_str("║              MASK PROVENANCE REPORT                          ║\n");
        report.push_str("╚══════════════════════════════════════════════════════════════╝\n\n");

        report.push_str(&format!("Dimensions: {}x{}\n", self.mask.nrows(), self.mask.ncols()));
        report.push_str(&format!(
            "Valid pixels: {} / {} ({:.1}%)\n\n",
            stats.valid_pixels, stats.total_pixels, stats.valid_percentage
        ));

        report.push_str("Processing History:\n");
        for (i, entry) in self.history.iter().enumerate() {
            report.push_str(&format!("  {}. {}\n", i + 1, entry));
        }

        report.push_str("\nInvalidity Reasons:\n");
        for (stage, count) in stats.reason_counts.iter() {
            let percentage = 100.0 * *count as f64 / stats.total_pixels as f64;
            report.push_str(&format!("  {:?}: {} pixels ({:.1}%)\n", stage, count, percentage));
        }

        report
    }
}

/// Mask statistics
#[derive(Debug, Clone)]
pub struct MaskStatistics {
    pub total_pixels: usize,
    pub valid_pixels: usize,
    pub invalid_pixels: usize,
    pub valid_percentage: f64,
    pub reason_counts: HashMap<MaskStage, usize>,
}

/// Convert u32 bit position to MaskStage
fn mask_stage_from_u32(val: u32) -> Option<MaskStage> {
    match val {
        0 => Some(MaskStage::SlcBorders),
        1 => Some(MaskStage::DemVoids),
        2 => Some(MaskStage::WaterMask),
        3 => Some(MaskStage::TerrainArtifacts),
        4 => Some(MaskStage::NoiseFloor),
        5 => Some(MaskStage::Deburst),
        6 => Some(MaskStage::Calibration),
        7 => Some(MaskStage::Multilook),
        8 => Some(MaskStage::TerrainCorrection),
        9 => Some(MaskStage::External),
        _ => None,
    }
}

/// DEM void handler
pub struct DemVoidHandler;

impl DemVoidHandler {
    /// Elevation limits for valid DEM data (meters)
    pub const MIN_VALID_ELEVATION: f32 = -500.0;  // Dead Sea minimum
    pub const MAX_VALID_ELEVATION: f32 = 9000.0;  // Above Everest
    
    /// Detect DEM voids (NaN, no-data values, extreme elevations)
    pub fn detect_voids(
        dem_data: &Array2<f32>,
        no_data_value: Option<f32>,
    ) -> Array2<u8> {
        let (height, width) = dem_data.dim();
        let mut void_mask = Array2::ones((height, width));

        for ((i, j), &elev) in dem_data.indexed_iter() {
            // Mark as void if:
            // 1. Not finite (NaN or Inf)
            let mut is_void = !elev.is_finite();
            
            // 2. Matches no-data value (if provided)
            if let Some(no_data) = no_data_value {
                is_void |= (elev - no_data).abs() < 0.1;
            }
            
            // 3. Extreme elevation (outside configurable bounds)
            is_void |= elev < Self::MIN_VALID_ELEVATION || elev > Self::MAX_VALID_ELEVATION;
            
            if is_void {
                void_mask[[i, j]] = 0;
            }
        }

        let void_count = void_mask.iter().filter(|&&x| x == 0).count();
        let void_percentage = 100.0 * void_count as f64 / void_mask.len() as f64;

        log::info!(
            "🗻 DEM void detection: {:.1}% voids detected",
            void_percentage
        );

        void_mask
    }

    /// Fill DEM voids using simple interpolation (for Newton-Raphson stability)
    pub fn fill_voids_simple(
        dem_data: &Array2<f32>,
        void_mask: &Array2<u8>,
    ) -> SarResult<Array2<f32>> {
        Self::fill_voids_with_fallback(dem_data, void_mask, 0.0)
    }
    
    /// Fill DEM voids with configurable fallback value
    pub fn fill_voids_with_fallback(
        dem_data: &Array2<f32>,
        void_mask: &Array2<u8>,
        fallback_elevation: f32,
    ) -> SarResult<Array2<f32>> {
        let (height, width) = dem_data.dim();
        let mut filled_dem = dem_data.clone();
        let mut no_neighbor_count = 0;

        // Simple nearest-neighbor fill for voids
        for i in 0..height {
            for j in 0..width {
                if void_mask[[i, j]] == 0 {
                    // Find nearest valid elevation
                    let mut search_radius = 1;
                    let mut found = false;

                    while search_radius <= 10 && !found {
                        for di in -(search_radius as isize)..=(search_radius as isize) {
                            for dj in -(search_radius as isize)..=(search_radius as isize) {
                                // Safe index computation with bounds checking
                                let ni_signed = i as isize + di;
                                let nj_signed = j as isize + dj;
                                
                                if ni_signed < 0 || ni_signed >= height as isize ||
                                   nj_signed < 0 || nj_signed >= width as isize {
                                    continue;
                                }
                                
                                let ni = ni_signed as usize;
                                let nj = nj_signed as usize;

                                if void_mask[[ni, nj]] == 1 {
                                    filled_dem[[i, j]] = dem_data[[ni, nj]];
                                    found = true;
                                    break;
                                }
                            }
                            if found {
                                break;
                            }
                        }
                        search_radius += 1;
                    }

                    // If no valid neighbor found, use configurable fallback
                    if !found {
                        filled_dem[[i, j]] = fallback_elevation;
                        no_neighbor_count += 1;
                    }
                }
            }
        }
        
        if no_neighbor_count > 0 {
            log::warn!(
                "🗻 DEM void fill: {} pixels had no valid neighbors within radius 10, using fallback elevation {:.1}m",
                no_neighbor_count, fallback_elevation
            );
        }

        Ok(filled_dem)
    }
}

/// Water/ocean mask integration
pub struct WaterMaskIntegration;

impl WaterMaskIntegration {
    /// Typical water threshold for C-band SAR (Sentinel-1)
    /// Values below this in dB are likely water/ocean
    pub const DEFAULT_WATER_THRESHOLD_DB: f32 = -22.0;
    
    /// Load water mask from auxiliary data (e.g., MODIS water mask)
    pub fn load_external_water_mask(
        _path: &str,
        height: usize,
        width: usize,
    ) -> SarResult<Array2<u8>> {
        // Placeholder: would load from MODIS, OSM, or other source
        log::info!("🌊 External water mask not applied (placeholder - returns all valid)");
        log::info!("   To enable: implement MODIS Water Mask or OSM water layer integration");
        Ok(Array2::ones((height, width)))
    }

    /// Detect water from SAR backscatter characteristics
    /// 
    /// # Arguments
    /// * `sigma0_db` - Sigma0 backscatter in dB
    /// * `threshold_db` - Threshold in dB (typically -22 dB for C-band)
    /// 
    /// # Notes
    /// Default threshold is sensor/mode-specific:
    /// - Sentinel-1 C-band: ~-22 dB (use DEFAULT_WATER_THRESHOLD_DB)
    /// - L-band (ALOS): ~-18 dB  
    /// - X-band: ~-25 dB
    pub fn detect_water_from_backscatter(
        sigma0_db: &Array2<f32>,
        threshold_db: f32,
    ) -> Array2<u8> {
        let (height, width) = sigma0_db.dim();
        let mut water_mask = Array2::ones((height, width));

        for ((i, j), &value) in sigma0_db.indexed_iter() {
            if value.is_finite() && value < threshold_db {
                water_mask[[i, j]] = 0; // Mark as water (invalid for land applications)
            }
        }

        let water_count = water_mask.iter().filter(|&&x| x == 0).count();
        let water_percentage = 100.0 * water_count as f64 / water_mask.len() as f64;

        log::info!(
            "🌊 Water detection (threshold={:.1} dB): {:.1}% water pixels detected",
            threshold_db, water_percentage
        );

        water_mask
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_provenance_mask_creation() {
        let mask = ProvenanceMask::new(100, 100);
        assert_eq!(mask.mask.dim(), (100, 100));
        assert_eq!(mask.mask.iter().filter(|&&x| x == 1).count(), 10000);
    }

    #[test]
    fn test_mask_propagation() {
        let mut mask = ProvenanceMask::new(100, 100);

        // Apply SLC borders mask
        let mut slc_mask = Array2::ones((100, 100));
        for i in 0..10 {
            for j in 0..100 {
                slc_mask[[i, j]] = 0; // Invalid top 10 rows
            }
        }

        mask.apply_stage_mask(MaskStage::SlcBorders, slc_mask)
            .unwrap();

        let stats = mask.statistics();
        assert_eq!(stats.valid_pixels, 9000);
        assert_eq!(stats.invalid_pixels, 1000);
        assert!((stats.valid_percentage - 90.0).abs() < 0.1);
    }

    #[test]
    fn test_dem_void_detection() {
        let mut dem = Array2::ones((50, 50)) * 100.0;

        // Add some voids
        dem[[10, 10]] = f32::NAN;
        dem[[20, 20]] = -32768.0; // No-data value
        dem[[30, 30]] = 10000.0; // Too high

        let void_mask = DemVoidHandler::detect_voids(&dem, Some(-32768.0));

        assert_eq!(void_mask[[10, 10]], 0);
        assert_eq!(void_mask[[20, 20]], 0);
        assert_eq!(void_mask[[30, 30]], 0);
        assert_eq!(void_mask[[0, 0]], 1); // Valid
    }

    #[test]
    fn test_mask_resampling() {
        let mut mask = ProvenanceMask::new(100, 100);

        // Mark some pixels invalid
        for i in 0..10 {
            for j in 0..100 {
                mask.mask[[i, j]] = 0;
            }
        }

        // Resample to 50x50
        let resampled = mask.resample_to(50, 50).unwrap();
        assert_eq!(resampled.mask.dim(), (50, 50));

        // First 5 rows should be invalid (conservative resampling)
        for i in 0..5 {
            assert_eq!(resampled.mask[[i, 0]], 0);
        }
    }
    
    #[test]
    fn test_bitfield_correctness() {
        let mut mask = ProvenanceMask::new(10, 10);
        
        // Apply multiple stage masks
        let mut slc_mask = Array2::ones((10, 10));
        slc_mask[[0, 0]] = 0;
        mask.apply_stage_mask(MaskStage::SlcBorders, slc_mask).unwrap();
        
        let mut dem_mask = Array2::ones((10, 10));
        dem_mask[[0, 0]] = 0; // Same pixel
        dem_mask[[1, 1]] = 0; // Different pixel
        mask.apply_stage_mask(MaskStage::DemVoids, dem_mask).unwrap();
        
        // Check bitfield for pixel [0,0] - should have both bits set
        let reason_00 = mask.reason_codes[[0, 0]];
        assert_eq!(reason_00 & (1u32 << (MaskStage::SlcBorders as u32)), 1u32 << (MaskStage::SlcBorders as u32));
        assert_eq!(reason_00 & (1u32 << (MaskStage::DemVoids as u32)), 1u32 << (MaskStage::DemVoids as u32));
        
        // Check bitfield for pixel [1,1] - should only have DemVoids bit
        let reason_11 = mask.reason_codes[[1, 1]];
        assert_eq!(reason_11 & (1u32 << (MaskStage::SlcBorders as u32)), 0);
        assert_eq!(reason_11 & (1u32 << (MaskStage::DemVoids as u32)), 1u32 << (MaskStage::DemVoids as u32));
        
        // Check valid pixel has no bits set
        assert_eq!(mask.reason_codes[[5, 5]], 0);
    }
    
    #[test]
    fn test_edge_resampling_3x3_to_1x1() {
        let mut mask = ProvenanceMask::new(3, 3);
        
        // Mark center pixel invalid
        mask.mask[[1, 1]] = 0;
        mask.reason_codes[[1, 1]] = 1u32 << (MaskStage::DemVoids as u32);
        
        // Resample to 1x1 - should be invalid (conservative)
        let resampled = mask.resample_to(1, 1).unwrap();
        assert_eq!(resampled.mask[[0, 0]], 0);
        assert_eq!(resampled.reason_codes[[0, 0]], 1u32 << (MaskStage::DemVoids as u32));
    }
    
    #[test]
    fn test_edge_resampling_1xn_to_1x1() {
        let mut mask = ProvenanceMask::new(1, 10);
        
        // Mark some pixels invalid
        mask.mask[[0, 0]] = 0;
        mask.mask[[0, 5]] = 0;
        mask.reason_codes[[0, 0]] = 1u32 << (MaskStage::SlcBorders as u32);
        mask.reason_codes[[0, 5]] = 1u32 << (MaskStage::Calibration as u32);
        
        // Resample to 1x1 - should be invalid with combined reasons
        let resampled = mask.resample_to(1, 1).unwrap();
        assert_eq!(resampled.mask[[0, 0]], 0);
        let expected_reasons = (1u32 << (MaskStage::SlcBorders as u32)) | 
                               (1u32 << (MaskStage::Calibration as u32));
        assert_eq!(resampled.reason_codes[[0, 0]], expected_reasons);
    }
    
    #[test]
    fn test_void_fill_at_borders() {
        // Create DEM with void at border
        let mut dem = Array2::ones((5, 5)) * 100.0;
        dem[[0, 0]] = f32::NAN; // Top-left corner
        
        let void_mask = DemVoidHandler::detect_voids(&dem, None);
        assert_eq!(void_mask[[0, 0]], 0);
        
        let filled = DemVoidHandler::fill_voids_simple(&dem, &void_mask).unwrap();
        // Should be filled from nearest neighbor
        assert!(filled[[0, 0]].is_finite());
        assert!((filled[[0, 0]] - 100.0).abs() < 1.0);
    }
    
    #[test]
    fn test_void_fill_no_neighbor() {
        // Create DEM with all voids
        let mut dem = Array2::from_elem((3, 3), f32::NAN);
        
        let void_mask = DemVoidHandler::detect_voids(&dem, None);
        assert_eq!(void_mask.iter().filter(|&&x| x == 0).count(), 9);
        
        // Fill with custom fallback
        let filled = DemVoidHandler::fill_voids_with_fallback(&dem, &void_mask, 500.0).unwrap();
        // All should be set to fallback
        for val in filled.iter() {
            assert_eq!(*val, 500.0);
        }
    }
    
    #[test]
    fn test_mask_stage_enum_discriminants() {
        // Verify explicit discriminants match expected bit positions
        assert_eq!(MaskStage::SlcBorders as u32, 0);
        assert_eq!(MaskStage::DemVoids as u32, 1);
        assert_eq!(MaskStage::WaterMask as u32, 2);
        assert_eq!(MaskStage::TerrainArtifacts as u32, 3);
        assert_eq!(MaskStage::NoiseFloor as u32, 4);
        assert_eq!(MaskStage::Deburst as u32, 5);
        assert_eq!(MaskStage::Calibration as u32, 6);
        assert_eq!(MaskStage::Multilook as u32, 7);
        assert_eq!(MaskStage::TerrainCorrection as u32, 8);
        assert_eq!(MaskStage::External as u32, 9);
    }
    
    #[test]
    fn test_water_detection_with_threshold() {
        let mut sigma0_db = Array2::ones((10, 10)) * -10.0; // Land
        
        // Add water pixels
        sigma0_db[[0, 0]] = -25.0;
        sigma0_db[[5, 5]] = -30.0;
        
        let water_mask = WaterMaskIntegration::detect_water_from_backscatter(
            &sigma0_db, 
            WaterMaskIntegration::DEFAULT_WATER_THRESHOLD_DB
        );
        
        // Water pixels should be marked invalid
        assert_eq!(water_mask[[0, 0]], 0);
        assert_eq!(water_mask[[5, 5]], 0);
        
        // Land pixels should be valid
        assert_eq!(water_mask[[1, 1]], 1);
    }
    
    #[test]
    fn test_apply_stage_mask_dimension_mismatch() {
        let mut mask = ProvenanceMask::new(10, 10);
        let wrong_size_mask = Array2::ones((5, 5));
        
        let result = mask.apply_stage_mask(MaskStage::DemVoids, wrong_size_mask);
        assert!(result.is_err());
        
        if let Err(SarError::InvalidInput(msg)) = result {
            assert!(msg.contains("dimension mismatch"));
        } else {
            panic!("Expected InvalidInput error for dimension mismatch");
        }
    }
}
