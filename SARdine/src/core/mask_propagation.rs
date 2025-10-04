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
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum MaskStage {
    /// Initial SLC invalid samples (firstValidSample/lastValidSample)
    SlcBorders,
    /// DEM void regions (NaN, no-data values)
    DemVoids,
    /// Water/ocean regions from auxiliary data
    WaterMask,
    /// Layover/shadow from terrain analysis  
    TerrainArtifacts,
    /// Noise floor detection
    NoiseFloor,
    /// Invalid after deburst (phase discontinuities)
    Deburst,
    /// Invalid after calibration (LUT edge effects)
    Calibration,
    /// Invalid after multilook (edge effects)
    Multilook,
    /// Invalid after terrain correction (geocoding failures)
    TerrainCorrection,
    /// User-provided external mask
    External,
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
            return Err(SarError::Processing(format!(
                "Mask dimension mismatch: {:?} vs {:?}",
                self.mask.dim(),
                stage_mask.dim()
            )));
        }

        // Logical AND: invalid in either mask → invalid in combined mask
        let mut invalid_count = 0;
        for ((i, j), &stage_val) in stage_mask.indexed_iter() {
            if stage_val == 0 {
                self.mask[[i, j]] = 0;
                self.reason_codes[[i, j]] |= stage as u32;
                invalid_count += 1;
            }
        }

        // Store stage mask for provenance
        self.stage_masks.insert(stage, stage_mask);

        let valid_percentage = 100.0 * (1.0 - invalid_count as f64 / self.mask.len() as f64);
        let history_entry = format!(
            "Applied {:?} mask: {:.1}% valid remaining",
            stage, valid_percentage
        );
        self.history.push(history_entry.clone());
        log::info!("  📋 {}", history_entry);

        Ok(())
    }

    /// Resample mask to new dimensions (used after multilook)
    pub fn resample_to(&self, new_height: usize, new_width: usize) -> SarResult<Self> {
        let (old_height, old_width) = self.mask.dim();

        let mut new_mask = Array2::zeros((new_height, new_width));
        let mut new_reason_codes = Array2::zeros((new_height, new_width));

        let height_ratio = old_height as f64 / new_height as f64;
        let width_ratio = old_width as f64 / new_width as f64;

        // Conservative resampling: invalid if any contributing pixel is invalid
        for i in 0..new_height {
            for j in 0..new_width {
                let old_i_start = (i as f64 * height_ratio).floor() as usize;
                let old_i_end =
                    ((i + 1) as f64 * height_ratio).ceil().min(old_height as f64) as usize;
                let old_j_start = (j as f64 * width_ratio).floor() as usize;
                let old_j_end =
                    ((j + 1) as f64 * width_ratio).ceil().min(old_width as f64) as usize;

                let mut all_valid = true;
                let mut combined_reasons = 0u32;

                for old_i in old_i_start..old_i_end {
                    for old_j in old_j_start..old_j_end {
                        if self.mask[[old_i, old_j]] == 0 {
                            all_valid = false;
                            combined_reasons |= self.reason_codes[[old_i, old_j]];
                        }
                    }
                }

                new_mask[[i, j]] = if all_valid { 1 } else { 0 };
                new_reason_codes[[i, j]] = combined_reasons;
            }
        }

        let mut history = self.history.clone();
        history.push(format!(
            "Resampled from {}x{} to {}x{}",
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
    /// Detect DEM voids (NaN, no-data values, extreme elevations)
    pub fn detect_voids(
        dem_data: &Array2<f32>,
        no_data_value: Option<f32>,
    ) -> Array2<u8> {
        let (height, width) = dem_data.dim();
        let mut void_mask = Array2::ones((height, width));

        let no_data = no_data_value.unwrap_or(-32768.0);

        for ((i, j), &elev) in dem_data.indexed_iter() {
            // Mark as void if:
            // 1. NaN
            // 2. Matches no-data value
            // 3. Extreme elevation (< -500m or > 9000m)
            if !elev.is_finite()
                || (elev - no_data).abs() < 0.1
                || elev < -500.0
                || elev > 9000.0
            {
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
        let (height, width) = dem_data.dim();
        let mut filled_dem = dem_data.clone();

        // Simple nearest-neighbor fill for voids
        for i in 0..height {
            for j in 0..width {
                if void_mask[[i, j]] == 0 {
                    // Find nearest valid elevation
                    let mut search_radius = 1;
                    let mut found = false;

                    while search_radius <= 10 && !found {
                        for di in -(search_radius as i32)..=(search_radius as i32) {
                            for dj in -(search_radius as i32)..=(search_radius as i32) {
                                let ni = (i as i32 + di) as usize;
                                let nj = (j as i32 + dj) as usize;

                                if ni < height && nj < width && void_mask[[ni, nj]] == 1 {
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

                    // If no valid neighbor found, use sea level
                    if !found {
                        filled_dem[[i, j]] = 0.0;
                    }
                }
            }
        }

        Ok(filled_dem)
    }
}

/// Water/ocean mask integration
pub struct WaterMaskIntegration;

impl WaterMaskIntegration {
    /// Load water mask from auxiliary data (e.g., MODIS water mask)
    pub fn load_external_water_mask(
        _path: &str,
        height: usize,
        width: usize,
    ) -> SarResult<Array2<u8>> {
        // Placeholder: would load from MODIS, OSM, or other source
        log::info!("🌊 Loading external water mask (placeholder)");
        Ok(Array2::ones((height, width)))
    }

    /// Detect water from SAR backscatter characteristics
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
            "🌊 Water detection: {:.1}% water pixels detected",
            water_percentage
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
}
