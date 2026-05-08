//! Mask propagation and seam-filling helpers for deburst outputs.
//!
//! This module provides systematic tracking of pixel validity through the processing pipeline.
//! Quality masks (shadow, layover, invalid pixels) must be propagated from raw data through
//! deburst → calibration → multilook → terrain correction to ensure users can identify
//! scientifically invalid pixels.
//!
//! # Pixel Validity States
//!
//! Pixels can be invalid for multiple reasons, tracked via bitflags:
//! - `VALID` (0x00): Pixel is scientifically valid
//! - `INVALID_SOURCE` (0x01): Source SLC data was NaN/Inf/zero
//! - `NO_COVERAGE` (0x02): No burst data covered this location (deburst gaps)
//! - `CALIBRATION_FAILED` (0x04): Calibration LUT lookup failed or out of bounds
//! - `NOISE_FLOOR` (0x08): Below noise floor after denoising
//! - `SHADOW` (0x10): Radar shadow detected (LIA > threshold)
//! - `LAYOVER` (0x20): Layover detected (surface facing away from radar)
//! - `RTC_CLAMPED_LOW` (0x40): RTC scale factor clamped at minimum
//! - `RTC_CLAMPED_HIGH` (0x80): RTC scale factor clamped at maximum

use ndarray::{Array2, Axis, Zip};
use rayon::prelude::*;

/// Pixel validity flags (bitwise OR to combine multiple reasons)
pub mod flags {
    pub const VALID: u8 = 0x00;
    pub const INVALID_SOURCE: u8 = 0x01;
    pub const NO_COVERAGE: u8 = 0x02;
    pub const CALIBRATION_FAILED: u8 = 0x04;
    pub const NOISE_FLOOR: u8 = 0x08;
    pub const SHADOW: u8 = 0x10;
    pub const LAYOVER: u8 = 0x20;
    pub const RTC_CLAMPED_LOW: u8 = 0x40;
    pub const RTC_CLAMPED_HIGH: u8 = 0x80;
}

/// Comprehensive pixel validity mask that propagates through the pipeline.
///
/// This struct tracks the validity state of each pixel through all processing stages.
/// Use `combine()` to merge masks from different processing stages.
#[derive(Debug, Clone)]
pub struct PixelValidityMask {
    /// Bitwise validity flags for each pixel
    pub flags: Array2<u8>,
    /// Processing stage that created/modified this mask
    pub stage: ProcessingStage,
    /// Statistics about the mask
    pub stats: MaskStats,
}

/// Processing stage identifier for mask provenance
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ProcessingStage {
    /// Initial mask from SLC data
    Source,
    /// After TOPSAR debursting
    Deburst,
    /// After radiometric calibration
    Calibration,
    /// After thermal noise removal
    Denoising,
    /// After multilooking
    Multilook,
    /// After terrain correction
    TerrainCorrection,
    /// Final output
    Output,
}

impl std::fmt::Display for ProcessingStage {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            ProcessingStage::Source => write!(f, "source"),
            ProcessingStage::Deburst => write!(f, "deburst"),
            ProcessingStage::Calibration => write!(f, "calibration"),
            ProcessingStage::Denoising => write!(f, "denoising"),
            ProcessingStage::Multilook => write!(f, "multilook"),
            ProcessingStage::TerrainCorrection => write!(f, "terrain_correction"),
            ProcessingStage::Output => write!(f, "output"),
        }
    }
}

/// Statistics about mask validity
#[derive(Debug, Clone, Default)]
pub struct MaskStats {
    pub total_pixels: usize,
    pub valid_pixels: usize,
    pub invalid_source: usize,
    pub no_coverage: usize,
    pub calibration_failed: usize,
    pub noise_floor: usize,
    pub shadow: usize,
    pub layover: usize,
    pub rtc_clamped_low: usize,
    pub rtc_clamped_high: usize,
}

impl MaskStats {
    /// Compute statistics from a flags array
    pub fn from_flags(flags: &Array2<u8>) -> Self {
        let mut stats = MaskStats {
            total_pixels: flags.len(),
            ..Default::default()
        };

        for &flag in flags.iter() {
            if flag == flags::VALID {
                stats.valid_pixels += 1;
            }
            if flag & flags::INVALID_SOURCE != 0 {
                stats.invalid_source += 1;
            }
            if flag & flags::NO_COVERAGE != 0 {
                stats.no_coverage += 1;
            }
            if flag & flags::CALIBRATION_FAILED != 0 {
                stats.calibration_failed += 1;
            }
            if flag & flags::NOISE_FLOOR != 0 {
                stats.noise_floor += 1;
            }
            if flag & flags::SHADOW != 0 {
                stats.shadow += 1;
            }
            if flag & flags::LAYOVER != 0 {
                stats.layover += 1;
            }
            if flag & flags::RTC_CLAMPED_LOW != 0 {
                stats.rtc_clamped_low += 1;
            }
            if flag & flags::RTC_CLAMPED_HIGH != 0 {
                stats.rtc_clamped_high += 1;
            }
        }

        stats
    }

    /// Get validity percentage
    pub fn validity_percent(&self) -> f64 {
        if self.total_pixels == 0 {
            0.0
        } else {
            100.0 * self.valid_pixels as f64 / self.total_pixels as f64
        }
    }

    /// Log summary
    pub fn log_summary(&self, stage: &str) {
        log::info!(
            "🎭 Mask stats [{}]: {:.2}% valid ({}/{} pixels)",
            stage,
            self.validity_percent(),
            self.valid_pixels,
            self.total_pixels
        );
        if self.invalid_source > 0 {
            log::info!("   - Invalid source: {} pixels", self.invalid_source);
        }
        if self.no_coverage > 0 {
            log::info!("   - No coverage: {} pixels", self.no_coverage);
        }
        if self.calibration_failed > 0 {
            log::info!(
                "   - Calibration failed: {} pixels",
                self.calibration_failed
            );
        }
        if self.noise_floor > 0 {
            log::info!("   - Below noise floor: {} pixels", self.noise_floor);
        }
        if self.shadow > 0 {
            log::info!("   - Radar shadow: {} pixels", self.shadow);
        }
        if self.layover > 0 {
            log::info!("   - Layover: {} pixels", self.layover);
        }
        if self.rtc_clamped_low > 0 {
            log::info!("   - RTC clamped low: {} pixels", self.rtc_clamped_low);
        }
        if self.rtc_clamped_high > 0 {
            log::info!("   - RTC clamped high: {} pixels", self.rtc_clamped_high);
        }
    }
}

impl PixelValidityMask {
    /// Create a new all-valid mask of the given dimensions
    pub fn new_valid(rows: usize, cols: usize, stage: ProcessingStage) -> Self {
        let flags = Array2::from_elem((rows, cols), flags::VALID);
        let stats = MaskStats {
            total_pixels: rows * cols,
            valid_pixels: rows * cols,
            ..Default::default()
        };
        Self {
            flags,
            stage,
            stats,
        }
    }

    /// Create mask from SLC complex data, marking NaN/Inf/zero as invalid
    pub fn from_slc_complex(data: &Array2<num_complex::Complex<f32>>) -> Self {
        let (rows, cols) = data.dim();
        let mut flags = Array2::from_elem((rows, cols), flags::VALID);

        Zip::from(&mut flags).and(data).for_each(|flag, &val| {
            if !val.re.is_finite() || !val.im.is_finite() || (val.re == 0.0 && val.im == 0.0) {
                *flag = flags::INVALID_SOURCE;
            }
        });

        let stats = MaskStats::from_flags(&flags);
        Self {
            flags,
            stage: ProcessingStage::Source,
            stats,
        }
    }

    /// Create mask from intensity/power data
    pub fn from_intensity(data: &Array2<f32>) -> Self {
        let (rows, cols) = data.dim();
        let mut flags = Array2::from_elem((rows, cols), flags::VALID);

        Zip::from(&mut flags).and(data).for_each(|flag, &val| {
            if !val.is_finite() || val <= 0.0 {
                *flag = flags::INVALID_SOURCE;
            }
        });

        let stats = MaskStats::from_flags(&flags);
        Self {
            flags,
            stage: ProcessingStage::Source,
            stats,
        }
    }

    /// Mark pixels as having no coverage (deburst gaps)
    pub fn mark_no_coverage(&mut self, hit_count: &Array2<u16>) {
        Zip::from(&mut self.flags)
            .and(hit_count)
            .for_each(|flag, &hits| {
                if hits == 0 {
                    *flag |= flags::NO_COVERAGE;
                }
            });
        self.stats = MaskStats::from_flags(&self.flags);
        self.stage = ProcessingStage::Deburst;
    }

    /// Mark pixels where calibration failed
    pub fn mark_calibration_failed(&mut self, calibration_valid: &Array2<bool>) {
        Zip::from(&mut self.flags)
            .and(calibration_valid)
            .for_each(|flag, &valid| {
                if !valid {
                    *flag |= flags::CALIBRATION_FAILED;
                }
            });
        self.stats = MaskStats::from_flags(&self.flags);
        self.stage = ProcessingStage::Calibration;
    }

    /// Mark pixels below noise floor
    pub fn mark_noise_floor(&mut self, denoised_data: &Array2<f32>, noise_lut: &Array2<f32>) {
        Zip::from(&mut self.flags)
            .and(denoised_data)
            .and(noise_lut)
            .for_each(|flag, &val, &noise| {
                // If denoised value is less than 10% of the noise level, mark as noise floor
                if val.is_finite() && noise.is_finite() && val < noise * 0.1 {
                    *flag |= flags::NOISE_FLOOR;
                }
            });
        self.stats = MaskStats::from_flags(&self.flags);
        self.stage = ProcessingStage::Denoising;
    }

    /// Mark shadow and layover pixels based on local incidence angle
    ///
    /// # Arguments
    /// * `cos_lia` - Cosine of local incidence angle for each pixel
    /// * `shadow_threshold` - cos(LIA) below this is shadow (default ~cos(85°) = 0.087)
    /// * `layover_threshold` - cos(LIA) indicates layover when surface faces away
    pub fn mark_shadow_layover(
        &mut self,
        cos_lia: &Array2<f32>,
        shadow_threshold: f32,
        layover_threshold: f32,
    ) {
        Zip::from(&mut self.flags)
            .and(cos_lia)
            .for_each(|flag, &cos_val| {
                if cos_val.is_finite() {
                    // Shadow: very grazing incidence (cos near 0)
                    if cos_val < shadow_threshold {
                        *flag |= flags::SHADOW;
                    }
                    // Layover: surface tilted past vertical relative to radar
                    // Indicated by negative cos_lia in some formulations
                    if cos_val < layover_threshold {
                        *flag |= flags::LAYOVER;
                    }
                }
            });
        self.stats = MaskStats::from_flags(&self.flags);
        self.stage = ProcessingStage::TerrainCorrection;
    }

    /// Mark pixels where RTC scale was clamped
    pub fn mark_rtc_clamped(&mut self, rtc_quality_flags: &Array2<u8>) {
        Zip::from(&mut self.flags)
            .and(rtc_quality_flags)
            .for_each(|flag, &rtc_flag| {
                if rtc_flag == 1 {
                    *flag |= flags::RTC_CLAMPED_LOW;
                } else if rtc_flag == 2 {
                    *flag |= flags::RTC_CLAMPED_HIGH;
                }
            });
        self.stats = MaskStats::from_flags(&self.flags);
    }

    /// Combine with another mask (bitwise OR of flags)
    pub fn combine(&mut self, other: &PixelValidityMask) {
        Zip::from(&mut self.flags)
            .and(&other.flags)
            .for_each(|flag, &other_flag| {
                *flag |= other_flag;
            });
        self.stats = MaskStats::from_flags(&self.flags);
    }

    /// Propagate mask through multilooking
    ///
    /// For each output pixel, if ANY input pixel in the block is invalid,
    /// the output pixel is marked with all the invalid flags from that block.
    pub fn propagate_through_multilook(&self, az_looks: usize, rg_looks: usize) -> Self {
        let (rows, cols) = self.flags.dim();
        let out_rows = (rows + az_looks - 1) / az_looks;
        let out_cols = (cols + rg_looks - 1) / rg_looks;

        let mut out_flags = Array2::from_elem((out_rows, out_cols), flags::VALID);

        // Parallel processing by output rows
        out_flags
            .axis_iter_mut(Axis(0))
            .into_par_iter()
            .enumerate()
            .for_each(|(out_row, mut row_view)| {
                let base_row = out_row * az_looks;

                for out_col in 0..out_cols {
                    let base_col = out_col * rg_looks;
                    let mut combined_flag = flags::VALID;

                    for az_offset in 0..az_looks {
                        let in_row = base_row + az_offset;
                        if in_row >= rows {
                            break;
                        }

                        for rg_offset in 0..rg_looks {
                            let in_col = base_col + rg_offset;
                            if in_col >= cols {
                                break;
                            }

                            combined_flag |= self.flags[[in_row, in_col]];
                        }
                    }

                    row_view[out_col] = combined_flag;
                }
            });

        let stats = MaskStats::from_flags(&out_flags);
        stats.log_summary("post-multilook");

        Self {
            flags: out_flags,
            stage: ProcessingStage::Multilook,
            stats,
        }
    }

    /// Get binary mask (1 = valid, 0 = invalid)
    pub fn to_binary_mask(&self) -> Array2<u8> {
        self.flags.mapv(|f| if f == flags::VALID { 1 } else { 0 })
    }

    /// Get boolean mask (true = valid)
    pub fn to_bool_mask(&self) -> Array2<bool> {
        self.flags.mapv(|f| f == flags::VALID)
    }

    /// Check if a specific flag is set for any pixel
    pub fn has_flag(&self, flag: u8) -> bool {
        self.flags.iter().any(|&f| f & flag != 0)
    }

    /// Count pixels with a specific flag
    pub fn count_flag(&self, flag: u8) -> usize {
        self.flags.iter().filter(|&&f| f & flag != 0).count()
    }

    /// Get dimensions
    pub fn dim(&self) -> (usize, usize) {
        self.flags.dim()
    }

    /// Update stage and recompute stats
    pub fn finalize(&mut self, stage: ProcessingStage) {
        self.stage = stage;
        self.stats = MaskStats::from_flags(&self.flags);
        self.stats.log_summary(&stage.to_string());
    }
}

/// Fill small gaps in the mask using morphological operations
pub fn fill_small_gaps(mask: &mut PixelValidityMask, max_gap_size: usize) {
    // Simple gap-filling: if a pixel is surrounded by valid pixels, mark it valid
    // This is a simplified version - production would use proper morphological closing

    let (rows, cols) = mask.flags.dim();
    if rows < 3 || cols < 3 {
        return;
    }

    let kernel_size = max_gap_size.min(5);
    let half = kernel_size / 2;

    let original = mask.flags.clone();

    for row in half..(rows - half) {
        for col in half..(cols - half) {
            if original[[row, col]] != flags::VALID {
                // Count valid neighbors
                let mut valid_count = 0;
                let mut total_count = 0;

                for dr in 0..kernel_size {
                    for dc in 0..kernel_size {
                        let r = row + dr - half;
                        let c = col + dc - half;
                        if r < rows && c < cols {
                            total_count += 1;
                            if original[[r, c]] == flags::VALID {
                                valid_count += 1;
                            }
                        }
                    }
                }

                // If >75% of neighbors are valid, fill this gap
                if valid_count * 4 > total_count * 3 {
                    mask.flags[[row, col]] = flags::VALID;
                }
            }
        }
    }

    mask.stats = MaskStats::from_flags(&mask.flags);
}

#[cfg(test)]
mod tests {
    use super::*;
    use ndarray::Array2;
    use num_complex::Complex;

    #[test]
    fn test_new_valid_mask() {
        let mask = PixelValidityMask::new_valid(100, 200, ProcessingStage::Source);
        assert_eq!(mask.dim(), (100, 200));
        assert_eq!(mask.stats.valid_pixels, 20000);
        assert_eq!(mask.stats.validity_percent(), 100.0);
    }

    #[test]
    fn test_from_slc_complex() {
        let mut data = Array2::from_elem((10, 10), Complex::new(1.0f32, 0.0));
        data[[5, 5]] = Complex::new(f32::NAN, 0.0);
        data[[3, 3]] = Complex::new(0.0, 0.0);

        let mask = PixelValidityMask::from_slc_complex(&data);
        assert_eq!(mask.stats.valid_pixels, 98);
        assert_eq!(mask.stats.invalid_source, 2);
    }

    #[test]
    fn test_propagate_through_multilook() {
        let mut mask = PixelValidityMask::new_valid(10, 10, ProcessingStage::Deburst);
        mask.flags[[0, 0]] = flags::INVALID_SOURCE;
        mask.flags[[2, 2]] = flags::NO_COVERAGE;

        let ml_mask = mask.propagate_through_multilook(2, 2);
        assert_eq!(ml_mask.dim(), (5, 5));
        // Pixel [0,0] in output should have INVALID_SOURCE flag from input [0,0]
        assert!(ml_mask.flags[[0, 0]] & flags::INVALID_SOURCE != 0);
        // Pixel [1,1] in output should have NO_COVERAGE flag from input [2,2]
        assert!(ml_mask.flags[[1, 1]] & flags::NO_COVERAGE != 0);
    }

    #[test]
    fn test_mark_shadow_layover() {
        let mut mask = PixelValidityMask::new_valid(5, 5, ProcessingStage::TerrainCorrection);
        let mut cos_lia = Array2::from_elem((5, 5), 0.5f32);
        cos_lia[[2, 2]] = 0.05; // Shadow (< 0.087 = cos(85°))
        cos_lia[[3, 3]] = -0.1; // Layover

        mask.mark_shadow_layover(&cos_lia, 0.087, 0.0);

        assert!(mask.flags[[2, 2]] & flags::SHADOW != 0);
        assert!(mask.flags[[3, 3]] & flags::LAYOVER != 0);
        assert_eq!(mask.flags[[0, 0]], flags::VALID);
    }
}
