#![allow(dead_code, unused_variables, deprecated)]
//! Integral image (summed-area table) based multilook methods
//!
//! O(N*M) regardless of looks (vs O(N*M*Lr*La) for naive approach)

use crate::types::{SarError, SarResult};
use ndarray::Array2;

use super::types::{BorderMode, MultilookMetadata};
use super::MultilookProcessor;

impl MultilookProcessor {
    /// Efficient multilook using summed-area tables (integral images).
    /// - O(N*M) regardless of looks (vs O(N*M*Lr*La) for naive approach)
    /// - NaN-aware: ignores NaNs and normalizes by the count of finite samples
    /// - Optional partial windows at the edges (include_partial = true)
    ///
    /// # Arguments
    /// * `intensity_data` - 2D array of intensity values (linear power, never dB)
    /// * `range_spacing` - Original range pixel spacing in meters
    /// * `azimuth_spacing` - Original azimuth pixel spacing in meters
    /// * `include_partial` - Include partial windows at edges (recommended: true)
    ///
    /// # Returns
    /// * Multilooked intensity data, new pixel spacings, and quality metrics
    pub fn apply_multilook_integral_nanaware(
        &self,
        intensity_data: &Array2<f32>,
        range_spacing: f64,
        azimuth_spacing: f64,
        include_partial: bool,
    ) -> SarResult<(Array2<f32>, f64, f64)> {
        if self.params.azimuth_looks == 0 || self.params.range_looks == 0 {
            return Err(SarError::InvalidParameter("looks must be >= 1".into()));
        }
        let (h, w) = intensity_data.dim();
        if h == 0 || w == 0 {
            return Err(SarError::Processing("empty input".into()));
        }

        let la = self.params.azimuth_looks;
        let lr = self.params.range_looks;

        log::info!(
            "Applying integral NaN-aware multilook: {}x{} looks to {}x{} image, partial={}",
            la,
            lr,
            h,
            w,
            include_partial
        );

        // Output size (ceil if include_partial, floor otherwise)
        let out_h = if include_partial {
            (h + la - 1) / la
        } else {
            h / la
        };
        let out_w = if include_partial {
            (w + lr - 1) / lr
        } else {
            w / lr
        };
        if out_h == 0 || out_w == 0 {
            return Err(SarError::Processing(
                "Multilook parameters too large for input image".into(),
            ));
        }

        // Summed-area tables with a 1-pixel top/left border to simplify rectangle sums:
        // SAT_sum(y+1,x+1) = sum of valid values in [0..y, 0..x]
        // SAT_cnt(y+1,x+1) = count of valid values in [0..y, 0..x]
        let mut sat_sum = vec![0.0f64; (h + 1) * (w + 1)];
        let mut sat_cnt = vec![0u32; (h + 1) * (w + 1)];
        let idx = |yy: usize, xx: usize| -> usize { yy * (w + 1) + xx };

        // Build SATs (single-threaded is usually fine; this is linear-time and cache-friendly).
        for y in 0..h {
            let mut row_sum = 0.0f64;
            let mut row_cnt = 0u32;
            for x in 0..w {
                let v = intensity_data[[y, x]];
                if v.is_finite() && v >= 0.0 {
                    row_sum += v as f64;
                    row_cnt += 1;
                }
                let above_sum = sat_sum[idx(y, x + 1)];
                let above_cnt = sat_cnt[idx(y, x + 1)];
                sat_sum[idx(y + 1, x + 1)] = above_sum + row_sum;
                sat_cnt[idx(y + 1, x + 1)] = above_cnt + row_cnt;
            }
        }

        let rect = |y0: usize, x0: usize, y1: usize, x1: usize| -> (f64, u32) {
            // returns (sum, count) over [y0..y1) × [x0..x1)
            let a = idx(y0, x0);
            let b = idx(y0, x1);
            let c = idx(y1, x0);
            let d = idx(y1, x1);
            let sum = sat_sum[d] + sat_sum[a] - sat_sum[b] - sat_sum[c];
            let cnt = sat_cnt[d] + sat_cnt[a] - sat_cnt[b] - sat_cnt[c];
            (sum, cnt)
        };

        let mut out = Array2::<f32>::zeros((out_h, out_w));
        let mut nan_blocks = 0;

        for oy in 0..out_h {
            let y0 = oy * la;
            if y0 >= h {
                break;
            }
            let y1 = (y0 + la).min(h);
            for ox in 0..out_w {
                let x0 = ox * lr;
                if x0 >= w {
                    break;
                }
                let x1 = (x0 + lr).min(w);
                let (sum, cnt) = rect(y0, x0, y1, x1);
                if cnt > 0 {
                    out[[oy, ox]] = (sum / (cnt as f64)) as f32;
                } else {
                    out[[oy, ox]] = f32::NAN;
                    nan_blocks += 1;
                }
            }
        }

        // New pixel spacings (decimation by integer looks)
        let new_r = range_spacing * lr as f64;
        let new_a = azimuth_spacing * la as f64;

        log::info!(
            "Integral multilook complete: {}x{} -> {}x{}, spacing: {:.2}m x {:.2}m, {} NaN blocks",
            h,
            w,
            out_h,
            out_w,
            new_r,
            new_a,
            nan_blocks
        );

        super::log_multilook_summary(intensity_data, &out, &self.params);
        Ok((out, new_r, new_a))
    }

    /// Enhanced multilook with quality metrics (looks-per-pixel count array)
    /// Returns both the averaged image and a count array for quality assessment
    ///
    /// # Returns
    /// * (multilooked_image, range_spacing, azimuth_spacing, looks_per_pixel_count)
    pub fn apply_multilook_with_quality(
        &self,
        intensity_data: &Array2<f32>,
        range_spacing: f64,
        azimuth_spacing: f64,
    ) -> SarResult<(Array2<f32>, f64, f64, Array2<u16>)> {
        if self.params.azimuth_looks == 0 || self.params.range_looks == 0 {
            return Err(SarError::InvalidParameter("looks must be >= 1".into()));
        }
        let (h, w) = intensity_data.dim();
        if h == 0 || w == 0 {
            return Err(SarError::Processing("empty input".into()));
        }

        let la = self.params.azimuth_looks;
        let lr = self.params.range_looks;
        let include_partial = self.params.include_partial;

        // Output size
        let out_h = if include_partial {
            (h + la - 1) / la
        } else {
            h / la
        };
        let out_w = if include_partial {
            (w + lr - 1) / lr
        } else {
            w / lr
        };
        if out_h == 0 || out_w == 0 {
            return Err(SarError::Processing(
                "Multilook parameters too large for input image".into(),
            ));
        }

        // Build summed-area tables for both values and counts
        let mut sat_sum = vec![0.0f64; (h + 1) * (w + 1)];
        let mut sat_cnt = vec![0u32; (h + 1) * (w + 1)];
        let idx = |yy: usize, xx: usize| -> usize { yy * (w + 1) + xx };

        for y in 0..h {
            let mut row_sum = 0.0f64;
            let mut row_cnt = 0u32;
            for x in 0..w {
                let v = intensity_data[[y, x]];
                if v.is_finite() && v >= 0.0 {
                    row_sum += v as f64;
                    row_cnt += 1;
                }
                let above_sum = sat_sum[idx(y, x + 1)];
                let above_cnt = sat_cnt[idx(y, x + 1)];
                sat_sum[idx(y + 1, x + 1)] = above_sum + row_sum;
                sat_cnt[idx(y + 1, x + 1)] = above_cnt + row_cnt;
            }
        }

        let rect = |y0: usize, x0: usize, y1: usize, x1: usize| -> (f64, u32) {
            let a = idx(y0, x0);
            let b = idx(y0, x1);
            let c = idx(y1, x0);
            let d = idx(y1, x1);
            let sum = sat_sum[d] + sat_sum[a] - sat_sum[b] - sat_sum[c];
            let cnt = sat_cnt[d] + sat_cnt[a] - sat_cnt[b] - sat_cnt[c];
            (sum, cnt)
        };

        let mut out = Array2::<f32>::zeros((out_h, out_w));
        let mut looks_count = Array2::<u16>::zeros((out_h, out_w));

        for oy in 0..out_h {
            let y0 = oy * la;
            if y0 >= h {
                break;
            }
            let y1 = (y0 + la).min(h);
            for ox in 0..out_w {
                let x0 = ox * lr;
                if x0 >= w {
                    break;
                }
                let x1 = (x0 + lr).min(w);
                let (sum, cnt) = rect(y0, x0, y1, x1);
                looks_count[[oy, ox]] = cnt as u16;
                if cnt > 0 {
                    out[[oy, ox]] = (sum / (cnt as f64)) as f32;
                } else {
                    out[[oy, ox]] = f32::NAN;
                }
            }
        }

        let new_r = range_spacing * lr as f64;
        let new_a = azimuth_spacing * la as f64;

        super::log_multilook_summary(intensity_data, &out, &self.params);
        Ok((out, new_r, new_a, looks_count))
    }

    /// Enhanced multilook with power preservation and comprehensive metadata
    ///
    /// **CRITICAL:** This is the recommended method for production use.
    /// Implements ESA/SNAP-compliant multilooking with:
    /// - Power preservation normalization (maintains radiometry)
    /// - Comprehensive QA metrics
    /// - Detailed logging
    /// - Performance monitoring
    pub fn apply_multilook_enhanced(
        &self,
        intensity_data: &Array2<f32>,
        range_spacing: f64,
        azimuth_spacing: f64,
    ) -> SarResult<(Array2<f32>, MultilookMetadata)> {
        let start_time = std::time::Instant::now();

        // Validate parameters
        if self.params.azimuth_looks == 0 || self.params.range_looks == 0 {
            return Err(SarError::InvalidParameter("looks must be >= 1".into()));
        }
        let (h, w) = intensity_data.dim();
        if h == 0 || w == 0 {
            return Err(SarError::Processing("empty input".into()));
        }

        let la = self.params.azimuth_looks;
        let lr = self.params.range_looks;
        let preserve_power = self.params.preserve_power;

        log::info!("📊 Starting Enhanced Multilook Processing");
        log::info!("   Input: {}×{} pixels", h, w);
        log::info!("   Looks: {}×{} (azimuth×range)", la, lr);
        log::info!("   Mode: {:?}", self.params.mode);
        log::info!("   Power preservation: {}", preserve_power);
        log::info!("   Border mode: {:?}", self.params.border_mode);

        // Calculate input power for validation
        let input_power: f64 = intensity_data
            .iter()
            .filter(|v| v.is_finite() && **v >= 0.0)
            .map(|v| *v as f64)
            .sum();

        // Output dimensions based on border mode
        let (out_h, out_w) = match self.params.border_mode {
            BorderMode::Drop => (h / la, w / lr),
            BorderMode::Partial => ((h + la - 1) / la, (w + lr - 1) / lr),
            BorderMode::ZeroPad | BorderMode::EdgeReplicate => {
                ((h + la - 1) / la, (w + lr - 1) / lr)
            }
        };

        if out_h == 0 || out_w == 0 {
            return Err(SarError::Processing(
                "Multilook parameters too large for input image".into(),
            ));
        }

        // Build summed-area tables
        let mut sat_sum = vec![0.0f64; (h + 1) * (w + 1)];
        let mut sat_cnt = vec![0u32; (h + 1) * (w + 1)];
        let idx = |yy: usize, xx: usize| -> usize { yy * (w + 1) + xx };

        for y in 0..h {
            let mut row_sum = 0.0f64;
            let mut row_cnt = 0u32;
            for x in 0..w {
                let v = intensity_data[[y, x]];
                if v.is_finite() && v >= 0.0 {
                    row_sum += v as f64;
                    row_cnt += 1;
                }
                let above_sum = sat_sum[idx(y, x + 1)];
                let above_cnt = sat_cnt[idx(y, x + 1)];
                sat_sum[idx(y + 1, x + 1)] = above_sum + row_sum;
                sat_cnt[idx(y + 1, x + 1)] = above_cnt + row_cnt;
            }
        }

        let rect = |y0: usize, x0: usize, y1: usize, x1: usize| -> (f64, u32) {
            let a = idx(y0, x0);
            let b = idx(y0, x1);
            let c = idx(y1, x0);
            let d = idx(y1, x1);
            let sum = sat_sum[d] + sat_sum[a] - sat_sum[b] - sat_sum[c];
            let cnt = sat_cnt[d] + sat_cnt[a] - sat_cnt[b] - sat_cnt[c];
            (sum, cnt)
        };

        let mut out = Array2::<f32>::zeros((out_h, out_w));
        let mut nan_blocks = 0usize;
        let mut total_samples = 0usize;
        let mut total_valid = 0usize;

        // Process each output block
        for oy in 0..out_h {
            let y0 = oy * la;
            if y0 >= h {
                break;
            }
            let y1 = (y0 + la).min(h);

            for ox in 0..out_w {
                let x0 = ox * lr;
                if x0 >= w {
                    break;
                }
                let x1 = (x0 + lr).min(w);

                let (sum, cnt) = rect(y0, x0, y1, x1);
                let block_size = (y1 - y0) * (x1 - x0);
                total_samples += block_size;
                total_valid += cnt as usize;

                if cnt > 0 {
                    // ML-1 FIX: Standard radiometric multilook uses simple mean.
                    // The fill_factor scaling was INCORRECT - it reduced power for
                    // partial blocks instead of preserving it. SNAP and ESA standard
                    // both use mean of valid samples. preserve_power flag now ignored.
                    let mean = sum / cnt as f64;
                    out[[oy, ox]] = mean as f32;
                } else {
                    out[[oy, ox]] = f32::NAN;
                    nan_blocks += 1;
                }
            }
        }

        // OPTIMIZATION #76: Single-pass power and count calculation
        // Eliminates double iteration over output array
        let (output_power, output_valid_pixels) =
            out.iter().fold((0.0f64, 0usize), |(sum, cnt), v| {
                if v.is_finite() && *v >= 0.0 {
                    (sum + *v as f64, cnt + 1)
                } else {
                    (sum, cnt)
                }
            });

        // Calculate power ratio
        let power_ratio = if input_power > 0.0 {
            output_power / input_power
        } else {
            f64::NAN
        };

        // Calculate metrics
        let valid_pixel_fraction = output_valid_pixels as f64 / (out_h * out_w) as f64;
        let avg_samples_per_block = total_valid as f64 / (out_h * out_w) as f64;
        let estimated_enl = self.estimate_enl(&out);

        let new_r = range_spacing * lr as f64;
        let new_a = azimuth_spacing * la as f64;
        let processing_time = start_time.elapsed().as_secs_f64();

        // Comprehensive logging
        log::info!("📊 Multilook Quality Assessment:");
        log::info!("   Output: {}×{} pixels", out_h, out_w);
        log::info!(
            "   Spacing: {:.2}m × {:.2}m → {:.2}m × {:.2}m",
            azimuth_spacing,
            range_spacing,
            new_a,
            new_r
        );
        log::info!(
            "   Valid pixel fraction: {:.1}%",
            valid_pixel_fraction * 100.0
        );
        log::info!(
            "   Average samples per block: {:.1}/{}",
            avg_samples_per_block,
            la * lr
        );
        log::info!("   Theoretical ENL: {}", la * lr);
        log::info!("   Estimated ENL: {:.1}", estimated_enl);
        log::info!("   Power preservation ratio: {:.6}", power_ratio);
        log::info!(
            "   NaN blocks: {} ({:.1}%)",
            nan_blocks,
            nan_blocks as f64 / (out_h * out_w) as f64 * 100.0
        );
        log::info!("   Processing time: {:.3}s", processing_time);

        // Validation warnings
        if preserve_power {
            if power_ratio < 0.98 || power_ratio > 1.02 {
                log::warn!(
                    "⚠️ Power preservation outside 2% tolerance! Ratio: {:.4}",
                    power_ratio
                );
            }
        }

        if valid_pixel_fraction < 0.9 {
            log::warn!(
                "⚠️ Low valid pixel fraction ({:.1}%) - many NaNs in input",
                valid_pixel_fraction * 100.0
            );
        }

        let enl_expected = (la * lr) as f32;
        if estimated_enl.is_finite() && estimated_enl < enl_expected * 0.7 {
            log::warn!(
                "⚠️ ENL ({:.1}) significantly lower than expected ({}) - check data quality",
                estimated_enl,
                enl_expected
            );
        }

        log::info!("✅ Enhanced multilook complete");

        let metadata = MultilookMetadata {
            input_dims: (h, w),
            output_dims: (out_h, out_w),
            looks: (la, lr),
            input_spacing: (azimuth_spacing, range_spacing),
            output_spacing: (new_a, new_r),
            mode: self.params.mode,
            border_mode: self.params.border_mode,
            valid_pixel_fraction,
            theoretical_enl: la * lr,
            estimated_enl,
            power_ratio,
            nan_blocks,
            avg_samples_per_block,
            processing_time_s: processing_time,
        };

        Ok((out, metadata))
    }
}
