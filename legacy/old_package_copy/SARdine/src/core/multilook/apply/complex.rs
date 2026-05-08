#![allow(dead_code, unused_variables, deprecated)]
//! Complex multilook for InSAR/PolSAR applications
//!
//! Coherent averaging of complex SLC data, preserving phase information.

use crate::types::{SarError, SarResult};
use ndarray::Array2;
use num_complex::Complex;

use super::types::{BorderMode, MultilookMetadata, MultilookMode};
use super::MultilookProcessor;

impl MultilookProcessor {
    /// Apply multilooking to complex SLC data (coherent averaging for InSAR/PolSAR)
    ///
    /// This method averages the complex I+jQ values directly, preserving phase information.
    /// Critical for interferometry, coherence estimation, and polarimetric processing.
    ///
    /// # Scientific Background
    /// Complex multilooking computes:
    /// ```text
    /// <I+jQ> = (1/N) Σ(I_k + jQ_k)
    /// ```
    /// NOT the intensity multilook <|I+jQ|²>.
    ///
    /// This is essential for:
    /// - **Interferometry**: preserves relative phase between master/slave
    /// - **Coherence**: requires complex averaging of both images before coherence estimation
    /// - **Polarimetry**: covariance/coherency matrix estimation requires complex averaging
    ///
    /// # Power Preservation
    /// For complex data, power preservation scales the averaged complex values by the fill factor:
    /// ```text
    /// <I+jQ>_out = <I+jQ>_mean × sqrt(fill_factor)
    /// ```
    /// This maintains the correct backscatter power while preserving phase.
    ///
    /// # References
    /// - Bamler & Hartl (1998): "Synthetic Aperture Radar Interferometry", Section 3.2
    /// - ESA S1-TN-ESA-GP-0028: "Sentinel-1 Processing" Section 4
    pub fn apply_multilook_complex(
        &self,
        complex_data: &Array2<Complex<f32>>,
        range_spacing: f64,
        azimuth_spacing: f64,
    ) -> SarResult<(Array2<Complex<f32>>, MultilookMetadata)> {
        let start_time = std::time::Instant::now();

        // Validate parameters
        if self.params.azimuth_looks == 0 || self.params.range_looks == 0 {
            return Err(SarError::InvalidParameter("looks must be >= 1".into()));
        }
        let (h, w) = complex_data.dim();
        if h == 0 || w == 0 {
            return Err(SarError::Processing("empty input".into()));
        }

        let la = self.params.azimuth_looks;
        let lr = self.params.range_looks;
        let preserve_power = self.params.preserve_power;

        log::info!("📊 Starting Enhanced Complex Multilook Processing");
        log::info!("   Input: {}×{} pixels", h, w);
        log::info!("   Looks: {}×{} (azimuth×range)", la, lr);
        log::info!("   Mode: Complex (coherent averaging for InSAR/PolSAR)");
        log::info!("   Power preservation: {}", preserve_power);
        log::info!("   Border mode: {:?}", self.params.border_mode);

        // Calculate input power for validation (from complex amplitude)
        let input_power: f64 = complex_data
            .iter()
            .filter(|z| z.re.is_finite() && z.im.is_finite())
            .map(|z| (z.re * z.re + z.im * z.im) as f64)
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

        // Build summed-area tables for real and imaginary parts
        let mut sat_sum_re = vec![0.0f64; (h + 1) * (w + 1)];
        let mut sat_sum_im = vec![0.0f64; (h + 1) * (w + 1)];
        let mut sat_cnt = vec![0u32; (h + 1) * (w + 1)];
        let idx = |yy: usize, xx: usize| -> usize { yy * (w + 1) + xx };

        for y in 0..h {
            let mut row_sum_re = 0.0f64;
            let mut row_sum_im = 0.0f64;
            let mut row_cnt = 0u32;
            for x in 0..w {
                let z = complex_data[[y, x]];
                if z.re.is_finite() && z.im.is_finite() {
                    row_sum_re += z.re as f64;
                    row_sum_im += z.im as f64;
                    row_cnt += 1;
                }
                let above_sum_re = sat_sum_re[idx(y, x + 1)];
                let above_sum_im = sat_sum_im[idx(y, x + 1)];
                let above_cnt = sat_cnt[idx(y, x + 1)];
                sat_sum_re[idx(y + 1, x + 1)] = above_sum_re + row_sum_re;
                sat_sum_im[idx(y + 1, x + 1)] = above_sum_im + row_sum_im;
                sat_cnt[idx(y + 1, x + 1)] = above_cnt + row_cnt;
            }
        }

        let rect = |y0: usize, x0: usize, y1: usize, x1: usize| -> ((f64, f64), u32) {
            let a = idx(y0, x0);
            let b = idx(y0, x1);
            let c = idx(y1, x0);
            let d = idx(y1, x1);
            let sum_re = sat_sum_re[d] + sat_sum_re[a] - sat_sum_re[b] - sat_sum_re[c];
            let sum_im = sat_sum_im[d] + sat_sum_im[a] - sat_sum_im[b] - sat_sum_im[c];
            let cnt = sat_cnt[d] + sat_cnt[a] - sat_cnt[b] - sat_cnt[c];
            ((sum_re, sum_im), cnt)
        };

        let mut out = Array2::<Complex<f32>>::zeros((out_h, out_w));
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

                let ((sum_re, sum_im), cnt) = rect(y0, x0, y1, x1);
                let block_size = (y1 - y0) * (x1 - x0);
                total_samples += block_size;
                total_valid += cnt as usize;

                if cnt > 0 {
                    // Complex average (coherent)
                    let mean_re = sum_re / cnt as f64;
                    let mean_im = sum_im / cnt as f64;

                    if preserve_power {
                        // SCIENTIFIC FIX: Power preservation for complex multilook
                        //
                        // For complex (coherent) multilook, averaging N samples reduces
                        // the expected intensity by factor N (for incoherent/speckle data).
                        //
                        // To preserve the mean intensity per pixel, we scale by √N where
                        // N is the number of valid samples actually averaged (cnt).
                        //
                        // This ensures that for random-phase (speckle) data:
                        //   - Expected output intensity ≈ expected input intensity
                        //   - Total power is preserved when accounting for fewer output pixels
                        //
                        // For coherent data (constant phase), this amplifies the signal,
                        // which is correct behavior - it maintains the power per pixel.
                        //
                        // Reference: Goodman, J.W. (1975) "Statistical Properties of Laser Speckle Patterns"
                        let scale = (cnt as f64).sqrt();
                        out[[oy, ox]] =
                            Complex::new((mean_re * scale) as f32, (mean_im * scale) as f32);
                    } else {
                        out[[oy, ox]] = Complex::new(mean_re as f32, mean_im as f32);
                    }
                } else {
                    out[[oy, ox]] = Complex::new(f32::NAN, f32::NAN);
                    nan_blocks += 1;
                }
            }
        }

        // Calculate output power (from multilooked complex amplitude)
        let output_power: f64 = out
            .iter()
            .filter(|z| z.re.is_finite() && z.im.is_finite())
            .map(|z| (z.re * z.re + z.im * z.im) as f64)
            .sum();
        let output_valid_pixels = out
            .iter()
            .filter(|z| z.re.is_finite() && z.im.is_finite())
            .count();

        // Calculate power ratio
        let power_ratio = if input_power > 0.0 {
            output_power / input_power
        } else {
            f64::NAN
        };

        // Calculate metrics
        let valid_pixel_fraction = output_valid_pixels as f64 / (out_h * out_w) as f64;
        let avg_samples_per_block = total_valid as f64 / (out_h * out_w) as f64;

        // For complex data, estimate ENL from intensity
        let intensity_out = out.map(|z| z.norm());
        let estimated_enl = self.estimate_enl(&intensity_out);

        let new_r = range_spacing * lr as f64;
        let new_a = azimuth_spacing * la as f64;
        let processing_time = start_time.elapsed().as_secs_f64();

        // Comprehensive logging
        log::info!("📊 Complex Multilook Quality Assessment:");
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
                log::warn!("   For complex data, this may indicate phase decorrelation or data quality issues");
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

        log::info!("✅ Enhanced complex multilook complete");

        let metadata = MultilookMetadata {
            input_dims: (h, w),
            output_dims: (out_h, out_w),
            looks: (la, lr),
            input_spacing: (azimuth_spacing, range_spacing),
            output_spacing: (new_a, new_r),
            mode: MultilookMode::Complex,
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
