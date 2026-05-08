#![allow(dead_code, unused_variables, deprecated)]
//! Multilook processing implementation
//!
//! This module provides various multilooking algorithms for SAR imagery:
//! - Basic multilook (parallel optimized)
//! - Filtered multilook (boxcar averaging)
//! - Separable multilook (1D cumulative sums)
//! - Integral image multilook (summed-area tables)
//! - Complex multilook (for InSAR/PolSAR)
//! - Enhanced multilook with quality metrics

#![allow(dead_code)]
#![allow(unused_variables)]

mod complex;
mod integral;
mod types;

pub use types::{BorderMode, MultilookMetadata, MultilookMode, MultilookParams};

use crate::types::{SarError, SarResult};
use ndarray::Array2;
use rayon::prelude::*;

/// Multilook processor for reducing speckle in SAR imagery
pub struct MultilookProcessor {
    params: MultilookParams,
}

impl MultilookProcessor {
    /// Create a new multilook processor
    pub fn new(params: MultilookParams) -> Self {
        Self { params }
    }

    /// Create processor with standard parameters
    pub fn standard() -> Self {
        Self::new(MultilookParams::default())
    }

    /// Apply multilooking to intensity data
    ///
    /// # Arguments
    /// * `intensity_data` - 2D array of intensity values (sigma0, beta0, etc.)
    /// * `range_spacing` - Original range pixel spacing in meters
    /// * `azimuth_spacing` - Original azimuth pixel spacing in meters
    ///
    /// # Returns
    /// * Multilooked intensity data and new pixel spacings
    pub fn apply_multilook(
        &self,
        intensity_data: &Array2<f32>,
        range_spacing: f64,
        azimuth_spacing: f64,
    ) -> SarResult<(Array2<f32>, f64, f64)> {
        // Parameter validation - prevent division by zero
        if self.params.azimuth_looks == 0 {
            return Err(SarError::InvalidParameter(
                "azimuth_looks must be >= 1".to_string(),
            ));
        }
        if self.params.range_looks == 0 {
            return Err(SarError::InvalidParameter(
                "range_looks must be >= 1".to_string(),
            ));
        }

        let (rows, cols) = intensity_data.dim();

        log::info!(
            "Applying multilook: {}x{} looks to {}x{} image",
            self.params.azimuth_looks,
            self.params.range_looks,
            rows,
            cols
        );

        // SCIENTIFIC WARNING: Multilooking should be applied to calibrated intensity data.
        // Applying multilook before calibration can produce incorrect results because:
        // 1. Calibration LUTs are defined at native resolution
        // 2. Interpolating LUTs to multilooked grid introduces errors
        // 3. The proper order is: deburst → calibrate → multilook
        //
        // We check for suspicious value ranges that might indicate uncalibrated data.
        let sample_max = intensity_data
            .iter()
            .filter(|v| v.is_finite())
            .take(1000)
            .fold(0.0f32, |max, &v| if v > max { v } else { max });
        if sample_max > 1e10 {
            log::warn!(
                "⚠️  MULTILOOK WARNING: Input data has very large values (max sample: {:.2e}). \
                This may indicate uncalibrated data. Multilooking should be applied AFTER \
                radiometric calibration for scientifically correct results.",
                sample_max
            );
        }

        // Calculate output dimensions
        // Use ceil-style division to retain partial blocks and avoid off-by-one drops
        let out_rows = (rows + self.params.azimuth_looks - 1) / self.params.azimuth_looks;
        let out_cols = (cols + self.params.range_looks - 1) / self.params.range_looks;

        if out_rows == 0 || out_cols == 0 {
            return Err(SarError::Processing(
                "Multilook parameters too large for input image".to_string(),
            ));
        }

        log::info!("Output dimensions: {}x{}", out_rows, out_cols);

        // Create output array
        let mut output = Array2::<f32>::zeros((out_rows, out_cols));
        let az_looks = self.params.azimuth_looks;
        let rg_looks = self.params.range_looks;
        let preserve_power = self.params.preserve_power;
        let input = intensity_data;

        if let Some(out_slice) = output.as_slice_mut() {
            debug_assert_eq!(out_slice.len() % out_cols, 0);
            out_slice
                .par_chunks_exact_mut(out_cols)
                .enumerate()
                .for_each(|(out_row, row_chunk)| {
                    let base_row = out_row * az_looks;
                    for (out_col, cell) in row_chunk.iter_mut().enumerate() {
                        let base_col = out_col * rg_looks;
                        let mut sum = 0.0f64; // Use f64 for precision in accumulation
                        let mut count = 0usize;
                        let block_size = az_looks * rg_looks;

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
                                let val = input[[in_row, in_col]];
                                if val.is_finite() {
                                    sum += val as f64; // Cast f32 to f64 for precision
                                    count += 1;
                                }
                            }
                        }

                        if count > 0 {
                            // Scientific standard: output = mean of valid samples
                            // ML-1 FIX: Removed fill_factor scaling which was incorrectly
                            // REDUCING power for partial blocks. Standard radiometric
                            // multilook always uses simple mean. preserve_power flag is
                            // now deprecated for intensity multilook (no effect).
                            *cell = (sum / count as f64) as f32; // Accumulate in f64, output f32
                        } else {
                            *cell = f32::NAN;
                        }
                    }
                });
        } else {
            // Fallback to sequential processing if the array is not contiguous
            let preserve_power = self.params.preserve_power;
            for out_row in 0..out_rows {
                let base_row = out_row * az_looks;
                for out_col in 0..out_cols {
                    let base_col = out_col * rg_looks;
                    let mut sum = 0.0f64; // Use f64 for precision in accumulation
                    let mut count = 0usize;
                    let block_size = az_looks * rg_looks;

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
                            let val = input[[in_row, in_col]];
                            if val.is_finite() {
                                sum += val as f64; // Cast f32 to f64 for precision
                                count += 1;
                            }
                        }
                    }

                    if count > 0 {
                        // Scientific standard: output = mean of valid samples
                        // ML-1 FIX: preserve_power flag deprecated for intensity multilook
                        output[[out_row, out_col]] = (sum / count as f64) as f32;
                    // Accumulate in f64, output f32
                    } else {
                        output[[out_row, out_col]] = f32::NAN;
                    }
                }
            }
        }

        // Calculate new pixel spacings
        let new_range_spacing = range_spacing * (self.params.range_looks as f64);
        let new_azimuth_spacing = azimuth_spacing * (self.params.azimuth_looks as f64);

        log::info!(
            "Multilooking complete: pixel spacing {}m x {}m -> {}m x {}m",
            range_spacing,
            azimuth_spacing,
            new_range_spacing,
            new_azimuth_spacing
        );

        log_multilook_summary(intensity_data, &output, &self.params);
        Ok((output, new_range_spacing, new_azimuth_spacing))
    }

    /// Apply multilooking with spatial filtering (boxcar average)
    pub fn apply_multilook_filtered(
        &self,
        intensity_data: &Array2<f32>,
        range_spacing: f64,
        azimuth_spacing: f64,
    ) -> SarResult<(Array2<f32>, f64, f64)> {
        if self.params.azimuth_looks == 0 {
            return Err(SarError::InvalidParameter(
                "azimuth_looks must be >= 1".to_string(),
            ));
        }
        if self.params.range_looks == 0 {
            return Err(SarError::InvalidParameter(
                "range_looks must be >= 1".to_string(),
            ));
        }

        let (rows, cols) = intensity_data.dim();

        log::info!(
            "Applying filtered multilook: {}x{} looks to {}x{} image",
            self.params.azimuth_looks,
            self.params.range_looks,
            rows,
            cols
        );

        let out_rows = (rows + self.params.azimuth_looks - 1) / self.params.azimuth_looks;
        let out_cols = (cols + self.params.range_looks - 1) / self.params.range_looks;

        if out_rows == 0 || out_cols == 0 {
            return Err(SarError::Processing(
                "Multilook parameters too large for input image".to_string(),
            ));
        }

        let mut output = Array2::<f32>::zeros((out_rows, out_cols));
        let az_looks = self.params.azimuth_looks;
        let rg_looks = self.params.range_looks;
        let input = intensity_data;

        if let Some(out_slice) = output.as_slice_mut() {
            debug_assert_eq!(out_slice.len() % out_cols, 0);
            out_slice
                .par_chunks_exact_mut(out_cols)
                .enumerate()
                .for_each(|(out_row, row_chunk)| {
                    let start_row = out_row * az_looks;
                    let end_row = ((out_row + 1) * az_looks).min(rows);
                    for (out_col, cell) in row_chunk.iter_mut().enumerate() {
                        let start_col = out_col * rg_looks;
                        let end_col = ((out_col + 1) * rg_looks).min(cols);
                        let mut sum = 0.0f64;
                        let mut count = 0usize;

                        for in_row in start_row..end_row {
                            for in_col in start_col..end_col {
                                let val = input[[in_row, in_col]];
                                if val.is_finite() {
                                    sum += val as f64;
                                    count += 1;
                                }
                            }
                        }

                        if count > 0 {
                            *cell = (sum / count as f64) as f32;
                        } else {
                            *cell = f32::NAN;
                        }
                    }
                });
        } else {
            for out_row in 0..out_rows {
                let start_row = out_row * az_looks;
                let end_row = ((out_row + 1) * az_looks).min(rows);
                for out_col in 0..out_cols {
                    let start_col = out_col * rg_looks;
                    let end_col = ((out_col + 1) * rg_looks).min(cols);
                    let mut sum = 0.0f64;
                    let mut count = 0usize;

                    for in_row in start_row..end_row {
                        for in_col in start_col..end_col {
                            let val = input[[in_row, in_col]];
                            if val.is_finite() {
                                sum += val as f64;
                                count += 1;
                            }
                        }
                    }

                    if count > 0 {
                        output[[out_row, out_col]] = (sum / count as f64) as f32;
                    } else {
                        output[[out_row, out_col]] = f32::NAN;
                    }
                }
            }
        }

        let new_range_spacing = range_spacing * (self.params.range_looks as f64);
        let new_azimuth_spacing = azimuth_spacing * (self.params.azimuth_looks as f64);

        log::info!(
            "Filtered multilooking complete: {}x{} -> {}x{}, spacing: {}m x {}m",
            rows,
            cols,
            out_rows,
            out_cols,
            new_range_spacing,
            new_azimuth_spacing
        );

        log_multilook_summary(intensity_data, &output, &self.params);
        Ok((output, new_range_spacing, new_azimuth_spacing))
    }

    /// High-performance separable boxcar multilook using 1D cumulative sums.
    pub fn apply_multilook_separable(
        &self,
        intensity_data: &Array2<f32>,
        range_spacing: f64,
        azimuth_spacing: f64,
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
        let include_partial = self.params.include_partial;
        let preserve_power = self.params.preserve_power;

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

        // First pass: horizontal sliding sum/count
        let mut horiz_sum = Array2::<f64>::zeros((h, out_w));
        let mut horiz_cnt = Array2::<u32>::zeros((h, out_w));
        for y in 0..h {
            let row = intensity_data.row(y);
            let mut win_sum = 0.0f64;
            let mut win_cnt = 0u32;
            for x in 0..w {
                let v = row[x];
                if v.is_finite() && v >= 0.0 {
                    win_sum += v as f64;
                    win_cnt += 1;
                }
                if x + 1 >= lr {
                    let out_x = x / lr;
                    horiz_sum[[y, out_x]] = win_sum;
                    horiz_cnt[[y, out_x]] = win_cnt;
                    let x_out = x + 1 - lr;
                    let v_out = row[x_out];
                    if v_out.is_finite() && v_out >= 0.0 {
                        win_sum -= v_out as f64;
                        win_cnt = win_cnt.saturating_sub(1);
                    }
                }
            }
        }

        // Second pass: vertical sliding on the horizontally aggregated rows
        // OPTIMIZATION #43: Parallelize across columns - each column is independent
        use rayon::prelude::*;
        let columns: Vec<Vec<f32>> = (0..out_w)
            .into_par_iter()
            .map(|x| {
                let mut col_result = vec![f32::NAN; out_h];
                let mut win_sum = 0.0f64;
                let mut win_cnt = 0u32;
                for y in 0..h {
                    win_sum += horiz_sum[[y, x]];
                    win_cnt += horiz_cnt[[y, x]];
                    if y + 1 >= la {
                        let out_y = y / la;
                        if win_cnt > 0 {
                            // Scientific standard: output = mean of valid samples
                            // ML-1 FIX: preserve_power flag deprecated for intensity multilook
                            col_result[out_y] = (win_sum / win_cnt as f64) as f32;
                        }

                        let y_out = y + 1 - la;
                        win_sum -= horiz_sum[[y_out, x]];
                        win_cnt = win_cnt.saturating_sub(horiz_cnt[[y_out, x]]);
                    }
                }
                col_result
            })
            .collect();

        // Assemble output from parallel column results
        let mut out = Array2::<f32>::from_elem((out_h, out_w), f32::NAN);
        for (x, col) in columns.into_iter().enumerate() {
            for (y, val) in col.into_iter().enumerate() {
                out[[y, x]] = val;
            }
        }

        let new_r = range_spacing * lr as f64;
        let new_a = azimuth_spacing * la as f64;
        log_multilook_summary(intensity_data, &out, &self.params);
        Ok((out, new_r, new_a))
    }

    /// Speckle-aware multilook using a simple Lee-sigma approach.
    pub fn apply_multilook_lee_sigma(
        &self,
        intensity_data: &Array2<f32>,
        range_spacing: f64,
        azimuth_spacing: f64,
        enl_prior: f32,
        sigma_factor: f32,
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
        let out_h = (h + la - 1) / la;
        let out_w = (w + lr - 1) / lr;

        let mut out = Array2::<f32>::zeros((out_h, out_w));

        for oy in 0..out_h {
            let y0 = oy * la;
            let y1 = (y0 + la).min(h);
            for ox in 0..out_w {
                let x0 = ox * lr;
                let x1 = (x0 + lr).min(w);

                let mut sum = 0.0f64;
                let mut sum_sq = 0.0f64;
                let mut n = 0usize;
                for yy in y0..y1 {
                    for xx in x0..x1 {
                        let v = intensity_data[[yy, xx]];
                        if v.is_finite() && v > 0.0 {
                            let vv = v as f64;
                            sum += vv;
                            sum_sq += vv * vv;
                            n += 1;
                        }
                    }
                }
                if n == 0 {
                    out[[oy, ox]] = f32::NAN;
                    continue;
                }

                let mean = sum / n as f64;
                let var = (sum_sq / n as f64) - mean * mean;
                let enl_est = if var > 0.0 && mean > 0.0 {
                    (mean * mean / var).max(1.0)
                } else {
                    enl_prior as f64
                };
                let enl_use = if enl_est.is_finite() {
                    enl_est as f32
                } else {
                    enl_prior
                };

                let sigma = (var.max(0.0)).sqrt();
                let cv = if mean > 0.0 { sigma / mean } else { 0.0 };
                let w = if cv <= 1.0 / (sigma_factor as f64 * enl_use as f64).sqrt() {
                    1.0
                } else {
                    1.0 / (1.0 + (cv * sigma_factor as f64))
                };

                // ML-1 FIX: Use simple mean, no fill_factor scaling
                let original = intensity_data[[y0, x0]].max(0.0) as f64;
                let value = w * mean + (1.0 - w) * original;
                out[[oy, ox]] = value as f32;
            }
        }

        let new_r = range_spacing * lr as f64;
        let new_a = azimuth_spacing * la as f64;
        log_multilook_summary(intensity_data, &out, &self.params);
        Ok((out, new_r, new_a))
    }

    /// Get the theoretical number of looks
    pub fn theoretical_looks(&self) -> usize {
        self.params.range_looks * self.params.azimuth_looks
    }

    /// Calculate equivalent number of looks (ENL) estimate
    pub fn estimate_enl(&self, data: &Array2<f32>) -> f32 {
        let mut m = 0.0f64;
        let mut s = 0.0f64;
        let mut n = 0usize;

        for &v in data.iter() {
            if v.is_finite() && v >= 0.0 {
                let vv = v as f64;
                m += vv;
                s += vv * vv;
                n += 1;
            }
        }
        if n < 2 {
            log::warn!(
                "🚨 SCIENTIFIC WARNING: Insufficient valid samples for ENL estimation (n={})",
                n
            );
            return f32::NAN;
        }

        let mean = m / n as f64;
        let var = (s - (m * m) / n as f64) / ((n - 1) as f64);

        if !mean.is_finite() || !var.is_finite() || mean <= 0.0 {
            log::warn!(
                "⚠️ Invalid statistics for ENL estimation: mean={:.6}, variance={:.6}",
                mean,
                var
            );
            return f32::NAN;
        }

        if var <= 1e-15 {
            log::info!(
                "📊 Uniform data detected (variance={:.2e}), returning high ENL",
                var
            );
            return 1e6;
        }

        let enl = (mean * mean) / var;
        if enl.is_finite() {
            enl.min(1e6) as f32
        } else {
            1e6
        }
    }
}

fn log_multilook_summary(input: &Array2<f32>, output: &Array2<f32>, params: &MultilookParams) {
    let total_out = output.len();
    if total_out == 0 {
        log::warn!("Multilook summary: empty output image");
        return;
    }

    let (input_sum, input_valid) = input.iter().fold((0.0f64, 0usize), |(sum, cnt), v| {
        if v.is_finite() && *v >= 0.0 {
            (sum + *v as f64, cnt + 1)
        } else {
            (sum, cnt)
        }
    });
    let (output_sum, output_valid) = output.iter().fold((0.0f64, 0usize), |(sum, cnt), v| {
        if v.is_finite() && *v >= 0.0 {
            (sum + *v as f64, cnt + 1)
        } else {
            (sum, cnt)
        }
    });

    let nan_blocks = total_out - output_valid;
    let valid_fraction = output_valid as f64 / total_out as f64;
    let input_mean = if input_valid > 0 {
        input_sum / input_valid as f64
    } else {
        f64::NAN
    };
    let output_mean = if output_valid > 0 {
        output_sum / output_valid as f64
    } else {
        f64::NAN
    };
    let power_ratio = if input_mean.is_finite() && input_mean > 0.0 {
        output_mean / input_mean
    } else {
        f64::NAN
    };

    log::info!(
        "Multilook summary: valid {:.1}% ({}/{}), NaN blocks {}, mean {:.4} → {:.4}",
        valid_fraction * 100.0,
        output_valid,
        total_out,
        nan_blocks,
        input_mean,
        output_mean
    );

    if params.preserve_power {
        if power_ratio.is_finite() {
            if power_ratio < 0.98 || power_ratio > 1.02 {
                log::warn!(
                    "Power preservation drift detected: ratio {:.4} (expected ~1.0 with preserve_power=true)",
                    power_ratio
                );
            } else {
                log::info!(
                    "Power ratio {:.4} within ±2% (preserve_power=true)",
                    power_ratio
                );
            }
        } else {
            log::warn!(
                "Power ratio unavailable (insufficient valid samples) with preserve_power=true"
            );
        }
    } else if power_ratio.is_finite() {
        log::debug!("Power ratio {:.4} (preserve_power=false)", power_ratio);
    }
}

#[cfg(test)]
mod tests;
