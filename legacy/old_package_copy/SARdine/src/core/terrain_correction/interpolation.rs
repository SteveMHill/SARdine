#![allow(dead_code, unused_variables)]
//! SAR image interpolation methods for terrain correction
//!
//! Provides multiple interpolation algorithms for resampling SAR data:
//! - Nearest neighbor (fastest, lowest quality)
//! - Bilinear (good balance of speed and quality)
//! - Bicubic (higher quality, smoother)
//! - Sinc with Hamming window (GAMMA standard)
//! - Lanczos-3 (highest quality)

use ndarray::Array2;
use std::sync::OnceLock;

use super::types::InterpolationMethod;
use super::TerrainCorrector;

// OPTIMIZATION #92: Pre-computed kernel lookup tables
// Sinc-Hamming and Lanczos kernels are expensive (sin, cos per pixel).
// Pre-compute for t in [-3.0, 3.0] at 0.001 resolution (6001 entries).
// Linear interpolation between table entries maintains accuracy.

const LUT_SIZE: usize = 6001;
const LUT_RANGE: f64 = 3.0;
const LUT_SCALE: f64 = (LUT_SIZE as f64 - 1.0) / (2.0 * LUT_RANGE);

/// Pre-computed windowed sinc (Hamming) kernel
static SINC_HAMMING_LUT: OnceLock<[f64; LUT_SIZE]> = OnceLock::new();

/// Pre-computed Lanczos-3 kernel
static LANCZOS3_LUT: OnceLock<[f64; LUT_SIZE]> = OnceLock::new();

/// Initialize sinc-Hamming lookup table
fn init_sinc_hamming_lut() -> [f64; LUT_SIZE] {
    let mut lut = [0.0f64; LUT_SIZE];
    for i in 0..LUT_SIZE {
        let t = (i as f64 / LUT_SCALE) - LUT_RANGE;
        lut[i] = if t.abs() < 1e-6 {
            1.0
        } else {
            let pi_t = std::f64::consts::PI * t;
            let sinc_val = pi_t.sin() / pi_t;
            let hamming = 0.54 + 0.46 * (std::f64::consts::PI * t / LUT_RANGE).cos();
            sinc_val * hamming
        };
    }
    lut
}

/// Initialize Lanczos-3 lookup table  
fn init_lanczos3_lut() -> [f64; LUT_SIZE] {
    let mut lut = [0.0f64; LUT_SIZE];
    let a = 3.0;
    for i in 0..LUT_SIZE {
        let t = (i as f64 / LUT_SCALE) - LUT_RANGE;
        lut[i] = if t.abs() < 1e-6 {
            1.0
        } else if t.abs() >= a {
            0.0
        } else {
            let pi_t = std::f64::consts::PI * t;
            let pi_t_a = std::f64::consts::PI * t / a;
            (pi_t.sin() / pi_t) * (pi_t_a.sin() / pi_t_a)
        };
    }
    lut
}

/// Fast LUT lookup with linear interpolation
#[inline(always)]
fn lut_lookup(lut: &[f64; LUT_SIZE], t: f64) -> f64 {
    // Clamp t to valid range
    let t_clamped = t.clamp(-LUT_RANGE, LUT_RANGE);
    let idx_f = (t_clamped + LUT_RANGE) * LUT_SCALE;
    let idx = idx_f as usize;
    let frac = idx_f - idx as f64;

    if idx + 1 < LUT_SIZE {
        // Linear interpolation between adjacent entries
        lut[idx] * (1.0 - frac) + lut[idx + 1] * frac
    } else {
        lut[idx]
    }
}

impl TerrainCorrector {
    /// Advanced SAR interpolation with multiple methods
    pub(super) fn interpolate_sar_value(
        &self,
        sar_image: &Array2<f32>,
        x: f64,
        y: f64,
        method: InterpolationMethod,
    ) -> f32 {
        match method {
            InterpolationMethod::Nearest => {
                // AUDIT FIX: Check for negative BEFORE casting to usize to prevent wrap-around
                let iy = y.round();
                let jx = x.round();
                if iy < 0.0 || jx < 0.0 {
                    return f32::NAN;
                }
                let i = iy as usize;
                let j = jx as usize;
                if i < sar_image.dim().0 && j < sar_image.dim().1 {
                    sar_image[[i, j]]
                } else {
                    f32::NAN
                }
            }
            InterpolationMethod::Bilinear => self.bilinear_interpolate_unified(sar_image, x, y),
            InterpolationMethod::Bicubic => self.bicubic_interpolate(sar_image, x, y),
            InterpolationMethod::Sinc => self.sinc_interpolate(sar_image, x, y),
            InterpolationMethod::Lanczos => self.lanczos_interpolate(sar_image, x, y),
        }
    }

    /// Bilinear interpolation for SAR data resampling
    ///
    /// CRITICAL: Uses consistent floor() for subpixel fraction to ensure
    /// reproducible results across all coordinate systems.
    ///
    /// # Arguments
    /// * `image` - Input SAR image
    /// * `x` - X coordinate (range direction, column index)
    /// * `y` - Y coordinate (azimuth direction, row index)
    ///
    /// # Returns
    /// * Interpolated value or NaN if out of bounds or any input sample is NaN
    pub(super) fn bilinear_interpolate_unified(&self, image: &Array2<f32>, x: f64, y: f64) -> f32 {
        let (height, width) = image.dim();

        // Early bounds check
        if x < 0.0 || y < 0.0 || x >= (width as f64) || y >= (height as f64) {
            return f32::NAN;
        }

        // Consistent floor() usage (never round for interpolation)
        let x1 = x.floor() as usize;
        let y1 = y.floor() as usize;

        // Ensure we have space for 2x2 interpolation kernel
        if x1 >= width.saturating_sub(1) || y1 >= height.saturating_sub(1) {
            // Edge case: use nearest neighbor
            let safe_x = x1.min(width - 1);
            let safe_y = y1.min(height - 1);
            return image[[safe_y, safe_x]];
        }

        // 2x2 interpolation kernel
        let x2 = x1 + 1;
        let y2 = y1 + 1;

        // Interpolation weights
        let dx = x - x1 as f64;
        let dy = y - y1 as f64;

        // Sample 2x2 kernel
        let v11 = image[[y1, x1]];
        let v12 = image[[y2, x1]];
        let v21 = image[[y1, x2]];
        let v22 = image[[y2, x2]];

        // Deterministic NaN propagation: if ANY sample is NaN, result is NaN
        if v11.is_nan() || v12.is_nan() || v21.is_nan() || v22.is_nan() {
            return f32::NAN;
        }

        // Bilinear interpolation using Horner's method for numerical stability
        let v1 = v11 as f64 + dx * (v21 as f64 - v11 as f64);
        let v2 = v12 as f64 + dx * (v22 as f64 - v12 as f64);
        let result = v1 + dy * (v2 - v1);

        result as f32
    }

    /// Bicubic interpolation for high-quality resampling
    ///
    /// Uses Keys' cubic kernel (a=-0.5) which is interpolating (passes through data points)
    /// and has good frequency response characteristics.
    ///
    /// SCIENTIFIC FIX: Only renormalize weights when pixels are missing (NaN or out-of-bounds).
    /// When all 16 pixels are valid, the Keys kernel already sums to 1.0 and renormalization
    /// would change the interpolation gain.
    pub(super) fn bicubic_interpolate(&self, sar_image: &Array2<f32>, x: f64, y: f64) -> f32 {
        let x0 = x.floor() as i32;
        let y0 = y.floor() as i32;
        let dx = x - x0 as f64;
        let dy = y - y0 as f64;

        // PERFORMANCE: Precompute all 8 kernel values (4 for x, 4 for y) once
        // instead of computing 32 times in the nested loop.
        // Keys' cubic kernel (a = -0.5, standard for image interpolation)
        #[inline(always)]
        fn cubic_kernel(t: f64) -> f64 {
            let t = t.abs();
            if t <= 1.0 {
                // Inner region: (a+2)|t|^3 - (a+3)|t|^2 + 1 with a=-0.5
                1.5 * t * t * t - 2.5 * t * t + 1.0
            } else if t <= 2.0 {
                // Outer region: a|t|^3 - 5a|t|^2 + 8a|t| - 4a with a=-0.5
                -0.5 * t * t * t + 2.5 * t * t - 4.0 * t + 2.0
            } else {
                0.0
            }
        }

        // Precompute kernel weights for all 4 offsets in each dimension
        let weights_x = [
            cubic_kernel(dx - (-1.0)), // j = -1
            cubic_kernel(dx - 0.0),    // j = 0
            cubic_kernel(dx - 1.0),    // j = 1
            cubic_kernel(dx - 2.0),    // j = 2
        ];
        let weights_y = [
            cubic_kernel(dy - (-1.0)), // i = -1
            cubic_kernel(dy - 0.0),    // i = 0
            cubic_kernel(dy - 1.0),    // i = 1
            cubic_kernel(dy - 2.0),    // i = 2
        ];

        let mut sum = 0.0;
        let mut weight_sum = 0.0;
        let mut total_weight = 0.0; // Track theoretical total weight
        let mut missing_pixels = false;

        // Sample 4x4 neighborhood using precomputed weights
        for (i_idx, i) in (-1..3).enumerate() {
            let weight_y = weights_y[i_idx];
            let yi = y0 + i;

            for (j_idx, j) in (-1..3).enumerate() {
                let weight_x = weights_x[j_idx];
                let xi = x0 + j;

                let weight = weight_x * weight_y;
                total_weight += weight;

                if yi >= 0
                    && yi < sar_image.dim().0 as i32
                    && xi >= 0
                    && xi < sar_image.dim().1 as i32
                {
                    let value = sar_image[[yi as usize, xi as usize]];
                    if !value.is_nan() {
                        sum += value as f64 * weight;
                        weight_sum += weight;
                    } else {
                        missing_pixels = true;
                    }
                } else {
                    missing_pixels = true;
                }
            }
        }

        if weight_sum.abs() < 1e-10 {
            return f32::NAN;
        }

        // SCIENTIFIC FIX: Only renormalize if pixels are missing
        // When all pixels are valid, Keys kernel already sums to ~1.0
        // Renormalizing a complete kernel changes the gain incorrectly
        if missing_pixels {
            // Renormalize to compensate for missing pixels
            (sum / weight_sum) as f32
        } else {
            // Use raw sum - kernel is already normalized
            sum as f32
        }
    }

    /// Sinc interpolation using windowed sinc function (GAMMA standard)
    /// OPTIMIZATION #92: Uses pre-computed LUT for 3-5x speedup
    pub(super) fn sinc_interpolate(&self, sar_image: &Array2<f32>, x: f64, y: f64) -> f32 {
        let x0 = x.floor() as i32;
        let y0 = y.floor() as i32;
        let dx = x - x0 as f64;
        let dy = y - y0 as f64;

        // Get or initialize the sinc-Hamming LUT
        let sinc_lut = SINC_HAMMING_LUT.get_or_init(init_sinc_hamming_lut);

        let mut sum = 0.0;
        let mut weight_sum = 0.0;

        // Sample 6x6 neighborhood for sinc interpolation
        for i in -2..4 {
            for j in -2..4 {
                let yi = y0 + i;
                let xi = x0 + j;

                if yi >= 0
                    && yi < sar_image.dim().0 as i32
                    && xi >= 0
                    && xi < sar_image.dim().1 as i32
                {
                    // Use LUT lookup instead of computing sinc
                    let weight_y = lut_lookup(sinc_lut, dy - i as f64);
                    let weight_x = lut_lookup(sinc_lut, dx - j as f64);
                    let weight = weight_x * weight_y;

                    let value = sar_image[[yi as usize, xi as usize]];
                    if !value.is_nan() && weight.abs() > 1e-6 {
                        sum += value as f64 * weight;
                        weight_sum += weight;
                    }
                }
            }
        }

        if weight_sum.abs() > 1e-6 {
            (sum / weight_sum) as f32
        } else {
            f32::NAN
        }
    }

    /// Lanczos interpolation using Lanczos kernel (high quality)
    /// OPTIMIZATION #92: Uses pre-computed LUT for 3-5x speedup
    pub(super) fn lanczos_interpolate(&self, sar_image: &Array2<f32>, x: f64, y: f64) -> f32 {
        let x0 = x.floor() as i32;
        let y0 = y.floor() as i32;
        let dx = x - x0 as f64;
        let dy = y - y0 as f64;

        // Get or initialize the Lanczos-3 LUT
        let lanczos_lut = LANCZOS3_LUT.get_or_init(init_lanczos3_lut);

        let mut sum = 0.0;
        let mut weight_sum = 0.0;

        // Sample 6x6 neighborhood for Lanczos-3
        for i in -2..4 {
            for j in -2..4 {
                let yi = y0 + i;
                let xi = x0 + j;

                if yi >= 0
                    && yi < sar_image.dim().0 as i32
                    && xi >= 0
                    && xi < sar_image.dim().1 as i32
                {
                    // Use LUT lookup instead of computing Lanczos
                    let weight_y = lut_lookup(lanczos_lut, dy - i as f64);
                    let weight_x = lut_lookup(lanczos_lut, dx - j as f64);
                    let weight = weight_x * weight_y;

                    let value = sar_image[[yi as usize, xi as usize]];
                    if !value.is_nan() && weight.abs() > 1e-6 {
                        sum += value as f64 * weight;
                        weight_sum += weight;
                    }
                }
            }
        }

        if weight_sum.abs() > 1e-6 {
            (sum / weight_sum) as f32
        } else {
            f32::NAN
        }
    }
}
