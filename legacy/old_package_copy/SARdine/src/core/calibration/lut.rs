#![allow(dead_code, unused_variables)]
use crate::core::calibration::CalibrationType;
use crate::types::{SarError, SarResult};
use ndarray::parallel::prelude::*;
use ndarray::Array2;

use super::model::CalibrationCoefficients;

fn interpolate_1d(pixels: &[usize], values: &[f32], query: usize) -> f32 {
    if pixels.is_empty() || values.is_empty() {
        return 1.0;
    }
    if query <= pixels[0] {
        return values[0];
    }
    if query >= *pixels.last().unwrap() {
        return *values.last().unwrap();
    }

    // Binary search for surrounding knots
    let mut lo = 0usize;
    let mut hi = pixels.len() - 1;
    while hi - lo > 1 {
        let mid = (hi + lo) / 2;
        if pixels[mid] <= query {
            lo = mid;
        } else {
            hi = mid;
        }
    }

    let x1 = pixels[lo] as f32;
    let x2 = pixels[hi] as f32;
    let y1 = values[lo];
    let y2 = values[hi];
    // SCIENTIFIC FIX (Jan 2026): Use practical tolerance for pixel coordinate comparison.
    // f32::EPSILON (~1.2e-7) is too strict for pixel coordinates (0-25000).
    // Use 1e-4 which is still much smaller than 1 pixel but handles floating point representation.
    if (x2 - x1).abs() < 1e-4 {
        y1
    } else {
        let t = ((query as f32) - x1) / (x2 - x1);
        y1 + t * (y2 - y1)
    }
}

fn interpolate_vector_value(
    vector: &super::model::CalibrationVector,
    pixel: usize,
    cal_type: CalibrationType,
) -> f32 {
    let result = match cal_type {
        CalibrationType::Sigma0 => interpolate_1d(&vector.pixels, &vector.sigma_nought, pixel),
        CalibrationType::Beta0 => {
            let beta_val = interpolate_1d(&vector.pixels, &vector.beta_nought, pixel);
            beta_val
        },
        CalibrationType::Gamma0 => interpolate_1d(&vector.pixels, &vector.gamma, pixel),
        CalibrationType::Dn => 1.0,
    };
    result
}

fn build_dense_gain_luts(
    coeffs: &CalibrationCoefficients,
    image_dims: (usize, usize),
) -> SarResult<(Array2<f32>, Array2<f32>, Array2<f32>)> {
    let (height, width) = image_dims;

    // Build context data needed for interpolation
    let vectors = &coeffs.vectors;
    if vectors.is_empty() {
        return Err(SarError::InvalidInput(
            "No calibration vectors present".into(),
        ));
    }

    // Build azimuth bracket indices and weights for each output line
    let az_brackets: Vec<(usize, usize, f32)> = (0..height)
        .map(|line_idx| {
            let query = line_idx;
            if query <= vectors[0].line as usize {
                (0, 0, 0.0)
            } else if query >= vectors.last().unwrap().line as usize {
                let last = vectors.len() - 1;
                (last, last, 0.0)
            } else {
                let mut lo = 0;
                let mut hi = vectors.len() - 1;
                while hi - lo > 1 {
                    let mid = (hi + lo) / 2;
                    if vectors[mid].line as usize <= query {
                        lo = mid;
                    } else {
                        hi = mid;
                    }
                }
                let line_lo = vectors[lo].line as f32;
                let line_hi = vectors[hi].line as f32;
                let w = if (line_hi - line_lo).abs() < f32::EPSILON {
                    0.0
                } else {
                    ((query as f32) - line_lo) / (line_hi - line_lo)
                };
                
                // Debug line 19 vector selection
                if line_idx == 19 {
                    log::trace!("LINE 19 VECTOR SELECTION: vector_lo={} (line={}), vector_hi={} (line={}), weight={:.6}", 
                        lo, vectors[lo].line, hi, vectors[hi].line, w);
                }
                
                (lo, hi, w)
            }
        })
        .collect();

    // Build range sample points
    let slc_pixels: Vec<usize> = (0..width).collect();

    let valid_ranges = coeffs.valid_sample_ranges.as_ref();

    // Precompute range-interpolated rows for each vector (all cal types in one sweep)
    let mut sigma_rows: Vec<Vec<f32>> = Vec::with_capacity(vectors.len());
    let mut beta_rows: Vec<Vec<f32>> = Vec::with_capacity(vectors.len());
    let mut gamma_rows: Vec<Vec<f32>> = Vec::with_capacity(vectors.len());

    for v in vectors.iter() {
        let mut row_sigma = vec![0.0f32; width];
        let mut row_beta = vec![0.0f32; width];
        let mut row_gamma = vec![0.0f32; width];

        for (col_idx, &slc_pixel) in slc_pixels.iter().enumerate() {
            let k_sigma = interpolate_vector_value(v, slc_pixel, CalibrationType::Sigma0);
            let k_beta = interpolate_vector_value(v, slc_pixel, CalibrationType::Beta0);
            let k_gamma = interpolate_vector_value(v, slc_pixel, CalibrationType::Gamma0);
            
            let to_gain = |k: f32| -> f32 {
                if k > 0.0 && k.is_finite() {
                    let denom = k * k;
                    if denom > f32::EPSILON {
                        1.0 / denom
                    } else {
                        0.0
                    }
                } else {
                    0.0
                }
            };

            row_sigma[col_idx] = to_gain(k_sigma);
            row_beta[col_idx] = to_gain(k_beta);
            row_gamma[col_idx] = to_gain(k_gamma);
        }

        sigma_rows.push(row_sigma);
        beta_rows.push(row_beta);
        gamma_rows.push(row_gamma);
    }

    // Apply absolute calibration constant once to all rows if present
    //
    // SCIENTIFIC NOTE (ESA Sentinel-1):
    // The absoluteCalibrationConstant from Sentinel-1 annotation XML is an additional
    // scaling factor that multiplies the calibration gain. It is typically 1.0 for SLC products.
    //
    // Formula: σ⁰ = |DN|² / K² × abs_const
    // where K is the per-pixel calibration constant from the LUT, and abs_const is this
    // absolute calibration constant (usually 1.0, may differ for different product types).
    //
    // This is applied AFTER computing gain = 1/K², so the final gain becomes: gain × abs_const
    if coeffs.abs_const != 0.0 && coeffs.abs_const.is_finite() {
        let abs_c = coeffs.abs_const as f32;
        for rows in [&mut sigma_rows, &mut beta_rows, &mut gamma_rows] {
            for row in rows.iter_mut() {
                for g in row.iter_mut() {
                    *g *= abs_c;
                }
            }
        }
    }

    let col_block: usize = std::env::var("SARDINE_CAL_LUT_COL_BLOCK")
        .ok()
        .and_then(|v| v.parse().ok())
        .filter(|&v| v > 0)
        .unwrap_or(256);

    let mut sigma = Array2::<f32>::zeros((height, width));
    let mut beta = Array2::<f32>::zeros((height, width));
    let mut gamma = Array2::<f32>::zeros((height, width));

    let valid_spans: Vec<(usize, usize)> = (0..height)
        .map(|row_idx| {
            valid_ranges
                .and_then(|vr| vr.ranges.get(row_idx).copied())
                .unwrap_or((0, width.saturating_sub(1)))
        })
        .collect();
    
    let row_indices: Vec<usize> = (0..height).collect();

    ndarray::Zip::from(sigma.outer_iter_mut())
        .and(beta.outer_iter_mut())
        .and(gamma.outer_iter_mut())
        .and(&az_brackets)
        .and(&valid_spans)
        .and(&row_indices)
        .into_par_iter()
        .for_each(
            |(mut s_row, mut b_row, mut g_row, &(lower, upper, w), &(valid_start, valid_end), &row_idx)| {
                let row_low_s = &sigma_rows[lower];
                let row_up_s = &sigma_rows[upper];
                let row_low_b = &beta_rows[lower];
                let row_up_b = &beta_rows[upper];
                let row_low_g = &gamma_rows[lower];
                let row_up_g = &gamma_rows[upper];

                let mut c = valid_start;
                // SCIENTIFIC FIX (Jan 2026): Use practical tolerance (1e-6) instead of f32::EPSILON.
                // f32::EPSILON (~1.2e-7) is appropriate for comparing values near 1.0, but for
                // interpolation weights (0.0 to 1.0), 1e-6 is more robust against floating point errors.
                if (w - 0.0).abs() < 1e-6 {
                    while c <= valid_end && c < width {
                        let end = (c + col_block).min(valid_end + 1).min(width);
                        for idx in c..end {
                            s_row[idx] = row_low_s[idx];
                            b_row[idx] = row_low_b[idx];
                            g_row[idx] = row_low_g[idx];
                        }
                        c = end;
                    }
                } else if (w - 1.0).abs() < 1e-6 {
                    while c <= valid_end && c < width {
                        let end = (c + col_block).min(valid_end + 1).min(width);
                        for idx in c..end {
                            s_row[idx] = row_up_s[idx];
                            b_row[idx] = row_up_b[idx];
                            g_row[idx] = row_up_g[idx];
                        }
                        c = end;
                    }
                } else {
                    let w0 = 1.0 - w;
                    while c <= valid_end && c < width {
                        let end = (c + col_block).min(valid_end + 1).min(width);
                        for idx in c..end {
                            s_row[idx] = w0 * row_low_s[idx] + w * row_up_s[idx];
                            b_row[idx] = w0 * row_low_b[idx] + w * row_up_b[idx];
                            g_row[idx] = w0 * row_low_g[idx] + w * row_up_g[idx];
                            

                        }
                        c = end;
                    }
                }

                if valid_start > 0 {
                    for idx in 0..valid_start {
                        s_row[idx] = 0.0;
                        b_row[idx] = 0.0;
                        g_row[idx] = 0.0;
                    }
                }
                if valid_end + 1 < width {
                    for idx in (valid_end + 1)..width {
                        s_row[idx] = 0.0;
                        b_row[idx] = 0.0;
                        g_row[idx] = 0.0;
                    }
                }
            },
        );

    Ok((sigma, beta, gamma))
}

pub fn precompute_calibration_lut(
    coeffs: &mut CalibrationCoefficients,
    image_dims: (usize, usize),
) -> SarResult<()> {
    let (sigma, beta, gamma) = build_dense_gain_luts(coeffs, image_dims)?;
    let dn = ndarray::Array2::<f32>::ones((image_dims.0, image_dims.1));

    coeffs.lut = Some(super::model::CalibrationLUT {
        sigma_values: sigma,
        beta_values: beta,
        gamma_values: gamma,
        dn_values: dn,
        is_precomputed: true,
    });
    Ok(())
}

/// Get calibration value at specific pixel coordinates from precomputed LUT.
///
/// Returns `f32::NAN` for out-of-LUT-bounds access (rather than silently applying
/// unity gain) so that downstream processing correctly marks edge pixels as no-data.
/// Returns `0.0` for pixels outside `ValidSampleRanges` (zero-fill from burst edge).
/// Returns `f32::NAN` when no LUT has been precomputed (caller should ensure the LUT
/// is built before invoking calibration).
pub fn get_calibration_value(
    coeffs: &CalibrationCoefficients,
    line: i32,
    pixel: usize,
    cal_type: CalibrationType,
) -> SarResult<f32> {
    if let Some(lut) = &coeffs.lut {
        let line_usize = line.max(0) as usize;
        if line_usize >= lut.sigma_values.shape()[0] || pixel >= lut.sigma_values.shape()[1] {
            // Out-of-bounds: return NaN so the pixel is flagged invalid rather than
            // receiving a meaningless unity gain.
            return Ok(f32::NAN);
        }
        if let Some(ranges) = &coeffs.valid_sample_ranges {
            if let Some((start, end)) = ranges.ranges.get(line_usize) {
                if pixel < *start || pixel > *end {
                    return Ok(0.0);
                }
            }
        }
        let val = match cal_type {
            CalibrationType::Sigma0 => lut.sigma_values[[line_usize, pixel]],
            CalibrationType::Beta0 => lut.beta_values[[line_usize, pixel]],
            CalibrationType::Gamma0 => lut.gamma_values[[line_usize, pixel]],
            CalibrationType::Dn => lut.dn_values[[line_usize, pixel]],
        };
        Ok(val)
    } else {
        // No LUT precomputed: returning NaN forces the pixel to be treated as
        // no-data rather than silently applying unity gain (0 dB), which would
        // produce radiometrically wrong output.
        log::warn!(
            "get_calibration_value called but no LUT is precomputed; \
             call precompute_calibration_lut() before calibration"
        );
        Ok(f32::NAN)
    }
}
