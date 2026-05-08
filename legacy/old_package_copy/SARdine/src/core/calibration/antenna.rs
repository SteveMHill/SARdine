#![allow(dead_code)]
use crate::types::{SarError, SarResult};
use serde::{Deserialize, Serialize};
use serde_json;
use std::fs;
use std::path::{Path, PathBuf};
use std::time::Instant;

use super::model::{AntennaPatternLUT, AntennaPatternVector, CalibrationCoefficients};

#[derive(Serialize, Deserialize)]
struct AntennaCache {
    swath: String,
    polarization: String,
    height: usize,
    width: usize,
    pattern_values: Vec<f32>,
}

fn antenna_cache_path(
    coeffs: &CalibrationCoefficients,
    image_dims: (usize, usize),
    cache_root: &Path,
) -> PathBuf {
    let (h, w) = image_dims;
    let fname = format!(
        "{}_{}_{}x{}.ant.json",
        coeffs.swath.to_lowercase(),
        coeffs.polarization.to_lowercase(),
        h,
        w
    );
    cache_root.join(fname)
}

fn try_load_antenna_cache(
    coeffs: &CalibrationCoefficients,
    image_dims: (usize, usize),
    cache_root: &Path,
) -> Option<AntennaPatternLUT> {
    let path = antenna_cache_path(coeffs, image_dims, cache_root);
    let data = fs::read(&path).ok()?;
    let cached: AntennaCache = serde_json::from_slice(&data).ok()?;
    if cached.height != image_dims.0 || cached.width != image_dims.1 {
        return None;
    }
    if cached.swath != coeffs.swath || cached.polarization != coeffs.polarization {
        return None;
    }
    if cached.pattern_values.len() != cached.height * cached.width {
        return None;
    }
    let pattern_values =
        ndarray::Array2::from_shape_vec((cached.height, cached.width), cached.pattern_values)
            .ok()?;
    Some(AntennaPatternLUT {
        pattern_values,
        is_precomputed: true,
    })
}

fn try_save_antenna_cache(
    coeffs: &CalibrationCoefficients,
    lut: &AntennaPatternLUT,
    cache_root: &Path,
) {
    let path = antenna_cache_path(coeffs, lut.pattern_values.dim(), cache_root);
    if let Some(parent) = path.parent() {
        let _ = fs::create_dir_all(parent);
    }
    let payload = AntennaCache {
        swath: coeffs.swath.clone(),
        polarization: coeffs.polarization.clone(),
        height: lut.pattern_values.nrows(),
        width: lut.pattern_values.ncols(),
        pattern_values: lut.pattern_values.iter().copied().collect(),
    };
    if let Ok(json) = serde_json::to_vec(&payload) {
        let _ = fs::write(path, json);
    }
}

// REMOVED: linspace_indices (Jan 2026 - only used by removed separable antenna pattern code)

fn interp_range(vec: &AntennaPatternVector, slc_pixel: usize) -> f32 {
    if vec.pixels.len() < 2 || vec.values.len() < 2 {
        return 1.0;
    }

    let clamp_val =
        |idx: isize| -> usize { idx.clamp(0, (vec.pixels.len() as isize) - 1) as usize };

    // Locate base segment
    let mut lo = 0usize;
    while lo + 1 < vec.pixels.len() && vec.pixels[lo + 1] <= slc_pixel {
        lo += 1;
    }
    let i1 = lo;
    let i2 = (i1 + 1).min(vec.pixels.len() - 1);
    let i0 = clamp_val(i1 as isize - 1);
    let i3 = clamp_val(i2 as isize + 1);

    let x1 = vec.pixels[i1] as f32;
    let x2 = vec.pixels[i2] as f32;
    if (x2 - x1).abs() < f32::EPSILON {
        return vec.values[i1] as f32;
    }
    let t = ((slc_pixel as f32) - x1) / (x2 - x1);

    catmull_rom(
        vec.values[i0] as f32,
        vec.values[i1] as f32,
        vec.values[i2] as f32,
        vec.values[i3] as f32,
        t.clamp(0.0, 1.0),
    )
    .max(1.0e-6)
}

fn evaluate_pattern(
    vectors: &[AntennaPatternVector],
    az_brackets: &[(usize, usize, f32)],
    slc_pixels: &[usize],
    row_idx: usize,
    col_idx: usize,
) -> f32 {
    let (i1, i2, t_az) = az_brackets[row_idx];
    let i0 = i1.saturating_sub(1);
    let i3 = (i2 + 1).min(vectors.len().saturating_sub(1));
    let slc_pixel = slc_pixels[col_idx];

    let p0 = interp_range(&vectors[i0], slc_pixel);
    let p1 = interp_range(&vectors[i1], slc_pixel);
    let p2 = interp_range(&vectors[i2], slc_pixel);
    let p3 = interp_range(&vectors[i3], slc_pixel);
    let interp = catmull_rom(p0, p1, p2, p3, t_az);
    (interp * interp).max(1.0e-6)
}

/// Parse antenna patterns from XML using the legacy robust parser (exposed via io_xml).
pub fn parse_antenna_pattern_from_xml(
    xml_content: &str,
    range_sampling_rate_hz: Option<f64>,
) -> SarResult<Vec<AntennaPatternVector>> {
    super::io_xml::parse_antenna_pattern_from_xml(xml_content, range_sampling_rate_hz)
}

/// Pre-compute antenna pattern LUT by forwarding to legacy implementation.
pub fn precompute_antenna_pattern_lut(
    coeffs: &mut CalibrationCoefficients,
    image_dims: (usize, usize),
) -> SarResult<()> {
    // Force the separable/range-only path to avoid building a dense 2D LUT.
    // This trims Step B cost and prevents corner clamp artifacts.
    let t_start = Instant::now();
    coeffs.antenna_pattern_lut = None; // Never carry forward a dense LUT

    let (h, w) = image_dims;
    let total = h
        .checked_mul(w)
        .ok_or_else(|| SarError::Processing("Antenna LUT dimension overflow".to_string()))?;
    let result = precompute_antenna_pattern_lut_separable(coeffs, image_dims);

    let mode = if coeffs.antenna_pattern_lut.is_some() {
        "dense"
    } else {
        "none"
    };
    let dense_bytes = coeffs
        .antenna_pattern_lut
        .as_ref()
        .map(|lut| lut.pattern_values.len() * std::mem::size_of::<f32>())
        .unwrap_or(0);
    let dense_mb = dense_bytes as f64 / (1024.0 * 1024.0);
    if dense_bytes > 0 {
        // Defensive: ensure we do not ship a dense LUT; drop it so apply path cannot touch 300M elements.
        log::warn!(
            "⚠️ Dense antenna LUT materialized unexpectedly ({:.1} MiB, {} px); dropping to avoid Step B blow-up",
            dense_mb,
            total
        );
        coeffs.antenna_pattern_lut = None;
    }

    log::info!(
        "📡 Antenna pattern precompute mode={}, dims={}x{}, dense={:.1} MiB ({}) px, elapsed={:.2?}",
        mode,
        h,
        w,
        dense_mb,
        total,
        t_start.elapsed()
    );

    result
}

/// Build a separable antenna pattern (azimuth × range) to avoid full 2D LUT allocation.
/// NOTE: Disabled for S1 IW mode — verified Jan 2026 that S1 calibration vectors
/// are constant in azimuth to within 0.001 dB and do NOT capture TOPSAR scalloping.
/// Scalloping is instead handled at the deburst stage via midpoint selection
/// (DeburstConfig::use_midpoint_selection = true), which always picks data from
/// the burst center where antenna gain is highest.
pub fn precompute_antenna_pattern_lut_separable(
    coeffs: &mut CalibrationCoefficients,
    image_dims: (usize, usize),
) -> SarResult<()> {
    let (height, width) = image_dims;

    // S1 SLC calibration LUTs are azimuth-invariant (sigma0 varies < 0.001 dB across bursts).
    // TOPSAR scalloping is handled by midpoint selection in the deburst stage.
    log::info!(
        "📡 Antenna pattern LUT skipped for {} {} ({}x{}) — scalloping handled by deburst midpoint selection",
        coeffs.swath,
        coeffs.polarization,
        height,
        width
    );

    coeffs.antenna_pattern_lut = None;

    Ok(())
}

#[inline]
fn catmull_rom(p0: f32, p1: f32, p2: f32, p3: f32, t: f32) -> f32 {
    let t2 = t * t;
    let t3 = t2 * t;
    0.5 * ((2.0 * p1)
        + (p2 - p0) * t
        + (2.0 * p0 - 5.0 * p1 + 4.0 * p2 - p3) * t2
        + (3.0 * p1 - p0 - 3.0 * p2 + p3) * t3)
}

fn smooth_gaussian(lut: &ndarray::Array2<f32>) -> ndarray::Array2<f32> {
    let kernel = [0.27901_f32, 0.44198_f32, 0.27901_f32];
    let (h, w) = lut.dim();
    let mut tmp = ndarray::Array2::<f32>::zeros((h, w));
    let mut out = ndarray::Array2::<f32>::zeros((h, w));

    // Horizontal pass
    for r in 0..h {
        for c in 0..w {
            let mut acc = 0.0;
            for k in 0..3 {
                let cc = c.saturating_sub(1).saturating_add(k).min(w - 1);
                acc += kernel[k] * lut[[r, cc]];
            }
            tmp[[r, c]] = acc;
        }
    }

    // Vertical pass
    for r in 0..h {
        for c in 0..w {
            let mut acc = 0.0;
            for k in 0..3 {
                let rr = r.saturating_sub(1).saturating_add(k).min(h - 1);
                acc += kernel[k] * tmp[[rr, c]];
            }
            out[[r, c]] = acc;
        }
    }

    out
}
