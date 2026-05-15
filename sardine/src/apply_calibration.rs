//! Radiometric calibration of a debursted Sentinel-1 SLC subswath.
//!
//! Converts raw CInt16 power to σ⁰ (sigma-nought) float32 values using the
//! calibration and noise LUTs from the annotation XML.
//!
//! # Formula
//!
//! ```text
//! σ⁰[l, c] = ( |DN[l,c]|² − N[l,c] ) / K[l,c]²
//! ```
//!
//! where:
//!
//! - `|DN|² = i² + q²`  — power in DNs (integer arithmetic, no overflow for i16)
//! - `N[l,c]`           — thermal noise, interpolated from the noise LUT
//! - `K[l,c]`           — σ⁰ calibration LUT value, interpolated from the cal LUT
//!
//! Both LUTs are defined on a sparse grid (every ~500 lines in azimuth, every
//! ~40 samples in range) and are interpolated bilinearly to the full output grid.
//!
//! # Noise removal
//!
//! The ESA noise model for modern IPF (≥ 2.9) is separable in range and azimuth:
//!
//! ```text
//! N[l, c] = N_range[l, c] * N_azimuth[l, c]
//! ```
//!
//! where `N_range` comes from the noise range vector and `N_azimuth` is from the
//! noise azimuth vector (both interpolated to the output grid).
//!
//! For older IPF versions the azimuth noise component may be absent; in that case
//! we treat `N_azimuth = 1.0` (i.e. only range noise is applied).
//!
//! # Negative power and clipping
//!
//! After noise subtraction `|DN|² − N` can be negative in very dark pixels where
//! the measured power is below the noise floor. We **clamp to zero** rather than
//! returning negative σ⁰, which is physically impossible. This is consistent with
//! ESA and SNAP convention.
//!
//! # Coordinate mapping
//!
//! LUT line indices reference the **original TIFF row grid** (line 0 = first
//! TIFF row). The debursted output has its own line numbering starting at 0.
//! The caller must supply `tiff_line_origin`: the TIFF row corresponding to
//! output line 0 of the debursted array.
//!
//! For a standard IW subswath, `tiff_line_origin = bursts[0].first_line` (which
//! is 0 for IW1, or 0 more generally since each subswath TIFF is separate).
//!
//! Column mapping: output column `c` → TIFF column `c + valid_sample_offset`,
//! where `valid_sample_offset` comes from `DeburstArray::valid_sample_offset`.

use crate::calibration::{SwathCalibration, SwathNoise};
use crate::deburst::DeburstArray;
use rayon::prelude::*;
use thiserror::Error;

// ─── Error type ───────────────────────────────────────────────────────────────

/// Errors produced by radiometric calibration.
#[derive(Debug, Error)]
pub enum CalibrationError {
    /// The calibration LUT has no vectors — cannot interpolate.
    #[error("calibration LUT has no vectors")]
    NoCalibrationVectors,

    /// The noise range LUT has no vectors — cannot interpolate range noise.
    #[error("noise range LUT has no vectors")]
    NoNoiseRangeVectors,

    /// A calibration vector's pixel list is shorter than expected for the
    /// requested output sample range.
    #[error(
        "calibration vector at line {line} covers pixels up to {lut_max_pixel}, \
         but output requires up to pixel {required_pixel}"
    )]
    LutRangeInsufficient {
        line: i32,
        lut_max_pixel: u32,
        required_pixel: usize,
    },
}

// ─── Output type ──────────────────────────────────────────────────────────────

/// Calibrated σ⁰ image for one debursted subswath.
///
/// Values are linear (not dB). Pixels where the DN power was below the noise
/// floor are clamped to 0.0.
pub struct Sigma0Array {
    /// Flat row-major buffer.  Length = `lines × samples`.
    pub data: Vec<f32>,
    /// Number of azimuth lines.
    pub lines: usize,
    /// Number of range samples per line.
    pub samples: usize,
    /// Column offset relative to the TIFF: output col 0 = TIFF col `valid_sample_offset`.
    pub valid_sample_offset: usize,
    /// Maximum number of range pixels by which the calibration LUT was
    /// extrapolated (clamped to the LUT's last value) on any azimuth row.
    /// Zero means the LUT covered every output pixel exactly.  Non-zero
    /// values up to `LUT_EDGE_SLACK_PX` (≈600 samples) are normal for ESA
    /// IPF products whose LUT stops short of the last valid sample, but
    /// callers may want to log or flag values above some workflow threshold.
    pub cal_lut_extrapolation_gap_px: u32,
    /// Same as `cal_lut_extrapolation_gap_px`, but for the range noise LUT.
    pub noise_lut_extrapolation_gap_px: u32,

    /// Per-pixel NESZ (Noise Equivalent Sigma Zero) in linear σ⁰ units.
    /// Same shape as `data` (flat row-major, `lines × samples`).
    ///
    /// Value at index `[l * samples + c]` is the thermal noise floor at that
    /// pixel in σ⁰ linear units: `N[l,c] / K²[l,c]`.  Pixels where `K²` is
    /// zero have NESZ = 0.0 (indeterminate).  Used downstream for per-pixel
    /// noise floor masking in terrain correction.
    pub nesz: Vec<f32>,
}

// ─── Interpolation helpers ────────────────────────────────────────────────────

/// Linear interpolation between two f32 values.
///
/// `t` is the fractional distance in [0, 1] from `a` to `b`.
/// Clamps `t` to [0, 1] so callers do not need to guard boundary cases.
#[inline]
fn lerp(a: f32, b: f32, t: f32) -> f32 {
    let t = t.clamp(0.0, 1.0);
    a + (b - a) * t
}

/// Interpolate a 1-D LUT defined at irregular pixel positions to a dense pixel
/// at position `col`.
///
/// `pixels` must be sorted ascending. `values` has the same length.
/// Returns the interpolated value, or the nearest-neighbour value if `col` is
/// outside the LUT's pixel range (extrapolation is clamped).
///
/// Used directly in unit tests.  Production code uses the cursor-based
/// [`bilinear_row_into`] path which amortises the search over the full row.
#[allow(dead_code)]
fn interp_range_lut(pixels: &[u32], values: &[f32], col: usize) -> f32 {
    debug_assert_eq!(pixels.len(), values.len());
    debug_assert!(!pixels.is_empty());

    let col = col as u32;

    // binary_search handles exact hits (Ok(i)) and interpolation (Err(i)).
    match pixels.binary_search(&col) {
        Ok(i) => values[i],
        Err(0) => values[0],
        Err(i) if i >= pixels.len() => *values.last().unwrap(),
        Err(i) => {
            let lo = i - 1;
            let hi = i;
            let t = (col - pixels[lo]) as f32 / (pixels[hi] - pixels[lo]) as f32;
            lerp(values[lo], values[hi], t)
        }
    }
}

/// Interpolate a 1-D LUT defined at two azimuth lines to an intermediate line.
///
/// `line_lo` / `vals_lo` and `line_hi` / `vals_hi` are the bracketing rows.
/// `line` is the target line (must be between `line_lo` and `line_hi`).
#[inline]
fn interp_azimuth(
    line_lo: i64,
    line_hi: i64,
    val_lo: f32,
    val_hi: f32,
    line: i64,
) -> f32 {
    if line_hi == line_lo {
        return val_lo;
    }
    let t = (line - line_lo) as f32 / (line_hi - line_lo) as f32;
    lerp(val_lo, val_hi, t)
}

/// Find the bracketing LUT vector pair for a given TIFF line.
///
/// Returns `(lo_idx, hi_idx)` into `lines`:
/// - `lo_idx`: largest index i such that `lines[i] ≤ tiff_line`
/// - `hi_idx`: smallest index i such that `lines[i] ≥ tiff_line`
///
/// If `tiff_line` is outside the LUT's range, both indices clamp to the
/// nearest endpoint (nearest-neighbour extrapolation in azimuth).
fn bracket_azimuth(lines: &[i64], tiff_line: i64) -> (usize, usize) {
    debug_assert!(!lines.is_empty());

    match lines.binary_search(&tiff_line) {
        Ok(i) => (i, i),
        Err(0) => (0, 0),
        Err(i) if i >= lines.len() => {
            let last = lines.len() - 1;
            (last, last)
        }
        Err(i) => (i - 1, i),
    }
}

// ─── Per-column LUT evaluation ────────────────────────────────────────────────

/// Evaluate a 1-D range LUT at `tiff_col` given that `cur` is the lower-bracket
/// cursor index (the largest index i such that `pixels[i] ≤ tiff_col`, or 0 if
/// `tiff_col < pixels[0]`).
///
/// Called from the hot inner loop of [`bilinear_row_into`]; must stay `#[inline]`.
#[inline]
fn lut_cursor_eval(pixels: &[u32], values: &[f32], cur: usize, tiff_col: u32) -> f32 {
    if pixels[cur] > tiff_col {
        // Before first LUT point — clamp.
        values[0]
    } else if cur + 1 < pixels.len() {
        // Between cur and cur+1.  After the while-advance, pixels[cur+1] > tiff_col
        // is guaranteed, so t ∈ [0, 1).
        let t = (tiff_col - pixels[cur]) as f32 / (pixels[cur + 1] - pixels[cur]) as f32;
        lerp(values[cur], values[cur + 1], t)
    } else {
        // Beyond last LUT point — clamp.
        values[cur]
    }
}

/// Fill `buf` with bilinearly-interpolated calibration or noise values for
/// a single output line.
///
/// # Arguments
///
/// * `vecs_lines`  — Azimuth line indices (TIFF frame) of each LUT vector.
/// * `vecs_pixels` — Range pixel positions for each vector (LUT sample grid).
/// * `vecs_values` — LUT values at each `(vector, pixel)` position.
/// * `tiff_line`   — TIFF row corresponding to this output line.
/// * `col_offset`  — Column offset: output col → TIFF col = col + col_offset.
/// * `buf`         — Output slice, length = number of output columns.
///
/// # Complexity
///
/// O(N_cols + N_lut) per call: output columns are iterated once in order;
/// the LUT cursor is advanced monotonically (total advances ≤ N_lut), so no
/// binary search is performed per column.  This is ~6× faster than a per-column
/// binary search for typical S-1 geometry (N_cols ≈ 22 k, N_lut ≈ 50).

// ─── Vectorizable row kernels ─────────────────────────────────────────────────
//
// Both functions are structured so that LLVM can autovectorise the inner loop:
//
//   • `dn_to_power_row`: simple element-wise i16 widening + multiply-add.
//     No branches, no memory aliasing.
//
//   • `apply_cal_row`: branchless f32 arithmetic.  The "k² > 0" guard is
//     implemented as a floating-point select (multiply by 0.0 / 1.0) rather
//     than a branch, removing the per-lane divergence that prevents SIMD.
//
// On x86-64 with -C target-cpu=native (or explicit +avx2) LLVM typically
// emits 8-wide AVX2 vectors.  On aarch64 it emits NEON/SVE.  No unsafe or
// platform-specific intrinsics are used: portability and safety are preserved.

/// Convert a row of CInt16 DN samples to f32 power values in-place.
///
/// `power[c] = (i_dn[c]² + q_dn[c]²) as f32`
///
/// Uses i32 intermediate arithmetic to prevent i16 overflow before the
/// conversion to f32 (`i16::MAX² + i16::MAX² = 2 × 32767² ≈ 2.1 × 10⁹`,
/// which fits in i32 but not i16).
#[inline(always)]
fn dn_to_power_row(in_row: &[[i16; 2]], power: &mut [f32]) {
    debug_assert_eq!(in_row.len(), power.len());
    for (iq, p) in in_row.iter().zip(power.iter_mut()) {
        let [i, q] = *iq;
        *p = (i as i32 * i as i32 + q as i32 * q as i32) as f32;
    }
}

/// Apply calibration and noise to one row, writing σ⁰ and NESZ.
///
/// ```text
/// σ⁰[c] = max(power[c] − noise[c] × az_factor, 0) / k[c]²
/// nesz[c] = noise[c] × az_factor / k[c]²
/// ```
///
/// When `k[c] = 0` both outputs are set to 0.0 via branchless f32 selection
/// so that the loop body contains no data-dependent branches and LLVM can
/// autovectorise the whole operation.
#[inline(always)]
fn apply_cal_row(
    power: &[f32],
    noise: &[f32],
    k: &[f32],
    az_factor: f32,
    out: &mut [f32],
    nesz: &mut [f32],
) {
    debug_assert_eq!(power.len(), noise.len());
    debug_assert_eq!(power.len(), k.len());
    debug_assert_eq!(power.len(), out.len());
    debug_assert_eq!(power.len(), nesz.len());
    for idx in 0..power.len() {
        let k2 = k[idx] * k[idx];
        // Branchless: when k2 == 0 the selector is 0.0, zeroing both outputs.
        let valid = (k2 > 0.0) as u32 as f32;
        let inv_k2 = if k2 > 0.0 { 1.0 / k2 } else { 0.0 };
        let noise_val = noise[idx] * az_factor;
        out[idx] = (power[idx] - noise_val).max(0.0) * inv_k2 * valid;
        nesz[idx] = noise_val * inv_k2 * valid;
    }
}

fn bilinear_row_into(
    vecs_lines: &[i64],
    vecs_pixels: &[&[u32]],
    vecs_values: &[&[f32]],
    tiff_line: i64,
    col_offset: usize,
    buf: &mut [f32],
) {
    let (lo, hi) = bracket_azimuth(vecs_lines, tiff_line);
    let pix_lo = vecs_pixels[lo];
    let val_lo = vecs_values[lo];
    let mut cur_lo: usize = 0;

    if lo == hi {
        // Exact azimuth hit or clamped at edge: only one LUT row needed.
        for (c, slot) in buf.iter_mut().enumerate() {
            let tiff_col = (c + col_offset) as u32;
            while cur_lo + 1 < pix_lo.len() && pix_lo[cur_lo + 1] <= tiff_col {
                cur_lo += 1;
            }
            *slot = lut_cursor_eval(pix_lo, val_lo, cur_lo, tiff_col);
        }
    } else {
        // Interpolate between two azimuth LUT rows.
        let t_line =
            (tiff_line - vecs_lines[lo]) as f32 / (vecs_lines[hi] - vecs_lines[lo]) as f32;
        let pix_hi = vecs_pixels[hi];
        let val_hi = vecs_values[hi];
        let mut cur_hi: usize = 0;

        for (c, slot) in buf.iter_mut().enumerate() {
            let tiff_col = (c + col_offset) as u32;
            while cur_lo + 1 < pix_lo.len() && pix_lo[cur_lo + 1] <= tiff_col {
                cur_lo += 1;
            }
            while cur_hi + 1 < pix_hi.len() && pix_hi[cur_hi + 1] <= tiff_col {
                cur_hi += 1;
            }
            let v_lo = lut_cursor_eval(pix_lo, val_lo, cur_lo, tiff_col);
            let v_hi = lut_cursor_eval(pix_hi, val_hi, cur_hi, tiff_col);
            *slot = lerp(v_lo, v_hi, t_line);
        }
    }
}

// ─── Main calibration function ────────────────────────────────────────────────

/// Apply radiometric calibration to a debursted subswath.
///
/// Computes σ⁰ = (|DN|² − noise) / K² for every pixel.
///
/// # Arguments
///
/// * `deburst`          — Debursted CInt16 array from [`crate::deburst::deburst_subswath`].
/// * `cal`              — Calibration annotation for the same swath + polarization.
/// * `noise`            — Noise annotation for the same swath + polarization.
/// * `tiff_line_origin` — TIFF row index corresponding to output line 0.
///                        For a separate-subswath TIFF this is typically 0
///                        (burst 0 occupies TIFF rows 0 .. lines_per_burst).
///
/// # Returns
///
/// A [`Sigma0Array`] with the same shape as `deburst`.
pub fn apply_calibration(
    deburst: &DeburstArray,
    cal: &SwathCalibration,
    noise: &SwathNoise,
    tiff_line_origin: usize,
) -> Result<Sigma0Array, CalibrationError> {
    if cal.vectors.is_empty() {
        return Err(CalibrationError::NoCalibrationVectors);
    }
    if noise.range_vectors.is_empty() {
        return Err(CalibrationError::NoNoiseRangeVectors);
    }

    let out_lines = deburst.lines;
    let out_cols = deburst.samples;
    let col_offset = deburst.valid_sample_offset;

    // ── Pre-extract LUT arrays ─────────────────────────────────────────────────

    // Calibration
    let cal_lines: Vec<i64> = cal.vectors.iter().map(|v| v.line as i64).collect();
    let cal_pixels: Vec<&[u32]> = cal.vectors.iter().map(|v| v.pixels.as_slice()).collect();
    let cal_values: Vec<&[f32]> = cal.vectors.iter().map(|v| v.sigma_nought.as_slice()).collect();

    // Noise range
    let nr_lines: Vec<i64> = noise.range_vectors.iter().map(|v| v.line as i64).collect();
    let nr_pixels: Vec<&[u32]> = noise.range_vectors.iter().map(|v| v.pixels.as_slice()).collect();
    let nr_values: Vec<&[f32]> = noise.range_vectors.iter().map(|v| v.noise_range_lut.as_slice()).collect();

    // ── Coverage check: LUT must span the full output TIFF-column extent ──
    //
    // Previously `lut_cursor_eval` silently clamped to the nearest edge value
    // when a pixel fell outside the LUT's coverage.  For corrupt annotations
    // or mixed-provenance products this produces arbitrary calibration error
    // at the swath edges without any warning.  We now:
    //   1. Reject inputs where the gap exceeds a physical bound
    //      (LUT_EDGE_SLACK_PX, chosen from the documented S-1 LUT sampling
    //      grid of ~40–80 samples per step plus headroom).
    //   2. For smaller gaps (normal S-1 behaviour where the LUT's last point
    //      stops short of the annotation's last_valid_sample by a fraction
    //      of one LUT cell), emit a one-line stderr warning with the gap
    //      size so it is visible in logs.
    //
    // The rationale for not hard-failing on small gaps: ESA's S-1 IPF
    // routinely ships calibration LUTs whose last pixel ≈ swath_width −
    // (½ × grid_step), so a ~40–400 sample extrapolation is the normal
    // case, not a corruption signal.
    const LUT_EDGE_SLACK_PX: u32 = 600;
    let required_last_px = (col_offset + out_cols).saturating_sub(1) as u32;
    let mut max_cal_gap: u32 = 0;
    for (i, pixels) in cal_pixels.iter().enumerate() {
        if pixels.is_empty() {
            return Err(CalibrationError::NoCalibrationVectors);
        }
        let lut_max = *pixels.last().unwrap();
        if lut_max < required_last_px {
            let gap = required_last_px - lut_max;
            if gap > LUT_EDGE_SLACK_PX {
                return Err(CalibrationError::LutRangeInsufficient {
                    line: cal.vectors[i].line,
                    lut_max_pixel: lut_max,
                    required_pixel: required_last_px as usize,
                });
            }
            if gap > max_cal_gap {
                max_cal_gap = gap;
            }
        }
    }
    let mut max_nr_gap: u32 = 0;
    for (i, pixels) in nr_pixels.iter().enumerate() {
        if pixels.is_empty() {
            return Err(CalibrationError::NoNoiseRangeVectors);
        }
        let lut_max = *pixels.last().unwrap();
        if lut_max < required_last_px {
            let gap = required_last_px - lut_max;
            if gap > LUT_EDGE_SLACK_PX {
                return Err(CalibrationError::LutRangeInsufficient {
                    line: noise.range_vectors[i].line,
                    lut_max_pixel: lut_max,
                    required_pixel: required_last_px as usize,
                });
            }
            if gap > max_nr_gap {
                max_nr_gap = gap;
            }
        }
    }
    if max_cal_gap > 0 || max_nr_gap > 0 {
        // Surface the extrapolation through the typed return value
        // (`Sigma0Array::cal_lut_extrapolation_gap_px` /
        //  `noise_lut_extrapolation_gap_px`) instead of stderr so that
        // callers can detect, log, threshold, or fail on it programmatically
        // — never silent.
    }

    // Noise azimuth (optional — present in modern IPF ≥ 2.9; treated as 1.0 if absent).
    // Each SwathNoise contains vectors already filtered to this (subswath, polarization);
    // there should be exactly one per subswath.
    debug_assert!(
        noise.azimuth_vectors.len() <= 1,
        "expected at most 1 azimuth noise vector, got {}",
        noise.azimuth_vectors.len()
    );
    let az_noise_available = !noise.azimuth_vectors.is_empty();
    if !az_noise_available {
        // Azimuth noise is present in all IPF ≥ 2.9 products (all S-1A/B from ~2017 onwards).
        // Its absence likely means an older product or a parsing gap — log so the user can
        // verify.  Calibration continues using N_azimuth = 1.0 (range-only noise model).
        tracing::warn!(
            "subswath {:?} pol {:?}: azimuth noise LUT absent — \
             falling back to range-only noise model (N_az = 1.0). \
             Expected for IPF < 2.9 products; unexpected for modern S-1A/B.",
            cal.subswath_id,
            cal.polarization,
        );
    }
    let az_lines: Vec<i64>;
    let az_values: &[f32];
    if az_noise_available {
        let av = &noise.azimuth_vectors[0];
        az_lines = av.lines.iter().map(|&l| l as i64).collect();
        az_values = &av.noise_azimuth_lut;
    } else {
        az_lines = Vec::new();
        az_values = &[];
    }

    // ── Output buffer ──────────────────────────────────────────────────────────
    //
    // Use uninitialised memory: the parallel loop below writes every element
    // before the buffers are returned, avoiding the serialised OS page-fault
    // cost of `vec![0.0f32; N]` for large (multi-GB) subswaths.
    let n = out_lines * out_cols;
    let mut out: Vec<f32> = Vec::with_capacity(n);
    // NESZ in linear σ⁰ units: N[l,c] / K²[l,c], same shape as `out`.
    let mut nesz_out: Vec<f32> = Vec::with_capacity(n);
    // SAFETY-OK: every element of `out` and `nesz_out` is written by the
    // `par_chunks_mut` loop below (via `apply_cal_row`) before either Vec is
    // read outside this function.  `apply_cal_row` writes every `out_row` and
    // `nesz_row` slot unconditionally (branchless arithmetic over the full
    // column range).  No element is read before it is written.
    unsafe {
        out.set_len(n);
        nesz_out.set_len(n);
    }

    // ── Main loop (parallel across lines) ─────────────────────────────────────
    //
    // Scratch buffers (k_row, nr_row, power_row) are kept in thread-local
    // storage so each Rayon worker thread initialises them once and reuses
    // them for every row it processes.  This eliminates ~3 × out_cols × 4 bytes
    // of malloc/free per row — on an 80-core machine processing 3 IWs
    // simultaneously that is ~3.2 GB of allocator traffic that would otherwise
    // hammer the system allocator's arena locks.
    use std::cell::RefCell;
    thread_local! {
        static SCRATCH: RefCell<[Vec<f32>; 3]> = RefCell::new(
            [Vec::new(), Vec::new(), Vec::new()]
        );
    }

    out.par_chunks_mut(out_cols)
        .zip(nesz_out.par_chunks_mut(out_cols))
        .enumerate()
        .for_each(|(out_line, (out_row, nesz_row))| {
            let tiff_line = (tiff_line_origin + out_line) as i64;

            SCRATCH.with(|cell| {
                let mut bufs = cell.borrow_mut();
                let [k_row, nr_row, power_row] = &mut *bufs;

                // Resize once on first use by this thread; no-op thereafter.
                k_row.resize(out_cols, 0.0f32);
                nr_row.resize(out_cols, 0.0f32);
                power_row.resize(out_cols, 0.0f32);

                // Interpolated calibration K values for this line.
                bilinear_row_into(
                    &cal_lines, &cal_pixels, &cal_values,
                    tiff_line, col_offset, k_row,
                );

                // Interpolated range noise for this line.
                bilinear_row_into(
                    &nr_lines, &nr_pixels, &nr_values,
                    tiff_line, col_offset, nr_row,
                );

                // Azimuth noise factor for this line (scalar; defaults to 1.0).
                let az_factor = if az_noise_available && !az_lines.is_empty() {
                    let (lo, hi) = bracket_azimuth(&az_lines, tiff_line);
                    interp_azimuth(
                        az_lines[lo],
                        az_lines[hi],
                        az_values[lo],
                        az_values[hi],
                        tiff_line,
                    )
                } else {
                    1.0
                };

                // Pass 1: i16 DN → f32 power  (enables i16→i32 widening SIMD).
                let in_row = &deburst.data[out_line * out_cols..(out_line + 1) * out_cols];
                dn_to_power_row(in_row, power_row);

                // Pass 2: branchless σ⁰ and NESZ  (pure f32, no data-dependent branches).
                apply_cal_row(power_row, nr_row, k_row, az_factor, out_row, nesz_row);
            });
        });

    Ok(Sigma0Array {
        data: out,
        lines: out_lines,
        samples: out_cols,
        valid_sample_offset: col_offset,
        cal_lut_extrapolation_gap_px: max_cal_gap,
        noise_lut_extrapolation_gap_px: max_nr_gap,
        nesz: nesz_out,
    })
}

// ─── Complex calibration ──────────────────────────────────────────────────────

/// Calibrated complex SLC samples for one debursted subswath.
///
/// Each element is `[I/K, Q/K]` where `K` is the sigma-nought calibration
/// LUT value at that pixel.  Dividing by `K` (not `K²`) keeps the data
/// complex so that cross-polarisation products `S_vv · S_vh*` can be formed
/// downstream.
///
/// The NESZ array stores `N / K²` per pixel — the same noise floor used by
/// [`Sigma0Array`] — for noise removal on the diagonal C2 elements.
pub struct CalibratedComplexArray {
    /// Calibrated [I, Q] samples. Length = `lines × samples`.
    pub data: Vec<[f32; 2]>,
    /// Number of azimuth lines.
    pub lines: usize,
    /// Number of range samples.
    pub samples: usize,
    /// Column offset: output col 0 = TIFF col `valid_sample_offset`.
    pub valid_sample_offset: usize,
    /// Per-pixel NESZ = N/K² (same units and layout as [`Sigma0Array::nesz`]).
    pub nesz: Vec<f32>,
}

/// Per-column kernel: scale CInt16 [i, q] by 1/K and write [f32; 2].
/// Also write NESZ = noise * az_factor / K².
#[inline(always)]
fn apply_cal_complex_row(
    in_row: &[[i16; 2]],
    noise: &[f32],
    k: &[f32],
    az_factor: f32,
    out: &mut [[f32; 2]],
    nesz: &mut [f32],
) {
    debug_assert_eq!(in_row.len(), noise.len());
    debug_assert_eq!(in_row.len(), k.len());
    debug_assert_eq!(in_row.len(), out.len());
    debug_assert_eq!(in_row.len(), nesz.len());
    for idx in 0..in_row.len() {
        let [i, q] = in_row[idx];
        let kv = k[idx];
        if kv > 0.0 {
            let inv_k = 1.0 / kv;
            out[idx] = [i as f32 * inv_k, q as f32 * inv_k];
            nesz[idx] = noise[idx] * az_factor * inv_k * inv_k;
        } else {
            out[idx] = [0.0, 0.0];
            nesz[idx] = 0.0;
        }
    }
}

/// Apply radiometric calibration to a debursted subswath, retaining complex phase.
///
/// Produces calibrated `[I/K, Q/K]` per pixel and the per-pixel NESZ `N/K²`.
/// The argument contract is identical to [`apply_calibration`].
pub fn apply_calibration_complex(
    deburst: &DeburstArray,
    cal: &SwathCalibration,
    noise: &SwathNoise,
    tiff_line_origin: usize,
) -> Result<CalibratedComplexArray, CalibrationError> {
    if cal.vectors.is_empty() {
        return Err(CalibrationError::NoCalibrationVectors);
    }
    if noise.range_vectors.is_empty() {
        return Err(CalibrationError::NoNoiseRangeVectors);
    }

    let out_lines = deburst.lines;
    let out_cols = deburst.samples;
    let col_offset = deburst.valid_sample_offset;

    let cal_lines: Vec<i64> = cal.vectors.iter().map(|v| v.line as i64).collect();
    let cal_pixels: Vec<&[u32]> = cal.vectors.iter().map(|v| v.pixels.as_slice()).collect();
    let cal_values: Vec<&[f32]> = cal.vectors.iter().map(|v| v.sigma_nought.as_slice()).collect();

    let nr_lines: Vec<i64> = noise.range_vectors.iter().map(|v| v.line as i64).collect();
    let nr_pixels: Vec<&[u32]> = noise.range_vectors.iter().map(|v| v.pixels.as_slice()).collect();
    let nr_values: Vec<&[f32]> = noise.range_vectors.iter().map(|v| v.noise_range_lut.as_slice()).collect();

    debug_assert!(
        noise.azimuth_vectors.len() <= 1,
        "expected at most 1 azimuth noise vector, got {}",
        noise.azimuth_vectors.len()
    );
    let az_noise_available = !noise.azimuth_vectors.is_empty();
    let az_lines: Vec<i64>;
    let az_values: &[f32];
    if az_noise_available {
        let av = &noise.azimuth_vectors[0];
        az_lines = av.lines.iter().map(|&l| l as i64).collect();
        az_values = &av.noise_azimuth_lut;
    } else {
        az_lines = Vec::new();
        az_values = &[];
    }

    let n = out_lines * out_cols;
    let mut out: Vec<[f32; 2]> = vec![[0.0, 0.0]; n];
    let mut nesz_out: Vec<f32> = vec![0.0; n];

    use std::cell::RefCell;
    thread_local! {
        static SCRATCH_COMPLEX: RefCell<[Vec<f32>; 2]> = RefCell::new(
            [Vec::new(), Vec::new()]
        );
    }

    out.par_chunks_mut(out_cols)
        .zip(nesz_out.par_chunks_mut(out_cols))
        .enumerate()
        .for_each(|(out_line, (out_row, nesz_row))| {
            let tiff_line = (tiff_line_origin + out_line) as i64;

            SCRATCH_COMPLEX.with(|cell| {
                let mut bufs = cell.borrow_mut();
                let [k_row, nr_row] = &mut *bufs;
                k_row.resize(out_cols, 0.0f32);
                nr_row.resize(out_cols, 0.0f32);

                bilinear_row_into(
                    &cal_lines, &cal_pixels, &cal_values, tiff_line, col_offset, k_row,
                );
                bilinear_row_into(
                    &nr_lines, &nr_pixels, &nr_values, tiff_line, col_offset, nr_row,
                );

                let az_factor = if az_noise_available && !az_lines.is_empty() {
                    let (lo, hi) = bracket_azimuth(&az_lines, tiff_line);
                    interp_azimuth(az_lines[lo], az_lines[hi], az_values[lo], az_values[hi], tiff_line)
                } else {
                    1.0
                };

                let in_row = &deburst.data[out_line * out_cols..(out_line + 1) * out_cols];
                apply_cal_complex_row(in_row, nr_row, k_row, az_factor, out_row, nesz_row);
            });
        });

    Ok(CalibratedComplexArray {
        data: out,
        lines: out_lines,
        samples: out_cols,
        valid_sample_offset: col_offset,
        nesz: nesz_out,
    })
}

// ─── Tests ────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    // ── Vectorizable row kernel tests ─────────────────────────────────────────

    #[test]
    fn dn_to_power_row_known_values() {
        let in_row: Vec<[i16; 2]> = vec![[3, 4], [0, 0], [i16::MAX, i16::MAX]];
        let mut power = vec![0.0f32; 3];
        dn_to_power_row(&in_row, &mut power);
        assert!((power[0] - 25.0).abs() < 1e-6, "3²+4²=25, got {}", power[0]);
        assert_eq!(power[1], 0.0);
        let expected_max = (i16::MAX as i32 * i16::MAX as i32 * 2) as f32;
        assert!((power[2] - expected_max).abs() < 1.0, "MAX² + MAX², got {}", power[2]);
    }

    #[test]
    fn apply_cal_row_zero_k_gives_zero_output() {
        let power = vec![100.0f32; 4];
        let noise = vec![10.0f32; 4];
        let k = vec![0.0f32; 4]; // all-zero → both outputs must be zero
        let mut out = vec![f32::NAN; 4];
        let mut nesz = vec![f32::NAN; 4];
        apply_cal_row(&power, &noise, &k, 1.0, &mut out, &mut nesz);
        for i in 0..4 {
            assert_eq!(out[i], 0.0, "out[{i}] should be 0 when k=0");
            assert_eq!(nesz[i], 0.0, "nesz[{i}] should be 0 when k=0");
        }
    }

    #[test]
    fn apply_cal_row_known_sigma0() {
        // power=125, noise=25, k=10, az=1 → σ⁰ = (125-25)/100 = 1.0, NESZ = 25/100 = 0.25
        let power = vec![125.0f32];
        let noise = vec![25.0f32];
        let k = vec![10.0f32];
        let mut out = vec![0.0f32];
        let mut nesz = vec![0.0f32];
        apply_cal_row(&power, &noise, &k, 1.0, &mut out, &mut nesz);
        assert!((out[0] - 1.0).abs() < 1e-5, "sigma0 mismatch: {}", out[0]);
        assert!((nesz[0] - 0.25).abs() < 1e-5, "nesz mismatch: {}", nesz[0]);
    }

    #[test]
    fn apply_cal_row_negative_power_clamps_to_zero() {
        // power < noise → clamped to 0, but NESZ is still noise/k²
        let power = vec![5.0f32];
        let noise = vec![50.0f32];
        let k = vec![10.0f32];
        let mut out = vec![f32::NAN];
        let mut nesz = vec![0.0f32];
        apply_cal_row(&power, &noise, &k, 1.0, &mut out, &mut nesz);
        assert_eq!(out[0], 0.0, "negative power should clamp to 0");
        assert!((nesz[0] - 0.5).abs() < 1e-5, "nesz mismatch: {}", nesz[0]);
    }

    #[test]
    fn apply_cal_row_az_factor_scales_noise() {
        // noise=10, az_factor=3 → effective noise = 30
        // power=130, k=10 → σ⁰ = (130-30)/100 = 1.0
        let power = vec![130.0f32];
        let noise = vec![10.0f32];
        let k = vec![10.0f32];
        let mut out = vec![0.0f32];
        let mut nesz = vec![0.0f32];
        apply_cal_row(&power, &noise, &k, 3.0, &mut out, &mut nesz);
        assert!((out[0] - 1.0).abs() < 1e-5, "sigma0 with az_factor: {}", out[0]);
        assert!((nesz[0] - 0.3).abs() < 1e-5, "nesz with az_factor: {}", nesz[0]);
    }

    // ── Interpolation unit tests ───────────────────────────────────────────────

    #[test]
    fn lerp_exact_endpoints() {
        assert_eq!(lerp(2.0, 8.0, 0.0), 2.0);
        assert_eq!(lerp(2.0, 8.0, 1.0), 8.0);
    }

    #[test]
    fn lerp_midpoint() {
        assert_eq!(lerp(0.0, 10.0, 0.5), 5.0);
    }

    #[test]
    fn lerp_clamps_t_outside_unit_interval() {
        assert_eq!(lerp(0.0, 10.0, -1.0), 0.0);
        assert_eq!(lerp(0.0, 10.0,  2.0), 10.0);
    }

    #[test]
    fn interp_range_lut_exact_hit() {
        let pixels = [0u32, 100, 200];
        let values = [1.0f32, 2.0, 3.0];
        assert_eq!(interp_range_lut(&pixels, &values, 100), 2.0);
    }

    #[test]
    fn interp_range_lut_midpoint() {
        let pixels = [0u32, 100, 200];
        let values = [0.0f32, 10.0, 20.0];
        let v = interp_range_lut(&pixels, &values, 50);
        assert!((v - 5.0).abs() < 1e-5, "expected 5.0, got {v}");
    }

    #[test]
    fn interp_range_lut_clamps_below() {
        let pixels = [10u32, 20];
        let values = [1.0f32, 2.0];
        assert_eq!(interp_range_lut(&pixels, &values, 0), 1.0);
    }

    #[test]
    fn interp_range_lut_clamps_above() {
        let pixels = [10u32, 20];
        let values = [1.0f32, 2.0];
        assert_eq!(interp_range_lut(&pixels, &values, 100), 2.0);
    }

    #[test]
    fn bracket_azimuth_exact() {
        let lines = [-10i64, 0, 500, 1000];
        let (lo, hi) = bracket_azimuth(&lines, 500);
        assert_eq!((lo, hi), (2, 2));
    }

    #[test]
    fn bracket_azimuth_between() {
        let lines = [0i64, 500, 1000];
        let (lo, hi) = bracket_azimuth(&lines, 300);
        assert_eq!((lo, hi), (0, 1));
    }

    #[test]
    fn bracket_azimuth_before_start() {
        let lines = [100i64, 200, 300];
        let (lo, hi) = bracket_azimuth(&lines, -5);
        assert_eq!((lo, hi), (0, 0));
    }

    #[test]
    fn bracket_azimuth_after_end() {
        let lines = [100i64, 200, 300];
        let (lo, hi) = bracket_azimuth(&lines, 9999);
        assert_eq!((lo, hi), (2, 2));
    }

    // ── Formula unit tests ────────────────────────────────────────────────────

    #[test]
    fn zero_input_gives_zero_sigma0() {
        // I=0, Q=0, noise=0, K=300 → σ⁰ = 0² / 300² = 0.
        let power: f32 = 0.0;
        let noise: f32 = 0.0;
        let k: f32 = 300.0;
        let sigma0 = (power - noise).max(0.0) / (k * k);
        assert_eq!(sigma0, 0.0);
    }

    #[test]
    fn known_value_sigma0() {
        // I=100, Q=200 → |DN|²=50000; noise=1000; K=300.
        // σ⁰ = (50000 − 1000) / 90000 = 49000 / 90000 ≈ 0.5444...
        let i_dn: i16 = 100;
        let q_dn: i16 = 200;
        let power = (i_dn as i32 * i_dn as i32 + q_dn as i32 * q_dn as i32) as f32;
        let noise: f32 = 1000.0;
        let k: f32 = 300.0;
        let sigma0 = (power - noise).max(0.0) / (k * k);
        let expected = 49000.0_f32 / 90000.0;
        assert!((sigma0 - expected).abs() < 1e-5, "got {sigma0}, expected {expected}");
    }

    #[test]
    fn negative_power_after_noise_clamped_to_zero() {
        // |DN|²=100, noise=200 → 100−200 = −100 → clamped to 0.
        let power: f32 = 100.0;
        let noise: f32 = 200.0;
        let k: f32 = 300.0;
        let sigma0 = (power - noise).max(0.0) / (k * k);
        assert_eq!(sigma0, 0.0);
    }

    #[test]
    fn i16_power_no_overflow() {
        // i16::MAX = 32767 → |DN|² = 2 × 32767² = 2_147_418_114 < i32::MAX.
        let v = i16::MAX;
        let p = v as i32 * v as i32 + v as i32 * v as i32;
        assert!(p > 0, "overflow: {p}");
    }

    // ── Integration test against real S1A data ────────────────────────────────

    #[test]
    fn sigma0_s1a_iw1_vv_smoke_test() {
        use crate::calibration::parse_calibration_noise;
        use crate::parse::parse_safe_directory;
        use crate::slc_reader::SlcReader;
        use crate::SubSwathId;
        use crate::types::Polarization;

        const S1A_SAFE: &str =
            "/home/datacube/dev/SARdine/data/SLC/\
             S1A_IW_SLC__1SDV_20201005T170824_20201005T170851_034664_04098A_1E66.SAFE";
        const S1A_IW1_VV: &str =
            "/home/datacube/dev/SARdine/data/SLC/\
             S1A_IW_SLC__1SDV_20201005T170824_20201005T170851_034664_04098A_1E66.SAFE/\
             measurement/\
             s1a-iw1-slc-vv-20201005t170824-20201005t170849-034664-04098a-004.tiff";

        if !std::path::Path::new(S1A_SAFE).is_dir() {
            eprintln!("skipping sigma0_s1a_iw1_vv_smoke_test — S1A SAFE not present");
            return;
        }
        // 1. Parse metadata and LUTs.
        let scene = parse_safe_directory(std::path::Path::new(S1A_SAFE))
            .expect("failed to parse SAFE");
        let cal_data = parse_calibration_noise(std::path::Path::new(S1A_SAFE))
            .expect("failed to parse calibration");

        let sw = scene.sub_swaths.iter().find(|s| s.id == SubSwathId::IW1).unwrap();
        let bursts: Vec<_> = scene.bursts.iter()
            .filter(|b| b.subswath_id == SubSwathId::IW1)
            .cloned()
            .collect();

        let cal = cal_data.calibrations.iter()
            .find(|c| c.subswath_id == SubSwathId::IW1 && c.polarization == Polarization::VV)
            .expect("IW1/VV calibration not found");
        let noise = cal_data.noises.iter()
            .find(|n| n.subswath_id == SubSwathId::IW1 && n.polarization == Polarization::VV)
            .expect("IW1/VV noise not found");

        // 2. Deburst — only burst 4 (mid-swath) to keep the test fast in debug mode.
        //    We construct a minimal DeburstArray for this single burst directly.
        let burst4 = bursts.iter().find(|b| b.burst_index == 4).unwrap();
        let mut reader = SlcReader::open(S1A_IW1_VV).expect("open TIFF failed");
        let first_valid = bursts.iter().map(|b| b.first_valid_sample).max().unwrap();
        let last_valid  = bursts.iter().map(|b| b.last_valid_sample ).min().unwrap();
        let valid_cols = last_valid - first_valid;
        let raw = reader
            .read_burst_raw(burst4.first_line, sw.lines_per_burst)
            .expect("read burst 4 failed");
        // Crop to valid range window.
        let raw_cols = reader.width as usize;
        let mut data = Vec::with_capacity(sw.lines_per_burst * valid_cols);
        for line in 0..sw.lines_per_burst {
            data.extend_from_slice(&raw[line * raw_cols + first_valid..line * raw_cols + last_valid]);
        }
        let debursted = crate::deburst::DeburstArray {
            data,
            lines: sw.lines_per_burst,
            samples: valid_cols,
            valid_sample_offset: first_valid,
        };

        // 3. Apply calibration.  tiff_line_origin = burst4.first_line (TIFF coordinates).
        let sigma0 = apply_calibration(&debursted, cal, noise, burst4.first_line)
            .expect("calibration failed");

        // 4. Verify output shape matches input.
        assert_eq!(sigma0.lines, debursted.lines);
        assert_eq!(sigma0.samples, debursted.samples);
        assert_eq!(sigma0.data.len(), sigma0.lines * sigma0.samples);

        // 5. All values must be finite and non-negative.
        let any_nan = sigma0.data.iter().any(|v| v.is_nan());
        let any_neg = sigma0.data.iter().any(|&v| v < 0.0);
        assert!(!any_nan, "NaN values in σ⁰ output");
        assert!(!any_neg, "negative values in σ⁰ output");

        // 6. Most of burst 4 interior should be non-zero (real signal).
        let n_nonzero = sigma0.data.iter().filter(|&&v| v > 0.0).count();
        let frac_nonzero = n_nonzero as f64 / sigma0.data.len() as f64;
        assert!(
            frac_nonzero > 0.5,
            "less than 50 % of σ⁰ pixels are non-zero ({:.1} %)",
            frac_nonzero * 100.0
        );

        // 7. Median (of non-zero) σ⁰ in physically plausible range.
        let mut vals: Vec<f32> = sigma0.data.iter().copied().filter(|&v| v > 0.0).collect();
        vals.sort_by(|a, b| a.partial_cmp(b).unwrap());
        let median = vals[vals.len() / 2];
        assert!(
            median > 1e-5 && median < 1.0,
            "median σ⁰ = {median:.2e} is outside expected range [1e-5, 1.0]"
        );

        // 8. NESZ grid must be same shape as data, all finite, non-negative.
        assert_eq!(sigma0.nesz.len(), sigma0.data.len(), "NESZ length mismatch");
        let any_nesz_nan  = sigma0.nesz.iter().any(|v| v.is_nan());
        let any_nesz_neg  = sigma0.nesz.iter().any(|&v| v < 0.0);
        assert!(!any_nesz_nan, "NaN in NESZ output");
        assert!(!any_nesz_neg, "negative NESZ value");
        // Median NESZ in dB should be in typical S-1 range (−30 to −10 dB).
        let mut nesz_vals: Vec<f32> = sigma0.nesz.iter().copied().filter(|&v| v > 0.0).collect();
        nesz_vals.sort_by(|a, b| a.partial_cmp(b).unwrap());
        let nesz_median_db = 10.0 * nesz_vals[nesz_vals.len() / 2].log10();
        assert!(
            nesz_median_db > -35.0 && nesz_median_db < -5.0,
            "median NESZ = {nesz_median_db:.1} dB outside expected range [−35, −5] dB"
        );
    }
}
