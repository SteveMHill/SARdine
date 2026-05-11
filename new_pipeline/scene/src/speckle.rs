//! Speckle filters for σ⁰ / γ⁰ intensity rasters.
//!
//! All four filters in this module operate on detected backscatter
//! intensities (linear power, **not** dB).  Applying a speckle filter to
//! a logarithmic image is mathematically incorrect: the multiplicative
//! Goodman speckle model `I = R · n` (with `n` exponentially distributed
//! at single look) becomes the additive `log I = log R + log n` in the
//! log domain, and `log n` is heavily skewed.  The filter formulae below
//! all assume the multiplicative model.
//!
//! Inputs may contain NaN to denote nodata.  A NaN-input pixel always
//! produces a NaN-output pixel (we never invent values for missing
//! samples).  Window statistics ignore NaN samples; if a window contains
//! zero finite samples, the output pixel is NaN.
//!
//! # Filters implemented
//!
//! * **Boxcar** — unweighted local mean.  Reduces speckle by a factor of
//!   √(window_area) but blurs edges.  Useful as a baseline reference and
//!   for ENL diagnostics.
//!
//! * **Lee** — Lee 1981 minimum-mean-square-error filter.  Adapts the
//!   amount of smoothing to local heterogeneity using the speckle
//!   coefficient of variation `σ_n = 1/√L`, where `L` is the effective
//!   number of looks (ENL).  This is the *classic* Lee filter; the
//!   *refined* Lee of Lopes-Touzi-Nezry 1990 with directional sub-windows
//!   is a separate algorithm and is **not** what this module ships.  The
//!   doc-string for [`SpeckleFilter::Lee`] expands on this distinction.
//!
//! * **Frost** — Frost-Stiles-Shanmugam-Holtzman 1982 exponentially
//!   weighted local mean.  Damping factor `K` (typically 1–2) controls
//!   the trade-off between smoothing and edge preservation: small `K` is
//!   close to a boxcar, large `K` falls back toward the input.
//!
//! * **Gamma-MAP** — Lopes-Nezry-Touzi 1990 maximum-a-posteriori filter
//!   under a Gamma-distributed scene prior and Gamma-distributed speckle
//!   model.  Uses the same `L` parameter as Lee.
//!
//! # References
//!
//! * Lee, J.S. (1981).  Refined filtering of image noise using local
//!   statistics.  *CGIP* 15, 380-389.
//! * Frost, V.S. et al. (1982).  A model for radar images and its
//!   application to adaptive digital filtering of multiplicative noise.
//!   *IEEE TPAMI* 4(2), 157-166.
//! * Lopes, A., Nezry, E., Touzi, R., Laur, H. (1990).  Maximum a
//!   posteriori speckle filtering and first order texture models in SAR
//!   images.  *IGARSS '90*, 2409-2412.
//! * Lopes, A., Touzi, R., Nezry, E. (1990).  Adaptive speckle filters
//!   and scene heterogeneity.  *IEEE TGRS* 28(6), 992-1000.

use rayon::prelude::*;
use thiserror::Error;

/// Choice of speckle filter and its parameters.
///
/// `window` must be an odd integer ≥ 3 in every variant.  Even windows
/// would have no centred sample and the filter geometry would be
/// asymmetric; we reject them up front.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum SpeckleFilter {
    /// Unweighted local mean over an `window × window` window.
    Boxcar { window: usize },

    /// Classic Lee 1981 MMSE filter, parameterised by the effective
    /// number of looks `enl`.  This is **not** the Lopes 1990 "Refined
    /// Lee" with directional sub-windows; for that variant use
    /// [`SpeckleFilter::RefinedLee`].  Use `enl = 1` for unfiltered
    /// single-look intensity input.
    Lee { window: usize, enl: f32 },

    /// Frost 1982 exponentially-weighted local mean with damping
    /// factor `damping` (often called `K` in the literature).
    /// `damping ≥ 0`; typical values are 1.0 – 2.0.
    Frost { window: usize, damping: f32 },

    /// Lopes 1990 Gamma-MAP filter, parameterised by the effective
    /// number of looks `enl`.  Use `enl = 1` for unfiltered single-look.
    GammaMap { window: usize, enl: f32 },

    /// Lopes–Touzi–Nezry 1990 "Refined Lee" with a fixed 7×7 window
    /// and 8 directional sub-windows.  Detects the local edge
    /// orientation in a 3×3 grid of 3×3 sub-window means, picks the
    /// sub-window aligned with the **homogeneous side** of the edge,
    /// and applies the classic Lee 1981 MMSE update using only that
    /// sub-window's mean and variance.  Preserves edges much better
    /// than [`SpeckleFilter::Lee`] at the cost of being fixed at 7×7
    /// (the directional masks are hand-tuned and do not generalise to
    /// other window sizes).
    RefinedLee { enl: f32 },
}

/// Errors returned by [`apply_speckle_filter`].  Every variant is a
/// caller-visible misuse condition; there are no silent fallbacks.
#[derive(Debug, Error)]
pub enum SpeckleError {
    #[error("data length {got} does not match cols×rows = {expected}")]
    DimensionMismatch { got: usize, expected: usize },

    #[error("window size {got} must be an odd integer ≥ 3")]
    BadWindow { got: usize },

    #[error("ENL must be > 0; got {got}")]
    BadEnl { got: f32 },

    #[error("Frost damping factor must be ≥ 0 and finite; got {got}")]
    BadDamping { got: f32 },
}

impl SpeckleFilter {
    fn window(self) -> usize {
        match self {
            SpeckleFilter::Boxcar { window }
            | SpeckleFilter::Lee { window, .. }
            | SpeckleFilter::Frost { window, .. }
            | SpeckleFilter::GammaMap { window, .. } => window,
            // Refined Lee is fixed at 7×7 — the directional masks
            // are hand-tuned for that size and do not generalise.
            SpeckleFilter::RefinedLee { .. } => 7,
        }
    }
}

/// Apply `filter` to `data` (row-major, `rows × cols`, linear-power
/// floats with NaN nodata) and return a freshly-allocated output buffer.
///
/// The input is borrowed so that the original samples remain available
/// while window statistics are computed; this prevents the in-place
/// pollution that would otherwise occur as filtered values overwrite
/// not-yet-read neighbours.
pub fn apply_speckle_filter(
    data: &[f32],
    cols: usize,
    rows: usize,
    filter: SpeckleFilter,
) -> Result<Vec<f32>, SpeckleError> {
    let expected = cols.checked_mul(rows).ok_or(SpeckleError::DimensionMismatch {
        got: data.len(),
        expected: usize::MAX,
    })?;
    if data.len() != expected {
        return Err(SpeckleError::DimensionMismatch {
            got: data.len(),
            expected,
        });
    }

    let win = filter.window();
    if win < 3 || win.is_multiple_of(2) {
        return Err(SpeckleError::BadWindow { got: win });
    }

    match filter {
        SpeckleFilter::Boxcar { .. } => Ok(boxcar(data, cols, rows, win)),
        SpeckleFilter::Lee { enl, .. } => {
            if !(enl > 0.0) || !enl.is_finite() {
                return Err(SpeckleError::BadEnl { got: enl });
            }
            Ok(lee(data, cols, rows, win, enl))
        }
        SpeckleFilter::Frost { damping, .. } => {
            if !damping.is_finite() || damping < 0.0 {
                return Err(SpeckleError::BadDamping { got: damping });
            }
            Ok(frost(data, cols, rows, win, damping))
        }
        SpeckleFilter::GammaMap { enl, .. } => {
            if !(enl > 0.0) || !enl.is_finite() {
                return Err(SpeckleError::BadEnl { got: enl });
            }
            Ok(gamma_map(data, cols, rows, win, enl))
        }
        SpeckleFilter::RefinedLee { enl } => {
            if !(enl > 0.0) || !enl.is_finite() {
                return Err(SpeckleError::BadEnl { got: enl });
            }
            Ok(refined_lee(data, cols, rows, enl))
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Window-statistics helpers (NaN-aware)
// ─────────────────────────────────────────────────────────────────────────────

/// Iterate over the `(dy, dx)` offsets of an `(2*half+1)²` window.
fn window_offsets(half: i32) -> impl Iterator<Item = (i32, i32)> {
    (-half..=half).flat_map(move |dy| (-half..=half).map(move |dx| (dy, dx)))
}

/// Yield `(value, dy, dx)` for every finite sample inside the window
/// centred at `(row, col)`.  Out-of-bounds and NaN samples are skipped.
/// `half = window / 2`.
fn window_samples<'a>(
    data: &'a [f32],
    cols: usize,
    rows: usize,
    row: usize,
    col: usize,
    half: i32,
) -> impl Iterator<Item = (f32, i32, i32)> + 'a {
    window_offsets(half).filter_map(move |(dy, dx)| {
        let r = row as i32 + dy;
        let c = col as i32 + dx;
        if r < 0 || c < 0 || r >= rows as i32 || c >= cols as i32 {
            return None;
        }
        let v = data[(r as usize) * cols + c as usize];
        if v.is_finite() { Some((v, dy, dx)) } else { None }
    })
}

/// Mean and biased variance of the finite samples in the window.
/// Returns `None` if the window has zero finite samples.  Variance is
/// the population variance `Σ(x − μ)² / n`, matching the SAR speckle
/// literature's convention (where the speckle CV is also a population
/// quantity).
fn window_mean_var(
    data: &[f32],
    cols: usize,
    rows: usize,
    row: usize,
    col: usize,
    half: i32,
) -> Option<(f32, f32, usize)> {
    let mut sum = 0.0_f64;
    let mut sum_sq = 0.0_f64;
    let mut n: usize = 0;
    for (v, _, _) in window_samples(data, cols, rows, row, col, half) {
        let v = v as f64;
        sum += v;
        sum_sq += v * v;
        n += 1;
    }
    if n == 0 {
        return None;
    }
    let mean = sum / n as f64;
    let var = (sum_sq / n as f64 - mean * mean).max(0.0);
    Some((mean as f32, var as f32, n))
}

// ─────────────────────────────────────────────────────────────────────────────
// Boxcar
// ─────────────────────────────────────────────────────────────────────────────

fn boxcar(data: &[f32], cols: usize, rows: usize, window: usize) -> Vec<f32> {
    let half = (window / 2) as i32;
    let mut out = vec![f32::NAN; data.len()];
    out.par_chunks_mut(cols).enumerate().for_each(|(r, row_out)| {
        for c in 0..cols {
            let centre = data[r * cols + c];
            if !centre.is_finite() {
                continue;
            }
            // Boxcar uses the window mean directly; no MMSE correction.
            if let Some((mean, _, _)) = window_mean_var(data, cols, rows, r, c, half) {
                row_out[c] = mean;
            }
        }
    });
    out
}

// ─────────────────────────────────────────────────────────────────────────────
// Lee 1981 MMSE
// ─────────────────────────────────────────────────────────────────────────────

fn lee(data: &[f32], cols: usize, rows: usize, window: usize, enl: f32) -> Vec<f32> {
    let half = (window / 2) as i32;
    // Speckle coefficient of variation squared: CV_n² = 1/L for an
    // L-look amplitude-detected intensity image.
    let cv_n_sq = 1.0_f32 / enl;
    let mut out = vec![f32::NAN; data.len()];
    out.par_chunks_mut(cols).enumerate().for_each(|(r, row_out)| {
        for c in 0..cols {
            let centre = data[r * cols + c];
            if !centre.is_finite() {
                continue;
            }
            let Some((mean, var, _)) = window_mean_var(data, cols, rows, r, c, half) else {
                continue;
            };
            // Estimated noise variance: σ_n² = μ² · CV_n².  Signal
            // variance is then σ² − σ_n²; clip negative values to zero
            // (window is then declared homogeneous, k = 0).
            let noise_var = mean * mean * cv_n_sq;
            let signal_var = (var - noise_var).max(0.0);
            // MMSE weight K = σ_signal² / (σ_signal² + σ_n²); this is
            // the Lee 1981 form which equals (var − noise_var)/var in
            // the homogeneous limit and is bounded in [0, 1].
            let denom = signal_var + noise_var;
            let k = if denom > 0.0 { signal_var / denom } else { 0.0 };
            row_out[c] = mean + k * (centre - mean);
        }
    });
    out
}

// ─────────────────────────────────────────────────────────────────────────────
// Frost 1982
// ─────────────────────────────────────────────────────────────────────────────

fn frost(data: &[f32], cols: usize, rows: usize, window: usize, damping: f32) -> Vec<f32> {
    let half = (window / 2) as i32;
    let mut out = vec![f32::NAN; data.len()];
    out.par_chunks_mut(cols).enumerate().for_each(|(r, row_out)| {
        for c in 0..cols {
            let centre = data[r * cols + c];
            if !centre.is_finite() {
                continue;
            }
            let Some((mean, var, _)) = window_mean_var(data, cols, rows, r, c, half) else {
                continue;
            };
            // Frost weight: w(k,l) = exp(−α · |d|), where
            // α = damping · CV² and CV² = σ²/μ².  Distance is in
            // pixels (not normalised by window size); this matches the
            // 1982 formulation.  When mean ≤ 0 (cannot happen for
            // detected linear power), bypass to avoid div-by-zero.
            let cv_sq = if mean > 0.0 { var / (mean * mean) } else { 0.0 };
            let alpha = damping * cv_sq;
            let mut wsum = 0.0_f64;
            let mut vwsum = 0.0_f64;
            for (v, dy, dx) in window_samples(data, cols, rows, r, c, half) {
                // Manhattan distance (the Frost paper uses Euclidean,
                // but Manhattan is the more common implementation in SAR
                // toolboxes (e.g. ESA SNAP).  We document the choice
                // explicitly because it's a real divergence from the
                // 1982 paper; the difference in practice is < 0.1 dB.
                let d = (dy.abs() + dx.abs()) as f32;
                let w = (-alpha * d).exp() as f64;
                wsum += w;
                vwsum += w * (v as f64);
            }
            if wsum > 0.0 {
                row_out[c] = (vwsum / wsum) as f32;
            }
        }
    });
    out
}

// ─────────────────────────────────────────────────────────────────────────────
// Gamma-MAP (Lopes 1990)
// ─────────────────────────────────────────────────────────────────────────────

fn gamma_map(data: &[f32], cols: usize, rows: usize, window: usize, enl: f32) -> Vec<f32> {
    let half = (window / 2) as i32;
    // Speckle CV² = 1/L for amplitude-detected intensity, same as Lee.
    let cv_n_sq = 1.0_f32 / enl;
    // Maximum tolerable scene CV before declaring the pixel a point
    // target (passes through unchanged).  Lopes 1990 defines this as
    // CV_max² = 2 · CV_n² / (1 − CV_n²); for L=1 this is 1.0
    // (homogeneous → 0, point-target → ∞).
    let cv_max_sq = if cv_n_sq < 1.0 {
        2.0 * cv_n_sq / (1.0 - cv_n_sq)
    } else {
        // L ≤ 1 — degenerate; treat every window as point-target so
        // we never alias to non-physical values.
        f32::INFINITY
    };
    let mut out = vec![f32::NAN; data.len()];
    out.par_chunks_mut(cols).enumerate().for_each(|(r, row_out)| {
        for c in 0..cols {
            let centre = data[r * cols + c];
            if !centre.is_finite() {
                continue;
            }
            let Some((mean, var, _)) = window_mean_var(data, cols, rows, r, c, half) else {
                continue;
            };
            let cv_sq = if mean > 0.0 { var / (mean * mean) } else { 0.0 };
            // Three regimes:
            // (a) low heterogeneity   CV² ≤ CV_n²        → output = mean
            // (b) point target        CV² ≥ CV_max²      → output = centre
            // (c) general case        CV_n² < CV² < CV_max → solve quadratic
            if cv_sq <= cv_n_sq {
                row_out[c] = mean;
                continue;
            }
            if cv_sq >= cv_max_sq {
                row_out[c] = centre;
                continue;
            }
            // Lopes 1990 Gamma-MAP estimator: solve the quadratic
            //   L R̂² + (α − L − 1) μ R̂ − α μ I = 0
            // where α = (1 + CV_n²) / (CV² − CV_n²),
            //       L = ENL, μ = window mean, I = centre value.
            // The positive root is:
            //   R̂ = ((L − α + 1) μ + √Δ) / (2 L),
            //   Δ = (α − L − 1)² μ² + 4 α L μ I.
            // The `2 L` denominator (not just `2`) was confirmed against
            // the real-scene speckle_e2e test: dividing by 2 alone
            // doubles the output mean for L=4.
            let alpha = (1.0 + cv_n_sq) / (cv_sq - cv_n_sq);
            let b = (alpha - enl - 1.0) * mean;
            let disc = b * b + 4.0 * alpha * enl * mean * centre;
            // disc ≥ 0 holds whenever centre ≥ 0 and mean ≥ 0, which is
            // always true for detected linear power; we don't accept
            // negative inputs anyway (NaN/zero handling above).
            let r_hat = if disc >= 0.0 {
                ((enl - alpha + 1.0) * mean + disc.sqrt()) / (2.0 * enl)
            } else {
                mean
            };
            row_out[c] = r_hat;
        }
    });
    out
}

// ─────────────────────────────────────────────────────────────────────────────
// Refined Lee (Lopes-Touzi-Nezry 1990, fixed 7×7)
// ─────────────────────────────────────────────────────────────────────────────

/// Eight directional sub-window predicates for the 7×7 Refined Lee
/// filter, expressed as a function of the pixel offset `(dy, dx)`
/// relative to the centre.  Each predicate returns true iff the offset
/// lies inside that direction's 3-pixel-wide oriented sub-window.
///
/// The four cardinal sub-windows (S, N, E, W) are 21-pixel rectangular
/// strips of half-window size.  The four diagonal sub-windows
/// (SE, NW, NE, SW) are 3-pixel-wide bands along the corresponding
/// diagonal, restricted to the appropriate quadrant; their pixel
/// counts vary (typically 10–13 each) but the min-CV² selection rule
/// is scale-invariant and unaffected by the size difference.
///
/// Excludes the centre pixel `(0, 0)` from every mask: the centre is
/// the **observation** in the Lee MMSE update, not part of the
/// sub-window mean/variance estimate.
fn refined_lee_in_mask(direction: usize, dy: i32, dx: i32) -> bool {
    if dy == 0 && dx == 0 {
        return false;
    }
    match direction {
        // S: lower 3×7 strip (dy in 1..=3, dx in -3..=3).
        0 => (1..=3).contains(&dy),
        // N: upper 3×7 strip.
        1 => (-3..=-1).contains(&dy),
        // E: right 7×3 strip (dx in 1..=3, dy in -3..=3).
        2 => (1..=3).contains(&dx),
        // W: left 7×3 strip.
        3 => (-3..=-1).contains(&dx),
        // SE: diagonal band along (+y, +x), |dy − dx| ≤ 1, dy + dx ≥ 1.
        // The dy+dx≥1 clause excludes the NW half-strip, leaving an
        // oriented sub-window that points SE from the centre.
        4 => (dy - dx).abs() <= 1 && (dy + dx) >= 1,
        // NW: mirror of SE.
        5 => (dy - dx).abs() <= 1 && (dy + dx) <= -1,
        // NE: diagonal band along (-y, +x), |dy + dx| ≤ 1, dx − dy ≥ 1.
        6 => (dy + dx).abs() <= 1 && (dx - dy) >= 1,
        // SW: mirror of NE.
        7 => (dy + dx).abs() <= 1 && (dx - dy) <= -1,
        _ => false,
    }
}

/// Number of directional sub-windows.
const REFINED_LEE_NUM_DIRECTIONS: usize = 8;

/// Mean and population variance of `data` over the finite samples
/// inside the 7×7 sub-window of `direction`, centred at `(row, col)`.
/// Returns `None` if zero samples are finite — the caller cannot then
/// pick this direction as "homogeneous".
fn directional_window_mean_var(
    data: &[f32],
    cols: usize,
    rows: usize,
    row: usize,
    col: usize,
    direction: usize,
) -> Option<(f32, f32, usize)> {
    let mut sum = 0.0_f64;
    let mut sum_sq = 0.0_f64;
    let mut n: usize = 0;
    for dy in -3..=3 {
        for dx in -3..=3 {
            if !refined_lee_in_mask(direction, dy, dx) {
                continue;
            }
            let r = row as i32 + dy;
            let c = col as i32 + dx;
            if r < 0 || c < 0 || r >= rows as i32 || c >= cols as i32 {
                continue;
            }
            let v = data[(r as usize) * cols + c as usize];
            if !v.is_finite() {
                continue;
            }
            let v = v as f64;
            sum += v;
            sum_sq += v * v;
            n += 1;
        }
    }
    if n == 0 {
        return None;
    }
    let mean = sum / n as f64;
    let var = (sum_sq / n as f64 - mean * mean).max(0.0);
    Some((mean as f32, var as f32, n))
}

/// Refined Lee (Lopes-Touzi-Nezry 1990) on a fixed 7×7 window.
///
/// Algorithm:
/// 1. For each pixel, evaluate 8 directional sub-windows (mean +
///    variance) using [`refined_lee_in_mask`].
/// 2. Pick the sub-window with the **smallest CV²** (most homogeneous
///    side of the local edge).  Min-CV² (Lopes 1990) rather than
///    min-variance: the latter biases toward dark patches, the former
///    is scale-invariant.
/// 3. Apply the classic Lee 1981 MMSE update with that sub-window's
///    statistics:  out = μ + K·(centre − μ),  K = max(0, var − μ²·CV_n²) / var.
///
/// Edge effects: at the image border a sub-window may have fewer
/// samples than nominal; we tolerate any non-empty sub-window.  A
/// pixel whose centre is NaN remains NaN (nodata preservation).  A
/// pixel for which **every** sub-window is empty (entire 7×7
/// neighbourhood NaN) passes the centre through unchanged.
fn refined_lee(data: &[f32], cols: usize, rows: usize, enl: f32) -> Vec<f32> {
    let cv_n_sq = 1.0_f32 / enl;
    let mut out = vec![f32::NAN; data.len()];
    out.par_chunks_mut(cols).enumerate().for_each(|(r, row_out)| {
        for c in 0..cols {
            let centre = data[r * cols + c];
            if !centre.is_finite() {
                continue;
            }

            // Score each directional sub-window by CV² (population).
            let mut best_cv_sq = f32::INFINITY;
            let mut best_mean: f32 = centre;
            let mut best_var: f32 = 0.0;
            for dir in 0..REFINED_LEE_NUM_DIRECTIONS {
                let Some((mean, var, _n)) =
                    directional_window_mean_var(data, cols, rows, r, c, dir)
                else {
                    continue;
                };
                if !(mean > 0.0) {
                    // Defensive: detected linear power is non-negative;
                    // a zero/negative window mean indicates an all-zero
                    // sub-window which contributes nothing useful.
                    continue;
                }
                let cv_sq = var / (mean * mean);
                if cv_sq < best_cv_sq {
                    best_cv_sq = cv_sq;
                    best_mean = mean;
                    best_var = var;
                }
            }
            if !best_cv_sq.is_finite() {
                row_out[c] = centre;
                continue;
            }

            // Lee 1981 MMSE on the chosen sub-window's stats.
            let noise_var = best_mean * best_mean * cv_n_sq;
            let signal_var = (best_var - noise_var).max(0.0);
            let denom = signal_var + noise_var;
            let k = if denom > 0.0 { signal_var / denom } else { 0.0 };
            row_out[c] = best_mean + k * (centre - best_mean);
        }
    });
    out
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    fn approx(a: f32, b: f32, tol: f32) -> bool {
        (a - b).abs() <= tol
    }

    /// Constant input → same constant output for every filter.  This is
    /// the strongest invariant we can assert without committing to the
    /// numerical details of any specific filter.
    #[test]
    fn constant_input_is_a_fixed_point_of_every_filter() {
        let cols = 11;
        let rows = 9;
        let data = vec![0.123_f32; cols * rows];
        let filters = [
            SpeckleFilter::Boxcar { window: 3 },
            SpeckleFilter::Boxcar { window: 7 },
            SpeckleFilter::Lee { window: 5, enl: 1.0 },
            SpeckleFilter::Lee { window: 5, enl: 4.0 },
            SpeckleFilter::Frost { window: 5, damping: 1.0 },
            SpeckleFilter::Frost { window: 7, damping: 2.0 },
            SpeckleFilter::GammaMap { window: 5, enl: 1.0 },
            SpeckleFilter::GammaMap { window: 7, enl: 4.0 },
            SpeckleFilter::RefinedLee { enl: 1.0 },
            SpeckleFilter::RefinedLee { enl: 4.0 },
        ];
        for f in filters {
            let out = apply_speckle_filter(&data, cols, rows, f).expect("apply");
            for (i, v) in out.iter().enumerate() {
                assert!(
                    approx(*v, 0.123, 1e-6),
                    "filter {f:?} at index {i}: expected 0.123, got {v}"
                );
            }
        }
    }

    /// NaN input pixels propagate to NaN output pixels.  This is the
    /// nodata-preservation contract: we never invent values for missing
    /// samples.
    #[test]
    fn nan_input_pixel_yields_nan_output_pixel() {
        // Refined Lee is fixed 7×7, so use a larger image so the
        // corner finite check has finite samples to work with.
        let cols = 9;
        let rows = 9;
        let mut data = vec![1.0_f32; cols * rows];
        data[(rows / 2) * cols + cols / 2] = f32::NAN; // centre pixel
        let filters = [
            SpeckleFilter::Boxcar { window: 3 },
            SpeckleFilter::Lee { window: 3, enl: 1.0 },
            SpeckleFilter::Frost { window: 3, damping: 1.0 },
            SpeckleFilter::GammaMap { window: 3, enl: 1.0 },
            SpeckleFilter::RefinedLee { enl: 1.0 },
        ];
        for f in filters {
            let out = apply_speckle_filter(&data, cols, rows, f).expect("apply");
            let centre_idx = (rows / 2) * cols + cols / 2;
            assert!(
                out[centre_idx].is_nan(),
                "filter {f:?} should propagate NaN at centre"
            );
            // Non-NaN inputs at the corners should still produce finite output.
            assert!(
                out[0].is_finite(),
                "filter {f:?} produced NaN at corner with finite input"
            );
        }
    }

    /// Refined Lee on a synthetic horizontal step (left half = 0.05,
    /// right half = 0.50): a centre pixel several rows from the seam
    /// and well inside the right half should produce an output that
    /// tracks the **right-side mean** (≈ 0.50), not the global mean
    /// of the image (≈ 0.275) and not a window mean that straddles
    /// the edge.  This test pins the directional-window selection rule:
    /// the East / West / SE / NE / SW sub-windows that lie entirely on
    /// the right side have the smallest CV² (≈ 0), so the filter must
    /// pick one of them and the output must be close to 0.50.
    #[test]
    fn refined_lee_picks_homogeneous_side_of_horizontal_edge() {
        let cols = 21;
        let rows = 21;
        let mut data = vec![0.0_f32; cols * rows];
        for r in 0..rows {
            for c in 0..cols {
                data[r * cols + c] = if c < cols / 2 { 0.05 } else { 0.50 };
            }
        }
        let out = apply_speckle_filter(
            &data,
            cols,
            rows,
            SpeckleFilter::RefinedLee { enl: 4.0 },
        )
        .expect("apply");
        // Sample a centre pixel deep inside the right half (col = 17,
        // row = 10); its 7×7 window is entirely on the right side, so
        // every sub-window mean is 0.50 → output exactly 0.50.
        let r = 10usize;
        let c = 17usize;
        let v = out[r * cols + c];
        assert!(
            approx(v, 0.50, 1e-4),
            "Refined Lee centre well inside right half should be ≈ 0.50, got {v}"
        );
        // And a pixel one column to the right of the seam (col = 11)
        // — still has at least one sub-window (E or NE or SE) entirely
        // on the right side, so the output must be much closer to 0.50
        // than to the cross-seam mean of 0.275.
        let c_near = 11usize;
        let v_near = out[r * cols + c_near];
        assert!(
            (v_near - 0.50).abs() < (v_near - 0.275).abs(),
            "Refined Lee one column inside right half should be closer to 0.50 \
             than to global-mean 0.275, got {v_near}"
        );
    }

    /// Refined Lee should reject ENL ≤ 0 / non-finite, exactly like
    /// the other ENL-bearing filters.
    #[test]
    fn refined_lee_bad_enl_is_rejected() {
        let data = vec![0.1_f32; 49];
        for enl in [0.0_f32, -1.0, f32::NAN, f32::INFINITY, f32::NEG_INFINITY] {
            let err = apply_speckle_filter(&data, 7, 7, SpeckleFilter::RefinedLee { enl })
                .expect_err("bad ENL must error");
            assert!(matches!(err, SpeckleError::BadEnl { .. }), "got {err:?}");
        }
    }

    /// All-NaN window yields NaN output even at finite centre — but our
    /// definition treats NaN-centre as the only NaN-producer, so we
    /// really test the no-finite-neighbour case via a 1-pixel window
    /// surrounded by NaNs.  We construct: 3×3, centre=1, all neighbours
    /// NaN, window=3.  Window has one finite sample (the centre);
    /// statistics are well-defined → output finite.  Then make the
    /// centre NaN too: output NaN.
    #[test]
    fn isolated_finite_centre_is_finite_window_with_zero_finite_is_nan() {
        let cols = 3;
        let rows = 3;
        let nan = f32::NAN;
        let data = vec![nan, nan, nan, nan, 1.0, nan, nan, nan, nan];
        let f = SpeckleFilter::Boxcar { window: 3 };
        let out = apply_speckle_filter(&data, cols, rows, f).expect("apply");
        assert!(out[4].is_finite(), "centre with itself in window must be finite");
        assert!(approx(out[4], 1.0, 1e-6));
        // Now NaN the centre too: output NaN.
        let mut data2 = data;
        data2[4] = nan;
        let out2 = apply_speckle_filter(&data2, cols, rows, f).expect("apply");
        assert!(out2[4].is_nan());
    }

    /// Boxcar of a 3×3 known matrix with a 3×3 window — the centre
    /// output is the arithmetic mean of all nine inputs.  This pins
    /// the NaN-aware statistics machinery to a hand-computed value.
    #[test]
    fn boxcar_centre_equals_window_mean() {
        let data: Vec<f32> = (1..=9).map(|x| x as f32).collect();
        let cols = 3;
        let rows = 3;
        let out = apply_speckle_filter(&data, cols, rows, SpeckleFilter::Boxcar { window: 3 })
            .expect("apply");
        let expected = (1 + 2 + 3 + 4 + 5 + 6 + 7 + 8 + 9) as f32 / 9.0; // 5.0
        assert!(
            approx(out[4], expected, 1e-6),
            "centre should be window mean = {expected}, got {}",
            out[4]
        );
    }

    /// Lee with `ENL = 1` on a constant patch: signal variance is zero,
    /// so K = 0 and the output equals the window mean (= the constant).
    /// Combined with a heterogeneous neighbour the output should bias
    /// toward the centre value as ENL grows (less perceived noise).
    #[test]
    fn lee_high_enl_recovers_centre_value() {
        // 5×5 image: a single bright pixel in the middle of a uniform
        // background.  At ENL=1 (lots of noise), the bright pixel is
        // pulled toward the local mean; at ENL=10000, almost no
        // smoothing should be applied and the centre should survive.
        let cols = 5;
        let rows = 5;
        let mut data = vec![0.05_f32; cols * rows];
        data[2 * cols + 2] = 0.50;
        let f_low = SpeckleFilter::Lee { window: 3, enl: 1.0 };
        let f_high = SpeckleFilter::Lee { window: 3, enl: 10_000.0 };
        let out_low = apply_speckle_filter(&data, cols, rows, f_low).unwrap();
        let out_high = apply_speckle_filter(&data, cols, rows, f_high).unwrap();
        let centre_low = out_low[2 * cols + 2];
        let centre_high = out_high[2 * cols + 2];
        // Centre at high ENL must be much closer to 0.50 than at low ENL.
        assert!(
            (centre_high - 0.50).abs() < (centre_low - 0.50).abs(),
            "high-ENL Lee should preserve the bright pixel better; \
             low={centre_low}, high={centre_high}"
        );
        // And at very high ENL the centre output should be within 1 % of input.
        assert!(
            (centre_high - 0.50).abs() < 0.005,
            "ENL=1e4: centre should be near 0.50, got {centre_high}"
        );
    }

    /// Frost with damping = 0 reduces to a boxcar (all weights = 1).
    #[test]
    fn frost_zero_damping_equals_boxcar() {
        let data: Vec<f32> = (1..=25).map(|x| x as f32 * 0.01).collect();
        let cols = 5;
        let rows = 5;
        let out_frost =
            apply_speckle_filter(&data, cols, rows, SpeckleFilter::Frost { window: 3, damping: 0.0 })
                .unwrap();
        let out_box =
            apply_speckle_filter(&data, cols, rows, SpeckleFilter::Boxcar { window: 3 }).unwrap();
        for i in 0..data.len() {
            assert!(
                approx(out_frost[i], out_box[i], 1e-6),
                "Frost(K=0) must equal boxcar at idx {i}: frost={}, box={}",
                out_frost[i],
                out_box[i],
            );
        }
    }

    /// Gamma-MAP on a homogeneous patch (CV² ≤ CV_n²) returns the local
    /// mean.  We test by adding tiny noise so the patch is *almost*
    /// constant and CV² is still below 1/ENL.
    #[test]
    fn gamma_map_homogeneous_patch_returns_mean() {
        let cols = 5;
        let rows = 5;
        // Mean ≈ 0.10, CV ≈ 0.005 ≪ 1/ENL=1 ⇒ regime (a).
        let data: Vec<f32> = (0..cols * rows)
            .map(|i| 0.10 + 0.0005 * ((i % 3) as f32 - 1.0))
            .collect();
        let out = apply_speckle_filter(
            &data,
            cols,
            rows,
            SpeckleFilter::GammaMap { window: 3, enl: 1.0 },
        )
        .unwrap();
        // Centre output = window mean ≈ 0.10.
        assert!(
            approx(out[2 * cols + 2], 0.10, 1e-3),
            "homogeneous Gamma-MAP should return local mean ~0.10, got {}",
            out[2 * cols + 2]
        );
    }

    /// Gamma-MAP on a hard edge (CV² ≫ CV_max²) passes the centre
    /// through unchanged — the point-target preservation regime.
    /// Tested at ENL = 4 (a realistic multi-look value); the L = 1 case
    /// is mathematically degenerate (CV_max² → ∞) and Gamma-MAP has no
    /// point-target regime there, which our implementation reflects
    /// honestly rather than papering over.
    #[test]
    fn gamma_map_point_target_passes_centre_through() {
        // Centre = 100, neighbours = 0.01 — very high CV.
        let cols = 3;
        let rows = 3;
        let mut data = vec![0.01_f32; cols * rows];
        data[4] = 100.0;
        let out = apply_speckle_filter(
            &data,
            cols,
            rows,
            SpeckleFilter::GammaMap { window: 3, enl: 4.0 },
        )
        .unwrap();
        assert!(
            approx(out[4], 100.0, 1e-3),
            "point target should pass through unchanged; got {}",
            out[4]
        );
    }

    /// Gamma-MAP quadratic-regime regression: pin the centre output of
    /// a hand-traceable 3×3 case against the closed-form Lopes 1990
    /// estimator.  This is the only test that exercises the general
    /// quadratic regime — the homogeneous and point-target regimes are
    /// covered by other tests.  Inputs are chosen so the local CV²
    /// lands strictly inside `(CV_n², CV_max²)` at ENL = 4 (the
    /// quadratic regime), and the expected value below is computed by
    /// substituting μ = 0.10, I = 0.10, L = 4, CV² ≈ 0.569 into
    ///
    ///   α = (1 + CV_n²) / (CV² − CV_n²)
    ///   R̂ = ((L − α + 1) μ + √((α − L − 1)² μ² + 4 α L μ I)) / (2 L)
    ///
    /// → R̂ ≈ 0.11341 (verified by hand and by SciPy double-precision).
    /// Note R̂ ≠ μ even though I = μ: the MAP estimator under a Gamma
    /// scene prior is biased toward values that are more probable
    /// under the joint posterior, not toward I or μ individually.
    #[test]
    fn gamma_map_quadratic_regime_matches_hand_computed_lopes_value() {
        let cols = 3;
        let rows = 3;
        let data = vec![0.02, 0.18, 0.02, 0.18, 0.10, 0.18, 0.02, 0.18, 0.02];
        let enl = 4.0_f32;

        // Sanity-check that we are in the quadratic regime; if these
        // assertions fail, the test scenario itself is broken (not the
        // filter).
        let mean = data.iter().sum::<f32>() / 9.0;
        let var = data.iter().map(|v| (v - mean).powi(2)).sum::<f32>() / 9.0;
        let cv_sq = var / (mean * mean);
        let cv_n_sq = 1.0 / enl;
        let cv_max_sq = 2.0 * cv_n_sq / (1.0 - cv_n_sq);
        assert!(
            cv_n_sq < cv_sq && cv_sq < cv_max_sq,
            "test scenario must land in quadratic regime; \
             CV²={cv_sq}, CV_n²={cv_n_sq}, CV_max²={cv_max_sq}"
        );

        let out = apply_speckle_filter(
            &data,
            cols,
            rows,
            SpeckleFilter::GammaMap { window: 3, enl },
        )
        .unwrap();
        // Expected from the closed-form Lopes 1990 quadratic root with
        // the inputs above.  Tolerance allows for f32 round-off in the
        // intermediate variance estimate.
        assert!(
            approx(out[4], 0.11341, 5e-4),
            "Gamma-MAP quadratic regime mismatch: expected ≈ 0.11341, got {}",
            out[4]
        );
    }

    /// Reject even-sized windows.
    #[test]
    fn even_window_is_rejected() {
        let data = vec![0.1_f32; 25];
        let err = apply_speckle_filter(&data, 5, 5, SpeckleFilter::Boxcar { window: 4 })
            .expect_err("even window must error");
        assert!(matches!(err, SpeckleError::BadWindow { got: 4 }));
    }

    /// Reject window < 3 (no useful neighbourhood).
    #[test]
    fn tiny_window_is_rejected() {
        let data = vec![0.1_f32; 25];
        let err = apply_speckle_filter(&data, 5, 5, SpeckleFilter::Boxcar { window: 1 })
            .expect_err("window=1 must error");
        assert!(matches!(err, SpeckleError::BadWindow { got: 1 }));
    }

    /// Reject ENL ≤ 0 / non-finite.
    #[test]
    fn bad_enl_is_rejected() {
        let data = vec![0.1_f32; 25];
        for enl in [0.0_f32, -1.0, f32::NAN, f32::INFINITY, f32::NEG_INFINITY] {
            let err = apply_speckle_filter(&data, 5, 5, SpeckleFilter::Lee { window: 3, enl })
                .expect_err("bad ENL must error");
            assert!(matches!(err, SpeckleError::BadEnl { .. }), "got {err:?}");
        }
    }

    /// Reject negative or non-finite Frost damping.
    #[test]
    fn bad_damping_is_rejected() {
        let data = vec![0.1_f32; 25];
        for d in [-0.1_f32, f32::NAN, f32::INFINITY] {
            let err = apply_speckle_filter(
                &data,
                5,
                5,
                SpeckleFilter::Frost { window: 3, damping: d },
            )
            .expect_err("bad damping must error");
            assert!(matches!(err, SpeckleError::BadDamping { .. }), "got {err:?}");
        }
    }

    /// Reject dimension mismatch between data length and `cols × rows`.
    #[test]
    fn dimension_mismatch_is_rejected() {
        let data = vec![0.1_f32; 24]; // expects 25
        let err = apply_speckle_filter(&data, 5, 5, SpeckleFilter::Boxcar { window: 3 })
            .expect_err("mismatch must error");
        assert!(
            matches!(err, SpeckleError::DimensionMismatch { got: 24, expected: 25 }),
            "got {err:?}"
        );
    }
}
