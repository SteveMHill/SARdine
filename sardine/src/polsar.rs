//! Dual-polarisation H/A/Alpha decomposition (Cloude–Pottier 1997).
//!
//! Operates on the 2×2 Hermitian covariance matrix C2 formed from
//! co-registered, calibrated VV and VH complex SLC samples.
//!
//! # Physical model
//!
//! For a dual-pol (VV + VH) acquisition the scattering vector is
//! `k = [S_vv, S_vh]ᵀ`.  The covariance matrix is
//!
//! ```text
//! C2 = ⟨k · k†⟩ = [[⟨|S_vv|²⟩,   ⟨S_vv · S_vh*⟩],
//!                   [⟨S_vh · S_vv*⟩, ⟨|S_vh|²⟩   ]]
//! ```
//!
//! where `⟨·⟩` denotes spatial (multilook) averaging over a rectangular window.
//!
//! # Outputs
//!
//! Eigendecomposing `C2 = U diag(λ₁,λ₂) U†` (λ₁ ≥ λ₂ ≥ 0) gives:
//!
//! * **H** (Entropy) ∈ [0, 1]:  `-Σ pᵢ log₂ pᵢ`  where `pᵢ = λᵢ/(λ₁+λ₂)`.
//!   0 = single pure scattering mechanism; 1 = completely random.
//! * **A** (Anisotropy) ∈ [0, 1]:  `(λ₁ − λ₂)/(λ₁ + λ₂)`.
//!   Measures contrast between the two mechanisms.
//! * **α** (mean alpha angle) ∈ [0°, 90°]:  `Σ pᵢ αᵢ`  where
//!   `αᵢ = arccos(|uᵢ[0]|)`.  0° ≈ surface scattering; 90° ≈ double-bounce.
//!
//! # Limitations
//!
//! Sentinel-1 IW operates in dual-pol mode (VV+VH or HH+HV).  Full-pol
//! (T3/C3) decomposition requires quad-pol capability which S-1 IW does not
//! provide.  The dual-pol H/A/Alpha here is a 2×2 restricted version; its
//! physical interpretation is more limited than the quad-pol variant.

use chrono::{DateTime, Utc};
use rayon::prelude::*;

// ─── Output geometry type ─────────────────────────────────────────────────────

/// Dual-polarisation 2×2 covariance matrix image in slant-range geometry.
///
/// Stores the four independent real components of the Hermitian C2 matrix
/// averaged over a multilook window.  The geometry fields mirror
/// [`crate::merge_subswaths::MergedSigma0`] so the arrays can be fed
/// directly into the terrain-correction engine.
pub struct MergedC2 {
    /// Upper-left diagonal: ⟨|S_vv|²⟩ − NESZ_vv.  Length = `lines × samples`.
    pub c11: Vec<f32>,
    /// Lower-right diagonal: ⟨|S_vh|²⟩ − NESZ_vh.  Length = `lines × samples`.
    pub c22: Vec<f32>,
    /// Real part of off-diagonal: Re⟨S_vv · S_vh*⟩.  Length = `lines × samples`.
    pub c12_re: Vec<f32>,
    /// Imaginary part of off-diagonal: Im⟨S_vv · S_vh*⟩.  Length = `lines × samples`.
    pub c12_im: Vec<f32>,

    /// Number of azimuth lines.
    pub lines: usize,
    /// Number of range samples.
    pub samples: usize,

    /// Two-way slant range time of output column 0, in seconds.
    pub near_slant_range_time_s: f64,
    /// Two-way range time per pixel, in seconds.
    pub range_pixel_spacing_s: f64,
    /// Range pixel spacing in metres.
    pub range_pixel_spacing_m: f64,

    /// Zero-Doppler UTC time of line 0.
    pub azimuth_start_time: DateTime<Utc>,
    /// Azimuth time interval between consecutive lines, in seconds.
    pub azimuth_time_interval_s: f64,
}

// ─── Core eigendecomposition ──────────────────────────────────────────────────

/// Per-pixel closed-form 2×2 Hermitian eigendecomposition.
///
/// The input Hermitian matrix is:
/// ```text
/// C2 = [[c11,          c12_re + i·c12_im],
///       [c12_re − i·c12_im,  c22        ]]
/// ```
///
/// # Returns
///
/// `(lambda1, lambda2, alpha1_deg, alpha2_deg)` with `lambda1 ≥ lambda2 ≥ 0`.
/// Alpha angles are in degrees, in `[0°, 90°]`.
///
/// Returns `(0, 0, 0, 0)` when total power is zero or negative.
pub fn c2_eigendecomp(c11: f32, c22: f32, c12_re: f32, c12_im: f32) -> (f32, f32, f32, f32) {
    let total = c11 + c22;
    if total <= 0.0 {
        return (0.0, 0.0, 0.0, 0.0);
    }

    // Discriminant D = (c11 − c22)² + 4|c12|²
    let diff = c11 - c22;
    let c12_mag_sq = c12_re * c12_re + c12_im * c12_im;
    let disc_sq = (diff as f64 * diff as f64) + 4.0 * c12_mag_sq as f64;
    let disc = disc_sq.sqrt() as f32;

    // Eigenvalues (clamped to ≥ 0 to absorb floating-point rounding).
    let lambda1 = ((total as f64 + disc_sq.sqrt()) * 0.5) as f32;
    let lambda2 = ((total as f64 - disc_sq.sqrt()) * 0.5) as f32;
    let lambda1 = lambda1.max(0.0);
    let lambda2 = lambda2.max(0.0);

    // Eigenvector for λ₁: u₁ ∝ [c12_re + i·c12_im,  λ₁ − c11]
    //   λ₁ − c11 = (total + disc)/2 − c11 = (disc − diff)/2
    let lam1_off = (disc - diff) * 0.5; // = λ₁ − c11
    let norm1_sq = c12_mag_sq + lam1_off * lam1_off;

    let alpha1_deg: f32;
    let alpha2_deg: f32;

    if norm1_sq < f32::EPSILON * lambda1 * lambda1 {
        // Near-degenerate case: either equal eigenvalues (arbitrary eigenvectors)
        // or diagonal matrix with c11 ≥ c22 (causing lam1_off ≈ 0).
        // Convention: if c11 dominates, u₁ ≡ [1,0] → α₁ = 0°.
        if c11 >= c22 {
            alpha1_deg = 0.0;
            alpha2_deg = 90.0;
        } else {
            // c22 > c11 with norm1_sq ≈ 0 can only happen when eigenvalues are equal.
            alpha1_deg = 45.0;
            alpha2_deg = 45.0;
        }
    } else {
        // General case: compute alpha angles from eigenvector component magnitudes.
        // |u₁[0]|² = |c12|² / norm1_sq
        // |u₂[0]|² = |u₁[1]|² = lam1_off² / norm1_sq  (orthogonal eigenvector)
        let u1_cos_sq = (c12_mag_sq / norm1_sq).clamp(0.0, 1.0);
        let u2_cos_sq = (lam1_off * lam1_off / norm1_sq).clamp(0.0, 1.0);
        alpha1_deg = u1_cos_sq.sqrt().acos().to_degrees();
        alpha2_deg = u2_cos_sq.sqrt().acos().to_degrees();
    }

    (lambda1, lambda2, alpha1_deg, alpha2_deg)
}

// ─── H/A/alpha computation ────────────────────────────────────────────────────

/// Compute H (Entropy), A (Anisotropy), and mean α angle from one C2 pixel.
///
/// Returns `(NaN, NaN, NaN)` when total power is zero or any input is NaN.
///
/// # Arguments
///
/// * `c11`    — Real VV power (diagonal, noise-subtracted).
/// * `c22`    — Real VH power (diagonal, noise-subtracted).
/// * `c12_re` — Real part of VV×VH* cross-product.
/// * `c12_im` — Imaginary part of VV×VH* cross-product.
#[inline]
pub fn compute_haa(c11: f32, c22: f32, c12_re: f32, c12_im: f32) -> (f32, f32, f32) {
    if c11.is_nan() || c22.is_nan() || c12_re.is_nan() || c12_im.is_nan() {
        return (f32::NAN, f32::NAN, f32::NAN);
    }

    let (lambda1, lambda2, alpha1_deg, alpha2_deg) = c2_eigendecomp(c11, c22, c12_re, c12_im);
    let total = lambda1 + lambda2;
    if total <= 0.0 {
        return (f32::NAN, f32::NAN, f32::NAN);
    }

    let p1 = lambda1 / total;
    let p2 = lambda2 / total;

    // Entropy H = −Σ pᵢ log₂ pᵢ  ∈ [0, 1] (log base 2 normalises to 1 for 2 eigenvalues).
    let h_term = |p: f32| if p > 0.0 { -p * p.log2() } else { 0.0 };
    let entropy = (h_term(p1) + h_term(p2)).clamp(0.0, 1.0);

    // Anisotropy A = (λ₁ − λ₂) / (λ₁ + λ₂)  ∈ [0, 1].
    let anisotropy = ((lambda1 - lambda2) / total).clamp(0.0, 1.0);

    // Mean alpha angle  ∈ [0°, 90°].
    let mean_alpha = (p1 * alpha1_deg + p2 * alpha2_deg).clamp(0.0, 90.0);

    (entropy, anisotropy, mean_alpha)
}

// ─── C2 multilook ─────────────────────────────────────────────────────────────

/// Spatially average a full-resolution C2 image by `rg_looks × az_looks`.
///
/// All four C2 components are averaged independently.  NaN pixels in any
/// component contribute nothing to the average; if all pixels in a window
/// cell are NaN, the output is NaN for that cell.
///
/// Updates the geometry fields accordingly:
/// * `lines` → `lines / az_looks`
/// * `samples` → `samples / rg_looks`
/// * `azimuth_time_interval_s` → `azimuth_time_interval_s * az_looks`
pub fn multilook_c2(c2: &MergedC2, rg_looks: usize, az_looks: usize) -> MergedC2 {
    assert!(rg_looks >= 1 && az_looks >= 1, "look counts must be ≥ 1");

    let out_lines = c2.lines / az_looks;
    let out_samples = c2.samples / rg_looks;
    let n = out_lines * out_samples;

    // Allocate output arrays.
    let mut out_c11 = vec![0.0f32; n];
    let mut out_c22 = vec![0.0f32; n];
    let mut out_c12_re = vec![0.0f32; n];
    let mut out_c12_im = vec![0.0f32; n];

    // Parallel over output lines.
    out_c11
        .par_chunks_mut(out_samples)
        .zip(out_c22.par_chunks_mut(out_samples))
        .zip(out_c12_re.par_chunks_mut(out_samples))
        .zip(out_c12_im.par_chunks_mut(out_samples))
        .enumerate()
        .for_each(|(out_l, (((row_c11, row_c22), row_c12_re), row_c12_im))| {
            for out_s in 0..out_samples {
                let mut sum_c11 = 0.0f64;
                let mut sum_c22 = 0.0f64;
                let mut sum_c12_re = 0.0f64;
                let mut sum_c12_im = 0.0f64;
                let mut count = 0u32;

                for al in 0..az_looks {
                    let in_l = out_l * az_looks + al;
                    for rl in 0..rg_looks {
                        let in_s = out_s * rg_looks + rl;
                        let idx = in_l * c2.samples + in_s;
                        let v = c2.c11[idx];
                        if v.is_nan() {
                            continue;
                        }
                        sum_c11 += v as f64;
                        sum_c22 += c2.c22[idx] as f64;
                        sum_c12_re += c2.c12_re[idx] as f64;
                        sum_c12_im += c2.c12_im[idx] as f64;
                        count += 1;
                    }
                }

                let out_idx = out_s; // chunk is already one row
                if count == 0 {
                    row_c11[out_idx] = f32::NAN;
                    row_c22[out_idx] = f32::NAN;
                    row_c12_re[out_idx] = f32::NAN;
                    row_c12_im[out_idx] = f32::NAN;
                } else {
                    let scale = 1.0 / count as f64;
                    row_c11[out_idx] = (sum_c11 * scale) as f32;
                    row_c22[out_idx] = (sum_c22 * scale) as f32;
                    row_c12_re[out_idx] = (sum_c12_re * scale) as f32;
                    row_c12_im[out_idx] = (sum_c12_im * scale) as f32;
                }
            }
        });

    MergedC2 {
        c11: out_c11,
        c22: out_c22,
        c12_re: out_c12_re,
        c12_im: out_c12_im,
        lines: out_lines,
        samples: out_samples,
        near_slant_range_time_s: c2.near_slant_range_time_s,
        range_pixel_spacing_s: c2.range_pixel_spacing_s * rg_looks as f64,
        range_pixel_spacing_m: c2.range_pixel_spacing_m * rg_looks as f64,
        azimuth_start_time: c2.azimuth_start_time,
        azimuth_time_interval_s: c2.azimuth_time_interval_s * az_looks as f64,
    }
}

/// Expand a MergedC2 into three flat arrays of H, A, and alpha (in degrees).
///
/// Applies [`compute_haa`] to every pixel in parallel.  The returned arrays
/// have the same flat row-major layout as the input C2 arrays.
pub fn compute_haa_image(c2: &MergedC2) -> (Vec<f32>, Vec<f32>, Vec<f32>) {
    let n = c2.lines * c2.samples;
    let mut h_buf = vec![f32::NAN; n];
    let mut a_buf = vec![f32::NAN; n];
    let mut alpha_buf = vec![f32::NAN; n];

    h_buf
        .par_iter_mut()
        .zip(a_buf.par_iter_mut())
        .zip(alpha_buf.par_iter_mut())
        .enumerate()
        .for_each(|(i, ((h, a), al))| {
            let (hv, av, alv) = compute_haa(c2.c11[i], c2.c22[i], c2.c12_re[i], c2.c12_im[i]);
            *h = hv;
            *a = av;
            *al = alv;
        });

    (h_buf, a_buf, alpha_buf)
}

// ─── Tests ────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    fn approx_eq(a: f32, b: f32, tol: f32) -> bool {
        if a.is_nan() && b.is_nan() {
            return true;
        }
        (a - b).abs() <= tol
    }

    // ── eigendecomp ──────────────────────────────────────────────────────────

    #[test]
    fn eigendecomp_diagonal_c11_dominant() {
        // c12 = 0, c11 > c22 → λ₁ = c11, α₁ = 0°, α₂ = 90°
        let (l1, l2, a1, a2) = c2_eigendecomp(4.0, 1.0, 0.0, 0.0);
        assert!(approx_eq(l1, 4.0, 1e-4), "lambda1={l1}");
        assert!(approx_eq(l2, 1.0, 1e-4), "lambda2={l2}");
        assert!(approx_eq(a1, 0.0, 1e-3), "alpha1={a1}");
        assert!(approx_eq(a2, 90.0, 1e-3), "alpha2={a2}");
    }

    #[test]
    fn eigendecomp_diagonal_c22_dominant() {
        // c12 = 0, c22 > c11 → λ₁ = c22, α₁ = 90° (VH-dominant)
        let (l1, l2, a1, a2) = c2_eigendecomp(1.0, 4.0, 0.0, 0.0);
        assert!(approx_eq(l1, 4.0, 1e-4), "lambda1={l1}");
        assert!(approx_eq(l2, 1.0, 1e-4), "lambda2={l2}");
        // VH dominant: u₁ = [0,1], so α₁ = arccos(0) = 90°
        assert!(approx_eq(a1, 90.0, 1e-3), "alpha1={a1}");
        assert!(approx_eq(a2, 0.0, 1e-3), "alpha2={a2}");
    }

    #[test]
    fn eigendecomp_equal_diagonal_no_offdiag() {
        // c11 == c22, c12 == 0 → equal eigenvalues
        let (l1, l2, a1, a2) = c2_eigendecomp(2.0, 2.0, 0.0, 0.0);
        assert!(approx_eq(l1, 2.0, 1e-4));
        assert!(approx_eq(l2, 2.0, 1e-4));
        // Eigenvectors arbitrary — convention gives 0°/90° (c11 >= c22)
        assert!(a1 >= 0.0 && a1 <= 90.0, "a1={a1}");
        assert!(a2 >= 0.0 && a2 <= 90.0, "a2={a2}");
    }

    #[test]
    fn eigendecomp_known_offdiag() {
        // C2 = [[2, 1+i], [1-i, 2]]
        // Eigenvalues: λ = 2 ± |1+i| = 2 ± √2 ≈ 3.414, 0.586
        let (l1, l2, _a1, _a2) = c2_eigendecomp(2.0, 2.0, 1.0, 1.0);
        let expected_l1 = 2.0 + 2.0_f32.sqrt();
        let expected_l2 = 2.0 - 2.0_f32.sqrt();
        assert!(approx_eq(l1, expected_l1, 1e-4), "l1={l1} expected={expected_l1}");
        assert!(approx_eq(l2, expected_l2, 1e-4), "l2={l2} expected={expected_l2}");
        // λ₂ ≥ 0
        assert!(l2 >= 0.0, "l2 must be non-negative, got {l2}");
    }

    #[test]
    fn eigendecomp_zero_power_returns_zeros() {
        let (l1, l2, a1, a2) = c2_eigendecomp(0.0, 0.0, 0.0, 0.0);
        assert_eq!(l1, 0.0);
        assert_eq!(l2, 0.0);
        assert_eq!(a1, 0.0);
        assert_eq!(a2, 0.0);
    }

    // ── compute_haa ──────────────────────────────────────────────────────────

    #[test]
    fn haa_single_mechanism_max_anisotropy() {
        // One eigenvalue dominates → A ≈ 1, H ≈ 0
        let (h, a, _alpha) = compute_haa(100.0, 0.001, 0.0, 0.0);
        assert!(a > 0.99, "A should be ≈1 for single dominant mechanism, got {a}");
        assert!(h < 0.02, "H should be ≈0 for single dominant mechanism, got {h}");
    }

    #[test]
    fn haa_equal_eigenvalues_max_entropy() {
        // Equal eigenvalues → H = 1, A = 0
        let (h, a, alpha) = compute_haa(1.0, 1.0, 0.0, 0.0);
        assert!(approx_eq(h, 1.0, 1e-5), "H should be 1.0 for equal eigenvalues, got {h}");
        assert!(approx_eq(a, 0.0, 1e-5), "A should be 0.0 for equal eigenvalues, got {a}");
        assert!(alpha >= 0.0 && alpha <= 90.0, "alpha out of range: {alpha}");
    }

    #[test]
    fn haa_output_ranges() {
        // Random-ish valid C2 matrices; verify all outputs are in range.
        let test_cases = [
            (3.0_f32, 1.0_f32, 0.5_f32, 0.5_f32),
            (10.0, 0.1, 1.0, -0.5),
            (5.0, 5.0, 2.0, 2.0),
            (1.0, 3.0, 0.0, 1.0),
            (0.5, 0.5, 0.3, 0.0),
        ];
        for (c11, c22, c12_re, c12_im) in test_cases {
            let (h, a, alpha) = compute_haa(c11, c22, c12_re, c12_im);
            assert!((0.0..=1.0).contains(&h), "H out of range: {h} for ({c11},{c22},{c12_re},{c12_im})");
            assert!((0.0..=1.0).contains(&a), "A out of range: {a}");
            assert!((0.0..=90.0).contains(&alpha), "alpha out of range: {alpha}");
        }
    }

    #[test]
    fn haa_nan_on_zero_power() {
        let (h, a, alpha) = compute_haa(0.0, 0.0, 0.0, 0.0);
        assert!(h.is_nan(), "H should be NaN for zero-power pixel");
        assert!(a.is_nan());
        assert!(alpha.is_nan());
    }

    #[test]
    fn haa_nan_input_propagates() {
        let (h, a, alpha) = compute_haa(f32::NAN, 1.0, 0.0, 0.0);
        assert!(h.is_nan());
        assert!(a.is_nan());
        assert!(alpha.is_nan());
    }

    // ── multilook_c2 ─────────────────────────────────────────────────────────

    #[test]
    fn multilook_c2_geometry_update() {
        use chrono::TimeZone;
        let epoch = Utc.with_ymd_and_hms(2019, 1, 23, 5, 33, 0).unwrap();
        let c2 = MergedC2 {
            c11: vec![1.0f32; 40 * 20],
            c22: vec![0.5f32; 40 * 20],
            c12_re: vec![0.1f32; 40 * 20],
            c12_im: vec![0.0f32; 40 * 20],
            lines: 40,
            samples: 20,
            near_slant_range_time_s: 0.005,
            range_pixel_spacing_s: 1e-8,
            range_pixel_spacing_m: 2.33,
            azimuth_start_time: epoch,
            azimuth_time_interval_s: 0.002,
        };
        let ml = multilook_c2(&c2, 4, 2);
        assert_eq!(ml.lines, 20);
        assert_eq!(ml.samples, 5);
        assert!((ml.azimuth_time_interval_s - 0.004).abs() < 1e-12);
        assert!((ml.range_pixel_spacing_m - 2.33 * 4.0).abs() < 1e-6);
        // All-uniform input → averages stay constant.
        assert!((ml.c11[0] - 1.0).abs() < 1e-5);
        assert!((ml.c22[0] - 0.5).abs() < 1e-5);
    }
}
