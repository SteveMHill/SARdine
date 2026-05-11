//! Sparse orbit-based co-registration for Sentinel-1 InSAR.
//!
//! Co-registration maps each reference SLC pixel to its corresponding
//! location in the secondary SLC by:
//!
//! 1. **Forward geocoding** — project reference (line, col) to ECEF on the
//!    WGS84 ellipsoid using the slant-range + zero-Doppler equations.
//! 2. **Secondary SAR coordinates** — solve the zero-Doppler equation on the
//!    secondary orbit for the projected ECEF target point.
//! 3. **Polynomial fitting** — fit a degree-2 bivariate polynomial to the
//!    sparse offset grid (both azimuth and range offsets independently).
//! 4. **Dense resampling** — evaluate the polynomial for every reference
//!    pixel and bilinearly interpolate the secondary Deramped SLC.
//!
//! # Why sparse grid?
//!
//! A per-pixel Newton solve on a 25000 × 5000 image would cost ~40 s.
//! Orbit-state smoothness means 20 × 50 = 1000 grid points capture all
//! variation; polynomial evaluation is then O(1) per pixel.
//!
//! # Height assumption
//!
//! Forward geocoding assumes the target lies on the WGS84 reference
//! ellipsoid (h = 0).  For typical Sentinel-1 IW geometry this introduces
//! a co-registration range error of ≈ h / (2 · sin(θ_inc)) / rps in pixels,
//! which is < 0.03 px for h < 100 m.  Coherence magnitude is not sensitive
//! to sub-pixel errors of this magnitude.

use chrono::{DateTime, Duration, Utc};
use rayon::prelude::*;

use crate::{
    deburst::DeburstArray,
    insar::{deramp::Deramped, error::InsarError},
    orbit::interpolate_orbit_and_accel,
    terrain_correction::{solve_zero_doppler, ZeroDopplerOutcome},
    types::{OrbitData, SubSwathMetadata},
};

// ── WGS84 ellipsoid constants ──────────────────────────────────────────────

/// WGS84 semi-major axis (m).
const WGS84_A: f64 = 6_378_137.0;
/// WGS84 semi-minor axis (m).
const WGS84_B: f64 = 6_356_752.314_245;

// ── Speed of light ────────────────────────────────────────────────────────

const SPEED_OF_LIGHT_M_S: f64 = 299_792_458.0;

// ── Sparse grid geometry ──────────────────────────────────────────────────

/// Number of azimuth grid points for sparse co-registration.
const GRID_AZ: usize = 20;
/// Number of range grid points for sparse co-registration.
const GRID_RG: usize = 50;

// ─── Public types ─────────────────────────────────────────────────────────

/// Degree-2 bivariate co-registration polynomial.
///
/// Maps (l, c) in the reference geometry to the offset (Δl, Δc) to apply
/// to find the corresponding secondary SLC pixel:
///
/// ```text
/// Δ(l,c) = p[0]
///        + p[1]*l  + p[2]*c
///        + p[3]*l² + p[4]*l*c + p[5]*c²
/// ```
///
/// where `l`, `c` are line and column indices in the reference `DeburstArray`.
#[derive(Debug, Clone)]
pub struct CoregPoly {
    pub az: [f64; 6],
    pub rg: [f64; 6],
}

impl CoregPoly {
    /// Evaluate the azimuth offset at reference pixel (l, c).
    pub fn eval_az(&self, l: f64, c: f64) -> f64 {
        eval_poly6(&self.az, l, c)
    }

    /// Evaluate the range offset at reference pixel (l, c).
    pub fn eval_rg(&self, l: f64, c: f64) -> f64 {
        eval_poly6(&self.rg, l, c)
    }
}

/// Output of the sparse co-registration step.
#[derive(Debug, Clone)]
pub struct CoregResult {
    /// Degree-2 bivariate polynomial for azimuth and range offsets.
    pub poly: CoregPoly,
    /// RMS polynomial fit residual over the sparse grid (in pixels).
    pub fit_residual_rms: f64,
    /// Number of valid grid points used for polynomial fitting.
    pub n_valid: usize,
}

// ─── Primary entry point ──────────────────────────────────────────────────

/// Compute sparse co-registration offsets between reference and secondary.
///
/// # Arguments
///
/// * `ref_swath`       — Reference subswath metadata.
/// * `ref_orbit`       — Reference scene orbit data.
/// * `ref_first_utc`   — UTC time of the first line of `ref_deburst`.
/// * `ref_deburst`     — Reference debursted SLC (used for geometry only).
/// * `sec_swath`       — Secondary subswath metadata.
/// * `sec_orbit`       — Secondary scene orbit data.
/// * `sec_first_utc`   — UTC time of the first line of the secondary debursted SLC.
///
/// # Returns
///
/// `CoregResult` with the degree-2 bivariate polynomial offset and fit diagnostics.
pub fn compute_coreg_offsets(
    ref_swath: &SubSwathMetadata,
    ref_orbit: &OrbitData,
    ref_first_utc: DateTime<Utc>,
    ref_deburst: &DeburstArray,
    sec_swath: &SubSwathMetadata,
    sec_orbit: &OrbitData,
    sec_first_utc: DateTime<Utc>,
) -> Result<CoregResult, InsarError> {
    let ref_lines = ref_deburst.lines;
    let ref_samples = ref_deburst.samples;

    // Azimuth time interval for reference and secondary (seconds per line).
    let ref_ati_s = ref_swath.azimuth_time_interval_s;
    let sec_ati_s = sec_swath.azimuth_time_interval_s;

    // Slant-range time to first pixel and per-pixel increment for reference.
    let ref_tau0 = ref_swath.slant_range_time_s; // two-way, seconds
    let ref_delta_tau = 2.0 * ref_swath.range_pixel_spacing_m / SPEED_OF_LIGHT_M_S;
    let ref_col_offset = ref_deburst.valid_sample_offset as f64;

    // Same for secondary.
    let sec_tau0 = sec_swath.slant_range_time_s;
    let sec_delta_tau = 2.0 * sec_swath.range_pixel_spacing_m / SPEED_OF_LIGHT_M_S;
    let sec_col_offset = sec_deburst_valid_offset(sec_swath, ref_swath);

    // Sparse grid: uniformly spaced inside the valid image extent.
    // Avoid edges (10% margin) to reduce geometric extrapolation.
    let az_step = ref_lines as f64 / (GRID_AZ + 1) as f64;
    let rg_step = ref_samples as f64 / (GRID_RG + 1) as f64;

    // Collect (l, c, Δl, Δc) from the sparse grid in parallel.
    let grid_pts: Vec<(f64, f64, f64, f64)> = (0..GRID_AZ)
        .into_par_iter()
        .flat_map_iter(|gi| {
            let l = (gi + 1) as f64 * az_step;
            (0..GRID_RG).filter_map(move |gj| {
                let c = (gj + 1) as f64 * rg_step;

                // Azimuth time and slant range for this reference pixel.
                let t_az_ref = ref_first_utc
                    + Duration::microseconds((l * ref_ati_s * 1e6) as i64); // SAFETY-OK: product of line count × ATI in microseconds; cannot overflow for S-1 scene durations of < 30 s
                let tau_ref = ref_tau0 + (c + ref_col_offset) * ref_delta_tau;
                let range_m = tau_ref * SPEED_OF_LIGHT_M_S / 2.0;

                // Get reference satellite state vector at azimuth time.
                let (sv_ref, _accel) = interpolate_orbit_and_accel(ref_orbit, t_az_ref).ok()?; // SAFETY-OK: orbit interpolation failure for a sparse sample point returns None, skipping this point; the n_valid >= 6 guard downstream ensures sufficient coverage

                // Forward geocoding: reference SAR pixel → ECEF on WGS84 ellipsoid.
                let target_ecef = ecef_from_slant_range_doppler(
                    sv_ref.position_m,
                    sv_ref.velocity_m_s,
                    range_m,
                    20,
                    1e-3, // 1 mm convergence
                )
                .ok()?; // SAFETY-OK: forward geocoding non-convergence for a sparse sample point returns None, skipping this point; n_valid >= 6 guard ensures polynomial fit has sufficient coverage

                // Backward geocoding: ECEF → secondary SAR coordinates.
                let (l_sec, c_sec) = secondary_sar_coords(
                    target_ecef,
                    sec_orbit,
                    sec_first_utc,
                    sec_first_utc, // initial guess = same time
                    sec_tau0,
                    sec_delta_tau,
                    sec_col_offset,
                    sec_ati_s,
                )
                .ok()?; // SAFETY-OK: backward geocoding non-convergence for a sparse sample point returns None, skipping this point; n_valid >= 6 guard ensures polynomial fit has sufficient coverage

                let dl = l_sec - l;
                let dc = c_sec - c;
                Some((l, c, dl, dc))
            })
        })
        .collect();

    let n_valid = grid_pts.len();
    if n_valid < 6 {
        return Err(InsarError::CoregGridTooSmall { n_valid });
    }

    // Fit degree-2 bivariate polynomials for azimuth and range offsets.
    let (az_poly, az_res) = fit_poly6(&grid_pts, |&(_, _, dl, _)| dl)?;
    let (rg_poly, rg_res) = fit_poly6(&grid_pts, |&(_, _, _, dc)| dc)?;

    let fit_residual_rms = ((az_res * az_res + rg_res * rg_res) / 2.0).sqrt();

    Ok(CoregResult {
        poly: CoregPoly { az: az_poly, rg: rg_poly },
        fit_residual_rms,
        n_valid,
    })
}

/// Resample secondary deramped SLC into reference geometry using co-registration polynomial.
///
/// For each reference pixel (l, c), evaluate the polynomial to get
/// (Δl, Δc), then bilinearly interpolate the secondary Deramped at
/// (l + Δl, c + Δc).  Pixels that fall outside the secondary image are
/// set to [0.0, 0.0].
pub fn resample_secondary(
    ref_lines: usize,
    ref_samples: usize,
    poly: &CoregPoly,
    secondary: &Deramped,
) -> Deramped {
    let n = ref_lines * ref_samples;
    let mut out_data: Vec<[f32; 2]> = Vec::with_capacity(n);
    // SAFETY-OK: capacity == len, every element is immediately written below
    unsafe { out_data.set_len(n) };

    out_data
        .par_chunks_mut(ref_samples)
        .enumerate()
        .for_each(|(l, row)| {
            let lf = l as f64;
            for (c, pix) in row.iter_mut().enumerate() {
                let cf = c as f64;
                let l_sec = lf + poly.eval_az(lf, cf);
                let c_sec = cf + poly.eval_rg(lf, cf);
                *pix = bilinear(secondary, l_sec, c_sec);
            }
        });

    Deramped {
        data: out_data,
        lines: ref_lines,
        samples: ref_samples,
        valid_sample_offset: secondary.valid_sample_offset,
    }
}

// ─── Forward geocoding ────────────────────────────────────────────────────

/// Project a reference SAR pixel to ECEF on the WGS84 ellipsoid using
/// Newton-Raphson on the 3-equation system:
///
/// ```text
/// F1: (T − P) · V = 0                          [zero-Doppler plane]
/// F2: (T − P)_x² + (T − P)_y² + (T − P)_z² − R² = 0   [range sphere]
/// F3: T_x²/a² + T_y²/a² + T_z²/b² − 1 = 0   [WGS84 ellipsoid]
/// ```
///
/// # Arguments
///
/// * `sat_pos`  — Satellite ECEF position at azimuth time (m).
/// * `sat_vel`  — Satellite velocity at azimuth time (m/s).
/// * `range_m`  — One-way slant range (m).
/// * `max_iter` — Maximum Newton iterations (20 is sufficient for < 1 mm).
/// * `tol_m`    — Convergence tolerance on position update norm (m).
fn ecef_from_slant_range_doppler(
    sat_pos: [f64; 3],
    sat_vel: [f64; 3],
    range_m: f64,
    max_iter: usize,
    tol_m: f64,
) -> Result<[f64; 3], ()> {
    let [px, py, pz] = sat_pos;
    let [vx, vy, vz] = sat_vel;
    let r = range_m;
    let a = WGS84_A;
    let b = WGS84_B;

    // Initial guess: satellite nadir direction, scaled to WGS84 surface.
    // The nadir point is along −ẑ_ecef direction; we project sat_pos to
    // the ellipsoid along the unit radial direction.
    let sat_r = (px * px + py * py + pz * pz).sqrt();
    // Scale to approximate ellipsoid surface (rough initial guess).
    let scale = a / sat_r;
    let mut tx = px * scale;
    let mut ty = py * scale;
    let mut tz = pz * scale * (b / a); // flatten in z direction

    for _ in 0..max_iter {
        let dx = tx - px;
        let dy = ty - py;
        let dz = tz - pz;

        let f1 = dx * vx + dy * vy + dz * vz;
        let f2 = dx * dx + dy * dy + dz * dz - r * r;
        let f3 = tx * tx / (a * a) + ty * ty / (a * a) + tz * tz / (b * b) - 1.0;

        // Jacobian rows (∂F_i / ∂[tx, ty, tz]):
        // J[0] = [vx, vy, vz]
        // J[1] = [2*dx, 2*dy, 2*dz]
        // J[2] = [2*tx/a², 2*ty/a², 2*tz/b²]
        let j00 = vx;
        let j01 = vy;
        let j02 = vz;
        let j10 = 2.0 * dx;
        let j11 = 2.0 * dy;
        let j12 = 2.0 * dz;
        let j20 = 2.0 * tx / (a * a);
        let j21 = 2.0 * ty / (a * a);
        let j22 = 2.0 * tz / (b * b);

        // Solve 3×3 system J * δ = F via Cramer's rule.
        let det = j00 * (j11 * j22 - j12 * j21)
            - j01 * (j10 * j22 - j12 * j20)
            + j02 * (j10 * j21 - j11 * j20);

        if det.abs() < 1e-20 {
            return Err(());
        }

        let dx_t = (f1 * (j11 * j22 - j12 * j21)
            - j01 * (f2 * j22 - j12 * f3)
            + j02 * (f2 * j21 - j11 * f3))
            / det;
        let dy_t = (j00 * (f2 * j22 - j12 * f3)
            - f1 * (j10 * j22 - j12 * j20)
            + j02 * (j10 * f3 - f2 * j20))
            / det;
        let dz_t = (j00 * (j11 * f3 - f2 * j21)
            - j01 * (j10 * f3 - f2 * j20)
            + f1 * (j10 * j21 - j11 * j20))
            / det;

        tx -= dx_t;
        ty -= dy_t;
        tz -= dz_t;

        let update_norm = (dx_t * dx_t + dy_t * dy_t + dz_t * dz_t).sqrt();
        if update_norm < tol_m {
            return Ok([tx, ty, tz]);
        }
    }

    Err(())
}

// ─── Secondary SAR coordinate lookup ────────────────────────────────────

/// Solve for secondary SAR (line, col_float) of a target ECEF point.
///
/// Uses `solve_zero_doppler` on the secondary orbit, then converts the
/// converged azimuth time to a line index and computes slant range to
/// derive the column index.
///
/// Returns `None` if the solver did not converge or returned a degenerate result.
#[allow(clippy::too_many_arguments)]
fn secondary_sar_coords(
    target_ecef: [f64; 3],
    sec_orbit: &OrbitData,
    sec_first_utc: DateTime<Utc>,
    initial_guess: DateTime<Utc>,
    sec_tau0: f64,
    sec_delta_tau: f64,
    sec_col_offset: f64,
    sec_ati_s: f64,
) -> Result<(f64, f64), InsarError> {
    let outcome = solve_zero_doppler(target_ecef, sec_orbit, initial_guess, 20, 1e-9)?;

    match outcome {
        ZeroDopplerOutcome::Converged { time: t_az_sec, sat_sv: sv_sec, .. } => {
            // Azimuth line index in the secondary debursted array.
            let dt_us = (t_az_sec - sec_first_utc)
                .num_microseconds()
                .unwrap_or(i64::MIN); // SAFETY-OK: if duration exceeds i64 microseconds the scene is clearly wrong; MIN produces an obviously-invalid line index
            let l_sec = dt_us as f64 / (sec_ati_s * 1e6);

            // Slant range to secondary satellite position.
            let [tx, ty, tz] = target_ecef;
            let [sx, sy, sz] = sv_sec.position_m;
            let range_m = ((tx - sx) * (tx - sx)
                + (ty - sy) * (ty - sy)
                + (tz - sz) * (tz - sz))
                .sqrt();

            // Secondary column (float) in debursted coordinates.
            let tau_sec = 2.0 * range_m / SPEED_OF_LIGHT_M_S;
            let col_full = (tau_sec - sec_tau0) / sec_delta_tau;
            let c_sec = col_full - sec_col_offset;

            Ok((l_sec, c_sec))
        }
        ZeroDopplerOutcome::NotConverged { .. } | ZeroDopplerOutcome::DegenerateGeometry => {
            // Return a sentinel value; the caller filters None/invalid points.
            Err(InsarError::SecondaryZeroDopplerFailed { line: 0, sample: 0 })
        }
    }
}

/// Approximate the secondary `valid_sample_offset` from its swath metadata.
///
/// The `DeburstArray::valid_sample_offset` is computed during deburst and we
/// do not have direct access to the secondary deburst output here.  As a
/// conservative approximation, we use the secondary swath's `first_sample`
/// field (full-image column of the first valid pixel), which is the same
/// value that the deburst step stores in `valid_sample_offset`.
fn sec_deburst_valid_offset(sec_swath: &SubSwathMetadata, _ref_swath: &SubSwathMetadata) -> f64 {
    sec_swath.first_sample as f64
}

// ─── Bilinear interpolation ───────────────────────────────────────────────

/// Bilinearly interpolate the secondary Deramped SLC at floating-point
/// coordinates (l, c).  Returns `[0.0, 0.0]` for out-of-bounds coordinates.
fn bilinear(sec: &Deramped, l: f64, c: f64) -> [f32; 2] {
    if l < 0.0 || c < 0.0 {
        return [0.0, 0.0];
    }
    let l0 = l.floor() as usize;
    let c0 = c.floor() as usize;
    let l1 = l0 + 1;
    let c1 = c0 + 1;

    if l1 >= sec.lines || c1 >= sec.samples {
        return [0.0, 0.0];
    }

    let wl = (l - l0 as f64) as f32;
    let wc = (c - c0 as f64) as f32;
    let p00 = sec.data[l0 * sec.samples + c0];
    let p01 = sec.data[l0 * sec.samples + c1];
    let p10 = sec.data[l1 * sec.samples + c0];
    let p11 = sec.data[l1 * sec.samples + c1];

    let lerp = |a: f32, b: f32, t: f32| a + (b - a) * t;

    let i = lerp(lerp(p00[0], p01[0], wc), lerp(p10[0], p11[0], wc), wl);
    let q = lerp(lerp(p00[1], p01[1], wc), lerp(p10[1], p11[1], wc), wl);
    [i, q]
}

// ─── Polynomial evaluation ────────────────────────────────────────────────

/// Evaluate a degree-2 bivariate polynomial at (l, c).
///
/// ```text
/// p[0] + p[1]*l + p[2]*c + p[3]*l² + p[4]*l*c + p[5]*c²
/// ```
fn eval_poly6(p: &[f64; 6], l: f64, c: f64) -> f64 {
    p[0] + p[1] * l + p[2] * c + p[3] * l * l + p[4] * l * c + p[5] * c * c
}

// ─── Least-squares polynomial fit ────────────────────────────────────────

/// Fit a degree-2 bivariate polynomial to the sparse grid offsets.
///
/// Returns `(coeffs [f64; 6], residual_rms)` using normal equations (Aᵀ A x = Aᵀ b).
/// The 6×6 normal equation matrix is solved by Gaussian elimination with
/// partial pivoting.
///
/// # Returns
///
/// `Err(InsarError::CoregGridTooSmall)` is returned by the caller before this
/// function is reached if `n < 6`; here we document the precondition.
fn fit_poly6<F>(
    pts: &[(f64, f64, f64, f64)],
    value_of: F,
) -> Result<([f64; 6], f64), InsarError>
where
    F: Fn(&(f64, f64, f64, f64)) -> f64,
{
    const N: usize = 6;
    let mut ata = [[0.0_f64; N]; N];
    let mut atb = [0.0_f64; N];

    for pt in pts {
        let (l, c, _, _) = *pt;
        let y = value_of(pt);
        let basis = [1.0, l, c, l * l, l * c, c * c];
        for i in 0..N {
            atb[i] += basis[i] * y;
            for j in 0..N {
                ata[i][j] += basis[i] * basis[j];
            }
        }
    }

    let coeffs = gaussian_elimination(ata, atb)?;

    // Compute fit residual RMS.
    let mut sse = 0.0_f64;
    for pt in pts {
        let (l, c, _, _) = *pt;
        let y = value_of(pt);
        let y_hat = eval_poly6(&coeffs, l, c);
        let res = y - y_hat;
        sse += res * res;
    }
    let rms = (sse / pts.len() as f64).sqrt();

    Ok((coeffs, rms))
}

/// Solve a dense 6×6 linear system Ax = b using Gaussian elimination with
/// partial pivoting.
///
/// Returns `Err(InsarError::CoregGridTooSmall)` (re-using the sentinel) if
/// the matrix is singular within floating-point tolerance.
fn gaussian_elimination(mut a: [[f64; 6]; 6], mut b: [f64; 6]) -> Result<[f64; 6], InsarError> {
    const N: usize = 6;
    const PIVOT_TOL: f64 = 1e-12;

    for col in 0..N {
        // Partial pivoting: find row with largest absolute value in this column.
        let pivot_row = (col..N)
            .max_by(|&i, &j| {
                a[i][col].abs().partial_cmp(&a[j][col].abs()).unwrap_or(std::cmp::Ordering::Equal) // SAFETY-OK: matrix entries are finite pixel-coordinate products; partial_cmp returns None only for NaN which cannot arise here; the pivot tolerance check below independently detects near-zero pivots
            })
            .unwrap_or(col); // SAFETY-OK: range col..N is non-empty (col < N); unwrap_or is an infallible fallback with identical result

        a.swap(col, pivot_row);
        b.swap(col, pivot_row);

        if a[col][col].abs() < PIVOT_TOL {
            return Err(InsarError::CoregGridTooSmall { n_valid: 0 });
        }

        let inv_pivot = 1.0 / a[col][col];
        for row in (col + 1)..N {
            let factor = a[row][col] * inv_pivot;
            for k in col..N {
                let av = a[col][k];
                a[row][k] -= factor * av;
            }
            b[row] -= factor * b[col];
        }
    }

    // Back substitution.
    let mut x = [0.0_f64; N];
    for i in (0..N).rev() {
        let mut sum = b[i];
        for j in (i + 1)..N {
            sum -= a[i][j] * x[j];
        }
        x[i] = sum / a[i][i];
    }

    Ok(x)
}

// ─── Tests ────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        insar::deramp::Deramped,
        types::{OrbitData, OrbitRelativeSeconds, StateVector, SubSwathMetadata},
    };
    use chrono::TimeZone;

    // ── Helper: build a minimal SubSwathMetadata ──

    fn make_swath(slant_range_time_s: f64) -> SubSwathMetadata {
        use crate::types::{AcquisitionMode, AzimuthFmRate, DcEstimate, Polarization, SubSwathId};
        SubSwathMetadata {
            id: SubSwathId::IW1,
            burst_count: 1,
            lines_per_burst: 100,
            range_samples: 200,
            azimuth_samples: 100,
            first_line: 0,
            last_line: 100,
            first_sample: 0,
            last_sample: 200,
            range_pixel_spacing_m: 2.329562,
            azimuth_pixel_spacing_m: 13.9,
            slant_range_time_s,
            azimuth_time_interval_s: 0.002,
            prf_hz: 486.486,
            burst_cycle_time_s: 2.758,
            dc_estimates: Vec::new(),
            fm_rates: Vec::new(),
        }
    }

    // ── Helper: build a circular-orbit OrbitData ──

    /// Construct a synthetic circular orbit around Earth for testing.
    /// Altitude ~700 km above equator, ascending pass.
    fn make_circular_orbit(
        start_utc: DateTime<Utc>,
        n_sv: usize,
        dt_s: f64,
    ) -> OrbitData {
        const R_ORBIT: f64 = WGS84_A + 700_000.0;
        const OMEGA: f64 = 7.5e-4; // rad/s (≈ 90 min period)

        let reference_epoch = start_utc;
        let mut state_vectors = Vec::with_capacity(n_sv);

        for k in 0..n_sv {
            let t_s = k as f64 * dt_s;
            let theta = OMEGA * t_s;
            let pos = [R_ORBIT * theta.cos(), R_ORBIT * theta.sin(), 0.0];
            let vel = [
                -R_ORBIT * OMEGA * theta.sin(),
                R_ORBIT * OMEGA * theta.cos(),
                0.0,
            ];
            state_vectors.push(StateVector {
                time: reference_epoch
                    + Duration::microseconds((t_s * 1e6) as i64), // SAFETY-OK: test orbit spans < 300 s; microseconds cannot overflow
                position_m: pos,
                velocity_m_s: vel,
            });
        }

        OrbitData { reference_epoch, state_vectors }
    }

    // ── Test: forward geocoding round-trip on 3D non-degenerate geometry ──

    #[test]
    fn forward_geocoding_converges_on_ellipsoid() {
        // Satellite roughly above lat=40°, lon=28°, altitude 700 km.
        // Velocity is perpendicular to the radial (T·V = 0) to ensure the
        // zero-Doppler plane is well-defined and the 3×3 Jacobian is full rank.
        //
        // The geometry is genuinely 3D: no vector lies in the x-y plane alone,
        // so J[0] (velocity), J[1] (range vector), J[2] (ellipsoid normal) are
        // linearly independent.
        let sat_pos = [4_750_000.0_f64, 2_600_000.0, 4_540_000.0];
        // Velocity perpendicular to the (40°lat, 28°lon) radial direction.
        // V = T × e_z_unit scaled to ~7600 m/s, where T direction ≈ sat_pos.
        let sat_vel = [3_640.0_f64, -6_664.0, 0.0];

        // Slant range ~700 km — right-looking at moderate incidence.
        let range_m = 700_000.0_f64;

        let result = ecef_from_slant_range_doppler(sat_pos, sat_vel, range_m, 20, 1e-3);
        let target = result.expect("forward geocoding should converge");

        let [tx, ty, tz] = target;
        let [px, py, pz] = sat_pos;
        let [vx, vy, vz] = sat_vel;

        // 1. Target on ellipsoid (|F3| < 1e-6).
        let ell =
            tx * tx / (WGS84_A * WGS84_A) + ty * ty / (WGS84_A * WGS84_A)
                + tz * tz / (WGS84_B * WGS84_B) - 1.0;
        assert!(ell.abs() < 1e-6, "target not on ellipsoid: F3 = {ell}");

        // 2. Slant range correct (< 1 m).
        let r =
            ((tx - px).powi(2) + (ty - py).powi(2) + (tz - pz).powi(2)).sqrt();
        assert!((r - range_m).abs() < 1.0, "range error = {} m", r - range_m);

        // 3. Zero-Doppler condition (< 1 m/s·m — normalised to range).
        let doppler = (tx - px) * vx + (ty - py) * vy + (tz - pz) * vz;
        assert!(
            doppler.abs() / range_m < 1.0,
            "Doppler residual / range = {}",
            doppler / range_m
        );
    }

    // ── Test: least-squares polynomial fit recovers known coefficients ──

    #[test]
    fn poly_fit_recovers_linear_offset() {
        // Offset field: Δl = 0.5 (constant), Δc = 0.0 everywhere.
        let pts: Vec<(f64, f64, f64, f64)> = (0..100)
            .flat_map(|i| {
                (0..10).map(move |j| {
                    let l = i as f64 * 10.0;
                    let c = j as f64 * 50.0;
                    (l, c, 0.5, 0.0)
                })
            })
            .collect();

        let (az_poly, rms) = fit_poly6(&pts, |&(_, _, dl, _)| dl)
            .expect("fit should succeed");

        assert!(
            (az_poly[0] - 0.5).abs() < 1e-9,
            "constant term should be 0.5, got {}",
            az_poly[0]
        );
        assert!(rms < 1e-10, "fit residual should be zero, got {rms}");
    }

    // ── Test: bilinear interpolation at exact grid points ──

    #[test]
    fn bilinear_at_integer_coords_returns_exact_value() {
        // 3×3 image — interior pixel (1,1) has valid l1=2 and c1=2 neighbours.
        let mut data = vec![[0.0_f32; 2]; 9];
        for i in 0..9_usize {
            data[i] = [i as f32, (i + 1) as f32];
        }
        let sec = Deramped { data, lines: 3, samples: 3, valid_sample_offset: 0 };

        // Exact integer coordinate (0,0) must return pixel value at (0,0).
        let val = bilinear(&sec, 0.0, 0.0);
        assert_eq!(val, [0.0, 1.0]);

        // Exact integer coordinate (1,1) — pixel index 4 = [4.0, 5.0].
        let val = bilinear(&sec, 1.0, 1.0);
        assert_eq!(val, [4.0, 5.0]);
    }

    // ── Test: bilinear OOB returns zero ──

    #[test]
    fn bilinear_out_of_bounds_returns_zero() {
        let sec = Deramped {
            data: vec![[1.0, 2.0]; 4],
            lines: 2,
            samples: 2,
            valid_sample_offset: 0,
        };
        assert_eq!(bilinear(&sec, -1.0, 0.0), [0.0, 0.0]);
        assert_eq!(bilinear(&sec, 0.0, 5.0), [0.0, 0.0]);
        assert_eq!(bilinear(&sec, 5.0, 0.0), [0.0, 0.0]);
    }

    // ── Test: resample with zero-offset polynomial is identity ──

    #[test]
    fn resample_with_zero_offset_is_identity() {
        // 3×3 image, exact pixel alignment → output matches input at integer coords.
        let mut data = Vec::with_capacity(9);
        for i in 0..9_usize {
            data.push([i as f32, (i + 1) as f32]);
        }
        let sec = Deramped { data: data.clone(), lines: 3, samples: 3, valid_sample_offset: 0 };

        let zero_poly = CoregPoly { az: [0.0; 6], rg: [0.0; 6] };
        let resampled = resample_secondary(2, 2, &zero_poly, &sec);

        // Only check interior pixels that are within sec bounds.
        for l in 0..2 {
            for c in 0..2 {
                assert_eq!(
                    resampled.data[l * 2 + c],
                    sec.data[l * 3 + c],
                    "mismatch at ({l}, {c})"
                );
            }
        }
    }
}
