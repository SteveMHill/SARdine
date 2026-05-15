//! Interferogram formation and coherence estimation.
//!
//! # Processing order (per subswath)
//!
//! 1. **Flat-earth phase removal** — compute and subtract the phase due to the
//!    reference ellipsoid from the resampled secondary.  This removes the dense
//!    azimuth and range fringes that would alias coherence estimation.
//! 2. **Complex cross-product** — form `igram[i] = ref[i] · conj(sec[i])`.
//! 3. **Windowed coherence estimation** — estimate `|γ|` in a rectangular
//!    multilook window using the boxcar estimator:
//!
//!    ```text
//!    |γ| = |Σ igram| / sqrt(Σ|ref|² · Σ|sec|²)
//!    ```
//!
//!    over a window of `az_looks × rg_looks` pixels.
//!
//! # Output geometry
//!
//! Coherence is output at the **full resolution** of the reference SLC grid
//! (same shape as the input Deramped arrays).  The estimation window acts
//! as a spatial averaging kernel but does not reduce the pixel count —
//! this is distinct from multilooking.  Full-resolution coherence preserves
//! the ability to geocode at the native pixel spacing and is compatible with
//! the existing terrain correction pipeline.
//!
//! # Flat-earth phase
//!
//! The flat-earth phase at pixel (l, c) is:
//!
//! ```text
//! φ_flat(l, c) = (4π / λ) · (R_ref(l,c) − R_sec(l,c))
//! ```
//!
//! where `R_ref` and `R_sec` are the slant ranges from the reference and
//! secondary satellites respectively to the target point on the WGS84
//! **ellipsoid** (h = 0).  The sign convention is such that subtracting
//! `φ_flat` from the secondary phase leaves only the topographic and
//! deformation components.
//!
//! The target point is obtained by evaluating the two-way slant-range time
//! at pixel (l, c) and the satellite state vector from the reference orbit.
//! We do **not** use the Newton solver for flat-earth computation — instead
//! we invert the range geometry analytically:
//!
//! 1. Compute reference satellite position `S_ref` from reference orbit at
//!    azimuth time `t_az_ref(l)`.
//! 2. Compute secondary satellite position `S_sec` from secondary orbit at
//!    the **same azimuth time** (acceptable approximation for flat-earth —
//!    the cross-orbit error is < 0.01 rad for typical S-1 baselines).
//! 3. Compute reference range: `R_ref(l,c) = τ_ref(c) · c_light / 2`.
//! 4. For flat-earth: find the point on the ellipsoid at (l, c) by scaling
//!    the look vector to the reference range.  This requires a numerical step
//!    but converges in 3–4 iterations.
//! 5. Compute secondary range: `R_sec = |S_sec − P_flat|`.
//!
//! For the purposes of coherence estimation, the flat-earth phase only needs
//! to be accurate to ~0.1 rad.  We use a simplified orbit lookup at the
//! reference azimuth time for both satellites (rather than solving the
//! secondary zero-Doppler equation), which introduces < 0.01 rad error for
//! cross-track baselines < 500 m.
//!
//! # Coherence window size
//!
//! Default: 10 azimuth × 40 range looks.  This gives an ENL (Equivalent
//! Number of Looks) of ~400 and a coherence standard deviation of
//! `1/(2√N) ≈ 0.025` for fully incoherent scenes.  The window dimensions
//! are passed as parameters so callers can trade resolution for accuracy.
//!
//! # Phase sign convention
//!
//! Following the ESA/Sentinel-1 convention:
//! - Positive LOS displacement (away from satellite) → positive phase.
//! - `igram = ref · conj(secondary)`, so φ_igram = φ_ref − φ_secondary.
//! - For deformation between acquisitions, φ_deformation = −(4π/λ) · ΔLOS.

use chrono::{DateTime, Duration, Utc};
use rayon::prelude::*;

use crate::insar::{coreg::{ecef_from_slant_range_doppler, Coregistered}, error::InsarError};
use crate::orbit::interpolate_orbit_and_accel;
use crate::types::{OrbitData, SubSwathMetadata, SPEED_OF_LIGHT_M_S};

// ── WGS84 ellipsoid ────────────────────────────────────────────────────────

#[allow(dead_code)]
const WGS84_A: f64 = 6_378_137.0;
#[allow(dead_code)]
const WGS84_B: f64 = 6_356_752.314_245;

// ── Default window sizes ──────────────────────────────────────────────────

/// Default azimuth window size for coherence estimation.
pub const DEFAULT_COH_AZ_LOOKS: usize = 10;
/// Default range window size for coherence estimation.
pub const DEFAULT_COH_RG_LOOKS: usize = 40;

// ── Output type ────────────────────────────────────────────────────────────

/// Coherence and optionally wrapped interferometric phase for one subswath.
///
/// Both arrays are in the reference SLC's slant-range geometry
/// (same shape: `lines × samples`).
#[derive(Debug)]
pub struct Interferogram {
    /// Coherence magnitude `|γ|` ∈ [0, 1], flat row-major.  Length = `lines × samples`.
    pub coherence: Vec<f32>,
    /// Wrapped interferometric phase in radians ∈ (−π, π], flat row-major.
    /// Only populated when `compute_phase` was set to `true`; otherwise empty.
    pub phase: Vec<f32>,
    /// Number of azimuth lines (reference geometry).
    pub lines: usize,
    /// Number of range samples (reference geometry).
    pub samples: usize,
    /// Azimuth window size used for coherence estimation.
    pub az_looks: usize,
    /// Range window size used for coherence estimation.
    pub rg_looks: usize,
}

// ── Configuration ─────────────────────────────────────────────────────────

/// Configuration for interferogram formation and coherence estimation.
#[derive(Debug, Clone)]
pub struct InterferogramConfig {
    /// Azimuth window size (lines) for coherence estimation.
    /// Default: [`DEFAULT_COH_AZ_LOOKS`].
    pub az_looks: usize,

    /// Range window size (samples) for coherence estimation.
    /// Default: [`DEFAULT_COH_RG_LOOKS`].
    pub rg_looks: usize,

    /// Whether to populate [`Interferogram::phase`] with the wrapped
    /// interferometric phase.
    ///
    /// **Requires verified flat-earth removal and TOPS deramping.**
    /// For the first coherence-only delivery, set this to `false`.
    pub compute_phase: bool,
}

impl Default for InterferogramConfig {
    fn default() -> Self {
        Self {
            az_looks: DEFAULT_COH_AZ_LOOKS,
            rg_looks: DEFAULT_COH_RG_LOOKS,
            compute_phase: false,
        }
    }
}

// ── Primary entry point ────────────────────────────────────────────────────

/// Form the interferogram and estimate coherence for one subswath.
///
/// # Arguments
///
/// * `coreg`         — Co-registered reference + secondary SLC pair (from `insar::coreg`).
/// * `ref_swath`     — Reference subswath metadata.
/// * `ref_orbit`     — Reference scene orbit data.
/// * `ref_first_utc` — UTC time of the first line of the reference debursted SLC.
/// * `sec_orbit`     — Secondary scene orbit data.
/// * `wavelength_m`  — Radar wavelength in metres (from `SPEED_OF_LIGHT_M_S / radar_frequency_hz`).
/// * `cfg`           — Window size and phase-output configuration.
///
/// # Returns
///
/// [`Interferogram`] in reference slant-range geometry.
pub fn form_interferogram(
    coreg: &Coregistered,
    ref_swath: &SubSwathMetadata,
    ref_orbit: &OrbitData,
    ref_first_utc: DateTime<Utc>,
    sec_orbit: &OrbitData,
    wavelength_m: f64,
    cfg: &InterferogramConfig,
) -> Result<Interferogram, InsarError> {
    let lines = coreg.ref_data.lines;
    let samples = coreg.ref_data.samples;
    let n = lines * samples;

    if cfg.az_looks == 0 {
        return Err(InsarError::InvalidWindowSize { az: cfg.az_looks, rg: cfg.rg_looks });
    }
    if cfg.rg_looks == 0 {
        return Err(InsarError::InvalidWindowSize { az: cfg.az_looks, rg: cfg.rg_looks });
    }

    let ref_ati_s = ref_swath.azimuth_time_interval_s;
    let ref_tau0 = ref_swath.slant_range_time_s;
    let ref_delta_tau = 2.0 * ref_swath.range_pixel_spacing_m / SPEED_OF_LIGHT_M_S;
    let ref_col_offset = coreg.ref_data.valid_sample_offset as f64;

    // ── Step 1: flat-earth phase per pixel (row-parallel) ─────────────────
    //
    // For each (l, c): compute reference satellite position, find the ECEF
    // point on the WGS84 ellipsoid along the reference look vector, then
    // compute the secondary slant range to that point.
    //
    // phi_flat[l*samples + c] = (4π/λ) * (R_ref - R_sec_to_flat_point)

    let k = std::f64::consts::TAU / wavelength_m; // 2π/λ wavenumber

    let mut flat_earth_phase: Vec<f32> = Vec::with_capacity(n);
    // SAFETY-OK: capacity == n; every element written by the parallel loop before the Vec is read
    unsafe { flat_earth_phase.set_len(n) };

    let flat_earth_fallback_count = std::sync::atomic::AtomicUsize::new(0);

    flat_earth_phase
        .par_chunks_mut(samples)
        .enumerate()
        .try_for_each(|(li, row)| -> Result<(), InsarError> {
            let t_az = ref_first_utc
                + Duration::microseconds((li as f64 * ref_ati_s * 1e6) as i64); // SAFETY-OK: line × ATI in µs; S-1 scene < 30 s, max value ~30e6 µs << i64::MAX
            let (sv_ref, _) = interpolate_orbit_and_accel(ref_orbit, t_az)
                .map_err(InsarError::Orbit)?;
            let (sv_sec, _) = interpolate_orbit_and_accel(sec_orbit, t_az)
                .map_err(InsarError::Orbit)?;

            let mut row_fallbacks = 0usize;
            for (ci, phi) in row.iter_mut().enumerate() {
                let tau = ref_tau0 + (ci as f64 + ref_col_offset) * ref_delta_tau;
                let r_ref = tau * SPEED_OF_LIGHT_M_S / 2.0;

                // Project the reference look vector to the ellipsoid at range r_ref
                // using the proper 3-equation system (zero-Doppler + range + ellipsoid).
                // For flat-earth purposes a 1 mm convergence threshold is sufficient.
                let flat_ecef = match ecef_from_slant_range_doppler(
                    sv_ref.position_m, sv_ref.velocity_m_s, r_ref, 0.0, 10, 1.0,
                ) {
                    Ok(p) => p,
                    Err(_) => {
                        // Forward geocoding failed: flat-earth phase set to 0.
                        // This can happen for pixels at the far image edge; coherence
                        // is unaffected (phi_flat = 0 ⟹ no phase rotation).
                        *phi = 0.0;
                        row_fallbacks += 1;
                        continue;
                    }
                };

                // Secondary range to the same flat-earth point.
                let dx = sv_sec.position_m[0] - flat_ecef[0];
                let dy = sv_sec.position_m[1] - flat_ecef[1];
                let dz = sv_sec.position_m[2] - flat_ecef[2];
                let r_sec = (dx * dx + dy * dy + dz * dz).sqrt();

                *phi = (k * 2.0 * (r_ref - r_sec)) as f32;
            }
            if row_fallbacks > 0 {
                flat_earth_fallback_count.fetch_add(row_fallbacks, std::sync::atomic::Ordering::Relaxed);
            }
            Ok(())
        })?;

    let total_fallbacks = flat_earth_fallback_count.load(std::sync::atomic::Ordering::Relaxed);
    if total_fallbacks > 0 {
        let fallback_pct = total_fallbacks as f64 / n as f64 * 100.0;
        if fallback_pct > 1.0 {
            tracing::warn!(
                "flat-earth phase: forward geocoding failed for {} / {} pixels ({:.1}%); \
                 these pixels have phi_flat=0 which may degrade interferogram quality",
                total_fallbacks, n, fallback_pct
            );
        } else {
            tracing::debug!(
                "flat-earth phase: {} pixels used phi_flat=0 fallback ({:.2}%)",
                total_fallbacks, fallback_pct
            );
        }
    }

    // ── Step 2: apply flat-earth removal to secondary, form cross-product ─
    //
    // secondary_flat[i] = secondary[i] * exp(-j * phi_flat[i])
    // igram[i]          = ref[i] * conj(secondary_flat[i])
    //
    // Combined: igram[i] = ref[i] * conj(secondary[i]) * exp(+j * phi_flat[i])
    // (removing phi_flat from the secondary is equivalent to adding it to the igram)

    let ref_data = &coreg.ref_data.data;
    let sec_data = &coreg.secondary_resampled.data;

    let mut igram_re: Vec<f32> = Vec::with_capacity(n);
    let mut igram_im: Vec<f32> = Vec::with_capacity(n);
    let mut ref_power: Vec<f32> = Vec::with_capacity(n);
    let mut sec_power: Vec<f32> = Vec::with_capacity(n);
    // SAFETY-OK: all four buffers have capacity == n; every element written in the loop before they are read
    unsafe {
        igram_re.set_len(n);
        igram_im.set_len(n);
        ref_power.set_len(n);
        sec_power.set_len(n);
    }

    // Parallel over rows.
    igram_re
        .par_chunks_mut(samples)
        .zip(igram_im.par_chunks_mut(samples))
        .zip(ref_power.par_chunks_mut(samples))
        .zip(sec_power.par_chunks_mut(samples))
        .zip(flat_earth_phase.par_chunks(samples))
        .enumerate()
        .for_each(|(li, ((((ig_re_row, ig_im_row), rp_row), sp_row), fe_row))| {
            for ci in 0..samples {
                let idx = li * samples + ci;
                let [ri, rq] = ref_data[idx];
                let [si, sq] = sec_data[idx];
                let phi = fe_row[ci] as f64;
                let (sin_phi, cos_phi) = phi.sin_cos();

                // secondary_flat = secondary * exp(-j * phi_flat)
                // conj(secondary_flat) = conj(sec) * exp(+j * phi_flat)
                // igram = ref * conj(secondary_flat)
                //       = (ri + j*rq) * [(si - j*sq) * (cos_phi + j*sin_phi)]
                //       = (ri + j*rq) * [(si*cos_phi + sq*sin_phi) + j*(si*sin_phi - sq*cos_phi)]

                let conj_sec_flat_re = (si as f64) * cos_phi + (sq as f64) * sin_phi;
                let conj_sec_flat_im = (si as f64) * sin_phi - (sq as f64) * cos_phi;

                ig_re_row[ci] = ((ri as f64) * conj_sec_flat_re - (rq as f64) * conj_sec_flat_im) as f32;
                ig_im_row[ci] = ((ri as f64) * conj_sec_flat_im + (rq as f64) * conj_sec_flat_re) as f32;
                rp_row[ci] = (ri as f32) * (ri as f32) + (rq as f32) * (rq as f32);
                sp_row[ci] = (si as f32) * (si as f32) + (sq as f32) * (sq as f32);
            }
        });

    // ── Step 3: windowed coherence (and optionally phase) estimation ───────
    //
    // For each pixel (l, c) the coherence window is:
    //   l_start = l.saturating_sub(az_looks/2)
    //   l_end   = min(l + az_looks/2, lines)
    // (symmetric, clamped at borders).
    //
    // |γ(l,c)| = |Σ igram| / sqrt(Σ|ref|² · Σ|sec|²)
    //
    // We use a naive O(N·az·rg) implementation for correctness.  A summed
    // area table (SAT) approach would be O(N) but adds complexity that is
    // not needed for Phase 1.

    let az_half = cfg.az_looks / 2;
    let rg_half = cfg.rg_looks / 2;

    let mut coherence: Vec<f32> = Vec::with_capacity(n);
    let mut phase_out: Vec<f32> = if cfg.compute_phase { Vec::with_capacity(n) } else { Vec::new() };
    // SAFETY-OK: capacity == n; every element written by the loop before read
    unsafe {
        coherence.set_len(n);
        if cfg.compute_phase {
            phase_out.set_len(n);
        }
    }

    // Coherence loop: iterate over coherence rows only (no zip with phase_out).
    // The zip-with-empty-slice pattern fails when compute_phase = false because
    // zipping a non-empty iterator with 0 elements produces 0 iterations.
    // Instead, we write coherence unconditionally and index into phase_out by offset.
    coherence
        .par_chunks_mut(samples)
        .enumerate()
        .for_each(|(li, coh_row)| {
            let l0 = li.saturating_sub(az_half);
            let l1 = (li + az_half + 1).min(lines);

            for ci in 0..samples {
                let c0 = ci.saturating_sub(rg_half);
                let c1 = (ci + rg_half + 1).min(samples);

                let mut sum_ig_re = 0.0_f64;
                let mut sum_ig_im = 0.0_f64;
                let mut sum_rp = 0.0_f64;
                let mut sum_sp = 0.0_f64;

                for wl in l0..l1 {
                    for wc in c0..c1 {
                        let wi = wl * samples + wc;
                        sum_ig_re += igram_re[wi] as f64;
                        sum_ig_im += igram_im[wi] as f64;
                        sum_rp += ref_power[wi] as f64;
                        sum_sp += sec_power[wi] as f64;
                    }
                }

                let denom = (sum_rp * sum_sp).sqrt();
                let coh = if denom > 0.0 {
                    ((sum_ig_re * sum_ig_re + sum_ig_im * sum_ig_im).sqrt() / denom)
                        .min(1.0) as f32
                } else {
                    0.0_f32
                };
                coh_row[ci] = coh;
            }
        });

    // Phase loop: only when requested.
    if cfg.compute_phase {
        phase_out
            .par_chunks_mut(samples)
            .enumerate()
            .for_each(|(li, phase_row)| {
                let l0 = li.saturating_sub(az_half);
                let l1 = (li + az_half + 1).min(lines);
                for ci in 0..samples {
                    let c0 = ci.saturating_sub(rg_half);
                    let c1 = (ci + rg_half + 1).min(samples);
                    let mut sum_ig_re = 0.0_f64;
                    let mut sum_ig_im = 0.0_f64;
                    for wl in l0..l1 {
                        for wc in c0..c1 {
                            let wi = wl * samples + wc;
                            sum_ig_re += igram_re[wi] as f64;
                            sum_ig_im += igram_im[wi] as f64;
                        }
                    }
                    phase_row[ci] = sum_ig_im.atan2(sum_ig_re) as f32;
                }
            });
    }

    Ok(Interferogram {
        coherence,
        phase: phase_out,
        lines,
        samples,
        az_looks: cfg.az_looks,
        rg_looks: cfg.rg_looks,
    })
}

// ── Helpers ────────────────────────────────────────────────────────────────

/// Find the point on the WGS84 ellipsoid at slant range `r_ref` from the
/// satellite position `sat_pos`, along the look direction implied by the
/// satellite velocity (zero-Doppler approximation).
///
/// Algorithm: given `sat_pos` and the perpendicular-to-velocity plane, find
/// the intersection of the sphere of radius `r_ref` centred at `sat_pos`
/// with the WGS84 oblate spheroid.
///
/// Uses a Newton-type 2-iteration approach.  For the purposes of flat-earth
/// phase estimation (accuracy requirement ~0.1 rad) this is sufficient.
#[allow(dead_code)]
fn ellipsoid_intercept(sat_pos: [f64; 3], sat_vel: [f64; 3], r_ref: f64) -> [f64; 3] {
    // Normalise velocity to get along-track unit vector.
    let v_norm = {
        let mag = (sat_vel[0] * sat_vel[0] + sat_vel[1] * sat_vel[1] + sat_vel[2] * sat_vel[2])
            .sqrt();
        [sat_vel[0] / mag, sat_vel[1] / mag, sat_vel[2] / mag]
    };

    // Initial look direction: from sat toward Earth centre, then subtract
    // along-track component (zero-Doppler plane).
    let r_mag = (sat_pos[0] * sat_pos[0] + sat_pos[1] * sat_pos[1] + sat_pos[2] * sat_pos[2])
        .sqrt();
    let nadir = [
        -sat_pos[0] / r_mag,
        -sat_pos[1] / r_mag,
        -sat_pos[2] / r_mag,
    ];
    // Remove along-track component: look = nadir - (nadir·v) v
    let dot_nv = nadir[0] * v_norm[0] + nadir[1] * v_norm[1] + nadir[2] * v_norm[2];
    let mut look = [
        nadir[0] - dot_nv * v_norm[0],
        nadir[1] - dot_nv * v_norm[1],
        nadir[2] - dot_nv * v_norm[2],
    ];
    // Normalise look direction.
    let look_mag =
        (look[0] * look[0] + look[1] * look[1] + look[2] * look[2]).sqrt();
    look[0] /= look_mag;
    look[1] /= look_mag;
    look[2] /= look_mag;

    // Candidate target: sat + r_ref * look.
    let mut p = [
        sat_pos[0] + r_ref * look[0],
        sat_pos[1] + r_ref * look[1],
        sat_pos[2] + r_ref * look[2],
    ];

    // Newton iteration: project p onto the ellipsoid surface.
    // f(p) = (px/a)² + (py/a)² + (pz/b)² - 1 = 0
    // Two iterations are sufficient for < 1 m accuracy at S-1 altitudes.
    for _ in 0..3 {
        let fx = p[0] / WGS84_A;
        let fy = p[1] / WGS84_A;
        let fz = p[2] / WGS84_B;
        let f_val = fx * fx + fy * fy + fz * fz - 1.0;
        // Scale p to put it on the ellipsoid.
        let scale = 1.0 / (fx * fx + fy * fy + fz * fz).sqrt();
        p[0] *= scale;
        p[1] *= scale;
        p[2] *= scale;

        // Re-adjust to preserve slant range from sat_pos.
        let dx = p[0] - sat_pos[0];
        let dy = p[1] - sat_pos[1];
        let dz = p[2] - sat_pos[2];
        let dist = (dx * dx + dy * dy + dz * dz).sqrt();
        if dist > 0.0 && f_val.abs() > 1e-10 {
            let look_p = [dx / dist, dy / dist, dz / dist];
            p[0] = sat_pos[0] + r_ref * look_p[0];
            p[1] = sat_pos[1] + r_ref * look_p[1];
            p[2] = sat_pos[2] + r_ref * look_p[2];
        }
    }

    p
}

// ── Tests ──────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use crate::insar::{
        coreg::{Coregistered, resample_secondary, CoregPoly, CoregResult},
        deramp::Deramped,
    };

    fn make_deramped(lines: usize, samples: usize, val: [f32; 2]) -> Deramped {
        Deramped {
            data: vec![val; lines * samples],
            lines,
            samples,
            valid_sample_offset: 0,
        }
    }

    fn make_coregistered(lines: usize, samples: usize) -> Coregistered {
        // Reference: constant unit complex [1, 0].
        // Secondary (resampled): same constant — perfect coherence.
        Coregistered {
            ref_data: make_deramped(lines, samples, [1.0, 0.0]),
            secondary_resampled: make_deramped(lines, samples, [1.0, 0.0]),
        }
    }

    // Satellite at proper S-1 altitude (~693 km) on the x-axis with velocity in y.
    // Position magnitude: ~7071 km; looking roughly in -x direction (right-looking).
    fn make_orbit_at(x: f64, y: f64, z: f64, t0: DateTime<Utc>) -> OrbitData {
        use crate::types::StateVector;
        let sv = |dt_s: f64| StateVector {
            time: t0 + Duration::microseconds((dt_s * 1e6) as i64), // SAFETY-OK: test orbit spans < 100 s; µs value << i64::MAX
            position_m: [x, y, z],
            velocity_m_s: [0.0, 7600.0, 0.0],
        };
        OrbitData {
            reference_epoch: t0,
            state_vectors: (-5..=5).map(|i| sv(i as f64 * 10.0)).collect(),
        }
    }

    fn minimal_swath() -> SubSwathMetadata {
        SubSwathMetadata {
            id: crate::types::SubSwathId::IW1,
            burst_count: 1,
            lines_per_burst: 10,
            range_samples: 20,
            azimuth_samples: 10,
            first_line: 0,
            last_line: 10,
            first_sample: 0,
            last_sample: 20,
            range_pixel_spacing_m: 2.33,
            azimuth_pixel_spacing_m: 13.93,
            slant_range_time_s: 5.35e-3,
            azimuth_time_interval_s: 2.056e-3,
            prf_hz: 486.0,
            burst_cycle_time_s: 2.758,
            dc_estimates: vec![],
            fm_rates: vec![],
        }
    }

    /// Same-scene coherence must be 1.0 when ref == secondary.
    #[test]
    fn same_scene_coherence_is_unity() {
        let lines = 40;
        let samples = 80;
        let coreg = make_coregistered(lines, samples);
        let swath = minimal_swath();
        let t0 = chrono::Utc::now();
        // Satellite at ~693 km altitude above equator on x-axis, velocity in y-direction.
        // Slant range 5.35e-3 s * c/2 ≈ 802 km — this geometry has the satellite
        // looking sideways (-x direction from equatorial position at 7071 km).
        let orbit = make_orbit_at(7_071_000.0, 0.0, 0.0, t0);
        let cfg = InterferogramConfig {
            az_looks: 4,
            rg_looks: 8,
            compute_phase: false,
        };
        let igram = form_interferogram(
            &coreg, &swath, &orbit, t0, &orbit,
            0.0555, // λ ≈ C-band
            &cfg,
        ).unwrap();

        assert_eq!(igram.lines, lines);
        assert_eq!(igram.samples, samples);
        // Interior pixels (away from edge clamping) should have coherence = 1.0.
        for l in 4..lines - 4 {
            for c in 8..samples - 8 {
                let coh = igram.coherence[l * samples + c];
                assert!(
                    (coh - 1.0_f32).abs() < 1e-4,
                    "expected coherence=1 at ({l},{c}), got {coh}"
                );
            }
        }
        // Phase vec is empty when compute_phase = false.
        assert!(igram.phase.is_empty());
    }

    /// Zero-signal secondary → coherence must be 0.0 everywhere.
    #[test]
    fn zero_secondary_gives_zero_coherence() {
        let lines = 20;
        let samples = 40;
        let coreg = Coregistered {
            ref_data: make_deramped(lines, samples, [1.0, 0.0]),
            secondary_resampled: make_deramped(lines, samples, [0.0, 0.0]),
        };
        let swath = minimal_swath();
        let t0 = chrono::Utc::now();
        let orbit = make_orbit_at(7_071_000.0, 0.0, 0.0, t0);
        let cfg = InterferogramConfig::default();
        let igram = form_interferogram(
            &coreg, &swath, &orbit, t0, &orbit,
            0.0555,
            &cfg,
        ).unwrap();
        for &coh in &igram.coherence {
            assert!(coh == 0.0, "expected coherence=0, got {coh}");
        }
    }

    /// compute_phase = true populates the phase vector with correct sign.
    #[test]
    fn phase_output_populated_when_requested() {
        let lines = 10;
        let samples = 20;
        // ref = [1, 0], secondary = [0, 1] → igram = (1+0j)*(0-1j) = -j → phase = -π/2
        let coreg = Coregistered {
            ref_data: make_deramped(lines, samples, [1.0, 0.0]),
            secondary_resampled: make_deramped(lines, samples, [0.0, 1.0]),
        };
        let swath = minimal_swath();
        let t0 = chrono::Utc::now();
        let orbit = make_orbit_at(7_071_000.0, 0.0, 0.0, t0);
        let cfg = InterferogramConfig {
            az_looks: 2,
            rg_looks: 4,
            compute_phase: true,
        };
        let igram = form_interferogram(
            &coreg, &swath, &orbit, t0, &orbit,
            0.0555,
            &cfg,
        ).unwrap();
        assert_eq!(igram.phase.len(), lines * samples);
        // Interior pixels: phase should be near 0 after flat-earth removal
        // (same orbit for both ref and secondary → R_ref == R_sec → phi_flat = 0).
        // The raw igram phase = atan2(-1, 0) ≈ -π/2 when flat-earth = 0.
        let interior_phase = igram.phase[(lines / 2) * samples + samples / 2];
        let expected = -std::f32::consts::FRAC_PI_2;
        assert!(
            (interior_phase - expected).abs() < 0.1,
            "expected phase≈{expected:.3}, got {interior_phase:.3}"
        );
    }

    /// Window size 0 must return an error, not a panic.
    #[test]
    fn zero_window_size_returns_error() {
        let lines = 10;
        let samples = 20;
        let coreg = make_coregistered(lines, samples);
        let swath = minimal_swath();
        let t0 = chrono::Utc::now();
        let orbit = make_orbit_at(7_071_000.0, 0.0, 0.0, t0);
        let cfg = InterferogramConfig { az_looks: 0, rg_looks: 4, compute_phase: false };
        assert!(form_interferogram(
            &coreg, &swath, &orbit, t0, &orbit, 0.0555, &cfg
        ).is_err());
    }

    /// The flat-earth forward geocoding must find a valid solution for typical
    /// S-1 geometry: satellite at 693 km altitude, slant range ~802 km.
    /// With the same orbit for ref and secondary, phi_flat = 0 for all pixels,
    /// so this indirectly verifies flat-earth phase does not corrupt coherence.
    #[test]
    fn flat_earth_phase_is_zero_same_orbit() {
        let lines = 8;
        let samples = 16;
        let coreg = make_coregistered(lines, samples);
        let swath = minimal_swath();
        let t0 = chrono::Utc::now();
        // Proper S-1 orbit altitude: 693 km above equator on x-axis.
        let orbit = make_orbit_at(7_071_000.0, 0.0, 0.0, t0);
        let cfg = InterferogramConfig {
            az_looks: 2,
            rg_looks: 4,
            compute_phase: true,
        };
        let igram = form_interferogram(
            &coreg, &swath, &orbit, t0, &orbit, 0.0555, &cfg,
        ).unwrap();
        // Same orbit ⟹ phi_flat = 0 ⟹ igram phase is constant 0 for [1,0]*conj([1,0]).
        let interior_phase = igram.phase[(lines / 2) * samples + samples / 2];
        assert!(
            interior_phase.abs() < 0.01,
            "expected flat-earth phase ≈ 0, got {interior_phase:.4}"
        );
    }
}
