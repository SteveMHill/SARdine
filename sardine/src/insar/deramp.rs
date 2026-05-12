//! TOPS azimuth deramping for Sentinel-1 IW SLC data.
//!
//! # Why deramping is required
//!
//! In TOPS (Terrain Observation with Progressive Scans) mode the radar beam
//! steers continuously from backward to forward across each burst.  This beam
//! steering imposes a time-varying Doppler centroid shift on every burst:
//!
//! ```text
//! f_DC(l) ≈ 0 at burst centre,   f_DC ≈ ±2300 Hz at burst edges (S1B IW1)
//! ```
//!
//! The instantaneous phase of each complex sample carries this steering phase:
//!
//! ```text
//! φ_TOPS(l, c) = 2π · f_DC(τ_c) · t_az(l)
//! ```
//!
//! where `τ_c` is the two-way slant-range time at column `c` and `t_az(l)` is
//! the azimuth time of line `l` relative to the first line of the debursted
//! array.
//!
//! For InSAR cross-multiplication `ref × conj(secondary)`, the steering phase
//! would appear as a residual ramp unless both scenes are processed with
//! identical Doppler centroid polynomials **and** are perfectly co-registered.
//! A 0.1-line azimuth co-registration error produces up to ±3 rad of phase
//! error at burst edges (see reevaluation notes).
//!
//! Deramping multiplies each sample by `exp(-jφ_TOPS)`, converting the burst
//! to an equivalent "zero-DC" signal.  After forming the interferogram, the
//! reramp step (multiplying by `exp(+jφ_ref_TOPS)`) restores the correct
//! phase for geocoding.
//!
//! **For coherence magnitude only** (Phase 1 deliverable) the ramp has no
//! effect — `|γ|` is phase-invariant.  However, deramping is required for
//! any wrapped-phase output, and it is implemented here so that Phase 4 can
//! gate wrapped-phase output behind a flag without restructuring this step.
//!
//! # DC polynomial evaluation
//!
//! Each [`DcEstimate`] covers approximately one burst cycle (~2.75 s).
//! For any given azimuth line `l`:
//! 1. Find the two nearest estimates bracketing `l` in azimuth time.
//! 2. Linearly interpolate the polynomial coefficients in time between them.
//! 3. Evaluate the interpolated polynomial at the pixel's slant-range time.
//!
//! The polynomial is:
//!
//! ```text
//! f_DC(τ) = a0 + a1·(τ − t0) + a2·(τ − t0)²
//! ```
//!
//! where `τ` is in seconds (slant-range time) and `t0` is per-estimate.
//!
//! # Output
//!
//! [`deramp_subswath`] returns a `Vec<[f32; 2]>` of the same shape as the
//! input `DeburstArray`.  Each element `[i_out, q_out]` is the complex product
//! of the original i16 sample with `exp(-jφ)`.  The `i16 → f32` promotion
//! happens here so downstream operations work in floating point.

use crate::deburst::DeburstArray;
use crate::types::{DcEstimate, SubSwathMetadata};

use super::error::InsarError;

use rayon::prelude::*;

// ── Public output type ─────────────────────────────────────────────────────

/// A deramped single-subswath SLC array in complex f32.
///
/// Same geometry as the input [`DeburstArray`]: `lines × samples` flat
/// row-major buffer with the same `valid_sample_offset`.
#[derive(Debug)]
pub struct Deramped {
    /// Flat row-major `[I, Q]` buffer.  Length = `lines × samples`.
    pub data: Vec<[f32; 2]>,
    /// Number of azimuth lines.
    pub lines: usize,
    /// Number of range samples per line (same as input `DeburstArray`).
    pub samples: usize,
    /// Same as the input `DeburstArray::valid_sample_offset`.
    pub valid_sample_offset: usize,
}

// ── Core function ──────────────────────────────────────────────────────────

/// Deramp a debursted IW SLC subswath.
///
/// Multiplies every complex sample by `exp(−j · 2π · f_DC(τ_c) · t_az(l))`
/// where:
///
/// - `f_DC(τ_c)` is the Doppler centroid (Hz) at the pixel's slant-range time,
///   obtained by temporally interpolating between the two nearest [`DcEstimate`]s.
/// - `t_az(l)` is the azimuth time of line `l` relative to the **first burst's
///   start time** stored in `burst_start_utc`.
///
/// # Arguments
///
/// * `deburst`         — Debursted SLC array (CInt16 pairs, row-major).
/// * `swath`           — Subswath metadata (ATI, slant-range time, pixel spacing).
/// * `burst_start_utc` — UTC azimuth time of the first line in `deburst`.
///                       Typically `bursts[first_burst_index].azimuth_time_utc`.
///
/// # Errors
///
/// - [`InsarError::NoDcEstimates`] if `swath.dc_estimates` is empty.
/// - [`InsarError::AzimuthTimeBeforeDcEstimates`] if a line's azimuth time
///   falls before the first DC estimate (should not occur for a valid product).
pub fn deramp_subswath(
    deburst: &DeburstArray,
    swath: &SubSwathMetadata,
    burst_start_utc: chrono::DateTime<chrono::Utc>,
) -> Result<Deramped, InsarError> {
    let dc_estimates = &swath.dc_estimates;
    if dc_estimates.is_empty() {
        return Err(InsarError::NoDcEstimates {
            swath: format!("{:?}", swath.id),
        });
    }

    let lines = deburst.lines;
    let samples = deburst.samples;
    let n_pixels = lines * samples;

    // Pre-compute the two-way slant-range time for every column.
    // τ_c = near_range_time + c * Δτ_r
    // where Δτ_r = 2 * range_pixel_spacing / c_light
    let delta_tau_r =
        2.0 * swath.range_pixel_spacing_m / crate::types::SPEED_OF_LIGHT_M_S;
    let near_tau = swath.slant_range_time_s;
    // valid_sample_offset shifts col 0 in the DeburstArray to col valid_sample_offset
    // in the full TIFF.  The slant-range time of DeburstArray col c is therefore
    // the time of TIFF col (c + valid_sample_offset):
    let col_offset = deburst.valid_sample_offset;

    let tau_col: Vec<f64> = (0..samples)
        .map(|c| near_tau + (c + col_offset) as f64 * delta_tau_r)
        .collect();

    // Derive the DC polynomial coefficients for each line by temporally
    // interpolating the two nearest estimates.
    //
    // We precompute the interpolated [a0, a1, a2, t0] for every LINE (not
    // per-pixel) because the polynomial only changes in the azimuth direction:
    // all columns in the same line share the same interpolated estimate.
    let ati_s = swath.azimuth_time_interval_s;

    // Convert DC estimate azimuth times to seconds since burst_start_utc for
    // fast arithmetic.
    let est_times_s: Vec<f64> = dc_estimates
        .iter()
        .map(|e| {
            (e.azimuth_time - burst_start_utc)
                .num_microseconds()
                .unwrap_or(0) as f64 // SAFETY-OK: chrono microseconds cannot overflow for intra-scene durations (< 30 s)
                / 1e6
        })
        .collect();

    // Compute the interpolated polynomial coefficients for each line.
    // Returns [a0_interp, a1_interp, a2_interp, t0_interp].
    let line_poly: Vec<[f64; 4]> = (0..lines)
        .map(|l| {
            let t_az = l as f64 * ati_s;
            interp_dc_poly(dc_estimates, &est_times_s, t_az)
        })
        .collect::<Result<Vec<_>, _>>()?;

    // Allocate output buffer (uninitialised; every element will be written).
    let mut out: Vec<[f32; 2]> = Vec::with_capacity(n_pixels);
    // SAFETY-OK: every element of `out` is written by the parallel loop below
    // (one write per input sample, branchless over the full lines × samples grid)
    // before `out` is returned or read outside this function.
    unsafe { out.set_len(n_pixels) };

    // Parallel over rows.  Each row gets its own polynomial [a0, a1, a2, t0].
    out.par_chunks_mut(samples)
        .zip(deburst.data.par_chunks(samples))
        .enumerate()
        .for_each(|(l, (out_row, in_row))| {
            let [a0, a1, a2, t0] = line_poly[l];
            let t_az = l as f64 * ati_s;

            for (c, (out_px, &in_px)) in out_row.iter_mut().zip(in_row.iter()).enumerate() {
                let tau = tau_col[c];
                let dt = tau - t0;
                let f_dc = a0 + a1 * dt + a2 * dt * dt;
                // Phase: φ = −2π · f_dc · t_az
                let phi = std::f64::consts::TAU * f_dc * t_az;
                // Complex multiply: (i + jq) · exp(−jφ) = (i·cos φ + q·sin φ) + j(q·cos φ − i·sin φ)
                // Note sign: exp(−jφ) = cos(−φ) + j·sin(−φ) = cos(φ) − j·sin(φ)
                let (sin_phi, cos_phi) = phi.sin_cos();
                let i_in = in_px[0] as f64;
                let q_in = in_px[1] as f64;
                *out_px = [
                    (i_in * cos_phi + q_in * sin_phi) as f32,
                    (q_in * cos_phi - i_in * sin_phi) as f32,
                ];
            }
        });

    Ok(Deramped {
        data: out,
        lines,
        samples,
        valid_sample_offset: deburst.valid_sample_offset,
    })
}

// ── Interpolation helpers ──────────────────────────────────────────────────

/// Linearly interpolate DC polynomial coefficients to azimuth time `t_az_s`.
///
/// `est_times_s` are the azimuth times of each estimate in seconds since
/// `burst_start_utc` (same origin as `t_az_s`).
///
/// Returns `[a0_interp, a1_interp, a2_interp, t0_interp]`.
///
/// Extrapolation rule:
/// - If `t_az_s` is before the first estimate, the first estimate is used
///   (no extrapolation — returns an error if the gap is large).
/// - If `t_az_s` is after the last estimate, the last estimate is used.
fn interp_dc_poly(
    estimates: &[DcEstimate],
    est_times_s: &[f64],
    t_az_s: f64,
) -> Result<[f64; 4], InsarError> {
    debug_assert_eq!(estimates.len(), est_times_s.len());

    // Find the bracketing pair (k, k+1) such that est_times_s[k] ≤ t_az_s < est_times_s[k+1].
    let n = estimates.len();

    if t_az_s < est_times_s[0] {
        // Extrapolation before first estimate.
        // Allow up to half a burst cycle (~1.4 s) before the first estimate
        // without error — this covers the partial first burst.
        let gap = est_times_s[0] - t_az_s;
        if gap > 2.0 {
            return Err(InsarError::AzimuthTimeBeforeDcEstimates {
                az_time_s: t_az_s,
                first_est_time_s: est_times_s[0],
            });
        }
        let e = &estimates[0];
        return Ok([e.data_dc_poly[0], e.data_dc_poly[1], e.data_dc_poly[2], e.t0_s]);
    }

    if t_az_s >= est_times_s[n - 1] {
        // Extrapolation past last estimate: use the last one.
        let e = &estimates[n - 1];
        return Ok([e.data_dc_poly[0], e.data_dc_poly[1], e.data_dc_poly[2], e.t0_s]);
    }

    // Binary search for the lower bracket index.
    let k = est_times_s
        .windows(2)
        .position(|w| t_az_s >= w[0] && t_az_s < w[1])
        .unwrap_or(n - 2); // SAFETY-OK: the early-return guards above ensure t_az_s is in [est_times_s[0], est_times_s[n-1]); a None result is therefore unreachable

    let t0 = est_times_s[k];
    let t1 = est_times_s[k + 1];
    let alpha = (t_az_s - t0) / (t1 - t0);

    let e0 = &estimates[k];
    let e1 = &estimates[k + 1];

    // Linearly interpolate each polynomial coefficient.
    let a0 = lerp(e0.data_dc_poly[0], e1.data_dc_poly[0], alpha);
    let a1 = lerp(e0.data_dc_poly[1], e1.data_dc_poly[1], alpha);
    let a2 = lerp(e0.data_dc_poly[2], e1.data_dc_poly[2], alpha);
    // Interpolate t0 as well so the polynomial evaluation point tracks the estimates.
    let t0_interp = lerp(e0.t0_s, e1.t0_s, alpha);

    Ok([a0, a1, a2, t0_interp])
}

#[inline]
fn lerp(a: f64, b: f64, t: f64) -> f64 {
    a + (b - a) * t
}

// ── Reramp ────────────────────────────────────────────────────────────────

/// Reramp a full-resolution interferometric phase array.
///
/// Adds the reference TOPS steering phase `+2π · f_DC(τ_c) · t_az(l)` to
/// every element of `phase` and wraps the result to (−π, π].
///
/// This is the inverse of the deramping applied to the reference SLC in
/// [`deramp_subswath`].  It must be called **after** coherence estimation
/// but **before** terrain-correction geocoding so that the geocoded wrapped
/// phase is free of TOPS-mode steering artefacts.
///
/// # Coherence output
///
/// Reramp is NOT required for coherence magnitude (`|γ|` is phase-invariant).
/// Only call this function when producing a wrapped-phase output.
///
/// # Arguments
///
/// * `phase`               — Wrapped interferometric phase (rad), modified in
///                           place.  Length must equal `lines × samples`.
/// * `lines`               — Number of azimuth lines in the reference SLC
///                           geometry.
/// * `samples`             — Number of range samples per line.
/// * `valid_sample_offset` — First valid range column in the full-TIFF
///                           geometry (same as [`Deramped::valid_sample_offset`]).
/// * `swath`               — Reference subswath metadata (DC estimates + timing).
/// * `burst_start_utc`     — UTC azimuth time of line 0 of `phase`.
///
/// # Errors
///
/// Returns [`InsarError::NoDcEstimates`] if `swath.dc_estimates` is empty.
pub fn reramp_interferogram_phase(
    phase: &mut [f32],
    lines: usize,
    samples: usize,
    valid_sample_offset: usize,
    swath: &SubSwathMetadata,
    burst_start_utc: chrono::DateTime<chrono::Utc>,
) -> Result<(), InsarError> {
    debug_assert_eq!(
        phase.len(),
        lines * samples,
        "phase buffer length {} ≠ lines({}) × samples({})",
        phase.len(),
        lines,
        samples
    );

    let dc_estimates = &swath.dc_estimates;
    if dc_estimates.is_empty() {
        return Err(InsarError::NoDcEstimates {
            swath: format!("{:?}", swath.id),
        });
    }

    let delta_tau_r =
        2.0 * swath.range_pixel_spacing_m / crate::types::SPEED_OF_LIGHT_M_S;
    let near_tau = swath.slant_range_time_s;
    let ati_s = swath.azimuth_time_interval_s;

    // Pre-compute slant-range time per column (identical to deramp_subswath).
    let tau_col: Vec<f64> = (0..samples)
        .map(|c| near_tau + (c + valid_sample_offset) as f64 * delta_tau_r)
        .collect();

    // Convert DC estimate azimuth times to seconds since burst_start_utc.
    let est_times_s: Vec<f64> = dc_estimates
        .iter()
        .map(|e| {
            (e.azimuth_time - burst_start_utc)
                .num_microseconds()
                .unwrap_or(0) as f64 // SAFETY-OK: intra-scene duration < 30 s; µs value << i64::MAX
                / 1e6
        })
        .collect();

    // Pre-compute interpolated DC polynomial per line.
    let line_poly: Vec<[f64; 4]> = (0..lines)
        .map(|l| {
            let t_az = l as f64 * ati_s;
            interp_dc_poly(dc_estimates, &est_times_s, t_az)
        })
        .collect::<Result<Vec<_>, _>>()?;

    // Parallel reramp: φ_out = wrap(φ_in + 2π·f_DC(τ_c)·t_az(l)).
    phase
        .par_chunks_mut(samples)
        .enumerate()
        .for_each(|(l, row)| {
            let [a0, a1, a2, t0] = line_poly[l];
            let t_az = l as f64 * ati_s;
            for (c, phi_out) in row.iter_mut().enumerate() {
                let tau = tau_col[c];
                let dt = tau - t0;
                let f_dc = a0 + a1 * dt + a2 * dt * dt;
                // Reramp: opposite sign to deramp (add back the steering phase).
                let phi_reramp = std::f64::consts::TAU * f_dc * t_az;
                // Wrap to (−π, π].
                let raw = *phi_out as f64 + phi_reramp;
                let wrapped = raw.rem_euclid(std::f64::consts::TAU);
                *phi_out = (if wrapped > std::f64::consts::PI {
                    wrapped - std::f64::consts::TAU
                } else {
                    wrapped
                }) as f32;
            }
        });

    Ok(())
}

// ── Unit tests ─────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use crate::types::{DcEstimate, SubSwathId, SubSwathMetadata};
    use chrono::{TimeZone, Utc};

    fn make_dc_estimate(az_time_s_from_epoch: f64, a0: f64) -> DcEstimate {
        let epoch = Utc.with_ymd_and_hms(2019, 1, 23, 5, 33, 48).unwrap();
        let az_time = epoch
            + chrono::Duration::microseconds((az_time_s_from_epoch * 1e6) as i64);
        DcEstimate {
            azimuth_time: az_time,
            t0_s: 5.354e-3,
            data_dc_poly: [a0, 0.0, 0.0],
        }
    }

    fn make_subswath(dc_estimates: Vec<DcEstimate>) -> SubSwathMetadata {
        SubSwathMetadata {
            id: SubSwathId::IW1,
            burst_count: 1,
            lines_per_burst: 10,
            range_samples: 5,
            azimuth_samples: 10,
            first_line: 0,
            last_line: 10,
            first_sample: 0,
            last_sample: 5,
            range_pixel_spacing_m: 2.329_562,
            azimuth_pixel_spacing_m: 13.93,
            slant_range_time_s: 5.354e-3,
            azimuth_time_interval_s: 2.055_556e-3,
            prf_hz: 1717.0,
            burst_cycle_time_s: 2.758,
            dc_estimates,
            fm_rates: Vec::new(),
        }
    }

    fn make_deburst(lines: usize, samples: usize) -> DeburstArray {
        // Fill with a constant non-zero complex sample [100, 0] (pure real).
        let data = vec![[100i16, 0i16]; lines * samples];
        DeburstArray {
            data,
            lines,
            samples,
            valid_sample_offset: 0,
        }
    }

    /// When f_DC = 0 for all pixels, the deramp is identity (exp(0) = 1).
    #[test]
    fn zero_dc_is_identity() {
        let epoch = Utc.with_ymd_and_hms(2019, 1, 23, 5, 33, 48).unwrap();
        let dc = vec![
            make_dc_estimate(0.0, 0.0),
            make_dc_estimate(3.0, 0.0),
        ];
        let swath = make_subswath(dc);
        let deburst = make_deburst(4, 3);

        let result = deramp_subswath(&deburst, &swath, epoch).unwrap();

        assert_eq!(result.lines, 4);
        assert_eq!(result.samples, 3);

        for &[i_out, q_out] in &result.data {
            // With f_DC = 0 the phase φ = 0, so exp(−jφ) = 1 + 0j.
            // Input [100, 0] → output [100, 0].
            assert!(
                (i_out - 100.0).abs() < 1e-3,
                "i_out = {i_out:.4} — expected 100"
            );
            assert!(
                q_out.abs() < 1e-3,
                "q_out = {q_out:.4} — expected 0"
            );
        }
    }

    /// Empty dc_estimates must return NoDcEstimates error.
    #[test]
    fn empty_dc_estimates_returns_error() {
        let epoch = Utc.with_ymd_and_hms(2019, 1, 23, 5, 33, 48).unwrap();
        let swath = make_subswath(Vec::new());
        let deburst = make_deburst(2, 2);

        let err = deramp_subswath(&deburst, &swath, epoch).unwrap_err();
        assert!(matches!(err, InsarError::NoDcEstimates { .. }));
    }

    /// Power is conserved by deramping: |I|² + |Q|² is unchanged.
    ///
    /// The deramp is a complex rotation, so it is isometric.
    #[test]
    fn power_is_conserved() {
        let epoch = Utc.with_ymd_and_hms(2019, 1, 23, 5, 33, 48).unwrap();
        // Non-zero DC: a0 = 100 Hz (will produce a non-trivial rotation)
        let dc = vec![
            make_dc_estimate(0.0, 100.0),
            make_dc_estimate(3.0, 100.0),
        ];
        let swath = make_subswath(dc);
        let deburst = make_deburst(5, 4);

        let result = deramp_subswath(&deburst, &swath, epoch).unwrap();

        let in_power = 100.0f32 * 100.0;
        for &[i_out, q_out] in &result.data {
            let out_power = i_out * i_out + q_out * q_out;
            assert!(
                (out_power - in_power).abs() < 0.1,
                "power not conserved: in={in_power:.1} out={out_power:.3}"
            );
        }
    }

    /// Verify temporal interpolation: at the midpoint between two estimates
    /// with different a0 values, the interpolated a0 should be the mean.
    #[test]
    fn dc_temporal_interpolation() {
        // Two estimates: a0=0 Hz at t=0, a0=200 Hz at t=2 s
        let epoch = Utc.with_ymd_and_hms(2019, 1, 23, 5, 33, 48).unwrap();
        let dc = vec![
            make_dc_estimate(0.0, 0.0),
            make_dc_estimate(2.0, 200.0),
        ];
        let est_times_s = vec![0.0, 2.0];

        // At t=1 s (midpoint): interpolated a0 should be 100 Hz
        let poly = interp_dc_poly(&dc, &est_times_s, 1.0).unwrap();
        assert!(
            (poly[0] - 100.0).abs() < 1e-9,
            "interp a0 at midpoint = {:.3} — expected 100.0",
            poly[0]
        );
    }

    /// A pure-real input [I, 0] deramped by 90° (quarter-cycle) becomes [0, -I].
    ///
    /// φ = −π/2 ⟹ exp(−jφ) = exp(+jπ/2) = cos(π/2) + j·sin(π/2) = j
    /// [I, 0] · j = [0, I] → but our convention: exp(-jφ) = cos(φ) − j·sin(φ)
    /// where φ = 2π·f_dc·t_az.
    ///
    /// To get φ = π/2: need f_dc · t_az = 0.25.
    /// With ATI = 2.055556e-3 s and f_dc = 121.556... Hz: f_dc * 1*ATI = 0.25.
    #[test]
    fn quarter_cycle_rotation() {
        let epoch = Utc.with_ymd_and_hms(2019, 1, 23, 5, 33, 48).unwrap();
        // f_dc = 0.25 / ATI so that at line 1: φ = 2π * 0.25 = π/2
        let ati = 2.055_556e-3;
        let f_dc_target = 0.25 / ati; // ≈ 121.556 Hz

        let dc = vec![
            make_dc_estimate(0.0, f_dc_target),
            make_dc_estimate(3.0, f_dc_target),
        ];
        let swath = make_subswath(dc);
        // Use a 2-line, 1-sample deburst so we can check line 1 specifically.
        let deburst = DeburstArray {
            data: vec![[100i16, 0i16]; 2],
            lines: 2,
            samples: 1,
            valid_sample_offset: 0,
        };

        let result = deramp_subswath(&deburst, &swath, epoch).unwrap();

        // Line 0: t_az = 0 → φ = 0 → identity → [100, 0]
        let [i0, q0] = result.data[0];
        assert!((i0 - 100.0).abs() < 0.1, "line 0 i_out = {i0:.3}");
        assert!(q0.abs() < 0.1, "line 0 q_out = {q0:.3}");

        // Line 1: t_az = ATI → φ = π/2
        // exp(−jπ/2) = cos(π/2) − j·sin(π/2) = 0 − j
        // [100 + 0j] × (0 − j) = 0 − 100j → [0, −100]
        let [i1, q1] = result.data[1];
        assert!(i1.abs() < 0.5, "line 1 i_out = {i1:.3} — expected ≈0");
        assert!(
            (q1 + 100.0).abs() < 0.5,
            "line 1 q_out = {q1:.3} — expected ≈−100"
        );
    }

    /// [`reramp_interferogram_phase`] is the exact inverse of the deramp.
    ///
    /// After deramp the phase of each sample equals `−2π·f_DC·t_az`.
    /// After reramp it should return to the original phase (0 for pure-real input).
    #[test]
    fn reramp_inverts_deramp() {
        let epoch = Utc.with_ymd_and_hms(2019, 1, 23, 5, 33, 48).unwrap();
        // Non-zero DC so the phase rotation is measurable.
        let dc = vec![
            make_dc_estimate(0.0, 80.0),
            make_dc_estimate(3.0, 80.0),
        ];
        let swath = make_subswath(dc);
        let deburst = make_deburst(5, 4);

        // Step 1: deramp.
        let deramped = deramp_subswath(&deburst, &swath, epoch).unwrap();

        // Step 2: extract phase from the deramped output.
        let mut phase: Vec<f32> = deramped
            .data
            .iter()
            .map(|&[i, q]| q.atan2(i))
            .collect();

        // Step 3: reramp — should recover the original phase (0 everywhere,
        // since the input is pure-real [100, 0]).
        reramp_interferogram_phase(
            &mut phase,
            deramped.lines,
            deramped.samples,
            deramped.valid_sample_offset,
            &swath,
            epoch,
        )
        .unwrap();

        for (idx, &p) in phase.iter().enumerate() {
            // Allow for f32 round-trip and sin/cos precision (~1e-4 rad).
            assert!(
                p.abs() < 2e-3,
                "reramp should restore phase to 0 at idx={idx}, got {p:.5} rad"
            );
        }
    }

    /// Reramping an all-zero phase array with zero DC returns all zeros.
    #[test]
    fn reramp_zero_dc_is_identity() {
        let epoch = Utc.with_ymd_and_hms(2019, 1, 23, 5, 33, 48).unwrap();
        let dc = vec![
            make_dc_estimate(0.0, 0.0),
            make_dc_estimate(3.0, 0.0),
        ];
        let swath = make_subswath(dc);
        let mut phase = vec![0.0_f32; 6];

        reramp_interferogram_phase(&mut phase, 2, 3, 0, &swath, epoch).unwrap();

        for &p in &phase {
            assert!(p.abs() < 1e-6, "expected 0, got {p:.6}");
        }
    }

    /// Empty DC estimates must return NoDcEstimates error.
    #[test]
    fn reramp_empty_dc_returns_error() {
        let epoch = Utc.with_ymd_and_hms(2019, 1, 23, 5, 33, 48).unwrap();
        let swath = make_subswath(vec![]);
        let mut phase = vec![0.0_f32; 4];
        let err = reramp_interferogram_phase(&mut phase, 2, 2, 0, &swath, epoch).unwrap_err();
        assert!(matches!(err, InsarError::NoDcEstimates { .. }));
    }
}
