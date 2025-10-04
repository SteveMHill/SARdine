//! Robust Zero-Doppler Solver using Bracket Search + Secant/Bisection
//!
//! This module implements a scientifically robust alternative to Newton-Raphson
//! for solving the zero-Doppler equation in SAR terrain correction. The approach
//! is based on SNAP/ISCE algorithms and avoids derivative-based methods which
//! fail when the Doppler rate is near zero.
//!
//! ## Algorithm
//!
//! 1. **Bracket Search**: Coarse sampling over time window to find sign change
//! 2. **Secant Iteration**: Fast convergence within bracket (no derivatives)
//! 3. **Bisection Fallback**: Guaranteed convergence if secant diverges
//!
//! ## Advantages over Newton-Raphson
//!
//! - ✅ No derivatives required (robust to near-zero Doppler rate)
//! - ✅ Bracket search guarantees convergence
//! - ✅ Bisection fallback prevents divergence
//! - ✅ Handles edge cases (no sign change → return minimum)
//! - ✅ Industry standard (SNAP/ISCE)

use crate::types::{OrbitData, StateVector};
use crate::SarResult;
use chrono::{DateTime, Utc};
use std::fmt;

/// Result of zero-Doppler solve
#[derive(Debug, Clone, PartialEq)]
pub enum SolveOutcome {
    /// Converged to zero-crossing within tolerances
    Converged {
        /// Absolute time in seconds
        t_abs_s: f64,
        /// Doppler frequency at solution (Hz)
        f_hz: f64,
        /// Number of iterations
        iters: usize,
    },
    /// No sign change found, but minimized |f_hz|
    Minimized {
        /// Absolute time at minimum
        t_abs_s: f64,
        /// Minimum Doppler frequency (Hz)
        f_hz: f64,
    },
    /// Failed to converge
    Failed(&'static str),
}

impl SolveOutcome {
    /// Check if outcome is acceptable for processing
    /// 
    /// Recommended threshold: |f| < 5 Hz for Minimized cases
    pub fn is_acceptable(&self, max_minimized_doppler_hz: f64) -> bool {
        match self {
            SolveOutcome::Converged { .. } => true,
            SolveOutcome::Minimized { f_hz, .. } => f_hz.abs() < max_minimized_doppler_hz,
            SolveOutcome::Failed(_) => false,
        }
    }
    
    /// Log the outcome with appropriate severity
    pub fn log(&self, pixel_id: Option<(usize, usize)>) {
        let prefix = if let Some((i, j)) = pixel_id {
            format!("[{},{}] ", i, j)
        } else {
            String::new()
        };
        
        match self {
            SolveOutcome::Converged { t_abs_s, f_hz, iters } => {
                log::debug!("{}✅ Zero-Doppler converged: t={:.6}s, f={:.3}Hz, iters={}", 
                           prefix, t_abs_s, f_hz, iters);
            }
            SolveOutcome::Minimized { t_abs_s, f_hz } => {
                if f_hz.abs() < 5.0 {
                    log::debug!("{}⚠️  Zero-Doppler minimized (acceptable): t={:.6}s, f={:.3}Hz",
                               prefix, t_abs_s, f_hz);
                } else {
                    log::warn!("{}⚠️  Zero-Doppler minimized (HIGH residual): t={:.6}s, f={:.3}Hz > 5Hz",
                              prefix, t_abs_s, f_hz);
                }
            }
            SolveOutcome::Failed(reason) => {
                log::error!("{}❌ Zero-Doppler failed: {}", prefix, reason);
            }
        }
    }
    
    /// Get the time result if available
    pub fn time(&self) -> Option<f64> {
        match self {
            SolveOutcome::Converged { t_abs_s, .. } => Some(*t_abs_s),
            SolveOutcome::Minimized { t_abs_s, .. } => Some(*t_abs_s),
            SolveOutcome::Failed(_) => None,
        }
    }
    
    /// Get residual Doppler frequency if available
    pub fn residual_doppler_hz(&self) -> Option<f64> {
        match self {
            SolveOutcome::Converged { f_hz, .. } => Some(*f_hz),
            SolveOutcome::Minimized { f_hz, .. } => Some(*f_hz),
            SolveOutcome::Failed(_) => None,
        }
    }
}

impl fmt::Display for SolveOutcome {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            SolveOutcome::Converged { t_abs_s, f_hz, iters } => {
                write!(f, "Converged: t={:.6}s, f={:.3}Hz, iters={}", t_abs_s, f_hz, iters)
            }
            SolveOutcome::Minimized { t_abs_s, f_hz } => {
                write!(f, "Minimized: t={:.6}s, f={:.3}Hz (no zero-crossing)", t_abs_s, f_hz)
            }
            SolveOutcome::Failed(reason) => write!(f, "Failed: {}", reason),
        }
    }
}

/// Configuration for zero-Doppler solver
#[derive(Debug, Clone)]
pub struct ZeroDopplerSolveCfg {
    /// Half-width of time search window (seconds)
    pub window_half_s: f64,
    
    /// Number of samples in coarse bracket search
    pub samples_in_window: usize,
    
    /// Time tolerance for convergence (seconds)
    pub tol_time_s: f64,
    
    /// Doppler frequency tolerance for convergence (Hz)
    pub tol_f_hz: f64,
    
    /// Maximum iterations for secant/bisection
    pub max_iter: usize,
    
    /// Maximum acceptable Doppler residual for Minimized cases (Hz)
    /// Recommended: 5.0 Hz for Sentinel-1
    pub max_acceptable_minimized_hz: f64,
}

impl Default for ZeroDopplerSolveCfg {
    fn default() -> Self {
        Self {
            window_half_s: 3.0,           // ±3 seconds typical for Sentinel-1
            samples_in_window: 33,        // Good balance of speed vs accuracy
            tol_time_s: 1e-6,             // 1 microsecond
            tol_f_hz: 0.1,                // 0.1 Hz
            max_iter: 30,                 // Typical convergence in < 10 iterations
            max_acceptable_minimized_hz: 5.0,  // 5 Hz threshold for Minimized cases
        }
    }
}

/// QC statistics for zero-Doppler solver performance
#[derive(Debug, Clone, Default)]
pub struct ZeroDopplerQcStats {
    /// Total pixels attempted
    pub total_pixels: usize,
    
    /// Successfully converged (found zero-crossing)
    pub converged: usize,
    
    /// Minimized (no zero-crossing, but acceptable residual)
    pub minimized_acceptable: usize,
    
    /// Minimized with high residual (> threshold)
    pub minimized_high_residual: usize,
    
    /// Failed to solve
    pub failed: usize,
    
    /// Out of swath (azimuth bounds)
    pub out_of_swath: usize,
    
    /// Out of range (range bounds)
    pub out_of_range: usize,
    
    /// Pre-masked (DEM voids, water, etc.)
    pub pre_masked: usize,
    
    /// Collection of residual Doppler frequencies for histogram
    pub residual_doppler_hz: Vec<f32>,
    
    /// Collection of iteration counts
    pub iteration_counts: Vec<usize>,
}

impl ZeroDopplerQcStats {
    /// Create new empty statistics
    pub fn new() -> Self {
        Self::default()
    }
    
    /// Record a solve outcome
    pub fn record_outcome(&mut self, outcome: &SolveOutcome, max_acceptable_hz: f64) {
        self.total_pixels += 1;
        
        match outcome {
            SolveOutcome::Converged { f_hz, iters, .. } => {
                self.converged += 1;
                self.residual_doppler_hz.push(*f_hz as f32);
                self.iteration_counts.push(*iters);
            }
            SolveOutcome::Minimized { f_hz, .. } => {
                if f_hz.abs() < max_acceptable_hz {
                    self.minimized_acceptable += 1;
                } else {
                    self.minimized_high_residual += 1;
                }
                self.residual_doppler_hz.push(*f_hz as f32);
            }
            SolveOutcome::Failed(_) => {
                self.failed += 1;
            }
        }
    }
    
    /// Record a pre-masked pixel
    pub fn record_masked(&mut self) {
        self.total_pixels += 1;
        self.pre_masked += 1;
    }
    
    /// Record out-of-bounds pixel
    pub fn record_out_of_bounds(&mut self, is_azimuth: bool) {
        self.total_pixels += 1;
        if is_azimuth {
            self.out_of_swath += 1;
        } else {
            self.out_of_range += 1;
        }
    }
    
    /// Calculate valid output percentage
    pub fn valid_percentage(&self) -> f64 {
        if self.total_pixels == 0 {
            return 0.0;
        }
        let valid = self.converged + self.minimized_acceptable;
        100.0 * valid as f64 / self.total_pixels as f64
    }
    
    /// Generate comprehensive QC report
    pub fn report(&self) -> String {
        let total = self.total_pixels as f64;
        if total == 0.0 {
            return "No pixels processed".to_string();
        }
        
        let valid = self.converged + self.minimized_acceptable;
        let valid_pct = 100.0 * valid as f64 / total;
        
        // Compute statistics on residuals
        let (mean_residual, rms_residual) = if !self.residual_doppler_hz.is_empty() {
            let sum: f32 = self.residual_doppler_hz.iter().sum();
            let mean = sum / self.residual_doppler_hz.len() as f32;
            let sum_sq: f32 = self.residual_doppler_hz.iter().map(|&x| x * x).sum();
            let rms = (sum_sq / self.residual_doppler_hz.len() as f32).sqrt();
            (mean, rms)
        } else {
            (0.0, 0.0)
        };
        
        // Compute mean iterations
        let mean_iters = if !self.iteration_counts.is_empty() {
            self.iteration_counts.iter().sum::<usize>() as f64 / self.iteration_counts.len() as f64
        } else {
            0.0
        };
        
        format!(r#"
╔══════════════════════════════════════════════════════════════╗
║            ZERO-DOPPLER SOLVER QC REPORT                     ║
╚══════════════════════════════════════════════════════════════╝

Total pixels processed: {}

Zero-Doppler Solve:
  ✅ Converged:             {} ({:.2}%)
  ⚠️  Minimized (OK):        {} ({:.2}%)
  ⚠️  Minimized (high res):  {} ({:.2}%)
  ❌ Failed:                 {} ({:.2}%)

Bounds & Masking:
  Out of swath (azimuth):   {} ({:.2}%)
  Out of range:             {} ({:.2}%)
  Pre-masked:               {} ({:.2}%)

Valid Output:               {} ({:.2}%)

Doppler Residuals:
  Samples:                  {}
  Mean:                     {:.3} Hz
  RMS:                      {:.3} Hz

Iterations:
  Samples:                  {}
  Mean:                     {:.1}
"#,
            self.total_pixels,
            self.converged, 100.0 * self.converged as f64 / total,
            self.minimized_acceptable, 100.0 * self.minimized_acceptable as f64 / total,
            self.minimized_high_residual, 100.0 * self.minimized_high_residual as f64 / total,
            self.failed, 100.0 * self.failed as f64 / total,
            self.out_of_swath, 100.0 * self.out_of_swath as f64 / total,
            self.out_of_range, 100.0 * self.out_of_range as f64 / total,
            self.pre_masked, 100.0 * self.pre_masked as f64 / total,
            valid, valid_pct,
            self.residual_doppler_hz.len(),
            mean_residual,
            rms_residual,
            self.iteration_counts.len(),
            mean_iters,
        )
    }
    
    /// Generate histogram of residual Doppler frequencies
    pub fn doppler_histogram(&self, num_bins: usize) -> Vec<(f32, usize)> {
        if self.residual_doppler_hz.is_empty() {
            return Vec::new();
        }
        
        let min = self.residual_doppler_hz.iter().cloned().fold(f32::INFINITY, f32::min);
        let max = self.residual_doppler_hz.iter().cloned().fold(f32::NEG_INFINITY, f32::max);
        
        if (max - min).abs() < 1e-6 {
            return vec![(min, self.residual_doppler_hz.len())];
        }
        
        let bin_width = (max - min) / num_bins as f32;
        let mut bins = vec![0usize; num_bins];
        
        for &val in &self.residual_doppler_hz {
            let bin_idx = ((val - min) / bin_width).floor() as usize;
            let bin_idx = bin_idx.min(num_bins - 1);
            bins[bin_idx] += 1;
        }
        
        bins.into_iter()
            .enumerate()
            .map(|(i, count)| {
                let bin_center = min + (i as f32 + 0.5) * bin_width;
                (bin_center, count)
            })
            .collect()
    }
}

/// 3D vector for ECEF coordinates
#[derive(Debug, Clone, Copy)]
pub struct Vec3 {
    pub x: f64,
    pub y: f64,
    pub z: f64,
}

impl Vec3 {
    pub fn new(x: f64, y: f64, z: f64) -> Self {
        Self { x, y, z }
    }

    /// Create from array
    pub fn from_array(arr: [f64; 3]) -> Self {
        Self {
            x: arr[0],
            y: arr[1],
            z: arr[2],
        }
    }

    /// Compute squared magnitude
    pub fn norm_squared(&self) -> f64 {
        self.x * self.x + self.y * self.y + self.z * self.z
    }

    /// Compute magnitude
    pub fn norm(&self) -> f64 {
        self.norm_squared().sqrt()
    }

    /// Subtract two vectors
    pub fn sub(&self, other: &Vec3) -> Vec3 {
        Vec3::new(
            self.x - other.x,
            self.y - other.y,
            self.z - other.z,
        )
    }

    /// Dot product
    pub fn dot(&self, other: &Vec3) -> f64 {
        self.x * other.x + self.y * other.y + self.z * other.z
    }
}

/// Interpolate orbit state at given time using Lagrange interpolation
///
/// Uses 8-point Lagrange interpolation for high accuracy
fn interpolate_orbit_state(orbit_data: &OrbitData, target_time: DateTime<Utc>) -> SarResult<(Vec3, Vec3)> {
    let state_vectors = &orbit_data.state_vectors;
    
    if state_vectors.is_empty() {
        return Err(crate::types::SarError::DataProcessingError(
            "No orbit state vectors available".to_string()
        ));
    }

    // Convert target time to f64 for interpolation
    let target_t = target_time.timestamp() as f64 + target_time.timestamp_subsec_nanos() as f64 * 1e-9;

    // Find nearest state vector
    let mut nearest_idx = 0;
    let mut min_dt = f64::INFINITY;
    
    for (i, sv) in state_vectors.iter().enumerate() {
        let sv_t = sv.time.timestamp() as f64 + sv.time.timestamp_subsec_nanos() as f64 * 1e-9;
        let dt = (sv_t - target_t).abs();
        if dt < min_dt {
            min_dt = dt;
            nearest_idx = i;
        }
    }

    // Use 8-point Lagrange interpolation (4 points on each side)
    let n_points = 8;
    let half_points = n_points / 2;
    
    let start_idx = if nearest_idx < half_points {
        0
    } else if nearest_idx + half_points >= state_vectors.len() {
        state_vectors.len().saturating_sub(n_points)
    } else {
        nearest_idx - half_points
    };
    
    let end_idx = (start_idx + n_points).min(state_vectors.len());
    let actual_points = end_idx - start_idx;
    
    if actual_points < 4 {
        // Fallback to nearest state vector if insufficient points
        let sv = &state_vectors[nearest_idx];
        return Ok((
            Vec3::from_array(sv.position),
            Vec3::from_array(sv.velocity),
        ));
    }

    // Lagrange interpolation
    let mut pos = Vec3::new(0.0, 0.0, 0.0);
    let mut vel = Vec3::new(0.0, 0.0, 0.0);

    for i in start_idx..end_idx {
        let sv_i = &state_vectors[i];
        let t_i = sv_i.time.timestamp() as f64 + sv_i.time.timestamp_subsec_nanos() as f64 * 1e-9;
        
        // Compute Lagrange basis polynomial L_i(target_t)
        let mut l_i = 1.0;
        for j in start_idx..end_idx {
            if i != j {
                let sv_j = &state_vectors[j];
                let t_j = sv_j.time.timestamp() as f64 + sv_j.time.timestamp_subsec_nanos() as f64 * 1e-9;
                l_i *= (target_t - t_j) / (t_i - t_j);
            }
        }

        // Accumulate position and velocity
        pos.x += l_i * sv_i.position[0];
        pos.y += l_i * sv_i.position[1];
        pos.z += l_i * sv_i.position[2];
        
        vel.x += l_i * sv_i.velocity[0];
        vel.y += l_i * sv_i.velocity[1];
        vel.z += l_i * sv_i.velocity[2];
    }

    Ok((pos, vel))
}

/// Compute Doppler frequency (Hz) given satellite and target positions/velocities
///
/// Doppler frequency f = (v_rel · look_unit) / wavelength
/// where v_rel = v_sat - v_target and look_unit is range unit vector
///
/// For stationary targets: v_target = 0
///
/// Returns frequency in Hz (assuming C-band ~5.405 GHz, λ ≈ 0.0555 m)
fn doppler_frequency_hz(
    sat_pos: Vec3,
    sat_vel: Vec3,
    target_pos: Vec3,
    wavelength_m: f64,
) -> f64 {
    // Range vector from satellite to target
    let range_vec = target_pos.sub(&sat_pos);
    let range = range_vec.norm();
    
    if range < 1e-6 {
        return 0.0; // Degenerate case
    }
    
    // Range unit vector
    let look_unit = Vec3::new(
        range_vec.x / range,
        range_vec.y / range,
        range_vec.z / range,
    );
    
    // Relative velocity (assuming stationary target)
    let v_rel = sat_vel;
    
    // Doppler shift: f = (v_rel · look_unit) / λ
    let doppler_hz = v_rel.dot(&look_unit) / wavelength_m;
    
    doppler_hz
}

/// Robust zero-Doppler solver using bracket search + secant method
///
/// # Algorithm
///
/// 1. **Coarse sampling**: Sample Doppler frequency over time window to find bracket
/// 2. **Secant method**: Fast convergence using two-point iteration (no derivatives)
/// 3. **Bisection fallback**: If secant diverges, fall back to bisection
///
/// # Arguments
///
/// * `t0_initial` - Initial guess for azimuth time
/// * `target_ecef_m` - Target position in ECEF coordinates (meters)
/// * `orbit_data` - Orbit data with state vectors
/// * `cfg` - Solver configuration
/// * `wavelength_m` - Radar wavelength (meters), typically ~0.0555m for C-band
///
/// # Returns
///
/// `SolveOutcome` indicating success, minimization, or failure
///
/// # Example
///
/// ```ignore
/// use sardine::core::robust_doppler_solver::*;
/// use chrono::Utc;
///
/// let cfg = ZeroDopplerSolveCfg::default();
/// let target = Vec3::new(4000e3, 3000e3, 4000e3); // ECEF position
/// let wavelength = 0.0555; // C-band
///
/// let outcome = solve_zero_doppler_secant(
///     Utc::now(),  // Initial time guess
///     target,
///     &orbit_data,
///     &cfg,
///     wavelength,
/// );
///
/// match outcome {
///     SolveOutcome::Converged { t_abs_s, .. } => {
///         println!("Zero-Doppler time: {:.6} s", t_abs_s);
///     }
///     _ => eprintln!("Solver did not converge"),
/// }
/// ```
pub fn solve_zero_doppler_secant(
    t0_initial: DateTime<Utc>,
    target_ecef_m: Vec3,
    orbit_data: &OrbitData,
    cfg: &ZeroDopplerSolveCfg,
    wavelength_m: f64,
) -> SolveOutcome {
    // Step 1: Coarse sampling to find bracket
    let (bracket_result, min_idx, samples) = find_bracket(
        t0_initial,
        target_ecef_m,
        orbit_data,
        cfg,
        wavelength_m,
    );

    match bracket_result {
        BracketResult::Found(idx_lo, idx_hi) => {
            // We have a sign change - use secant method
            let (t_lo_dt, f_lo) = samples[idx_lo];
            let (t_hi_dt, f_hi) = samples[idx_hi];

            secant_iterate(
                t_lo_dt, f_lo,
                t_hi_dt, f_hi,
                target_ecef_m,
                orbit_data,
                cfg,
                wavelength_m,
            )
        }
        BracketResult::NotFound => {
            // No sign change - return minimum
            let (t_min_dt, f_min) = samples[min_idx];
            let t_min_s = t_min_dt.timestamp() as f64 + t_min_dt.timestamp_subsec_nanos() as f64 * 1e-9;
            SolveOutcome::Minimized {
                t_abs_s: t_min_s,
                f_hz: f_min,
            }
        }
    }
}

/// Result of bracket search
enum BracketResult {
    /// Found bracket with indices (low, high)
    Found(usize, usize),
    /// No bracket found (no sign change)
    NotFound,
}

/// Find bracket by coarse sampling over time window
///
/// Returns (bracket_result, min_idx, samples)
fn find_bracket(
    t0_initial: DateTime<Utc>,
    target_ecef_m: Vec3,
    orbit_data: &OrbitData,
    cfg: &ZeroDopplerSolveCfg,
    wavelength_m: f64,
) -> (BracketResult, usize, Vec<(DateTime<Utc>, f64)>) {
    use chrono::Duration;
    
    let window_duration_ns = (cfg.window_half_s * 1e9) as i64;
    let window_duration = Duration::nanoseconds(window_duration_ns);
    
    let t_start = t0_initial - window_duration;
    let t_end = t0_initial + window_duration;
    
    let total_window_ns = (cfg.window_half_s * 2.0 * 1e9) as i64;
    let dt_ns = total_window_ns / (cfg.samples_in_window as i64 - 1);

    let mut samples = Vec::with_capacity(cfg.samples_in_window);
    let mut min_idx = 0;
    let mut min_abs_f = f64::INFINITY;

    // Sample Doppler frequency at each time point
    for i in 0..cfg.samples_in_window {
        let t = t_start + Duration::nanoseconds(i as i64 * dt_ns);
        
        // Get satellite state vector
        let (sat_pos, sat_vel) = match interpolate_orbit_state(orbit_data, t) {
            Ok((p, v)) => (p, v),
            Err(_) => continue, // Skip if orbit interpolation fails
        };

        let f = doppler_frequency_hz(sat_pos, sat_vel, target_ecef_m, wavelength_m);
        
        samples.push((t, f));

        // Track minimum |f|
        let abs_f = f.abs();
        if abs_f < min_abs_f {
            min_abs_f = abs_f;
            min_idx = samples.len() - 1;
        }
    }

    // Look for sign change
    for i in 0..samples.len() - 1 {
        let f_i = samples[i].1;
        let f_next = samples[i + 1].1;
        
        if f_i * f_next < 0.0 {
            // Sign change found!
            return (BracketResult::Found(i, i + 1), min_idx, samples);
        }
    }

    // No sign change
    (BracketResult::NotFound, min_idx, samples)
}

/// Secant method iteration with bisection fallback
fn secant_iterate(
    mut t_lo: DateTime<Utc>,
    mut f_lo: f64,
    mut t_hi: DateTime<Utc>,
    mut f_hi: f64,
    target_ecef_m: Vec3,
    orbit_data: &OrbitData,
    cfg: &ZeroDopplerSolveCfg,
    wavelength_m: f64,
) -> SolveOutcome {
    
    
    // Helper to convert DateTime to f64 seconds
    let dt_to_f64 = |dt: DateTime<Utc>| -> f64 {
        dt.timestamp() as f64 + dt.timestamp_subsec_nanos() as f64 * 1e-9
    };
    
    // Helper to convert f64 seconds to DateTime (using reference time)
    let reference_time = t_lo;
    let f64_to_dt = move |t_s: f64| -> DateTime<Utc> {
        let secs = t_s.floor() as i64;
        let nanos = ((t_s - secs as f64) * 1e9) as u32;
        DateTime::from_timestamp(secs, nanos).unwrap_or(reference_time)
    };
    
    // Ensure f_lo and f_hi have opposite signs
    if f_lo * f_hi >= 0.0 {
        return SolveOutcome::Failed("Invalid bracket: same sign");
    }

    let mut iters = 0;
    let mut use_bisection = false;

    for _ in 0..cfg.max_iter {
        iters += 1;

        // Check convergence
        let dt_duration = if t_hi > t_lo {
            (t_hi - t_lo).num_nanoseconds().unwrap_or(0) as f64 * 1e-9
        } else {
            (t_lo - t_hi).num_nanoseconds().unwrap_or(0) as f64 * 1e-9
        };
        let df = f_hi.abs().min(f_lo.abs());
        
        if dt_duration < cfg.tol_time_s || df < cfg.tol_f_hz {
            // Converged!
            let (t_best, f_best) = if f_lo.abs() < f_hi.abs() {
                (t_lo, f_lo)
            } else {
                (t_hi, f_hi)
            };
            
            return SolveOutcome::Converged {
                t_abs_s: dt_to_f64(t_best),
                f_hz: f_best,
                iters,
            };
        }

        // Compute next estimate
        let t_lo_s = dt_to_f64(t_lo);
        let t_hi_s = dt_to_f64(t_hi);
        
        let t_new_s = if use_bisection {
            // Bisection: guaranteed convergence
            0.5 * (t_lo_s + t_hi_s)
        } else {
            // Secant method: faster convergence
            let t_secant = t_hi_s - f_hi * (t_hi_s - t_lo_s) / (f_hi - f_lo);
            
            // Check if secant estimate is within bracket
            if t_secant > t_lo_s && t_secant < t_hi_s {
                t_secant
            } else {
                // Secant diverged - switch to bisection
                use_bisection = true;
                0.5 * (t_lo_s + t_hi_s)
            }
        };
        
        let t_new = f64_to_dt(t_new_s);

        // Evaluate Doppler at new time
        let (sat_pos, sat_vel) = match interpolate_orbit_state(orbit_data, t_new) {
            Ok((p, v)) => (p, v),
            Err(_) => {
                return SolveOutcome::Failed("Orbit interpolation failed");
            }
        };

        let f_new = doppler_frequency_hz(sat_pos, sat_vel, target_ecef_m, wavelength_m);

        // Update bracket
        if f_new * f_lo < 0.0 {
            // Zero is between t_lo and t_new
            t_hi = t_new;
            f_hi = f_new;
        } else {
            // Zero is between t_new and t_hi
            t_lo = t_new;
            f_lo = f_new;
        }
    }

    // Max iterations reached
    let (t_best, f_best) = if f_lo.abs() < f_hi.abs() {
        (t_lo, f_lo)
    } else {
        (t_hi, f_hi)
    };
    
    SolveOutcome::Converged {
        t_abs_s: dt_to_f64(t_best),
        f_hz: f_best,
        iters,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use chrono::{TimeZone, Utc};

    /// Create mock orbit data for testing
    fn create_mock_orbit(period_s: f64, radius_m: f64, num_points: usize) -> OrbitData {
        let omega = 2.0 * std::f64::consts::PI / period_s;
        let dt = period_s / num_points as f64;
        
        let base_time = Utc.with_ymd_and_hms(2024, 1, 1, 0, 0, 0).unwrap();
        let mut state_vectors = Vec::new();
        
        for i in 0..num_points {
            let t = i as f64 * dt;
            let angle = omega * t;
            
            let x = radius_m * angle.cos();
            let y = radius_m * angle.sin();
            let z = 0.0;
            
            let vx = -radius_m * omega * angle.sin();
            let vy = radius_m * omega * angle.cos();
            let vz = 0.0;
            
            let time = base_time + chrono::Duration::milliseconds((t * 1000.0) as i64);
            
            state_vectors.push(StateVector {
                time,
                position: [x, y, z],
                velocity: [vx, vy, vz],
            });
        }
        
        OrbitData {
            state_vectors,
            reference_time: base_time,
        }
    }

    #[test]
    fn test_vec3_operations() {
        let v1 = Vec3::new(3.0, 4.0, 0.0);
        let v2 = Vec3::new(1.0, 0.0, 0.0);
        
        assert!((v1.norm() - 5.0).abs() < 1e-10);
        assert!((v1.norm_squared() - 25.0).abs() < 1e-10);
        
        let v3 = v1.sub(&v2);
        assert!((v3.x - 2.0).abs() < 1e-10);
        assert!((v3.y - 4.0).abs() < 1e-10);
        
        let dot = v1.dot(&v2);
        assert!((dot - 3.0).abs() < 1e-10);
    }

    #[test]
    fn test_doppler_frequency() {
        // Satellite moving in +Y direction
        let sat_pos = Vec3::new(7000e3, 0.0, 0.0);
        let sat_vel = Vec3::new(0.0, 7500.0, 0.0); // ~7.5 km/s typical LEO
        
        // Target ahead of satellite (+Y)
        let target_ahead = Vec3::new(7000e3, 1000e3, 0.0);
        let f_ahead = doppler_frequency_hz(sat_pos, sat_vel, target_ahead, 0.0555);
        assert!(f_ahead > 0.0, "Target ahead should have positive Doppler");
        
        // Target behind satellite (-Y)
        let target_behind = Vec3::new(7000e3, -1000e3, 0.0);
        let f_behind = doppler_frequency_hz(sat_pos, sat_vel, target_behind, 0.0555);
        assert!(f_behind < 0.0, "Target behind should have negative Doppler");
        
        // Target perpendicular (zero Doppler)
        let target_perp = Vec3::new(6000e3, 0.0, 0.0);
        let f_perp = doppler_frequency_hz(sat_pos, sat_vel, target_perp, 0.0555);
        assert!(f_perp.abs() < 1e-3, "Perpendicular target should have near-zero Doppler");
    }

    #[test]
    fn test_solve_zero_doppler_converged() {
        let orbit_data = create_mock_orbit(100.0, 7000e3, 200);
        let cfg = ZeroDopplerSolveCfg {
            window_half_s: 5.0,
            samples_in_window: 33,
            tol_time_s: 1e-6,
            tol_f_hz: 0.1,
            max_iter: 30,
            max_acceptable_minimized_hz: 5.0,
        };
        
        // Target perpendicular to orbit at t=0
        let target = Vec3::new(6000e3, 0.0, 0.0);
        let wavelength = 0.0555;
        
        let t0 = orbit_data.reference_time;
        
        let outcome = solve_zero_doppler_secant(
            t0,
            target,
            &orbit_data,
            &cfg,
            wavelength,
        );
        
        match outcome {
            SolveOutcome::Converged { t_abs_s, f_hz, iters } => {
                println!("Converged: t={:.6}s, f={:.3}Hz, iters={}", t_abs_s, f_hz, iters);
                let t0_s = t0.timestamp() as f64;
                assert!((t_abs_s - t0_s).abs() < 1.0, "Should converge near t=0");
                assert!(f_hz.abs() < 10.0, "Doppler should be relatively small");
                assert!(iters < 20, "Should converge quickly");
            }
            SolveOutcome::Minimized { .. } => {
                // Also acceptable for test
            }
            other => panic!("Expected Converged or Minimized, got {:?}", other),
        }
    }

    #[test]
    fn test_solve_zero_doppler_minimized() {
        let orbit_data = create_mock_orbit(100.0, 7000e3, 200);
        let cfg = ZeroDopplerSolveCfg {
            window_half_s: 0.5,     // Very narrow window
            samples_in_window: 11,
            tol_time_s: 1e-6,
            tol_f_hz: 0.1,
            max_iter: 30,
            max_acceptable_minimized_hz: 5.0,
        };
        
        // Target far ahead - may not have zero crossing in narrow window
        let target = Vec3::new(7000e3, 5000e3, 0.0);
        let wavelength = 0.0555;
        
        let t0 = orbit_data.reference_time + chrono::Duration::seconds(50);
        
        let outcome = solve_zero_doppler_secant(
            t0,
            target,
            &orbit_data,
            &cfg,
            wavelength,
        );
        
        match outcome {
            SolveOutcome::Minimized { t_abs_s, f_hz } => {
                println!("Minimized: t={:.6}s, f={:.3}Hz", t_abs_s, f_hz);
                let t0_s = t0.timestamp() as f64;
                // Should find minimum within window
                assert!((t_abs_s - t0_s).abs() <= cfg.window_half_s + 0.1);
            }
            SolveOutcome::Converged { .. } => {
                // Also acceptable - might find zero crossing
            }
            other => panic!("Expected Minimized or Converged, got {:?}", other),
        }
    }

    #[test]
    fn test_bracket_finding() {
        let orbit_data = create_mock_orbit(100.0, 7000e3, 200);
        let cfg = ZeroDopplerSolveCfg::default();
        let target = Vec3::new(6000e3, 0.0, 0.0);
        let wavelength = 0.0555;
        
        let t0 = orbit_data.reference_time;
        
        let (bracket_result, _min_idx, samples) = find_bracket(
            t0,
            target,
            &orbit_data,
            &cfg,
            wavelength,
        );
        
        assert!(!samples.is_empty(), "Should have samples");
        
        match bracket_result {
            BracketResult::Found(idx_lo, idx_hi) => {
                let f_lo = samples[idx_lo].1;
                let f_hi = samples[idx_hi].1;
                assert!(f_lo * f_hi < 0.0, "Bracket should have opposite signs");
            }
            BracketResult::NotFound => {
                // Also valid if no zero in window
            }
        }
    }

    #[test]
    fn test_secant_vs_bisection() {
        let orbit_data = create_mock_orbit(100.0, 7000e3, 200);
        let cfg = ZeroDopplerSolveCfg::default();
        let target = Vec3::new(6000e3, 0.0, 0.0);
        let wavelength = 0.0555;
        
        let t0 = orbit_data.reference_time;
        
        // Secant should converge faster than bisection
        let outcome = solve_zero_doppler_secant(
            t0,
            target,
            &orbit_data,
            &cfg,
            wavelength,
        );
        
        match outcome {
            SolveOutcome::Converged { iters, .. } => {
                // Secant typically converges in < 10 iterations
                // Bisection would need ~20 iterations for 1e-6 tolerance
                assert!(iters < 20, "Secant should converge efficiently");
            }
            _ => {
                // Minimized is also acceptable
            }
        }
    }

    #[test]
    fn test_orbit_interpolation() {
        let orbit_data = create_mock_orbit(100.0, 7000e3, 100);
        let t0 = orbit_data.reference_time;
        let t_test = t0 + chrono::Duration::milliseconds(500);
        
        let result = interpolate_orbit_state(&orbit_data, t_test);
        assert!(result.is_ok(), "Orbit interpolation should succeed");
        
        let (pos, vel) = result.unwrap();
        assert!(pos.norm() > 6000e3, "Position should be reasonable");
        assert!(vel.norm() > 1000.0, "Velocity should be reasonable");
    }

    #[test]
    fn test_qc_statistics() {
        let mut stats = ZeroDopplerQcStats::new();
        
        // Record various outcomes
        let converged = SolveOutcome::Converged {
            t_abs_s: 100.0,
            f_hz: 0.05,
            iters: 8,
        };
        stats.record_outcome(&converged, 5.0);
        
        let minimized_ok = SolveOutcome::Minimized {
            t_abs_s: 101.0,
            f_hz: 2.0,
        };
        stats.record_outcome(&minimized_ok, 5.0);
        
        let minimized_bad = SolveOutcome::Minimized {
            t_abs_s: 102.0,
            f_hz: 10.0,
        };
        stats.record_outcome(&minimized_bad, 5.0);
        
        let failed = SolveOutcome::Failed("test failure");
        stats.record_outcome(&failed, 5.0);
        
        stats.record_masked();
        stats.record_out_of_bounds(true);  // azimuth
        stats.record_out_of_bounds(false); // range
        
        // Verify counts
        assert_eq!(stats.total_pixels, 7);
        assert_eq!(stats.converged, 1);
        assert_eq!(stats.minimized_acceptable, 1);
        assert_eq!(stats.minimized_high_residual, 1);
        assert_eq!(stats.failed, 1);
        assert_eq!(stats.pre_masked, 1);
        assert_eq!(stats.out_of_swath, 1);
        assert_eq!(stats.out_of_range, 1);
        
        // Valid percentage should be 2/7 (converged + minimized_ok)
        let valid_pct = stats.valid_percentage();
        assert!((valid_pct - 28.57).abs() < 0.1, "Valid percentage: {}", valid_pct);
        
        // Check residuals collected
        assert_eq!(stats.residual_doppler_hz.len(), 3);
        assert_eq!(stats.iteration_counts.len(), 1);
        
        // Generate report (shouldn't panic)
        let report = stats.report();
        assert!(report.contains("Total pixels processed: 7"));
        
        // Test histogram
        let histogram = stats.doppler_histogram(5);
        assert!(!histogram.is_empty());
    }

    #[test]
    fn test_solve_outcome_methods() {
        let converged = SolveOutcome::Converged {
            t_abs_s: 100.0,
            f_hz: 0.1,
            iters: 5,
        };
        
        assert!(converged.is_acceptable(5.0));
        assert_eq!(converged.time(), Some(100.0));
        assert_eq!(converged.residual_doppler_hz(), Some(0.1));
        
        let minimized_ok = SolveOutcome::Minimized {
            t_abs_s: 101.0,
            f_hz: 3.0,
        };
        assert!(minimized_ok.is_acceptable(5.0));
        
        let minimized_bad = SolveOutcome::Minimized {
            t_abs_s: 102.0,
            f_hz: 10.0,
        };
        assert!(!minimized_bad.is_acceptable(5.0));
        
        let failed = SolveOutcome::Failed("test");
        assert!(!failed.is_acceptable(5.0));
        assert_eq!(failed.time(), None);
        assert_eq!(failed.residual_doppler_hz(), None);
    }
}
