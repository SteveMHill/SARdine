//! Unified Zero-Doppler Solver for Range-Doppler Terrain Correction
//!
//! This module consolidates multiple redundant zero-Doppler solvers into one
//! configurable, robust implementation based on Cumming & Wong (2005).
//!
//! # Algorithm
//! The zero-Doppler condition finds the azimuth time when the Doppler frequency
//! between satellite and ground target is zero (perpendicular geometry).
//!
//! f_D = -2 * (r⃗ · v⃗) / (λ * |r⃗|) = 0
//!
//! # Implementation
//! Uses bracketed Newton-Raphson with analytic derivative and bisection fallback.
//! Supports optional seeding for faster convergence in dense grid processing.

use super::{OrbitData, RangeDopplerParams, TerrainCorrector};

/// Configuration for zero-Doppler solver behavior
#[derive(Debug, Clone)]
pub struct ZeroDopplerOptions {
    /// Maximum iterations before giving up
    pub max_iterations: usize,
    /// Doppler convergence threshold in Hz (typically 1e-3)
    pub doppler_tolerance_hz: f64,
    /// Time step convergence threshold in seconds
    pub time_tolerance_s: f64,
    /// Maximum allowed time step per iteration (prevents divergence)
    pub max_time_step_s: f64,
    /// Optional seed time for faster convergence (orbit-relative seconds)
    pub seed_time: Option<f64>,
    /// Whether to use bisection fallback when Newton fails
    pub use_bisection_fallback: bool,
    /// Enable detailed logging for debugging
    pub verbose: bool,
}

impl Default for ZeroDopplerOptions {
    fn default() -> Self {
        Self {
            max_iterations: 30,
            doppler_tolerance_hz: 1e-3, // 1 mHz
            time_tolerance_s: 1e-6,     // 1 microsecond
            max_time_step_s: 5.0,       // 5 seconds max step
            seed_time: None,
            use_bisection_fallback: true,
            verbose: false,
        }
    }
}

impl ZeroDopplerOptions {
    /// Create options optimized for tie-point grid (coarse, fast)
    pub fn for_tie_point_grid() -> Self {
        Self {
            max_iterations: 20,
            doppler_tolerance_hz: 1e-2, // 10 mHz - slightly relaxed
            time_tolerance_s: 1e-5,
            max_time_step_s: 5.0,
            seed_time: None,
            use_bisection_fallback: true,
            verbose: false,
        }
    }

    /// Create options for seeded refinement (fast convergence with good initial guess)
    pub fn for_seeded_refinement(seed: f64) -> Self {
        Self {
            max_iterations: 5, // Should converge in 1-2 iterations
            doppler_tolerance_hz: 1e-3,
            time_tolerance_s: 1e-6,
            max_time_step_s: 2.0, // Smaller steps since we're close
            seed_time: Some(seed),
            use_bisection_fallback: false, // Seed should be good enough
            verbose: false,
        }
    }

    /// Create options for high-precision processing
    pub fn high_precision() -> Self {
        Self {
            max_iterations: 50,
            doppler_tolerance_hz: 1e-4, // 0.1 mHz
            time_tolerance_s: 1e-7,
            max_time_step_s: 3.0,
            seed_time: None,
            use_bisection_fallback: true,
            verbose: false,
        }
    }
}

/// Result of zero-Doppler solving
#[derive(Debug, Clone)]
pub struct ZeroDopplerResult {
    /// Azimuth time in orbit-relative seconds
    pub time_rel_s: f64,
    /// Final Doppler frequency achieved (should be near zero)
    pub doppler_hz: f64,
    /// Number of iterations used
    pub iterations: usize,
    /// Whether the solution converged
    pub converged: bool,
}

impl TerrainCorrector {
    /// Unified zero-Doppler solver with configurable behavior
    ///
    /// # Arguments
    /// * `target_ecef` - Target point in ECEF coordinates [x, y, z] meters
    /// * `orbit_data` - Satellite orbit state vectors
    /// * `params` - Range-Doppler processing parameters
    /// * `options` - Solver configuration
    ///
    /// # Returns
    /// `Some(time)` - Orbit-relative azimuth time in seconds when Doppler is zero
    /// `None` - If no solution found within constraints
    pub fn solve_zero_doppler(
        &self,
        target_ecef: &[f64; 3],
        orbit_data: &OrbitData,
        params: &RangeDopplerParams,
        options: &ZeroDopplerOptions,
    ) -> Option<ZeroDopplerResult> {
        if orbit_data.state_vectors.is_empty() {
            log::debug!("solve_zero_doppler: No orbit state vectors");
            return None;
        }

        // Get time bounds (orbit-relative)
        // AUDIT FIX: Safe access to orbit vectors (already checked is_empty above)
        let orbit_ref_epoch = crate::types::datetime_to_utc_seconds(orbit_data.reference_time);
        let t_start = crate::types::datetime_to_utc_seconds(orbit_data.state_vectors[0].time)
            - orbit_ref_epoch;
        let last_sv = orbit_data.state_vectors.last()?; // Returns None if empty
        let t_end = crate::types::datetime_to_utc_seconds(last_sv.time) - orbit_ref_epoch;

        // Constrain to product window if available
        let (search_lo, search_hi) = if params.product_duration > 0.0 {
            let product_end = params.product_start_rel_s + params.product_duration;
            let margin = 0.5; // Small margin for edge cases
            (
                (params.product_start_rel_s - margin).max(t_start),
                (product_end + margin).min(t_end),
            )
        } else {
            (t_start, t_end)
        };

        if search_hi <= search_lo {
            log::debug!(
                "solve_zero_doppler: Invalid search window [{:.3}, {:.3}]",
                search_lo,
                search_hi
            );
            return None;
        }

        // Determine initial guess
        let initial_time = match options.seed_time {
            Some(seed) => seed.clamp(search_lo, search_hi),
            None => {
                // Default: use product mid-time or search mid-time
                if params.product_duration > 0.0 {
                    params.product_start_rel_s + params.product_duration / 2.0
                } else {
                    (search_lo + search_hi) / 2.0
                }
            }
        };

        // Try Newton-Raphson first
        if let Some(result) = self.newton_raphson_solve(
            target_ecef,
            orbit_data,
            params,
            initial_time,
            search_lo,
            search_hi,
            options,
        ) {
            return Some(result);
        }

        // Fallback to bracketed bisection if enabled
        if options.use_bisection_fallback {
            return self.bracketed_bisection_solve(
                target_ecef,
                orbit_data,
                params,
                search_lo,
                search_hi,
                options,
            );
        }

        None
    }

    /// Newton-Raphson solver with analytic derivative for zero-Doppler time calculation
    ///
    /// ## Mathematical Basis
    /// Solves the zero-Doppler condition: f_D(t) = -2·(r⃗(t)·v⃗(t)) / (λ·|r⃗(t)|) = 0
    ///
    /// Where:
    /// - r⃗(t) = target_ecef - sat_pos(t): Range vector from satellite to target
    /// - v⃗(t): Satellite velocity vector at time t
    /// - λ: Radar wavelength
    /// - f_D: Doppler frequency (Hz)
    ///
    /// ## Newton-Raphson Iteration
    /// Iterative update: t_{n+1} = t_n - f(t_n) / f'(t_n)
    ///
    /// Where:
    /// - f(t) = doppler_freq: The zero-Doppler condition function
    /// - f'(t) = doppler_deriv: Analytic derivative of Doppler frequency
    ///
    /// ## Convergence Properties
    /// - **Quadratic convergence**: When well-conditioned, error decreases quadratically
    /// - **Typical iterations**: 3-5 for most cases (well-conditioned)
    /// - **Convergence criteria**: |f_D| < tolerance (typically 1e-3 Hz) OR |Δt| < time_tolerance
    ///
    /// ## Analytic Derivative Derivation
    /// f_D = -2·(r⃗·v⃗) / (λ·|r⃗|)
    ///
    /// Using chain rule and assuming constant target position:
    /// df_D/dt = -2/λ · d/dt[(r⃗·v⃗) / |r⃗|]
    ///
    /// After simplification:
    /// df_D/dt = -2/λ · [(R_dot² - |v⃗|²) / |r⃗|]
    ///
    /// Where R_dot = (r⃗·v⃗) / |r⃗| is the range rate.
    ///
    /// Reference: Press et al. (2007) "Numerical Recipes", Cumming & Wong (2005)
    fn newton_raphson_solve(
        &self,
        target_ecef: &[f64; 3],
        orbit_data: &OrbitData,
        params: &RangeDopplerParams,
        initial_time: f64,
        t_lo: f64,
        t_hi: f64,
        options: &ZeroDopplerOptions,
    ) -> Option<ZeroDopplerResult> {
        let orbit_ref_epoch = crate::types::datetime_to_utc_seconds(orbit_data.reference_time);
        let mut t = initial_time;

        for iteration in 0..options.max_iterations {
            // Bounds check: ensure time stays within valid search window
            if t < t_lo || t > t_hi {
                if options.verbose {
                    log::debug!(
                        "Newton-Raphson: time {:.6}s outside bounds [{:.6}, {:.6}]",
                        t,
                        t_lo,
                        t_hi
                    );
                }
                return None;
            }

            // Step 1: Interpolate satellite state (position + velocity) at current time
            // Uses cubic spline interpolation for high precision
            let absolute_time = t + orbit_ref_epoch;
            let (sat_pos, sat_vel) =
                match self.scientific_orbit_interpolation(orbit_data, absolute_time) {
                    Ok(state) => state,
                    Err(_) => return None,
                };

            // Step 2: Calculate range vector r⃗ = target - satellite
            let range_vec = [
                target_ecef[0] - sat_pos.x,
                target_ecef[1] - sat_pos.y,
                target_ecef[2] - sat_pos.z,
            ];
            let range_mag =
                (range_vec[0].powi(2) + range_vec[1].powi(2) + range_vec[2].powi(2)).sqrt();

            // Sanity check: range should be realistic (1000m to 10,000km)
            if range_mag < 1000.0 || range_mag > 1e7 {
                return None; // Unrealistic range - likely numerical error
            }

            // Step 3: Calculate Doppler frequency f_D = -2·(r⃗·v⃗) / (λ·|r⃗|)
            // This is the function we're solving: f_D(t) = 0
            let range_dot_vel =
                range_vec[0] * sat_vel.x + range_vec[1] * sat_vel.y + range_vec[2] * sat_vel.z;
            let range_rate = range_dot_vel / range_mag; // R_dot = (r⃗·v⃗) / |r⃗|
            let doppler_freq = -2.0 * range_rate / params.wavelength;

            // Step 4: Check convergence - zero-Doppler condition satisfied?
            // Convergence criterion: |f_D| < tolerance (typically 1 mHz = 1e-3 Hz)
            if doppler_freq.abs() < options.doppler_tolerance_hz {
                return Some(ZeroDopplerResult {
                    time_rel_s: t,
                    doppler_hz: doppler_freq,
                    iterations: iteration + 1,
                    converged: true,
                });
            }

            // Step 5: Calculate analytic derivative f'(t) = df_D/dt
            //
            // Derivation:
            //   f_D = -2·R_dot / λ  where R_dot = (r⃗·v⃗) / |r⃗|
            //   df_D/dt = -2/λ · dR_dot/dt
            //
            // Using: dR_dot/dt = (R_dot² - |v⃗|²) / |r⃗|
            // (This comes from differentiating R_dot = (r⃗·v⃗) / |r⃗| with respect to time,
            //  assuming constant target position and using satellite acceleration)
            let vel_mag_sq = sat_vel.x.powi(2) + sat_vel.y.powi(2) + sat_vel.z.powi(2);
            let range_rate_sq = range_rate * range_rate;
            let range_rate_deriv = (range_rate_sq - vel_mag_sq) / range_mag;
            let doppler_deriv = -2.0 * range_rate_deriv / params.wavelength;

            // Step 6: Newton-Raphson update: t_{n+1} = t_n - f(t_n) / f'(t_n)
            // Check for zero derivative (can't make progress)
            if doppler_deriv.abs() < 1e-12 {
                if options.verbose {
                    log::debug!(
                        "Newton-Raphson: derivative too small at iteration {}",
                        iteration
                    );
                }
                return None; // Can't make progress - likely at local extremum or numerical issue
            }

            // Compute time step: Δt = -f(t) / f'(t)
            let time_step = -doppler_freq / doppler_deriv;

            // Clamp step size to prevent divergence (max 5 seconds per iteration)
            let clamped_step = time_step.clamp(-options.max_time_step_s, options.max_time_step_s);

            // Alternative convergence check: time step is negligible
            if clamped_step.abs() < options.time_tolerance_s {
                return Some(ZeroDopplerResult {
                    time_rel_s: t,
                    doppler_hz: doppler_freq,
                    iterations: iteration + 1,
                    converged: true,
                });
            }

            // Update time for next iteration
            t += clamped_step;
        }

        // Did not converge within max_iterations
        None
    }

    /// Bracketed bisection solver with sign-change detection
    fn bracketed_bisection_solve(
        &self,
        target_ecef: &[f64; 3],
        orbit_data: &OrbitData,
        params: &RangeDopplerParams,
        t_lo: f64,
        t_hi: f64,
        options: &ZeroDopplerOptions,
    ) -> Option<ZeroDopplerResult> {
        const BRACKET_SAMPLES: usize = 32;

        // Sample to find sign change
        let delta = (t_hi - t_lo) / (BRACKET_SAMPLES as f64);
        let mut bracket: Option<(f64, f64, f64, f64)> = None;
        let mut prev: Option<(f64, f64)> = None;
        let mut best: Option<(f64, f64)> = None;

        for i in 0..=BRACKET_SAMPLES {
            let t = t_lo + delta * (i as f64);
            if let Some(f) = self.doppler_frequency_at(target_ecef, t, orbit_data, params) {
                // Track best (smallest absolute Doppler)
                if best.map_or(true, |(_, bf)| f.abs() < bf.abs()) {
                    best = Some((t, f));
                }

                // Check for sign change
                if let Some((pt, pf)) = prev {
                    if pf.signum() != f.signum() {
                        bracket = Some((pt, pf, t, f));
                        break;
                    }
                }
                prev = Some((t, f));
            }
        }

        // If we found a bracket, refine with bisection
        if let Some((mut a, mut fa, mut b, mut _fb)) = bracket {
            for iteration in 0..options.max_iterations {
                if (b - a).abs() < options.time_tolerance_s {
                    let mid = (a + b) / 2.0;
                    let f_mid = self.doppler_frequency_at(target_ecef, mid, orbit_data, params)?;
                    return Some(ZeroDopplerResult {
                        time_rel_s: mid,
                        doppler_hz: f_mid,
                        iterations: iteration + 1,
                        converged: f_mid.abs() < options.doppler_tolerance_hz,
                    });
                }

                let mid = (a + b) / 2.0;
                let f_mid = self.doppler_frequency_at(target_ecef, mid, orbit_data, params)?;

                if f_mid.abs() < options.doppler_tolerance_hz {
                    return Some(ZeroDopplerResult {
                        time_rel_s: mid,
                        doppler_hz: f_mid,
                        iterations: iteration + 1,
                        converged: true,
                    });
                }

                if f_mid.signum() == fa.signum() {
                    a = mid;
                    fa = f_mid;
                } else {
                    b = mid;
                    _fb = f_mid;
                }
            }
        }

        // Return best found only if it meets a reasonable tolerance (10x nominal).
        // Returning non-converged results with converged=true would corrupt downstream geocoding.
        best.and_then(|(t, f)| {
            let relaxed_tol = options.doppler_tolerance_hz * 10.0;
            if f.abs() < relaxed_tol {
                Some(ZeroDopplerResult {
                    time_rel_s: t,
                    doppler_hz: f,
                    iterations: options.max_iterations,
                    converged: f.abs() < options.doppler_tolerance_hz,
                })
            } else {
                log::debug!(
                    "Bisection fallback: best Doppler {:.2} Hz exceeds relaxed tolerance {:.2} Hz — returning None",
                    f, relaxed_tol
                );
                None
            }
        })
    }

    /// Convenience method: solve with default options
    pub fn solve_zero_doppler_default(
        &self,
        target_ecef: &[f64; 3],
        orbit_data: &OrbitData,
        params: &RangeDopplerParams,
    ) -> Option<f64> {
        let result = self.solve_zero_doppler(
            target_ecef,
            orbit_data,
            params,
            &ZeroDopplerOptions::default(),
        );

        // INSTRUMENTATION #5: Zero-Doppler solver result validation
        static ZERO_DOPPLER_CHECK_LOGGED: std::sync::atomic::AtomicBool =
            std::sync::atomic::AtomicBool::new(false);
        if let Some(ref res) = result {
            if res.converged
                && !ZERO_DOPPLER_CHECK_LOGGED.load(std::sync::atomic::Ordering::Relaxed)
            {
                let orbit_ref_epoch =
                    crate::types::datetime_to_utc_seconds(orbit_data.reference_time);
                let absolute_time = res.time_rel_s + orbit_ref_epoch;
                let product_start_abs = params.product_start_absolute();
                let product_end_abs = product_start_abs + params.product_duration;

                log::debug!("🔍 ZERO-DOPPLER CHECK #5:");
                log::debug!("  Solver returned time_rel: {:.6}s", res.time_rel_s);
                log::debug!("  Absolute time: {:.6}s", absolute_time);
                log::debug!(
                    "  Product window: [{:.6}, {:.6}]s",
                    product_start_abs,
                    product_end_abs
                );
                log::debug!("  Doppler frequency: {:.3} Hz", res.doppler_hz);
                log::debug!("  Iterations: {}", res.iterations);

                if absolute_time < product_start_abs || absolute_time > product_end_abs {
                    log::debug!("  ⚠️ WARNING: Zero-Doppler time outside product window!");
                } else {
                    log::debug!("  ✅ Zero-Doppler time within product window");
                }

                // Validate time is orbit-relative
                if res.time_rel_s < 0.0 || res.time_rel_s >= 3600.0 {
                    log::debug!("  ⚠️ WARNING: Zero-Doppler time_rel should be orbit-relative (0-3600s), got {}", res.time_rel_s);
                } else {
                    log::debug!("  ✅ Zero-Doppler time_rel in valid orbit-relative range");
                }

                ZERO_DOPPLER_CHECK_LOGGED.store(true, std::sync::atomic::Ordering::Relaxed);
            }
        }

        result.filter(|r| r.converged).map(|r| r.time_rel_s)
    }

    /// Convenience method: solve with seed for fast refinement
    pub fn solve_zero_doppler_seeded(
        &self,
        target_ecef: &[f64; 3],
        orbit_data: &OrbitData,
        params: &RangeDopplerParams,
        seed_time: f64,
    ) -> Option<f64> {
        let options = ZeroDopplerOptions::for_seeded_refinement(seed_time);
        self.solve_zero_doppler(target_ecef, orbit_data, params, &options)
            .filter(|r| r.converged)
            .map(|r| r.time_rel_s)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_zero_doppler_options_defaults() {
        let opts = ZeroDopplerOptions::default();
        assert_eq!(opts.max_iterations, 30);
        assert!((opts.doppler_tolerance_hz - 1e-3).abs() < 1e-9);
        assert!(opts.use_bisection_fallback);
    }

    #[test]
    fn test_zero_doppler_options_presets() {
        let tie_point = ZeroDopplerOptions::for_tie_point_grid();
        assert_eq!(tie_point.max_iterations, 20);

        let seeded = ZeroDopplerOptions::for_seeded_refinement(100.0);
        assert_eq!(seeded.max_iterations, 5);
        assert_eq!(seeded.seed_time, Some(100.0));

        let precise = ZeroDopplerOptions::high_precision();
        assert_eq!(precise.max_iterations, 50);
    }

    /// Test zero-Doppler solver convergence validation
    ///
    /// Validates that the solver:
    /// 1. Converges to zero Doppler frequency (|f_D| < tolerance)
    /// 2. Achieves convergence within max_iterations
    /// 3. Returns valid orbit-relative time within product window
    #[test]
    fn test_zero_doppler_solver_convergence() {
        use super::super::*;
        use crate::types::{OrbitData, StateVector};
        use chrono::{TimeZone, Utc};

        // Create a simple orbit: circular orbit at 700 km altitude
        let reference_time = Utc.timestamp_opt(1_600_000_000, 0).unwrap();
        let orbit_radius = 7_000_000.0; // 700 km altitude
        let angular_velocity = 0.0011; // rad/s (approximate for Sentinel-1)

        // Create 5 state vectors for cubic spline interpolation
        let mut state_vectors = Vec::new();
        for i in 0..5 {
            let t_offset = i as f64 * 10.0; // 10 seconds apart
            let time = reference_time + chrono::Duration::seconds(t_offset as i64);
            let angle = angular_velocity * t_offset;

            // Circular orbit in XY plane
            let position = [orbit_radius * angle.cos(), orbit_radius * angle.sin(), 0.0];
            let velocity = [
                -orbit_radius * angular_velocity * angle.sin(),
                orbit_radius * angular_velocity * angle.cos(),
                0.0,
            ];

            state_vectors.push(StateVector {
                time,
                position,
                velocity,
            });
        }

        let orbit_data = OrbitData {
            reference_time,
            state_vectors,
        };

        // Create target on Earth surface (directly below satellite at t=20s)
        let target_time = 20.0;
        let target_angle = angular_velocity * target_time;
        let earth_radius = 6_371_000.0; // Earth radius
        let target_ecef = [
            earth_radius * target_angle.cos(),
            earth_radius * target_angle.sin(),
            0.0,
        ];

        // Create Range-Doppler parameters
        #[allow(deprecated)]
        let params = RangeDopplerParams {
            range_pixel_spacing: 2.33,
            azimuth_pixel_spacing: 0.000588,
            slant_range_time: 0.00015,
            prf: 1700.0,
            azimuth_time_interval: 0.0018,
            wavelength: 0.055465763, // Sentinel-1 C-band
            speed_of_light: 299_792_458.0,
            orbit_ref_epoch_utc: crate::types::datetime_to_utc_seconds(reference_time),
            product_start_rel_s: 0.0,
            product_start_time_abs: crate::types::datetime_to_utc_seconds(reference_time),
            product_stop_time_abs: crate::types::datetime_to_utc_seconds(reference_time) + 50.0,
            product_duration: 50.0,
            total_azimuth_lines: Some(10000),
            doppler_centroid: None,
            first_valid_line: None,
            last_valid_line: None,
            first_valid_sample: None,
            last_valid_sample: None,
            range_multilook_factor: 1.0,
            azimuth_multilook_factor: 1.0,
            range_multilook_safe: 1.0,
            azimuth_multilook_safe: 1.0,
            subswaths: std::collections::HashMap::new(),
            burst_timings: Vec::new(),
            burst_segments: Vec::new(),
            reference_incidence_angle_deg: None,
            incidence_angle_near_deg: None,
            incidence_angle_far_deg: None,
            total_range_samples: None,
        };
        let mut params = params;
        params.compute_safe_multilook_factors();
        let params = params;

        // Create a dummy corrector (we only need the solver)
        let dem = ndarray::Array2::<f32>::zeros((2, 2));
        let dem_transform = GeoTransform {
            top_left_x: 0.0,
            pixel_width: 1.0,
            rotation_x: 0.0,
            top_left_y: 0.0,
            rotation_y: 0.0,
            pixel_height: -1.0,
        };
        let mut corrector = TerrainCorrector::new(dem, dem_transform, -32768.0, 4326, 4326, 20.0);
        corrector.set_orbit_data(orbit_data.clone());

        // Test solver convergence
        let options = ZeroDopplerOptions::default();
        let result = corrector.solve_zero_doppler(&target_ecef, &orbit_data, &params, &options);

        assert!(result.is_some(), "Solver should converge for valid target");
        let result = result.unwrap();

        // Validate convergence
        assert!(result.converged, "Solver should report convergence");
        assert!(
            result.doppler_hz.abs() < options.doppler_tolerance_hz,
            "Doppler frequency should be within tolerance: |f_D| = {:.6} Hz < {:.6} Hz",
            result.doppler_hz.abs(),
            options.doppler_tolerance_hz
        );
        assert!(
            result.iterations <= options.max_iterations,
            "Should converge within max_iterations: {} <= {}",
            result.iterations,
            options.max_iterations
        );

        // Validate time is within product window
        assert!(
            result.time_rel_s >= params.product_start_rel_s - 1.0
                && result.time_rel_s <= params.product_start_rel_s + params.product_duration + 1.0,
            "Time should be within product window: {:.3}s in [{:.3}, {:.3}]",
            result.time_rel_s,
            params.product_start_rel_s,
            params.product_start_rel_s + params.product_duration
        );
    }
}
