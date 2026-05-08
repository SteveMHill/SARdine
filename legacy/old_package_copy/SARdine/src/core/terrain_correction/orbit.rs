#![allow(dead_code, unused_variables)]
use std::collections::HashMap;

use crate::types::{OrbitData, SarError, SarResult};

use super::types::{Position3D, Vector3, Velocity3D};

/// High-performance orbit interpolation cache
///
/// # Deprecation notice
/// This struct is superseded by [`OrbitSplineCache`], which uses Hermite cubic
/// interpolation and is the actively maintained interpolation path.
/// `OrbitCache` exists only as a placeholder for a future SoA SIMD layout and
/// should not be used in production code.
#[derive(Debug)]
pub struct OrbitCache {
    /// Pre-computed interpolation coefficients
    interpolation_cache: HashMap<(i64, i64), (Vector3, Vector3)>, // (time_key, precision) -> (position, velocity)
    /// Lookup table for time to index mapping
    time_index_lut: Vec<(f64, usize)>,
    /// SIMD-optimized X-position buffer: each [f64;4] holds X of 4 consecutive state vectors
    positions_x_simd: Vec<[f64; 4]>,
    /// SIMD-optimized Y-position buffer
    positions_y_simd: Vec<[f64; 4]>,
    /// SIMD-optimized Z-position buffer
    positions_z_simd: Vec<[f64; 4]>,
    /// SIMD-optimized X-velocity buffer
    velocities_x_simd: Vec<[f64; 4]>,
    /// SIMD-optimized Y-velocity buffer
    velocities_y_simd: Vec<[f64; 4]>,
    /// SIMD-optimized Z-velocity buffer
    velocities_z_simd: Vec<[f64; 4]>,
}

impl OrbitCache {
    pub fn new(orbit_data: &OrbitData) -> Self {
        let mut cache = Self {
            interpolation_cache: HashMap::new(),
            time_index_lut: Vec::new(),
            positions_x_simd: Vec::new(),
            positions_y_simd: Vec::new(),
            positions_z_simd: Vec::new(),
            velocities_x_simd: Vec::new(),
            velocities_y_simd: Vec::new(),
            velocities_z_simd: Vec::new(),
        };
        cache.precompute_orbit_data(orbit_data);
        cache
    }

    fn precompute_orbit_data(&mut self, orbit_data: &OrbitData) {
        // Build time index lookup table
        for (i, sv) in orbit_data.state_vectors.iter().enumerate() {
            let time_seconds =
                sv.time.timestamp() as f64 + sv.time.timestamp_subsec_nanos() as f64 * 1e-9;
            self.time_index_lut.push((time_seconds, i));
        }
        self.time_index_lut.sort_by(|a, b| a.0.total_cmp(&b.0));

        // Pre-pack data for SIMD operations using structure-of-arrays layout.
        // Each [f64; 4] chunk holds one coordinate of 4 consecutive state vectors.
        for chunk in orbit_data.state_vectors.chunks(4) {
            let mut px = [0.0; 4];
            let mut py = [0.0; 4];
            let mut pz = [0.0; 4];
            let mut vx = [0.0; 4];
            let mut vy = [0.0; 4];
            let mut vz = [0.0; 4];

            for (i, sv) in chunk.iter().enumerate() {
                px[i] = sv.position[0];
                py[i] = sv.position[1];
                pz[i] = sv.position[2];
                vx[i] = sv.velocity[0];
                vy[i] = sv.velocity[1];
                vz[i] = sv.velocity[2];
            }
            self.positions_x_simd.push(px);
            self.positions_y_simd.push(py);
            self.positions_z_simd.push(pz);
            self.velocities_x_simd.push(vx);
            self.velocities_y_simd.push(vy);
            self.velocities_z_simd.push(vz);
        }
    }
}

/// Phase 3.2A: Orbit spline cache using Hermite cubic interpolation
/// Precomputes spline coefficients once, enabling ~100 ns orbit queries (vs 2-5 μs)
pub struct OrbitSplineCache {
    /// Knot times (Unix seconds absolute)
    times: Vec<f64>,

    /// Position coefficients: flat array of [a,b,c,d] for each axis and segment
    /// Layout: segment 0 [x_a, x_b, x_c, x_d, y_a, y_b, y_c, y_d, z_a, z_b, z_c, z_d],
    ///         segment 1 [...], ...
    /// Total: (n_state_vectors - 1) * 12 coefficients
    pos_coeffs: Vec<f64>,

    /// Velocity coefficients: flat array of [a,b,c,d] for each axis and segment
    /// Layout: same as pos_coeffs
    vel_coeffs: Vec<f64>,
}

impl OrbitSplineCache {
    /// Build spline cache from orbit data
    /// Uses Hermite cubic splines with position and velocity constraints at knots
    pub fn from_orbit_data(orbit_data: &OrbitData) -> SarResult<Self> {
        let n = orbit_data.state_vectors.len();

        if n < 2 {
            return Err(SarError::Processing(
                "Need at least 2 orbit state vectors for spline interpolation".to_string(),
            ));
        }

        let mut times = Vec::with_capacity(n);

        // Extract times and positions/velocities
        for sv in &orbit_data.state_vectors {
            times.push(crate::types::datetime_to_utc_seconds(sv.time));
        }

        // Fit Hermite cubics for each segment
        let n_segments = n - 1;
        let mut pos_coeffs = Vec::with_capacity(n_segments * 12);
        let mut vel_coeffs = Vec::with_capacity(n_segments * 12);

        for i in 0..n_segments {
            let sv0 = &orbit_data.state_vectors[i];
            let sv1 = &orbit_data.state_vectors[i + 1];
            let h = times[i + 1] - times[i];

            // Hermite cubic for each axis (x, y, z)
            // Given: p(t0) = p0, p(t1) = p1, p'(t0) = v0, p'(t1) = v1
            // Solve: p(t) = a + b*(t-t0) + c*(t-t0)² + d*(t-t0)³

            for axis in 0..3 {
                let p0 = sv0.position[axis];
                let p1 = sv1.position[axis];
                let v0 = sv0.velocity[axis];
                let v1 = sv1.velocity[axis];

                // Hermite basis coefficients
                let a = p0;
                let b = v0;
                let c = (3.0 * (p1 - p0) / h - 2.0 * v0 - v1) / h;
                let d = (2.0 * (p0 - p1) / h + v0 + v1) / (h * h);

                pos_coeffs.extend_from_slice(&[a, b, c, d]);

                // Velocity is derivative: v(t) = b + 2c*(t-t0) + 3d*(t-t0)²
                // For Hermite, velocity also uses cubic (or we can use linear approximation)
                // Here we use the derivative of position spline
                vel_coeffs.extend_from_slice(&[b, 2.0 * c, 3.0 * d, 0.0]);
            }
        }

        log::info!(
            "Phase 3.2: Orbit spline cache built ({} state vectors, {} segments, {} coeffs)",
            n,
            n_segments,
            pos_coeffs.len()
        );

        Ok(OrbitSplineCache {
            times,
            pos_coeffs,
            vel_coeffs,
        })
    }

    /// Fast orbit interpolation using precomputed splines (~100 ns vs 2-5 μs)
    #[inline(always)]
    pub fn interpolate(&self, t_abs: f64) -> SarResult<(Position3D, Velocity3D)> {
        // Binary search for segment containing t_abs
        let idx = match self.times.binary_search_by(|probe| {
            probe
                .partial_cmp(&t_abs)
                .unwrap_or(std::cmp::Ordering::Equal)
        }) {
            Ok(i) => {
                // Exact match - use this segment (or previous if last knot)
                if i >= self.times.len() - 1 {
                    i - 1
                } else {
                    i
                }
            }
            Err(i) => {
                if i == 0 || i >= self.times.len() {
                    return Err(SarError::Processing(format!(
                        "Time {} outside orbit range [{:.6}, {:.6}]",
                        t_abs,
                        self.times[0],
                        self.times[self.times.len() - 1]
                    )));
                }
                i - 1
            }
        };

        let t0 = self.times[idx];
        let dt = t_abs - t0;

        // Evaluate position: p(t) = a + b*dt + c*dt² + d*dt³ (Horner's method)
        // Evaluate velocity: v(t) = b + 2c*dt + 3d*dt²

        let base = idx * 12; // 12 coeffs per segment (3 axes × 4 coeffs)

        let mut pos = [0.0; 3];
        let mut vel = [0.0; 3];

        for axis in 0..3 {
            let offset = base + axis * 4;
            let a = self.pos_coeffs[offset];
            let b = self.pos_coeffs[offset + 1];
            let c = self.pos_coeffs[offset + 2];
            let d = self.pos_coeffs[offset + 3];

            // Horner's method for position
            pos[axis] = a + dt * (b + dt * (c + dt * d));

            // Derivative for velocity
            vel[axis] = b + dt * (2.0 * c + dt * 3.0 * d);
        }

        Ok((
            Position3D {
                x: pos[0],
                y: pos[1],
                z: pos[2],
            },
            Velocity3D {
                x: vel[0],
                y: vel[1],
                z: vel[2],
            },
        ))
    }
}
