#![allow(dead_code, unused_variables)]
use chrono::{DateTime, Utc};

use crate::constants::physical::geodetic::{
    WGS84_ECCENTRICITY_SQUARED, WGS84_SEMI_MAJOR_AXIS_M, WGS84_SEMI_MINOR_AXIS_M,
};
use crate::types::{BurstOrbitData, OrbitData, SarError, SarResult, StateVector};
use log;

/// Convert DateTime<Utc> to f64 seconds with nanosecond precision
#[inline]
pub(crate) fn datetime_to_secs_precise(dt: DateTime<Utc>) -> f64 {
    let nanos = dt.timestamp_nanos_opt().expect("timestamp out of range") as f64;
    nanos * 1e-9
}

/// Ensure target time lies within OSV coverage
pub(crate) fn validate_orbit_coverage(
    orbit: &OrbitData,
    target_time: DateTime<Utc>,
) -> SarResult<()> {
    if orbit.state_vectors.is_empty() {
        return Err(SarError::Processing(
            "Orbit file parsed but contains no state vectors".to_string(),
        ));
    }
    let first = orbit
        .state_vectors
        .first()
        .map(|sv| sv.time)
        .expect("non-empty checked");
    let last = orbit
        .state_vectors
        .last()
        .map(|sv| sv.time)
        .expect("non-empty checked");
    if target_time < first || target_time > last {
        return Err(SarError::Processing(format!(
            "Orbit coverage mismatch: target {} outside [{}, {}]",
            target_time, first, last
        )));
    }
    Ok(())
}

pub fn interpolate_position(orbit: &OrbitData, target_time: DateTime<Utc>) -> SarResult<[f64; 3]> {
    if orbit.state_vectors.is_empty() {
        return Err(SarError::Processing(
            "No state vectors in orbit data".to_string(),
        ));
    }

    validate_orbit_coverage(orbit, target_time)?;

    let target_timestamp = datetime_to_secs_precise(target_time);
    let selected_svs = select_interp_vectors(&orbit.state_vectors, target_timestamp)?;
    lagrange_interpolate(&selected_svs, target_timestamp)
}

pub fn interpolate_velocity(orbit: &OrbitData, target_time: DateTime<Utc>) -> SarResult<[f64; 3]> {
    if orbit.state_vectors.is_empty() {
        return Err(SarError::Processing(
            "No state vectors in orbit data".to_string(),
        ));
    }

    validate_orbit_coverage(orbit, target_time)?;

    let target_timestamp = datetime_to_secs_precise(target_time);
    let selected_svs = select_interp_vectors(&orbit.state_vectors, target_timestamp)?;
    lagrange_interpolate_velocity(&selected_svs, target_timestamp)
}

/// Select up to four interpolation vectors around a target timestamp
pub(crate) fn select_interp_vectors(
    state_vectors: &[StateVector],
    target_timestamp: f64,
) -> SarResult<Vec<&StateVector>> {
    if state_vectors.is_empty() {
        return Err(SarError::Processing(
            "No state vectors available".to_string(),
        ));
    }

    let closest_idx = binary_search_closest_time(state_vectors, target_timestamp);

    let num_points = std::cmp::min(4, state_vectors.len());
    if num_points < 2 {
        log::warn!("Not enough state vectors for interpolation, using nearest");
        return Ok(vec![&state_vectors[closest_idx]]);
    }

    let start_idx = if closest_idx >= num_points / 2 {
        std::cmp::min(
            closest_idx - num_points / 2,
            state_vectors.len() - num_points,
        )
    } else {
        0
    };

    let selected_svs: Vec<&StateVector> = (start_idx..start_idx + num_points)
        .map(|i| &state_vectors[i])
        .collect();

    Ok(selected_svs)
}

/// Binary search to find the state vector closest in time to target
pub(crate) fn binary_search_closest_time(
    state_vectors: &[StateVector],
    target_timestamp: f64,
) -> usize {
    if state_vectors.is_empty() {
        return 0;
    }

    let mut left = 0;
    let mut right = state_vectors.len() - 1;

    while left < right {
        let mid = left + (right - left) / 2;
        let mid_timestamp = datetime_to_secs_precise(state_vectors[mid].time);

        if mid_timestamp < target_timestamp {
            left = mid + 1;
        } else {
            right = mid;
        }
    }

    // Check if left-1 is closer than left
    if left > 0 {
        let left_dist =
            (datetime_to_secs_precise(state_vectors[left].time) - target_timestamp).abs();
        let left_minus_1_dist =
            (datetime_to_secs_precise(state_vectors[left - 1].time) - target_timestamp).abs();

        if left_minus_1_dist < left_dist {
            return left - 1;
        }
    }

    left
}

/// Lagrange polynomial interpolation for orbit position
pub(crate) fn lagrange_interpolate(
    state_vectors: &[&StateVector],
    target_time: f64,
) -> SarResult<[f64; 3]> {
    let n = state_vectors.len();
    let mut result = [0.0; 3];

    for (i, sv_i) in state_vectors.iter().enumerate().take(n) {
        let ti = datetime_to_secs_precise(sv_i.time);
        let mut li = 1.0;

        for (j, sv_j) in state_vectors.iter().enumerate().take(n) {
            if i != j {
                let tj = datetime_to_secs_precise(sv_j.time);
                li *= (target_time - tj) / (ti - tj);
            }
        }

        for (coord, res) in result.iter_mut().enumerate() {
            *res += li * sv_i.position[coord];
        }
    }

    Ok(result)
}

/// Lagrange interpolation for velocity
pub(crate) fn lagrange_interpolate_velocity(
    state_vectors: &[&StateVector],
    target_time: f64,
) -> SarResult<[f64; 3]> {
    let n = state_vectors.len();
    let mut result = [0.0; 3];

    for (i, sv_i) in state_vectors.iter().enumerate().take(n) {
        let ti = datetime_to_secs_precise(sv_i.time);
        let mut li = 1.0;

        for (j, sv_j) in state_vectors.iter().enumerate().take(n) {
            if i != j {
                let tj = datetime_to_secs_precise(sv_j.time);
                li *= (target_time - tj) / (ti - tj);
            }
        }

        for (coord, res) in result.iter_mut().enumerate() {
            *res += li * sv_i.velocity[coord];
        }
    }

    Ok(result)
}

/// Get satellite position and velocity for a burst at specific azimuth times
pub fn interpolate_burst_orbit(
    orbit: &OrbitData,
    burst_start_time: DateTime<Utc>,
    azimuth_time_interval: f64,
    num_azimuth_lines: usize,
) -> SarResult<BurstOrbitData> {
    if orbit.state_vectors.is_empty() {
        return Err(SarError::Processing(
            "No state vectors in orbit data".to_string(),
        ));
    }

    // Validate coverage for the full burst time span before interpolating
    let burst_start_timestamp = datetime_to_secs_precise(burst_start_time);
    let burst_end_timestamp =
        burst_start_timestamp + (num_azimuth_lines as f64 * azimuth_time_interval);
    let orbit_first = datetime_to_secs_precise(
        orbit
            .state_vectors
            .first()
            .map(|sv| sv.time)
            .expect("non-empty"),
    );
    let orbit_last = datetime_to_secs_precise(
        orbit
            .state_vectors
            .last()
            .map(|sv| sv.time)
            .expect("non-empty"),
    );
    if burst_start_timestamp < orbit_first || burst_end_timestamp > orbit_last {
        return Err(SarError::Processing(format!(
            "Orbit coverage mismatch: [{:.3}, {:.3}]s outside [{:.3}, {:.3}]s",
            burst_start_timestamp, burst_end_timestamp, orbit_first, orbit_last
        )));
    }

    // Pre-sort state vectors by time for efficient binary search
    let mut sorted_state_vectors = orbit.state_vectors.clone();
    sorted_state_vectors.sort_by(|a, b| a.time.cmp(&b.time));

    let mut positions = Vec::with_capacity(num_azimuth_lines);
    let mut velocities = Vec::with_capacity(num_azimuth_lines);
    let mut attitudes = Vec::with_capacity(num_azimuth_lines);
    let mut geodetic = Vec::with_capacity(num_azimuth_lines);
    let mut azimuth_times = Vec::with_capacity(num_azimuth_lines);

    // Find the range of state vectors we'll need for this entire burst
    let start_search_idx = binary_search_closest_time(&sorted_state_vectors, burst_start_timestamp);
    let end_search_idx = binary_search_closest_time(&sorted_state_vectors, burst_end_timestamp);

    // Expand search range to ensure we have enough vectors for interpolation
    let search_start = start_search_idx.saturating_sub(2);
    let search_end = std::cmp::min(end_search_idx + 3, sorted_state_vectors.len());

    let relevant_vectors = &sorted_state_vectors[search_start..search_end];

    // Quick residual sanity check for the interpolator; fallback to a linear scheme if it blows up
    let (pos_rms, vel_rms) = compute_interp_residuals(relevant_vectors)?;
    let bad_fit = !pos_rms.is_finite() || !vel_rms.is_finite() || pos_rms > 50.0 || vel_rms > 5.0;
    if bad_fit {
        log::warn!(
            "⚠️  Orbit interpolation residuals high (pos RMS {:.3} m, vel RMS {:.3} m/s); falling back to linear interpolation.",
            pos_rms,
            vel_rms
        );
    } else {
        log::info!(
            "Orbit interpolation residuals: position RMS {:.3} m, velocity RMS {:.3} m/s",
            pos_rms,
            vel_rms
        );
    }

    // Compute azimuth times in integer nanoseconds to avoid accumulation drift
    let dt_ns_per_line = (azimuth_time_interval * 1e9).round() as i128;
    for line_idx in 0..num_azimuth_lines {
        let t_ns = (line_idx as i128) * dt_ns_per_line;
        let azimuth_time = burst_start_time + chrono::Duration::nanoseconds(t_ns as i64);

        let target_timestamp = datetime_to_secs_precise(azimuth_time);
        let (position, velocity) = if bad_fit {
            linear_interpolate(relevant_vectors, target_timestamp)?
        } else {
            let selected_svs = select_interp_vectors(relevant_vectors, target_timestamp)?;
            let pos = lagrange_interpolate(&selected_svs, target_timestamp)?;
            let vel = lagrange_interpolate_velocity(&selected_svs, target_timestamp)?;
            (pos, vel)
        };

        positions.push(position);
        velocities.push(velocity);
        attitudes.push(derive_lvlh_quaternion(position, velocity));
        geodetic.push(ecef_to_geodetic(position));
        azimuth_times.push(azimuth_time);
    }

    Ok(BurstOrbitData {
        positions,
        velocities,
        attitudes: Some(attitudes),
        geodetic: Some(geodetic),
        azimuth_times,
        burst_start_time,
        azimuth_time_interval,
    })
}

pub fn calculate_doppler_centroid(
    satellite_velocity: [f64; 3],
    look_direction: [f64; 3],
    wavelength: f64,
) -> f64 {
    let velocity_dot_look = satellite_velocity[0] * look_direction[0]
        + satellite_velocity[1] * look_direction[1]
        + satellite_velocity[2] * look_direction[2];

    2.0 * velocity_dot_look / wavelength
}

fn compute_interp_residuals(vectors: &[StateVector]) -> SarResult<(f64, f64)> {
    if vectors.len() < 2 {
        return Ok((0.0, 0.0));
    }

    let mut pos_err_sq = 0.0;
    let mut vel_err_sq = 0.0;
    let mut count = 0.0;

    for sv in vectors {
        let ts = datetime_to_secs_precise(sv.time);
        let selected = select_interp_vectors(vectors, ts)?;
        let pred_pos = lagrange_interpolate(&selected, ts)?;
        let pred_vel = lagrange_interpolate_velocity(&selected, ts)?;

        let dp = [
            pred_pos[0] - sv.position[0],
            pred_pos[1] - sv.position[1],
            pred_pos[2] - sv.position[2],
        ];
        let dv = [
            pred_vel[0] - sv.velocity[0],
            pred_vel[1] - sv.velocity[1],
            pred_vel[2] - sv.velocity[2],
        ];

        pos_err_sq += dp[0] * dp[0] + dp[1] * dp[1] + dp[2] * dp[2];
        vel_err_sq += dv[0] * dv[0] + dv[1] * dv[1] + dv[2] * dv[2];
        count += 1.0;
    }

    if count == 0.0 {
        return Ok((0.0, 0.0));
    }

    Ok(((pos_err_sq / count).sqrt(), (vel_err_sq / count).sqrt()))
}

fn linear_interpolate(vectors: &[StateVector], target_ts: f64) -> SarResult<([f64; 3], [f64; 3])> {
    if vectors.is_empty() {
        return Err(SarError::Processing(
            "No state vectors available for linear interpolation".to_string(),
        ));
    }

    // Find bracketing vectors
    let idx = binary_search_closest_time(vectors, target_ts);
    let (a, b) = if idx + 1 < vectors.len() {
        (idx, idx + 1)
    } else if idx > 0 {
        (idx - 1, idx)
    } else {
        (idx, idx)
    };

    let ta = datetime_to_secs_precise(vectors[a].time);
    let tb = datetime_to_secs_precise(vectors[b].time);
    let alpha = if tb > ta {
        ((target_ts - ta) / (tb - ta)).clamp(0.0, 1.0)
    } else {
        0.0
    };

    let lerp = |a: f64, b: f64| a + alpha * (b - a);

    let mut pos = [0.0; 3];
    let mut vel = [0.0; 3];
    for i in 0..3 {
        pos[i] = lerp(vectors[a].position[i], vectors[b].position[i]);
        vel[i] = lerp(vectors[a].velocity[i], vectors[b].velocity[i]);
    }

    Ok((pos, vel))
}

fn derive_lvlh_quaternion(position_ecef: [f64; 3], velocity_ecef: [f64; 3]) -> [f64; 4] {
    // Build an LVLH-style frame from position/velocity and convert to quaternion (ECEF->LVLH)
    let r = position_ecef;
    let v = velocity_ecef;

    let norm = |v: [f64; 3]| -> f64 { (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt() };
    let cross = |a: [f64; 3], b: [f64; 3]| -> [f64; 3] {
        [
            a[1] * b[2] - a[2] * b[1],
            a[2] * b[0] - a[0] * b[2],
            a[0] * b[1] - a[1] * b[0],
        ]
    };

    let mut z_b = [-r[0], -r[1], -r[2]]; // toward Earth center
    let z_norm = norm(z_b);
    if z_norm > 0.0 {
        z_b.iter_mut().for_each(|c| *c /= z_norm);
    }

    let h = cross(r, v);
    let mut y_b = [-h[0], -h[1], -h[2]]; // along -orbit normal
    let y_norm = norm(y_b);
    if y_norm > 0.0 {
        y_b.iter_mut().for_each(|c| *c /= y_norm);
    }

    let mut x_b = cross(y_b, z_b);
    let x_norm = norm(x_b);
    if x_norm > 0.0 {
        x_b.iter_mut().for_each(|c| *c /= x_norm);
    }

    // Rotation matrix from ECEF to LVLH (columns are basis vectors)
    let rmat = [
        [x_b[0], y_b[0], z_b[0]],
        [x_b[1], y_b[1], z_b[1]],
        [x_b[2], y_b[2], z_b[2]],
    ];

    // Convert rotation matrix to quaternion (w, x, y, z)
    let trace = rmat[0][0] + rmat[1][1] + rmat[2][2];
    let (w, x, y, z) = if trace > 0.0 {
        let s = (trace + 1.0).sqrt() * 2.0;
        let w = 0.25 * s;
        let x = (rmat[2][1] - rmat[1][2]) / s;
        let y = (rmat[0][2] - rmat[2][0]) / s;
        let z = (rmat[1][0] - rmat[0][1]) / s;
        (w, x, y, z)
    } else if rmat[0][0] > rmat[1][1] && rmat[0][0] > rmat[2][2] {
        let s = (1.0 + rmat[0][0] - rmat[1][1] - rmat[2][2]).sqrt() * 2.0;
        let w = (rmat[2][1] - rmat[1][2]) / s;
        let x = 0.25 * s;
        let y = (rmat[0][1] + rmat[1][0]) / s;
        let z = (rmat[0][2] + rmat[2][0]) / s;
        (w, x, y, z)
    } else if rmat[1][1] > rmat[2][2] {
        let s = (1.0 + rmat[1][1] - rmat[0][0] - rmat[2][2]).sqrt() * 2.0;
        let w = (rmat[0][2] - rmat[2][0]) / s;
        let x = (rmat[0][1] + rmat[1][0]) / s;
        let y = 0.25 * s;
        let z = (rmat[1][2] + rmat[2][1]) / s;
        (w, x, y, z)
    } else {
        let s = (1.0 + rmat[2][2] - rmat[0][0] - rmat[1][1]).sqrt() * 2.0;
        let w = (rmat[1][0] - rmat[0][1]) / s;
        let x = (rmat[0][2] + rmat[2][0]) / s;
        let y = (rmat[1][2] + rmat[2][1]) / s;
        let z = 0.25 * s;
        (w, x, y, z)
    };

    let norm_q = (w * w + x * x + y * y + z * z).sqrt();
    if norm_q > 0.0 {
        [w / norm_q, x / norm_q, y / norm_q, z / norm_q]
    } else {
        [1.0, 0.0, 0.0, 0.0]
    }
}

fn ecef_to_geodetic(ecef: [f64; 3]) -> [f64; 3] {
    // Returns [lat(rad), lon(rad), height(m)]
    let x = ecef[0];
    let y = ecef[1];
    let z = ecef[2];

    let lon = y.atan2(x);
    let a = WGS84_SEMI_MAJOR_AXIS_M;
    let b = WGS84_SEMI_MINOR_AXIS_M;
    let e2 = WGS84_ECCENTRICITY_SQUARED;

    let p = (x * x + y * y).sqrt();
    let mut lat = z.atan2(p * (1.0 - e2));
    for _ in 0..5 {
        let sin_lat = lat.sin();
        let n = a / (1.0 - e2 * sin_lat * sin_lat).sqrt();
        lat = (z + e2 * n * sin_lat).atan2(p);
    }
    let sin_lat = lat.sin();
    let n = a / (1.0 - e2 * sin_lat * sin_lat).sqrt();
    let height = p / lat.cos() - n;

    [lat, lon, height]
}
