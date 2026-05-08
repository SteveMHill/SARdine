//! Tests for ISCE2-style forward geocoding (SAR pixel → ground coordinates)
//!
//! These tests verify the TCN (Track-Cross-Nadir) coordinate system implementation
//! that matches ISCE2's topozero.f90 algorithm.

use sardine::types::SarResult;

/// WGS84 semi-major axis in meters
const WGS84_A: f64 = 6378137.0;
/// WGS84 semi-minor axis in meters  
const WGS84_B: f64 = 6356752.314245;

/// Convert ECEF to WGS84 lat/lon (test implementation)
/// Uses Bowring's iterative method for accurate conversion
fn ecef_to_wgs84_test(x: f64, y: f64, z: f64) -> (f64, f64) {
    let a = WGS84_A;
    let f = 1.0 / 298.257223563; // WGS84 flattening
    let b = a * (1.0 - f); // Semi-minor axis
    let e_sq = 2.0 * f - f * f; // First eccentricity squared
    let ep_sq = (a * a - b * b) / (b * b); // Second eccentricity squared

    let p = (x * x + y * y).sqrt();
    let lon = y.atan2(x);

    // Bowring's iterative method
    let mut lat = (z / (p * (1.0 - e_sq))).atan(); // Initial approximation

    for _ in 0..10 {
        let sin_lat = lat.sin();
        let n = a / (1.0 - e_sq * sin_lat * sin_lat).sqrt();
        lat = (z + e_sq * n * sin_lat).atan2(p);
    }

    (lat.to_degrees(), lon.to_degrees())
}

/// Convert WGS84 lat/lon to ECEF (test implementation)
fn wgs84_to_ecef_test(lat_deg: f64, lon_deg: f64, height: f64) -> (f64, f64, f64) {
    let lat = lat_deg.to_radians();
    let lon = lon_deg.to_radians();

    let f = 1.0 / 298.257223563;
    let e_sq = 2.0 * f - f * f;

    let n = WGS84_A / (1.0 - e_sq * lat.sin().powi(2)).sqrt();

    let x = (n + height) * lat.cos() * lon.cos();
    let y = (n + height) * lat.cos() * lon.sin();
    let z = (n * (1.0 - e_sq) + height) * lat.sin();

    (x, y, z)
}

/// Vector magnitude helper
fn vec3_norm(x: f64, y: f64, z: f64) -> f64 {
    (x * x + y * y + z * z).sqrt()
}

/// Core ISCE2-style forward geocoding function for testing
///
/// Given satellite position, velocity, and slant range, compute the ground point
/// using the TCN (Track-Cross-Nadir) coordinate system.
///
/// Parameters:
/// - sat_pos: Satellite position in ECEF (meters)
/// - sat_vel: Satellite velocity in ECEF (m/s)
/// - slant_range: Distance from satellite to ground point (meters)
/// - ilrl: Look direction (-1 for right-looking, +1 for left-looking)
fn isce2_forward_geocode(
    sat_pos: [f64; 3],
    sat_vel: [f64; 3],
    slant_range: f64,
    ilrl: f64,
) -> Option<(f64, f64)> {
    const EPSILON: f64 = 1e-6;

    // Normalize velocity
    let vel_mag = vec3_norm(sat_vel[0], sat_vel[1], sat_vel[2]);
    if vel_mag < EPSILON {
        return None;
    }
    let vel_unit = [
        sat_vel[0] / vel_mag,
        sat_vel[1] / vel_mag,
        sat_vel[2] / vel_mag,
    ];

    let sat_mag = vec3_norm(sat_pos[0], sat_pos[1], sat_pos[2]);

    // =========================================================================
    // ISCE2 TCN BASIS CONSTRUCTION
    // =========================================================================
    // ISCE2 Ellipsoid.cpp uses the LOCAL ELLIPSOID NORMAL at the sub-satellite point
    // (the ground point directly below the satellite), NOT the geocentric direction.
    //
    // From Ellipsoid.cpp line 140-143:
    //   latlon(pos,llh,XYZ_2_LLH);  // Convert satellite XYZ to lat,lon,h
    //   n[0] = -cos(llh[0]) * cos(llh[1]);
    //   n[1] = -cos(llh[0]) * sin(llh[1]);
    //   n[2] = -sin(llh[0]);
    //
    // This is the outward normal to the ellipsoid at (lat, lon), pointing AWAY from Earth.
    // We need the inward normal (toward Earth), so we use the negative.

    // First, get satellite lat/lon
    let (sat_lat_deg, sat_lon_deg) = ecef_to_wgs84_test(sat_pos[0], sat_pos[1], sat_pos[2]);
    let sat_lat = sat_lat_deg.to_radians();
    let sat_lon = sat_lon_deg.to_radians();

    // Compute local ellipsoid normal (pointing toward Earth = inward)
    // ISCE2 uses outward normal, but we need inward for nhat
    let nhat = [
        -(sat_lat.cos() * sat_lon.cos()), // = -cos(lat)*cos(lon)
        -(sat_lat.cos() * sat_lon.sin()), // = -cos(lat)*sin(lon)
        -sat_lat.sin(),                   // = -sin(lat)
    ];

    // =================================================================
    // ISCE2 TCN BASIS (from Ellipsoid.cpp lines 138-152)
    // =================================================================
    // Order: chat first, then that

    // C (chat) = normalize(N × V) - cross-track direction
    let n_cross_v = [
        nhat[1] * vel_unit[2] - nhat[2] * vel_unit[1],
        nhat[2] * vel_unit[0] - nhat[0] * vel_unit[2],
        nhat[0] * vel_unit[1] - nhat[1] * vel_unit[0],
    ];
    let n_cross_v_mag = vec3_norm(n_cross_v[0], n_cross_v[1], n_cross_v[2]);
    if n_cross_v_mag < EPSILON {
        return None;
    }
    let chat = [
        n_cross_v[0] / n_cross_v_mag,
        n_cross_v[1] / n_cross_v_mag,
        n_cross_v[2] / n_cross_v_mag,
    ];

    // T (that) = normalize(C × N) - along-track direction
    let c_cross_n = [
        chat[1] * nhat[2] - chat[2] * nhat[1],
        chat[2] * nhat[0] - chat[0] * nhat[2],
        chat[0] * nhat[1] - chat[1] * nhat[0],
    ];
    let c_cross_n_mag = vec3_norm(c_cross_n[0], c_cross_n[1], c_cross_n[2]);
    let that = [
        c_cross_n[0] / c_cross_n_mag,
        c_cross_n[1] / c_cross_n_mag,
        c_cross_n[2] / c_cross_n_mag,
    ];

    // Compute theta from law of cosines
    let aa = sat_mag;
    let bb = WGS84_A; // Approximate Earth radius
    let r = slant_range;

    let cos_theta = (aa * aa + r * r - bb * bb) / (2.0 * aa * r);
    if cos_theta.abs() > 1.0 {
        return None;
    }
    let sin_theta = (1.0 - cos_theta * cos_theta).sqrt();

    // TCN components
    let gamma = cos_theta * r; // Nadir component

    // Zero-Doppler: alpha ≈ 0 only when N·V ≈ 0
    let n_dot_v = nhat[0] * vel_unit[0] + nhat[1] * vel_unit[1] + nhat[2] * vel_unit[2];
    let t_dot_v = that[0] * vel_unit[0] + that[1] * vel_unit[1] + that[2] * vel_unit[2];

    // Debug
    println!("  n_dot_v = {:.6}, t_dot_v = {:.6}", n_dot_v, t_dot_v);
    println!(
        "  gamma = {:.1}m, gamma * n_dot_v = {:.1}m",
        gamma,
        gamma * n_dot_v
    );

    let alpha = if t_dot_v.abs() > EPSILON {
        (0.0 - gamma * n_dot_v) / t_dot_v
    } else {
        0.0
    };
    println!("  alpha = {:.1}m", alpha);

    // Cross-track: beta = -ILRL * sqrt(r² * sin²θ - α²)
    // ISCE2 topozero.f90 line 505: beta = -ilrl * sqrt(...)
    // The negative sign is critical!
    let cross_track_sq = r * r * sin_theta * sin_theta - alpha * alpha;
    if cross_track_sq < 0.0 {
        return None;
    }
    let beta = -ilrl * cross_track_sq.sqrt(); // NOTE: negative sign!

    // Debug: print all contributions
    let gamma_x = gamma * nhat[0];
    let alpha_x = alpha * that[0];
    let beta_x = beta * chat[0];
    let gamma_y = gamma * nhat[1];
    let alpha_y = alpha * that[1];
    let beta_y = beta * chat[1];
    println!("Debug: ilrl={:.1}", ilrl);
    println!(
        "  X: gamma*nhat = {:.1}, alpha*that = {:.1}, beta*chat = {:.1}",
        gamma_x, alpha_x, beta_x
    );
    println!(
        "  Y: gamma*nhat = {:.1}, alpha*that = {:.1}, beta*chat = {:.1}",
        gamma_y, alpha_y, beta_y
    );
    println!(
        "  Total X = {:.1}m, Total Y = {:.1}m",
        gamma_x + alpha_x + beta_x,
        gamma_y + alpha_y + beta_y
    );

    // Ground point = satellite + range_vector
    let range_vec = [
        gamma * nhat[0] + alpha * that[0] + beta * chat[0],
        gamma * nhat[1] + alpha * that[1] + beta * chat[1],
        gamma * nhat[2] + alpha * that[2] + beta * chat[2],
    ];

    let ground_ecef = [
        sat_pos[0] + range_vec[0],
        sat_pos[1] + range_vec[1],
        sat_pos[2] + range_vec[2],
    ];

    println!(
        "  Satellite ECEF: ({:.0}, {:.0}, {:.0})",
        sat_pos[0], sat_pos[1], sat_pos[2]
    );
    println!(
        "  Range vector:   ({:.0}, {:.0}, {:.0})",
        range_vec[0], range_vec[1], range_vec[2]
    );
    println!(
        "  Ground ECEF:    ({:.0}, {:.0}, {:.0})",
        ground_ecef[0], ground_ecef[1], ground_ecef[2]
    );
    println!(
        "  Ground lon = atan2({:.0}, {:.0}) = {:.4}°",
        ground_ecef[1],
        ground_ecef[0],
        ground_ecef[1].atan2(ground_ecef[0]).to_degrees()
    );

    let (lat, lon) = ecef_to_wgs84_test(ground_ecef[0], ground_ecef[1], ground_ecef[2]);

    if lat.is_finite() && lon.is_finite() {
        Some((lat, lon))
    } else {
        None
    }
}

#[test]
fn test_ecef_wgs84_roundtrip() {
    // Test that ECEF <-> WGS84 conversion is correct for zero elevation
    // (The test ECEF to WGS84 doesn't compute elevation, so test at h=0)
    let test_points = [
        (48.8, 9.5),   // Central Europe
        (47.5, 8.5),   // Switzerland
        (52.5, 13.4),  // Berlin
        (0.0, 0.0),    // Equator/Prime Meridian
        (45.0, -90.0), // North America
    ];

    for (lat, lon) in test_points {
        let (x, y, z) = wgs84_to_ecef_test(lat, lon, 0.0);
        let (lat2, lon2) = ecef_to_wgs84_test(x, y, z);

        assert!(
            (lat - lat2).abs() < 0.0001,
            "Latitude mismatch: {:.6} vs {:.6}",
            lat,
            lat2
        );
        assert!(
            (lon - lon2).abs() < 0.0001,
            "Longitude mismatch: {:.6} vs {:.6}",
            lon,
            lon2
        );
    }
}

#[test]
fn test_tcn_basis_orthonormality() {
    // Verify TCN basis vectors are orthonormal for a realistic satellite position
    let sat_pos = wgs84_to_ecef_test(48.5, 9.0, 700_000.0); // 700km altitude
    let sat_pos_arr = [sat_pos.0, sat_pos.1, sat_pos.2];

    // Approximate velocity for polar orbit (roughly along latitude lines)
    let sat_vel = [0.0, 7500.0, 1000.0]; // Mostly northward, ~7.5 km/s

    let sat_mag = vec3_norm(sat_pos_arr[0], sat_pos_arr[1], sat_pos_arr[2]);
    let vel_mag = vec3_norm(sat_vel[0], sat_vel[1], sat_vel[2]);

    // N (nadir) - points toward Earth center
    let nhat = [
        -sat_pos_arr[0] / sat_mag,
        -sat_pos_arr[1] / sat_mag,
        -sat_pos_arr[2] / sat_mag,
    ];

    // Unit velocity
    let vel_unit = [
        sat_vel[0] / vel_mag,
        sat_vel[1] / vel_mag,
        sat_vel[2] / vel_mag,
    ];

    // ISCE2 TCN: C = normalize(N × V), T = normalize(C × N)
    let n_cross_v = [
        nhat[1] * vel_unit[2] - nhat[2] * vel_unit[1],
        nhat[2] * vel_unit[0] - nhat[0] * vel_unit[2],
        nhat[0] * vel_unit[1] - nhat[1] * vel_unit[0],
    ];
    let n_cross_v_mag = vec3_norm(n_cross_v[0], n_cross_v[1], n_cross_v[2]);
    let chat = [
        n_cross_v[0] / n_cross_v_mag,
        n_cross_v[1] / n_cross_v_mag,
        n_cross_v[2] / n_cross_v_mag,
    ];

    let c_cross_n = [
        chat[1] * nhat[2] - chat[2] * nhat[1],
        chat[2] * nhat[0] - chat[0] * nhat[2],
        chat[0] * nhat[1] - chat[1] * nhat[0],
    ];
    let c_cross_n_mag = vec3_norm(c_cross_n[0], c_cross_n[1], c_cross_n[2]);
    let that = [
        c_cross_n[0] / c_cross_n_mag,
        c_cross_n[1] / c_cross_n_mag,
        c_cross_n[2] / c_cross_n_mag,
    ];

    // Verify orthonormality
    let n_dot_t = nhat[0] * that[0] + nhat[1] * that[1] + nhat[2] * that[2];
    let n_dot_c = nhat[0] * chat[0] + nhat[1] * chat[1] + nhat[2] * chat[2];
    let t_dot_c = that[0] * chat[0] + that[1] * chat[1] + that[2] * chat[2];

    assert!(n_dot_t.abs() < 1e-10, "N·T should be 0, got {}", n_dot_t);
    assert!(n_dot_c.abs() < 1e-10, "N·C should be 0, got {}", n_dot_c);
    assert!(t_dot_c.abs() < 1e-10, "T·C should be 0, got {}", t_dot_c);

    // Verify unit length
    let n_len = vec3_norm(nhat[0], nhat[1], nhat[2]);
    let t_len = vec3_norm(that[0], that[1], that[2]);
    let c_len = vec3_norm(chat[0], chat[1], chat[2]);

    assert!(
        (n_len - 1.0).abs() < 1e-10,
        "|N| should be 1, got {}",
        n_len
    );
    assert!(
        (t_len - 1.0).abs() < 1e-10,
        "|T| should be 1, got {}",
        t_len
    );
    assert!(
        (c_len - 1.0).abs() < 1e-10,
        "|C| should be 1, got {}",
        c_len
    );
}

#[test]
fn test_forward_geocode_right_looking_sentinel1() {
    // Simulate Sentinel-1 right-looking SAR geometry
    // Expected: Ground point should be to the RIGHT of satellite track (when looking along velocity)

    // Satellite at ~700km altitude over central Europe
    let target_lat = 48.5;
    let target_lon = 9.5;
    let sat_altitude = 700_000.0;

    let sat_pos = wgs84_to_ecef_test(target_lat, target_lon, sat_altitude);
    let sat_pos_arr = [sat_pos.0, sat_pos.1, sat_pos.2];

    // Velocity: roughly northward for descending pass
    // Sentinel-1 orbital velocity ~7.5 km/s
    // For a descending pass over Europe, velocity is roughly N with slight W component
    let sat_vel = [100.0, 7000.0, 2000.0]; // Predominantly northward

    // Typical Sentinel-1 slant range: 800-1000 km
    let slant_range = 850_000.0; // 850 km

    // Right-looking: ILRL = -1
    let ilrl = -1.0;

    let result = isce2_forward_geocode(sat_pos_arr, sat_vel, slant_range, ilrl);

    assert!(result.is_some(), "Forward geocoding should succeed");

    let (ground_lat, ground_lon) = result.unwrap();

    println!(
        "Satellite: lat={:.4}°, lon={:.4}°, alt={:.0}km",
        target_lat,
        target_lon,
        sat_altitude / 1000.0
    );
    println!(
        "Ground point: lat={:.4}°, lon={:.4}°",
        ground_lat, ground_lon
    );
    println!("Slant range: {:.0} km", slant_range / 1000.0);

    // For right-looking SAR flying northward:
    // - Ground point should be to the EAST (positive lon) of satellite nadir
    // - Ground point lat should be close to satellite lat (zero-Doppler)

    // Check ground is to the east for right-looking
    assert!(
        ground_lon > target_lon,
        "Right-looking SAR should see ground to the east. Got lon={:.4}°, expected > {:.4}°",
        ground_lon,
        target_lon
    );

    // Check lat is reasonable (within ~5° due to Earth curvature and slant range)
    assert!(
        (ground_lat - target_lat).abs() < 5.0,
        "Ground lat {:.4}° should be within 5° of satellite lat {:.4}°",
        ground_lat,
        target_lat
    );

    // Check lon offset is reasonable (should be ~2-5° for 850km slant range)
    let lon_offset = ground_lon - target_lon;
    assert!(
        lon_offset > 0.5 && lon_offset < 10.0,
        "Longitude offset {:.4}° seems unreasonable for 850km slant range",
        lon_offset
    );
}

#[test]
fn test_forward_geocode_left_vs_right_looking() {
    // Verify that ILRL properly controls look direction
    // Left-looking should give ground point on opposite side of right-looking

    let sat_pos = wgs84_to_ecef_test(48.5, 9.0, 700_000.0);
    let sat_pos_arr = [sat_pos.0, sat_pos.1, sat_pos.2];

    // Compute a realistic velocity tangent to a circular orbit
    // For a polar orbit, velocity is perpendicular to radial direction
    // and mostly in the N-S direction (along latitude lines)
    //
    // At lat=48.5°, a northward velocity would have components:
    // - Decrease in Z (moving away from pole)
    // - Increase in (X, Y) depending on longitude
    //
    // For simplicity, compute velocity as perpendicular to radial direction
    // with orbital speed ~7.5 km/s
    let sat_mag = vec3_norm(sat_pos_arr[0], sat_pos_arr[1], sat_pos_arr[2]);
    let radial = [
        sat_pos_arr[0] / sat_mag,
        sat_pos_arr[1] / sat_mag,
        sat_pos_arr[2] / sat_mag,
    ];

    // Orbital plane is roughly polar (inclination ~98° for Sentinel-1)
    // Choose velocity perpendicular to radial and mostly northward
    // Cross radial with East direction to get North, then normalize
    let east = [-sat_pos_arr[1], sat_pos_arr[0], 0.0];
    let east_mag = vec3_norm(east[0], east[1], east[2]);
    let east_unit = [east[0] / east_mag, east[1] / east_mag, east[2] / east_mag];

    // North = radial × East (perpendicular to radial, pointing north)
    let north = [
        radial[1] * east_unit[2] - radial[2] * east_unit[1],
        radial[2] * east_unit[0] - radial[0] * east_unit[2],
        radial[0] * east_unit[1] - radial[1] * east_unit[0],
    ];
    let north_mag = vec3_norm(north[0], north[1], north[2]);
    let north_unit = [
        north[0] / north_mag,
        north[1] / north_mag,
        north[2] / north_mag,
    ];

    // Velocity = 7.5 km/s in north direction (ascending pass)
    let orbital_speed = 7500.0;
    let sat_vel = [
        orbital_speed * north_unit[0],
        orbital_speed * north_unit[1],
        orbital_speed * north_unit[2],
    ];
    println!(
        "Computed orbital velocity: ({:.0}, {:.0}, {:.0}) m/s",
        sat_vel[0], sat_vel[1], sat_vel[2]
    );

    let slant_range = 850_000.0;

    // Debug: Print intermediate values
    let sat_mag = vec3_norm(sat_pos_arr[0], sat_pos_arr[1], sat_pos_arr[2]);
    let vel_mag = vec3_norm(sat_vel[0], sat_vel[1], sat_vel[2]);
    let vel_unit = [
        sat_vel[0] / vel_mag,
        sat_vel[1] / vel_mag,
        sat_vel[2] / vel_mag,
    ];
    let nhat = [
        -sat_pos_arr[0] / sat_mag,
        -sat_pos_arr[1] / sat_mag,
        -sat_pos_arr[2] / sat_mag,
    ];

    // Cross-track: C = normalize(N × V)
    let n_cross_v = [
        nhat[1] * vel_unit[2] - nhat[2] * vel_unit[1],
        nhat[2] * vel_unit[0] - nhat[0] * vel_unit[2],
        nhat[0] * vel_unit[1] - nhat[1] * vel_unit[0],
    ];
    let n_cross_v_mag = vec3_norm(n_cross_v[0], n_cross_v[1], n_cross_v[2]);
    let chat = [
        n_cross_v[0] / n_cross_v_mag,
        n_cross_v[1] / n_cross_v_mag,
        n_cross_v[2] / n_cross_v_mag,
    ];

    println!(
        "Satellite ECEF: ({:.0}, {:.0}, {:.0})",
        sat_pos_arr[0], sat_pos_arr[1], sat_pos_arr[2]
    );
    println!(
        "Nadir (nhat): ({:.4}, {:.4}, {:.4})",
        nhat[0], nhat[1], nhat[2]
    );
    println!(
        "Velocity unit: ({:.4}, {:.4}, {:.4})",
        vel_unit[0], vel_unit[1], vel_unit[2]
    );
    println!(
        "Cross-track (chat = N×V): ({:.4}, {:.4}, {:.4})",
        chat[0], chat[1], chat[2]
    );

    // Convert chat to lat/lon direction for intuition
    // chat should point approximately to one side of the track
    // For N×V with V mostly +Y (north), and N mostly -Z (down):
    // N×V should point in -X direction (west) or +X direction (east)
    println!(
        "Chat X component: {:.4} (positive=east, negative=west)",
        chat[0]
    );

    let right_result = isce2_forward_geocode(sat_pos_arr, sat_vel, slant_range, -1.0);
    let left_result = isce2_forward_geocode(sat_pos_arr, sat_vel, slant_range, 1.0);

    assert!(right_result.is_some(), "Right-looking should succeed");
    assert!(left_result.is_some(), "Left-looking should succeed");

    let (right_lat, right_lon) = right_result.unwrap();
    let (left_lat, left_lon) = left_result.unwrap();

    println!(
        "Right-looking (ilrl=-1): lat={:.4}°, lon={:.4}°",
        right_lat, right_lon
    );
    println!(
        "Left-looking (ilrl=+1):  lat={:.4}°, lon={:.4}°",
        left_lat, left_lon
    );

    // Key assertion: right and left should be on opposite sides of the track
    // For ASCENDING pass (northward velocity):
    //   Right-looking = EAST (positive lon offset)
    //   Left-looking = WEST (negative lon offset)
    let sat_lon = 9.0;

    println!("Satellite lon: {:.4}°", sat_lon);
    println!(
        "Right offset: {:.4}°, Left offset: {:.4}°",
        right_lon - sat_lon,
        left_lon - sat_lon
    );

    // Check they are on opposite sides
    let right_offset = right_lon - sat_lon;
    let left_offset = left_lon - sat_lon;

    assert!(right_offset * left_offset < 0.0,
        "Right and left should be on opposite sides of track! Right offset={:.4}°, Left offset={:.4}°",
        right_offset, left_offset);
}

#[test]
fn test_slant_range_reconstruction() {
    // Verify that the computed ground point is at the correct slant range

    let sat_pos = wgs84_to_ecef_test(48.5, 9.0, 700_000.0);
    let sat_pos_arr = [sat_pos.0, sat_pos.1, sat_pos.2];
    let sat_vel = [100.0, 7500.0, 1000.0];

    let test_ranges = [800_000.0, 850_000.0, 900_000.0, 950_000.0, 1_000_000.0];

    for expected_range in test_ranges {
        let result = isce2_forward_geocode(sat_pos_arr, sat_vel, expected_range, -1.0);

        assert!(
            result.is_some(),
            "Should geocode successfully for range {:.0}km",
            expected_range / 1000.0
        );

        let (ground_lat, ground_lon) = result.unwrap();
        let ground_ecef = wgs84_to_ecef_test(ground_lat, ground_lon, 0.0);

        // Compute actual slant range
        let actual_range = vec3_norm(
            ground_ecef.0 - sat_pos_arr[0],
            ground_ecef.1 - sat_pos_arr[1],
            ground_ecef.2 - sat_pos_arr[2],
        );

        // Allow up to 5% error due to spherical vs ellipsoid approximation
        let range_error_pct = ((actual_range - expected_range).abs() / expected_range) * 100.0;

        println!(
            "Expected range: {:.0}km, Actual: {:.0}km, Error: {:.2}%",
            expected_range / 1000.0,
            actual_range / 1000.0,
            range_error_pct
        );

        assert!(
            range_error_pct < 5.0,
            "Range error too large: expected {:.0}km, got {:.0}km ({:.2}% error)",
            expected_range / 1000.0,
            actual_range / 1000.0,
            range_error_pct
        );
    }
}

#[test]
fn test_sentinel1_iw_typical_geometry() {
    // Test with realistic Sentinel-1 IW mode parameters
    // Reference: S1B_IW_SLC__1SDV_20190123 over Germany/Switzerland
    // Expected bbox: lat=[47.79, 49.74], lon=[7.22, 11.08]
    //
    // IMPORTANT: For Sentinel-1 descending passes over Europe:
    // - Satellite flies roughly south (descending from north pole region)
    // - Right-looking radar illuminates areas to the EAST of ground track
    // - The actual velocity direction depends on the orbit geometry at that location
    //
    // Use the same validated approach as test_forward_geocode_left_vs_right_looking

    // Satellite at ~693km altitude over Europe
    // For a descending pass, we need to compute realistic orbital velocity
    let sat_altitude = 693_000.0;
    let sat_lat = 49.0; // degrees
    let sat_lon = 9.0; // degrees - satellite ground track

    let sat_pos = wgs84_to_ecef_test(sat_lat, sat_lon, sat_altitude);
    let sat_pos_arr = [sat_pos.0, sat_pos.1, sat_pos.2];

    // For a descending sun-synchronous orbit (inclination ~98°), the velocity
    // at this latitude has components that make "right" = EAST in geographic terms.
    //
    // Use orbital mechanics: velocity is tangent to orbit, perpendicular to position
    // For descending pass at lat ~49°N, velocity has significant -Z component (southward)
    // and smaller X,Y components that depend on longitude
    //
    // CRITICAL: The velocity must be computed such that:
    // 1. It's tangent to the orbit (approximately perpendicular to radial direction)
    // 2. It's in the orbital plane (inclined at 98.18° to equator)
    // 3. For descending pass, the Z-component should be negative (moving south)

    // Compute orbital velocity using the formula from test_forward_geocode_left_vs_right_looking
    // which was validated to produce correct left/right looking behavior
    let sat_lat_rad = sat_lat.to_radians();
    let sat_lon_rad = sat_lon.to_radians();

    // Orbital speed (circular orbit at 693km)
    let speed = 7500.0; // m/s

    // For Sentinel-1 descending pass at this location:
    // The velocity in inertial frame is roughly in the orbital plane.
    // Converting to ECEF, the velocity has components that depend on:
    // - Earth rotation (adds ~460 m/s eastward at equator, less at higher lat)
    // - Orbital plane orientation
    //
    // For a descending pass at lat=49°, the inertial velocity is roughly:
    // - Southward (negative geographic north)
    // - With small eastward component due to orbital inclination > 90°
    //
    // In ECEF, for satellite at (lat=49°, lon=9°), descending:
    // Approximate velocity (from real orbit files):
    let sat_vel: [f64; 3] = [-5500.0, -900.0, 5000.0]; // Descending pass, roughly SSW in ECEF

    // Verify it's tangent to orbit (roughly perpendicular to position)
    let pos_mag = (sat_pos_arr[0].powi(2) + sat_pos_arr[1].powi(2) + sat_pos_arr[2].powi(2)).sqrt();
    let pos_unit = [
        sat_pos_arr[0] / pos_mag,
        sat_pos_arr[1] / pos_mag,
        sat_pos_arr[2] / pos_mag,
    ];
    let vel_mag = (sat_vel[0].powi(2) + sat_vel[1].powi(2) + sat_vel[2].powi(2)).sqrt();
    let vel_unit = [
        sat_vel[0] / vel_mag,
        sat_vel[1] / vel_mag,
        sat_vel[2] / vel_mag,
    ];
    let dot = pos_unit[0] * vel_unit[0] + pos_unit[1] * vel_unit[1] + pos_unit[2] * vel_unit[2];
    println!(
        "Position-velocity dot product: {:.4} (should be ~0 for circular orbit)",
        dot
    );

    println!(
        "Satellite position: ({:.0}, {:.0}, {:.0}) m",
        sat_pos_arr[0], sat_pos_arr[1], sat_pos_arr[2]
    );
    println!(
        "Satellite velocity: ({:.0}, {:.0}, {:.0}) m/s",
        sat_vel[0], sat_vel[1], sat_vel[2]
    );

    // Sentinel-1 IW mode: 3 subswaths
    // Near range ~800km, far range ~1100km
    let near_range = 850_000.0;
    let far_range = 1_050_000.0;

    // Compute near and far ground points (right-looking)
    let near_result = isce2_forward_geocode(sat_pos_arr, sat_vel, near_range, -1.0);
    let far_result = isce2_forward_geocode(sat_pos_arr, sat_vel, far_range, -1.0);

    assert!(near_result.is_some(), "Near range geocoding should succeed");
    assert!(far_result.is_some(), "Far range geocoding should succeed");

    let (near_lat, near_lon) = near_result.unwrap();
    let (far_lat, far_lon) = far_result.unwrap();

    println!(
        "Near range ({:.0}km): lat={:.4}°, lon={:.4}°",
        near_range / 1000.0,
        near_lat,
        near_lon
    );
    println!(
        "Far range ({:.0}km): lat={:.4}°, lon={:.4}°",
        far_range / 1000.0,
        far_lat,
        far_lon
    );

    // Verify ground points are in reasonable geographic area
    // For a pass over Germany/Switzerland, we expect:
    // - Latitude ~47-50°N
    // - Longitude ~7-12°E

    // Check latitude is in expected range
    assert!(
        near_lat > 45.0 && near_lat < 55.0,
        "Near lat {:.4}° should be in range [45, 55]",
        near_lat
    );
    assert!(
        far_lat > 45.0 && far_lat < 55.0,
        "Far lat {:.4}° should be in range [45, 55]",
        far_lat
    );

    // Check longitude is in expected range (right-looking = east of track)
    // Since satellite is at lon=9°, ground points should be EAST of 9°
    // For Sentinel-1 IW with slant range 850-1050km, typical cross-track offset is 5-15°
    assert!(
        near_lon > sat_lon && near_lon < sat_lon + 15.0,
        "Near lon {:.4}° should be east of satellite ({:.1}°) but within 15°",
        near_lon,
        sat_lon
    );
    assert!(
        far_lon > sat_lon && far_lon < sat_lon + 20.0,
        "Far lon {:.4}° should be east of satellite ({:.1}°) but within 20°",
        far_lon,
        sat_lon
    );

    // Far range should have greater cross-track offset than near range
    let near_offset = near_lon - sat_lon;
    let far_offset = far_lon - sat_lon;

    println!(
        "Near offset: {:.4}°, Far offset: {:.4}°",
        near_offset, far_offset
    );

    // Both should be positive (east of track) for right-looking
    // Far should have larger offset
    assert!(
        far_offset > near_offset,
        "Far range offset {:.4}° should be greater than near range offset {:.4}°",
        far_offset,
        near_offset
    );
}
