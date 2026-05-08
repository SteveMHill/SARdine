//! WGS84 geodetic ↔ ECEF coordinate conversions.
//!
//! Implements the standard closed-form geodetic-to-ECEF transform and
//! Bowring's iterative ECEF-to-geodetic inverse.
//!
//! All angles are in **radians** internally and **degrees** at the public API boundary.

use std::f64::consts::PI;

// ── WGS84 ellipsoid constants ──────────────────────────────────────

/// Semi-major axis (equatorial radius) in metres.
pub const WGS84_A: f64 = 6_378_137.0;

/// Flattening.
pub const WGS84_F: f64 = 1.0 / 298.257_223_563;

/// Semi-minor axis (polar radius) in metres.
pub const WGS84_B: f64 = WGS84_A * (1.0 - WGS84_F);

/// First eccentricity squared: e² = 2f - f².
pub const WGS84_E2: f64 = 2.0 * WGS84_F - WGS84_F * WGS84_F;

/// Second eccentricity squared: e'² = e² / (1 - e²).
pub const WGS84_EP2: f64 = WGS84_E2 / (1.0 - WGS84_E2);

// ── Public API ─────────────────────────────────────────────────────

/// Convert geodetic coordinates (degrees, metres) to ECEF (metres).
///
/// # Arguments
/// * `lat_deg` — WGS84 latitude in degrees (−90 to +90).
/// * `lon_deg` — WGS84 longitude in degrees (−180 to +180).
/// * `height_m` — Height above WGS84 ellipsoid in metres.
///
/// # Returns
/// `[X, Y, Z]` in metres (ECEF, Earth-fixed).
pub fn geodetic_to_ecef(lat_deg: f64, lon_deg: f64, height_m: f64) -> [f64; 3] {
    let lat = lat_deg * PI / 180.0;
    let lon = lon_deg * PI / 180.0;

    let sin_lat = lat.sin();
    let cos_lat = lat.cos();
    let sin_lon = lon.sin();
    let cos_lon = lon.cos();

    // Prime vertical radius of curvature
    let n = WGS84_A / (1.0 - WGS84_E2 * sin_lat * sin_lat).sqrt();

    let x = (n + height_m) * cos_lat * cos_lon;
    let y = (n + height_m) * cos_lat * sin_lon;
    let z = (n * (1.0 - WGS84_E2) + height_m) * sin_lat;

    [x, y, z]
}

/// Convert ECEF coordinates (metres) to geodetic (degrees, metres).
///
/// Uses Bowring's iterative method (converges to < 1 mm after 3 iterations).
///
/// # Returns
/// `(latitude_deg, longitude_deg, height_m)`.
pub fn ecef_to_geodetic(ecef: [f64; 3]) -> (f64, f64, f64) {
    let [x, y, z] = ecef;

    let lon = y.atan2(x);
    let p = (x * x + y * y).sqrt();

    // Initial guess (Bowring)
    let theta = (z * WGS84_A).atan2(p * WGS84_B);
    let mut lat = (z + WGS84_EP2 * WGS84_B * theta.sin().powi(3))
        .atan2(p - WGS84_E2 * WGS84_A * theta.cos().powi(3));

    // Iterate for convergence
    for _ in 0..5 {
        let sin_lat = lat.sin();
        let n = WGS84_A / (1.0 - WGS84_E2 * sin_lat * sin_lat).sqrt();
        let new_lat = (z + WGS84_E2 * n * sin_lat).atan2(p);
        if (new_lat - lat).abs() < 1e-12 {
            break;
        }
        lat = new_lat;
    }

    let sin_lat = lat.sin();
    let n = WGS84_A / (1.0 - WGS84_E2 * sin_lat * sin_lat).sqrt();
    let height = if lat.cos().abs() > 1e-10 {
        p / lat.cos() - n
    } else {
        z / lat.sin() - n * (1.0 - WGS84_E2)
    };

    (lat * 180.0 / PI, lon * 180.0 / PI, height)
}

/// Compute the dot product of two 3-vectors.
#[inline]
pub fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

/// Compute the Euclidean norm of a 3-vector.
#[inline]
pub fn norm3(v: [f64; 3]) -> f64 {
    dot3(v, v).sqrt()
}

/// Subtract two 3-vectors: a - b.
#[inline]
pub fn sub3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

/// Cross product of two 3-vectors.
#[inline]
pub fn cross3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

/// Normalize a 3-vector. Returns zero vector if input norm is zero.
#[inline]
pub fn normalize3(v: [f64; 3]) -> [f64; 3] {
    let n = norm3(v);
    if n == 0.0 {
        return [0.0, 0.0, 0.0];
    }
    [v[0] / n, v[1] / n, v[2] / n]
}

// ── Tests ──────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    /// Round-trip: geodetic → ECEF → geodetic should recover original.
    #[test]
    fn round_trip_equator() {
        let lat = 0.0;
        let lon = 0.0;
        let h = 0.0;
        let ecef = geodetic_to_ecef(lat, lon, h);
        // Equator should give X ≈ a, Y = 0, Z = 0
        assert!((ecef[0] - WGS84_A).abs() < 1e-3);
        assert!(ecef[1].abs() < 1e-3);
        assert!(ecef[2].abs() < 1e-3);

        let (lat2, lon2, h2) = ecef_to_geodetic(ecef);
        assert!((lat2 - lat).abs() < 1e-10);
        assert!((lon2 - lon).abs() < 1e-10);
        assert!((h2 - h).abs() < 1e-3);
    }

    #[test]
    fn round_trip_pole() {
        let lat = 90.0;
        let lon = 0.0;
        let h = 0.0;
        let ecef = geodetic_to_ecef(lat, lon, h);
        // North pole: X,Y ≈ 0, Z ≈ b
        assert!(ecef[0].abs() < 1e-3);
        assert!(ecef[1].abs() < 1e-3);
        assert!((ecef[2] - WGS84_B).abs() < 1e-3);

        let (lat2, _lon2, h2) = ecef_to_geodetic(ecef);
        assert!((lat2 - 90.0).abs() < 1e-10);
        assert!(h2.abs() < 1e-3);
    }

    #[test]
    fn round_trip_mid_latitude_with_altitude() {
        let lat = 50.5;
        let lon = 10.0;
        let h = 300.0;
        let ecef = geodetic_to_ecef(lat, lon, h);
        let (lat2, lon2, h2) = ecef_to_geodetic(ecef);
        assert!(
            (lat2 - lat).abs() < 1e-9,
            "lat: {} vs {}",
            lat2,
            lat
        );
        assert!(
            (lon2 - lon).abs() < 1e-9,
            "lon: {} vs {}",
            lon2,
            lon
        );
        assert!(
            (h2 - h).abs() < 1e-3,
            "height: {} vs {}",
            h2,
            h
        );
    }

    /// Known point: London (51.5°N, 0°W, 0m).
    /// ECEF values from standard references.
    #[test]
    fn london_ecef() {
        let ecef = geodetic_to_ecef(51.5, 0.0, 0.0);
        // Expected approximately: X ≈ 3978.7 km, Y ≈ 0, Z ≈ 4968.4 km
        assert!((ecef[0] - 3_978_700.0).abs() < 2000.0);
        assert!(ecef[1].abs() < 1.0);
        assert!((ecef[2] - 4_968_400.0).abs() < 2000.0);
    }

    /// Sentinel-1 typical orbit altitude: ~693 km above equator.
    #[test]
    fn satellite_altitude_round_trip() {
        let alt = 693_000.0;
        let ecef = geodetic_to_ecef(0.0, 45.0, alt);
        let (lat2, lon2, h2) = ecef_to_geodetic(ecef);
        assert!((lat2).abs() < 1e-9);
        assert!((lon2 - 45.0).abs() < 1e-9);
        assert!((h2 - alt).abs() < 0.01);
    }

    #[test]
    fn vector_ops() {
        let a = [1.0, 2.0, 3.0];
        let b = [4.0, 5.0, 6.0];
        assert_eq!(dot3(a, b), 32.0);
        assert!((norm3(a) - (14.0_f64).sqrt()).abs() < 1e-12);
        let d = sub3(b, a);
        assert_eq!(d, [3.0, 3.0, 3.0]);
        let c = cross3([1.0, 0.0, 0.0], [0.0, 1.0, 0.0]);
        assert_eq!(c, [0.0, 0.0, 1.0]);
    }
}
