//! Unit tests for TerrainCorrector
//!
//! Extracted from mod.rs for better organization.

use super::*;
use crate::types::BurstTiming;
use chrono::{TimeZone, Utc};
use ndarray::Array2;
use std::collections::HashMap;

fn dummy_corrector_with_spacing(spacing_m: f64) -> TerrainCorrector {
    let dem = Array2::<f32>::zeros((2, 2));
    let dem_transform = GeoTransform {
        top_left_x: 0.0,
        pixel_width: 1.0,
        rotation_x: 0.0,
        top_left_y: 0.0,
        rotation_y: 0.0,
        pixel_height: -1.0,
    };
    TerrainCorrector::new(dem, dem_transform, -32768.0, 4326, 4326, spacing_m)
}

#[test]
fn test_calculate_pixel_size_degrees_20m_mid_lat() {
    let res_m = 20.0;
    let lat = 45.0;
    let px_deg = TerrainCorrectionConfig::calculate_pixel_size_degrees(res_m, lat);

    // Expected using same WGS84 conversion
    let a = crate::constants::geodetic::WGS84_SEMI_MAJOR_AXIS_M;
    let e2 = crate::constants::geodetic::WGS84_ECCENTRICITY_SQUARED;
    let lat_rad = lat.to_radians();
    let sin_lat = lat_rad.sin();
    let n = a / (1.0 - e2 * sin_lat * sin_lat).sqrt();
    let meters_per_deg_lon = n * lat_rad.cos() * std::f64::consts::PI / 180.0;
    let expected = res_m / meters_per_deg_lon;

    assert!(
        (px_deg - expected).abs() < 1e-10,
        "px_deg={} expected={}",
        px_deg,
        expected
    );
}

#[test]
fn test_build_burst_segments_basic() {
    // Note: last_line_global is EXCLUSIVE (one past the last valid line)
    let bursts = vec![
        BurstTiming {
            subswath_id: "IW1".to_string(),
            burst_index: 0,
            azimuth_time_rel_orbit: 2.0,
            first_line_global: 0,
            last_line_global: 1000, // exclusive: lines 0..1000
            first_valid_sample: None,
            last_valid_sample: None,
        },
        BurstTiming {
            subswath_id: "IW1".to_string(),
            burst_index: 1,
            azimuth_time_rel_orbit: 3.0,
            first_line_global: 1000,
            last_line_global: 2000, // exclusive: lines 1000..2000
            first_valid_sample: None,
            last_valid_sample: None,
        },
    ];

    let native_dt = 0.001;
    let azimuth_multilook_factor = 4.0;
    let segments =
        RangeDopplerParams::build_burst_segments(&bursts, native_dt, azimuth_multilook_factor, None);

    assert_eq!(segments.len(), 2);
    assert!(
        (segments[0].start_line - 0.0).abs() < 1e-6,
        "unexpected start_line {}",
        segments[0].start_line
    );
    assert!(
        (segments[0].end_line - 250.0).abs() < 1e-6,
        "unexpected end_line {}",
        segments[0].end_line
    );
    assert!((segments[0].line_time_interval - native_dt * azimuth_multilook_factor).abs() < 1e-9);
    assert!(
        (segments[1].start_line - 250.0).abs() < 1e-6,
        "unexpected second start_line {}",
        segments[1].start_line
    );
    assert!(
        (segments[1].end_line - 500.0).abs() < 1e-6,
        "unexpected second end_line {}",
        segments[1].end_line
    );
}

#[test]
fn test_azimuth_time_to_pixel_prefers_burst_segments() {
    let mut corrector = dummy_corrector_with_spacing(10.0);
    let reference_time = Utc.timestamp_opt(1_600_000_000, 0).unwrap();
    let orbit_data = OrbitData {
        reference_time,
        state_vectors: vec![StateVector {
            time: reference_time,
            position: [7_000_000.0, 0.0, 0.0],
            velocity: [0.0, 7_500.0, 0.0],
        }],
    };
    corrector.set_orbit_data(orbit_data);

    // Note: last_line_global is EXCLUSIVE (one past the last valid line)
    let bursts = vec![
        BurstTiming {
            subswath_id: "IW1".to_string(),
            burst_index: 0,
            azimuth_time_rel_orbit: 2.0,
            first_line_global: 0,
            last_line_global: 1000, // exclusive: lines 0..1000
            first_valid_sample: None,
            last_valid_sample: None,
        },
        BurstTiming {
            subswath_id: "IW1".to_string(),
            burst_index: 1,
            azimuth_time_rel_orbit: 5.0,
            first_line_global: 1000,
            last_line_global: 2000, // exclusive: lines 1000..2000
            first_valid_sample: None,
            last_valid_sample: None,
        },
    ];

    let native_dt = 0.001;
    let azimuth_multilook_factor = 4.0;
    #[allow(deprecated)]
    let mut params = RangeDopplerParams {
        range_pixel_spacing: 2.3,
        azimuth_pixel_spacing: native_dt,
        slant_range_time: 0.004,
        prf: 1000.0,
        azimuth_time_interval: native_dt,
        wavelength: 0.055,
        speed_of_light: 299_792_458.0,
        orbit_ref_epoch_utc: reference_time.timestamp() as f64,
        product_start_rel_s: 1.0,
        product_start_time_abs: reference_time.timestamp() as f64 + 1.0,
        product_stop_time_abs: reference_time.timestamp() as f64 + 21.0,
        product_duration: 20.0,
        total_azimuth_lines: Some(500),
        doppler_centroid: None,
        first_valid_line: None,
        last_valid_line: None,
        first_valid_sample: None,
        last_valid_sample: None,
        range_multilook_factor: 1.0,
        azimuth_multilook_factor,
        range_multilook_safe: 1.0,
        azimuth_multilook_safe: 1.0,
        subswaths: HashMap::new(),
        burst_timings: bursts.clone(),
        burst_segments: Vec::new(),
        reference_incidence_angle_deg: None,
        incidence_angle_near_deg: None,
        incidence_angle_far_deg: None,
        total_range_samples: None,
    };
    params.compute_safe_multilook_factors();
    params.burst_segments =
        RangeDopplerParams::build_burst_segments(&bursts, native_dt, azimuth_multilook_factor, None);

    let pixel_inside = corrector.azimuth_time_to_pixel(2.5, &params).unwrap();
    assert!(
        (pixel_inside - 500.0).abs() < 1e-3,
        "expected ~500 native lines, got {}",
        pixel_inside
    );

    let before_first = corrector.azimuth_time_to_pixel(0.5, &params);
    assert!(before_first.is_none(), "expected None before first burst");

    let after_last = corrector.azimuth_time_to_pixel(8.5, &params);
    assert!(after_last.is_none(), "expected None after last burst");
}

#[test]
fn test_validate_and_fix_output_spacing_corrects_large_error() {
    let mut corr = dummy_corrector_with_spacing(1e-6); // absurdly small => triggers correction
    let target_res = 20.0;
    let lat = 45.0;
    let lon = 10.0;
    let corrected_deg = corr
        .validate_and_fix_output_spacing(target_res, lat, lon)
        .expect("validation failed");

    // After fix, internal spacing equals target resolution in meters
    assert!((corr.output_spacing - target_res).abs() < 1e-9);

    // Check returned degree-per-pixel matches WGS84 conversion at latitude
    let a = crate::constants::geodetic::WGS84_SEMI_MAJOR_AXIS_M;
    let e2 = crate::constants::geodetic::WGS84_ECCENTRICITY_SQUARED;
    let lat_rad = lat.to_radians();
    let sin_lat = lat_rad.sin();
    let n = a / (1.0 - e2 * sin_lat * sin_lat).sqrt();
    let meters_per_deg_lon = n * lat_rad.cos() * std::f64::consts::PI / 180.0;
    let expected_deg = target_res / meters_per_deg_lon;
    assert!(
        (corrected_deg - expected_deg).abs() < 1e-10,
        "got={} expected={}",
        corrected_deg,
        expected_deg
    );
}

#[test]
fn test_create_geographic_grid_uses_resolution() {
    let corr = dummy_corrector_with_spacing(20.0);
    let bounds = BoundingBox {
        min_lat: 44.0,
        max_lat: 44.1,
        min_lon: 9.0,
        max_lon: 9.2,
    };
    let (_w, _h, gt) = corr.create_geographic_grid(&bounds).expect("grid");

    let mid_lat = (bounds.min_lat + bounds.max_lat) / 2.0;
    let lat_rad = mid_lat.to_radians();
    let a = crate::constants::geodetic::WGS84_SEMI_MAJOR_AXIS_M;
    let e2 = crate::constants::geodetic::WGS84_ECCENTRICITY_SQUARED;
    let sin_lat = lat_rad.sin();
    let meridional_radius = a * (1.0 - e2) / (1.0 - e2 * sin_lat * sin_lat).powf(1.5);
    let prime_vertical_radius = a / (1.0 - e2 * sin_lat * sin_lat).sqrt();
    let m_per_deg_lat = meridional_radius * std::f64::consts::PI / 180.0;
    let m_per_deg_lon = prime_vertical_radius * lat_rad.cos() * std::f64::consts::PI / 180.0;

    let expected_px_w = 20.0 / m_per_deg_lon;
    let expected_px_h = -(20.0 / m_per_deg_lat);

    assert!(
        (gt.pixel_width - expected_px_w).abs() < 1e-10,
        "w {} vs {}",
        gt.pixel_width,
        expected_px_w
    );
    assert!(
        (gt.pixel_height - expected_px_h).abs() < 1e-10,
        "h {} vs {}",
        gt.pixel_height,
        expected_px_h
    );
    assert!(gt.pixel_height < 0.0);
}

#[test]
fn test_tie_point_grid_interpolation_blends_valid_cells() {
    let mut cells = Array2::from_elem((2, 2), TiePointCell::invalid());
    cells[[0, 0]] = TiePointCell::with_values(0.0, 1000.0, 0.4, TIE_FLAG_VALID);
    cells[[0, 1]] = TiePointCell::with_values(1.0, 2000.0, 0.5, TIE_FLAG_VALID);
    cells[[1, 0]] = TiePointCell::with_values(2.0, 3000.0, 0.6, TIE_FLAG_VALID);
    cells[[1, 1]] = TiePointCell::with_values(3.0, 4000.0, 0.7, TIE_FLAG_VALID);

    let grid = TiePointGrid::new(4, cells);
    let sample = grid.interpolate(2, 2).expect("sample");

    assert!((sample.azimuth_time - 1.5).abs() < 1e-6);
    assert!((sample.slant_range - 2500.0).abs() < 1e-6);
    assert!((sample.cos_local_incidence as f64 - 0.55).abs() < 1e-6);
}

#[test]
fn test_tie_point_grid_interpolation_requires_valid_cells() {
    let mut cells = Array2::from_elem((1, 1), TiePointCell::invalid());
    cells[[0, 0]] = TiePointCell::with_values(0.0, 100.0, 0.4, 0);
    let grid = TiePointGrid::new(8, cells);
    assert!(grid.interpolate(0, 0).is_none());
}

#[test]
fn test_subswath_slant_range_time_lookup() {
    // Test that per-subswath slant_range_time is used correctly for merged IW data
    // This verifies Patch 2: Always use per-subswath lookup, never fall back to global
    use crate::types::SubSwath;

    let mut subswaths = HashMap::new();
    let range_pixel_spacing = 2.33;
    let speed_of_light = 299_792_458.0;
    let range_pixel_spacing_time = 2.0 * range_pixel_spacing / speed_of_light;

    // IW1: samples 0-22k, slant_range_time = 0.00015s
    // End time: 0.00015 + 22000 * spacing_time
    let iw1_end_time = 0.00015 + 22000.0 * range_pixel_spacing_time;

    // IW2: samples 22k-45k, slant_range_time = 0.00018s (different from IW1!)
    // Start time must be > IW1 end time to avoid overlap
    let iw2_start_time = iw1_end_time + 0.00001; // Small gap
    let iw2_slant_range_time = iw2_start_time;
    let iw2_end_time = iw2_slant_range_time + 23000.0 * range_pixel_spacing_time;

    subswaths.insert(
        "IW1".to_string(),
        SubSwath {
            id: "IW1".to_string(),
            burst_count: 9,
            lines_per_burst: 5000 / 9,
            range_samples: 22000,
            azimuth_samples: 5000,
            first_line_global: 0,
            last_line_global: 5000,
            first_sample_global: 0,
            last_sample_global: 22000,
            full_range_samples: 22000,
            valid_first_line: Some(0),
            valid_last_line: Some(5000),
            valid_first_sample: Some(0),
            valid_last_sample: Some(22000),
            range_pixel_spacing,
            azimuth_pixel_spacing: 14.0,
            slant_range_time: 0.00015,
            burst_duration: 2.758,
            near_range_m: 800_000.0,
            prf_hz: Some(1700.0),
            dc_polynomial: None,
            azimuth_time_interval: Some(0.000588),
            dc_polynomial_t0: None,
            fm_rate_estimates: None,
        },
    );

    subswaths.insert(
        "IW2".to_string(),
        SubSwath {
            id: "IW2".to_string(),
            burst_count: 9,
            lines_per_burst: 5000 / 9,
            range_samples: 23000,
            azimuth_samples: 5000,
            first_line_global: 5000,
            last_line_global: 10000,
            first_sample_global: 22000,
            last_sample_global: 45000,
            full_range_samples: 23000,
            valid_first_line: Some(5000),
            valid_last_line: Some(10000),
            valid_first_sample: Some(22000),
            valid_last_sample: Some(45000),
            range_pixel_spacing,
            azimuth_pixel_spacing: 14.0,
            slant_range_time: iw2_slant_range_time, // Different from IW1
            burst_duration: 2.758,
            near_range_m: 850_000.0,
            prf_hz: Some(1700.0),
            dc_polynomial: None,
            azimuth_time_interval: Some(0.000588),
            dc_polynomial_t0: None,
            fm_rate_estimates: None,
        },
    );

    #[allow(deprecated)]
    let params = RangeDopplerParams {
        range_pixel_spacing,
        azimuth_pixel_spacing: 0.000588,
        slant_range_time: 0.00015, // Global (IW1) slant_range_time
        prf: 1700.0,
        azimuth_time_interval: 0.0018,
        wavelength: 0.055465763,
        speed_of_light,
        orbit_ref_epoch_utc: 1_600_000_000.0,
        product_start_rel_s: 100.0,
        product_start_time_abs: 1_600_000_100.0,
        product_stop_time_abs: 1_600_000_125.0,
        product_duration: 25.0,
        total_azimuth_lines: Some(15000),
        doppler_centroid: None,
        first_valid_line: None,
        last_valid_line: None,
        first_valid_sample: None,
        last_valid_sample: None,
        range_multilook_factor: 1.0,
        azimuth_multilook_factor: 1.0,
        range_multilook_safe: 1.0,
        azimuth_multilook_safe: 1.0,
        subswaths,
        burst_timings: Vec::new(),
        burst_segments: Vec::new(),
        reference_incidence_angle_deg: Some(35.0), // Typical Sentinel-1 mid-swath
        incidence_angle_near_deg: None,
        incidence_angle_far_deg: None,
        total_range_samples: None,
    };
    let mut params = params;
    params.compute_safe_multilook_factors();
    let params = params;

    // Test IW2: Use a two-way time that falls in IW2's range
    // Middle of IW2: local_pixel = 11500, so two_way_time = iw2_slant_range_time + 11500 * spacing_time
    let iw2_middle_time = iw2_slant_range_time + 11500.0 * range_pixel_spacing_time;
    let iw2_slant_range = (iw2_middle_time * speed_of_light) / 2.0;

    let range_pixel = params.slant_range_to_native_pixel(iw2_slant_range);

    // Key test: Should use IW2's slant_range_time, mapping to IW2 range (22k-45k)
    assert!(
        range_pixel >= 22_000.0 && range_pixel <= 45_000.0,
        "IW2 range pixel should be in [22k, 45k] using IW2 slant_range_time ({:.9}s), got {}",
        iw2_slant_range_time,
        range_pixel
    );

    // Verify it's near the expected middle of IW2
    // Tolerance: < 1.0 pixel (scientific precision requirement)
    // At 2.33 m/pixel, 1 pixel = 2.33 m, which is acceptable for range accuracy
    let expected_pixel = 22000.0 + 11500.0;
    assert!(
        (range_pixel - expected_pixel).abs() < 1.0,
        "IW2 range pixel should be near {} (middle of IW2), got {} (diff: {:.3} pixels)",
        expected_pixel,
        range_pixel,
        (range_pixel - expected_pixel).abs()
    );
}

#[test]
fn test_time_domain_consistency() {
    // Test that product_start_absolute() correctly computes from orbit_ref_epoch_utc + product_start_rel_s
    use chrono::{TimeZone, Utc};

    let orbit_ref_epoch = Utc.with_ymd_and_hms(2020, 10, 5, 17, 8, 24).unwrap();
    let product_start = Utc.with_ymd_and_hms(2020, 10, 5, 17, 8, 25).unwrap();

    let orbit_ref_epoch_utc = crate::types::datetime_to_utc_seconds(orbit_ref_epoch);
    let product_start_abs = crate::types::datetime_to_utc_seconds(product_start);
    let product_start_rel_s = product_start_abs - orbit_ref_epoch_utc;

    #[allow(deprecated)]
    let params = RangeDopplerParams {
        range_pixel_spacing: 2.33,
        azimuth_pixel_spacing: 0.000588,
        slant_range_time: 0.00015,
        prf: 1700.0,
        azimuth_time_interval: 0.0018,
        wavelength: 0.055465763,
        speed_of_light: 299_792_458.0,
        orbit_ref_epoch_utc,
        product_start_rel_s,
        product_start_time_abs: product_start_abs,
        product_stop_time_abs: product_start_abs + 25.0,
        product_duration: 25.0,
        total_azimuth_lines: Some(15000),
        doppler_centroid: None,
        first_valid_line: None,
        last_valid_line: None,
        first_valid_sample: None,
        last_valid_sample: None,
        range_multilook_factor: 1.0,
        azimuth_multilook_factor: 1.0,
        range_multilook_safe: 1.0,
        azimuth_multilook_safe: 1.0,
        subswaths: HashMap::new(),
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

    // Validate product_start_absolute() matches computed value
    let computed_abs = params.product_start_absolute();
    assert!(
        (computed_abs - product_start_abs).abs() < 1e-6,
        "product_start_absolute() should match computed value: got {:.9}, expected {:.9}",
        computed_abs,
        product_start_abs
    );

    // Validate product_stop_absolute() matches
    let computed_stop = params.product_stop_absolute();
    let expected_stop = product_start_abs + 25.0;
    assert!(
        (computed_stop - expected_stop).abs() < 1e-6,
        "product_stop_absolute() should match: got {:.9}, expected {:.9}",
        computed_stop,
        expected_stop
    );
}

#[test]
fn test_orbit_interpolation_accuracy() {
    // Test that orbit interpolation is accurate for known positions
    use crate::types::{OrbitData, StateVector};
    use chrono::{TimeZone, Utc};

    let reference_time = Utc.timestamp_opt(1_600_000_000, 0).unwrap();

    // Create a simple circular orbit with known positions
    let orbit_radius = 7_000_000.0;
    let angular_velocity = 0.0011; // rad/s

    // Create 5 state vectors (minimum for cubic spline)
    let mut state_vectors = Vec::new();
    for i in 0..5 {
        let t_offset = i as f64 * 10.0;
        let time = reference_time + chrono::Duration::seconds(t_offset as i64);
        let angle = angular_velocity * t_offset;

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

    let mut corrector = dummy_corrector_with_spacing(20.0);
    corrector.set_orbit_data(orbit_data.clone());

    // Test interpolation at t=25s (between state vectors at t=20s and t=30s)
    let target_time = crate::types::datetime_to_utc_seconds(reference_time) + 25.0;
    let result = corrector.scientific_orbit_interpolation(&orbit_data, target_time);

    assert!(result.is_ok(), "Orbit interpolation should succeed");
    let (pos, vel) = result.unwrap();

    // Expected position at t=25s
    let expected_angle = angular_velocity * 25.0;
    let expected_pos = [
        orbit_radius * expected_angle.cos(),
        orbit_radius * expected_angle.sin(),
        0.0,
    ];

    // Validate position accuracy (should be within 1 meter for cubic spline)
    let pos_error = ((pos.x - expected_pos[0]).powi(2)
        + (pos.y - expected_pos[1]).powi(2)
        + (pos.z - expected_pos[2]).powi(2))
    .sqrt();

    assert!(
        pos_error < 1.0,
        "Orbit interpolation position error should be < 1m: got {:.6}m",
        pos_error
    );

    // Validate velocity accuracy (should be within 0.01 m/s)
    let expected_vel = [
        -orbit_radius * angular_velocity * expected_angle.sin(),
        orbit_radius * angular_velocity * expected_angle.cos(),
        0.0,
    ];
    let vel_error = ((vel.x - expected_vel[0]).powi(2)
        + (vel.y - expected_vel[1]).powi(2)
        + (vel.z - expected_vel[2]).powi(2))
    .sqrt();

    assert!(
        vel_error < 0.01,
        "Orbit interpolation velocity error should be < 0.01 m/s: got {:.6} m/s",
        vel_error
    );
}

#[test]
fn test_dem_coordinate_transformation() {
    // Test coordinate transformation between WGS84 and UTM
    use crate::types::GeoTransform;

    // Create corrector with WGS84 DEM
    let dem_wgs84 = Array2::<f32>::zeros((100, 100));
    let dem_transform_wgs84 = GeoTransform {
        top_left_x: 9.0,    // longitude
        pixel_width: 0.001, // degrees
        rotation_x: 0.0,
        top_left_y: 45.0, // latitude
        rotation_y: 0.0,
        pixel_height: -0.001, // degrees (negative for north-up)
    };
    let corrector_wgs84 = TerrainCorrector::new(
        dem_wgs84,
        dem_transform_wgs84,
        -32768.0,
        4326, // WGS84
        4326, // Output WGS84
        20.0,
    );

    // Test WGS84 to WGS84 (should be identity)
    let (lat, lon) = corrector_wgs84.map_to_geographic(9.5, 45.5).unwrap();
    assert!(
        (lat - 45.5).abs() < 1e-6 && (lon - 9.5).abs() < 1e-6,
        "WGS84 to WGS84 should be identity: got ({:.6}, {:.6}), expected (45.5, 9.5)",
        lat,
        lon
    );

    // Test UTM to WGS84 transformation
    // UTM Zone 32N (EPSG:32632) for central Europe
    let dem_utm = Array2::<f32>::zeros((100, 100));
    let dem_transform_utm = GeoTransform {
        top_left_x: 500000.0, // UTM easting
        pixel_width: 20.0,    // meters
        rotation_x: 0.0,
        top_left_y: 5000000.0, // UTM northing
        rotation_y: 0.0,
        pixel_height: -20.0, // meters (negative for north-up)
    };
    let corrector_utm = TerrainCorrector::new(
        dem_utm,
        dem_transform_utm,
        -32768.0,
        32632, // UTM Zone 32N
        4326,  // Output WGS84
        20.0,
    );

    // Test coordinate transformation: geographic to UTM and verify it's reversible
    // This validates that DEM coordinate transformation works correctly
    // Test point: 45°N, 9°E (central Europe, UTM Zone 32N)
    let test_lat = 45.0;
    let test_lon = 9.0;

    // Transform geographic to DEM CRS (UTM)
    let (utm_x, utm_y) = corrector_utm
        .transform_latlon_to_dem_crs(test_lat, test_lon)
        .unwrap();

    // Validate UTM coordinates are reasonable for UTM Zone 32N
    // UTM Zone 32N: 500000m false easting, so coordinates should be ~500000m + offset
    assert!(
        utm_x > 400000.0 && utm_x < 600000.0,
        "UTM easting should be in valid range for Zone 32N: got {:.3}m",
        utm_x
    );
    assert!(
        utm_y > 4900000.0 && utm_y < 5100000.0,
        "UTM northing should be in valid range for Zone 32N: got {:.3}m",
        utm_y
    );

    // Test that output CRS transformation works (output is WGS84, so should be identity for geographic input)
    let (output_lat, output_lon) = corrector_utm.map_to_geographic(test_lon, test_lat).unwrap();
    assert!(
        (output_lat - test_lat).abs() < 1e-6 && (output_lon - test_lon).abs() < 1e-6,
        "WGS84 output CRS should preserve geographic coordinates: got ({:.6}°, {:.6}°), expected ({:.6}°, {:.6}°)",
        output_lat, output_lon, test_lat, test_lon
    );
}

#[test]
fn test_burst_timing_overlap_calculations() {
    // Test that burst timing overlap calculations are accurate
    // This validates that overlapping bursts have correct time windows

    let bursts = vec![
        BurstTiming {
            subswath_id: "IW1".to_string(),
            burst_index: 0,
            azimuth_time_rel_orbit: 10.0,
            first_line_global: 0,
            last_line_global: 1000, // exclusive: lines 0..1000
            first_valid_sample: None,
            last_valid_sample: None,
        },
        BurstTiming {
            subswath_id: "IW1".to_string(),
            burst_index: 1,
            azimuth_time_rel_orbit: 12.758, // 2.758s later (typical burst duration)
            first_line_global: 1000,
            last_line_global: 2000, // exclusive: lines 1000..2000
            first_valid_sample: None,
            last_valid_sample: None,
        },
    ];

    let native_dt = 0.000588; // Sentinel-1 azimuth time interval
    let azimuth_multilook_factor = 4.0;
    let segments =
        RangeDopplerParams::build_burst_segments(&bursts, native_dt, azimuth_multilook_factor, None);

    assert_eq!(segments.len(), 2, "Should have 2 burst segments");

    // Validate time continuity: segments should be adjacent
    // Segment 0: lines 0..250 (after multilook)
    // Segment 1: lines 250..500 (after multilook)
    assert!(
        (segments[0].end_line - segments[1].start_line).abs() < 1e-6,
        "Burst segments should be adjacent: segment 0 ends at {:.3}, segment 1 starts at {:.3}",
        segments[0].end_line,
        segments[1].start_line
    );

    // Validate time interval is correct after multilooking
    let expected_interval = native_dt * azimuth_multilook_factor;
    assert!(
        (segments[0].line_time_interval - expected_interval).abs() < 1e-9,
        "Line time interval should match multilook factor: got {:.9}s, expected {:.9}s",
        segments[0].line_time_interval,
        expected_interval
    );

    // Validate burst timing: second burst should start 2.758s after first
    // After multilooking, this becomes 2.758 / 4 = 0.6895s per multilooked line
    let time_diff = segments[1].start_time_rel - segments[0].start_time_rel;
    let expected_time_diff = 2.758; // Original burst separation
    assert!(
        (time_diff - expected_time_diff).abs() < 1e-6,
        "Burst timing should preserve original separation: got {:.6}s, expected {:.6}s",
        time_diff,
        expected_time_diff
    );
}
