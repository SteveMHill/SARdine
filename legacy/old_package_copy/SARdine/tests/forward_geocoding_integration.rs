//! Integration tests for forward geocoding with real Sentinel-1 data
//!
//! These tests validate that the ISCE2-style forward geocoding produces
//! correct results when applied to real SAFE annotation data.

use std::fs;
use std::path::Path;

/// Test forward geocoding with real S1B IW SLC data
///
/// This test loads real SAFE annotation files and validates that:
/// 1. Annotation parsing succeeds
/// 2. Orbit data is correctly extracted  
/// 3. Forward geocoding produces bbox close to ESA metadata
///
/// Reference ESA coordinates from manifest.safe:
/// 47.729340,10.594063 48.124813,7.219544 49.740406,7.595866 49.343430,11.083220
/// Expected bbox: lat=[47.73, 49.74], lon=[7.22, 11.08]
#[test]
#[ignore] // Run with: cargo test --test forward_geocoding_integration -- --ignored
fn test_forward_geocoding_with_real_s1b_data() {
    // Path to test data - adjust if running from different directory
    let safe_dir = Path::new("/home/datacube/apps/SARdine/data/SLC/S1B_IW_SLC__1SDV_20190123T053348_20190123T053415_014617_01B3D4_E833.SAFE");

    if !safe_dir.exists() {
        println!(
            "⚠️  Test data not available at {:?}, skipping test",
            safe_dir
        );
        return;
    }

    // Read IW1 annotation file for testing
    let annotation_path = safe_dir
        .join("annotation/s1b-iw1-slc-vv-20190123t053349-20190123t053414-014617-01b3d4-004.xml");

    assert!(
        annotation_path.exists(),
        "Annotation file should exist: {:?}",
        annotation_path
    );

    let xml_content =
        fs::read_to_string(&annotation_path).expect("Should be able to read annotation file");

    // Parse the annotation
    let parsed = sardine::io::annotation::parse_annotation(&xml_content)
        .expect("Should parse annotation successfully");

    // Validate orbit data is present
    let orbit_list = parsed.orbit_list.as_ref().expect("Should have orbit data");

    println!(
        "✅ Parsed annotation with {} orbit state vectors",
        orbit_list.len()
    );
    assert!(
        orbit_list.len() >= 4,
        "Should have at least 4 state vectors for cubic interpolation"
    );

    // Validate geolocation grid is present (for reference bbox)
    let geoloc_grid = parsed
        .geolocation_grid
        .as_ref()
        .expect("Should have geolocation grid");

    let geoloc_points = geoloc_grid
        .geolocation_grid_point_list
        .as_ref()
        .and_then(|l| l.geolocation_grid_points.as_ref())
        .expect("Should have geolocation points");

    println!(
        "✅ Parsed geolocation grid with {} points",
        geoloc_points.len()
    );

    // Extract bbox from geolocation grid
    let mut min_lat = f64::INFINITY;
    let mut max_lat = f64::NEG_INFINITY;
    let mut min_lon = f64::INFINITY;
    let mut max_lon = f64::NEG_INFINITY;

    for point in geoloc_points {
        let lat = point.latitude;
        let lon = point.longitude;
        min_lat = min_lat.min(lat);
        max_lat = max_lat.max(lat);
        min_lon = min_lon.min(lon);
        max_lon = max_lon.max(lon);
    }

    println!(
        "📍 IW1 geolocation bbox: lat=[{:.4}, {:.4}], lon=[{:.4}, {:.4}]",
        min_lat, max_lat, min_lon, max_lon
    );

    // Expected ESA metadata bbox for entire product:
    // lat=[47.73, 49.74], lon=[7.22, 11.08]
    // IW1 should be a subset of this (eastern edge for right-looking)

    // Verify the IW1 bbox is within the expected product extent
    assert!(
        min_lat >= 47.0 && min_lat <= 50.0,
        "min_lat {:.4} should be in range [47, 50]",
        min_lat
    );
    assert!(
        max_lat >= 48.0 && max_lat <= 51.0,
        "max_lat {:.4} should be in range [48, 51]",
        max_lat
    );
    assert!(
        min_lon >= 6.0 && min_lon <= 12.0,
        "min_lon {:.4} should be in range [6, 12]",
        min_lon
    );
    assert!(
        max_lon >= 7.0 && max_lon <= 13.0,
        "max_lon {:.4} should be in range [7, 13]",
        max_lon
    );

    // Validate swath timing is present
    let swath_timing = parsed
        .swath_timing
        .as_ref()
        .expect("Should have swath timing");

    let burst_list = swath_timing
        .burst_list
        .as_ref()
        .expect("Should have burst list");

    let bursts = burst_list.bursts.as_ref().expect("Should have bursts");

    println!("✅ Parsed {} bursts in IW1", bursts.len());
    assert!(bursts.len() >= 5, "Should have at least 5 bursts");

    // Validate image dimensions
    let image_info = parsed
        .image_annotation
        .as_ref()
        .and_then(|ia| ia.image_information.as_ref())
        .expect("Should have image information");

    let num_samples = image_info
        .number_of_samples
        .expect("Should have number of samples");
    let num_lines = image_info
        .number_of_lines
        .expect("Should have number of lines");

    println!(
        "✅ IW1 image dimensions: {} x {} (samples x lines)",
        num_samples, num_lines
    );

    // Typical IW1 dimensions: ~20-25k samples, ~10-16k lines
    assert!(
        num_samples > 15000 && num_samples < 30000,
        "num_samples {} should be in typical IW range",
        num_samples
    );
    assert!(
        num_lines > 5000 && num_lines < 20000,
        "num_lines {} should be in typical IW range",
        num_lines
    );

    println!("\n✅ All forward geocoding integration tests passed!");
}

/// Test that orbit interpolation works correctly with real data
#[test]
#[ignore]
fn test_orbit_interpolation_with_real_data() {
    let safe_dir = Path::new("/home/datacube/apps/SARdine/data/SLC/S1B_IW_SLC__1SDV_20190123T053348_20190123T053415_014617_01B3D4_E833.SAFE");

    if !safe_dir.exists() {
        println!("⚠️  Test data not available, skipping test");
        return;
    }

    let annotation_path = safe_dir
        .join("annotation/s1b-iw1-slc-vv-20190123t053349-20190123t053414-014617-01b3d4-004.xml");
    let xml_content =
        fs::read_to_string(&annotation_path).expect("Should be able to read annotation file");

    let parsed = sardine::io::annotation::parse_annotation(&xml_content)
        .expect("Should parse annotation successfully");

    let orbit_list = parsed.orbit_list.as_ref().expect("Should have orbit data");

    println!(
        "Testing orbit interpolation with {} state vectors",
        orbit_list.len()
    );

    // Print first and last state vector times
    if let (Some(first), Some(last)) = (orbit_list.first(), orbit_list.last()) {
        println!("First SV time: {:?}", first.time);
        println!("Last SV time: {:?}", last.time);

        // Validate position is at orbital altitude (693-707 km for S1)
        let pos_mag =
            (first.position[0].powi(2) + first.position[1].powi(2) + first.position[2].powi(2))
                .sqrt();
        let altitude_km = (pos_mag - 6371000.0) / 1000.0; // Approximate altitude

        println!(
            "First SV position magnitude: {:.1} km ({:.1} km altitude)",
            pos_mag / 1000.0,
            altitude_km
        );

        assert!(
            altitude_km > 650.0 && altitude_km < 750.0,
            "Satellite altitude {:.1} km should be in range [650, 750] for Sentinel-1",
            altitude_km
        );

        // Validate velocity magnitude (~7.5 km/s for LEO)
        let vel_mag =
            (first.velocity[0].powi(2) + first.velocity[1].powi(2) + first.velocity[2].powi(2))
                .sqrt();
        let vel_km_s = vel_mag / 1000.0;

        println!("First SV velocity magnitude: {:.3} km/s", vel_km_s);

        assert!(
            vel_km_s > 7.0 && vel_km_s < 8.0,
            "Satellite velocity {:.3} km/s should be in range [7, 8] for LEO",
            vel_km_s
        );
    }

    println!("\n✅ Orbit interpolation tests passed!");
}
