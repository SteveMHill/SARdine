//! Tests for SLC reader module

use super::*;
use chrono::{TimeZone, Utc};

#[test]
fn test_burst_range_window_constant_across_bursts() {
    // This test verifies Fix A2: burst range window does NOT advance per burst
    // In TOPS/IW, bursts tile in azimuth only; range window is constant.

    // Create a mock subswath with known geometry
    let subswath = crate::types::SubSwath {
        id: "IW1".to_string(),
        burst_count: 3,
        lines_per_burst: 1000,
        first_line_global: 0,
        last_line_global: 2999,    // 3000 lines total
        first_sample_global: 100,  // Constant range start
        last_sample_global: 25000, // Constant range end
        azimuth_samples: 3000,
        range_samples: 24901,
        full_range_samples: 24901,
        slant_range_time: 0.005,
        range_pixel_spacing: 2.3,
        azimuth_pixel_spacing: 14.0,
        burst_duration: 2.76,
        near_range_m: 800000.0,
        azimuth_time_interval: Some(0.002056),
        prf_hz: Some(486.0),
        dc_polynomial: None,
        dc_polynomial_t0: None,
        valid_first_line: None,
        valid_last_line: None,
        valid_first_sample: None,
        valid_last_sample: None,
        fm_rate_estimates: None,
    };

    // Create mock annotation with 3 bursts
    let mock_bursts = vec![
        create_mock_burst(0, "2020-01-01T12:00:00.000000Z"),
        create_mock_burst(1, "2020-01-01T12:00:02.000000Z"),
        create_mock_burst(2, "2020-01-01T12:00:04.000000Z"),
    ];

    let annotation = create_mock_annotation(mock_bursts, 1000); // 1000 lines per burst

    let orbit_epoch = Some(Utc.with_ymd_and_hms(2020, 1, 1, 12, 0, 0).unwrap());

    let records = burst_records::build_burst_records_for_subswath_fixed(
        &annotation,
        "IW1",
        &subswath,
        orbit_epoch,
    );

    assert_eq!(records.len(), 3, "Should have 3 bursts");

    // Verify ALL bursts have the SAME range window
    for (i, rec) in records.iter().enumerate() {
        assert_eq!(
            rec.start_sample_global, 100,
            "Burst {} should have constant start_sample_global=100, got {}",
            i, rec.start_sample_global
        );
        assert_eq!(
            rec.end_sample_global, 25000,
            "Burst {} should have constant end_sample_global=25000, got {}",
            i, rec.end_sample_global
        );
    }

    // Verify azimuth lines DO advance per burst
    assert_eq!(records[0].first_line_global, 0);
    assert_eq!(records[0].last_line_global, 999);
    assert_eq!(records[1].first_line_global, 1000);
    assert_eq!(records[1].last_line_global, 1999);
    assert_eq!(records[2].first_line_global, 2000);
    assert_eq!(records[2].last_line_global, 2999);
}

#[test]
fn test_orbit_relative_timing_preferred() {
    let subswath = crate::types::SubSwath {
        id: "IW1".to_string(),
        burst_count: 1,
        lines_per_burst: 1000,
        first_line_global: 0,
        last_line_global: 999,
        first_sample_global: 0,
        last_sample_global: 1000,
        azimuth_samples: 1000,
        range_samples: 1001,
        full_range_samples: 1001,
        slant_range_time: 0.005,
        range_pixel_spacing: 2.3,
        azimuth_pixel_spacing: 14.0,
        burst_duration: 2.76,
        near_range_m: 800000.0,
        azimuth_time_interval: Some(0.002056),
        prf_hz: Some(486.0),
        dc_polynomial: None,
        dc_polynomial_t0: None,
        valid_first_line: None,
        valid_last_line: None,
        valid_first_sample: None,
        valid_last_sample: None,
        fm_rate_estimates: None,
    };

    // Burst at 12:00:05 with orbit epoch at 12:00:00
    // Expected relative time: 5.5 seconds
    let mock_bursts = vec![create_mock_burst(0, "2020-01-01T12:00:05.500000Z")];

    let annotation = create_mock_annotation(mock_bursts, 1000);
    let orbit_epoch = Some(Utc.with_ymd_and_hms(2020, 1, 1, 12, 0, 0).unwrap());

    let records = burst_records::build_burst_records_for_subswath_fixed(
        &annotation,
        "IW1",
        &subswath,
        orbit_epoch,
    );

    assert_eq!(records.len(), 1);

    let rel_time = records[0]
        .azimuth_time_rel_orbit
        .expect("Should have relative time");
    assert!(
        (rel_time - 5.5).abs() < 1e-6,
        "Expected 5.5s relative to orbit epoch, got {}",
        rel_time
    );
}

// Helper functions to create mock annotation structures
fn create_mock_burst(idx: usize, azimuth_time: &str) -> crate::io::annotation::Burst {
    crate::io::annotation::Burst {
        azimuth_time: Some(azimuth_time.to_string()),
        sensing_time: None,
        azimuth_anx_time: (idx as f64) * 2.0, // Fallback timing
        byte_offset: Some(idx as u64 * 1000000),
        first_valid_sample: vec![10, 10, 10],
        last_valid_sample: vec![990, 990, 990],
    }
}

fn create_mock_annotation(
    bursts: Vec<crate::io::annotation::Burst>,
    lines_per_burst: u32,
) -> crate::io::annotation::ProductRoot {
    crate::io::annotation::ProductRoot {
        ads_header: None,
        general_annotation: None,
        image_annotation: None,
        doppler_centroid: None,
        antenna_pattern: None,
        swath_timing: Some(crate::io::annotation::SwathTiming {
            lines_per_burst: Some(lines_per_burst),
            samples_per_burst: Some(24901),
            burst_list: Some(crate::io::annotation::BurstList {
                count: Some(bursts.len() as u32),
                bursts: Some(bursts),
            }),
        }),
        geolocation_grid: None,
        coordinate_conversion: None,
        quality_information: None,
        orbit_list: None,
    }
}
