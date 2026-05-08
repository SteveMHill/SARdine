use crate::core::deburst::geometry::{build_line_timing_with_offset, eval_dc_fm, BurstInfo};
use crate::core::deburst::quality;
use crate::core::deburst::burst_ops::{compute_pairwise_weights, overlap_weight, valid_window, w_cos2};
use crate::core::deburst::{DeburstConfig, TopSarDeburstProcessor};
use crate::types::SarComplex;
use ndarray::Array2;

#[test]
fn topsar_deburst_processor_runs() {
    let burst_info = vec![BurstInfo {
        burst_id: 0,
        start_line: 0,
        end_line: 499,
        start_sample: 0,
        end_sample: 999,
        azimuth_time: "2020-01-03T17:08:16.618328".to_string(),
        sensing_time: "2020-01-03T17:08:17.623236".to_string(),
        first_valid_sample: vec![100; 500],
        last_valid_sample: vec![900; 500],
        byte_offset: 109035,
        azimuth_fm_rate: 2000.0,
        azimuth_steering_rate: 0.0015,
        slant_range_time: 0.006,
        doppler_centroid: 0.0,
        azimuth_bandwidth: 320.0,
        range_sampling_rate: 64000000.0,
        range_pixel_spacing: 2.329562,
        azimuth_pixel_spacing: 14.059906,
        azimuth_time_interval: 0.00021,
        dc_polynomial: vec![0.0, 0.0],
        fm_polynomial: vec![2000.0, 0.0],
        dc_range_poly: None,
        fm_range_poly: None,
        dc_polynomial_t0: None,
        fm_polynomial_t0: None,
        burst_reference_time_seconds: None,
    }];

    let config = DeburstConfig::default();
    let processor = TopSarDeburstProcessor::new(burst_info, config, 7500.0);
    let test_data = Array2::zeros((1000, 1000));

    let result = processor.deburst_topsar(&test_data);
    assert!(result.is_ok());

    let enhanced_result = processor.deburst_topsar_enhanced(&test_data);
    assert!(enhanced_result.is_ok());
    let result = enhanced_result.unwrap();
    assert!(result.power_ratio >= 0.0);
    assert_eq!(result.image.dim().0, 500);
}

#[test]
fn enhanced_timing_functions() {
    let timings = build_line_timing_with_offset(10, 0.001, 0.0);
    assert_eq!(timings.len(), 10);
    assert!((timings[0].t_az + 0.0045).abs() < 1e-6);

    let dc_poly = vec![100.0, 0.0];
    let fm_poly = vec![2000.0, 0.0];
    let (dc, fm) = eval_dc_fm(0.001, &dc_poly, &fm_poly, None, None);
    assert!((dc - 100.0).abs() < 1e-6);
    assert!((fm - 2000.0).abs() < 1e-6);

    let first_valid = vec![10, 20, 30];
    let last_valid = vec![90, 80, 70];
    let (start, end) = valid_window(1, &first_valid, &last_valid, 100);
    assert_eq!(start, 20);
    assert_eq!(end, 81);
}

#[test]
fn complementary_blending_weights() {
    assert!((w_cos2(0.0) - 1.0).abs() < 1e-6);
    assert!((w_cos2(1.0) - 0.0).abs() < 1e-6);
    assert!((w_cos2(0.5) - 0.5).abs() < 1e-6);

    let w1 = overlap_weight(10, 100, 20, true);
    let w2 = overlap_weight(10, 100, 20, false);
    assert!((w1 + w2 - 1.0).abs() < 1e-6);
}

#[test]
fn pairwise_complementary_weights() {
    for overlap_len in [10, 20, 50, 100] {
        for line_in_overlap in 0..overlap_len {
            let (w_current, w_next) = compute_pairwise_weights(overlap_len, line_in_overlap);
            assert!((w_current + w_next - 1.0).abs() < 1e-6);
            assert!(w_current >= 0.0 && w_current <= 1.0);
            assert!(w_next >= 0.0 && w_next <= 1.0);
        }
    }

    let (w_start_current, w_start_next) = compute_pairwise_weights(100, 0);
    assert!((w_start_current - 1.0).abs() < 1e-6);
    assert!((w_start_next - 0.0).abs() < 1e-6);

    let (w_end_current, w_end_next) = compute_pairwise_weights(100, 99);
    assert!((w_end_current - 0.0).abs() < 1e-6);
    assert!((w_end_next - 1.0).abs() < 1e-6);
}

#[test]
fn burst_power_diagnostics_are_sane() {
    let burst_info = vec![
        BurstInfo {
            burst_id: 0,
            start_line: 0,
            end_line: 99,
            start_sample: 0,
            end_sample: 99,
            azimuth_time: "2020-01-03T17:08:16.618328".to_string(),
            sensing_time: "2020-01-03T17:08:17.623236".to_string(),
            first_valid_sample: vec![10; 100],
            last_valid_sample: vec![90; 100],
            byte_offset: 0,
            azimuth_fm_rate: 2000.0,
            azimuth_steering_rate: 0.0015,
            slant_range_time: 0.006,
            doppler_centroid: 0.0,
            azimuth_bandwidth: 320.0,
            range_sampling_rate: 64000000.0,
            range_pixel_spacing: 2.329562,
            azimuth_pixel_spacing: 14.059906,
            azimuth_time_interval: 0.00021,
            dc_polynomial: vec![0.0],
            fm_polynomial: vec![2000.0],
            dc_range_poly: None,
            fm_range_poly: None,
            dc_polynomial_t0: None,
            fm_polynomial_t0: None,
            burst_reference_time_seconds: None,
        },
        BurstInfo {
            burst_id: 1,
            start_line: 100,
            end_line: 199,
            start_sample: 0,
            end_sample: 99,
            azimuth_time: "2020-01-03T17:08:17.618328".to_string(),
            sensing_time: "2020-01-03T17:08:18.623236".to_string(),
            first_valid_sample: vec![10; 100],
            last_valid_sample: vec![90; 100],
            byte_offset: 10000,
            azimuth_fm_rate: 2000.0,
            azimuth_steering_rate: 0.0015,
            slant_range_time: 0.006,
            doppler_centroid: 0.0,
            azimuth_bandwidth: 320.0,
            range_sampling_rate: 64000000.0,
            range_pixel_spacing: 2.329562,
            azimuth_pixel_spacing: 14.059906,
            azimuth_time_interval: 0.00021,
            dc_polynomial: vec![0.0],
            fm_polynomial: vec![2000.0],
            dc_range_poly: None,
            fm_range_poly: None,
            dc_polynomial_t0: None,
            fm_polynomial_t0: None,
            burst_reference_time_seconds: None,
        },
    ];

    let mut test_data = Array2::zeros((200, 100));
    for i in 0..200 {
        for j in 0..100 {
            test_data[[i, j]] = SarComplex::new(1.0, 1.0);
        }
    }

    let burst_info_clone = burst_info.clone();
    let _processor = TopSarDeburstProcessor::new(burst_info, DeburstConfig::default(), 7500.0);

    let diagnostics = quality::calculate_burst_power_diagnostics(&burst_info_clone, &test_data);
    assert_eq!(diagnostics.len(), 2);
    for (_, _, mean_power) in &diagnostics {
        assert!((*mean_power - 2.0).abs() < 0.1);
    }
}

#[test]
fn range_dependent_deramp_default_enabled() {
    let config = DeburstConfig::default();
    assert!(config.use_range_dependent_deramp);
}
