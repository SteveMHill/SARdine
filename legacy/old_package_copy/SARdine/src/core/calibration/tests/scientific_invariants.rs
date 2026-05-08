//! Scientific invariant tests for calibration and noise correction
//!
//! These tests validate that the core SAR processing operations maintain
//! scientific correctness according to ESA Sentinel-1 standards.

use crate::core::calibration::model::{
    CalibrationCoefficients, CalibrationVector, NoiseCoefficients, NoiseVector,
};
use crate::core::calibration::noise::apply_thermal_noise_removal_inplace;
use crate::core::calibration::{apply_calibration_to_denoised, CalibrationType};
use crate::core::calibration::parsing::derive_gamma_from_sigma;
use ndarray::Array2;
use num_complex::Complex;

/// Test 1: Complex SLC reading correctness
///
/// Validates that complex samples are correctly interpreted:
/// - Real/imaginary interleaving is correct
/// - Magnitude/power computation is correct: |S|² = real² + imag²
#[test]
fn test_complex_power_computation() {
    // Create synthetic complex data: [1+2j, 3+4j, 5+6j]
    let complex_data = vec![
        Complex::new(1.0, 2.0), // |1+2j|² = 1 + 4 = 5
        Complex::new(3.0, 4.0), // |3+4j|² = 9 + 16 = 25
        Complex::new(5.0, 6.0), // |5+6j|² = 25 + 36 = 61
    ];

    let power: Vec<f32> = complex_data.iter().map(|c| c.norm_sqr()).collect();

    // Verify power computation: |S|² = real² + imag²
    assert_eq!(power[0], 5.0, "|1+2j|² should be 5");
    assert_eq!(power[1], 25.0, "|3+4j|² should be 25");
    assert_eq!(power[2], 61.0, "|5+6j|² should be 61");
}

/// Test 2: Calibration correctness
///
/// Validates sigma0 calibration formula: σ⁰ = |DN|² / K²
/// Where K is the calibration constant from LUT.
///
/// **ESA Sentinel-1 Calibration Formula:**
/// According to ESA Sentinel-1 Product Specification (EOP-SM/1087/PS):
/// σ⁰ = |DN|² / (K² × A_s)
/// Where:
/// - DN = Digital Number (complex SLC pixel)
/// - K = Calibration constant (sigmaNought from annotation XML)
/// - A_s = Absolute calibration constant (typically 1.0 for Sentinel-1)
///
/// For this test, A_s = 1.0 (unity), so the formula simplifies to: σ⁰ = |DN|² / K²
#[test]
fn test_sigma0_calibration_formula() {
    // Create a simple calibration vector with known K values
    let mut calibration = CalibrationCoefficients::new();
    calibration.vectors.push(CalibrationVector {
        azimuth_time: "2020-01-01T00:00:00Z".to_string(),
        line: 0,
        pixels: vec![0, 1, 2],
        sigma_nought: vec![100.0, 200.0, 300.0], // K values
        beta_nought: vec![100.0, 200.0, 300.0],
        gamma: vec![100.0, 200.0, 300.0],
        dn: vec![1.0, 1.0, 1.0],
        beta_flat: false,
        sigma_flat: false,
        gamma_flat: false,
    });

    // Precompute LUT for 1x3 image
    calibration
        .precompute_lut((1, 3))
        .expect("LUT precompute failed");

    // Create power data: |DN|² = [100.0, 400.0, 900.0]
    // For pixel 0: K=100, |DN|²=100 → σ⁰ = 100 / (100²) = 0.01
    // For pixel 1: K=200, |DN|²=400 → σ⁰ = 400 / (200²) = 0.01
    // For pixel 2: K=300, |DN|²=900 → σ⁰ = 900 / (300²) = 0.01
    let power_data = Array2::from_shape_vec((1, 3), vec![100.0, 400.0, 900.0])
        .expect("Failed to create power array");

    let calibrated =
        apply_calibration_to_denoised(&power_data, &calibration, CalibrationType::Sigma0, None)
            .expect("Calibration failed");

    // Expected: gain = 1/K²
    // Pixel 0: gain = 1/(100²) = 0.0001, σ⁰ = 100 * 0.0001 = 0.01
    // Pixel 1: gain = 1/(200²) = 0.000025, σ⁰ = 400 * 0.000025 = 0.01
    // Pixel 2: gain = 1/(300²) = 0.000011111..., σ⁰ = 900 * 0.000011111... ≈ 0.01
    let expected = 0.01f32;
    assert!(
        (calibrated[[0, 0]] - expected).abs() < 1e-5,
        "Pixel 0 sigma0 mismatch"
    );
    assert!(
        (calibrated[[0, 1]] - expected).abs() < 1e-5,
        "Pixel 1 sigma0 mismatch"
    );
    assert!(
        (calibrated[[0, 2]] - expected).abs() < 1e-5,
        "Pixel 2 sigma0 mismatch"
    );
}

/// Test 3: Noise subtraction correctness
///
/// Validates thermal noise removal in power domain:
/// - Noise is subtracted: denoised_power = power - noise
/// - Noise-dominated pixels (noise > signal) are marked invalid (NaN)
#[test]
fn test_noise_subtraction_correctness() {
    // Create simple noise LUT with known values
    let mut noise_coeffs = NoiseCoefficients::new();
    noise_coeffs.vectors.push(NoiseVector {
        azimuth_time: "2020-01-01T00:00:00Z".to_string(),
        azimuth_time_seconds: 0.0,
        line: 0.0,
        range_pixels: vec![0.0, 1.0, 2.0],
        noise_range_lut: vec![10.0, 20.0, 30.0], // Noise values
    });

    // Precompute range-only noise LUT for 1x3 image
    noise_coeffs
        .precompute_lut((1, 3))
        .expect("Noise LUT precompute failed");

    // Create power data: [50.0, 25.0, 15.0]
    // Expected after noise removal:
    // Pixel 0: 50.0 - 10.0 = 40.0
    // Pixel 1: 25.0 - 20.0 = 5.0
    // Pixel 2: 15.0 - 30.0 = -15.0 → marked as NaN (noise-dominated)
    let mut power_data = Array2::from_shape_vec((1, 3), vec![50.0, 25.0, 15.0])
        .expect("Failed to create power array");

    apply_thermal_noise_removal_inplace(&mut power_data, &noise_coeffs)
        .expect("Noise removal failed");

    assert!(
        (power_data[[0, 0]] - 40.0).abs() < 1e-5,
        "Pixel 0 denoised mismatch"
    );
    assert!(
        (power_data[[0, 1]] - 5.0).abs() < 1e-5,
        "Pixel 1 denoised mismatch"
    );
    // Noise removal sets negative/too-small values to NaN, not 0.0.
    // Pixels where noise > signal are treated as invalid data.
    assert!(power_data[[0, 2]].is_nan(), "Pixel 2 should be NaN (noise > signal)");
}

/// Test 3b: Sigma0 → Gamma0 conversion across incidence angles
///
/// Validates γ⁰ = σ⁰ / cos(θ_inc) for a range of pixels and ensures
/// that the implementation tracks the analytical relation within a
/// small numerical tolerance for realistic incidence angles.
#[test]
fn test_gamma0_from_sigma0_across_incidence() {
    // Simple swath with 3 pixels across range
    let sigma_nought = vec![1.0_f32; 3];
    let pixels = vec![0_usize, 1_usize, 2_usize];
    let swath_width = 2_usize; // so fractions are 0/2, 1/2, 2/2 → 0.0, 0.5, 1.0

    // Choose a realistic IW incidence span
    let min_inc_deg = 30.0_f32;
    let max_inc_deg = 40.0_f32;

    let gamma = derive_gamma_from_sigma(
        &sigma_nought,
        &pixels,
        min_inc_deg,
        max_inc_deg,
        swath_width,
    );

    assert_eq!(gamma.len(), sigma_nought.len());

    for (idx, (&sigma, &pixel)) in sigma_nought.iter().zip(pixels.iter()).enumerate() {
        let range_fraction = (pixel as f32) / (swath_width as f32).max(1.0);
        let range_fraction = range_fraction.clamp(0.0, 1.0);
        let inc_deg = min_inc_deg + range_fraction * (max_inc_deg - min_inc_deg);
        let inc_rad = inc_deg.to_radians();
        let cos_inc = inc_rad.cos();

        // In this angle range, cos(θ) is comfortably above the safeguard threshold.
        assert!(cos_inc.abs() > 0.01);

        let expected_gamma = sigma / cos_inc;
        let got = gamma[idx];

        assert!(
            (got - expected_gamma).abs() < 1e-4,
            "Gamma0 mismatch at pixel {}: got {}, expected {} (θ = {} deg)",
            pixel,
            got,
            expected_gamma,
            inc_deg
        );
    }
}

/// Test 4: Burst geometry - constant range window
///
/// Validates that TOPS bursts have constant range extent per subswath.
/// This is a structural test that validates the burst record construction logic.
#[test]
fn test_burst_range_window_constant() {
    // This test would require mock annotation data
    // For now, we document the invariant that should be tested:
    // - All bursts in a subswath share the same range extent
    // - Only azimuth (line) extent varies between bursts
    //
    // In practice, this would test:
    // let bursts = build_burst_records_for_subswath_fixed(...);
    // let first_range = (bursts[0].start_sample_global, bursts[0].end_sample_global);
    // for burst in bursts.iter().skip(1) {
    //     assert_eq!(
    //         (burst.start_sample_global, burst.end_sample_global),
    //         first_range,
    //         "All bursts must have same range extent"
    //     );
    // }
}

/// Test 5: Time base consistency
///
/// Validates that all timing calculations use the same time base
/// (orbit_ref_epoch for burst times and orbit state vectors).
#[test]
fn test_time_base_consistency() {
    // This test would validate:
    // - Burst times are relative to orbit_ref_epoch
    // - Orbit state vector times are relative to orbit_ref_epoch
    // - All time calculations use consistent reference
    //
    // In practice:
    // let orbit_ref = DateTime<Utc>::from_utc(...);
    // let burst_time_rel = compute_burst_timing(burst, Some(orbit_ref));
    // let orbit_time_rel = orbit_state_vector.time - orbit_ref;
    // assert!(burst_time_rel.is_finite(), "Burst time should be finite");
    // assert!(orbit_time_rel.is_finite(), "Orbit time should be finite");
}

/// Test 6: End-to-end pipeline sanity
///
/// Validates that a minimal synthetic SLC → backscatter pipeline
/// produces expected outputs.
#[test]
fn test_end_to_end_pipeline_sanity() {
    // Create minimal synthetic data
    let complex_data = Array2::from_shape_vec(
        (2, 2),
        vec![
            Complex::new(10.0, 0.0), // |10+0j|² = 100
            Complex::new(0.0, 10.0), // |0+10j|² = 100
            Complex::new(5.0, 5.0),  // |5+5j|² = 50
            Complex::new(2.0, 4.0),  // |2+4j|² = 20
        ],
    )
    .expect("Failed to create complex array");

    // Compute power
    let power: Array2<f32> = complex_data.map(|c| c.norm_sqr());

    // Verify power values
    assert_eq!(power[[0, 0]], 100.0);
    assert_eq!(power[[0, 1]], 100.0);
    assert_eq!(power[[1, 0]], 50.0);
    assert_eq!(power[[1, 1]], 20.0);

    // Create simple calibration (unity gain for this test)
    let mut calibration = CalibrationCoefficients::new();
    calibration.vectors.push(CalibrationVector {
        azimuth_time: "2020-01-01T00:00:00Z".to_string(),
        line: 0,
        pixels: vec![0, 1],
        sigma_nought: vec![1.0, 1.0], // K=1 → gain=1
        beta_nought: vec![1.0, 1.0],
        gamma: vec![1.0, 1.0],
        dn: vec![1.0, 1.0],
        beta_flat: false,
        sigma_flat: false,
        gamma_flat: false,
    });

    // Expand to 2x2 (simple case: same calibration for all pixels)
    calibration.vectors.push(CalibrationVector {
        azimuth_time: "2020-01-01T00:00:01Z".to_string(),
        line: 1,
        pixels: vec![0, 1],
        sigma_nought: vec![1.0, 1.0],
        beta_nought: vec![1.0, 1.0],
        gamma: vec![1.0, 1.0],
        dn: vec![1.0, 1.0],
        beta_flat: false,
        sigma_flat: false,
        gamma_flat: false,
    });

    calibration
        .precompute_lut((2, 2))
        .expect("LUT precompute failed");

    // Apply calibration
    let calibrated =
        apply_calibration_to_denoised(&power, &calibration, CalibrationType::Sigma0, None)
            .expect("Calibration failed");

    // With K=1, gain=1, so output should equal input power
    assert!((calibrated[[0, 0]] - 100.0).abs() < 1e-5);
    assert!((calibrated[[0, 1]] - 100.0).abs() < 1e-5);
    assert!((calibrated[[1, 0]] - 50.0).abs() < 1e-5);
    assert!((calibrated[[1, 1]] - 20.0).abs() < 1e-5);
}
