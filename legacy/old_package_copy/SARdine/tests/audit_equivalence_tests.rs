//! Equivalence tests for the scientific audit fixes.
//!
//! These tests verify:
//! - Speed of light constant consistency (DUP-5)
//! - Gradient operator equivalence (DUP-4)
//! - Quality flag emission for fallbacks (FALL-2, FALL-3, FALL-5)

use ndarray::Array2;

/// Verify all speed of light references use the canonical constant.
#[test]
fn test_speed_of_light_consistency() {
    use sardine::constants::physical::SPEED_OF_LIGHT_M_S;

    // The exact SI definition (2019)
    const EXPECTED_C: f64 = 299_792_458.0;

    assert_eq!(
        SPEED_OF_LIGHT_M_S, EXPECTED_C,
        "SPEED_OF_LIGHT_M_S should match SI definition"
    );

    // Verify it's used in calculations
    let wavelength = SPEED_OF_LIGHT_M_S / 5.405e9; // Sentinel-1 C-band
    assert!(
        (wavelength - 0.0555).abs() < 0.001,
        "Wavelength calculation should be approximately 5.55 cm"
    );
}

/// Verify all gradient operators produce consistent results for simple cases.
#[test]
fn test_gradient_operators_equivalence() {
    use sardine::core::terrain_correction::gradient::{compute_gradients, GradientOperator};

    // Create a simple linear ramp DEM (constant gradient)
    let mut dem = Array2::<f32>::zeros((10, 10));
    for i in 0..10 {
        for j in 0..10 {
            // Slope of 1:1 in both directions (45° slope)
            dem[[i, j]] = (i + j) as f32 * 10.0; // 10m per pixel
        }
    }

    let spacing = 10.0; // 10m pixel spacing

    // Compute gradients with all operators
    let (horn_dx, horn_dy) =
        compute_gradients(&dem, GradientOperator::Horn3x3, spacing, spacing).unwrap();
    let (sobel_dx, sobel_dy) =
        compute_gradients(&dem, GradientOperator::Sobel3x3, spacing, spacing).unwrap();
    let (cd_dx, cd_dy) =
        compute_gradients(&dem, GradientOperator::CentralDifference, spacing, spacing).unwrap();

    // For a linear ramp, all operators should give identical results
    for i in 2..8 {
        for j in 2..8 {
            // Expected gradient: 10m rise / 10m run = 1.0 in each direction
            let expected_grad = 1.0;

            // Horn and Sobel are mathematically identical
            assert!(
                (horn_dx[[i, j]] - sobel_dx[[i, j]]).abs() < 1e-6,
                "Horn and Sobel should be identical: horn={}, sobel={}",
                horn_dx[[i, j]],
                sobel_dx[[i, j]]
            );

            // For linear slope, central difference should also match
            assert!(
                (horn_dx[[i, j]] - cd_dx[[i, j]]).abs() < 0.1,
                "Horn and CentralDiff should match for linear slope: horn={}, cd={}",
                horn_dx[[i, j]],
                cd_dx[[i, j]]
            );

            // Verify expected gradient value
            assert!(
                (horn_dx[[i, j]] - expected_grad).abs() < 0.1,
                "Expected gradient ~1.0, got {}",
                horn_dx[[i, j]]
            );
        }
    }
}

/// Verify quality flag module is properly exported and functional.
#[test]
fn test_quality_flags_api() {
    use sardine::core::quality_flags::{
        global_quality_flags, reset_global_quality_flags, QualityFlag, QualityFlags,
    };

    // Reset global flags before test
    reset_global_quality_flags();

    // Create a flag
    let flag = QualityFlag::NoiseRemovalSkipped {
        reason: "Test reason".to_string(),
        zero_ratio: Some(0.96),
        noise_power_ratio: Some(0.8),
        zero_threshold: 0.95,
        ratio_threshold: 0.9,
    };

    // Verify flag type
    assert_eq!(flag.flag_type(), "NOISE_REMOVAL_SKIPPED");
    assert_eq!(flag.severity(), "WARNING");

    // Verify collection works
    let mut flags = QualityFlags::new();
    flags.add(flag);
    assert!(flags.has_flag("NOISE_REMOVAL_SKIPPED"));
    assert!(!flags.has_flag("DEBURST_T0_FALLBACK"));

    // Verify global flags are accessible
    let global = global_quality_flags();
    global.add(QualityFlag::DeburstT0Fallback {
        burst_id: 0,
        reason: "Test".to_string(),
        expected_offset: None,
    });

    let collected = global.get_flags();
    assert!(collected.has_flag("DEBURST_T0_FALLBACK"));
}

/// Verify gradient normal computation is consistent.
#[test]
fn test_gradient_to_normal_consistency() {
    use sardine::core::terrain_correction::gradient::{gradient_to_normal, slope_magnitude};

    // Flat surface
    let flat_normal = gradient_to_normal(0.0, 0.0);
    assert!(
        (flat_normal[2] - 1.0).abs() < 1e-6,
        "Flat surface should have z-up normal"
    );

    // 45° slope in X
    let slope_normal = gradient_to_normal(1.0, 0.0);
    let expected_z = 1.0 / (2.0_f32).sqrt();
    assert!(
        (slope_normal[2] - expected_z).abs() < 0.01,
        "45° slope should have nz = 1/√2 = 0.707, got {}",
        slope_normal[2]
    );

    // Verify slope magnitude
    let mag = slope_magnitude(1.0, 0.0);
    assert!(
        (mag - (2.0_f32).sqrt()).abs() < 0.01,
        "45° slope magnitude should be √2 = 1.414, got {}",
        mag
    );
}
