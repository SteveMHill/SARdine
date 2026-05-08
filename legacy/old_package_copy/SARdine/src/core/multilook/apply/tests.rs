//! Tests for multilook processing

use super::*;
use ndarray::Array;
use num_complex::Complex;

#[test]
fn test_multilook_basic() {
    // Create a 4x4 test image
    let data = Array::from_shape_vec(
        (4, 4),
        vec![
            1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0,
        ],
    )
    .unwrap();

    let params = MultilookParams {
        range_looks: 2,
        azimuth_looks: 2,
        mode: MultilookMode::Intensity,
        preserve_power: false, // Use simple average for this test
        border_mode: BorderMode::Partial,
        include_partial: true,
    };

    let processor = MultilookProcessor::new(params);
    let (result, new_range, new_azimuth) = processor.apply_multilook(&data, 10.0, 10.0).unwrap();

    // Should be 2x2 output
    assert_eq!(result.dim(), (2, 2));

    // Check values (should be averages of 2x2 blocks)
    assert_eq!(result[[0, 0]], 3.5); // (1+2+5+6)/4
    assert_eq!(result[[0, 1]], 5.5); // (3+4+7+8)/4
    assert_eq!(result[[1, 0]], 11.5); // (9+10+13+14)/4
    assert_eq!(result[[1, 1]], 13.5); // (11+12+15+16)/4

    // Check new spacings
    assert_eq!(new_range, 20.0);
    assert_eq!(new_azimuth, 20.0);
}

#[test]
fn test_multilook_asymmetric() {
    // Test with different looks in each direction
    let data = Array::from_shape_vec((2, 4), vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0]).unwrap();

    let params = MultilookParams {
        range_looks: 2,
        azimuth_looks: 1,
        mode: MultilookMode::Intensity,
        preserve_power: false,
        border_mode: BorderMode::Partial,
        include_partial: true,
    };

    let processor = MultilookProcessor::new(params);
    let (result, _, _) = processor.apply_multilook(&data, 10.0, 10.0).unwrap();

    // Should be 2x2 output
    assert_eq!(result.dim(), (2, 2));

    // Check values
    assert_eq!(result[[0, 0]], 1.5); // (1+2)/2
    assert_eq!(result[[0, 1]], 3.5); // (3+4)/2
    assert_eq!(result[[1, 0]], 5.5); // (5+6)/2
    assert_eq!(result[[1, 1]], 7.5); // (7+8)/2
}

#[test]
fn test_enl_calculation() {
    // Create uniform data (should have high ENL)
    let uniform_data = Array2::<f32>::from_elem((100, 100), 1.0);
    let processor = MultilookProcessor::standard();
    let enl = processor.estimate_enl(&uniform_data);

    // Uniform data should have very high ENL (near infinity, but numerically limited)
    assert!(enl > 1000.0);
}

#[test]
fn test_complex_multilook_basic() {
    // Create a 4x4 complex test image with known phase pattern
    let mut data = Array2::<Complex<f32>>::zeros((4, 4));
    for i in 0..4 {
        for j in 0..4 {
            // Create a complex pattern: amplitude varies with position, constant phase
            let amp = (i * 4 + j + 1) as f32;
            let phase = std::f32::consts::PI / 4.0; // 45 degrees
            data[[i, j]] = Complex::new(amp * phase.cos(), amp * phase.sin());
        }
    }

    let params = MultilookParams {
        range_looks: 2,
        azimuth_looks: 2,
        mode: MultilookMode::Complex,
        preserve_power: false, // Use simple average for this test
        border_mode: BorderMode::Partial,
        include_partial: true,
    };

    let processor = MultilookProcessor::new(params);
    let (result, metadata) = processor
        .apply_multilook_complex(&data, 10.0, 10.0)
        .unwrap();

    // Should be 2x2 output
    assert_eq!(result.dim(), (2, 2));
    assert_eq!(metadata.output_dims, (2, 2));
    assert_eq!(metadata.mode, MultilookMode::Complex);

    // Check that phase is preserved (all samples have same phase)
    for i in 0..2 {
        for j in 0..2 {
            let z = result[[i, j]];
            let phase = z.arg();
            // Should be approximately PI/4 (allowing for numerical precision)
            assert!((phase - std::f32::consts::PI / 4.0).abs() < 0.01);
        }
    }

    // Check that amplitudes are averaged correctly
    // Block [0,0] should average samples (1,2,5,6)
    let expected_amp_00 = (1.0 + 2.0 + 5.0 + 6.0) / 4.0;
    let actual_amp_00 = result[[0, 0]].norm();
    assert!((actual_amp_00 - expected_amp_00).abs() < 0.01);
}

#[test]
fn test_complex_multilook_power_preservation() {
    // Test power preservation for complex multilook
    //
    // SCIENTIFIC CONTEXT:
    // For fully developed speckle (random phases), averaging N complex samples:
    //   - Reduces amplitude by 1/√N (coherent sum of random phasors)
    //   - Reduces intensity by 1/N
    //
    // With preserve_power=true, we scale by √N to restore intensity per pixel.
    // This means power_ratio ≈ 1.0 (total power preserved).
    //
    // For coherent data (constant phase), averaging preserves the signal,
    // and scaling by √N amplifies it by √N. This is correct behavior because
    // it preserves the expected power for incoherent targets in the scene.

    let mut data = Array2::<Complex<f32>>::zeros((4, 4));
    for i in 0..4 {
        for j in 0..4 {
            let amp = 2.0;
            let phase = (i as f32 + j as f32) * 0.1;
            data[[i, j]] = Complex::new(amp * phase.cos(), amp * phase.sin());
        }
    }

    let params = MultilookParams {
        range_looks: 2,
        azimuth_looks: 2,
        mode: MultilookMode::Complex,
        preserve_power: true,
        border_mode: BorderMode::Drop, // Use only complete blocks
        include_partial: false,
    };

    let processor = MultilookProcessor::new(params);
    let (result, metadata) = processor
        .apply_multilook_complex(&data, 10.0, 10.0)
        .unwrap();

    println!("Power ratio: {}", metadata.power_ratio);

    // With preserve_power=true, total power should be approximately preserved.
    // For this test with semi-coherent data (slowly varying phase), the power ratio
    // will be higher than 1.0 because the coherent signal is amplified by √N.
    // For truly random phase data, power_ratio ≈ 1.0.
    //
    // We allow a wide range here because the test data has coherent phase structure,
    // which gets amplified. The key check is that it's NOT 0.25 (the wrong behavior).
    assert!(
        metadata.power_ratio > 0.8,
        "Power ratio {} too low - power preservation not working (should be >= 0.8)",
        metadata.power_ratio
    );

    // Theoretical ENL should match looks
    assert_eq!(metadata.theoretical_enl, 4);

    // Output spacing should be double
    assert_eq!(metadata.output_spacing, (20.0, 20.0));
}

#[test]
fn test_complex_multilook_with_nans() {
    // Test NaN handling in complex multilook
    let mut data = Array2::<Complex<f32>>::zeros((4, 4));
    for i in 0..4 {
        for j in 0..4 {
            if i == 0 && j == 0 {
                data[[i, j]] = Complex::new(f32::NAN, f32::NAN);
            } else {
                data[[i, j]] = Complex::new(1.0, 0.5);
            }
        }
    }

    let params = MultilookParams {
        range_looks: 2,
        azimuth_looks: 2,
        mode: MultilookMode::Complex,
        preserve_power: false,
        border_mode: BorderMode::Partial,
        include_partial: true,
    };

    let processor = MultilookProcessor::new(params);
    let (result, metadata) = processor
        .apply_multilook_complex(&data, 10.0, 10.0)
        .unwrap();

    // Should handle NaN and produce valid output
    assert!(result[[0, 0]].re.is_finite());
    assert!(result[[0, 0]].im.is_finite());

    // Valid pixel fraction should reflect the NaN
    assert!(metadata.valid_pixel_fraction > 0.9);
}

// =============================================================================
// ML-1 FIX VERIFICATION TESTS
// =============================================================================
// These tests verify the scientifically correct multilook behavior after the
// ML-1 fix. The key requirement: intensity multilook output = mean of valid
// samples with NO fill_factor scaling for partial blocks.

#[test]
fn test_ml1_full_block_constant_field() {
    // ML-1 Test 1: Full block with constant values
    // Expected: output = input value (mean of identical values)
    let data = Array2::<f32>::from_elem((4, 4), 1.0);

    let params = MultilookParams {
        range_looks: 2,
        azimuth_looks: 2,
        mode: MultilookMode::Intensity,
        preserve_power: true, // Should have no effect after fix
        border_mode: BorderMode::Partial,
        include_partial: true,
    };

    let processor = MultilookProcessor::new(params);
    let (result, _, _) = processor.apply_multilook(&data, 10.0, 10.0).unwrap();

    // All output pixels should equal exactly 1.0 (mean of [1,1,1,1] = 1)
    for i in 0..2 {
        for j in 0..2 {
            assert!(
                (result[[i, j]] - 1.0).abs() < 1e-6,
                "Full block at [{},{}]: expected 1.0, got {}",
                i,
                j,
                result[[i, j]]
            );
        }
    }
}

#[test]
fn test_ml1_partial_block_no_bias() {
    // ML-1 Test 2: Partial block with NaN - CRITICAL regression test
    // This was the bug: preserve_power=true would multiply mean by fill_factor,
    // resulting in 0.75 instead of 1.0 for a 3/4-filled block.

    // Create 4x4 with one NaN in each 2x2 block corner
    let mut data = Array2::<f32>::from_elem((4, 4), 1.0);
    data[[0, 0]] = f32::NAN; // Block [0,0] has 3 valid pixels
    data[[1, 3]] = f32::NAN; // Block [0,1] has 3 valid pixels
    data[[2, 0]] = f32::NAN; // Block [1,0] has 3 valid pixels
    data[[3, 3]] = f32::NAN; // Block [1,1] has 3 valid pixels

    let params = MultilookParams {
        range_looks: 2,
        azimuth_looks: 2,
        mode: MultilookMode::Intensity,
        preserve_power: true, // Would cause bug before fix
        border_mode: BorderMode::Partial,
        include_partial: true,
    };

    let processor = MultilookProcessor::new(params);
    let (result, _, _) = processor.apply_multilook(&data, 10.0, 10.0).unwrap();

    // AFTER ML-1 FIX: All output should be 1.0 (mean of [1,1,1] = 1)
    // BEFORE FIX (BUG): Output was 0.75 (mean * 3/4 fill_factor)
    for i in 0..2 {
        for j in 0..2 {
            assert!(
                (result[[i, j]] - 1.0).abs() < 1e-6,
                "ML-1 regression: partial block at [{},{}]: expected 1.0 (mean), got {} (was 0.75 before fix)",
                i, j, result[[i, j]]
            );
        }
    }
}

#[test]
fn test_ml1_preserve_power_flag_no_effect() {
    // ML-1 Test 3: preserve_power flag should have NO effect for intensity multilook
    let mut data = Array2::<f32>::from_elem((4, 4), 2.0);
    data[[0, 0]] = f32::NAN; // Partial block

    let params_with = MultilookParams {
        range_looks: 2,
        azimuth_looks: 2,
        mode: MultilookMode::Intensity,
        preserve_power: true,
        border_mode: BorderMode::Partial,
        include_partial: true,
    };

    let params_without = MultilookParams {
        range_looks: 2,
        azimuth_looks: 2,
        mode: MultilookMode::Intensity,
        preserve_power: false,
        border_mode: BorderMode::Partial,
        include_partial: true,
    };

    let processor_with = MultilookProcessor::new(params_with);
    let processor_without = MultilookProcessor::new(params_without);

    let (result_with, _, _) = processor_with.apply_multilook(&data, 10.0, 10.0).unwrap();
    let (result_without, _, _) = processor_without
        .apply_multilook(&data, 10.0, 10.0)
        .unwrap();

    // Both should produce identical results
    for i in 0..2 {
        for j in 0..2 {
            let diff = (result_with[[i, j]] - result_without[[i, j]]).abs();
            assert!(
                diff < 1e-6,
                "preserve_power should have no effect: [{},{}] with={} without={} diff={}",
                i,
                j,
                result_with[[i, j]],
                result_without[[i, j]],
                diff
            );
        }
    }
}

#[test]
fn test_ml1_output_equals_mean_property() {
    // ML-1 Test 4: Property test - output always equals mean of valid inputs
    use std::collections::HashMap;

    // Create random-ish test data with some NaNs
    let mut data = Array2::<f32>::zeros((6, 8));
    let mut rng_state = 12345u64;
    for i in 0..6 {
        for j in 0..8 {
            rng_state = rng_state.wrapping_mul(6364136223846793005).wrapping_add(1);
            let val = (rng_state % 1000) as f32 / 100.0; // 0.0 to 9.99
            if rng_state % 7 == 0 {
                data[[i, j]] = f32::NAN;
            } else {
                data[[i, j]] = val;
            }
        }
    }

    let params = MultilookParams {
        range_looks: 2,
        azimuth_looks: 3,
        mode: MultilookMode::Intensity,
        preserve_power: true,
        border_mode: BorderMode::Partial,
        include_partial: true,
    };

    let processor = MultilookProcessor::new(params);
    let (result, _, _) = processor.apply_multilook(&data, 10.0, 10.0).unwrap();

    // Manually compute expected means for each block
    let out_rows = (6 + 3 - 1) / 3; // 2
    let out_cols = (8 + 2 - 1) / 2; // 4

    for out_row in 0..out_rows {
        for out_col in 0..out_cols {
            let mut sum = 0.0f32;
            let mut count = 0;

            for az in 0..3 {
                let in_row = out_row * 3 + az;
                if in_row >= 6 {
                    break;
                }
                for rg in 0..2 {
                    let in_col = out_col * 2 + rg;
                    if in_col >= 8 {
                        break;
                    }
                    let val = data[[in_row, in_col]];
                    if val.is_finite() {
                        sum += val;
                        count += 1;
                    }
                }
            }

            let expected = if count > 0 {
                sum / count as f32
            } else {
                f32::NAN
            };
            let actual = result[[out_row, out_col]];

            if expected.is_finite() {
                assert!(
                    (actual - expected).abs() < 1e-5,
                    "Block [{},{}]: expected mean={}, got {}",
                    out_row,
                    out_col,
                    expected,
                    actual
                );
            } else {
                assert!(
                    actual.is_nan(),
                    "Block [{},{}]: expected NaN",
                    out_row,
                    out_col
                );
            }
        }
    }
}

#[test]
fn test_ml1_separable_method_consistency() {
    // ML-1 Test 5: Separable method should match standard method
    let mut data = Array2::<f32>::from_elem((6, 8), 5.0);
    data[[0, 0]] = f32::NAN;
    data[[2, 4]] = f32::NAN;
    data[[5, 7]] = f32::NAN;

    let params = MultilookParams {
        range_looks: 2,
        azimuth_looks: 2,
        mode: MultilookMode::Intensity,
        preserve_power: true,
        border_mode: BorderMode::Partial,
        include_partial: true,
    };

    let processor = MultilookProcessor::new(params);
    let (result_standard, _, _) = processor.apply_multilook(&data, 10.0, 10.0).unwrap();
    let (result_separable, _, _) = processor
        .apply_multilook_separable(&data, 10.0, 10.0)
        .unwrap();

    // Both methods should produce same output within numerical tolerance
    assert_eq!(result_standard.dim(), result_separable.dim());
    for i in 0..result_standard.nrows() {
        for j in 0..result_standard.ncols() {
            let std = result_standard[[i, j]];
            let sep = result_separable[[i, j]];
            if std.is_finite() && sep.is_finite() {
                let diff = (std - sep).abs();
                assert!(
                    diff < 1e-4,
                    "Method mismatch at [{},{}]: standard={}, separable={}, diff={}",
                    i,
                    j,
                    std,
                    sep,
                    diff
                );
            }
        }
    }
}

#[test]
fn test_ml1_edge_blocks_no_power_reduction() {
    // ML-1 Test 6: Edge blocks (partial due to image size) should not have reduced power
    // Create 5x5 image with 2x2 multilook -> 3x3 output with edge partials
    let data = Array2::<f32>::from_elem((5, 5), 4.0);

    let params = MultilookParams {
        range_looks: 2,
        azimuth_looks: 2,
        mode: MultilookMode::Intensity,
        preserve_power: true,
        border_mode: BorderMode::Partial,
        include_partial: true,
    };

    let processor = MultilookProcessor::new(params);
    let (result, _, _) = processor.apply_multilook(&data, 10.0, 10.0).unwrap();

    // Output should be 3x3
    assert_eq!(result.dim(), (3, 3));

    // ALL blocks should have value 4.0 (the input constant)
    // - Interior block [0,0]: 4 samples -> mean = 4.0
    // - Edge block [0,2]: 2 samples -> mean = 4.0 (NOT 2.0 with old bug)
    // - Edge block [2,0]: 2 samples -> mean = 4.0
    // - Corner block [2,2]: 1 sample -> mean = 4.0 (NOT 1.0 with old bug)
    for i in 0..3 {
        for j in 0..3 {
            assert!(
                (result[[i, j]] - 4.0).abs() < 1e-6,
                "Edge block [{},{}] power reduction: expected 4.0, got {}",
                i,
                j,
                result[[i, j]]
            );
        }
    }
}
