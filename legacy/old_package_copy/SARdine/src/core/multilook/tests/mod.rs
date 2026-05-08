use super::*;
use ndarray::{Array, Array2};

#[test]
fn includes_partial_blocks_by_default() {
    let data = Array2::<f32>::from_elem((101, 103), 1.0);
    let processor = MultilookProcessor::new(MultilookParams {
        range_looks: 8,
        azimuth_looks: 6,
        ..MultilookParams::default()
    });

    let (result, _, _) = processor
        .apply_multilook(&data, 10.0, 10.0)
        .expect("multilook should succeed");

    // ceil-style output size with partial blocks retained
    assert_eq!(result.dim(), (17, 13));
    assert!(result.iter().all(|v| (v - 1.0).abs() < 1e-6));
}

#[test]
fn rejects_zero_looks() {
    let data = Array2::<f32>::ones((4, 4));
    let processor = MultilookProcessor::new(MultilookParams {
        range_looks: 0,
        ..MultilookParams::default()
    });

    let result = processor.apply_multilook(&data, 10.0, 10.0);
    assert!(result.is_err());
}

#[test]
fn nan_only_block_is_marked_nan() {
    let mut data = Array2::<f32>::ones((4, 4));
    for i in 0..2 {
        for j in 0..2 {
            data[[i, j]] = f32::NAN;
        }
    }

    let processor = MultilookProcessor::new(MultilookParams {
        range_looks: 2,
        azimuth_looks: 2,
        ..MultilookParams::default()
    });
    let (result, _, _) = processor
        .apply_multilook(&data, 10.0, 10.0)
        .expect("multilook should succeed");

    assert!(result[[0, 0]].is_nan());
    assert!(result[[0, 1]].is_finite());
    assert!(result[[1, 0]].is_finite());
    assert!(result[[1, 1]].is_finite());
}

#[test]
fn filtered_matches_standard_output() {
    let data = Array::from_shape_vec((6, 8), (0..48).map(|x| x as f32).collect()).unwrap();

    let processor = MultilookProcessor::new(MultilookParams {
        range_looks: 3,
        azimuth_looks: 2,
        ..MultilookParams::default()
    });

    let (standard, r1, a1) = processor
        .apply_multilook(&data, 10.0, 20.0)
        .expect("standard multilook should succeed");
    let (filtered, r2, a2) = processor
        .apply_multilook_filtered(&data, 10.0, 20.0)
        .expect("filtered multilook should succeed");

    assert_eq!(standard.dim(), filtered.dim());
    assert_eq!(r1, r2);
    assert_eq!(a1, a2);

    for (a, b) in standard.iter().zip(filtered.iter()) {
        assert!((a - b).abs() < 1e-5);
    }
}

#[test]
fn enl_on_uniform_is_high() {
    let data = Array2::<f32>::from_elem((20, 20), 5.0);
    let processor = MultilookProcessor::standard();
    let enl = processor.estimate_enl(&data);

    assert!(enl > 1e5);
}

#[test]
fn integral_nanaware_respects_include_partial() {
    let data = Array2::<f32>::ones((5, 5));
    let processor = MultilookProcessor::new(MultilookParams {
        range_looks: 2,
        azimuth_looks: 2,
        ..MultilookParams::default()
    });

    let (compact, _, _) = processor
        .apply_multilook_integral_nanaware(&data, 10.0, 10.0, false)
        .expect("compact multilook should succeed");
    assert_eq!(compact.dim(), (2, 2));

    let (partial, _, _) = processor
        .apply_multilook_integral_nanaware(&data, 10.0, 10.0, true)
        .expect("partial multilook should succeed");
    assert_eq!(partial.dim(), (3, 3));
}
