use super::{SpeckleFilter, SpeckleFilterParams, SpeckleFilterType};
use ndarray::Array2;

fn small_window_params() -> SpeckleFilterParams {
    SpeckleFilterParams {
        window_size: 3,
        ..SpeckleFilterParams::default()
    }
}

#[test]
fn validate_rejects_invalid_window() {
    let params = SpeckleFilterParams {
        window_size: 0,
        ..SpeckleFilterParams::default()
    };
    assert!(params.validate_window().is_err());

    let params = SpeckleFilterParams {
        window_size: 4,
        ..SpeckleFilterParams::default()
    };
    assert!(params.validate_window().is_err());
}

#[test]
fn mean_filter_preserves_shape() {
    let image = Array2::from_elem((6, 6), 10.0f32);
    let filter = SpeckleFilter::with_params(small_window_params());

    let result = filter
        .apply_filter(&image, SpeckleFilterType::Mean)
        .expect("mean filter should succeed");

    assert_eq!(result.dim(), image.dim());
}

#[test]
fn tiled_processing_matches_standard_for_mean() {
    let mut image = Array2::zeros((8, 8));
    for i in 0..8 {
        for j in 0..8 {
            image[[i, j]] = (i * 8 + j) as f32;
        }
    }

    let filter = SpeckleFilter::with_params(small_window_params());

    let standard = filter
        .apply_filter(&image, SpeckleFilterType::Mean)
        .expect("standard mean filter should succeed");

    // Use a small tile size to exercise the tiled path
    let tiled = filter
        .apply_filter_tiled(&image, SpeckleFilterType::Mean, Some(4))
        .expect("tiled mean filter should succeed");

    assert_eq!(standard.dim(), tiled.dim());

    for ((a, b), idx) in standard
        .iter()
        .zip(tiled.iter())
        .zip(standard.indexed_iter().map(|(idx, _)| idx))
    {
        let diff = (a - b).abs();
        assert!(
            diff < 1e-5,
            "difference {} at {:?} (standard {}, tiled {})",
            diff,
            idx,
            a,
            b
        );
    }
}
