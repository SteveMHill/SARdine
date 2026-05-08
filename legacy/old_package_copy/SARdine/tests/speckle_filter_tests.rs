use ndarray::Array2;
use sardine::core::precision_standards::DeterministicRng;
use sardine::core::speckle_filter::{SpeckleFilter, SpeckleFilterParams, SpeckleFilterType};
use sardine::types::SarResult;

fn make_positive_image(height: usize, width: usize, seed: u64) -> Array2<f32> {
    let mut rng = DeterministicRng::with_seed(seed);
    let mut data = Array2::zeros((height, width));
    for pixel in data.iter_mut() {
        // Keep values strictly positive to avoid log-domain corner cases.
        *pixel = 1.0 + rng.gen_f32();
    }
    data
}

fn variance(image: &Array2<f32>) -> f32 {
    let mut sum = 0.0f32;
    let mut sum_sq = 0.0f32;
    let mut count = 0f32;
    for &v in image.iter() {
        if v.is_finite() && v >= 0.0 {
            sum += v;
            sum_sq += v * v;
            count += 1.0;
        }
    }
    if count <= 1.0 {
        return 0.0;
    }
    let mean = sum / count;
    (sum_sq / count - mean * mean).max(0.0)
}

#[test]
fn lee_filter_reduces_variance() -> SarResult<()> {
    let image = make_positive_image(64, 64, 0xA5A5_1234);
    let params = SpeckleFilterParams {
        window_size: 5,
        num_looks: 4.0,
        ..Default::default()
    };
    let filter = SpeckleFilter::with_params(params);

    let before = variance(&image);
    let after = variance(&filter.apply_filter(&image, SpeckleFilterType::Lee)?);

    assert!(
        after <= before,
        "Lee filter should not increase variance (before: {before}, after: {after})"
    );
    Ok(())
}

#[test]
fn tiled_and_direct_gamma_map_match() -> SarResult<()> {
    let image = make_positive_image(48, 48, 0xBEEF_BEEF);
    let params = SpeckleFilterParams {
        window_size: 5,
        num_looks: 3.0,
        ..Default::default()
    };
    let filter = SpeckleFilter::with_params(params);

    let direct = filter.apply_filter(&image, SpeckleFilterType::GammaMAP)?;
    let tiled = filter.apply_filter_tiled(&image, SpeckleFilterType::GammaMAP, Some(16))?;

    let mut max_diff = 0.0f32;
    for (a, b) in direct.iter().zip(tiled.iter()) {
        max_diff = max_diff.max((a - b).abs());
    }

    assert!(
        max_diff < 1e-3,
        "Tiled Gamma-MAP output diverged (max diff: {max_diff})"
    );
    Ok(())
}

#[test]
fn filters_do_not_emit_negative_or_nan() -> SarResult<()> {
    let image = make_positive_image(32, 32, 0xC0FF_EE00);
    let params = SpeckleFilterParams {
        window_size: 3,
        num_looks: 2.0,
        damping_factor: 0.8,
        ..Default::default()
    };
    let filter = SpeckleFilter::with_params(params);

    for filter_type in [
        SpeckleFilterType::Mean,
        SpeckleFilterType::Median,
        SpeckleFilterType::Frost,
        SpeckleFilterType::RefinedLee,
    ] {
        let filtered = filter.apply_filter(&image, filter_type)?;
        assert!(filtered.iter().all(|v| v.is_finite() && *v >= 0.0));
    }

    Ok(())
}
