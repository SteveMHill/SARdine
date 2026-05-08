use ndarray::Array2;

/// Simple mask helper: valid if finite and non-negative.
#[inline]
pub fn is_valid(v: f32) -> bool {
    v.is_finite() && v >= 0.0
}

/// Apply a boolean mask to zero-out invalid pixels (returns new array).
pub fn apply_mask(image: &Array2<f32>, mask: &Array2<bool>) -> Array2<f32> {
    let mut out = image.clone();
    for ((y, x), v) in out.indexed_iter_mut() {
        if !mask[[y, x]] {
            *v = 0.0;
        }
    }
    out
}
