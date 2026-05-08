//! Multilook weighting helpers (box/feather placeholders).

/// Uniform weight (boxcar) for multilook averaging.
#[inline]
pub fn box_weight() -> f32 {
    1.0
}

/// Placeholder for future windowed weights (e.g., Hamming/Hann).
#[inline]
pub fn sample_weight(_row: usize, _col: usize) -> f32 {
    1.0
}
