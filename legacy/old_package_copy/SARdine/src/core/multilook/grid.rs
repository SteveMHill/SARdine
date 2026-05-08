//! Multilook grid helpers: output dimensions and spacing calculations.
use crate::types::{SarError, SarResult};

/// Compute multilook output dimensions using ceil-style division when `include_partial`.
pub fn compute_output_dims(
    input_rows: usize,
    input_cols: usize,
    azimuth_looks: usize,
    range_looks: usize,
    include_partial: bool,
) -> SarResult<(usize, usize)> {
    if azimuth_looks == 0 || range_looks == 0 {
        return Err(SarError::InvalidParameter("looks must be >= 1".into()));
    }
    if input_rows == 0 || input_cols == 0 {
        return Err(SarError::Processing("empty input".into()));
    }

    let out_rows = if include_partial {
        (input_rows + azimuth_looks - 1) / azimuth_looks
    } else {
        input_rows / azimuth_looks
    };
    let out_cols = if include_partial {
        (input_cols + range_looks - 1) / range_looks
    } else {
        input_cols / range_looks
    };

    if out_rows == 0 || out_cols == 0 {
        return Err(SarError::Processing(
            "Multilook parameters too large for input image".into(),
        ));
    }

    Ok((out_rows, out_cols))
}

/// Compute updated pixel spacing after multilooking.
pub fn compute_output_spacing(
    range_spacing: f64,
    azimuth_spacing: f64,
    range_looks: usize,
    azimuth_looks: usize,
) -> (f64, f64) {
    let new_range = range_spacing * range_looks as f64;
    let new_azimuth = azimuth_spacing * azimuth_looks as f64;
    (new_range, new_azimuth)
}
