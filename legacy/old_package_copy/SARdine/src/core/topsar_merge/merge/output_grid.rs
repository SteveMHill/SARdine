#![allow(dead_code)]
//! Output grid computation for TOPSAR merge.
//!
//! Determines the final merged image dimensions and allocates output arrays.

use ndarray::Array2;
use num_complex::Complex32;

use crate::types::{SarResult, SubSwath};

use super::types::{AzimuthTimingModel, OutputGrid};
use super::assert_subswath_geometry;

/// Default satellite velocity (m/s) when not available from metadata
const DEFAULT_SATELLITE_VELOCITY_MPS: f64 = 7_500.0;

/// Compute the output grid dimensions from subswath extents.
pub fn compute_output_grid(
    subswaths: &[SubSwath],
    azimuth_timing: AzimuthTimingModel,
    azimuth_index_origin: usize,
) -> SarResult<OutputGrid> {
    if subswaths.is_empty() {
        return Err(crate::types::SarError::Processing(
            "No subswaths provided for output grid computation".into(),
        ));
    }

    #[cfg(debug_assertions)]
    for sw in subswaths {
        assert_subswath_geometry(sw);
    }

    // Compute range extent across all swaths
    let mut min_range_start = usize::MAX;
    let mut max_range_end = 0usize;
    let mut min_azimuth_start = usize::MAX;
    let mut max_azimuth_end = 0usize;

    for sw in subswaths {
        // Use global coordinates for range
        let range_start = sw.first_sample_global;
        let range_end = sw.last_sample_global;
        // Use global coordinates for azimuth
        let az_start = sw.first_line_global;
        let az_end = sw.last_line_global;

        min_range_start = min_range_start.min(range_start);
        max_range_end = max_range_end.max(range_end);
        min_azimuth_start = min_azimuth_start.min(az_start);
        max_azimuth_end = max_azimuth_end.max(az_end);
    }

    let range_samples = max_range_end.saturating_sub(min_range_start);
    let azimuth_samples = max_azimuth_end.saturating_sub(azimuth_index_origin);

    // Use first subswath for pixel spacing (assumes consistent spacing across swaths)
    let first = &subswaths[0];
    let range_pixel_spacing = first.range_pixel_spacing;

    // Compute azimuth pixel spacing from timing and velocity
    // SCIENTIFIC FIX: Require PRF from metadata - hardcoded fallback produces incorrect timing
    let prf = first.prf_hz.ok_or_else(|| {
        crate::types::SarError::MissingMetadata(
            "PRF (prf_hz) is required for azimuth timing computation. \
            Hardcoded fallback values produce incorrect geocoding results."
                .into(),
        )
    })?;
    let dt = first.azimuth_time_interval.unwrap_or_else(|| {
        log::debug!(
            "Using dt = 1/PRF = {} seconds (azimuth_time_interval not set)",
            1.0 / prf
        );
        1.0 / prf
    });
    let velocity = DEFAULT_SATELLITE_VELOCITY_MPS;
    let azimuth_pixel_spacing = dt * velocity;

    let range_time_start = first.slant_range_time;
    let azimuth_time_start = azimuth_timing.reference_azimuth_time;

    Ok(OutputGrid {
        range_samples,
        azimuth_samples,
        range_pixel_spacing,
        azimuth_pixel_spacing,
        range_time_start,
        azimuth_time_start,
        azimuth_timing,
    })
}

/// Allocate output intensity array initialized to zero.
pub fn allocate_intensity_output(output_grid: &OutputGrid) -> Array2<f32> {
    Array2::zeros((output_grid.azimuth_samples, output_grid.range_samples))
}

/// Allocate output complex array initialized to zero.
pub fn allocate_complex_output(output_grid: &OutputGrid) -> Array2<Complex32> {
    Array2::zeros((output_grid.azimuth_samples, output_grid.range_samples))
}

/// Allocate weight sum accumulator.
pub fn allocate_weight_sum(output_grid: &OutputGrid) -> Array2<f32> {
    Array2::zeros((output_grid.azimuth_samples, output_grid.range_samples))
}

/// Allocate hit count array for coverage tracking.
pub fn allocate_hitcount(output_grid: &OutputGrid) -> Array2<f32> {
    Array2::zeros((output_grid.azimuth_samples, output_grid.range_samples))
}

/// Allocate uncovered mask (0 = covered, 1 = uncovered).
pub fn allocate_uncovered_mask(output_grid: &OutputGrid) -> Array2<u8> {
    Array2::ones((output_grid.azimuth_samples, output_grid.range_samples))
}

/// Container for all output arrays needed during merge.
pub struct MergeOutputBuffers {
    pub intensity: Array2<f32>,
    pub complex: Option<Array2<Complex32>>,
    pub weight_sum: Array2<f32>,
    pub hitcount: Array2<f32>,
    pub uncovered_mask: Array2<u8>,
}

impl MergeOutputBuffers {
    /// Allocate all output buffers for the given output grid.
    pub fn allocate(output_grid: &OutputGrid, include_complex: bool) -> Self {
        Self {
            intensity: allocate_intensity_output(output_grid),
            complex: if include_complex {
                Some(allocate_complex_output(output_grid))
            } else {
                None
            },
            weight_sum: allocate_weight_sum(output_grid),
            hitcount: allocate_hitcount(output_grid),
            uncovered_mask: allocate_uncovered_mask(output_grid),
        }
    }

    /// Get output dimensions (rows, cols).
    pub fn shape(&self) -> (usize, usize) {
        self.intensity.dim()
    }
}

/// Validate that output grid dimensions are sensible.
pub fn validate_output_grid(output_grid: &OutputGrid) -> SarResult<()> {
    if output_grid.azimuth_samples == 0 {
        return Err(crate::types::SarError::Processing(
            "Output grid has zero azimuth samples".into(),
        ));
    }
    if output_grid.range_samples == 0 {
        return Err(crate::types::SarError::Processing(
            "Output grid has zero range samples".into(),
        ));
    }
    if output_grid.range_pixel_spacing <= 0.0 {
        return Err(crate::types::SarError::Processing(
            "Invalid range pixel spacing".into(),
        ));
    }
    if output_grid.azimuth_pixel_spacing <= 0.0 {
        return Err(crate::types::SarError::Processing(
            "Invalid azimuth pixel spacing".into(),
        ));
    }

    // Sanity check on dimensions
    const MAX_DIMENSION: usize = 100_000;
    if output_grid.azimuth_samples > MAX_DIMENSION || output_grid.range_samples > MAX_DIMENSION {
        log::warn!(
            "⚠️  Output grid dimensions very large: {}x{} - ensure this is expected",
            output_grid.azimuth_samples,
            output_grid.range_samples
        );
    }

    Ok(())
}

/// Compute per-subswath column offset in the output grid.
///
/// Returns a map from subswath ID to the column offset where that swath starts.
pub fn compute_swath_offsets(subswaths: &[SubSwath]) -> std::collections::HashMap<String, usize> {
    #[cfg(debug_assertions)]
    for sw in subswaths {
        assert_subswath_geometry(sw);
    }

    subswaths
        .iter()
        .map(|sw| (sw.id.clone(), sw.first_sample_global))
        .collect()
}

/// Compute per-subswath azimuth row offset in the output grid.
pub fn compute_swath_azimuth_offsets(
    subswaths: &[SubSwath],
    azimuth_index_origin: usize,
) -> std::collections::HashMap<String, usize> {
    #[cfg(debug_assertions)]
    for sw in subswaths {
        assert_subswath_geometry(sw);
    }

    subswaths
        .iter()
        .map(|sw| {
            let offset = sw.first_line_global.saturating_sub(azimuth_index_origin);
            (sw.id.clone(), offset)
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::types::SubSwath;

    fn mock_subswath(id: &str, range_start: usize, range_samples: usize) -> SubSwath {
        SubSwath {
            id: id.to_string(),
            burst_count: 8,
            lines_per_burst: 125,
            range_samples,
            azimuth_samples: 1000,
            first_line_global: 0,
            last_line_global: 1000,
            first_sample_global: range_start,
            last_sample_global: range_start + range_samples,
            full_range_samples: range_samples,
            valid_first_line: Some(0),
            valid_last_line: Some(1000),
            valid_first_sample: Some(range_start),
            valid_last_sample: Some(range_start + range_samples),
            range_pixel_spacing: 2.3,
            azimuth_pixel_spacing: 14.0,
            slant_range_time: 0.005,
            burst_duration: 2.75,
            near_range_m: 800_000.0,
            prf_hz: Some(1679.902),
            dc_polynomial: None,
            azimuth_time_interval: Some(0.002),
            dc_polynomial_t0: None,
            fm_rate_estimates: None,
        }
    }

    #[test]
    fn test_output_grid_dimensions() {
        let subswaths = vec![
            mock_subswath("IW1", 0, 10000),
            mock_subswath("IW2", 8000, 10000),
            mock_subswath("IW3", 16000, 10000),
        ];

        let timing = AzimuthTimingModel {
            prf: 1679.902,
            azimuth_time_interval: 0.002,
            burst_timing: vec![],
            reference_azimuth_time: 0.0,
        };

        let grid = compute_output_grid(&subswaths, timing, 0).unwrap();

        // With the mock setup, range should span from 0 to 26000
        assert!(grid.range_samples > 0);
        assert!(grid.azimuth_samples > 0);
    }

    #[test]
    fn test_swath_offsets() {
        let subswaths = vec![
            mock_subswath("IW1", 0, 10000),
            mock_subswath("IW2", 8000, 10000),
            mock_subswath("IW3", 16000, 10000),
        ];

        let offsets = compute_swath_offsets(&subswaths);

        assert_eq!(*offsets.get("IW1").unwrap(), 0);
        assert_eq!(*offsets.get("IW2").unwrap(), 8000);
        assert_eq!(*offsets.get("IW3").unwrap(), 16000);
    }
}
