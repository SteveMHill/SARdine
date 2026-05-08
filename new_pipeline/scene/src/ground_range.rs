//! Ground range projection and multilooking for merged σ⁰ images.
//!
//! Converts a slant-range merged image to a ground-range image with
//! configurable output pixel spacing (default 10 m, matching ESA GRD High-Res).
//!
//! # Algorithm
//!
//! 1. Interpolate incidence angle at each merged range column from the
//!    sparse geolocation grid (typically 10 azimuth × 21 range per subswath).
//! 2. Compute cumulative ground range distance: each slant-range pixel maps to
//!    `slant_range_spacing / sin(θ)` in ground range.
//! 3. Resample onto a regular ground-range grid via linear interpolation.
//! 4. Average (multilook) in azimuth to achieve approximately square pixels.
//!
//! # SAR-specific notes
//!
//! - Incidence angle varies ~30° (IW1 near range) to ~46° (IW3 far range).
//! - Ground range spacing = slant range spacing / sin(incidence angle).
//! - At 30° incidence, ground pixel is 2× slant pixel; at 46°, ~1.4×.
//! - Multilooking in azimuth is performed AFTER ground range projection.
//! - This is NOT geocoding (no map projection). The output is still in
//!   radar geometry (azimuth × ground range).

use crate::merge_subswaths::MergedSigma0;
use crate::types::{GeolocationGridPoint, SubSwathId};

// ═══════════════════════════════════════════════════════════════════════
// Error type
// ═══════════════════════════════════════════════════════════════════════

#[derive(Debug, thiserror::Error)]
pub enum GroundRangeError {
    #[error("no geolocation grid points for {0}")]
    NoGridPoints(SubSwathId),

    #[error("merged image is empty ({0} lines × {1} samples)")]
    EmptyImage(usize, usize),

    #[error("target pixel spacing must be positive, got {0}")]
    InvalidSpacing(f64),

    #[error("input azimuth pixel spacing must be positive and finite, got {0}")]
    InvalidAzimuthSpacing(f64),

    #[error("incidence angle interpolation failed at column {col}: {detail}")]
    InterpolationFailed { col: usize, detail: String },
}

// ═══════════════════════════════════════════════════════════════════════
// Output type
// ═══════════════════════════════════════════════════════════════════════

/// Ground-range projected and multilooked σ⁰ image.
#[derive(Debug)]
pub struct GroundRangeImage {
    /// Flat row-major buffer. Length = `lines × samples`.
    pub data: Vec<f32>,

    /// Number of azimuth lines in the output.
    pub lines: usize,

    /// Number of ground-range samples in the output.
    pub samples: usize,

    /// Ground-range pixel spacing in metres.
    pub range_pixel_spacing_m: f64,

    /// Azimuth pixel spacing in metres.
    pub azimuth_pixel_spacing_m: f64,

    /// Number of slant-range pixels averaged per ground-range pixel (range multilook factor).
    pub range_looks: usize,

    /// Number of azimuth lines averaged per output line (azimuth multilook factor).
    pub azimuth_looks: usize,
}

// ═══════════════════════════════════════════════════════════════════════
// Incidence angle interpolation
// ═══════════════════════════════════════════════════════════════════════

/// A 1D incidence-angle profile evaluated at every column of the merged image.
///
/// The profile is computed at a representative azimuth line (typically mid-scene)
/// and is assumed constant across azimuth lines. This is valid because the
/// incidence angle is a geometric property that depends on range, not azimuth
/// (to first order).
pub(crate) struct IncidenceProfile {
    /// Incidence angle in radians for each column 0..num_samples.
    pub angles_rad: Vec<f64>,
}

/// Build a per-column incidence angle profile for the merged image.
///
/// Uses bilinear interpolation on the sparse geolocation grids from the
/// three IW subswaths. The grids are stitched by slant range time, matching
/// the merge offset logic in `merge_subswaths`.
///
/// # Arguments
///
/// * `grids` - Per-subswath geolocation grids, sorted by SubSwathId.
/// * `merged` - The merged σ⁰ image (for dimensions and range timing).
pub(crate) fn build_incidence_profile(
    grids: &[(SubSwathId, Vec<GeolocationGridPoint>)],
    merged: &MergedSigma0,
) -> Result<IncidenceProfile, GroundRangeError> {
    // Strategy: for each merged column, determine which subswath it belongs to
    // based on slant range time, then interpolate incidence angle from that
    // subswath's geolocation grid.
    //
    // We use a mid-scene azimuth line for interpolation. The incidence angle
    // variation in azimuth is negligible for a single frame (<0.01°).

    let num_samples = merged.samples;
    let mut angles_rad = Vec::with_capacity(num_samples);

    // Build per-subswath range-only interpolators at mid-azimuth
    let mut swath_interps: Vec<RangeInterpolator> = Vec::new();

    for (id, grid_points) in grids {
        if grid_points.is_empty() {
            return Err(GroundRangeError::NoGridPoints(*id));
        }
        let interp = RangeInterpolator::from_grid(grid_points)?;
        swath_interps.push(interp);
    }

    // Merge the subswath interpolators into one covering the full range
    // Sort by slant range time
    swath_interps.sort_by(|a, b| {
        a.min_slant_time
            .partial_cmp(&b.min_slant_time)
            .unwrap_or(std::cmp::Ordering::Equal) // SAFETY-OK: NaN-safe sort comparator (NaN treated as equal); slant times come from validated geolocation grid
    });

    for col in 0..num_samples {
        let slant_time =
            merged.near_slant_range_time_s + col as f64 * merged.range_pixel_spacing_s;

        // Find the best interpolator for this slant range time
        let angle_deg = interpolate_at_slant_time(&swath_interps, slant_time)
            .ok_or_else(|| GroundRangeError::InterpolationFailed {
                col,
                detail: format!(
                    "slant time {:.6e} s outside all grid ranges",
                    slant_time
                ),
            })?;

        angles_rad.push(angle_deg.to_radians());
    }

    Ok(IncidenceProfile { angles_rad })
}

/// Interpolator for one subswath: maps slant range time → incidence angle
/// using the mid-azimuth row of the geolocation grid.
struct RangeInterpolator {
    /// Slant range times of the grid points (sorted, ascending).
    slant_times: Vec<f64>,
    /// Incidence angles in degrees at corresponding slant range times.
    angles_deg: Vec<f64>,
    /// Min/max slant range time for this subswath.
    min_slant_time: f64,
    max_slant_time: f64,
}

impl RangeInterpolator {
    fn from_grid(grid_points: &[GeolocationGridPoint]) -> Result<Self, GroundRangeError> {
        // Find unique line indices (azimuth rows)
        let mut lines: Vec<u32> = grid_points.iter().map(|p| p.line).collect();
        lines.sort();
        lines.dedup();

        // Pick the middle azimuth row
        let mid_line = lines[lines.len() / 2];

        // Extract the range profile at mid-azimuth
        let mut range_points: Vec<(f64, f64)> = grid_points
            .iter()
            .filter(|p| p.line == mid_line)
            .map(|p| (p.slant_range_time_s, p.incidence_angle_deg))
            .collect();

        range_points.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal)); // SAFETY-OK: NaN-safe sort comparator (NaN treated as equal); slant times come from validated geolocation grid

        let slant_times: Vec<f64> = range_points.iter().map(|&(t, _)| t).collect();
        let angles_deg: Vec<f64> = range_points.iter().map(|&(_, a)| a).collect();

        let min_slant_time = slant_times[0];
        let max_slant_time = *slant_times.last().unwrap();

        Ok(RangeInterpolator {
            slant_times,
            angles_deg,
            min_slant_time,
            max_slant_time,
        })
    }

    /// Linearly interpolate incidence angle at a given slant range time.
    /// Returns `None` if outside the interpolator's range.
    fn interpolate(&self, slant_time: f64) -> Option<f64> {
        if slant_time < self.min_slant_time || slant_time > self.max_slant_time {
            return None;
        }

        // Binary search for the bracketing interval
        let idx = match self
            .slant_times
            .binary_search_by(|t| t.partial_cmp(&slant_time).unwrap())
        {
            Ok(i) => return Some(self.angles_deg[i]), // exact match
            Err(i) => i,
        };

        if idx == 0 {
            return Some(self.angles_deg[0]);
        }
        if idx >= self.slant_times.len() {
            return Some(*self.angles_deg.last().unwrap());
        }

        // Linear interpolation
        let t0 = self.slant_times[idx - 1];
        let t1 = self.slant_times[idx];
        let a0 = self.angles_deg[idx - 1];
        let a1 = self.angles_deg[idx];
        let frac = (slant_time - t0) / (t1 - t0);
        Some(a0 + frac * (a1 - a0))
    }
}

/// Find the best interpolator for a given slant range time and interpolate.
///
/// When multiple subswaths overlap (they do at the seams), use the one whose
/// grid center is closest. When none overlap, extrapolate from the nearest.
fn interpolate_at_slant_time(interps: &[RangeInterpolator], slant_time: f64) -> Option<f64> {
    // Try direct interpolation first
    let mut best: Option<(f64, f64)> = None; // (distance to center, angle)
    for interp in interps {
        if let Some(angle) = interp.interpolate(slant_time) {
            let center = (interp.min_slant_time + interp.max_slant_time) / 2.0;
            let dist = (slant_time - center).abs();
            if best.is_none() || dist < best.unwrap().0 {
                best = Some((dist, angle));
            }
        }
    }
    if let Some((_, angle)) = best {
        return Some(angle);
    }

    // Extrapolate from nearest edge (only for small gap regions at merge seams)
    let mut nearest: Option<(f64, &RangeInterpolator)> = None;
    for interp in interps {
        let dist = if slant_time < interp.min_slant_time {
            interp.min_slant_time - slant_time
        } else {
            slant_time - interp.max_slant_time
        };
        if nearest.is_none() || dist < nearest.unwrap().0 {
            nearest = Some((dist, interp));
        }
    }

    // Only allow small extrapolation (< 1% of swath width)
    if let Some((dist, interp)) = nearest {
        let swath_width = interp.max_slant_time - interp.min_slant_time;
        if dist < swath_width * 0.01 {
            // Clamp to edge
            let clamped = slant_time.clamp(interp.min_slant_time, interp.max_slant_time);
            return interp.interpolate(clamped);
        }
    }

    None
}

// ═══════════════════════════════════════════════════════════════════════
// Ground range projection + multilooking
// ═══════════════════════════════════════════════════════════════════════

/// Project a merged slant-range σ⁰ image to ground range and multilook.
///
/// # Arguments
///
/// * `merged` - Merged slant-range σ⁰ image from [`crate::merge_subswaths::merge_subswaths`].
/// * `grids`  - Per-subswath geolocation grids from [`crate::parse::parse_geolocation_grids`].
/// * `target_spacing_m` - Desired output pixel spacing in metres (both range and azimuth).
/// * `azimuth_pixel_spacing_m` - Azimuth pixel spacing of the input slant-range image
///   in metres.  Must be taken from
///   [`crate::types::SubSwathMetadata::azimuth_pixel_spacing_m`] of a representative
///   subswath (all IW subswaths share this value to within ≪ 0.1 %).  Must be positive.
///
/// # Algorithm
///
/// 1. Build a per-column incidence angle profile from the geolocation grids.
/// 2. Compute cumulative ground range distance for each slant-range column.
/// 3. Determine the range multilook factor and azimuth multilook factor to
///    achieve approximately `target_spacing_m` in both dimensions.
/// 4. Average (multilook) the image to produce the output.
pub fn to_ground_range(
    merged: &MergedSigma0,
    grids: &[(SubSwathId, Vec<GeolocationGridPoint>)],
    target_spacing_m: f64,
    azimuth_pixel_spacing_m: f64,
) -> Result<GroundRangeImage, GroundRangeError> {
    if merged.lines == 0 || merged.samples == 0 {
        return Err(GroundRangeError::EmptyImage(merged.lines, merged.samples));
    }
    if target_spacing_m <= 0.0 {
        return Err(GroundRangeError::InvalidSpacing(target_spacing_m));
    }
    if !azimuth_pixel_spacing_m.is_finite() || azimuth_pixel_spacing_m <= 0.0 {
        return Err(GroundRangeError::InvalidAzimuthSpacing(
            azimuth_pixel_spacing_m,
        ));
    }

    // 1. Incidence angle profile
    let profile = build_incidence_profile(grids, merged)?;

    // 2. Compute ground range spacing per slant-range column
    //    ground_range_spacing = slant_range_spacing_m / sin(θ)
    let slant_spacing_m = merged.range_pixel_spacing_m;

    let ground_spacings: Vec<f64> = profile
        .angles_rad
        .iter()
        .map(|&theta| slant_spacing_m / theta.sin())
        .collect();

    // 3. Cumulative ground range distance from column 0
    let mut cum_ground_range = Vec::with_capacity(merged.samples);
    let mut acc = 0.0_f64;
    cum_ground_range.push(0.0);
    for i in 1..merged.samples {
        // Average ground spacing between adjacent columns
        acc += (ground_spacings[i - 1] + ground_spacings[i]) / 2.0;
        cum_ground_range.push(acc);
    }

    let total_ground_range = *cum_ground_range.last().unwrap();

    // 4. Define the regular output ground-range grid
    let out_range_samples = (total_ground_range / target_spacing_m).floor() as usize;
    if out_range_samples == 0 {
        return Err(GroundRangeError::InvalidSpacing(target_spacing_m));
    }

    // 5. Resample each azimuth line from slant range to regular ground range
    //    using linear interpolation in the cumulative ground range table.
    //
    //    For each output column j, the output ground range is j * target_spacing_m.
    //    Find the two slant-range columns that bracket that ground range distance
    //    and linearly interpolate the σ⁰ value.
    let mut projected = vec![f32::NAN; merged.lines * out_range_samples];

    // Pre-compute the mapping: for each output column, find the slant-range interval
    let mapping = build_ground_to_slant_mapping(&cum_ground_range, target_spacing_m, out_range_samples);

    for line in 0..merged.lines {
        let in_row_start = line * merged.samples;
        let out_row_start = line * out_range_samples;

        for (j, map) in mapping.iter().enumerate() {
            let val = match *map {
                ColumnMapping::Exact(idx) => merged.data[in_row_start + idx],
                ColumnMapping::Interpolate { left, right, frac } => {
                    let v0 = merged.data[in_row_start + left];
                    let v1 = merged.data[in_row_start + right];
                    if v0.is_nan() || v1.is_nan() {
                        f32::NAN
                    } else {
                        v0 + (v1 - v0) * frac as f32
                    }
                }
            };
            projected[out_row_start + j] = val;
        }
    }

    // 6. Determine azimuth multilook factor from the caller-provided input
    //    azimuth pixel spacing.  The merged slant-range image preserves the SLC
    //    line spacing, so the input azimuth spacing equals the subswath's
    //    `azimuth_pixel_spacing_m` from the product annotation.
    let azimuth_looks =
        (target_spacing_m / azimuth_pixel_spacing_m).round().max(1.0) as usize;

    // 7. Azimuth multilook: average `azimuth_looks` consecutive lines
    let out_lines = merged.lines / azimuth_looks;
    let mut data = vec![f32::NAN; out_lines * out_range_samples];

    for out_line in 0..out_lines {
        let az_start = out_line * azimuth_looks;
        let az_end = az_start + azimuth_looks;
        let out_row_start = out_line * out_range_samples;

        for col in 0..out_range_samples {
            let mut sum = 0.0_f64;
            let mut count = 0u32;
            for line in az_start..az_end {
                let val = projected[line * out_range_samples + col];
                if !val.is_nan() {
                    sum += val as f64;
                    count += 1;
                }
            }
            data[out_row_start + col] = if count > 0 {
                (sum / count as f64) as f32
            } else {
                f32::NAN
            };
        }
    }

    Ok(GroundRangeImage {
        data,
        lines: out_lines,
        samples: out_range_samples,
        range_pixel_spacing_m: target_spacing_m,
        azimuth_pixel_spacing_m: azimuth_pixel_spacing_m * azimuth_looks as f64,
        range_looks: 1, // ground range resampling, not pixel averaging
        azimuth_looks,
    })
}

// ═══════════════════════════════════════════════════════════════════════
// Column mapping helpers
// ═══════════════════════════════════════════════════════════════════════

#[derive(Debug, Clone, Copy)]
enum ColumnMapping {
    Exact(usize),
    Interpolate {
        left: usize,
        right: usize,
        frac: f64,
    },
}

/// Pre-compute the slant-range column mapping for each output ground-range column.
fn build_ground_to_slant_mapping(
    cum_ground_range: &[f64],
    target_spacing: f64,
    out_samples: usize,
) -> Vec<ColumnMapping> {
    let mut mapping = Vec::with_capacity(out_samples);
    let mut search_start = 0;

    for j in 0..out_samples {
        let target_gr = j as f64 * target_spacing;

        // Binary search in cum_ground_range for the bracketing interval
        // (use search_start to avoid redundant scanning; cum_ground_range is monotonic)
        let idx = match cum_ground_range[search_start..]
            .binary_search_by(|r| r.partial_cmp(&target_gr).unwrap())
        {
            Ok(i) => {
                mapping.push(ColumnMapping::Exact(search_start + i));
                search_start = search_start + i;
                continue;
            }
            Err(i) => search_start + i,
        };

        if idx == 0 {
            mapping.push(ColumnMapping::Exact(0));
        } else if idx >= cum_ground_range.len() {
            mapping.push(ColumnMapping::Exact(cum_ground_range.len() - 1));
        } else {
            let r0 = cum_ground_range[idx - 1];
            let r1 = cum_ground_range[idx];
            let frac = (target_gr - r0) / (r1 - r0);
            mapping.push(ColumnMapping::Interpolate {
                left: idx - 1,
                right: idx,
                frac,
            });
            search_start = idx - 1;
        }
    }

    mapping
}

// ═══════════════════════════════════════════════════════════════════════
// Tests
// ═══════════════════════════════════════════════════════════════════════

#[cfg(test)]
mod tests {
    use super::*;

    /// Create a simple synthetic incidence angle profile for testing.
    fn make_test_profile(num_cols: usize, min_deg: f64, max_deg: f64) -> IncidenceProfile {
        let angles_rad: Vec<f64> = (0..num_cols)
            .map(|i| {
                let frac = i as f64 / (num_cols - 1) as f64;
                (min_deg + frac * (max_deg - min_deg)).to_radians()
            })
            .collect();
        IncidenceProfile { angles_rad }
    }

    #[test]
    fn test_ground_spacing_increases_with_steeper_incidence() {
        // At 30° incidence: ground spacing = slant / sin(30°) = 2× slant
        // At 45° incidence: ground spacing = slant / sin(45°) ≈ 1.41× slant
        let slant_spacing = 2.33;
        let gs_30 = slant_spacing / (30.0_f64).to_radians().sin();
        let gs_45 = slant_spacing / (45.0_f64).to_radians().sin();

        assert!((gs_30 - 4.66).abs() < 0.01, "gs_30 = {}", gs_30);
        assert!((gs_45 - 3.295).abs() < 0.01, "gs_45 = {}", gs_45);
        assert!(gs_30 > gs_45, "ground spacing should decrease with increasing incidence");
    }

    #[test]
    fn test_column_mapping_monotonic() {
        // Synthetic ground range: linearly increasing
        let n = 1000;
        let cum: Vec<f64> = (0..n).map(|i| i as f64 * 3.5).collect(); // 3.5m per column
        let target = 10.0; // 10m output spacing
        let out_samples = 349; // 1000 * 3.5 / 10 ≈ 350

        let mapping = build_ground_to_slant_mapping(&cum, target, out_samples);
        assert_eq!(mapping.len(), out_samples);

        // Check monotonicity: each output column should map to same or later input column
        let mut prev_idx = 0;
        for m in &mapping {
            let idx = match m {
                ColumnMapping::Exact(i) => *i,
                ColumnMapping::Interpolate { left, .. } => *left,
            };
            assert!(idx >= prev_idx, "mapping not monotonic");
            prev_idx = idx;
        }
    }

    #[test]
    fn test_column_mapping_covers_range() {
        let n = 500;
        let cum: Vec<f64> = (0..n).map(|i| i as f64 * 4.0).collect();
        let target = 10.0;
        let out_samples = (499.0_f64 * 4.0 / 10.0).floor() as usize;

        let mapping = build_ground_to_slant_mapping(&cum, target, out_samples);

        // First output column should map near column 0
        match mapping[0] {
            ColumnMapping::Exact(i) => assert_eq!(i, 0),
            ColumnMapping::Interpolate { left, .. } => assert_eq!(left, 0),
        }
    }

    #[test]
    fn test_range_interpolator_basic() {
        let points = vec![
            GeolocationGridPoint {
                azimuth_time_utc: chrono::Utc::now(),
                slant_range_time_s: 5.0e-3,
                line: 100,
                pixel: 0,
                latitude_deg: 50.0,
                longitude_deg: 8.0,
                height_m: 100.0,
                incidence_angle_deg: 30.0,
                elevation_angle_deg: 27.0,
            },
            GeolocationGridPoint {
                azimuth_time_utc: chrono::Utc::now(),
                slant_range_time_s: 6.0e-3,
                line: 100,
                pixel: 1000,
                latitude_deg: 50.0,
                longitude_deg: 9.0,
                height_m: 100.0,
                incidence_angle_deg: 36.0,
                elevation_angle_deg: 33.0,
            },
        ];

        let interp = RangeInterpolator::from_grid(&points).unwrap();

        // At midpoint: (5.0 + 6.0) / 2 = 5.5e-3 → expect 33°
        let angle = interp.interpolate(5.5e-3).unwrap();
        assert!((angle - 33.0).abs() < 0.01, "midpoint angle = {}", angle);

        // At start
        let angle = interp.interpolate(5.0e-3).unwrap();
        assert!((angle - 30.0).abs() < 0.01);

        // At end
        let angle = interp.interpolate(6.0e-3).unwrap();
        assert!((angle - 36.0).abs() < 0.01);

        // Outside range: should return None
        assert!(interp.interpolate(4.9e-3).is_none());
        assert!(interp.interpolate(6.1e-3).is_none());
    }

    #[test]
    fn test_to_ground_range_synthetic() {
        // Create a small synthetic merged image with uniform σ⁰ = 0.1
        let lines = 100;
        let samples = 200;
        let slant_spacing_m = 2.33;
        let near_time = 5.0e-3;
        let c = 299_792_458.0;
        let rps_s = slant_spacing_m * 2.0 / c;

        let merged = MergedSigma0 {
            data: vec![0.1_f32; lines * samples],
            nesz: vec![0.0f32; lines * samples],
            lines,
            samples,
            near_slant_range_time_s: near_time,
            range_pixel_spacing_s: rps_s,
            range_pixel_spacing_m: slant_spacing_m,
            cal_lut_extrapolation_gap_px: 0,
            noise_lut_extrapolation_gap_px: 0,
            azimuth_start_time: chrono::Utc::now(),
        };

        // Create a simple geolocation grid for one subswath
        // Incidence angle varies linearly from 30° to 36°
        let mut grid_points = Vec::new();
        for &line in &[0u32, 50, 99] {
            for i in 0..5 {
                let pixel = i * 50;
                let frac = pixel as f64 / 199.0;
                let inc = 30.0 + frac * 6.0;
                let srt = near_time + pixel as f64 * rps_s;
                grid_points.push(GeolocationGridPoint {
                    azimuth_time_utc: chrono::Utc::now(),
                    slant_range_time_s: srt,
                    line,
                    pixel,
                    latitude_deg: 50.0,
                    longitude_deg: 8.0,
                    height_m: 100.0,
                    incidence_angle_deg: inc,
                    elevation_angle_deg: inc - 3.0,
                });
            }
        }

        let grids = vec![(SubSwathId::IW1, grid_points)];
        let result = to_ground_range(&merged, &grids, 10.0, 14.08).unwrap();

        // Output should have fewer range samples (ground range > slant at <45°)
        assert!(result.samples < samples, "expected fewer ground range samples");
        assert!(result.samples > 0);
        assert!(result.lines > 0);
        assert_eq!(result.range_pixel_spacing_m, 10.0);

        // All values should be close to 0.1 (uniform input → uniform output)
        let valid: Vec<f32> = result.data.iter().copied().filter(|v| !v.is_nan()).collect();
        assert!(!valid.is_empty(), "should have valid pixels");
        let mean = valid.iter().map(|&v| v as f64).sum::<f64>() / valid.len() as f64;
        assert!(
            (mean - 0.1).abs() < 0.001,
            "mean should be ~0.1, got {}",
            mean
        );
    }

    #[test]
    fn test_to_ground_range_rejects_invalid_input() {
        let merged = MergedSigma0 {
            data: vec![],
            nesz: vec![],
            lines: 0,
            samples: 0,
            near_slant_range_time_s: 5e-3,
            range_pixel_spacing_s: 1e-8,
            range_pixel_spacing_m: 2.33,
            cal_lut_extrapolation_gap_px: 0,
            noise_lut_extrapolation_gap_px: 0,
            azimuth_start_time: chrono::Utc::now(),
        };
        let grids = vec![];
        assert!(to_ground_range(&merged, &grids, 10.0, 14.08).is_err());
    }

    #[test]
    fn test_to_ground_range_rejects_negative_spacing() {
        let merged = MergedSigma0 {
            data: vec![1.0; 100],
            nesz: vec![0.0f32; 100],
            lines: 10,
            samples: 10,
            near_slant_range_time_s: 5e-3,
            range_pixel_spacing_s: 1e-8,
            range_pixel_spacing_m: 2.33,
            cal_lut_extrapolation_gap_px: 0,
            noise_lut_extrapolation_gap_px: 0,
            azimuth_start_time: chrono::Utc::now(),
        };
        let grids = vec![];
        assert!(to_ground_range(&merged, &grids, -10.0, 14.08).is_err());
    }

    #[test]
    fn test_to_ground_range_rejects_invalid_azimuth_spacing() {
        let merged = MergedSigma0 {
            data: vec![1.0; 100],
            nesz: vec![0.0f32; 100],
            lines: 10,
            samples: 10,
            near_slant_range_time_s: 5e-3,
            range_pixel_spacing_s: 1e-8,
            range_pixel_spacing_m: 2.33,
            cal_lut_extrapolation_gap_px: 0,
            noise_lut_extrapolation_gap_px: 0,
            azimuth_start_time: chrono::Utc::now(),
        };
        let grids = vec![];
        assert!(to_ground_range(&merged, &grids, 10.0, 0.0).is_err());
        assert!(to_ground_range(&merged, &grids, 10.0, -1.0).is_err());
        assert!(to_ground_range(&merged, &grids, 10.0, f64::NAN).is_err());
    }
}
