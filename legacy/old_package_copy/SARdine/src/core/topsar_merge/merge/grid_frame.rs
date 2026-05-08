#![allow(dead_code, unused_variables)]
//! GridFrame and IndexConvention: enforce a single coordinate system across merge.
//!
//! # Coordinate Frame Conventions
//!
//! All indices in the TOPSAR merge pipeline use **exclusive-end** convention:
//! - `first_line..last_line` means lines `[first_line, last_line)`
//! - `first_sample..last_sample` means samples `[first_sample, last_sample)`
//!
//! This matches Rust's standard range semantics and prevents off-by-one errors.
//!
//! # Frame Types
//!
//! - **Global**: Absolute indices after origin normalization (starts at 0)
//! - **Local**: Indices relative to a specific subswath's first_line/first_sample
//! - **Overlap**: Indices within an overlap region's local coordinate system

use crate::types::{SarError, SarResult, SubSwath};

use super::assert_subswath_geometry;

/// Index convention marker - all extents use exclusive end
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum IndexConvention {
    /// Start is inclusive, end is exclusive: [start, end)
    #[default]
    ExclusiveEnd,
}

/// Grid frame definition with explicit origin and convention
#[derive(Debug, Clone)]
pub struct GridFrame {
    /// Azimuth origin that was subtracted from raw annotation indices
    pub az_origin: usize,
    /// Range origin that was subtracted from raw annotation indices
    pub rg_origin: usize,
    /// Index convention (always ExclusiveEnd)
    pub convention: IndexConvention,
}

impl Default for GridFrame {
    fn default() -> Self {
        Self {
            az_origin: 0,
            rg_origin: 0,
            convention: IndexConvention::ExclusiveEnd,
        }
    }
}

impl GridFrame {
    /// Create a new grid frame with specified origins
    pub fn new(az_origin: usize, rg_origin: usize) -> Self {
        Self {
            az_origin,
            rg_origin,
            convention: IndexConvention::ExclusiveEnd,
        }
    }

    /// Convert global (normalized) coordinates to local swath coordinates
    ///
    /// Returns None if the global coordinate falls outside the swath's extent
    #[inline]
    pub fn global_to_local_swath(
        &self,
        swath: &SubSwath,
        global_row: usize,
        global_col: usize,
    ) -> Option<(usize, usize)> {
        #[cfg(debug_assertions)]
        assert_subswath_geometry(swath);

        // Check bounds in global coordinates
        if global_row < swath.first_line_global || global_row >= swath.last_line_global {
            return None;
        }
        if global_col < swath.first_sample_global || global_col >= swath.last_sample_global {
            return None;
        }

        // Convert to local coordinates
        let local_row = global_row.checked_sub(swath.first_line_global)?;
        let local_col = global_col.checked_sub(swath.first_sample_global)?;

        // Verify local coordinates are within array bounds
        if local_row >= swath.azimuth_samples || local_col >= swath.range_samples {
            return None;
        }

        Some((local_row, local_col))
    }

    /// Convert local swath coordinates to global (normalized) coordinates
    #[inline]
    pub fn swath_local_to_global(
        &self,
        swath: &SubSwath,
        local_row: usize,
        local_col: usize,
    ) -> Option<(usize, usize)> {
        #[cfg(debug_assertions)]
        assert_subswath_geometry(swath);

        // Verify local coordinates are within swath bounds
        if local_row >= swath.azimuth_samples || local_col >= swath.range_samples {
            return None;
        }

        let global_row = swath.first_line_global.checked_add(local_row)?;
        let global_col = swath.first_sample_global.checked_add(local_col)?;

        Some((global_row, global_col))
    }

    /// Check if a global coordinate is within a swath's valid region
    #[inline]
    pub fn is_in_valid_region(
        &self,
        swath: &SubSwath,
        global_row: usize,
        global_col: usize,
    ) -> bool {
        #[cfg(debug_assertions)]
        assert_subswath_geometry(swath);

        let valid_first_line = swath.valid_first_line.unwrap_or(swath.first_line_global);
        let valid_last_line = swath.valid_last_line.unwrap_or(swath.last_line_global);
        let valid_first_sample = swath
            .valid_first_sample
            .unwrap_or(swath.first_sample_global);
        let valid_last_sample = swath.valid_last_sample.unwrap_or(swath.last_sample_global);

        global_row >= valid_first_line
            && global_row < valid_last_line
            && global_col >= valid_first_sample
            && global_col < valid_last_sample
    }
}

/// Extent in a single dimension using exclusive-end convention
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Extent {
    /// Start index (inclusive)
    pub start: usize,
    /// End index (exclusive)
    pub end: usize,
}

impl Extent {
    /// Create a new extent
    #[inline]
    pub fn new(start: usize, end: usize) -> Self {
        debug_assert!(end >= start, "Extent end must be >= start");
        Self { start, end }
    }

    /// Length of this extent
    #[inline]
    pub fn len(&self) -> usize {
        self.end.saturating_sub(self.start)
    }

    /// Check if extent is empty
    #[inline]
    pub fn is_empty(&self) -> bool {
        self.end <= self.start
    }

    /// Compute intersection with another extent
    #[inline]
    pub fn intersect(&self, other: &Extent) -> Option<Extent> {
        let start = self.start.max(other.start);
        let end = self.end.min(other.end);
        if end > start {
            Some(Extent { start, end })
        } else {
            None
        }
    }

    /// Check if this extent contains a specific index
    #[inline]
    pub fn contains(&self, idx: usize) -> bool {
        idx >= self.start && idx < self.end
    }

    /// Translate extent by an offset (for coordinate frame conversion)
    #[inline]
    pub fn translate(&self, offset: isize) -> SarResult<Extent> {
        let new_start = if offset >= 0 {
            self.start.checked_add(offset as usize)
        } else {
            self.start.checked_sub((-offset) as usize)
        };
        let new_end = if offset >= 0 {
            self.end.checked_add(offset as usize)
        } else {
            self.end.checked_sub((-offset) as usize)
        };

        match (new_start, new_end) {
            (Some(s), Some(e)) if e >= s => Ok(Extent::new(s, e)),
            _ => Err(SarError::Processing(format!(
                "Extent translation by {} would produce invalid range",
                offset
            ))),
        }
    }
}

/// 2D extent for azimuth × range regions
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Extent2D {
    pub azimuth: Extent,
    pub range: Extent,
}

impl Extent2D {
    /// Create a new 2D extent
    pub fn new(az_start: usize, az_end: usize, rg_start: usize, rg_end: usize) -> Self {
        Self {
            azimuth: Extent::new(az_start, az_end),
            range: Extent::new(rg_start, rg_end),
        }
    }

    /// Check if both dimensions are non-empty
    pub fn is_valid(&self) -> bool {
        !self.azimuth.is_empty() && !self.range.is_empty()
    }

    /// Get dimensions (height, width)
    pub fn dims(&self) -> (usize, usize) {
        (self.azimuth.len(), self.range.len())
    }

    /// Compute 2D intersection
    pub fn intersect(&self, other: &Extent2D) -> Option<Extent2D> {
        let az = self.azimuth.intersect(&other.azimuth)?;
        let rg = self.range.intersect(&other.range)?;
        Some(Extent2D {
            azimuth: az,
            range: rg,
        })
    }
}

/// Normalize subswath metadata to use exclusive-end indices
///
/// Many annotation sources store last_* as inclusive; this function
/// converts to the exclusive-end convention used throughout merge.
///
/// NOTE: azimuth_samples is per-burst, so total azimuth = burst_count * azimuth_samples
pub fn normalize_to_exclusive_end(sw: &mut SubSwath) {
    // For azimuth: use total lines (burst_count * azimuth_samples)
    let total_azimuth_lines = sw.burst_count * sw.azimuth_samples;
    let inclusive_last_line = sw
        .first_line_global
        .saturating_add(total_azimuth_lines.saturating_sub(1));
    if sw.last_line_global == inclusive_last_line && total_azimuth_lines > 0 {
        sw.last_line_global = sw.last_line_global.saturating_add(1);
    }

    // For range: use range_samples directly
    let inclusive_last_sample = sw
        .first_sample_global
        .saturating_add(sw.range_samples.saturating_sub(1));
    if sw.last_sample_global == inclusive_last_sample && sw.range_samples > 0 {
        sw.last_sample_global = sw.last_sample_global.saturating_add(1);
    }

    // Also normalize valid_* fields
    if let Some(v_last) = sw.valid_last_line {
        if v_last == inclusive_last_line && total_azimuth_lines > 0 {
            sw.valid_last_line = Some(v_last.saturating_add(1));
        }
    }

    if let Some(v_last) = sw.valid_last_sample {
        if v_last == inclusive_last_sample && sw.range_samples > 0 {
            sw.valid_last_sample = Some(v_last.saturating_add(1));
        }
    }

    #[cfg(debug_assertions)]
    {
        if sw.burst_count > 0 && sw.lines_per_burst > 0 && sw.range_samples > 0 {
            assert_subswath_geometry(sw);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn test_swath() -> SubSwath {
        SubSwath {
            id: "IW1".to_string(),
            burst_count: 1,
            lines_per_burst: 50,
            range_samples: 100,
            azimuth_samples: 50,
            first_line_global: 10,
            last_line_global: 60, // exclusive
            first_sample_global: 20,
            last_sample_global: 120, // exclusive
            full_range_samples: 100,
            valid_first_line: Some(10),
            valid_last_line: Some(60),
            valid_first_sample: Some(20),
            valid_last_sample: Some(120),
            range_pixel_spacing: 2.3,
            azimuth_pixel_spacing: 14.0,
            slant_range_time: 0.005,
            burst_duration: 2.7,
            near_range_m: 800000.0,
            prf_hz: Some(1680.0),
            dc_polynomial: Some(vec![0.0, 0.0, 0.0]),
            azimuth_time_interval: Some(0.000595),
            dc_polynomial_t0: Some(0.0),
            fm_rate_estimates: None,
        }
    }

    #[test]
    fn global_to_local_in_bounds() {
        let frame = GridFrame::default();
        let sw = test_swath();

        // First pixel
        let result = frame.global_to_local_swath(&sw, 10, 20);
        assert_eq!(result, Some((0, 0)));

        // Last pixel (exclusive, so -1)
        let result = frame.global_to_local_swath(&sw, 59, 119);
        assert_eq!(result, Some((49, 99)));
    }

    #[test]
    fn global_to_local_out_of_bounds() {
        let frame = GridFrame::default();
        let sw = test_swath();

        // Before swath
        assert!(frame.global_to_local_swath(&sw, 9, 20).is_none());
        assert!(frame.global_to_local_swath(&sw, 10, 19).is_none());

        // After swath (at exclusive end)
        assert!(frame.global_to_local_swath(&sw, 60, 20).is_none());
        assert!(frame.global_to_local_swath(&sw, 10, 120).is_none());
    }

    #[test]
    fn local_to_global() {
        let frame = GridFrame::default();
        let sw = test_swath();

        let result = frame.swath_local_to_global(&sw, 0, 0);
        assert_eq!(result, Some((10, 20)));

        let result = frame.swath_local_to_global(&sw, 49, 99);
        assert_eq!(result, Some((59, 119)));

        // Out of bounds
        assert!(frame.swath_local_to_global(&sw, 50, 0).is_none());
    }

    #[test]
    fn extent_intersection() {
        let a = Extent::new(10, 30);
        let b = Extent::new(20, 40);

        let intersect = a.intersect(&b);
        assert_eq!(intersect, Some(Extent::new(20, 30)));

        // No overlap
        let c = Extent::new(40, 50);
        assert!(a.intersect(&c).is_none());
    }

    #[test]
    fn extent_translate() {
        let e = Extent::new(10, 20);

        let shifted_pos = e.translate(5).unwrap();
        assert_eq!(shifted_pos, Extent::new(15, 25));

        let shifted_neg = e.translate(-5).unwrap();
        assert_eq!(shifted_neg, Extent::new(5, 15));
    }

    #[test]
    fn normalize_inclusive_to_exclusive() {
        let mut sw = test_swath();
        // Simulate inclusive last values
        sw.last_line_global = sw.first_line_global + sw.azimuth_samples - 1;
        sw.last_sample_global = sw.first_sample_global + sw.range_samples - 1;

        normalize_to_exclusive_end(&mut sw);

        assert_eq!(sw.last_line_global, 60);
        assert_eq!(sw.last_sample_global, 120);
    }
}
