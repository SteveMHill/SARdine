//! Multi-TIFF SLC reader for slice-assembled Sentinel-1 scenes.
//!
//! A single subswath in a slice-assembled scene spans N consecutive TIFF files
//! (one per `.SAFE` slice).  [`MultiSlcReader`] wraps N [`SlcReader`]s behind
//! a single [`BurstReader`] interface, translating logical line numbers (as
//! stored in assembled [`BurstEntry::first_line`][crate::types::BurstEntry])
//! into per-TIFF offsets.
//!
//! # Line numbering contract
//!
//! In an assembled scene produced by [`crate::slice_assembly::assemble_slices`]:
//!
//! ```text
//! burst.first_line = slice_line_offsets[slice_index] + burst_local_first_line
//! ```
//!
//! where `slice_line_offsets[k]` = Σ `azimuth_samples` for slices 0..k.
//!
//! [`MultiSlcReader`] reconstructs those offsets from the physical TIFF heights
//! of each opened reader (`offset[k] = Σ height(readers[0..k])`).  The two
//! calculations are equivalent because `azimuth_samples = n_bursts * lines_per_burst`
//! and the TIFF height of a subswath is exactly `n_bursts * lines_per_burst`.
//!
//! # Safety contract
//!
//! Every burst read must be entirely within one TIFF file.  A read that would
//! cross a slice boundary is rejected with [`SlcReadError::CrossSliceBoundary`].
//! In a correctly assembled scene this never occurs, so any such error indicates
//! corrupted assembled burst geometry.

use crate::slc_reader::{BurstReader, SlcReadError, SlcReader};

/// Multi-TIFF reader for one (subswath, polarization) channel of a
/// slice-assembled Sentinel-1 scene.
///
/// Construct with [`MultiSlcReader::from_readers`] after opening one
/// [`SlcReader`] per slice.  The slices must be provided in the same temporal
/// order as the [`crate::slice_assembly::AssembledScene::safe_paths`] list.
#[derive(Debug)]
pub struct MultiSlcReader {
    /// One reader per slice, in slice order.
    readers: Vec<SlcReader>,
    /// `slice_line_offsets[k]` = sum of `readers[0..k].height()`.
    /// Length = `readers.len()`.  Entry 0 is always 0.
    slice_line_offsets: Vec<usize>,
    /// Sum of all TIFF heights (logical height of the assembled channel).
    logical_height: u32,
    /// Range samples per line (equal for all slices; validated on construction).
    width: u32,
}

// ─── Error type ───────────────────────────────────────────────────────────────

/// Errors produced exclusively during [`MultiSlcReader`] construction.
///
/// Read-time errors use the existing [`SlcReadError`] variants.
#[derive(Debug, thiserror::Error)]
pub enum MultiSlcOpenError {
    /// No readers were provided.
    #[error("MultiSlcReader requires at least one slice reader")]
    EmptyInput,

    /// Two slices have different TIFF widths; all slices of the same sub-swath
    /// must have the same number of range samples.
    #[error("slice {index} width {found} != slice 0 width {expected}")]
    WidthMismatch {
        index: usize,
        expected: u32,
        found: u32,
    },
}

// ─── Construction ─────────────────────────────────────────────────────────────

impl MultiSlcReader {
    /// Build a [`MultiSlcReader`] from pre-opened slice readers.
    ///
    /// # Arguments
    ///
    /// * `readers` — One open [`SlcReader`] per slice, in temporal order
    ///   (earliest slice first).  All must have the same TIFF width.
    ///
    /// # Errors
    ///
    /// Returns [`MultiSlcOpenError`] if:
    /// - `readers` is empty.
    /// - Any reader's width differs from `readers[0].width`.
    pub fn from_readers(readers: Vec<SlcReader>) -> Result<Self, MultiSlcOpenError> {
        if readers.is_empty() {
            return Err(MultiSlcOpenError::EmptyInput);
        }

        let width = readers[0].width;
        for (i, r) in readers.iter().enumerate().skip(1) {
            if r.width != width {
                return Err(MultiSlcOpenError::WidthMismatch {
                    index: i,
                    expected: width,
                    found: r.width,
                });
            }
        }

        // Build cumulative line offsets.  offset[k] = Σ height(readers[0..k]).
        let mut slice_line_offsets = Vec::with_capacity(readers.len());
        let mut cumulative: usize = 0;
        for r in &readers {
            slice_line_offsets.push(cumulative);
            cumulative = cumulative
                .checked_add(r.height as usize)
                .expect("total TIFF height overflows usize"); // SAFETY-OK: S-1 IW stacks are O(10^5) lines; usize overflow is physically impossible
        }

        let logical_height = cumulative as u32;

        Ok(Self {
            readers,
            slice_line_offsets,
            logical_height,
            width,
        })
    }
}

// ─── BurstReader impl ─────────────────────────────────────────────────────────

impl BurstReader for MultiSlcReader {
    fn width(&self) -> u32 {
        self.width
    }

    fn height(&self) -> u32 {
        self.logical_height
    }

    /// Read `line_count` azimuth lines from the assembled (logical) line space.
    ///
    /// # Dispatch
    ///
    /// 1. Binary-search `slice_line_offsets` to find the owning slice `k` such
    ///    that `offset[k] <= logical_first < offset[k+1]`.
    /// 2. Translate: `local_line = logical_first - offset[k]`.
    /// 3. Guard: the read must not cross into slice `k+1`.
    /// 4. Delegate to `readers[k].read_burst_raw(local_line, line_count)`.
    ///
    /// # Errors
    ///
    /// - [`SlcReadError::LineOutOfBounds`] if `logical_first + line_count > logical_height`.
    /// - [`SlcReadError::CrossSliceBoundary`] if the read would span a slice boundary.
    /// - [`SlcReadError::Io`] on file read failure from the owning slice reader.
    fn read_burst_raw(
        &mut self,
        logical_first: usize,
        line_count: usize,
    ) -> Result<Vec<[i16; 2]>, SlcReadError> {
        let logical_end = logical_first
            .checked_add(line_count)
            .filter(|&e| e <= self.logical_height as usize)
            .ok_or(SlcReadError::LineOutOfBounds {
                first: logical_first,
                end: logical_first.saturating_add(line_count),
                height: self.logical_height,
            })?;

        // Binary-search for the owning slice: largest k where offset[k] <= logical_first.
        let slice_idx = self
            .slice_line_offsets
            .partition_point(|&off| off <= logical_first)
            .saturating_sub(1); // SAFETY-OK: partition_point returns a value in [0, len]; saturating_sub(1) is 0 when logical_first < offset[0], which cannot happen since offset[0] == 0 and logical_first is non-negative

        let slice_start = self.slice_line_offsets[slice_idx];
        let slice_end = self
            .slice_line_offsets
            .get(slice_idx + 1)
            .copied()
            .unwrap_or(self.logical_height as usize); // SAFETY-OK: last slice has no next entry; logical_height is the correct upper bound

        // Guard: the read must not cross into the next slice.
        if logical_end > slice_end {
            return Err(SlcReadError::CrossSliceBoundary {
                first: logical_first,
                end: logical_end,
                slice_end,
            });
        }

        let local_first = logical_first - slice_start;
        self.readers[slice_idx].read_burst_raw(local_first, line_count)
    }
}

// ─── Tests ────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    /// Build a minimal in-memory-like synthetic `SlcReader` is not feasible
    /// since SlcReader requires a real file.  These tests use the on-disk S1A
    /// measurement TIFF (IW1/VV) as a single real reader, then split it in two
    /// halves to simulate a two-slice assembly and verify the dispatch logic.

    const S1A_IW1_VV: &str =
        "/home/datacube/dev/SARdine/data/SLC/\
         S1A_IW_SLC__1SDV_20201005T170824_20201005T170851_034664_04098A_1E66.SAFE/\
         measurement/\
         s1a-iw1-slc-vv-20201005t170824-20201005t170849-034664-04098a-004.tiff";

    fn open_iw1() -> SlcReader {
        SlcReader::open(S1A_IW1_VV).expect("open failed")
    }

    fn skip_if_no_tiff() -> bool {
        !std::path::Path::new(S1A_IW1_VV).exists()
    }

    // ── Construction ─────────────────────────────────────────────────────────

    #[test]
    fn empty_input_rejected() {
        let err = MultiSlcReader::from_readers(vec![]).unwrap_err();
        assert!(matches!(err, MultiSlcOpenError::EmptyInput));
    }

    #[test]
    fn single_reader_passthrough_width_height() {
        if skip_if_no_tiff() {
            return;
        }
        let r = open_iw1();
        let expected_w = r.width;
        let expected_h = r.height;
        let mut mr = MultiSlcReader::from_readers(vec![r]).unwrap();
        assert_eq!(mr.width(), expected_w);
        assert_eq!(mr.height(), expected_h);
        // Reading first 10 lines via MultiSlcReader must equal direct SlcReader read.
        let mut direct = open_iw1();
        let via_multi = mr.read_burst_raw(0, 10).unwrap();
        let via_direct = direct.read_burst_raw(0, 10).unwrap();
        assert_eq!(via_multi, via_direct, "data mismatch on first 10 lines");
    }

    #[test]
    fn two_readers_height_is_sum() {
        if skip_if_no_tiff() {
            return;
        }
        let r0 = open_iw1();
        let r1 = open_iw1();
        let h_each = r0.height;
        let mr = MultiSlcReader::from_readers(vec![r0, r1]).unwrap();
        assert_eq!(mr.height(), h_each * 2);
    }

    #[test]
    fn line_out_of_bounds_rejected() {
        if skip_if_no_tiff() {
            return;
        }
        let r = open_iw1();
        let h = r.height as usize;
        let mut mr = MultiSlcReader::from_readers(vec![r]).unwrap();
        let err = mr.read_burst_raw(h, 1).unwrap_err();
        assert!(
            matches!(err, SlcReadError::LineOutOfBounds { .. }),
            "expected LineOutOfBounds, got: {err}"
        );
    }

    /// With two identical readers stacked: reading lines in the second slice
    /// (at logical offset `height`) must return the same data as reading the
    /// same lines directly from the first TIFF file (because both readers point
    /// to the same file in this test fixture).
    #[test]
    fn second_slice_dispatches_correctly() {
        if skip_if_no_tiff() {
            return;
        }
        let r0 = open_iw1();
        let r1 = open_iw1();
        let h = r0.height as usize;
        let mut mr = MultiSlcReader::from_readers(vec![r0, r1]).unwrap();

        // Read 5 lines starting at logical line `h` (== first line of slice 1).
        // Since both readers are the same TIFF, this should equal lines [0,5).
        let via_multi = mr.read_burst_raw(h, 5).unwrap();
        let mut direct = open_iw1();
        let via_direct = direct.read_burst_raw(0, 5).unwrap();
        assert_eq!(via_multi, via_direct, "dispatch to second reader failed");
    }

    /// A read that straddles the slice boundary must be rejected.
    #[test]
    fn cross_boundary_read_rejected() {
        if skip_if_no_tiff() {
            return;
        }
        let r0 = open_iw1();
        let r1 = open_iw1();
        let h = r0.height as usize;
        let mut mr = MultiSlcReader::from_readers(vec![r0, r1]).unwrap();

        // Read starting 2 lines before the boundary, requesting 5 lines:
        // [h-2, h+3) crosses the boundary at h.
        let err = mr.read_burst_raw(h - 2, 5).unwrap_err();
        assert!(
            matches!(err, SlcReadError::CrossSliceBoundary { .. }),
            "expected CrossSliceBoundary, got: {err}"
        );
    }
}
