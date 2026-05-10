//! IW subswath merge for Sentinel-1 TOPS σ⁰ products.
//!
//! Combines the calibrated σ⁰ arrays from the three IW subswaths (IW1, IW2, IW3)
//! into a single contiguous slant-range image.
//!
//! # Performance
//!
//! The merge allocates two output buffers totalling several GB.  To avoid
//! serialising millions of OS page faults on a single thread, both buffers are
//! left uninitialised after allocation and the fill+copy is done in parallel
//! with Rayon.  Every element of every output line is written by the parallel
//! loop before the buffers are returned, so the `unsafe { set_len }` calls are
//! sound.
//!
//! # Range geometry
//!
//! All three IW subswaths share **identical range pixel spacing** (confirmed on
//! S1A/S1B real data: 2.329562 m).  This means the merge is pure integer-offset
//! arithmetic: no range resampling is needed.
//!
//! The three subswath grids are position-shifted copies of the same pixel spacing,
//! but their near-range slant times are not aligned to an integer number of pixels
//! relative to each other.  We round the inter-subswath offset to the nearest
//! integer pixel.  The resulting positioning error is < 0.5 pixel (≈ 1.2 m in
//! slant range) — well below any incoherent σ⁰ spatial requirement.
//!
//! # Overlap handling
//!
//! Adjacent subswaths overlap by ~2100–2200 pixels in slant range.  We apply a
//! **midpoint hard cut**: each output column is assigned to exactly one subswath,
//! with the seam placed at the centre of the overlap region.  No blending is
//! applied — blending raw σ⁰ values from different squint angles is physically
//! incorrect for slant-range imagery.
//!
//! # Azimuth alignment
//!
//! After deburst, the three subswaths have slightly different line counts
//! (because `linesPerBurst` and `overlapLines` differ).  This module clips all
//! inputs to the **minimum line count** across all input swaths.  This is valid
//! for incoherent σ⁰ because the subswaths are approximately co-registered in
//! azimuth (azimuth start-time differences are < 1 line ≈ 14 ms).
//!
//! # NaN fill
//!
//! Output pixels where no input subswath provides valid data (gap or seam
//! outside all swath extents) are filled with `f32::NAN`.

use crate::apply_calibration::Sigma0Array;
use crate::types::SubSwathMetadata;
use chrono::{DateTime, Utc};
use rayon::prelude::*;
use thiserror::Error;

// ─── constants ────────────────────────────────────────────────────────────────

const SPEED_OF_LIGHT_M_S: f64 = 299_792_458.0;

// ─── Error type ───────────────────────────────────────────────────────────────

/// Errors produced by [`merge_subswaths`].
#[derive(Debug, Error)]
pub enum MergeError {
    /// At least one subswath input is required.
    #[error("at least 1 subswath input is required, got 0")]
    EmptyInputs,

    /// The range pixel spacings of the input swaths are not consistent.
    /// All S-1 IW subswaths should share the same value.
    #[error(
        "range pixel spacing inconsistency: swath 0 = {first:.6e} m, swath {idx} = {found:.6e} m"
    )]
    InconsistentPixelSpacing {
        first: f64,
        idx: usize,
        found: f64,
    },

    /// A Sigma0Array is empty (zero samples).
    #[error("subswath {idx} has zero range samples")]
    EmptySwath { idx: usize },
}

// ─── Types ────────────────────────────────────────────────────────────────────

/// One input to the merge: a calibrated σ⁰ array paired with its subswath metadata.
pub struct SwathInput<'a> {
    /// Calibrated σ⁰ output from [`crate::apply_calibration::apply_calibration`].
    pub sigma0: &'a Sigma0Array,
    /// Annotation metadata for this subswath (provides slant range time +
    /// pixel spacing).
    pub swath: &'a SubSwathMetadata,
    /// Zero-Doppler UTC time of debursted line 0 for *this* subswath.
    ///
    /// For a Sentinel-1 IW SLC debursted in the standard way, this is the
    /// annotation `azimuthTime` of `bursts[0]` (S-1 convention: the
    /// zero-Doppler sensing time of the first line of the first burst of
    /// the subswath).  The caller must look this up on the per-subswath
    /// burst list of the owning [`SceneMetadata`] and pass it in
    /// explicitly, rather than the merge function reaching through
    /// `swath` for it — `SubSwathMetadata` itself does not carry bursts.
    pub azimuth_start_time: DateTime<Utc>,
}

/// Merged σ⁰ image covering all input subswaths in slant range.
#[derive(Debug)]
pub struct MergedSigma0 {
    /// Flat row-major buffer.  Length = `lines × samples`.
    ///
    /// Pixels with no valid input swath coverage are `f32::NAN`.
    pub data: Vec<f32>,

    /// Number of azimuth lines (= minimum line count across all input swaths).
    pub lines: usize,

    /// Number of range samples in the merged output.
    pub samples: usize,

    /// Two-way slant range time of output column 0, in seconds.
    ///
    /// Column `c` has slant time `near_slant_range_time_s + c * range_pixel_spacing_s`.
    pub near_slant_range_time_s: f64,

    /// Two-way range time per pixel, in seconds.  Equal to
    /// `range_pixel_spacing_m * 2 / c`.
    pub range_pixel_spacing_s: f64,

    /// Range pixel spacing in metres (slant range, inherited from inputs).
    pub range_pixel_spacing_m: f64,

    /// Maximum calibration LUT extrapolation gap (in range pixels) observed
    /// across all merged subswaths.  Propagated from
    /// [`Sigma0Array::cal_lut_extrapolation_gap_px`] so the caller can audit
    /// LUT coverage of the merged product without losing information at the
    /// merge boundary.
    pub cal_lut_extrapolation_gap_px: u32,

    /// Maximum range-noise LUT extrapolation gap (in range pixels) observed
    /// across all merged subswaths.
    pub noise_lut_extrapolation_gap_px: u32,

    /// Per-pixel NESZ in linear σ⁰ units, merged from input subswaths.
    ///
    /// Same shape as `data` (flat row-major, `lines × samples`).
    /// NaN where no input swath covers the pixel (seam gaps, swath edges).
    pub nesz: Vec<f32>,

    /// Zero-Doppler UTC time of merged line 0.
    ///
    /// Set to the **earliest** `bursts[0].azimuth_time_utc` among the input
    /// subswaths: debursted line 0 of each input corresponds to its first
    /// burst's line 0 (annotation `azimuthTime`), and by convention the
    /// merged image's line-0 clock is anchored to the earliest such time so
    /// that downstream per-pixel geocoding references the physically
    /// earliest sample in the mosaic.
    ///
    /// Subswaths whose `bursts[0].azimuth_time_utc` is later than this
    /// anchor are treated as line-aligned to within ≤ 1 azimuth line
    /// (TOPS burst-cycle synchronisation across IW subswaths).  The per-
    /// pixel zero-Doppler Newton solve downstream absorbs this sub-line
    /// offset.
    pub azimuth_start_time: DateTime<Utc>,
}

// ─── Main function ────────────────────────────────────────────────────────────

/// Single-subswath passthrough that **moves** `sigma0` data buffers into
/// [`MergedSigma0`] without cloning.
///
/// Callers that hold an owned [`Sigma0Array`] and process only one subswath
/// (e.g. `--iw IW2` mode) should prefer this over [`merge_subswaths`] to
/// avoid O(n) copies of `data` and `nesz`.
pub fn merge_single_subswath_owned(
    sigma0: crate::apply_calibration::Sigma0Array,
    swath: &SubSwathMetadata,
    azimuth_start_time: DateTime<Utc>,
) -> Result<MergedSigma0, MergeError> {
    if sigma0.samples == 0 {
        return Err(MergeError::EmptySwath { idx: 0 });
    }
    let rps_m = swath.range_pixel_spacing_m;
    let dt = rps_m * 2.0 / SPEED_OF_LIGHT_M_S;
    let near_slant_range_time_s =
        swath.slant_range_time_s + sigma0.valid_sample_offset as f64 * dt;
    Ok(MergedSigma0 {
        data: sigma0.data,
        nesz: sigma0.nesz,
        lines: sigma0.lines,
        samples: sigma0.samples,
        near_slant_range_time_s,
        range_pixel_spacing_s: dt,
        range_pixel_spacing_m: rps_m,
        cal_lut_extrapolation_gap_px: sigma0.cal_lut_extrapolation_gap_px,
        noise_lut_extrapolation_gap_px: sigma0.noise_lut_extrapolation_gap_px,
        azimuth_start_time,
    })
}

/// Merge calibrated σ⁰ arrays from adjacent IW subswaths into one slant-range image.
///
/// `inputs` must be ordered by increasing near-range slant time (i.e. IW1, IW2, IW3
/// for a standard VV/VH product).  Passing swaths out of order produces an incorrect
/// seam placement but is not detected at runtime.
///
/// # Panics
///
/// Does not panic; all validation is returned as [`MergeError`].
pub fn merge_subswaths(inputs: &[SwathInput<'_>]) -> Result<MergedSigma0, MergeError> {
    if inputs.is_empty() {
        return Err(MergeError::EmptyInputs);
    }

    // ── Single-subswath passthrough (no overlap or seam computation needed) ────
    if inputs.len() == 1 {
        let inp = &inputs[0];
        if inp.sigma0.samples == 0 {
            return Err(MergeError::EmptySwath { idx: 0 });
        }
        let rps_m = inp.swath.range_pixel_spacing_m;
        let dt = rps_m * 2.0 / SPEED_OF_LIGHT_M_S;
        let near_slant_range_time_s =
            inp.swath.slant_range_time_s + inp.sigma0.valid_sample_offset as f64 * dt;
        return Ok(MergedSigma0 {
            data: inp.sigma0.data.clone(),
            nesz: inp.sigma0.nesz.clone(),
            lines: inp.sigma0.lines,
            samples: inp.sigma0.samples,
            near_slant_range_time_s,
            range_pixel_spacing_s: dt,
            range_pixel_spacing_m: rps_m,
            cal_lut_extrapolation_gap_px: inp.sigma0.cal_lut_extrapolation_gap_px,
            noise_lut_extrapolation_gap_px: inp.sigma0.noise_lut_extrapolation_gap_px,
            azimuth_start_time: inp.azimuth_start_time,
        });
    }

    // ── Validate consistent pixel spacing ─────────────────────────────────────

    let rps_m = inputs[0].swath.range_pixel_spacing_m;
    for (idx, sw) in inputs.iter().enumerate().skip(1) {
        let rps = sw.swath.range_pixel_spacing_m;
        // Tolerance: 1 mm (hardware precision >> this).
        if (rps - rps_m).abs() > 1e-3 {
            return Err(MergeError::InconsistentPixelSpacing {
                first: rps_m,
                idx,
                found: rps,
            });
        }
    }

    // ── Validate non-empty inputs ──────────────────────────────────────────────

    for (idx, inp) in inputs.iter().enumerate() {
        if inp.sigma0.samples == 0 {
            return Err(MergeError::EmptySwath { idx });
        }
    }

    // ── Geometry calculation ───────────────────────────────────────────────────

    // Two-way range time step [s] per pixel.
    let dt = rps_m * 2.0 / SPEED_OF_LIGHT_M_S;

    // Slant time of sigma0 column 0 for each input swath.
    // sigma0 column 0 = TIFF column valid_sample_offset, so:
    //   t_sw_start = t0_sw + valid_sample_offset * dt
    let sw_start_times: Vec<f64> = inputs
        .iter()
        .map(|inp| {
            inp.swath.slant_range_time_s
                + inp.sigma0.valid_sample_offset as f64 * dt
        })
        .collect();

    // Overall merged output starts at the earliest (leftmost) subswath valid near edge.
    let t_out_start = sw_start_times[0];

    // Output column offset for each input swath:
    //   out_col(in_col) = round((t_sw_start - t_out_start) / dt) + in_col
    // So the starting output column for each swath (sigma0 col 0) is:
    let sw_out_offsets: Vec<i64> = sw_start_times
        .iter()
        .map(|&t| ((t - t_out_start) / dt).round() as i64)
        .collect();

    // Far (exclusive) output column for each swath:
    //   sw_out_end = sw_out_offset + sigma0.samples
    let sw_out_ends: Vec<i64> = inputs
        .iter()
        .zip(sw_out_offsets.iter())
        .map(|(inp, &off)| off + inp.sigma0.samples as i64)
        .collect();

    // Total output column count.
    let out_cols = *sw_out_ends.last().unwrap() as usize;

    // ── Seam positions (midpoint of overlap between adjacent swaths) ───────────
    //
    // The seam is the first output column assigned to the right swath.
    // In the overlap region [IW(i)_far, IW(i+1)_near) in time, both swaths
    // provide valid data; the seam is placed at the midpoint.
    //
    // For non-overlapping swaths (positive gap), the seam is at sw_out_end[i].
    let mut seams: Vec<usize> = Vec::with_capacity(inputs.len() - 1);
    for i in 0..inputs.len() - 1 {
        let left_end  = sw_out_ends[i];       // exclusive end of left swath
        let right_start = sw_out_offsets[i + 1]; // start of right swath
        // midpoint (rounded toward left half to keep integer arithmetic clean)
        let seam = ((left_end + right_start) / 2) as usize;
        seams.push(seam);
    }

    // ── Azimuth: clip to minimum line count ────────────────────────────────────

    let out_lines = inputs.iter().map(|inp| inp.sigma0.lines).min().unwrap();

    // ── Per-swath column fill ranges (precomputed once, used by every line) ───
    //
    // For swath `sw_idx`:
    //   owned output columns: [col_start, col_end)
    //   corresponding input columns: [in_start, in_start + n_copy)
    //
    // `valid_out_start/end` clips the owned range to where the swath actually
    // provides data (handles edge cases where a seam falls outside the swath's
    // valid sample window; in normal IW processing this never happens but we
    // check for correctness).
    struct SwathFillRange {
        /// First owned output column (fill NaN from previous fill_end up to here).
        col_start: usize,
        /// One-past-last owned output column (fill NaN from fill_end to next col_start).
        col_end: usize,
        /// First output column that has valid input data from this swath (≥ col_start).
        valid_out_start: usize,
        /// One-past-last output column with valid input data (≤ col_end).
        valid_out_end: usize,
        /// Input sigma0 column corresponding to `valid_out_start`.
        in_start: usize,
    }

    let fill_ranges: Vec<SwathFillRange> = inputs
        .iter()
        .enumerate()
        .map(|(sw_idx, inp)| {
            let col_start: usize = if sw_idx == 0 { 0 } else { seams[sw_idx - 1] };
            let col_end: usize = if sw_idx == inputs.len() - 1 {
                out_cols
            } else {
                seams[sw_idx]
            };
            let off = sw_out_offsets[sw_idx] as usize; // all offsets are ≥ 0 (IW1 anchors t_out_start)
            let valid_out_start = col_start.max(off);
            let valid_out_end = col_end.min(off.saturating_add(inp.sigma0.samples));
            let in_start = valid_out_start.saturating_sub(off);
            SwathFillRange { col_start, col_end, valid_out_start, valid_out_end, in_start }
        })
        .collect();

    // ── Allocate output buffers ────────────────────────────────────────────────
    //
    // Use uninitialised memory: the parallel fill loop below writes every
    // element before the buffers are returned.  This avoids the serialised
    // OS page-fault cost of `vec![f32::NAN; N]` for several-GB allocations.

    let n = out_lines * out_cols;
    let mut data: Vec<f32> = Vec::with_capacity(n);
    let mut nesz: Vec<f32> = Vec::with_capacity(n);
    // SAFETY-OK: every element of `data` and `nesz` is written by the
    // `par_chunks_mut` loop below (each line is either copy_from_slice'd from
    // a swath or fill'd with NaN for gap regions) before either Vec is read
    // outside this function.  No element is read before it is written.
    unsafe {
        data.set_len(n);
        nesz.set_len(n);
    }

    // ── Fill output lines in parallel ──────────────────────────────────────────
    //
    // Each line is independent: it reads from fixed sigma0 slices and writes to
    // a disjoint output chunk.  `par_chunks_mut` gives each Rayon worker a
    // disjoint slice; workers only read from `inputs` (shared immutable refs).

    data.par_chunks_mut(out_cols)
        .zip(nesz.par_chunks_mut(out_cols))
        .enumerate()
        .for_each(|(out_line, (out_row, nesz_row))| {
            let mut fill_end = 0usize; // tracks how far we have initialised

            for (sw_idx, inp) in inputs.iter().enumerate() {
                let fr = &fill_ranges[sw_idx];

                // Fill any gap before this swath's valid data with NaN.
                // In standard IW processing fr.valid_out_start == fr.col_start,
                // so this range is empty and no NaN fill is emitted.
                if fill_end < fr.valid_out_start {
                    out_row[fill_end..fr.valid_out_start].fill(f32::NAN);
                    nesz_row[fill_end..fr.valid_out_start].fill(f32::NAN);
                }

                if fr.valid_out_start < fr.valid_out_end {
                    let n_copy = fr.valid_out_end - fr.valid_out_start;
                    let base = out_line * inp.sigma0.samples + fr.in_start;
                    out_row[fr.valid_out_start..fr.valid_out_end]
                        .copy_from_slice(&inp.sigma0.data[base..base + n_copy]);
                    nesz_row[fr.valid_out_start..fr.valid_out_end]
                        .copy_from_slice(&inp.sigma0.nesz[base..base + n_copy]);
                }

                // Advance past the owned range (NaN in gap after valid data).
                let after_valid = fr.valid_out_end;
                if after_valid < fr.col_end {
                    out_row[after_valid..fr.col_end].fill(f32::NAN);
                    nesz_row[after_valid..fr.col_end].fill(f32::NAN);
                }

                fill_end = fr.col_end;
            }

            // Fill any trailing gap (e.g. non-contiguous swaths, final swath short).
            if fill_end < out_cols {
                out_row[fill_end..out_cols].fill(f32::NAN);
                nesz_row[fill_end..out_cols].fill(f32::NAN);
            }
        });

    Ok(MergedSigma0 {
        data,
        nesz,
        lines: out_lines,
        samples: out_cols,
        near_slant_range_time_s: t_out_start,
        range_pixel_spacing_s: dt,
        range_pixel_spacing_m: rps_m,
        cal_lut_extrapolation_gap_px: inputs
            .iter()
            .map(|s| s.sigma0.cal_lut_extrapolation_gap_px)
            .max()
            .unwrap_or(0), // SAFETY-OK: single-swath early return and EmptyInputs check guarantee inputs.len() >= 2 at this point; .max() over a non-empty iter is always Some; unwrap_or(0) is unreachable
        noise_lut_extrapolation_gap_px: inputs
            .iter()
            .map(|s| s.sigma0.noise_lut_extrapolation_gap_px)
            .max()
            .unwrap_or(0), // SAFETY-OK: same precondition as cal_lut_extrapolation_gap_px above
        azimuth_start_time: inputs
            .iter()
            .map(|s| s.azimuth_start_time)
            .min()
            .unwrap_or_else(|| inputs[0].azimuth_start_time), // SAFETY-OK: single-swath early return and EmptyInputs check guarantee inputs.len() >= 2 at this point; .min() always returns Some; this branch is unreachable at runtime
    })
}

// ─── Tests ────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use crate::apply_calibration::Sigma0Array;
    use crate::types::{SubSwathId, SubSwathMetadata};

    // ── Helper: build a minimal SubSwathMetadata ──────────────────────────────

    fn make_swath(t0_s: f64, rps_m: f64) -> SubSwathMetadata {
        SubSwathMetadata {
            id: SubSwathId::IW1,
            burst_count: 1,
            lines_per_burst: 1,
            range_samples: 10,
            azimuth_samples: 1,
            first_line: 0,
            last_line: 1,
            first_sample: 0,
            last_sample: 10,
            range_pixel_spacing_m: rps_m,
            azimuth_pixel_spacing_m: 14.0,
            slant_range_time_s: t0_s,
            azimuth_time_interval_s: 2e-3,
            prf_hz: 1717.13,
            burst_cycle_time_s: 2.758_091,
            dc_estimates: Vec::new(),
            fm_rates: Vec::new(),
        }
    }

    fn make_sigma0(data: Vec<f32>, samples: usize, valid_offset: usize) -> Sigma0Array {
        let lines = if samples == 0 { 0 } else { data.len() / samples };
        let nesz_len = lines * samples;
        Sigma0Array {
            data,
            lines,
            samples,
            valid_sample_offset: valid_offset,
            cal_lut_extrapolation_gap_px: 0,
            noise_lut_extrapolation_gap_px: 0,
            nesz: vec![0.0f32; nesz_len],
        }
    }

    // ── Error cases ───────────────────────────────────────────────────────────

    #[test]
    fn empty_inputs_rejected() {
        let result = merge_subswaths(&[]);
        assert!(matches!(result, Err(MergeError::EmptyInputs)));
    }

    #[test]
    fn single_swath_passthrough() {
        // A single subswath with 2 samples and 1 line should pass through directly.
        let swath = make_swath(0.005_000, 2.329_562);
        let s = make_sigma0(vec![1.0, 2.0], 2, 0);
        let result = merge_subswaths(&[SwathInput {
            sigma0: &s,
            swath: &swath,
            azimuth_start_time: chrono::Utc::now(),
        }]);
        let merged = result.expect("single swath should produce Ok");
        assert_eq!(merged.samples, 2);
        assert_eq!(merged.lines, 1);
        assert_eq!(merged.data, vec![1.0, 2.0]);
        assert_eq!(merged.range_pixel_spacing_m, 2.329_562);
    }

    #[test]
    fn inconsistent_pixel_spacing_rejected() {
        let sw1 = make_swath(0.005, 2.33);
        let sw2 = make_swath(0.006, 2.50);  // different spacing
        let s1 = make_sigma0(vec![1.0; 10], 10, 0);
        let s2 = make_sigma0(vec![1.0; 10], 10, 0);
        let result = merge_subswaths(&[
            SwathInput { sigma0: &s1, swath: &sw1, azimuth_start_time: chrono::Utc::now() },
            SwathInput { sigma0: &s2, swath: &sw2, azimuth_start_time: chrono::Utc::now() },
        ]);
        assert!(matches!(result, Err(MergeError::InconsistentPixelSpacing { .. })));
    }

    #[test]
    fn empty_swath_rejected() {
        let sw1 = make_swath(0.005, 2.33);
        let sw2 = make_swath(0.006, 2.33);
        let s1 = make_sigma0(vec![], 0, 0);
        let s2 = make_sigma0(vec![1.0; 10], 10, 0);
        let result = merge_subswaths(&[
            SwathInput { sigma0: &s1, swath: &sw1, azimuth_start_time: chrono::Utc::now() },
            SwathInput { sigma0: &s2, swath: &sw2, azimuth_start_time: chrono::Utc::now() },
        ]);
        assert!(matches!(result, Err(MergeError::EmptySwath { idx: 0 })));
    }

    // ── Geometry unit tests ───────────────────────────────────────────────────

    /// Two adjacent swaths, no overlap, back-to-back: output should be
    /// simple concatenation with no NaN.
    #[test]
    fn two_swaths_no_overlap_concatenated() {
        let rps = 2.329562_f64;
        let dt = rps * 2.0 / SPEED_OF_LIGHT_M_S;
        // sw1: t0 = 0, 10 samples  → covers output cols 0..10
        // sw2: t0 = 10 * dt (exact), 10 samples → covers output cols 10..20
        let t0_sw2 = 10.0 * dt;
        let sw1 = make_swath(0.0, rps);
        let sw2 = make_swath(t0_sw2, rps);
        let s1 = make_sigma0((1..=10).map(|v| v as f32).collect(), 10, 0);
        let s2 = make_sigma0((11..=20).map(|v| v as f32).collect(), 10, 0);

        let merged = merge_subswaths(&[
            SwathInput { sigma0: &s1, swath: &sw1, azimuth_start_time: chrono::Utc::now() },
            SwathInput { sigma0: &s2, swath: &sw2, azimuth_start_time: chrono::Utc::now() },
        ])
        .unwrap();

        assert_eq!(merged.samples, 20);
        assert_eq!(merged.lines, 1);
        let any_nan = merged.data.iter().any(|v| v.is_nan());
        assert!(!any_nan, "unexpected NaN in no-gap merge");
        // IW1 owns cols 0..10, IW2 owns cols 10..20
        assert_eq!(merged.data[0], 1.0);
        assert_eq!(merged.data[9], 10.0);
        assert_eq!(merged.data[10], 11.0);
        assert_eq!(merged.data[19], 20.0);
    }

    /// Two swaths with 4-pixel overlap: seam should be at midpoint (col 8 of a
    /// 20-sample output where left ends at col 12 and right starts at col 8).
    #[test]
    fn seam_at_midpoint_of_overlap() {
        let rps = 2.329562_f64;
        let dt = rps * 2.0 / SPEED_OF_LIGHT_M_S;
        // sw1: 12 samples starting at t=0  → tiff cols 0..12, valid [0..12]
        // sw2: 12 samples starting 8*dt later → tiff cols 8..20
        // overlap between sw1 col 8..12 (4 pixels)
        // seam = midpoint of (sw1_end=12, sw2_start=8) = (12+8)/2 = 10
        let t0_sw2 = 8.0 * dt;
        let sw1 = make_swath(0.0, rps);
        let sw2 = make_swath(t0_sw2, rps);
        // sw1 data: all 1.0; sw2 data: all 2.0
        let s1 = make_sigma0(vec![1.0; 12], 12, 0);
        let s2 = make_sigma0(vec![2.0; 12], 12, 0);

        let merged = merge_subswaths(&[
            SwathInput { sigma0: &s1, swath: &sw1, azimuth_start_time: chrono::Utc::now() },
            SwathInput { sigma0: &s2, swath: &sw2, azimuth_start_time: chrono::Utc::now() },
        ])
        .unwrap();

        // output spans t=0 to t=19*dt → 20 cols
        assert_eq!(merged.samples, 20);
        // col 0..9: sw1 → 1.0
        assert_eq!(merged.data[9], 1.0, "col 9 should come from sw1");
        // col 10..19: sw2 → 2.0  (note: sw2 starts at output col 8,
        //   so output col 10 = sw2 input col 10-8=2, value=2.0)
        assert_eq!(merged.data[10], 2.0, "col 10 should come from sw2");
    }

    /// Output lines clipped to minimum when inputs differ.
    #[test]
    fn azimuth_clipped_to_minimum_lines() {
        let rps = 2.329562_f64;
        let dt = rps * 2.0 / SPEED_OF_LIGHT_M_S;
        let sw1 = make_swath(0.0, rps);
        let sw2 = make_swath(10.0 * dt, rps);
        // sw1: 5 lines, sw2: 3 lines
        let s1 = make_sigma0(vec![1.0; 10 * 5], 10, 0);
        let s2 = make_sigma0(vec![2.0; 10 * 3], 10, 0);

        let merged = merge_subswaths(&[
            SwathInput { sigma0: &s1, swath: &sw1, azimuth_start_time: chrono::Utc::now() },
            SwathInput { sigma0: &s2, swath: &sw2, azimuth_start_time: chrono::Utc::now() },
        ])
        .unwrap();

        assert_eq!(merged.lines, 3);
    }

    /// Calibration / noise LUT extrapolation gaps must propagate to the
    /// merged product as the per-input maximum, so downstream callers can
    /// audit LUT coverage end-to-end without losing information at the
    /// merge boundary.
    #[test]
    fn lut_extrapolation_gaps_propagate_as_max() {
        let rps = 2.329562_f64;
        let dt = rps * 2.0 / SPEED_OF_LIGHT_M_S;
        let sw1 = make_swath(0.0, rps);
        let sw2 = make_swath(10.0 * dt, rps);
        let mut s1 = make_sigma0(vec![1.0; 10 * 3], 10, 0);
        let mut s2 = make_sigma0(vec![2.0; 10 * 3], 10, 0);
        s1.cal_lut_extrapolation_gap_px = 17;
        s1.noise_lut_extrapolation_gap_px = 0;
        s2.cal_lut_extrapolation_gap_px = 5;
        s2.noise_lut_extrapolation_gap_px = 42;

        let merged = merge_subswaths(&[
            SwathInput { sigma0: &s1, swath: &sw1, azimuth_start_time: chrono::Utc::now() },
            SwathInput { sigma0: &s2, swath: &sw2, azimuth_start_time: chrono::Utc::now() },
        ])
        .unwrap();

        assert_eq!(merged.cal_lut_extrapolation_gap_px, 17);
        assert_eq!(merged.noise_lut_extrapolation_gap_px, 42);
    }

    /// valid_sample_offset shifts the sigma0 column mapping correctly.
    #[test]
    fn valid_sample_offset_applied_correctly() {
        // sw1: t0=0, 20 samples TIFF, but valid window is col 5..15 (10 valid cols)
        // sw2: starts at t0 = TIFF col 10 * dt from sw1 origin, 10 samples, offset 0
        // The first sigma0 col of sw1 maps to TIFF col 5, slant time 5*dt.
        let rps = 2.329562_f64;
        let dt = rps * 2.0 / SPEED_OF_LIGHT_M_S;
        let sw1 = make_swath(0.0, rps);                // sw1 t0 = 0
        let sw2 = make_swath(10.0 * dt, rps);          // sw2 t0 = 10*dt

        // sw1 sigma0: 10 samples with valid_sample_offset=5
        //   → covers slant times [5*dt, 14*dt]
        let s1 = make_sigma0((1..=10).map(|v| v as f32).collect(), 10, 5);
        // sw2 sigma0: 10 samples with valid_sample_offset=0
        //   → covers slant times [10*dt, 19*dt]
        let s2 = make_sigma0((101..=110).map(|v| v as f32).collect(), 10, 0);

        let merged = merge_subswaths(&[
            SwathInput { sigma0: &s1, swath: &sw1, azimuth_start_time: chrono::Utc::now() },
            SwathInput { sigma0: &s2, swath: &sw2, azimuth_start_time: chrono::Utc::now() },
        ])
        .unwrap();

        // Output starts at sw1 valid near edge: t=5*dt → 15 output cols (5..19)
        assert_eq!(merged.samples, 15, "unexpected output width");

        // The seam is at midpoint of (sw1_end=10, sw2_start=5) → (10+5)/2=7
        // in OUTPUT column space:
        //   sw1_end_outcol = 10 - 5 = 5 (sw1 has 10 samples, offset 5, starts at out_col 0, ends at 10)
        // Wait, let me trace this carefully:
        // t_out_start = t0_sw1 + valid_offset_sw1 * dt = 0 + 5*dt = 5*dt
        // sw_out_offsets[0] = round((5*dt - 5*dt) / dt) = 0  ← sw1 starts at out_col 0
        // sw_out_offsets[1] = round((10*dt + 0*dt - 5*dt) / dt) = round(5) = 5  ← sw2 starts at out_col 5
        // sw_out_ends[0] = 0 + 10 = 10
        // sw_out_ends[1] = 5 + 10 = 15
        // seam = (10 + 5) / 2 = 7  (integer division of 15)
        // sw1 owns cols [0..7), sw2 owns cols [7..15)
        assert_eq!(merged.data[0], 1.0, "out_col 0 → sw1 in_col 0 = 1.0");
        assert_eq!(merged.data[6], 7.0, "out_col 6 → sw1 in_col 6 = 7.0");
        // out_col 7 → sw2 in_col = 7-5 = 2 → value 103.0
        assert_eq!(merged.data[7], 103.0, "out_col 7 → sw2 in_col 2 = 103.0");
        assert_eq!(merged.data[14], 110.0, "last col → sw2 in_col 9 = 110.0");
    }

    /// Near-range time of output is correctly set.
    #[test]
    fn near_range_time_matches_sw1_valid_near() {
        let rps = 2.329562_f64;
        let dt = rps * 2.0 / SPEED_OF_LIGHT_M_S;
        let sw1 = make_swath(0.005, rps);
        let sw2 = make_swath(0.005 + 15.0 * dt, rps);
        let s1 = make_sigma0(vec![1.0; 10], 10, 3);
        let s2 = make_sigma0(vec![2.0; 10], 10, 0);

        let merged = merge_subswaths(&[
            SwathInput { sigma0: &s1, swath: &sw1, azimuth_start_time: chrono::Utc::now() },
            SwathInput { sigma0: &s2, swath: &sw2, azimuth_start_time: chrono::Utc::now() },
        ])
        .unwrap();

        let expected_t0 = 0.005 + 3.0 * dt;
        let diff = (merged.near_slant_range_time_s - expected_t0).abs();
        assert!(diff < 1e-15, "near_slant_range_time_s off by {diff:.2e}");
    }

    // ── Integration test against real S1A data ────────────────────────────────

    #[test]
    fn merge_s1a_iw_vv_smoke_test() {
        use crate::apply_calibration::apply_calibration;
        use crate::calibration::parse_calibration_noise;
        use crate::parse::parse_safe_directory;
        use crate::slc_reader::SlcReader;
        use crate::types::Polarization;

        const S1A_SAFE: &str =
            "/home/datacube/dev/SARdine/data/SLC/\
             S1A_IW_SLC__1SDV_20201005T170824_20201005T170851_034664_04098A_1E66.SAFE";
        const TIFF_IW: [&str; 3] = [
            "/home/datacube/dev/SARdine/data/SLC/\
             S1A_IW_SLC__1SDV_20201005T170824_20201005T170851_034664_04098A_1E66.SAFE/\
             measurement/s1a-iw1-slc-vv-20201005t170824-20201005t170849-034664-04098a-004.tiff",
            "/home/datacube/dev/SARdine/data/SLC/\
             S1A_IW_SLC__1SDV_20201005T170824_20201005T170851_034664_04098A_1E66.SAFE/\
             measurement/s1a-iw2-slc-vv-20201005t170824-20201005t170850-034664-04098a-005.tiff",
            "/home/datacube/dev/SARdine/data/SLC/\
             S1A_IW_SLC__1SDV_20201005T170824_20201005T170851_034664_04098A_1E66.SAFE/\
             measurement/s1a-iw3-slc-vv-20201005t170825-20201005t170851-034664-04098a-006.tiff",
        ];
        use crate::SubSwathId;
        const IW_IDS: [SubSwathId; 3] = [SubSwathId::IW1, SubSwathId::IW2, SubSwathId::IW3];

        if !std::path::Path::new(S1A_SAFE).is_dir() {
            eprintln!("skipping merge_s1a_iw_vv_smoke_test — S1A SAFE not present");
            return;
        }
        let scene = parse_safe_directory(std::path::Path::new(S1A_SAFE))
            .expect("parse SAFE");
        let cal_data = parse_calibration_noise(std::path::Path::new(S1A_SAFE))
            .expect("parse calibration");

        // Process each subswath — only burst 4 (mid-scene) to keep the test
        // quick in debug mode.
        let mut sigma0_arrays: Vec<Sigma0Array> = Vec::new();
        let mut swath_metas: Vec<&crate::types::SubSwathMetadata> = Vec::new();
        let mut az_start_times: Vec<DateTime<Utc>> = Vec::new();

        for (iw_idx, &iw_id) in IW_IDS.iter().enumerate() {
            let sw = scene.sub_swaths.iter().find(|s| s.id == iw_id).unwrap();
            let bursts: Vec<_> = scene.bursts.iter()
                .filter(|b| b.subswath_id == iw_id)
                .cloned()
                .collect();
            let burst4 = bursts.iter().find(|b| b.burst_index == 4).unwrap();

            let cal = cal_data.calibrations.iter()
                .find(|c| c.subswath_id == iw_id && c.polarization == Polarization::VV)
                .unwrap();
            let noise = cal_data.noises.iter()
                .find(|n| n.subswath_id == iw_id && n.polarization == Polarization::VV)
                .unwrap();

            // Read + crop burst 4 to valid range window (same as smoke test in G).
            let first_valid = bursts.iter().map(|b| b.first_valid_sample).max().unwrap();
            let last_valid  = bursts.iter().map(|b| b.last_valid_sample ).min().unwrap();
            let valid_cols  = last_valid - first_valid;

            let mut reader = SlcReader::open(TIFF_IW[iw_idx]).unwrap();
            let raw_cols = reader.width as usize;
            let raw = reader.read_burst_raw(burst4.first_line, sw.lines_per_burst).unwrap();

            let mut raw_data = Vec::with_capacity(sw.lines_per_burst * valid_cols);
            for line in 0..sw.lines_per_burst {
                raw_data.extend_from_slice(
                    &raw[line * raw_cols + first_valid..line * raw_cols + last_valid],
                );
            }
            let debursted = crate::deburst::DeburstArray {
                data: raw_data,
                lines: sw.lines_per_burst,
                samples: valid_cols,
                valid_sample_offset: first_valid,
            };

            let sigma0 = apply_calibration(&debursted, cal, noise, burst4.first_line)
                .expect("calibration failed");

            sigma0_arrays.push(sigma0);
            swath_metas.push(sw);
            az_start_times.push(burst4.azimuth_time_utc);
        }

        // Merge.
        let inputs: Vec<SwathInput<'_>> = sigma0_arrays.iter()
            .zip(swath_metas.iter())
            .zip(az_start_times.iter())
            .map(|((s, sw), &t0)| SwathInput {
                sigma0: s,
                swath: sw,
                azimuth_start_time: t0,
            })
            .collect();

        let merged = merge_subswaths(&inputs).expect("merge failed");

        // Shape sanity.
        assert_eq!(merged.data.len(), merged.lines * merged.samples);
        // Output should be wider than any individual subswath.
        let max_sw_width = sigma0_arrays.iter().map(|s| s.samples).max().unwrap();
        assert!(
            merged.samples > max_sw_width,
            "merged width {} not wider than largest subswath {}",
            merged.samples, max_sw_width
        );

        // Near-range time should equal IW1's valid start.
        let dt = swath_metas[0].range_pixel_spacing_m * 2.0 / SPEED_OF_LIGHT_M_S;
        let expected_near = swath_metas[0].slant_range_time_s
            + sigma0_arrays[0].valid_sample_offset as f64 * dt;
        let diff = (merged.near_slant_range_time_s - expected_near).abs();
        assert!(diff < 1e-12, "near range time off by {diff:.2e} s");

        // Most pixels should be finite and non-negative where covered.
        let n_valid = merged.data.iter().filter(|v| !v.is_nan()).count();
        let frac_valid = n_valid as f64 / merged.data.len() as f64;
        assert!(
            frac_valid > 0.85,
            "less than 85% of merged pixels are valid: {:.1}%",
            frac_valid * 100.0
        );
        let any_neg = merged.data.iter().any(|&v| !v.is_nan() && v < 0.0);
        assert!(!any_neg, "negative σ⁰ in merged output");

        // Range pixel spacing should be preserved.
        assert!((merged.range_pixel_spacing_m - 2.329562).abs() < 1e-4);
    }
}
