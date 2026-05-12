//! TOPS-mode deburst for Sentinel-1 IW/EW single-look complex data.
//!
//! Merges the per-burst SLC strips from one subswath into a continuous
//! 2D complex array using **azimuth midpoint selection**.
//!
//! # Algorithm
//!
//! Each TOPS burst overlaps its neighbour by `overlap_lines` azimuth lines.
//! Midpoint selection assigns each output line to exactly one burst:
//!
//! ```text
//! Burst 0   → local lines [0,           lps − ⌊ov/2⌋)
//! Burst k   → local lines [⌈ov/2⌉,     lps − ⌊ov/2⌋)
//! Burst N−1 → local lines [⌈ov/2⌉,     lps            )
//! ```
//!
//! where `ov = overlap_lines` for the adjacent pair and `lps = lines_per_burst`.
//!
//! The total output line count equals roughly `N × (lps − ov)`.
//!
//! # Why midpoint selection, not cos²-blending?
//!
//! Cosine blending operates on complex amplitudes in the accumulator:
//!
//! ```text
//! pixel_out = w₁·DN₁ + w₂·DN₂·e^{jΔφ}
//! ```
//!
//! The TOPSAR azimuth FM ramp introduces a non-zero `Δφ` across the overlap
//! region.  Without removing that ramp first (deramp), the complex-weighted
//! sum partially cancels and the power in the blend region is
//! `|w₁·DN₁ + w₂·DN₂·e^{jΔφ}|²` ≠ `w₁|DN₁|² + w₂|DN₂|²`.
//!
//! For incoherent backscatter (σ⁰ = |DN|² / K²) the final product is power,
//! so the residual phase in a single burst never matters.  But the
//! *blend* step introduces spurious phase-cancellation artefacts unless
//! deramp is applied first.  Midpoint selection avoids all of this: every
//! output sample comes from exactly one burst.
//!
//! # Range validity
//!
//! `first_valid_sample` and `last_valid_sample` from the annotation XML are
//! constant within each burst (verified on S1A + S1B real data).  The output
//! range window is the **intersection** (most conservative) of valid extents
//! across all bursts, guaranteeing that every output column is backed by real
//! signal in every contributing burst.

use crate::slc_reader::{BurstReader, SlcReadError};
use crate::{BurstEntry, SubSwathMetadata};
use thiserror::Error;

// ─── Error type ───────────────────────────────────────────────────────────────

/// Errors produced by deburst operations.
#[derive(Debug, Error)]
pub enum DeburstError {
    /// I/O or TIFF format error from the SLC reader.
    #[error("SLC read error: {0}")]
    SlcRead(#[from] SlcReadError),

    /// Burst list is empty.
    #[error("burst list is empty")]
    NoBursts,

    /// Burst count given to the function does not match the subswath declaration.
    #[error("subswath declares {declared} bursts but {provided} were given")]
    BurstCountMismatch { declared: usize, provided: usize },

    /// The TIFF dimensions are too small for the number of bursts × lines.
    #[error(
        "TIFF height {tiff_height} is too small for {n_bursts} bursts × \
         {lines_per_burst} lines = {needed}"
    )]
    TiffTooSmall {
        tiff_height: u32,
        n_bursts: usize,
        lines_per_burst: usize,
        needed: usize,
    },

    /// The intersection of valid range windows across all bursts is empty.
    #[error(
        "no valid range window after intersecting all burst valid-sample extents: \
         first_valid={first_valid}, last_valid={last_valid}"
    )]
    EmptyValidWindow {
        first_valid: usize,
        last_valid: usize,
    },

    /// The computed overlap between two adjacent bursts is non-positive.
    ///
    /// This indicates inconsistent burst timing or an incorrect
    /// `azimuth_time_interval_s` value.
    #[error(
        "non-positive overlap between burst {k} and burst {k1}: \
         computed overlap = {overlap_lines} lines; \
         check burst azimuth times and azimuth_time_interval_s"
    )]
    InvalidOverlap {
        k: usize,
        k1: usize,
        overlap_lines: i64,
    },

    /// The computed overlap falls outside the plausible physical range for
    /// Sentinel-1 IW TOPS.  Hitting this means either the annotation's
    /// azimuth timing is corrupt or the product is not an IW TOPS SLC.
    #[error(
        "implausible overlap between burst {k} and burst {k1}: \
         {overlap_lines} lines is outside [{min}, {max}] (S1 IW typical ≈ 150); \
         check burst timing and product mode"
    )]
    ImplausibleOverlap {
        k: usize,
        k1: usize,
        overlap_lines: usize,
        min: usize,
        max: usize,
    },
}

// ─── Output type ──────────────────────────────────────────────────────────────

/// A debursted single-subswath SLC array.
///
/// Contains raw CInt16 samples for one subswath after merging all bursts
/// by azimuth midpoint selection.  No calibration or noise removal is applied.
///
/// # Coordinate convention
///
/// - `data[line * samples + col]` = `[I, Q]`.
/// - Lines increase in azimuth; line 0 is the first contributing line of burst 0.
/// - Samples increase in slant range; sample 0 corresponds to
///   `valid_sample_offset` in the original TIFF column.  Use
///   `col + valid_sample_offset` to recover the TIFF column index.
#[derive(Debug)]
pub struct DeburstArray {
    /// Flat row-major buffer.  Length = `lines × samples`.
    pub data: Vec<[i16; 2]>,
    /// Number of azimuth lines in the output.
    pub lines: usize,
    /// Number of range samples per line (width of valid window).
    pub samples: usize,
    /// Column offset: output column 0 maps to TIFF column `valid_sample_offset`.
    pub valid_sample_offset: usize,
}

// ─── Core functions ───────────────────────────────────────────────────────────

/// Compute the overlap in azimuth lines between two adjacent bursts.
///
/// Derived from burst start times:
///
/// ```text
/// step_lines = round((t_{k+1} − t_k) / ati_s)
/// overlap    = lines_per_burst − step_lines
/// ```
///
/// # Arguments
///
/// * `burst_k`       — Earlier burst (lower `burst_index`).
/// * `burst_k1`      — Later burst (higher `burst_index`).
/// * `ati_s`         — Azimuth time interval in seconds per SLC line
///                     (from `SubSwathMetadata::azimuth_time_interval_s`).
/// * `lines_per_burst` — Lines per burst (from `SubSwathMetadata::lines_per_burst`).
///
/// # Errors
///
/// Returns [`DeburstError::InvalidOverlap`] if the computed overlap is ≤ 0,
/// which indicates an impossible timing relationship.
pub fn compute_overlap_lines(
    burst_k: &BurstEntry,
    burst_k1: &BurstEntry,
    ati_s: f64,
    lines_per_burst: usize,
) -> Result<usize, DeburstError> {
    // Sentinel-1 IW TOPS burst overlap is physically bounded: typical values
    // are 140–180 lines.  We widen the accepted band to [50, 400] to tolerate
    // IPF revisions and other TOPS modes (e.g. EW) without silently masking a
    // corrupt annotation.  Hitting the bounds is a loud failure, not a warning.
    const MIN_PLAUSIBLE_OVERLAP: usize = 50;
    const MAX_PLAUSIBLE_OVERLAP: usize = 400;

    let dt_us = (burst_k1.azimuth_time_utc - burst_k.azimuth_time_utc)
        .num_microseconds()
        .unwrap_or(0); // SAFETY-OK: chrono microseconds cannot overflow for inter-burst deltas (≈3 s)
    let dt_s = dt_us as f64 * 1e-6;
    let step_lines = (dt_s / ati_s).round() as i64;
    let overlap = lines_per_burst as i64 - step_lines;

    if overlap <= 0 {
        return Err(DeburstError::InvalidOverlap {
            k: burst_k.burst_index,
            k1: burst_k1.burst_index,
            overlap_lines: overlap,
        });
    }
    let overlap = overlap as usize;
    if overlap < MIN_PLAUSIBLE_OVERLAP || overlap > MAX_PLAUSIBLE_OVERLAP {
        return Err(DeburstError::ImplausibleOverlap {
            k: burst_k.burst_index,
            k1: burst_k1.burst_index,
            overlap_lines: overlap,
            min: MIN_PLAUSIBLE_OVERLAP,
            max: MAX_PLAUSIBLE_OVERLAP,
        });
    }
    Ok(overlap)
}

/// Deburst one subswath into a continuous complex array.
///
/// Reads raw CInt16 data from `reader`, applies midpoint selection in overlap
/// zones, crops to the valid range window, and returns the merged array.
///
/// # Arguments
///
/// * `reader` — Open [`SlcReader`] for the measurement TIFF corresponding to
///              `sw` and the chosen polarization.
/// * `sw`     — Subswath metadata (from `SceneMetadata::sub_swaths`).
/// * `bursts` — Burst entries for this subswath only
///              (`burst.subswath_id == sw.id`), sorted by `burst_index`
///              ascending.  Length must equal `sw.burst_count`.
///
/// # Output
///
/// Returns a [`DeburstArray`] with raw CInt16 samples.  No calibration,
/// noise removal, or phase correction is applied.
pub fn deburst_subswath<R: BurstReader>(
    reader: &mut R,
    sw: &SubSwathMetadata,
    bursts: &[BurstEntry],
) -> Result<DeburstArray, DeburstError> {
    // ── Precondition checks ───────────────────────────────────────────────────
    if bursts.is_empty() {
        return Err(DeburstError::NoBursts);
    }
    if bursts.len() != sw.burst_count {
        return Err(DeburstError::BurstCountMismatch {
            declared: sw.burst_count,
            provided: bursts.len(),
        });
    }

    let n = bursts.len();
    let lps = sw.lines_per_burst;
    let ati = sw.azimuth_time_interval_s;

    // ── Validate TIFF height covers all expected burst lines ──────────────────
    let needed = n * lps;
    if (reader.height() as usize) < needed {
        return Err(DeburstError::TiffTooSmall {
            tiff_height: reader.height(),
            n_bursts: n,
            lines_per_burst: lps,
            needed,
        });
    }

    // ── Valid range window: intersection across all bursts ────────────────────
    // max of firsts → most restrictive start
    // min of lasts  → most restrictive end
    let first_valid = bursts.iter().map(|b| b.first_valid_sample).max().unwrap();
    let last_valid = bursts.iter().map(|b| b.last_valid_sample).min().unwrap();
    if first_valid >= last_valid {
        return Err(DeburstError::EmptyValidWindow {
            first_valid,
            last_valid,
        });
    }
    let valid_cols = last_valid - first_valid;

    // ── Compute pairwise overlaps ─────────────────────────────────────────────
    // overlaps[k] = overlap between burst k and burst k+1, in lines.
    let mut overlaps: Vec<usize> = Vec::with_capacity(n.saturating_sub(1));
    for k in 0..n.saturating_sub(1) {
        overlaps.push(compute_overlap_lines(&bursts[k], &bursts[k + 1], ati, lps)?);
    }

    // ── Midpoint selection: per-burst [local_first, local_end) ───────────────
    //
    // For adjacent bursts k and k+1 with overlap ov:
    //   burst k   yields up to:   lps − ⌊ov/2⌋  (drops trailing ⌊ov/2⌋ lines)
    //   burst k+1 starts from:    ⌈ov/2⌉        (skips leading ⌈ov/2⌉ lines)
    //
    // This choice ensures skip_start(k+1) + keep_end_trim(k) = ov exactly,
    // so no line from the overlap is emitted twice or dropped.
    let contributions: Vec<(usize, usize)> = (0..n)
        .map(|k| {
            let local_first = if k == 0 {
                0
            } else {
                let ov = overlaps[k - 1];
                (ov + 1) / 2 // ⌈ov/2⌉
            };
            let local_end = if k == n - 1 {
                lps
            } else {
                let ov = overlaps[k];
                lps - ov / 2 // lps − ⌊ov/2⌋
            };
            (local_first, local_end)
        })
        .collect();

    // ── Allocate output ───────────────────────────────────────────────────────
    let total_lines: usize = contributions.iter().map(|(s, e)| e - s).sum();
    let mut out = vec![[0i16; 2]; total_lines * valid_cols];

    // ── Fill output ───────────────────────────────────────────────────────────
    let raw_cols = reader.width() as usize;
    let mut out_line = 0usize;

    for (k, burst) in bursts.iter().enumerate() {
        let (local_first, local_end) = contributions[k];
        let n_lines = local_end - local_first;
        if n_lines == 0 {
            continue;
        }

        // Absolute TIFF row for the first contributing line of this burst.
        // burst.first_line = burst_index * lines_per_burst (set by the parser).
        let tiff_first = burst.first_line + local_first;
        let raw = reader.read_burst_raw(tiff_first, n_lines)?;

        // Slice out [first_valid, last_valid) from each raw row.
        for line in 0..n_lines {
            let src = &raw[line * raw_cols + first_valid..line * raw_cols + last_valid];
            let dst = &mut out[(out_line + line) * valid_cols..(out_line + line + 1) * valid_cols];
            dst.copy_from_slice(src);
        }
        out_line += n_lines;
    }

    debug_assert_eq!(out_line, total_lines, "line count mismatch after filling output");

    Ok(DeburstArray {
        data: out,
        lines: total_lines,
        samples: valid_cols,
        valid_sample_offset: first_valid,
    })
}

// ─── Tests ────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use crate::parse::parse_safe_directory;
    use crate::slc_reader::SlcReader;
    use crate::{SceneMetadata, SubSwathId};

    // ── Test data paths ───────────────────────────────────────────────────────

    const S1A_SAFE: &str =
        "/home/datacube/dev/SARdine/data/SLC/\
         S1A_IW_SLC__1SDV_20201005T170824_20201005T170851_034664_04098A_1E66.SAFE";

    const S1A_IW1_VV: &str =
        "/home/datacube/dev/SARdine/data/SLC/\
         S1A_IW_SLC__1SDV_20201005T170824_20201005T170851_034664_04098A_1E66.SAFE/\
         measurement/\
         s1a-iw1-slc-vv-20201005t170824-20201005t170849-034664-04098a-004.tiff";

    fn s1a_fixtures_present() -> bool {
        std::path::Path::new(S1A_SAFE).is_dir()
    }

    fn s1a_scene() -> SceneMetadata {
        parse_safe_directory(std::path::Path::new(S1A_SAFE)).expect("failed to parse S1A SAFE")
    }

    fn iw1_bursts(scene: &SceneMetadata) -> Vec<BurstEntry> {
        scene
            .bursts
            .iter()
            .filter(|b| b.subswath_id == SubSwathId::IW1)
            .cloned()
            .collect()
    }

    fn iw1_meta(scene: &SceneMetadata) -> &crate::SubSwathMetadata {
        scene
            .sub_swaths
            .iter()
            .find(|s| s.id == SubSwathId::IW1)
            .expect("IW1 subswath not found")
    }

    // ── Unit test: overlap computation ────────────────────────────────────────

    #[test]
    fn overlap_s1a_iw1_is_in_expected_range() {
        if !s1a_fixtures_present() { eprintln!("skipping — S1A SAFE not present"); return; }
        // S1A IW TOPS bursts overlap by roughly 150–165 lines.
        // This is derived from burst timing, not a hardcoded constant.
        let scene = s1a_scene();
        let bursts = iw1_bursts(&scene);
        let sw = iw1_meta(&scene);

        assert!(bursts.len() >= 2, "need at least 2 bursts");

        let ov = compute_overlap_lines(&bursts[0], &bursts[1], sw.azimuth_time_interval_s, sw.lines_per_burst)
            .expect("overlap computation failed");

        // Empirical range for S1A IW: ~150–165 lines.
        // If this fails the azimuth_time_interval_s or burst times are wrong.
        assert!(
            ov >= 140 && ov <= 200,
            "overlap {ov} outside expected range 140–200 lines"
        );
    }

    #[test]
    fn overlap_consistent_across_all_burst_pairs() {
        if !s1a_fixtures_present() { eprintln!("skipping — S1A SAFE not present"); return; }
        // All adjacent burst pairs within a subswath should have nearly the
        // same overlap (they differ by at most 1 line due to rounding).
        let scene = s1a_scene();
        let bursts = iw1_bursts(&scene);
        let sw = iw1_meta(&scene);

        let overlaps: Vec<usize> = bursts
            .windows(2)
            .map(|w| {
                compute_overlap_lines(&w[0], &w[1], sw.azimuth_time_interval_s, sw.lines_per_burst)
                    .expect("overlap failed")
            })
            .collect();

        let min = *overlaps.iter().min().unwrap();
        let max = *overlaps.iter().max().unwrap();
        // Overlap should vary by at most 2 lines across all burst pairs
        // (rounding from burst timing quantisation).
        assert!(
            max - min <= 2,
            "overlap varies by more than 2 lines across burst pairs: {:?}",
            overlaps
        );
    }

    #[test]
    fn overlap_invalid_returns_error() {
        if !s1a_fixtures_present() { eprintln!("skipping — S1A SAFE not present"); return; }
        let scene = s1a_scene();
        let bursts = iw1_bursts(&scene);
        let sw = iw1_meta(&scene);

        // Supply lines_per_burst = 1.  The real inter-burst azimuth step is
        // ~1340 lines, so overlap = 1 − 1340 ≪ 0 → InvalidOverlap.
        let err = compute_overlap_lines(
            &bursts[0],
            &bursts[1],
            sw.azimuth_time_interval_s,
            1, // lines_per_burst = 1 is far smaller than the inter-burst step
        )
        .unwrap_err();
        assert!(
            matches!(err, DeburstError::InvalidOverlap { .. }),
            "expected InvalidOverlap, got {err}"
        );
    }

    #[test]
    fn overlap_implausibly_large_returns_error() {
        // Construct a fake burst pair with an azimuth interval short enough
        // that the computed overlap exceeds the plausible-range ceiling.
        // dt = 0.1 s, ati = 2e-3 s → step_lines = 50.  With lines_per_burst
        // = 1000 this gives overlap = 950, far above the 400-line ceiling.
        use chrono::TimeZone;
        use chrono::Utc;

        let t0 = Utc.with_ymd_and_hms(2020, 10, 5, 17, 8, 24).unwrap();
        let t1 = t0 + chrono::Duration::milliseconds(100);
        let b0 = BurstEntry {
            subswath_id: crate::types::SubSwathId::IW1,
            burst_index: 0,
            azimuth_time_utc: t0,
            first_line: 0,
            last_line: 1000,
            first_valid_sample: 0,
            last_valid_sample: 100,
            slice_index: 0,
        };
        let b1 = BurstEntry {
            subswath_id: crate::types::SubSwathId::IW1,
            burst_index: 1,
            azimuth_time_utc: t1,
            first_line: 1000,
            last_line: 2000,
            first_valid_sample: 0,
            last_valid_sample: 100,
            slice_index: 0,
        };

        let err = compute_overlap_lines(&b0, &b1, 2e-3, 1000).unwrap_err();
        assert!(
            matches!(err, DeburstError::ImplausibleOverlap { .. }),
            "expected ImplausibleOverlap, got {err}"
        );
    }

    // ── Integration test: output dimensions ───────────────────────────────────

    #[test]
    fn deburst_s1a_iw1_vv_output_shape() {
        if !s1a_fixtures_present() { eprintln!("skipping — S1A SAFE not present"); return; }
        let scene = s1a_scene();
        let bursts = iw1_bursts(&scene);
        let sw = iw1_meta(&scene);
        let mut reader = SlcReader::open(S1A_IW1_VV).expect("open failed");

        let result = deburst_subswath(&mut reader, sw, &bursts).expect("deburst failed");

        // Lines: should be N × (lps − overlap) ± small rounding.
        // For IW1 with 9 bursts, lps=1498, overlap≈158:
        //   ≈ 9 × 1340 = 12060 lines (± a few from rounding).
        // We check a generous range rather than a hardcoded value.
        assert!(
            result.lines > 11500 && result.lines < 13000,
            "debursted line count {} is outside the expected range 11500–13000",
            result.lines
        );

        // Samples: valid range window; must be positive and less than TIFF width.
        assert!(result.samples > 0, "zero samples in output");
        assert!(
            result.samples < reader.width() as usize,
            "samples {} >= TIFF width {}",
            result.samples,
            reader.width()
        );

        // Data buffer length must match declared dimensions.
        assert_eq!(
            result.data.len(),
            result.lines * result.samples,
            "data buffer length mismatch"
        );
    }

    #[test]
    fn deburst_valid_sample_offset_matches_burst_first_valid() {
        if !s1a_fixtures_present() { eprintln!("skipping — S1A SAFE not present"); return; }
        // valid_sample_offset should equal max(burst.first_valid_sample) across bursts.
        let scene = s1a_scene();
        let bursts = iw1_bursts(&scene);
        let sw = iw1_meta(&scene);
        let mut reader = SlcReader::open(S1A_IW1_VV).expect("open failed");

        let result = deburst_subswath(&mut reader, sw, &bursts).expect("deburst failed");

        let expected_offset = bursts.iter().map(|b| b.first_valid_sample).max().unwrap();
        assert_eq!(
            result.valid_sample_offset, expected_offset,
            "valid_sample_offset should equal max(burst.first_valid_sample)"
        );
    }

    // ── Integration test: seam continuity ─────────────────────────────────────

    #[test]
    fn no_all_zero_rows_at_burst_zero_one_seam() {
        if !s1a_fixtures_present() { eprintln!("skipping — S1A SAFE not present"); return; }
        // The seam between burst 0 and burst 1 must not produce all-zero rows.
        // If midpoint selection picks the wrong local line index, it can read
        // from the leading invalid guard lines of the next burst (all zero).
        let scene = s1a_scene();
        let bursts = iw1_bursts(&scene);
        let sw = iw1_meta(&scene);
        let mut reader = SlcReader::open(S1A_IW1_VV).expect("open failed");

        let result = deburst_subswath(&mut reader, sw, &bursts).expect("deburst failed");

        // Find seam: sum of lines contributed by burst 0.
        let ov0 = compute_overlap_lines(
            &bursts[0],
            &bursts[1],
            sw.azimuth_time_interval_s,
            sw.lines_per_burst,
        )
        .unwrap();
        let burst0_lines = sw.lines_per_burst - ov0 / 2; // lps − ⌊ov/2⌋
        let seam = burst0_lines; // first line from burst 1

        // Check 10 lines centred on the seam.
        let check_start = seam.saturating_sub(5);
        let check_end = (seam + 5).min(result.lines);

        for line in check_start..check_end {
            let row = &result.data[line * result.samples..(line + 1) * result.samples];
            let any_nonzero = row.iter().any(|[i, q]| *i != 0 || *q != 0);
            assert!(
                any_nonzero,
                "all-zero row at output line {line} near burst-0/burst-1 seam (seam={seam})"
            );
        }
    }

    #[test]
    fn interior_of_each_burst_region_has_signal() {
        if !s1a_fixtures_present() { eprintln!("skipping — S1A SAFE not present"); return; }
        // Check that the middle 100 lines of each burst's contribution contain
        // at least one non-zero pixel.  Guards against reading from the wrong
        // burst or inverting the local line offset.
        let scene = s1a_scene();
        let bursts = iw1_bursts(&scene);
        let sw = iw1_meta(&scene);
        let mut reader = SlcReader::open(S1A_IW1_VV).expect("open failed");

        let result = deburst_subswath(&mut reader, sw, &bursts).expect("deburst failed");

        // Rebuild contributions list (same logic as the main function).
        let overlaps: Vec<usize> = bursts
            .windows(2)
            .map(|w| {
                compute_overlap_lines(&w[0], &w[1], sw.azimuth_time_interval_s, sw.lines_per_burst)
                    .unwrap()
            })
            .collect();

        let n = bursts.len();
        let lps = sw.lines_per_burst;
        let mut out_start = 0usize;

        for k in 0..n {
            let local_first = if k == 0 { 0 } else { (overlaps[k - 1] + 1) / 2 };
            let local_end = if k == n - 1 { lps } else { lps - overlaps[k] / 2 };
            let contrib = local_end - local_first;

            // Middle 100 lines of this burst's contribution.
            let mid = out_start + contrib / 2;
            let lo = mid.saturating_sub(50);
            let hi = (mid + 50).min(out_start + contrib);

            let any_nonzero = result.data[lo * result.samples..hi * result.samples]
                .iter()
                .any(|[i, q]| *i != 0 || *q != 0);
            assert!(
                any_nonzero,
                "burst {k} interior (output lines {lo}–{hi}) is all zeros"
            );

            out_start += contrib;
        }
    }

    // ── Error path tests ──────────────────────────────────────────────────────

    #[test]
    fn deburst_error_no_bursts() {
        if !s1a_fixtures_present() { eprintln!("skipping — S1A SAFE not present"); return; }
        let scene = s1a_scene();
        let sw = iw1_meta(&scene);
        let mut reader = SlcReader::open(S1A_IW1_VV).expect("open failed");

        let err = deburst_subswath(&mut reader, sw, &[]).unwrap_err();
        assert!(matches!(err, DeburstError::NoBursts), "expected NoBursts, got {err}");
    }

    #[test]
    fn deburst_error_count_mismatch() {
        if !s1a_fixtures_present() { eprintln!("skipping — S1A SAFE not present"); return; }
        let scene = s1a_scene();
        let bursts = iw1_bursts(&scene);
        let sw = iw1_meta(&scene);
        let mut reader = SlcReader::open(S1A_IW1_VV).expect("open failed");

        // Provide only 3 bursts for a 9-burst subswath.
        let err = deburst_subswath(&mut reader, sw, &bursts[..3]).unwrap_err();
        assert!(
            matches!(err, DeburstError::BurstCountMismatch { declared: 9, provided: 3 }),
            "expected BurstCountMismatch, got {err}"
        );
    }
}
