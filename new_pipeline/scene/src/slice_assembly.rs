//! Slice assembly: merge multiple consecutive Sentinel-1 IW SLC `.SAFE`
//! products from the same orbit pass into a single [`AssembledScene`].
//!
//! # Overview
//!
//! Sentinel-1 IW SLC data is distributed as time-contiguous *slices* (also
//! called *frames*), each covering roughly 250 km along-track.  A full swath
//! segment is the concatenation of consecutive slices from the same orbit pass.
//!
//! This module implements **Phase 1** of slice assembly — metadata merging:
//!
//! 1. Parse each `.SAFE` independently via [`crate::parse::parse_safe_directory`].
//! 2. Validate that all slices come from the same acquisition (mission, mode,
//!    frequency, sub-swath geometry, burst timing continuity at seams).
//! 3. Merge validated single-slice metadata into one [`SceneMetadata`] with
//!    globally-renumbered bursts, a composite bounding box, and a merged orbit.
//!
//! The returned [`AssembledScene`] carries the merged `SceneMetadata` plus the
//! ordered list of SAFE directory paths.  Each [`BurstEntry`] in the merged
//! scene records its `slice_index` so that Phase 2 (multi-TIFF SLC reader
//! dispatch) can route reads to the correct TIFF file.
//!
//! # Correctness constraints
//!
//! - `BurstEntry::first_line` in the assembled scene is a **logical** azimuth
//!   offset across the stacked slices:
//!   `global_offset_for_slice_k + burst_index_within_slice * lines_per_burst`.
//!   The Phase 2 reader must translate this to a TIFF-local line by subtracting
//!   the cumulative line offset for the owning slice.
//!
//! - Burst timing at slice seams is validated to within
//!   [`BURST_CYCLE_TOLERANCE`] of the expected burst cycle time.  A gap outside
//!   this range indicates non-consecutive slices or a corrupt product.
//!
//! - The assembled orbit is the union of all slice annotation orbit state
//!   vectors, deduplicated within [`ORBIT_DEDUP_TOLERANCE_S`] and sorted by
//!   time.  When a precise orbit (POEORB) is applied later it replaces this
//!   annotation orbit entirely, so no slice-specific handling is needed.

use std::collections::BTreeSet;
use std::path::{Path, PathBuf};

use crate::types::{
    AcquisitionMode, BoundingBox, BurstEntry, Mission, OrbitData, Polarization, SceneMetadata,
    StateVector, SubSwathId, SubSwathMetadata,
};

// ── Constants ──────────────────────────────────────────────────────────────

/// Fraction of `burst_cycle_time_s` by which the inter-slice burst gap may
/// deviate before being rejected.
///
/// Sentinel-1 burst cycle times vary by ~0.4 ms across products; 10 % provides
/// a wide enough band without masking truly non-consecutive slices or products
/// with severe timing anomalies.
const BURST_CYCLE_TOLERANCE: f64 = 0.10;

/// Two state vectors whose timestamps differ by less than this (in seconds) are
/// considered duplicates; only the first is kept.
const ORBIT_DEDUP_TOLERANCE_S: f64 = 0.001;

/// Maximum relative difference tolerated between the same scalar physical
/// parameter in two different slices (e.g. `slant_range_time_s`).
const PARAM_REL_TOLERANCE: f64 = 1e-6;

// ── Error type ──────────────────────────────────────────────────────────────

/// Errors produced by [`assemble_slices`].
///
/// Every variant identifies exactly which slice (by 0-based index) is involved
/// and what the concrete mismatch is.  There are no silent fallbacks.
#[derive(Debug, thiserror::Error)]
pub enum SliceAssemblyError {
    #[error("no SAFE paths provided; need at least one")]
    EmptyInput,

    #[error("slice {index}: {source}")]
    Parse {
        index: usize,
        #[source]
        source: crate::parse::ParseError,
    },

    #[error("slice {index} mission {found:?} ≠ slice 0 mission {expected:?}")]
    MissionMismatch {
        index: usize,
        expected: Mission,
        found: Mission,
    },

    #[error("slice {index} acquisition mode {found:?} ≠ slice 0 mode {expected:?}")]
    ModeMismatch {
        index: usize,
        expected: AcquisitionMode,
        found: AcquisitionMode,
    },

    #[error("slice {index} subswath ID set differs from slice 0")]
    SubswathSetMismatch { index: usize },

    #[error(
        "slice {index} radar_frequency_hz={found:.3} differs from slice 0 \
         ({expected:.3}) by {delta:.3} Hz (rel tol {tol:.1e})"
    )]
    RadarFrequencyMismatch {
        index: usize,
        expected: f64,
        found: f64,
        delta: f64,
        tol: f64,
    },

    #[error(
        "slice {index} range_sampling_rate_hz={found:.3} differs from slice 0 \
         ({expected:.3}) by {delta:.3} Hz (rel tol {tol:.1e})"
    )]
    RangeSamplingRateMismatch {
        index: usize,
        expected: f64,
        found: f64,
        delta: f64,
        tol: f64,
    },

    #[error(
        "slices are not in temporal order: slice {index} start_time {start} is \
         not after slice {prev} stop_time {prev_stop}"
    )]
    TemporalOrderViolation {
        prev: usize,
        index: usize,
        prev_stop: chrono::DateTime<chrono::Utc>,
        start: chrono::DateTime<chrono::Utc>,
    },

    #[error(
        "slice {index} subswath {subswath}: lines_per_burst={found} ≠ \
         slice 0 value {expected}"
    )]
    LinesPerBurstMismatch {
        index: usize,
        subswath: SubSwathId,
        expected: usize,
        found: usize,
    },

    #[error(
        "slice {index} subswath {subswath}: range_samples={found} ≠ \
         slice 0 value {expected}"
    )]
    RangeSamplesMismatch {
        index: usize,
        subswath: SubSwathId,
        expected: usize,
        found: usize,
    },

    #[error(
        "slice {index} subswath {subswath}: {param}={found:.8} differs from \
         slice 0 ({expected:.8}) by more than rel tol {tol:.1e}"
    )]
    SubswathParamMismatch {
        index: usize,
        subswath: SubSwathId,
        param: &'static str,
        expected: f64,
        found: f64,
        tol: f64,
    },

    #[error(
        "slice {index} subswath {subswath}: burst timing gap at seam = {dt_s:.6} s, \
         expected ≈ {expected_s:.6} s (tolerance ±{tolerance_pct:.0}%)"
    )]
    BurstTimingSeamViolation {
        index: usize,
        subswath: SubSwathId,
        dt_s: f64,
        expected_s: f64,
        tolerance_pct: f64,
    },

    #[error("assembled scene failed internal validation: {0}")]
    AssembledSceneInvalid(#[from] crate::validate::ValidationErrors),
}

// ── Output type ─────────────────────────────────────────────────────────────

/// Output of [`assemble_slices`].
///
/// `scene.bursts[i].slice_index` maps each burst back to the originating SAFE
/// directory: `safe_paths[scene.bursts[i].slice_index]`.
pub struct AssembledScene {
    /// Merged, validated scene metadata covering all input slices.
    ///
    /// `burst.first_line` values are **logical** offsets across the stacked
    /// slices, not TIFF-local offsets.  See module-level docs.
    pub scene: SceneMetadata,

    /// SAFE directory paths in slice order.
    ///
    /// `safe_paths[i]` is the source SAFE for all bursts whose
    /// `slice_index == i`.
    pub safe_paths: Vec<PathBuf>,
}

// ── Public API ───────────────────────────────────────────────────────────────

/// Assemble N consecutive Sentinel-1 IW SLC slices into one [`AssembledScene`].
///
/// # Arguments
///
/// * `safe_paths` — Ordered list of `.SAFE` directory paths, ascending in
///   time (first → last along the orbit pass).  A single-element slice is
///   valid and produces an `AssembledScene` equivalent to a direct
///   `parse_safe_directory` call (all bursts get `slice_index = 0`).
///
/// # Errors
///
/// Returns [`SliceAssemblyError`] if:
///
/// - An individual SAFE fails to parse.
/// - Slices are not from the same mission / mode.
/// - Slices are not in temporal order.
/// - Sub-swath geometry (lines_per_burst, range_samples, physical params)
///   differs across slices.
/// - The burst timing gap at any slice seam is outside the expected
///   burst-cycle window (i.e. slices are not truly consecutive).
pub fn assemble_slices(safe_paths: &[&Path]) -> Result<AssembledScene, SliceAssemblyError> {
    if safe_paths.is_empty() {
        return Err(SliceAssemblyError::EmptyInput);
    }

    // 1. Parse every SAFE independently.
    let scenes: Vec<SceneMetadata> = safe_paths
        .iter()
        .enumerate()
        .map(|(i, &p)| {
            crate::parse::parse_safe_directory(p).map_err(|e| SliceAssemblyError::Parse {
                index: i,
                source: e,
            })
        })
        .collect::<Result<Vec<_>, _>>()?;

    // 2. Validate cross-slice compatibility.
    validate_compatibility(&scenes)?;

    // 3. Merge into a single SceneMetadata.
    let scene = merge_scenes(&scenes)?;

    Ok(AssembledScene {
        scene,
        safe_paths: safe_paths.iter().map(|&p| p.to_path_buf()).collect(),
    })
}

// ── Validation ───────────────────────────────────────────────────────────────

fn validate_compatibility(scenes: &[SceneMetadata]) -> Result<(), SliceAssemblyError> {
    let ref0 = &scenes[0];

    // ── Per-slice scalar checks ───────────────────────────────────────────────
    for (i, sc) in scenes.iter().enumerate().skip(1) {
        if sc.mission != ref0.mission {
            return Err(SliceAssemblyError::MissionMismatch {
                index: i,
                expected: ref0.mission,
                found: sc.mission,
            });
        }
        if sc.acquisition_mode != ref0.acquisition_mode {
            return Err(SliceAssemblyError::ModeMismatch {
                index: i,
                expected: ref0.acquisition_mode,
                found: sc.acquisition_mode,
            });
        }

        check_rel_eq(ref0.radar_frequency_hz, sc.radar_frequency_hz, PARAM_REL_TOLERANCE)
            .map_err(|(exp, fnd, dlt, tol)| SliceAssemblyError::RadarFrequencyMismatch {
                index: i,
                expected: exp,
                found: fnd,
                delta: dlt,
                tol,
            })?;

        check_rel_eq(
            ref0.range_sampling_rate_hz,
            sc.range_sampling_rate_hz,
            PARAM_REL_TOLERANCE,
        )
        .map_err(|(exp, fnd, dlt, tol)| SliceAssemblyError::RangeSamplingRateMismatch {
            index: i,
            expected: exp,
            found: fnd,
            delta: dlt,
            tol,
        })?;

        // Sub-swath ID sets must be identical.
        let ids0: BTreeSet<SubSwathId> = ref0.sub_swaths.iter().map(|s| s.id).collect();
        let ids_i: BTreeSet<SubSwathId> = sc.sub_swaths.iter().map(|s| s.id).collect();
        if ids0 != ids_i {
            return Err(SliceAssemblyError::SubswathSetMismatch { index: i });
        }

        // Per-subswath geometry.
        for sw0 in &ref0.sub_swaths {
            // Safety: existence of sw_i in sc is guaranteed by the set-equality check above.
            let sw_i = sc
                .sub_swaths
                .iter()
                .find(|s| s.id == sw0.id)
                .expect("subswath IDs validated equal above"); // SAFETY-OK: SubswathSetMismatch guard above ensures sw_i exists

            if sw_i.lines_per_burst != sw0.lines_per_burst {
                return Err(SliceAssemblyError::LinesPerBurstMismatch {
                    index: i,
                    subswath: sw0.id,
                    expected: sw0.lines_per_burst,
                    found: sw_i.lines_per_burst,
                });
            }
            if sw_i.range_samples != sw0.range_samples {
                return Err(SliceAssemblyError::RangeSamplesMismatch {
                    index: i,
                    subswath: sw0.id,
                    expected: sw0.range_samples,
                    found: sw_i.range_samples,
                });
            }

            for (name, v0, vi) in [
                ("slant_range_time_s", sw0.slant_range_time_s, sw_i.slant_range_time_s),
                (
                    "azimuth_time_interval_s",
                    sw0.azimuth_time_interval_s,
                    sw_i.azimuth_time_interval_s,
                ),
                (
                    "range_pixel_spacing_m",
                    sw0.range_pixel_spacing_m,
                    sw_i.range_pixel_spacing_m,
                ),
            ] {
                check_rel_eq(v0, vi, PARAM_REL_TOLERANCE).map_err(
                    |(exp, fnd, _dlt, tol)| SliceAssemblyError::SubswathParamMismatch {
                        index: i,
                        subswath: sw0.id,
                        param: name,
                        expected: exp,
                        found: fnd,
                        tol,
                    },
                )?;
            }
        }
    }

    // ── Temporal ordering ─────────────────────────────────────────────────────
    for (i, pair) in scenes.windows(2).enumerate() {
        let (prev, curr) = (&pair[0], &pair[1]);
        if curr.start_time <= prev.stop_time {
            return Err(SliceAssemblyError::TemporalOrderViolation {
                prev: i,
                index: i + 1,
                prev_stop: prev.stop_time,
                start: curr.start_time,
            });
        }
    }

    // ── Burst timing at seams ─────────────────────────────────────────────────
    //
    // Between the last burst of slice k and the first burst of slice k+1, the
    // azimuth-time gap must be approximately one burst cycle: the satellite
    // images continuously across slice boundaries so the burst cadence should
    // not change.  A gap outside [1-tol, 1+tol] × burst_cycle indicates
    // non-consecutive slices.
    for (k, pair) in scenes.windows(2).enumerate() {
        let (prev_scene, next_scene) = (&pair[0], &pair[1]);
        for sw in &prev_scene.sub_swaths {
            let burst_cycle = sw.burst_cycle_time_s;

            let last_prev = prev_scene
                .bursts
                .iter()
                .filter(|b| b.subswath_id == sw.id)
                .last();
            let first_next = next_scene
                .bursts
                .iter()
                .find(|b| b.subswath_id == sw.id);

            let (last_prev, first_next) = match (last_prev, first_next) {
                (Some(a), Some(b)) => (a, b),
                _ => continue, // No bursts in one of the slices; skip this subswath.
            };

            let dt_us = (first_next.azimuth_time_utc - last_prev.azimuth_time_utc)
                .num_microseconds()
                .unwrap_or(0); // SAFETY-OK: chrono microseconds cannot overflow for inter-burst deltas (≈3 s)
            let dt_s = dt_us as f64 * 1e-6;

            let lo = burst_cycle * (1.0 - BURST_CYCLE_TOLERANCE);
            let hi = burst_cycle * (1.0 + BURST_CYCLE_TOLERANCE);
            if dt_s < lo || dt_s > hi {
                return Err(SliceAssemblyError::BurstTimingSeamViolation {
                    index: k + 1,
                    subswath: sw.id,
                    dt_s,
                    expected_s: burst_cycle,
                    tolerance_pct: BURST_CYCLE_TOLERANCE * 100.0,
                });
            }
        }
    }

    Ok(())
}

// ── Merge ─────────────────────────────────────────────────────────────────────

fn merge_scenes(scenes: &[SceneMetadata]) -> Result<SceneMetadata, SliceAssemblyError> {
    let ref0 = &scenes[0];

    // Composite product_id: "first..last" for multi-slice, original for single.
    let product_id = if scenes.len() == 1 {
        ref0.product_id.clone()
    } else {
        // Use 30-char prefix of each to keep the ID manageable.
        let first = ref0.product_id.chars().take(30).collect::<String>();
        let last = scenes
            .last() // SAFETY-OK: scenes.len() >= 1 (EmptyInput guard above)
            .expect("non-empty")
            .product_id
            .chars()
            .take(30)
            .collect::<String>();
        format!("{first}..{last}")
    };

    // Polarizations: intersection (only process channels present in all slices).
    let polarizations: Vec<Polarization> = ref0
        .polarizations
        .iter()
        .copied()
        .filter(|p| scenes.iter().all(|s| s.polarizations.contains(p)))
        .collect();

    // Bounding box: union of all scene footprints.
    let bounding_box = scenes.iter().fold(
        BoundingBox {
            min_lat_deg: f64::INFINITY,
            max_lat_deg: f64::NEG_INFINITY,
            min_lon_deg: f64::INFINITY,
            max_lon_deg: f64::NEG_INFINITY,
        },
        |acc, s| BoundingBox {
            min_lat_deg: acc.min_lat_deg.min(s.bounding_box.min_lat_deg),
            max_lat_deg: acc.max_lat_deg.max(s.bounding_box.max_lat_deg),
            min_lon_deg: acc.min_lon_deg.min(s.bounding_box.min_lon_deg),
            max_lon_deg: acc.max_lon_deg.max(s.bounding_box.max_lon_deg),
        },
    );

    // Temporal window: span all slices.
    let start_time = scenes
        .iter()
        .map(|s| s.start_time)
        .min()
        .expect("non-empty"); // SAFETY-OK: scenes.len() >= 1
    let stop_time = scenes
        .iter()
        .map(|s| s.stop_time)
        .max()
        .expect("non-empty"); // SAFETY-OK: scenes.len() >= 1

    let orbit = merge_orbits(scenes);
    let sub_swaths = merge_subswaths_meta(scenes);
    let bursts = merge_bursts(scenes);

    let scene = SceneMetadata {
        product_id,
        mission: ref0.mission,
        acquisition_mode: ref0.acquisition_mode,
        polarizations,
        start_time,
        stop_time,
        radar_frequency_hz: ref0.radar_frequency_hz,
        range_sampling_rate_hz: ref0.range_sampling_rate_hz,
        bounding_box,
        sub_swaths,
        bursts,
        orbit,
    };

    scene.validated().map_err(SliceAssemblyError::AssembledSceneInvalid)
}

/// Build the merged orbit: union of all state vectors, deduplicated, sorted.
///
/// This produces the best possible annotation orbit coverage across the full
/// assembled time window.  When a precise orbit (POEORB) is applied later,
/// `apply_precise_orbit` replaces the entire orbit, so this merge only matters
/// if the caller uses annotation orbit mode.
fn merge_orbits(scenes: &[SceneMetadata]) -> OrbitData {
    let mut all_svs: Vec<StateVector> = scenes
        .iter()
        .flat_map(|s| s.orbit.state_vectors.iter().cloned())
        .collect();

    all_svs.sort_by_key(|sv| sv.time);

    // Deduplicate: keep first of each group whose timestamps are within
    // ORBIT_DEDUP_TOLERANCE_S of each other.
    let mut deduped: Vec<StateVector> = Vec::with_capacity(all_svs.len());
    for sv in all_svs {
        let is_dup = deduped.last().map(|last| {
            let dt_us = (sv.time - last.time)
                .num_microseconds()
                .unwrap_or(i64::MAX); // SAFETY-OK: chrono microseconds for typical orbit SV spacing (≤ hours) cannot overflow
            (dt_us as f64).abs() * 1e-6 < ORBIT_DEDUP_TOLERANCE_S
        });
        if is_dup != Some(true) {
            deduped.push(sv);
        }
    }

    let reference_epoch = deduped
        .first()
        .map(|sv| sv.time)
        .unwrap_or(scenes[0].orbit.reference_epoch); // SAFETY-OK: scenes.len() >= 1; fallback only if no SVs survived dedup (cannot happen with valid input)

    OrbitData {
        reference_epoch,
        state_vectors: deduped,
    }
}

/// Build one `SubSwathMetadata` per subswath spanning all input slices.
///
/// Physical parameters (pixel spacing, timing) are taken from slice 0;
/// geometric extent (`azimuth_samples`, `burst_count`) is summed across all.
fn merge_subswaths_meta(scenes: &[SceneMetadata]) -> Vec<SubSwathMetadata> {
    scenes[0]
        .sub_swaths
        .iter()
        .map(|sw0| {
            let (total_azimuth, total_bursts) =
                scenes.iter().fold((0usize, 0usize), |(az, bc), s| {
                    let sw = s
                        .sub_swaths
                        .iter()
                        .find(|s| s.id == sw0.id)
                        .expect("subswath IDs validated equal"); // SAFETY-OK: SubswathSetMismatch guard ensures all scenes have the same subswaths
                    (az + sw.azimuth_samples, bc + sw.burst_count)
                });

            SubSwathMetadata {
                id: sw0.id,
                burst_count: total_bursts,
                lines_per_burst: sw0.lines_per_burst,
                range_samples: sw0.range_samples,
                azimuth_samples: total_azimuth,
                first_line: 0,
                last_line: total_azimuth,
                first_sample: sw0.first_sample,
                last_sample: sw0.last_sample,
                range_pixel_spacing_m: sw0.range_pixel_spacing_m,
                azimuth_pixel_spacing_m: sw0.azimuth_pixel_spacing_m,
                slant_range_time_s: sw0.slant_range_time_s,
                azimuth_time_interval_s: sw0.azimuth_time_interval_s,
                prf_hz: sw0.prf_hz,
                burst_cycle_time_s: sw0.burst_cycle_time_s,
            }
        })
        .collect()
}

/// Build the globally-renumbered burst list.
///
/// For each subswath, iterates slices in order and assigns:
///
/// - `burst_index`: monotonically increasing across slices
/// - `first_line` / `last_line`: logical offset = cumulative azimuth lines from
///   earlier slices + local offset within the current slice's TIFF
/// - `slice_index`: the owning slice (for Phase 2 multi-TIFF dispatch)
fn merge_bursts(scenes: &[SceneMetadata]) -> Vec<BurstEntry> {
    let subswath_ids: Vec<SubSwathId> = scenes[0].sub_swaths.iter().map(|s| s.id).collect();
    let mut all_bursts: Vec<BurstEntry> = Vec::new();

    for &sw_id in &subswath_ids {
        let mut global_burst_index: usize = 0;
        let mut global_line_offset: usize = 0;

        for (slice_idx, scene) in scenes.iter().enumerate() {
            let sw = scene
                .sub_swaths
                .iter()
                .find(|s| s.id == sw_id)
                .expect("subswath IDs validated equal"); // SAFETY-OK: same as merge_subswaths_meta

            let slice_bursts: Vec<&BurstEntry> = scene
                .bursts
                .iter()
                .filter(|b| b.subswath_id == sw_id)
                .collect();

            for burst in &slice_bursts {
                // burst.first_line within the slice's TIFF equals
                // burst.burst_index * lines_per_burst (set by the parser).
                let local_first = burst.burst_index * sw.lines_per_burst;

                all_bursts.push(BurstEntry {
                    subswath_id: sw_id,
                    burst_index: global_burst_index,
                    azimuth_time_utc: burst.azimuth_time_utc,
                    first_line: global_line_offset + local_first,
                    last_line: global_line_offset + local_first + sw.lines_per_burst,
                    first_valid_sample: burst.first_valid_sample,
                    last_valid_sample: burst.last_valid_sample,
                    slice_index: slice_idx,
                });
                global_burst_index += 1;
            }

            global_line_offset += sw.azimuth_samples;
        }
    }

    // Sort by (subswath, burst_index): matches the convention in parse/mod.rs.
    all_bursts.sort_by(|a, b| {
        a.subswath_id
            .cmp(&b.subswath_id)
            .then(a.burst_index.cmp(&b.burst_index))
    });

    all_bursts
}

// ── Helpers ──────────────────────────────────────────────────────────────────

/// Return `Ok(())` when `|found - expected| / |expected| ≤ tol`.
/// Return `Err((expected, found, abs_delta, tol))` otherwise.
fn check_rel_eq(expected: f64, found: f64, tol: f64) -> Result<(), (f64, f64, f64, f64)> {
    let delta = (found - expected).abs();
    let rel = if expected.abs() > 0.0 {
        delta / expected.abs()
    } else {
        delta
    };
    if rel <= tol {
        Ok(())
    } else {
        Err((expected, found, delta, tol))
    }
}

// ── Unit tests ────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use chrono::{TimeZone, Utc};

    // ── Synthetic metadata builders ───────────────────────────────────────────

    fn make_orbit(epoch: chrono::DateTime<Utc>) -> crate::types::OrbitData {
        // 17 state vectors at 10-second intervals — passes the MIN_ORBIT_VECTORS check.
        let state_vectors = (0..17u32)
            .map(|i| crate::types::StateVector {
                time: epoch + chrono::Duration::seconds(i64::from(i) * 10),
                position_m: [4_000_000.0 + f64::from(i) * 1000.0, 5_000_000.0, 3_000_000.0],
                velocity_m_s: [100.0, 200.0, 7500.0],
            })
            .collect();
        crate::types::OrbitData {
            reference_epoch: epoch,
            state_vectors,
        }
    }

    /// Build a single-subswath (IW1 only) `SceneMetadata` with `n_bursts`
    /// bursts whose first burst starts at `t0_ms` milliseconds after the fixed
    /// base epoch (2020-10-05 17:08:24 UTC).
    ///
    /// `burst_cycle_ms` is the inter-burst interval in milliseconds (default
    /// realistic S-1 value: 2758 ms).
    fn make_scene(n_bursts: usize, t0_ms: i64, burst_cycle_ms: i64) -> SceneMetadata {
        const LINES_PER_BURST: usize = 1526;
        const RANGE_SAMPLES: usize = 21608;

        let base = Utc.with_ymd_and_hms(2020, 10, 5, 17, 8, 24).unwrap();
        let orbit_epoch = Utc.with_ymd_and_hms(2020, 10, 5, 17, 6, 0).unwrap();

        let bursts: Vec<BurstEntry> = (0..n_bursts)
            .map(|i| BurstEntry {
                subswath_id: SubSwathId::IW1,
                burst_index: i,
                azimuth_time_utc: base
                    + chrono::Duration::milliseconds(t0_ms + i as i64 * burst_cycle_ms),
                first_line: i * LINES_PER_BURST,
                last_line: (i + 1) * LINES_PER_BURST,
                first_valid_sample: 100,
                last_valid_sample: 21600,
                slice_index: 0,
            })
            .collect();

        let azimuth_samples = n_bursts * LINES_PER_BURST;
        let start_time = bursts[0].azimuth_time_utc;
        let stop_time = *bursts
            .last() // SAFETY-OK: n_bursts > 0 (tests always pass ≥ 1)
            .map(|b| &b.azimuth_time_utc)
            .expect("n_bursts > 0 in test helper");

        SceneMetadata {
            product_id: format!("SYNTHETIC_T{t0_ms}_N{n_bursts}"),
            mission: crate::types::Mission::S1A,
            acquisition_mode: crate::types::AcquisitionMode::IW,
            polarizations: vec![crate::types::Polarization::VV],
            start_time,
            stop_time,
            radar_frequency_hz: 5.405_000_454_334_349e9,
            range_sampling_rate_hz: 64.345_238_095_238_1e6,
            bounding_box: crate::types::BoundingBox {
                min_lat_deg: 48.0 + t0_ms as f64 * 1e-6,
                max_lat_deg: 50.0 + t0_ms as f64 * 1e-6,
                min_lon_deg: 11.0,
                max_lon_deg: 13.0,
            },
            sub_swaths: vec![crate::types::SubSwathMetadata {
                id: SubSwathId::IW1,
                burst_count: n_bursts,
                lines_per_burst: LINES_PER_BURST,
                range_samples: RANGE_SAMPLES,
                azimuth_samples,
                first_line: 0,
                last_line: azimuth_samples,
                first_sample: 0,
                last_sample: RANGE_SAMPLES,
                range_pixel_spacing_m: 2.329_562,
                azimuth_pixel_spacing_m: 13.968,
                slant_range_time_s: 0.005_387,
                azimuth_time_interval_s: 2.055_556e-3,
                prf_hz: 1717.0,
                burst_cycle_time_s: burst_cycle_ms as f64 * 1e-3,
            }],
            bursts,
            orbit: make_orbit(orbit_epoch),
        }
    }

    // ── Compatibility check tests ─────────────────────────────────────────────

    /// Two identical scenes pass compatibility.
    #[test]
    fn validate_two_identical_scenes_passes() {
        let s0 = make_scene(9, 0, 2758);
        let s1 = make_scene(9, 9 * 2758, 2758); // starts right after s0
        validate_compatibility(&[s0, s1]).unwrap();
    }

    #[test]
    fn validate_mission_mismatch_rejected() {
        let s0 = make_scene(9, 0, 2758);
        let mut s1 = make_scene(9, 9 * 2758, 2758);
        s1.mission = crate::types::Mission::S1B;
        let err = validate_compatibility(&[s0, s1]).unwrap_err();
        assert!(
            matches!(err, SliceAssemblyError::MissionMismatch { index: 1, .. }),
            "expected MissionMismatch, got: {err}"
        );
    }

    #[test]
    fn validate_mode_mismatch_rejected() {
        let s0 = make_scene(9, 0, 2758);
        let mut s1 = make_scene(9, 9 * 2758, 2758);
        // There is no EW subswath in the synthetic scene so just change the mode field.
        s1.acquisition_mode = crate::types::AcquisitionMode::EW;
        let err = validate_compatibility(&[s0, s1]).unwrap_err();
        assert!(
            matches!(err, SliceAssemblyError::ModeMismatch { index: 1, .. }),
            "expected ModeMismatch, got: {err}"
        );
    }

    #[test]
    fn validate_temporal_order_violation_rejected() {
        // s0 starts later than s1 → wrong order.
        let s0 = make_scene(9, 9 * 2758, 2758);
        let s1 = make_scene(9, 0, 2758); // starts before s0 stop
        let err = validate_compatibility(&[s0, s1]).unwrap_err();
        assert!(
            matches!(err, SliceAssemblyError::TemporalOrderViolation { .. }),
            "expected TemporalOrderViolation, got: {err}"
        );
    }

    #[test]
    fn validate_lines_per_burst_mismatch_rejected() {
        let s0 = make_scene(9, 0, 2758);
        let mut s1 = make_scene(9, 9 * 2758, 2758);
        s1.sub_swaths[0].lines_per_burst = 1500; // differs from s0's 1526
        let err = validate_compatibility(&[s0, s1]).unwrap_err();
        assert!(
            matches!(err, SliceAssemblyError::LinesPerBurstMismatch { index: 1, .. }),
            "expected LinesPerBurstMismatch, got: {err}"
        );
    }

    #[test]
    fn validate_burst_timing_seam_too_large_rejected() {
        let s0 = make_scene(9, 0, 2758);
        // s1 starts 5 × burst_cycle after s0 ends → gap = 5 * 2758 ms, way outside tolerance.
        let s1 = make_scene(9, 50 * 2758, 2758);
        let err = validate_compatibility(&[s0, s1]).unwrap_err();
        assert!(
            matches!(err, SliceAssemblyError::BurstTimingSeamViolation { .. }),
            "expected BurstTimingSeamViolation, got: {err}"
        );
    }

    // ── Merge correctness tests ───────────────────────────────────────────────

    /// Single-slice assembly is a pass-through: burst count and indices unchanged.
    #[test]
    fn merge_single_slice_is_passthrough() {
        let s0 = make_scene(9, 0, 2758);
        let scenes = [s0];
        let merged = merge_scenes(&scenes).unwrap();
        assert_eq!(merged.bursts.len(), 9);
        assert_eq!(merged.sub_swaths[0].burst_count, 9);
        // All bursts have slice_index 0.
        assert!(merged.bursts.iter().all(|b| b.slice_index == 0));
        // first_line of burst i = i * lines_per_burst (unchanged from single-SAFE).
        for (i, b) in merged.bursts.iter().enumerate() {
            assert_eq!(b.first_line, i * 1526, "burst {i} first_line");
            assert_eq!(b.burst_index, i, "burst {i} burst_index");
        }
    }

    /// Two-slice assembly: burst count doubles, global line offsets are correct,
    /// slice_index is assigned correctly.
    #[test]
    fn merge_two_slices_burst_count_and_offsets() {
        const N: usize = 9;
        const LPB: usize = 1526;
        const CYCLE_MS: i64 = 2758;

        let s0 = make_scene(N, 0, CYCLE_MS);
        let s1 = make_scene(N, N as i64 * CYCLE_MS, CYCLE_MS);
        let merged = merge_scenes(&[s0, s1]).unwrap();

        // Burst count.
        assert_eq!(merged.bursts.len(), 2 * N);
        assert_eq!(merged.sub_swaths[0].burst_count, 2 * N);
        assert_eq!(merged.sub_swaths[0].azimuth_samples, 2 * N * LPB);

        // First N bursts come from slice 0; next N from slice 1.
        for (i, b) in merged.bursts.iter().enumerate() {
            assert_eq!(b.burst_index, i, "burst_index at position {i}");
            if i < N {
                assert_eq!(b.slice_index, 0, "slice_index at {i}");
                assert_eq!(b.first_line, i * LPB, "first_line at {i}");
            } else {
                assert_eq!(b.slice_index, 1, "slice_index at {i}");
                // Slice 1 line offset = N * LPB (azimuth_samples of slice 0).
                let local_idx = i - N;
                let expected_first = N * LPB + local_idx * LPB;
                assert_eq!(b.first_line, expected_first, "first_line at {i}");
            }
        }
    }

    /// Bounding box of assembled scene is the union of input bounding boxes.
    #[test]
    fn merge_bounding_box_is_union() {
        let s0 = make_scene(9, 0, 2758); // lat [48.0, 50.0]
        let s1 = make_scene(9, 9 * 2758, 2758); // lat [48 + ε, 50 + ε] (ε small)
        let merged = merge_scenes(&[s0.clone(), s1]).unwrap();
        assert!(merged.bounding_box.min_lat_deg <= s0.bounding_box.min_lat_deg + 1e-9);
        assert!(merged.bounding_box.max_lat_deg >= s0.bounding_box.max_lat_deg - 1e-9);
    }

    /// Merged scene passes SceneMetadata::validated() (no internal inconsistencies).
    #[test]
    fn merge_two_slices_passes_validation() {
        let s0 = make_scene(9, 0, 2758);
        let s1 = make_scene(9, 9 * 2758, 2758);
        // merge_scenes already calls validated() internally; just check no error.
        merge_scenes(&[s0, s1]).unwrap();
    }

    /// Empty input to validate_compatibility is a no-op (single-element is trivially valid).
    #[test]
    fn validate_single_scene_passes() {
        let s0 = make_scene(9, 0, 2758);
        validate_compatibility(&[s0]).unwrap();
    }

    /// Orbit union: merged orbit has at least as many vectors as the largest input orbit.
    #[test]
    fn merge_orbit_has_at_least_as_many_vectors_as_largest() {
        let s0 = make_scene(9, 0, 2758);
        let s1 = make_scene(9, 9 * 2758, 2758);
        let max_sv = s0.orbit.state_vectors.len().max(s1.orbit.state_vectors.len());
        let merged = merge_scenes(&[s0, s1]).unwrap();
        assert!(merged.orbit.state_vectors.len() >= max_sv);
    }
}
