# TOPSAR Merge Refactor Plan

Target structure (matches your proposal):
```
src/core/topsar/
  mod.rs
  subswath_model.rs
  alignment.rs
  merge.rs
  overlap_weights.rs
  steering.rs
  diagnostics.rs
  tests/
    alignment_tests.rs
    merge_tests.rs
```

## Mapping from current code (single `topsar_merge.rs` ~3.9k LOC)
- **Types/metadata**: `SubswathAlignment`, `OutputGrid`, `OverlapRegion/Quality`, `MergePlan`, `MergeRowSegment`, `MergeWeight`, `MergeParameters`, `QualityControl`, `MergedSwathData`, `QualityResults`, `ProcessingMetadata`, `PerformanceMetrics`, `AzimuthTimingModel`, `BurstTiming`, `DopplerPolynomial` → move to `subsasth_model.rs` (SubSwath + MergeGrid + timing) and `alignment.rs` (alignment types) and `merge.rs` (MergePlan/segments/weights). Keep `MergeGrid` as new struct to replace ad-hoc fields.
- **Overlap weighting**: routines that build per-overlap weights (`compute_overlap_gains`, tapering/feathering constants) → `overlap_weights.rs`.
- **Alignment**: functions turning subswath metadata/DC polynomials into azimuth/range offsets and a common merge grid; also `has_required_subswaths` check → `alignment.rs`.
- **Merge orchestration**: main merge paths (intensity, complex, hitcount/uncovered masks, seam filling) → `merge.rs`. Keep public entry `merge_iw_subswaths` here and re-export from `mod.rs`.
- **Steering**: azimuth steering angle/rate helpers, phase adjustments tied to steering and Doppler → `steering.rs`.
- **Diagnostics**: radiometric/phase checks in overlaps, debug profiles, warning aggregation → `diagnostics.rs`.
- **Tests**: lift overlap/merge-specific unit tests from bottom of `topsar_merge.rs` into `tests/alignment_tests.rs` and `tests/merge_tests.rs` (preserve fixtures and expected values).

## Migration steps (minimize churn)
1) **Scaffold modules**: create `src/core/topsar/mod.rs` with `pub mod subswath_model; alignment; merge; overlap_weights; steering; diagnostics;` and `pub use merge::merge_iw_subswaths;`. Add `pub use` re-exports for key types to keep external API stable.
2) **Move types first**: relocate type definitions into `subsasth_model.rs` and `merge.rs` (for plan/weights). Update `topsar_merge.rs` to `pub use crate::core::topsar::*;` so callers continue compiling. (We already started a `topsar_merge_types.rs`; replace it with the new structure.)
3) **Extract pure helpers**: move overlap weighting and alignment helpers into their modules with `pub(crate)` visibility. Wire `merge.rs` to import them. Keep function signatures identical initially.
4) **Split merge logic**: move the main merge functions into `merge.rs`; leave a thin wrapper in the old file to delegate (or deprecate/remove the old file once all call sites point to `core::topsar`).
5) **Tests**: move/adjust existing inline tests into `src/core/topsar/tests/*.rs`; ensure `#[cfg(test)] mod tests;` in `topsar/mod.rs`.
6) **Clean up legacy**: remove `topsar_merge.rs` after consumers import from `core::topsar`; update `core/mod.rs` to expose `pub mod topsar;` and keep a compatibility `pub use topsar::merge_iw_subswaths;`.

## Notes / gotchas
- **SubSwath source**: currently `types::SubSwath`. Decide whether to redefine a slimmer `topsar::SubSwath` and convert in `io::annotation::derive`, or keep the existing type but centralize construction helpers in `subsasth_model.rs`.
- **Constants**: move IW defaults/thresholds into `mod.rs` or `alignment.rs` (public constants) to avoid scattering magic numbers.
- **DC/FM providers**: keep `DcFmRateProvider` trait imported in `subsasth_model.rs` to avoid circular deps. Alignment/steering functions should take providers as params rather than reach into global state.
- **Parallelism**: keep rayon usage in `merge.rs`; helper modules should stay `std`-only.
- **Diagnostics outputs**: consider optional logging hooks (feature-flag?) so core algorithms stay clean.

## Immediate next action (suggested)
- Create `src/core/topsar/mod.rs` + empty module files with re-exports.
- Move the type definitions from `topsar_merge.rs` into `subsasth_model.rs` and `merge.rs` (for plan/weights). Keep `merge_iw_subswaths` as the single public entry in `merge.rs` and re-export it.
- Add a temporary `pub use topsar::*;` in `core/mod.rs` to keep the API stable while we migrate the implementation.
