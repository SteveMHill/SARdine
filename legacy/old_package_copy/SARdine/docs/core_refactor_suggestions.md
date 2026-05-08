# Core Refactor Suggestions

Generated on 2025-12-10 after a quick pass over `src/core`. Focus: break down oversized modules, surface redundancy, and sketch a target layout that keeps calibration/deburst clean while shrinking ~47k lines in core.

## Quick Size Scan (top offenders)
- `terrain_correction.rs` (~9k LOC)
- `topsar_merge.rs` (~4k)
- `deburst/iw_deburst.rs` (~2.4k) and `deburst/ew_deburst.rs` (~2.0k) plus shared helpers
- `speckle_filter.rs` (~1.9k)
- `terrain_flatten.rs` (~1.7k) vs `scientific_terrain_flatten.rs` (~1.1k)
- `calibrate/parsing.rs` (~1.6k) and `calibrate/processing.rs` (~1.0k) (legacy vs new `calibration/`)
- `multilook.rs` (~1.4k)
- `advanced_masking.rs` (~1.3k)
- `deburst_calibrate_fused.rs` (~1.25k)

## Suggested Modularization Targets
- **Terrain correction (`terrain_correction.rs` ~9k)**
  - Split into `terrain_correction/{io.rs, grids.rs, geocoding.rs, radiometry.rs, stitching.rs, qc.rs}`.
  - Isolate DEM handling and geoid lookup into `io.rs`; move resampling/interpolation kernels to `grids.rs`; keep orthorectification math in `geocoding.rs`; stripe equalization / power preservation in `radiometry.rs`; stitching/mosaicking into `stitching.rs`; assertions + metrics into `qc.rs`.
- **TOPSAR merge (`topsar_merge.rs` ~4k)**
  - Extract burst alignment and overlap weighting into `topsar_merge/alignment.rs`.
  - Move spectral weighting kernels into `topsar_merge/kernels.rs`.
  - Keep high-level driver/orchestration in `mod.rs`; share helpers with `topsar_merge_optimized.rs` to avoid drift.
- **Deburst**
  - Directory already exists; factor shared code between `iw_deburst.rs` and `ew_deburst.rs` into `deburst/common/{timing.rs, doppler.rs, interpolation.rs, masking.rs}`.
  - Move mask propagation (currently duplicated) into `deburst/mask.rs` and have both IW/EW call it.
  - Consider keeping burst parsing/IO in `deburst/io.rs` and math kernels in `deburst/kernels.rs`.
- **Calibration**
  - Treat `calibrate/` as legacy shim only; move remaining parsing/processing logic into `calibration/{io_xml.rs, apply.rs, noise.rs, antenna.rs, lut.rs, model.rs}` (already present) and deprecate direct `calibrate::*` uses.
  - Add a small `calibration/bridge.rs` if any legacy call sites need compatibility wrappers, instead of duplicating conversions.
- **Terrain flattening**
  - Merge or clearly separate `terrain_flatten.rs` (legacy) vs `scientific_terrain_flatten.rs` (modern). If both needed, isolate shared math in `terrain_flatten/common.rs` and keep two thin strategy files.
- **Speckle filter (`speckle_filter.rs` ~1.9k)**
  - Split into `speckle_filter/{lee.rs,frost.rs,filters.rs}` and keep orchestrator in `mod.rs`.
- **Validation and QA**
  - `validation_gates.rs` and `quality_assessment.rs` can move per-check logic into `validation/{power.rs,coverage.rs,noise.rs}` with `mod.rs` coordinating.
- **Multilook & masking**
  - Break `multilook.rs` into `multilook/{kernels.rs,weights.rs,io.rs}`.
  - `advanced_masking.rs` could expose per-mask modules: `water.rs`, `terrain.rs`, `speckle.rs`.

## Redundant / Deprecated / Duplicated Areas Noted
- **Calibration duplication**: Legacy `calibrate/parsing.rs` and `calibrate/processing.rs` overlap with the newer `calibration/` module. Recommend marking `calibrate` as deprecated, keeping only a compatibility wrapper, and migrating call sites to `calibration::*`.
- **Terrain flattening**: Two implementations (`terrain_flatten.rs` vs `scientific_terrain_flatten.rs`) without clear separation; either consolidate or document strategy differences and share math helpers.
- **Deburst mask propagation**: Similar logic appears in `deburst/mask_propagation.rs` and in main deburst files; consolidate in one helper module.
- **TOPSAR merge duplication**: `topsar_merge.rs` and `topsar_merge_optimized.rs` share weighting/alignment code; extract shared kernels to avoid divergence.
- **Deprecated flags/comments**: Several `#[deprecated]` helpers remain in deburst timing code; can be moved to a `legacy.rs` or removed after confirming no call sites.

## Suggested Directory Shape (illustrative)
```
src/core/
  calibration/ (new home; shim in calibrate/ only)
  deburst/
    common/{timing.rs,doppler.rs,mask.rs,interp.rs}
    iw.rs
    ew.rs
    io.rs
    kernels.rs
  terrain_correction/{io.rs,grids.rs,geocode.rs,radiometry.rs,stitch.rs,qc.rs}
  topsar_merge/{alignment.rs,kernels.rs,driver.rs}
  terrain_flatten/{legacy.rs,scientific.rs,common.rs}
  speckle_filter/{lee.rs,frost.rs,kernels.rs,mod.rs}
  validation/{power.rs,coverage.rs,noise.rs,mod.rs}
  multilook/{kernels.rs,weights.rs,io.rs,mod.rs}
```

## Process Notes
- Prioritize carving out pure helpers (math kernels, interpolation, weighting) first; they are low-risk moves that shrink monoliths quickly.
- Keep public APIs stable via re-exports in `mod.rs` to avoid downstream churn during the split.
- Add smoke tests around extracted kernels before relocating large orchestration code.
