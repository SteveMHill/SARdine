# SARdine Rebuild Progress

*Last updated: April 26, 2026*

## Current state in one sentence

End-to-end pipeline is complete and empirically validated against ASF RTC10
GAMMA at +0.016 dB median linear bias (S1B, after 10×10 multilook). **301
unit tests + 1 guard integration test** pass. The pipeline reads a full
S1A/S1B IW SLC product, applies a precise POEORB orbit, debursts,
calibrates (σ⁰ + per-pixel NESZ), merges all three IW subswaths, geocodes
via backward Range-Doppler into a selectable output CRS (EPSG:4326 or UTM
326XX/327XX via `proj4rs`), optionally applies Small (2011) terrain
flattening (γ⁰), optional speckle filtering (boxcar / Lee 1981 / Frost /
Gamma-MAP / Refined Lee), and writes a self-contained GeoTIFF. The same
library entry points are exposed as a Python extension
(`sardine-py`, PyO3 + maturin) with `process`, `grd`, `fetch_orbit`,
`download_slc`, `fetch_geoid`, and `features`.

## Quick-start for a new session

```
cd /home/datacube/dev/SARdine/new_pipeline/scene
cargo test                                  # all 301 unit + 1 guard test must pass
cargo build --release

# Inspect a SAFE
cargo run --release -- inspect /path/to/S1B.SAFE

# Full pipeline (POEORB required, or SARDINE_ALLOW_ANNOTATION_ORBIT=1)
cargo run --release --features geoid-fetch -- process \
    --safe   /path/to/S1B.SAFE \
    --orbit  /path/to/S1B_POEORB.EOF \
    --dem    /path/to/dem_tiles_dir/ \
    --output sardine_out.tiff \
    --polarization VV \
    --geoid auto \
    --crs auto                              # optional; default EPSG:4326
```

Python extension (built once, then importable):

```
cd /home/datacube/dev/SARdine/new_pipeline/sardine-py
maturin develop --release --features fetch-all
python3 -c "import sardine; print(sardine.features())"
```

Test data:
- `data/SLC/S1B_IW_SLC__1SDV_20190123T053348…SAFE` (POEORB present in `legacy/old_package_copy/SARdine/orbit_cache/`)
- `data/SLC/S1A_IW_SLC__1SDV_20201005T170824…SAFE` (no POEORB; needs `SARDINE_ALLOW_ANNOTATION_ORBIT=1` or `--features orbit-fetch`)
- DEM: `data/dem/srtm1/` (SRTM-1 `.hgt` tiles)
- Reference: `data/ASF/S1B_IW_20190123T053348_DVP_RTC10_G_gpufem_6EC9/`

See [HANDOVER.md §3](HANDOVER.md) for the full quick-start and last verified run.

---

## `sardine-scene` crate — `new_pipeline/scene/`

Rust crate. **301 unit tests + 1 guard integration test**, no GDAL dependency.

**Cargo.toml deps:** `chrono 0.4`, `quick-xml 0.36`, `serde 1`, `thiserror 1`,
`proj4rs 0.1`, `rayon 1`, `tempfile 3` (dev). Optional: `reqwest`, `zip`,
`tokio` behind `orbit-fetch` / `geoid-fetch` / `slc-fetch` features.

The Python extension lives in a sibling crate `new_pipeline/sardine-py/`
(PyO3 0.22, abi3-py38), built with `maturin develop`. It depends on
`sardine-scene` and re-exports library entry points to Python.

### Modules

| Module | Purpose | Tests | Step |
|--------|---------|-------|------|
| `types.rs` | All domain types: Mission, AcquisitionMode, Polarization, SubSwathId, OrbitData, BoundingBox, SubSwathMetadata, BurstEntry, SceneMetadata | — | A |
| `validate.rs` | 24 structural/scientific invariant checks, collected (not first-error) | 23 | A |
| `parse/xml.rs` | Serde XML structs for S1 annotation format (private) | — | B |
| `parse/mod.rs` | `parse_safe_directory()` — reads annotation XMLs, builds validated `SceneMetadata` | 11 | B |
| `orbit.rs` | EOF orbit parser + `apply_precise_orbit()` — clips, validates altitude/velocity | 11 | C |
| `calibration.rs` | `parse_calibration_noise()` — reads cal + noise XMLs, structured LUT data | 9 | D |
| `slc_reader.rs` | Pure-Rust CInt16 TIFF reader (no GDAL); per-burst `read_burst_raw()` | 13 | E |
| `deburst.rs` | TOPS midpoint-selection deburst; `deburst_subswath()` merges 9 bursts | 9 | F + G |
| `apply_calibration.rs` | σ⁰ = (\|DN\|² − N) / K²; cursor-based bilinear LUT interpolation | 15 | G |
| `merge_subswaths.rs` | Integer-offset IW1+IW2+IW3 merge; midpoint hard-cut seam | 9 | H |
| `ground_range.rs` | Ground range projection + multilooking; bilinear incidence interpolation from geolocation grid | 7 | I |
| `lib.rs` | Module declarations, re-exports | — | — |

#### Key design decisions

- **No optional fields** for scientifically required values
- **Units in field names**: `_m`, `_s`, `_hz`, `_deg` suffixes
- **Exclusive-end bounds**: `[first, last)` everywhere (matches Rust ranges)
- **Burst times as absolute UTC**, converted to orbit-relative seconds on demand via `BurstEntry::azimuth_time_rel(epoch)` — avoids epoch-coupling bug
- **Explicit failures**: no silent fallbacks, no fake data paths

## Verified assumptions (cumulative)

| Fact | Status |
|------|--------|
| Calibration LUTs are linear (not dB), sigma0 ~280–340 | ✅ Confirmed |
| `absoluteCalibrationConstant` already included in LUTs | ✅ Confirmed |
| Correct formula: `σ⁰ = (\|DN\|² − N) / K²` (noise subtraction included) | ✅ Confirmed |
| Legacy formula `(DN²/K²)×A_s` double-applies the constant | ✅ Confirmed wrong |
| LUT pixel sampling: mostly every 40 range samples, shorter last step | ✅ Confirmed |
| Noise range LUT: positive linear power | ✅ Confirmed |
| Noise azimuth LUT: multiplicative factor near 1.0 | ✅ Confirmed |
| S1A `absoluteCalibrationConstant` = 1.0 | ✅ Confirmed |
| Valid sample constancy within burst | ✅ Confirmed |
| `azimuth_time_interval` = 1/azFreq ≠ 1/PRF | ✅ Confirmed |
| `radarFrequency` and `rangeSamplingRate` identical across all annotations in a product | ✅ Confirmed |
| PRF differs per swath (IW1≠IW2≠IW3) | ✅ Confirmed |
| `burstCount × linesPerBurst == numberOfLines` | ✅ Confirmed |
| Both S1A and S1B use modern noise format (range + azimuth) | ✅ Confirmed |
| 6 annotations, 6 calibrations, 6 noises per DV product | ✅ Confirmed |
| All 3 VV subswaths have identical range pixel spacing (2.329562 m) → integer-offset merge | ✅ Confirmed |
| Merged σ⁰ mean ≈ 0.118 linear (≈ −9 dB VV) over S1A vegetated land scene | ✅ Confirmed |
| Near range time IW1 = 5.356 ms; IW1 shorter in azimuth (12233) than IW2 (12246) / IW3 (12251) | ✅ Confirmed |
| CInt16 interleaved I/Q TIFF is readable without GDAL using a pure-Rust strip reader | ✅ Confirmed |
| Geolocation grid: 210 points per subswath (10 azimuth × 21 range), incidence 30°–46° | ✅ Confirmed |
| Incidence angle increases monotonically in range within each subswath | ✅ Confirmed |
| IW1 near-range ~30°, IW3 far-range ~46° — matches expected C-band IW geometry | ✅ Confirmed |
| Ground range projection preserves mean σ⁰ (0.118 → 0.121, <3% change from resampling) | ✅ Confirmed |
| **S1B γ⁰ radiometric bias vs ASF RTC10 GAMMA** — see detailed note below | ✅ Investigated |

### Radiometric comparison vs ASF RTC10 GAMMA (S1B, VV) — investigation complete

**Scene**: S1B_IW_SLC__1SDV_20190123T053348, **ASF ref**: S1B_IW_20190123T053348_DVP_RTC10_G_gpufem_6EC9

**Apparent bias (dB domain)**: −1.28 dB (sardine lower than ASF).

**Root cause**: This is **not a calibration error**. It is a Jensen's inequality / logarithm-domain
statistical artifact arising from comparing:
- sardine: unfiltered (ENL ≈ 1, speckle fully present)
- ASF RTC10: Enhanced Lee 7×7 speckle filter applied (effective ENL ≈ 30–50)

For single-look speckle (exponential distribution), `E[log x] ≈ log(E[x]) − 2.507 dB`.
For a heavily filtered product this correction is near 0. The difference explains the
−1.28 dB apparent offset with no algorithmic change needed.

**Proof** (`scripts/filter_bias_check.py`): boxcar-filtering sardine before computing the
dB-domain mean eliminates the offset systematically:

| Sardine filter | dB-domain bias | Linear-domain bias |
|:--------------:|:--------------:|:------------------:|
| 1×1 (raw)      | −1.279 dB      | **+0.466 dB**      |
| 3×3            | −0.023 dB      | +0.466 dB          |
| 7×7            | +0.446 dB      | +0.466 dB          |
| 11×11          | +0.630 dB      | +0.466 dB          |

The **linear-domain mean bias is +0.466 dB** (sardine ≈ 11% above ASF), constant across all
filter sizes. This is the true inter-implementation calibration difference. At 7×7 (matching
ASF's filter size) the dB bias converges to +0.446 dB ≈ +0.466 dB, confirming the ENL
explanation.

**Conclusion**: The +0.466 dB linear residual is within the typical ±1 dB inter-processor
tolerance for RTC products with different DEMs (SRTM1 vs Copernicus GLO-30), different
speckle-filter order (sardine: no filter; ASF GAMMA: filter before geocoding), and different
terrain-flattening implementations. The sardine calibration and flattening formulas are correct.
No Rust code change is required for this issue.

## Bugs found and fixed (cumulative)

1. **`burst_duration_s` was actually burst cycle time** — renamed to `burst_cycle_time_s`
2. **Epoch-coupling bug** — burst times changed to `DateTime<Utc>`, orbit-relative seconds computed on demand
3. **Silent fallback on invalid bursts** — `first_nonneg` returning `None` now returns explicit `ParseError`
4. **BoundingBox anti-meridian** — documented limitation
5. **Missing burst-count cross-validation** — added check that burst entries match `SubSwathMetadata.burst_count`
6. **Missing burst monotonicity** — added per-subswath time-ordering check

## End-to-end pipeline run (verified)

The full pipeline (deburst → calibrate → merge → terrain-correct → flatten →
dB → GeoTIFF) is wired into the `process` subcommand of the `sardine` CLI
(`bin/sardine.rs`) and into the two examples
`examples/dump_merged_sigma0.rs` and `examples/dump_s1b_tc.rs`.

Last verified end-to-end run (S1B, VV, 2019-01-23, 30 threads):
- Output: `sardine_s1b_vv_30threads.tiff` (2.8 GiB Float32 EPSG:4326)
- 522 600 024 valid pixels, 2 446 flat-masked, 0 DEM-missing, 0 non-converged
- Wall-clock: ~30 min on 30 threads
- Radiometric validation vs ASF RTC10 GAMMA: **+0.016 dB median linear bias**
  after 10×10 multilook (`scripts/multilook_compare.py`)
- Seam continuity (`scripts/seam_continuity.py`): IW1/IW2 = −0.055 dB,
  IW2/IW3 = −0.129 dB (both PASS, < 0.2 dB threshold)

See [HANDOVER.md §7](HANDOVER.md) for the full validation write-up,
including the explanation of the −1.28 dB single-look dB-domain offset as a
Jensen's-inequality / ENL-mismatch artifact (not a calibration bug).

---

## Stage status (current)

All eleven processing stages from the original roadmap are complete. The
full table now lives in [architecture_overview.md §Processing stages](architecture_overview.md).
A short version:

| Stage | Status |
|-------|--------|
| Parse + validate metadata | ✅ |
| Apply precise orbit (POEORB) | ✅ |
| Parse cal + noise LUTs | ✅ |
| Per-subswath deburst (intensity, midpoint) | ✅ |
| σ⁰ calibration + per-pixel NESZ | ✅ |
| IW1+IW2+IW3 merge (integer offset, midpoint seam) | ✅ |
| DEM mosaic (SRTM-1 + GLO-30) + pre-flight coverage check | ✅ |
| Backward Range-Doppler geocoding | ✅ |
| Small (2011) terrain flattening (σ⁰ → γ⁰) | ✅ |
| Per-pixel NESZ noise-floor masking (optional) | ✅ |
| dB + GeoTIFF export | ✅ |
| Ground range projection (`sardine grd`) | ✅ |
| Speckle filter (5 kinds, on linear power) | ✅ |
| Selectable output CRS (EPSG:4326 + UTM via `proj4rs`) | ✅ |
| LIA / shadow-layover mask raster outputs (`--write-lia`, `--write-mask`) | ✅ |
| Provenance JSON sidecar | ✅ |
| Library entry points (`sardine_scene::run`) | ✅ |
| Python bindings (`sardine-py`, PyO3 + maturin) | ✅ |
| BigTIFF output (>4 GiB) | ✅ (`needs_bigtiff()` + `write_bigtiff_raw_inner()` wired into all three write paths in `export.rs`) |
| Coherence / polarimetric workflows | ❌ (intensity-only deburst by design) |

---

## Test data locations

| Data | Path |
|------|------|
| S1A SAFE | `data/SLC/S1A_IW_SLC__1SDV_20201005T170824_20201005T170851_034664_04098A_1E66.SAFE` |
| S1B SAFE | `data/SLC/S1B_IW_SLC__1SDV_20190123T053348_20190123T053415_014617_01B3D4_E833.SAFE` |
| S1A POEORB | `legacy/old_package_copy/SARdine/orbit_cache/S1A_IW_SLC__1SDV_20201005T170824_20201005T170851_034664_04098A_1E66_POEORB.EOF` |
| S1B POEORB | `legacy/old_package_copy/SARdine/orbit_cache/S1B_IW_SLC__1SDV_20190123T053348_20190123T053415_014617_01B3D4_E833_POEORB.EOF` |
| ASF reference (RTC10 GAMMA) | `data/ASF/S1B_IW_20190123T053348_DVP_RTC10_G_gpufem_6EC9/` |
| SRTM-1 tiles | `data/dem/srtm1/` |

---

## 12. Suggested Next Steps (April 2026 review)

The pipeline is functionally complete and radiometrically validated. The
following items came out of earlier critical reviews and are listed in
rough order of payoff. They are **not** invitations to start new feature
branches — see also "Do NOT build yet" at the end.

### Done since the previous revision

- **Speckle filtering** — 5 filters (boxcar, Lee 1981, Frost, Gamma-MAP,
  Refined Lee) on the linear-power buffer; wired into `sardine process`
  and `sardine grd` and recorded in `provenance.json`.
- **Selectable output CRS** — `OutputCrs` enum + `Projector` (proj4rs)
  integrated into the geocoding loop; `--crs <SPEC>` and
  `--pixel-spacing-m` CLI flags; EPSG:4326 + UTM 326XX/327XX + `auto`.
- **Python bindings** — `sardine-py` PyO3 crate built with maturin;
  exposes `process`, `grd`, `fetch_orbit`, `download_slc`, `fetch_geoid`,
  `features`. Releases the GIL via `py.allow_threads`.
- **LIA / shadow-layover mask sidecars** — `--write-lia`, `--write-mask`.
- **Default geoid** is no longer silently `Zero` — `--geoid` is required
  with `auto`, an explicit path, or `zero` (the last opts in to the
  ±80 m bias explicitly).
- **Pre-flight DEM coverage check** — `DemMosaic::covers_bbox` runs in
  `run.rs::prepare_merged_scene` before any heavy work.
- **Library extraction** — pipeline orchestration moved out of
  `bin/sardine.rs` into `run.rs` (`run_process`, `run_grd`, helpers).
  The CLI is a thin shim; the same code path drives the Python extension.
- **Radiometric regression test** (`tests/regression_s1b_20190123.rs`) —
  `#[ignore]`-gated, `--features geoid-fetch`.  Runs the full pipeline on
  the S1B 2019-01-23 VV scene and compares against the ASF RTC10-GAMMA
  reference.  Asserts **linear-domain mean bias ≤ ±1.5 dB** (filter-invariant).
  Measured: +0.34 dB (655k valid pairs, 103s runtime).  The dB-domain median
  (~−1 dB) is printed for information but not asserted — it reflects ENL
  mismatch (sardine unfiltered vs ASF 7×7 Enhanced Lee), not a calibration
  error (see `docs/PROGRESS.md §filter_bias_check`).  Also fixed
  `resolve_crs` to accept `"wgs84"` as documented (previously only
  `"EPSG:4326"` worked).  Also fixed `tiff::Decoder` to use
  `Limits::unlimited()` — default 256 MiB limit is too small for the
  ~1.7 GB sardine output TIFF at 0.0001° spacing.
  `features()`, all error paths (bogus SAFE/DEM/token → `RuntimeError`),
  and every new kwarg (`mode`, `speckle_order`, `iw`, `burst_range`). All pass.
  `terrain_correction.rs`, `slc_fetch.rs`, `orbit_fetch.rs`, and
  `geoid_fetch.rs` replaced with `tracing::info!`, `tracing::warn!`, or
  `tracing::debug!`.  `bin/sardine.rs` initialises a `tracing-subscriber`
  with `EnvFilter`; users can set `RUST_LOG=sardine_scene=debug`.
- **BigTIFF output** — `needs_bigtiff()` + `write_bigtiff_raw_inner()` wired
  into all three write paths in `export.rs`; the 4 GiB pre-flight gate
  was removed.

### Done since the previous revision (cont.)

- **COG output** (tiled + overviews) — `write_cog_with_crs` in `export.rs`
  builds a full overview pyramid via `downsample_2x` and writes a
  spec-compliant COG with all overview IFDs before the full-res IFD.
  `sardine process --cog` uses this pure-Rust writer; `sardine grd --cog`
  still delegates to `gdal_translate`.
- **`ground_range.rs`** — wired into `sardine grd` via `run::run_grd` →
  `ground_range::to_ground_range`. Not orphaned.
- **Dual-pol** — `--polarization VV+VH` / `dual` / `both` supported in the
  `process` and `grd` subcommands; the pipeline runs once per polarization
  and inserts `_VV` / `_VH` into the output filename.

### Should fix next

1. **`Stage`/`Pipeline` trait** so a third-party crate can insert a
   stage — only after a second concrete implementation exists.
2. **Geometric accuracy validation** — sub-pixel absolute geolocation has
   not been confirmed. Visual GIS overlay (QGIS at 1:50 000) or a corner
   reflector / point target comparison is needed (see HANDOVER.md §8.1).

### Validation scripts added

| Script | Reference | Status |
|--------|-----------|--------|
| `scripts/compare_asf.py` | ASF RTC10 GAMMA (EPSG:32632, 10 m) | ✅ used in regression |
| `scripts/compare_mpc.py` | MPC Sentinel-1 RTC, Catalyst Earth processor (EPSG:32632, 10 m, γ⁰ linear, nodata=−32768) | ✅ added — streams windowed COG via MPC SAS token; downloads once to `data/MPC/` cache |

### Nice to have later

- EGM2008 geoid (currently EGM96 only)
- More CRSs in `output_crs.rs` (currently EPSG:4326 + UTM 326XX/327XX only;
  `proj4rs` supports more but each needs an explicit constructor)
- `RadarImage` trait so an SLC / coherence / polarimetric branch can reuse
  the geocoding code
- Cubic / sinc resampling kernels (today bilinear only)
- Batch / multi-scene driver
- SIMD inner loops in calibration + geocoding

### Do NOT build yet

- **Coherence workflow** — current deburst is intensity-midpoint; coherent
  workflows need deramp + reramp + complex overlap blending. Not a feature
  flag.
- **Quad-pol H-α decomposition** — S1 IW is dual-pol (VV+VH or HH+HV)
  only; the math does not apply.
- **Generic plugin framework** — the legacy package tried this and is the
  reason for the rebuild. Add the trait when there is a second
  implementation, not before.

For each item, the rule from [AGENTS.md](../AGENTS.md) applies: inspect
first, summarize, list assumptions, list risks, propose smallest safe next
step, only then edit code. No silent fallbacks; every new failure mode is a
typed error variant; PR self-checklist must be filled in.
