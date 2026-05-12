# SARdine Architecture Overview

## Status

*Last updated: April 26, 2026*

Sentinel-1 IW SLC backscatter processing pipeline. The active codebase is
two Rust crates: [`sardine/`](../sardine/)
(`sardine` v0.1.0, library + CLI) and
[`sardine-py/`](../sardine-py/) (`sardine-py`
v0.1.0, PyO3 extension built with maturin). The legacy Python+Rust package
under `legacy/` is **reference-only** — see
[LEGACY_STATUS.md](../LEGACY_STATUS.md) and [AGENTS.md](../AGENTS.md).

**Test status:** `cargo test` (scene crate) → **220 unit tests + 1 guard integration test**, all passing.
End-to-end pipeline runs and is empirically validated against ASF RTC10 GAMMA
at **+0.016 dB median linear bias** after 10×10 multilook (S1B,
2019-01-23). See [HANDOVER.md §7](HANDOVER.md) for validation details.

A mandatory forbidden-pattern guard (`scripts/check_no_silent_fallbacks.sh`,
enforced by `cargo test`) prevents silent fallbacks and hardcoded Sentinel-1
constants from re-entering the codebase. Every legitimate exception carries a
same-line `// SAFETY-OK: <reason>` annotation.

## Core principle
Legacy code is reference material, not truth. See [AGENTS.md](../AGENTS.md).

## What this package currently is

A **Rust library + CLI + Python extension** that reads a Sentinel-1 IW SLC
`.SAFE` product, a POEORB precise orbit file, and SRTM-1 / Copernicus GLO-30
DEM tiles, and writes a Float32 GeoTIFF of σ⁰ or terrain-flattened γ⁰ in dB.
Output CRS is selectable via `--crs` (default EPSG:4326; UTM zones supported
via `proj4rs` — `EPSG:326{01..60}` north, `EPSG:327{01..60}` south, or
`auto`/`utm` to pick the centroid zone). One CLI binary (`sardine`) with
subcommands `process`, `grd`, `inspect`, `fetch-orbit` (feature-gated),
`download-slc` (feature-gated). The `grd` subcommand produces a
non-georeferenced ground-range σ⁰ TIFF in radar geometry (no DEM, no
flattening) — radar geometry cannot truthfully advertise EPSG:4326 so the
GeoTIFF tags are deliberately omitted.

The Python extension (`sardine-py`, built via `maturin develop`) re-exports
the library entry points as `sardine.process`, `sardine.grd`,
`sardine.fetch_orbit`, `sardine.download_slc`, `sardine.fetch_geoid`, and
`sardine.features()`. All functions release the GIL via `py.allow_threads`
and surface failures as `RuntimeError`.

Alongside every output raster the CLI writes a `.provenance.json` sidecar
recording inputs, parameters, and outcome counts. Suppress with
`--no-provenance`. The schema is polymorphic on `mode` (`tc` or `grd`):
TC-only fields and GRD-only fields are mutually exclusive.

Not yet covered: coherence and polarimetric workflows (deburst is
intensity-only by design), multi-scene / pair / stack drivers, BigTIFF
output (>4 GiB classic-TIFF cap is enforced pre-flight), and CRSs other
than EPSG:4326 / UTM.

## Processing stages (verified against `bin/sardine.rs`)

| # | Stage | Status | Module |
|---|-------|--------|--------|
| 1 | Parse SAFE metadata | ✅ | `parse/mod.rs` |
| 2 | Apply precise orbit (POEORB; refuses extrapolation) | ✅ | `orbit.rs` |
| 3 | Parse calibration + noise LUTs | ✅ | `calibration.rs` |
| 4 | Parse geolocation grids | ✅ | `parse/mod.rs` |
| 5 | Per-subswath SLC read + TOPS midpoint deburst | ✅ (parallel) | `slc_reader.rs`, `deburst.rs` |
| 6 | Radiometric calibration σ⁰ = (\|DN\|² − N) / K² + per-pixel NESZ | ✅ | `apply_calibration.rs` |
| 7 | Merge IW1+IW2+IW3 (integer-offset, midpoint hard-cut seam) | ✅ | `merge_subswaths.rs` |
| 8 | DEM mosaic (SRTM-1 `.hgt` or GLO-30 GeoTIFF) | ✅ | `dem.rs` |
| 9 | Backward Range-Doppler geocoding (Newton zero-Doppler, parallel rows) | ✅ | `terrain_correction.rs` |
| 10 | Terrain flattening σ⁰ → γ⁰ (Small 2011 area-projection) | ✅ | `terrain_correction.rs` |
| 11 | Per-pixel NESZ noise-floor masking (optional) | ✅ | `terrain_correction.rs` |
| 12 | dB conversion + GeoTIFF export (classic TIFF ≤ 4 GiB, pre-flight size check) | ✅ | `export.rs` |
| 13 | Ground-range projection + multilook (slant→ground), no-CRS TIFF writer | ✅ (CLI: `sardine grd`) | `ground_range.rs`, `export.rs::write_tiff_no_crs` |
| 14 | LIA + shadow/layover mask raster sidecars (`.lia.tif`, `.mask.tif`) | ✅ (CLI: `--write-lia`, `--write-mask`) | `terrain_correction.rs`, `export.rs` |
| 15 | Provenance JSON sidecar (`.provenance.json`, schema v1) | ✅ | `provenance.rs` |
| 16 | Speckle filter on linear-power buffer (boxcar / Lee 1981 / Frost / Gamma-MAP / Refined Lee) | ✅ (CLI: `--speckle KIND`) | `speckle.rs` |
| 17 | Selectable output CRS (EPSG:4326, UTM 326XX/327XX, `auto`) via `proj4rs` | ✅ (CLI: `--crs`, `--pixel-spacing-m`) | `output_crs.rs`, `terrain_correction.rs` |
| 18 | Library entry points (`sardine::run::{run_process, run_grd, …}`) | ✅ | `run.rs` |
| 19 | Python bindings (`sardine.process` / `grd` / `fetch_orbit` / `download_slc` / `fetch_geoid` / `features`) | ✅ (`sardine-py` crate, maturin) | `sardine-py/src/lib.rs` |
| — | Coherence / polarimetric workflows | ❌ Foundations not present (intensity-only deburst) | — |
| — | BigTIFF (>4 GiB) output | ❌ Not implemented; classic-TIFF size enforced pre-flight | — |

## Module architecture (implemented)

All in `sardine/src/` (line counts current as of 2026-04-24):

```
types.rs              domain types, no I/O
validate.rs           24 invariant checks; collects all errors
parse/mod.rs          parse_safe_directory, parse_calibration_noise, parse_geolocation_grids
parse/xml.rs          serde XML structs (private)
orbit.rs              EOF parser, 8-pt Lagrange interp, refuses extrapolation
orbit_fetch.rs        POEORB downloader (feature `orbit-fetch`)
geodesy.rs            WGS84 ↔ ECEF, Bowring inverse, vector algebra
calibration.rs        cal + noise XML parser
slc_reader.rs         pure-Rust CInt16 TIFF reader (no GDAL)
deburst.rs            midpoint-selection deburst; intensity-only
apply_calibration.rs  σ⁰ + per-pixel NESZ; bilinear LUT
merge_subswaths.rs    integer-offset slant-range merge; midpoint seam
dem.rs                SRTM-1 + GLO-30 loaders; DemMosaic; DemSource trait
geoid.rs              GeoidModel::{Zero, Egm96}; bilinear with lon-wrap
geoid_fetch.rs        EGM96 downloader from PROJ CDN (feature `geoid-fetch`)
terrain_correction.rs backward geocoding + Small 2011 flattening; OutputCrs-aware
ground_range.rs       slant→ground projection + multilook (wired by `sardine grd`)
speckle.rs            Boxcar / Lee 1981 / Frost / Gamma-MAP / Refined Lee on linear power
output_crs.rs         OutputCrs enum + Projector (proj4rs); EPSG:4326 + UTM 326XX/327XX + auto
export.rs             to_db_inplace + write_geotiff (CRS-aware) + write_tiff_no_crs (radar geometry); 4 GiB pre-flight
provenance.rs         Provenance JSON sidecar (schema v1, polymorphic on TC vs GRD)
slc_fetch.rs          ASF DAAC SLC downloader (feature `slc-fetch`)
run.rs                Library entry points: run_process, run_grd, prepare_merged_scene, resolve_crs/geoid/speckle, …
bin/sardine.rs        Thin clap shim over run.rs (process, grd, inspect, fetch-orbit, download-slc)
```

Python extension crate `sardine-py/`:

```
src/lib.rs                       PyO3 wrappers: process, grd, fetch_orbit, download_slc, fetch_geoid, features
python/sardine/__init__.py       Re-exports from compiled `_sardine` module
pyproject.toml                   maturin config (module-name = sardine._sardine, abi3-py38)
Cargo.toml                       Standalone crate; passes through fetch features to sardine
```

Guard: `scripts/check_no_silent_fallbacks.sh` + `tests/no_silent_fallbacks.rs`
Examples: `examples/dump_merged_sigma0.rs`, `examples/dump_s1b_tc.rs`
(both require `--orbit` POEORB or `SARDINE_ALLOW_ANNOTATION_ORBIT=1`).

## Known architectural blockers (must be resolved before adding new branches)

1. **Pipeline is monomorphic over `MergedSigma0`.** No `RadarImage` trait. A
   complex / polarimetric / coherence branch cannot reuse the geocoding code
   without inflating that struct or duplicating the chain.
2. **Deburst is intensity-only (midpoint selection).** Forecloses coherence
   and polarimetric workflows by design — those need deramp + reramp +
   coherent overlap blending on complex SLC.
3. **Classic TIFF output cap (≤ 4 GiB).** `export.rs` enforces a pre-flight
   size check; large UTM scenes at 10 m spacing approach the limit. Adding
   BigTIFF requires a rewrite of `write_geotiff` (or a switch to
   `tiff = "…"` with the `bigtiff` feature) plus equivalent updates to the
   no-CRS writer.

Resolved since the previous revision:

* **Output CRS is selectable.** `output_crs.rs` defines `OutputCrs::{Wgs84,
  Utm{epsg, zone, hemisphere}}`; `Projector` (proj4rs) projects per-pixel
  inside the geocoding loop. `terrain_correction.rs` builds the output grid
  in CRS units (`pixel_spacing_m` for UTM, degrees for WGS84) and
  `write_geotiff` embeds the right tags. Driven by `--crs <SPEC>` and
  `--pixel-spacing-m`.
* **Python bindings.** `sardine-py` exposes `process`, `grd`, `fetch_orbit`,
  `download_slc`, `fetch_geoid`, and `features` via PyO3 (abi3-py38). All
  functions release the GIL and surface failures as `RuntimeError`.
* **Library entry points.** Pipeline orchestration moved out of
  `bin/sardine.rs` into `run.rs` (`ProcessOptions`, `GrdOptions`,
  `run_process`, `run_grd`, `prepare_merged_scene`, helpers). The clap
  binary is now a thin shim and the same code path is used by the Python
  extension.
* `speckle.rs` ships five filters (boxcar, Lee 1981, Frost, Gamma-MAP, Refined Lee)
  operating on the linear-power buffer; wired into both `sardine process`
  (post-geocoding, pre-dB) and `sardine grd` (post-ground-range, pre-write)
  via `--speckle`, `--speckle-window`, `--enl`, `--frost-damping`.  The
  selection is recorded in `provenance.json`.
* `GeoidModel::Zero` no longer defaults — `--geoid` is required, with
  `auto`, an explicit path, or `zero` (the last requires opting into the
  ±80 m bias explicitly).
* `ground_range.rs` is wired into the CLI as `sardine grd` and writes a
  non-georeferenced TIFF via `export::write_tiff_no_crs`.
* LIA and shadow/layover mask sidecars are writable via `--write-lia` /
  `--write-mask`.
* Pre-flight DEM coverage check (`DemMosaic::covers_bbox`) runs in
  `run.rs::prepare_merged_scene` before the first heavy step.

## Roadmap pointer

For the prioritised list of what to do next, see
[HANDOVER.md §16 — Next Steps (April 2026 review)](HANDOVER.md) and
[PROGRESS.md §12 — Suggested Next Steps](PROGRESS.md).