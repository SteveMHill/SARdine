# SARdine

Sentinel-1 IW SLC backscatter processing pipeline in Rust.

> **Status (April 24, 2026):** pre-1.0. The end-to-end pipeline is
> implemented and empirically validated against ASF RTC10 GAMMA at
> **+0.016 dB median linear bias** on S1B (10×10 multilook). 165 tests
> pass (`cargo test`).

The active codebase is a single Rust crate at
[`new_pipeline/scene/`](new_pipeline/scene/) (`sardine-scene` v0.1.0). The
Python-and-Rust package under [`legacy/`](legacy/) is **reference-only** —
see [LEGACY_STATUS.md](LEGACY_STATUS.md) and [AGENTS.md](AGENTS.md).

## What works

- Full pipeline: parse SAFE → POEORB orbit → calibration + per-pixel NESZ
  → TOPS midpoint deburst → IW1+IW2+IW3 merge → DEM mosaic
  (SRTM-1 or Copernicus GLO-30) → backward Range-Doppler geocoding
  → Small (2011) terrain flattening (σ⁰ → γ⁰) → dB → GeoTIFF
- Single-polarization (`VV` or `VH`, one per run)
- Output: Float32 GeoTIFF, **EPSG:4326 only**, classic TIFF (≤4 GiB)
- Optional EGM96 geoid (feature `geoid-fetch`); precise-orbit auto-fetch
  (feature `orbit-fetch`); SLC download from ASF (feature `slc-fetch`)

## What is **not** in the active code

- No Python bindings
- No speckle filter
- No output CRS other than EPSG:4326
- No COG / BigTIFF output
- No dual-pol single-call processing
- No coherence / polarimetric workflow (intensity-only deburst)
- No batch / multi-scene driver
- No layover/shadow or LIA mask raster outputs (computed, then discarded)

## Quick start

```sh
cd new_pipeline/scene
cargo test                              # 165 tests must pass
cargo build --release

cargo run --release -- inspect /path/to/S1B.SAFE

cargo run --release --features geoid-fetch -- process \
    --safe   /path/to/S1B.SAFE \
    --orbit  /path/to/POEORB.EOF \
    --dem    /path/to/dem_tiles_dir/ \
    --output sardine_out.tiff \
    --polarization VV \
    --geoid auto                          # 'auto' fetches EGM96 (needs --features geoid-fetch)
```

## Documentation

- [docs/HANDOVER.md](docs/HANDOVER.md) — full developer handover
- [docs/architecture_overview.md](docs/architecture_overview.md) — module map
- [docs/PROGRESS.md](docs/PROGRESS.md) — current state and prioritised next steps
- [AGENTS.md](AGENTS.md) — non-negotiable working rules for contributors and AI agents
- [LEGACY_STATUS.md](LEGACY_STATUS.md) — why the legacy package is reference-only

## License

See repository.
