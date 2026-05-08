# New Pipeline Area

Contract-driven rebuild area for new verified modules.
Only narrow, testable, justified implementations are added here.

## Current state (April 24, 2026)

The active codebase is a **single Rust crate**: [`scene/`](scene/)
(`sardine-scene` v0.1.0). It implements the full Sentinel-1 IW SLC
backscatter pipeline: parse → orbit → calibrate → deburst → merge → DEM →
backward Range-Doppler geocode → Small (2011) terrain flatten → dB → GeoTIFF.

- **Tests:** 168 unit + 1 forbidden-pattern guard integration test, all passing
  (`cargo test`).
- **Validation:** +0.016 dB median linear bias against ASF RTC10 GAMMA on
  S1B (2019-01-23) after 10×10 multilook.
- **Output:** Float32 GeoTIFF, EPSG:4326, classic TIFF (≤4 GiB).

## Quick start

```sh
cd scene
cargo test                      # 169 tests must pass
cargo build --release

cargo run --release -- inspect /path/to/S1B.SAFE

cargo run --release --features geoid-fetch -- process \
    --safe   /path/to/S1B.SAFE \
    --orbit  /path/to/POEORB.EOF \
    --dem    /path/to/dem_tiles_dir/ \
    --output sardine_out.tiff \
    --polarization VV \
    --geoid auto
```

## What this is **not** (yet)

- No Python bindings (the crate is Rust-only)
- No speckle filter
- No output CRS choice (EPSG:4326 only)
- No COG / BigTIFF output
- No coherence / polarimetric branch (intensity-only deburst)
- No batch / multi-scene driver

See [../docs/HANDOVER.md](../docs/HANDOVER.md) for the full handover
document, [../docs/architecture_overview.md](../docs/architecture_overview.md)
for the module map, and [../docs/PROGRESS.md](../docs/PROGRESS.md) for
the prioritised next-steps list. Working rules for any contributor or
agent are in [../AGENTS.md](../AGENTS.md).