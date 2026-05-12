# SARdine

Sentinel-1 IW SLC backscatter processing pipeline in Rust.

> **Status (May 2026):** pre-1.0. The end-to-end pipeline is implemented and
> empirically validated against ASF RTC10 GAMMA at **+0.016 dB median linear
> bias** on S1B (10×10 multilook). **358 unit tests + 1 guard integration
> test** pass (`cargo test`).

The active codebase is a single Rust crate at
[`new_pipeline/scene/`](new_pipeline/scene/) (`sardine-scene` v0.1.0). Python
bindings live at [`new_pipeline/sardine-py/`](new_pipeline/sardine-py/)
(`sardine` PyO3 wheel, built with maturin). The Python-and-Rust package under
[`legacy/`](legacy/) is **reference-only** — see
[LEGACY_STATUS.md](LEGACY_STATUS.md) and [AGENTS.md](AGENTS.md).

## What works

- **Full terrain-correction pipeline** (`sardine process`): parse SAFE →
  POEORB precise orbit → calibration + per-pixel NESZ → TOPS midpoint deburst
  → IW1+IW2+IW3 merge → DEM mosaic (SRTM-1 or Copernicus GLO-30) → backward
  Range-Doppler geocoding → Small (2011) terrain flattening (σ⁰ → γ⁰) → dB
  → GeoTIFF
- **Ground-range projection** (`sardine grd`): deburst + calibrate + project
  to ground range, no DEM required
- Single-polarization per run (`--polarization VV` or `VH`)
- Output CRS: EPSG:4326, UTM 326XX/327XX, or `auto` (UTM zone from scene
  centre); selectable pixel spacing (`--pixel-spacing-m`)
- Output formats: stripped Float32 GeoTIFF, tiled Cloud-Optimised GeoTIFF
  (`--cog`), BigTIFF for outputs >4 GiB
- Speckle filtering on linear power before dB conversion: `none`, `boxcar`,
  `lee`, `frost`, `gamma_map`, `refined_lee` (`--speckle KIND --speckle-window N`)
- LIA and shadow/layover mask raster sidecars (`--write-lia`, `--write-mask`)
- Provenance JSON sidecar written next to every output GeoTIFF
- Per-pixel NESZ noise-floor masking (`--noise-floor-margin-db`)
- EGM96 geoid correction (`--geoid auto|<path>|zero`); `--geoid` is **required**
  (no silent default)
- Pre-flight DEM coverage check before heavy processing starts
- Auto-fetch POEORB from ESA STEP / AWS S3 (`--features orbit-fetch`)
- SLC download from ASF Vertex (`--features slc-fetch`)
- EGM96 geoid auto-fetch and cache (`--features geoid-fetch`)
- Python bindings: `sardine.process()`, `sardine.grd()`, `sardine.fetch_orbit()`,
  `sardine.download_slc()`, `sardine.fetch_geoid()`, `sardine.features()`
- No GDAL dependency in `sardine-scene`; pure-Rust TIFF reader and writer

## What is **not** in the active code

- No single CLI call for dual-pol VV+VH (run twice, once per polarization)
- No coherence / polarimetric workflow (intensity-only deburst by design)
- No batch / multi-scene mosaic driver
- No anti-meridian handling (scenes crossing ±180° produce incorrect bboxes)
- IW mode only (EW and StripMap acquisition modes are not supported)

## Quick start

### Rust CLI

```sh
cd new_pipeline/scene
cargo test                              # 358 unit + 1 guard test must pass
cargo build --release

# Inspect a SAFE product
cargo run --release -- inspect /path/to/S1B.SAFE

# Full terrain-correction pipeline
cargo run --release --features orbit-fetch,dem-fetch,geoid-fetch -- process \
    --safe         /path/to/S1B.SAFE \
    --dem          /path/to/dem_tiles_dir/ \
    --output       sardine_out.tiff \
    --polarization VV \
    --geoid        auto \
    --crs          auto \
    --cog                               # write Cloud-Optimised GeoTIFF

# Auto-fetch orbit (no --orbit flag needed with --features orbit-fetch):
#   tries AWS S3 s1-orbits first, falls back to ESA STEP
# --geoid auto fetches and caches EGM96 (needs --features geoid-fetch)
# Set RUST_LOG=sardine_scene=debug for detailed progress
```

### Python bindings

```sh
cd new_pipeline/sardine-py
python3 -m venv .venv && source .venv/bin/activate
pip install maturin pytest
maturin develop --features slc-fetch
python3 -c "import sardine; print(sardine.features())"
```

## Documentation

- [docs/HANDOVER.md](docs/HANDOVER.md) — full developer handover (current state,
  validation results, open items)
- [docs/PROGRESS.md](docs/PROGRESS.md) — step-by-step build log and next steps
- [docs/architecture_overview.md](docs/architecture_overview.md) — module map
- [AGENTS.md](AGENTS.md) — non-negotiable working rules for contributors and AI agents
- [LEGACY_STATUS.md](LEGACY_STATUS.md) — why the legacy package is reference-only

## License

See repository.
