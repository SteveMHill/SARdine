# SARdine

**Sentinel-1 IW SLC backscatter processing pipeline in pure Rust.**

> ⚠️ **Prototype — not production-ready.**
> SARdine is research and reference software. It is not validated for
> operational use, does not carry a stability guarantee, and has known
> limitations (see [below](#limitations)). Use results at your own risk and
> always verify outputs against a trusted reference product.

[![CI](https://github.com/SteveMHill/SARdine/actions/workflows/ci.yml/badge.svg)](https://github.com/SteveMHill/SARdine/actions/workflows/ci.yml)

---

## What is SARdine?

SARdine reads a Sentinel-1 IW SLC product (`.SAFE` directory) and produces a
calibrated, terrain-corrected, geocoded backscatter raster (σ⁰ or γ⁰) with no
external dependencies beyond the Rust toolchain. There is no GDAL dependency;
TIFF reading and writing are handled in pure Rust.

The pipeline implements:

1. **Parse** — annotation XMLs, calibration XMLs, noise XMLs
2. **Precise orbit** — apply POEORB/RESORB (auto-fetch from ESA STEP / AWS S3
   with `--features orbit-fetch`)
3. **Calibration** — σ⁰ = (|DN|² − N) / K² using the LUT-based formula from
   the Sentinel-1 Product Specification; per-pixel NESZ noise-floor correction
4. **TOPS deburst** — midpoint-selection intensity deburst across all 9 bursts
   per IW subswath
5. **IW merge** — pixel-exact integer-offset merge of IW1 + IW2 + IW3
6. **Terrain correction** — backward Range-Doppler geocoding with Newton
   zero-Doppler solver; EGM96 or EGM2008 geoid undulation (auto-fetch with
   `--features geoid-fetch`); DEM from SRTM-1 tiles or Copernicus GLO-30
   (auto-fetch with `--features dem-fetch`)
7. **Terrain flattening** — Small (2011) γ⁰ correction (optional, `--flatten`)
8. **Speckle filtering** — boxcar, Lee (1981), Frost (1982), Gamma-MAP (Lopes
   1990), Refined Lee (Lopes–Touzi–Nezry 1990); all behind a `SpeckleKernel`
   trait for third-party extensions
9. **Export** — Float32 GeoTIFF or Cloud-Optimised GeoTIFF (`--cog`),
   BigTIFF auto-selected for outputs > 4 GiB; EPSG:4326, UTM, or `auto`
10. **Provenance** — `.provenance.json` and `.stac.json` sidecars next to
    every output

**Python bindings** (PyO3, `sardine` wheel, `sardine-py/`):
`sardine.process()`, `sardine.grd()`, `sardine.fetch_orbit()`,
`sardine.download_slc()`, `sardine.fetch_geoid()`, `sardine.features()`.

---

## Quick start

### Requirements

- Rust stable ≥ 1.80 (`rustup update stable`)
- A Sentinel-1 IW SLC product (`.SAFE` directory), freely available from
  [ASF Vertex](https://search.asf.alaska.edu/) or
  [Copernicus Dataspace](https://dataspace.copernicus.eu/)

### Build and test

```sh
cd sardine
cargo test          # 372 unit + 1 guard test must pass
cargo build --release
```

### Inspect a SAFE product

```sh
cargo run --release -- inspect /path/to/S1B_IW_SLC__1SDV_....SAFE
```

### Full terrain-correction pipeline

```sh
cargo run --release --features orbit-fetch,dem-fetch,geoid-fetch -- process \
    --safe         /path/to/S1B.SAFE \
    --dem          /path/to/dem_tiles_dir/ \
    --output       output_sigma0.tif \
    --polarization VV \
    --geoid        auto \
    --crs          auto \
    --cog
```

Key flags:

| Flag | Default | Notes |
|------|---------|-------|
| `--polarization` | *(required)* | `VV`, `VH`, `VV+VH` (dual-pol writes two files) |
| `--geoid` | *(required)* | `auto` (fetch+cache EGM96), `egm2008`, or path to grid file |
| `--crs` | `epsg:4326` | `auto` (UTM from scene centre), `utm32n`, … |
| `--pixel-spacing-m` | 10.0 | Output pixel spacing in metres |
| `--flatten` | off | Apply Small (2011) terrain flattening (γ⁰) |
| `--speckle` | `none` | `boxcar`, `lee`, `frost`, `gamma_map`, `refined_lee` |
| `--resampling` | `bilinear` | `bicubic`, `lanczos3` |
| `--cog` | off | Write Cloud-Optimised GeoTIFF |
| `--write-lia` | off | Write local incidence angle sidecar |
| `--write-mask` | off | Write shadow/layover mask sidecar |
| `--noise-floor-margin-db` | off | Mask pixels within N dB of NESZ |

Auto-fetch behaviour (requires the corresponding feature flag at compile time):

- **`--features orbit-fetch`** — `--orbit` can be omitted; POEORB fetched from
  AWS S3 `s1-orbits`, fallback to ESA STEP
- **`--features dem-fetch`** — `--dem` can be omitted; Copernicus GLO-30 tiles
  fetched and cached
- **`--features geoid-fetch`** — `--geoid auto` downloads and caches the EGM96
  grid (≈ 7 MB)

Set `RUST_LOG=sardine=info` for a progress log.

### Python bindings

```sh
cd sardine-py
python3 -m venv .venv && source .venv/bin/activate
pip install maturin
maturin develop --release --features slc-fetch
```

```python
import sardine

result = sardine.process(
    safe_path="/path/to/S1B.SAFE",
    output_path="output.tif",
    polarization="VV",
    geoid="auto",
    crs="auto",
    cog=True,
)
print(result)
```

---

## Repository layout

```
sardine/            sardine Rust crate (the pipeline)
sardine-py/         Python bindings (PyO3 / maturin)
sardine-validate/   End-to-end radiometric regression harness
docs/               Architecture notes, handover doc, progress log
scripts/            Validation and diagnostic Python scripts
data/               Test data (not committed — see .gitignore)
```

---

## Limitations

- **Prototype only** — outputs have not been independently validated beyond the
  single S1B scene referenced above
- IW acquisition mode only (EW, StripMap not supported)
- Anti-meridian scenes (crossing ±180°) produce incorrect bounding boxes
- No coherence or polarimetric decomposition (intensity-only deburst)
- No batch/multi-scene mosaic driver

---

## License

Apache License 2.0 — see [`LICENSE`](LICENSE) and [`NOTICE`](NOTICE).
