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

SARdine reads a Sentinel-1 IW SLC product (`.SAFE` directory) and produces
calibrated, terrain-corrected, geocoded output rasters with no external
dependencies beyond the Rust toolchain. There is no GDAL dependency; TIFF
reading and writing are handled in pure Rust.

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
9. **PolSAR decomposition** — Cloude–Pottier H/A/Alpha from the dual-pol
   covariance matrix C2; geocoded entropy, anisotropy, and mean scattering
   angle bands (`sardine polsar`)
10. **InSAR coherence** — per-IW boxcar coherence estimator and wrapped phase
    for an SLC pair (`sardine insar`)
11. **Export** — Float32 GeoTIFF or Cloud-Optimised GeoTIFF (`--cog`),
    BigTIFF auto-selected for outputs > 4 GiB; EPSG:4326, UTM, or `auto`
12. **Provenance** — `.provenance.json` and `.stac.json` sidecars next to
    every output

**Output modes** (`sardine process --mode <MODE>` or dedicated subcommand):

| Mode | Subcommand | Description |
|------|-----------|-------------|
| `rtc` | `process` | Terrain-corrected σ⁰/γ⁰ — default |
| `nrb` | `process --mode nrb` | RTC + mandatory LIA, quality mask, STAC (NRB profile) |
| `grd` | `grd` | Ground-range radar geometry (no geocoding) |
| `polsar` | `polsar` | H/A/Alpha Cloude–Pottier decomposition (dual-pol) |
| — | `insar` | InSAR coherence + wrapped phase for an SLC pair |

**Python bindings** (PyO3, `sardine` wheel, `sardine-py/`):
`sardine.process()`, `sardine.grd()`, `sardine.insar()`, `sardine.polsar()`,
`sardine.fetch_orbit()`, `sardine.download_slc()`, `sardine.fetch_geoid()`,
`sardine.features()`.

---

## Quick start

### Requirements

- Rust stable ≥ 1.80 (`rustup update stable`)
- A Sentinel-1 IW SLC product (`.SAFE` directory), freely available from
  [ASF Vertex](https://search.asf.alaska.edu/) or
  [Copernicus Dataspace](https://dataspace.copernicus.eu/)

### Build and test

```sh
cargo test -p sardine --lib   # 399 unit tests must pass
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

### PolSAR decomposition (H/A/Alpha)

```sh
cargo run --release --features orbit-fetch,dem-fetch,geoid-fetch -- polsar \
    --safe              /path/to/S1B_IW_SLC__1SDV_....SAFE \
    --output            output_polsar \
    --geoid             auto \
    --polarization      VV+VH \
    --multilook-range   3 \
    --crs               auto
# writes output_polsar_H.tif, output_polsar_A.tif, output_polsar_alpha.tif
# plus .polsar.provenance.json and .polsar.stac.json sidecars
```

### InSAR coherence

```sh
cargo run --release --features orbit-fetch,dem-fetch,geoid-fetch -- insar \
    --reference  /path/to/S1B_reference.SAFE \
    --secondary  /path/to/S1B_secondary.SAFE \
    --output     coherence_out \
    --geoid      auto \
    --az-looks   10 \
    --rg-looks   2
# writes coherence_out_iw{1,2,3}_coherence.tif + sidecars
```

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
maturin develop --release --features orbit-fetch,dem-fetch,geoid-fetch,slc-fetch
```

```python
import sardine

# Terrain-corrected backscatter
sardine.process(
    safe="/path/to/S1B.SAFE",
    dem="/path/to/dem_tiles/",
    output="sigma0.tif",
    geoid="auto",
    crs="auto",
    cog=True,
)

# H/A/Alpha PolSAR decomposition (writes _H.tif, _A.tif, _alpha.tif)
sardine.polsar(
    safe="/path/to/S1B_IW_SLC__1SDV_....SAFE",
    output="polsar_out",
    geoid="auto",
    polarization="VV+VH",
    multilook_range=3,
    crs="auto",
)

# InSAR coherence
sardine.insar(
    reference="/path/to/S1B_ref.SAFE",
    secondary="/path/to/S1B_sec.SAFE",
    output="coherence_out",
    geoid="auto",
    az_looks=10,
    rg_looks=2,
)

# Check which auto-fetch features are compiled in
print(sardine.features())
# {'geoid_fetch': True, 'orbit_fetch': True, 'slc_fetch': True}
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
- InSAR phase correctness requires validated flat-earth deramping (coherence
  is reliable; treat wrapped phase as experimental)
- No multi-scene mosaic / burst-stitch driver

---

## License

Apache License 2.0 — see [`LICENSE`](LICENSE) and [`NOTICE`](NOTICE).
