git clone https://github.com/SteveMHill/SARdine.git
# SARdine - Experimental SAR Processing Package

## 🛰️ SARdine

### Experimental SAR Processing Toolkit

```
⚠️  EXPERIMENTAL - NOT PRODUCTION READY  ⚠️
```

SARdine is an *experimental* Synthetic Aperture Radar (SAR) processing playground. It is designed for research, learning, and prototyping—not for mission-critical production use. Expect sharp edges, rapid iteration, and breaking changes.

- **Status:** Development / Research Tool
- **Stability:** Experimental
- **Support:** Community-driven exploration only

At its core, SARdine couples a Rust processing engine with Python orchestration to help you experiment with Sentinel-1 backscatter workflows, tweak algorithms, and observe impacts quickly.

---

## 🎯 Quick Start

Run the backscatter pipeline end-to-end:

```shell
python -m sardine.cli backscatter path/to/S1A_IW_SLC__*.SAFE ./output/
```

More variations for experimentation:

```shell
# Default VV processing
python -m sardine.cli backscatter S1A_IW_SLC__*.SAFE ./my_output/

# VH polarization with a different speckle filter
python -m sardine.cli backscatter S1A_IW_SLC__*.SAFE ./output/ \
    --polarization VH --speckle-filter lee

# Higher nominal ground resolution
python -m sardine.cli backscatter S1A_IW_SLC__*.SAFE ./output/ \
    --resolution 5 --filter-window 5

# Faster prototyping (skip terrain flattening & geocoding)
python -m sardine.cli backscatter S1A_IW_SLC__*.SAFE ./output/ \
    --no-geocode --no-terrain-flatten
```

---

## � Package Structure (high level)

```
SARdine/
├── backscatter_cli.py*            # Convenience wrapper around the CLI module
├── examples/                      # Usage demonstrations & quick sanity checks
├── tests/                         # Pytest suite for regression scenarios
├── data/                          # Example Sentinel-1 inputs (not tracked)
├── docs/                          # Documentation & experiment notes
│   ├── user-guide/
│   ├── scientific/
│   └── development/
├── SARdine/                       # Core library (Rust crate + Python bindings)
│   ├── src/                       # Rust implementation
│   ├── python/                    # Processors, CLI entry points, utilities
│   └── pyproject.toml             # Build metadata for the Python package
└── archive/                       # Legacy experiments kept for reference
```

\* `backscatter_cli.py` is optional. The canonical launcher is `python -m sardine.cli backscatter`, but a lightweight script can mirror those arguments if you prefer.

---

## 🧪 Testing (research-focused)

```shell
# High-level orbit + metadata regression
python -m pytest tests/test_orbit_debug.py

# Enhanced deburst & merge validation
python -m pytest tests/test_enhanced_deburst.py

# Strict metadata parsing checks
python -m pytest tests/test_subswaths_strict.py
```

Feel free to craft additional pytest modules in `tests/` to capture new hypotheses or datasets.

---

## � Processing Steps Implemented

1. **SLC Reading** – parse SAFE archives and prime cache metadata
2. **Precise Orbit Application** – download & integrate POE orbits
3. **IW Subswath Handling** – manage burst timing per polarization
4. **TOPSAR Deburst** – stitch bursts with orbit-aware phasing
5. **Radiometric Calibration** – compute σ⁰ / β⁰ / γ⁰ backscatter
6. **Multilooking** – configurable looks for speckle reduction
7. **Speckle Filtering** – Enhanced Lee and other research filters
8. **Terrain Flattening** – γ⁰ normalization for steep terrain
9. **Terrain Correction & Geocoding** – DEM-driven range–Doppler solver
10. **Product Export** – GeoTIFF, NumPy, and log summaries

Every stage exposes hooks for parameter tuning so you can observe how algorithm tweaks ripple through the pipeline.

---

## 📖 Documentation Highlights

- `docs/user-guide/README.md` – orientation & usage walkthroughs
- `docs/scientific/SCIENTIFIC_AUDIT_REPORT.md` – algorithm references and validation notes
- `docs/scientific/COMPLETE_14_STEP_PIPELINE_SUCCESS_REPORT.md` – narrative of the full IW pipeline

These documents evolve alongside the research effort; treat them as living notes.

---

## � Requirements

- Python **3.8+** (tested most on 3.10)
- Rust toolchain (for building the core crate)
- GDAL / PROJ libraries available on your system
- Sentinel-1 SLC scenes (SAFE directories or ZIP archives)

Install the editable package after cloning:

```shell
cd SARdine/SARdine
pip install -e .
```

Rebuild the native extension when iterating on Rust code:

```shell
cargo build --release
```

---

## ⚠️ Important Notes

- **Not production ready:** Interfaces, file formats, and algorithms may change abruptly.
- **Limited validation:** Most experiments target a handful of Sentinel-1 IW acquisitions.
- **No warranty:** Always verify outputs against trusted tools (e.g., ESA SNAP, GAMMA).
- **Use at your own risk:** Suitable for labs, prototypes, and academic exploration only.

---

## 🤝 Contributing

We welcome contributions that improve the experimental toolkit:

- Hardening or profiling specific pipeline stages
- Documenting new experiments and data findings
- Expanding test coverage for additional Sentinel-1 scenes
- Porting algorithms or suggesting alternative approaches

Open a pull request with a clear description, reproduction steps, and (when possible) sample data references.

---

## 📄 License

Released under the **MIT License** – see [`LICENSE`](../LICENSE) for details.

---

## 🏷️ Versioning

Current status: **Experimental / Development**. Semantic versioning will begin once the API surface stabilizes toward a production-ready release.

---

**Reminder:** SARdine is for learning and experimentation. For production SAR processing, rely on established suites like ESA SNAP, GAMMA, or commercial offerings.
