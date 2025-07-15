# SARius: A Fast, Modular Sentinel-1 Backscatter Processor

## 📦 Project Overview

SARius is an open-source Python+Rust library for processing Sentinel-1 SLC data into calibrated, terrain-corrected backscatter products (Gamma0). It is designed as a modern alternative to ESA SNAP and GAMMA, with a focus on speed, modularity, and full transparency.

---

## 🎯 Goals

- Provide a modular open-source pipeline for Sentinel-1 IW SLC processing
- Replace slow and closed-source tools (SNAP, GAMMA)
- Expose a clean Python API powered by a fast Rust backend
- Output high-quality Gamma0 VV/VH GeoTIFFs ready for analysis

---

## 📐 Pipeline Steps

1. Read SLC
2. Apply orbit file
3. IW split
4. Deburst
5. Radiometric calibration (|SLC|^2 \* beta0)
6. Multilooking
7. Terrain flattening (gamma0 = sigma0 / cos(theta\_lia))
8. Speckle filtering
9. Terrain correction (map projection)
10. Output (VV, VH Gamma0 GeoTIFFs)

---

## 🛠 Draft Directory Layout

```
SARius/
├── Cargo.toml              # Rust config
├── pyproject.toml          # Python-Rust binding
├── src/
│   ├── lib.rs              # PyO3 interface
│   ├── types.rs            # Core structs and types
│   ├── io/
│   │   ├── slc_reader.rs
│   │   ├── orbit.rs
│   │   ├── annotation.rs
│   │   └── dem.rs
│   ├── core/
│   │   ├── deburst.rs
│   │   ├── calibrate.rs
│   │   ├── multilook.rs
│   │   ├── terrain_flatten.rs
│   │   ├── speckle.rs
│   │   └── terrain_correction.rs
├── py/
│   └── api.rs              # PyO3 bindings to Rust core
└── examples/
    └── process_slc.py
```

---

## 🗺 Roadmap

### Week 1: Core Infrastructure

-

### Week 2: Core Processing

-

### Week 3: Terrain Integration

-

### Week 4: Refinement

-

---

## 📚 References & Sources

### 📄 ESA Documentation

- [Sentinel-1 Level-1 Product Definition](https://sentinel.esa.int/documents/247904/685163/Sentinel-1-Level-1-Product-Definition)
- [Detailed Algorithm Definition (DAD)](https://sentinel.esa.int/documents/247904/349490/Sentinel-1-Level-1-Detailed-Algorithm-Definition)
- [Sentinel-1 Toolbox (SNAP) Documentation](https://step.esa.int/main/toolboxes/snap/)

### 🔓 Open Source Projects

- [pyroSAR](https://github.com/johntruckenbrodt/pyroSAR)
- [ISCE2](https://github.com/isce-framework/isce2)
- [sarsen](https://github.com/Open-EO/sarsen)
- [MintPy](https://github.com/insarlab/MintPy)
- [GMTSAR](https://topex.ucsd.edu/gmtsar/)

### 📦 Tools & Crates

- [GDAL Rust Bindings](https://docs.rs/gdal)
- [Quick-XML for Rust](https://docs.rs/quick-xml)
- [ndarray + num-complex](https://docs.rs/ndarray)
- [PyO3 + maturin](https://pyo3.rs)

---

## 🤝 Contributions Welcome

If you're interested in helping build an open, fast, modern SAR processor, contributions are welcome. Let's replace black boxes with transparent, reproducible science.

License: MIT Author: Steven Hill and contributors

