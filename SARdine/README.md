# SARdine: SAR Processing Development Playground

<div align="center">
  <img src="logo.png" alt="SARdine Logo" width="300"/>
  
  [![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
  [![Rust](https://img.shields.io/badge/rust-1.70+-blue.svg)](https://www.rust-lang.org)
  [![Python](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org)
</div>

## **DEVELOPMENT STATUS - NOT PRODUCTION READY**

**This repository is a development playground for testing and experimenting with SAR processing algorithms. This package is NOT production-ready and should only be used for research, learning, and experimentation purposes.**

## Overview

SARdine is an experimental SAR data processing library for Sentinel-1 data, implemented in Rust with Python bindings. It aims to provide a complete 14-step processing pipeline from SLC to analysis-ready products, designed as a learning platform for SAR data processing workflows.

## Acknowledgments & Inspiration

This project is heavily inspired by and builds upon the excellent work of several established SAR processing packages:

- **[pyroSAR](https://github.com/johntruckenbrodt/pyroSAR)** - A Python framework for large-scale SAR satellite data processing
- **[sarsen](https://github.com/bopen/sarsen)** - Sentinel-1 SAR geocoding with xarray and rasterio  
- **[OpenSARToolkit](https://github.com/ESA-PhiLab/OpenSarToolkit)** - ESA's open-source toolkit for SAR data processing
- **[SNAP](https://step.esa.int/main/toolboxes/snap/)** - ESA's Sentinel Application Platform
- **[GAMMA](https://www.gamma-rs.ch/)** - Professional SAR processing software

We acknowledge and appreciate the foundational work these projects have contributed to the SAR processing community.

## Experimental Features

- **High-Performance**: Rust backend with parallel processing experiments
- **Python Integration**: Python API development with numpy compatibility  
- **Processing Pipeline**: Experimental 14-step SLC to backscatter workflow
- **Multi-Swath**: Testing IW swaths (IW1, IW2, IW3) processing
- **Dual-Polarization**: VV + VH processing experiments
- **Terrain Processing**: DEM integration experiments
- **Radiometric Calibration**: σ⁰, β⁰, γ⁰ calibration testing
- **Speckle Filtering**: Filter algorithm implementations
- **Output Formats**: GeoTIFF export testing

## Experimental 14-Step Pipeline

### **Implementation Status (Work in Progress)**

1. **SLC Reading & Metadata Extraction** - File parsing experiments
2. **Orbit Data Application** - Orbit integration testing
3. **IW Split** - Multi-swath data separation (experimental)
4. **Debursting** - Burst concatenation development
5. **Radiometric Calibration** - Multi-polarization calibration testing
6. **IW Merge** - Swath combination experiments
7. **Multilooking** - Speckle reduction (configurable)
8. **Terrain Flattening** - Topographic normalization testing
9. **Speckle Filtering** - Noise reduction algorithms (Lee filter)
10. **Terrain Correction** - Geocoding experiments
11. **Masking** - Quality control development
12. **dB Conversion** - Logarithmic scaling testing
13. **Export** - GeoTIFF export testing
14. **Metadata Generation** - Processing documentation experiments

### **Development Status**
- **Processing Pipeline**: Under active development
- **Output Quality**: Experimental implementations
- **Testing**: Limited validation with sample data
- **Documentation**: Work in progress

## **Important Disclaimers**

- **NOT FOR PRODUCTION USE**: This is experimental software for research and learning
- **Data Loss Risk**: Always backup your original data before processing
- **Quality Not Guaranteed**: Output products may contain errors or artifacts
- **API Instability**: Function signatures and behavior may change frequently
- **Limited Support**: Community-driven development with no warranty

## Experimental Installation

### Prerequisites

- **Rust**: 1.70+ ([Install Rust](https://rustup.rs/))
- **Python**: 3.8+ with pip
- **System Dependencies**: 
  ```bash
  # Ubuntu/Debian
  sudo apt-get install build-essential pkg-config libssl-dev libgdal-dev
  
  # macOS (with Homebrew)
  brew install rust gdal
  
  # Windows
  # Install Visual Studio Build Tools and GDAL
  ```

### Installation

```bash
# Clone the repository
git clone https://github.com/SteveMHill/SARdine.git
cd SARdine

# Install Python dependencies
pip install maturin numpy rasterio pyproj

# Build and install (experimental build)
maturin develop --release

# Verify installation
python -c "import sardine; print('SARdine experimental build installed')"
```

## Experimental Python API

```python
import sardine
import numpy as np

# EXPERIMENTAL: Basic terrain flattening
sigma0 = np.array(...)  # Your SAR backscatter data
dem = np.array(...)     # Digital elevation model

# Apply experimental terrain flattening
result = sardine.apply_terrain_flattening(sigma0, dem)
gamma0 = result['gamma0']  # Terrain-flattened backscatter
incidence_angles = result['incidence_angles']

# EXPERIMENTAL: Speckle filtering
filtered = sardine.apply_speckle_filter_numpy(
    sigma0, "enhanced_lee", window_size=7, num_looks=1.0
)
```

## Experimental Performance

Performance testing is ongoing:

- **Terrain Flattening**: Under development and testing
- **Speckle Filtering**: Experimental multi-core processing
- **Memory Usage**: Optimization in progress
- **Scalability**: CPU utilization experiments

## Development Examples

See the `examples/` directory for experimental usage examples:

- `complete_backscatter_pipeline.py` - Experimental SLC processing
- `complete_terrain_correction_workflow.py` - Terrain correction testing
- `complete_speckle_filtering_workflow.py` - Speckle filtering experiments
- `production_backscatter_processor.py` - Development workflow testing

## Experimental Testing

```bash
# Run experimental tests  
python tests/test_basic.py

# Run development tests
python tests/test_pipeline.py
python tests/test_terrain_flattening.py

# Run Rust unit tests
cargo test --release
```

## Contributing to the Playground

We welcome contributions to this experimental project:

1. Fork the repository
2. Create a feature branch (`git checkout -b experiment/new-feature`)
3. Test your experimental changes
4. Commit your changes (`git commit -m 'Add experimental feature'`)
5. Push to the branch (`git push origin experiment/new-feature`)
6. Open a Pull Request

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Acknowledgments

- **pyroSAR, sarsen, and OpenSARToolkit** for foundational SAR processing concepts
- ESA for Sentinel-1 data and comprehensive documentation
- The Rust and Python communities for excellent development tooling
- GDAL/OGR for robust geospatial data handling capabilities

---

**SARdine**: Experimental SAR processing playground

*Remember: This is experimental software. Use at your own risk and always validate results with established tools.*
