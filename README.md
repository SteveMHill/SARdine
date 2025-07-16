# SARdine: A Fast, Complete Sentinel-1 SAR Processor 

SARdine is a modern, production-ready SAR data processing library for Sentinel-1 SLC products, implemented in Rust with Python bindings. It provides a complete processing pipeline from raw SLC to analysis-ready data products.

## 🌟 Key Features

- 🚀 **High-Performance Rust Backend**: Optimized core processing (>1.5M pixels/second)
- 🐍 **Intuitive Python API**: Easy-to-use interface with full feature access
- 📱 **Complete CLI Toolset**: Command-line tools for production batch processing
- 🛰️ **Automatic Orbit Handling**: Download and apply precise orbit files
- 📡 **IW Split & Deburst**: Extract sub-swaths and create seamless images
- 📊 **Radiometric Calibration**: Sigma0/Beta0/Gamma0 with bilinear interpolation
- 🏔️ **Automatic DEM & Terrain Flattening**: AWS SRTM/Copernicus DEM with topographic correction
- 🎯 **Advanced Speckle Filtering**: 8 filter types (Lee, Enhanced Lee, Frost, Gamma MAP, etc.)
- 📈 **Adaptive Multilooking**: Smart spatial averaging with noise estimation
- 🧪 **Production Ready**: Comprehensive testing and error handling

## 🎯 Complete SAR Processing Pipeline

**SLC Reading** → **Orbit Application** → **Debursting** → **Calibration** → **Terrain Flattening** → **Speckle Filtering** → **Multilooking** → **Analysis-Ready Products**

### ✅ Fully Implemented Components
- **Orbit Management**: Automatic download, validation, and precise interpolation
- **SLC I/O**: Efficient ZIP archive and metadata handling
- **IW Processing**: Sub-swath extraction and burst concatenation
- **Radiometric Calibration**: Complete σ⁰/β⁰/γ⁰ calibration pipeline
- **DEM Integration**: Automatic SRTM/Copernicus DEM download from AWS
- **Terrain Flattening**: Local incidence angle correction with DEM
- **Speckle Filtering**: 8 advanced algorithms with adaptive selection
- **Multilooking**: ENL-adaptive spatial averaging
- **Python/CLI APIs**: Complete interface coverage

## Installation

### Development Installation

```bash
# Clone the repository
git clone https://github.com/SteveMHill/SARdine.git
cd SARdine

# Build and install
cargo build --release
pip install -e .
```

## 🚀 Quick Start

### Command Line Interface

```bash
# Complete processing pipeline
sardine info S1A_IW_SLC__1SDV_*.zip                    # Product information
sardine orbit S1A_IW_SLC__1SDV_*.zip                   # Download orbit files
sardine deburst S1A_IW_SLC__1SDV_*.zip --subswath IW1  # Deburst processing
sardine calibrate S1A_IW_SLC__1SDV_*.zip --type sigma0 # Radiometric calibration

# Advanced processing with DEM and speckle filtering
sardine dem-download --bbox "37.0,38.0,-122.5,-121.5" # Download DEM tiles
sardine speckle-filter data.npy --filter lee --looks 4 # Apply speckle filter
sardine estimate-nlooks data.npy                       # Estimate noise levels

# One-step terrain flattening (automatic DEM + orbit + calibration)
sardine terrain S1A_IW_SLC__1SDV_*.zip --polarization VV --range-looks 4 --azimuth-looks 1
```

### Python API

```python
import sardine

# Product information and orbit handling
info = sardine.get_product_info("S1A_IW_SLC__1SDV_*.zip")
orbit_data = sardine.download_orbit_file("S1A_IW_SLC__1SDV_*.zip", "./orbit_files/")

# Complete calibration workflow
reader = sardine.SlcReader("S1A_IW_SLC__1SDV_*.zip")
calibrated_data = reader.calibrate_and_multilook("VV", "sigma0", range_looks=3, azimuth_looks=3)

# Advanced speckle filtering
import numpy as np
filtered_data = sardine.apply_speckle_filter(calibrated_data, "enhanced_lee", num_looks=4)
noise_level = sardine.estimate_num_looks(calibrated_data)

# Complete terrain flattening pipeline (DEM + orbit + calibration)
reader.set_orbit_data(orbit_data)
gamma0_data, mask = reader.calibrate_multilook_and_flatten_auto_dem(
    "VV", "sigma0", range_looks=4, azimuth_looks=1, dem_cache_dir="./dem_cache"
)

# DEM download and processing
from sardine.types import BoundingBox
bbox = BoundingBox(min_lat=37.0, max_lat=38.0, min_lon=-122.5, max_lon=-121.5)
dem_files = sardine.download_srtm_tiles(bbox, "./dem_cache")
```

## 📊 Performance & Quality

- **Processing Speed**: >1.5M pixels/second for speckle filtering
- **Memory Efficiency**: In-place processing for large datasets  
- **Noise Reduction**: 1.5x-3x improvement typical with speckle filtering
- **DEM Coverage**: Global SRTM and Copernicus DEM support
- **Accuracy**: Maintains scientific accuracy while optimizing performance

## 📁 Repository Structure

```
SARdine/
├── src/                    # Rust source code
│   ├── core/              # Core processing algorithms
│   │   ├── calibrate.rs   # Radiometric calibration
│   │   ├── deburst.rs     # Burst concatenation
│   │   ├── speckle_filter.rs # Speckle filtering algorithms
│   │   └── multilook.rs   # Spatial averaging
│   ├── io/                # Input/output handling
│   │   ├── dem.rs         # DEM download and processing
│   │   ├── orbit.rs       # Orbit file management
│   │   └── slc_reader.rs  # SLC data reading
│   └── lib.rs             # Python bindings
├── python/sardine/        # Python API
├── examples/              # Usage examples
├── docs/                  # Documentation
│   └── implementation/    # Technical implementation guides
├── tests/                 # Test suite
└── data/                  # Test data directory
```

## 📚 Documentation & Examples

- **[docs/](docs/)** - Complete technical documentation
- **[examples/](examples/)** - Working code examples for all features  
- **[data/](data/)** - Test data setup instructions
- **API Documentation**: Generate with `cargo doc --open`

## 🧪 Testing

```bash
# Rust unit tests
cargo test

# Python integration tests  
python -m pytest tests/

# Run example workflows
python examples/complete_speckle_filtering_workflow.py
```

## 🤝 Contributing

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/amazing-feature`)
3. Make your changes and add tests
4. Update documentation if needed
5. Submit a pull request

See [docs/](docs/) for implementation details and [examples/](examples/) for usage patterns.

## 📄 License

MIT License - see LICENSE file for details.

## 🙏 Acknowledgments

SARdine provides a modern, high-performance alternative to traditional SAR processing tools while maintaining scientific accuracy. Developed for the SAR remote sensing community.

---

**🌟 Star this repository if SARdine helps your research!**

*For questions, issues, or contributions: [GitHub Issues](https://github.com/SteveMHill/SARdine/issues)*
