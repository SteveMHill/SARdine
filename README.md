# SARdine: High-Performance SAR Processing Library

<div align="center">
  <img src="logo.png" alt="SARdine Logo" width="300"/>
  
  [![CI](https://img.shields.io/github/workflow/status/your-username/SARdine/CI)](https://github.com/your-username/SARdine/actions)
  [![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
  [![Rust](https://img.shields.io/badge/rust-1.70+-blue.svg)](https://www.rust-lang.org)
  [![Python](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org)
</div>

## Overview

SARdine is a modern, high-performance SAR data processing library for Sentinel-1 data, implemented in Rust with Python bindings. It provides a complete processing pipeline from SLC to analysis-ready products, optimized for both research and operational use.

## 🌟 Key Features

- 🚀 **High-Performance**: Rust backend with parallel processing and memory optimization
- 🐍 **Python Integration**: Easy-to-use Python API with numpy compatibility
- 📡 **Complete Pipeline**: SLC to backscatter products in one library
- 🛰️ **Orbit Handling**: Automatic precise orbit file download and application
- 🏔️ **Terrain Processing**: DEM integration with terrain correction and flattening
- 📊 **Radiometric Calibration**: σ⁰, β⁰, γ⁰ calibration with quality control
- 🎯 **Speckle Filtering**: Advanced filtering algorithms (Lee, Enhanced Lee, Gamma MAP)
- ⚡ **Memory Efficient**: Chunked processing for large datasets

## 🎯 Core Processing Components

### ✅ Implemented
- **SLC Processing**: Reading, debursting, and concatenation
- **Orbit Correction**: Precise orbit file application
- **Radiometric Calibration**: Multi-polarization calibration to backscatter
- **Terrain Processing**: DEM-based terrain correction and flattening  
- **Speckle Filtering**: Multiple adaptive filters with quality metrics
- **Geocoding**: Map projection and coordinate transformation
- **Export**: GeoTIFF and other format support

### Performance Highlights
- **Terrain Flattening**: 9M+ pixels/second average throughput
- **Speckle Filtering**: Parallel processing with coefficient of variation improvements
- **Memory Usage**: Chunked processing prevents memory overflow
- **Scaling**: Multi-core CPU utilization with Rayon

## 🚀 Installation

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
git clone https://github.com/your-username/SARdine.git
cd SARdine

# Install Python dependencies
pip install maturin numpy rasterio pyproj

# Build and install
maturin develop --release

# Verify installation
python -c "import sardine; print('SARdine installed successfully')"
```

## 🐍 Python API Quick Start

```python
import sardine
import numpy as np

# Basic terrain flattening
sigma0 = np.array(...)  # Your SAR backscatter data
dem = np.array(...)     # Digital elevation model

# Apply terrain flattening
result = sardine.apply_terrain_flattening(sigma0, dem)
gamma0 = result['gamma0']  # Terrain-flattened backscatter
incidence_angles = result['incidence_angles']

# Speckle filtering
filtered = sardine.apply_speckle_filter_numpy(
    sigma0, "enhanced_lee", window_size=7, num_looks=1.0
)
```

## 📊 Performance

SARdine is optimized for high-throughput processing:

- **Terrain Flattening**: 9M+ pixels/second average
- **Speckle Filtering**: Multi-core parallel processing
- **Memory Efficient**: Chunked processing for large datasets
- **Scalable**: Automatic CPU core utilization

## 📁 Examples

See the `examples/` directory for comprehensive usage examples:

- `complete_backscatter_pipeline.py` - Full SLC to backscatter processing
- `complete_terrain_correction_workflow.py` - Terrain correction with DEM
- `complete_speckle_filtering_workflow.py` - Advanced speckle filtering
- `production_backscatter_processor.py` - Production-ready processing

## 🧪 Testing

```bash
# Run basic tests
python tests/test_basic.py

# Run comprehensive tests
python tests/test_pipeline.py
python tests/test_terrain_flattening.py

# Run Rust tests
cargo test --release
```

## 🤝 Contributing

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/amazing-feature`)
3. Commit your changes (`git commit -m 'Add amazing feature'`)
4. Push to the branch (`git push origin feature/amazing-feature`)
5. Open a Pull Request

## 📄 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## 🙏 Acknowledgments

- ESA for Sentinel-1 data and documentation
- The Rust and Python communities for excellent tooling
- GDAL/OGR for geospatial data handling

---

**SARdine**: High-performance SAR processing made simple 🚀
