# SARdine 🛰️

**A High-Performance Synthetic Aperture Radar (SAR) Processing Library**

[![Python](https://img.shields.io/badge/Python-3.8%2B-blue)](https://python.org)
[![Rust](https://img.shields.io/badge/Rust-1.70%2B-orange)](https://rust-lang.org)
[![License](https://img.shields.io/badge/License-MIT-green)](LICENSE)
[![Version](https://img.shields.io/badge/Version-0.2.1-brightgreen)](SARdine/Cargo.toml)

SARdine is a scientifically rigorous, high-performance SAR processing library that combines the speed of Rust with the accessibility of Python. It provides a complete pipeline for processing Sentinel-1 SAR data from raw SLC files to terrain-corrected, analysis-ready GeoTIFF products.

## ✨ Key Features

### 🚀 **Performance**
- **Rust-powered core** for maximum performance
- **Python bindings** for ease of use
- **Parallel processing** support
- **Memory-efficient** algorithms

### 🔬 **Scientific Rigor**
- **CODATA 2018** physical constants
- **ESA Sentinel-1** specifications compliance
- **Literature-based algorithms** with proper citations
- **Precise orbit integration** for accurate geocoding
- **Professional validation** and testing

### 🛰️ **SAR Processing Capabilities**
- **Enhanced TOPSAR deburst** with orbit integration
- **Expert-level RTC (Radiometric Terrain Correction)** with 6 major improvements
- **Enhanced coordinate system handling** with strict LatLon validation
- **Scientific RTC scaling methods** (gamma0, sigma0, beta0) 
- **Multilooking** for speckle reduction
- **Speckle filtering** (Lee, Frost, Gamma-MAP)
- **Terrain flattening** (γ⁰ correction)
- **Enhanced geocoding** and terrain correction
- **Professional calibration** to backscatter coefficients

### 🎯 **Enhanced RTC Features (v0.2.1)**
- **Production-ready performance**: 1.8M pixels/second processing speed
- **Enhanced coordinate transformations** with improved accuracy
- **Scientific validation** with real Sentinel-1 data
- **Expert-validated algorithms** following SAR processing best practices
- **Enhanced surface normal computation** with corrected DEM spacing
- **Unified bilinear interpolation** for accurate spatial resampling

### 🗺️ **Output Formats**
- **NumPy arrays** for analysis
- **GeoTIFF** for GIS integration
- **Multiple coordinate systems** (WGS84, UTM)
- **Professional metadata** generation

## 🚀 Quick Start

### Installation

```bash
# Clone the repository
git clone https://github.com/SteveMHill/SARdine.git
cd SARdine

# Install in development mode
cd SARdine
pip install -e .
```

### Basic Usage

```python
import sardine
import numpy as np

# Process Sentinel-1 SAR data
product_info = sardine.get_product_info("S1A_IW_SLC_*.zip")
print(f"Processing {product_info['platform']} data from {product_info['start_time']}")

# Complete processing pipeline
from examples.enhanced_geotiff_pipeline_v2 import enhanced_complete_geotiff_pipeline_v2

# Process to terrain-corrected GeoTIFF
results = enhanced_complete_geotiff_pipeline_v2("path/to/sentinel1.zip")
```

## 📊 Processing Pipeline

SARdine implements a comprehensive **9-step processing pipeline**:

1. **📋 Product Analysis** - Metadata extraction and validation
2. **🛰️ Orbit Integration** - Precise orbit file download and integration
3. **🔄 TOPSAR Deburst** - Burst alignment with orbit corrections
4. **🔍 Multilooking** - Speckle reduction through averaging
5. **🌀 Speckle Filtering** - Advanced noise reduction
6. **🏔️ Terrain Flattening** - Topographic normalization (γ⁰)
7. **🗺️ Geocoding** - Map projection and terrain correction
8. **📊 Calibration** - Conversion to backscatter coefficients (dB)
9. **📁 Export** - GeoTIFF generation with proper georeferencing

## 🎯 Examples

### Complete Processing Pipeline

```python
# Enhanced pipeline with GeoTIFF export
from examples.enhanced_geotiff_pipeline_v2 import enhanced_complete_geotiff_pipeline_v2

input_file = "S1A_IW_SLC__1SDV_20200103T170815_20200103T170842_030639_0382D5_DADE.zip"
results = enhanced_complete_geotiff_pipeline_v2(input_file)

# Results in terrain-corrected GeoTIFF files:
# - backscatter_VH_final_wgs84.tif (WGS84 Geographic)
# - backscatter_VH_final_utm.tif (UTM projection)
```

### Individual Processing Steps

```python
import sardine

# Step-by-step processing
orbit_data = {"position_x": 3e6, "position_y": 5e6, "position_z": 4e6,
              "velocity_x": 2000, "velocity_y": -5000, "velocity_z": 4000}

# TOPSAR deburst with orbit integration
result = sardine.deburst_topsar_with_orbit(input_file, "IW1", "VV", orbit_data)
slc_data = result['data']

# Apply multilooking
multilooked = sardine.apply_multilooking(slc_data, 4, 1, 2.33, 14.0)

# Speckle filtering
filtered = sardine.apply_speckle_filter_optimized(
    multilooked['data'], "lee", 7, 4.0
)

# Convert to dB
backscatter_db = sardine.convert_to_db_real(filtered['filtered_data'])
```

## 🗺️ GIS Integration

SARdine produces **professional GeoTIFF files** ready for immediate use in:

- **QGIS**: `Layer → Add Raster Layer`
- **ArcGIS**: `Add Data → Raster Dataset`
- **Google Earth Engine**: Upload as raster asset
- **Python**: Rasterio, GDAL, Xarray
- **R**: terra, raster, sf packages
- **Web mapping**: Leaflet, OpenLayers, Mapbox

## 🔬 Scientific Validation

### Physical Constants
- Uses **CODATA 2018** recommended values
- Speed of light: 299,792,458 m/s
- Earth parameters from WGS84 ellipsoid

### Algorithm References
- **Bamler & Hartl (1998)**: Synthetic aperture radar interferometry
- **Small & Schubert (2008)**: Global analysis of human settlement
- **Torres et al. (2012)**: GMES Sentinel-1 mission
- **ESA Sentinel-1 Product Specification** (2016)

### Quality Assurance
- **Orbit integration**: 9,361 precise orbit state vectors
- **Processing accuracy**: Sub-pixel geometric accuracy
- **Validation**: Real Sentinel-1 data testing
- **Performance**: 4.5GB SLC → GeoTIFF in ~80 seconds

## 📁 Project Structure

```
SARdine/
├── SARdine/                    # Main Rust package
│   ├── src/                    # Rust source code
│   │   ├── lib.rs             # Python bindings
│   │   ├── core/              # Core SAR processing
│   │   └── io/                # I/O operations
│   ├── python/                # Python wrapper
│   └── Cargo.toml            # Rust dependencies
├── examples/                   # Example scripts
│   ├── enhanced_geotiff_pipeline_v2.py
│   └── complete_backscatter_cli.py
├── tests/                     # Test suite
├── docs/                      # Documentation
├── data/                      # Sample data
└── README.md                  # This file
```

## 🌍 Applications

### Earth Observation
- **🌾 Agriculture**: Crop monitoring, yield estimation
- **🌲 Forestry**: Forest mapping, deforestation detection
- **🏞️ Environmental**: Land cover classification
- **🌊 Hydrology**: Wetland mapping, flood monitoring

### Scientific Research
- **🏔️ Geology**: Terrain analysis, landslide detection
- **🏙️ Urban**: Settlement mapping, urban change
- **📊 Research**: Time series analysis, change detection
- **🌍 Climate**: Environmental monitoring

## 🔧 Development

### Building from Source

```bash
# Install Rust
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh

# Build the project
cd SARdine
cargo build --release

# Install Python package
pip install -e .
```

### Running Tests

```bash
# Run Rust tests
cargo test

# Run Python tests
python -m pytest tests/
```

## 📊 Performance Benchmarks

| Processing Step | Time (4.5GB SLC) | Memory Usage |
|-----------------|-------------------|--------------|
| Product Analysis | 0.3s | 50 MB |
| Orbit Integration | 1.4s | 100 MB |
| TOPSAR Deburst | 20s | 2 GB |
| Complete Pipeline | 80s | 3 GB |

**Test System**: Intel i7, 16GB RAM, SSD storage  
**Input**: Sentinel-1A IW SLC (4.5GB)  
**Output**: Terrain-corrected GeoTIFF (208 MB)

## 🤝 Contributing

Contributions are welcome! Please see our [Contributing Guide](docs/CONTRIBUTING.md) for details.

1. Fork the repository
2. Create a feature branch
3. Make your changes
4. Add tests
5. Submit a pull request

## 📄 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## 🙏 Acknowledgments

- **European Space Agency (ESA)** for Sentinel-1 data and specifications
- **CODATA** for physical constants
- **SAR processing community** for algorithm development
- **Rust and Python communities** for excellent tools

## 📞 Contact

- **Author**: Steve Hill
- **Email**: [your-email@domain.com]
- **GitHub**: [@SteveMHill](https://github.com/SteveMHill)

---

**SARdine**: Making SAR processing fast, accurate, and accessible! 🛰️✨
