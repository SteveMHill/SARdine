# SARdine: A Fast, Modular Sentinel-1 SAR Data Processor

SARdine is a modern, high-performance SAR data processing library for Sentinel-1 SLC products, implemented in Rust with Python bindings. It provides fast, reliable processing of complex SAR workflows including orbit correction, IW split, deburst, and radiometric calibration.

## Features

- 🚀 **Fast Rust Backend**: High-performance core processing in Rust
- 🐍 **Python API**: Easy-to-use Python interface with full feature access
- 📱 **CLI Interface**: Command-line tools for batch processing
- 🛰️ **Orbit Handling**: Automatic orbit file download and precise state vector interpolation
- 📡 **IW Split**: Extract individual sub-swaths from IW SLC products
- 🔗 **Deburst**: Seamless burst concatenation with precise timing
- 📊 **Radiometric Calibration**: Sigma0/Beta0/Gamma0 calibration with bilinear interpolation
- 🧪 **Comprehensive Testing**: Rust unit tests and Python integration tests

## Current Implementation Status

### ✅ Completed Components
- **Orbit File Management**: Download, validation, and interpolation
- **SLC Reader**: ZIP archive handling and metadata extraction
- **IW Split**: Sub-swath extraction from SLC products
- **Deburst**: Burst concatenation with precise geolocation
- **Calibration Vector Parsing**: Robust XML parsing and bilinear interpolation
- **Python API**: Complete exposure of Rust functionality
- **CLI Interface**: Command-line access to all processing steps

### 🔄 In Progress
- End-to-end workflow integration
- Performance optimization
- Error handling improvements

## Installation

### Development Installation

```bash
# Clone the repository
git clone https://github.com/yourusername/sardine.git
cd sardine

# Build Rust components
cargo build --release

# Install Python package in development mode
pip install -e .
```

## Usage

### Command Line Interface

```bash
# Get product information
sardine info path/to/S1A_IW_SLC__1SDV_*.zip

# Download orbit files
sardine orbit path/to/S1A_IW_SLC__1SDV_*.zip --output ./orbit_files/

# Extract sub-swath
sardine iw-split path/to/S1A_IW_SLC__1SDV_*.zip --subswath IW1 --output ./iw1_data/

# Deburst processing
sardine deburst path/to/S1A_IW_SLC__1SDV_*.zip --subswath IW1 --output ./deburst_data/

# Radiometric calibration
sardine calibrate path/to/S1A_IW_SLC__1SDV_*.zip --subswath IW1 --calibration-type sigma0 --output ./calibrated_data/
```

### Python API

```python
import sardine

# Product information
info = sardine.get_product_info("path/to/S1A_IW_SLC__1SDV_*.zip")
print(f"Product: {info['product_id']}")
print(f"Acquisition: {info['start_time']} - {info['stop_time']}")

# Orbit processing
orbit_data = sardine.download_orbit_file("path/to/S1A_IW_SLC__1SDV_*.zip", "./orbit_files/")

# IW split
iw_data = sardine.split_iw_subswath("path/to/S1A_IW_SLC__1SDV_*.zip", "IW1")

# Deburst
deburst_data = sardine.deburst_subswath("path/to/S1A_IW_SLC__1SDV_*.zip", "IW1")

# Calibration
calibrated_data = sardine.calibrate_data("path/to/S1A_IW_SLC__1SDV_*.zip", "IW1", "sigma0")
```

## Examples

See the `examples/` directory for complete workflow demonstrations:

- `complete_orbit_workflow.py` - Orbit file handling
- `complete_iw_split_workflow.py` - Sub-swath extraction
- `complete_deburst_workflow.py` - Deburst processing
- `complete_calibration_workflow.py` - Radiometric calibration
- `python_orbit_api.py` - Python API usage

## Architecture

SARdine follows a modular architecture:

```
src/
├── lib.rs              # Python bindings (PyO3)
├── types.rs            # Core data structures
├── core/
│   ├── calibrate.rs    # Radiometric calibration
│   └── deburst.rs      # Deburst processing
└── io/
    ├── annotation.rs   # XML annotation parsing
    ├── orbit.rs        # Orbit file handling
    ├── slc_reader.rs   # SLC data reading
    └── dem.rs          # DEM integration
```

## Dependencies

### Rust Dependencies
- `gdal` - Geospatial data processing
- `pyo3` - Python bindings
- `zip` - Archive handling
- `quick-xml` - XML parsing
- `reqwest` - HTTP requests
- `chrono` - Date/time handling

### Python Dependencies
- `numpy` - Numerical computing
- `rasterio` - Raster I/O
- `click` - CLI framework

## Testing

### Rust Tests
```bash
cargo test
```

### Python Tests
```bash
python -m pytest tests/
```

## Contributing

1. Fork the repository
2. Create a feature branch
3. Make your changes
4. Add tests for new functionality
5. Submit a pull request

## Documentation

- `ORBIT_IMPLEMENTATION.md` - Detailed orbit handling documentation
- `DEBURST_IMPLEMENTATION.md` - Deburst algorithm documentation
- `sar_pipeline_readme.md` - Pipeline overview and roadmap

## License

MIT License - see LICENSE file for details.

## Acknowledgments

This project was developed to provide a fast, modern alternative to traditional SAR processing tools while maintaining scientific accuracy and reliability.
