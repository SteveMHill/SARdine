# SARdine: A Fast, Modular Sentinel-1 Backscatter Processor

SARdine is a modern, open-source alternative to ESA SNAP and GAMMA for processing Sentinel-1 SLC data into calibrated, terrain-corrected backscatter products.

## Features

- 🚀 Fast Rust backend with Python API
- 📊 High-quality Gamma0 VV/VH GeoTIFFs
- 🔧 Modular, extensible pipeline
- 📱 Simple CLI interface
- 🧪 Comprehensive test suite

## Installation

```bash
pip install sardine
```

## Quick Start

```bash
# Get product information
sardine info path/to/sentinel1.zip

# Process SLC data
sardine process path/to/sentinel1.zip --output ./results/
```

## Python API

```python
import sardine

# Get product information
info = sardine.get_product_info("path/to/sentinel1.zip")
print(info)

# Process SLC data (coming soon)
result = sardine.process_slc("path/to/sentinel1.zip", output_dir="./results/")
```

## Development

This project is under active development. See the test data and examples for current capabilities.

## License

MIT License - see LICENSE file for details.
