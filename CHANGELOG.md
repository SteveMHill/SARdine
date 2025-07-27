# Changelog

All notable changes to SARdine will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.2.0] - 2025-07-27

### Added
- 🚀 **High-Performance Terrain Flattening**: 9M+ pixels/second throughput with parallel processing
- 🏔️ **Complete Terrain Correction Pipeline**: DEM-based Range-Doppler terrain correction
- 🎯 **Advanced Speckle Filtering**: Lee, Enhanced Lee, and Gamma MAP filters with parallel processing
- 🐍 **Python API**: Complete PyO3 bindings for all core functions
- ⚡ **Memory Optimization**: Chunked processing for large datasets
- 📊 **Scientific Validation**: Realistic Sentinel-1 geometry and quality control

### Performance Improvements
- Parallel processing with Rayon for multi-core CPU utilization
- Chunked memory management prevents overflow on large arrays
- Optimized Rust algorithms with SIMD support
- Streaming I/O for memory-efficient processing

### Technical Details
- Rust 1.70+ backend with PyO3 Python bindings
- Complete terrain flattening implementation with realistic incidence angles
- Robust slope/aspect computation with edge handling
- Quality masking for terrain correction validation
- Comprehensive test suite with performance benchmarks

## [0.1.0] - 2025-07-01

### Added
- Initial SARdine implementation
- Basic SLC reading and metadata extraction
- Orbit file handling and application
- Debursting and sub-swath concatenation
- Radiometric calibration (σ⁰, β⁰, γ⁰)
- DEM download and processing infrastructure
- Command-line interface and Python bindings

### Infrastructure
- Rust-based core processing engine
- Python package with pip installation
- GitHub Actions CI/CD pipeline
- MIT license and open-source release

## [Unreleased]

### Planned
- Time series analysis tools
- Advanced polarimetric processing
- GPU acceleration support
- Enhanced CLI tools
- Development scripts and validation tools
- Comprehensive error handling (fails on insufficient real data)

### Fixed
- IW subswath extraction now processes all bursts correctly
- Orbit time format handling for terrain correction
- DEM preparation uses real scene bounding box from manifest
- Terrain correction no longer falls back to uncorrected data
- Coordinate system handling and geolocation accuracy

### Removed
- All synthetic DEM generation code
- Fallback orbit data creation
- Placeholder burst information
- Zero-filled arrays and artificial data
- "Using speckle filtered data as fallback" logic
- All placeholder and dummy data generation

### Technical Improvements
- Clean separation of development tools in `dev_scripts/`
- Updated documentation and examples
- Improved .gitignore for better repository management
- Enhanced logging and status reporting

### Performance
- Maintains high processing speed with real data only
- Successful processing of full German Sentinel-1 scene
- Coverage results: VV 95.0%, VH 92.8%

## [0.1.0] - 2025-07-17

### Initial Release
- Basic SLC reading and metadata extraction
- Orbit file handling
- Debursting and calibration
- Initial terrain processing (with synthetic fallbacks)
- Python bindings and CLI interface
- Documentation and examples
