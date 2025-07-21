# Changelog

All notable changes to SARdine will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.2.0] - 2025-07-21

### Major Changes
- **🎯 Removed ALL synthetic data and fallback mechanisms** - SARdine now uses only real data
- **🏔️ Real DEM integration** - Uses only real SRTM tiles, no synthetic DEMs
- **🛰️ Real orbit data** - Extracts actual orbit state vectors from SLC files
- **📡 Complete IW processing** - All 6 IW subswaths (VV+VH IW1,IW2,IW3) properly extracted and processed

### Added
- Complete terrain correction pipeline with real DEM data
- Automatic SRTM tile download and mosaicking
- Masking and quality assessment workflows
- Geocoded product export (Cloud-Optimized GeoTIFF)
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
