# SARdine Documentation

Welcome to the SARdine documentation! This directory contains comprehensive guides for understanding and extending SARdine's SAR processing capabilities.

## Implementation Guides

The `implementation/` directory contains detailed technical documentation for each major component:

### Core Processing Components
- **[Deburst Implementation](implementation/DEBURST_IMPLEMENTATION.md)** - Burst boundary removal and seamless merging
- **[Calibration & Orbit](implementation/ORBIT_IMPLEMENTATION.md)** - Radiometric calibration and precise orbit application
- **[Multilook Implementation](implementation/MULTILOOK_IMPLEMENTATION.md)** - Spatial averaging for speckle reduction

### Advanced Processing
- **[Speckle Filtering](implementation/SPECKLE_FILTERING_IMPLEMENTATION.md)** - Advanced speckle noise reduction algorithms
- **[Terrain Flattening](implementation/TERRAIN_FLATTENING_COMPLETE.md)** - Topographic normalization and DEM processing
- **[DEM Implementation](implementation/DEM_IMPLEMENTATION_SUMMARY.md)** - Automatic DEM acquisition and processing

## Quick Start

For getting started with SARdine, see the main [README.md](../README.md) in the root directory.

## Examples

Practical usage examples can be found in the [examples/](../examples/) directory.

## API Reference

- **Rust API**: Generated automatically with `cargo doc --open`
- **Python API**: See docstrings in the [python/sardine/](../python/sardine/) module

## Contributing

When adding new features, please:
1. Add implementation documentation to `implementation/`
2. Include working examples in `examples/`
3. Update the main README if the public API changes
4. Add tests to verify functionality

---

*For questions or contributions, please see the main repository at [GitHub](https://github.com/SteveMHill/SARdine).*
