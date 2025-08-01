# SARdine - Experimental SAR Processing Package

<div align="center">
  <h2>🛰️ SARdine</h2>
  <p><em>Experimental SAR Processing Toolkit</em></p>
</div>

## ⚠️ **EXPERIMENTAL - NOT PRODUCTION READY** ⚠️

**This package is experimental and designed for research, learning, and testing purposes only. It is NOT production-ready and should not be used in operational environments.**

**Status**: Development/Research Tool  
**Stability**: Experimental  
**Support**: Community/Educational Use Only

---

A streamlined SAR processing package with essential tools for Sentinel-1 data processing, designed for learning and experimentation with SAR processing algorithms.

## 🎯 Quick Start

### Process Sentinel-1 backscatter with one command:

```bash
python backscatter_cli.py input_S1A_*.zip ./output/
```

### Basic Usage Examples:

```bash
# Default VV processing
python backscatter_cli.py S1A_IW_SLC_*.zip ./my_output/

# VH polarization with custom speckle filter
python backscatter_cli.py S1A_*.zip ./output/ --polarization VH --speckle-filter lee

# High-resolution processing
python backscatter_cli.py S1A_*.zip ./output/ --resolution 5 --filter-window 5

# Quick processing (skip geocoding)
python backscatter_cli.py S1A_*.zip ./output/ --no-geocode --no-terrain-flatten
```

## 📁 Package Structure

```
SARdine/
├── backscatter_cli.py          # Main CLI for SAR processing
├── examples.py                 # Usage examples and quick tests
├── test_14_step_pipeline.py    # Complete processing pipeline test
├── test_real_data_only.py      # Real Sentinel-1 data test
├── test_direct_dem.py          # DEM processing test
├── data/                       # Input data directory
├── docs/                       # Documentation
│   ├── USER_GUIDE.md
│   ├── SCIENTIFIC_AUDIT_REPORT.md
│   └── COMPLETE_14_STEP_PIPELINE_SUCCESS_REPORT.md
├── SARdine/                    # Core library (Rust + Python)
└── archive/                    # Legacy files (for reference)
```

## 🧪 Testing

Run the essential tests to verify functionality:

```bash
# Test with real Sentinel-1 data (most comprehensive)
python test_real_data_only.py

# Test complete 14-step pipeline
python test_14_step_pipeline.py

# Test DEM processing
python test_direct_dem.py

# Quick examples and functionality check
python examples.py
```

## 📊 Processing Steps

The pipeline implements these key SAR processing steps:

1. **SLC Reading** - Load Sentinel-1 SLC data
2. **Orbit Correction** - Apply precise orbit files
3. **IW Split** - Split interferometric wide swath data
4. **Deburst** - Remove TOPSAR burst boundaries
5. **Radiometric Calibration** - Convert to σ⁰ backscatter
6. **Multilooking** - Reduce speckle noise
7. **Speckle Filtering** - Advanced noise reduction
8. **Terrain Correction** - Geometric correction using DEM
9. **Geocoding** - Project to geographic coordinates

## 📖 Documentation

- **User Guide**: `docs/USER_GUIDE.md` - Complete usage instructions
- **Scientific Report**: `docs/SCIENTIFIC_AUDIT_REPORT.md` - Algorithm validation
- **Pipeline Report**: `docs/COMPLETE_14_STEP_PIPELINE_SUCCESS_REPORT.md` - Implementation details

## 🔧 Requirements

- Python 3.8+
- Rust (for core processing)
- GDAL
- Sentinel-1 SLC data

## ⚠️ Important Notes

- **Not for production use** - This is an experimental research tool
- **Limited testing** - Primarily tested with specific Sentinel-1 datasets
- **Active development** - API and functionality may change
- **No warranty** - Use at your own risk for research purposes

## 🤝 Contributing

This is a research/learning project. Contributions welcome for:
- Bug fixes
- Algorithm improvements
- Documentation
- Additional test cases

## 📄 License

MIT License - See SARdine/LICENSE for details.

## 🏷️ Version

Experimental/Development - Not versioned for production use.

---

**Remember**: This package is for learning and experimentation only. For production SAR processing, use established tools like SNAP, GAMMA, or commercial solutions.
