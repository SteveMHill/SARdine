# SARdine: Enhanced Sentinel-1 SAR Processor

A fast, modular, and scientifically accurate Sentinel-1 backscatter processor with enhanced capabilities:

## Key Features

- **Enhanced SLC Reader**: Dual ZIP/SAFE format support with real metadata extraction
- **Real Parameter Processing**: XML-based parameter parsing replacing hardcoded values
- **TOPSAR Processing**: Fixed deburst implementation using real burst parameters
- **Scientific Accuracy**: Research-grade SAR processor suitable for scientific applications
- **Python Integration**: Complete PyO3 bindings for 14-step SAR processing pipeline

## Recent Major Improvements

### Enhanced SLC Reader (August 2025)
- ZIP and SAFE format support with automatic detection
- Real metadata extraction from annotation XML files
- Universal file access methods for both formats
- Comprehensive parameter extraction (pixel spacing, timing, geographic bounds)

### Scientific Accuracy Improvements
- Real burst parameters from annotation XML (not hardcoded)
- XML-based calibration coefficients
- Actual pixel spacing from product metadata
- Real orbit data integration
- Proper coordinate system handling

### Processing Pipeline
Complete 14-step SAR processing pipeline:
1. SLC Data Reading with real metadata
2. Precise Orbit File application
3. IW Split using real geometry
4. TOPSAR Deburst with real parameters
5. Radiometric Calibration
6. IW Subswath Merging
7. Multilooking
8. Terrain Flattening
9. Speckle Filtering
10. Terrain Correction
11. Advanced Masking
12. dB Conversion
13. GeoTIFF Export
14. Quality Assessment

## Installation

```bash
# Build from source
maturin develop

# Or build wheel
maturin build
```

## Testing

Run comprehensive tests:
```bash
python test_enhanced_sardine.py
```

## License

MIT License - See LICENSE file for details.
