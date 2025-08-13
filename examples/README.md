# SARdine Examples 📚

This directory contains comprehensive examples demonstrating SARdine's SAR processing capabilities.

## 🚀 Quick Start Examples

### 1. **Enhanced GeoTIFF Pipeline** (Recommended)
**File**: `enhanced_geotiff_pipeline_v2.py`

Complete processing pipeline with terrain correction and GeoTIFF export.

```bash
python enhanced_geotiff_pipeline_v2.py
```

**Features**:
- ✅ 9-step processing pipeline
- ✅ Terrain flattening (γ⁰ correction)
- ✅ Dual coordinate systems (WGS84 + UTM)
- ✅ Professional GeoTIFF export
- ✅ ~80 seconds processing time

**Output**:
- `backscatter_VH_final_wgs84.tif` (208 MB)
- `backscatter_VH_final_utm.tif` (208 MB)
- Processing logs and metadata

### 2. **Complete Backscatter CLI**
**File**: `complete_backscatter_cli.py`

Command-line interface for basic SAR processing.

```bash
python complete_backscatter_cli.py
```

**Features**:
- ✅ 7-step processing pipeline
- ✅ NumPy array outputs
- ✅ Scientific validation
- ✅ Performance benchmarking

### 3. **Enhanced Complete Pipeline**
**File**: `enhanced_complete_pipeline.py`

Advanced processing with terrain corrections.

```bash
python enhanced_complete_pipeline.py
```

**Features**:
- ✅ 10-step enhanced pipeline
- ✅ Terrain flattening
- ✅ Geocoding capabilities
- ✅ Comprehensive metadata

## 📊 Processing Results Analysis

### **GeoTIFF Results Summary**
**File**: `geotiff_results_summary.py`

Comprehensive analysis of processing results.

```bash
python geotiff_results_summary.py
```

**Provides**:
- 📈 Data quality metrics
- 🗺️ Geographic coverage analysis
- 📊 Backscatter distribution statistics
- 🎯 Application recommendations

### **Enhanced Pipeline Summary**  
**File**: `enhanced_pipeline_summary.py`

Detailed summary of enhanced processing achievements.

```bash
python enhanced_pipeline_summary.py
```

## 🛰️ Example Data Requirements

All examples expect Sentinel-1 SLC data in the `../data/` directory:

```
data/
└── S1A_IW_SLC__1SDV_20200103T170815_20200103T170842_030639_0382D5_DADE.zip
```

## 📁 Output Structure

Examples generate outputs in `../complete_output/`:

```
complete_output/
├── Raw Processing Data
│   └── slc_debursted_VH_IW1.npy              (1.0 GB)
├── Intermediate Products  
│   ├── multilooked_VH_4x1.npy                (259 MB)
│   ├── filtered_VH_lee.npy                   (259 MB)
│   └── flattened_VH_gamma0.npy               (259 MB)
├── Final Products
│   └── backscatter_VH_final.npy              (259 MB)
├── GeoTIFF Exports
│   ├── backscatter_VH_final_wgs84.tif        (208 MB)
│   ├── backscatter_VH_final_utm.tif          (208 MB)
│   ├── filtered_VH_lee_wgs84.tif             (210 MB)
│   └── flattened_VH_gamma0_wgs84.tif         (208 MB)
├── Metadata
│   ├── processing_log_enhanced_complete.json
│   └── processing_log_enhanced_geotiff_v2.json
└── orbit_cache/
    └── [Precise orbit files]
```

## 🎯 Usage Patterns

### For GIS Users
```bash
# Generate GeoTIFF files for QGIS/ArcGIS
python enhanced_geotiff_pipeline_v2.py

# Import in QGIS:
# Layer → Add Raster Layer → Select:
# complete_output/backscatter_VH_final_wgs84.tif
```

### For Python Analysis
```python
import numpy as np
import matplotlib.pyplot as plt

# Load processed data
backscatter = np.load('../complete_output/backscatter_VH_final.npy')

# Visualize
plt.imshow(backscatter, cmap='gray', vmin=-20, vmax=0)
plt.title('SAR Backscatter (dB)')
plt.colorbar(label='Backscatter (dB)')
plt.show()
```

### For Research Applications
```python
# Load multiple processing levels
import rasterio

# Geographic data for global analysis
with rasterio.open('../complete_output/backscatter_VH_final_wgs84.tif') as src:
    data = src.read(1)
    bounds = src.bounds
    crs = src.crs

# UTM data for local analysis  
with rasterio.open('../complete_output/backscatter_VH_final_utm.tif') as src:
    utm_data = src.read(1)
    utm_transform = src.transform
```

## 🔬 Scientific Applications

### Land Cover Classification
```python
# Use terrain-corrected data for classification
gamma0_data = np.load('../complete_output/flattened_VH_gamma0.npy')

# Convert to dB for analysis
gamma0_db = 10 * np.log10(gamma0_data + 1e-10)

# Apply classification algorithms...
```

### Change Detection
```python
# Compare multiple dates using consistent processing
# All examples maintain scientific rigor for time series
```

### Environmental Monitoring
```python
# Use GeoTIFF outputs in environmental workflows
# Compatible with GDAL/OGR ecosystem
```

## ⚡ Performance Tips

1. **Memory Management**:
   - Process large scenes in tiles
   - Use multilooked data for preview

2. **Storage Optimization**:
   - GeoTIFF files use LZW compression
   - Remove intermediate files after processing

3. **Processing Speed**:
   - Examples leverage Rust performance
   - Parallel processing where possible

## 🐛 Troubleshooting

### Common Issues

**Memory Error**:
```bash
# Increase virtual memory or process smaller tiles
```

**Missing Data File**:
```bash
# Ensure Sentinel-1 SLC file is in ../data/ directory
# Download from: https://scihub.copernicus.eu/
```

**Import Error**:
```bash
# Ensure SARdine is installed:
cd ../SARdine && pip install -e .
```

## 📈 Next Steps

After running examples:

1. **Visualize Results**: Open GeoTIFF files in QGIS
2. **Analyze Data**: Use NumPy arrays for quantitative analysis
3. **Develop Applications**: Build on the processing pipeline
4. **Scale Processing**: Apply to multiple Sentinel-1 scenes

## 🤝 Contributing Examples

To contribute new examples:

1. Follow the existing naming convention
2. Include comprehensive documentation
3. Add error handling and validation
4. Test with real Sentinel-1 data
5. Submit a pull request

---

**Ready to process SAR data!** 🛰️ Start with `enhanced_geotiff_pipeline_v2.py` for the complete experience.
