# SARdine User Guide: Easy SAR Processing for Scientists
## Complete Guide to Research-Ready Backscatter Products

### 🚀 Quick Start (< 5 minutes)

**For users who want immediate results:**

```bash
# Download this repository
git clone https://github.com/SARdine/SARdine
cd SARdine

# Install dependencies
pip install -r requirements.txt

# Process your data with one command
python quick_sar.py your_sentinel1_file.zip

# Done! Check the sardine_output_* directory for results
```

### 📋 What You Get

SARdine produces **analysis-ready** SAR backscatter products that are:
- ✅ **Scientifically validated** - All algorithms verified against literature
- ✅ **Research-grade quality** - Suitable for peer-reviewed publications
- ✅ **Fully documented** - Complete metadata and provenance tracking
- ✅ **Standard formats** - GeoTIFF with complete georeferencing

### 🎯 Processing Options

#### Option 1: Quick Start (Recommended for beginners)
```bash
python quick_sar.py input.zip
```
**Automatically applies best-practice settings for research-quality results.**

#### Option 2: User-Friendly CLI (Recommended for most users)
```bash
python sardine_user_friendly_cli.py input.zip output/ --polarization VV --terrain-flatten --speckle-filter enhanced_lee
```

#### Option 3: Advanced Control (For experts)
```bash
python sardine_user_friendly_cli.py input.zip output/ \\
  --polarization VV \\
  --range-looks 2 \\
  --azimuth-looks 1 \\
  --terrain-flatten \\
  --speckle-filter refined_lee \\
  --filter-window 7 \\
  --geocode \\
  --resolution 10 \\
  --quality-report
```

### 🔬 Scientific Processing Steps

SARdine implements the complete **14-step SAR processing pipeline**:

1. **📖 Read SLC data** - Complex I+Q data from Sentinel-1
2. **🛰️ Orbit correction** - Precise orbit files for accurate positioning  
3. **📡 IW split** - Separate sub-swaths (IW1, IW2, IW3)
4. **💥 TOPSAR deburst** - Remove burst artifacts and azimuth deramp
5. **⚡ Radiometric calibration** - Convert to σ⁰ backscatter (ESA standard)
6. **🔗 IW merge** - Seamless wide-swath combination
7. **👁️ Multilooking** - Reduce speckle noise through spatial averaging
8. **🏔️ Terrain flattening** - Convert σ⁰ → γ⁰ for topographic normalization
9. **🌊 Speckle filtering** - Advanced noise reduction (6 algorithms available)
10. **🗺️ Terrain correction** - Geocoding to map coordinates (optional)
11. **🎭 Quality masking** - Multi-dimensional quality assessment
12. **📊 dB conversion** - Convert to logarithmic scale
13. **💾 Export products** - Standard GeoTIFF format
14. **📋 Generate metadata** - Complete provenance and quality reports

### 🌊 Speckle Filtering Options

Choose the best algorithm for your application:

- **`none`** - No filtering (fastest, preserves all details)
- **`lee`** - Classic Lee filter (good balance, widely used)
- **`enhanced_lee`** - **Recommended for most users** (better edge preservation)
- **`gamma_map`** - MAP estimation (good for heterogeneous scenes)
- **`lee_sigma`** - Statistical outlier rejection (robust)
- **`frost`** - Exponential damping (smooth results)
- **`refined_lee`** - Most advanced (best quality, slower)

### 🏔️ Terrain Processing

**Terrain flattening** is highly recommended for research:

```bash
--terrain-flatten  # Converts σ⁰ → γ⁰ for topographic normalization
```

**Benefits:**
- ✅ Removes topographic effects from backscatter
- ✅ Enables comparison across different terrain types
- ✅ Essential for vegetation and land cover analysis
- ✅ Uses automatic SRTM DEM (no manual downloads needed)

### 📊 Output Products

Each processing run produces:

```
output_directory/
├── sardine_VV_backscatter.tif    # Main backscatter product (dB)
├── processing_metadata.json      # Complete processing parameters
├── quality_report.txt           # Quality assessment summary
└── dem_cache/                   # Cached DEM tiles (reusable)
```

**Main Product Features:**
- 📊 **Format**: GeoTIFF with complete georeferencing
- 📏 **Units**: dB (decibels) for backscatter
- 🌍 **Coordinates**: WGS84 geographic (EPSG:4326) or UTM
- 📈 **Dynamic range**: Typically -30 to +10 dB
- ✅ **Quality masked**: Invalid pixels removed

### 💡 Usage Examples

#### Vegetation Monitoring
```bash
python sardine_user_friendly_cli.py input.zip vegetation_study/ \\
  --polarization VH \\
  --terrain-flatten \\
  --speckle-filter enhanced_lee \\
  --resolution 20
```

#### Urban Analysis
```bash
python sardine_user_friendly_cli.py input.zip urban_study/ \\
  --polarization VV \\
  --speckle-filter lee \\
  --geocode \\
  --resolution 10
```

#### Water Body Mapping
```bash
python sardine_user_friendly_cli.py input.zip water_study/ \\
  --polarization VV \\
  --terrain-flatten \\
  --speckle-filter gamma_map \\
  --quality-report
```

#### Change Detection Time Series
```bash
# Process multiple dates with identical parameters
for file in S1A_*.zip; do
  python sardine_user_friendly_cli.py "$file" "timeseries/$(basename $file .zip)/" \\
    --polarization VV \\
    --terrain-flatten \\
    --speckle-filter enhanced_lee \\
    --resolution 20
done
```

### 📚 Quality Assessment

SARdine provides comprehensive quality metrics:

**Automatic Quality Checks:**
- ✅ **Valid pixel percentage** - Data coverage assessment
- ✅ **SNR analysis** - Signal-to-noise ratio evaluation
- ✅ **Geometric accuracy** - Spatial positioning validation
- ✅ **Radiometric quality** - Calibration accuracy check
- ✅ **Statistical validation** - Data distribution analysis

**Quality Report Contents:**
```
📊 Valid pixels: 94.5%
📈 Overall quality score: 0.89
🎯 Geometric accuracy: <10m CE90
📡 Radiometric accuracy: ±0.8 dB
🌊 Speckle reduction: 94.2% MSE improvement
```

### 🔧 Troubleshooting

**Common Issues and Solutions:**

**"File not found" error:**
```bash
# Check file path and extension
ls -la your_file.zip
python sardine_user_friendly_cli.py /full/path/to/file.zip output/
```

**"No polarization available" error:**
```bash
# Check available polarizations first
python sardine_user_friendly_cli.py input.zip output/ --polarization VH
# Try VV if VH not available
python sardine_user_friendly_cli.py input.zip output/ --polarization VV
```

**Memory or disk space issues:**
```bash
# Use higher multilooking to reduce memory usage
python sardine_user_friendly_cli.py input.zip output/ \\
  --range-looks 8 --azimuth-looks 2
```

**Processing too slow:**
```bash
# Skip optional steps for faster processing
python sardine_user_friendly_cli.py input.zip output/ \\
  --polarization VV  # Skip terrain-flatten and speckle-filter
```

### 📖 Scientific References

When using SARdine in research, please cite:

1. **SARdine**: [Citation to be added]
2. **Terrain flattening**: Small, D. (2011). Flattening Gamma: Radiometric Terrain Correction for SAR Imagery. IEEE TGRS.
3. **Speckle filtering**: Lee, J.S. (1980). Digital image enhancement and noise filtering by use of local statistics. IEEE TPAMI.
4. **SAR processing**: Cumming, I.G. & Wong, F.H. (2005). Digital Processing of Synthetic Aperture Radar Data. Artech House.

### 💻 System Requirements

**Minimum:**
- Python 3.8+
- 8 GB RAM
- 20 GB free disk space
- Internet connection (for DEM downloads)

**Recommended:**
- Python 3.10+
- 16 GB RAM
- 100 GB free disk space (for caching)
- Multi-core CPU

### 🆘 Getting Help

1. **Check this user guide** - Most questions answered here
2. **Read quality reports** - Processing statistics and warnings
3. **Check example outputs** - Compare with provided examples
4. **GitHub Issues** - Report bugs and request features
5. **Documentation** - Full technical documentation available

### 🎓 Learning SAR Processing

**New to SAR?** Check these resources:
- ESA Sentinel-1 User Handbook
- NASA SAR Handbook
- "Digital Processing of Synthetic Aperture Radar Data" (Cumming & Wong)
- SARdine technical documentation

### 🚀 Advanced Features

**For expert users:**
- Custom DEM input support
- Batch processing scripts
- Integration with other tools
- Custom quality thresholds
- Algorithm parameter tuning

**See the full technical documentation for details.**

---

## 🎉 You're Ready to Start!

SARdine makes SAR processing accessible while maintaining scientific rigor. Whether you're a student learning SAR or a researcher needing production-quality results, SARdine provides the tools you need.

**Start with the Quick Start option and explore from there!**
