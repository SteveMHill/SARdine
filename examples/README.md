# SARdine Examples

This directory contains practical examples demonstrating SARdine's capabilities for Sentinel-1 SAR processing.

## Complete Workflows

These examples show end-to-end processing pipelines:

### Basic Processing
- **[complete_calibration_workflow.py](complete_calibration_workflow.py)** - Radiometric calibration from raw SLC to sigma0
- **[complete_deburst_workflow.py](complete_deburst_workflow.py)** - Remove burst boundaries and create seamless images
- **[complete_iw_split_workflow.py](complete_iw_split_workflow.py)** - Split interferometric wide swath data by sub-swath

### Advanced Processing
- **[complete_speckle_filtering_workflow.py](complete_speckle_filtering_workflow.py)** - Apply speckle filters for noise reduction
- **[complete_terrain_flattening_workflow.py](complete_terrain_flattening_workflow.py)** - Topographic normalization using DEM

### Orbit Processing
- **[complete_orbit_workflow.py](complete_orbit_workflow.py)** - Download and apply precise orbit files
- **[orbit_file_example.py](orbit_file_example.py)** - Work with orbit data directly
- **[python_orbit_api.py](python_orbit_api.py)** - Python API for orbit operations

## Usage

Each example can be run independently:

```bash
# Basic usage
python examples/complete_calibration_workflow.py

# With custom data
python examples/complete_speckle_filtering_workflow.py --input your_data.zip
```

## Requirements

- SARdine installed (`pip install .` from repository root)
- Test data in `data/` directory (see [data/README.md](../data/README.md))
- Internet connection for DEM/orbit downloads

## Typical Processing Chain

For a complete SAR processing workflow:

1. **SLC Reading & Splitting** → `complete_iw_split_workflow.py`
2. **Orbit Application** → `complete_orbit_workflow.py` 
3. **Debursting** → `complete_deburst_workflow.py`
4. **Calibration** → `complete_calibration_workflow.py`
5. **Terrain Flattening** → `complete_terrain_flattening_workflow.py`
6. **Speckle Filtering** → `complete_speckle_filtering_workflow.py`

## Performance Tips

- Use SSD storage for large SLC files
- Enable multi-threading for processing large scenes
- Cache DEM tiles to avoid repeated downloads
- Use appropriate speckle filter based on your application

---

*For implementation details, see [docs/](../docs/) directory.*
