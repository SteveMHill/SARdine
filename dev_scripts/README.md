# Development Scripts

This directory contains development and debugging tools for SARdine.

## Scripts

### Debugging Tools
- **`debug_orbit_times.py`** - Debug orbit time format and extraction issues
- **`debug_terrain_correction_params.py`** - Debug terrain correction parameters and DEM issues
- **`compare_coordinates.py`** - Compare coordinate systems and transformations

### Validation Tools
- **`validate_no_synthetic_data.py`** - Validate that no synthetic/fallback data is used
- **`test_mask_fix.py`** - Test mask generation and validity

## Usage

These scripts are for development and debugging purposes only. They are not part of the main SARdine pipeline.

Example:
```bash
# Run validation
python3 dev_scripts/validate_no_synthetic_data.py

# Debug orbit times
python3 dev_scripts/debug_orbit_times.py

# Debug terrain correction
python3 dev_scripts/debug_terrain_correction_params.py
```

## Output

Scripts may create temporary output files in this directory. These are ignored by git and cleaned up automatically.
