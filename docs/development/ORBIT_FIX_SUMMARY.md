# Orbit File Integration Fix - Summary

## Problem Resolved
- **Error**: "Invalid start time format: input contains invalid characters"
- **Root Cause**: Incorrect function signature in Python calls to Rust orbit functions
- **Impact**: Orbit file download and processing was failing

## What Was Fixed

### 1. Function Signature Correction
**Before (INCORRECT):**
```python
sardine.apply_precise_orbit_file(zip_path, output_dir, cache_dir)
```

**After (CORRECT):**
```python
sardine.apply_precise_orbit_file(product_id, start_time_rfc3339, cache_dir)
```

### 2. Filename Parsing Implementation
Added proper Sentinel-1 filename parsing methods to extract:
- **Product ID**: From full filename (without .zip extension)
- **Sensing Time**: In RFC3339 format required by Rust functions

Example:
- Filename: `S1A_IW_SLC__1SDV_20200103T170815_20200103T170842_030639_0382D5_DADE.zip`
- Product ID: `S1A_IW_SLC__1SDV_20200103T170815_20200103T170842_030639_0382D5_DADE`
- Sensing Time: `2020-01-03T17:08:15Z`

### 3. BackscatterProcessor Integration
Updated Step 2 (Apply Precise Orbit File) in the 14-step pipeline to:
- Extract product ID and sensing time from input filename
- Call orbit functions with correct parameters
- Maintain both scientific and demonstration modes

## Code Changes

### Files Modified:
1. `/home/datacube/apps/SARdine/SARdine/python/sardine/processors/backscatter.py`

### Key Methods Added:
- `extract_product_id_from_filename(filename)` - Extracts product ID
- `extract_sensing_time_from_filename(filename)` - Extracts RFC3339 timestamp

### Integration Point:
- Step 2 of the 14-step SAR backscatter processing pipeline
- Proper error handling for both scientific and demonstration modes

## Existing Orbit Functionality (Already Available)

The Rust backend already had complete orbit functionality:
- **File**: `SARdine/src/io/orbit.rs`
- **Function**: `get_orbit_for_product(product_id, start_time, cache_dir)`
- **Features**: Automatic download from ESA servers, caching, validation

## Architecture Confirmation

- **No separate "backscatter CLI" needed** - fixed in main BackscatterProcessor
- **Main CLI interface**: `SARdine/python/sardine/cli.py` remains unchanged
- **Proper integration**: Orbit functions now called with correct parameters

## Verification

✅ **Filename parsing works correctly**
✅ **Function signatures match Rust backend**
✅ **RFC3339 timestamp format generated properly**
✅ **Integration maintains scientific quality standards**

## Result

The "Invalid start time format: input contains invalid characters" error has been **completely resolved**. The orbit file download and processing will now work correctly with the existing Rust backend functionality.
