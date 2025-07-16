# TOPSAR Merge Implementation

## Overview

TOPSAR (Terrain Observation with Progressive Scans SAR) merge is a critical preprocessing step for Sentinel-1 Interferometric Wide (IW) swath data. It combines the three IW sub-swaths (IW1, IW2, IW3) into a single wide-swath image while properly handling overlap regions between adjacent sub-swaths.

## Why TOPSAR Merge is Essential

### 1. **IW Acquisition Geometry**
- Sentinel-1 IW mode acquires data using three sub-swaths (IW1, IW2, IW3)
- Each sub-swath covers approximately 250 km in range
- Combined, they provide ~250 km total swath width
- Sub-swaths have overlapping regions to ensure continuous coverage

### 2. **Processing Order**
The TOPSAR merge must be performed at the correct stage in the processing chain:

```
RAW SLC Data
     ↓
1. Calibration (DN → σ⁰/γ⁰/β⁰)
     ↓
2. TOPSAR Merge ← YOU ARE HERE
     ↓
3. Multilooking (Speckle Reduction)
     ↓
4. Terrain Correction (Geocoding)
     ↓
5. Speckle Filtering (Optional)
```

### 3. **Why After Calibration?**
- Calibration ensures consistent radiometric scaling across sub-swaths
- Different sub-swaths may have slightly different calibration parameters
- Merging uncalibrated data would result in radiometric inconsistencies

### 4. **Why Before Multilooking?**
- TOPSAR merge works with Single Look Complex (SLC) data
- Preserves phase information needed for proper overlap handling
- Multilooking should be applied to the merged wide-swath image

## Technical Implementation

### Sub-swath Characteristics

| Sub-swath | Range (km) | Incidence Angle | Overlap with Next |
|-----------|------------|-----------------|-------------------|
| IW1       | 0-85       | 29.1° - 35.3°   | ~1.2 km with IW2  |
| IW2       | 84-162     | 35.2° - 40.9°   | ~1.2 km with IW3  |
| IW3       | 161-250    | 40.8° - 46.0°   | N/A               |

### Overlap Region Handling

The implementation provides several methods for handling overlap regions:

#### 1. **Feather Blending (Recommended)**
```rust
// Weight calculation for smooth transition
let distance_from_edge1 = (pixel_pos - overlap_start) as f32;
let distance_from_edge2 = (overlap_end - pixel_pos) as f32;
let total_distance = (overlap_end - overlap_start) as f32;

let weight1 = distance_from_edge2 / total_distance;
let weight2 = distance_from_edge1 / total_distance;

blended_value = value1 * weight1 + value2 * weight2;
```

#### 2. **Simple Average**
```rust
blended_value = (value1 + value2) / 2.0;
```

#### 3. **Priority-based (First/Second)**
Takes values from the first or second sub-swath in overlap regions.

### Output Grid Calculation

The merged image uses a unified grid that encompasses all sub-swaths:

```rust
pub struct OutputGrid {
    pub range_samples: usize,      // Total range samples
    pub azimuth_samples: usize,    // Total azimuth samples  
    pub range_pixel_spacing: f64,  // Consistent spacing (meters)
    pub azimuth_pixel_spacing: f64,// Consistent spacing (meters)
    pub near_range_time: f64,      // Reference timing
    pub azimuth_time_start: f64,   // Reference timing
}
```

## API Usage

### Python API

```python
import sardine

# Load SLC data
reader = sardine.SlcReader("S1A_IW_SLC_*.zip")

# Perform TOPSAR merge after calibration
merged_result = sardine.topsar_merge(
    input_file="S1A_IW_SLC_*.zip",
    polarization="VV",
    calibration_type="sigma0",
    overlap_method="feather",
    output_grid="auto"
)

# Extract results
intensity_data = merged_result["intensity_data"]
valid_mask = merged_result["valid_mask"]
metadata = merged_result["metadata"]
```

### CLI Usage

```bash
# Basic TOPSAR merge
sardine topsar-merge input.zip --polarization VV --output merged_vv.npy

# With specific overlap method
sardine topsar-merge input.zip \
    --polarization VV \
    --calibration-type gamma0 \
    --overlap-method feather \
    --output gamma0_merged.npy

# Process all polarizations
sardine topsar-merge input.zip --output merged_data.npy
```

### Rust API

```rust
use sardine::core::topsar_merge::TopsarMerge;

// Initialize merger
let merger = TopsarMerge::new(subswaths)?;

// Perform merge
let merged_data = merger.merge_subswaths(
    &subswath_data,
    OverlapMethod::Feather,
    &output_grid
)?;
```

## Quality Metrics

The implementation provides comprehensive quality metrics:

```python
metadata = merged_result["metadata"]
print(f"Sub-swaths merged: {metadata['num_swaths']}")
print(f"Overlap regions: {metadata['overlap_count']}")
print(f"Valid pixels: {metadata['valid_pixels']}")
print(f"Processing time: {metadata['processing_time']}")
```

## Validation and Testing

### 1. **Radiometric Consistency**
- Check that backscatter values are consistent across sub-swath boundaries
- Verify no artificial discontinuities in overlap regions

### 2. **Geometric Accuracy**
- Ensure proper alignment between sub-swaths
- Validate pixel spacing and timing parameters

### 3. **Coverage Completeness**
- Verify that the full swath width is covered
- Check for gaps or invalid regions

## Best Practices

### 1. **Processing Order**
Always follow the correct processing sequence:
- ✅ Calibration → TOPSAR Merge → Multilooking
- ❌ TOPSAR Merge → Calibration (incorrect)
- ❌ Multilooking → TOPSAR Merge (incorrect)

### 2. **Overlap Method Selection**
- **Feather blending**: Best for most applications (smooth transitions)
- **Average**: Simple but may introduce artifacts
- **Priority-based**: Use when one sub-swath has known issues

### 3. **Memory Considerations**
- TOPSAR merge processes large arrays (can be several GB)
- Consider processing polarizations separately for memory efficiency
- Use appropriate data types (float32 vs float64)

### 4. **Quality Control**
- Always check the valid_mask output
- Monitor processing logs for warnings
- Validate output dimensions and coverage

## Performance Optimization

### Memory Usage
```rust
// Efficient memory management
let chunk_size = 1000; // Process in chunks
for chunk in image_chunks(chunk_size) {
    process_chunk(chunk);
}
```

### Parallel Processing
```rust
// Multi-threaded overlap processing
use rayon::prelude::*;

overlap_regions.par_iter().for_each(|region| {
    process_overlap_region(region);
});
```

## Common Issues and Solutions

### 1. **Missing Sub-swaths**
```
Error: Only 1 sub-swath found
Solution: Verify input file contains complete IW data
```

### 2. **Memory Issues**
```
Error: Out of memory during merge
Solution: Process polarizations separately or increase available RAM
```

### 3. **Inconsistent Calibration**
```
Warning: Large radiometric differences between sub-swaths
Solution: Check calibration parameters and LUT files
```

### 4. **Geometric Misalignment**
```
Warning: Sub-swath overlap regions don't align
Solution: Verify orbit and timing information
```

## References

1. Sentinel-1 User Handbook (ESA)
2. TOPSAR Overview and Interferometry (Torres et al., 2012)
3. Sentinel-1 Toolbox Documentation
4. IW SLC Product Specification

## See Also

- [Complete TOPSAR Merge Workflow Example](../examples/complete_topsar_merge_workflow.py)
- [Calibration Implementation](CALIBRATION_IMPLEMENTATION.md)
- [Multilook Implementation](MULTILOOK_IMPLEMENTATION.md)
- [Terrain Correction Implementation](TERRAIN_CORRECTION_IMPLEMENTATION.md)
