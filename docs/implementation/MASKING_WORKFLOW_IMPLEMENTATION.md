# Masking Workflow Implementation

## Overview

The masking workflow provides post-processing quality control for terrain-corrected SAR gamma0 data. It combines multiple validity criteria to identify and mask unreliable pixels, improving the quality and interpretability of SAR backscatter products.

## Implementation Details

### Core Components

#### 1. Local Incidence Angle (LIA) Computation

The local incidence angle represents the angle between the radar look vector and the local terrain surface normal. It's crucial for:
- Identifying areas with steep slopes that may have geometric distortions
- Understanding radiometric variations due to topography
- Quality assessment of terrain correction

**Surface Normal Calculation:**
```rust
// Compute surface gradients using central differences
let dz_dx = (dem[row, col+1] - dem[row, col-1]) / (2 * pixel_spacing);
let dz_dy = (dem[row+1, col] - dem[row-1, col]) / (2 * pixel_spacing);

// Surface normal vector (normalized)
let normal = normalize([-dz_dx, -dz_dy, pixel_spacing]);
```

**Local Incidence Angle:**
```rust
// cos(θ_lia) = |look_vector · surface_normal|
let cos_lia = abs(look_vector.dot(surface_normal));
```

#### 2. Multi-Criteria Masking

The workflow applies multiple validity criteria:

**Gamma0 Range Check:**
- Minimum threshold: -50 dB (default)
- Maximum threshold: +10 dB (default)
- Removes unrealistic backscatter values
- Filters outliers and processing artifacts

**DEM Validity Check:**
- Elevation threshold: -100 m (default)
- Removes areas below reasonable sea level
- Filters invalid DEM values (NaN, extreme values)

**Local Incidence Angle Check:**
- Cosine threshold: 0.1 (≈84°, default)
- Removes very steep slopes prone to geometric distortion
- Maintains areas with reliable terrain correction

#### 3. Mask Combination

All individual masks are combined using logical AND:
```rust
combined_mask = gamma0_mask && dem_mask && lia_mask;
```

### Rust Implementation

#### Core Structures

```rust
#[derive(Debug, Clone)]
pub struct MaskingWorkflow {
    pub lia_threshold: f64,     // cos(θ) threshold
    pub dem_threshold: f64,     // elevation threshold (m)
    pub gamma0_min: f32,        // minimum γ⁰ (dB)
    pub gamma0_max: f32,        // maximum γ⁰ (dB)
}

#[derive(Debug, Clone)]
pub struct MaskResult {
    pub combined_mask: Array2<bool>,    // final validity mask
    pub lia_cosine: Array2<f32>,        // LIA cosine values
    pub gamma0_mask: Array2<bool>,      // γ⁰ validity
    pub dem_mask: Array2<bool>,         // DEM validity
    pub lia_mask: Array2<bool>,         // LIA validity
    pub valid_pixels: usize,            // statistics
    pub total_pixels: usize,
    pub coverage_percent: f64,
}
```

#### Key Methods

**Surface Normal Computation:**
```rust
pub fn compute_surface_normal(&self, dem: &Array2<f32>, row: usize, col: usize) -> SurfaceNormal
```

**Local Incidence Angle:**
```rust
pub fn compute_local_incidence_angle(&self, normal: &SurfaceNormal, look_vector: &Vector3) -> f64
```

**Masking Workflow:**
```rust
pub fn apply_masking_workflow(
    &self,
    gamma0_data: &Array2<f32>,
    dem_array: &Array2<f32>,
    workflow: &MaskingWorkflow,
) -> Result<MaskResult, TerrainCorrectionError>
```

**Mask Application:**
```rust
pub fn apply_mask_to_gamma0(
    gamma0_data: &Array2<f32>,
    mask: &Array2<bool>,
    fill_value: f32,
) -> Array2<f32>
```

### Python API

#### Classes

```python
# Workflow configuration
workflow = sardine.create_masking_workflow(
    lia_threshold=0.1,      # cos(84°)
    dem_threshold=-100.0,   # meters
    gamma0_min=-50.0,       # dB
    gamma0_max=10.0         # dB
)

# Apply masking
result = sardine.apply_masking_workflow(
    corrector,    # TerrainCorrector instance
    gamma0_data,  # numpy array
    dem_data,     # numpy array
    workflow      # MaskingWorkflow instance
)

# Access results
combined_mask = result.get_combined_mask()
lia_cosine = result.get_lia_cosine()
coverage = result.coverage_percent

# Apply mask
masked_gamma0 = sardine.apply_mask_to_gamma0(
    gamma0_data, combined_mask, fill_value=np.nan
)
```

#### Statistics and Quality Metrics

```python
print(f"Total pixels: {result.total_pixels}")
print(f"Valid pixels: {result.valid_pixels}")
print(f"Coverage: {result.coverage_percent:.1f}%")

# Individual mask statistics
gamma0_valid = np.sum(result.get_gamma0_mask())
dem_valid = np.sum(result.get_dem_mask())
lia_valid = np.sum(result.get_lia_mask())
```

### CLI Interface

```bash
# Basic masking workflow
sardine mask input_gamma0.tif input_dem.tif output_masked.tif

# Custom thresholds
sardine mask input_gamma0.tif input_dem.tif output_masked.tif \
    --lia-threshold 0.15 \
    --dem-threshold -50.0 \
    --gamma0-min -40.0 \
    --gamma0-max 5.0 \
    --fill-value -999.0

# Save additional outputs
sardine mask input_gamma0.tif input_dem.tif output_masked.tif \
    --save-masks \    # Save individual mask components
    --save-lia        # Save LIA cosine values
```

## Configuration Guidelines

### Threshold Selection

**Local Incidence Angle (cos θ):**
- 0.05 (87°): Very permissive, includes steep slopes
- 0.1 (84°): Standard threshold for most applications
- 0.2 (78°): Conservative, excludes moderate slopes
- 0.3 (73°): Strict, only gentle terrain

**DEM Thresholds:**
- -200 m: Very permissive (deep ocean)
- -100 m: Standard (reasonable sea level)
- 0 m: Conservative (land only)
- 100 m: Strict (elevated terrain only)

**Gamma0 Range:**
- Conservative: [-50, +20] dB (very wide range)
- Standard: [-35, +5] dB (typical land surfaces)
- Strict: [-25, 0] dB (natural surfaces only)

### Application-Specific Settings

**Forest Monitoring:**
```python
workflow = sardine.create_masking_workflow(
    lia_threshold=0.15,     # Moderate slopes OK
    dem_threshold=0.0,      # Land surfaces only
    gamma0_min=-30.0,       # Forest range
    gamma0_max=-5.0
)
```

**Urban Analysis:**
```python
workflow = sardine.create_masking_workflow(
    lia_threshold=0.2,      # Gentle slopes preferred
    dem_threshold=0.0,      # Above sea level
    gamma0_min=-20.0,       # Include bright targets
    gamma0_max=10.0
)
```

**Coastal/Marine:**
```python
workflow = sardine.create_masking_workflow(
    lia_threshold=0.1,      # Standard slopes
    dem_threshold=-50.0,    # Include shallow water
    gamma0_min=-40.0,       # Include water surfaces
    gamma0_max=0.0
)
```

## Integration Workflow

### 1. Standard Processing Chain

```python
# 1. Terrain correction
corrector = sardine.create_terrain_corrector(...)
gamma0_corrected = sardine.terrain_correction(...)

# 2. Masking workflow
workflow = sardine.create_masking_workflow(...)
mask_result = sardine.apply_masking_workflow(
    corrector, gamma0_corrected, dem_data, workflow
)

# 3. Apply mask
gamma0_masked = sardine.apply_mask_to_gamma0(
    gamma0_corrected, 
    mask_result.get_combined_mask(),
    fill_value=np.nan
)

# 4. Quality assessment
print(f"Data coverage: {mask_result.coverage_percent:.1f}%")
```

### 2. Batch Processing

```python
def process_scene_with_masking(slc_file, orbit_file, dem_file, output_dir):
    # ... [standard processing steps] ...
    
    # Apply masking
    workflow = sardine.create_masking_workflow()
    mask_result = sardine.apply_masking_workflow(
        corrector, gamma0_data, dem_data, workflow
    )
    
    # Save products
    save_geotiff(gamma0_masked, f"{output_dir}/gamma0_masked.tif")
    save_geotiff(mask_result.get_combined_mask(), f"{output_dir}/validity_mask.tif")
    save_geotiff(mask_result.get_lia_cosine(), f"{output_dir}/lia_cosine.tif")
    
    # Return quality metrics
    return {
        'coverage_percent': mask_result.coverage_percent,
        'valid_pixels': mask_result.valid_pixels,
        'total_pixels': mask_result.total_pixels
    }
```

## Performance Considerations

### Memory Usage
- Masks use boolean arrays (1 byte per pixel)
- LIA computation requires additional float array
- Memory usage ≈ 6 bytes per pixel for all outputs

### Computational Complexity
- Surface normal: O(n) with boundary handling
- LIA computation: O(n) dot products
- Mask combination: O(n) logical operations
- Overall: Linear in image size

### Optimization Tips
- Process in tiles for large images
- Reuse TerrainCorrector for multiple scenes
- Consider downsampling for initial assessment

## Output Products

### 1. Masked Gamma0
- Primary product with invalid pixels removed
- Fill value configurable (NaN, -999, etc.)
- Preserves original data range and units

### 2. Combined Validity Mask
- Boolean array (true = valid pixel)
- Combines all quality criteria
- Use for downstream analysis

### 3. Individual Masks
- Gamma0 validity mask
- DEM validity mask  
- LIA validity mask
- Useful for diagnostic analysis

### 4. Local Incidence Angle
- Cosine values [0, 1]
- Angle values [0°, 90°]
- Important for radiometric analysis

## Quality Assessment

### Coverage Metrics
```python
# Overall coverage
coverage = mask_result.coverage_percent

# Per-criterion statistics
gamma0_coverage = np.sum(mask_result.get_gamma0_mask()) / mask_result.total_pixels * 100
dem_coverage = np.sum(mask_result.get_dem_mask()) / mask_result.total_pixels * 100
lia_coverage = np.sum(mask_result.get_lia_mask()) / mask_result.total_pixels * 100
```

### Spatial Patterns
- Identify systematic masking patterns
- Check for edge effects
- Validate against known terrain features

### Threshold Sensitivity
- Test multiple threshold combinations
- Analyze coverage vs. quality trade-offs
- Optimize for specific applications

## Future Enhancements

### 1. Advanced LIA Computation
- Use actual SAR geometry instead of vertical look
- Account for orbit position and look angle
- Implement range-dependent corrections

### 2. Additional Quality Criteria
- Layover/shadow detection
- Multi-temporal consistency
- Speckle noise assessment

### 3. Adaptive Thresholding
- Terrain-dependent thresholds
- Statistical outlier detection
- Machine learning-based classification

### 4. Performance Optimization
- GPU acceleration for large scenes
- Parallel processing capabilities
- Memory-efficient streaming

## References

1. Small, D., et al. (2022). "Guide to Sentinel-1 Geocoding", ESA Technical Note
2. Kellndorfer, J., et al. (2019). "SAR Data Processing and Analysis"
3. Flores-Anderson, A., et al. (2019). "Utility of SAR for Land Cover Mapping"

## See Also

- [Terrain Correction Implementation](TERRAIN_CORRECTION_IMPLEMENTATION.md)
- [Complete Masking Workflow Example](../examples/complete_masking_workflow.py)
- [Python API Reference](../README.md#python-api)
- [CLI Reference](../README.md#command-line-interface)
