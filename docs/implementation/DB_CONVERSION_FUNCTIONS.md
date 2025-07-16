# dB Conversion Functions in SARdine

## Overview

SARdine now includes built-in functions for converting between linear and decibel (dB) scales, which are commonly used in SAR data analysis and visualization.

## Functions

### `linear_to_db(data)`

Converts linear values to decibels (dB) using the formula: `dB = 10 * log10(linear)`

**Parameters:**
- `data`: 2D numpy array of float64 values in linear scale

**Returns:**
- 2D numpy array of float64 values in dB scale

**Special value handling:**
- Values ≤ 0 are converted to `-inf` (negative infinity)
- `inf` values remain `inf`
- `NaN` values remain `NaN`

### `db_to_linear(data)`

Converts decibel (dB) values to linear scale using the formula: `linear = 10^(dB/10)`

**Parameters:**
- `data`: 2D numpy array of float64 values in dB scale

**Returns:**
- 2D numpy array of float64 values in linear scale

**Special value handling:**
- `-inf` values are converted to 0
- `inf` values remain `inf`
- `NaN` values remain `NaN`

## Usage Examples

### Basic Usage

```python
import numpy as np
import sardine

# Create linear data
linear_data = np.array([[1.0, 10.0, 100.0], [0.1, 0.01, 0.001]])

# Convert to dB
db_data = sardine.linear_to_db(linear_data)
print(db_data)
# Output: [[0. 10. 20.] [-10. -20. -30.]]

# Convert back to linear
restored_data = sardine.db_to_linear(db_data)
print(restored_data)
# Output: [[1.e+00 1.e+01 1.e+02] [1.e-01 1.e-02 1.e-03]]
```

### SAR Backscatter Analysis

```python
import sardine
import numpy as np

# Typical gamma0 backscatter values (linear scale)
gamma0_linear = np.array([
    [0.1, 0.05, 0.2],   # Water, forest, agricultural
    [0.8, 0.3, 0.15]    # Urban, bare soil, vegetation
])

# Convert to dB for analysis and visualization
gamma0_db = sardine.linear_to_db(gamma0_linear)
print("Gamma0 in dB:")
print(gamma0_db)

# Typical ranges:
# Water: < -15 dB
# Vegetation: -15 to -8 dB  
# Urban: > -8 dB

# Land cover classification
water_mask = gamma0_db < -15
vegetation_mask = (gamma0_db >= -15) & (gamma0_db < -8)
urban_mask = gamma0_db >= -8
```

### Integration with Masking Workflow

```python
# After terrain correction and masking
corrected_gamma0_linear = terrain_corrector.process(sar_data, dem_data)
mask_result = masking_workflow.apply_mask(corrected_gamma0_linear, ...)

# Convert to dB for visualization
gamma0_db = sardine.linear_to_db(corrected_gamma0_linear)

# Apply mask in dB domain
gamma0_db_masked = gamma0_db.copy()
gamma0_db_masked[~mask_result.combined_mask] = np.nan

# Statistics in dB (excluding masked pixels)
valid_pixels = ~np.isnan(gamma0_db_masked) & ~np.isinf(gamma0_db_masked)
mean_db = np.mean(gamma0_db_masked[valid_pixels])
std_db = np.std(gamma0_db_masked[valid_pixels])
```

## Performance Notes

- Both functions are implemented in Rust for optimal performance
- Operations are vectorized and use SIMD instructions when available
- Memory efficient with minimal copying
- Suitable for large SAR datasets

## Typical SAR Value Ranges

### Linear Scale (Gamma0)
- **Water**: 0.01 - 0.1
- **Vegetation**: 0.05 - 0.3
- **Urban**: 0.2 - 2.0
- **Snow/Ice**: 0.001 - 0.05

### dB Scale (Gamma0)
- **Water**: -20 to -10 dB
- **Vegetation**: -13 to -5 dB
- **Urban**: -7 to +3 dB
- **Snow/Ice**: -30 to -13 dB

## Error Handling

The functions handle edge cases gracefully:

```python
# Zero and negative values
data = np.array([[0.0, -1.0], [1e-10, np.inf]])
db_data = sardine.linear_to_db(data)
# Result: [[-inf, -inf], [-100., inf]]

# NaN propagation
data_with_nan = np.array([[1.0, np.nan], [10.0, 100.0]])
db_data = sardine.linear_to_db(data_with_nan)
# Result: [[0., nan], [10., 20.]]
```

## See Also

- `examples/test_db_conversion.py` - Basic function testing
- `examples/db_conversion_workflow.py` - Complete workflow example
- `examples/complete_masking_workflow.py` - Integration with masking
