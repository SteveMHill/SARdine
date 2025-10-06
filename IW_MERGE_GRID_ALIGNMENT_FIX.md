# IW Merge Grid Alignment Fix

## Problem
Range-Doppler solver now works (burst time fix successful!), but IW merge fails with:
```
CRITICAL: Range grid misalignment 44079.0 samples exceeds threshold (5.0)
```

## Root Cause
The `first_sample_global` field in `SubSwath` is being set to **burst-local coordinates** (first valid sample within burst, typically 0-100) instead of **swath-global coordinates** (range pixel offset in unified grid).

### Current Buggy Logic (annotation.rs lines 2177):
```rust
first_sample_global = derived_first_sample;  // ❌ WRONG: Local burst boundary
```

- `derived_first_sample` comes from `geom.first_valid_sample` 
- This is the first VALID sample WITHIN the burst (after invalid edge samples)
- For IW1/IW2/IW3, this is typically 0, 0, 0 or small values
- Result: All subswaths appear to start at same range → 44k sample misalignment when compared to physics-based calculation

### Correct Physics-Based Approach
For Sentinel-1 IW mode:
- IW1 has shortest slant range (closest to satellite)
- IW2 has medium slant range  
- IW3 has longest slant range (farthest from satellite)

The range offset between subswaths should be calculated from **slant range time differences**:

```
Δn_r = (τ_swath - τ_ref) × (c/2) / range_pixel_spacing
```

Where:
- `τ_swath` = `slant_range_time` for current swath
- `τ_ref` = `slant_range_time` for reference swath (typically IW1)
- `c` = speed of light (299,792,458 m/s)
- Factor of 2 because radar uses two-way travel time

## Solution

### Option 1: Calculate `first_sample_global` from slant range (RECOMMENDED)
In `annotation.rs`, after extracting all subswaths, calculate relative offsets:

```rust
// After creating all subswaths, set first_sample_global from slant range
let ref_srt = subswaths.values()
    .next()
    .map(|sw| sw.slant_range_time)
    .unwrap_or(0.0);

for swath in subswaths.values_mut() {
    let delta_tau = swath.slant_range_time - ref_srt;
    let delta_range_m = delta_tau * (SPEED_OF_LIGHT_M_S / 2.0);
    let range_offset_pixels = (delta_range_m / swath.range_pixel_spacing).round() as usize;
    swath.first_sample_global = range_offset_pixels;
}
```

### Option 2: Relax validation threshold (TEMPORARY WORKAROUND)
If annotation doesn't provide grid info, relax the 5-sample threshold:
```rust
if max_range_misalignment > 100.0 {  // Was 5.0
    // Only fail for truly enormous misalignments
}
```

## Expected Outcome
After fix:
- ✅ IW1: `first_sample_global = 0` (reference)
- ✅ IW2: `first_sample_global ≈ 5000-10000` (medium range offset)
- ✅ IW3: `first_sample_global ≈ 10000-20000` (far range offset)
- ✅ Range misalignment < 1 sample (physics-based calculation matches annotation)
- ✅ IW merge succeeds

## Priority
**HIGH** - This blocks terrain-corrected backscatter output after Range-Doppler fix

## References
- Sentinel-1 SAR User Guide: TOPS IW mode swath geometry
- ESA annotation XML: `slantRangeTime` field provides physical basis for grid alignment
