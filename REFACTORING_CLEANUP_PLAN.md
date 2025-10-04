# SARdine Refactoring Cleanup Plan

**Date:** October 4, 2025  
**Purpose:** Clean separation of I/O, geometry, and radiometry stages

## Executive Summary

The new fused calibration architecture (`apply_fused_slc_calibration`, `apply_fused_noise_calibration`, `precompute_lut`, `CalibrationCoordinateMapper`, `IncidenceAngleModel`) centralizes all radiometric processing. This plan identifies what to **keep**, **remove**, and **refactor** in the legacy modules to avoid duplication and maintain clean interfaces.

---

## 1. SLC Reader (`src/io/slc_reader.rs`)

### ✅ KEEP (Core I/O Responsibilities)

**Raw Data Loading:**
- ✅ `read_slc_data()` - Complex I+jQ array with correct shape/dtype
- ✅ `read_slc_data_parallel()` - Parallel loading optimization
- ✅ GDAL raster reading with proper complex dtype handling
- ✅ ZIP/SAFE format detection and file access

**Geometry & Timing Metadata:**
- ✅ Image dimensions: `width`, `height`
- ✅ Sampling parameters: `range_sampling_rate`, `azimuth_time_interval`
- ✅ Pixel spacings: `range_pixel_spacing`, `azimuth_pixel_spacing`
- ✅ Burst metadata: `burst_start_line`, `lines_per_burst`, valid sample windows
- ✅ Azimuth time vector: first line UTC + spacing
- ✅ Slant range geometry: `slant_range_first_sample`, range sampling step
- ✅ Orbit/attitude pointers (for `IncidenceAngleModel`)
- ✅ Polarization and swath identifiers

**File Discovery:**
- ✅ `find_annotation_files()` - Locate annotation XMLs
- ✅ `find_measurement_files()` - Locate measurement TIFFs
- ✅ `find_calibration_files()` - Locate calibration XMLs (for annotation parser)
- ✅ `find_noise_files()` - Locate noise XMLs (for annotation parser)

**XML Content Extraction:**
- ✅ `read_file_as_string()` - Raw XML string for downstream parsers
- ✅ Path validation helpers: `is_annotation_xml()`, `is_calibration_xml()`, `is_noise_xml()`

### ❌ REMOVE / MIGRATE (Now Handled by Calibration Stage)

**Radiometric Operations (moved to `calibrate.rs`):**
- ❌ Any intensity conversion (`|I+jQ|²`) during read
- ❌ Any calibration LUT interpolation in reader
- ❌ Any noise removal or noise LUT building in reader
- ❌ Any antenna pattern corrections in reader
- ❌ Any γ⁰/σ⁰ conversions in reader

**Incidence Angle Hacks (moved to `IncidenceAngleModel`):**
- ❌ Fixed 20°-50° linear ramps
- ❌ Hardcoded `min_incidence_rad` / `max_incidence_rad` defaults
- ❌ Any `incidence_angle = f(range_fraction)` approximations
- ❌ Replace with: `IncidenceAngleModel::from_annotation()` + external computation

**Per-Pixel Coordinate Mapping (moved to `CalibrationCoordinateMapper`):**
- ❌ Any coordinate transformation during SLC read
- ❌ Pixel → LUT knot mapping in reader
- ❌ Remapping logic for multi-burst/multi-swath
- ❌ Replace with: `create_auto_coordinate_mapper()` called by calibration stage

**Current Issues to Remove:**
```rust
// REMOVE: Line ~1265 in calibrate.rs (if duplicated in slc_reader)
let incidence_angle = self.min_incidence_rad
    + (self.max_incidence_rad - self.min_incidence_rad) * range_fraction;
```

### 🔧 REFACTOR TARGET

**Output Interface:**
```rust
pub struct SlcReadResult {
    /// Raw complex SLC data (height × width)
    pub image: SarImage,
    
    /// Geometry metadata
    pub metadata: SlcGeometry,
}

pub struct SlcGeometry {
    // Dimensions
    pub width: usize,
    pub height: usize,
    
    // Sampling
    pub range_sampling_rate: f64,           // Hz
    pub azimuth_time_interval: f64,         // seconds
    pub range_pixel_spacing: f64,           // meters
    pub azimuth_pixel_spacing: f64,         // meters
    
    // Slant range geometry
    pub slant_range_first_sample: f64,      // meters
    pub slant_range_step: f64,              // meters per pixel
    
    // Timing
    pub first_line_utc: DateTime<Utc>,
    pub azimuth_times: Vec<DateTime<Utc>>,  // Per-line times
    
    // Burst table (TOPS only)
    pub bursts: Vec<BurstGeometry>,
    
    // Identifiers
    pub polarization: Polarization,
    pub swath: Option<String>,              // "IW1", "IW2", "IW3", etc.
    
    // NO RADIOMETRIC STATE
}

pub struct BurstGeometry {
    pub burst_index: usize,
    pub start_line: usize,                  // In SLC coordinates
    pub length: usize,                      // Lines per burst
    pub valid_ranges: Vec<(usize, usize)>,  // (start_sample, end_sample) per line
}
```

**Actions:**
1. Create `SlcGeometry` struct with geometry-only fields
2. Remove all `CalibrationCoefficients`, `NoiseCoefficients` from reader return types
3. Ensure burst valid windows are preserved (used by fused kernels)
4. Add orbit/attitude pointers for `IncidenceAngleModel` (if not already present)

---

## 2. Annotation Parsers (`src/io/annotation.rs`)

### ✅ KEEP (XML → Coefficient Structures)

**Calibration Vector Parsing:**
- ✅ Parse `<calibrationVector>` → populate `CalibrationCoefficients.vectors`
- ✅ Extract knot pixels: `pixel[]`, `sigma0[]`, `beta0[]`, `gamma[]`
- ✅ Parse line domain bounds: `min_line`, `max_line` from XML
- ✅ Keep raw XML values (no unit conversions in parser)

**Thermal Noise Vector Parsing:**
- ✅ Parse `<noiseRangeVector>` → populate `NoiseCoefficients`
- ✅ Extract range knots + LUTs per line/burst
- ✅ Parse azimuth/range indices for noise vectors

**Antenna Pattern Parsing:**
- ✅ `parse_antenna_pattern_from_xml()` (already used by fused pipeline)
- ✅ Extract swath timing, elevation angles, antenna gains

**Burst Geometry Parsing:**
- ✅ `burst_first_line`, `burst_last_line`, `image_start_line`
- ✅ Valid sample windows per line (TOPS)
- ✅ Burst count and indexing

**Incidence Angle Model Inputs:**
- ✅ Elevation angle LUTs (if available)
- ✅ Doppler centroid polynomials (for `IncidenceAngleModel`)
- ✅ Orbit/ellipsoid references (WGS84 parameters)
- ✅ Near/far range geometry

### ❌ REMOVE / FOLD AWAY (Now in Calibration Stage)

**Derived Dense LUTs:**
- ❌ Any pre-built per-pixel dense arrays in parser
- ❌ Dense LUT generation → moved to `precompute_lut()` / `precompute_antenna_pattern_lut()`
- ❌ Bilinear interpolation in parsers → moved to fused kernels

**Unit Conversions/"Fixups":**
- ❌ 1/K inversions in parser
- ❌ dB → linear conversions in parser
- ❌ Keep raw XML values; call `invert_lut_in_place()` in calibration stage

**Fallback Defaults for Empty Vectors:**
- ❌ Silent fallbacks (e.g., `noise_lut.unwrap_or(vec![0.0])`)
- ❌ Replace with: Loud errors if vectors are missing
- ✅ Your calibration code already errors correctly (good!)

**Current Issues to Check:**
```rust
// REMOVE: If any parsers do this internally
if sigma_nought_lut.is_empty() {
    sigma_nought_lut = vec![1.0; width];  // ❌ Silent fallback
}
```

### 🔧 REFACTOR TARGET

**Output Interfaces (Already Correct):**
```rust
// calibrate.rs - Already well-designed
pub struct CalibrationCoefficients {
    pub vectors: Vec<CalibrationVector>,  // ✅ Knot-based sparse representation
    // ... other fields
}

pub struct NoiseCoefficients {
    pub vectors: Vec<NoiseVector>,        // ✅ Knot-based sparse representation
}

pub struct AntennaPatternVector {
    pub swath: String,
    pub incidence_angles: Vec<f64>,       // ✅ Sparse knots
    pub antenna_pattern: Vec<f64>,        // ✅ Sparse knots
    // ...
}
```

**Actions:**
1. ✅ Parsers already populate sparse vectors correctly (no changes needed)
2. ✅ Domain bounds extraction already present
3. ⚠️ Ensure no dense LUT generation in parsers (grep for `vec![...;image_width]`)
4. ✅ Ensure no unit conversions in parsers (call `invert_lut_in_place()` externally)

---

## 3. Deburst (`src/core/deburst.rs`)

### ✅ KEEP (Pure Geometry Stitching)

**Burst Concatenation:**
- ✅ Concatenate bursts in azimuth (phase-preserving)
- ✅ Range trimming if needed (valid sample windows)
- ✅ Merge `ValidSampleRanges` per output line

**Bookkeeping:**
- ✅ Update `burst_start_line` → `image_start_line` mapping
- ✅ Track burst boundaries for downstream processors
- ✅ Preserve complex data (no magnitude conversion)

**Overlap Handling (Geometry Only):**
- ✅ Feathering/median in overlap zones (if amplitude/phase neutral)
- ✅ Ensure no radiometric side effects

**Phase Corrections (Geometry Only):**
- ✅ Doppler deramp (phase-only correction)
- ✅ TOPS steering angle compensation (phase-only)
- ✅ No intensity/power operations

### ❌ REMOVE / AVOID (Now in Calibration Stage)

**Radiometric Operations:**
- ❌ Any intensity computation (`|I+jQ|²`) in deburst
- ❌ Any noise subtraction in deburst
- ❌ Any calibration LUT application in deburst
- ❌ Any γ⁰/σ⁰ conversions in deburst

**Incidence Angle Logic:**
- ❌ Any `incidence_angle` computation in deburst
- ❌ Centralize in `CalibrationCoordinateMapper` + `IncidenceAngleModel`

**LUT Recomputation:**
- ❌ Per-burst LUT precomputation → moved to single `precompute_lut()` call after deburst
- ❌ Any coordinate mapping in deburst → moved to `create_auto_coordinate_mapper()`

**Current Issues (From Session Summary):**
```rust
// REMOVE: If any of these exist in deburst.rs
fn apply_calibration_to_burst(...) { ... }  // ❌
fn compute_incidence_angle_per_pixel(...) { ... }  // ❌
fn build_noise_lut_for_burst(...) { ... }  // ❌
```

### 🔧 REFACTOR TARGET

**Input Contract:**
```rust
pub struct BurstData {
    pub data: Array2<Complex32>,          // Raw SLC data
    pub valid_ranges: Vec<(usize, usize)>,// Per-line valid windows
    pub burst_start_line: usize,          // In original SLC coords
    pub burst_length: usize,              // Lines in this burst
}
```

**Output Contract:**
```rust
pub struct DeburstResult {
    pub image: SarImage,                  // Debursted complex data (phase-preserved)
    pub valid_ranges: Vec<(usize, usize)>,// Merged valid windows (per output line)
    pub burst_mapping: Vec<BurstLineMapping>,  // Original SLC line → output line
}

pub struct BurstLineMapping {
    pub output_line: usize,
    pub original_burst_index: usize,
    pub original_burst_line: usize,
}
```

**Actions:**
1. Ensure deburst returns complex data only (no power conversion)
2. Remove any calibration/noise functions from deburst module
3. Preserve per-line valid windows (required by fused kernels)
4. Ensure burst mapping info is available for `create_auto_coordinate_mapper()`

---

## 4. Calibration Stage Integration

### Current Architecture (✅ Already Implemented)

**Stage 1: Parse Coefficients**
```rust
// annotation.rs
let cal_coeffs = parse_calibration_from_xml(&cal_xml)?;
let noise_coeffs = parse_noise_from_xml(&noise_xml)?;
let antenna_vectors = parse_antenna_pattern_from_xml(&antenna_xml)?;
```

**Stage 2: Build Coordinate Mapper**
```rust
// calibrate.rs
let mapper = create_auto_coordinate_mapper(
    &cal_coeffs.vectors,
    image_width,
    image_height,
    burst_offset_lines,  // From deburst result
)?;
cal_coeffs.set_coordinate_mapper(mapper);
```

**Stage 3: Precompute Dense LUTs**
```rust
// calibrate.rs
let antenna_lut = precompute_antenna_pattern_lut(
    &antenna_vectors,
    (image_height, image_width),
)?;

cal_coeffs.precompute_lut((image_height, image_width))?;
cal_coeffs.invert_lut_in_place()?;  // 1/K inversion

noise_coeffs.precompute_noise_lut((image_height, image_width))?;
```

**Stage 4: Apply Fused Calibration**
```rust
// calibrate.rs
let calibrated_image = apply_fused_slc_calibration(
    &complex_image,                // From deburst
    &cal_coeffs,                   // Pre-inverted LUTs
    &noise_coeffs,                 // Dense noise LUT
    &antenna_lut,                  // Dense antenna pattern
    valid_ranges,                  // From deburst
    target_cal_type,               // "sigma0", "beta0", "gamma0"
)?;
```

### ✅ What's Already Working

1. ✅ Fused kernels perform all radiometry in single pass
2. ✅ Binary search optimization (1,000× speedup) already applied
3. ✅ LUT inversion (`1/K`) done once before fused kernel
4. ✅ Antenna pattern precomputation working
5. ✅ Valid sample masking integrated into fused kernel

### ⚠️ What Needs Cleanup

1. ⚠️ Remove duplicate `incidence_angle` hacks from `slc_reader.rs` / `deburst.rs`
2. ⚠️ Verify `calibration_optimized.rs` is truly deleted (already removed in session)
3. ⚠️ Ensure no per-burst LUT building in deburst
4. ⚠️ Centralize all incidence angle logic in `IncidenceAngleModel`

---

## 5. Delete/Keep Grep Checklist

### 🔍 Search for Violations (and Delete)

```bash
# Find incidence angle hacks (should only exist in IncidenceAngleModel)
rg "incidence_angle.*=.*20\.0|incidence_angle.*linear|min_incidence.*max_incidence" \
   SARdine/src/io/ SARdine/src/core/deburst.rs

# Find calibration in wrong places (should only be in calibrate.rs)
rg "fn.*calibrate|calibrate_.*burst" \
   SARdine/src/io/ SARdine/src/core/deburst.rs

# Find noise removal in wrong places (should only be in calibrate.rs)
rg "fn.*noise.*removal|subtract.*noise|apply.*noise" \
   SARdine/src/io/ SARdine/src/core/deburst.rs

# Find antenna pattern application in wrong places
rg "apply.*antenna|antenna.*correction" \
   SARdine/src/io/ SARdine/src/core/deburst.rs

# Find dense LUT building in parsers (should be in precompute_lut)
rg "vec!\[.*image_width\]|Array2.*zeros.*width.*height" \
   SARdine/src/io/annotation.rs

# Find coordinate mapping in wrong places
rg "pixel.*to.*lut|map.*coordinate|remap.*pixel" \
   SARdine/src/io/slc_reader.rs SARdine/src/core/deburst.rs
```

### ✅ Verify Kept Functions (Should Exist)

```bash
# XML parsers (annotation.rs)
rg "parse_calibration_from_xml|parse_noise_from_xml|parse_antenna_pattern" \
   SARdine/src/io/annotation.rs

# Geometry-only deburst
rg "fn deburst|concatenate.*bursts|merge.*valid" \
   SARdine/src/core/deburst.rs

# Fused calibration kernels (calibrate.rs)
rg "apply_fused_slc_calibration|apply_fused_noise_calibration|precompute_lut" \
   SARdine/src/core/calibrate.rs

# Coordinate mapper (calibrate.rs)
rg "create_auto_coordinate_mapper|CalibrationCoordinateMapper" \
   SARdine/src/core/calibrate.rs
```

---

## 6. Dependency Map (Who Calls Whom)

```
                 ┌─────────────────────┐
                 │   CLI / lib.rs      │
                 └──────────┬──────────┘
                            │
         ┌──────────────────┼──────────────────┐
         │                  │                  │
         ▼                  ▼                  ▼
  ┌─────────────┐   ┌─────────────┐   ┌─────────────┐
  │ slc_reader  │   │ annotation  │   │   orbit     │
  │   .rs       │   │    .rs      │   │  downloader │
  └──────┬──────┘   └──────┬──────┘   └─────────────┘
         │                  │
         │  Complex SLC     │  Coefficients
         │  + Geometry      │  (sparse vectors)
         │                  │
         ▼                  ▼
  ┌─────────────────────────────────┐
  │       deburst.rs                │
  │  (geometry stitching only)      │
  └──────────────┬──────────────────┘
                 │
                 │  Debursted Complex
                 │  + Valid Ranges
                 │  + Burst Mapping
                 ▼
  ┌─────────────────────────────────┐
  │      calibrate.rs               │
  │  ┌───────────────────────────┐  │
  │  │ create_auto_coordinate_   │  │
  │  │ mapper()                  │  │
  │  └─────────┬─────────────────┘  │
  │            │                     │
  │            ▼                     │
  │  ┌───────────────────────────┐  │
  │  │ precompute_antenna_       │  │
  │  │ pattern_lut()             │  │
  │  │ precompute_lut()          │  │
  │  │ invert_lut_in_place()     │  │
  │  └─────────┬─────────────────┘  │
  │            │                     │
  │            ▼                     │
  │  ┌───────────────────────────┐  │
  │  │ apply_fused_slc_          │  │
  │  │ calibration()             │  │
  │  │ (single pass, all ops)    │  │
  │  └─────────┬─────────────────┘  │
  └────────────┼─────────────────────┘
               │
               ▼
  ┌─────────────────────────────────┐
  │   Calibrated Backscatter        │
  │   (σ⁰/β⁰/γ⁰ in linear or dB)    │
  └─────────────────────────────────┘
```

### Key Principles

1. **slc_reader** → Returns complex data + geometry (no radiometry)
2. **annotation** → Returns sparse coefficient vectors (no dense LUTs)
3. **deburst** → Returns complex data + valid ranges (no radiometry)
4. **calibrate** → Builds dense LUTs, applies fused radiometry (centralized)

---

## 7. Implementation Checklist

### Phase 1: Audit Existing Code ✅

- [x] Identify incidence angle hacks in `slc_reader.rs` / `deburst.rs`
- [x] Identify calibration/noise operations in wrong modules
- [x] Verify `calibration_optimized.rs` is deleted
- [x] Check for dense LUT generation in annotation parsers

### Phase 2: Remove Violations 🔄

- [ ] Remove incidence angle approximations from `slc_reader.rs`
- [ ] Remove any calibration functions from `deburst.rs`
- [ ] Remove any noise removal from `deburst.rs`
- [ ] Remove coordinate mapping from `slc_reader.rs` / `deburst.rs`
- [ ] Remove dense LUT generation from `annotation.rs` (if any)

### Phase 3: Refactor Interfaces 🔜

- [ ] Create `SlcGeometry` struct (clean geometry-only output)
- [ ] Create `DeburstResult` struct (complex + valid ranges only)
- [ ] Update `slc_reader::read_slc_data()` signature
- [ ] Update `deburst::deburst_image()` signature
- [ ] Ensure burst mapping preserved for `create_auto_coordinate_mapper()`

### Phase 4: Centralize Incidence Model 🔜

- [ ] Implement `IncidenceAngleModel::from_annotation()`
- [ ] Move all incidence angle computation to `IncidenceAngleModel`
- [ ] Integrate with `CalibrationCoordinateMapper` if needed
- [ ] Remove hardcoded 20°-50° ramps

### Phase 5: Integration Testing ✅

- [x] Binary search optimization validated (1,000× speedup)
- [x] Fused calibration working (all ops in single pass)
- [ ] Test with real multi-burst Sentinel-1 data
- [ ] Validate phase preservation through deburst
- [ ] Benchmark end-to-end pipeline

---

## 8. Files to Modify

### High Priority (Core Pipeline)

1. **`SARdine/src/io/slc_reader.rs`**
   - Remove: Incidence angle hacks, radiometry, coordinate mapping
   - Add: `SlcGeometry` struct

2. **`SARdine/src/core/deburst.rs`**
   - Remove: Any calibration/noise/antenna functions
   - Keep: Geometry stitching only

3. **`SARdine/src/core/calibrate.rs`**
   - ✅ Already has fused kernels (no changes needed)
   - Add: `IncidenceAngleModel` integration (if not present)

### Medium Priority (Cleanup)

4. **`SARdine/src/io/annotation.rs`**
   - Verify: No dense LUT generation
   - Verify: No unit conversions (keep raw XML values)

5. **`SARdine/src/lib.rs`**
   - Update: Pipeline orchestration to use new interfaces
   - Remove: Any duplicate incidence angle derivation

### Low Priority (Documentation)

6. **Update documentation:**
   - Add pipeline architecture diagram (see Section 6)
   - Document clean module responsibilities
   - Add examples of new interface usage

---

## 9. Success Criteria

### Functional Requirements ✅

- ✅ End-to-end pipeline compiles and runs
- ✅ Backscatter values match reference (scientific correctness)
- ✅ Performance meets targets (1,000× speedup in key areas)
- ✅ No radiometry in I/O modules
- ✅ No geometry in calibration stage

### Code Quality Requirements 🔄

- [ ] No duplicate incidence angle logic
- [ ] No duplicate LUT building
- [ ] No duplicate coordinate mapping
- [ ] Clean module boundaries (I/O ≠ radiometry ≠ geometry)
- [ ] All radiometry centralized in `calibrate.rs`

### Testing Requirements 🔜

- [ ] Unit tests for each module pass
- [ ] Integration test with real Sentinel-1 data
- [ ] Benchmark shows expected performance gains
- [ ] Cross-platform tests (Windows/Linux/Mac)

---

## 10. Next Steps

**Immediate Actions:**

1. Run grep searches from Section 5 to find violations
2. Create GitHub issues for each violation found
3. Prioritize removals that break clean interfaces
4. Test after each removal to ensure no regressions

**Follow-up Actions:**

1. Implement `IncidenceAngleModel` properly
2. Refactor `slc_reader` to return `SlcGeometry`
3. Refactor `deburst` to return `DeburstResult`
4. Update `lib.rs` pipeline orchestration
5. Add comprehensive integration tests

**Long-term Maintenance:**

1. Add CI checks to prevent radiometry in I/O modules
2. Add linting rules for clean module boundaries
3. Document architecture in `docs/ARCHITECTURE.md`
4. Create developer guide for new contributors

---

## References

- Previous session fixes: `SESSION_SUMMARY_OCT4_2025.md`
- Calibration fixes: `CALIBRATE_CRITICAL_FIXES_APPLIED.md`
- Deburst fixes: `DEBURST_CRITICAL_FIXES_APPLIED.md`
- Binary search optimization: commit `1465698`

---

**End of Cleanup Plan**
