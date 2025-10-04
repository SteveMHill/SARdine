# Refactoring Violations Found

**Date:** October 4, 2025  
**Script:** `find_violations.sh`  
**Status:** 3 violation categories found

---

## Summary

✅ **4 checks passed** (no violations)
❌ **3 checks failed** (violations found)

### Passed Checks ✅

1. ✅ **Incidence Angle Hacks** - None found in slc_reader.rs or deburst.rs
2. ✅ **Noise Removal** - Correctly only in calibrate.rs
3. ✅ **Dense LUT Building** - None in annotation parsers
4. ✅ **Coordinate Mapping** - Correctly only in calibrate.rs

### Failed Checks ❌

1. ❌ **Calibration Functions in I/O** (1 violation)
2. ❌ **Antenna Pattern Application** (7 references in wrong places)
3. ❌ **Power/Intensity Conversion** (11 instances in wrong places)

---

## Detailed Findings

### ❌ Violation 1: Calibration Function in slc_reader.rs

**File:** `src/io/slc_reader.rs`  
**Line:** 2646  
**Issue:** `calibrate_and_multilook()` method exists in reader

```rust
// src/io/slc_reader.rs:2646
pub fn calibrate_and_multilook(
```

**Problem:** This is a high-level processing function that combines multiple stages. It should not be in the I/O module.

**Solution:**
- Move to `lib.rs` as a pipeline orchestration function
- Or move to a new `pipeline.rs` module
- slc_reader should only return raw complex data + geometry

**Priority:** HIGH - This violates clean module boundaries

---

### ❌ Violation 2: Antenna Pattern References in deburst.rs

**File:** `src/core/deburst.rs`  
**Lines:** 660, 676, 726, 799, 800, 802  
**Issue:** Antenna pattern correction config exists but not implemented

```rust
// src/core/deburst.rs:660
/// Apply azimuth antenna pattern correction (ESA best practice)

// src/core/deburst.rs:676
pub antenna_pattern_correction: bool,

// src/core/deburst.rs:726
antenna_pattern_correction: false, // NOT IMPLEMENTED YET

// src/core/deburst.rs:799-802
if self.config.antenna_pattern_correction {
    // Validation warning
    log::warn!("Azimuth antenna pattern LUT correction is not yet applied");
}
```

**Also in:** `src/io/slc_reader.rs:2248`

**Problem:** 
- Config flag exists but feature not implemented
- Antenna pattern correction should happen in calibration stage, not deburst
- Creates confusion about where radiometry happens

**Solution:**
- **Option A (Recommended):** Remove all antenna_pattern_correction references from deburst
  - Antenna pattern already handled by `precompute_antenna_pattern_lut()` in calibrate.rs
  - Deburst should be geometry-only (phase-preserving)
  
- **Option B (If azimuth-specific pattern needed):**
  - Document that azimuth antenna pattern is applied in calibration stage
  - Remove from deburst config entirely
  - Add to calibration config if needed

**Priority:** MEDIUM - Not causing bugs (feature disabled), but violates architecture

---

### ❌ Violation 3: Power/Intensity Conversion in Wrong Places

**File:** `src/core/deburst.rs`  
**Lines:** 2262, 2338, 2365, 2366, 2371, 2408, 2417, 2418  
**Issue:** Intensity computation in helper functions

```rust
// src/core/deburst.rs:2262
let intensity_data = if band_count == 1 {

// src/core/deburst.rs:2338
let mut intensity_array = Array2::zeros((height, width));

// src/core/deburst.rs:2365-2366
let intensity = real_f32 * real_f32 + imag_f32 * imag_f32;
intensity_array[[row, col]] = intensity;

// Similar pattern at lines 2408, 2417-2418
```

**File:** `src/io/slc_reader.rs`  
**Line:** 2670  
**Issue:** Intensity data passed to calibration

```rust
// src/io/slc_reader.rs:2670
let intensity_data = processor.calibrate(&deburst_data)?;
```

**Also in:** `src/io/dem.rs:1042`, `src/io/annotation.rs:2508` (false positives - gradient/velocity magnitude, not radiometric)

**Problem:**
- Deburst should preserve complex data (phase information)
- Power/intensity conversion (`|I+jQ|²`) should only happen in fused calibration kernel
- Current code may be computing power twice (once in deburst, once in calibration)

**Solution:**

1. **Check if these are test/debug functions:**
   - If yes, mark with `#[cfg(test)]` or move to tests module
   - If helper functions for visualization, rename to clarify (e.g., `debug_intensity_for_visualization`)

2. **If these are production code paths:**
   - Refactor to keep complex data through deburst
   - Remove intensity conversion from deburst
   - Let `apply_fused_slc_calibration()` handle power conversion

3. **Specific action for line 2670:**
   - Change from: `let intensity_data = processor.calibrate(&deburst_data)?;`
   - Change to: `let calibrated_data = apply_fused_slc_calibration(&complex_data, ...)?;`

**Priority:** HIGH - This affects correctness (phase information loss) and performance (duplicate work)

---

## False Positives (Acceptable)

### ✅ `src/io/dem.rs:1042` - Gradient magnitude
```rust
let gradient_magnitude = (dz_dx * dz_dx + dz_dy * dz_dy).sqrt();
```
**Verdict:** OK - This is terrain slope calculation, not radiometric intensity

### ✅ `src/io/annotation.rs:2508` - Velocity magnitude
```rust
let velocity_magnitude = (vx * vx + vy * vy + vz * vz).sqrt();
```
**Verdict:** OK - This is orbit velocity, not radiometric intensity

---

## Action Plan

### Immediate (High Priority)

1. **Remove `calibrate_and_multilook()` from slc_reader.rs**
   ```bash
   # Check what it does
   rg -A 20 "pub fn calibrate_and_multilook" SARdine/src/io/slc_reader.rs
   
   # Move to lib.rs or mark deprecated
   ```

2. **Fix intensity conversion in deburst.rs**
   ```bash
   # Check context of intensity computations
   rg -B 5 -A 10 "intensity_data|intensity_array" SARdine/src/core/deburst.rs
   
   # Determine if test code or production
   # Refactor to preserve complex data
   ```

### Medium Priority

3. **Remove antenna_pattern_correction from deburst config**
   ```bash
   # Remove config field
   sed -i '/antenna_pattern_correction/d' SARdine/src/core/deburst.rs
   
   # Remove from slc_reader.rs as well
   sed -i '/antenna_pattern_correction/d' SARdine/src/io/slc_reader.rs
   
   # Test that nothing breaks
   cargo test
   ```

### After Fixes

4. **Re-run violation finder**
   ```bash
   ./find_violations.sh
   ```

5. **Verify with real data test**
   ```bash
   # Test that complex data flows through correctly
   # Test that calibration produces correct results
   ```

---

## Investigation Needed

### Questions to Answer

1. **Is `calibrate_and_multilook()` used anywhere?**
   ```bash
   rg "calibrate_and_multilook" SARdine/
   ```

2. **Are the intensity conversions in deburst actually used?**
   ```bash
   # Check if they're in main code path or just debug/test
   rg -B 10 "intensity_data.*if band_count" SARdine/src/core/deburst.rs
   ```

3. **What does the intensity_data flow look like?**
   ```bash
   # Trace from deburst output to calibration input
   rg "intensity.*=" SARdine/src/lib.rs
   ```

---

## Expected Outcomes After Fixes

### Clean Architecture ✅
- slc_reader returns: Complex SLC + Geometry
- deburst returns: Complex image + Valid ranges
- calibrate performs: All radiometry in fused kernel

### Performance ✅
- No duplicate power computation
- Single pass through pixel data
- Maintain 1,000× speedup from binary search

### Correctness ✅
- Phase information preserved until calibration
- No premature intensity conversion
- Fused kernel handles all radiometry

---

## References

- **Cleanup Plan:** `REFACTORING_CLEANUP_PLAN.md`
- **Session Summary:** `SESSION_SUMMARY_OCT4_2025.md`
- **Violation Script:** `find_violations.sh`
- **This Report:** commit `38cd37c`

---

## Next Steps

1. Investigate `calibrate_and_multilook()` usage
2. Check deburst intensity conversion context
3. Create GitHub issues for each violation
4. Implement fixes one by one
5. Re-run `find_violations.sh` after each fix
6. Update this document with resolution status

---

**Status:** Violations identified, awaiting investigation and fixes
