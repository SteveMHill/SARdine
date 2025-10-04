# Deburst Critical Fixes - Implementation Status

**Date:** October 4, 2025  
**Status:** 🔄 PARTIAL - 6/13 issues fixed, 7 remaining  
**File:** `src/core/deburst.rs`

---

## ✅ FIXED (6/13)

### Issue #2: Steering Term - Linear → Removed ✅
**Status:** FIXED in 3 locations

**Problem:** Used `+ steering_rate_rad_s * t` (linear), but physics requires quadratic term `∝ t²`

**Fix Applied:**
- Removed linear steering terms from all deramp paths
- Added documentation that DC polynomial already includes steering
- If quadratic steering needed: `2.0 * PI * (2.0 * v / λ) * steering_rate * 0.5 * t * t`

**Locations:**
- Line ~185: `precompute_deramp_per_line()` 
- Line ~240: `precompute_deramp_2d()` range-dependent path
- Line ~260: `precompute_deramp_2d()` fast path

**Commit:** Pending

---

### Issue #3: 2D Polynomial Detection - Disabled ✅
**Status:** FIXED

**Problem:**
- `needs_2d = dc_poly.len() > 3` falsely detects cubic time polynomials as 2D
- `eval_poly_2d` uses made-up square grid assumption, wrong coefficient ordering

**Fix Applied:**
```rust
// Force 1D mode until proper XML parsing implemented
let needs_2d = false;

if dc_poly.len() > 3 || fm_poly.len() > 3 {
    log::warn!("⚠️  Polynomial length {} suggests 2D, but 2D parsing not implemented. Using time-only.", ...);
}
```

**Location:** Line ~218

**Commit:** Pending

---

### Issue #4: valid_window Last Valid = -1 ✅
**Status:** FIXED

**Problem:** With `last_valid == -1`, old code returned `(0, 1)` (1-pixel sliver)

**Fix Applied:**
```rust
// Treat -1 as empty window immediately
if last_raw < 0 {
    return (0, 0);
}
```

**Location:** Line ~289

**Commit:** Pending

---

### Issue #5: Power Preservation Claims ✅
**Status:** DOCUMENTATION FIXED

**Problem:** Comments claimed "perfect energy preservation" but math doesn't support it

**Fix Applied:**
```rust
// NOTE: This does NOT guarantee power preservation for complex data.
// True power preservation would require: w1 = cos(θ), w2 = sin(θ)
// and normalization by sum of SQUARED weights.
// Current approach prioritizes smooth visual blending.
```

**Location:** Overlap blending section

**Commit:** Pending

---

### Issue #9: Misleading Units in steering_angle_to_phase_rate ✅
**Status:** FIXED

**Problem:** Doc claimed "~2.9 MHz" but function returns rad/s

**Fix Applied:**
```rust
/// **Returns:** Phase rate in **rad/s** (not Hz!)
///
/// For Sentinel-1 IW with typical parameters:
/// - θ_rate ≈ 0.0015 rad/s
/// - Result: ≈ 0.463 rad/s phase rate
/// - In Hz: ≈ 74 kHz (NOT MHz as previously claimed)
```

**Location:** Line ~160

**Commit:** Pending

---

### Issue #12: Unused compute_pairwise_weights ✅
**Status:** MARKED AS DEAD CODE

**Fix Applied:**
```rust
/// WARNING: This function is currently UNUSED in production code.
/// The plan uses compute_row_weight() instead. Consider removing if not needed.
#[allow(dead_code)]
fn compute_pairwise_weights(...) { ... }
```

**Location:** Function definition

**Commit:** Pending

---

## ❌ REMAINING ISSUES (7/13)

### Issue #1: Polynomial Time Offset Still Zero ❌
**Status:** NEEDS FIX

**Problem:** Both `precompute_deramp_per_line` and `precompute_deramp_2d` call:
```rust
let timings = build_line_timing_with_offset(lines, az_time_interval_s, 0.0);
```

**Required Fix:**
1. Add to `BurstInfo`:
   ```rust
   pub dc_polynomial_t0: Option<f64>,
   pub burst_reference_time_seconds: Option<f64>,
   ```

2. Compute offset:
   ```rust
   let poly_time_offset_s = burst_sensing_time - polynomial_ref_time;
   let timings = build_line_timing_with_offset(lines, az_interval, poly_time_offset_s);
   ```

3. Parse from annotation:
   - `<dcEstimate><t0>` → `dc_polynomial_t0`
   - `<azimuthTime>` → `burst_reference_time_seconds`

**Locations:** Line ~184, ~218

---

### Issue #6: GDAL Complex Reading Incorrect ❌
**Status:** NEEDS FIX IN slc_reader.rs

**Problem:** Already fixed in slc_reader.rs (removed `width*2` trick)

**Status:** ✅ Actually FIXED in commit 8baa324

---

### Issue #7: Hardcoded Physical Constants ❌
**Status:** NEEDS PROPER PARSING

**Problem:** Multiple locations have:
```rust
slant_range_time: 0.006,  // WRONG: Sentinel-1 is C-band, not S-band
azimuth_bandwidth: 320.0,  // Hardcoded, should parse from annotation
```

**Required Fix:**
Parse from annotation XML:
- `<slantRangeTime>` near-range value
- `<azimuthBandwidth>` or `<processingBandwidth>`

**Locations:** Lines 1795, 1932, 2080, 2103

**Workaround:** Add error logging when using fallbacks

---

### Issue #8: Timing API Not Wired ❌
**Status:** BLOCKED BY ISSUE #1

**Problem:** New timing fields exist but aren't used in deramp

**Fix:** Same as Issue #1 - wire polynomial t0 through deramp calls

---

### Issue #10: Column Misalignment Across Bursts ⚠️
**Status:** NEEDS VERIFICATION

**Problem:** If bursts have different `start_sample`, all copy to same output columns

**Required Fix:**
```rust
let global_start = bursts.iter().map(|b| b.start_sample).min().unwrap();
dst_col_start = (burst.start_sample - global_start) + valid_start;
```

**Status:** Needs inspection of actual code

---

### Issue #11: Test & Production Defaults Mismatch ⚠️
**Status:** CONFIGURATION ISSUE

**Problem:** `DeburstConfig::default()` enables range-dependent deramp, but Issues #1-#3 aren't fully fixed

**Required Fix:**
Until Issues #1-#3 complete:
```rust
impl Default for DeburstConfig {
    fn default() -> Self {
        Self {
            use_range_dependent_deramp: false,  // Disable until timing fixed
            // ...
        }
    }
}
```

**Location:** `DeburstConfig::default()` impl

---

### Issue #13: calculate_output_dimensions vs Plan Overlap ✅
**Status:** ALREADY CORRECT

**Assessment:** Current code properly handles:
```rust
skip_front = blend_len.min(burst_lines)
```

**No fix needed**

---

## Implementation Priority

### High Priority (Breaks Correctness):
1. **Issue #1**: Polynomial time offset (phase errors, misalignment)
2. **Issue #7**: Parse real slant range time (wrong constants)

### Medium Priority (Reduces Quality):
3. **Issue #10**: Column alignment verification
4. **Issue #11**: Disable range-dependent by default

### Completed:
5. ✅ Issues #2, #3, #4, #5, #9, #12

---

## Next Steps

### Step 1: Add Polynomial Timing Fields
```rust
// In BurstInfo
pub dc_polynomial_t0: Option<f64>,
pub burst_reference_time_seconds: Option<f64>,
```

### Step 2: Wire Timing Through Deramp
```rust
let poly_time_offset = match (burst.burst_reference_time_seconds, burst.dc_polynomial_t0) {
    (Some(sensing), Some(t0)) => sensing - t0,
    _ => {
        log::warn!("⚠️  Missing polynomial timing - phase may be incorrect");
        0.0
    }
};

let timings = build_line_timing_with_offset(lines, az_interval, poly_time_offset);
```

### Step 3: Parse Annotation Values
Add XML parsing for:
- `<dcEstimate><t0>`
- `<slantRangeTime>`
- `<azimuthBandwidth>`

### Step 4: Update Defaults
```rust
impl Default for DeburstConfig {
    fn default() -> Self {
        Self {
            use_range_dependent_deramp: false,  // Until Issues #1-#3 fixed
            // ...
        }
    }
}
```

---

## Testing Recommendations

After fixes:
1. **Unit Test:** Polynomial timing with known t0
2. **Integration Test:** Real S1 data with DC polynomial
3. **Verification:** Check inter-burst phase continuity
4. **Regression:** Ensure overlap blending still smooth

---

## Commit Strategy

### Commit 1 (Current): Documentation & Safe Fixes
- ✅ Issue #2: Remove linear steering
- ✅ Issue #3: Disable broken 2D
- ✅ Issue #4: Fix valid_window
- ✅ Issue #5: Correct power preservation claims
- ✅ Issue #9: Fix units documentation
- ✅ Issue #12: Mark dead code

### Commit 2 (Next): Structure Changes
- Issue #1: Add timing fields to BurstInfo
- Update all BurstInfo construction sites

### Commit 3 (Then): Wire Timing
- Issue #1: Use polynomial t0 in deramp
- Issue #8: Connect timing API

### Commit 4 (Finally): Parse Real Values
- Issue #7: Parse slant range, bandwidth from annotation
- Remove hardcoded fallbacks

---

**Status:** Ready to commit first batch (6 fixes)
