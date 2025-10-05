# Doppler Polynomial Parsing Fix - Implementation Summary

**Date**: October 5, 2025  
**Status**: ✅ IMPLEMENTED  
**Impact**: CRITICAL - Prevents silent data corruption from zero polynomial fallback

---

## Problem Statement

### Root Cause
The `SinglePassXmlParser` in `metadata_parser.rs` had `process_timing_element` "abbreviated for brevity", so Sentinel-1 Doppler centroid polynomials from `<dopplerCentroid><dcEstimateList><dcEstimate><dataDcPolynomial>` were **never parsed** into `CompactTimingData.dc_polynomials`.

### Observed Symptoms
- Deburst fallback to zero DC polynomial (logged as warning)
- 6-7% power loss due to phase misalignment
- 5-10% uncovered pixels in seam regions
- Broken burst alignment despite "valid" processing

### Why It Wasn't Caught Earlier
1. **Silent Fallback**: Deburst used `vec![0.0, 0.0]` instead of erroring out
2. **Validator Weakness**: `validate_doppler_polynomials` only checked SubSwath presence, not actual polynomial data
3. **Logging Confusion**: Warning logged but processing continued

---

## Solution Design - Minimal Surgical Fix

### 1. Parse `<dataDcPolynomial>` Correctly (metadata_parser.rs)

**File**: `src/core/metadata_parser.rs`

#### A. Extended `handle_start_element` to Track Doppler/FM Context

```rust
// In SinglePassXmlParser::handle_start_element(...)
match tag_name {
    // ...existing...
    "dopplerCentroid" | "dcEstimateList" | "dcEstimate" | "dataDcPolynomial" => {
        // Track path for Doppler polynomial parsing
        self.parsing_state.accumulated_text.clear();
        if self.timing_data.is_none() {
            self.timing_data = Some(CompactTimingData::new());
        }
    }
    "fmRate" | "fmRateList" | "fmRateEstimate" | "dataFmratePolynomial" | "dataAzimuthFmRatePolynomial" => {
        // Track path for FM rate polynomial parsing
        self.parsing_state.accumulated_text.clear();
        if self.timing_data.is_none() {
            self.timing_data = Some(CompactTimingData::new());
        }
    }
    _ => {}
}
```

**Why**: Ensures `timing_data` is initialized and `accumulated_text` is cleared before parsing polynomial content.

#### B. Implemented Full `process_timing_element`

```rust
fn process_timing_element(&mut self, tag_name: &str, text: &str) -> SarResult<()> {
    let path = self.current_path.join("/");
    
    if let Some(ref mut td) = self.timing_data {
        // Doppler centroid polynomial (Sentinel-1: dataDcPolynomial)
        if tag_name == "dataDcPolynomial" && 
           (path.ends_with("dopplerCentroid/dcEstimateList/dcEstimate/dataDcPolynomial") ||
            path.contains("dcEstimate/dataDcPolynomial")) {
            let coeffs = self.parse_whitespace_separated_f64(text)?;
            if !coeffs.is_empty() {
                td.dc_polynomials.push(coeffs);
                log::debug!("Parsed DC polynomial with {} coefficients", td.dc_polynomials.last().unwrap().len());
            }
        }
        
        // Alternative DC polynomial tags (vendor variants)
        if (tag_name == "dcPolynomial" || tag_name == "dopplerCentroidCoefficients") &&
           path.contains("doppler") {
            let coeffs = self.parse_whitespace_separated_f64(text)?;
            if !coeffs.is_empty() {
                td.dc_polynomials.push(coeffs);
                log::debug!("Parsed DC polynomial (variant tag) with {} coefficients", td.dc_polynomials.last().unwrap().len());
            }
        }
        
        // FM rate polynomial (Sentinel-1: dataFmratePolynomial or dataAzimuthFmRatePolynomial)
        if (tag_name == "dataFmratePolynomial" || tag_name == "dataAzimuthFmRatePolynomial" || tag_name == "fmRatePolynomial") &&
           (path.contains("fmRate") || path.contains("azimuthFmRate")) {
            let coeffs = self.parse_whitespace_separated_f64(text)?;
            if !coeffs.is_empty() {
                td.fm_polynomials.push(coeffs);
                log::debug!("Parsed FM rate polynomial with {} coefficients", td.fm_polynomials.last().unwrap().len());
            }
        }
        
        // Optional: PRF per-burst if present
        if tag_name == "prf" && path.contains("burst") {
            let prf_val = text.parse::<f64>().map_err(|e|
                SarError::InvalidMetadata(format!("Invalid PRF '{}': {}", text, e)))?;
            td.prf_values.push(Hertz::new(prf_val));
        }
    }
    
    Ok(())
}
```

**Key Features**:
- **Path-aware**: Matches exact XML hierarchy to avoid false positives
- **Vendor variants**: Handles `dataDcPolynomial`, `dcPolynomial`, `dopplerCentroidCoefficients`
- **FM rate variants**: Handles `dataFmratePolynomial`, `dataAzimuthFmRatePolynomial`, `fmRatePolynomial`
- **Diagnostic logging**: Logs each polynomial as it's parsed

#### C. Added Polynomial Count Validation in `finalize_parsing`

```rust
fn finalize_parsing(&mut self) -> SarResult<()> {
    // ...existing bounds computation...
    
    // Validate timing data consistency
    if let Some(ref mut td) = self.timing_data {
        log::info!(
            "Timing parse complete: bursts={}, dc_polys={}, fm_polys={}, prf_values={}",
            td.num_bursts, td.dc_polynomials.len(), td.fm_polynomials.len(), td.prf_values.len()
        );
        
        // If burst count is known, enforce consistency
        if td.num_bursts > 0 {
            if !td.dc_polynomials.is_empty() && td.dc_polynomials.len() != td.num_bursts {
                return Err(SarError::InvalidMetadata(format!(
                    "DC polynomials count ({}) != num_bursts ({})",
                    td.dc_polynomials.len(), td.num_bursts
                )));
            }
            if !td.fm_polynomials.is_empty() && td.fm_polynomials.len() != td.num_bursts {
                return Err(SarError::InvalidMetadata(format!(
                    "FM polynomials count ({}) != num_bursts ({})",
                    td.fm_polynomials.len(), td.num_bursts
                )));
            }
        } else if !td.dc_polynomials.is_empty() {
            // Infer burst count from polynomial count if not explicitly set
            td.num_bursts = td.dc_polynomials.len();
            log::debug!("Inferred num_bursts={} from DC polynomial count", td.num_bursts);
        }
    }
    
    Ok(())
}
```

**Why**: Catches count mismatches early before they cause runtime failures in deburst.

---

### 2. Make Strict Validation Truly Strict (metadata_strictness.rs)

**File**: `src/core/metadata_strictness.rs`

#### Enhanced `validate_doppler_polynomials`

**Before**:
```rust
// Check for Doppler data presence
let has_doppler_data = subswath.burst_count > 0;
if !has_doppler_data {
    critical_missing.push(format!("doppler_centroid_{}", swath_id));
    errors.push(format!("CRITICAL: Doppler centroid polynomial missing for subswath {}", swath_id));
}
```

**After**:
```rust
// STRICT: Check for actual DC polynomial presence in SubSwath
if subswath.dc_polynomial.is_none() || subswath.dc_polynomial.as_ref().map_or(true, |p| p.is_empty()) {
    critical_missing.push(format!("dc_polynomial_{}", swath_id));
    errors.push(format!(
        "CRITICAL: Doppler centroid polynomial missing or empty for subswath {}. \
         Deburst cannot proceed without DC polynomials to prevent phase misalignment.",
        swath_id
    ));
} else if let Some(dc_poly) = &subswath.dc_polynomial {
    // Validate all coefficients are finite
    for (idx, &coeff) in dc_poly.iter().enumerate() {
        if !coeff.is_finite() {
            errors.push(format!(
                "CRITICAL: DC polynomial coefficient[{}] is non-finite for {}: {}",
                idx, swath_id, coeff
            ));
        }
    }
    
    // Evaluate DC polynomial at mid-burst to check Nyquist constraint
    let mid_time = subswath.burst_duration / 2.0;
    let dc_value = dc_poly.iter().enumerate()
        .fold(0.0, |acc, (i, &c)| acc + c * mid_time.powi(i as i32));
    
    let nyquist_freq = subswath.prf_hz.unwrap_or(0.0) / 2.0;
    if dc_value.abs() > nyquist_freq {
        errors.push(format!(
            "CRITICAL: DC value {:.1} Hz at mid-burst exceeds Nyquist frequency {:.1} Hz for {}",
            dc_value, nyquist_freq, swath_id
        ));
    }
    
    log::debug!("DC polynomial for {} has degree {} (coeffs: {:?})", 
               swath_id, dc_poly.len().saturating_sub(1), dc_poly);
}
```

**Key Improvements**:
1. **Checks actual polynomial data**, not just proxy flags
2. **Validates coefficient finiteness** (catches NaN/Inf)
3. **Evaluates polynomial at mid-burst** to check Nyquist constraint
4. **Detailed diagnostic logging** for debugging

---

### 3. Remove Zero Fallback in Deburst (deburst.rs)

**File**: `src/core/deburst.rs`

**Before**:
```rust
let dc_polynomial = Self::extract_dc_estimates_from_annotation(annotation_data)
    .unwrap_or_else(|| {
        log::error!("❌ CRITICAL: DC polynomial not found in annotation - this will cause azimuth misalignment");
        log::error!("   Expected: <dopplerCentroid><dcEstimateList><dcEstimate><dataDcPolynomial>");
        log::error!("   Impact: 6-7% power loss, 5-10% uncovered pixels, broken burst alignment");
        log::error!("   Using zero fallback (NOT RECOMMENDED for production)");
        vec![0.0, 0.0]
    });
```

**After**:
```rust
let dc_polynomial = Self::extract_dc_estimates_from_annotation(annotation_data)
    .ok_or_else(|| {
        SarError::Processing(
            "CRITICAL: DC polynomial not found in annotation. \
             Expected: <dopplerCentroid><dcEstimateList><dcEstimate><dataDcPolynomial>. \
             Impact: Deburst cannot proceed without DC polynomials to prevent phase misalignment, \
             6-7% power loss, and 5-10% uncovered pixels. \
             ABORTING to prevent silent data corruption.".to_string()
        )
    })?;

log::info!("✅ Parsed DC polynomial with {} coefficients: {:?}", 
           dc_polynomial.len(), dc_polynomial);
```

**Why**: Fail fast instead of silent corruption. Consistent with "no silent fallbacks" requirement.

---

### 4. Add Diagnostic Logging (deburst.rs)

#### A. Per-Burst Polynomial Logging

```rust
for (burst_idx, burst) in self.burst_info.iter().enumerate() {
    // ...existing code...
    
    // Diagnostic logging: Show which polynomial is being used
    log::debug!(
        "Deburst: using DC deg={} FM deg={} for burst {} (lines={}, samples={})",
        burst.dc_polynomial.len().saturating_sub(1),
        burst.fm_polynomial.len().saturating_sub(1),
        burst_idx,
        burst_lines,
        burst_samples
    );
    
    // ...deramp computation...
}
```

**Benefit**: Quickly verify that non-zero polynomials are being used at runtime.

---

## Common Pitfalls Addressed

### 1. Namespace Stripping
✅ **Status**: Handled  
`quick-xml` gives local name already; no prefix matching needed.

### 2. Case Sensitivity
✅ **Status**: Handled  
Tag matches are exact; no `.to_lowercase()` used.

### 3. Path Tracking
✅ **Status**: Handled  
- `current_path.push(tag)` on `Event::Start`
- `current_path.pop()` on `Event::End`
- `accumulated_text.clear()` on each `Event::Start`

### 4. Multiple Products (IW1/2/3)
✅ **Status**: Handled  
`td.dc_polynomials.push(coeffs)` appends, doesn't overwrite.

### 5. Variant Tag Names
✅ **Status**: Handled  
Matches:
- `dataDcPolynomial`, `dcPolynomial`, `dopplerCentroidCoefficients`
- `dataFmratePolynomial`, `dataAzimuthFmRatePolynomial`, `fmRatePolynomial`

---

## Testing Strategy

### Quick Validation Test

**Print first 80 chars of any `Event::Text` whose path contains "doppler":**

```rust
if path.contains("doppler") {
    log::debug!("Doppler text: {}", &text[..text.len().min(80)]);
}
```

If you see **numeric strings** near scientific notation values → tags are correct.  
If you see **nothing** → path matching is wrong.

### Integration Test

1. Run backscatter pipeline on real Sentinel-1 IW data
2. Check logs for:
   - `"Timing parse complete: bursts=X, dc_polys=X, fm_polys=X"`
   - `"Parsed DC polynomial with N coefficients: [...]"`
   - `"Deburst: using DC deg=N FM deg=M for burst X"`
3. Verify NO "Using zero fallback" messages
4. Verify power preservation > 99%

---

## Verification Checklist

### Build & Compilation
- [ ] `cargo build --release` succeeds
- [ ] No new warnings introduced
- [ ] All existing tests pass

### Parsing Validation
- [ ] `SinglePassXmlParser` extracts DC polynomials
- [ ] Polynomial count matches burst count
- [ ] Coefficients are finite
- [ ] Log shows `dc_polys=N` where N > 0

### Strict Validation
- [ ] `validate_doppler_polynomials` checks actual polynomial data
- [ ] Empty polynomials trigger CRITICAL error
- [ ] Non-finite coefficients detected
- [ ] Nyquist constraint enforced

### Deburst Behavior
- [ ] NO "Using zero fallback" log messages
- [ ] Fails fast with clear error if polynomials missing
- [ ] Per-burst diagnostic shows `DC deg=N` where N > 0
- [ ] Power preservation > 99%

---

## Impact Assessment

### Before Fix
- ❌ Zero polynomial fallback → 6-7% power loss
- ❌ 5-10% uncovered pixels in seams
- ❌ Silent data corruption
- ❌ Broken burst alignment

### After Fix
- ✅ Actual Doppler polynomials parsed from annotation
- ✅ Fail fast if polynomials missing (prevents corruption)
- ✅ Polynomial count validated against burst count
- ✅ Nyquist constraint enforced
- ✅ Comprehensive diagnostic logging
- ✅ Expected power preservation > 99%

---

## Next Steps

1. **Extract `range_sampling_rate` from annotation XML**  
   Currently set to `None` with TODO comment at line ~1429 in `slc_reader.rs`.  
   Parse from `<product><imageAnnotation><productInformation><rangeSamplingRate>`.

2. **Add comprehensive validation test suite**  
   Test polynomial parsing with various real-world annotation files.

3. **Test backscatter pipeline with real data**  
   Process full Sentinel-1 IW scene and verify:
   - Polynomial count matches burst count
   - Power preservation > 99%
   - No uncovered pixels
   - Terrain correction coordinates stay within bounds

4. **Document tag variant handling**  
   Add to user documentation which XML tag variants are supported.

---

## Files Modified

1. **`src/core/metadata_parser.rs`** (3 changes)
   - Extended `handle_start_element` to track Doppler/FM context
   - Implemented full `process_timing_element` with polynomial parsing
   - Added polynomial count validation in `finalize_parsing`

2. **`src/core/metadata_strictness.rs`** (1 change)
   - Enhanced `validate_doppler_polynomials` to check actual polynomial data

3. **`src/core/deburst.rs`** (2 changes)
   - Removed zero fallback, replaced with hard error
   - Added per-burst diagnostic logging

---

## References

- **Original Issue**: "Coordinates exceed bounds (azimuth 134,000 vs max 50,000)"
- **Root Cause Analysis**: `process_timing_element` abbreviated, never parsed polynomials
- **Design Principle**: "No silent fallbacks" from strict validation requirements
- **Sentinel-1 XML Structure**: `<dopplerCentroid><dcEstimateList><dcEstimate><dataDcPolynomial>`

---

**STATUS**: ✅ READY FOR BUILD & TEST
