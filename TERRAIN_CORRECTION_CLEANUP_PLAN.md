# Terrain Correction Module Cleanup Plan

**Date:** October 5, 2025  
**Module:** `src/core/terrain_correction.rs`  
**Current State:** Functional but with excessive debug logging and deprecated field usage

## Issues Identified

### 1. Excessive Debug Logging (High Priority)
**Problem:** Many `log::error!()` calls used for routine diagnostic information

**Examples:**
```rust
log::error!("🔍 RD TRANSFORM START: lat={:.6}, lon={:.6}, elevation={:.1}", ...);
log::error!("🔍 RD TRANSFORM: About to call latlon_to_ecef");
log::error!("🔍 CLAMP DEBUG #1: best_time.clamp({:.1}, {:.1})", ...);
```

**Impact:**
- Misleading error logs (not actual errors)
- Log spam during normal operation
- Harder to identify real errors

**Fix:**
- Change routine diagnostics to `log::debug!()` or `log::trace!()`
- Keep `log::error!()` only for actual error conditions
- Use `log::warn!()` for suspicious but recoverable conditions

### 2. Deprecated Field Usage (Medium Priority)
**Problem:** Using deprecated `product_start_time_abs` and `product_stop_time_abs` fields

**Locations:**
- Line 2656: `params.product_start_time_abs`
- Line 2815: `params.product_start_time_abs`
- Multiple other locations throughout

**Impact:**
- 19 deprecation warnings in build
- Future incompatibility risk
- Confusing time base mixing

**Fix:**
- Migrate to `orbit_ref_epoch_utc + product_start_rel_s`
- Add helper methods to encapsulate the conversion
- Document time base usage clearly

### 3. Complex Clamping Logic (Low Priority)
**Problem:** Custom `diag_clamp()` function with complex inversion detection

**Current Implementation:**
```rust
fn diag_clamp(value: f64, min: f64, max: f64, label: &str) -> f64 {
    if min > max {
        if std::env::var("SARDINE_STRICT_BOUNDS").is_ok() {
            panic!("STRICT_CLAMP: inversion...");
        }
        log::error!("🚫 diag_clamp inversion...");
        return value.clamp(max, min); // Swap!
    }
    value.clamp(min, max)
}
```

**Issues:**
- Auto-swapping min/max can hide bugs
- Panic on SARDINE_STRICT_BOUNDS is too aggressive
- Doesn't help identify root cause

**Fix:**
- Add debug assertions to catch inversions in development
- Log warning but don't auto-swap (fail fast)
- Add context about where the inversion originated

## Implementation Phases

### Phase 1: Debug Logging Cleanup (This Session)

**Target Lines:**
- 2372-2375: RD TRANSFORM START → `log::debug!()`
- 2378: "About to call latlon_to_ecef" → `log::trace!()`
- 2381: ECEF result → `log::trace!()`
- 2471: COORDINATE DEBUG → `log::debug!()`
- 2527: PIXEL CALC DEBUG → `log::debug!()`
- All "CLAMP DEBUG #N" → `log::trace!()`

**Estimated Impact:** ~30-40 log level changes

### Phase 2: Deprecated Field Migration (Next Session)

**Strategy:**
1. Add helper methods:
   ```rust
   impl RangeDopplerParams {
       pub fn product_start_time_abs(&self) -> f64 {
           self.orbit_ref_epoch_utc + self.product_start_rel_s
       }
       
       pub fn product_stop_time_abs(&self) -> f64 {
           self.product_start_time_abs() + self.product_duration
       }
   }
   ```

2. Replace direct field access with method calls
3. Remove #[allow(deprecated)] attributes
4. Verify all tests pass

**Estimated Impact:** 19 deprecation warnings → 0

### Phase 3: Clamping Logic Hardening (Future)

**Strategy:**
1. Add debug_assert! for min <= max
2. Remove auto-swapping behavior
3. Add proper error context
4. Document expected ranges

**Example:**
```rust
fn diag_clamp(value: f64, min: f64, max: f64, label: &str) -> f64 {
    debug_assert!(
        min <= max,
        "Clamp inversion at {}: min={} > max={} (value={})",
        label, min, max, value
    );
    
    if min > max {
        log::error!(
            "🚫 Clamp bounds inverted at {}: min={:.6} > max={:.6} (value={:.6}). \
             This is a bug in the calling code. Returning unclamped value.",
            label, min, max, value
        );
        return value; // Fail loudly, don't hide the bug
    }
    
    value.clamp(min, max)
}
```

## Expected Outcomes

### Phase 1 Completion
- ✅ Clean error logs (only real errors)
- ✅ Useful debug logs (can be enabled when needed)
- ✅ Trace logs for deep debugging
- ✅ Easier to diagnose real issues

### Phase 2 Completion
- ✅ Zero deprecation warnings
- ✅ Consistent time base usage
- ✅ Future-proof API usage
- ✅ Clear time reference documentation

### Phase 3 Completion
- ✅ Bugs caught early (debug assertions)
- ✅ No silent bug hiding (auto-swapping removed)
- ✅ Better error context
- ✅ Documented constraints

## Testing Strategy

After each phase:
1. Run `cargo test --release` (verify no regressions)
2. Run `cargo clippy` (check for new warnings)
3. Test with real Sentinel-1 data
4. Verify terrain correction output quality unchanged

## Breaking Changes

**None expected** - These are internal implementation improvements that maintain API compatibility and output quality.

## Timeline

- **Phase 1:** 1-2 hours (this session)
- **Phase 2:** 2-3 hours (next session)
- **Phase 3:** 1-2 hours (future enhancement)

Total estimated effort: 4-7 hours across 2-3 sessions
