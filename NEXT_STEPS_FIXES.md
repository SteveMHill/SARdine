# Next Steps: Violation Fixes

**Date:** October 4, 2025  
**Status:** Investigation Complete → Ready for Implementation  
**Reference:** `VIOLATION_INVESTIGATION_REPORT.md`

---

## ✅ Investigation Complete

All 3 violations have been investigated and classified:

1. ✅ **calibrate_and_multilook()** - Confirmed unused, ready to remove
2. ✅ **antenna_pattern_correction** - Not implemented, ready to remove  
3. ✅ **Intensity conversion** - Standalone helper, NOT in main deburst path

---

## 🎯 Phase 1: Quick Wins (30 minutes - DO THIS NOW)

### Fix 1: Remove `calibrate_and_multilook()` from slc_reader.rs

**File:** `SARdine/src/io/slc_reader.rs`  
**Lines:** 2643-2693 (51 lines to delete)

**Command:**
```bash
cd /home/datacube/apps/SARdine/SARdine

# Delete the function
# Lines 2643-2693 in src/io/slc_reader.rs
```

**Verification:**
```bash
# Compile
cargo build --release

# Should succeed with no errors
```

**Commit Message:**
```
refactor: Remove unused calibrate_and_multilook from slc_reader

This function violated module boundaries (I/O orchestrating processing).
- No callers found in codebase (verified via grep)
- Used old CalibrationProcessor API (not fused kernels)
- If needed, pipeline orchestration should be in lib.rs

Part of refactoring cleanup (see REFACTORING_CLEANUP_PLAN.md)
```

---

### Fix 2: Remove `antenna_pattern_correction` from deburst config

**File:** `SARdine/src/core/deburst.rs`  
**Lines to modify:** 4 locations

**Changes:**

1. **Remove documentation comment** (lines 660-675):
```rust
// DELETE these 16 lines of documentation
```

2. **Remove struct field** (line 676):
```rust
// DELETE this line:
pub antenna_pattern_correction: bool,
```

3. **Remove from Default impl** (line 726):
```rust
// DELETE this line:
antenna_pattern_correction: false, // NOT IMPLEMENTED YET...
```

4. **Remove validation check** (lines 799-803):
```rust
// DELETE this entire if block:
if self.config.antenna_pattern_correction {
    log::warn!("   Azimuth antenna pattern LUT correction is not yet applied");
}
```

**Verification:**
```bash
cargo build --release
cargo test
```

**Commit Message:**
```
refactor: Remove unimplemented antenna_pattern_correction from deburst

This config field was never implemented and creates confusion:
- Default was false (not implemented)
- Antenna pattern correction properly implemented in calibrate.rs
- Deburst should remain geometry-only (phase-preserving)

Antenna pattern is applied via precompute_antenna_pattern_lut() in
the calibration stage, which is the correct architectural location.

Part of refactoring cleanup (see REFACTORING_CLEANUP_PLAN.md)
```

---

## 📊 Phase 2: Verify (10 minutes - DO AFTER PHASE 1)

### Verification Steps

**Step 1: Run violation scanner**
```bash
cd /home/datacube/apps/SARdine
./find_violations.sh
```

**Expected Output:**
```
==================================================
🔍 SEARCHING FOR REFACTORING VIOLATIONS
==================================================

1. Incidence Angle Hacks
   ✅ Clean

2. Calibration Functions
   ✅ Clean  <-- Should be fixed now

3. Noise Removal
   ✅ Clean

4. Antenna Pattern Application
   ⚠️  (May still show extract_subswath_from_zip) <-- Expected, not main path

5. Dense LUT Building
   ✅ Clean

6. Coordinate Mapping
   ✅ Clean

7. Power/Intensity Conversion
   ⚠️  (May still show extract_subswath_from_zip) <-- Expected, not main path

==================================================
Status: 2/3 violations fixed (extract_subswath_from_zip is helper)
```

**Step 2: Run tests**
```bash
cd /home/datacube/apps/SARdine/SARdine
cargo test 2>&1 | tee test_results.txt
```

**Step 3: Check compilation**
```bash
cargo build --release
```

---

## 🔬 Phase 3: Decision on `extract_subswath_from_zip()` (Later This Week)

### Context

**Function:** `extract_subswath_from_zip()` at line 2168  
**Purpose:** Standalone helper to extract subswath data from ZIP  
**Issue:** Returns intensity (power) instead of complex data

### Options

**Option A: Mark as Test Helper** (Recommended if only used in tests)
```rust
#[cfg(test)]
pub fn extract_subswath_from_zip(...) -> SarResult<Array2<f32>> {
    // Existing implementation
}
```

**Option B: Refactor to Return Complex** (If used in production)
```rust
pub fn extract_subswath_complex_from_zip(...) -> SarResult<Array2<Complex32>> {
    // Modify convert_cint16_to_intensity → convert_cint16_to_complex
}
```

**Option C: Remove** (If obsolete)
```rust
// Delete function entirely if not used
```

### Investigation Needed

```bash
# Check if used anywhere
rg "extract_subswath_from_zip" SARdine/

# Check git history
cd SARdine
git log -p --all -S "extract_subswath_from_zip" | head -50

# If only in tests, use Option A
# If in production, use Option B
# If unused, use Option C
```

---

## 📈 Success Metrics

### After Phase 1
- [ ] 2 violations fixed (calibrate_and_multilook, antenna_pattern_correction)
- [ ] `cargo build --release` succeeds
- [ ] `cargo test` passes
- [ ] `find_violations.sh` shows improvement

### After Phase 2
- [ ] Verification complete
- [ ] Test results documented
- [ ] Ready for Phase 3 decision

### After Phase 3 (Later)
- [ ] Decision on extract_subswath_from_zip made
- [ ] All violations fixed or justified
- [ ] Documentation updated
- [ ] `find_violations.sh` clean (or acceptable)

---

## 🚀 Start Here

### Immediate Actions (Right Now)

```bash
# 1. Navigate to repo
cd /home/datacube/apps/SARdine

# 2. Create feature branch (recommended)
git checkout -b refactor/remove-violations

# 3. Open slc_reader.rs and delete lines 2643-2693

# 4. Open deburst.rs and remove 4 locations (see Fix 2 above)

# 5. Compile
cd SARdine
cargo build --release

# 6. Test
cargo test

# 7. Commit
git add src/io/slc_reader.rs src/core/deburst.rs
git commit -m "refactor: Remove unused calibrate_and_multilook and antenna_pattern_correction

See VIOLATION_INVESTIGATION_REPORT.md for details"

# 8. Run scanner
cd ..
./find_violations.sh

# 9. If all good, merge to main
git checkout main
git merge refactor/remove-violations
```

---

## 📚 Documentation Trail

1. **Start:** `REFACTORING_INDEX.md` - Master index
2. **Plan:** `REFACTORING_CLEANUP_PLAN.md` - Detailed cleanup plan
3. **Scan:** `find_violations.sh` - Automated violation scanner
4. **Results:** `VIOLATIONS_FOUND.md` - Initial scan results
5. **Investigation:** `VIOLATION_INVESTIGATION_REPORT.md` - Detailed investigation ← YOU ARE HERE
6. **Next Steps:** This document (`NEXT_STEPS_FIXES.md`)

---

## ⏱️ Time Estimate

- **Phase 1 (Fixes):** 30 minutes
- **Phase 2 (Verify):** 10 minutes
- **Phase 3 (Decision):** 1-2 hours (later this week)
- **Total Today:** 40 minutes to significantly improve codebase

---

## 🎯 Goal

**Immediate:** Remove 2 violations today (40 minutes)  
**This Week:** Complete Phase 3 decision  
**This Month:** Full refactoring complete with integration tests

---

## 💡 Pro Tips

1. **Work on a branch** - Safer for experimentation
2. **Commit after each fix** - Easy to revert if needed
3. **Run tests frequently** - Catch issues early
4. **Check git blame** - Understand why code was added
5. **Ask questions** - Better to clarify than guess

---

**Status:** Ready to execute Phase 1  
**Confidence:** HIGH (investigation thorough, fixes straightforward)  
**Risk:** LOW (no callers, not implemented, safe to remove)

**GO FOR IT! 🚀**
