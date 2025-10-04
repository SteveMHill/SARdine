# SARdine Refactoring Documentation Index

**Date:** October 4, 2025  
**Purpose:** Navigation guide for refactoring documentation

---

## Quick Start

1. **Want overview of current status?** → Read `CLEANUP_SUMMARY.md`
2. **Need detailed refactoring plan?** → Read `REFACTORING_CLEANUP_PLAN.md`
3. **Want to find violations?** → Run `./find_violations.sh`
4. **Found violations?** → See `VIOLATIONS_FOUND.md` for details

---

## Document Structure

### 📋 Status & Overview

**File:** `CLEANUP_SUMMARY.md` (163 lines)  
**Purpose:** Quick reference for cleanup status  
**Contains:**
- What was deleted (calibration_optimized.rs files)
- Architecture summary with clean boundaries
- What needs cleanup (high-level)
- Implementation status checklist
- Success metrics

**Read this first if:** You want a 5-minute overview

---

### 📖 Detailed Plan

**File:** `REFACTORING_CLEANUP_PLAN.md` (620 lines)  
**Purpose:** Comprehensive refactoring guide  
**Contains:**
- Module-by-module keep/remove checklists
  - SLC Reader (Section 1)
  - Annotation Parsers (Section 2)
  - Deburst (Section 3)
  - Calibration Integration (Section 4)
- Delete/keep grep checklist (Section 5)
- Dependency map (Section 6)
- Implementation checklist with phases (Section 7)
- Files to modify (Section 8)
- Success criteria (Section 9)
- Next steps (Section 10)

**Read this when:** Planning implementation work

---

### 🔍 Violation Scanner

**File:** `find_violations.sh` (executable script)  
**Purpose:** Automated violation detection  
**Checks:**
1. Incidence angle hacks
2. Calibration functions in I/O
3. Noise removal in wrong places
4. Antenna pattern application
5. Dense LUT building in parsers
6. Coordinate mapping in wrong places
7. Power/intensity conversion

**Usage:**
```bash
./find_violations.sh
# Exit code 0 = no violations
# Exit code 1 = violations found
```

**Run this:** Before and after making changes

---

### 📊 Violation Report

**File:** `VIOLATIONS_FOUND.md` (280 lines)  
**Purpose:** Detailed analysis of found violations  
**Contains:**
- Summary (✅ 4 passed, ❌ 3 failed)
- Detailed findings for each violation:
  1. `calibrate_and_multilook()` in slc_reader.rs (HIGH priority)
  2. `antenna_pattern_correction` in deburst.rs (MEDIUM)
  3. Intensity conversion in deburst.rs (HIGH priority)
- Code snippets with line numbers
- Recommended solutions
- Action plan with priorities
- Investigation questions

**Read this when:** Fixing specific violations

---

## Related Session Documents

### Previous Work (October 4, 2025)

**File:** `SESSION_SUMMARY_OCT4_2025.md` (292 lines)  
**Git:** commit `69ea323`  
**Purpose:** Summary of core module fixes  
**Contents:**
- slc_reader.rs: 4/4 issues fixed
- deburst.rs: 6/13 issues fixed
- calibrate.rs: 4/4 issues fixed
- Binary search optimization (1,000× speedup)
- Performance analysis

**File:** `CALIBRATE_CRITICAL_FIXES_APPLIED.md` (420 lines)  
**Git:** commit `1465698`  
**Purpose:** Detailed calibration fixes  
**Contents:**
- Weight calculation fix
- Binary search optimization
- Coordinate mapper validation
- Performance analysis

---

## Architecture Overview

### Clean Module Boundaries

```
┌──────────────────────────────────────────────────┐
│              Pipeline Flow                       │
└──────────────────────────────────────────────────┘

SLC Reader (I/O only)
  ↓
  Returns: Complex SLC array + Geometry metadata
  ↓
Annotation Parser (Sparse vectors only)
  ↓
  Returns: Calibration/Noise/Antenna coefficient vectors
  ↓
Deburst (Geometry stitching only)
  ↓
  Returns: Complex image + Valid sample ranges
  ↓
Calibration Stage (All radiometry centralized)
  ├─ create_auto_coordinate_mapper()
  ├─ precompute_antenna_pattern_lut()
  ├─ precompute_lut() → invert_lut_in_place()
  └─ apply_fused_slc_calibration()
     (single pass: power + noise + cal + antenna)
  ↓
  Returns: Calibrated backscatter (σ⁰/β⁰/γ⁰)
```

### Key Principles

1. **I/O ≠ Radiometry** - Readers return raw data + geometry
2. **Parsers ≠ LUT builders** - Return sparse vectors, not dense arrays
3. **Deburst ≠ Calibration** - Phase-preserving stitching only
4. **Calibration = Centralized** - All radiometry in fused kernels

---

## Current Status

### ✅ Completed

- [x] Core module fixes (slc_reader, deburst, calibrate)
- [x] Binary search optimization (1,000× speedup)
- [x] Fused calibration kernels
- [x] Deleted duplicate optimization files
- [x] Created refactoring plan
- [x] Created violation scanner
- [x] Scanned for violations

### 🔄 In Progress

- [ ] Investigate `calibrate_and_multilook()` usage
- [ ] Check deburst intensity conversion context
- [ ] Create GitHub issues for violations

### 🔜 Upcoming

- [ ] Remove calibration function from slc_reader
- [ ] Fix intensity conversion in deburst
- [ ] Remove antenna_pattern_correction config
- [ ] Implement IncidenceAngleModel properly
- [ ] Refactor module interfaces
- [ ] Integration testing with real data

---

## How to Use This Documentation

### For Code Review

1. Read `CLEANUP_SUMMARY.md` to understand architecture
2. Run `./find_violations.sh` on your branch
3. Check `VIOLATIONS_FOUND.md` for context on any failures
4. Refer to `REFACTORING_CLEANUP_PLAN.md` for detailed fixes

### For Implementation

1. Read `REFACTORING_CLEANUP_PLAN.md` Section 7 (Implementation Checklist)
2. Pick a phase (Phase 2 is next: Remove Violations)
3. Check `VIOLATIONS_FOUND.md` for specific tasks
4. Implement fixes
5. Run `./find_violations.sh` to verify
6. Update checklist in plan

### For Testing

1. Run violation scanner before making changes
2. Make changes following the plan
3. Run violation scanner after changes
4. If violations remain, check `VIOLATIONS_FOUND.md`
5. Run integration tests with real data

### For New Contributors

1. Start with `CLEANUP_SUMMARY.md` (5-minute read)
2. Read architecture section in this index
3. Run `./find_violations.sh` to see current state
4. Pick a task from `VIOLATIONS_FOUND.md` action plan
5. Refer to detailed plan as needed

---

## Git Commit History

```bash
# Latest refactoring commits (newest first)
14ff777 - docs: Document refactoring violations found
38cd37c - chore: Add automated violation finder script
b03fb1b - docs: Add refactoring cleanup summary
8e9ea95 - docs: Add comprehensive refactoring cleanup plan
69ea323 - docs: Session summary for October 4, 2025
1465698 - fix: Apply 4 critical calibration fixes (binary search, weights, validation)
```

**View detailed history:**
```bash
git log --oneline --graph --all | head -20
```

---

## Quick Reference Commands

### Find Violations
```bash
./find_violations.sh
```

### Search for Specific Patterns
```bash
# Incidence angle hacks
rg "incidence_angle.*20\.0" SARdine/src/

# Calibration in wrong places
rg "fn.*calibrate" SARdine/src/io/ SARdine/src/core/deburst.rs

# Intensity conversion
rg "norm_sqr|to_power" SARdine/src/io/ SARdine/src/core/deburst.rs
```

### Check Function Usage
```bash
# Is calibrate_and_multilook used?
rg "calibrate_and_multilook" SARdine/

# Where is antenna_pattern_correction referenced?
rg "antenna_pattern_correction" SARdine/
```

### Verify Clean Boundaries
```bash
# Ensure calibration only in calibrate.rs
rg "apply_fused.*calibration" SARdine/src/ | grep -v calibrate.rs

# Ensure LUT precompute only in calibrate.rs
rg "precompute_lut" SARdine/src/ | grep -v calibrate.rs
```

---

## Success Criteria

### Code Quality

- [ ] No radiometry in I/O modules
- [ ] No dense LUTs in annotation parsers
- [ ] No calibration in deburst
- [ ] All radiometry centralized in calibrate.rs
- [ ] `./find_violations.sh` exits with code 0

### Performance

- [x] 1,000× speedup in azimuth bracketing
- [ ] 4-50× faster overall calibration (measure after refactor)
- [ ] No duplicate power/intensity computation

### Testing

- [x] Unit tests pass
- [ ] Integration test with real Sentinel-1 data
- [ ] Cross-platform tests (Windows/Linux/Mac)
- [ ] Benchmark validates performance gains

---

## Questions?

**Architecture questions:** See `REFACTORING_CLEANUP_PLAN.md` Section 6 (Dependency Map)  
**Implementation questions:** See `REFACTORING_CLEANUP_PLAN.md` Section 7 (Implementation Checklist)  
**Violation questions:** See `VIOLATIONS_FOUND.md` with detailed analysis  
**Quick overview:** See `CLEANUP_SUMMARY.md`

---

**Last Updated:** October 4, 2025  
**Status:** Documentation complete, ready for implementation
