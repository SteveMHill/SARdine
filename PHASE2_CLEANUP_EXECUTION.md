# Phase 2 Cleanup Execution Report

## Date: October 4, 2025
## Status: 🔄 IN PROGRESS

---

## Task List

### ✅ Part 1: Replace deprecated datetime_to_seconds() calls (COMPLETED)
- Updated 19 calls to use `crate::types::datetime_to_utc_seconds()`
- Files: `terrain_correction.rs`, `context_extraction.rs`
- Committed: 1fe8d86

### 🔄 Part 2: Update extract_range_doppler_params() callers (IN PROGRESS)
- **Issue**: Tests call `extract_range_doppler_params()` without orbit_vectors parameter
- **Files**: `tests/annotation_parser.rs` (2 calls)
- **Solution**: Add mock orbit vectors or make orbit_vectors optional with deprecation

### 📋 Part 3: Remove deprecated function definitions
- See CODE_CLEANUP_ANALYSIS.md for ~1,000 lines to remove
- Target files:
  - `src/core/terrain_correction.rs` (~500 lines)
  - `src/core/deburst.rs` (~200-300 lines)

### 📋 Part 4: Remove deprecated shim modules
- `src/core/deburst_optimized.rs` (244 lines)
- `src/core/optimized_calibration.rs` (8 lines)

### 📋 Part 5: orbit.rs time precision fixes
- Replace millisecond conversions with nanosecond precision
- Update chrono idioms (from_utc_datetime)
- Implement barycentric interpolation
- Fix binary search for nanosecond precision

### 📋 Part 6: dem.rs correctness fixes
- Fix .gz decompression (Cursor wrapper)
- SRTM tile size detection (1201 vs 3601)
- Array2::ones replacement
- Mosaic dimension fixes

---

## Detailed Execution Log

### Part 2: extract_range_doppler_params() Update Strategy

**Problem**: Tests don't have orbit data, but new signature requires it.

**Options**:
1. ❌ Make orbit_vectors optional - breaks type safety
2. ✅ Create minimal mock orbit vectors in tests
3. ✅ Add a test-only method that uses dummy values

**Chosen Solution**: Create mock orbit vectors in test

