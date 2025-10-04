# Refactoring Cleanup Summary

**Date:** October 4, 2025  
**Status:** Planning Complete ✅ | Implementation Ready 🔄

## What Was Done

### 1. Files Deleted ✅
- `/home/datacube/apps/SARdine/SARdine/src/core/optimized_calibration.rs` (removed)
- `/home/datacube/apps/SARdine/SARdine/src/core/calibration_optimized.rs` (removed)

These files contained duplicate/experimental optimization code that has been superseded by the fused calibration architecture in `calibrate.rs`.

### 2. Documentation Created ✅
- `REFACTORING_CLEANUP_PLAN.md` - Comprehensive 620-line cleanup guide

## Architecture Summary

### Clean Module Boundaries

```
┌─────────────────┐
│  slc_reader.rs  │  → Returns: Complex SLC + Geometry (NO radiometry)
└────────┬────────┘
         │
         ▼
┌─────────────────┐
│  annotation.rs  │  → Returns: Sparse coefficient vectors (NO dense LUTs)
└────────┬────────┘
         │
         ▼
┌─────────────────┐
│   deburst.rs    │  → Returns: Complex image + Valid ranges (NO radiometry)
└────────┬────────┘
         │
         ▼
┌─────────────────┐
│  calibrate.rs   │  → Builds dense LUTs, applies fused radiometry (CENTRALIZED)
└─────────────────┘
```

### Key Principles

1. **I/O modules don't do radiometry** - No calibration, noise removal, or antenna corrections in readers
2. **Parsers return sparse data** - No dense per-pixel arrays in annotation parsers
3. **Deburst is geometry-only** - Phase-preserving stitching, no intensity operations
4. **Calibration is centralized** - All radiometry happens in fused kernels

## What Needs Cleanup

### High Priority Violations to Remove

1. **Incidence angle hacks** in `slc_reader.rs` / `deburst.rs`
   - Fixed 20°-50° linear ramps
   - Hardcoded min/max incidence angles
   - → Move to `IncidenceAngleModel`

2. **Radiometric operations** in wrong modules
   - Any `calibrate_*()` functions in deburst
   - Any `noise_removal()` in deburst
   - Any coordinate mapping in slc_reader
   - → Move to `calibrate.rs` fused kernels

3. **Dense LUT building** in annotation parsers
   - Any `vec![...;image_width]` in parsers
   - Any pre-interpolated arrays
   - → Move to `precompute_lut()` in calibrate.rs

### Search Commands (From Cleanup Plan)

```bash
# Find incidence angle hacks
rg "incidence_angle.*=.*20\.0|min_incidence.*max_incidence" \
   SARdine/src/io/ SARdine/src/core/deburst.rs

# Find calibration in wrong places
rg "fn.*calibrate|calibrate_.*burst" \
   SARdine/src/io/ SARdine/src/core/deburst.rs

# Find noise removal in wrong places  
rg "fn.*noise.*removal|subtract.*noise" \
   SARdine/src/io/ SARdine/src/core/deburst.rs

# Find dense LUTs in parsers
rg "vec!\[.*image_width\]|Array2.*zeros.*width" \
   SARdine/src/io/annotation.rs
```

## Implementation Status

### ✅ Already Working (From Previous Session)

- Binary search optimization (1,000× speedup) ✅
- Fused calibration kernels (single pass) ✅
- LUT inversion (1/K) done once before kernel ✅
- Antenna pattern precomputation ✅
- Valid sample masking in fused kernel ✅

### 🔄 Needs Implementation

- [ ] Remove incidence angle hacks from slc_reader/deburst
- [ ] Remove calibration functions from deburst
- [ ] Create `SlcGeometry` struct (clean interface)
- [ ] Create `DeburstResult` struct (complex + valid ranges only)
- [ ] Implement `IncidenceAngleModel::from_annotation()`
- [ ] Update pipeline orchestration in `lib.rs`

### 🔜 Future Enhancements

- [ ] Add CI checks for module boundaries
- [ ] Create `docs/ARCHITECTURE.md`
- [ ] Add comprehensive integration tests
- [ ] Benchmark with real Sentinel-1 multi-burst data

## Next Steps

### Immediate (Today)

1. Run grep searches to find violations
2. Create GitHub issues for each violation
3. Start with highest priority removals

### This Week

1. Remove incidence angle hacks
2. Remove radiometry from deburst
3. Refactor slc_reader interface
4. Test with real data

### This Month

1. Implement `IncidenceAngleModel` properly
2. Complete integration testing
3. Document architecture
4. Add CI checks

## Success Metrics

### Performance ✅
- 1,000-14,000× speedup in azimuth bracketing (achieved)
- 4-50× faster overall calibration (expected)

### Code Quality 🔄
- Zero radiometry in I/O modules (pending cleanup)
- Zero duplicate LUT building (pending cleanup)
- Single source of truth for incidence angles (pending)

### Testing 🔜
- All unit tests pass (current: passing)
- Integration test with real data (pending)
- Cross-platform tests (pending)

## References

- **Detailed Plan:** `REFACTORING_CLEANUP_PLAN.md` (620 lines)
- **Previous Fixes:** `SESSION_SUMMARY_OCT4_2025.md`
- **Calibration Fixes:** `CALIBRATE_CRITICAL_FIXES_APPLIED.md`
- **Binary Search:** commit `1465698`
- **This Summary:** commit `8e9ea95`

---

**Status:** Ready to begin systematic cleanup following the plan in `REFACTORING_CLEANUP_PLAN.md`
