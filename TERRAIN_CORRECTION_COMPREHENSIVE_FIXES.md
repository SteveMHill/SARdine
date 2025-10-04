# Terrain Correction: Complete Fix List & Code Cleanup

## Executive Summary

Based on comprehensive analysis of `terrain_correction.rs` (8,390 lines), this document provides:
1. **9 Critical Fixes** for correctness (time bases, precision, DEM indexing, orbit)
2. **Removal of 12+ duplicated/heuristic functions** (reduce from ~8.4K to ~6.5K lines)
3. **Logging hygiene** (demote 19+ diagnostic log::error! calls)

---

## Part A: CRITICAL CORRECTNESS FIXES (Must Apply)

### Fix 1: Preserve Sub-Pixel Precision End-to-End

**Location:** Line ~3161  
**Problem:** `fast_range_doppler_calculation` wrapper casts `(usize, usize)` back to `(f64, f64)`, quantizing before resampling  
**Fix:**

```rust
fn fast_range_doppler_calculation(
    &self,
    lat: f64,
    lon: f64,
    elevation: f32,
    orbit_data: &OrbitData,
    params: &RangeDopplerParams,
    _sar_nrows: usize,
    _sar_ncols: usize,
) -> Option<(f64, f64)> {
    // FIXED: Return floats directly, no premature rounding
    self.scientific_range_doppler_transformation(
        lat, lon, elevation as f64, orbit_data, params,
    )
}
```

**Impact:** Prevents 0.5-pixel quantization error before bilinear interpolation

---

### Fix 2: Use Actual Azimuth Line Time (Not Just PRF)

**Locations:**
- `doppler_to_azimuth_pixel_fast` (Line ~3250)
- `latlon_to_sar_pixel_direct` (Line ~5023)
- `azimuth_time_to_pixel` (Line ~6130)
- Any azimuth time → pixel conversion

**Problem:** Using `1.0 / PRF` assumes uniform line spacing; TOPS burst-merged products need `azimuth_time_interval`  
**Fix:**

```rust
fn doppler_to_azimuth_pixel_fast(
    &self,
    _doppler_freq: f64,
    azimuth_time: f64,
    params: &RangeDopplerParams,
) -> SarResult<f64> {
    // FIXED: Use actual line time interval, fallback to PRF only if missing
    let mut dt = params.azimuth_time_interval;
    if !dt.is_finite() || dt <= 0.0 {
        dt = 1.0 / params.prf;
    }
    Ok(azimuth_time / dt)
}
```

**In `latlon_to_sar_pixel_direct`** (Line ~5023):

```rust
// FIXED: Use actual line time
let mut dt_az = params.azimuth_time_interval;
if !dt_az.is_finite() || dt_az <= 0.0 {
    dt_az = 1.0 / params.prf;
}

// Convert orbit-relative → absolute → product-relative
let orbit_ref_epoch = crate::types::datetime_to_utc_seconds(orbit_data.reference_time);
let zero_doppler_abs = zero_doppler_time + orbit_ref_epoch;
let azimuth_pixel = (zero_doppler_abs - params.product_start_time_abs) / dt_az;
```

**Impact:** Fixes off-by-many-lines errors in TOPS/IW products

---

### Fix 3: Remove Accidental EPSG:4326 Branch Inside `create_projected_grid`

**Location:** Line ~4078 (inside `create_projected_grid`)  
**Problem:** Already routed via `create_output_grid` → `create_geographic_grid` for 4326; internal 4326 branch does degree math in "projected" code path  
**Fix:**

```rust
fn create_projected_grid(
    &self,
    bounds: &BoundingBox,
) -> SarResult<(usize, usize, GeoTransform)> {
    // Convert geographic bounds to projected coordinates
    let (proj_min_x, proj_min_y) =
        self.geographic_to_projected(bounds.min_lon, bounds.min_lat)?;
    let (proj_max_x, proj_max_y) =
        self.geographic_to_projected(bounds.max_lon, bounds.max_lat)?;

    let x_extent = proj_max_x - proj_min_x;
    let y_extent = proj_max_y - proj_min_y;

    // Grid dimensions in meters (NO degree logic here)
    let original_width  = (x_extent / self.output_spacing).ceil() as usize;
    let original_height = (y_extent / self.output_spacing).ceil() as usize;

    let max_dimension = self.config.max_output_dimension;
    let width  = original_width.min(max_dimension);
    let height = original_height.min(max_dimension);

    if width != original_width || height != original_height {
        log::warn!("Output dimensions clamped: {}x{} -> {}x{}", 
                   original_width, original_height, width, height);
    }

    // Geotransform in meters (north-up)
    let transform = GeoTransform {
        top_left_x: proj_min_x,
        pixel_width:  self.output_spacing,
        rotation_x:   0.0,
        top_left_y:   proj_max_y,
        rotation_y:   0.0,
        pixel_height: -self.output_spacing, // north-up
    };

    Ok((width, height, transform))
}
```

**Remove:** Lines ~4140-4160 (the internal `if self.output_crs == 4326` block with degree calculations)

**Impact:** Prevents mixed meter/degree math in UTM path

---

### Fix 4: DEM Indexing with Negative `pixel_height`

**Location:** `dem_lookup_with_indices` (Line ~4491 or similar, search for "inv_pixel_height")  
**Problem:** Uses `1.0 / pixel_height` and casts to `usize` before clamp; breaks for north-up rasters (negative height)  
**Fix:**

```rust
fn dem_lookup_with_indices(&self, dem_x_coord: f64, dem_y_coord: f64) -> Option<DemLookupSample> {
    // FIXED: Division + floor() for both positive and negative pixel_height
    let dem_x = (dem_x_coord - self.dem_transform.top_left_x) / self.dem_transform.pixel_width;
    let dem_y = (dem_y_coord - self.dem_transform.top_left_y) / self.dem_transform.pixel_height;

    let dem_col_i = dem_x.floor() as isize;
    let dem_row_i = dem_y.floor() as isize;

    let (dem_h, dem_w) = self.dem.dim();

    // Clamp before casting to usize
    if dem_row_i < 0 || dem_col_i < 0 || 
       dem_row_i as usize >= dem_h - 1 || dem_col_i as usize >= dem_w - 1 {
        return None;
    }

    let dem_row = dem_row_i as usize;
    let dem_col = dem_col_i as usize;

    // ... rest of bilinear logic
}
```

**Impact:** Fixes coordinate inversions/NaN seams for typical north-up DEMs

---

### Fix 5: Orbit Interpolation Coverage & Finiteness Guard

**Location:** `scientific_orbit_interpolation` (Line ~3101)  
**Current:** Already enforces cubic-only ✅  
**Add:**

```rust
fn scientific_orbit_interpolation(
    &self,
    orbit_data: &OrbitData,
    time_seconds: f64,
) -> SarResult<(Position3D, Velocity3D)> {
    use crate::io::orbit::OrbitReader;
    use chrono::DateTime;

    // NEW: Check temporal coverage
    let first = crate::types::datetime_to_utc_seconds(
        orbit_data.state_vectors.first().unwrap().time
    );
    let last = crate::types::datetime_to_utc_seconds(
        orbit_data.state_vectors.last().unwrap().time
    );
    if time_seconds < first - 5.0 || time_seconds > last + 5.0 {
        return Err(SarError::DataProcessingError(format!(
            "Interpolation time {:.3}s outside orbit coverage [{:.3}, {:.3}]",
            time_seconds, first, last
        )));
    }

    // Convert to DateTime
    let time = DateTime::from_timestamp(
        time_seconds as i64, 
        (time_seconds.fract() * 1e9) as u32
    ).ok_or_else(|| {
        SarError::DataProcessingError(format!("Invalid timestamp: {}", time_seconds))
    })?;

    // Cubic interpolation (already enforced)
    if orbit_data.state_vectors.len() >= 4 {
        let position = OrbitReader::interpolate_position(orbit_data, time)?;
        let velocity = OrbitReader::interpolate_velocity(orbit_data, time)?;

        // NEW: Finite checks
        if !position[0].is_finite() || !position[1].is_finite() || !position[2].is_finite() {
            return Err(SarError::DataProcessingError(
                "Orbit interpolation returned non-finite position".to_string()
            ));
        }
        if !velocity[0].is_finite() || !velocity[1].is_finite() || !velocity[2].is_finite() {
            return Err(SarError::DataProcessingError(
                "Orbit interpolation returned non-finite velocity".to_string()
            ));
        }

        return Ok((
            Position3D { x: position[0], y: position[1], z: position[2] },
            Velocity3D { x: velocity[0], y: velocity[1], z: velocity[2] },
        ));
    } else {
        return Err(SarError::DataProcessingError(format!(
            "Insufficient orbit state vectors: {} (need ≥4 for cubic spline)",
            orbit_data.state_vectors.len()
        )));
    }
}
```

**Impact:** Prevents extrapolation/NaN returns from orbit interp

---

### Fix 6: Kill Heuristic Azimuth Approximations (or Fence Them)

**Functions to Remove or Gate:**

1. **`latlon_to_sar_pixel_optimized`** (Line ~3175)  
   - Uses `best_state_idx * 6235.0` magic scaling  
   - **Action:** Delete or wrap in `#[cfg(feature = "diag_fallbacks")]` + never call in production

2. **`latlon_to_sar_pixel`** (Line ~5089)  
   - Hardcoded 800–1500 km slant ranges, 80 s window normalization  
   - **Action:** Delete or feature-gate

3. **`interpolate_satellite_state`** (Line ~3665 or similar)  
   - Index-based, not time-based interpolation  
   - **Action:** Delete (already have `scientific_orbit_interpolation`)

4. **`find_azimuth_time_with_lut`** (Line ~3625)  
   - Returns fraction of LUT length, not physical time  
   - **Action:** Delete

**If keeping for diagnostics:**

```rust
/// WARNING: Diagnostic-only fallback. Not used in scientific mode.
/// Uses crude approximations and should not be enabled in production.
#[cfg(feature = "diag_fallbacks")]
fn latlon_to_sar_pixel_optimized(...) { ... }
```

**Production path:** Always use `scientific_range_doppler_transformation` → Newton/secant zero-Doppler

**Impact:** Eliminates non-physical geometric errors

---

### Fix 7: Make `latlon_to_sar_pixel_direct` Consistent with Time Bases

**Location:** Line ~5023  
**Problem:** Multiplies "zero_doppler_time * PRF" without clear epoch handling  
**Fix:** (Already shown in Fix 2; consolidate here)

```rust
fn latlon_to_sar_pixel_direct(
    &self,
    lat: f64,
    lon: f64,
    elevation: f64,
    orbit_data: &OrbitData,
    params: &RangeDopplerParams,
) -> Option<(f64, f64)> {
    let target_ecef = self.latlon_to_ecef(lat, lon, elevation);

    // Zero-Doppler time (orbit-relative)
    let zero_doppler_time = match self.newton_raphson_zero_doppler(&target_ecef, orbit_data, params) {
        Ok(t) => t,
        Err(_) => return None,
    };

    // Convert orbit-relative → absolute → product-relative
    let orbit_ref_epoch = crate::types::datetime_to_utc_seconds(orbit_data.reference_time);
    let zero_doppler_abs = zero_doppler_time + orbit_ref_epoch;

    // Get satellite position at zero-Doppler time
    let (sat_pos, _) = match self.scientific_orbit_interpolation(orbit_data, zero_doppler_abs) {
        Ok(pv) => pv,
        Err(_) => return None,
    };

    // Slant range
    let range_vec = [
        target_ecef[0] - sat_pos.x,
        target_ecef[1] - sat_pos.y,
        target_ecef[2] - sat_pos.z,
    ];
    let slant_range = (range_vec[0].powi(2) + range_vec[1].powi(2) + range_vec[2].powi(2)).sqrt();

    // Range pixel (two-way timing)
    let two_way_travel_time = 2.0 * slant_range / params.speed_of_light;
    let dt_range = 2.0 * params.range_pixel_spacing / params.speed_of_light;
    let range_pixel = (two_way_travel_time - params.slant_range_time) / dt_range;

    // Azimuth pixel using product start and actual line time
    let mut dt_az = params.azimuth_time_interval;
    if !dt_az.is_finite() || dt_az <= 0.0 {
        dt_az = 1.0 / params.prf;
    }
    let azimuth_pixel = (zero_doppler_abs - params.product_start_time_abs) / dt_az;

    if range_pixel.is_finite() && azimuth_pixel.is_finite() {
        Some((range_pixel, azimuth_pixel))
    } else {
        None
    }
}
```

**Impact:** Aligns time-base handling with NR solver

---

### Fix 8: Bilinear at Borders (Optional but Cleaner)

**Location:** `bilinear_interpolate_unified` (Line ~3533)  
**Current:** Falls back to NN on outermost edge  
**Option A (recommended):**

```rust
fn bilinear_interpolate_unified(&self, image: &Array2<f32>, x: f64, y: f64) -> f32 {
    let (height, width) = image.dim();

    if x < 0.0 || y < 0.0 || x >= width as f64 || y >= height as f64 {
        return f32::NAN;
    }

    let x1 = x.floor() as usize;
    let y1 = y.floor() as usize;

    // FIXED: Clamp x2/y2 for true bilinear at borders
    let x2 = (x1 + 1).min(width  - 1);
    let y2 = (y1 + 1).min(height - 1);

    let dx = x - x1 as f64;
    let dy = y - y1 as f64;

    let v11 = image[[y1, x1]];
    let v12 = image[[y2, x1]];
    let v21 = image[[y1, x2]];
    let v22 = image[[y2, x2]];

    // NaN propagation
    if v11.is_nan() || v12.is_nan() || v21.is_nan() || v22.is_nan() {
        return f32::NAN;
    }

    // Bilinear
    let v1 = v11 as f64 + dx * (v21 as f64 - v11 as f64);
    let v2 = v12 as f64 + dx * (v22 as f64 - v12 as f64);
    (v1 + dy * (v2 - v1)) as f32
}
```

**Option B:** Keep existing NN fallback; document it clearly

**Impact:** Prevents 0.5-1 pixel NN artifacts at image edges

---

### Fix 9: Logging Level Hygiene

**Problem:** 19+ diagnostic `log::error!()` calls on happy paths (e.g., "CLAMP TRAP", "RD TRANSFORM START")  
**Fix:** Demote to `debug!` or `trace!`

**Examples:**

- Line ~2389: `log::error!("🔍 RD TRANSFORM START...")` → `log::trace!(...)`
- Line ~2597: `log::error!("🔍 CLAMP DEBUG #1-7...")` → `log::trace!(...)`
- Line ~2838: `log::error!("🔍 ELEVATION BOUNDS...")` → `log::debug!(...)`
- Line ~2956: `log::error!("🔍 DEM STATS...")` → `log::debug!(...)`
- Line ~3021: `log::error!("🔍 All derivative attempts failed...")` → `log::warn!(...)`

**Keep `error!` only for:**
- Return `Err(...)`
- Drop pixel / fail-fast conditions
- User-facing failures

**Impact:** Production logs show only real errors

---

## Part B: REMOVE DUPLICATED/UNNECESSARY CODE

### Group 1: Interpolation Duplicates

**Remove:**
1. **`bilinear_interpolate`** (Line ~3584) — superseded by `bilinear_interpolate_unified`
2. **`bilinear_interpolate_fast`** (Line ~2223) — wrapper that just calls `unified`; inline or delete

**Action:**
- Delete both functions
- All callers should use `bilinear_interpolate_unified` directly
- **Savings:** ~100 lines

---

### Group 2: Heuristic/Geometric Azimuth Fallbacks

**Remove:**
1. **`latlon_to_sar_pixel_optimized`** (Line ~3175) — index scaling + 6235 magic
2. **`latlon_to_sar_pixel`** (Line ~5089) — 800-1500 km / 80 s normalization
3. **`interpolate_satellite_state`** (Line ~3665) — index-based, not time
4. **`find_azimuth_time_with_lut`** (Line ~3625) — LUT fraction

**Action:** Delete all four + associated orbit LUT helpers if unused elsewhere  
**Savings:** ~500-600 lines

---

### Group 3: Row/Tile Chunk Duplicates

**Remove:**
1. **`process_row_chunk_optimized`** (Line ~2229) — just calls `process_tile_chunk_optimized` with full width

**Action:** Delete proxy; callers use tile version directly  
**Savings:** ~10-20 lines

---

### Group 4: Batch DEM Lookup Duplication

**Functions:**
- `get_elevations_batch_optimized` (sorted)
- `get_elevations_batch` (parallel)

**Action:**
- Keep **one**: the sorted/optimized version
- Make `get_elevations_batch` either alias it or delete
- **Savings:** ~50-100 lines

---

### Group 5: Legacy/Stubs

**Remove:**
1. **`geographic_to_utm`** (legacy) → keep only `enhanced_geographic_to_utm`
2. **`future_gdal_transform`** stub (if not TODO)

**Savings:** ~50 lines

---

### Group 6: Pixel Conversion Variants

**Keep one authoritative:**
- `latlon_to_sar_pixel_direct` (after Fix 7 applied)

**Remove:**
- `latlon_to_sar_pixel` (heuristic version)
- `latlon_to_sar_pixel_optimized` (heuristic)
- `latlon_to_sar_pixel_with_lut` (if redundant after scientific path)
- `latlon_to_sar_pixel_bounded` (if just a wrapper)

**Savings:** ~400-500 lines

---

### Group 7: Orbit Lookup Helpers

**If not needed after deleting heuristic functions:**
- `find_nearest_orbit_state_fast` (Line ~3239)
- `compute_spatial_hash` (Line ~3231)
- Any orbit LUT building code (e.g., `build_orbit_lookup_table_optimized`)

**Savings:** ~200-300 lines

---

### Group 8: Doppler Quick Converter

**`doppler_to_azimuth_pixel_fast`** (Line ~3250):
- If only used by heuristic paths → delete
- If used by tie-grid → keep but apply Fix 2 (use `azimuth_time_interval`)

---

## Part C: TIE-POINT GRID SPECIFICS

**No removal needed; fixes:**

1. **Grid stride:** Already clamps to ≥4 ✅
2. **Flagging:** `TIE_FLAG_VALID` logic correct ✅
3. **Interpolation:** Uses proper bilinear ✅

**Optional enhancement:**
- Expose `tie_point_stride` via `TerrainCorrectionConfig` for tuning (already exists)

---

## Part D: VALIDATION & LIMITS

**Already good:**
- `validate_processing_parameters` extracts frequency from metadata ✅
- `max_realistic_azimuth` derived from `total_azimuth_lines` or duration×PRF ✅
- `max_valid_range_pixel` from config ✅

**No changes needed** (unless you want to make `max_valid_range_pixel` auto-populate from annotation)

---

## Part E: HOUSEKEEPING / API POLISH

### 1. Deprecations Cleanup

**In `RangeDopplerParams`:**
- `product_start_time_abs` / `product_stop_time_abs` already marked `#[deprecated]` ✅
- Usage: Convert remaining uses to `orbit_ref_epoch_utc + product_start_rel_s`

### 2. Single Public Elevation Accessor

**Expose:** `get_elevation_at_latlon_fast` (bilinear)  
**Make private:** `get_elevation_at_latlon`, `get_elevation_at_latlon_optimized` (or delete if redundant)

### 3. `map_to_geographic` Return Order

**Already documented:** Returns `(lat, lon)` ✅  
**Action:** Add doc comment if missing

### 4. `create_output_grid` / Pixel Height

**Already correct:** `pixel_height < 0.0` for north-up ✅

---

## Part F: NICE-TO-HAVE (OPTIONAL)

1. **One-time warning if `product_start_time_abs` used:**

```rust
static WARN_ONCE: Once = Once::new();
WARN_ONCE.call_once(|| {
    log::warn!("Using deprecated product_start_time_abs; prefer product_start_rel_s + orbit_ref_epoch_utc");
});
```

2. **Cubic spline note for 4–6 vectors:**

```rust
if orbit_data.state_vectors.len() < 7 {
    log::info!("Cubic spline with {} state vectors – ensure vectors well-sampled to reduce endpoint noise", 
               orbit_data.state_vectors.len());
}
```

---

## SUMMARY: CODE SIZE REDUCTION

| Category | Lines Removed | Lines Fixed | Net Change |
|----------|---------------|-------------|------------|
| Interpolation duplicates | ~100 | — | -100 |
| Heuristic azimuth fallbacks | ~600 | — | -600 |
| Row/tile duplicates | ~20 | — | -20 |
| Batch DEM duplicates | ~100 | — | -100 |
| Legacy/stubs | ~50 | — | -50 |
| Pixel conversion variants | ~500 | ~50 (keep one) | -450 |
| Orbit lookup helpers | ~300 | — | -300 |
| **Subtotal (removals)** | **~1,670** | — | **-1,670** |
| **Fixes (9 items)** | — | ~200 | +200 |
| **Logging hygiene** | — | ~50 | ~0 (change level) |
| **TOTAL NET REDUCTION** | | | **~1,470 lines** |

**Before:** 8,390 lines  
**After:** ~6,920 lines (–17.5%)

---

## ACTION PLAN

### Phase 1: Apply Critical Fixes (Day 1)
1. Fix 1: Sub-pixel precision wrapper
2. Fix 2: Azimuth time interval (3 locations)
3. Fix 4: DEM negative pixel_height
4. Fix 5: Orbit coverage guard
5. Fix 9: Logging hygiene (demote 19 calls)

**Validation:** Run existing test suite; no regressions

### Phase 2: Remove Heuristic Functions (Day 2)
1. Delete `latlon_to_sar_pixel_optimized`, `latlon_to_sar_pixel`, `interpolate_satellite_state`, `find_azimuth_time_with_lut`
2. Delete orbit LUT helpers if orphaned
3. Grep for remaining calls; replace with `scientific_range_doppler_transformation`

**Validation:** Cargo build succeeds; no orphaned references

### Phase 3: Consolidate Duplicates (Day 2)
1. Delete `bilinear_interpolate`, `bilinear_interpolate_fast`
2. Delete `process_row_chunk_optimized` (tile proxy)
3. Consolidate batch DEM lookups
4. Delete legacy `geographic_to_utm`

**Validation:** Test suite passes; performance benchmarks (should improve)

### Phase 4: Remove Projected Grid 4326 Branch (Day 3)
1. Fix 3: Delete internal EPSG:4326 block in `create_projected_grid`
2. Validate UTM processing remains correct

### Phase 5: Polish & Document (Day 3)
1. Fix 7: `latlon_to_sar_pixel_direct` time bases
2. Fix 8: Bilinear border handling (optional)
3. Add doc comments for single elevation accessor
4. Add deprecation one-time warnings

**Final validation:** Full test suite + real Sentinel-1 scene

---

## TESTING STRATEGY

### Unit Tests
- Time base conversions (orbit_ref → absolute → product_rel)
- DEM indexing with negative pixel_height
- Bilinear interpolation at borders

### Integration Tests
- End-to-end terrain correction with TOPS IW product
- Compare output against pre-fix results (should match or improve)
- Check for NaN seams in output (should disappear)

### Performance Tests
- Benchmark before/after code removal (expect 5-10% speedup)
- Memory profile (expect slight reduction)

---

## RISK MITIGATION

1. **Branch before changes:** `git checkout -b terrain-correction-fixes`
2. **Incremental commits:** One fix per commit with descriptive message
3. **Test after each phase:** Don't proceed if tests fail
4. **Keep deleted code in git history:** Can resurrect if needed
5. **Document rationale:** This markdown serves as reference

---

## COMPLETION CHECKLIST

- [ ] Phase 1: Critical fixes applied (5 items)
- [ ] Phase 2: Heuristic functions removed (~600 lines)
- [ ] Phase 3: Duplicates consolidated (~400 lines)
- [ ] Phase 4: Projected grid fix (~50 lines)
- [ ] Phase 5: Polish & docs (~50 lines)
- [ ] All tests passing
- [ ] Performance benchmarks run
- [ ] Real Sentinel-1 validation complete
- [ ] Code review complete
- [ ] Documentation updated

**Estimated total effort:** 3 days  
**Expected improvement:** –17.5% code, +correctness, +maintainability
