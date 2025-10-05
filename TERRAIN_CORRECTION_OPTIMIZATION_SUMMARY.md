# Terrain Correction Optimization Summary

## Critical Performance Fixes Applied

### 1. **Epoch Alignment Fix** ✅ (COMPLETE)
**Impact**: Eliminates huge azimuth index errors (~1.7 billion → reasonable values)

**Problem**: 
- Newton-Raphson solver initialized `best_time` using mixed time references
- Used: `product_start_time_abs - orbit_ref_epoch` (mixing absolute Unix timestamps)
- Result: Incorrect initial guess → excessive iterations, extreme clamping needed

**Solution** (`terrain_correction.rs` line 2557):
```rust
// BEFORE (WRONG):
let product_mid_time = params.product_start_time_abs + (params.product_duration / 2.0);
let mut best_time = product_mid_time - orbit_ref_epoch; // Mixed epochs!

// AFTER (CORRECT):
let mut best_time = params.product_start_rel_s + (params.product_duration / 2.0);
// product_start_rel_s is already orbit-relative!
```

**Expected Results**:
- Azimuth indices: 0 to ~num_lines (instead of 1.7 billion)
- Newton-Raphson convergence: 3-5 iterations (instead of 10-50)
- No extreme clamping (25 seconds → reasonable bounds)
- Doppler derivatives correctly map to azimuth positions

---

### 2. **Analytic Doppler Derivative** ✅ (COMPLETE)
**Impact**: 3× reduction in orbit evaluations per iteration

**Problem**:
- Finite differences required 2 extra orbit interpolations per NR iteration
- `doppler(t+Δ)` and `doppler(t-Δ)` → 3 orbit queries total

**Solution** (`terrain_correction.rs` lines 2705-2710):
```rust
// Analytic derivative from geometry:
// ∂f_d/∂t = -2/λ * [ |v⃗|² - (R_dot)² ] / |r⃗|
//
// where:
// - |v⃗|² = sat_vel · sat_vel (velocity magnitude squared)
// - R_dot = (r⃗ · v⃗) / |r⃗| (radial rate)
// - Assumes a⃗ ≈ 0 (satellite in free-fall orbit)

let vel_magnitude_sq = sat_vel.x * sat_vel.x + sat_vel.y * sat_vel.y + sat_vel.z * sat_vel.z;
let doppler_derivative = -2.0 / params.wavelength * 
    (vel_magnitude_sq / range_magnitude - range_rate * range_rate / range_magnitude);
```

**Fallback**: Finite difference only if analytic derivative fails (rare)

**Expected Results**:
- Orbit evaluations: 1 per iteration (down from 3)
- Wall-clock speedup: ~2-3× on Newton-Raphson solve
- Same accuracy as finite differences

---

### 3. **Early Stop Conditions** ✅ (COMPLETE)
**Impact**: Stop unnecessary iterations once converged

**Implementation** (`terrain_correction.rs` lines 2605-2608, 2681-2684):
```rust
const DOPPLER_CONVERGENCE_HZ: f64 = 1e-3;  // 1 mHz Doppler tolerance
const TIME_CONVERGENCE_S: f64 = 1e-4;       // 0.1 ms time step tolerance
const MAX_ITERATIONS: usize = 8;             // Hard cap (3-5 is typical with good guess)

// Check Doppler convergence
if doppler_freq.abs() < DOPPLER_CONVERGENCE_HZ {
    return Ok(azimuth_time);
}

// Check time step convergence
if time_update.abs() < TIME_CONVERGENCE_S {
    return Ok(azimuth_time);
}
```

**Expected Results**:
- Median iterations: ≤5 (down from 10-17)
- Early exit when either Doppler or time step converges
- Reduced wasted computation

---

### 4. **Debug Logging Cleanup** ✅ (COMPLETE)
**Impact**: Eliminate logging overhead in hot loops

**Changes**:
- Removed per-iteration `log::debug!()` calls
- Only log on first iteration (diagnostics)
- Removed redundant clamp debug messages
- Eliminated duplicate trace calls

**Before**: 10+ log calls per pixel in Newton-Raphson loop
**After**: 0-1 log calls per pixel (only on failures)

**Expected**: 5-15% speedup from reduced I/O overhead

---

## Performance Targets (After Phase 1)

| Metric | Before | Target After Phase 1 | Notes |
|--------|--------|----------------------|-------|
| **NR Iterations (median)** | 10-17 | ≤5 | Epoch fix + early stop |
| **NR Iterations (max)** | 50 | 8 | Hard cap |
| **Orbit Evals/Iteration** | 3 | 1 | Analytic derivative |
| **Doppler Convergence** | 1e-3 Hz | 1e-3 Hz | Same |
| **Time Step Convergence** | N/A | 1e-4 s | New early stop |
| **Logging Overhead** | High | Minimal | Debug removed from loops |

---

## Phase 2 Optimizations (TODO - High Impact)

### A. **Orbit Interpolation Caching**
Pre-spline orbit across entire scene, cache coefficients per knot interval:
```rust
struct OrbitSpline {
    knots: Vec<f64>,           // Time points
    position_coeffs: Vec<[f64; 4]>,  // Cubic coeffs per interval [a, b, c, d]
    velocity_coeffs: Vec<[f64; 4]>,
}

// Evaluate: Horner's rule (no re-fitting)
fn evaluate_position(&self, t: f64) -> [f64; 3] {
    let idx = find_interval(t);
    let dt = t - self.knots[idx];
    let [a, b, c, d] = self.position_coeffs[idx];
    // x(t) = a + b*dt + c*dt² + d*dt³ (repeat for y, z)
}
```
**Expected**: 5-10× faster orbit queries

### B. **Grid Seeding (Coarse→Fine)**
Solve NR on coarse grid (every 16×16 pixels), bilinear interpolate for neighbors:
```rust
// Coarse solve: 16×16 spacing
for (i, j) in coarse_grid {
    solution[i][j] = newton_raphson_full(i, j);
}

// Fine refine: 0-2 iterations per pixel
for (i, j) in fine_grid {
    let seed = bilinear_interp(coarse_solution, i, j);
    solution[i][j] = newton_raphson_refine(i, j, seed, max_its=2);
}
```
**Expected**: 3-10× wall-clock speedup

### C. **Batch Processing with SIMD**
Evaluate orbit spline for 4096 times vectorized:
```rust
use std::simd::f64x8;

fn eval_orbit_batch(spline: &OrbitSpline, times: &[f64]) -> Vec<[f64; 3]> {
    times.chunks(8).map(|chunk| {
        // Vectorized Horner eval with f64x8
    }).collect()
}
```
**Expected**: 2-4× with AVX2/AVX-512

---

## Phase 3 Optimizations (TODO - Medium Impact)

### D. **DEM Tile Caching**
- Memory-map DEM, keep 3×3 tile window in RAM
- Precompute `sin φ, cos φ, sin λ, cos λ` per tile corner
- Bilinear DEM interpolation (not bicubic)

### E. **Parallel Tile Processing**
```rust
use rayon::prelude::*;

tiles.par_chunks(512 * 1024).for_each(|tile| {
    // Process independent tiles in parallel
    // Orbit/DEM reads stay hot in L2/L3
});
```

### F. **Compile Optimizations**
```bash
cargo build --release \
    -C opt-level=3 \
    -C target-cpu=native \
    -C codegen-units=1 \
    -C lto=thin
```

---

## Verification Commands

### Test Epoch Fix + Optimizations
```bash
cd /home/datacube/apps/SARdine/SARdine
RUST_LOG=info python3 -m sardine.cli backscatter \
  ../data/S1A_IW_SLC__1SDV_20201230T165244_20201230T165311_035918_0434F0_6788.SAFE \
  /tmp/test_optimized \
  --resolution 30 \
  --multilook 2 2 \
  --polarization VV
```

### Look For in Logs
```
✅ GOOD SIGNS:
- No "🚨 EPOCH MISMATCH" warnings
- azimuth_time_from_start: 0-60s (not 1.7 billion)
- "Newton-Raphson converged in 3 iterations" (not 10+)
- No "clamping time" warnings on every pixel
- Faster overall processing

❌ BAD SIGNS:
- "EPOCH MISMATCH: azimuth_time_from_start=..." >60s
- "Newton-Raphson did not converge within 8 iterations"
- "clamping time 25.666 to ..." (extreme clamping)
```

---

## Implementation Status

| Optimization | Status | Impact | Effort |
|--------------|--------|--------|--------|
| **Epoch alignment** | ✅ DONE | CRITICAL | Low |
| **Analytic derivative** | ✅ DONE | HIGH | Low |
| **Early stop conditions** | ✅ DONE | MEDIUM | Low |
| **Debug logging cleanup** | ✅ DONE | LOW | Low |
| **Orbit spline caching** | ⏳ TODO | HIGH | Medium |
| **Grid seeding** | ⏳ TODO | VERY HIGH | Medium |
| **Batch SIMD processing** | ⏳ TODO | MEDIUM | High |
| **DEM tile caching** | ⏳ TODO | MEDIUM | Medium |
| **Parallel tiles** | ⏳ TODO | HIGH | Low |

---

## References
- Schaub & Junkins (2003), *Analytical Mechanics of Space Systems*
- Curlander & McDonough (1991), *Synthetic Aperture Radar: Systems and Signal Processing*
- ESA Sentinel-1 Product Specification (S1-RS-MDA-52-7440)

---

## Credits
Optimization strategy based on expert SAR processing best practices:
1. Fix time-base alignment (epoch mix-up)
2. Eliminate unnecessary math (analytic derivatives)
3. Early stopping (Doppler + time convergence)
4. Batch work (future: orbit splines, grid seeding, SIMD)
