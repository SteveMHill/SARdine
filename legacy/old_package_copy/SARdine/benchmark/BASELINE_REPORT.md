# SARdine Performance Baseline Report

**Date:** 2026-01-09  
**Scene:** S1B_IW_SLC__1SDV_20190123T053348_20190123T053415_014617_01B3D4_E833.SAFE  
**Configuration:** `--no-geocode --no-terrain-flatten` (focused on core processing)

## Machine Configuration

| Property | Value |
|----------|-------|
| CPU | Intel Xeon Silver 4316 @ 2.30GHz |
| CPU Cores | 80 (logical) |
| RAM | 1.4 TB |
| OS | Linux 5.15.0-157-generic |
| Python | 3.10.12 |
| SARdine | 0.2.1 |
| Rayon Threads | 8 |

## Benchmark Configuration

| Parameter | Value |
|-----------|-------|
| Polarization | VV |
| Multilook | 2×2 |
| Geocode | Disabled |
| Terrain Flatten | Disabled |
| Speckle Filter | enhanced_lee |
| Thread Count | 8 |

## Summary Results

| Metric | Value |
|--------|-------|
| **Total Wall Time** | **348.4s** (~5.8 min) |
| Logged Step Time | 348.4s |
| Unaccounted Time | ~84s (startup, shutdown, I/O) |

## Per-Stage Breakdown

| Step | Name | Duration | % of Total |
|------|------|----------|------------|
| 1 | Read Metadata & Files | 0.1s | 0.0% |
| 2 | Apply Precise Orbit File | 2.0s | 0.6% |
| 3 | IW Split | 47.6s | 13.7% |
| 4 | **Deburst & Radiometric Calibration** | **122.9s** | **35.3%** |
| 5 | Radiometric Calibration | 0.0s | 0.0% |
| 6 | **Expert IW Merge** | **115.1s** | **33.0%** |
| 7 | **Multilooking** | **35.9s** | **10.3%** |
| 8 | Terrain Flattening | 0.0s | 0.0% |
| 9 | Speckle Filtering | 8.8s | 2.5% |
| 10 | Terrain Correction | 0.0s | 0.0% |
| 11 | Mask Invalid Areas | 6.2s | 1.8% |
| 12 | Convert to dB | 9.8s | 2.8% |
| **Total** | - | **348.4s** | **100%** |

## Hotspot Analysis

### Top 3 Hotspots (92% of total time)

```
1. Step 4: Deburst & Radiometric Calibration - 122.9s (35.3%)
2. Step 6: Expert IW Merge                   - 115.1s (33.0%)
3. Step 3: IW Split                          -  47.6s (13.7%)
```

### Detailed Breakdown from Rust Logs

**Deburst + Calibration (per subswath):**
- IW1: deburst 11.4s, calibration 29.7s (total 41.1s)
- IW2: deburst 10.6s, calibration 32.3s (total 42.9s)
- IW3: deburst 9.1s, calibration 14.1s (total 23.2s)

**Calibration Sub-steps (from Rust timing):**
- Step A: complex→power (SIMD): 0.5-0.7s per subswath
- Step B: calibration coefficient prep: 5-7s per subswath ⚠️ **SLOW**
- Step C: thermal noise removal: 0.8-1.5s per subswath
- Step D: calibration apply: 0.3-0.5s per subswath

**Merge:**
- TOPSAR merge core: 97.4s
- Destriping filter: contributes to merge time

**I/O Operations (SLC Read):**
- IW1: 11.1s for 13518×21728 pixels
- IW2: 9.3s for 13626×25594 pixels  
- IW3: 12.9s for 13671×24680 pixels
- Total I/O: ~33s (overlapped with processing)

## Bound Type Analysis

| Hotspot | Likely Bound | Evidence |
|---------|--------------|----------|
| Deburst | Compute + Memory | Complex deramp + FFT operations |
| Calibration Prep (Step B) | I/O + Compute | "calibration coefficient prep" taking 5-7s |
| IW Merge | Memory | Large array manipulation, ~1GB merged data |
| Multilook | Compute | Array averaging, but 35.9s seems high |

## Key Observations

1. **DEM download dominates "IW Split" step**: The 47.6s for Step 3 is actually DEM tile download (~45s) plus geometry validation (~2.6s). This is I/O-bound and can be parallelized or cached.

2. **Calibration prep (Step B) is unexpectedly slow**: 5-7s per subswath just for "coefficient prep" - this should be sub-second for LUT interpolation. The `precompute_separable_calibration_lut_fast` function is already parallelized with Rayon but is doing redundant work.

3. **Merge is NOT parallelized**: The `execute_merge_plan` function processes rows sequentially. 97s for TOPSAR merge on 12471×67735 output (~844M pixels) = ~9M pixels/second. Parallelizing this could give 4-8x speedup.

4. **Operations are sequential**: Each subswath is processed one at a time. Parallel subswath processing could save ~40% (process IW1/IW2/IW3 in parallel).

5. **Destriping filter adds overhead**: 452413 pixels filtered at 400 boundary columns - this seems like an O(n) operation that could be optimized.

## Profiling Deep-Dive

### Step 3 (IW Split) - 47.6s breakdown:
| Operation | Time |
|-----------|------|
| DEM tile download (15 tiles) | ~45s |
| DEM mosaic creation | ~2s |
| Geometry validation | <1s |

This is NOT a compute bottleneck - it's network I/O. Can be eliminated with cached DEM.

### Step 4 (Deburst+Calibration) - 122.9s breakdown:
| Subswath | Deburst | Calibration Prep | Noise | Apply Cal | Total |
|----------|---------|------------------|-------|-----------|-------|
| IW1 | 11.4s | 5.2s | 1.2s | 0.3s | 18.1s |
| IW2 | 10.6s | 5.4s | 0.8s | 0.3s | 17.1s |
| IW3 | 9.1s | 6.7s | 1.5s | 0.5s | 17.8s |
| **Sequential Total** | **31.1s** | **17.3s** | **3.5s** | **1.1s** | **53.0s** |

Gap analysis: 122.9s - 53.0s = **69.9s unaccounted** (likely SLC I/O and overhead)

### Step 6 (Merge) - 115.1s breakdown:
| Operation | Time |
|-----------|------|
| TOPSAR merge core | 97.4s |
| Destriping filter | ~17.7s |
| **Total** | **115.1s** |

Merge is NOT parallelized - `execute_merge_plan` processes rows sequentially.

## Revised Optimization Priorities

| Priority | Optimization | Est. Speedup | Effort | Risk |
|----------|--------------|--------------|--------|------|
| **P0** | Cache DEM tiles | **30-40% of Step 3** | Low | None |
| **P1** | Parallelize merge rows | **50-75% of Step 6** | Medium | Low |
| **P2** | Parallel subswath deburst+calibrate | **30-40% of Step 4** | Medium | Low |
| **P3** | Optimize calibration LUT prep | **15-20% of Step 4** | Low | Low |
| **P4** | Skip destriping if unnecessary | **15% of Step 6** | Low | Low |

## Reproducibility

```bash
cd /home/datacube/apps/SARdine/SARdine
export RAYON_NUM_THREADS=8
export RUST_LOG=info
python3 -m sardine.cli backscatter \
    "/home/datacube/apps/SARdine/data/SLC/S1B_IW_SLC__1SDV_20190123T053348_20190123T053415_014617_01B3D4_E833.SAFE" \
    ./benchmark_output \
    --polarization VV \
    --multilook 2 2 \
    --num-threads 8 \
    --no-geocode \
    --no-terrain-flatten
```

## Next Steps

1. **Profile Step B (calibration prep)** - Why is LUT interpolation taking 5-7s?
2. **Profile merge** - Is it memory bandwidth limited or compute limited?
3. **Implement parallel subswath processing** - Lowest-hanging fruit for ~40% speedup
4. **Create equivalence tests** - Capture golden outputs before any optimization
