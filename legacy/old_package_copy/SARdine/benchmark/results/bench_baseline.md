# SARdine Baseline Benchmark Report

**Generated:** 2026-01-09T15:30:54.769109

## Machine Configuration

| Property | Value |
|----------|-------|
| CPU | Intel(R) Xeon(R) Silver 4316 CPU @ 2.30GHz |
| CPU Cores | 80 |
| RAM | 1452.8 GB |
| OS | Linux 5.15.0-157-generic |
| Python | 3.10.12 |
| SARdine | 0.2.1 |

## Benchmark Configuration

| Parameter | Value |
|-----------|-------|
| Scene | `S1B_IW_SLC__1SDV_20190123T053348_20190123T053415_014617_01B3D4_E833.SAFE` |
| Polarization | VV |
| Multilook | (2, 2) |
| Threads | 8 |
| Geocode | False |
| Terrain Flatten | False |
| Resolution | 10.0 m |

## Summary Results

| Metric | Value |
|--------|-------|
| **Total Wall Time** | **393.70s** |
| I/O Time | 2.20s |
| Compute Time | 311.10s |
| Peak Memory | 14.8 MB |
| Success | ✅ Yes |

## Per-Stage Breakdown

| Step | Name | Duration (s) | % of Total |
|------|------|-------------|------------|
| 1 | Read Metadata & Files | 0.10 | 0.0% |
| 2 | Apply Precise Orbit File | 2.10 | 0.7% |
| 3 | IW Split | 0.00 | 0.0% |
| 4 | Deburst & Radiometric Calibration | 137.00 | 43.7% |
| 5 | Radiometric Calibration | 0.00 | 0.0% |
| 6 | Expert IW Merge | 120.10 | 38.3% |
| 7 | Multilooking | 34.80 | 11.1% |
| 9 | Speckle Filtering | 7.40 | 2.4% |
| 11 | Mask Invalid Areas | 6.00 | 1.9% |
| 12 | Convert to dB | 5.80 | 1.9% |
| **Total** | - | **313.30** | **100%** |

## Hotspot Analysis

Top 5 slowest stages:

1. **Step 4: Deburst & Radiometric Calibration** - 137.00s (34.8%)
2. **Step 6: Expert IW Merge** - 120.10s (30.5%)
3. **Step 7: Multilooking** - 34.80s (8.8%)
4. **Step 9: Speckle Filtering** - 7.40s (1.9%)
5. **Step 11: Mask Invalid Areas** - 6.00s (1.5%)


## Reproducibility Notes

To reproduce this benchmark:

```bash
cd /home/datacube/apps/SARdine/SARdine
export RAYON_NUM_THREADS=8
export RUST_LOG=info
python -m sardine.cli backscatter \
    "/home/datacube/apps/SARdine/data/SLC/S1B_IW_SLC__1SDV_20190123T053348_20190123T053415_014617_01B3D4_E833.SAFE" \
    ./benchmark_output \
    --polarization VV \
    --multilook 2 2 \
    --num-threads 8 \
    --no-geocode \
    --no-terrain-flatten
```

## Raw Data

The complete benchmark data is saved to `benchmark_results.json`.
