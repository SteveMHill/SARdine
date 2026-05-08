#!/usr/bin/env python3
"""Diagnose burst boundary dark seams in debursted/calibrated output."""

import sys
import numpy as np
import rasterio

def analyze_subswath(path, name):
    src = rasterio.open(path)
    print(f"\n{'='*60}")
    print(f"{name}: {src.height} x {src.width}")
    
    # Read a vertical strip at image center
    mid_col = src.width // 2
    data = src.read(1, window=rasterio.windows.Window(mid_col - 5, 0, 10, src.height))
    # Average across the 10 columns for stability
    col_avg = np.nanmean(data, axis=1)
    
    valid = np.isfinite(col_avg) & (col_avg > 0)
    db = np.full_like(col_avg, np.nan)
    db[valid] = 10 * np.log10(col_avg[valid])
    
    print(f"Total lines: {len(db)}, valid: {valid.sum()}")
    
    # Compute running mean (window=50 lines)
    win = 50
    db_clean = db.copy()
    db_clean[~np.isfinite(db_clean)] = np.nan
    
    # Find burst boundaries by looking for periodically spaced power dips
    # Compute local mean in sliding window
    local_means = []
    for i in range(len(db)):
        lo = max(0, i - win//2)
        hi = min(len(db), i + win//2)
        chunk = db[lo:hi]
        v = chunk[np.isfinite(chunk)]
        local_means.append(np.mean(v) if len(v) > 0 else np.nan)
    local_means = np.array(local_means)
    
    # Find significant dips (> 2 dB below local neighborhood)
    dip_rows = []
    for i in range(win, len(db) - win):
        if not np.isfinite(db[i]):
            continue
        # Compare this pixel to surrounding +/- 100 lines (excluding +-20 lines around it)
        surr_lo = db[max(0, i-150):max(0, i-30)]
        surr_hi = db[min(len(db), i+30):min(len(db), i+150)]
        surr = np.concatenate([surr_lo, surr_hi])
        surr = surr[np.isfinite(surr)]
        if len(surr) < 10:
            continue
        surr_mean = np.mean(surr)
        if db[i] < surr_mean - 3.0:  # More than 3 dB below surroundings
            dip_rows.append((i, db[i], surr_mean, surr_mean - db[i]))
    
    # Cluster nearby dips
    if dip_rows:
        clusters = []
        current = [dip_rows[0]]
        for d in dip_rows[1:]:
            if d[0] - current[-1][0] < 30:
                current.append(d)
            else:
                clusters.append(current)
                current = [d]
        clusters.append(current)
        
        print(f"\nFound {len(clusters)} dark bands (>3 dB below surroundings):")
        spacings = []
        prev_center = None
        for cl in clusters:
            rows = [d[0] for d in cl]
            depths = [d[3] for d in cl]
            center = int(np.mean(rows))
            max_depth = max(depths)
            span = max(rows) - min(rows) + 1
            spacing_str = ""
            if prev_center is not None:
                s = center - prev_center
                spacings.append(s)
                spacing_str = f"  (spacing: {s} lines)"
            prev_center = center
            print(f"  Row {center:5d} (width ~{span:3d} lines): max dip = {max_depth:.1f} dB{spacing_str}")
        
        if spacings:
            print(f"\n  Mean burst spacing: {np.mean(spacings):.0f} lines")
            print(f"  Std burst spacing: {np.std(spacings):.0f} lines")
    else:
        print("\nNo significant dark bands found! (threshold: 3 dB)")
    
    # Also check row-by-row for zero/NaN bands
    all_data = src.read(1, window=rasterio.windows.Window(src.width//4, 0, src.width//2, src.height))
    row_valid_frac = np.array([np.sum(np.isfinite(all_data[r]) & (all_data[r] > 0)) / all_data.shape[1] for r in range(all_data.shape[0])])
    low_coverage = np.where(row_valid_frac < 0.5)[0]
    if len(low_coverage) > 0:
        # Cluster
        gaps = np.diff(low_coverage)
        breaks = np.where(gaps > 5)[0]
        segments = np.split(low_coverage, breaks + 1)
        print(f"\n  Found {len(segments)} low-coverage row bands (< 50% valid):")
        for seg in segments[:20]:
            print(f"    Rows {seg[0]}-{seg[-1]} ({len(seg)} rows, coverage {row_valid_frac[seg[0]]:.2f}-{row_valid_frac[seg[-1]]:.2f})")
    
    src.close()

# Analyze per-subswath outputs
base = 'outputs/validation_run/S1A_IW_SLC__1SDV_20201005T170824_20201005T170851_034664_04098A_1E66'
for sw in ['iw1', 'iw2', 'iw3']:
    path = f'{base}_{sw}_calibrated.tif'
    try:
        analyze_subswath(path, sw.upper())
    except Exception as e:
        print(f"{sw}: ERROR: {e}")

# Also check merged
analyze_subswath(f'{base}_merged_calibrated.tif', 'MERGED')
