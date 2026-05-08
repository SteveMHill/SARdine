#!/usr/bin/env python3
"""Precise burst boundary detection by looking for row-to-row power discontinuities."""

import rasterio, numpy as np

base_old = 'outputs/validation_run/S1A_IW_SLC__1SDV_20201005T170824_20201005T170851_034664_04098A_1E66'
base_new = 'outputs/validation_run_v2/S1A_IW_SLC__1SDV_20201005T170824_20201005T170851_034664_04098A_1E66'

def find_burst_boundaries(data_path, label):
    """Find actual burst boundaries and measure power dips."""
    src = rasterio.open(data_path)
    h, w = src.height, src.width
    
    # Read a wide strip
    strip_w = min(500, w // 4)
    data = src.read(1, window=rasterio.windows.Window(w//2 - strip_w//2, 0, strip_w, h))
    row_power = np.nanmean(data, axis=1)
    valid = np.isfinite(row_power) & (row_power > 0)
    row_db = np.full(h, np.nan)
    row_db[valid] = 10 * np.log10(row_power[valid])
    
    # Smooth to get the envelope
    kernel = 201
    padded = np.pad(row_db, kernel//2, mode='edge')
    padded[~np.isfinite(padded)] = np.nanmedian(row_db[valid])
    smoothed = np.convolve(padded, np.ones(kernel)/kernel, 'valid')[:h]
    
    # Deviation from smooth envelope
    deviation = row_db - smoothed
    
    # Find the actual burst scalloping pattern: look for periodic dips
    # with spacing close to expected burst spacing
    expected_spacing = h / 9  # ~1365-1382 lines
    
    # Count rows that deviate more than thresholds
    severe_dip_rows = np.sum(np.isfinite(deviation) & (deviation < -3.0))
    moderate_dip_rows = np.sum(np.isfinite(deviation) & (deviation < -2.0))
    mild_dip_rows = np.sum(np.isfinite(deviation) & (deviation < -1.0))
    
    print(f"\n  {label}")
    print(f"    Rows below envelope by >3 dB: {severe_dip_rows}")
    print(f"    Rows below envelope by >2 dB: {moderate_dip_rows}")
    print(f"    Rows below envelope by >1 dB: {mild_dip_rows}")
    
    src.close()
    return severe_dip_rows, moderate_dip_rows, mild_dip_rows

print("BURST BOUNDARY DARK BAND COMPARISON")
print("="*60)
print("Counting rows that dip below the local smoothed envelope")
print("(lower = better, means fewer dark bands)")

for sw in ['iw1', 'iw2', 'iw3']:
    print(f"\n{'='*40}")
    print(f"  {sw.upper()}")
    print(f"{'='*40}")
    old_s, old_m, old_l = find_burst_boundaries(f'{base_old}_{sw}_calibrated.tif', "BEFORE (cos2)")
    new_s, new_m, new_l = find_burst_boundaries(f'{base_new}_{sw}_calibrated.tif', "AFTER (midpoint)")
    
    print(f"\n  Improvement:")
    for thresh, o, n in [("3 dB", old_s, new_s), ("2 dB", old_m, new_m), ("1 dB", old_l, new_l)]:
        change = n - o
        pct = (change / max(o, 1)) * 100
        print(f"    >{thresh} dip rows: {o} -> {n} ({change:+d}, {pct:+.0f}%)")

# Also check merged
print(f"\n{'='*40}")
print(f"  MERGED")
print(f"{'='*40}")
old_s, old_m, old_l = find_burst_boundaries(f'{base_old}_merged_calibrated.tif', "BEFORE (cos2)")
new_s, new_m, new_l = find_burst_boundaries(f'{base_new}_merged_calibrated.tif', "AFTER (midpoint)")
print(f"\n  Improvement:")
for thresh, o, n in [("3 dB", old_s, new_s), ("2 dB", old_m, new_m), ("1 dB", old_l, new_l)]:
    change = n - o
    pct = (change / max(o, 1)) * 100
    print(f"    >{thresh} dip rows: {o} -> {n} ({change:+d}, {pct:+.0f}%)")
