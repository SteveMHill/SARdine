#!/usr/bin/env python3
"""Verify if a constant shift can align SARdine with ASF reference."""
import numpy as np
import rasterio
from scipy.ndimage import shift as ndshift
from scipy.signal import fftconvolve

ASF_PATH = "../data/ASF/S1B_IW_20190123T053348_DVP_RTC10_G_gpufem_6EC9/S1B_IW_20190123T053348_DVP_RTC10_G_gpufem_6EC9_VV.tif"
SARDINE_PATH = "outputs/s1b_gamma0_eval/backscatter_VV_linear.tif"

with rasterio.open(ASF_PATH) as asf_ds, rasterio.open(SARDINE_PATH) as sar_ds:
    # Find overlap in geographic coordinates
    ax0, ay0, ax1, ay1 = asf_ds.bounds
    sx0, sy0, sx1, sy1 = sar_ds.bounds
    
    ox0 = max(ax0, sx0)
    oy0 = max(ay0, sy0)
    ox1 = min(ax1, sx1)
    oy1 = min(ay1, sy1)
    
    # Read overlapping windows
    asf_win = rasterio.windows.from_bounds(ox0, oy0, ox1, oy1, asf_ds.transform)
    sar_win = rasterio.windows.from_bounds(ox0, oy0, ox1, oy1, sar_ds.transform)
    
    asf_data = asf_ds.read(1, window=asf_win)
    sar_data = sar_ds.read(1, window=sar_win)

print(f"ASF shape: {asf_data.shape}, SARdine shape: {sar_data.shape}")

# Match shapes
min_h = min(asf_data.shape[0], sar_data.shape[0])
min_w = min(asf_data.shape[1], sar_data.shape[1])
asf_data = asf_data[:min_h, :min_w]
sar_data = sar_data[:min_h, :min_w]

# Convert to dB
asf_db = np.where(asf_data > 0, 10 * np.log10(asf_data), np.nan)
sar_db = np.where(sar_data > 0, 10 * np.log10(sar_data), np.nan)

# Compute correlation in blocks for different shifts
print("\n=== Testing shifts ===")
# Use downsampled data for faster cross-correlation
ds = 4  # downsample factor
asf_small = asf_db[::ds, ::ds]
sar_small = sar_db[::ds, ::ds]

# Replace NaN with 0 for correlation
asf_clean = np.where(np.isfinite(asf_small), asf_small, 0)
sar_clean = np.where(np.isfinite(sar_small), sar_small, 0)
mask = np.isfinite(asf_small) & np.isfinite(sar_small)

# Phase correlation for precise shift
from numpy.fft import fft2, ifft2
f1 = fft2(asf_clean)
f2 = fft2(sar_clean)
cross = f1 * np.conj(f2)
cross /= np.abs(cross) + 1e-10
corr = np.abs(ifft2(cross))

# Find peak
peak = np.unravel_index(np.argmax(corr), corr.shape)
dy, dx = peak[0], peak[1]
if dy > corr.shape[0] // 2:
    dy -= corr.shape[0]
if dx > corr.shape[1] // 2:
    dx -= corr.shape[1]
print(f"Phase correlation peak (downsampled): dy={dy}, dx={dx}")
print(f"Phase correlation peak (full res):    dy={dy*ds}, dx={dx*ds}")
print(f"In meters: dy={dy*ds*10}m, dx={dx*ds*10}m")
print(f"Peak value: {corr.max():.4f}")

# Now test correlation at multiple shift candidates
best_r = -1
best_shift = (0, 0)
shifts_to_test = [
    (0, 0),
    (dy*ds, dx*ds),
    (-814, -227),  # from original deep_comparison
]

# Also do template matching on a central 500x500 block
cy, cx = min_h // 2, min_w // 2
bh, bw = 500, 500
search_r = 1200  # search radius in pixels

asf_template = asf_db[cy-bh//2:cy+bh//2, cx-bw//2:cx+bw//2]
valid_template = np.isfinite(asf_template).sum() / asf_template.size
print(f"\nTemplate from center: {asf_template.shape}, {valid_template:.1%} valid")

if valid_template > 0.5:
    # Search in SARdine image
    sy_lo = max(0, cy - bh//2 - search_r)
    sy_hi = min(min_h, cy + bh//2 + search_r)
    sx_lo = max(0, cx - bw//2 - search_r)
    sx_hi = min(min_w, cx + bw//2 + search_r)
    sar_search = sar_db[sy_lo:sy_hi, sx_lo:sx_hi]
    
    # Clean for correlation
    tmpl = np.where(np.isfinite(asf_template), asf_template, 0)
    srch = np.where(np.isfinite(sar_search), sar_search, 0)
    
    # Normalized cross-correlation via FFT
    tmpl_norm = tmpl - tmpl.mean()
    srch_norm = srch - srch.mean()
    
    corr_map = fftconvolve(srch_norm, tmpl_norm[::-1, ::-1], mode='valid')
    # Normalize
    from scipy.ndimage import uniform_filter
    srch_sq = uniform_filter(srch_norm**2, size=tmpl_norm.shape)
    local_var = srch_sq[bh//2:srch_sq.shape[0]-bh//2+1, bw//2:srch_sq.shape[1]-bw//2+1]
    
    tmpl_energy = np.sqrt(np.sum(tmpl_norm**2))
    
    min_shape = min(corr_map.shape[0], local_var.shape[0]), min(corr_map.shape[1], local_var.shape[1])
    corr_map = corr_map[:min_shape[0], :min_shape[1]]
    local_var = local_var[:min_shape[0], :min_shape[1]]
    
    # Peak location
    peak_idx = np.unravel_index(np.argmax(corr_map), corr_map.shape)
    # Shift relative to expected center position
    expected_y = cy - bh//2 - sy_lo
    expected_x = cx - bw//2 - sx_lo
    shift_y = peak_idx[0] - expected_y
    shift_x = peak_idx[1] - expected_x
    print(f"Template match shift: dy={shift_y}, dx={shift_x}")
    print(f"In meters: dy={shift_y*10}m, dx={shift_x*10}m")
    shifts_to_test.append((shift_y, shift_x))

# Test each shift by computing correlation
print("\n=== Correlation at each shift ===")
for dy_test, dx_test in shifts_to_test:
    # Crop to overlapping region after shift
    if dy_test >= 0:
        a_ylo, s_ylo = dy_test, 0
    else:
        a_ylo, s_ylo = 0, -dy_test
    if dx_test >= 0:
        a_xlo, s_xlo = dx_test, 0
    else:
        a_xlo, s_xlo = 0, -dx_test
    
    h = min(min_h - a_ylo, min_h - s_ylo)
    w = min(min_w - a_xlo, min_w - s_xlo)
    
    if h <= 0 or w <= 0:
        print(f"  shift=({dy_test},{dx_test}): invalid overlap")
        continue
    
    a_crop = asf_db[int(a_ylo):int(a_ylo+h), int(a_xlo):int(a_xlo+w)]
    s_crop = sar_db[int(s_ylo):int(s_ylo+h), int(s_xlo):int(s_xlo+w)]
    
    valid = np.isfinite(a_crop) & np.isfinite(s_crop)
    n_valid = valid.sum()
    if n_valid < 1000:
        print(f"  shift=({dy_test},{dx_test}): too few valid ({n_valid})")
        continue
    
    a_vals = a_crop[valid]
    s_vals = s_crop[valid]
    
    r = np.corrcoef(a_vals, s_vals)[0, 1]
    bias = np.mean(s_vals - a_vals)
    rmse = np.sqrt(np.mean((s_vals - a_vals)**2))
    print(f"  shift=({dy_test:+5d},{dx_test:+5d}): r={r:.4f}, bias={bias:+.2f} dB, RMSE={rmse:.2f} dB, n={n_valid}")
    
    if r > best_r:
        best_r = r
        best_shift = (dy_test, dx_test)

print(f"\nBest shift: {best_shift}, r={best_r:.4f}")

# Fine search around best shift
print("\n=== Fine search around best shift ===")
for ddy in range(-50, 51, 10):
    for ddx in range(-50, 51, 10):
        dy_test = best_shift[0] + ddy
        dx_test = best_shift[1] + ddx
        
        if dy_test >= 0:
            a_ylo, s_ylo = dy_test, 0
        else:
            a_ylo, s_ylo = 0, -dy_test
        if dx_test >= 0:
            a_xlo, s_xlo = dx_test, 0
        else:
            a_xlo, s_xlo = 0, -dx_test
        
        h = min(min_h - a_ylo, min_h - s_ylo)
        w = min(min_w - a_xlo, min_w - s_xlo)
        
        if h <= 0 or w <= 0:
            continue
        
        # Subsample for speed
        a_crop = asf_db[int(a_ylo):int(a_ylo+h):2, int(a_xlo):int(a_xlo+w):2]
        s_crop = sar_db[int(s_ylo):int(s_ylo+h):2, int(s_xlo):int(s_xlo+w):2]
        
        valid = np.isfinite(a_crop) & np.isfinite(s_crop)
        n_valid = valid.sum()
        if n_valid < 1000:
            continue
        
        r = np.corrcoef(a_crop[valid], s_crop[valid])[0, 1]
        if r > best_r:
            best_r = r
            best_shift = (dy_test, dx_test)

print(f"Refined best shift: {best_shift}, r={best_r:.4f}")
print(f"In meters: dy={best_shift[0]*10}m, dx={best_shift[1]*10}m")
