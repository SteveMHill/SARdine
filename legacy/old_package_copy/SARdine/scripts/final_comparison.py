#!/usr/bin/env python3
"""Final A/B comparison: original cos² blending vs corrected midpoint selection."""

import rasterio, numpy as np

base_old = 'outputs/validation_run/S1A_IW_SLC__1SDV_20201005T170824_20201005T170851_034664_04098A_1E66'
base_new = 'outputs/validation_run_v3/S1A_IW_SLC__1SDV_20201005T170824_20201005T170851_034664_04098A_1E66'

print("A/B COMPARISON: cos² blending vs midpoint selection (corrected)")
print("="*65)

# 1. Overall statistics
for label, base in [("ORIGINAL (cos2)", base_old), ("NEW (midpoint v3)", base_new)]:
    src = rasterio.open(f'{base}_merged_calibrated.tif')
    w = rasterio.windows.Window(src.width//4, src.height//4, src.width//2, src.height//2)
    d = src.read(1, window=w)
    v = d[np.isfinite(d) & (d > 0)]
    db = 10 * np.log10(v)
    print(f"  {label}: mean={db.mean():.2f} dB, med={np.median(db):.2f}, std={db.std():.2f}, valid={100*v.size/d.size:.1f}%")
    src.close()

# 2. Per-subswath burst boundary analysis
print(f"\n{'='*65}")
print("BURST BOUNDARY DIP ANALYSIS (per-subswath)")
print("="*65)

for sw in ['iw1', 'iw2', 'iw3']:
    print(f"\n  {sw.upper()}:")
    for label, base in [("cos2", base_old), ("midpt", base_new)]:
        src = rasterio.open(f'{base}_{sw}_calibrated.tif')
        h, w = src.height, src.width
        
        strip_w = min(1000, w // 3)
        data = src.read(1, window=rasterio.windows.Window(w//2 - strip_w//2, 0, strip_w, h))
        row_power = np.nanmean(data, axis=1)
        valid = np.isfinite(row_power) & (row_power > 0)
        row_db = np.full(h, np.nan)
        row_db[valid] = 10 * np.log10(row_power[valid])
        
        # Envelope with wide kernel
        kernel = 201
        padded = np.pad(row_db, kernel//2, mode='edge')
        padded[~np.isfinite(padded)] = np.nanmedian(row_db[valid])
        env = np.convolve(padded, np.ones(kernel)/kernel, 'valid')[:h]
        dev = row_db - env
        
        severe = np.sum(np.isfinite(dev) & (dev < -3.0))
        moderate = np.sum(np.isfinite(dev) & (dev < -2.0))
        
        # Spectral analysis at burst frequency
        dev_filled = dev.copy()
        dev_filled[~np.isfinite(dev_filled)] = 0
        spectrum = np.abs(np.fft.rfft(dev_filled))
        burst_spacing = h / 9.0
        burst_freq = 1.0 / burst_spacing
        target_idx = int(burst_freq * h)
        lo = max(1, int(target_idx * 0.8))
        hi_ = min(len(spectrum)-1, int(target_idx * 1.2))
        burst_power = np.mean(spectrum[lo:hi_+1]**2)
        broad_lo = max(1, int(3/burst_spacing))
        broad_hi = min(len(spectrum)-1, h//4)
        broad_power = np.mean(spectrum[broad_lo:broad_hi]**2) if broad_hi > broad_lo else 1e-10
        burst_snr = 10*np.log10(burst_power / max(broad_power, 1e-20))
        
        print(f"    {label}: >3dB dips={severe:3d}, >2dB dips={moderate:3d}, "
              f"burst_freq_SNR={burst_snr:+.1f} dB, row_std={np.nanstd(row_db):.2f} dB")
        src.close()

# 3. Direct comparison at known burst boundary position
print(f"\n{'='*65}")
print("DIRECT BURST BOUNDARY POWER COMPARISON (IW2, wider strip)")
print("="*65)

for label, base in [("cos2", base_old), ("midpt", base_new)]:
    src = rasterio.open(f'{base}_iw2_calibrated.tif')
    h, w = src.height, src.width
    
    # Read VERY wide strip for stable estimates 
    data = src.read(1, window=rasterio.windows.Window(w//4, 0, w//2, h))
    row_power = np.nanmean(data, axis=1)
    valid = np.isfinite(row_power) & (row_power > 0)
    row_db = np.full(h, np.nan)
    row_db[valid] = 10 * np.log10(row_power[valid])
    
    burst_spacing = h / 9
    print(f"\n  {label}:")
    dips = []
    for bi in range(1, 9):
        br = int(bi * burst_spacing)
        # Power at boundary (±10 lines) vs surroundings (±75 lines, excl ±20)
        bv = np.nanmean(row_db[max(0,br-10):min(h,br+10)])
        sl = np.nanmean(row_db[max(0,br-75):max(0,br-20)])
        sh = np.nanmean(row_db[min(h,br+20):min(h,br+75)])
        sm = (sl + sh) / 2
        delta = bv - sm
        dips.append(delta)
        print(f"    B{bi} (row {br}): {delta:+.1f} dB")
    
    dips = np.array(dips)
    print(f"    Mean: {np.mean(dips):+.2f} dB, Worst: {np.min(dips):+.1f} dB, "
          f"Dips>1dB: {np.sum(dips < -1.0)}/{len(dips)}")
    src.close()
