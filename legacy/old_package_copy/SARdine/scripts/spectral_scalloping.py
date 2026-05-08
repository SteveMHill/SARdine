#!/usr/bin/env python3
"""Spectral analysis: detect periodic scalloping at burst frequency."""

import rasterio, numpy as np

base_old = 'outputs/validation_run/S1A_IW_SLC__1SDV_20201005T170824_20201005T170851_034664_04098A_1E66'
base_new = 'outputs/validation_run_v2/S1A_IW_SLC__1SDV_20201005T170824_20201005T170851_034664_04098A_1E66'

print("SPECTRAL SCALLOPING ANALYSIS")
print("Looking for periodic power modulation at burst frequency")
print("="*60)

for sw in ['iw1', 'iw2', 'iw3']:
    print(f"\n{sw.upper()}:")
    
    for label, base in [("cos2 blend", base_old), ("midpoint", base_new)]:
        src = rasterio.open(f'{base}_{sw}_calibrated.tif')
        h, w = src.height, src.width
        
        # Average across many range columns for stable azimuth profile
        strip_w = min(1000, w // 3)
        data = src.read(1, window=rasterio.windows.Window(w//2 - strip_w//2, 0, strip_w, h))
        row_power = np.nanmean(data, axis=1)
        valid = np.isfinite(row_power) & (row_power > 0)
        row_db = np.full(h, np.nan)
        row_db[valid] = 10 * np.log10(row_power[valid])
        
        # Fill NaN for FFT
        row_db_filled = row_db.copy()
        nanmask = ~np.isfinite(row_db_filled)
        if nanmask.any():
            row_db_filled[nanmask] = np.nanmedian(row_db)
        
        # Remove trend (detrend)
        x = np.arange(h)
        coeffs = np.polyfit(x, row_db_filled, 2)
        trend = np.polyval(coeffs, x)
        detrended = row_db_filled - trend
        
        # FFT
        spectrum = np.abs(np.fft.rfft(detrended))
        freqs = np.fft.rfftfreq(h)
        
        # Expected burst frequency: 9 bursts over h lines => fundamental at 9/h
        # But we want the modulation frequency which is at the burst spacing
        burst_spacing = h / 9.0
        burst_freq = 1.0 / burst_spacing  # cycles per line
        
        # Look at spectrum near burst frequency and its harmonics
        spec_db = 20 * np.log10(spectrum + 1e-10)
        
        # Find peak near burst frequency (within ±20%)
        target_idx = int(burst_freq * h)
        search_lo = max(1, int(target_idx * 0.8))
        search_hi = min(len(spectrum) - 1, int(target_idx * 1.2))
        
        burst_band_power = np.mean(spectrum[search_lo:search_hi+1]**2)
        
        # Compare to broadband power (excluding DC and very low frequencies)
        broadband_lo = max(1, int(3 / burst_spacing))  # exclude first 3 harmonics
        broadband_hi = min(len(spectrum) - 1, h // 4)
        broadband_power = np.mean(spectrum[broadband_lo:broadband_hi]**2) if broadband_hi > broadband_lo else 1e-10
        
        burst_to_broad = 10 * np.log10(burst_band_power / max(broadband_power, 1e-20))
        
        # Also check the actual peak near burst frequency
        peak_idx = search_lo + np.argmax(spectrum[search_lo:search_hi+1])
        peak_freq = freqs[peak_idx]
        peak_period = 1.0 / peak_freq if peak_freq > 0 else float('inf')
        
        print(f"  {label:15s}: burst_freq_power={burst_to_broad:+.1f} dB above broadband, "
              f"peak_period={peak_period:.0f} lines (expected ~{burst_spacing:.0f})")
        
        src.close()

# Now check the most meaningful metric: direct row-to-row power discontinuity
print(f"\n{'='*60}")
print("ROW-TO-ROW POWER JUMPS")
print("(measures sharp transitions - burst boundary signature)")
print("="*60)

for sw in ['iw2']:  # Focus on IW2 (middle, most representative)
    print(f"\n{sw.upper()}:")
    
    for label, base in [("cos2 blend", base_old), ("midpoint", base_new)]:
        src = rasterio.open(f'{base}_{sw}_calibrated.tif')
        h, w = src.height, src.width
        
        # Wide strip for stability
        strip_w = min(2000, w // 2)
        data = src.read(1, window=rasterio.windows.Window(w//2 - strip_w//2, 0, strip_w, h))
        row_power = np.nanmean(data, axis=1)
        valid = np.isfinite(row_power) & (row_power > 0)
        row_db = np.full(h, np.nan)
        row_db[valid] = 10 * np.log10(row_power[valid])
        
        # Compute row-to-row dB jumps
        jumps = np.abs(np.diff(row_db))
        jumps = jumps[np.isfinite(jumps)]
        
        print(f"  {label:15s}: row-to-row jumps: "
              f"mean={np.mean(jumps):.3f} dB, "
              f"p95={np.percentile(jumps, 95):.3f} dB, "
              f"p99={np.percentile(jumps, 99):.3f} dB, "
              f"max={np.max(jumps):.3f} dB, "
              f"jumps>1dB: {np.sum(jumps > 1)}/{len(jumps)}")
        
        src.close()
