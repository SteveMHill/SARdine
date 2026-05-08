#!/usr/bin/env python3
"""Compare burst boundary quality: old (cos2 blending) vs new (midpoint selection)."""

import rasterio, numpy as np

base_old = 'outputs/validation_run/S1A_IW_SLC__1SDV_20201005T170824_20201005T170851_034664_04098A_1E66'
base_new = 'outputs/validation_run_v2/S1A_IW_SLC__1SDV_20201005T170824_20201005T170851_034664_04098A_1E66'

for label, base in [("BEFORE (cos2 blending)", base_old), ("AFTER (midpoint selection)", base_new)]:
    print(f"\n{'='*60}")
    print(f"  {label}")
    print(f"{'='*60}")
    
    for sw in ['iw1', 'iw2', 'iw3']:
        src = rasterio.open(f'{base}_{sw}_calibrated.tif')
        h, w = src.height, src.width
        
        strip_w = 200
        data = src.read(1, window=rasterio.windows.Window(w//2 - strip_w//2, 0, strip_w, h))
        row_power = np.nanmean(data, axis=1)
        valid = np.isfinite(row_power) & (row_power > 0)
        row_db = np.full(h, np.nan)
        row_db[valid] = 10 * np.log10(row_power[valid])
        
        print(f"\n  {sw.upper()}: {h} x {w}")
        print(f"  Row power std: {np.nanstd(row_db):.2f} dB")
        
        burst_spacing = h / 9
        dips = []
        for bi in range(1, 9):
            boundary_row = int(bi * burst_spacing)
            lo = max(0, boundary_row - 75)
            hi = min(h, boundary_row + 75)
            clo = max(0, boundary_row - 10)
            chi = min(h, boundary_row + 10)
            
            bv = np.nanmean(row_db[clo:chi])
            sl = np.nanmean(row_db[lo:clo])
            sh = np.nanmean(row_db[chi:hi])
            sm = (sl + sh) / 2
            delta = bv - sm
            dips.append(delta)
            flag = " <-- DIP!" if delta < -1.0 else ""
            print(f"    Boundary {bi} (row {boundary_row}): delta={delta:+.1f} dB{flag}")
        
        dips = np.array(dips)
        print(f"  Summary: mean={np.mean(dips):+.2f} dB, worst={np.min(dips):+.1f} dB")
        src.close()
