#!/usr/bin/env python3
"""Template matching to find the actual offset between ASF and SARdine geocoded products."""
import rasterio
import numpy as np
from scipy.signal import fftconvolve

asf_path = "/home/datacube/apps/SARdine/data/ASF/S1B_IW_20190123T053348_DVP_RTC10_G_gpufem_6EC9/S1B_IW_20190123T053348_DVP_RTC10_G_gpufem_6EC9_VV.tif"
sardine_path = "outputs/s1b_gamma0_eval/backscatter_VV_linear.tif"

with rasterio.open(asf_path) as a, rasterio.open(sardine_path) as s:
    # Read center 4000×4000 patches from both using rasterio coordinate mapping
    cx = (s.bounds.left + s.bounds.right) / 2
    cy = (s.bounds.bottom + s.bounds.top) / 2
    
    # Convert dB and standardize for better cross-correlation
    a_data = a.read(1)
    s_data = s.read(1)
    
    # ASF center area (use rasterio index for coordinate mapping)
    ar, ac = a.index(cx, cy)
    sr, sc = s.index(cx, cy)
    
    print(f"Center: ({cx:.0f}, {cy:.0f})")
    print(f"ASF pixel: ({ar}, {ac})")
    print(f"SARdine pixel: ({sr}, {sc})")
    
    # --- Approach 1: Large-scale cross-correlation ---
    # Take a 200×200 template from ASF, search in 2000×2000 area of SARdine
    template_size = 200
    search_size = 2000
    
    # ASF template (centered on coordinate center)
    ts = template_size // 2
    a_tmpl = a_data[ar-ts:ar+ts, ac-ts:ac+ts].astype(float).copy()
    a_tmpl[~np.isfinite(a_tmpl) | (a_tmpl <= 0)] = np.nan
    
    # Convert to dB for better feature matching
    a_tmpl_db = 10 * np.log10(np.clip(a_tmpl, 1e-10, None))
    a_tmpl_db[np.isnan(a_tmpl)] = 0
    a_tmpl_db -= np.mean(a_tmpl_db)
    a_tmpl_db /= np.std(a_tmpl_db) + 1e-10
    
    # SARdine search area
    ss = search_size // 2
    s_search = s_data[max(0,sr-ss):min(s_data.shape[0],sr+ss), 
                       max(0,sc-ss):min(s_data.shape[1],sc+ss)].astype(float).copy()
    s_search[~np.isfinite(s_search) | (s_search <= 0)] = np.nan
    s_search_db = 10 * np.log10(np.clip(s_search, 1e-10, None))
    s_search_db[np.isnan(s_search)] = 0
    s_search_db -= np.mean(s_search_db)
    s_search_db /= np.std(s_search_db) + 1e-10
    
    print(f"\nASF template: {a_tmpl_db.shape}, SARdine search: {s_search_db.shape}")
    
    # Cross-correlate
    cc = fftconvolve(s_search_db, a_tmpl_db[::-1, ::-1], mode='same')
    peak = np.unravel_index(np.argmax(cc), cc.shape)
    peak_val = cc[peak[0], peak[1]]
    
    # The expected match position (if perfectly aligned) would be at the center
    expected_r = ss
    expected_c = ss
    
    shift_r = peak[0] - expected_r
    shift_c = peak[1] - expected_c
    
    print(f"\n--- Template matching results ---")
    print(f"Peak at: ({peak[0]}, {peak[1]}), expected: ({expected_r}, {expected_c})")
    print(f"Shift: dy={shift_r} px ({shift_r*10:.0f}m), dx={shift_c} px ({shift_c*10:.0f}m)")
    print(f"Peak correlation value: {peak_val:.2f}")
    
    # --- Approach 2: Try multiple template positions ---
    print(f"\n--- Multi-point template matching ---")
    offsets = [(-3000, -3000), (-3000, 3000), (3000, -3000), (3000, 3000), (0, 0)]
    
    for dy_off, dx_off in offsets:
        try:
            # ASF template position
            ar2 = ar + dy_off
            ac2 = ac + dx_off
            if ar2-ts < 0 or ar2+ts >= a_data.shape[0] or ac2-ts < 0 or ac2+ts >= a_data.shape[1]:
                continue
                
            a_t = a_data[ar2-ts:ar2+ts, ac2-ts:ac2+ts].astype(float).copy()
            a_t[~np.isfinite(a_t) | (a_t <= 0)] = np.nan
            if np.isnan(a_t).sum() > a_t.size * 0.3:
                continue
            a_t_db = 10 * np.log10(np.clip(a_t, 1e-10, None))
            a_t_db[np.isnan(a_t)] = 0
            a_t_db -= np.mean(a_t_db)
            a_t_db /= np.std(a_t_db) + 1e-10
            
            # SARdine search (larger to find possible shifts)
            sr2 = sr + dy_off
            sc2 = sc + dx_off
            s_s = s_data[max(0,sr2-ss):min(s_data.shape[0],sr2+ss),
                         max(0,sc2-ss):min(s_data.shape[1],sc2+ss)].astype(float).copy()
            s_s[~np.isfinite(s_s) | (s_s <= 0)] = np.nan
            if np.isnan(s_s).sum() > s_s.size * 0.5:
                continue
            s_s_db = 10 * np.log10(np.clip(s_s, 1e-10, None))
            s_s_db[np.isnan(s_s)] = 0
            s_s_db -= np.mean(s_s_db)
            s_s_db /= np.std(s_s_db) + 1e-10
            
            cc2 = fftconvolve(s_s_db, a_t_db[::-1, ::-1], mode='same')
            pk = np.unravel_index(np.argmax(cc2), cc2.shape)
            sy2 = pk[0] - s_s.shape[0]//2
            sx2 = pk[1] - s_s.shape[1]//2
            
            # Map ASF geo coordinate at template center
            geo_x = a.transform.c + a.transform.a * ac2
            geo_y = a.transform.f + a.transform.e * ar2
            
            print(f"  Template@({geo_x:.0f},{geo_y:.0f}): shift=({sx2},{sy2})px, peak={cc2[pk[0],pk[1]]:.2f}")
        except Exception as e:
            print(f"  Template@offset({dy_off},{dx_off}): {e}")
    
    # --- Approach 3: Check if SARdine data might be in radar geometry ---
    # If geocoding just stretched the raw SAR image to fill the bbox, 
    # the row power profile would show a specific pattern
    print("\n--- Check for radar geometry remnants ---")
    # In radar geometry, near range (left) has different characteristics than far range (right)
    # Also, different subswaths would show up as distinct zones
    col_means = []
    for c in range(0, s_data.shape[1], 100):
        col = s_data[:, c]
        valid = (col > 0) & np.isfinite(col)
        if valid.sum() > 100:
            col_means.append(10*np.log10(np.clip(np.mean(col[valid]), 1e-10, None)))
        else:
            col_means.append(np.nan)
    col_means = np.array(col_means)
    valid_cm = np.isfinite(col_means)
    if valid_cm.sum() > 10:
        trend = np.polyfit(np.where(valid_cm)[0], col_means[valid_cm], 1)
        print(f"  Column-wise trend: {trend[0]:.4f} dB/100px (flat → geocoded correctly)")
        print(f"  Column mean range: {np.nanmin(col_means):.1f} to {np.nanmax(col_means):.1f} dB")
        print(f"  Column mean std: {np.nanstd(col_means):.2f} dB")
        
        # In radar geometry, we'd expect a smooth gradient from near to far range
        # In geocoded data, the trend should be small
        if abs(trend[0]) > 0.1:
            print(f"  ⚠️  Strong across-track trend suggests incomplete geocoding")
        else:
            print(f"  ✅ No strong across-track trend")

    # --- Approach 4: Detailed check of transform handling ---
    print("\n--- Transform verification ---")
    # Check both transforms align visually
    print(f"  ASF top-left:     ({a.bounds.left:.2f}, {a.bounds.top:.2f})")
    print(f"  SARdine top-left: ({s.bounds.left:.2f}, {s.bounds.top:.2f})")
    print(f"  Expected offset:  ({s.bounds.left-a.bounds.left:.2f}, {s.bounds.top-a.bounds.top:.2f})")
    print(f"  In pixels:        ({(s.bounds.left-a.bounds.left)/10:.1f}, {(a.bounds.top-s.bounds.top)/10:.1f})")

    # Check if SARdine grid snaps to ASF grid
    x_offset = (s.bounds.left - a.bounds.left) % 10
    y_offset = (s.bounds.top - a.bounds.top) % 10
    print(f"  Grid alignment:   x_mod10={x_offset:.2f}, y_mod10={y_offset:.2f}")
    if x_offset > 0.5 or y_offset > 0.5:
        print(f"  ⚠️  Pixel grids not aligned")
    else:
        print(f"  ✅ Pixel grids approximately aligned")
