#!/usr/bin/env python3
"""Check if SARdine geocoding is fundamentally correct by:
1. Comparing valid-data mask shapes
2. Looking for distinctive features (water bodies) at same coordinates
3. Checking if SARdine has spatial structure (autocorrelation)
"""
import rasterio
import numpy as np

asf_path = "/home/datacube/apps/SARdine/data/ASF/S1B_IW_20190123T053348_DVP_RTC10_G_gpufem_6EC9/S1B_IW_20190123T053348_DVP_RTC10_G_gpufem_6EC9_VV.tif"
sardine_path = "outputs/s1b_gamma0_eval/backscatter_VV_linear.tif"

with rasterio.open(asf_path) as a, rasterio.open(sardine_path) as s:
    a_data = a.read(1)
    s_data = s.read(1)
    
    # --- Test 1: Valid mask shape correlation ---
    # Downsample to 100m resolution for mask comparison
    print("=== TEST 1: Valid mask correlation (100m blocks) ===")
    # Use the SARdine extent as common grid
    cx = (s.bounds.left + s.bounds.right) / 2
    cy = (s.bounds.bottom + s.bounds.top) / 2
    
    # Create coarse grids at 100m
    step = 10  # 10 pixels = 100m
    grid_h = s_data.shape[0] // step
    grid_w = s_data.shape[1] // step
    
    s_valid_coarse = np.zeros((grid_h, grid_w), dtype=float)
    a_valid_coarse = np.zeros((grid_h, grid_w), dtype=float)
    
    for iy in range(grid_h):
        for ix in range(grid_w):
            # SARdine block validity
            block_s = s_data[iy*step:(iy+1)*step, ix*step:(ix+1)*step]
            s_valid_coarse[iy, ix] = np.mean((block_s > 0) & np.isfinite(block_s))
            
            # Map to ASF coordinates
            # SARdine pixel center to geo coordinates
            geo_x = s.transform.c + s.transform.a * (ix * step + step/2)
            geo_y = s.transform.f + s.transform.e * (iy * step + step/2)
            
            # Geo coordinates to ASF pixel
            a_col = (geo_x - a.transform.c) / a.transform.a
            a_row = (geo_y - a.transform.f) / a.transform.e
            
            a_r0 = int(a_row)
            a_c0 = int(a_col)
            if 0 <= a_r0 < a_data.shape[0]-step and 0 <= a_c0 < a_data.shape[1]-step:
                block_a = a_data[a_r0:a_r0+step, a_c0:a_c0+step]
                a_valid_coarse[iy, ix] = np.mean((block_a > 0) & np.isfinite(block_a))
    
    # Flatten and compare valid fractions
    mask_corr = np.corrcoef(s_valid_coarse.ravel(), a_valid_coarse.ravel())[0, 1]
    print(f"  Coarse valid-mask correlation: r = {mask_corr:.4f}")
    
    # IoU of the main valid area
    s_binary = s_valid_coarse > 0.5
    a_binary = a_valid_coarse > 0.5
    intersection = np.sum(s_binary & a_binary)
    union = np.sum(s_binary | a_binary)
    iou = intersection / union if union > 0 else 0
    print(f"  Valid mask IoU: {iou:.4f}")
    print(f"  SARdine valid blocks: {s_binary.sum()}, ASF valid blocks: {a_binary.sum()}")

    # --- Test 2: Spatial autocorrelation (lag-1) in each product ---
    print("\n=== TEST 2: Spatial autocorrelation (lag-1 pixel) ===")
    # Use center 2000×2000 patch
    sy, sx = s_data.shape[0]//2, s_data.shape[1]//2
    sz = 1000
    s_patch = s_data[sy-sz:sy+sz, sx-sz:sx+sz].astype(float).copy()
    s_patch[s_patch <= 0] = np.nan
    
    # x-lag autocorrelation
    p1 = s_patch[:, :-1]
    p2 = s_patch[:, 1:]
    valid = np.isfinite(p1) & np.isfinite(p2)
    if valid.sum() > 1000:
        r_lag_x = np.corrcoef(p1[valid], p2[valid])[0, 1]
        print(f"  SARdine lag-1 x autocorrelation: r = {r_lag_x:.4f} (expect >0.8 for real terrain)")
    
    # y-lag autocorrelation
    p1 = s_patch[:-1, :]
    p2 = s_patch[1:, :]
    valid = np.isfinite(p1) & np.isfinite(p2)
    if valid.sum() > 1000:
        r_lag_y = np.corrcoef(p1[valid], p2[valid])[0, 1]
        print(f"  SARdine lag-1 y autocorrelation: r = {r_lag_y:.4f} (expect >0.8 for real terrain)")
    
    # --- Test 3: Water body detection at same coordinates ---
    print("\n=== TEST 3: Feature detection (low backscatter regions = water) ===")
    # Find large water body in ASF (block mean < -18 dB in 100m blocks)
    water_threshold_db = -15.0
    water_threshold_lin = 10**(water_threshold_db/10)
    
    water_count = 0
    match_count = 0
    n_tested = 0
    
    # Test 1000 random blocks
    rng = np.random.RandomState(42)
    for _ in range(2000):
        iy = rng.randint(0, grid_h)
        ix = rng.randint(0, grid_w)
        
        if a_valid_coarse[iy, ix] < 0.8 or s_valid_coarse[iy, ix] < 0.8:
            continue
        
        # Get ASF mean backscatter for this block
        geo_x = s.transform.c + s.transform.a * (ix * step + step/2)
        geo_y = s.transform.f + s.transform.e * (iy * step + step/2)
        
        a_col = int((geo_x - a.transform.c) / a.transform.a)
        a_row = int((geo_y - a.transform.f) / a.transform.e)
        
        if 0 <= a_row < a_data.shape[0]-step and 0 <= a_col < a_data.shape[1]-step:
            block_a = a_data[a_row:a_row+step, a_col:a_col+step]
            block_s = s_data[iy*step:(iy+1)*step, ix*step:(ix+1)*step]
            
            va = (block_a > 0) & np.isfinite(block_a)
            vs = (block_s > 0) & np.isfinite(block_s)
            
            if va.sum() > 50 and vs.sum() > 50:
                mean_a = np.mean(block_a[va])
                mean_s = np.mean(block_s[vs])
                n_tested += 1
                
                if mean_a < water_threshold_lin:
                    water_count += 1
                    if mean_s < water_threshold_lin * 3:  # Allow 5dB tolerance
                        match_count += 1
    
    print(f"  Tested blocks: {n_tested}")
    print(f"  ASF water blocks (< {water_threshold_db} dB): {water_count}")
    if water_count > 0:
        print(f"  Water blocks that also appear dark in SARdine: {match_count} ({100*match_count/water_count:.1f}%)")
    
    # --- Test 4: Block-mean correlation (alternative to pixel-level) ---
    print("\n=== TEST 4: Block-mean backscatter correlation (100m blocks) ===")
    block_size = 100
    n = 0
    asf_vals = []
    sar_vals = []
    
    for iy in range(0, s_data.shape[0] - block_size, block_size):
        for ix in range(0, s_data.shape[1] - block_size, block_size):
            # SARdine block
            block_s = s_data[iy:iy+block_size, ix:ix+block_size]
            vs = (block_s > 0) & np.isfinite(block_s)
            if vs.sum() < block_size * block_size * 0.5:
                continue
            
            # Map center of SARdine block to ASF coordinates
            geo_x = s.transform.c + s.transform.a * (ix + block_size/2)
            geo_y = s.transform.f + s.transform.e * (iy + block_size/2)
            a_col = int((geo_x - a.transform.c) / a.transform.a)
            a_row = int((geo_y - a.transform.f) / a.transform.e)
            
            if 0 <= a_row < a_data.shape[0]-block_size and 0 <= a_col < a_data.shape[1]-block_size:
                block_a = a_data[a_row:a_row+block_size, a_col:a_col+block_size]
                va = (block_a > 0) & np.isfinite(block_a)
                if va.sum() < block_size * block_size * 0.5:
                    continue
                
                asf_vals.append(10*np.log10(np.clip(np.mean(block_a[va]), 1e-10, None)))
                sar_vals.append(10*np.log10(np.clip(np.mean(block_s[vs]), 1e-10, None)))
                n += 1
    
    if n > 10:
        asf_arr = np.array(asf_vals)
        sar_arr = np.array(sar_vals)
        r = np.corrcoef(asf_arr, sar_arr)[0, 1]
        bias = np.mean(sar_arr - asf_arr)
        print(f"  Valid blocks: {n}")
        print(f"  Block correlation: r = {r:.4f}")
        print(f"  Block bias: {bias:+.2f} dB")
        print(f"  ASF block range: {asf_arr.min():.1f} to {asf_arr.max():.1f} dB")  
        print(f"  SARdine block range: {sar_arr.min():.1f} to {sar_arr.max():.1f} dB")
    else:
        print(f"  Only {n} valid blocks found")

    # --- Test 5: Check SARdine output for geocoding artifacts ---
    print("\n=== TEST 5: Geocoding artifact check ===")
    # Check for grid/tile boundaries in the SARdine output
    # Look for discontinuities every 64/128/256 pixels (common tile sizes)
    for tile_size in [64, 128, 256, 512]:
        # Compare row means at tile boundaries vs. elsewhere
        row_means = np.nanmean(np.where(s_data > 0, s_data, np.nan), axis=1)
        valid_means = ~np.isnan(row_means)
        if valid_means.sum() < 100:
            continue
        boundary_rows = np.arange(tile_size, s_data.shape[0], tile_size)
        boundary_rows = boundary_rows[boundary_rows < len(row_means)]
        boundary_rows = boundary_rows[valid_means[boundary_rows]]
        
        if len(boundary_rows) > 5:
            boundary_vals = row_means[boundary_rows]
            all_vals = row_means[valid_means]
            boundary_mean = np.mean(boundary_vals)
            all_mean = np.mean(all_vals)
            boundary_std = np.std(boundary_vals)
            all_std = np.std(all_vals)
            diff = abs(boundary_mean - all_mean) / (all_std + 1e-10)
            if diff > 0.5:
                print(f"  ⚠️  Tile-{tile_size} boundaries differ: {diff:.2f} sigma")
            else:
                print(f"  ✅ Tile-{tile_size}: no boundary artifacts")
        
    print("\nDone.")
