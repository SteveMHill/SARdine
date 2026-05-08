#!/usr/bin/env python3
"""Diagnose spatial alignment between SARdine and ASF products."""
import rasterio
import numpy as np

asf_path = "/home/datacube/apps/SARdine/data/ASF/S1B_IW_20190123T053348_DVP_RTC10_G_gpufem_6EC9/S1B_IW_20190123T053348_DVP_RTC10_G_gpufem_6EC9_VV.tif"
sardine_path = "outputs/s1b_gamma0_eval/backscatter_VV_linear.tif"

with rasterio.open(asf_path) as a, rasterio.open(sardine_path) as s:
    print("ASF transform:", a.transform)
    print("SARdine transform:", s.transform)
    print()

    # Center of SARdine extent
    cx = (s.bounds.left + s.bounds.right) / 2
    cy = (s.bounds.bottom + s.bounds.top) / 2
    print(f"Center point: ({cx:.0f}, {cy:.0f})")

    a_row, a_col = a.index(cx, cy)
    s_row, s_col = s.index(cx, cy)
    a_data = a.read(1)
    s_data = s.read(1)
    print(f"ASF pixel [{a_row},{a_col}]: {a_data[a_row,a_col]:.6f}  dB={10*np.log10(max(a_data[a_row,a_col],1e-10)):.2f}")
    print(f"SARdine pixel [{s_row},{s_col}]: {s_data[s_row,s_col]:.6f}  dB={10*np.log10(max(s_data[s_row,s_col],1e-10)):.2f}")

    print("\nASF 5x5 patch (dB):")
    for r in range(a_row-2, a_row+3):
        print("  " + " ".join(f"{10*np.log10(max(a_data[r,c],1e-10)):7.2f}" for c in range(a_col-2, a_col+3)))
    print("SARdine 5x5 patch (dB):")
    for r in range(s_row-2, s_row+3):
        print("  " + " ".join(f"{10*np.log10(max(s_data[r,c],1e-10)):7.2f}" for c in range(s_col-2, s_col+3)))

    s_flat = s_data[s_data > 0]
    a_flat = a_data[a_data > 0]
    print(f"\nSARdine range: {s_flat.min():.6f} to {s_flat.max():.4f}, P1={np.percentile(s_flat,1):.6f} P50={np.percentile(s_flat,50):.6f} P99={np.percentile(s_flat,99):.4f}")
    print(f"ASF range:     {a_flat.min():.8f} to {a_flat.max():.4f}, P1={np.percentile(a_flat,1):.6f} P50={np.percentile(a_flat,50):.6f} P99={np.percentile(a_flat,99):.4f}")

    # Check pixel-by-pixel correlation at the geocoded centre
    # Use rasterio index to ensure coordinate-based alignment
    print("\n--- Coordinate-based spot checks ---")
    test_points = [
        (cx, cy),
        (cx - 10000, cy),
        (cx + 10000, cy),
        (cx, cy - 10000),
        (cx, cy + 10000),
    ]
    for px, py in test_points:
        ar, ac = a.index(px, py)
        sr, sc = s.index(px, py)
        av = a_data[ar, ac] if 0 <= ar < a_data.shape[0] and 0 <= ac < a_data.shape[1] else 0
        sv = s_data[sr, sc] if 0 <= sr < s_data.shape[0] and 0 <= sc < s_data.shape[1] else 0
        a_db = 10*np.log10(max(av,1e-10))
        s_db = 10*np.log10(max(sv,1e-10))
        print(f"  ({px:.0f},{py:.0f}): ASF={a_db:+7.2f} dB  SARdine={s_db:+7.2f} dB  diff={s_db-a_db:+.2f} dB")

    # Check 1000x1000 spatial correlation using rasterio-based alignment
    print("\n--- Small-patch spatial correlation (200x200 around centre) ---")
    sz = 100
    a_patch = a_data[a_row-sz:a_row+sz, a_col-sz:a_col+sz].astype(float)
    s_patch = s_data[s_row-sz:s_row+sz, s_col-sz:s_col+sz].astype(float)
    a_patch[a_patch <= 0] = np.nan
    s_patch[s_patch <= 0] = np.nan
    a_db_p = 10 * np.log10(a_patch)
    s_db_p = 10 * np.log10(s_patch)
    mask = np.isfinite(a_db_p) & np.isfinite(s_db_p)
    if mask.sum() > 100:
        r = np.corrcoef(a_db_p[mask], s_db_p[mask])[0, 1]
        print(f"  Correlation (coordinate-aligned): r = {r:.4f}  (N={mask.sum()})")
    else:
        print(f"  Insufficient overlap at centre ({mask.sum()} valid)")

    # Also try phase correlation on a larger patch to estimate true shift
    print("\n--- Phase correlation shift estimate ---")
    sz2 = 500
    # Read patches at same geographic position
    a_p2 = a_data[a_row-sz2:a_row+sz2, a_col-sz2:a_col+sz2].astype(float).copy()
    s_p2 = s_data[s_row-sz2:s_row+sz2, s_col-sz2:s_col+sz2].astype(float).copy()
    a_p2[~np.isfinite(a_p2) | (a_p2 <= 0)] = 0
    s_p2[~np.isfinite(s_p2) | (s_p2 <= 0)] = 0
    # Phase correlation (more robust than cross-correlation)
    from numpy.fft import fft2, ifft2, fftshift
    F = fft2(a_p2)
    G = fft2(s_p2)
    R = F * np.conj(G)
    denom = np.abs(R)
    denom[denom == 0] = 1
    R /= denom
    cc = fftshift(np.real(ifft2(R)))
    peak = np.unravel_index(np.argmax(cc), cc.shape)
    dy = peak[0] - sz2
    dx = peak[1] - sz2
    print(f"  Phase correlation peak: dy={dy} px ({dy*10:.0f}m), dx={dx} px ({dx*10:.0f}m)")
    print(f"  Peak value: {cc[peak[0], peak[1]]:.4f}")

    # Check the SARdine output metadata JSON if exists
    import json, os
    json_path = "outputs/s1b_gamma0_eval/backscatter_VV_linear.json"
    if os.path.exists(json_path):
        with open(json_path) as f:
            meta = json.load(f)
        print(f"\nSARdine metadata JSON keys: {list(meta.keys())}")
        for k in ["calibration_type", "output_scale", "terrain_correction", "rtc_applied"]:
            if k in meta:
                print(f"  {k}: {meta[k]}")
