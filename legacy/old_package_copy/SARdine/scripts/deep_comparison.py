#!/usr/bin/env python3
"""
Deep comparison: SARdine output vs ASF GAMMA RTC reference.

Both products:
  - Same scene: S1B_IW_SLC 2019-01-23T05:33:48
  - Same CRS: EPSG:32632
  - Same resolution: 10m
  - gamma0, linear power, VV polarization
  
ASF product: GAMMA software, RTC, gamma0, power, filtered
SARdine:     sardine pipeline, RTC (area projection), gamma0, power, enhanced_lee filtered
"""
import os
import numpy as np
import rasterio
from rasterio.warp import reproject, Resampling
from scipy.ndimage import median_filter

# Paths
asf_path = "/home/datacube/apps/SARdine/data/ASF/S1B_IW_20190123T053348_DVP_RTC10_G_gpufem_6EC9/S1B_IW_20190123T053348_DVP_RTC10_G_gpufem_6EC9_VV.tif"
sardine_path = "/home/datacube/apps/SARdine/SARdine/outputs/s1b_gamma0_eval/backscatter_VV_linear.tif"

print("=" * 70)
print("DEEP COMPARISON: SARdine vs ASF GAMMA RTC")
print("=" * 70)

# ─── 1. Load both products ───
print("\n1. PRODUCT METADATA")
print("-" * 50)

with rasterio.open(asf_path) as ds:
    asf_meta = {
        "width": ds.width, "height": ds.height,
        "crs": str(ds.crs), "res": ds.res,
        "bounds": ds.bounds, "transform": ds.transform,
    }
    asf_data = ds.read(1)
    asf_profile = ds.profile.copy()
    print(f"  ASF:     {ds.width}x{ds.height}, {ds.crs}, {ds.res[0]}m, dtype={ds.dtypes[0]}")
    print(f"           bounds: ({ds.bounds.left:.0f}, {ds.bounds.bottom:.0f}) to ({ds.bounds.right:.0f}, {ds.bounds.top:.0f})")

with rasterio.open(sardine_path) as ds:
    sardine_meta = {
        "width": ds.width, "height": ds.height,
        "crs": str(ds.crs), "res": ds.res,
        "bounds": ds.bounds, "transform": ds.transform,
    }
    sardine_data = ds.read(1)
    sardine_profile = ds.profile.copy()
    print(f"  SARdine: {ds.width}x{ds.height}, {ds.crs}, {ds.res[0]}m, dtype={ds.dtypes[0]}")
    print(f"           bounds: ({ds.bounds.left:.0f}, {ds.bounds.bottom:.0f}) to ({ds.bounds.right:.0f}, {ds.bounds.top:.0f})")

# ─── 2. Co-register to common grid ───
print("\n2. CO-REGISTRATION")
print("-" * 50)

# Find overlapping extent
left = max(asf_meta["bounds"].left, sardine_meta["bounds"].left)
bottom = max(asf_meta["bounds"].bottom, sardine_meta["bounds"].bottom)
right = min(asf_meta["bounds"].right, sardine_meta["bounds"].right)
top = min(asf_meta["bounds"].top, sardine_meta["bounds"].top)

print(f"  Overlap bbox: ({left:.0f}, {bottom:.0f}) to ({right:.0f}, {top:.0f})")
print(f"  Overlap size: {(right-left)/10:.0f} x {(top-bottom)/10:.0f} pixels at 10m")

# Convert ASF bbox to pixel coords
def bbox_to_window(transform, left, bottom, right, top):
    """Convert bbox to rasterio-style (row_off, col_off, height, width)."""
    col_off = int((left - transform.c) / transform.a)
    row_off = int((top - transform.f) / transform.e)
    width = int((right - left) / transform.a)
    height = int((bottom - top) / transform.e)
    return max(0, row_off), max(0, col_off), abs(height), abs(width)

asf_r, asf_c, asf_h, asf_w = bbox_to_window(asf_meta["transform"], left, bottom, right, top)
sar_r, sar_c, sar_h, sar_w = bbox_to_window(sardine_meta["transform"], left, bottom, right, top)

# Use minimum extent
h = min(asf_h, sar_h, asf_data.shape[0] - asf_r, sardine_data.shape[0] - sar_r)
w = min(asf_w, sar_w, asf_data.shape[1] - asf_c, sardine_data.shape[1] - sar_c)

print(f"  ASF window:     row={asf_r}, col={asf_c}, h={h}, w={w}")
print(f"  SARdine window: row={sar_r}, col={sar_c}, h={h}, w={w}")

asf_crop = asf_data[asf_r:asf_r+h, asf_c:asf_c+w]
sar_crop = sardine_data[sar_r:sar_r+h, sar_c:sar_c+w]

# Valid mask: both products have data
asf_valid = (asf_crop > 0) & np.isfinite(asf_crop)
sar_valid = (sar_crop > 0) & np.isfinite(sar_crop)
both_valid = asf_valid & sar_valid

print(f"  Co-registered area: {w}x{h} pixels")
print(f"  ASF valid:     {asf_valid.sum():,} ({100*asf_valid.sum()/(h*w):.1f}%)")
print(f"  SARdine valid: {sar_valid.sum():,} ({100*sar_valid.sum()/(h*w):.1f}%)")
print(f"  Both valid:    {both_valid.sum():,} ({100*both_valid.sum()/(h*w):.1f}%)")

if both_valid.sum() < 1000:
    print("ERROR: Insufficient overlap! Check bounds alignment.")
    exit(1)

# ─── 3. Global radiometric comparison ───
print("\n3. RADIOMETRIC COMPARISON (co-registered overlap)")
print("-" * 50)

asf_v = asf_crop[both_valid]
sar_v = sar_crop[both_valid]

# Convert to dB for comparison
asf_db = 10 * np.log10(np.clip(asf_v, 1e-10, None))
sar_db = 10 * np.log10(np.clip(sar_v, 1e-10, None))

diff_db = sar_db - asf_db  # positive = SARdine brighter

print(f"  ASF gamma0:     mean={asf_db.mean():.2f} dB, std={asf_db.std():.2f} dB")
print(f"  SARdine gamma0: mean={sar_db.mean():.2f} dB, std={sar_db.std():.2f} dB")
print(f"")
print(f"  Difference (SARdine - ASF):")
print(f"    Mean bias:  {diff_db.mean():+.3f} dB")
print(f"    Std:        {diff_db.std():.3f} dB")
print(f"    Median:     {np.median(diff_db):+.3f} dB")
print(f"    P5/P95:     {np.percentile(diff_db,5):+.3f} / {np.percentile(diff_db,95):+.3f} dB")
print(f"    MAE:        {np.abs(diff_db).mean():.3f} dB")
print(f"    RMSE:       {np.sqrt((diff_db**2).mean()):.3f} dB")

# Linear ratio
ratio = sar_v / asf_v
print(f"\n  Linear ratio (SARdine/ASF):")
print(f"    Mean:   {ratio.mean():.4f}")
print(f"    Median: {np.median(ratio):.4f}")
print(f"    Std:    {ratio.std():.4f}")

# ─── 4. Spatial correlation ───
print("\n4. SPATIAL CORRELATION")
print("-" * 50)

# Downsample for correlation (100x100 blocks)
block = 100
bh, bw = h // block, w // block
asf_blocked = np.full((bh, bw), np.nan)
sar_blocked = np.full((bh, bw), np.nan)

for i in range(bh):
    for j in range(bw):
        patch_asf = asf_crop[i*block:(i+1)*block, j*block:(j+1)*block]
        patch_sar = sar_crop[i*block:(i+1)*block, j*block:(j+1)*block]
        mask = (patch_asf > 0) & np.isfinite(patch_asf) & (patch_sar > 0) & np.isfinite(patch_sar)
        if mask.sum() > block * block * 0.5:
            asf_blocked[i, j] = np.mean(patch_asf[mask])
            sar_blocked[i, j] = np.mean(patch_sar[mask])

valid_blocks = np.isfinite(asf_blocked) & np.isfinite(sar_blocked)
if valid_blocks.sum() > 10:
    asf_b = 10 * np.log10(np.clip(asf_blocked[valid_blocks], 1e-10, None))
    sar_b = 10 * np.log10(np.clip(sar_blocked[valid_blocks], 1e-10, None))
    corr = np.corrcoef(asf_b, sar_b)[0, 1]
    print(f"  Block correlation (100x100 px blocks): r = {corr:.6f}")
    print(f"  Valid blocks: {valid_blocks.sum()}")
    block_diff = sar_b - asf_b
    print(f"  Block diff: mean={block_diff.mean():+.3f} dB, std={block_diff.std():.3f} dB")
else:
    print("  Insufficient valid blocks for correlation")

# ─── 5. Histogram comparison ───
print("\n5. HISTOGRAM COMPARISON (dB)")
print("-" * 50)

bins = np.arange(-30, 10, 0.5)
asf_hist, _ = np.histogram(asf_db, bins=bins)
sar_hist, _ = np.histogram(sar_db, bins=bins)

# Normalize
asf_hist_n = asf_hist / asf_hist.sum()
sar_hist_n = sar_hist / sar_hist.sum()

# KL divergence (symmetrized)
eps = 1e-10
kl_ab = np.sum(asf_hist_n * np.log((asf_hist_n + eps) / (sar_hist_n + eps)))
kl_ba = np.sum(sar_hist_n * np.log((sar_hist_n + eps) / (asf_hist_n + eps)))
jsd = 0.5 * kl_ab + 0.5 * kl_ba

# Histogram intersection
hist_intersection = np.minimum(asf_hist_n, sar_hist_n).sum()

print(f"  Histogram intersection: {hist_intersection:.4f} (1.0 = identical)")
print(f"  Jensen-Shannon divergence: {jsd:.6f} (0.0 = identical)")

# Show peak positions
asf_peak = bins[np.argmax(asf_hist_n)]
sar_peak = bins[np.argmax(sar_hist_n)]
print(f"  ASF histogram peak:     {asf_peak:.1f} dB")
print(f"  SARdine histogram peak: {sar_peak:.1f} dB")

# ─── 6. Burst seam check in SARdine ───
print("\n6. BURST SEAM ARTIFACT CHECK (SARdine output)")
print("-" * 50)

row_power_sar = np.nanmean(sar_crop, axis=1)
valid_rows = np.isfinite(row_power_sar) & (row_power_sar > 0)
if valid_rows.sum() > 100:
    rp = 10 * np.log10(np.clip(row_power_sar[valid_rows], 1e-10, None))
    med = median_filter(rp, size=51)
    dev = rp - med
    dips_2 = np.sum(dev < -2.0)
    dips_3 = np.sum(dev < -3.0)
    print(f"  Rows with >2dB dip: {dips_2}")
    print(f"  Rows with >3dB dip: {dips_3}")
    print(f"  Row power std: {np.std(dev):.3f} dB")
    print(f"  Row power range: {dev.min():.2f} to {dev.max():.2f} dB")

# Same check for ASF (should be clean)
row_power_asf = np.nanmean(asf_crop, axis=1)
valid_rows_asf = np.isfinite(row_power_asf) & (row_power_asf > 0)
if valid_rows_asf.sum() > 100:
    rp_a = 10 * np.log10(np.clip(row_power_asf[valid_rows_asf], 1e-10, None))
    med_a = median_filter(rp_a, size=51)
    dev_a = rp_a - med_a
    dips_2a = np.sum(dev_a < -2.0)
    print(f"  ASF reference >2dB dips: {dips_2a} (should be ~0)")

# ─── 7. Spatial shift estimation ───
print("\n7. SPATIAL SHIFT ESTIMATION (cross-correlation)")
print("-" * 50)

# Use central 2000x2000 patch for shift detection
cy, cx = h // 2, w // 2
ps = 1000  # half-size
patch_asf = asf_crop[cy-ps:cy+ps, cx-ps:cx+ps].copy()
patch_sar = sar_crop[cy-ps:cy+ps, cx-ps:cx+ps].copy()

# Replace invalid with 0
patch_asf[~np.isfinite(patch_asf) | (patch_asf <= 0)] = 0
patch_sar[~np.isfinite(patch_sar) | (patch_sar <= 0)] = 0

if patch_asf.sum() > 0 and patch_sar.sum() > 0:
    # Normalize
    patch_asf = (patch_asf - patch_asf.mean()) / (patch_asf.std() + 1e-10)
    patch_sar = (patch_sar - patch_sar.mean()) / (patch_sar.std() + 1e-10)
    
    from numpy.fft import fft2, ifft2, fftshift
    cc = fftshift(np.abs(ifft2(fft2(patch_asf) * np.conj(fft2(patch_sar)))))
    peak = np.unravel_index(np.argmax(cc), cc.shape)
    shift_y = peak[0] - ps
    shift_x = peak[1] - ps
    print(f"  Estimated shift: dy={shift_y} px ({shift_y*10:.0f}m), dx={shift_x} px ({shift_x*10:.0f}m)")
    if abs(shift_y) <= 2 and abs(shift_x) <= 2:
        print(f"  ✅ Sub-pixel or negligible shift — good geolocation accuracy")
    elif abs(shift_y) <= 5 and abs(shift_x) <= 5:
        print(f"  ⚠️  Small shift — acceptable for 10m products")
    else:
        print(f"  ❌ Significant shift — possible geolocation issue")

# ─── 8. Summary ───
print("\n" + "=" * 70)
print("SUMMARY: SARdine vs ASF GAMMA RTC Quality Assessment")
print("=" * 70)

bias = diff_db.mean()
rmse = np.sqrt((diff_db**2).mean())

print(f"  Radiometric bias:       {bias:+.3f} dB {'✅ <0.5dB' if abs(bias) < 0.5 else '⚠️  >0.5dB' if abs(bias) < 1.0 else '❌ >1dB'}")
print(f"  Radiometric RMSE:       {rmse:.3f} dB {'✅ <2dB' if rmse < 2.0 else '⚠️  >2dB' if rmse < 3.0 else '❌ >3dB'}")
print(f"  Spatial correlation:    r={corr:.4f} {'✅ >0.99' if corr > 0.99 else '⚠️  >0.95' if corr > 0.95 else '❌ <0.95'}")
print(f"  Histogram similarity:   {hist_intersection:.4f} {'✅ >0.90' if hist_intersection > 0.90 else '⚠️  >0.80' if hist_intersection > 0.80 else '❌ <0.80'}")
print(f"  Burst seam artifacts:   {dips_2} rows >2dB {'✅ none' if dips_2 == 0 else '⚠️  few' if dips_2 < 10 else '❌ many'}")
shift_ok = abs(shift_y) <= 2 and abs(shift_x) <= 2
print(f"  Geolocation shift:      ({shift_x},{shift_y}) px {'✅ <20m' if shift_ok else '⚠️  >20m'}")

# Overall grade
issues = 0
if abs(bias) >= 1.0: issues += 2
elif abs(bias) >= 0.5: issues += 1
if rmse >= 3.0: issues += 2
elif rmse >= 2.0: issues += 1
if corr < 0.95: issues += 2
elif corr < 0.99: issues += 1
if dips_2 > 10: issues += 1

if issues == 0:
    print(f"\n  OVERALL: ✅ EXCELLENT — SARdine output matches ASF reference closely")
elif issues <= 2:
    print(f"\n  OVERALL: ⚠️  GOOD — Minor differences from ASF reference")
else:
    print(f"\n  OVERALL: ❌ NEEDS WORK — Significant differences from ASF reference")
