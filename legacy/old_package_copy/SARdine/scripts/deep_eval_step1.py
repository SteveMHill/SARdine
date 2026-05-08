#!/usr/bin/env python3
"""
Deep evaluation: Process S1B scene and compare against ASF RTC reference.

ASF reference: S1B_IW_20190123T053348 RTC10 gamma0 VV (linear power, 10m)
SLC input:     S1B_IW_SLC__1SDV_20190123T053348_20190123T053415_014617_01B3D4_E833.SAFE
"""
import os
import sys
import numpy as np

# Step 1: Inspect ASF reference
print("=" * 70)
print("STEP 1: ASF REFERENCE DATA INSPECTION")
print("=" * 70)

import rasterio

asf_dir = "/home/datacube/apps/SARdine/data/ASF/S1B_IW_20190123T053348_DVP_RTC10_G_gpufem_6EC9"
asf_vv = os.path.join(asf_dir, "S1B_IW_20190123T053348_DVP_RTC10_G_gpufem_6EC9_VV.tif")

with rasterio.open(asf_vv) as ds:
    print(f"  File: {os.path.basename(asf_vv)}")
    print(f"  Size: {ds.width} x {ds.height}")
    print(f"  CRS: {ds.crs}")
    print(f"  Pixel size: {ds.res}")
    print(f"  Bounds: {ds.bounds}")
    print(f"  Dtype: {ds.dtypes[0]}")
    print(f"  NoData: {ds.nodata}")
    asf_data = ds.read(1)
    asf_bounds = ds.bounds
    asf_crs = ds.crs
    asf_transform = ds.transform
    asf_res = ds.res

nodata_val = 0.0
valid = (asf_data != nodata_val) & np.isfinite(asf_data) & (asf_data > 0)
vd = asf_data[valid]
print(f"  Valid pixels: {valid.sum():,} / {asf_data.size:,} ({100*valid.sum()/asf_data.size:.1f}%)")
print(f"  Range: {vd.min():.8f} to {vd.max():.4f}")
print(f"  Mean: {vd.mean():.6f}")

# Check if linear power or dB
if vd.mean() < 10 and vd.min() > 0:
    print("  Scale: LINEAR POWER (gamma0)")
    db = 10 * np.log10(vd)
    print(f"  In dB: mean={db.mean():.2f}, std={db.std():.2f}")
    print(f"  P1={np.percentile(db,1):.2f}, P5={np.percentile(db,5):.2f}, P50={np.percentile(db,50):.2f}, P95={np.percentile(db,95):.2f}, P99={np.percentile(db,99):.2f}")
else:
    print("  Scale: dB or unknown")

# Check the README for product details
readme = os.path.join(asf_dir, "S1B_IW_20190123T053348_DVP_RTC10_G_gpufem_6EC9.README.md.txt")
if os.path.exists(readme):
    print(f"\n  README snippet:")
    with open(readme) as f:
        for i, line in enumerate(f):
            if i < 30:
                print(f"    {line.rstrip()}")

# Step 2: Check SLC input
print("\n" + "=" * 70)
print("STEP 2: S1B SLC INPUT CHECK")
print("=" * 70)

slc_path = "/home/datacube/apps/SARdine/data/SLC/S1B_IW_SLC__1SDV_20190123T053348_20190123T053415_014617_01B3D4_E833.SAFE"
if os.path.isdir(slc_path):
    print(f"  SLC directory exists: {os.path.basename(slc_path)}")
    # Check for annotation files
    ann_dir = os.path.join(slc_path, "annotation")
    if os.path.isdir(ann_dir):
        anns = os.listdir(ann_dir)
        print(f"  Annotation files: {len(anns)}")
        for a in sorted(anns)[:6]:
            print(f"    {a}")
    meas_dir = os.path.join(slc_path, "measurement")
    if os.path.isdir(meas_dir):
        meas = os.listdir(meas_dir)
        print(f"  Measurement files: {len(meas)}")
        for m in sorted(meas):
            fpath = os.path.join(meas_dir, m)
            sz = os.path.getsize(fpath)
            print(f"    {m} ({sz/1e9:.2f} GB)")
else:
    print(f"  ERROR: SLC not found at {slc_path}")
    sys.exit(1)

print("\nReady to process. Run the full pipeline on this scene.")
