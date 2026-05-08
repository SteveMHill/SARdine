#!/usr/bin/env python3
"""Check the geocoded output for remaining artifacts and quality issues."""
import rasterio
import numpy as np
from scipy.ndimage import median_filter

tif = "outputs/full_pipeline/backscatter_VV_final.tif"

with rasterio.open(tif) as ds:
    data = ds.read(1)
    print(f"Size: {ds.width} x {ds.height}")
    print(f"CRS: {ds.crs}")
    print(f"Pixel size: {ds.res}")

# Row-average power profile to detect artifacts
row_power = np.nanmean(data, axis=1)
valid_rows = np.isfinite(row_power)

med = median_filter(row_power[valid_rows], size=51)
deviation = row_power[valid_rows] - med

dips_2db = np.sum(deviation < -2.0)
dips_3db = np.sum(deviation < -3.0)
dips_5db = np.sum(deviation < -5.0)

print(f"\n--- Burst Seam Check ---")
print(f"Total valid rows: {valid_rows.sum()}")
print(f"Rows with >2dB dip: {dips_2db}")
print(f"Rows with >3dB dip: {dips_3db}")
print(f"Rows with >5dB dip: {dips_5db}")
print(f"Row power std: {np.std(deviation):.3f} dB")
print(f"Row power range: {deviation.min():.2f} to {deviation.max():.2f} dB")

# NaN column check
col_valid = np.sum(np.isfinite(data), axis=0)
nan_cols = np.sum(col_valid == 0)
print(f"\n--- Coverage ---")
print(f"Data shape: {data.shape}")
print(f"All-NaN columns: {nan_cols}")

first100 = np.isfinite(data[:100]).mean()
last100 = np.isfinite(data[-100:]).mean()
mid100 = np.isfinite(data[data.shape[0]//2-50:data.shape[0]//2+50]).mean()
print(f"First 100 rows coverage: {first100*100:.1f}%")
print(f"Middle 100 rows coverage: {mid100*100:.1f}%")
print(f"Last 100 rows coverage: {last100*100:.1f}%")

# Overall stats
valid = np.isfinite(data)
vdata = data[valid]
print(f"\n--- Radiometry ---")
print(f"Valid pixels: {valid.sum():,} / {data.size:,} ({100*valid.sum()/data.size:.1f}%)")
print(f"Range: {vdata.min():.2f} to {vdata.max():.2f} dB")
print(f"Mean: {vdata.mean():.2f} +/- {vdata.std():.2f} dB")
print(f"P1/P5/P50/P95/P99: {np.percentile(vdata,1):.2f} / {np.percentile(vdata,5):.2f} / {np.percentile(vdata,50):.2f} / {np.percentile(vdata,95):.2f} / {np.percentile(vdata,99):.2f} dB")

# Check for ASF reference comparison
import os
asf_dir = "../../data/ASF/S1B_IW_20190123T053348_DVP_RTC10_G_gpufem_6EC9"
if os.path.isdir(asf_dir):
    print(f"\n--- ASF Reference Available ---")
    for f in os.listdir(asf_dir):
        if f.endswith('.tif'):
            print(f"  {f}")

# Range-Doppler failure rate
rd_failures = 1.5  # from log
print(f"\n--- Known Issues ---")
print(f"Range-Doppler failures: ~{rd_failures}% (edge pixels outside SAR image)")
print(f"Coverage: 60.9% fill in bounding rectangle (normal for IW parallelogram)")
if dips_2db > 0:
    print(f"WARNING: {dips_2db} rows with >2dB dips - may indicate residual burst seams")
else:
    print("OK: No significant burst seam artifacts detected")
