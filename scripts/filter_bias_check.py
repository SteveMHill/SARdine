"""
Apply a 7x7 boxcar filter to sardine linear output, then compare mean dB with ASF.
If the -1.28 dB bias drops to near zero after filtering, that confirms:
  - the bias is purely ENL/speckle-statistics artifact
  - sardine calibration is correct

Also measures bias as a function of box size to characterise ENL sensitivity.
"""

import numpy as np
import rasterio
from rasterio.warp import reproject, Resampling
from scipy.ndimage import uniform_filter

SARDINE = "sardine_s1b_vv_30threads.tiff"
ASF_VV  = (
    "data/ASF/S1B_IW_20190123T053348_DVP_RTC10_G_gpufem_6EC9/"
    "S1B_IW_20190123T053348_DVP_RTC10_G_gpufem_6EC9_VV.tif"
)

print("Reading sardine (dB) → linear ...", flush=True)
with rasterio.open(SARDINE) as ds:
    sardine_db    = ds.read(1)
    sardine_tf    = ds.transform
    sardine_crs   = ds.crs
    sardine_shape = (ds.height, ds.width)

sardine_linear = np.where(np.isfinite(sardine_db), 10.0 ** (sardine_db / 10.0), 0.0)
sardine_valid  = np.isfinite(sardine_db)

print("Reprojecting ASF ...", flush=True)
with rasterio.open(ASF_VV) as ds:
    asf_raw  = ds.read(1).astype(np.float32)
    asf_raw[asf_raw <= 0] = np.nan
    asf_crs  = ds.crs
    asf_tf   = ds.transform

asf_reproj = np.full(sardine_shape, np.nan, dtype=np.float32)
reproject(
    source=asf_raw,
    destination=asf_reproj,
    src_transform=asf_tf, src_crs=asf_crs,
    dst_transform=sardine_tf, dst_crs=sardine_crs,
    resampling=Resampling.bilinear,
    src_nodata=np.nan, dst_nodata=np.nan,
)

joint_base = sardine_valid & np.isfinite(asf_reproj)
print(f"Base joint pixels: {joint_base.sum():,}")
print()

# Compute ASF dB (filtered product, reference)
asf_db = np.where(np.isfinite(asf_reproj), 10.0 * np.log10(asf_reproj), np.nan)

print("Boxcar filter sizes vs mean dB bias:")
print(f"  {'Filter':>12s}  {'bias(dB-dB)':>12s}  {'bias(lin mean)':>14s}  {'n_valid':>10s}")
print(f"  {'-'*12}  {'-'*12}  {'-'*14}  {'-'*10}")

for size in [1, 3, 5, 7, 9, 11, 15, 21, 31]:
    # Filter sardine LINEAR, then convert to dB
    valid_float = sardine_valid.astype(np.float32)
    filt_sum    = uniform_filter(sardine_linear, size=size, mode='constant', cval=0.0)
    filt_cnt    = uniform_filter(valid_float,    size=size, mode='constant', cval=0.0)

    # Avoid division by zero; only keep pixels where filter window was >50% valid
    filt_valid  = filt_cnt > 0.5
    filt_linear = np.where(filt_valid, filt_sum / np.maximum(filt_cnt, 1e-9), np.nan)
    filt_db     = np.where(filt_valid & (filt_linear > 0),
                           10.0 * np.log10(filt_linear), np.nan)

    joint = filt_valid & np.isfinite(filt_db) & np.isfinite(asf_db)
    n = joint.sum()
    if n < 1000:
        print(f"  {size:>4d}×{size:<4d}      (insufficient overlap: {n})")
        continue

    s = filt_db[joint]
    a = asf_db[joint]
    mean_db_bias = float(np.mean(s - a))

    # Linear-domain mean bias for this filtered sardine
    s_lin = filt_linear[joint]
    a_lin = asf_reproj[joint]
    lin_bias = 10.0 * np.log10(np.mean(s_lin) / np.mean(a_lin))

    print(f"  {size:>4d}×{size:<4d}      {mean_db_bias:+11.3f}  {lin_bias:+13.3f}  {n:>10,}")

print()
print("Note: ASF uses Enhanced Lee 7×7 (adaptive, not boxcar).")
print("      Boxcar 7×7 gives equivalent ENL as lower bound.")
print("      When bias approaches linear-domain bias (~+0.47 dB),")
print("      the ENL explanation is confirmed.")
