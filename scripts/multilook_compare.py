"""
Multi-look comparison: suppress speckle to isolate true calibration bias.
Block-averages sardine and ASF to 10×10 pixels (~110 m) before comparing.
"""
import numpy as np
import rasterio
from rasterio.warp import reproject, Resampling
import os

SARDINE = os.environ.get(
    "SARDINE",
    "/home/datacube/dev/SARdine/sardine_s1b_tc_db.tiff",
)
ASF_VV = (
    "/home/datacube/dev/SARdine/data/ASF/"
    "S1B_IW_20190123T053348_DVP_RTC10_G_gpufem_6EC9/"
    "S1B_IW_20190123T053348_DVP_RTC10_G_gpufem_6EC9_VV.tif"
)
BLK = 10  # pixels per side for multi-looking

print("Reading sardine ...", flush=True)
with rasterio.open(SARDINE) as ds:
    sard_db    = ds.read(1)
    sard_t     = ds.transform
    sard_crs   = ds.crs
    sard_shape = (ds.height, ds.width)

print("Reading + reprojecting ASF ...", flush=True)
with rasterio.open(ASF_VV) as ds:
    asf_raw = ds.read(1).astype(np.float32)
    asf_raw[asf_raw == 0] = np.nan
    asf_repr = np.full(sard_shape, np.nan, dtype=np.float32)
    reproject(
        asf_raw, asf_repr,
        src_transform=ds.transform, src_crs=ds.crs,
        dst_transform=sard_t, dst_crs=sard_crs,
        resampling=Resampling.bilinear,
        src_nodata=np.nan, dst_nodata=np.nan,
    )

# Convert sardine dB → linear for block-averaging
sard_lin = np.where(np.isfinite(sard_db), 10.0 ** (sard_db / 10.0), np.nan)

print(f"Block-averaging {BLK}×{BLK} ...", flush=True)
h_blk = sard_shape[0] // BLK
w_blk = sard_shape[1] // BLK

# Reshape-and-mean trick (fast, avoids Python loops)
# Trim to exact multiple of BLK
sl = sard_lin[:h_blk*BLK, :w_blk*BLK]
ar = asf_repr [:h_blk*BLK, :w_blk*BLK]

# nan-aware block mean: reshape to (h_blk, BLK, w_blk, BLK) then nanmean
sl_r = sl.reshape(h_blk, BLK, w_blk, BLK)
ar_r = ar.reshape(h_blk, BLK, w_blk, BLK)

import warnings
with warnings.catch_warnings():
    warnings.simplefilter("ignore", RuntimeWarning)
    s_ml = np.nanmean(sl_r, axis=(1, 3))
    a_ml = np.nanmean(ar_r, axis=(1, 3))

    # Count valid pixels per block; require >= 50% valid
    s_cnt = np.sum(np.isfinite(sl_r), axis=(1, 3))
    a_cnt = np.sum(np.isfinite(ar_r), axis=(1, 3))
    min_valid = BLK * BLK // 2
    s_ml[s_cnt < min_valid] = np.nan
    a_ml[a_cnt < min_valid] = np.nan

    ml_joint = np.isfinite(s_ml) & np.isfinite(a_ml) & (a_ml > 0) & (s_ml > 0)
    s_ml_db = np.where(ml_joint, 10.0 * np.log10(s_ml), np.nan)
    a_ml_db = np.where(ml_joint, 10.0 * np.log10(a_ml), np.nan)
    diff_ml  = s_ml_db[ml_joint] - a_ml_db[ml_joint]

print(f"\nMulti-looked ({BLK}×{BLK}, ~{BLK*11:.0f} m) comparison: n={ml_joint.sum():,}")
print(f"  Sardine median: {np.nanmedian(s_ml_db):.3f} dB")
print(f"  ASF    median:  {np.nanmedian(a_ml_db):.3f} dB")
print(f"  Mean bias:      {diff_ml.mean():+.4f} dB")
print(f"  Median bias:    {np.median(diff_ml):+.4f} dB")
print(f"  Std dev:        {diff_ml.std():.4f} dB  (single-look was 5.37 dB)")
print(f"  RMSE:           {np.sqrt((diff_ml**2).mean()):.4f} dB")

ratio = s_ml[ml_joint] / a_ml[ml_joint]
print(f"\n  Linear ratio (sardine/ASF): mean={ratio.mean():.4f}  median={np.median(ratio):.4f}")
print(f"  Ratio in dB: {10.0*np.log10(np.median(ratio)):+.4f} dB")

# Histogram of multi-looked differences
print("\n── Multi-looked bias histogram ────────────────────────────────────")
bins   = np.arange(-3.0, 3.1, 0.5)
counts, edges = np.histogram(diff_ml, bins=bins)
total  = counts.sum()
peak   = counts.max()
for i, c in enumerate(counts):
    bar = '█' * int(40 * c / (peak + 1e-9))
    lo, hi = edges[i], edges[i+1]
    print(f"  [{lo:+.1f},{hi:+.1f})  {100*c/total:5.1f}%  {bar}")

# Percentiles of |error|
print("\n── |error| percentiles (multi-looked) ─────────────────────────────")
for p in [50, 68, 90, 95]:
    print(f"  p{p:2d}: |error| <= {np.percentile(np.abs(diff_ml), p):.3f} dB")
