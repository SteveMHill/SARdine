"""
Measure the spatial offset between sardine and ASF by cross-correlating
a 512x512 chip from a high-contrast area (near image centre).
Reports offset in pixels and metres.
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

print("Reading sardine ...")
with rasterio.open(SARDINE) as ds:
    sardine_data = ds.read(1)
    sard_t       = ds.transform
    sard_crs     = ds.crs
    sard_shape   = (ds.height, ds.width)

print("Reading + reprojecting ASF ...")
with rasterio.open(ASF_VV) as ds:
    asf_raw = ds.read(1).astype(np.float32)
    asf_raw[asf_raw == 0] = np.nan
    asf_reproj = np.full(sard_shape, np.nan, dtype=np.float32)
    reproject(asf_raw, asf_reproj,
              src_transform=ds.transform, src_crs=ds.crs,
              dst_transform=sard_t, dst_crs=sard_crs,
              resampling=Resampling.bilinear,
              src_nodata=np.nan, dst_nodata=np.nan)

# Try multiple chips at different scene locations to get robust offset
px_deg = abs(sard_t.e)
h, w = sard_shape
results = []

for frac_r, frac_c, label in [
    (0.35, 0.5, "upper-centre"),
    (0.5,  0.5, "centre"),
    (0.65, 0.5, "lower-centre"),
    (0.5,  0.35, "centre-left"),
    (0.5,  0.65, "centre-right"),
]:
    crow = int(h * frac_r)
    ccol = int(w * frac_c)
    R = 256
    if crow - R < 0 or crow + R > h or ccol - R < 0 or ccol + R > w:
        continue
    chip_s = sardine_data[crow-R:crow+R, ccol-R:ccol+R]
    chip_a = asf_reproj  [crow-R:crow+R, ccol-R:ccol+R]
    valid = np.isfinite(chip_s) & np.isfinite(chip_a)
    if valid.mean() < 0.4:
        print(f"  {label}: only {100*valid.mean():.0f}% valid, skipping")
        continue

    # fill NaN with mean, normalise
    cs = np.where(valid, chip_s, float(np.nanmean(chip_s)))
    ca = np.where(valid, chip_a, float(np.nanmean(chip_a)))
    cs = (cs - cs.mean()) / (cs.std() + 1e-9)
    ca = (ca - ca.mean()) / (ca.std() + 1e-9)

    # FFT cross-correlation (128×128 inner chip for speed)
    R2 = 64
    cs2 = cs[R-R2:R+R2, R-R2:R+R2]
    ca2 = ca[R-R2:R+R2, R-R2:R+R2]
    C = np.fft.ifft2(np.fft.fft2(cs2) * np.conj(np.fft.fft2(ca2))).real
    C = np.fft.fftshift(C)
    peak = np.unravel_index(C.argmax(), C.shape)
    dy = peak[0] - R2
    dx = peak[1] - R2

    dy_m = dy * px_deg * 110540.0
    dx_m = dx * px_deg * 111320.0 * np.cos(np.radians(48.5))
    print(f"  {label}: row={dy:+d} col={dx:+d} px  "
          f"({dy_m:+.0f} m N/S, {dx_m:+.0f} m E/W, "
          f"total {np.hypot(dy_m,dx_m):.0f} m)")
    results.append((dy, dx))

if results:
    med_dy = float(np.median([r[0] for r in results]))
    med_dx = float(np.median([r[1] for r in results]))
    dy_m = med_dy * px_deg * 110540.0
    dx_m = med_dx * px_deg * 111320.0 * np.cos(np.radians(48.5))
    print()
    print(f"Median offset: row={med_dy:+.1f} col={med_dx:+.1f} px")
    print(f"  lat: {dy_m:+.1f} m  (sardine {'S' if dy_m>0 else 'N'} of ASF)")
    print(f"  lon: {dx_m:+.1f} m  (sardine {'E' if dx_m>0 else 'W'} of ASF)")
    print(f"  total: {np.hypot(dy_m,dx_m):.1f} m")
