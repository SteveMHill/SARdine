"""Compare sardine linear sigma0 vs ASF linear gamma0 at a specific flat pixel."""
import math
import numpy as np
import rasterio
from rasterio.warp import transform as rio_transform

NOFLATTEN = "/home/datacube/dev/SARdine/sardine_s1b_vv_noflatten.tiff"
ASF_VV = (
    "/home/datacube/dev/SARdine/data/ASF/"
    "S1B_IW_20190123T053348_DVP_RTC10_G_gpufem_6EC9/"
    "S1B_IW_20190123T053348_DVP_RTC10_G_gpufem_6EC9_VV.tif"
)
ASF_INC = (
    "/home/datacube/dev/SARdine/data/ASF/"
    "S1B_IW_20190123T053348_DVP_RTC10_G_gpufem_6EC9/"
    "S1B_IW_20190123T053348_DVP_RTC10_G_gpufem_6EC9_inc_map.tif"
)
LAT, LON = 49.3, 9.5
HALF = 25

print(f"Comparing at lat={LAT}, lon={LON}, ±{HALF} pixels", flush=True)

with rasterio.open(NOFLATTEN) as ds:
    row, col = ds.index(LON, LAT)
    w = rasterio.windows.Window(col - HALF, row - HALF, HALF * 2, HALF * 2)
    s_db = ds.read(1, window=w).astype(np.float32)

s_lin = 10 ** (s_db / 10)
valid_s = s_lin[np.isfinite(s_lin) & (s_lin > 0)]
s_med = float(np.median(valid_s))
print(f"  Sardine σ⁰ (no-flatten, linear): median = {s_med:.6f}  n = {len(valid_s)}")

# Reproject WGS84 → UTM32N for ASF lookup
xs, ys = rio_transform("EPSG:4326", "EPSG:32632", [LON], [LAT])
x_utm, y_utm = xs[0], ys[0]

with rasterio.open(ASF_VV) as ds:
    row, col = ds.index(x_utm, y_utm)
    w = rasterio.windows.Window(col - HALF, row - HALF, HALF * 2, HALF * 2)
    a_lin = ds.read(1, window=w).astype(np.float32)

valid_a = a_lin[(a_lin > 0) & np.isfinite(a_lin)]
a_med = float(np.median(valid_a))
print(f"  ASF γ⁰ (linear):                 median = {a_med:.6f}  n = {len(valid_a)}")

# Read local incidence angle at same location
with rasterio.open(ASF_INC) as ds:
    row, col = ds.index(x_utm, y_utm)
    w = rasterio.windows.Window(col - HALF, row - HALF, HALF * 2, HALF * 2)
    inc_patch = ds.read(1, window=w).astype(np.float32)

valid_inc = inc_patch[(inc_patch > 0) & np.isfinite(inc_patch)]
local_inc_rad = float(np.median(valid_inc))
local_inc_deg = math.degrees(local_inc_rad)
cos_local = math.cos(local_inc_rad)
print(f"  Local incidence (ASF inc_map):   {local_inc_deg:.2f}° → cos = {cos_local:.4f}")
print()

ratio = s_med / a_med
ratio_db = 10 * math.log10(ratio)
print(f"  Measured ratio σ⁰/γ⁰ = {ratio:.4f}  ({ratio_db:+.4f} dB)")
print(f"  cos(local_inc)        = {cos_local:.4f}  ({10*math.log10(cos_local):+.4f} dB)")
print(f"  Residual after cos    = {ratio_db - 10*math.log10(cos_local):+.4f} dB")
print()

# If sardine correctly computes sigma0 and ASF is gamma0:
# expected ratio = sigma0/gamma0 = cos(theta)
# Any deviation from cos(theta) is an upstream calibration error.
expected = cos_local
deviation_db = ratio_db - 10 * math.log10(expected)
print(f"  DIAGNOSIS:")
print(f"    If sardine σ⁰ is correct: expected ratio = cos({local_inc_deg:.1f}°) = {expected:.4f}")
print(f"    Deviation from expected:  {deviation_db:+.4f} dB")
if abs(deviation_db) < 0.2:
    print("    → Sardine σ⁰ is correct. The γ⁰ bias is purely from the cos(θ)")
    print("      flattening, which is already applied. The -1.28 dB mean bias is")
    print("      explained by a scene-wide mean incidence angle mismatch between")
    print("      the spatial grids (sardine vs ASF cover different pixels after reproject).")
else:
    print(f"    → Extra calibration offset of {deviation_db:+.4f} dB unexplained by cos(θ).")
