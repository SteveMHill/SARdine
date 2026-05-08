"""Bias profile along a latitude strip to see within-IW vs between-IW variation."""
import numpy as np
import rasterio
import math
from rasterio.warp import transform as rio_transform

NOFLATTEN = "sardine_s1b_vv_noflatten.tiff"
FLATTEN   = "sardine_s1b_vv_30threads.tiff"
ASF_VV = (
    "data/ASF/S1B_IW_20190123T053348_DVP_RTC10_G_gpufem_6EC9/"
    "S1B_IW_20190123T053348_DVP_RTC10_G_gpufem_6EC9_VV.tif"
)
ASF_INC = (
    "data/ASF/S1B_IW_20190123T053348_DVP_RTC10_G_gpufem_6EC9/"
    "S1B_IW_20190123T053348_DVP_RTC10_G_gpufem_6EC9_inc_map.tif"
)

LAT = 49.2
LONS = np.arange(7.5, 11.05, 0.05)
WIN = 30  # half-window in pixels

print(f"{'lon':>6}  {'inc_deg':>7}  {'sard_s0':>9}  {'asf_s0':>9}  "
      f"{'s0_bias_dB':>10}  {'sard_g0':>9}  {'asf_g0':>9}  {'g0_bias_dB':>10}")

for lon in LONS:
    xs, ys = rio_transform("EPSG:4326", "EPSG:32632", [lon], [LAT])
    x_utm, y_utm = xs[0], ys[0]

    # ASF incidence angle
    with rasterio.open(ASF_INC) as ds:
        r, c = ds.index(x_utm, y_utm)
        if r < WIN or c < WIN or r > ds.height - WIN or c > ds.width - WIN:
            continue
        inc_arr = ds.read(1, window=rasterio.windows.Window(c-WIN, r-WIN, 2*WIN, 2*WIN)).astype(np.float32)
    inc_valid = inc_arr[inc_arr > 0]
    if len(inc_valid) == 0:
        continue
    inc = float(np.median(inc_valid))
    cos_inc = math.cos(inc)

    # ASF gamma0
    with rasterio.open(ASF_VV) as ds:
        r2, c2 = ds.index(x_utm, y_utm)
        if r2 < WIN or c2 < WIN or r2 > ds.height - WIN or c2 > ds.width - WIN:
            continue
        a = ds.read(1, window=rasterio.windows.Window(c2-WIN, r2-WIN, 2*WIN, 2*WIN)).astype(np.float32)
    a_valid = a[a > 0]
    if len(a_valid) == 0:
        continue
    asf_g0 = float(np.median(a_valid))
    asf_s0 = asf_g0 * cos_inc

    # Sardine no-flatten sigma0 (dB output → linear)
    with rasterio.open(NOFLATTEN) as ds:
        r3, c3 = ds.index(lon, LAT)
        if r3 < WIN or c3 < WIN or r3 > ds.height - WIN or c3 > ds.width - WIN:
            continue
        s_db = ds.read(1, window=rasterio.windows.Window(c3-WIN, r3-WIN, 2*WIN, 2*WIN)).astype(np.float32)
    s_lin = 10 ** (s_db / 10)
    s_valid = s_lin[np.isfinite(s_lin) & (s_lin > 0)]
    if len(s_valid) == 0:
        continue
    sard_s0 = float(np.median(s_valid))

    # Sardine flatten gamma0 (dB output → linear)
    with rasterio.open(FLATTEN) as ds:
        r4, c4 = ds.index(lon, LAT)
        if r4 < WIN or c4 < WIN or r4 > ds.height - WIN or c4 > ds.width - WIN:
            continue
        g_db = ds.read(1, window=rasterio.windows.Window(c4-WIN, r4-WIN, 2*WIN, 2*WIN)).astype(np.float32)
    g_lin = 10 ** (g_db / 10)
    g_valid = g_lin[np.isfinite(g_lin) & (g_lin > 0)]
    if len(g_valid) == 0:
        continue
    sard_g0 = float(np.median(g_valid))

    s0_bias = 10 * math.log10(sard_s0 / asf_s0) if asf_s0 > 0 else float("nan")
    g0_bias = 10 * math.log10(sard_g0 / asf_g0) if asf_g0 > 0 else float("nan")

    print(
        f"{lon:6.2f}  {math.degrees(inc):7.2f}  "
        f"{sard_s0:9.5f}  {asf_s0:9.5f}  {s0_bias:+10.4f}  "
        f"{sard_g0:9.5f}  {asf_g0:9.5f}  {g0_bias:+10.4f}"
    )
