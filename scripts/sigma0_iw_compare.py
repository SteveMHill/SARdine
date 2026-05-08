"""Compare sardine sigma0 vs ASF sigma0 at three IW zones."""
import math
import numpy as np
import rasterio
from rasterio.warp import transform as rio_transform

NOFLATTEN = "sardine_s1b_vv_noflatten.tiff"
ASF_VV = (
    "data/ASF/S1B_IW_20190123T053348_DVP_RTC10_G_gpufem_6EC9/"
    "S1B_IW_20190123T053348_DVP_RTC10_G_gpufem_6EC9_VV.tif"
)
ASF_INC = (
    "data/ASF/S1B_IW_20190123T053348_DVP_RTC10_G_gpufem_6EC9/"
    "S1B_IW_20190123T053348_DVP_RTC10_G_gpufem_6EC9_inc_map.tif"
)

for lat, lon, label in [
    (49.3, 9.0, "IW1"),
    (49.3, 10.1, "IW2"),
    (49.3, 11.0, "IW3"),
]:
    xs, ys = rio_transform("EPSG:4326", "EPSG:32632", [lon], [lat])
    x_utm, y_utm = xs[0], ys[0]

    with rasterio.open(ASF_INC) as ds:
        r, c = ds.index(x_utm, y_utm)
        inc_p = ds.read(1, window=rasterio.windows.Window(c - 5, r - 5, 10, 10)).astype(np.float32)
    inc = float(np.nanmedian(inc_p[inc_p > 0]))
    cos_inc = math.cos(inc)

    with rasterio.open(ASF_VV) as ds:
        r, c = ds.index(x_utm, y_utm)
        a = ds.read(1, window=rasterio.windows.Window(c - 50, r - 50, 100, 100)).astype(np.float32)
    a = a[a > 0]
    asf_g0 = float(np.median(a))
    asf_s0 = asf_g0 * cos_inc

    with rasterio.open(NOFLATTEN) as ds:
        r, c = ds.index(lon, lat)
        s = ds.read(1, window=rasterio.windows.Window(c - 50, r - 50, 100, 100)).astype(np.float32)
    s_lin = 10 ** (s / 10)
    s_lin = s_lin[np.isfinite(s_lin) & (s_lin > 0)]
    sardine_s0 = float(np.median(s_lin))

    ratio = sardine_s0 / asf_s0
    print(
        f"{label}: inc={math.degrees(inc):.1f}°  "
        f"sardine_σ⁰={sardine_s0:.5f}  asf_σ⁰={asf_s0:.5f}  "
        f"ratio={ratio:.4f} ({10*math.log10(ratio):+.3f} dB)"
    )
