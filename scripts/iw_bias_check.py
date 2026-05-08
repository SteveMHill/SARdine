"""Compare sardine gamma0 vs ASF gamma0 at three IW incidence zones."""
import math
import numpy as np
import rasterio
from rasterio.warp import transform as rio_transform

ASF_VV = (
    "data/ASF/S1B_IW_20190123T053348_DVP_RTC10_G_gpufem_6EC9/"
    "S1B_IW_20190123T053348_DVP_RTC10_G_gpufem_6EC9_VV.tif"
)
ASF_INC = (
    "data/ASF/S1B_IW_20190123T053348_DVP_RTC10_G_gpufem_6EC9/"
    "S1B_IW_20190123T053348_DVP_RTC10_G_gpufem_6EC9_inc_map.tif"
)
FLATTEN = "sardine_s1b_vv_30threads.tiff"

HALF = 50

for lat, lon, label in [
    (49.3, 9.0,  "IW1 ~36°"),
    (49.3, 10.1, "IW2 ~40°"),
    (49.3, 11.0, "IW3 ~43°"),
]:
    xs, ys = rio_transform("EPSG:4326", "EPSG:32632", [lon], [lat])
    x_utm, y_utm = xs[0], ys[0]

    with rasterio.open(ASF_INC) as ds:
        r, c = ds.index(x_utm, y_utm)
        inc_patch = ds.read(
            1, window=rasterio.windows.Window(c - 5, r - 5, 10, 10)
        ).astype(np.float32)
    inc_rad = float(np.nanmedian(inc_patch[inc_patch > 0]))

    with rasterio.open(ASF_VV) as ds:
        r, c = ds.index(x_utm, y_utm)
        a_patch = ds.read(
            1, window=rasterio.windows.Window(c - HALF, r - HALF, HALF * 2, HALF * 2)
        ).astype(np.float32)
    a_patch = a_patch[a_patch > 0]
    a_med = float(np.median(a_patch))

    with rasterio.open(FLATTEN) as ds:
        r, c = ds.index(lon, lat)
        s_patch = ds.read(
            1, window=rasterio.windows.Window(c - HALF, r - HALF, HALF * 2, HALF * 2)
        ).astype(np.float32)
    s_lin = 10 ** (s_patch / 10)
    s_lin = s_lin[np.isfinite(s_lin) & (s_lin > 0)]
    s_med = float(np.median(s_lin))

    ratio_db = 10 * math.log10(s_med / a_med)
    print(
        f"{label}: inc={math.degrees(inc_rad):.1f}°  "
        f"sardine_γ⁰={s_med:.5f}  asf_γ⁰={a_med:.5f}  "
        f"ratio={ratio_db:+.3f} dB"
    )
