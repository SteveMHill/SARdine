#!/usr/bin/env python3
"""
Spatial alignment check: Rhine river profile method.

Extracts a west-east backscatter profile at a fixed latitude that crosses
the Rhine river, from both sardine and ASF.  The Rhine appears as a sharp
dark feature (~−25 dB) flanked by bright urban areas.  If sardine is
spatially aligned, the dark trough should occur at the same longitude in
both products.

Pass threshold: river-minimum longitude agrees within ±3 pixels (~0.0003 deg ≈ 25 m).
"""
import os, sys
import numpy as np
import rasterio
from rasterio.warp import reproject, Resampling

_REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SARDINE = os.environ.get("SARDINE", os.path.join(_REPO, "sardine_s1b_tc_db.tiff"))
ASF_VV  = os.path.join(
    _REPO,
    "data/ASF/S1B_IW_20190123T053348_DVP_RTC10_G_gpufem_6EC9/"
    "S1B_IW_20190123T053348_DVP_RTC10_G_gpufem_6EC9_VV.tif"
)

# Rhine near Rastatt/Karlsruhe: lat ≈ 48.85°N, Rhine at lon ≈ 8.15°E
PROFILE_LAT = 48.85

# Also check Rhine near Breisach: lat ≈ 48.03°N, Rhine at lon ≈ 7.58°E
PROFILE_LAT2 = 48.03


def profile_at_lat(sardine_path, asf_path, lat_deg, smooth_km=1.0):
    """
    Return (lons_sar, vals_sar_db, lons_asf_db, vals_asf_db) along a latitude.

    Both profiles are extracted at the sardine grid pixel spacing (0.0001 deg).
    ASF is reprojected to the sardine grid first.
    smooth_km: Gaussian smoothing to reduce speckle before alignment.
    """
    with rasterio.open(sardine_path) as ds:
        tf   = ds.transform
        crs  = ds.crs
        # Row index for this latitude
        row = int((lat_deg - tf.f) / tf.e)
        if row < 0 or row >= ds.height:
            return None, None, None, None
        # Read ±5 rows and average (5 × 11m = 55m along-track)
        r0 = max(0, row - 5)
        r1 = min(ds.height, row + 6)
        strip = ds.read(1, window=rasterio.windows.Window(0, r0, ds.width, r1 - r0))
        strip = strip.astype(np.float32)
        strip[strip < -9000] = np.nan
        vals_sar = np.nanmean(strip, axis=0)
        lons_sar = tf.c + np.arange(ds.width) * tf.a
        sar_tf   = tf
        sar_crs  = crs
        sar_w    = ds.width
        sar_h    = ds.height

    # Reproject ASF to sardine grid
    with rasterio.open(asf_path) as ds:
        asf_raw = ds.read(1).astype(np.float32)
        asf_raw[asf_raw <= 0] = np.nan
        # dB conversion
        asf_db = np.where(np.isfinite(asf_raw), 10.0 * np.log10(asf_raw), np.nan)

        asf_reproj = np.full((sar_h, sar_w), np.nan, dtype=np.float32)
        reproject(asf_db, asf_reproj,
                  src_transform=ds.transform, src_crs=ds.crs,
                  dst_transform=sar_tf, dst_crs=sar_crs,
                  resampling=Resampling.bilinear,
                  src_nodata=np.nan, dst_nodata=np.nan)

    strip_asf = asf_reproj[r0:r1, :]
    vals_asf  = np.nanmean(strip_asf, axis=0)

    # Apply 1-km smoothing to suppress speckle (100 pixels at 0.0001 deg)
    from scipy.ndimage import uniform_filter
    smooth_px = max(1, int(smooth_km * 1000 / 11.1))
    mask_sar = np.isfinite(vals_sar)
    mask_asf = np.isfinite(vals_asf)
    vals_sar_s = np.where(mask_sar, uniform_filter(np.where(mask_sar, vals_sar, 0.0), smooth_px)
                          / uniform_filter(mask_sar.astype(float), smooth_px).clip(0.01), np.nan)
    vals_asf_s = np.where(mask_asf, uniform_filter(np.where(mask_asf, vals_asf, 0.0), smooth_px)
                          / uniform_filter(mask_asf.astype(float), smooth_px).clip(0.01), np.nan)

    return lons_sar, vals_sar_s, vals_asf_s


def find_minimum_lon(lons, vals, lon_min, lon_max):
    """Return longitude of minimum value within [lon_min, lon_max]."""
    mask = (lons >= lon_min) & (lons <= lon_max) & np.isfinite(vals)
    if not mask.any():
        return None
    idx = np.where(mask, vals, np.inf).argmin()
    return lons[idx]


print(f"Sardine: {SARDINE}")
print(f"ASF VV:  {ASF_VV}")
print()

PAX_DEG = 0.0001   # pixel size

for lat, label, lon_search_lo, lon_search_hi, expected_lon in [
    (PROFILE_LAT,  "Rastatt/Karlsruhe (lat=48.85°N)", 7.9, 8.35, 8.14),
    (PROFILE_LAT2, "Breisach (lat=48.03°N)",           7.4, 7.9,  7.58),
]:
    print(f"Profile at lat = {lat}°N — {label}")
    lons, vals_sar, vals_asf = profile_at_lat(SARDINE, ASF_VV, lat)
    if lons is None:
        print("  ERROR: latitude outside sardine extent")
        continue

    lon_min_sar = find_minimum_lon(lons, vals_sar, lon_search_lo, lon_search_hi)
    lon_min_asf = find_minimum_lon(lons, vals_asf, lon_search_lo, lon_search_hi)

    if lon_min_sar is None:
        print(f"  WARNING: No sardine data in search window [{lon_search_lo}, {lon_search_hi}]")
    else:
        sar_db = vals_sar[np.argmin(np.abs(lons - lon_min_sar))]
        print(f"  Sardine  dark-water minimum: lon = {lon_min_sar:.4f}°E   val = {sar_db:.1f} dB")

    if lon_min_asf is None:
        print(f"  WARNING: No ASF data in search window [{lon_search_lo}, {lon_search_hi}]")
    else:
        asf_db = vals_asf[np.argmin(np.abs(lons - lon_min_asf))]
        print(f"  ASF RTC  dark-water minimum: lon = {lon_min_asf:.4f}°E   val = {asf_db:.1f} dB")

    if lon_min_sar is not None and lon_min_asf is not None:
        delta_deg = lon_min_sar - lon_min_asf
        delta_px  = delta_deg / PAX_DEG
        delta_m   = abs(delta_deg) * 111_320 * np.cos(np.radians(lat))
        PASS = delta_m < 40.0
        verdict = "PASS" if PASS else "FAIL"
        print(f"  Offset:  Δlon = {delta_deg:+.4f}° = {delta_px:+.1f} px = {delta_m:.0f} m   [{verdict}]")
    print()

print("Done.")
