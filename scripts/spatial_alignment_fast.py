#!/usr/bin/env python3
"""
Fast spatial alignment check: windowed reading, no full-image reproject.

Reads only a narrow latitude strip from each product, finds dark-water
features (River Rhine, Bodensee) in the 1-km-smoothed profiles, and
reports the longitude offset between sardine and ASF.

Run:
  python3 scripts/spatial_alignment_fast.py
or:
  SARDINE=/path/to/output.tiff python3 scripts/spatial_alignment_fast.py
"""
import os, sys
import numpy as np
import rasterio
from rasterio.windows import from_bounds
from rasterio.warp import transform_bounds, reproject, Resampling
from scipy.ndimage import uniform_filter1d

_REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SARDINE = os.environ.get("SARDINE", os.path.join(_REPO, "sardine_s1b_tc_db.tiff"))
ASF_VV  = os.path.join(
    _REPO,
    "data/ASF/S1B_IW_20190123T053348_DVP_RTC10_G_gpufem_6EC9/"
    "S1B_IW_20190123T053348_DVP_RTC10_G_gpufem_6EC9_VV.tif"
)

# (lat_deg, lon_search_lo, lon_search_hi, expected_water_lon, label)
PROFILES = [
    # Rhine near Rastatt: water at ~8.09°E at 48.85°N
    (48.85, 7.9, 8.35, 8.09, "Rhine @ Rastatt/KA (48.85°N)"),
    # Rhine near Basel: water at ~7.57°E at 47.82°N
    (47.82, 7.35, 7.80, 7.57, "Rhine @ Basel (47.82°N)"),
    # Main river east of Frankfurt: water at ~8.76°E at 50.10°N
    # (outside sardine coverage – skip if sardine stops at 49.74°N)
]

STRIP_HEIGHT = 11   # rows to average for the profile (±5 rows)
SMOOTH_PX    = 90   # 1-km smoothing (90 × 11.1m ≈ 1 km)
PASS_M       = 200  # pass threshold (metres)


def sardine_profile(path, lat_deg, lon_lo, lon_hi):
    """Return (lons, vals_db) for a latitude strip in the sardine TIFF."""
    with rasterio.open(path) as ds:
        tf  = ds.transform
        row_centre = int((lat_deg - tf.f) / tf.e)
        if row_centre < 0 or row_centre >= ds.height:
            return None, None
        r0 = max(0, row_centre - STRIP_HEIGHT // 2)
        r1 = min(ds.height, row_centre + STRIP_HEIGHT // 2 + 1)
        col_lo = max(0, int((lon_lo - tf.c) / tf.a))
        col_hi = min(ds.width, int((lon_hi - tf.c) / tf.a) + 1)
        if col_lo >= col_hi:
            return None, None
        win = rasterio.windows.Window(col_lo, r0, col_hi - col_lo, r1 - r0)
        strip = ds.read(1, window=win).astype(np.float32)
        strip[strip < -9000] = np.nan
        vals = np.nanmean(strip, axis=0)
        lons = tf.c + (col_lo + np.arange(strip.shape[1])) * tf.a
    return lons, vals


def asf_profile(path, lat_deg, lon_lo, lon_hi, target_lons):
    """Return vals_db sampled at target_lons by reading a small window from ASF."""
    with rasterio.open(path) as ds:
        tf_asf  = ds.transform
        crs_asf = ds.crs
        crs_sar = rasterio.CRS.from_epsg(4326)

        # Transform the sardine latitude strip bounds into ASF CRS
        pad = 0.05  # degrees padding
        bounds_4326 = (lon_lo - pad, lat_deg - 0.05, lon_hi + pad, lat_deg + 0.05)
        bounds_asf  = transform_bounds(crs_sar, crs_asf, *bounds_4326)
        win = from_bounds(*bounds_asf, transform=tf_asf)
        win = win.intersection(rasterio.windows.Window(0, 0, ds.width, ds.height))
        strip_asf = ds.read(1, window=win).astype(np.float32)
        strip_asf[strip_asf <= 0] = np.nan
        # dB
        strip_asf_db = np.where(np.isfinite(strip_asf), 10 * np.log10(strip_asf), np.nan)
        win_tf = ds.window_transform(win)

        # For each target longitude, find the ASF column; average ±5 rows
        n_rows = strip_asf.shape[0]
        row_mid = n_rows // 2
        r0 = max(0, row_mid - STRIP_HEIGHT // 2)
        r1 = min(n_rows, row_mid + STRIP_HEIGHT // 2 + 1)

        vals = np.full(len(target_lons), np.nan, dtype=np.float32)
        for i, lon in enumerate(target_lons):
            # We need to map (lat_deg, lon) → (row, col) in ASF CRS
            # Use rasterio.transform.rowcol on win_tf
            from rasterio.transform import rowcol
            from pyproj import Transformer
            t = Transformer.from_crs(crs_sar, crs_asf, always_xy=True)
            x, y = t.transform(lon, lat_deg)
            row_f, col_f = rowcol(win_tf, x, y, op=float)
            col_i = int(round(col_f))
            if col_i < 0 or col_i >= strip_asf.shape[1]:
                continue
            col_block = strip_asf_db[r0:r1, max(0, col_i-1):col_i+2]
            v = col_block[np.isfinite(col_block)]
            if v.size > 0:
                vals[i] = v.mean()
    return vals


print(f"Sardine: {SARDINE}")
print(f"ASF VV:  {ASF_VV}")
print()

results = []
for lat, lon_lo, lon_hi, expected, label in PROFILES:
    print(f"Profile: {label}")
    lons, vals_sar = sardine_profile(SARDINE, lat, lon_lo, lon_hi)
    if lons is None:
        print("  SKIP: latitude outside sardine coverage")
        print()
        continue

    # Smooth to 1 km to suppress speckle
    valid = np.isfinite(vals_sar)
    vals_sar_s = np.where(valid,
        uniform_filter1d(np.where(valid, vals_sar, 0.0), SMOOTH_PX)
        / uniform_filter1d(valid.astype(float), SMOOTH_PX).clip(0.01),
        np.nan)

    # Find minimum in search window
    mask = (lons >= lon_lo) & (lons <= lon_hi) & np.isfinite(vals_sar_s)
    if not mask.any():
        print("  WARNING: no valid sardine data in search window")
        print()
        continue
    sar_min_idx = np.where(mask, vals_sar_s, np.inf).argmin()
    sar_min_lon = lons[sar_min_idx]
    sar_min_db  = vals_sar_s[sar_min_idx]
    print(f"  Sardine min: lon = {sar_min_lon:.4f}°  val = {sar_min_db:.1f} dB")

    # Sample ASF at the same longitudes
    vals_asf = asf_profile(ASF_VV, lat, lon_lo, lon_hi, lons)
    valid_a = np.isfinite(vals_asf)
    if not valid_a.any():
        print("  WARNING: no valid ASF data in search window")
        print()
        continue
    vals_asf_s = np.where(valid_a,
        uniform_filter1d(np.where(valid_a, vals_asf, 0.0), SMOOTH_PX)
        / uniform_filter1d(valid_a.astype(float), SMOOTH_PX).clip(0.01),
        np.nan)

    mask_a = (lons >= lon_lo) & (lons <= lon_hi) & np.isfinite(vals_asf_s)
    if not mask_a.any():
        print("  WARNING: no valid ASF data after smoothing")
        print()
        continue
    asf_min_idx = np.where(mask_a, vals_asf_s, np.inf).argmin()
    asf_min_lon = lons[asf_min_idx]
    asf_min_db  = vals_asf_s[asf_min_idx]
    print(f"  ASF RTC min: lon = {asf_min_lon:.4f}°  val = {asf_min_db:.1f} dB")

    delta_deg = sar_min_lon - asf_min_lon
    delta_px  = delta_deg / 0.0001
    delta_m   = abs(delta_deg) * 111_320 * np.cos(np.radians(lat))
    PASS = delta_m < PASS_M
    verdict = "PASS" if PASS else "UNCERTAIN"
    print(f"  Offset:      Δlon = {delta_deg:+.4f}° = {delta_px:+.0f} px = {delta_m:.0f} m  [{verdict}]")
    if abs(sar_min_db) < 5 or abs(asf_min_db) < 5:
        print("  NOTE: minimum value > -5 dB — likely NOT open water; treat result as unreliable")
    results.append((label, delta_m, PASS))
    print()

if results:
    print("=== Summary ===")
    for label, delta_m, passed in results:
        print(f"  {'PASS' if passed else 'UNCERTAIN':9s}  {delta_m:6.0f} m  {label}")
    print()
    if all(p for _, _, p in results):
        print("VERDICT: PASS — all water-feature offsets within threshold")
    else:
        print("VERDICT: UNCERTAIN — manual visual verification recommended")
        print("  Load sardine_s1b_tc_db.tiff and ASF VV.tif in QGIS,")
        print("  align to EPSG:4326, overlay, and visually check river edges.")
print("Done.")
