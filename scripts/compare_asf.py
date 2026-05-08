"""
Compare sardine VV dB output against ASF RTC10 GAMMA reference.

ASF product: S1B_IW_20190123T053348_DVP_RTC10_G_gpufem_6EC9
  - CRS: EPSG:32632 (UTM 32N), 10 m resolution
  - Radiometry: linear power (sigma0/gamma0), nodata=0.0

Sardine output: sardine_s1b_vv_db.tiff
  - CRS: EPSG:4326, 0.0001 deg resolution
  - Radiometry: dB (10 * log10(gamma0)), NaN = nodata

Method:
  1. Reproject ASF linear → sardine grid (bilinear)
  2. Convert ASF linear → dB
  3. Compute bias, RMSE, percentiles on joint valid mask
"""

import sys
import os
import numpy as np
import rasterio
from rasterio.warp import reproject, Resampling

SARDINE = os.environ.get(
    "SARDINE",
    "/home/datacube/dev/SARdine/sardine_s1b_tc_db.tiff",
)
ASF_VV  = (
    "/home/datacube/dev/SARdine/data/ASF/"
    "S1B_IW_20190123T053348_DVP_RTC10_G_gpufem_6EC9/"
    "S1B_IW_20190123T053348_DVP_RTC10_G_gpufem_6EC9_VV.tif"
)

print("Reading sardine output ...", flush=True)
with rasterio.open(SARDINE) as ds:
    sardine_data      = ds.read(1)
    sardine_transform = ds.transform
    sardine_crs       = ds.crs
    sardine_shape     = (ds.height, ds.width)
    sardine_bounds    = ds.bounds

sardine_valid = np.isfinite(sardine_data)
print(f"Sardine: {sardine_data.shape[1]}x{sardine_data.shape[0]}, "
      f"valid={sardine_valid.sum():,} ({100*sardine_valid.mean():.1f}%)")
print(f"  dB range:  {np.nanmin(sardine_data):.2f} .. {np.nanmax(sardine_data):.2f}")
print(f"  median dB: {np.nanmedian(sardine_data):.3f}")

print("\nReading ASF RTC10 and reprojecting ...", flush=True)
with rasterio.open(ASF_VV) as ds:
    asf_raw       = ds.read(1).astype(np.float32)
    asf_crs       = ds.crs
    asf_transform = ds.transform
    asf_nodata    = ds.nodata   # 0.0

# Mask nodata before reprojection to avoid bilinear interpolation artefacts
asf_raw[asf_raw == asf_nodata] = np.nan

asf_reproj = np.full(sardine_shape, np.nan, dtype=np.float32)
reproject(
    source=asf_raw,
    destination=asf_reproj,
    src_transform=asf_transform,
    src_crs=asf_crs,
    dst_transform=sardine_transform,
    dst_crs=sardine_crs,
    resampling=Resampling.bilinear,
    src_nodata=np.nan,
    dst_nodata=np.nan,
)

# Convert linear power → dB
with np.errstate(divide='ignore', invalid='ignore'):
    asf_db = np.where(asf_reproj > 0, 10.0 * np.log10(asf_reproj), np.nan)

asf_valid = np.isfinite(asf_db)
print(f"ASF (reprojected): valid={asf_valid.sum():,} ({100*asf_valid.mean():.1f}%)")
print(f"  dB range:  {np.nanmin(asf_db):.2f} .. {np.nanmax(asf_db):.2f}")
print(f"  median dB: {np.nanmedian(asf_db):.3f}")

# Joint valid mask
joint = sardine_valid & asf_valid
n_joint = joint.sum()
print(f"\nJoint valid pixels: {n_joint:,}")
if n_joint == 0:
    print("ERROR: no overlap — check extents/CRS", file=sys.stderr)
    sys.exit(1)

s = sardine_data[joint]
a = asf_db[joint]
diff = s - a

print("\n══════════════════════════════════════════════════════")
print("  Radiometric comparison: sardine − ASF RTC10 GAMMA   ")
print("══════════════════════════════════════════════════════")
print(f"  Mean bias:    {diff.mean():+.4f} dB")
print(f"  Median bias:  {np.median(diff):+.4f} dB")
print(f"  Std dev:      {diff.std():.4f} dB")
print(f"  RMSE:         {np.sqrt((diff**2).mean()):.4f} dB")

print("\n── |error| percentiles ──────────────────────────────")
for p in [50, 68, 90, 95, 99]:
    print(f"  p{p:2d}: |error| <= {np.percentile(np.abs(diff), p):.3f} dB")

print("\n── Bias histogram (sardine − ASF) ───────────────────")
bins  = np.arange(-5, 5.1, 0.5)
counts, edges = np.histogram(diff, bins=bins)
total = counts.sum()
peak  = counts.max()
for i, c in enumerate(counts):
    bar = '█' * int(40 * c / peak)
    lo, hi = edges[i], edges[i+1]
    print(f"  [{lo:+.1f},{hi:+.1f})  {100*c/total:5.1f}%  {bar}")

# Spatial bias map (block-average to ~1 deg tiles for quick read)
print("\n── Spatial bias by latitude band (2° bins) ──────────")
rows, cols = np.where(joint)
lats_all = sardine_transform.f + rows * sardine_transform.e  # north-up
bin_edges = np.arange(47, 51, 2)
for i in range(len(bin_edges)-1):
    mask = (lats_all >= bin_edges[i]) & (lats_all < bin_edges[i+1])
    if mask.sum() > 0:
        b = diff[mask]
        print(f"  lat [{bin_edges[i]:.0f},{bin_edges[i+1]:.0f}): "
              f"n={mask.sum():,}  bias={b.mean():+.3f}  std={b.std():.3f}")

# ── Flat-terrain bias isolation ───────────────────────────────────────────────
# Mosaic the available SRTM1 tiles, compute terrain slope, reproject to the
# sardine grid, then repeat the bias comparison restricted to flat pixels
# (slope < SLOPE_THRESHOLD_DEG).  The result tells us:
#   flat bias ≈ 0 dB  → residual is DEM-source artefact (SRTM vs GLO-30)
#   flat bias ≈ full  → algorithmic error independent of terrain slope

import glob
import math
import os
from rasterio.merge import merge as rasterio_merge

SLOPE_THRESHOLD_DEG = 2.0
SRTM_DIR = "/home/datacube/dev/SARdine/data/dem/srtm1"

print(f"\n── Flat-terrain bias isolation (slope < {SLOPE_THRESHOLD_DEG}°) ────────")

hgt_paths = sorted(glob.glob(os.path.join(SRTM_DIR, "*.hgt")))
if not hgt_paths:
    print("  ERROR: no .hgt files found in", SRTM_DIR)
else:
    # Select tiles that overlap the sardine bounding box
    overlap_paths = []
    for p in hgt_paths:
        with rasterio.open(p) as ds:
            b = ds.bounds
        if (b.left < sardine_bounds.right and b.right > sardine_bounds.left and
                b.bottom < sardine_bounds.top and b.top > sardine_bounds.bottom):
            overlap_paths.append(p)

    print(f"  SRTM tiles overlapping scene: {len(overlap_paths)}")
    if len(overlap_paths) == 0:
        print("  WARNING: no overlapping tiles – skipping flat-terrain test")
    else:
        # Mosaic tiles
        datasets = [rasterio.open(p) for p in overlap_paths]
        dem_mosaic, dem_transform_m = rasterio_merge(datasets, nodata=-32768)
        for ds in datasets:
            ds.close()

        dem_elev = dem_mosaic[0].astype(np.float32)
        dem_elev[dem_elev == -32768] = np.nan

        # Compute slope in degrees using numpy gradient.
        # SRTM is in EPSG:4326; pixel spacing converted to metres at mosaic centre.
        n_rows_dem, n_cols_dem = dem_elev.shape
        lat_centre = dem_transform_m.f + (n_rows_dem / 2.0) * dem_transform_m.e
        dy_m = abs(dem_transform_m.e) * 110540.0          # deg → m, latitude
        dx_m = abs(dem_transform_m.a) * 111320.0 * math.cos(math.radians(lat_centre))

        # gradient returns (d/dy_row, d/dx_col) with units m/m when spacings given
        dz_dy, dz_dx = np.gradient(dem_elev, dy_m, dx_m)
        slope_deg = np.degrees(np.arctan(np.sqrt(dz_dx**2 + dz_dy**2)))
        # NaN elevations propagate to NaN slope
        slope_deg[~np.isfinite(dem_elev)] = np.nan

        # Reproject slope to sardine grid (EPSG:4326, same CRS – just resample)
        dem_crs = rasterio.crs.CRS.from_epsg(4326)
        slope_reproj = np.full(sardine_shape, np.nan, dtype=np.float32)
        reproject(
            source=slope_deg.astype(np.float32),
            destination=slope_reproj,
            src_transform=dem_transform_m,
            src_crs=dem_crs,
            dst_transform=sardine_transform,
            dst_crs=sardine_crs,
            resampling=Resampling.bilinear,
            src_nodata=np.nan,
            dst_nodata=np.nan,
        )

        slope_known  = np.isfinite(slope_reproj)
        flat_mask    = slope_known & (slope_reproj <  SLOPE_THRESHOLD_DEG)
        steep_mask   = slope_known & (slope_reproj >= SLOPE_THRESHOLD_DEG)
        flat_joint   = flat_mask  & joint
        steep_joint  = steep_mask & joint
        n_flat       = int(flat_joint.sum())
        n_steep      = int(steep_joint.sum())

        print(f"  Pixels with slope known (in output grid): {slope_known.sum():>12,}")
        print(f"  Flat  (slope < {SLOPE_THRESHOLD_DEG}°, in joint mask):  {n_flat:>12,}")
        print(f"  Steep (slope >= {SLOPE_THRESHOLD_DEG}°, in joint mask): {n_steep:>12,}")

        if n_flat > 0:
            diff_flat = sardine_data[flat_joint] - asf_db[flat_joint]
            print(f"\n  Flat-terrain:  mean bias = {diff_flat.mean():+.4f} dB  "
                  f"std = {diff_flat.std():.4f} dB  n={n_flat:,}")
        if n_steep > 0:
            diff_steep = sardine_data[steep_joint] - asf_db[steep_joint]
            print(f"  Steep-terrain: mean bias = {diff_steep.mean():+.4f} dB  "
                  f"std = {diff_steep.std():.4f} dB  n={n_steep:,}")

        print()
        if n_flat > 0:
            flat_bias = float(diff_flat.mean())
            all_bias  = float(diff.mean())
            # NOTE: single-look dB comparisons have ~5.4 dB std from speckle statistics
            # regardless of terrain type; flat vs steep std will always be similar.
            # The VERDICT here is therefore unreliable without multi-looking.
            # Use multilook_compare.py for a calibration bias free of speckle artefacts.
            if abs(flat_bias) < 0.5:
                verdict = ("flat-terrain bias ≈ 0 dB → likely OK; run multilook_compare.py "
                           "for speckle-free calibration check")
            elif abs(flat_bias - all_bias) < 0.3:
                verdict = (f"flat-terrain bias ({flat_bias:+.3f} dB) ≈ all-pixel bias "
                           f"({all_bias:+.3f} dB) — dominated by single-look speckle; "
                           f"run multilook_compare.py for true calibration check")
            else:
                verdict = (f"flat={flat_bias:+.3f} dB  all={all_bias:+.3f} dB — "
                           "terrain-dependent; check DEM source or flattening")
            print(f"  VERDICT: {verdict}")
        else:
            print("  VERDICT: insufficient flat-terrain pixels in overlap area")

print("\nDone.")
