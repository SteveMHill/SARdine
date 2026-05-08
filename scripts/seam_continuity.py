"""
IW1/IW2 and IW2/IW3 seam continuity diagnostic.

Strategy
--------
1. Read sardine γ⁰ output (dB, EPSG:4326).
2. Convert to linear power.
3. Locate the two seam columns (IW1/IW2, IW2/IW3) by scanning for narrow
   longitude bands with abruptly lower valid-pixel density — seam columns can
   have gaps where neither swath reached.  If no NaN gap exists, scan for the
   steepest step in a smoothed column-mean profile.
4. For each seam: compare mean linear backscatter 20 columns to the left vs
   20 columns to the right, using only rows where both sides are valid.
5. Report the step in dB and print a local longitude profile to show the
   transition region.

The sardine output uses 0.0001° pixel spacing; 20 columns ≈ 220 m — wide
enough to average out speckle but narrow enough to stay within one sub-swath.
"""

import sys
import numpy as np
import rasterio

SARDINE_FLATTEN = "sardine_s1b_vv_30threads.tiff"
WINDOW_COLS = 30           # columns each side of seam for step estimate
SMOOTH_COLS = 100          # boxcar width for seam finder (suppress speckle)
PROFILE_HALF_WIDTH = 50    # columns printed around each seam in lon profile

print("Reading sardine γ⁰ (dB) ...", flush=True)
with rasterio.open(SARDINE_FLATTEN) as ds:
    db_data   = ds.read(1)                # float32, NaN = nodata
    tf        = ds.transform              # EPSG:4326, 0.0001° pixels
    n_rows, n_cols = db_data.shape

lon_origin = tf.c                         # west edge of pixel col 0
lon_step   = tf.a                         # positive eastward (0.0001°)
lat_origin = tf.f                         # north edge of pixel row 0
lat_step   = tf.e                         # negative southward

print(f"  Image size: {n_cols} cols × {n_rows} rows")
print(f"  Lon range: {lon_origin:.4f}° .. {lon_origin + n_cols*lon_step:.4f}°")
print(f"  Lat range: {lat_origin + n_rows*lat_step:.4f}° .. {lat_origin:.4f}°")

# Convert dB → linear; keep NaN for nodata
lin = np.where(np.isfinite(db_data), 10.0 ** (db_data / 10.0), np.nan)

# ── Step 1: column-mean profile (linear) ──────────────────────────────────────
# Using nanmean so nodata columns get NaN column mean.
print("Computing column-mean profile ...", flush=True)
col_mean_lin = np.nanmean(lin, axis=0)          # shape (n_cols,)
valid_count  = np.sum(np.isfinite(lin), axis=0) # valid pixels per column

# ── Step 2: locate seam columns ───────────────────────────────────────────────
# Smooth column mean with a wide boxcar, then find the two steepest slope
# transitions in the resulting profile.

from scipy.ndimage import uniform_filter1d

# Find column range with any valid data (trim all-NaN edge columns)
has_valid = valid_count > 0
valid_col_indices = np.where(has_valid)[0]
if len(valid_col_indices) < 100:
    print("ERROR: fewer than 100 valid columns — is the right tiff loaded?", file=sys.stderr)
    sys.exit(1)
first_valid_col = int(valid_col_indices[0])
last_valid_col  = int(valid_col_indices[-1])
print(f"  Valid column range: [{first_valid_col}, {last_valid_col}]  ({last_valid_col - first_valid_col + 1} cols)")

# Fill all-NaN columns with linear interpolation so the smoothing doesn't
# propagate NaN into the interior valid region.
col_mean_filled = col_mean_lin.copy()
if not np.isfinite(col_mean_filled).all():
    xp = np.where(np.isfinite(col_mean_filled))[0]
    fp = col_mean_filled[xp]
    col_mean_filled = np.interp(np.arange(n_cols), xp, fp)

smooth = uniform_filter1d(col_mean_filled, size=SMOOTH_COLS, mode='nearest')

# First derivative of smoothed profile
diff1 = np.gradient(smooth)

# Restrict search to the valid column range with a small inset margin.
margin = (last_valid_col - first_valid_col) // 8
search = np.abs(diff1).copy()
search[:first_valid_col + margin] = 0
search[last_valid_col  - margin:] = 0

# Divide the valid span into thirds; look for the IW1/IW2 seam in the left
# two-thirds and the IW2/IW3 seam in the right two-thirds.
span = last_valid_col - first_valid_col
left_boundary  = first_valid_col + (span * 2) // 3
right_boundary = first_valid_col + span // 3

search_left  = search.copy(); search_left[left_boundary:] = 0
search_right = search.copy(); search_right[:right_boundary] = 0

seam1_col = int(np.argmax(search_left))
seam2_col = int(np.argmax(search_right))

seam1_lon = lon_origin + seam1_col * lon_step
seam2_lon = lon_origin + seam2_col * lon_step

print(f"\nDetected seam columns (derivative-based, may be inaccurate):")
print(f"  raw seam1: col={seam1_col}  lon={seam1_lon:.4f}°E")
print(f"  raw seam2: col={seam2_col}  lon={seam2_lon:.4f}°E")

# ── Robust seam finding via ASF incidence-angle map ───────────────────────────
# The ASF inc_map.tif is in EPSG:32632 (UTM 32N).  Within each IW subswath,
# incidence angle increases smoothly with range.  At the IW1→IW2 and IW2→IW3
# boundaries the incidence angle jumps *back* (large negative eastward step).
# These two negative steps are the cleanest locator for the subswath seams.

import os, rasterio.crs

ASF_INC = (
    "data/ASF/S1B_IW_20190123T053348_DVP_RTC10_G_gpufem_6EC9/"
    "S1B_IW_20190123T053348_DVP_RTC10_G_gpufem_6EC9_inc_map.tif"
)

if os.path.exists(ASF_INC):
    print("\nLocating seam longitudes from ASF inc_map.tif ...", flush=True)
    with rasterio.open(ASF_INC) as ds:
        inc_raw    = ds.read(1).astype(np.float32)   # degrees, nodata≈0
        inc_tf     = ds.transform
        inc_crs    = ds.crs
        inc_nodata = ds.nodata or 0.0

    inc_raw[inc_raw == inc_nodata] = np.nan
    inc_raw[inc_raw <= 0]          = np.nan

    # Column mean of incidence angle
    inc_col_mean = np.nanmean(inc_raw, axis=0)        # (n_inc_cols,)
    inc_valid    = np.isfinite(inc_col_mean)
    inc_vcols    = np.where(inc_valid)[0]
    inc_first    = int(inc_vcols[0])
    inc_last     = int(inc_vcols[-1])

    # Fill NaN then smooth (IW transitions are sharp; 100-col smooth still
    # leaves clear negative spikes).
    inc_filled = inc_col_mean.copy()
    xp = np.where(np.isfinite(inc_filled))[0]
    inc_filled = np.interp(np.arange(len(inc_filled)), xp, inc_filled[xp])
    inc_smooth = uniform_filter1d(inc_filled, size=30, mode='nearest')
    inc_diff   = np.gradient(inc_smooth)              # d(inc)/d(col)

    # Restrict to interior valid columns
    margin_inc = (inc_last - inc_first) // 10
    inc_diff[:inc_first + margin_inc] = 0
    inc_diff[inc_last  - margin_inc:] = 0

    # Find two most-negative steps (IW far→near transitions)
    inc_span = inc_last - inc_first
    inc_diff_left  = inc_diff.copy(); inc_diff_left[inc_first + inc_span * 2 // 3:] = 0
    inc_diff_right = inc_diff.copy(); inc_diff_right[:inc_first + inc_span // 3] = 0

    inc_seam1_col = int(np.argmin(inc_diff_left))
    inc_seam2_col = int(np.argmin(inc_diff_right))

    # Convert UTM column → longitude (approximate, use col-centre easting)
    def utm_col_to_lon(col, tf, utm_crs):
        easting  = tf.c + (col + 0.5) * tf.a
        northing = tf.f + (tf.e * (inc_raw.shape[0] // 2))
        from pyproj import Transformer
        t = Transformer.from_crs(utm_crs, "EPSG:4326", always_xy=True)
        lon, _ = t.transform(easting, northing)
        return lon

    seam1_lon_geo = utm_col_to_lon(inc_seam1_col, inc_tf, inc_crs)
    seam2_lon_geo = utm_col_to_lon(inc_seam2_col, inc_tf, inc_crs)

    # Map back to sardine image columns
    seam1_col = int(round((seam1_lon_geo - lon_origin) / lon_step))
    seam2_col = int(round((seam2_lon_geo - lon_origin) / lon_step))
    seam1_lon = lon_origin + seam1_col * lon_step
    seam2_lon = lon_origin + seam2_col * lon_step

    print(f"  IW1/IW2 seam: inc_col={inc_seam1_col}  lon={seam1_lon:.4f}°E  sardine_col={seam1_col}")
    print(f"  IW2/IW3 seam: inc_col={inc_seam2_col}  lon={seam2_lon:.4f}°E  sardine_col={seam2_col}")
else:
    print(f"\nWARN: ASF inc_map not found at {ASF_INC} — using derivative-based seam estimate")

print(f"\nFinal seam columns:")
print(f"  IW1/IW2 seam: col={seam1_col}  lon={seam1_lon:.4f}°E")
print(f"  IW2/IW3 seam: col={seam2_col}  lon={seam2_lon:.4f}°E")

# ── Step 3: per-seam step measurement ─────────────────────────────────────────

def measure_seam_step(lin_img, seam_col, window_cols, label):
    """Measure dB step across seam_col using window_cols on each side."""
    n_rows, n_cols = lin_img.shape
    left_start  = max(0, seam_col - window_cols)
    left_end    = seam_col
    right_start = seam_col
    right_end   = min(n_cols, seam_col + window_cols)

    left_strip  = lin_img[:, left_start:left_end]
    right_strip = lin_img[:, right_start:right_end]

    # Only use rows where both sides have at least one valid pixel
    left_row_mean  = np.nanmean(left_strip,  axis=1)  # (n_rows,)
    right_row_mean = np.nanmean(right_strip, axis=1)

    both_valid = np.isfinite(left_row_mean) & np.isfinite(right_row_mean)
    n_pairs = both_valid.sum()

    if n_pairs < 100:
        print(f"\n  {label}: insufficient row pairs ({n_pairs}) — seam may be in a gap")
        return

    L = left_row_mean[both_valid]
    R = right_row_mean[both_valid]

    row_diff_db = 10.0 * np.log10(R / L)
    step_mean_db   = 10.0 * np.log10(np.mean(R) / np.mean(L))
    step_median_db = float(np.median(row_diff_db))
    std_db         = float(np.std(row_diff_db))

    print(f"\n{'─'*60}")
    print(f"  {label}  (col {seam_col}, lon {lon_origin + seam_col*lon_step:.4f}°E)")
    print(f"{'─'*60}")
    print(f"  Valid row pairs: {n_pairs:,} / {n_rows:,}")
    print(f"  Left  strip cols [{left_start},{left_end}): mean = {10*np.log10(np.mean(L)):+.3f} dB")
    print(f"  Right strip cols [{right_start},{right_end}): mean = {10*np.log10(np.mean(R)):+.3f} dB")
    print(f"  Step mean   (linear means):   {step_mean_db:+.3f} dB  [sensitive to bright targets]")
    print(f"  Step median (row-wise): {step_median_db:+.3f} dB  [robust metric for speckled imagery]")
    print(f"  Step std (row-wise):    {std_db:.3f} dB")
    # Use median as the authoritative metric — mean is skewed by isolated bright targets
    if abs(step_median_db) < 0.1:
        verdict = "PASS — median step < 0.1 dB; seam is radiometrically continuous"
    elif abs(step_median_db) < 0.5:
        verdict = "WARN — median step 0.1–0.5 dB; minor seam may be visible at high contrast"
    else:
        verdict = "FAIL — median step > 0.5 dB; seam will be visible in output"
    print(f"  Verdict: {verdict}")

    return step_db


print()
step1 = measure_seam_step(lin, seam1_col, WINDOW_COLS, "IW1/IW2 seam")
step2 = measure_seam_step(lin, seam2_col, WINDOW_COLS, "IW2/IW3 seam")

# ── Step 4: print lon profile around each seam ────────────────────────────────

def print_lon_profile(col_mean_lin, seam_col, half_width, lon_origin, lon_step, label):
    c0 = max(0, seam_col - half_width)
    c1 = min(len(col_mean_lin), seam_col + half_width + 1)
    vals_db = 10.0 * np.log10(np.maximum(col_mean_lin[c0:c1], 1e-20))
    ref_db  = float(np.nanmedian(vals_db))
    peak    = float(np.nanmax(np.abs(vals_db - ref_db)))
    if peak == 0:
        peak = 1.0

    print(f"\n── {label} profile  (col_mean γ⁰ dB, centred at col {seam_col}) ──")
    for i, col in enumerate(range(c0, c1)):
        v = vals_db[i]
        bar_len = int(20 * (v - ref_db + peak) / (2 * peak))
        bar_len = max(0, min(40, bar_len))
        marker = "◄seam" if col == seam_col else ""
        lon = lon_origin + col * lon_step
        print(f"  col={col:5d}  lon={lon:.4f}°  {v:+.2f} dB  {'█'*bar_len} {marker}")

print_lon_profile(col_mean_filled, seam1_col, PROFILE_HALF_WIDTH, lon_origin, lon_step, "IW1/IW2")
print_lon_profile(col_mean_filled, seam2_col, PROFILE_HALF_WIDTH, lon_origin, lon_step, "IW2/IW3")

print("\nDone.")
