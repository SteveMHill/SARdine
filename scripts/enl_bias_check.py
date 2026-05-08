"""
Compare sardine vs ASF using MEAN LINEAR (not dB) to remove the ENL/logarithm bias.

For SAR speckle: E[log(x)] < log(E[x]) by a factor that depends on ENL.
If sardine ENL << ASF ENL, the mean of sardine_dB will be biased low compared
to mean of asf_dB even if the TRUE mean backscatter is identical.

This script computes 10*log10(mean_linear) for both products — which gives
the TRUE mean difference without ENL artefacts.
"""

import numpy as np
import rasterio
from rasterio.warp import reproject, Resampling

SARDINE = "sardine_s1b_vv_30threads.tiff"
ASF_VV  = (
    "data/ASF/S1B_IW_20190123T053348_DVP_RTC10_G_gpufem_6EC9/"
    "S1B_IW_20190123T053348_DVP_RTC10_G_gpufem_6EC9_VV.tif"
)

print("Reading sardine (dB) and converting to linear ...", flush=True)
with rasterio.open(SARDINE) as ds:
    sardine_db    = ds.read(1)
    sardine_tf    = ds.transform
    sardine_crs   = ds.crs
    sardine_shape = (ds.height, ds.width)

sardine_linear = np.where(np.isfinite(sardine_db), 10 ** (sardine_db / 10), np.nan)

print("Reading and reprojecting ASF (already linear) ...", flush=True)
with rasterio.open(ASF_VV) as ds:
    asf_raw   = ds.read(1).astype(np.float32)
    asf_raw[asf_raw <= 0] = np.nan
    asf_crs   = ds.crs
    asf_tf    = ds.transform

asf_reproj = np.full(sardine_shape, np.nan, dtype=np.float32)
reproject(
    source=asf_raw,
    destination=asf_reproj,
    src_transform=asf_tf, src_crs=asf_crs,
    dst_transform=sardine_tf, dst_crs=sardine_crs,
    resampling=Resampling.bilinear,
    src_nodata=np.nan, dst_nodata=np.nan,
)

# Joint valid mask
joint = np.isfinite(sardine_linear) & np.isfinite(asf_reproj)
n = joint.sum()
print(f"Joint valid pixels: {n:,}")

s_lin = sardine_linear[joint]
a_lin = asf_reproj[joint]

# Method 1: mean of pixel-wise dB difference (includes ENL bias)
s_db = sardine_db[np.isfinite(sardine_db) & joint]
a_db = 10 * np.log10(a_lin)
mean_diff_db = float(np.mean(s_db - a_db))

# Method 2: dB of mean linear (TRUE mean calibration difference, no ENL bias)
mean_sardine_db = 10 * np.log10(np.mean(s_lin))
mean_asf_db     = 10 * np.log10(np.mean(a_lin))
true_bias_db    = mean_sardine_db - mean_asf_db

# Method 3: median linear
med_sardine     = float(np.median(s_lin))
med_asf         = float(np.median(a_lin))
median_bias_db  = 10 * np.log10(med_sardine / med_asf)

# ENL estimate from sardine variance
# For L-look SAR: Var[x] = E[x]^2 / L  → L = E[x]^2 / Var[x]
enl_sardine = float(np.mean(s_lin)**2 / np.var(s_lin))
enl_asf     = float(np.mean(a_lin)**2 / np.var(a_lin))

import math

print()
print("══════════════════════════════════════════════════════════════")
print("  ENL analysis and logarithm bias correction                  ")
print("══════════════════════════════════════════════════════════════")
print(f"  Sardine ENL estimate:     {enl_sardine:.2f}")
print(f"  ASF ENL estimate:         {enl_asf:.2f}")

# digamma-based log bias: E[dB(x)] = dB(E[x]) + 10/ln(10) * (psi(L) - log(L))
# digamma approximation: psi(L) ≈ log(L) - 1/(2L) for large L
def psi_approx(L):
    """Digamma approximation (harmonic series, good for L>=1)."""
    s = 0.0
    for k in range(1, int(L) + 1):
        s += 1.0 / k
    return s - 0.5772156649  # subtract Euler–Mascheroni

log_bias_sardine = (psi_approx(enl_sardine) - math.log(enl_sardine)) * 10 / math.log(10)
log_bias_asf     = (psi_approx(enl_asf) - math.log(enl_asf)) * 10 / math.log(10)

print()
print(f"  Log-bias (Jensen) sardine: {log_bias_sardine:+.3f} dB (E[dB] - dB(E[linear]))")
print(f"  Log-bias (Jensen) ASF:     {log_bias_asf:+.3f} dB")
print(f"  Expected bias from ENL:    {log_bias_sardine - log_bias_asf:+.3f} dB")
print()
print("  Method 1: mean(sardine_dB - asf_dB) [includes ENL bias]")
print(f"    = {mean_diff_db:+.4f} dB")
print()
print("  Method 2: 10*log10(mean_sardine_linear / mean_asf_linear) [TRUE bias]")
print(f"    sardine mean = {mean_sardine_db:+.4f} dBγ⁰")
print(f"    ASF mean     = {mean_asf_db:+.4f} dBγ⁰")
print(f"    TRUE bias    = {true_bias_db:+.4f} dB")
print()
print("  Method 3: 10*log10(median_sardine_linear / median_asf_linear)")
print(f"    = {median_bias_db:+.4f} dB")
print()
print("══════════════════════════════════════════════════════════════")
print("  Verdict:")
expected_method1 = true_bias_db + (log_bias_sardine - log_bias_asf)
print(f"    True calibration bias = {true_bias_db:+.3f} dB")
print(f"    ENL-driven dB bias    = {log_bias_sardine - log_bias_asf:+.3f} dB")
print(f"    Predicted method1     = {expected_method1:+.3f} dB (check vs {mean_diff_db:+.3f})")
