"""
Three-way comparison: sardine (flatten) vs sardine (no-flatten) vs ASF γ⁰.

Scenario A: flattening doing nothing  → flatten_bias ≈ noflatten_bias
Scenario B: flattening works, σ⁰ low  → noflatten_bias ≈ 2 × flatten_bias
"""
import math
import numpy as np
import rasterio
from rasterio.warp import reproject, Resampling

FLATTEN   = "/home/datacube/dev/SARdine/sardine_s1b_vv_30threads.tiff"
NOFLATTEN = "/home/datacube/dev/SARdine/sardine_s1b_vv_noflatten.tiff"
ASF_VV    = (
    "/home/datacube/dev/SARdine/data/ASF/"
    "S1B_IW_20190123T053348_DVP_RTC10_G_gpufem_6EC9/"
    "S1B_IW_20190123T053348_DVP_RTC10_G_gpufem_6EC9_VV.tif"
)

print("Reading sardine outputs ...", flush=True)
with rasterio.open(FLATTEN) as ds:
    flat_data = ds.read(1)
    t = ds.transform
    crs = ds.crs
    shape = (ds.height, ds.width)

with rasterio.open(NOFLATTEN) as ds:
    noflat_data = ds.read(1)

print("Reading + reprojecting ASF ...", flush=True)
with rasterio.open(ASF_VV) as ds:
    asf_raw = ds.read(1).astype(np.float32)
    asf_raw[asf_raw == ds.nodata] = np.nan

asf_reproj = np.full(shape, np.nan, dtype=np.float32)
with rasterio.open(ASF_VV) as ds:
    reproject(
        source=asf_raw,
        destination=asf_reproj,
        src_transform=ds.transform,
        src_crs=ds.crs,
        dst_transform=t,
        dst_crs=crs,
        resampling=Resampling.bilinear,
        src_nodata=np.nan,
        dst_nodata=np.nan,
    )

with np.errstate(divide='ignore', invalid='ignore'):
    asf_db = np.where(asf_reproj > 0, 10.0 * np.log10(asf_reproj), np.nan)

joint = np.isfinite(flat_data) & np.isfinite(noflat_data) & np.isfinite(asf_db)
print(f"Joint valid pixels: {joint.sum():,}", flush=True)

f  = flat_data[joint]
nf = noflat_data[joint]
a  = asf_db[joint]

diff_fn = f - nf          # flatten − no-flatten
diff_fa = f - a           # flatten − ASF
diff_nfa = nf - a         # no-flatten − ASF

print()
print("══════════════════════════════════════════════════════")
print("  flatten − no_flatten  (should be +1.12 dB if flatten works)")
print(f"  mean = {diff_fn.mean():+.4f} dB   std = {diff_fn.std():.4f} dB")
print()
print("  flatten − ASF         (was -1.28 dB overall, -1.13 dB flat)")
print(f"  mean = {diff_fa.mean():+.4f} dB   std = {diff_fa.std():.4f} dB")
print()
print("  no_flatten − ASF")
print(f"  mean = {diff_nfa.mean():+.4f} dB   std = {diff_nfa.std():.4f} dB")
print("══════════════════════════════════════════════════════")
print()

fn_mean  = float(diff_fn.mean())
nfa_mean = float(diff_nfa.mean())
fa_mean  = float(diff_fa.mean())

if abs(fn_mean) < 0.1:
    print("VERDICT: Scenario A — terrain flattening is doing NOTHING.")
    print("  The flatten / no-flatten outputs are identical.")
    print("  The -1.12 dB bias vs ASF is because sardine writes σ⁰ not γ⁰.")
    print("  Fix: check why sigma0 /= w is not executing or w ≈ 1.")
elif abs(nfa_mean - 2 * fa_mean) < 0.3:
    print("VERDICT: Scenario B — flattening works, but sardine σ⁰ is biased low.")
    print(f"  noflatten-ASF ≈ {nfa_mean:+.3f} dB ≈ 2× (flatten-ASF {fa_mean:+.3f} dB).")
    print("  Possible cause: calibration LUT scale error or missing absolute constant.")
else:
    print(f"VERDICT: mixed — flatten-noflatten={fn_mean:+.3f} dB, "
          f"flatten-ASF={fa_mean:+.3f} dB, noflatten-ASF={nfa_mean:+.3f} dB.")
    print("  Investigate further.")
