"""Quick incidence-angle + area diagnostic for the ASF reference product."""
import math
import numpy as np
import rasterio
from rasterio.windows import Window

BASE = (
    "/home/datacube/dev/SARdine/data/ASF/"
    "S1B_IW_20190123T053348_DVP_RTC10_G_gpufem_6EC9/"
    "S1B_IW_20190123T053348_DVP_RTC10_G_gpufem_6EC9"
)

# Read a decimated overview (every 20th row/col) to stay in RAM
STRIDE = 20
print("Reading incidence angle map (decimated) ...")
with rasterio.open(BASE + "_inc_map.tif") as ds:
    rows = list(range(0, ds.height, STRIDE))
    cols = list(range(0, ds.width, STRIDE))
    inc_vals = []
    for r in rows:
        win = Window(0, r, ds.width, 1)
        row_data = ds.read(1, window=win)[0]
        inc_vals.append(row_data[::STRIDE])
    inc_arr = np.concatenate(inc_vals).astype(np.float32)

valid_inc = inc_arr[(inc_arr > 0) & np.isfinite(inc_arr)]
mn_rad = float(np.mean(valid_inc))
mn_deg = math.degrees(mn_rad)
print(f"  n_valid = {len(valid_inc):,}")
print(f"  min={math.degrees(float(np.min(valid_inc))):.2f}°  max={math.degrees(float(np.max(valid_inc))):.2f}°  mean={mn_deg:.2f}°")
print(f"  cos(mean) = {math.cos(mn_rad):.4f}")
print(f"  10*log10(cos(mean)) = {10*math.log10(math.cos(mn_rad)):.4f} dB")

print()
print("Reading scattering area map (decimated) ...")
with rasterio.open(BASE + "_area.tif") as ds:
    area_vals = []
    for r in rows:
        win = Window(0, r, ds.width, 1)
        row_data = ds.read(1, window=win)[0]
        area_vals.append(row_data[::STRIDE])
    area_arr = np.concatenate(area_vals).astype(np.float32)

valid_area = area_arr[(area_arr > 0) & np.isfinite(area_arr)]
med_area = float(np.median(valid_area))
print(f"  n_valid = {len(valid_area):,}")
print(f"  median area = {med_area:.2f} m^2  (10m×10m flat = 100 m^2)")
cos_implied = 100.0 / med_area
if 0 < cos_implied <= 1:
    print(f"  cos_implied (100/median) = {cos_implied:.4f}  θ = {math.degrees(math.acos(cos_implied)):.2f}°")
    print(f"  10*log10(cos_implied) = {10*math.log10(cos_implied):.4f} dB")
else:
    print(f"  cos_implied = {cos_implied:.4f} (out of range — area != flat)")

print()
print("Done.")
