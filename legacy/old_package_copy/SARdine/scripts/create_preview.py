#!/usr/bin/env python3
"""Create a proper preview PNG from the fixed merged GeoTIFF."""

import rasterio
import numpy as np
from PIL import Image
from pathlib import Path

output_dir = Path("/home/datacube/apps/SARdine/SARdine/outputs/early_pipeline_test")

# Read the multilooked dB GeoTIFF
tif_path = output_dir / "04_multilooked_VV_dB.tif"
print(f"Reading: {tif_path}")

with rasterio.open(tif_path) as src:
    data = src.read(1)
    print(f"Shape: {data.shape} (height x width)")

# Data is in dB, typically -30 to +10 dB range
# Normalize to 0-255 for display
vmin, vmax = -25, 5  # Typical SAR dB range for good visualization
print(f"Normalizing dB values from [{vmin}, {vmax}] to [0, 255]")

# Clip and normalize
normalized = np.clip(data, vmin, vmax)
normalized = (normalized - vmin) / (vmax - vmin) * 255
normalized = normalized.astype(np.uint8)

# Handle NaN values (set to black)
normalized[~np.isfinite(data)] = 0

# Create PIL image (grayscale)
img = Image.fromarray(normalized, mode='L')
print(f"Image size: {img.size} (width x height)")

# For SAR, the data is typically displayed with:
# - X axis = Range (cross-track) - this is the width
# - Y axis = Azimuth (along-track) - this is the height

# The image should be taller than wide for a single swath,
# but merged IW image is wider because 3 subswaths are side-by-side.

# Save full resolution
output_png = output_dir / "04_multilooked_VV_preview_fixed.png"
img.save(output_png)
print(f"Saved full resolution: {output_png}")

# Create a downsampled version for easier viewing
max_width = 2000  # Reduce to reasonable size
if img.width > max_width:
    scale = max_width / img.width
    new_size = (int(img.width * scale), int(img.height * scale))
    img_small = img.resize(new_size, Image.Resampling.LANCZOS)
    small_path = output_dir / "04_multilooked_VV_preview_small.png"
    img_small.save(small_path)
    print(f"Saved downsampled ({new_size[0]}x{new_size[1]}): {small_path}")

print("\n=== Preview Statistics ===")
print(f"Original multilooked dimensions: {data.shape}")
print(f"Valid pixels: {np.sum(np.isfinite(data))}/{data.size} ({100*np.sum(np.isfinite(data))/data.size:.1f}%)")
print(f"dB range: [{np.nanmin(data):.1f}, {np.nanmax(data):.1f}]")
