#!/usr/bin/env python3
"""Create visualization with proper aspect ratio and annotated artifacts."""

import rasterio
import numpy as np
from PIL import Image, ImageDraw, ImageFont
from pathlib import Path

output_dir = Path("/home/datacube/apps/SARdine/SARdine/outputs/early_pipeline_test")

print("Creating annotated visualization...")

# Load multilooked dB data
with rasterio.open(output_dir / "04_multilooked_VV_dB.tif") as src:
    data = src.read(1)
    print(f"Multilooked shape: {data.shape}")

# Normalize to 0-255
vmin, vmax = -25, 5
normalized = np.clip(data, vmin, vmax)
normalized = (normalized - vmin) / (vmax - vmin) * 255
normalized = normalized.astype(np.uint8)
normalized[~np.isfinite(data)] = 0

# Create image with proper SAR orientation
# SAR convention: azimuth (flight direction) is vertical, range is horizontal
# Our data has rows=azimuth, cols=range, which is correct
# But for display, we might want to rotate or adjust aspect

# Original dimensions
h, w = data.shape
print(f"Original: {h} rows (azimuth) × {w} cols (range)")

# SAR pixel spacings (approximate for S1 IW after 3x3 multilook)
# Original: ~2.3m range × ~14m azimuth
# After 3x3 multilook: ~7m range × ~42m azimuth
range_spacing = 7.0  # meters
azimuth_spacing = 42.0  # meters

# To display with correct aspect ratio, we need to scale
# If we want square pixels in display, we need to stretch azimuth
aspect_correction = azimuth_spacing / range_spacing  # ~6

print(f"Aspect correction factor: {aspect_correction:.1f}")
print(f"  (azimuth pixels need to be {aspect_correction:.1f}x taller)")

# Create image with aspect-corrected dimensions
# Stretch height by aspect ratio
display_height = int(h * aspect_correction)
display_width = w

# But that would make the image huge, so let's scale down
max_dimension = 3000  # max pixels in any dimension
scale = min(max_dimension / display_height, max_dimension / display_width)

final_height = int(display_height * scale)
final_width = int(display_width * scale)

print(f"Display dimensions: {final_height} × {final_width} (aspect-corrected)")

# Create PIL image and resize
img = Image.fromarray(normalized, mode='L')
img_resized = img.resize((final_width, final_height), Image.Resampling.LANCZOS)

# Convert to RGB for annotations
img_rgb = img_resized.convert('RGB')
draw = ImageDraw.Draw(img_rgb)

# Try to get a font for annotations
try:
    font = ImageFont.truetype("/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf", 16)
    small_font = ImageFont.truetype("/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf", 12)
except:
    font = ImageFont.load_default()
    small_font = font

# Calculate positions for overlap regions (scaled)
# Original overlap positions (from analysis):
# IW1-IW2: cols 21015-21570
# IW2-IW3: cols 44826-45103
col_scale = final_width / w
row_scale = final_height / h

# Mark IW1-IW2 overlap
iw12_start = int(21015 / 3 * col_scale)  # /3 for multilook
iw12_end = int(21570 / 3 * col_scale)
draw.rectangle([iw12_start, 0, iw12_end, final_height], outline='red', width=2)
draw.text((iw12_start + 5, 10), "IW1-IW2 overlap", fill='red', font=small_font)

# Mark IW2-IW3 overlap
iw23_start = int(44826 / 3 * col_scale)
iw23_end = int(45103 / 3 * col_scale)
draw.rectangle([iw23_start, 0, iw23_end, final_height], outline='green', width=2)
draw.text((iw23_start + 5, 10), "IW2-IW3 overlap", fill='green', font=small_font)

# Mark burst boundaries (approximately every 1/9 of height)
burst_height = final_height // 9
for i in range(1, 9):
    y = i * burst_height
    draw.line([(0, y), (100, y)], fill='yellow', width=1)
    draw.text((5, y - 15), f"Burst {i}/{i+1}", fill='yellow', font=small_font)

# Add title and legend
draw.text((final_width - 200, 10), "SAR Artifacts:", fill='white', font=font)
draw.text((final_width - 200, 30), "• Red: IW overlaps", fill='red', font=small_font)
draw.text((final_width - 200, 45), "• Yellow: Burst boundaries", fill='yellow', font=small_font)
draw.text((final_width - 200, 60), "• Dark H-lines: Burst gaps", fill='cyan', font=small_font)

# Add scale info
draw.text((10, final_height - 30), f"Aspect-corrected ({azimuth_spacing:.0f}m × {range_spacing:.0f}m pixels)", fill='white', font=small_font)

# Save
output_path = output_dir / "04_multilooked_VV_annotated.png"
img_rgb.save(output_path)
print(f"Saved: {output_path}")

# Also save a version without aspect correction for comparison
img_original = Image.fromarray(normalized, mode='L')
scale2 = min(3000 / h, 3000 / w)
img_original_scaled = img_original.resize((int(w * scale2), int(h * scale2)), Image.Resampling.LANCZOS)
output_path2 = output_dir / "04_multilooked_VV_original_aspect.png"
img_original_scaled.save(output_path2)
print(f"Saved (original aspect): {output_path2}")

print("\n=== SUMMARY OF ISSUES ===")
print("""
1. STRETCHED APPEARANCE:
   - SAR pixels are non-square: ~7m range × ~42m azimuth (after multilook)
   - The aspect-corrected image compensates for this
   - Real ground coverage should look more square

2. DARK HORIZONTAL STRIPES (Burst Seams):
   - Gaps at burst boundaries in deburst output
   - Rows like 2769-2797 (29 lines), 4123-4173 (51 lines) have missing data
   - This is a deburst algorithm issue - bursts not properly joined

3. BRIGHT VERTICAL STRIPES (Overlap Seams):
   - At IW1-IW2 and IW2-IW3 boundaries
   - Radiometric inconsistency between subswaths
   - The blend is averaging different calibration levels

4. LAST BURST DISLOCATION:
   - Last burst in each subswath has 77-94 pixel range offset
   - Valid sample range starts later than previous bursts
   - This is a geometry/timing issue in burst metadata
""")
