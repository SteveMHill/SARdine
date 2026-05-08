#!/usr/bin/env python3
"""Analyze the merged GeoTIFF to understand orientation and overlap issues."""

import rasterio
import numpy as np
from pathlib import Path

output_dir = Path("/home/datacube/apps/SARdine/SARdine/outputs/early_pipeline_test")

# Check the merged GeoTIFF
print("=" * 60)
print("MERGED GEOTIFF ANALYSIS")
print("=" * 60)

with rasterio.open(output_dir / '03_merged_VV_linear.tif') as src:
    print(f'\nShape: {src.count} band(s), {src.height} rows (height), {src.width} cols (width)')
    print(f'Transform: {src.transform}')
    data = src.read(1)
    print(f'NumPy array shape: {data.shape}')
    print(f'Data range: {np.nanmin(data):.4e} to {np.nanmax(data):.4e}')
    
    # In SAR convention:
    # - Rows = azimuth (along-track, should be tall)
    # - Columns = range (across-track, should be narrower)
    # So expected: height >> width for a typical S1 swath
    
    aspect_ratio = src.width / src.height
    print(f'\nAspect ratio (width/height): {aspect_ratio:.2f}')
    if aspect_ratio > 1:
        print("⚠️  Image is WIDER than TALL - unusual for SAR!")
        print("   SAR convention: rows=azimuth (tall), cols=range (narrow)")
        print("   This suggests the data may be transposed.")
    else:
        print("✓ Image is TALLER than WIDE - expected for SAR")
    
    # Analyze for bright columns/rows (potential overlap regions)
    print("\n" + "-" * 40)
    print("OVERLAP ANALYSIS (looking for bright stripes)")
    print("-" * 40)
    
    # Column-wise analysis (if image is transposed, overlaps would be columns)
    col_means = np.nanmean(data, axis=0)
    col_p50 = np.nanpercentile(col_means, 50)
    col_p99 = np.nanpercentile(col_means, 99)
    print(f'\nColumn means: median={col_p50:.4e}, 99th pctl={col_p99:.4e}')
    
    # Find exceptionally bright columns
    threshold = col_p50 * 3  # 3x median
    bright_cols = np.where(col_means > threshold)[0]
    if len(bright_cols) > 0:
        # Group consecutive columns
        groups = []
        start = bright_cols[0]
        end = bright_cols[0]
        for i in range(1, len(bright_cols)):
            if bright_cols[i] == end + 1:
                end = bright_cols[i]
            else:
                groups.append((start, end))
                start = bright_cols[i]
                end = bright_cols[i]
        groups.append((start, end))
        
        print(f"\n⚠️  Found {len(groups)} bright column region(s) (>3x median):")
        for start, end in groups[:5]:  # Show first 5
            width = end - start + 1
            print(f"   Columns {start}-{end} (width: {width} pixels)")
    
    # Row-wise analysis
    row_means = np.nanmean(data, axis=1)
    row_p50 = np.nanpercentile(row_means, 50)
    row_p99 = np.nanpercentile(row_means, 99)
    print(f'\nRow means: median={row_p50:.4e}, 99th pctl={row_p99:.4e}')
    
    bright_rows = np.where(row_means > row_p50 * 3)[0]
    if len(bright_rows) > 0:
        # Group consecutive rows
        groups = []
        start = bright_rows[0]
        end = bright_rows[0]
        for i in range(1, len(bright_rows)):
            if bright_rows[i] == end + 1:
                end = bright_rows[i]
            else:
                groups.append((start, end))
                start = bright_rows[i]
                end = bright_rows[i]
        groups.append((start, end))
        
        print(f"\n⚠️  Found {len(groups)} bright row region(s) (>3x median):")
        for start, end in groups[:5]:
            height = end - start + 1
            print(f"   Rows {start}-{end} (height: {height} pixels)")

# Check individual subswaths
print("\n" + "=" * 60)
print("INDIVIDUAL SUBSWATH DIMENSIONS")
print("=" * 60)

for iw in [1, 2, 3]:
    tif_path = output_dir / f'02_calibrated_IW{iw}_VV.tif'
    if tif_path.exists():
        with rasterio.open(tif_path) as src:
            print(f'\nIW{iw}: {src.height} rows × {src.width} cols')
            data = src.read(1)
            print(f'     Data range: {np.nanmin(data):.4e} to {np.nanmax(data):.4e}')
    else:
        print(f'\nIW{iw}: file not found')

# Expected subswath arrangement
print("\n" + "-" * 40)
print("EXPECTED S1 IW GEOMETRY:")
print("-" * 40)
print("""
Sentinel-1 IW mode:
- 3 subswaths: IW1 (near range) | IW2 (mid range) | IW3 (far range)
- Each subswath is tall (azimuth) × narrow (range)
- Overlaps occur in RANGE direction between IW1-IW2 and IW2-IW3
- After merging: combined image should be tall × wide(r)

If image appears wide × short, the axes are likely transposed.
If bright stripes appear vertically, overlaps weren't properly blended.
""")
