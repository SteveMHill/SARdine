#!/usr/bin/env python3
"""Analyze artifacts in merged SAR data: burst seams, overlap stripes, and geometry."""

import rasterio
import numpy as np
from pathlib import Path
import matplotlib
matplotlib.use('Agg')  # Non-interactive backend
import matplotlib.pyplot as plt

output_dir = Path("/home/datacube/apps/SARdine/SARdine/outputs/early_pipeline_test")

print("=" * 70)
print("ARTIFACT ANALYSIS")
print("=" * 70)

# Load the merged dB data
with rasterio.open(output_dir / "03_merged_VV_dB.tif") as src:
    merged_db = src.read(1)
    print(f"\nMerged image: {merged_db.shape} (rows=azimuth, cols=range)")

# =============================================================================
# 1. ANALYZE HORIZONTAL STRIPES (burst boundaries)
# =============================================================================
print("\n" + "=" * 50)
print("1. HORIZONTAL STRIPES (Burst Boundaries)")
print("=" * 50)

# Calculate row-wise statistics to find discontinuities
row_means = np.nanmean(merged_db, axis=1)
row_stds = np.nanstd(merged_db, axis=1)

# Find sudden jumps in row means (potential burst boundaries)
row_diffs = np.diff(row_means)
threshold = np.nanstd(row_diffs) * 3

jump_rows = np.where(np.abs(row_diffs) > threshold)[0]
print(f"Found {len(jump_rows)} potential burst boundary rows (>3σ jumps)")
if len(jump_rows) > 0:
    print(f"First 10 jump locations: {jump_rows[:10]}")
    
    # Expected burst boundaries: ~1370 lines per burst for 9 bursts
    expected_burst_height = merged_db.shape[0] / 9
    print(f"Expected burst height: ~{expected_burst_height:.0f} lines")
    
    # Check spacing between jumps
    if len(jump_rows) > 1:
        spacings = np.diff(jump_rows)
        print(f"Spacing between jumps: min={spacings.min()}, max={spacings.max()}, median={np.median(spacings):.0f}")

# Find dark rows (potential gaps)
dark_threshold = np.nanpercentile(row_means, 5)  # 5th percentile
dark_rows = np.where(row_means < dark_threshold)[0]
print(f"\nDark rows (<5th percentile): {len(dark_rows)} rows")
if len(dark_rows) > 0:
    # Group consecutive dark rows
    groups = []
    start = dark_rows[0]
    end = dark_rows[0]
    for i in range(1, len(dark_rows)):
        if dark_rows[i] == end + 1:
            end = dark_rows[i]
        else:
            groups.append((start, end))
            start = dark_rows[i]
            end = dark_rows[i]
    groups.append((start, end))
    
    print(f"Dark row groups (potential burst gaps):")
    for start, end in groups[:15]:
        print(f"  Rows {start}-{end} (height: {end-start+1})")

# =============================================================================
# 2. ANALYZE VERTICAL STRIPES (subswath overlaps)
# =============================================================================
print("\n" + "=" * 50)
print("2. VERTICAL STRIPES (Subswath Overlaps)")
print("=" * 50)

# Calculate column-wise statistics
col_means = np.nanmean(merged_db, axis=0)
col_stds = np.nanstd(merged_db, axis=0)

# Find bright columns
bright_threshold = np.nanpercentile(col_means, 98)  # 98th percentile
bright_cols = np.where(col_means > bright_threshold)[0]

print(f"Bright columns (>98th percentile): {len(bright_cols)} columns")

# Expected overlap regions based on subswath widths
# IW1: 21571, IW2: 25434, IW3: 24531
# Overlaps should be around column 21571 and 21571+25434-overlap
print("\nLooking for overlap regions...")

# Find contiguous bright column regions
if len(bright_cols) > 0:
    groups = []
    start = bright_cols[0]
    end = bright_cols[0]
    for i in range(1, len(bright_cols)):
        if bright_cols[i] <= end + 5:  # Allow small gaps
            end = bright_cols[i]
        else:
            groups.append((start, end))
            start = bright_cols[i]
            end = bright_cols[i]
    groups.append((start, end))
    
    # Filter to significant groups
    significant_groups = [(s, e) for s, e in groups if e - s > 50]
    print(f"Significant bright column regions (>50 cols wide):")
    for start, end in significant_groups:
        width = end - start + 1
        mean_val = np.nanmean(col_means[start:end+1])
        print(f"  Columns {start}-{end} (width: {width}, mean dB: {mean_val:.1f})")

# =============================================================================
# 3. CHECK INDIVIDUAL SUBSWATHS FOR BURST ARTIFACTS
# =============================================================================
print("\n" + "=" * 50)
print("3. INDIVIDUAL SUBSWATH BURST ANALYSIS")
print("=" * 50)

for iw in [1, 2, 3]:
    tif_path = output_dir / f"01_deburst_IW{iw}_VV.tif"
    if tif_path.exists():
        with rasterio.open(tif_path) as src:
            data = src.read(1)
            
        # Convert intensity to dB for analysis
        with np.errstate(divide='ignore', invalid='ignore'):
            db_data = 10 * np.log10(np.maximum(data, 1e-10))
        
        print(f"\nIW{iw}: {data.shape}")
        
        # Row statistics
        row_means = np.nanmean(db_data, axis=1)
        row_diffs = np.diff(row_means)
        
        # Find burst boundaries (sudden jumps)
        threshold = np.nanstd(row_diffs) * 2.5
        jumps = np.where(np.abs(row_diffs) > threshold)[0]
        
        if len(jumps) > 0:
            print(f"  Potential burst boundaries at rows: {jumps[:10]}")
            if len(jumps) > 1:
                spacings = np.diff(jumps)
                print(f"  Spacing: median={np.median(spacings):.0f}, std={np.std(spacings):.1f}")
        
        # Check for dark rows (gaps)
        dark_thresh = np.nanpercentile(row_means, 3)
        dark_rows = np.where(row_means < dark_thresh)[0]
        if len(dark_rows) > 0:
            print(f"  Dark rows: {len(dark_rows)} (potential gaps)")

# =============================================================================
# 4. LAST BURST DISLOCATION CHECK
# =============================================================================
print("\n" + "=" * 50)
print("4. LAST BURST GEOMETRY CHECK")
print("=" * 50)

# Check the bottom portion of each subswath
for iw in [1, 2, 3]:
    tif_path = output_dir / f"02_calibrated_IW{iw}_VV.tif"
    if tif_path.exists():
        with rasterio.open(tif_path) as src:
            data = src.read(1)
        
        rows = data.shape[0]
        burst_height = rows // 9
        
        # Get last burst (bottom ~1/9 of image)
        last_burst_start = rows - burst_height
        last_burst = data[last_burst_start:, :]
        second_last = data[last_burst_start-burst_height:last_burst_start, :]
        
        # Compare statistics
        last_valid = np.sum(np.isfinite(last_burst) & (last_burst > 0))
        second_valid = np.sum(np.isfinite(second_last) & (second_last > 0))
        
        print(f"\nIW{iw}:")
        print(f"  Total rows: {rows}, estimated burst height: {burst_height}")
        print(f"  Second-to-last burst valid pixels: {second_valid}")
        print(f"  Last burst valid pixels: {last_valid}")
        
        if last_valid < second_valid * 0.8:
            print(f"  ⚠️  Last burst has {100*(1-last_valid/second_valid):.0f}% fewer valid pixels!")
        
        # Check for horizontal offset in last burst
        # Compare valid sample ranges
        second_valid_cols = np.where(np.any(np.isfinite(second_last) & (second_last > 0), axis=0))[0]
        last_valid_cols = np.where(np.any(np.isfinite(last_burst) & (last_burst > 0), axis=0))[0]
        
        if len(second_valid_cols) > 0 and len(last_valid_cols) > 0:
            second_start, second_end = second_valid_cols[0], second_valid_cols[-1]
            last_start, last_end = last_valid_cols[0], last_valid_cols[-1]
            
            if abs(second_start - last_start) > 10 or abs(second_end - last_end) > 10:
                print(f"  ⚠️  Range offset detected!")
                print(f"      Second-last burst: cols {second_start}-{second_end}")
                print(f"      Last burst: cols {last_start}-{last_end}")
                print(f"      Offset: start={last_start-second_start}, end={last_end-second_end}")

# =============================================================================
# 5. CREATE DIAGNOSTIC PLOTS
# =============================================================================
print("\n" + "=" * 50)
print("5. CREATING DIAGNOSTIC PLOTS")
print("=" * 50)

fig, axes = plt.subplots(2, 2, figsize=(16, 12))

# Plot 1: Row means (horizontal profile)
ax1 = axes[0, 0]
with rasterio.open(output_dir / "03_merged_VV_dB.tif") as src:
    merged_db = src.read(1)
row_means = np.nanmean(merged_db, axis=1)
ax1.plot(row_means, linewidth=0.5)
ax1.set_xlabel('Row (azimuth)')
ax1.set_ylabel('Mean dB')
ax1.set_title('Row-wise Mean (shows burst boundaries as dips)')
ax1.grid(True, alpha=0.3)

# Mark expected burst boundaries
burst_height = len(row_means) // 9
for i in range(1, 9):
    ax1.axvline(i * burst_height, color='r', linestyle='--', alpha=0.5, label='Expected burst boundary' if i==1 else '')
ax1.legend()

# Plot 2: Column means (vertical profile)
ax2 = axes[0, 1]
col_means = np.nanmean(merged_db, axis=0)
ax2.plot(col_means, linewidth=0.5)
ax2.set_xlabel('Column (range)')
ax2.set_ylabel('Mean dB')
ax2.set_title('Column-wise Mean (shows subswath overlaps as peaks)')
ax2.grid(True, alpha=0.3)

# Mark expected subswath boundaries (approximate)
# Based on IW1=21571, IW2=25434, IW3=24531 with overlaps
ax2.axvline(21000, color='r', linestyle='--', alpha=0.5, label='IW1-IW2 overlap region')
ax2.axvline(45000, color='g', linestyle='--', alpha=0.5, label='IW2-IW3 overlap region')
ax2.legend()

# Plot 3: Sample rows showing burst seams
ax3 = axes[1, 0]
# Take a vertical slice through the middle
mid_col = merged_db.shape[1] // 2
vertical_slice = merged_db[:, mid_col-100:mid_col+100]
vertical_mean = np.nanmean(vertical_slice, axis=1)
ax3.plot(vertical_mean, linewidth=0.5)
ax3.set_xlabel('Row (azimuth)')
ax3.set_ylabel('Mean dB (central 200 cols)')
ax3.set_title('Vertical Slice Through Center - Burst Seam Visibility')
ax3.grid(True, alpha=0.3)

# Plot 4: Sample columns showing overlap seams
ax4 = axes[1, 1]
# Take horizontal slices through regions with overlap artifacts
mid_row = merged_db.shape[0] // 2
horizontal_slice = merged_db[mid_row-100:mid_row+100, :]
horizontal_mean = np.nanmean(horizontal_slice, axis=0)
ax4.plot(horizontal_mean, linewidth=0.5)
ax4.set_xlabel('Column (range)')
ax4.set_ylabel('Mean dB (central 200 rows)')
ax4.set_title('Horizontal Slice Through Center - Overlap Seam Visibility')
ax4.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig(output_dir / 'artifact_analysis.png', dpi=150)
print(f"Saved diagnostic plot: {output_dir / 'artifact_analysis.png'}")

print("\n" + "=" * 70)
print("ANALYSIS COMPLETE")
print("=" * 70)
