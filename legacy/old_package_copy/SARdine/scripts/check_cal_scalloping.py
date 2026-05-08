#!/usr/bin/env python3
"""Check if calibration vectors capture azimuth variation (scalloping)."""

import numpy as np
from lxml import etree
import os, glob

SAFE = '/home/datacube/apps/SARdine/data/SLC/S1A_IW_SLC__1SDV_20201005T170824_20201005T170851_034664_04098A_1E66.SAFE'
cal_dir = os.path.join(SAFE, 'annotation', 'calibration')
cal_file = sorted(glob.glob(os.path.join(cal_dir, 'calibration-*iw2*vv*.xml')))[0]

tree = etree.parse(cal_file)
root = tree.getroot()
vectors = root.findall('.//calibrationVector')

print(f"File: {os.path.basename(cal_file)}")
print(f"Vectors: {len(vectors)}\n")

# Extract sigma0 values at a fixed range pixel across all vectors
# Use pixel index 325 (middle of range)
ref_pixel_idx = 325

prev_sigma = None
for i, v in enumerate(vectors):
    line = int(v.findtext('line', '0'))
    pixels_str = v.findtext('pixel', '')
    sigma_str = v.findtext('sigmaNought', '')
    
    pixels = [int(x) for x in pixels_str.split()]
    sigmas = [float(x) for x in sigma_str.split()]
    
    # Get sigma at reference pixel
    if ref_pixel_idx < len(sigmas):
        sig = sigmas[ref_pixel_idx]
        pix = pixels[ref_pixel_idx]
        
        delta_str = ""
        if prev_sigma is not None:
            delta = sig - prev_sigma
            pct = (delta / prev_sigma * 100) if prev_sigma != 0 else 0
            delta_str = f"  delta={delta:+.6f} ({pct:+.3f}%)"
        prev_sigma = sig
        
        # Show gain (1/K^2)
        gain_db = -20 * np.log10(sig) if sig > 0 else float('nan')
        print(f"  Vector {i:2d} line={line:6d}: sigma0[pix={pix}] = {sig:.6f} (gain={gain_db:.2f} dB){delta_str}")

# Now check: within each burst, how much does sigma0 vary?
print("\n\nPER-BURST SIGMA0 VARIATION (azimuth scalloping check)")
print("="*60)

# IW2 bursts: 9 bursts, ~1510 lines each. Lines 0..13589
# Burst boundaries approximately at 0, 1510, 3020, 4530, ...
burst_lines = 1510
n_bursts = 9

all_lines = [int(v.findtext('line', '0')) for v in vectors]

for bi in range(n_bursts):
    burst_start = bi * burst_lines
    burst_end = burst_start + burst_lines
    
    # Find vectors within this burst
    burst_vectors = [(l, v) for l, v in zip(all_lines, vectors) if burst_start <= l < burst_end]
    
    if burst_vectors:
        sigs_at_ref = []
        for l, v in burst_vectors:
            sigma_str = v.findtext('sigmaNought', '')
            sigmas = [float(x) for x in sigma_str.split()]
            if ref_pixel_idx < len(sigmas):
                sigs_at_ref.append((l, sigmas[ref_pixel_idx]))
        
        if len(sigs_at_ref) >= 2:
            vals = [s for _, s in sigs_at_ref]
            gains = [-20 * np.log10(v) for v in vals if v > 0]
            gain_range = max(gains) - min(gains) if gains else 0
            print(f"  Burst {bi} (lines {burst_start}-{burst_end}): "
                  f"{len(sigs_at_ref)} vectors, "
                  f"sigma0 range: [{min(vals):.4f}, {max(vals):.4f}], "
                  f"gain range: {gain_range:.3f} dB")
        elif sigs_at_ref:
            print(f"  Burst {bi}: only 1 vector")
    else:
        print(f"  Burst {bi}: no vectors")

# Compare sigma0 at burst centers vs burst edges
print("\n\nBURST CENTER vs EDGE COMPARISON")
print("="*60)

lines_arr = np.array(all_lines)
sigs_arr = []
for v in vectors:
    sigma_str = v.findtext('sigmaNought', '')
    sigmas = [float(x) for x in sigma_str.split()]
    sigs_arr.append(sigmas[ref_pixel_idx] if ref_pixel_idx < len(sigmas) else 0)
sigs_arr = np.array(sigs_arr)

# Identify vectors near burst centers and edges
for bi in range(n_bursts):
    burst_center = bi * burst_lines + burst_lines // 2
    burst_edge_start = bi * burst_lines
    burst_edge_end = (bi + 1) * burst_lines
    
    # Closest vector to center
    center_idx = np.argmin(np.abs(lines_arr - burst_center))
    center_sig = sigs_arr[center_idx]
    center_gain = -20 * np.log10(center_sig) if center_sig > 0 else float('nan')
    
    # Closest vector to edge
    edge_dists = np.minimum(np.abs(lines_arr - burst_edge_start), np.abs(lines_arr - burst_edge_end))
    edge_idx = np.argmin(edge_dists)
    edge_sig = sigs_arr[edge_idx]
    edge_gain = -20 * np.log10(edge_sig) if edge_sig > 0 else float('nan')
    
    print(f"  Burst {bi}: center_gain={center_gain:.2f} dB (line {lines_arr[center_idx]}), "
          f"edge_gain={edge_gain:.2f} dB (line {lines_arr[edge_idx]}), "
          f"diff={center_gain - edge_gain:+.3f} dB")
