#!/usr/bin/env python3
"""Deep diagnosis of burst scalloping and calibration LUT alignment."""

import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'python'))

os.environ['SARDINE_SERDE_ONLY'] = '1'
os.environ['SARDINE_REQUIRE_SUBSWATHS'] = '1'

import sardine
import numpy as np
import rasterio
from lxml import etree

SAFE = '/home/datacube/apps/SARdine/data/SLC/S1A_IW_SLC__1SDV_20201005T170824_20201005T170851_034664_04098A_1E66.SAFE'

# 1. Get burst metadata
reader = sardine.create_cached_slc_reader(SAFE)

# Get burst metadata from annotation XML directly
print("="*70)
print("BURST METADATA (from annotation XML)")
print("="*70)

import glob
ann_dir = os.path.join(SAFE, 'annotation')
ann_files = sorted(glob.glob(os.path.join(ann_dir, 's1a-iw*-slc-vv-*.xml')))
for af in ann_files:
    fname = os.path.basename(af)
    tree = etree.parse(af)
    root = tree.getroot()
    bursts = root.findall('.//burst')
    sw_name = fname.split('-')[1].upper()
    n_lines = root.findtext('.//numberOfLines', '?')
    n_samples = root.findtext('.//numberOfSamples', '?')
    print(f"\n{sw_name} ({fname}): {n_lines} x {n_samples}, {len(bursts)} bursts")
    for i, b in enumerate(bursts):
        fl = b.findtext('firstValidSample', '')
        ll = b.findtext('lastValidSample', '')
        az_time = b.findtext('azimuthTime', '')
        # Count valid samples
        fvs = [int(x) for x in fl.split()] if fl else []
        lvs = [int(x) for x in ll.split()] if ll else []
        # firstValidSample/lastValidSample are per-line within the burst
        n_burst_lines = len(fvs) if fvs else 0
        valid_lines = sum(1 for f in fvs if f >= 0) if fvs else 0
        print(f"    Burst {i}: {n_burst_lines} lines, {valid_lines} valid, time={az_time}")

# 2. Parse calibration annotation XML to check vector density
print("\n" + "="*70)
print("CALIBRATION VECTOR DENSITY")
print("="*70)

import glob
import zipfile

# Find calibration files
annotation_dir = os.path.join(SAFE, 'annotation', 'calibration')
cal_files = sorted(glob.glob(os.path.join(annotation_dir, 'calibration-*.xml')))
print(f"Found {len(cal_files)} calibration files")

for cf in cal_files:
    fname = os.path.basename(cf)
    # Extract subswath and pol from filename
    parts = fname.lower()
    tree = etree.parse(cf)
    root = tree.getroot()
    
    cal_vectors = root.findall('.//calibrationVector')
    print(f"\n  {fname}:")
    print(f"    Total vectors: {len(cal_vectors)}")
    
    if cal_vectors:
        lines = []
        for cv in cal_vectors:
            line = int(cv.findtext('line', '0'))
            lines.append(line)
            # Check number of pixel values in first vector
            if cv == cal_vectors[0]:
                pixels = cv.findtext('pixel', '')
                n_pixels = len(pixels.split()) if pixels else 0
                sigma = cv.findtext('sigmaNought', '')
                n_sigma = len(sigma.split()) if sigma else 0
                print(f"    Pixels per vector: {n_pixels}, sigma0 values: {n_sigma}")
        
        lines = np.array(lines)
        spacings = np.diff(lines)
        print(f"    Lines range: {lines[0]} - {lines[-1]}")
        print(f"    Spacings: min={spacings.min()}, max={spacings.max()}, mean={spacings.mean():.1f}")
        
        # Show all vector lines
        if len(lines) <= 30:
            print(f"    All vector lines: {lines.tolist()}")

# 3. Analyze scalloping in debursted output
print("\n" + "="*70)
print("SCALLOPING ANALYSIS (per-subswath)")
print("="*70)

base = 'outputs/validation_run/S1A_IW_SLC__1SDV_20201005T170824_20201005T170851_034664_04098A_1E66'

for sw in ['iw1', 'iw2', 'iw3']:
    src = rasterio.open(f'{base}_{sw}_calibrated.tif')
    h, w = src.height, src.width
    
    # Read a wide swath at center for averaging
    strip_w = min(200, w // 4)
    data = src.read(1, window=rasterio.windows.Window(w//2 - strip_w//2, 0, strip_w, h))
    
    # Average across range for each azimuth line
    row_power = np.nanmean(data, axis=1)
    valid = np.isfinite(row_power) & (row_power > 0)
    row_db = np.full(h, np.nan)
    row_db[valid] = 10 * np.log10(row_power[valid])
    
    # Smooth heavily to see the scalloping envelope
    kernel = 51
    valid_db = row_db.copy()
    valid_db[~np.isfinite(valid_db)] = np.nanmedian(row_db[np.isfinite(row_db)])  # fill gaps
    smoothed = np.convolve(valid_db, np.ones(kernel)/kernel, mode='same')
    
    # Compute expected burst periodicity
    # IW mode: ~1350-1500 lines per burst, with overlap
    expected_period = h / 9  # roughly
    
    # Look for periodic modulation by computing autocorrelation
    detrended = valid_db - smoothed
    detrended[~np.isfinite(detrended)] = 0
    autocorr = np.correlate(detrended[:5000], detrended[:5000], mode='full')
    autocorr = autocorr[len(autocorr)//2:]  # positive lags only
    autocorr = autocorr / autocorr[0]  # normalize
    
    # Find peaks in autocorrelation (indicating periodic pattern)
    peaks = []
    for i in range(200, min(3000, len(autocorr) - 1)):
        if autocorr[i] > autocorr[i-1] and autocorr[i] > autocorr[i+1] and autocorr[i] > 0.05:
            peaks.append((i, autocorr[i]))
    
    print(f"\n{sw.upper()}: {h} x {w}")
    print(f"  Row power range: {np.nanmin(row_db):.1f} to {np.nanmax(row_db):.1f} dB")
    print(f"  Mean: {np.nanmean(row_db):.1f} dB, Std: {np.nanstd(row_db):.2f} dB")
    print(f"  Smoothed range: {np.min(smoothed):.1f} to {np.max(smoothed):.1f} dB (envelope variation: {np.max(smoothed)-np.min(smoothed):.1f} dB)")
    
    if peaks:
        print(f"  Autocorrelation peaks (periodic modulation):")
        for lag, corr in sorted(peaks, key=lambda x: -x[1])[:5]:
            print(f"    lag={lag} lines, correlation={corr:.3f}")
    else:
        print(f"  No significant periodic autocorrelation found")
    
    # Check specific burst boundaries (based on burst spacing)
    n_bursts = 9
    burst_spacing = h / n_bursts
    print(f"  Expected burst spacing: {burst_spacing:.0f} lines")
    print(f"  Power at expected burst boundaries:")
    for bi in range(1, n_bursts):
        boundary_row = int(bi * burst_spacing)
        lo = max(0, boundary_row - 75)
        hi = min(h, boundary_row + 75)
        center_lo = max(0, boundary_row - 10)
        center_hi = min(h, boundary_row + 10)
        
        # Power at boundary vs surroundings
        boundary_val = np.nanmean(row_db[center_lo:center_hi])
        surround_lo = np.nanmean(row_db[lo:center_lo])
        surround_hi = np.nanmean(row_db[center_hi:hi])
        surround_mean = (surround_lo + surround_hi) / 2
        delta = boundary_val - surround_mean
        flag = " <-- DIP!" if delta < -1.0 else ""
        print(f"    Boundary {bi} (row {boundary_row}): {boundary_val:.1f} dB, surroundings: {surround_mean:.1f} dB, delta: {delta:+.1f} dB{flag}")
    
    src.close()
