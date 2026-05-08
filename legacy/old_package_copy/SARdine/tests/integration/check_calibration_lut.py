#!/usr/bin/env python3
"""Check calibration LUT values from Sentinel-1 data"""

import sys
sys.path.insert(0, '/home/datacube/apps/SARdine/SARdine/python')

import sardine
import numpy as np

# Path to SLC data
slc_path = "/home/datacube/apps/SARdine/data/SLC/S1B_IW_SLC__1SDV_20190123T053348_20190123T053415_014617_01B3D4_E833.SAFE"

print("🔍 Creating cached reader...")
reader = sardine.create_cached_slc_reader(slc_path)

print("📖 Getting metadata...")
metadata = reader.get_cached_metadata()

print("\n📊 Metadata keys:")
for key in metadata.keys():
    print(f"  - {key}")

# Get calibration vectors if available
if 'calibration_vectors' in metadata:
    cal_vec = metadata['calibration_vectors']
    print(f"\n🎯 Calibration vectors type: {type(cal_vec)}")
    print(f"Keys: {cal_vec.keys() if hasattr(cal_vec, 'keys') else 'N/A'}")

# Try to get calibration LUTs
print("\n🔧 Attempting to load calibration LUTs...")
try:
    # This calls the Rust function that prepares calibration LUTs
    result = sardine.prepare_calibration_luts(
        slc_path,
        polarization='VV',
        subswath='IW1',
        calibration_type='sigma0',
        enable_noise_removal=True,
    )
    print(f"✅ Got calibration LUTs: {type(result)}")
    print(f"Keys: {result.keys() if hasattr(result, 'keys') else 'N/A'}")
    
    # Check the values
    if isinstance(result, dict):
        for key, val in result.items():
            if isinstance(val, np.ndarray):
                print(f"\n📈 {key}:")
                print(f"  Shape: {val.shape}")
                print(f"  Min: {val.min():.6e}")
                print(f"  Max: {val.max():.6e}")
                print(f"  Mean: {val.mean():.6e}")
                print(f"  Median: {np.median(val):.6e}")
                # Sample a few values
                flat = val.flatten()
                mid_idx = len(flat) // 2
                print(f"  Sample values (mid): {flat[mid_idx:mid_idx+5]}")
except Exception as e:
    print(f"❌ Error loading calibration LUTs: {e}")
    import traceback
    traceback.print_exc()

# Let's also check what the XML says
print("\n\n🗂️  Checking raw XML calibration values...")
import os
import xml.etree.ElementTree as ET

# Find calibration XML file
cal_dir = os.path.join(slc_path, "annotation", "calibration")
cal_files = [f for f in os.listdir(cal_dir) if 'vv' in f.lower()]
if cal_files:
    cal_file = os.path.join(cal_dir, cal_files[0])
    print(f"📄 Reading: {cal_file}")
    
    tree = ET.parse(cal_file)
    root = tree.getroot()
    
    # Get absoluteCalibrationConstant
    abs_cal = root.find('.//absoluteCalibrationConstant')
    if abs_cal is not None:
        print(f"🔢 absoluteCalibrationConstant: {abs_cal.text}")
    
    # Get first calibrationVector
    cal_vec = root.find('.//calibrationVector')
    if cal_vec is not None:
        sigma = cal_vec.find('sigmaNought')
        beta = cal_vec.find('betaNought')
        gamma = cal_vec.find('gamma')
        
        if sigma is not None:
            sigma_vals = [float(x) for x in sigma.text.strip().split()]
            print(f"\n📊 Raw XML sigmaNought values:")
            print(f"  Count: {len(sigma_vals)}")
            print(f"  Min: {min(sigma_vals):.6e}")
            print(f"  Max: {max(sigma_vals):.6e}")
            print(f"  Mean: {np.mean(sigma_vals):.6e}")
            print(f"  Sample: {sigma_vals[len(sigma_vals)//2:len(sigma_vals)//2+5]}")
            
            # Check if these look inverted
            print(f"\n🤔 Analysis:")
            mean_val = np.mean(sigma_vals)
            if mean_val > 100:
                print(f"  ❌ These are large (~{mean_val:.1f}), suggesting they are K values (not inverted)")
                print(f"  Expected after inversion: ~{1.0/mean_val:.6e}")
            elif mean_val < 0.01:
                print(f"  ✅ These are small (~{mean_val:.6e}), suggesting they are already 1/K (inverted)")
            else:
                print(f"  ⚠️  Unclear - mean is {mean_val:.6e}")
