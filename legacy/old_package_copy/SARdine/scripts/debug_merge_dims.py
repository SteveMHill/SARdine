#!/usr/bin/env python3
"""Debug the merge dimension issue."""

import os
import sys
import numpy as np

# Set strict parsing
os.environ['SARDINE_SERDE_ONLY'] = '1'
os.environ['SARDINE_REQUIRE_SUBSWATHS'] = '1'
os.environ['RUST_LOG'] = 'debug'

from pathlib import Path

SAFE_PATH = Path("/home/datacube/apps/SARdine/data/SLC/S1A_IW_SLC__1SDV_20201005T170824_20201005T170851_034664_04098A_1E66.SAFE")
POLARIZATION = "VV"

import sardine

print("=" * 70)
print("DEBUG: Merge Dimension Analysis")
print("=" * 70)

# Create reader
reader = sardine.create_cached_slc_reader(str(SAFE_PATH))
metadata = reader.get_cached_metadata()
subswaths_by_pol = reader.get_all_iw_subswaths()
subswaths_dict = subswaths_by_pol[POLARIZATION]

print("\n=== Subswath Metadata from Cached Reader ===")
for sw_name, sw_info in sorted(subswaths_dict.items()):
    print(f"\n{sw_name}:")
    print(f"  first_line_global: {sw_info.get('first_line_global', 'N/A')}")
    print(f"  last_line_global: {sw_info.get('last_line_global', 'N/A')}")
    print(f"  first_sample_global: {sw_info.get('first_sample_global', 'N/A')}")
    print(f"  last_sample_global: {sw_info.get('last_sample_global', 'N/A')}")
    print(f"  range_samples: {sw_info.get('range_samples', 'N/A')}")
    print(f"  azimuth_samples: {sw_info.get('azimuth_samples', 'N/A')}")
    print(f"  burst_count: {sw_info.get('burst_count', 'N/A')}")

# Run deburst to get actual data dimensions
print("\n=== Deburst Actual Data Dimensions ===")
deburst_data = {}
for sw_name in ['IW1', 'IW2', 'IW3']:
    result = sardine.deburst_topsar_cached(reader, sw_name, POLARIZATION)
    if isinstance(result, dict):
        print(f"{sw_name} result keys: {list(result.keys())}")
        if 'data' in result:
            data = result['data']
            print(f"{sw_name}: {data.shape} (rows={data.shape[0]}, cols={data.shape[1]})")
            deburst_data[sw_name] = data
        elif result.get('status') == 'error':
            print(f"{sw_name}: ERROR - {result.get('message', 'Unknown')}")
        else:
            # Look for any numpy array
            for k, v in result.items():
                if isinstance(v, np.ndarray):
                    print(f"{sw_name}: Found data in '{k}' - shape={v.shape}")
                    deburst_data[sw_name] = v
                    break
            else:
                print(f"{sw_name}: No numpy array found in dict")
    else:
        print(f"{sw_name}: Unexpected type {type(result)}")

# Calibrate
print("\n=== Calibrated Data Dimensions ===")
cal_data = {}
for sw_name, data in deburst_data.items():
    cal_result = sardine.radiometric_calibration(
        str(SAFE_PATH), sw_name, POLARIZATION, 'sigma0', data
    )
    if isinstance(cal_result, dict) and 'data' in cal_result:
        cdata = cal_result['data']
        print(f"{sw_name}: {cdata.shape}")
        cal_data[sw_name] = cdata

# Now merge with verbose logging
print("\n=== Merge Input Summary ===")
for sw_name, data in cal_data.items():
    print(f"{sw_name}: shape={data.shape}, dtype={data.dtype}")

print("\n=== Calling topsar_merge_cached ===")
merge_result = sardine.topsar_merge_cached(cal_data, POLARIZATION, reader)

if isinstance(merge_result, dict):
    print(f"Result keys: {list(merge_result.keys())}")
    if 'data' in merge_result:
        merged = merge_result['data']
        print(f"Merged shape: {merged.shape}")
        print(f"Expected: (~12300 rows, ~68800 cols)")
        print(f"Got: ({merged.shape[0]} rows, {merged.shape[1]} cols)")
        
        if merged.shape[0] < 5000:
            print("\n⚠️  SEVERE ISSUE: Row count way too low!")
            print("This indicates the merge is either:")
            print("  1. Truncating data based on incorrect metadata")
            print("  2. Using overlap extent instead of full azimuth extent")
            print("  3. Misinterpreting global vs local coordinates")
