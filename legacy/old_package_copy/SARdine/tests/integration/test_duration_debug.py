#!/usr/bin/env python3
"""Quick test script to debug the product_duration issue in terrain_correction."""

import os
import sys

# Set environment variables BEFORE importing sardine
os.environ["SARDINE_SERDE_ONLY"] = "1"
os.environ["SARDINE_REQUIRE_SUBSWATHS"] = "1"
os.environ["RUST_LOG"] = "info"

import sardine
import numpy as np
import json

# Test SAFE path
safe_path = "/home/datacube/apps/SARdine/data/SLC/S1A_IW_SLC__1SDV_20201005T170824_20201005T170851_034664_04098A_1E66.SAFE"
output_dir = "/home/datacube/apps/SARdine/tmp_out_duration_test"
os.makedirs(output_dir, exist_ok=True)

print("=" * 70)
print("Testing product_duration metadata handoff")
print("=" * 70)

# Create cached reader
print("\n1. Creating cached SLC reader...")
reader = sardine.create_cached_slc_reader(safe_path)
metadata = reader.get_cached_metadata()

print(f"\n2. Cached metadata keys ({len(metadata)} items):")
timing_keys = [k for k in metadata.keys() if 'duration' in k.lower() or 'time' in k.lower() or 'start' in k.lower() or 'stop' in k.lower()]
for key in sorted(timing_keys):
    value = metadata.get(key)
    if isinstance(value, (int, float)):
        print(f"   {key}: {value}")
    elif isinstance(value, str) and len(value) < 100:
        print(f"   {key}: {value}")

# Get the subswath metadata
print("\n3. Getting IW subswaths...")
subswaths = reader.get_all_iw_subswaths()
print(f"   Found {len(subswaths)} subswaths")

for sw_name, sw_data in sorted(subswaths.items()):
    if isinstance(sw_data, dict):
        prf = sw_data.get('prf_hz')
        lines = sw_data.get('num_lines')
        azimuth_dt = sw_data.get('azimuth_time_interval')
        print(f"   {sw_name}: PRF={prf} Hz, lines={lines}, azimuth_dt={azimuth_dt}")
        if prf and lines:
            duration_from_lines = float(lines) / float(prf)
            print(f"      → duration from lines = {duration_from_lines:.3f}s")

# Extract key timing values
print("\n4. Constructing real_metadata dict (like backscatter.py does)...")

# Get orbit vectors
orbit_str = metadata.get('orbit_state_vectors_json')
if orbit_str:
    import json
    orbits = json.loads(orbit_str)
    if orbits:
        first_orbit = orbits[0]
        print(f"   First orbit vector time: {first_orbit.get('time')}")
        print(f"   First orbit vector UTC_s (approx): {first_orbit.get('time_utc_seconds', 'N/A')}")

product_start_time = metadata.get('product_start_time')
product_stop_time = metadata.get('product_stop_time')
print(f"\n   product_start_time (str): {product_start_time}")
print(f"   product_stop_time (str): {product_stop_time}")

# Try to get the numeric values that backscatter.py would compute
# This is where the bug likely is
number_of_lines = metadata.get('number_of_lines')
print(f"\n   number_of_lines: {number_of_lines}")

# Simulate what the Rust side would receive
real_metadata = {
    'product_duration': 25.0,  # Expected: ~25 seconds
    'product_stop_rel_s': 91.0,  # Expected: product_start_rel_s + ~25s
    'product_start_time_abs': 1601917704.0,  # Example
    'product_stop_time_abs': 1601917729.0,  # Example
    'total_azimuth_lines': number_of_lines,  # BUG: If this gets confused with duration!
    'number_of_lines': number_of_lines,
}

print("\n5. Checking for confusion between lines and duration...")
print(f"   number_of_lines = {number_of_lines}")
print(f"   If confused, product_duration would be ~{number_of_lines} instead of ~25s")

# KEY CHECK: Does any metadata field contain the wrong value?
for key in ['product_duration', 'product_stop_time_abs', 'product_stop_rel_s']:
    val = metadata.get(key)
    if val is not None:
        if isinstance(val, (int, float)) and 13000 < val < 14000:
            print(f"\n   ⚠️  WARNING: {key} = {val} looks like number_of_lines!")
        else:
            print(f"   {key} = {val}")

print("\n" + "=" * 70)
print("Done - check output above for anomalies")
print("=" * 70)
