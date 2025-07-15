#!/usr/bin/env python3
"""
Example: Using the Python API for Sentinel-1 Orbit File Handling

This example demonstrates the complete orbit file handling workflow
through the Python API, including status checking, downloading,
and accessing state vector data.
"""

import os
import sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import sardine
import tempfile


def main():
    print("🛰️  SARdine Python Orbit API Example")
    print("=" * 50)
    
    # Path to the example SLC file
    slc_path = "data/S1A_IW_SLC__1SDV_20200103T170815_20200103T170842_030639_0382D5_DADE.zip"
    
    if not os.path.exists(slc_path):
        print(f"❌ SLC file not found: {slc_path}")
        print("Please make sure the test data is available.")
        return
    
    try:
        # 1. Create SLC reader
        print(f"📖 Loading SLC file: {os.path.basename(slc_path)}")
        reader = sardine.SlcReader(slc_path)
        print("✅ SLC reader created successfully")
        
        # 2. Get basic metadata
        print("\n📋 Getting metadata...")
        try:
            metadata = reader.get_metadata("VV")
            print(f"   Product ID: {metadata.product_id}")
            print(f"   Mission: {metadata.mission}")
            print(f"   Start time: {metadata.start_time}")
            print(f"   Polarizations: {metadata.polarizations}")
        except Exception as e:
            print(f"   Warning: Could not get metadata - {e}")
        
        # 3. Check orbit status
        print("\n🔍 Checking orbit status...")
        with tempfile.TemporaryDirectory() as temp_cache:
            status = reader.check_orbit_status(temp_cache)
            print(f"   Status: {status}")
            
            # 4. Download orbit files if needed
            print("\n⬇️  Downloading orbit files...")
            try:
                downloaded_files = reader.download_orbit_files(temp_cache)
                print(f"   Downloaded {len(downloaded_files)} orbit file(s):")
                for file_path in downloaded_files:
                    print(f"     - {os.path.basename(file_path)}")
            except Exception as e:
                print(f"   ⚠️  Download failed: {e}")
                print("   (This may be expected if offline or if files don't exist)")
                return
            
            # 5. Get orbit data
            print("\n📡 Getting orbit data...")
            try:
                orbit_data = reader.get_orbit_data(temp_cache)
                print(f"   Orbit data: {orbit_data}")
                print(f"   Number of state vectors: {len(orbit_data)}")
                print(f"   Reference time: {orbit_data.reference_time}")
                
                # 6. Access state vectors
                print("\n🧭 Examining state vectors...")
                state_vectors = orbit_data.state_vectors
                if state_vectors:
                    print(f"   Total state vectors: {len(state_vectors)}")
                    
                    # Show first few state vectors
                    for i, sv in enumerate(state_vectors[:3]):
                        print(f"   Vector {i+1}:")
                        print(f"     Time: {sv.time}")
                        print(f"     Position (m): {sv.position}")
                        print(f"     Velocity (m/s): {sv.velocity}")
                    
                    if len(state_vectors) > 3:
                        print(f"   ... and {len(state_vectors) - 3} more vectors")
                
            except Exception as e:
                print(f"   ❌ Failed to get orbit data: {e}")
            
            # 7. Get satellite position at specific pixel
            print("\n🎯 Getting satellite position at specific pixel...")
            try:
                # Try to get position for a pixel in the middle of the image
                azimuth_line = 100
                range_sample = 500
                pos_vel = reader.get_satellite_position_at_pixel('VV', azimuth_line, range_sample, temp_cache)
                
                print(f"   Position and velocity at pixel ({azimuth_line}, {range_sample}):")
                print(f"     Position (m): ({pos_vel[0]:.2f}, {pos_vel[1]:.2f}, {pos_vel[2]:.2f})")
                print(f"     Velocity (m/s): ({pos_vel[3]:.2f}, {pos_vel[4]:.2f}, {pos_vel[5]:.2f})")
                
                # Calculate distance from Earth center
                import math
                earth_distance = math.sqrt(pos_vel[0]**2 + pos_vel[1]**2 + pos_vel[2]**2)
                print(f"     Distance from Earth center: {earth_distance/1000:.2f} km")
                
            except Exception as e:
                print(f"   ⚠️  Failed to get satellite position: {e}")
    
    except Exception as e:
        print(f"❌ Error: {e}")
        return
    
    print("\n🎉 Orbit API example completed successfully!")
    print("\nKey features demonstrated:")
    print("  ✓ Orbit status checking")
    print("  ✓ Orbit file downloading")
    print("  ✓ State vector access")
    print("  ✓ Satellite position interpolation")


if __name__ == "__main__":
    main()
