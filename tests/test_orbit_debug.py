#!/usr/bin/env python3
"""
Debug orbit processing issue
"""

import sardine
import sys
from pathlib import Path

def test_orbit_processing():
    """Test orbit processing step by step"""
    
    input_file = './data/S1A_IW_SLC__1SDV_20200103T170815_20200103T170842_030639_0382D5_DADE.zip'
    
    print("🔍 DEBUGGING ORBIT PROCESSING")
    print("=" * 50)
    
    # Check if input file exists
    if not Path(input_file).exists():
        print(f"❌ Input file not found: {input_file}")
        return
    else:
        print(f"✅ Input file found: {input_file}")
    
    try:
        # Test 1: Product info
        print("\n1. Testing product info extraction...")
        product_info = sardine.get_product_info(input_file)
        print(f"   ✅ Product info: {len(product_info)} fields")
        for key, value in list(product_info.items())[:5]:  # Show first 5 fields
            print(f"      {key}: {value}")
        
        # Test 2: Orbit file processing
        print("\n2. Testing orbit file processing...")
        
        # Extract parameters from filename
        filename = Path(input_file).name
        product_id = filename[:-4] if filename.endswith('.zip') else filename
        print(f"   Product ID: {product_id}")
        
        # Parse sensing time
        parts = filename.split('_')
        start_time_str = parts[5]  # Should be '20200103T170815'
        print(f"   Start time string: {start_time_str}")
        
        from datetime import datetime
        sensing_time = datetime.strptime(start_time_str, '%Y%m%dT%H%M%S')
        start_time_rfc3339 = sensing_time.isoformat() + 'Z'
        print(f"   RFC3339 time: {start_time_rfc3339}")
        
        # Test orbit processing
        orbit_result = sardine.apply_precise_orbit_file(
            product_id, 
            start_time_rfc3339, 
            '/tmp/orbit_cache'
        )
        
        print(f"   ✅ Orbit processing successful!")
        print(f"   Result keys: {list(orbit_result.keys())}")
        
        if 'result' in orbit_result:
            result = orbit_result['result']
            print(f"   Status: {result.get('status', 'unknown')}")
            print(f"   Orbit vectors: {result.get('orbit_vectors_count', 0)}")
            print(f"   Orbit type: {result.get('orbit_type', 'unknown')}")
            print(f"   Reference time: {result.get('reference_time', 'unknown')}")
            
            # Check for the "0 orbit state vectors" issue
            vector_count = result.get('orbit_vectors_count', 0)
            if vector_count == 0:
                print(f"   ❌ FOUND THE ISSUE: 0 orbit state vectors!")
            elif vector_count < 10:
                print(f"   ⚠️  WARNING: Only {vector_count} orbit vectors (minimum 10 required)")
            else:
                print(f"   ✅ Sufficient orbit vectors: {vector_count}")
        
        # Test 3: Check if there are any pipeline components that reset orbit data
        print("\n3. Testing IW split (next pipeline step)...")
        try:
            iw_result = sardine.iw_split_with_real_data(input_file, 'VV', 'IW1')
            print(f"   ✅ IW split successful")
            print(f"   Result keys: {list(iw_result.keys()) if isinstance(iw_result, dict) else 'Non-dict result'}")
        except Exception as e:
            print(f"   ❌ IW split failed: {e}")
            
        # Test 4: Check deburst processing
        print("\n4. Testing deburst processing...")
        try:
            deburst_result = sardine.deburst_topsar(input_file, 'IW1', 'VV')
            print(f"   ✅ Deburst successful")
            print(f"   Result keys: {list(deburst_result.keys()) if isinstance(deburst_result, dict) else 'Non-dict result'}")
        except Exception as e:
            print(f"   ❌ Deburst failed: {e}")
            
    except Exception as e:
        print(f"❌ Error during testing: {e}")
        import traceback
        traceback.print_exc()
        return
    
    print("\n" + "=" * 50)
    print("🎯 DEBUG COMPLETE")

if __name__ == "__main__":
    test_orbit_processing()
