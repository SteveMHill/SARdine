#!/usr/bin/env python3
"""
Quick test of terrain correction fixes.
"""

import numpy as np
import sardine
import os

def test_terrain_fixes():
    """Test if terrain correction fixes resolved the all-zero issue"""
    
    print("🧪 Testing Terrain Correction Fixes")
    
    try:
        # Create small test data
        test_data = np.random.rand(100, 100).astype(np.float32) * 1000.0 + 100.0
        test_bbox = [9.5, 48.5, 9.7, 48.7]  # Small area covering our scene
        
        print(f"Input data: {test_data.shape}, range: {np.min(test_data):.1f} to {np.max(test_data):.1f}")
        print(f"Test bbox: {test_bbox}")
        
        # Load orbit data
        orbit_result = sardine.load_orbit_file(
            '/home/datacube/apps/SARdine/orbit_files/S1A_IW_SLC__1SDV_20240720T052639_20240720T052707_054840_06ADB3_57D2_clean_POEORB.EOF'
        )
        
        annotation_path = '/tmp/s1a-iw1-slc-vv-20240720t052639-20240720t052707-054840-06adb3-004.xml'
        
        if not os.path.exists(annotation_path):
            print(f"❌ Annotation file not found: {annotation_path}")
            return False
        
        # Test terrain correction
        print("🌍 Testing terrain correction...")
        result = sardine.apply_terrain_correction_with_real_orbits(
            test_data,
            test_bbox,
            orbit_result,
            annotation_path,
            '/tmp/sardine_cache',
            50.0  # 50m resolution for faster processing
        )
        
        if isinstance(result, dict) and 'data' in result:
            output_data = result['data']
            
            print(f"\n✅ Terrain correction completed!")
            print(f"Output shape: {output_data.shape}")
            print(f"Output range: {np.min(output_data):.3f} to {np.max(output_data):.3f}")
            
            finite_count = np.sum(np.isfinite(output_data))
            nonzero_count = np.count_nonzero(output_data)
            
            print(f"Finite values: {finite_count}/{output_data.size}")
            print(f"Non-zero values: {nonzero_count}/{output_data.size}")
            
            if 'valid_pixel_percentage' in result:
                print(f"Valid percentage: {result['valid_pixel_percentage']:.1f}%")
            
            # Check if we resolved the all-zero issue
            if finite_count > 0 and nonzero_count > 0:
                print("🎉 SUCCESS: Terrain correction now preserves data values!")
                return True
            elif finite_count > 0:
                print("⚠️  PARTIAL: Finite values but all zeros")
                return False  
            else:
                print("❌ FAILED: Still all NaN/infinite values")
                return False
                
        else:
            print(f"❌ Invalid result format: {type(result)}")
            return False
            
    except Exception as e:
        print(f"❌ Error: {e}")
        import traceback
        traceback.print_exc()
        return False

if __name__ == '__main__':
    success = test_terrain_fixes()
    print(f"\n{'✅ FIXES SUCCESSFUL' if success else '❌ FIXES NOT COMPLETE'}")
    
    if success:
        print("🚀 Ready to test full backscatter pipeline!")
    else:
        print("🔧 More debugging needed...")
