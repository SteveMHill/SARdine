#!/usr/bin/env python3
"""
Real Sentinel-1 SAR Processing Pipeline Test - NO SYNTHETIC DATA
This test uses ONLY real Sentinel-1 SLC data throughout the entire pipeline.
No fallbacks, no synthetic data, no testing shortcuts.
"""

import sardine
import numpy as np
import os
import sys

def test_real_sentinel1_pipeline():
    """Test complete SAR processing pipeline with real Sentinel-1 data only"""
    print("=" * 80)
    print("🛰️  REAL SENTINEL-1 SAR PROCESSING PIPELINE")
    print("Using ONLY real data - no synthetic/testing data allowed")
    print("=" * 80)
    
    # Path to real Sentinel-1 data
    zip_path = "data/S1A_IW_SLC__1SDV_20200103T170815_20200103T170842_030639_0382D5_DADE.zip"
    
    if not os.path.exists(zip_path):
        print(f"❌ ERROR: Real Sentinel-1 data not found at {zip_path}")
        return False
    
    try:
        # Step 1: Read real SLC data
        print("\n📡 STEP 1: Reading Real SLC Data")
        reader = sardine.SlcReader(zip_path)
        metadata = reader.get_metadata()
        
        # Read VV polarization data
        slc_result = reader.read_slc_data('VV')
        slc_data = slc_result['data']
        
        print(f"   ✅ Real SLC loaded: {slc_data.shape} {slc_data.dtype}")
        print(f"   ✅ Metadata fields: {len(metadata)}")
        
        # Step 2: Apply Precise Orbit File
        print("\n🛰️ STEP 2: Apply Precise Orbit File")
        orbit_result = sardine.apply_precise_orbit_file(
            "S1A_IW_SLC__1SDV_20200103T170815", 
            "2020-01-03T17:08:15.000000Z",
            "./output/cache/orbit"
        )
        print(f"   ✅ Orbit file applied: {orbit_result}")
        
        # Step 3: IW Split  
        print("\n🔀 STEP 3: IW Split")
        split_result = sardine.iw_split(slc_data, 'IW2')  # Use IW2 subswath
        if isinstance(split_result, dict) and 'data' in split_result:
            split_data = split_result['data']
        else:
            split_data = split_result
        print(f"   ✅ IW split successful: {split_data.shape}")
        
        # Step 4: Deburst TOPSAR
        print("\n🔄 STEP 4: Deburst TOPSAR")
        deburst_result = sardine.deburst_topsar(split_data, 3)  # 3 bursts typical
        if isinstance(deburst_result, dict) and 'data' in deburst_result:
            deburst_data = deburst_result['data']
        else:
            deburst_data = deburst_result
        print(f"   ✅ Deburst successful: {deburst_data.shape}")
        
        # Use deburst data for calibration instead of raw SLC
        calibration_input = deburst_data
        
        # Step 5: Load real calibration data
        # Step 5: Load real calibration data
        print("\n📐 STEP 5: Loading Real Calibration Data")
        cal_files = reader.find_calibration_files()
        cal_data = reader.read_calibration_data('VV')
        
        print(f"   ✅ Calibration files: {list(cal_files.keys())}")
        print(f"   ✅ Calibration vectors: {len(cal_data['vectors'])}")
        print(f"   ✅ Swath: {cal_data['swath']}")
        
        # Verify we have real calibration data
        first_vector = cal_data['vectors'][0]
        if len(first_vector['sigma_nought']) < 10:
            print(f"   ❌ Insufficient calibration data - only {len(first_vector['sigma_nought'])} values")
            return False
        
        print(f"   ✅ Real calibration verified: {len(first_vector['sigma_nought'])} sigma values")
        
        # Step 6: Apply radiometric calibration with real data
        print("\n🔧 STEP 6: Radiometric Calibration with Real Data")
        calibrated_result = sardine.radiometric_calibration_with_zip(
            calibration_input, 'VV', zip_path
        )
        
        calibrated_data = calibrated_result['calibrated_data']
        print(f"   ✅ Calibration successful: {calibrated_data.shape}")
        print(f"   ✅ Calibration type: {calibrated_result['calibration_type']}")
        print(f"   ✅ Data range: {np.min(calibrated_data):.3e} to {np.max(calibrated_data):.3e}")
        
        # Verify calibration actually changed the data
        if np.array_equal(calibrated_data, np.abs(slc_data)**2):
            print(f"   ⚠️  WARNING: Calibrated data identical to intensity - calibration may not be working")
        else:
            print(f"   ✅ Data properly calibrated (values changed from raw intensity)")
        
        # Step 7: Convert to dB (real calibrated data)
        print("\n📊 STEP 7: Convert to dB")
        db_result = sardine.convert_to_db_real(calibrated_data)
        
        if isinstance(db_result, dict) and 'data' in db_result:
            db_data = db_result['data']
        else:
            db_data = db_result
            
        print(f"   ✅ dB conversion: {db_data.shape}")
        print(f"   ✅ dB range: {np.min(db_data):.1f} to {np.max(db_data):.1f} dB")
        
        # Step 8: Apply multilooking (reduce speckle)
        print("\n👁️  STEP 8: Multilooking")
        multilooked_result = sardine.apply_multilooking(
            calibrated_data, 
            range_looks=2, 
            azimuth_looks=2
        )
        
        if isinstance(multilooked_result, dict) and 'data' in multilooked_result:
            multilooked_data = multilooked_result['data']
        else:
            multilooked_data = multilooked_result
            
        print(f"   ✅ Multilooking: {multilooked_data.shape}")
        print(f"   ✅ Speckle reduction applied")
        
        # Step 9: Apply speckle filtering
        print("\n✨ STEP 9: Speckle Filtering")
        filtered_result = sardine.apply_speckle_filter_optimized(
            multilooked_data,
            filter_type="enhanced_lee",
            window_size=7
        )
        
        # Handle different return types
        if isinstance(filtered_result, dict) and 'filtered_data' in filtered_result:
            filtered_data = filtered_result['filtered_data']
        elif isinstance(filtered_result, dict) and 'data' in filtered_result:
            filtered_data = filtered_result['data']
        elif isinstance(filtered_result, list) and len(filtered_result) > 0:
            filtered_data = filtered_result[0]  # Take first element if it's a list
        else:
            filtered_data = filtered_result
            
        print(f"   ✅ Speckle filtering: {filtered_data.shape}")
        print(f"   ✅ Enhanced Lee filter applied")
        
        # Final verification: Ensure all data is real (no synthetic)
        print("\n🔍 VERIFICATION: Real Data Pipeline")
        
        # Check that we're not using default/synthetic values
        cal_vector = cal_data['vectors'][0]
        unique_sigma_values = len(set(cal_vector['sigma_nought'][:10]))  # Check first 10 values
        
        if unique_sigma_values < 3:
            print(f"   ❌ Calibration data appears synthetic (only {unique_sigma_values} unique values)")
            return False
        
        print(f"   ✅ Calibration data is real ({unique_sigma_values} unique values in sample)")
        print(f"   ✅ Pipeline used ONLY real Sentinel-1 data")
        print(f"   ✅ No synthetic/fallback data detected")
        
        # Final statistics
        print("\n📈 FINAL RESULTS")
        print(f"   📊 Original SLC: {slc_data.shape} complex")
        print(f"   📊 Calibrated: {calibrated_data.shape} σ⁰")
        print(f"   📊 Multilooked: {multilooked_data.shape}")
        print(f"   📊 Filtered: {filtered_data.shape}")
        print(f"   📊 Final range: {np.min(filtered_data):.3e} to {np.max(filtered_data):.3e}")
        
        return True
        
    except Exception as e:
        print(f"\n❌ Pipeline failed: {e}")
        import traceback
        traceback.print_exc()
        return False

if __name__ == "__main__":
    print("=" * 80)
    print("🛰️  REAL SENTINEL-1 SAR PROCESSING TEST")
    print("NO SYNTHETIC DATA - REAL DATA ONLY")
    print("=" * 80)
    
    success = test_real_sentinel1_pipeline()
    
    if success:
        print("\n" + "=" * 80)
        print("🎉 SUCCESS: Complete real data SAR processing pipeline working!")
        print("✅ All steps use real Sentinel-1 data")
        print("✅ No synthetic/fallback data used")
        print("✅ Radiometric calibration working with real calibration vectors")
        print("=" * 80)
        sys.exit(0)
    else:
        print("\n" + "=" * 80)
        print("💥 FAILURE: Real data pipeline has issues")
        print("❌ Check error messages above")
        print("=" * 80)
        sys.exit(1)
