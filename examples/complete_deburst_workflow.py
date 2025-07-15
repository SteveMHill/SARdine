#!/usr/bin/env python3
"""
Complete deburst workflow example for SARdine.

This script demonstrates the deburst step in the SAR processing pipeline:
1. Load Sentinel-1 IW SLC data
2. Extract burst information from annotation XML
3. Deburst each polarization
4. Export results
"""

import sys
import time
from pathlib import Path


def run_complete_deburst_workflow():
    """Run a complete deburst workflow."""
    
    slc_path = "/home/datacube/SARdine/data/S1A_IW_SLC__1SDV_20200103T170815_20200103T170842_030639_0382D5_DADE.zip"
    
    if not Path(slc_path).exists():
        print(f"❌ Sample SLC file not found: {slc_path}")
        return False
    
    print("🚀 SARdine Complete Deburst Workflow")
    print("=" * 60)
    print(f"📁 Input: {Path(slc_path).name}")
    
    try:
        # Import sardine
        import sardine
        
        # Initialize reader
        reader = sardine.SlcReader(slc_path)
        
        # Verify this is IW mode
        is_iw = reader.is_iw_mode()
        print(f"📡 IW mode: {'✅' if is_iw else '❌'}")
        
        if not is_iw:
            print("⚠️  Warning: This workflow is optimized for IW mode data")
        
        # Get available polarizations
        annotation_files = reader.find_annotation_files()
        available_pols = list(annotation_files.keys())
        print(f"📡 Available polarizations: {', '.join(available_pols)}")
        
        # Step 1: IW Split Analysis (prerequisite for deburst)
        print(f"\n📋 Step 1: IW Split Analysis")
        print("-" * 40)
        
        all_subswaths = reader.get_all_iw_subswaths()
        
        for pol, subswaths in all_subswaths.items():
            print(f"\n📡 {pol} Polarization:")
            for swath_id, swath in subswaths.items():
                print(f"  🎯 {swath_id}: {swath.burst_count} bursts, "
                      f"{swath.azimuth_samples:,} x {swath.range_samples:,}")
        
        # Step 2: Deburst Processing
        print(f"\n🔄 Step 2: Deburst Processing")
        print("-" * 40)
        
        deburst_results = {}
        total_start = time.time()
        
        for pol in available_pols:
            print(f"\n🔄 Processing {pol} polarization...")
            pol_start = time.time()
            
            try:
                # Perform deburst
                deburst_data, (rows, cols) = reader.deburst_slc(pol)
                pol_end = time.time()
                
                processing_time = pol_end - pol_start
                data_size_gb = (rows * cols * 8) / (1024**3)  # Complex64 = 8 bytes
                
                deburst_results[pol] = {
                    'dimensions': (rows, cols),
                    'processing_time': processing_time,
                    'data_size_gb': data_size_gb
                }
                
                print(f"✅ {pol} deburst completed:")
                print(f"   • Output dimensions: {rows:,} x {cols:,}")
                print(f"   • Data size: {data_size_gb:.2f} GB")
                print(f"   • Processing time: {processing_time:.1f} seconds")
                print(f"   • Processing rate: {data_size_gb/processing_time:.2f} GB/s")
                
                # Quick data validation
                if len(deburst_data) > 0 and len(deburst_data[0]) > 0:
                    sample_pixel = deburst_data[rows//2][cols//2]
                    amplitude = (sample_pixel[0]**2 + sample_pixel[1]**2)**0.5
                    print(f"   • Sample amplitude: {amplitude:.3f}")
                
            except Exception as e:
                print(f"❌ Error processing {pol}: {e}")
                return False
        
        total_end = time.time()
        
        # Step 3: Results Summary
        print(f"\n📊 Step 3: Results Summary")
        print("-" * 40)
        
        total_time = total_end - total_start
        total_data_gb = sum(r['data_size_gb'] for r in deburst_results.values())
        
        print(f"🎯 Deburst processing completed!")
        print(f"   • Total polarizations: {len(deburst_results)}")
        print(f"   • Total data processed: {total_data_gb:.2f} GB")
        print(f"   • Total processing time: {total_time:.1f} seconds")
        print(f"   • Average processing rate: {total_data_gb/total_time:.2f} GB/s")
        
        # Individual polarization summary
        for pol, result in deburst_results.items():
            dims = result['dimensions']
            print(f"   • {pol}: {dims[0]:,} x {dims[1]:,} "
                  f"({result['data_size_gb']:.2f} GB in {result['processing_time']:.1f}s)")
        
        # Step 4: Next Steps
        print(f"\n🚀 Step 4: Next Steps in SAR Pipeline")
        print("-" * 40)
        print("✅ Completed: SLC Read → Orbit Application → IW Split → Deburst")
        print("🔄 Next steps:")
        print("   1. Radiometric Calibration (SLC → Sigma0/Beta0)")
        print("   2. Multilooking (reduce speckle)")
        print("   3. Terrain Flattening (Sigma0 → Gamma0)")
        print("   4. Speckle Filtering")
        print("   5. Terrain Correction (geocoding)")
        print("   6. Output Generation (GeoTIFF)")
        
        print(f"\n🎉 Deburst workflow completed successfully!")
        return True
        
    except ImportError as e:
        print(f"❌ Import error: {e}")
        print("Make sure SARdine is properly installed: pip install -e .")
        return False
    except Exception as e:
        print(f"❌ Error during workflow: {e}")
        import traceback
        traceback.print_exc()
        return False


if __name__ == "__main__":
    success = run_complete_deburst_workflow()
    sys.exit(0 if success else 1)
