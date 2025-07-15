#!/usr/bin/env python3
"""
Complete IW Split Example for SARdine

This example demonstrates the full IW split workflow:
1. Load Sentinel-1 SLC product
2. Check IW mode and analyze sub-swaths
3. Extract detailed sub-swath information
4. Show practical usage for next pipeline steps

The IW split step is crucial for SAR processing as it:
- Identifies individual sub-swaths (IW1, IW2, IW3) within the SLC
- Extracts geometric and sampling parameters for each sub-swath
- Provides burst information needed for debursting
- Prepares data for independent processing of each sub-swath
"""

import sys
import json
from pathlib import Path

# Add the python directory to sys.path to import sardine
sys.path.insert(0, str(Path(__file__).parent.parent / "python"))

import sardine

def complete_iw_split_workflow():
    """Complete workflow demonstrating IW split functionality."""
    
    # Path to the sample SLC file
    slc_path = "/home/datacube/SARdine/data/S1A_IW_SLC__1SDV_20200103T170815_20200103T170842_030639_0382D5_DADE.zip"
    
    if not Path(slc_path).exists():
        print(f"❌ Sample SLC file not found: {slc_path}")
        return False
    
    print("🚀 Complete IW Split Workflow")
    print("=" * 60)
    
    try:
        # Step 1: Initialize SLC Reader
        print(f"📂 Step 1: Loading SLC file")
        print(f"   File: {Path(slc_path).name}")
        reader = sardine.SlcReader(slc_path)
        
        # Step 2: Verify IW mode
        print(f"\n🔍 Step 2: Verifying acquisition mode")
        is_iw = reader.is_iw_mode()
        print(f"   IW Mode: {'✅ Confirmed' if is_iw else '❌ Not IW mode'}")
        
        if not is_iw:
            print("⚠️  This workflow requires an IW mode product")
            return False
        
        # Step 3: Analyze available polarizations
        print(f"\n📊 Step 3: Analyzing available polarizations")
        annotations = reader.find_annotation_files()
        polarizations = list(annotations.keys())
        print(f"   Available polarizations: {polarizations}")
        print(f"   Annotation files:")
        for pol, file_path in annotations.items():
            filename = Path(file_path).name
            print(f"     {pol}: {filename}")
        
        # Step 4: Extract sub-swaths for each polarization
        print(f"\n🎯 Step 4: Extracting IW sub-swaths")
        all_subswaths = reader.get_all_iw_subswaths()
        
        if not all_subswaths:
            print("❌ No sub-swaths extracted")
            return False
        
        total_subswaths = sum(len(sw) for sw in all_subswaths.values())
        print(f"   Successfully extracted {total_subswaths} sub-swaths")
        
        # Step 5: Detailed analysis of each sub-swath
        print(f"\n📈 Step 5: Detailed sub-swath analysis")
        
        for pol_str, subswaths in all_subswaths.items():
            print(f"\n   📡 Polarization: {pol_str}")
            print(f"   {'─' * 40}")
            
            for swath_id, swath in subswaths.items():
                print(f"     🎯 Sub-swath: {swath_id}")
                print(f"        • Identifier: {swath.id}")
                print(f"        • Burst count: {swath.burst_count}")
                print(f"        • Dimensions: {swath.range_samples} × {swath.azimuth_samples} pixels")
                print(f"        • Range pixel spacing: {swath.range_pixel_spacing:.3f} m")
                print(f"        • Azimuth pixel spacing: {swath.azimuth_pixel_spacing:.3f} m")
                print(f"        • Slant range time: {swath.slant_range_time:.6f} s")
                print(f"        • Burst duration: {swath.burst_duration:.6f} s")
                
                # Calculate derived parameters
                range_extent = swath.range_samples * swath.range_pixel_spacing / 1000  # km
                azimuth_extent = swath.azimuth_samples * swath.azimuth_pixel_spacing / 1000  # km
                total_pixels = swath.range_samples * swath.azimuth_samples
                data_size_mb = total_pixels * 8 / (1024 * 1024)  # Assuming complex data (8 bytes)
                
                print(f"        • Coverage: {range_extent:.2f} × {azimuth_extent:.2f} km")
                print(f"        • Data size: ~{data_size_mb:.1f} MB (complex)")
        
        # Step 6: Demonstrate individual sub-swath extraction
        print(f"\n🔬 Step 6: Individual sub-swath extraction test")
        for pol_str in all_subswaths.keys():
            individual_swaths = reader.extract_iw_subswaths(pol_str)
            print(f"   {pol_str}: ✅ Extracted {len(individual_swaths)} sub-swath(s)")
        
        # Step 7: Generate processing recommendations
        print(f"\n💡 Step 7: Processing recommendations")
        
        # Find the sub-swath with most bursts (typically central sub-swath)
        max_bursts = 0
        best_swath = None
        best_pol = None
        
        for pol_str, subswaths in all_subswaths.items():
            for swath_id, swath in subswaths.items():
                if swath.burst_count > max_bursts:
                    max_bursts = swath.burst_count
                    best_swath = swath_id
                    best_pol = pol_str
        
        print(f"   📍 Recommended primary sub-swath: {best_swath} ({best_pol})")
        print(f"      Reason: Highest burst count ({max_bursts} bursts)")
        
        # Calculate total data volume
        total_pixels = sum(
            sum(swath.range_samples * swath.azimuth_samples for swath in subswaths.values())
            for subswaths in all_subswaths.values()
        )
        total_size_gb = total_pixels * 8 / (1024**3)  # GB
        print(f"   💾 Total data volume: ~{total_size_gb:.2f} GB (all sub-swaths)")
        
        # Step 8: Export results for further processing
        print(f"\n💾 Step 8: Exporting results")
        
        # Create comprehensive results dictionary
        results = {
            "product_info": {
                "filename": Path(slc_path).name,
                "mode": "IW",
                "polarizations": polarizations,
                "total_subswaths": total_subswaths
            },
            "subswaths": {}
        }
        
        for pol_str, subswaths in all_subswaths.items():
            results["subswaths"][pol_str] = {}
            for swath_id, swath in subswaths.items():
                results["subswaths"][pol_str][swath_id] = {
                    "id": swath.id,
                    "burst_count": swath.burst_count,
                    "range_samples": swath.range_samples,
                    "azimuth_samples": swath.azimuth_samples,
                    "range_pixel_spacing": swath.range_pixel_spacing,
                    "azimuth_pixel_spacing": swath.azimuth_pixel_spacing,
                    "slant_range_time": swath.slant_range_time,
                    "burst_duration": swath.burst_duration,
                    "range_extent_km": swath.range_samples * swath.range_pixel_spacing / 1000,
                    "azimuth_extent_km": swath.azimuth_samples * swath.azimuth_pixel_spacing / 1000,
                    "total_pixels": swath.range_samples * swath.azimuth_samples
                }
        
        # Save to file
        output_file = "/tmp/complete_iw_split_results.json"
        with open(output_file, 'w') as f:
            json.dump(results, f, indent=2)
        
        print(f"   ✅ Results saved to: {output_file}")
        
        # Step 9: Next steps guidance
        print(f"\n🎯 Step 9: Next steps in SAR processing pipeline")
        print(f"   After IW split, the next pipeline steps are:")
        print(f"   1. 🎯 Burst extraction - Extract individual bursts from each sub-swath")
        print(f"   2. 🔄 Debursting - Remove burst boundaries and create continuous images")
        print(f"   3. 📏 Radiometric calibration - Convert to calibrated backscatter")
        print(f"   4. 🔍 Multilooking - Reduce speckle and resolution")
        print(f"   5. 🌍 Terrain correction - Geocode to map projection")
        
        print(f"\n✅ IW Split workflow completed successfully!")
        
        return True
        
    except Exception as e:
        print(f"❌ Error during IW split workflow: {e}")
        import traceback
        traceback.print_exc()
        return False


def main():
    """Main function."""
    print("🌟 SARdine: Complete IW Split Workflow Example")
    print("=" * 60)
    
    success = complete_iw_split_workflow()
    
    if success:
        print("\n🎉 Workflow completed successfully!")
        print("\n📚 Key takeaways:")
        print("   • IW split successfully extracts sub-swath information")
        print("   • Each sub-swath has unique geometry and sampling parameters")
        print("   • Results can be exported for use in subsequent processing steps")
        print("   • The next step is burst extraction and debursting")
        return 0
    else:
        print("\n💥 Workflow failed!")
        return 1


if __name__ == "__main__":
    sys.exit(main())
