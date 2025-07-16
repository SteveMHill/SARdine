#!/usr/bin/env python3
"""
Complete TOPSAR Merge Workflow Example

This example demonstrates the complete processing chain for Sentinel-1 IW data:
1. Calibration (sigma0, beta0, or gamma0)
2. TOPSAR merge (combine IW1, IW2, IW3 sub-swaths)
3. Multilooking (speckle reduction)
4. Optional: Terrain correction (geocoding)

The TOPSAR merge step is critical for IW data as it:
- Combines the three IW sub-swaths (IW1, IW2, IW3) into a single wide-swath image
- Handles overlap regions between adjacent sub-swaths with proper blending
- Should be performed after calibration but before multilooking
- Preserves radiometric accuracy across the full swath width

Usage:
    python complete_topsar_merge_workflow.py <input_slc.zip> [output_directory]
"""

import sys
import time
import numpy as np
from pathlib import Path
import sardine


def print_processing_step(step_name: str, description: str):
    """Print a formatted processing step header."""
    print(f"\n{'='*80}")
    print(f"🔄 STEP: {step_name}")
    print(f"📋 {description}")
    print('='*80)


def print_data_summary(data: np.ndarray, name: str, units: str = ""):
    """Print summary statistics for data array."""
    print(f"\n📊 {name} Summary:")
    print(f"   • Shape: {data.shape}")
    print(f"   • Data type: {data.dtype}")
    print(f"   • Min: {np.min(data):.2e} {units}")
    print(f"   • Max: {np.max(data):.2e} {units}")
    print(f"   • Mean: {np.mean(data):.2e} {units}")
    print(f"   • Std: {np.std(data):.2e} {units}")
    print(f"   • Size: {data.size:,} pixels ({data.nbytes / 1024**2:.1f} MB)")


def complete_topsar_merge_workflow(input_file: str, output_dir: str = None):
    """
    Run the complete TOPSAR merge workflow.
    
    Args:
        input_file: Path to Sentinel-1 SLC ZIP file
        output_dir: Output directory (default: same as input)
    """
    
    print("🛰️  SARdine - Complete TOPSAR Merge Workflow")
    print("=" * 80)
    
    input_path = Path(input_file)
    if not input_path.exists():
        raise FileNotFoundError(f"Input file not found: {input_path}")
    
    # Set up output directory
    if output_dir is None:
        output_dir = input_path.parent / f"{input_path.stem}_topsar_processed"
    else:
        output_dir = Path(output_dir)
    
    output_dir.mkdir(exist_ok=True)
    print(f"📁 Output directory: {output_dir}")
    
    # Initialize SLC reader
    print(f"\n🔍 Initializing SLC reader for: {input_path.name}")
    reader = sardine.SlcReader(str(input_path))
    
    # Get product information
    annotation_files = reader.find_annotation_files()
    available_polarizations = list(annotation_files.keys())
    print(f"📡 Available polarizations: {', '.join(available_polarizations)}")
    
    # Process each polarization
    for pol in available_polarizations:
        print(f"\n" + "🎯" * 20 + f" PROCESSING {pol} " + "🎯" * 20)
        
        # Step 1: Get sub-swath information
        print_processing_step("1. SUB-SWATH ANALYSIS", 
                            "Analyzing IW sub-swaths and overlap regions")
        
        subswaths = reader.get_subswath_info(pol)
        print(f"📊 Sub-swaths found: {len(subswaths)}")
        
        for i, swath in enumerate(subswaths):
            print(f"   • {swath['swath_id']}: {swath['num_bursts']} bursts, "
                  f"{swath['samples_per_burst']} samples/burst")
        
        if len(subswaths) < 2:
            print("⚠️  Warning: Less than 2 sub-swaths found, skipping TOPSAR merge")
            continue
        
        # Step 2: Calibration
        print_processing_step("2. RADIOMETRIC CALIBRATION", 
                            "Converting DN values to sigma0 backscatter")
        
        start_time = time.time()
        
        # Perform calibration (sigma0 is most common for backscatter analysis)
        calibrated_data, (rows, cols) = reader.calibrate_slc(pol, "sigma0")
        calibrated_array = np.array(calibrated_data)
        
        calibration_time = time.time() - start_time
        print(f"✅ Calibration completed in {calibration_time:.2f} seconds")
        print_data_summary(calibrated_array, "Calibrated Data", "m²/m²")
        
        # Save calibrated data
        calib_output = output_dir / f"{input_path.stem}_{pol.lower()}_calibrated_sigma0.npy"
        np.save(calib_output, calibrated_array.astype(np.float32))
        print(f"💾 Calibrated data saved: {calib_output}")
        
        # Step 3: TOPSAR Merge
        print_processing_step("3. TOPSAR MERGE", 
                            "Merging IW sub-swaths with overlap blending")
        
        start_time = time.time()
        
        # Perform TOPSAR merge with feathering in overlap regions
        merged_result = sardine.topsar_merge(
            str(input_path),
            pol,
            "sigma0",  # Use same calibration as above
            "feather",  # Smooth blending in overlap regions
            "auto"     # Automatic output grid determination
        )
        
        merge_time = time.time() - start_time
        
        # Extract results
        merged_data = merged_result["intensity_data"]
        valid_mask = merged_result["valid_mask"]
        metadata = merged_result["metadata"]
        
        print(f"✅ TOPSAR merge completed in {merge_time:.2f} seconds")
        print(f"📊 Merge Statistics:")
        print(f"   • Sub-swaths merged: {metadata['num_swaths']}")
        print(f"   • Overlap regions: {metadata['overlap_count']}")
        print(f"   • Valid pixels: {np.sum(valid_mask):,} / {valid_mask.size:,} "
              f"({100 * np.sum(valid_mask) / valid_mask.size:.1f}%)")
        
        print_data_summary(merged_data, "Merged Data", "m²/m²")
        
        # Save merged data
        merged_output = output_dir / f"{input_path.stem}_{pol.lower()}_merged_sigma0.npy"
        mask_output = output_dir / f"{input_path.stem}_{pol.lower()}_merged_mask.npy"
        
        np.save(merged_output, merged_data.astype(np.float32))
        np.save(mask_output, valid_mask)
        
        print(f"💾 Merged data saved: {merged_output}")
        print(f"💾 Valid mask saved: {mask_output}")
        
        # Step 4: Multilooking (optional speckle reduction)
        print_processing_step("4. MULTILOOKING", 
                            "Reducing speckle noise through spatial averaging")
        
        start_time = time.time()
        
        # Apply multilooking (4x1 is common for Sentinel-1 IW)
        range_looks = 4
        azimuth_looks = 1
        
        # Simple multilooking implementation
        h, w = merged_data.shape
        ml_height = h // azimuth_looks
        ml_width = w // range_looks
        
        # Reshape and average
        reshaped = merged_data[:ml_height*azimuth_looks, :ml_width*range_looks]
        reshaped = reshaped.reshape(ml_height, azimuth_looks, ml_width, range_looks)
        multilooked = np.mean(reshaped, axis=(1, 3))
        
        # Also process the mask
        mask_reshaped = valid_mask[:ml_height*azimuth_looks, :ml_width*range_looks]
        mask_reshaped = mask_reshaped.reshape(ml_height, azimuth_looks, ml_width, range_looks)
        ml_mask = np.mean(mask_reshaped, axis=(1, 3)) > 0.5  # At least 50% valid
        
        multilook_time = time.time() - start_time
        
        print(f"✅ Multilooking completed in {multilook_time:.2f} seconds")
        print(f"📊 Multilook parameters: {range_looks}x{azimuth_looks} looks")
        print_data_summary(multilooked, "Multilooked Data", "m²/m²")
        
        # Save multilooked data
        ml_output = output_dir / f"{input_path.stem}_{pol.lower()}_multilooked_{range_looks}x{azimuth_looks}.npy"
        ml_mask_output = output_dir / f"{input_path.stem}_{pol.lower()}_multilooked_mask.npy"
        
        np.save(ml_output, multilooked.astype(np.float32))
        np.save(ml_mask_output, ml_mask)
        
        print(f"💾 Multilooked data saved: {ml_output}")
        print(f"💾 Multilooked mask saved: {ml_mask_output}")
        
        # Step 5: Convert to dB (optional)
        print_processing_step("5. dB CONVERSION", 
                            "Converting linear values to dB scale")
        
        # Convert to dB scale (standard for visualization and analysis)
        db_data = 10 * np.log10(np.maximum(multilooked, 1e-10))  # Avoid log(0)
        
        print_data_summary(db_data, "dB Data", "dB")
        
        # Save dB data
        db_output = output_dir / f"{input_path.stem}_{pol.lower()}_multilooked_db.npy"
        np.save(db_output, db_data.astype(np.float32))
        print(f"💾 dB data saved: {db_output}")
        
        # Summary for this polarization
        total_time = calibration_time + merge_time + multilook_time
        print(f"\n📈 Processing Summary for {pol}:")
        print(f"   • Calibration: {calibration_time:.2f}s")
        print(f"   • TOPSAR merge: {merge_time:.2f}s") 
        print(f"   • Multilooking: {multilook_time:.2f}s")
        print(f"   • Total: {total_time:.2f}s")
        
        # Data size summary
        original_size = calibrated_array.nbytes / 1024**2
        merged_size = merged_data.nbytes / 1024**2
        final_size = multilooked.nbytes / 1024**2
        
        print(f"   • Original size: {original_size:.1f} MB")
        print(f"   • Merged size: {merged_size:.1f} MB")
        print(f"   • Final size: {final_size:.1f} MB (compression: {final_size/original_size:.1%})")
    
    print(f"\n🎉 TOPSAR merge workflow completed successfully!")
    print(f"📁 All outputs saved to: {output_dir}")
    print(f"\n📋 Processing Order Summary:")
    print(f"   1. ✅ Radiometric Calibration (DN → σ⁰)")
    print(f"   2. ✅ TOPSAR Merge (IW1+IW2+IW3 → Full Swath)")
    print(f"   3. ✅ Multilooking (Speckle Reduction)")
    print(f"   4. ✅ dB Conversion (Linear → Logarithmic)")
    print(f"\n💡 Next steps:")
    print(f"   • Apply terrain correction (geocoding) for geographic projection")
    print(f"   • Use speckle filtering for further noise reduction")
    print(f"   • Perform change detection or classification analysis")


def main():
    """Main entry point."""
    if len(sys.argv) < 2:
        print("Usage: python complete_topsar_merge_workflow.py <input_slc.zip> [output_directory]")
        print("\nExample:")
        print("  python complete_topsar_merge_workflow.py data/S1A_IW_SLC__1SDV_20200103T170815_20200103T170842_030639_0382D5_DADE.zip")
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_dir = sys.argv[2] if len(sys.argv) > 2 else None
    
    try:
        complete_topsar_merge_workflow(input_file, output_dir)
    except Exception as e:
        print(f"❌ Error: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()
