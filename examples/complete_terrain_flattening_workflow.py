#!/usr/bin/env python3
"""
Complete Terrain Flattening Workflow Example

This example demonstrates the automatic DEM-based terrain flattening (gamma0) 
processing pipeline in SARdine. It shows how to:

1. Process Sentinel-1 SLC data with automatic DEM download/preparation
2. Apply radiometric calibration, multilooking, and terrain flattening in one step
3. Handle orbit data for improved geometric accuracy
4. Export results for further analysis

Features demonstrated:
- Automatic SRTM DEM download and preparation
- Per-pixel local incidence angle computation
- Cosine division terrain flattening method
- Configurable multilooking parameters
- Error handling and logging
"""

import sys
import numpy as np
from pathlib import Path

# Add the Python module to the path
sys.path.insert(0, str(Path(__file__).parent.parent / "python"))

import sardine


def terrain_flattening_workflow():
    """Complete terrain flattening workflow example."""
    
    # Input data paths
    slc_file = "../data/S1A_IW_SLC__1SDV_20200103T170815_20200103T170842_030639_0382D5_DADE.zip"
    orbit_file = "../orbit_test/S1A_IW_SLC__1SDV_20200103T170815_20200103T170842_030639_0382D5_DADE_RESORB.EOF"
    
    print("🏔️  SARdine Terrain Flattening Workflow")
    print("=" * 50)
    
    # Check if input files exist
    if not Path(slc_file).exists():
        print(f"❌ SLC file not found: {slc_file}")
        print("   Please ensure the test data is available")
        return False
    
    try:
        # Step 1: Open SLC reader
        print(f"📁 Opening SLC file: {Path(slc_file).name}")
        reader = sardine.SlcReader(slc_file)
        
        # Step 2: Load orbit data for improved accuracy
        if Path(orbit_file).exists():
            print(f"📡 Loading orbit data: {Path(orbit_file).name}")
            orbit_data = sardine.load_orbit_file(orbit_file)
            reader.set_orbit_data(orbit_data)
            print("✅ Orbit data loaded successfully")
        else:
            print(f"⚠️  Orbit file not found: {orbit_file}")
            print("   Processing will continue with reduced accuracy")
        
        # Step 3: Configure processing parameters
        polarization = "VV"
        calibration_type = "Sigma0"  # Input calibration before terrain flattening
        range_looks = 4
        azimuth_looks = 1
        dem_cache_dir = "./dem_cache"
        
        print(f"\n🔧 Processing Configuration:")
        print(f"   Polarization: {polarization}")
        print(f"   Calibration: {calibration_type}")
        print(f"   Multilooking: {range_looks} x {azimuth_looks}")
        print(f"   DEM cache: {dem_cache_dir}")
        
        # Step 4: Process with automatic DEM preparation and terrain flattening
        print(f"\n🚀 Starting terrain flattening workflow...")
        print("   This includes:")
        print("   • Automatic DEM download/preparation")
        print("   • Radiometric calibration")
        print("   • Multilooking")
        print("   • DEM resampling to SAR geometry")
        print("   • Slope/aspect computation")
        print("   • Local incidence angle calculation")
        print("   • Terrain flattening (gamma0 = sigma0 / cos(θ_lia))")
        
        # This method handles the complete workflow automatically
        gamma0_data, incidence_angles, range_spacing, azimuth_spacing = reader.calibrate_multilook_and_flatten_auto_dem(
            polarization,
            calibration_type,
            range_looks,
            azimuth_looks,
            dem_cache_dir
        )
        
        print(f"✅ Terrain flattening completed!")
        
        # Step 5: Analyze results
        print(f"\n📊 Results Summary:")
        print(f"   Gamma0 shape: {gamma0_data.shape}")
        print(f"   Pixel spacing: {range_spacing:.2f}m x {azimuth_spacing:.2f}m")
        print(f"   Gamma0 range: {float(gamma0_data.min()):.6f} - {float(gamma0_data.max()):.6f}")
        print(f"   Incidence angle range: {float(incidence_angles.min()):.2f}° - {float(incidence_angles.max()):.2f}°")
        
        # Calculate statistics
        valid_gamma0 = gamma0_data[~np.isnan(gamma0_data)]
        valid_angles = incidence_angles[~np.isnan(incidence_angles)]
        
        print(f"\n📈 Statistics (valid pixels only):")
        print(f"   Valid pixels: {len(valid_gamma0):,} / {gamma0_data.size:,} ({100*len(valid_gamma0)/gamma0_data.size:.1f}%)")
        print(f"   Gamma0 mean: {valid_gamma0.mean():.6f}")
        print(f"   Gamma0 std:  {valid_gamma0.std():.6f}")
        print(f"   Angle mean:  {valid_angles.mean():.2f}°")
        print(f"   Angle std:   {valid_angles.std():.2f}°")
        
        # Step 6: Save results
        output_dir = Path("terrain_flattening_results")
        output_dir.mkdir(exist_ok=True)
        
        # Save gamma0 data
        gamma0_file = output_dir / f"gamma0_{polarization}_{range_looks}x{azimuth_looks}.npy"
        gamma0_data.save(str(gamma0_file))
        print(f"💾 Gamma0 saved: {gamma0_file}")
        
        # Save incidence angles
        angles_file = output_dir / f"incidence_angles_{polarization}_{range_looks}x{azimuth_looks}.npy"
        incidence_angles.save(str(angles_file))
        print(f"💾 Incidence angles saved: {angles_file}")
        
        # Save in dB scale
        gamma0_db = 10 * np.log10(np.maximum(valid_gamma0, 1e-10))
        gamma0_db_full = np.full_like(gamma0_data, np.nan)
        gamma0_db_full[~np.isnan(gamma0_data)] = gamma0_db
        
        gamma0_db_file = output_dir / f"gamma0_db_{polarization}_{range_looks}x{azimuth_looks}.npy"
        np.save(str(gamma0_db_file), gamma0_db_full)
        print(f"💾 Gamma0 (dB) saved: {gamma0_db_file}")
        
        # Step 7: Generate processing summary
        summary_file = output_dir / "processing_summary.txt"
        with open(summary_file, 'w') as f:
            f.write("SARdine Terrain Flattening Results\n")
            f.write("="*40 + "\n\n")
            f.write(f"Input SLC: {Path(slc_file).name}\n")
            f.write(f"Orbit file: {Path(orbit_file).name if Path(orbit_file).exists() else 'Not available'}\n")
            f.write(f"Polarization: {polarization}\n")
            f.write(f"Multilooking: {range_looks} x {azimuth_looks}\n")
            f.write(f"Output shape: {gamma0_data.shape}\n")
            f.write(f"Pixel spacing: {range_spacing:.2f}m x {azimuth_spacing:.2f}m\n")
            f.write(f"Valid pixels: {len(valid_gamma0):,} / {gamma0_data.size:,} ({100*len(valid_gamma0)/gamma0_data.size:.1f}%)\n")
            f.write(f"Gamma0 range: {float(gamma0_data.min()):.6f} - {float(gamma0_data.max()):.6f}\n")
            f.write(f"Gamma0 mean ± std: {valid_gamma0.mean():.6f} ± {valid_gamma0.std():.6f}\n")
            f.write(f"Incidence angle range: {float(incidence_angles.min()):.2f}° - {float(incidence_angles.max()):.2f}°\n")
            f.write(f"Incidence angle mean ± std: {valid_angles.mean():.2f}° ± {valid_angles.std():.2f}°\n")
        
        print(f"📄 Processing summary saved: {summary_file}")
        
        print(f"\n✅ Terrain flattening workflow completed successfully!")
        print(f"📁 Results saved in: {output_dir}")
        print(f"\n🔬 Next steps:")
        print(f"   • Convert to GeoTIFF for GIS analysis")
        print(f"   • Apply additional filtering or classification")
        print(f"   • Compare with other terrain correction methods")
        print(f"   • Validate results with ground truth data")
        
        return True
        
    except Exception as e:
        print(f"❌ Error during terrain flattening: {e}")
        return False


def main():
    """Main function."""
    print("SARdine Terrain Flattening Example")
    print("This example demonstrates automatic DEM-based terrain flattening")
    print("using the cosine division method with per-pixel local incidence angles.\n")
    
    success = terrain_flattening_workflow()
    
    if success:
        print("\n🎉 Example completed successfully!")
        return 0
    else:
        print("\n❌ Example failed!")
        return 1


if __name__ == "__main__":
    sys.exit(main())
