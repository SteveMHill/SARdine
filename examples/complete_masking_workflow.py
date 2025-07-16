#!/usr/bin/env python3
"""
Complete Masking Workflow Example for SARdine

This example demonstrates how to apply a masking workflow to terrain-corrected
SAR gamma0 data using local incidence angle, DEM validity, and gamma0 value thresholds.

The masking workflow:
1. Computes local incidence angle cosine from DEM surface normals
2. Identifies invalid pixels based on multiple criteria:
   - Local incidence angle threshold (steep slopes)
   - DEM validity (below sea level, NaN values)
   - Gamma0 value range (outliers, invalid backscatter)
3. Combines all masks into a single validity mask
4. Applies the mask to gamma0 data with configurable fill values
5. Optionally saves individual mask components and LIA values

Prerequisites:
- Terrain-corrected gamma0 GeoTIFF file
- Corresponding DEM file with same projection and extent
- Optional: rasterio for advanced GeoTIFF handling
"""

import sys
import numpy as np
import sardine

def main():
    # Input files - update these paths
    gamma0_file = "test_terrain_correction/gamma0_terrain_corrected.tif"
    dem_file = "test_terrain_correction/dem_resampled.tif"
    
    # Output files
    masked_gamma0_file = "test_terrain_correction/gamma0_masked.tif"
    
    print("SARdine Masking Workflow Example")
    print("=" * 50)
    
    try:
        # For this example, we'll create synthetic data
        print("Creating synthetic test data...")
        
        # Generate synthetic gamma0 data (terrain-corrected backscatter in dB)
        height, width = 100, 100
        np.random.seed(42)
        
        # Realistic gamma0 values: -25 to -5 dB for most terrain
        gamma0_data = np.random.normal(-15.0, 5.0, (height, width)).astype(np.float32)
        
        # Add some outliers and invalid values
        gamma0_data[10:20, 10:20] = -60.0  # Very low values (water/shadow)
        gamma0_data[80:90, 80:90] = 15.0   # Very high values (urban/bright targets)
        gamma0_data[45:55, 45:55] = np.nan # Invalid values
        
        # Generate synthetic DEM data (elevation in meters)
        x = np.linspace(-2, 2, width)
        y = np.linspace(-2, 2, height)
        X, Y = np.meshgrid(x, y)
        
        # Create terrain with varying slopes
        dem_data = (100 * np.sin(X) * np.cos(Y) + 
                   50 * X**2 + 30 * Y**2 + 200).astype(np.float32)
        
        # Add some water areas (below sea level)
        dem_data[20:30, 70:80] = -50.0  # Below sea level
        dem_data[60:70, 20:30] = np.nan # Invalid DEM values
        
        # Create a simple terrain corrector for the masking workflow
        # For this demo, we'll use a placeholder file and simplified workflow
        print(f"Creating simplified masking workflow demo...")
        
        # Since the actual terrain corrector requires a real DEM file,
        # we'll demonstrate the workflow with simulated results
        corrector_info = "demo_terrain_corrector"
        
        # Create masking workflow with different thresholds
        print("\nConfiguring masking workflow...")
        
        # Conservative masking (keep more pixels)
        conservative_workflow = sardine.create_masking_workflow(
            lia_threshold=0.05,    # cos(87°) - very steep slopes only
            dem_threshold=-200.0,  # Very low elevation threshold
            gamma0_min=-45.0,      # Wide gamma0 range
            gamma0_max=20.0
        )
        
        # Standard masking (recommended settings)
        standard_workflow = sardine.create_masking_workflow(
            lia_threshold=0.1,     # cos(84°) - moderate slope threshold
            dem_threshold=-100.0,  # Below sea level threshold
            gamma0_min=-35.0,      # Typical gamma0 range
            gamma0_max=5.0
        )
        
        # Strict masking (remove more pixels)
        strict_workflow = sardine.create_masking_workflow(
            lia_threshold=0.2,     # cos(78°) - gentler slopes only
            dem_threshold=0.0,     # Above sea level only
            gamma0_min=-30.0,      # Narrow gamma0 range
            gamma0_max=0.0
        )
        
        # Apply different masking workflows
        workflows = {
            "conservative": conservative_workflow,
            "standard": standard_workflow,
            "strict": strict_workflow
        }
        
        results = {}
        
        for name, workflow in workflows.items():
            print(f"\nApplying {name} masking workflow...")
            
            # Apply masking workflow
            mask_result = sardine.apply_masking_workflow(
                corrector_info, gamma0_data, dem_data, workflow
            )
            
            results[name] = mask_result
            
            print(f"  Total pixels: {mask_result.total_pixels}")
            print(f"  Valid pixels: {mask_result.valid_pixels}")
            print(f"  Coverage: {mask_result.coverage_percent:.1f}%")
            
            # Get individual mask statistics
            combined_mask = mask_result.get_combined_mask()
            gamma0_mask = mask_result.get_gamma0_mask()
            dem_mask = mask_result.get_dem_mask()
            lia_mask = mask_result.get_lia_mask()
            
            # Convert to numpy arrays for analysis (our simplified implementation returns lists)
            combined_mask = np.array(combined_mask)
            gamma0_mask = np.array(gamma0_mask)
            dem_mask = np.array(dem_mask)
            lia_mask = np.array(lia_mask)
            
            print(f"  Gamma0 valid: {np.sum(gamma0_mask)}/{gamma0_mask.size} "
                  f"({np.sum(gamma0_mask)/gamma0_mask.size*100:.1f}%)")
            print(f"  DEM valid: {np.sum(dem_mask)}/{dem_mask.size} "
                  f"({np.sum(dem_mask)/dem_mask.size*100:.1f}%)")
            print(f"  LIA valid: {np.sum(lia_mask)}/{lia_mask.size} "
                  f"({np.sum(lia_mask)/lia_mask.size*100:.1f}%)")
        
        # Apply standard masking to create output
        print("\nApplying standard masking to gamma0 data...")
        
        standard_result = results["standard"]
        combined_mask = np.array(standard_result.get_combined_mask())
        
        # Apply mask with NaN fill value
        masked_gamma0_nan = sardine.apply_mask_to_gamma0(
            gamma0_data.tolist(), combined_mask.tolist(), float('nan')
        )
        
        # Apply mask with custom fill value
        masked_gamma0_fill = sardine.apply_mask_to_gamma0(
            gamma0_data.tolist(), combined_mask.tolist(), -999.0
        )
        
        print(f"Masked data shape: {gamma0_data.shape}")
        print(f"Valid pixels in combined mask: {np.sum(combined_mask)}")
        print(f"Invalid pixels: {np.sum(~combined_mask)}")
        
        # Analyze local incidence angle
        lia_cosine = np.array(standard_result.get_lia_cosine())
        valid_lia = lia_cosine[~np.isnan(lia_cosine)]
        
        if len(valid_lia) > 0:
            print(f"\nLocal Incidence Angle Analysis:")
            print(f"  LIA cosine range: [{np.min(valid_lia):.3f}, {np.max(valid_lia):.3f}]")
            print(f"  LIA angle range: [{np.degrees(np.arccos(np.max(valid_lia))):.1f}°, "
                  f"{np.degrees(np.arccos(np.min(valid_lia))):.1f}°]")
            print(f"  Mean LIA cosine: {np.mean(valid_lia):.3f}")
        
        # Demo: Quality assessment
        print("\nQuality Assessment:")
        
        # Since our simplified implementation doesn't apply actual masking,
        # we'll demonstrate with synthetic statistics
        print(f"  Simulated gamma0 statistics:")
        print(f"    Mean: {np.mean(gamma0_data):.2f} dB")
        print(f"    Std:  {np.std(gamma0_data):.2f} dB")
        print(f"    Range: [{np.min(gamma0_data):.2f}, {np.max(gamma0_data):.2f}] dB")
        
        # Compare coverage between workflows
        print(f"\nCoverage Comparison:")
        for name, result in results.items():
            print(f"  {name.capitalize():12}: {result.coverage_percent:5.1f}%")
        
        # Demo: Save results (if rasterio is available)
        try:
            import rasterio
            from rasterio.transform import Affine
            
            print(f"\nSaving masked results...")
            
            # Create geotransform for output
            geotransform = [-1.0, 0.01, 0.0, 1.0, 0.0, -0.01]  # 0.01 degree pixels
            transform = Affine.from_gdal(*geotransform)
            
            # Save masked gamma0
            profile = {
                'driver': 'GTiff',
                'width': width,
                'height': height,
                'count': 1,
                'dtype': 'float32',
                'crs': 'EPSG:4326',
                'transform': transform,
                'nodata': np.nan
            }
            
            # Save with different masking levels
            for name, result in results.items():
                output_file = f"test_terrain_correction/gamma0_masked_{name}.tif"
                mask = np.array(result.get_combined_mask())
                
                # Simulate masked data by applying NaN where mask is False
                masked_data = gamma0_data.copy()
                masked_data[~mask] = np.nan
                
                with rasterio.open(output_file, 'w', **profile) as dst:
                    dst.write(masked_data, 1)
                print(f"  Saved {name} masked gamma0: {output_file}")
            
            # Save mask components for standard workflow
            mask_profile = profile.copy()
            mask_profile.update({'dtype': 'uint8', 'nodata': 255})
            
            masks = {
                'combined': np.array(standard_result.get_combined_mask()),
                'gamma0': np.array(standard_result.get_gamma0_mask()),
                'dem': np.array(standard_result.get_dem_mask()),
                'lia': np.array(standard_result.get_lia_mask())
            }
            
            for mask_name, mask_data in masks.items():
                mask_file = f"test_terrain_correction/mask_{mask_name}.tif"
                with rasterio.open(mask_file, 'w', **mask_profile) as dst:
                    dst.write((mask_data * 255).astype(np.uint8), 1)
                print(f"  Saved {mask_name} mask: {mask_file}")
            
            # Save LIA cosine values
            lia_file = "test_terrain_correction/lia_cosine.tif"
            with rasterio.open(lia_file, 'w', **profile) as dst:
                dst.write(np.array(standard_result.get_lia_cosine()), 1)
            print(f"  Saved LIA cosine: {lia_file}")
            
        except ImportError:
            print("  rasterio not available - skipping file output")
            print("  Install rasterio to save GeoTIFF results")
        
        print("\nMasking workflow completed successfully!")
        print("\nNext steps:")
        print("1. Integrate masking into your processing pipeline")
        print("2. Adjust thresholds based on your data and requirements")
        print("3. Consider local incidence angle for radiometric correction")
        print("4. Use masks for quality assessment and filtering")
        
    except Exception as e:
        print(f"Error in masking workflow: {e}")
        import traceback
        traceback.print_exc()
        return 1
    
    return 0

if __name__ == "__main__":
    sys.exit(main())
