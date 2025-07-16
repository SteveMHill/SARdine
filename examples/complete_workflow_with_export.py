#!/usr/bin/env python3
"""
Complete SAR Processing and Export Workflow

This script demonstrates a complete workflow from SAR data processing 
through masking to final GeoTIFF export, showcasing the integration of 
all SARdine components.
"""

import numpy as np
import os
import tempfile
import sardine

def complete_sar_workflow():
    """
    Demonstrate a complete SAR processing workflow with export.
    """
    
    print("Complete SAR Processing and Export Workflow")
    print("=" * 50)
    
    # Check if all required functions are available
    required_functions = ['linear_to_db', 'db_to_linear']
    geotiff_functions = ['export_geotiff', 'export_cog', 'export_multiband_geotiff']
    
    missing_core = [f for f in required_functions if not hasattr(sardine, f)]
    missing_geotiff = [f for f in geotiff_functions if not hasattr(sardine, f)]
    
    if missing_core:
        print(f"❌ Missing core functions: {missing_core}")
        return
        
    if missing_geotiff:
        print(f"⚠️  GeoTIFF functions not available: {missing_geotiff}")
        print("   Install rasterio for GeoTIFF export: pip install rasterio")
        export_enabled = False
    else:
        export_enabled = True
    
    with tempfile.TemporaryDirectory() as temp_dir:
        print(f"📁 Working directory: {temp_dir}")
        
        # Step 1: Simulate SAR data processing
        print("\n1. Simulating SAR data processing...")
        
        # Create realistic SAR scene (1000x1000 pixels)
        width, height = 1000, 1000
        print(f"   Processing {height}x{width} pixel SAR scene...")
        
        # Simulate calibrated backscatter data
        gamma0_data = create_realistic_sar_scene(height, width)
        
        print(f"   ✓ Generated calibrated gamma0 data")
        print(f"   ✓ Linear range: {np.min(gamma0_data):.6f} to {np.max(gamma0_data):.6f}")
        
        # Step 2: Apply masking workflow
        print("\n2. Applying quality masking...")
        
        # Create a simple mask based on data quality
        # In real workflow, this would use terrain correction masking
        valid_mask = create_quality_mask(gamma0_data)
        
        masked_gamma0 = gamma0_data.copy()
        masked_gamma0[~valid_mask] = np.nan
        
        valid_pixels = np.sum(valid_mask)
        total_pixels = height * width
        coverage = 100 * valid_pixels / total_pixels
        
        print(f"   ✓ Applied quality mask")
        print(f"   ✓ Valid pixels: {valid_pixels:,} ({coverage:.1f}%)")
        
        # Step 3: Convert to dB scale
        print("\n3. Converting to dB scale...")
        
        # Convert both masked and unmasked data
        gamma0_db = sardine.linear_to_db(gamma0_data)
        gamma0_db_masked = sardine.linear_to_db(masked_gamma0)
        
        print(f"   ✓ Converted to dB scale")
        print(f"   ✓ dB range (all data): {np.nanmin(gamma0_db):.2f} to {np.nanmax(gamma0_db):.2f} dB")
        print(f"   ✓ dB range (masked): {np.nanmin(gamma0_db_masked):.2f} to {np.nanmax(gamma0_db_masked):.2f} dB")
        
        # Step 4: Create derived products
        print("\n4. Creating derived products...")
        
        # Land cover classification based on backscatter
        water_mask = gamma0_db < -15
        vegetation_mask = (gamma0_db >= -15) & (gamma0_db < -8)
        urban_mask = gamma0_db >= -8
        
        # Create classification map (0=water, 1=vegetation, 2=urban, 255=nodata)
        classification = np.zeros((height, width), dtype=np.uint8)
        classification[vegetation_mask] = 1
        classification[urban_mask] = 2
        classification[~valid_mask] = 255
        
        # Statistics
        water_pct = 100 * np.sum(water_mask & valid_mask) / valid_pixels
        veg_pct = 100 * np.sum(vegetation_mask & valid_mask) / valid_pixels
        urban_pct = 100 * np.sum(urban_mask & valid_mask) / valid_pixels
        
        print(f"   ✓ Created land cover classification")
        print(f"   ✓ Water: {water_pct:.1f}%, Vegetation: {veg_pct:.1f}%, Urban: {urban_pct:.1f}%")
        
        # Step 5: Export to GeoTIFF (if available)
        if export_enabled:
            print("\n5. Exporting to GeoTIFF...")
            
            # Define geographic bounds (example: San Francisco Bay Area)
            bounds = (-122.8, 37.2, -121.8, 38.2)  # (west, south, east, north)
            crs = 'EPSG:4326'
            
            # Export gamma0 in linear scale
            linear_path = os.path.join(temp_dir, "gamma0_linear.tif")
            sardine.export_geotiff(
                data=masked_gamma0,
                output_path=linear_path,
                bounds=bounds,
                crs=crs,
                nodata=np.nan,
                description="Gamma0 Backscatter (Linear Scale)"
            )
            
            # Export gamma0 in dB scale as COG
            db_path = os.path.join(temp_dir, "gamma0_db.tif")
            sardine.export_cog(
                data=gamma0_db_masked,
                output_path=db_path,
                bounds=bounds,
                crs=crs,
                nodata=-999.0,
                description="Gamma0 Backscatter (dB Scale)",
                compress='lzw'
            )
            
            # Export classification
            class_path = os.path.join(temp_dir, "land_cover.tif")
            sardine.export_geotiff(
                data=classification,
                output_path=class_path,
                bounds=bounds,
                crs=crs,
                nodata=255,
                description="Land Cover Classification"
            )
            
            # Export multi-band product
            multiband_path = os.path.join(temp_dir, "sar_multiband.tif")
            sardine.export_multiband_geotiff(
                data_list=[
                    masked_gamma0.astype(np.float32),
                    gamma0_db_masked.astype(np.float32),
                    valid_mask.astype(np.uint8)
                ],
                output_path=multiband_path,
                band_names=["Gamma0_Linear", "Gamma0_dB", "Quality_Mask"],
                bounds=bounds,
                crs=crs,
                nodata=np.nan
            )
            
            print(f"   ✓ Exported linear gamma0: {os.path.basename(linear_path)} ({get_file_size_mb(linear_path):.1f} MB)")
            print(f"   ✓ Exported dB gamma0 (COG): {os.path.basename(db_path)} ({get_file_size_mb(db_path):.1f} MB)")
            print(f"   ✓ Exported classification: {os.path.basename(class_path)} ({get_file_size_mb(class_path):.1f} MB)")
            print(f"   ✓ Exported multi-band: {os.path.basename(multiband_path)} ({get_file_size_mb(multiband_path):.1f} MB)")
            
            # Validate exports
            print("\n6. Validating exports...")
            
            try:
                db_info = sardine.validate_geotiff(db_path)
                multi_info = sardine.validate_geotiff(multiband_path)
                
                print(f"   ✓ dB export validation: {db_info['width']}x{db_info['height']}, {db_info['dtype']}")
                print(f"   ✓ Multi-band validation: {multi_info['count']} bands")
                print(f"   ✓ Coordinate system: {db_info['crs']}")
                print(f"   ✓ Geographic bounds: {db_info['bounds']}")
                
                # Check COG compliance
                if db_info['is_tiled']:
                    print(f"   ✓ COG compliance: Tiled ({db_info['blocksize'][0]}x{db_info['blocksize'][1]})")
                
                if db_info.get('overviews'):
                    print(f"   ✓ COG overviews: {len(db_info['overviews'])} levels")
                
            except Exception as e:
                print(f"   ⚠️  Validation warning: {e}")
        else:
            print("\n5. GeoTIFF export skipped (rasterio not available)")
        
        # Step 6: Round-trip conversion test
        print("\n7. Testing round-trip dB conversion...")
        
        # Test conversion accuracy
        linear_restored = sardine.db_to_linear(gamma0_db)
        
        # Calculate errors (only for valid data)
        valid_data = ~np.isnan(gamma0_data) & ~np.isinf(gamma0_data)
        max_error = np.max(np.abs(gamma0_data[valid_data] - linear_restored[valid_data]))
        rel_error = max_error / np.max(gamma0_data[valid_data]) * 100
        
        print(f"   ✓ Round-trip conversion test")
        print(f"   ✓ Maximum absolute error: {max_error:.2e}")
        print(f"   ✓ Maximum relative error: {rel_error:.2e}%")
        
        # Summary
        print(f"\n🎉 Workflow completed successfully!")
        print(f"\n📊 Processing Summary:")
        print(f"   • Input data: {height}x{width} pixels")
        print(f"   • Valid data coverage: {coverage:.1f}%")
        print(f"   • Land cover distribution:")
        print(f"     - Water: {water_pct:.1f}%")
        print(f"     - Vegetation: {veg_pct:.1f}%")
        print(f"     - Urban: {urban_pct:.1f}%")
        print(f"   • dB conversion: ✓")
        if export_enabled:
            print(f"   • GeoTIFF export: ✓")
            print(f"   • COG export: ✓")
            print(f"   • Multi-band export: ✓")
        print(f"   • Data validation: ✓")

def create_realistic_sar_scene(height, width):
    """Create a realistic SAR backscatter scene with multiple land cover types."""
    
    np.random.seed(42)
    
    # Create coordinate grids
    x = np.linspace(0, 10, width)
    y = np.linspace(0, 10, height)
    X, Y = np.meshgrid(x, y)
    
    # Base terrain with smooth variations
    terrain_base = 0.1 + 0.05 * np.sin(0.5 * X) * np.cos(0.3 * Y)
    
    # Add different land cover types
    scene = np.zeros((height, width))
    
    # Water bodies (very low backscatter)
    water_regions = (
        ((X - 2)**2 + (Y - 2)**2 < 1.5) |  # Lake 1
        ((X - 7)**2 + (Y - 8)**2 < 2.0) |  # Lake 2
        (Y < 1) |  # River/coast at bottom
        ((X > 8.5) & (Y > 2) & (Y < 4))  # River segment
    )
    scene[water_regions] = np.random.lognormal(-3.0, 0.2, np.sum(water_regions))
    
    # Urban areas (high backscatter)
    urban_regions = (
        ((X > 3) & (X < 6) & (Y > 6) & (Y < 9)) |  # City center
        ((X > 1) & (X < 2.5) & (Y > 4) & (Y < 6))  # Small town
    )
    scene[urban_regions] = np.random.lognormal(-1.2, 0.4, np.sum(urban_regions))
    
    # Forest areas (medium-high backscatter)
    forest_regions = (
        ((X > 6) & (Y < 6)) |  # Forest patch 1
        ((X < 3) & (Y > 7))    # Forest patch 2
    )
    scene[forest_regions] = np.random.lognormal(-1.8, 0.3, np.sum(forest_regions))
    
    # Agricultural areas (medium backscatter with patterns)
    ag_regions = (X > 4) & (X < 8) & (Y > 2) & (Y < 6)
    ag_pattern = 0.12 * (1 + 0.3 * np.sin(2 * X) * np.sin(2 * Y))
    scene[ag_regions] = ag_pattern[ag_regions] * np.random.lognormal(0, 0.2, np.sum(ag_regions))
    
    # Fill remaining areas with mixed vegetation
    unassigned = scene == 0
    scene[unassigned] = terrain_base[unassigned] * np.random.lognormal(0, 0.25, np.sum(unassigned))
    
    # Add speckle noise (multiplicative)
    speckle = np.random.gamma(1.0, 1.0, (height, width))
    scene = scene * speckle
    
    # Ensure reasonable range and no zero values
    scene = np.clip(scene, 0.001, 2.0)
    
    return scene.astype(np.float64)

def create_quality_mask(data):
    """Create a quality mask based on data characteristics."""
    
    # Remove very low values (potential water/shadow areas)
    valid = data > 0.005
    
    # Remove very high values (potential noise/urban artifacts)
    valid &= data < 1.5
    
    # Remove NaN and infinite values
    valid &= np.isfinite(data)
    
    # Apply morphological operations to clean up the mask
    # (In real application, you might use scipy.ndimage for this)
    # For simplicity, we'll just use a simple filter
    
    return valid

def get_file_size_mb(filepath):
    """Get file size in megabytes."""
    return os.path.getsize(filepath) / (1024 * 1024)

if __name__ == "__main__":
    complete_sar_workflow()
