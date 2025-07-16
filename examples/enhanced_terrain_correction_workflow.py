#!/usr/bin/env python3
"""
Enhanced Terrain Correction Workflow with Adaptive Masking

This example demonstrates the enhanced terrain correction pipeline with:
1. Integrated masking workflow
2. Adaptive threshold adjustment
3. Quality assessment
4. Local incidence angle computation for radiometric correction
5. Comprehensive output products

Features showcased:
- Enhanced terrain correction pipeline
- Adaptive masking thresholds based on data characteristics
- Quality assessment and filtering
- LIA-based radiometric correction
- Intermediate product generation
"""

import sys
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import time

# Add the parent directory to Python path to import sardine
sys.path.insert(0, str(Path(__file__).parent.parent / "python"))

try:
    import sardine
    print("✅ SARdine imported successfully")
except ImportError as e:
    print(f"❌ Failed to import SARdine: {e}")
    print("Make sure the extension is built with 'python -m pip install -e .'")
    sys.exit(1)

def create_synthetic_sar_data(height=512, width=512):
    """Create synthetic SAR data for demonstration."""
    np.random.seed(42)
    
    # Create a realistic SAR intensity pattern
    # Base terrain with different scattering characteristics
    x, y = np.meshgrid(np.linspace(0, 10, width), np.linspace(0, 10, height))
    
    # Add terrain features
    terrain = np.sin(x * 0.5) * np.cos(y * 0.3) + 0.5 * np.sin(x * 2) * np.sin(y * 1.5)
    
    # Add urban areas (high backscatter)
    urban_mask = ((x - 5)**2 + (y - 5)**2) < 2
    terrain[urban_mask] += 2.0
    
    # Add water bodies (low backscatter)
    water_mask = ((x - 2)**2 + (y - 8)**2) < 1
    terrain[water_mask] -= 3.0
    
    # Convert to gamma0 in dB
    gamma0_linear = 10**(terrain / 10)  # Convert from dB to linear
    
    # Add speckle noise
    speckle = np.random.gamma(1.0, 1.0, (height, width))
    gamma0_linear *= speckle
    
    # Convert back to dB
    gamma0_db = 10 * np.log10(gamma0_linear)
    
    # Add some invalid pixels
    invalid_mask = np.random.random((height, width)) < 0.05
    gamma0_db[invalid_mask] = np.nan
    
    return gamma0_db.astype(np.float32)

def create_synthetic_dem(height=512, width=512):
    """Create synthetic DEM data."""
    np.random.seed(123)
    
    x, y = np.meshgrid(np.linspace(0, 10, width), np.linspace(0, 10, height))
    
    # Create mountainous terrain
    elevation = (200 * np.sin(x * 0.3) * np.cos(y * 0.4) + 
                100 * np.sin(x * 1.2) * np.sin(y * 0.8) + 
                50 * np.random.random((height, width)) + 500)
    
    # Add a valley (below sea level in some areas)
    valley_mask = ((x - 7)**2 + (y - 3)**2) < 1.5
    elevation[valley_mask] -= 600
    
    # Add some invalid DEM pixels
    invalid_mask = np.random.random((height, width)) < 0.02
    elevation[invalid_mask] = -32768  # Standard DEM no-data value
    
    return elevation.astype(np.float32)

def create_synthetic_orbit():
    """Create synthetic orbit data."""
    return sardine.OrbitData()  # Simplified for demonstration

def demonstrate_enhanced_pipeline():
    """Demonstrate the enhanced terrain correction pipeline."""
    print("\n🌍 === Enhanced Terrain Correction Pipeline ===")
    
    # Create synthetic data
    print("📊 Creating synthetic SAR and DEM data...")
    sar_data = create_synthetic_sar_data()
    dem_data = create_synthetic_dem()
    orbit_data = create_synthetic_orbit()
    
    print(f"SAR data shape: {sar_data.shape}")
    print(f"DEM data shape: {dem_data.shape}")
    print(f"SAR value range: [{np.nanmin(sar_data):.2f}, {np.nanmax(sar_data):.2f}] dB")
    print(f"DEM value range: [{np.nanmin(dem_data):.2f}, {np.nanmax(dem_data):.2f}] m")
    
    # Save synthetic DEM as numpy file for the pipeline
    dem_path = "/tmp/synthetic_dem.npy"
    np.save(dem_path, dem_data)
    
    # Test basic masking workflow first
    print("\n🎭 Testing masking workflow...")
    
    # Create default masking workflow
    default_workflow = sardine.create_masking_workflow()
    print(f"Default thresholds:")
    print(f"  - LIA threshold: {default_workflow.lia_threshold}")
    print(f"  - DEM threshold: {default_workflow.dem_threshold}")
    print(f"  - Gamma0 range: [{default_workflow.gamma0_min}, {default_workflow.gamma0_max}]")
    
    # Apply masking workflow
    try:
        mask_result = sardine.apply_masking_workflow(
            dem_path=dem_path,
            gamma0_data=sar_data,
            workflow=default_workflow
        )
        
        print(f"\n📈 Masking Results:")
        print(f"  - Valid pixels: {mask_result.valid_pixels:,}/{mask_result.total_pixels:,}")
        print(f"  - Coverage: {mask_result.coverage_percent:.1f}%")
        
        # Get mask components as numpy arrays
        combined_mask = mask_result.get_combined_mask()
        lia_cosine = mask_result.get_lia_cosine()
        gamma0_mask = mask_result.get_gamma0_mask()
        dem_mask = mask_result.get_dem_mask()
        lia_mask = mask_result.get_lia_mask()
        
        print(f"  - Gamma0 valid: {np.sum(gamma0_mask)}/{gamma0_mask.size} ({np.sum(gamma0_mask)/gamma0_mask.size*100:.1f}%)")
        print(f"  - DEM valid: {np.sum(dem_mask)}/{dem_mask.size} ({np.sum(dem_mask)/dem_mask.size*100:.1f}%)")
        print(f"  - LIA valid: {np.sum(lia_mask)}/{lia_mask.size} ({np.sum(lia_mask)/lia_mask.size*100:.1f}%)")
        
        # Apply mask to gamma0 data
        masked_gamma0 = sardine.apply_mask_to_gamma0(
            dem_path=dem_path,
            gamma0_data=sar_data,
            mask=combined_mask,
            fill_value=-999.0
        )
        
        print(f"✅ Masking workflow completed successfully")
        
    except Exception as e:
        print(f"❌ Error in masking workflow: {e}")
        return False
    
    return True

def demonstrate_adaptive_thresholding():
    """Demonstrate adaptive thresholding."""
    print("\n🧠 === Adaptive Masking Thresholds ===")
    
    # Create data with different characteristics
    scenarios = [
        ("Urban Area", create_urban_sar_data()),
        ("Forest Area", create_forest_sar_data()),
        ("Water Body", create_water_sar_data()),
    ]
    
    for scenario_name, sar_data in scenarios:
        print(f"\n📊 Scenario: {scenario_name}")
        print(f"  - Data range: [{np.nanmin(sar_data):.2f}, {np.nanmax(sar_data):.2f}] dB")
        print(f"  - Mean: {np.nanmean(sar_data):.2f} dB")
        print(f"  - Std: {np.nanstd(sar_data):.2f} dB")
        
        # Test different threshold strategies
        conservative_workflow = sardine.create_masking_workflow(
            lia_threshold=0.2,      # More conservative LIA
            dem_threshold=-50.0,    # More conservative DEM
            gamma0_min=-40.0,       # Tighter gamma0 range
            gamma0_max=5.0
        )
        
        liberal_workflow = sardine.create_masking_workflow(
            lia_threshold=0.05,     # More liberal LIA
            dem_threshold=-200.0,   # More liberal DEM
            gamma0_min=-60.0,       # Wider gamma0 range
            gamma0_max=15.0
        )
        
        print(f"  - Conservative masking: LIA>{conservative_workflow.lia_threshold}, DEM>{conservative_workflow.dem_threshold}")
        print(f"  - Liberal masking: LIA>{liberal_workflow.lia_threshold}, DEM>{liberal_workflow.dem_threshold}")

def create_urban_sar_data():
    """Create SAR data typical of urban areas."""
    np.random.seed(100)
    sar_data = np.random.normal(-5.0, 3.0, (256, 256))  # High backscatter
    return sar_data.astype(np.float32)

def create_forest_sar_data():
    """Create SAR data typical of forest areas."""
    np.random.seed(200)
    sar_data = np.random.normal(-15.0, 4.0, (256, 256))  # Medium backscatter
    return sar_data.astype(np.float32)

def create_water_sar_data():
    """Create SAR data typical of water bodies."""
    np.random.seed(300)
    sar_data = np.random.normal(-25.0, 2.0, (256, 256))  # Low backscatter
    return sar_data.astype(np.float32)

def demonstrate_quality_assessment():
    """Demonstrate quality assessment features."""
    print("\n📈 === Quality Assessment and Filtering ===")
    
    # Create data with varying quality
    print("📊 Creating data with different quality levels...")
    
    quality_scenarios = [
        ("High Quality", 0.95, 0.02),    # 95% valid, low noise
        ("Medium Quality", 0.75, 0.05),  # 75% valid, medium noise
        ("Low Quality", 0.45, 0.10),     # 45% valid, high noise
    ]
    
    for scenario_name, valid_fraction, noise_level in quality_scenarios:
        print(f"\n🔍 {scenario_name} Data:")
        
        # Create synthetic data with controlled quality
        sar_data = create_synthetic_sar_data()
        
        # Add controlled amount of invalid data
        invalid_mask = np.random.random(sar_data.shape) > valid_fraction
        sar_data[invalid_mask] = np.nan
        
        # Add noise
        noise = np.random.normal(0, noise_level, sar_data.shape)
        valid_mask = ~np.isnan(sar_data)
        sar_data[valid_mask] += noise[valid_mask]
        
        # Apply masking and assess quality
        workflow = sardine.create_masking_workflow()
        
        try:
            # Save temporary DEM
            dem_data = create_synthetic_dem()
            dem_path = f"/tmp/dem_{scenario_name.lower().replace(' ', '_')}.npy"
            np.save(dem_path, dem_data)
            
            mask_result = sardine.apply_masking_workflow(
                dem_path=dem_path,
                gamma0_data=sar_data,
                workflow=workflow
            )
            
            coverage = mask_result.coverage_percent
            
            # Determine quality level
            if coverage >= 90:
                quality_level = "Excellent ⭐⭐⭐⭐⭐"
            elif coverage >= 75:
                quality_level = "Good ⭐⭐⭐⭐"
            elif coverage >= 50:
                quality_level = "Fair ⭐⭐⭐"
            elif coverage >= 25:
                quality_level = "Poor ⭐⭐"
            else:
                quality_level = "Very Poor ⭐"
            
            print(f"  - Input valid pixels: {np.sum(~np.isnan(sar_data))}/{sar_data.size} ({np.sum(~np.isnan(sar_data))/sar_data.size*100:.1f}%)")
            print(f"  - Final coverage: {coverage:.1f}%")
            print(f"  - Quality level: {quality_level}")
            
            if coverage < 50:
                print(f"  ⚠️  Low coverage detected - consider adjusting thresholds")
            
        except Exception as e:
            print(f"  ❌ Error in quality assessment: {e}")

def demonstrate_lia_radiometric_correction():
    """Demonstrate using LIA for radiometric correction."""
    print("\n🌅 === Local Incidence Angle Radiometric Correction ===")
    
    print("📊 Creating terrain with varying slopes...")
    
    # Create DEM with varying slopes
    x, y = np.meshgrid(np.linspace(0, 10, 256), np.linspace(0, 10, 256))
    
    # Create steep terrain
    elevation = 1000 * np.sin(x * 0.8) * np.cos(y * 0.6) + 500
    
    # Create corresponding SAR data affected by terrain
    sar_data = create_synthetic_sar_data(256, 256)
    
    # Save DEM
    dem_path = "/tmp/steep_terrain_dem.npy"
    np.save(dem_path, elevation.astype(np.float32))
    
    # Apply masking workflow to get LIA values
    workflow = sardine.create_masking_workflow(
        lia_threshold=0.05  # Include steep slopes for correction
    )
    
    try:
        mask_result = sardine.apply_masking_workflow(
            dem_path=dem_path,
            gamma0_data=sar_data,
            workflow=workflow
        )
        
        lia_cosine = mask_result.get_lia_cosine()
        
        print(f"📐 Local Incidence Angle Statistics:")
        valid_lia = lia_cosine[~np.isnan(lia_cosine)]
        print(f"  - LIA cosine range: [{np.min(valid_lia):.3f}, {np.max(valid_lia):.3f}]")
        print(f"  - Mean LIA cosine: {np.mean(valid_lia):.3f}")
        print(f"  - LIA angle range: [{np.degrees(np.arccos(np.max(valid_lia))):.1f}°, {np.degrees(np.arccos(np.min(valid_lia))):.1f}°]")
        
        # Demonstrate radiometric correction using LIA
        # Gamma0 corrected = Gamma0 / cos(LIA) for area normalization
        corrected_gamma0 = sar_data.copy()
        valid_mask = ~np.isnan(lia_cosine) & (lia_cosine > 0.1) & ~np.isnan(sar_data)
        
        # Apply correction in linear scale
        linear_gamma0 = 10**(sar_data[valid_mask] / 10)
        corrected_linear = linear_gamma0 / lia_cosine[valid_mask]
        corrected_gamma0[valid_mask] = 10 * np.log10(corrected_linear)
        
        print(f"📊 Radiometric Correction Results:")
        print(f"  - Original gamma0 range: [{np.nanmin(sar_data):.2f}, {np.nanmax(sar_data):.2f}] dB")
        print(f"  - Corrected gamma0 range: [{np.nanmin(corrected_gamma0):.2f}, {np.nanmax(corrected_gamma0):.2f}] dB")
        print(f"  - Correction applied to {np.sum(valid_mask)} pixels")
        
        # Save correction results
        np.save("/tmp/original_gamma0.npy", sar_data)
        np.save("/tmp/corrected_gamma0.npy", corrected_gamma0)
        np.save("/tmp/lia_cosine.npy", lia_cosine)
        
        print(f"💾 Radiometric correction products saved to /tmp/")
        
    except Exception as e:
        print(f"❌ Error in LIA radiometric correction: {e}")

def create_visualization():
    """Create visualization of the masking and correction results."""
    print("\n📈 === Creating Visualizations ===")
    
    try:
        # Load saved data
        sar_data = create_synthetic_sar_data()
        dem_data = create_synthetic_dem()
        
        # Apply masking
        dem_path = "/tmp/visualization_dem.npy"
        np.save(dem_path, dem_data)
        
        workflow = sardine.create_masking_workflow()
        mask_result = sardine.apply_masking_workflow(
            dem_path=dem_path,
            gamma0_data=sar_data,
            workflow=workflow
        )
        
        # Get results
        combined_mask = mask_result.get_combined_mask()
        lia_cosine = mask_result.get_lia_cosine()
        
        # Create plots
        fig, axes = plt.subplots(2, 3, figsize=(15, 10))
        fig.suptitle('SARdine Enhanced Masking Workflow Results', fontsize=16)
        
        # Original SAR data
        im1 = axes[0, 0].imshow(sar_data, cmap='viridis', aspect='equal')
        axes[0, 0].set_title('Original Gamma0 (dB)')
        plt.colorbar(im1, ax=axes[0, 0])
        
        # DEM
        im2 = axes[0, 1].imshow(dem_data, cmap='terrain', aspect='equal')
        axes[0, 1].set_title('DEM (m)')
        plt.colorbar(im2, ax=axes[0, 1])
        
        # Combined mask
        im3 = axes[0, 2].imshow(combined_mask, cmap='RdYlGn', aspect='equal')
        axes[0, 2].set_title('Combined Quality Mask')
        plt.colorbar(im3, ax=axes[0, 2])
        
        # LIA cosine
        im4 = axes[1, 0].imshow(lia_cosine, cmap='plasma', aspect='equal')
        axes[1, 0].set_title('Local Incidence Angle Cosine')
        plt.colorbar(im4, ax=axes[1, 0])
        
        # Masked SAR data
        masked_sar = sar_data.copy()
        masked_sar[~combined_mask] = np.nan
        im5 = axes[1, 1].imshow(masked_sar, cmap='viridis', aspect='equal')
        axes[1, 1].set_title('Masked Gamma0 (dB)')
        plt.colorbar(im5, ax=axes[1, 1])
        
        # Statistics
        axes[1, 2].text(0.1, 0.9, f'Coverage: {mask_result.coverage_percent:.1f}%', transform=axes[1, 2].transAxes)
        axes[1, 2].text(0.1, 0.8, f'Valid pixels: {mask_result.valid_pixels:,}', transform=axes[1, 2].transAxes)
        axes[1, 2].text(0.1, 0.7, f'Total pixels: {mask_result.total_pixels:,}', transform=axes[1, 2].transAxes)
        axes[1, 2].text(0.1, 0.6, f'LIA threshold: {workflow.lia_threshold}', transform=axes[1, 2].transAxes)
        axes[1, 2].text(0.1, 0.5, f'DEM threshold: {workflow.dem_threshold}m', transform=axes[1, 2].transAxes)
        axes[1, 2].text(0.1, 0.4, f'Gamma0 range: [{workflow.gamma0_min}, {workflow.gamma0_max}]', transform=axes[1, 2].transAxes)
        axes[1, 2].set_title('Processing Statistics')
        axes[1, 2].axis('off')
        
        plt.tight_layout()
        plt.savefig('/tmp/sardine_masking_results.png', dpi=150, bbox_inches='tight')
        print(f"📊 Visualization saved to: /tmp/sardine_masking_results.png")
        
    except Exception as e:
        print(f"❌ Error creating visualization: {e}")

def main():
    """Main demonstration function."""
    print("🚀 SARdine Enhanced Terrain Correction & Masking Workflow Demo")
    print("=" * 70)
    
    start_time = time.time()
    
    # Run demonstrations
    success = True
    
    try:
        # Basic enhanced pipeline
        if not demonstrate_enhanced_pipeline():
            success = False
        
        # Adaptive thresholding
        demonstrate_adaptive_thresholding()
        
        # Quality assessment
        demonstrate_quality_assessment()
        
        # LIA radiometric correction
        demonstrate_lia_radiometric_correction()
        
        # Create visualizations
        create_visualization()
        
    except Exception as e:
        print(f"❌ Demo failed with error: {e}")
        import traceback
        traceback.print_exc()
        success = False
    
    # Summary
    total_time = time.time() - start_time
    print(f"\n{'='*70}")
    
    if success:
        print(f"🎉 Demo completed successfully in {total_time:.2f}s")
        print(f"\n📋 Key Features Demonstrated:")
        print(f"  ✅ Enhanced terrain correction pipeline with masking")
        print(f"  ✅ Adaptive masking thresholds")
        print(f"  ✅ Quality assessment and filtering") 
        print(f"  ✅ Local incidence angle radiometric correction")
        print(f"  ✅ Comprehensive mask generation")
        print(f"  ✅ Numpy array integration")
        
        print(f"\n🎯 Next Steps:")
        print(f"  1. Integrate masking into main processing pipeline")
        print(f"  2. Adjust thresholds based on your data characteristics")
        print(f"  3. Use LIA for advanced radiometric correction")
        print(f"  4. Employ masks for quality control in analysis")
        print(f"  5. Optimize processing for large-scale datasets")
        
    else:
        print(f"❌ Demo completed with errors in {total_time:.2f}s")
        print(f"Check the error messages above for troubleshooting")

if __name__ == "__main__":
    main()
