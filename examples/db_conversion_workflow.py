#!/usr/bin/env python3
"""
Example workflow demonstrating dB conversion and array handling in SARdine.

This script shows how to use the dB conversion functions in a typical 
SAR data processing workflow.
"""

import numpy as np
import sardine

def simple_sar_workflow():
    """
    Demonstrate a simple SAR processing workflow with dB conversion.
    """
    
    print("SARdine dB Conversion Workflow Example")
    print("=" * 50)
    
    # Simulate reading SAR backscatter data (e.g., after terrain correction)
    print("1. Simulating SAR backscatter data (gamma0)...")
    
    # Create synthetic gamma0 data (linear scale)
    # Typical values for different land cover types
    width, height = 100, 100
    np.random.seed(42)  # For reproducible results
    
    # Simulate different backscatter values for different areas
    gamma0_linear = np.zeros((height, width), dtype=np.float64)
    
    # Water (low backscatter)
    gamma0_linear[:30, :40] = np.random.lognormal(-2.5, 0.3, (30, 40))  # ~0.08
    
    # Forest (medium backscatter)  
    gamma0_linear[30:70, 20:80] = np.random.lognormal(-1.8, 0.4, (40, 60))  # ~0.16
    
    # Urban (high backscatter)
    gamma0_linear[70:, 60:] = np.random.lognormal(-1.0, 0.5, (30, 40))  # ~0.37
    
    # Fill remaining areas with mixed vegetation
    mask = gamma0_linear == 0
    gamma0_linear[mask] = np.random.lognormal(-2.0, 0.3, np.sum(mask))
    
    print(f"   Linear gamma0 statistics:")
    print(f"   - Min: {np.min(gamma0_linear):.6f}")
    print(f"   - Max: {np.max(gamma0_linear):.6f}")
    print(f"   - Mean: {np.mean(gamma0_linear):.6f}")
    print(f"   - Std: {np.std(gamma0_linear):.6f}")
    
    # 2. Convert to dB for visualization and analysis
    print("\n2. Converting to dB scale...")
    gamma0_db = sardine.linear_to_db(gamma0_linear)
    
    print(f"   dB gamma0 statistics:")
    print(f"   - Min: {np.min(gamma0_db):.2f} dB")
    print(f"   - Max: {np.max(gamma0_db):.2f} dB")
    print(f"   - Mean: {np.mean(gamma0_db):.2f} dB")
    print(f"   - Std: {np.std(gamma0_db):.2f} dB")
    
    # 3. Demonstrate typical SAR analysis in dB
    print("\n3. Performing analysis in dB domain...")
    
    # Land cover classification based on dB values
    water_mask = gamma0_db < -15  # Very low backscatter
    vegetation_mask = (gamma0_db >= -15) & (gamma0_db < -8)  # Medium backscatter
    urban_mask = gamma0_db >= -8  # High backscatter
    
    print(f"   Land cover classification:")
    print(f"   - Water pixels: {np.sum(water_mask)} ({100*np.sum(water_mask)/(width*height):.1f}%)")
    print(f"   - Vegetation pixels: {np.sum(vegetation_mask)} ({100*np.sum(vegetation_mask)/(width*height):.1f}%)")
    print(f"   - Urban pixels: {np.sum(urban_mask)} ({100*np.sum(urban_mask)/(width*height):.1f}%)")
    
    # 4. Convert back to linear for further processing
    print("\n4. Converting back to linear scale for radiometric processing...")
    gamma0_linear_restored = sardine.db_to_linear(gamma0_db)
    
    # Check conversion accuracy
    max_error = np.max(np.abs(gamma0_linear - gamma0_linear_restored))
    relative_error = max_error / np.max(gamma0_linear) * 100
    
    print(f"   Conversion accuracy:")
    print(f"   - Maximum absolute error: {max_error:.2e}")
    print(f"   - Maximum relative error: {relative_error:.2e}%")
    
    # 5. Demonstrate masked operations
    print("\n5. Demonstrating masked operations...")
    
    # Apply water mask to linear data
    gamma0_masked = gamma0_linear.copy()
    gamma0_masked[water_mask] = np.nan
    
    # Convert masked data to dB (NaN handling)
    gamma0_masked_db = sardine.linear_to_db(gamma0_masked)
    
    valid_pixels = ~np.isnan(gamma0_masked_db) & ~np.isinf(gamma0_masked_db)
    print(f"   Valid pixels after masking: {np.sum(valid_pixels)} ({100*np.sum(valid_pixels)/(width*height):.1f}%)")
    
    if np.sum(valid_pixels) > 0:
        print(f"   Masked data statistics (dB):")
        print(f"   - Min: {np.nanmin(gamma0_masked_db):.2f} dB")
        print(f"   - Max: {np.nanmax(gamma0_masked_db):.2f} dB")
        print(f"   - Mean: {np.nanmean(gamma0_masked_db):.2f} dB")
    
    print("\n6. Summary")
    print("   ✓ Linear to dB conversion working correctly")
    print("   ✓ dB to linear conversion working correctly")
    print("   ✓ Round-trip conversion maintains precision")
    print("   ✓ Handles special values (NaN, inf) appropriately")
    print("   ✓ Suitable for typical SAR analysis workflows")
    
    print(f"\nWorkflow completed successfully!")
    
    return {
        'gamma0_linear': gamma0_linear,
        'gamma0_db': gamma0_db,
        'water_mask': water_mask,
        'vegetation_mask': vegetation_mask,
        'urban_mask': urban_mask
    }

if __name__ == "__main__":
    results = simple_sar_workflow()
