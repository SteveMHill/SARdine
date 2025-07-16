#!/usr/bin/env python3
"""
Simple masking workflow demonstration

This script demonstrates the basic masking workflow functionality
without requiring actual DEM files.
"""

import sys
import numpy as np
from pathlib import Path

# Add the parent directory to Python path to import sardine
sys.path.insert(0, str(Path(__file__).parent.parent / "python"))

try:
    import sardine
    print("✅ SARdine imported successfully")
except ImportError as e:
    print(f"❌ Failed to import SARdine: {e}")
    sys.exit(1)

def test_workflow_creation():
    """Test creation of different masking workflows"""
    print("\n🔧 Testing masking workflow creation...")
    
    # Test default workflow
    default_workflow = sardine.create_masking_workflow()
    print(f"✅ Default workflow created:")
    print(f"  - LIA threshold: {default_workflow.lia_threshold}")
    print(f"  - DEM threshold: {default_workflow.dem_threshold}")
    print(f"  - Gamma0 range: [{default_workflow.gamma0_min}, {default_workflow.gamma0_max}]")
    
    # Test custom workflows
    conservative = sardine.create_masking_workflow(
        lia_threshold=0.2,
        dem_threshold=-50.0,
        gamma0_min=-40.0,
        gamma0_max=5.0
    )
    print(f"✅ Conservative workflow created:")
    print(f"  - LIA threshold: {conservative.lia_threshold}")
    print(f"  - DEM threshold: {conservative.dem_threshold}")
    print(f"  - Gamma0 range: [{conservative.gamma0_min}, {conservative.gamma0_max}]")
    
    liberal = sardine.create_masking_workflow(
        lia_threshold=0.05,
        dem_threshold=-200.0,
        gamma0_min=-60.0,
        gamma0_max=15.0
    )
    print(f"✅ Liberal workflow created:")
    print(f"  - LIA threshold: {liberal.lia_threshold}")
    print(f"  - DEM threshold: {liberal.dem_threshold}")
    print(f"  - Gamma0 range: [{liberal.gamma0_min}, {liberal.gamma0_max}]")
    
    return True

def test_workflow_parameters():
    """Test different workflow parameter combinations"""
    print("\n📊 Testing workflow parameter variations...")
    
    scenarios = [
        ("Urban Processing", {"lia_threshold": 0.15, "dem_threshold": -50.0, "gamma0_min": -30.0, "gamma0_max": 10.0}),
        ("Forest Processing", {"lia_threshold": 0.1, "dem_threshold": -100.0, "gamma0_min": -45.0, "gamma0_max": 5.0}),
        ("Water Processing", {"lia_threshold": 0.05, "dem_threshold": -200.0, "gamma0_min": -60.0, "gamma0_max": 0.0}),
        ("Mountain Processing", {"lia_threshold": 0.2, "dem_threshold": 0.0, "gamma0_min": -40.0, "gamma0_max": 15.0}),
    ]
    
    for scenario_name, params in scenarios:
        workflow = sardine.create_masking_workflow(**params)
        print(f"✅ {scenario_name}:")
        print(f"  - LIA ≥ {workflow.lia_threshold} (angle ≤ {np.degrees(np.arccos(workflow.lia_threshold)):.1f}°)")
        print(f"  - DEM > {workflow.dem_threshold}m")
        print(f"  - γ₀ ∈ [{workflow.gamma0_min}, {workflow.gamma0_max}] dB")
    
    return True

def demonstrate_threshold_effects():
    """Demonstrate the effects of different thresholds"""
    print("\n🎯 Demonstrating threshold effects...")
    
    # Simulate different data characteristics
    data_types = {
        "Flat Terrain": {"mean_lia": 0.9, "std_lia": 0.1, "mean_gamma0": -12, "std_gamma0": 3},
        "Hilly Terrain": {"mean_lia": 0.5, "std_lia": 0.3, "mean_gamma0": -15, "std_gamma0": 5},
        "Steep Mountains": {"mean_lia": 0.2, "std_lia": 0.2, "mean_gamma0": -18, "std_gamma0": 8},
        "Mixed Terrain": {"mean_lia": 0.4, "std_lia": 0.4, "mean_gamma0": -14, "std_gamma0": 6},
    }
    
    thresholds = [0.05, 0.1, 0.2, 0.3]
    
    for terrain_type, characteristics in data_types.items():
        print(f"\n📍 {terrain_type}:")
        print(f"  Characteristics: LIA~{characteristics['mean_lia']:.1f}±{characteristics['std_lia']:.1f}, "
              f"γ₀~{characteristics['mean_gamma0']}±{characteristics['std_gamma0']} dB")
        
        # Simulate data
        size = 1000
        lia_values = np.clip(np.random.normal(characteristics['mean_lia'], characteristics['std_lia'], size), 0, 1)
        gamma0_values = np.random.normal(characteristics['mean_gamma0'], characteristics['std_gamma0'], size)
        
        for threshold in thresholds:
            # Calculate what percentage would pass different thresholds
            lia_pass_rate = np.sum(lia_values >= threshold) / size * 100
            gamma0_pass_rate = np.sum((gamma0_values >= -50) & (gamma0_values <= 10)) / size * 100
            combined_pass_rate = min(lia_pass_rate, gamma0_pass_rate)  # Simplified estimation
            
            angle_threshold = np.degrees(np.arccos(threshold))
            print(f"    LIA ≥ {threshold:.2f} (≤{angle_threshold:.0f}°): {lia_pass_rate:.1f}% pass")
    
    return True

def calculate_quality_metrics():
    """Calculate and display quality assessment metrics"""
    print("\n📈 Quality Assessment Metrics...")
    
    coverage_levels = [95, 85, 75, 65, 50, 35, 20]
    
    print("Coverage Level → Quality Rating:")
    for coverage in coverage_levels:
        if coverage >= 90:
            quality = "Excellent ⭐⭐⭐⭐⭐"
        elif coverage >= 75:
            quality = "Good ⭐⭐⭐⭐"
        elif coverage >= 50:
            quality = "Fair ⭐⭐⭐"
        elif coverage >= 25:
            quality = "Poor ⭐⭐"
        else:
            quality = "Very Poor ⭐"
        
        print(f"  {coverage:2d}% → {quality}")
    
    print(f"\n💡 Recommendations:")
    print(f"  • Coverage ≥ 90%: Excellent for analysis")
    print(f"  • Coverage 75-89%: Good for most applications")
    print(f"  • Coverage 50-74%: Consider threshold adjustment")
    print(f"  • Coverage < 50%: Review processing parameters")
    
    return True

def demonstrate_radiometric_correction():
    """Demonstrate LIA-based radiometric correction principles"""
    print("\n🌅 Local Incidence Angle Radiometric Correction...")
    
    # Simulate different terrain slopes and their LIA values
    slopes = np.array([0, 10, 20, 30, 40, 50, 60])  # degrees
    look_angle = 35  # Typical SAR look angle in degrees
    
    print(f"SAR Look Angle: {look_angle}°")
    print(f"Terrain Slope → Local Incidence Angle → Correction Factor:")
    
    for slope in slopes:
        # Simplified calculation (actual calculation is more complex)
        lia = abs(look_angle - slope)  # Simplified
        lia_rad = np.radians(lia)
        lia_cosine = np.cos(lia_rad)
        
        # Correction factor for area normalization
        correction_factor = 1.0 / lia_cosine if lia_cosine > 0.1 else "N/A (too steep)"
        
        if isinstance(correction_factor, float):
            correction_db = 10 * np.log10(correction_factor)
            print(f"  {slope:2d}° → {lia:4.1f}° → {lia_cosine:.3f} → {correction_factor:.2f} ({correction_db:+.1f} dB)")
        else:
            print(f"  {slope:2d}° → {lia:4.1f}° → {lia_cosine:.3f} → {correction_factor}")
    
    print(f"\n📝 Notes:")
    print(f"  • Correction = 1/cos(LIA) for area normalization")
    print(f"  • Steep slopes (LIA > ~84°) may need different treatment")
    print(f"  • Correction reduces terrain-induced radiometric variations")
    
    return True

def main():
    """Main demonstration function"""
    print("🚀 SARdine Enhanced Masking Workflow Demonstration")
    print("=" * 60)
    
    success = True
    
    try:
        # Test basic functionality
        if not test_workflow_creation():
            success = False
        
        if not test_workflow_parameters():
            success = False
        
        if not demonstrate_threshold_effects():
            success = False
        
        if not calculate_quality_metrics():
            success = False
        
        if not demonstrate_radiometric_correction():
            success = False
        
    except Exception as e:
        print(f"❌ Demo failed with error: {e}")
        import traceback
        traceback.print_exc()
        success = False
    
    # Summary
    print(f"\n{'='*60}")
    if success:
        print("🎉 Demonstration completed successfully!")
        print(f"\n📋 Features Demonstrated:")
        print(f"  ✅ Masking workflow creation and configuration")
        print(f"  ✅ Parameter variation for different scenarios")
        print(f"  ✅ Threshold effect analysis")
        print(f"  ✅ Quality assessment framework")
        print(f"  ✅ LIA-based radiometric correction principles")
        
        print(f"\n🎯 Next Steps for Full Implementation:")
        print(f"  1. Integrate with actual DEM files")
        print(f"  2. Process real SAR data")
        print(f"  3. Use enhanced terrain correction pipeline")
        print(f"  4. Apply adaptive thresholding")
        print(f"  5. Generate visualization outputs")
        
    else:
        print("❌ Demonstration completed with errors")
    
    return 0 if success else 1

if __name__ == "__main__":
    exit(main())
