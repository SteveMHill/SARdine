#!/usr/bin/env python3
"""
SARdine Usage Examples
Quick examples for common SAR processing tasks
"""

import os
import subprocess
import sys

def print_examples():
    """Print usage examples"""
    
    print("=" * 80)
    print("🛰️  SARdine SAR Processing - Usage Examples")
    print("=" * 80)
    
    print("\n📁 1. BASIC BACKSCATTER PROCESSING")
    print("   Process Sentinel-1 SLC to backscatter with default settings:")
    print("   python backscatter_cli.py data/S1A_*.zip ./output/")
    
    print("\n📡 2. DUAL POLARIZATION PROCESSING")
    print("   Process both VV and VH:")
    print("   python backscatter_cli.py data/S1A_*.zip ./output_vv/ --polarization VV")
    print("   python backscatter_cli.py data/S1A_*.zip ./output_vh/ --polarization VH")
    
    print("\n🔧 3. CUSTOM SPECKLE FILTERING")
    print("   Use different speckle filters:")
    print("   python backscatter_cli.py data/S1A_*.zip ./output/ --speckle-filter lee")
    print("   python backscatter_cli.py data/S1A_*.zip ./output/ --speckle-filter gamma_map")
    
    print("\n👁️  4. CUSTOM MULTILOOKING")
    print("   Adjust multilook factors:")
    print("   python backscatter_cli.py data/S1A_*.zip ./output/ --multilook 3 3")
    print("   python backscatter_cli.py data/S1A_*.zip ./output/ --multilook 1 1  # No multilooking")
    
    print("\n🗺️  5. HIGH-RESOLUTION PROCESSING")
    print("   Generate high-resolution products:")
    print("   python backscatter_cli.py data/S1A_*.zip ./output/ --resolution 5 --filter-window 5")
    
    print("\n⚡ 6. QUICK PROCESSING (NO GEOCODING)")
    print("   Fast processing without geocoding:")
    print("   python backscatter_cli.py data/S1A_*.zip ./output/ --no-geocode --no-terrain-flatten")
    
    print("\n🧪 7. VALIDATION TESTS")
    print("   Run validation tests:")
    print("   python test_14_step_pipeline.py")
    print("   python test_real_data_only.py")
    
    print("\n📊 8. COMPREHENSIVE PROCESSING")
    print("   Full processing with all options:")
    print("   python backscatter_cli.py data/S1A_*.zip ./output/ \\")
    print("     --polarization VV \\")
    print("     --speckle-filter enhanced_lee \\")
    print("     --filter-window 7 \\")
    print("     --multilook 2 2 \\")
    print("     --resolution 10")
    
    print("\n📋 OUTPUT FILES GENERATED:")
    print("   - backscatter_VV.tif      # Main GeoTIFF product")
    print("   - metadata.json           # Processing metadata")
    print("   - metadata.xml            # XML metadata")
    print("   - quality_report.json     # Quality assessment")
    print("   - processing_log.json     # Detailed processing log")
    
    print("\n🚀 REQUIREMENTS:")
    print("   - Sentinel-1 SLC ZIP file in data/ directory")
    print("   - SARdine Python package installed")
    print("   - Output directory (will be created automatically)")
    
    print("\n" + "=" * 80)

def check_requirements():
    """Check if requirements are met"""
    
    # Check if SARdine is available
    try:
        import sardine
        print("✅ SARdine package available")
    except ImportError:
        print("❌ SARdine package not found. Please install first.")
        return False
    
    # Check for data directory
    if os.path.exists("data/"):
        print("✅ Data directory found")
        
        # Look for SLC files
        slc_files = [f for f in os.listdir("data/") if f.endswith('.zip')]
        if slc_files:
            print(f"✅ Found {len(slc_files)} SLC files:")
            for f in slc_files[:3]:  # Show first 3
                print(f"   - {f}")
            if len(slc_files) > 3:
                print(f"   ... and {len(slc_files) - 3} more")
        else:
            print("⚠️  No SLC ZIP files found in data/ directory")
            print("   Place Sentinel-1 SLC ZIP files in data/ directory")
    else:
        print("⚠️  Data directory not found. Create it and add SLC files.")
    
    return True

def run_basic_example():
    """Run a basic example if data is available"""
    
    if not os.path.exists("data/"):
        print("❌ Cannot run example: data/ directory not found")
        return
    
    slc_files = [f for f in os.listdir("data/") if f.endswith('.zip')]
    if not slc_files:
        print("❌ Cannot run example: No SLC files found in data/")
        return
    
    # Use first SLC file
    slc_file = f"data/{slc_files[0]}"
    output_dir = "./example_output/"
    
    print(f"\n🚀 Running basic example with: {slc_files[0]}")
    print(f"   Output directory: {output_dir}")
    
    cmd = [
        "python", "backscatter_cli.py",
        slc_file, output_dir,
        "--no-geocode",  # For faster example
        "--no-terrain-flatten"
    ]
    
    print(f"   Command: {' '.join(cmd)}")
    print("\n" + "="*50)
    
    try:
        subprocess.run(cmd, check=True)
        print("\n✅ Example completed successfully!")
        print(f"   Check outputs in: {output_dir}")
    except subprocess.CalledProcessError as e:
        print(f"\n❌ Example failed with error code: {e.returncode}")
    except FileNotFoundError:
        print("\n❌ Could not find backscatter_cli.py")

if __name__ == "__main__":
    if len(sys.argv) > 1:
        if sys.argv[1] == "run":
            check_requirements()
            run_basic_example()
        elif sys.argv[1] == "check":
            check_requirements()
        else:
            print("Usage: python examples.py [run|check]")
            print("  run   - Run basic example")
            print("  check - Check requirements")
    else:
        print_examples()
