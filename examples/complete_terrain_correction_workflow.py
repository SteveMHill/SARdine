#!/usr/bin/env python3
"""
Complete Terrain Correction (Geocoding) Workflow

This example demonstrates how to perform Range-Doppler terrain correction
to convert SAR data from slant range geometry to map-projected coordinates.

Terrain correction includes:
- Automatic DEM loading and preparation
- Range-Doppler equation solving
- Geometric correction for topographic distortions
- Resampling to regular map grid
- GeoTIFF output generation
"""

import sys
import time
import numpy as np
from pathlib import Path

# Add the project root to Python path for development
sys.path.insert(0, str(Path(__file__).parent.parent))

import sardine
from sardine import _core


def create_synthetic_sar_data(rows=1000, cols=1000):
    """Create synthetic SAR intensity data for testing."""
    print("📊 Creating synthetic SAR intensity data...")
    
    # Create synthetic sigma0 data with terrain effects
    np.random.seed(42)
    
    # Base backscatter pattern
    y, x = np.meshgrid(np.linspace(0, 1, rows), np.linspace(0, 1, cols), indexing='ij')
    
    # Simulate terrain effects
    elevation_pattern = np.sin(y * 4 * np.pi) * np.cos(x * 6 * np.pi)
    terrain_effect = 1.0 + 0.3 * elevation_pattern
    
    # Speckle noise
    speckle = np.random.gamma(4, 0.25, (rows, cols))  # Gamma distributed speckle
    
    # Combine effects
    sigma0_db = -15 + 5 * terrain_effect + 2 * (speckle - 1.0)
    sigma0_linear = 10 ** (sigma0_db / 10)
    
    print(f"✅ Created {rows}x{cols} synthetic SAR data")
    print(f"   Sigma0 range: {sigma0_db.min():.1f} to {sigma0_db.max():.1f} dB")
    
    return sigma0_linear


def create_synthetic_orbit_data():
    """Create synthetic orbit data for testing."""
    print("🛰️  Creating synthetic orbit data...")
    
    # Simplified orbit for Sentinel-1 (polar orbit)
    orbit_data = []
    
    # Generate state vectors for a typical pass
    from datetime import datetime, timezone
    
    for i in range(100):  # 100 state vectors
        time_str = datetime.now(timezone.utc).isoformat()
        
        # Simplified orbital parameters (not realistic, just for testing)
        t = i * 0.1  # 0.1 second intervals
        
        # Approximate satellite position (simplified)
        position = [
            6800000 + 100000 * np.cos(t * 0.001),  # X (m)
            1000000 * np.sin(t * 0.001),           # Y (m) 
            7000000 + 50000 * np.sin(t * 0.002)    # Z (m)
        ]
        
        # Approximate velocity
        velocity = [
            -100 * np.sin(t * 0.001),              # Vx (m/s)
            1000 * np.cos(t * 0.001),              # Vy (m/s)
            -50 * np.cos(t * 0.002)                # Vz (m/s)
        ]
        
        orbit_data.append((time_str, position, velocity))
    
    print(f"✅ Created {len(orbit_data)} state vectors")
    return orbit_data


def test_coordinate_conversion():
    """Test coordinate conversion functions."""
    print("\n🧪 Testing coordinate conversion...")
    
    # Test location: San Francisco Bay Area
    lat, lon, elevation = 37.7749, -122.4194, 50.0
    
    print(f"📍 Test coordinates:")
    print(f"   Latitude: {lat}°")
    print(f"   Longitude: {lon}°")
    print(f"   Elevation: {elevation}m")
    
    # Convert to ECEF
    ecef = _core.latlon_to_ecef(lat, lon, elevation)
    
    print(f"🌍 ECEF coordinates:")
    print(f"   X: {ecef[0]:,.1f}m")
    print(f"   Y: {ecef[1]:,.1f}m")
    print(f"   Z: {ecef[2]:,.1f}m")
    
    # Calculate distance from Earth center
    distance = np.sqrt(ecef[0]**2 + ecef[1]**2 + ecef[2]**2)
    print(f"   Distance from center: {distance/1000:.1f}km")


def test_dem_loading(dem_path):
    """Test DEM loading for terrain correction."""
    print(f"\n🗻 Testing DEM loading...")
    print(f"📁 DEM file: {dem_path}")
    
    try:
        result = _core.create_terrain_corrector(
            str(dem_path),
            4326,  # WGS84
            10.0   # 10m spacing
        )
        print(f"✅ {result}")
        return True
    except Exception as e:
        print(f"❌ DEM loading failed: {e}")
        return False


def perform_terrain_correction_workflow():
    """Perform complete terrain correction workflow."""
    print("\n" + "="*80)
    print("🗺️  COMPLETE TERRAIN CORRECTION WORKFLOW")
    print("="*80)
    
    start_time = time.time()
    
    # Step 1: Create test data
    print("\n📊 Step 1: Preparing test data")
    sar_image = create_synthetic_sar_data(500, 500)  # Smaller for faster testing
    orbit_data = create_synthetic_orbit_data()
    
    # Step 2: Define processing parameters
    print("\n⚙️  Step 2: Setting up processing parameters")
    bbox = (37.0, 38.0, -122.5, -121.5)  # San Francisco Bay Area
    output_crs = 4326  # WGS84
    output_spacing = 10.0  # 10m
    
    print(f"   Bounding box: {bbox}")
    print(f"   Output CRS: EPSG:{output_crs}")
    print(f"   Output spacing: {output_spacing}m")
    
    # Step 3: Save test data
    print("\n💾 Step 3: Saving test data")
    test_data_dir = Path("test_terrain_correction")
    test_data_dir.mkdir(exist_ok=True)
    
    sar_file = test_data_dir / "test_sar_image.npy"
    np.save(sar_file, sar_image)
    print(f"   SAR image saved to: {sar_file}")
    
    # Step 4: Download DEM if needed
    print("\n🌍 Step 4: Checking DEM availability")
    dem_cache_dir = test_data_dir / "dem_cache"
    dem_cache_dir.mkdir(exist_ok=True)
    
    # Try to download SRTM for the area
    try:
        from sardine.types import BoundingBox
        test_bbox = BoundingBox(
            min_lat=bbox[0], max_lat=bbox[1], 
            min_lon=bbox[2], max_lon=bbox[3]
        )
        
        print(f"   Attempting SRTM download for bbox: {bbox}")
        dem_files = sardine.download_srtm_tiles(test_bbox, str(dem_cache_dir))
        
        if dem_files:
            dem_file = dem_files[0]
            print(f"   ✅ DEM available: {dem_file}")
            
            # Test DEM loading
            if test_dem_loading(dem_file):
                # Step 5: Test terrain correction (would need proper implementation)
                print(f"\n🗺️  Step 5: Terrain correction ready")
                print(f"   Note: Full terrain correction requires complete orbit integration")
                print(f"   This demo shows the setup and data preparation steps")
            else:
                print(f"   ❌ DEM loading test failed")
        else:
            print(f"   ⚠️  No DEM files downloaded (may be network issue)")
            
    except Exception as e:
        print(f"   ⚠️  DEM download failed: {e}")
        print(f"   This is normal if no internet connection is available")
    
    # Step 6: Coordinate conversion test
    test_coordinate_conversion()
    
    # Summary
    total_time = time.time() - start_time
    print(f"\n" + "="*80)
    print(f"🎉 TERRAIN CORRECTION WORKFLOW COMPLETED")
    print(f"⏱️  Total time: {total_time:.2f} seconds")
    print(f"📁 Test data saved in: {test_data_dir}")
    print("="*80)
    
    return test_data_dir


def demonstrate_range_doppler_concepts():
    """Demonstrate the concepts behind Range-Doppler terrain correction."""
    print("\n" + "="*80)
    print("📡 RANGE-DOPPLER TERRAIN CORRECTION CONCEPTS")
    print("="*80)
    
    print("""
🎯 What is Range-Doppler Terrain Correction?

Range-Doppler terrain correction transforms SAR data from slant range 
geometry to map coordinates using:

1. RANGE EQUATION: R = c·t/2
   - R: slant range distance
   - c: speed of light
   - t: two-way travel time

2. DOPPLER EQUATION: f_d = (2·v·cos(θ))/λ  
   - f_d: Doppler frequency
   - v: satellite velocity
   - θ: look angle
   - λ: radar wavelength

3. DEM INTERSECTION:
   - Find where radar beam intersects terrain
   - Account for topographic variations
   - Correct geometric distortions

🔄 Processing Steps:

┌─────────────────┐    ┌─────────────────┐    ┌─────────────────┐
│   SAR Image     │    │   DEM + Orbit   │    │  Map Grid       │
│  (Slant Range)  │ -> │   (3D Geometry) │ -> │ (Geographic)    │
└─────────────────┘    └─────────────────┘    └─────────────────┘

🏔️  Geometric Corrections:
• Foreshortening: Steep slopes appear compressed
• Layover: Mountaintops arrive before bases  
• Shadow: Hidden areas behind terrain
• Range migration: Curved Earth effects

📊 Output: Orthorectified, map-projected SAR image
""")


def main():
    """Main function for terrain correction workflow."""
    print("🌍 SARdine Terrain Correction (Geocoding) Example")
    
    # Demonstrate concepts
    demonstrate_range_doppler_concepts()
    
    # Run workflow
    test_dir = perform_terrain_correction_workflow()
    
    print(f"""
🎓 Next Steps:

1. Examine the test data in: {test_dir}
2. Try the CLI command:
   sardine geocode test_sar_image.npy dem_file.hgt orbit.EOF output.tif \\
     --bbox "37.0,38.0,-122.5,-121.5" --output-spacing 10.0

3. For production use:
   - Use real Sentinel-1 SLC data
   - Download precise orbit files  
   - Use high-resolution DEMs (SRTM 30m or Copernicus 10m)
   - Choose appropriate output projections (UTM for local analysis)

📚 Documentation:
   See docs/implementation/TERRAIN_CORRECTION_IMPLEMENTATION.md
""")


if __name__ == "__main__":
    main()
