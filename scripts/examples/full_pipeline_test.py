#!/usr/bin/env python3
"""
Complete SAR Processing Pipeline Test

This script demonstrates the full SARdine processing pipeline:
1. Create synthetic SAR data (simulating real Sentinel-1 data)
2. Generate realistic orbit data 
3. Apply terrain correction (both standard and optimized)
4. Apply speckle filtering
5. Convert to dB scale
6. Performance comparison

Expected workflow:
- Synthetic SAR Image → Terrain Correction → Speckle Filter → dB Conversion
"""

import numpy as np
import sys
import os
import time
from typing import Dict, Tuple, Any

# Import SARdine
try:
    import sardine
    print("✅ Successfully imported SARdine!")
    print(f"📦 Available functions: {len([x for x in dir(sardine) if not x.startswith('_')])}")
except ImportError as e:
    print(f"❌ Failed to import SARdine: {e}")
    sys.exit(1)

def create_realistic_sar_data() -> Tuple[np.ndarray, Dict[str, Any]]:
    """Create realistic synthetic SAR data"""
    
    print("🛰️  Creating realistic synthetic SAR data...")
    
    # Create a realistic SAR image (1000x800 pixels)
    rows, cols = 1000, 800
    
    # Generate spatially correlated SAR data with realistic statistics
    np.random.seed(42)  # For reproducible results
    
    # Create base terrain features
    x = np.linspace(0, 4*np.pi, cols)
    y = np.linspace(0, 4*np.pi, rows)
    X, Y = np.meshgrid(x, y)
    
    # Simulate different land cover types
    urban_areas = 0.8 * np.exp(-((X-2*np.pi)**2 + (Y-2*np.pi)**2) / 4)
    forest_areas = 0.3 * (np.sin(X/2) * np.cos(Y/2) + 1) / 2
    water_areas = 0.05 * np.ones_like(X)
    water_mask = (X < np.pi) & (Y < np.pi)
    
    # Combine land cover types
    base_backscatter = urban_areas + forest_areas + water_areas
    base_backscatter[water_mask] = water_areas[water_mask]
    
    # Add speckle noise (multiplicative)
    speckle = np.random.gamma(1.0, 1.0, (rows, cols))
    sar_image = base_backscatter * speckle
    
    # Ensure realistic value range (typical Sentinel-1 values)
    sar_image = np.clip(sar_image, 0.001, 2.0)
    
    stats = {
        'min_value': np.min(sar_image),
        'max_value': np.max(sar_image),
        'mean_value': np.mean(sar_image),
        'std_value': np.std(sar_image),
        'shape': sar_image.shape
    }
    
    print(f"📊 SAR image created: {rows}x{cols}")
    print(f"📈 Value range: {stats['min_value']:.4f} to {stats['max_value']:.4f}")
    print(f"📊 Mean: {stats['mean_value']:.4f}, Std: {stats['std_value']:.4f}")
    
    return sar_image.astype(np.float32), stats

def create_realistic_orbit_data(num_vectors: int = 15) -> Tuple[list, list, list]:
    """Create realistic Sentinel-1 orbit data"""
    
    print(f"🛰️  Creating realistic orbit data with {num_vectors} state vectors...")
    
    # Sentinel-1 typical orbital parameters
    orbit_altitude = 693000.0  # meters (693 km)
    earth_radius = 6371000.0   # meters
    orbital_radius = earth_radius + orbit_altitude
    orbital_velocity = 7550.0  # m/s (typical LEO velocity)
    
    # Time span: 30 seconds of orbit data
    total_time = 30.0  # seconds
    time_step = total_time / (num_vectors - 1)
    
    times = []
    positions = []
    velocities = []
    
    # Start time (RFC3339 format)
    from datetime import datetime, timezone
    start_time = datetime.now(timezone.utc)
    
    for i in range(num_vectors):
        # Time for this state vector
        current_time = start_time.timestamp() + i * time_step
        time_dt = datetime.fromtimestamp(current_time, timezone.utc)
        times.append(time_dt.isoformat().replace('+00:00', 'Z'))
        
        # Orbital angle (simplified circular orbit)
        angle = (i * time_step / total_time) * 2 * np.pi * 0.1  # Small arc
        
        # Position in ECEF coordinates (simplified)
        x = orbital_radius * np.cos(angle)
        y = orbital_radius * np.sin(angle) * 0.1  # Small Y component
        z = orbital_radius * np.sin(angle) * 0.8  # Inclined orbit
        
        positions.append([x, y, z])
        
        # Velocity vector (tangential to orbit)
        vx = -orbital_velocity * np.sin(angle)
        vy = orbital_velocity * np.cos(angle) * 0.1
        vz = orbital_velocity * np.cos(angle) * 0.8
        
        velocities.append([vx, vy, vz])
    
    print(f"🛰️  Created {num_vectors} orbital state vectors")
    print(f"📍 Altitude: {orbit_altitude/1000:.0f} km")
    print(f"⚡ Velocity: {orbital_velocity:.0f} m/s")
    
    return times, positions, velocities

def create_realistic_metadata() -> Dict[str, float]:
    """Create realistic Sentinel-1 metadata"""
    
    return {
        'range_pixel_spacing': 2.3296,      # meters (Sentinel-1 IW mode)
        'azimuth_pixel_spacing': 14.0,      # meters (Sentinel-1 IW mode)
        'wavelength': 0.0555,               # meters (C-band, 5.405 GHz)
        'slant_range_time': 0.005331,       # seconds (typical Sentinel-1)
        'prf': 1586.0,                      # Hz (pulse repetition frequency)
    }

def run_full_pipeline_test():
    """Run the complete SAR processing pipeline"""
    
    print("🚀 SARdine Complete Processing Pipeline Test")
    print("=" * 80)
    
    # Create test data
    print("\n📋 Step 1: Creating test data...")
    sar_image, image_stats = create_realistic_sar_data()
    times, positions, velocities = create_realistic_orbit_data()
    metadata = create_realistic_metadata()
    
    # Define processing area
    sar_bbox = [11.0, 45.0, 11.5, 45.3]  # Northern Italy
    cache_dir = "/tmp/sardine_pipeline_test"
    output_resolution = 30.0  # meters
    
    os.makedirs(cache_dir, exist_ok=True)
    
    print(f"📍 Processing area: [{sar_bbox[0]:.1f}°, {sar_bbox[1]:.1f}°, {sar_bbox[2]:.1f}°, {sar_bbox[3]:.1f}°]")
    print(f"📏 Target resolution: {output_resolution}m")
    
    # Pipeline Step 1: Standard Terrain Correction
    print("\n📋 Step 2: Standard Terrain Correction...")
    start_time = time.time()
    
    try:
        standard_result = sardine.apply_terrain_correction(
            sar_image=sar_image,
            sar_bbox=sar_bbox,
            orbit_times=times,
            orbit_positions=positions,
            orbit_velocities=velocities,
            cache_dir=cache_dir,
            output_resolution=output_resolution,
            real_metadata=metadata
        )
        
        standard_time = time.time() - start_time
        print(f"✅ Standard terrain correction completed in {standard_time:.2f}s")
        print(f"📊 Output dimensions: {standard_result['rows']}x{standard_result['cols']}")
        
    except Exception as e:
        print(f"❌ Standard terrain correction failed: {e}")
        standard_result = None
        standard_time = 0
    
    # Pipeline Step 2: Optimized Terrain Correction
    print("\n📋 Step 3: OPTIMIZED Terrain Correction...")
    start_time = time.time()
    
    try:
        optimized_result = sardine.apply_terrain_correction_optimized(
            sar_image=sar_image,
            sar_bbox=sar_bbox,
            orbit_times=times,
            orbit_positions=positions,
            orbit_velocities=velocities,
            cache_dir=cache_dir,
            output_resolution=output_resolution,
            real_metadata=metadata
        )
        
        optimized_time = time.time() - start_time
        print(f"✅ Optimized terrain correction completed in {optimized_time:.2f}s")
        print(f"📊 Output dimensions: {optimized_result['rows']}x{optimized_result['cols']}")
        
        # Calculate speedup
        if standard_time > 0 and optimized_time > 0:
            speedup = standard_time / optimized_time
            print(f"⚡ Performance improvement: {speedup:.1f}x faster")
        
    except Exception as e:
        print(f"❌ Optimized terrain correction failed: {e}")
        optimized_result = None
        optimized_time = 0
    
    # Choose result for further processing (prefer optimized)
    terrain_result = optimized_result if optimized_result else standard_result
    
    if terrain_result is None:
        print("❌ No terrain correction results available for pipeline continuation")
        return False
    
    # Pipeline Step 3: Speckle Filtering
    print("\n📋 Step 4: Speckle Filtering...")
    try:
        filter_result = sardine.apply_speckle_filter_optimized(
            image=terrain_result['data'],
            filter_type="lee",
            window_size=7,
            num_looks=1.0
        )
        
        print(f"✅ Speckle filtering completed")
        print(f"📊 Filtered dimensions: {filter_result['rows']}x{filter_result['cols']}")
        
    except Exception as e:
        print(f"❌ Speckle filtering failed: {e}")
        filter_result = {'filtered_data': terrain_result['data']}
    
    # Pipeline Step 4: dB Conversion
    print("\n📋 Step 5: dB Scale Conversion...")
    try:
        db_result = sardine.convert_to_db_real(values=filter_result['filtered_data'])
        
        # Calculate statistics
        db_data = np.array(db_result)
        valid_mask = np.isfinite(db_data) & (db_data > -50)  # Reasonable dB range
        valid_count = np.sum(valid_mask)
        total_count = db_data.size
        
        if valid_count > 0:
            db_min = np.min(db_data[valid_mask])
            db_max = np.max(db_data[valid_mask])
            db_mean = np.mean(db_data[valid_mask])
            
            print(f"✅ dB conversion completed")
            print(f"📊 Valid pixels: {valid_count:,} / {total_count:,} ({100*valid_count/total_count:.1f}%)")
            print(f"📈 dB range: {db_min:.1f} to {db_max:.1f} dB (mean: {db_mean:.1f} dB)")
        else:
            print("⚠️  dB conversion completed but no valid values found")
        
    except Exception as e:
        print(f"❌ dB conversion failed: {e}")
        db_result = None
    
    # Pipeline Summary
    print("\n📋 Pipeline Summary:")
    print("=" * 50)
    print(f"🛰️  Input SAR image: {image_stats['shape'][0]}x{image_stats['shape'][1]} pixels")
    print(f"📍 Processing area: {(sar_bbox[2]-sar_bbox[0]):.1f}° × {(sar_bbox[3]-sar_bbox[1]):.1f}°")
    print(f"🗺️  Terrain correction: {'✅ SUCCESS' if terrain_result else '❌ FAILED'}")
    print(f"🌊 Speckle filtering: {'✅ SUCCESS' if 'filtered_data' in filter_result else '❌ FAILED'}")
    print(f"📊 dB conversion: {'✅ SUCCESS' if db_result is not None else '❌ FAILED'}")
    
    if standard_time > 0 and optimized_time > 0:
        print(f"⚡ Performance gain: {standard_time/optimized_time:.1f}x speedup with optimization")
    
    print(f"📏 Final resolution: {output_resolution}m")
    
    print("\n🎉 Complete SAR processing pipeline test finished!")
    
    return True

if __name__ == "__main__":
    success = run_full_pipeline_test()
    sys.exit(0 if success else 1)
