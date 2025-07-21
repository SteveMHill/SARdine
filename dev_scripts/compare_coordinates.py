#!/usr/bin/env python3
"""
Compare coordinate extraction methods
"""

import sys
import os
import numpy as np
from pathlib import Path

# Add SARdine Python bindings to path
sys.path.insert(0, "/home/datacube/SARdine/python")

import sardine

def compare_coordinate_extraction():
    """Compare different ways of extracting scene coordinates"""
    
    slc_path = "/home/datacube/SARdine/data/S1A_IW_SLC__1SDV_20200103T170815_20200103T170842_030639_0382D5_DADE.zip"
    
    print("=== Official Manifest Footprint ===")
    print("Coordinates: 51.492371,8.047090 51.896465,11.731247 50.275928,12.124238 49.874306,8.568501")
    print("This gives a bounding box:")
    
    # Parse the official coordinates
    coords = [
        (51.492371, 8.047090),
        (51.896465, 11.731247), 
        (50.275928, 12.124238),
        (49.874306, 8.568501)
    ]
    lats = [c[0] for c in coords]
    lons = [c[1] for c in coords]
    official_bbox = [min(lons), min(lats), max(lons), max(lats)]  # [W, S, E, N]
    
    print(f"Official bbox [W,S,E,N]: {official_bbox}")
    print(f"Longitude range: {min(lons):.6f} to {max(lons):.6f} ({max(lons)-min(lons):.6f}°)")
    print(f"Latitude range: {min(lats):.6f} to {max(lats):.6f} ({max(lats)-min(lats):.6f}°)")
    
    print("\n=== SARdine Product Info ===")
    try:
        product_info = sardine.get_product_info(slc_path)
        print(f"Product info: {product_info}")
        
        # Check if product_info has bbox information
        if hasattr(product_info, 'bbox'):
            sardine_bbox = product_info.bbox
            print(f"SARdine bbox: {sardine_bbox}")
        elif hasattr(product_info, 'footprint'):
            sardine_footprint = product_info.footprint
            print(f"SARdine footprint: {sardine_footprint}")
        elif hasattr(product_info, 'coordinates'):
            sardine_coords = product_info.coordinates
            print(f"SARdine coordinates: {sardine_coords}")
        else:
            print("No coordinate information found in product_info")
            # List all attributes
            attrs = [attr for attr in dir(product_info) if not attr.startswith('_')]
            print(f"Available attributes: {attrs}")
            
    except Exception as e:
        print(f"Error getting product info: {e}")
    
    print("\n=== Our XML Annotation Extraction ===")
    # This is what our pipeline currently does
    import zipfile
    import xml.etree.ElementTree as ET
    
    with zipfile.ZipFile(slc_path, 'r') as z:
        annotation_files = [f for f in z.namelist() 
                           if 'annotation/' in f and f.endswith('.xml') 
                           and 'calibration' not in f and 'noise' not in f]
        
        print(f"Found {len(annotation_files)} annotation files:")
        for ann_file in annotation_files:
            print(f"  {ann_file}")
        
        # Test with the first file
        ann_file = annotation_files[0]
        print(f"\nExtracting from: {ann_file}")
        
        with z.open(ann_file) as xml_file:
            content = xml_file.read().decode('utf-8')
        
        root = ET.fromstring(content)
        
        coordinates = []
        
        # Find all latitude/longitude pairs in geolocationGrid
        for elem in root.iter():
            if 'geolocationGrid' in elem.tag.lower() or 'geolocation' in elem.tag.lower():
                for child in elem.iter():
                    if child.tag == 'latitude' or child.tag == 'longitude':
                        try:
                            coord_value = float(child.text)
                            coordinates.append((child.tag, coord_value))
                        except (ValueError, TypeError):
                            pass
        
        lats = [val for tag, val in coordinates if tag == 'latitude']
        lons = [val for tag, val in coordinates if tag == 'longitude']
        
        if lats and lons:
            annotation_bbox = [min(lons), min(lats), max(lons), max(lats)]
            print(f"Annotation bbox [W,S,E,N]: {annotation_bbox}")
            print(f"Longitude range: {min(lons):.6f} to {max(lons):.6f} ({max(lons)-min(lons):.6f}°)")
            print(f"Latitude range: {min(lats):.6f} to {max(lats):.6f} ({max(lats)-min(lats):.6f}°)")
            
            # Compare with official
            lon_diff = abs((max(lons) - min(lons)) - (official_bbox[2] - official_bbox[0]))
            lat_diff = abs((max(lats) - min(lats)) - (official_bbox[3] - official_bbox[1]))
            
            print(f"\nComparison with official:")
            print(f"  Longitude extent difference: {lon_diff:.6f}° ({'LARGE DIFFERENCE' if lon_diff > 1.0 else 'small difference'})")
            print(f"  Latitude extent difference: {lat_diff:.6f}° ({'LARGE DIFFERENCE' if lat_diff > 0.5 else 'small difference'})")
            
            if lon_diff > 1.0:
                print("  ⚠️  ANNOTATION EXTRACTION IS MISSING SIGNIFICANT LONGITUDE COVERAGE")
                print("  This explains why DEM and SAR scene don't overlap!")
    
    print("\n=== DEM with Official Coordinates ===")
    try:
        print("Testing DEM preparation with official coordinates...")
        official_bbox_tuple = tuple(official_bbox)
        dem_result = sardine.prepare_dem_for_scene(
            bbox=official_bbox_tuple,
            cache_dir='/tmp/sardine_dem_cache',
            output_resolution=20.0
        )
        
        if isinstance(dem_result, tuple) and len(dem_result) == 2:
            dem_data, dem_geotransform = dem_result
            print(f"DEM prepared successfully with official coordinates:")
            print(f"  DEM shape: {dem_data.shape}")
            print(f"  DEM geotransform: {dem_geotransform}")
            
            # Check DEM coverage
            dem_min_x = dem_geotransform[0]
            dem_max_x = dem_geotransform[0] + dem_data.shape[1] * dem_geotransform[1]
            dem_max_y = dem_geotransform[3]
            dem_min_y = dem_geotransform[3] + dem_data.shape[0] * dem_geotransform[5]
            
            print(f"  DEM coverage: X={dem_min_x:.6f} to {dem_max_x:.6f}, Y={dem_min_y:.6f} to {dem_max_y:.6f}")
            
            # Check overlap with official SAR scene
            overlap_x = max(0, min(dem_max_x, official_bbox[2]) - max(dem_min_x, official_bbox[0]))
            overlap_y = max(0, min(dem_max_y, official_bbox[3]) - max(dem_min_y, official_bbox[1]))
            
            if overlap_x > 0 and overlap_y > 0:
                print(f"  ✅ DEM overlaps with official SAR scene: {overlap_x:.6f}° x {overlap_y:.6f}°")
            else:
                print(f"  ❌ Still no overlap with official coordinates")
                
    except Exception as e:
        print(f"DEM preparation with official coordinates failed: {e}")

if __name__ == "__main__":
    compare_coordinate_extraction()
