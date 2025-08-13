"""
GeoTIFF export and STAC metadata generation for SARdine
"""

import json
import os
import uuid
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, Any, Optional, Tuple, Union
import numpy as np


def export_to_geotiff(
    data: np.ndarray,
    output_path: Union[str, Path],
    geotransform: Tuple[float, float, float, float, float, float],
    crs: str = "EPSG:4326",
    nodata_value: Optional[float] = None,
    compress: str = "lzw",
    tiled: bool = True,
    overview_levels: Optional[list] = None,
    metadata: Optional[Dict[str, Any]] = None
) -> Path:
    """
    Export numpy array to Cloud Optimized GeoTIFF
    
    Args:
        data: 2D numpy array to export
        output_path: Path for output GeoTIFF file
        geotransform: GDAL geotransform (ulx, xres, xskew, uly, yskew, yres)
        crs: Coordinate reference system (default: EPSG:4326)
        nodata_value: Value to use for nodata pixels
        compress: Compression method (lzw, deflate, jpeg, etc.)
        tiled: Whether to create tiled GeoTIFF
        overview_levels: Overview levels for COG (default: [2, 4, 8, 16])
        metadata: Additional metadata to write to GeoTIFF
        
    Returns:
        Path to created GeoTIFF file
    """
    try:
        import rasterio
        from rasterio.crs import CRS
        from rasterio.enums import Resampling
        from rasterio.transform import from_bounds
    except ImportError:
        raise ImportError("rasterio package required for GeoTIFF export. Install with: pip install rasterio")
    
    output_path = Path(output_path)
    
    # Ensure data is 2D
    if data.ndim != 2:
        raise ValueError(f"Data must be 2D array, got {data.ndim}D")
    
    # Set default overview levels
    if overview_levels is None:
        overview_levels = [2, 4, 8, 16]
    
    # Handle nodata
    if nodata_value is None:
        if data.dtype == np.float32 or data.dtype == np.float64:
            nodata_value = np.nan
        else:
            nodata_value = 0
    
    # Create transform from geotransform
    transform = rasterio.Affine.from_gdal(*geotransform)
    
    # COG creation profile
    profile = {
        'driver': 'GTiff',
        'height': data.shape[0],
        'width': data.shape[1],
        'count': 1,
        'dtype': data.dtype,
        'crs': CRS.from_string(crs),
        'transform': transform,
        'compress': compress,
        'tiled': tiled,
        'blockxsize': 512,
        'blockysize': 512,
        'interleave': 'band',
    }
    
    if nodata_value is not None and not np.isnan(nodata_value):
        profile['nodata'] = nodata_value
    
    # Write the GeoTIFF
    with rasterio.open(output_path, 'w', **profile) as dst:
        dst.write(data, 1)
        
        # Add metadata
        if metadata:
            dst.update_tags(**metadata)
        
        # Build overviews for COG
        if overview_levels:
            dst.build_overviews(overview_levels, Resampling.average)
            dst.update_tags(ns='rio_overview', resampling='average')
    
    print(f"✅ GeoTIFF exported: {output_path}")
    print(f"   📊 Size: {data.shape[1]} x {data.shape[0]}")
    print(f"   🗺️  CRS: {crs}")
    print(f"   💾 Compression: {compress}")
    print(f"   📈 Overviews: {overview_levels}")
    
    return output_path


def generate_stac_metadata(
    geotiff_path: Union[str, Path],
    sar_metadata: Dict[str, Any],
    output_path: Optional[Union[str, Path]] = None,
    collection_id: str = "sardine-backscatter",
    description: Optional[str] = None
) -> Path:
    """
    Generate STAC (SpatioTemporal Asset Catalog) metadata for SAR backscatter product
    
    Args:
        geotiff_path: Path to the GeoTIFF file
        sar_metadata: SAR processing metadata from SARdine
        output_path: Path for output STAC JSON (default: same as GeoTIFF with .json)
        collection_id: STAC collection identifier
        description: Product description
        
    Returns:
        Path to created STAC JSON file
    """
    try:
        import rasterio
        from rasterio.crs import CRS
    except ImportError:
        raise ImportError("rasterio package required for STAC metadata. Install with: pip install rasterio")
    
    geotiff_path = Path(geotiff_path)
    
    if output_path is None:
        output_path = geotiff_path.with_suffix('.json')
    else:
        output_path = Path(output_path)
    
    # Read GeoTIFF metadata
    with rasterio.open(geotiff_path) as src:
        bounds = src.bounds
        crs = src.crs
        transform = src.transform
        shape = src.shape
        
        # Convert bounds to WGS84 if needed
        if crs != CRS.from_epsg(4326):
            from rasterio.warp import transform_bounds
            bounds = transform_bounds(crs, CRS.from_epsg(4326), *bounds)
    
    # Extract key information from SAR metadata
    platform = sar_metadata.get('platform', 'unknown')
    acquisition_time = sar_metadata.get('acquisition_start_time', datetime.now(timezone.utc).isoformat())
    polarization = sar_metadata.get('polarization', 'unknown')
    orbit_direction = sar_metadata.get('orbit_direction', 'unknown')
    processing_level = sar_metadata.get('processing_level', 'L1')
    
    # Create unique ID
    item_id = f"{platform.lower()}_{Path(geotiff_path).stem}_{uuid.uuid4().hex[:8]}"
    
    # Default description
    if description is None:
        description = f"SAR backscatter product from {platform} {polarization} polarization"
    
    # STAC Item structure (v1.0.0)
    stac_item = {
        "stac_version": "1.0.0",
        "stac_extensions": [
            "https://stac-extensions.github.io/sar/v1.0.0/schema.json",
            "https://stac-extensions.github.io/projection/v1.0.0/schema.json"
        ],
        "type": "Feature",
        "id": item_id,
        "collection": collection_id,
        "geometry": {
            "type": "Polygon",
            "coordinates": [[
                [bounds[0], bounds[1]],  # min_x, min_y
                [bounds[2], bounds[1]],  # max_x, min_y
                [bounds[2], bounds[3]],  # max_x, max_y
                [bounds[0], bounds[3]],  # min_x, max_y
                [bounds[0], bounds[1]]   # close polygon
            ]]
        },
        "bbox": [bounds[0], bounds[1], bounds[2], bounds[3]],
        "properties": {
            "datetime": acquisition_time,
            "title": f"SAR Backscatter - {platform} {polarization}",
            "description": description,
            "platform": platform.lower(),
            "instruments": ["c-sar"] if platform.upper().startswith('S1') else ["sar"],
            "constellation": "sentinel-1" if platform.upper().startswith('S1') else platform.lower(),
            "mission": "sentinel-1" if platform.upper().startswith('S1') else platform.lower(),
            "processing:level": processing_level,
            "processing:software": {
                "sardine": {
                    "version": "0.2.0",
                    "repository": "https://github.com/SteveMHill/SARdine"
                }
            },
            # SAR extension properties
            "sar:instrument_mode": sar_metadata.get('acquisition_mode', 'IW'),
            "sar:frequency_band": "C",
            "sar:polarizations": [polarization] if polarization != 'unknown' else [],
            "sar:product_type": "backscatter",
            "sar:resolution_range": float(sar_metadata.get('range_pixel_spacing', 0)),
            "sar:resolution_azimuth": float(sar_metadata.get('azimuth_pixel_spacing', 0)),
            "sar:pixel_spacing_range": float(sar_metadata.get('range_pixel_spacing', 0)),
            "sar:pixel_spacing_azimuth": float(sar_metadata.get('azimuth_pixel_spacing', 0)),
            "sar:orbit_state": orbit_direction.lower() if orbit_direction != 'unknown' else None,
            # Projection extension properties
            "proj:epsg": int(crs.to_epsg()) if crs.to_epsg() else None,
            "proj:wkt2": crs.to_wkt() if crs else None,
            "proj:transform": list(transform)[:6],
            "proj:shape": [shape[0], shape[1]]
        },
        "assets": {
            "backscatter": {
                "href": f"./{geotiff_path.name}",
                "type": "image/tiff; application=geotiff; profile=cloud-optimized",
                "title": "SAR Backscatter",
                "description": f"Terrain-corrected SAR backscatter in {polarization} polarization",
                "roles": ["data"],
                "raster:bands": [{
                    "name": f"backscatter_{polarization}",
                    "description": f"SAR backscatter coefficient {polarization} polarization",
                    "data_type": "float32",
                    "unit": "dB",
                    "nodata": "nan"
                }]
            }
        },
        "links": [
            {
                "rel": "self",
                "type": "application/json",
                "href": f"./{output_path.name}"
            },
            {
                "rel": "collection",
                "type": "application/json",
                "href": f"./collection.json"
            }
        ]
    }
    
    # Add additional SAR metadata if available
    if 'orbit_number' in sar_metadata:
        stac_item['properties']['sar:absolute_orbit'] = int(sar_metadata['orbit_number'])
    
    if 'relative_orbit_number' in sar_metadata:
        stac_item['properties']['sar:relative_orbit'] = int(sar_metadata['relative_orbit_number'])
    
    if 'look_direction' in sar_metadata:
        stac_item['properties']['sar:looks_equivalent_number'] = float(sar_metadata.get('multilook_range', 1)) * float(sar_metadata.get('multilook_azimuth', 1))
    
    if 'incidence_angle_min' in sar_metadata and 'incidence_angle_max' in sar_metadata:
        stac_item['properties']['sar:incidence_angle'] = {
            "minimum": float(sar_metadata['incidence_angle_min']),
            "maximum": float(sar_metadata['incidence_angle_max'])
        }
    
    # Write STAC JSON
    with open(output_path, 'w') as f:
        json.dump(stac_item, f, indent=2)
    
    print(f"✅ STAC metadata generated: {output_path}")
    print(f"   🆔 Item ID: {item_id}")
    print(f"   🗂️  Collection: {collection_id}")
    print(f"   🌍 Extent: {bounds}")
    print(f"   📡 Platform: {platform}")
    print(f"   📊 Polarization: {polarization}")
    
    return output_path


def create_cog_with_stac(
    data: np.ndarray,
    output_dir: Union[str, Path],
    filename_base: str,
    geotransform: Tuple[float, float, float, float, float, float],
    sar_metadata: Dict[str, Any],
    crs: str = "EPSG:4326",
    compress: str = "lzw"
) -> Tuple[Path, Path]:
    """
    Create Cloud Optimized GeoTIFF with STAC metadata
    
    Args:
        data: 2D numpy array to export
        output_dir: Output directory
        filename_base: Base filename (without extension)
        geotransform: GDAL geotransform
        sar_metadata: SAR processing metadata
        crs: Coordinate reference system
        compress: Compression method
        
    Returns:
        Tuple of (geotiff_path, stac_path)
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Create GeoTIFF
    geotiff_path = output_dir / f"{filename_base}.tif"
    
    # Add processing metadata to GeoTIFF tags
    processing_metadata = {
        'SARDINE_VERSION': '0.2.0',
        'PROCESSING_DATE': datetime.now(timezone.utc).isoformat(),
        'PLATFORM': sar_metadata.get('platform', 'unknown'),
        'POLARIZATION': sar_metadata.get('polarization', 'unknown'),
        'PROCESSING_LEVEL': sar_metadata.get('processing_level', 'L1'),
        'ORBIT_DIRECTION': sar_metadata.get('orbit_direction', 'unknown')
    }
    
    export_to_geotiff(
        data=data,
        output_path=geotiff_path,
        geotransform=geotransform,
        crs=crs,
        compress=compress,
        metadata=processing_metadata
    )
    
    # Create STAC metadata
    stac_path = generate_stac_metadata(
        geotiff_path=geotiff_path,
        sar_metadata=sar_metadata,
        description=f"SAR backscatter product processed with SARdine v0.2.0"
    )
    
    return geotiff_path, stac_path
