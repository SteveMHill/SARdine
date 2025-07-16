"""
GeoTIFF export utilities for SARdine.

Simple utilities for exporting processed SAR data to GeoTIFF format.
These functions use rasterio for compatibility and ease of use.
"""

import numpy as np
try:
    import rasterio
    from rasterio.transform import from_bounds
    from rasterio.crs import CRS
    RASTERIO_AVAILABLE = True
except ImportError:
    RASTERIO_AVAILABLE = False

def export_geotiff(data, output_path, bounds=None, crs='EPSG:4326', nodata=None, 
                   compress='lzw', tiled=True, description=None):
    """
    Export a 2D array to GeoTIFF format.
    
    Parameters
    ----------
    data : numpy.ndarray
        2D array to export
    output_path : str
        Output file path
    bounds : tuple, optional
        Bounding box as (west, south, east, north)
    crs : str, optional
        Coordinate reference system (default: 'EPSG:4326')
    nodata : float, optional
        NoData value
    compress : str, optional
        Compression method ('lzw', 'deflate', 'none')
    tiled : bool, optional
        Whether to create a tiled GeoTIFF
    description : str, optional
        Band description
        
    Returns
    -------
    str
        Path to created file
    """
    if not RASTERIO_AVAILABLE:
        raise ImportError("rasterio is required for GeoTIFF export. Install with: pip install rasterio")
    
    if data.ndim != 2:
        raise ValueError("Input data must be 2D")
    
    height, width = data.shape
    
    # Create transform from bounds or use identity
    if bounds is not None:
        west, south, east, north = bounds
        transform = from_bounds(west, south, east, north, width, height)
    else:
        transform = rasterio.transform.from_origin(0, height, 1, 1)
    
    # Setup profile
    profile = {
        'driver': 'GTiff',
        'dtype': data.dtype,
        'width': width,
        'height': height,
        'count': 1,
        'crs': CRS.from_string(crs),
        'transform': transform,
        'compress': compress,
        'tiled': tiled,
        'blockxsize': 512,
        'blockysize': 512,
    }
    
    if nodata is not None:
        profile['nodata'] = nodata
    
    # Write the GeoTIFF
    with rasterio.open(output_path, 'w', **profile) as dst:
        dst.write(data, 1)
        if description:
            dst.set_band_description(1, description)
    
    return output_path

def export_cog(data, output_path, bounds=None, crs='EPSG:4326', nodata=None,
               compress='lzw', description=None, overviews=True):
    """
    Export a 2D array to Cloud Optimized GeoTIFF (COG) format.
    
    Parameters
    ----------
    data : numpy.ndarray
        2D array to export
    output_path : str
        Output file path
    bounds : tuple, optional
        Bounding box as (west, south, east, north)
    crs : str, optional
        Coordinate reference system (default: 'EPSG:4326')
    nodata : float, optional
        NoData value
    compress : str, optional
        Compression method ('lzw', 'deflate', 'none')
    description : str, optional
        Band description
    overviews : bool, optional
        Whether to build overviews
        
    Returns
    -------
    str
        Path to created file
    """
    if not RASTERIO_AVAILABLE:
        raise ImportError("rasterio is required for COG export. Install with: pip install rasterio")
    
    if data.ndim != 2:
        raise ValueError("Input data must be 2D")
    
    height, width = data.shape
    
    # Create transform from bounds or use identity
    if bounds is not None:
        west, south, east, north = bounds
        transform = from_bounds(west, south, east, north, width, height)
    else:
        transform = rasterio.transform.from_origin(0, height, 1, 1)
    
    # COG-optimized profile
    profile = {
        'driver': 'GTiff',
        'dtype': data.dtype,
        'width': width,
        'height': height,
        'count': 1,
        'crs': CRS.from_string(crs),
        'transform': transform,
        'compress': compress,
        'tiled': True,
        'blockxsize': 512,
        'blockysize': 512,
        'interleave': 'pixel',
    }
    
    if nodata is not None:
        profile['nodata'] = nodata
    
    # Write the COG
    with rasterio.open(output_path, 'w', **profile) as dst:
        dst.write(data, 1)
        if description:
            dst.set_band_description(1, description)
        
        # Build overviews for COG compliance
        if overviews:
            overview_factors = [2, 4, 8, 16]
            dst.build_overviews(overview_factors, resampling=rasterio.enums.Resampling.average)
            dst.update_tags(ns='gdal', TILED='YES')
    
    return output_path

def export_multiband_geotiff(data_list, output_path, band_names=None, bounds=None, 
                            crs='EPSG:4326', nodata=None, compress='lzw', tiled=True):
    """
    Export multiple 2D arrays as a multi-band GeoTIFF.
    
    Parameters
    ----------
    data_list : list of numpy.ndarray
        List of 2D arrays (one per band)
    output_path : str
        Output file path
    band_names : list of str, optional
        Names for each band
    bounds : tuple, optional
        Bounding box as (west, south, east, north)
    crs : str, optional
        Coordinate reference system (default: 'EPSG:4326')
    nodata : float, optional
        NoData value
    compress : str, optional
        Compression method ('lzw', 'deflate', 'none')
    tiled : bool, optional
        Whether to create a tiled GeoTIFF
        
    Returns
    -------
    str
        Path to created file
    """
    if not RASTERIO_AVAILABLE:
        raise ImportError("rasterio is required for GeoTIFF export. Install with: pip install rasterio")
    
    if not data_list:
        raise ValueError("data_list cannot be empty")
    
    # Check all arrays have same dimensions
    first_shape = data_list[0].shape
    if not all(arr.shape == first_shape for arr in data_list):
        raise ValueError("All arrays must have the same dimensions")
    
    if any(arr.ndim != 2 for arr in data_list):
        raise ValueError("All arrays must be 2D")
    
    height, width = first_shape
    num_bands = len(data_list)
    
    # Create transform from bounds or use identity
    if bounds is not None:
        west, south, east, north = bounds
        transform = from_bounds(west, south, east, north, width, height)
    else:
        transform = rasterio.transform.from_origin(0, height, 1, 1)
    
    # Setup profile
    profile = {
        'driver': 'GTiff',
        'dtype': data_list[0].dtype,
        'width': width,
        'height': height,
        'count': num_bands,
        'crs': CRS.from_string(crs),
        'transform': transform,
        'compress': compress,
        'tiled': tiled,
        'blockxsize': 512,
        'blockysize': 512,
    }
    
    if nodata is not None:
        profile['nodata'] = nodata
    
    # Write the multi-band GeoTIFF
    with rasterio.open(output_path, 'w', **profile) as dst:
        for i, data in enumerate(data_list):
            dst.write(data, i + 1)
            
            # Set band description if provided
            if band_names and i < len(band_names):
                dst.set_band_description(i + 1, band_names[i])
    
    return output_path

def validate_geotiff(file_path):
    """
    Validate and get information about a GeoTIFF file.
    
    Parameters
    ----------
    file_path : str
        Path to GeoTIFF file
        
    Returns
    -------
    dict
        Information about the GeoTIFF file
    """
    if not RASTERIO_AVAILABLE:
        raise ImportError("rasterio is required for GeoTIFF validation. Install with: pip install rasterio")
    
    info = {}
    
    with rasterio.open(file_path) as src:
        info['width'] = src.width
        info['height'] = src.height
        info['count'] = src.count
        info['dtype'] = str(src.dtypes[0])
        info['crs'] = str(src.crs) if src.crs else None
        info['bounds'] = src.bounds
        info['transform'] = src.transform
        info['nodata'] = src.nodata
        
        # Band information
        bands = []
        for i in range(1, src.count + 1):
            band_info = {
                'band': i,
                'dtype': str(src.dtypes[i-1]),
                'nodata': src.nodatavals[i-1],
                'description': src.descriptions[i-1] or f'Band {i}'
            }
            bands.append(band_info)
        
        info['bands'] = bands
        
        # Check if it's a valid COG
        info['is_tiled'] = src.profile.get('tiled', False)
        info['blocksize'] = (src.profile.get('blockxsize', 0), src.profile.get('blockysize', 0))
        info['compression'] = src.profile.get('compress', 'none')
        
        # Overview information
        overviews = []
        for i in range(1, src.count + 1):
            band_overviews = src.overviews(i)
            if band_overviews:
                overviews.append({f'band_{i}': band_overviews})
        info['overviews'] = overviews
    
    return info
