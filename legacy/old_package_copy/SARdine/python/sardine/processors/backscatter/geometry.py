"""Geospatial helper utilities for the backscatter processor."""

from __future__ import annotations

import math
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np


WGS84_SEMI_MAJOR_AXIS_M = 6_378_137.0
WGS84_ECCENTRICITY_SQUARED = 6.694_379_990_14e-3


def meters_per_degree(latitude_deg: float) -> Tuple[float, float]:
    """Return meridional and prime-vertical arc lengths (meters per degree)."""

    lat_rad = math.radians(latitude_deg)
    sin_lat = math.sin(lat_rad)
    cos_lat = math.cos(lat_rad)

    denom = math.sqrt(1.0 - WGS84_ECCENTRICITY_SQUARED * sin_lat * sin_lat)
    if denom == 0.0:
        raise ValueError("Invalid latitude for WGS84 radius computation")

    prime_vertical = WGS84_SEMI_MAJOR_AXIS_M / denom
    meters_per_degree_lon = prime_vertical * cos_lat * math.pi / 180.0

    meridional_radius = (
        WGS84_SEMI_MAJOR_AXIS_M
        * (1.0 - WGS84_ECCENTRICITY_SQUARED)
        / (1.0 - WGS84_ECCENTRICITY_SQUARED * sin_lat * sin_lat) ** 1.5
    )
    meters_per_degree_lat = meridional_radius * math.pi / 180.0

    if meters_per_degree_lat <= 0.0 or meters_per_degree_lon <= 0.0:
        raise ValueError("Derived meters-per-degree values must be positive")

    return meters_per_degree_lat, meters_per_degree_lon


def refine_geocoding_bbox(
    bbox: Sequence[float],
    rows: int,
    cols: int,
    range_spacing_m: float,
    azimuth_spacing_m: float,
) -> Optional[Tuple[List[float], Dict[str, float]]]:
    """Refine geocoding bbox based on actual SAR image dimensions.
    
    The metadata footprint may be larger than the actual SAR coverage (especially
    for merged IW data where timing metadata covers the full acquisition window).
    This function computes the expected ground extent from image dimensions and
    pixel spacing, then shrinks the bbox if needed to match the actual coverage.
    """

    if rows <= 0 or cols <= 0:
        return None

    if not np.isfinite(range_spacing_m) or range_spacing_m <= 0.0:
        return None

    if not np.isfinite(azimuth_spacing_m) or azimuth_spacing_m <= 0.0:
        return None

    try:
        min_lon, min_lat, max_lon, max_lat = [float(value) for value in bbox]
    except (TypeError, ValueError):
        return None

    if min_lat >= max_lat or min_lon >= max_lon:
        return None

    original_lat_extent = max_lat - min_lat
    original_lon_extent = max_lon - min_lon

    # Compute expected ground extent from image dimensions and spacing
    center_lat = (min_lat + max_lat) / 2.0
    try:
        m_per_deg_lat, m_per_deg_lon = meters_per_degree(center_lat)
    except (ValueError, ZeroDivisionError):
        m_per_deg_lat = 111320.0  # fallback
        m_per_deg_lon = 111320.0 * math.cos(math.radians(center_lat))

    # Expected extent based on image size (azimuth = lat direction for descending pass)
    expected_azimuth_extent_m = rows * azimuth_spacing_m
    expected_range_extent_m = cols * range_spacing_m
    
    # Convert to degrees (approximate - SAR geometry is complex)
    expected_lat_extent_deg = expected_azimuth_extent_m / m_per_deg_lat
    expected_lon_extent_deg = expected_range_extent_m / m_per_deg_lon

    # Add generous margin (20%) for edge effects, interpolation, and SAR geometry
    # The SAR footprint doesn't map perfectly to a rectangular bbox
    margin_factor = 1.20
    expected_lat_extent_deg *= margin_factor
    expected_lon_extent_deg *= margin_factor

    # Only shrink if metadata bbox is MUCH larger than expected
    # (more than 50% larger indicates metadata footprint exceeds actual coverage significantly)
    # FIXED: Increased threshold from 1.2 to 1.5 to reduce false positives
    shrink_threshold = 1.5
    
    # FIXED: Maximum allowed shrink to prevent over-aggressive clipping
    # Never shrink more than 30% - if the difference is larger, the spacing metadata is likely wrong
    max_shrink_fraction = 0.30
    
    new_lat_extent = original_lat_extent
    new_lon_extent = original_lon_extent
    new_min_lat = min_lat
    new_max_lat = max_lat
    new_min_lon = min_lon
    new_max_lon = max_lon

    if original_lat_extent > expected_lat_extent_deg * shrink_threshold:
        # Metadata lat extent is too large - shrink it, but cap the shrink
        proposed_extent = expected_lat_extent_deg
        min_allowed_extent = original_lat_extent * (1.0 - max_shrink_fraction)
        new_lat_extent = max(proposed_extent, min_allowed_extent)
        lat_center = (min_lat + max_lat) / 2.0
        new_min_lat = lat_center - new_lat_extent / 2.0
        new_max_lat = lat_center + new_lat_extent / 2.0

    if original_lon_extent > expected_lon_extent_deg * shrink_threshold:
        # Metadata lon extent is too large - shrink it, but cap the shrink
        proposed_extent = expected_lon_extent_deg
        min_allowed_extent = original_lon_extent * (1.0 - max_shrink_fraction)
        new_lon_extent = max(proposed_extent, min_allowed_extent)
        lon_center = (min_lon + max_lon) / 2.0
        new_min_lon = lon_center - new_lon_extent / 2.0
        new_max_lon = lon_center + new_lon_extent / 2.0

    shrink_lat_pct = (1.0 - new_lat_extent / original_lat_extent) if original_lat_extent > 0 else 0.0
    shrink_lon_pct = (1.0 - new_lon_extent / original_lon_extent) if original_lon_extent > 0 else 0.0

    metrics: Dict[str, float] = {
        "original_lat_extent_deg": original_lat_extent,
        "original_lon_extent_deg": original_lon_extent,
        "expected_lat_extent_deg": expected_lat_extent_deg,
        "expected_lon_extent_deg": expected_lon_extent_deg,
        "new_lat_extent_deg": new_lat_extent,
        "new_lon_extent_deg": new_lon_extent,
        "shrink_lat_pct": shrink_lat_pct,
        "shrink_lon_pct": shrink_lon_pct,
    }

    refined_bbox: List[float] = [
        new_min_lon,
        new_min_lat,
        new_max_lon,
        new_max_lat,
    ]
    return refined_bbox, metrics


def compute_bbox_from_transform(
    transform: Sequence[float], 
    height: int, 
    width: int,
    epsg: Optional[int] = None
):
    """Derive bounding box from a GDAL-style transform.
    
    Args:
        transform: GDAL-style geotransform (6 elements)
        height: Image height in pixels
        width: Image width in pixels
        epsg: Optional EPSG code. If None, assumes geographic (EPSG:4326).
              If UTM (32601-32660 or 32701-32760), returns UTM coordinates.
              For geographic (4326), returns (min_lon, min_lat, max_lon, max_lat).
              For UTM, returns (min_easting, min_northing, max_easting, max_northing).
    
    Returns:
        Tuple of (min_x, min_y, max_x, max_y) where x/y are lon/lat for geographic
        or easting/northing for UTM. Returns None if transform is invalid.
    """

    if transform is None or len(transform) != 6:
        return None
    origin_x, pixel_x, rot_x, origin_y, rot_y, pixel_y = (
        float(transform[0]),
        float(transform[1]),
        float(transform[2]),
        float(transform[3]),
        float(transform[4]),
        float(transform[5]),
    )
    corner_x = origin_x + pixel_x * width + rot_x * height
    corner_y = origin_y + rot_y * width + pixel_y * height
    min_x = min(origin_x, corner_x)
    max_x = max(origin_x, corner_x)
    min_y = min(origin_y, corner_y)
    max_y = max(origin_y, corner_y)
    
    # For geographic coordinates, return as (min_lon, min_lat, max_lon, max_lat)
    # For UTM, return as (min_easting, min_northing, max_easting, max_northing)
    # Note: The function signature historically returns (min_lon, min_lat, max_lon, max_lat)
    # but we preserve this for backward compatibility even for UTM coordinates
    return (min_x, min_y, max_x, max_y)


def ensure_north_up_georeferencing(
    data: np.ndarray,
    transform: Sequence[float],
    metadata: Optional[Dict[str, float]] = None,
):
    """Ensure geotransform aligns with north-up, east-positive orientation."""

    if transform is None:
        raise ValueError("Terrain correction result did not include a geotransform")

    if not isinstance(data, np.ndarray) or data.ndim != 2:
        raise ValueError("Terrain correction data must be a 2D numpy array")

    transform_list = [float(value) for value in transform]
    if len(transform_list) != 6:
        raise ValueError("GeoTransform must contain exactly six elements")

    height, width = data.shape
    adjusted_data = data

    if transform_list[5] > 0:
        adjusted_data = np.flipud(adjusted_data)
        transform_list[3] += transform_list[5] * (height - 1)
        transform_list[5] = -transform_list[5]
        print("   🔃 Adjusted GeoTIFF orientation: flipped vertically for north-up alignment")

    if transform_list[1] < 0:
        adjusted_data = np.fliplr(adjusted_data)
        transform_list[0] += transform_list[1] * (width - 1)
        transform_list[1] = abs(transform_list[1])
        print("   🔃 Adjusted GeoTIFF orientation: flipped horizontally for east-positive alignment")

    if isinstance(metadata, dict):
        try:
            expected_min_lat = float(metadata.get('min_latitude'))
            expected_max_lat = float(metadata.get('max_latitude'))
            expected_min_lon = float(metadata.get('min_longitude'))
            expected_max_lon = float(metadata.get('max_longitude'))
        except (TypeError, ValueError):
            expected_min_lat = expected_max_lat = expected_min_lon = expected_max_lon = None

        lat_top = transform_list[3]
        lat_bottom = lat_top + transform_list[5] * (height - 1)
        lon_left = transform_list[0]
        lon_right = lon_left + transform_list[1] * (width - 1)

        if expected_min_lat is not None and expected_max_lat is not None:
            top_closer_to_min = abs(lat_top - expected_min_lat) < abs(lat_top - expected_max_lat)
            bottom_closer_to_max = abs(lat_bottom - expected_max_lat) < abs(lat_bottom - expected_min_lat)
            if top_closer_to_min and bottom_closer_to_max:
                pixel_height = transform_list[5]
                adjusted_data = np.flipud(adjusted_data)
                transform_list[3] = lat_bottom
                transform_list[5] = -abs(pixel_height)
                print("   🔃 Adjusted GeoTIFF orientation: flipped vertically to align with scene metadata")

    return transform_list, adjusted_data


def crop_geocoded_output(data: np.ndarray, transform: Sequence[float], bbox):
    """Trim large empty borders around the geocoded raster to boost valid coverage."""

    if not isinstance(data, np.ndarray) or data.size == 0:
        return None
    finite_mask = np.isfinite(data)
    valid_pixels = int(finite_mask.sum())
    if valid_pixels == 0:
        return None

    total_pixels = data.size
    valid_pct = 100.0 * valid_pixels / float(total_pixels)
    if valid_pct >= 75.0:
        return None

    rows_with_data = np.any(finite_mask, axis=1)
    cols_with_data = np.any(finite_mask, axis=0)
    if not rows_with_data.any() or not cols_with_data.any():
        return None

    row_indices = np.where(rows_with_data)[0]
    col_indices = np.where(cols_with_data)[0]
    row_start = int(row_indices[0])
    row_end = int(row_indices[-1]) + 1
    col_start = int(col_indices[0])
    col_end = int(col_indices[-1]) + 1

    if row_end - row_start < 2 or col_end - col_start < 2:
        return None

    cropped = data[row_start:row_end, col_start:col_end]
    cropped_mask = finite_mask[row_start:row_end, col_start:col_end]
    cropped_valid = int(cropped_mask.sum())
    cropped_pct = 100.0 * cropped_valid / float(cropped.size) if cropped.size else 0.0

    new_transform = list(transform)
    new_transform[0] = (
        float(transform[0])
        + float(transform[1]) * col_start
        + float(transform[2]) * row_start
    )
    new_transform[3] = (
        float(transform[3])
        + float(transform[4]) * col_start
        + float(transform[5]) * row_start
    )

    # Note: crop_geocoded_output doesn't have CRS info, so we pass None (assumes geographic)
    # This is acceptable since cropping typically happens before CRS-specific validation
    new_bbox = compute_bbox_from_transform(new_transform, cropped.shape[0], cropped.shape[1], epsg=None)
    if new_bbox is None:
        new_bbox = bbox

    return {
        "data": cropped,
        "transform": tuple(new_transform),
        "bbox": tuple(new_bbox) if new_bbox else bbox,
        "valid_pct": cropped_pct,
    }
