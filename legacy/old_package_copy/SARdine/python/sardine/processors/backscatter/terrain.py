"""
Terrain flattening and terrain correction helpers for backscatter processing.
"""

import logging
import time
from pathlib import Path
from typing import Optional
import numpy as np
import sardine

from .metadata_normalize import (
    normalize_metadata_for_rust,
    validate_metadata_completeness,
    validate_orbit_vector_count,
    MINIMUM_ORBIT_VECTORS,
)
from .geometry import (
    compute_bbox_from_transform, 
    crop_geocoded_output, 
    refine_geocoding_bbox,
    ensure_north_up_georeferencing,
)

logger = logging.getLogger(__name__)

# Default DEM resolution in meters when not specified
DEFAULT_DEM_RESOLUTION_M = 30.0
# Earth constants for coordinate conversion
METERS_PER_DEG_LAT_BASE = 111_132.954
METERS_PER_DEG_LON_BASE = 111_320.0

# Valid ranges for geographic coordinates
VALID_LAT_RANGE = (-90.0, 90.0)
VALID_LON_RANGE = (-180.0, 180.0)

# Default multilook factor when not provided
DEFAULT_MULTILOOK_FACTOR = 1.0

# Coverage validation thresholds for terrain correction output
# These prevent silent failures where geocoding produces no/few valid pixels
# For IW TOPSAR products, the SAR swath is a parallelogram projected onto a
# rectangular output grid, so ~55-65% fill is normal. The critical threshold
# catches genuine failures (wrong timing, missing DEM, corrupt orbit).
CRITICAL_COVERAGE_THRESHOLD = 45.0  # Below 45%: fail hard (likely DEM gap or timing error)
WARNING_COVERAGE_THRESHOLD = 75.0  # Below 75%: warn loudly (partial failure)
DEM_VOID_WARNING_THRESHOLD = 20.0  # Above 20% voids: warn about DEM quality


def _validate_bbox(
    bbox: tuple, 
    source: str = "unknown",
    epsg: Optional[int] = None
) -> tuple[bool, str]:
    """Validate a geographic bounding box.
    
    Args:
        bbox: Bounding box tuple (min_lon, min_lat, max_lon, max_lat)
        source: Description of where the bbox came from (for error messages)
        
    Returns:
        Tuple of (is_valid, error_message)
    """
    if bbox is None:
        return False, f"Bounding box from {source} is None"
    
    if not isinstance(bbox, (tuple, list)):
        return False, f"Bounding box from {source} is not a tuple/list: {type(bbox)}"
    
    if len(bbox) != 4:
        return False, f"Bounding box from {source} has {len(bbox)} elements, expected 4"
    
    try:
        min_lon, min_lat, max_lon, max_lat = [float(v) for v in bbox]
    except (TypeError, ValueError) as e:
        return False, f"Bounding box from {source} contains non-numeric values: {e}"
    
    # Check for finite values
    if not all(np.isfinite([min_lon, min_lat, max_lon, max_lat])):
        return False, f"Bounding box from {source} contains non-finite values: {bbox}"
    
    # Determine if coordinates are geographic or projected (UTM)
    is_geographic = epsg is None or epsg == 4326
    is_utm = epsg is not None and ((32601 <= epsg <= 32660) or (32701 <= epsg <= 32760))
    
    # Only validate lat/lon ranges for geographic coordinates
    # For UTM coordinates, skip lat/lon validation (coordinates are in meters)
    if is_geographic:
        # Validate latitude range
        if min_lat < VALID_LAT_RANGE[0] or max_lat > VALID_LAT_RANGE[1]:
            return False, (
                f"Bounding box from {source} has invalid latitude: "
                f"[{min_lat}, {max_lat}] not in {VALID_LAT_RANGE}"
            )
        
        # Validate longitude range
        if min_lon < VALID_LON_RANGE[0] or max_lon > VALID_LON_RANGE[1]:
            return False, (
                f"Bounding box from {source} has invalid longitude: "
                f"[{min_lon}, {max_lon}] not in {VALID_LON_RANGE}"
            )
    elif is_utm:
        # For UTM, validate that coordinates are reasonable easting/northing values
        # UTM easting: typically 166000 to 834000 (with false easting of 500000)
        # UTM northing: 0 to 10000000 (Northern) or 0 to 10000000 (Southern, but with false northing)
        # We'll just check they're positive and reasonable (not negative, not too large)
        if min_lon < 0 or max_lon > 10000000 or min_lat < 0 or max_lat > 10000000:
            # This is a warning, not an error - UTM coordinates can vary
            # We'll still accept them but log a warning
            pass  # Accept UTM coordinates even if outside typical ranges
    
    # Check min < max
    if min_lat >= max_lat:
        return False, (
            f"Bounding box from {source} has min_lat >= max_lat: "
            f"{min_lat} >= {max_lat}"
        )
    
    if min_lon >= max_lon:
        return False, (
            f"Bounding box from {source} has min_lon >= max_lon: "
            f"{min_lon} >= {max_lon}"
        )
    
    return True, ""


def _apply_north_up_orientation(
    data: np.ndarray,
    transform: list,
) -> tuple[np.ndarray, list, Optional[str]]:
    """Apply north-up, east-positive orientation to geocoded data.
    
    Standard GeoTIFF convention expects:
    - transform[5] < 0 (negative y-spacing for north-up)
    - transform[1] > 0 (positive x-spacing for east-positive)
    
    If the geotransform violates these conventions, flip the data and adjust
    the transform accordingly.
    
    Args:
        data: 2D geocoded data array
        transform: 6-element GDAL geotransform
        
    Returns:
        Tuple of (adjusted_data, adjusted_transform, message)
        message is None if no adjustment was needed
    """
    if transform is None or len(transform) != 6:
        return data, transform, None
    
    transform = list(transform)  # Make a copy
    adjusted_data = data
    messages = []
    height, width = data.shape
    
    # Check for south-up orientation (positive y-spacing)
    if transform[5] > 0:
        adjusted_data = np.flipud(adjusted_data)
        # Adjust origin to top-left of flipped image
        transform[3] = transform[3] + transform[5] * (height - 1)
        transform[5] = -transform[5]
        messages.append("🔃 Flipped vertically for north-up alignment")
    
    # Check for west-positive orientation (negative x-spacing)
    if transform[1] < 0:
        adjusted_data = np.fliplr(adjusted_data)
        # Adjust origin to left side of flipped image
        transform[0] = transform[0] + transform[1] * (width - 1)
        transform[1] = abs(transform[1])
        messages.append("🔃 Flipped horizontally for east-positive alignment")
    
    if messages:
        return adjusted_data, transform, "; ".join(messages)
    return adjusted_data, transform, None


def _validate_dem_voids(dem_data: np.ndarray, stage: str = "DEM loading") -> None:
    """Validate DEM data for voids/missing values.
    
    Args:
        dem_data: DEM elevation data array
        stage: Description of where this is called from (for error messages)
        
    Raises:
        RuntimeError: If DEM has excessive voids (>50%)
    """
    if dem_data is None or dem_data.size == 0:
        raise RuntimeError(f"DEM data is empty or None during {stage}")
    
    # Count invalid values (NaN, inf, or sentinel void fill values).
    # Use explicit nodata checks rather than an elevation threshold so that
    # valid low-elevation terrain (Dead Sea ~-431 m, Danakil ~-125 m) is not
    # incorrectly flagged.  The -9000 floor covers any value below the deepest
    # point on Earth (~-10,994 m Mariana Trench floor) while safely above the
    # common integer nodata sentinels -32768 and -9999.
    invalid_mask = (
        ~np.isfinite(dem_data) |
        (dem_data == -32768) |   # Standard SRTM/ASTER integer nodata
        (dem_data == -9999) |    # Common floating-point nodata sentinel
        (dem_data < -9000)       # Below any real Earth surface elevation
    )
    
    void_count = np.sum(invalid_mask)
    total_pixels = dem_data.size
    void_percentage = (void_count / total_pixels) * 100.0
    
    if void_percentage > 50.0:
        raise RuntimeError(
            f"DEM has excessive voids during {stage}: "
            f"{void_percentage:.1f}% invalid ({void_count}/{total_pixels} pixels). "
            f"This indicates DEM coverage issues."
        )
    elif void_percentage > DEM_VOID_WARNING_THRESHOLD:
        logger.warning(
            f"⚠️  DEM has significant voids during {stage}: "
            f"{void_percentage:.1f}% invalid ({void_count}/{total_pixels} pixels)"
        )


def _extract_orbit_data(cached_meta: dict, processor) -> tuple[list, list, list, str]:
    """Extract orbit data from cached metadata with validation.
    
    PRIORITY ORDER (tries sources in sequence until sufficient vectors found):
    0. PRECISE ORBIT FILE via sardine.load_orbit_file() (BEST - ~9000+ vectors over 24h)
       - Downloaded .EOF file from Copernicus/AWS
       - Covers full orbital pass, not just scene acquisition
       - Provides best interpolation accuracy for geocoding
    
    1. Processor's cached precise orbit (from apply_precise_orbit step)
       - Stored in processor._cached_orbit_times/positions/velocities
       - Same as Priority 0 but already loaded in memory
       - Typical count: 9000+ vectors
    
    2. Processor's _precise_orbit_records (structured records)  
       - List of dicts with 'time', 'position', 'velocity' keys
       - Alternative storage format for precise orbits
    
    3. cached_meta["orbit_state_vectors"] (SAFE annotation - WARNING: insufficient!)
       - Direct from product annotation XML
       - Typical count: ~17 vectors over ~25s acquisition
       - May not cover scene edges adequately
       - Triggers warning if count < 50
    
    4. Flat keys in cached_meta (legacy fallback)
       - orbit_times, orbit_positions, orbit_velocities as top-level keys
       - Last resort, rarely used
    
    The precise orbit file typically has 9000+ state vectors (10s intervals over 24h),
    while SAFE annotation only has ~17 vectors (~10s intervals around acquisition).
    Using the precise orbit is CRITICAL for accurate geocoding.
    
    Args:
        cached_meta: Cached metadata dictionary from reader
        processor: BackscatterProcessor instance for warnings
        
    Returns:
        Tuple of (orbit_times, orbit_positions, orbit_velocities, status_message)
        
    Raises:
        RuntimeError: If insufficient orbit data is available
    """
    orbit_times = []
    orbit_positions = []
    orbit_velocities = []
    extraction_source = None
    
    if not cached_meta or not isinstance(cached_meta, dict):
        raise RuntimeError(
            "Cannot extract orbit data: cached metadata is empty or invalid"
        )
    
    # PRIORITY 0: Load PRECISE ORBIT FILE directly (BEST SOURCE - ~9000+ vectors)
    # This bypasses all intermediate caches and loads the actual .EOF file
    if processor is not None:
        precise_orbit_path = getattr(processor, "precise_orbit_path", None)
        if precise_orbit_path and Path(precise_orbit_path).exists():
            try:
                import sardine
                if hasattr(sardine, "load_orbit_file"):
                    orbit_data = sardine.load_orbit_file(str(precise_orbit_path))
                    if isinstance(orbit_data, dict) and orbit_data.get("osv_times"):
                        orbit_times = list(orbit_data["osv_times"])
                        orbit_positions = [list(pos) for pos in orbit_data["osv_positions"]]
                        orbit_velocities = [list(vel) for vel in orbit_data["osv_velocities"]]
                        extraction_source = f"precise_orbit_file:{Path(precise_orbit_path).name}"
                        logger.info(
                            f"✅ Loaded {len(orbit_times)} precise orbit vectors from EOF file: "
                            f"{Path(precise_orbit_path).name}"
                        )
            except Exception as e:
                logger.warning(f"Failed to load precise orbit file {precise_orbit_path}: {e}")
                orbit_times = []
                orbit_positions = []
                orbit_velocities = []
    
    # PRIORITY 1: Check processor's cached precise orbit data (from apply_precise_orbit step)
    # This is the second best source - precise orbit file typically has 9000+ state vectors
    if len(orbit_times) < MINIMUM_ORBIT_VECTORS and processor is not None:
        cached_times = getattr(processor, "_cached_orbit_times", None)
        cached_positions = getattr(processor, "_cached_orbit_positions", None)
        cached_velocities = getattr(processor, "_cached_orbit_velocities", None)
        
        if cached_times and cached_positions and len(cached_times) >= MINIMUM_ORBIT_VECTORS:
            try:
                from datetime import datetime
                for i, t in enumerate(cached_times):
                    # Convert to ISO8601 with 'Z' suffix (required by Rust)
                    if isinstance(t, datetime):
                        # Convert datetime to string with 'Z' suffix
                        time_str = t.strftime("%Y-%m-%dT%H:%M:%S.%fZ") if t.tzinfo is None else t.astimezone(__import__('datetime').timezone.utc).strftime("%Y-%m-%dT%H:%M:%S.%fZ")
                        orbit_times.append(time_str)
                    elif isinstance(t, str):
                        # If already string, ensure it ends with 'Z' not '+00:00'
                        time_str = t.replace('+00:00', 'Z') if '+00:00' in t else t
                        orbit_times.append(time_str)
                    else:
                        orbit_times.append(str(t))
                    
                    orbit_positions.append([float(p) for p in cached_positions[i]])
                    if cached_velocities and i < len(cached_velocities):
                        orbit_velocities.append([float(v) for v in cached_velocities[i]])
                    else:
                        orbit_velocities.append([0.0, 0.0, 0.0])
                extraction_source = "processor_cached_precise_orbit"
                logger.info(f"✅ Using {len(orbit_times)} precise orbit vectors from processor cache")
            except (IndexError, TypeError, ValueError) as e:
                logger.warning(f"Failed to extract from processor cached orbit: {e}")
                orbit_times = []
                orbit_positions = []
                orbit_velocities = []
    
    # PRIORITY 2: Check processor's _precise_orbit_records (structured)
    if len(orbit_times) < MINIMUM_ORBIT_VECTORS and processor is not None:
        precise_records = getattr(processor, "_precise_orbit_records", None)
        if precise_records and isinstance(precise_records, list) and len(precise_records) >= MINIMUM_ORBIT_VECTORS:
            orbit_times = []
            orbit_positions = []
            orbit_velocities = []
            try:
                from datetime import datetime
                for rec in precise_records:
                    t = rec.get("time") or rec.get("utc_time")
                    pos = rec.get("position")
                    vel = rec.get("velocity")
                    if t and pos:
                        # Convert to ISO8601 with 'Z' suffix (required by Rust)
                        if isinstance(t, datetime):
                            time_str = t.strftime("%Y-%m-%dT%H:%M:%S.%fZ") if t.tzinfo is None else t.astimezone(__import__('datetime').timezone.utc).strftime("%Y-%m-%dT%H:%M:%S.%fZ")
                            orbit_times.append(time_str)
                        elif isinstance(t, str):
                            time_str = t.replace('+00:00', 'Z') if '+00:00' in t else t
                            orbit_times.append(time_str)
                        else:
                            orbit_times.append(str(t))
                        
                        orbit_positions.append([float(pos[0]), float(pos[1]), float(pos[2])])
                        if vel:
                            orbit_velocities.append([float(vel[0]), float(vel[1]), float(vel[2])])
                        else:
                            orbit_velocities.append([0.0, 0.0, 0.0])
                extraction_source = "processor_precise_orbit_records"
                logger.info(f"✅ Using {len(orbit_times)} precise orbit vectors from processor records")
            except (IndexError, TypeError, ValueError) as e:
                logger.warning(f"Failed to extract from processor precise_orbit_records: {e}")
                orbit_times = []
                orbit_positions = []
                orbit_velocities = []
    
    # PRIORITY 3: Try structured orbit_state_vectors list from cached_meta (SAFE annotation)
    # WARNING: This typically has only ~17 vectors which may not cover the scene adequately.
    # Annotation vectors are sufficient for basic geocoding but degrade accuracy at scene edges.
    # 
    # In strict orbit mode, fail if precise orbit is not available.
    strict_orbit_mode = getattr(processor, "strict_orbit_mode", False)
    if len(orbit_times) < MINIMUM_ORBIT_VECTORS:
        if strict_orbit_mode:
            raise RuntimeError(
                f"Strict orbit mode enabled: Precise orbit file required but not available.\n"
                f"Only {len(orbit_times)} orbit vectors found (minimum {MINIMUM_ORBIT_VECTORS} required).\n"
                f"Annotation vectors (~17) are insufficient for high-precision terrain correction.\n"
                f"Please ensure precise orbit file is available or disable strict_orbit_mode."
            )
        orbit_vectors = cached_meta.get("orbit_state_vectors")
    
        if orbit_vectors and isinstance(orbit_vectors, list):
            if extraction_source is None:
                extraction_source = "orbit_state_vectors"
            orbit_times = []
            orbit_positions = []
            orbit_velocities = []
            for vec in orbit_vectors:
                if isinstance(vec, dict):
                    t = vec.get("time") or vec.get("utc_time")
                    pos = vec.get("position")
                    vel = vec.get("velocity")
                    if t and pos and vel:
                        try:
                            from datetime import datetime
                            # Convert to ISO8601 with 'Z' suffix (required by Rust)
                            if isinstance(t, datetime):
                                time_str = t.strftime("%Y-%m-%dT%H:%M:%S.%fZ") if t.tzinfo is None else t.astimezone(__import__('datetime').timezone.utc).strftime("%Y-%m-%dT%H:%M:%S.%fZ")
                                orbit_times.append(time_str)
                            elif isinstance(t, str):
                                time_str = t.replace('+00:00', 'Z') if '+00:00' in t else t
                                orbit_times.append(time_str)
                            else:
                                orbit_times.append(str(t))
                            
                            orbit_positions.append([float(pos[0]), float(pos[1]), float(pos[2])])
                            orbit_velocities.append([float(vel[0]), float(vel[1]), float(vel[2])])
                        except (IndexError, TypeError, ValueError) as e:
                            logger.warning(f"Skipping malformed orbit vector: {e}")
            
            # CRITICAL WARNING: SAFE annotation orbit may be insufficient for accurate geocoding
            if len(orbit_times) < 50:
                logger.warning(
                    f"⚠️  Using only {len(orbit_times)} orbit vectors from SAFE annotation. "
                    "Precise orbit file has ~9000+ vectors. Geocoding accuracy may be degraded!"
                )

    # Fallback: reconstruct from flat keys (orbit_state_vector_N_x_position, etc.)
    if len(orbit_times) < MINIMUM_ORBIT_VECTORS:
        if extraction_source:
            logger.warning(
                f"Only {len(orbit_times)} vectors from {extraction_source}, "
                f"trying flat key reconstruction"
            )
        extraction_source = "flat_keys"
        
        from datetime import datetime, timezone

        # Find all orbit vector indices
        osv_indices = set()
        for key in cached_meta.keys():
            if key.startswith("orbit_state_vector_") and "_x_position" in key:
                try:
                    idx = int(key.split("_")[3])
                    osv_indices.add(idx)
                except (ValueError, IndexError):
                    pass

        for idx in sorted(osv_indices):
            prefix = f"orbit_state_vector_{idx}_"
            t = cached_meta.get(f"{prefix}time")
            x_pos = cached_meta.get(f"{prefix}x_position")
            y_pos = cached_meta.get(f"{prefix}y_position")
            z_pos = cached_meta.get(f"{prefix}z_position")
            x_vel = cached_meta.get(f"{prefix}x_velocity")
            y_vel = cached_meta.get(f"{prefix}y_velocity")
            z_vel = cached_meta.get(f"{prefix}z_velocity")

            if t and x_pos is not None and y_pos is not None and z_pos is not None:
                # Convert Unix timestamp to ISO format if needed
                t_str = str(t)
                try:
                    ts_val = float(t)
                    # If it looks like a Unix timestamp (> year 2000 in seconds), convert it
                    if ts_val > 946684800:  # Jan 1, 2000 as Unix timestamp
                        dt = datetime.fromtimestamp(ts_val, tz=timezone.utc)
                        t_str = dt.strftime("%Y-%m-%dT%H:%M:%S.%f")[:-3] + "Z"
                except (ValueError, TypeError):
                    pass
                
                try:
                    orbit_times.append(t_str)
                    orbit_positions.append([float(x_pos), float(y_pos), float(z_pos)])
                    # SCIENTIFIC FIX: Require velocity data - zero velocities cause Doppler errors
                    if x_vel is not None and y_vel is not None and z_vel is not None:
                        orbit_velocities.append([float(x_vel), float(y_vel), float(z_vel)])
                    else:
                        raise ValueError(
                            f"Orbit state vector {idx} missing velocity components. "
                            "Velocity data is required for accurate Doppler calculations. "
                            "Check that precise orbit file (.EOF) was applied correctly."
                        )
                except (TypeError, ValueError) as e:
                    logger.warning(f"Skipping malformed flat orbit vector {idx}: {e}")

    # Validate extracted data
    is_valid, vector_count, orbit_msg = validate_orbit_vector_count(
        orbit_times=orbit_times,
        orbit_positions=orbit_positions,
        minimum_vectors=MINIMUM_ORBIT_VECTORS,
    )
    
    if not is_valid:
        raise RuntimeError(
            f"Insufficient orbit data for terrain correction: {orbit_msg}. "
            f"Extracted {vector_count} vectors from {extraction_source}, "
            f"need at least {MINIMUM_ORBIT_VECTORS}."
        )
    
    status_msg = f"Extracted {vector_count} orbit vectors from {extraction_source}"
    return orbit_times, orbit_positions, orbit_velocities, status_msg


def _extract_and_validate_bbox(processor, cached_meta: dict) -> tuple:
    """Extract and validate bounding box from multiple sources.
    
    Args:
        processor: BackscatterProcessor instance
        cached_meta: Cached metadata dictionary
        
    Returns:
        Validated bounding box tuple (min_lon, min_lat, max_lon, max_lat)
        
    Raises:
        RuntimeError: If no valid bounding box can be found
    """
    bbox_sources = []
    validated_bbox = None
    
    # Source 1: processor.product_bbox
    product_bbox = getattr(processor, "product_bbox", None)
    if product_bbox:
        is_valid, error = _validate_bbox(product_bbox, "processor.product_bbox")
        if is_valid:
            validated_bbox = product_bbox
        else:
            bbox_sources.append(("processor.product_bbox", error))
    
    # Source 2: processor.metadata
    if validated_bbox is None:
        meta = getattr(processor, "metadata", None)
        if isinstance(meta, dict):
            try:
                meta_bbox = (
                    float(meta.get("min_longitude")),
                    float(meta.get("min_latitude")),
                    float(meta.get("max_longitude")),
                    float(meta.get("max_latitude")),
                )
                is_valid, error = _validate_bbox(meta_bbox, "processor.metadata")
                if is_valid:
                    validated_bbox = meta_bbox
                else:
                    bbox_sources.append(("processor.metadata", error))
            except (TypeError, ValueError) as e:
                bbox_sources.append(("processor.metadata", f"extraction failed: {e}"))
    
    # Source 3: cached_meta
    if validated_bbox is None and cached_meta and isinstance(cached_meta, dict):
        try:
            cached_bbox = (
                float(cached_meta.get("min_longitude")),
                float(cached_meta.get("min_latitude")),
                float(cached_meta.get("max_longitude")),
                float(cached_meta.get("max_latitude")),
            )
            is_valid, error = _validate_bbox(cached_bbox, "cached_metadata")
            if is_valid:
                validated_bbox = cached_bbox
            else:
                bbox_sources.append(("cached_metadata", error))
        except (TypeError, ValueError) as e:
            bbox_sources.append(("cached_metadata", f"extraction failed: {e}"))
    
    # All sources failed
    if validated_bbox is None:
        error_details = "; ".join([f"{src}: {err}" for src, err in bbox_sources])
        raise RuntimeError(
            f"No valid bounding box available for terrain correction. "
            f"Tried: {error_details}"
        )
    
    # FIXED: Apply bbox refinement if enabled and we have image dimensions
    # This prevents metadata footprint from being much larger than actual coverage
    refine_enabled = getattr(processor, 'options', {}).get('refine_geocoding_bbox', True)
    working_data = getattr(processor, '_working_data', None)
    
    if refine_enabled and working_data is not None and isinstance(working_data, np.ndarray):
        rows, cols = working_data.shape[:2]
        
        # Get pixel spacing from metadata
        range_spacing = None
        azimuth_spacing = None
        
        meta = getattr(processor, "metadata", {}) or {}
        range_spacing = meta.get("range_pixel_spacing") or meta.get("native_range_pixel_spacing")
        azimuth_spacing = meta.get("azimuth_pixel_spacing") or meta.get("native_azimuth_pixel_spacing")
        
        # CRITICAL FIX: Account for multilook factors when computing expected extent
        # The spacing from metadata is typically NATIVE spacing, but rows/cols are MULTILOOKED.
        # We must multiply spacing by multilook factor to get the actual ground spacing.
        range_ml = getattr(processor, "_range_multilook_factor", None) or 1.0
        azimuth_ml = getattr(processor, "_azimuth_multilook_factor", None) or 1.0
        
        if range_spacing and azimuth_spacing:
            try:
                # Use multilooked spacing = native spacing * multilook factor
                effective_range_spacing = float(range_spacing) * float(range_ml)
                effective_azimuth_spacing = float(azimuth_spacing) * float(azimuth_ml)
                
                refined = refine_geocoding_bbox(
                    validated_bbox,
                    rows=rows,
                    cols=cols,
                    range_spacing_m=effective_range_spacing,
                    azimuth_spacing_m=effective_azimuth_spacing,
                )
                if refined is not None:
                    refined_bbox, metrics = refined
                    shrink_lat = metrics.get("shrink_lat_pct", 0) * 100
                    shrink_lon = metrics.get("shrink_lon_pct", 0) * 100
                    if shrink_lat > 0 or shrink_lon > 0:
                        logger.info(
                            f"Refined geocoding bbox: lat shrink {shrink_lat:.1f}%, "
                            f"lon shrink {shrink_lon:.1f}%"
                        )
                        validated_bbox = tuple(refined_bbox)
            except Exception as e:
                logger.debug(f"Bbox refinement skipped: {e}")
    
    return validated_bbox


def _compute_sar_footprint_bbox(
    processor,
    working_data: np.ndarray,
    cached_meta: dict,
    reader,
    orbit_times: list,
    orbit_positions: list,
    orbit_velocities: list,
    real_metadata: dict,
    burst_timing_json: Optional[str],
    margin_percent: float = 2.0,
) -> tuple:
    """
    Compute bounding box from SAR image footprint geometry.
    
    This is the CORRECT approach: derive bbox from actual SAR image corners
    using Range-Doppler geocoding, NOT from annotation/metadata bbox.
    
    Args:
        processor: BackscatterProcessor instance
        working_data: SAR image array (post-merge, post-multilook)
        cached_meta: Cached metadata dictionary
        reader: SLC reader instance
        orbit_times: Orbit state vector times
        orbit_positions: Orbit positions (ECEF)
        orbit_velocities: Orbit velocities (ECEF)
        real_metadata: Metadata dict for Range-Doppler params
        burst_timing_json: Burst timing data (JSON string)
        margin_percent: Margin to add around footprint (default: 2%)
        
    Returns:
        Tuple of (min_lon, min_lat, max_lon, max_lat)
        
    Raises:
        RuntimeError: If bbox computation fails
    """
    try:
        logger.info("🎯 Computing SAR footprint bbox from image corners (SCIENTIFIC MODE)")
        logger.info(f"   Image shape: {working_data.shape}")
        logger.info(f"   Margin: {margin_percent:.1f}%")
        
        # CRITICAL FIX: Extract only numeric metadata for Rust function
        # The Rust function signature is: real_metadata: HashMap<String, f64>
        # It cannot accept strings, dicts, or lists - only float values
        numeric_metadata = {}
        for key, value in real_metadata.items():
            # Skip non-numeric values (subswaths, burst_timings, etc.)
            if isinstance(value, (int, float)):
                numeric_metadata[key] = float(value)
            elif isinstance(value, str):
                # Try to parse strings that might be numbers
                try:
                    numeric_metadata[key] = float(value)
                except (ValueError, TypeError):
                    # Skip non-numeric strings
                    pass
        
        logger.debug(f"   Numeric metadata keys passed to Rust: {sorted(numeric_metadata.keys())}")
        
        # Check if strict science mode requires accurate bbox
        strict_science = bool(getattr(processor, 'strict_science', False))
        
        # TEMPORARY LIMITATION: Rust bbox computation lacks subswath geometry (hardcoded empty HashMap)
        # This causes incorrect corner geocoding for merged IW data.
        # In strict science mode, we fail hard rather than silently using oversized bbox.
        if strict_science:
            raise RuntimeError(
                "STRICT SCIENCE MODE FAILURE: SAR footprint bbox computation not yet supported.\n"
                "The Rust implementation currently lacks per-subswath geometry needed for accurate\n"
                "corner geocoding of merged IW TOPSAR data. Using the oversized metadata bbox\n"
                "would introduce systematic geocoding errors of 100s of meters at swath edges.\n\n"
                "To proceed:\n"
                "  1. Set processor.strict_science = False to allow metadata bbox fallback, OR\n"
                "  2. Implement proper subswath geometry passing to Rust bbox function\n\n"
                "This error prevents scientifically invalid output in strict mode."
            )
        
        # In non-strict mode, fall back to metadata bbox with warning
        logger.warning(
            "⚠️  SAR footprint bbox computation skipped: Rust function lacks subswath geometry.\n"
            "   This causes incorrect corner geocoding for merged IW data.\n"
            "   Falling back to metadata bbox (oversized but safer for now).\n"
            "   Enable strict_science mode to prevent this fallback."
        )
        return None
        
        # DISABLED until Rust subswath support is added:
        # Call Rust function to compute bbox from SAR image corners
        result = sardine.compute_sar_footprint_bbox(
            working_data,
            orbit_times,
            orbit_positions,
            orbit_velocities,
            numeric_metadata,  # ✅ FIXED: Only numeric values
            reader,
            burst_timing_json,
            margin_percent,
        )
        
        # Extract bbox from result dict
        bbox = tuple(result["bbox"])  # [min_lon, min_lat, max_lon, max_lat]
        bbox_source = result["bbox_source"]
        footprint_area_km2 = result["footprint_area_km2"]
        
        logger.info(f"✅ SAR footprint bbox computed:")
        logger.info(f"   Lat: [{bbox[1]:.6f}, {bbox[3]:.6f}]")
        logger.info(f"   Lon: [{bbox[0]:.6f}, {bbox[2]:.6f}]")
        logger.info(f"   Area: {footprint_area_km2:.1f} km²")
        logger.info(f"   Source: {bbox_source}")
        
        # Store provenance info in processor for metadata export
        if hasattr(processor, '_bbox_provenance'):
            processor._bbox_provenance.update({
                'bbox_source': bbox_source,
                'bbox_footprint': bbox,
                'footprint_area_km2': footprint_area_km2,
                'margin_percent': margin_percent,
                'bbox_without_margin': result.get('bbox_without_margin'),
            })
        else:
            processor._bbox_provenance = {
                'bbox_source': bbox_source,
                'bbox_footprint': bbox,
                'footprint_area_km2': footprint_area_km2,
                'margin_percent': margin_percent,
                'bbox_without_margin': result.get('bbox_without_margin'),
            }
        
        return bbox
        
    except Exception as e:
        logger.error(f"❌ Failed to compute SAR footprint bbox: {e}")
        raise RuntimeError(
            f"Failed to compute SAR footprint bbox from image corners: {e}\n"
            f"This is required for correct terrain correction grid sizing.\n"
            f"Check orbit data coverage and metadata completeness."
        ) from e


def _validate_multilook_factors(
    processor, 
    working_data: np.ndarray
) -> tuple[float, float]:
    """Validate multilook factors match actual image dimensions.
    
    Args:
        processor: BackscatterProcessor instance
        working_data: Current working data array
        
    Returns:
        Tuple of (range_multilook_factor, azimuth_multilook_factor)
    """
    range_ml = getattr(processor, "_range_multilook_factor", None)
    azimuth_ml = getattr(processor, "_azimuth_multilook_factor", None)
    
    # Get native dimensions if available
    native_range = getattr(processor, "_native_range_samples", None)
    native_azimuth = getattr(processor, "_native_azimuth_lines", None)
    
    actual_azimuth, actual_range = working_data.shape
    
    # Validate range multilook factor
    if range_ml is None or range_ml <= 0:
        if native_range and native_range > 0 and actual_range > 0:
            inferred_range_ml = native_range / actual_range
            logger.warning(
                f"Missing range_multilook_factor; inferred {inferred_range_ml:.2f} "
                f"from dimensions ({native_range} / {actual_range})"
            )
            range_ml = inferred_range_ml
        else:
            logger.warning(
                f"Missing range_multilook_factor and cannot infer from dimensions. "
                f"Using default {DEFAULT_MULTILOOK_FACTOR}. "
                f"This may cause coordinate scaling errors."
            )
            range_ml = DEFAULT_MULTILOOK_FACTOR
    else:
        # Validate factor matches dimensions
        if native_range and native_range > 0:
            expected_range = native_range / range_ml
            diff = abs(expected_range - actual_range)
            if diff > 2:  # Allow 2 pixel tolerance
                logger.warning(
                    f"⚠️  Range multilook factor mismatch: factor={range_ml}, "
                    f"native={native_range}, expected_multilooked={expected_range:.0f}, "
                    f"actual={actual_range} (diff={diff:.0f})"
                )
    
    # Validate azimuth multilook factor
    if azimuth_ml is None or azimuth_ml <= 0:
        if native_azimuth and native_azimuth > 0 and actual_azimuth > 0:
            inferred_azimuth_ml = native_azimuth / actual_azimuth
            logger.warning(
                f"Missing azimuth_multilook_factor; inferred {inferred_azimuth_ml:.2f} "
                f"from dimensions ({native_azimuth} / {actual_azimuth})"
            )
            azimuth_ml = inferred_azimuth_ml
        else:
            logger.warning(
                f"Missing azimuth_multilook_factor and cannot infer from dimensions. "
                f"Using default {DEFAULT_MULTILOOK_FACTOR}. "
                f"This may cause coordinate scaling errors."
            )
            azimuth_ml = DEFAULT_MULTILOOK_FACTOR
    else:
        # Validate factor matches dimensions
        if native_azimuth and native_azimuth > 0:
            expected_azimuth = native_azimuth / azimuth_ml
            diff = abs(expected_azimuth - actual_azimuth)
            if diff > 2:  # Allow 2 pixel tolerance
                logger.warning(
                    f"⚠️  Azimuth multilook factor mismatch: factor={azimuth_ml}, "
                    f"native={native_azimuth}, expected_multilooked={expected_azimuth:.0f}, "
                    f"actual={actual_azimuth} (diff={diff:.0f})"
                )
    
    return float(range_ml), float(azimuth_ml)


def _resample_dem_to_shape(dem: np.ndarray, target_shape: tuple[int, int]) -> np.ndarray:
    """Resample DEM to the target image shape.

    Uses scipy.ndimage.zoom when available for speed; falls back to a
    separable np.interp implementation to avoid extra dependencies.
    
    Args:
        dem: Input DEM array
        target_shape: Target (rows, cols) shape
        
    Returns:
        Resampled DEM array matching target_shape
        
    Raises:
        ValueError: If input DEM is invalid or target shape is invalid
    """
    if dem is None or dem.size == 0:
        raise ValueError("Cannot resample empty or None DEM")
    
    if len(target_shape) != 2 or target_shape[0] <= 0 or target_shape[1] <= 0:
        raise ValueError(f"Invalid target shape: {target_shape}")

    if dem.shape == target_shape:
        return dem
    
    # Check for upsampling case (DEM lower resolution than SAR) - this is problematic
    src_rows, src_cols = dem.shape
    tgt_rows, tgt_cols = target_shape
    is_upsampling = (tgt_rows > src_rows) or (tgt_cols > src_cols)
    
    if is_upsampling:
        upsample_ratio = max(tgt_rows / src_rows, tgt_cols / src_cols)
        if upsample_ratio > 2.0:
            logger.warning(
                f"⚠️  SCIENTIFIC WARNING: DEM being upsampled by {upsample_ratio:.1f}× for terrain flattening! "
                f"DEM resolution ({src_rows}×{src_cols}) is coarser than SAR image ({tgt_rows}×{tgt_cols}). "
                f"This smooths DEM gradients and can cause inaccurate cos(θ_local) correction. "
                f"Consider using a higher resolution DEM (e.g., Copernicus DEM 10m) or reducing multilook factors."
            )

    try:  # Prefer fast C-backed zoom if SciPy is present
        import scipy.ndimage  # type: ignore

        zoom_factors = (target_shape[0] / dem.shape[0], target_shape[1] / dem.shape[1])
        
        # Use cubic interpolation (order=3) for gradient-preserving resampling
        # Bilinear (order=1) smooths DEM gradients which degrades terrain flattening accuracy
        # Cubic better preserves local slope information critical for cos(θ_local)
        interpolation_order = 3 if is_upsampling else 1  # Only use cubic when needed
        
        result = scipy.ndimage.zoom(dem, zoom_factors, order=interpolation_order)
        
        # Validate output dimensions
        if result.shape != target_shape:
            logger.warning(
                f"DEM resample shape mismatch: expected {target_shape}, got {result.shape}. "
                f"Trimming/padding to match."
            )
            # Trim or pad to exact target shape
            final = np.zeros(target_shape, dtype=result.dtype)
            copy_rows = min(result.shape[0], target_shape[0])
            copy_cols = min(result.shape[1], target_shape[1])
            final[:copy_rows, :copy_cols] = result[:copy_rows, :copy_cols]
            return final
            
        return result
    except ImportError:
        logger.warning(
            "SciPy not available for DEM resampling; using slower numpy fallback. "
            "Install scipy for faster processing: pip install scipy"
        )
    except Exception as e:
        logger.warning(
            f"SciPy DEM resampling failed: {e}. Falling back to numpy interpolation."
        )
    
    # Fallback: separable bilinear interpolation using numpy only
    src_rows, src_cols = dem.shape
    tgt_rows, tgt_cols = target_shape

    col_positions = np.linspace(0, src_cols - 1, tgt_cols, dtype=np.float32)
    row_positions = np.linspace(0, src_rows - 1, tgt_rows, dtype=np.float32)

    # Resample columns for each source row
    resampled_cols = np.empty((src_rows, tgt_cols), dtype=np.float32)
    src_col_indices = np.arange(src_cols, dtype=np.float32)
    for i in range(src_rows):
        resampled_cols[i, :] = np.interp(col_positions, src_col_indices, dem[i, :])

    # Resample rows for each target column
    resampled = np.empty((tgt_rows, tgt_cols), dtype=np.float32)
    src_row_indices = np.arange(src_rows, dtype=np.float32)
    for j in range(tgt_cols):
        resampled[:, j] = np.interp(row_positions, src_row_indices, resampled_cols[:, j])

    return resampled


def _calculate_effective_spacing(
    dem_bbox: tuple, rows: int, cols: int, default_resolution: float
) -> float:
    """Calculate effective DEM pixel spacing in meters from bbox and dimensions.
    
    Args:
        dem_bbox: Bounding box (min_lon, min_lat, max_lon, max_lat)
        rows: Number of rows in the image
        cols: Number of columns in the image
        default_resolution: Fallback resolution if calculation fails
        
    Returns:
        Effective spacing in meters
    """
    if not dem_bbox or len(dem_bbox) != 4:
        logger.debug(
            f"Invalid or missing DEM bbox; using default spacing {default_resolution}m"
        )
        return float(default_resolution)
    
    try:
        min_lon, min_lat, max_lon, max_lat = dem_bbox
        
        # Validate bbox values
        if not all(np.isfinite([min_lon, min_lat, max_lon, max_lat])):
            raise ValueError("Non-finite values in bbox")
        
        if min_lat >= max_lat or min_lon >= max_lon:
            raise ValueError(f"Invalid bbox extent: {dem_bbox}")
        
        mean_lat = 0.5 * (min_lat + max_lat)
        lat_rad = np.deg2rad(mean_lat)
        
        # Calculate meters per degree at this latitude
        m_per_deg_lat = (
            METERS_PER_DEG_LAT_BASE 
            - 559.822 * np.cos(2.0 * lat_rad) 
            + 1.175 * np.cos(4.0 * lat_rad)
        )
        m_per_deg_lon = METERS_PER_DEG_LON_BASE * np.cos(lat_rad)
        
        # Calculate pixel spacing
        dy_m = abs(max_lat - min_lat) * m_per_deg_lat / max(rows, 1)
        dx_m = abs(max_lon - min_lon) * m_per_deg_lon / max(cols, 1)
        
        effective_spacing = float(0.5 * (dx_m + dy_m))
        
        # Validate result
        if not np.isfinite(effective_spacing) or effective_spacing <= 0:
            raise ValueError(f"Invalid computed spacing: {effective_spacing}")
        
        # Sanity check: spacing should be reasonable (1m to 1000m)
        if effective_spacing < 1.0 or effective_spacing > 1000.0:
            logger.warning(
                f"Computed DEM spacing {effective_spacing:.1f}m outside expected range [1, 1000]m. "
                f"Using default {default_resolution}m instead."
            )
            return float(default_resolution)
        
        return effective_spacing
        
    except Exception as e:
        logger.warning(
            f"Failed to calculate effective DEM spacing: {e}. "
            f"Using default {default_resolution}m."
        )
        return float(default_resolution)


def _get_primary_subswath(processor) -> str:
    """Extract primary subswath ID with validation.
    
    Args:
        processor: BackscatterProcessor instance
        
    Returns:
        Primary subswath ID (e.g., "IW1")
        
    Raises:
        RuntimeError: If no valid subswath can be determined
    """
    # First try explicit primary subswath
    primary_subswath = getattr(processor, "_primary_subswath", None)
    if primary_subswath and isinstance(primary_subswath, str):
        return primary_subswath
    
    # Try calibrated subswaths
    calibrated_subswaths = getattr(processor, "_calibrated_subswaths", None)
    if calibrated_subswaths and isinstance(calibrated_subswaths, dict) and calibrated_subswaths:
        primary_subswath = next(iter(calibrated_subswaths.keys()), None)
        if primary_subswath:
            logger.debug(f"Using first calibrated subswath as primary: {primary_subswath}")
            return primary_subswath
    
    # Try used subswaths for merge
    used_subswaths = getattr(processor, "_used_subswaths_for_merge", None)
    if used_subswaths and len(used_subswaths) > 0:
        primary_subswath = used_subswaths[0]
        if primary_subswath:
            logger.debug(f"Using first merge subswath as primary: {primary_subswath}")
            return primary_subswath
    
    raise RuntimeError(
        "Cannot determine primary subswath for terrain flattening. "
        "No calibrated subswaths or merge subswaths available. "
        "This indicates earlier pipeline stages failed silently."
    )


def _get_terrain_flatten_function():
    """Get the terrain flattening function from sardine bindings.
    
    Returns:
        Callable terrain flattening function
        
    Raises:
        RuntimeError: If no terrain flattening function is available
    """
    # Try primary function name first
    terrain_flatten_fn = getattr(sardine, "terrain_flattening", None)
    if terrain_flatten_fn is not None and callable(terrain_flatten_fn):
        return terrain_flatten_fn
    
    # Try alternative function name
    terrain_flatten_fn = getattr(sardine, "apply_scientific_terrain_flattening", None)
    if terrain_flatten_fn is not None and callable(terrain_flatten_fn):
        return terrain_flatten_fn
    
    # List available functions for debugging
    available_funcs = [
        name for name in dir(sardine) 
        if "terrain" in name.lower() and callable(getattr(sardine, name, None))
    ]
    
    raise RuntimeError(
        f"Terrain flattening function not available in sardine bindings. "
        f"Expected 'terrain_flattening' or 'apply_scientific_terrain_flattening'. "
        f"Available terrain-related functions: {available_funcs}"
    )


def run_terrain_flattening(processor) -> None:
    """Apply DEM-based terrain flattening before speckle filtering.
    
        DEPRECATED: RTC is now integrated directly into the geocoding step
        (terrain_correction_with_rtc in the Rust core). This legacy, pre-
        geocoding terrain flattening should not be combined with RTC.

        Behaviour (Jan 2026):
        - If processor.rtc_mode is not "none": this step is SKIPPED to avoid
            double terrain correction. The correct path is calibrated σ⁰ →
            terrain_correction_with_rtc → γ⁰_tc.
        - If processor.rtc_mode is "none" but terrain_flatten=True: we allow
            this legacy path for backwards compatibility, producing γ⁰_tc on the
            SLC grid prior to geocoding. New workflows should prefer RTC.
    """
    processor.announce_step(8, "Terrain Flattening", "Applying DEM-based terrain flattening")
    step_start = time.time()

    working_data = getattr(processor, "_working_data", None)

    # DEPRECATED: Skip if RTC is enabled in terrain correction (new default behavior)
    # RTC is now integrated into geocoding for scientific correctness
    rtc_mode = getattr(processor, "rtc_mode", "area")
    if rtc_mode and rtc_mode.lower() != "none":
        step_duration = time.time() - step_start
        processor.log_step(
            8,
            "Terrain Flattening",
            "skipped",
            f"RTC is integrated into geocoding (rtc_mode={rtc_mode})",
            step_duration,
        )
        print(f"   ℹ️  Legacy terrain flattening skipped - RTC now integrated into geocoding")
        print(f"   ℹ️  Using rtc_mode='{rtc_mode}' during terrain correction (Small 2011 method)")
        return

    if not processor.terrain_flatten:
        step_duration = time.time() - step_start
        processor.log_step(8, "Terrain Flattening", "skipped", "Disabled by user", step_duration)
        return

    try:
        if not isinstance(working_data, np.ndarray):
            raise ValueError("Working data is not a valid numpy array")
        
        if working_data.ndim != 2:
            raise ValueError(
                f"Working data has invalid dimensions: {working_data.ndim}D (expected 2D)"
            )
        
        if working_data.size == 0:
            raise ValueError("Working data is empty")

        rows, cols = working_data.shape

        if not getattr(processor, "metadata", None) or not isinstance(processor.metadata, dict):
            raise ValueError("No valid metadata available for terrain flattening")

        subswath_specific_bbox = None
        requested_subswaths = processor.options.get("subswaths") if hasattr(processor, "options") else None
        if requested_subswaths and hasattr(processor, "reader") and processor.reader:
            try:
                # CRITICAL FIX: Use processor.polarization (the actual processing polarization)
                # instead of options.get("polarizations") which may default to ["VV"] incorrectly
                pol = getattr(processor, "polarization", None) or "VV"
                bbox_dict = processor.reader.get_subswath_bounding_box(list(requested_subswaths), pol)
                if bbox_dict and isinstance(bbox_dict, dict) and len(bbox_dict) == 4:
                    subswath_specific_bbox = (
                        float(bbox_dict.get("min_lon")),
                        float(bbox_dict.get("min_lat")),
                        float(bbox_dict.get("max_lon")),
                        float(bbox_dict.get("max_lat")),
                    )
                    print(f"   ✅ Using subswath-specific bbox for DEM: {subswath_specific_bbox}")
            except Exception as bbox_exc:
                logger.warning(f"Could not retrieve subswath bbox for DEM cropping: {bbox_exc}")

        dem_bbox = subswath_specific_bbox or getattr(processor, "product_bbox", None)
        dem_resolution = (
            processor.dem_resolution 
            if processor.dem_resolution is not None 
            else DEFAULT_DEM_RESOLUTION_M
        )
        cosine_clip = processor.cosine_clip_threshold

        dem_data, dem_transform, dem_crs = processor._load_dem_for_scene(dem_bbox, dem_resolution)
        
        if dem_data is None or dem_data.size == 0:
            raise RuntimeError(
                "DEM loading failed - no elevation data available for terrain flattening"
            )

        # Validate DEM void percentage
        _validate_dem_voids(dem_data, stage="terrain flattening")

        # Resample DEM to match the working image shape expected by the terrain flattener
        dem_data_resampled = _resample_dem_to_shape(
            np.asarray(dem_data, dtype=np.float32), (rows, cols)
        )
        
        # Validate resampled DEM
        if dem_data_resampled.shape != (rows, cols):
            raise RuntimeError(
                f"DEM resampling failed: expected shape ({rows}, {cols}), "
                f"got {dem_data_resampled.shape}"
            )

        # Calculate effective DEM pixel spacing in meters
        effective_spacing_m = _calculate_effective_spacing(
            dem_bbox, rows, cols, dem_resolution
        )
        logger.debug(f"Effective DEM spacing: {effective_spacing_m:.2f}m")

        # Get primary subswath with validation
        primary_subswath = _get_primary_subswath(processor)
        logger.debug(f"Using primary subswath: {primary_subswath}")

        # Get terrain flattening function with clear error
        terrain_flatten_fn = _get_terrain_flatten_function()

        tf_result = terrain_flatten_fn(
            working_data,
            dem_data_resampled,
            str(processor.input_path),
            primary_subswath,
            None,
            None,
            float(effective_spacing_m),
        )

        # Validate terrain flattening result with clear error messages
        if tf_result is None:
            raise RuntimeError(
                "Terrain flattening returned None - this indicates a failure in the Rust bindings"
            )
        
        flattened_data = None
        if isinstance(tf_result, dict):
            # Try primary key first
            flattened_data = tf_result.get("data")
            if flattened_data is None:
                # Try alternative key
                flattened_data = tf_result.get("flattened_data")
            if flattened_data is None:
                available_keys = list(tf_result.keys())
                raise RuntimeError(
                    f"Terrain flattening result dict missing data. "
                    f"Expected 'data' or 'flattened_data' key. "
                    f"Available keys: {available_keys}"
                )
        elif isinstance(tf_result, np.ndarray):
            flattened_data = tf_result
        else:
            raise RuntimeError(
                f"Unexpected terrain flattening result type: {type(tf_result)}. "
                f"Expected dict or numpy.ndarray."
            )

        if not isinstance(flattened_data, np.ndarray):
            raise RuntimeError(
                f"Terrain flattening data is not a numpy array: {type(flattened_data)}"
            )
        
        if flattened_data.ndim != 2:
            raise RuntimeError(
                f"Terrain flattening produced {flattened_data.ndim}D array, expected 2D"
            )
        
        if flattened_data.size == 0:
            raise RuntimeError("Terrain flattening produced empty array")
        
        # Validate output shape matches input
        if flattened_data.shape != working_data.shape:
            logger.warning(
                f"Terrain flattening changed data shape: {working_data.shape} -> {flattened_data.shape}"
            )

        flattened_array = np.asarray(flattened_data, dtype=np.float32)
        
        # Diagnostic logging: After terrain flattening
        processor._log_diagnostic_statistics(
            "After Terrain Flattening",
            flattened_array,
            context=f"calibration_type={processor.calibration_type}"
        )
        
        processor._working_data = flattened_array

        step_duration = time.time() - step_start
        processor.log_step(
            8,
            "Terrain Flattening",
            "success",
            f"Applied terrain flattening ({rows}x{cols} -> {flattened_data.shape})",
            step_duration,
        )
        processor._record_stage_timing("Terrain Flattening", step_duration)

    except Exception as exc:
        step_duration = time.time() - step_start
        processor.log_step(8, "Terrain Flattening", "error", f"Terrain flattening failed: {exc}", step_duration)
        raise


def run_terrain_correction(processor) -> None:
    """Execute scientific Range-Doppler terrain correction."""
    from sardine.processors.backscatter.processor import STEP_NUMBERS
    step_number = STEP_NUMBERS["Terrain Correction"]
    processor.announce_step(step_number, "Terrain Correction", "Range-Doppler geocoding to analysis-ready grid")
    step_start = time.time()

    working_data = getattr(processor, "_working_data", None)

    try:
        if not processor.geocode:
            step_duration = time.time() - step_start
            processor.log_step(step_number, "Terrain Correction", "skipped", "Geocoding disabled by user", step_duration)
            return

        if not isinstance(working_data, np.ndarray):
            raise RuntimeError("SCIENTIFIC MODE FAILURE: Working data is not a valid numpy array")

        if not isinstance(processor.metadata, dict):
            raise RuntimeError("SCIENTIFIC MODE FAILURE: Missing metadata required for terrain correction")

        # Try to get orbit data from reader's cached metadata
        reader = getattr(processor, "reader", None)
        if reader is None:
            raise RuntimeError("SCIENTIFIC MODE FAILURE: No SLC reader available for terrain correction")

        cached_meta = None
        try:
            cached_meta = reader.get_cached_metadata()
        except Exception as e:
            raise RuntimeError(
                f"Failed to read cached metadata from SLC reader: {e}. "
                "This is required for terrain correction. Check that orbit files were applied."
            )

        if not cached_meta or not isinstance(cached_meta, dict):
            raise RuntimeError(
                "Cached metadata is empty or invalid. "
                "This indicates orbit file application failed or metadata was not cached."
            )

        # Extract orbit data with validation using helper function
        orbit_times, orbit_positions, orbit_velocities, orbit_msg = _extract_orbit_data(
            cached_meta, processor
        )
        print(f"   ✅ {orbit_msg}")

        # TC-1 FIX: Extract metadata bbox for fallback/comparison
        # Extract metadata bbox first - it's used as fallback if SAR footprint computation fails
        try:
            metadata_bbox = _extract_and_validate_bbox(processor, cached_meta)
            print(f"   📋 Metadata bbox extracted: {metadata_bbox}")
        except RuntimeError as metadata_error:
            step_duration = time.time() - step_start
            processor.log_step(
                step_number,
                "Terrain Correction",
                "error",
                f"Could not extract metadata bbox: {metadata_error}",
                step_duration,
            )
            raise RuntimeError(
                f"Terrain correction cannot proceed without valid bounding box: {metadata_error}"
            )

        # Validate and get multilook factors using helper function
        range_ml, azimuth_ml = _validate_multilook_factors(processor, working_data)
        
        # Get DEM cache directory
        dem_cache_dir = getattr(processor, "dem_cache_dir", None) or "/tmp/sardine_dem_cache"

        # Build real_metadata dict with required SLC parameters
        real_metadata = {}
        meta = processor.metadata
        for key in [
            "azimuth_time_interval",
            "range_pixel_spacing",
            "native_range_pixel_spacing",
            "native_azimuth_pixel_spacing",  # Required for terrain correction
            "azimuth_pixel_spacing",         # Fallback name
            "slant_range_time",
            "radar_frequency",
            "wavelength",                    # Required for Range-Doppler
            "first_line_utc_seconds",
            "number_of_lines",
            "number_of_samples",
            "total_azimuth_lines",           # Used for azimuth bounds
            "near_range",
            "far_range",
            "prf",
            "product_start_time_abs",        # Timing parameters
            "product_stop_time_abs",
            "product_duration",
            "orbit_ref_epoch_utc",
            "product_start_rel_s",
            "incidence_angle_mid_swath",     # Reference incidence angle for RTC (Small 2011)
        ]:
            if key in meta:
                try:
                    real_metadata[key] = float(meta[key])
                except (TypeError, ValueError):
                    pass
            elif cached_meta and key in cached_meta:
                try:
                    real_metadata[key] = float(cached_meta[key])
                except (TypeError, ValueError):
                    pass

        # Add validated multilook factors to metadata for proper coordinate scaling
        # The Rust terrain correction needs these to scale SAR coordinates from native to multilooked dimensions
        real_metadata["range_multilook_factor"] = range_ml
        real_metadata["azimuth_multilook_factor"] = azimuth_ml
        
        print(
            f"   📐 Multilook factors for terrain correction: range={range_ml}, "
            f"azimuth={azimuth_ml}"
        )

        # *** CRITICAL FIX: Pass target output CRS to Rust terrain correction ***
        # Without this, Rust defaults to EPSG:4326 (geographic coordinates) instead of
        # the user-requested projection (e.g., EPSG:32632 for UTM zone 32N).
        # This caused the bug where GeoTIFF had geographic coordinates but was labeled as UTM.
        if processor.output_epsg is not None:
            real_metadata["target_output_epsg"] = float(processor.output_epsg)
            print(f"   🗺️  Target output CRS: EPSG:{processor.output_epsg}")
        else:
            # Fallback: resolve output EPSG from scene metadata if not explicitly set
            if hasattr(processor, "_resolve_output_epsg"):
                resolved_epsg = processor._resolve_output_epsg()
                real_metadata["target_output_epsg"] = float(resolved_epsg)
                print(f"   🗺️  Resolved output CRS: EPSG:{resolved_epsg}")
            else:
                print("   ⚠️  No output EPSG specified, Rust will default to EPSG:4326 (geographic)")

        # *** CRITICAL FIX: Override total_azimuth_lines with ACTUAL image dimensions ***
        # The annotation metadata has pre-deburst/pre-clipped values (e.g., 13672), but the
        # actual merged/debursted/multilooked image may have fewer lines (e.g., 12471).
        # Range-Doppler needs the ACTUAL multilooked image dimensions, not the annotation values.
        actual_azimuth_lines, actual_range_samples = working_data.shape
        real_metadata["total_azimuth_lines"] = float(actual_azimuth_lines)
        real_metadata["number_of_lines"] = float(actual_azimuth_lines)
        real_metadata["number_of_samples"] = float(actual_range_samples)
        
        # Also compute the native (pre-multilook) dimensions for bounds checking
        native_azimuth_lines = actual_azimuth_lines * azimuth_ml
        native_range_samples = actual_range_samples * range_ml
        
        print(
            f"   📊 Image dimensions for terrain correction: multilooked={actual_azimuth_lines}x{actual_range_samples}, "
            f"native={int(native_azimuth_lines)}x{int(native_range_samples)}"
        )

        # Add subswath geometry from geocoding cache if available
        geocoding_cache = getattr(processor, "_geocoding_metadata_cache", None)
        if geocoding_cache:
            # The cache may have these under different keys
            subswath_data_all = (
                geocoding_cache.get("subswaths") or
                geocoding_cache.get("subswath_metadata") or
                geocoding_cache.get("subswath_geometry")
            )
            if subswath_data_all:
                # CRITICAL FIX: subswath_data_all has structure {pol: {swath_id: {...}}}
                # We need to extract the subswaths for the current polarization
                # The structure from get_all_iw_subswaths() is: {polarization: {swath_id: subswath_dict}}
                if isinstance(subswath_data_all, dict):
                    # Find the entry for the current polarization
                    # The processor has the current polarization being processed
                    current_pol = getattr(processor, 'polarization', None) or metadata.get('polarization', 'VV')
                    
                    # Try exact match first
                    subswath_data = subswath_data_all.get(current_pol)
                    if not subswath_data:
                        # Try with "Polarization::" prefix (Rust Debug format)
                        pol_key = f"Polarization::{current_pol}"
                        subswath_data = subswath_data_all.get(pol_key)
                    if not subswath_data:
                        # Try all keys that contain the polarization name
                        for key in subswath_data_all.keys():
                            if current_pol in str(key):
                                subswath_data = subswath_data_all[key]
                                break
                    
                    if subswath_data and isinstance(subswath_data, dict):
                        print(f"   📊 Using subswath geometry for {current_pol}: {len(subswath_data)} subswaths")
                        real_metadata["subswaths"] = subswath_data
                    else:
                        print(f"   ⚠️  Could not find subswath data for polarization {current_pol}")
                        print(f"       Available keys: {list(subswath_data_all.keys())}")
                else:
                    # Fallback: if it's not nested by polarization, use as-is
                    real_metadata["subswaths"] = subswath_data_all
            
            # Also merge any cached metadata
            cached_from_cache = geocoding_cache.get("cached_metadata")
            if cached_from_cache and isinstance(cached_from_cache, dict):
                for k, v in cached_from_cache.items():
                    if k not in real_metadata:
                        real_metadata[k] = v

        # *** NORMALIZE METADATA FOR RUST BRIDGE ***
        # This ensures all key names match what the Rust code expects
        is_multilooked = range_ml > 1.0 or azimuth_ml > 1.0
        range_looks = int(range_ml)
        azimuth_looks = int(azimuth_ml)
        
        real_metadata = normalize_metadata_for_rust(
            real_metadata,
            is_multilooked=is_multilooked,
            range_looks=range_looks,
            azimuth_looks=azimuth_looks,
            log_changes=True,
        )

        # Validate required keys using the normalization module
        is_valid, missing_keys = validate_metadata_completeness(real_metadata)
        if not is_valid and missing_keys:
            # Check if it's just optional keys or truly required ones
            critical_missing = [k for k in missing_keys if k in [
                "range_pixel_spacing", "native_range_pixel_spacing", "azimuth_time_interval"
            ]]
            if critical_missing:
                # SCIENTIFIC FIX: This is a FATAL error, not a skip condition
                # Missing critical metadata means terrain correction would be scientifically invalid
                step_duration = time.time() - step_start
                processor.log_step(
                    step_number,
                    "Terrain Correction",
                    "error",
                    f"FATAL: Missing required metadata keys for terrain correction: {critical_missing}",
                    step_duration,
                )
                raise RuntimeError(
                    f"SCIENTIFIC MODE FAILURE: Cannot perform terrain correction.\n"
                    f"Missing required metadata keys: {critical_missing}\n"
                    f"This indicates SAFE annotation parsing failed to extract essential parameters.\n"
                    f"Check that the SLC product has valid annotation XMLs.\n"
                    f"Required keys: range_pixel_spacing, native_range_pixel_spacing, azimuth_time_interval"
                )
            else:
                # Non-critical keys missing, log warning but continue
                print(f"   ⚠️  Optional metadata keys missing (may affect quality): {missing_keys}")

        # Fail fast if subswaths are missing/empty – Range-Doppler needs per-subswath timing/geometry
        subswaths = real_metadata.get("subswaths")
        if not subswaths:
            step_duration = time.time() - step_start
            processor.log_step(
                step_number,
                "Terrain Correction",
                "error",
                "FATAL: Missing subswath metadata (subswaths is empty)",
                step_duration,
            )
            raise RuntimeError(
                "SCIENTIFIC MODE FAILURE: No subswath metadata available.\n"
                "Range-Doppler geocoding requires per-subswath timing and geometry for IW TOPSAR.\n"
                "Ensure SARDINE_REQUIRE_SUBSWATHS=1 is set before import and that SAFE parsing succeeded."
            )

        # Additional physics validation for timing consistency
        from .metadata_normalize import (
            validate_timing_consistency,
            validate_geometry_metadata,
            log_metadata_summary,
        )
        
        timing_ok, timing_warnings = validate_timing_consistency(real_metadata)
        for warning in timing_warnings:
            logger.debug(f"Timing validation: {warning}")
        
        geom_ok, geom_issues = validate_geometry_metadata(real_metadata)
        for issue in geom_issues:
            if "Missing" in issue:
                logger.warning(f"Geometry validation: {issue}")
            else:
                logger.debug(f"Geometry validation: {issue}")
        
        # Log metadata summary for debugging (only at debug level)
        log_metadata_summary(real_metadata, "terrain_correction")

        # === Extract burst timing if available ===
        # Burst timing key PRIORITY ORDER (tries in sequence until data found):
        # 1. "burst_timing_json" - JSON string format (how Rust stores it, needs deserialization)
        # 2. "burst_timings" - Direct list format (most common from cached reader)
        # 3. "burst_timing_records" - Alternative list format (some SAFE parsers use this)
        # 4. "burst_records" - Generic key (fallback)
        # 5. processor.burst_timing_records - Processor instance attribute (last resort)
        # 
        # Each burst timing record must contain these keys for terrain correction:
        # - subswath_id (str): Subswath identifier (IW1, IW2, IW3, etc.)
        # - burst_index (int): Burst number within subswath (0-based)
        # - azimuth_time_rel_orbit (float): CRITICAL! Time in seconds relative to orbit epoch
        #   Must be non-zero (except potentially first burst). Zero value indicates deserialization bug.
        # - first_line_global, last_line_global (int): Line extents in merged debursted grid
        import json
        burst_timing_json = None
        burst_timings = None
        
        # Try multiple possible keys for burst timing data
        if cached_meta:
            # Priority 1: Check for JSON string (how Rust stores it)
            raw_json = cached_meta.get("burst_timing_json")
            if raw_json:
                try:
                    burst_timings = json.loads(raw_json)
                    print(f"   📋 Parsed {len(burst_timings)} burst timing records from burst_timing_json")
                except Exception as e:
                    print(f"   ⚠️  Failed to parse burst_timing_json: {e}")
            
            # Priority 2-4: Fallback to list-based keys
            if not burst_timings:
                burst_timings = cached_meta.get("burst_timings")  # Priority 2 (most common)
            if not burst_timings:
                burst_timings = cached_meta.get("burst_timing_records")  # Priority 3
            if not burst_timings:
                burst_timings = cached_meta.get("burst_records")  # Priority 4
        
        # Priority 5: Check processor instance attribute as last resort
        if not burst_timings and hasattr(processor, 'burst_timing_records') and processor.burst_timing_records:
            burst_timings = processor.burst_timing_records
            print(f"   📋 Using {len(burst_timings)} burst timing records from processor")
        
        # Log burst timing status and serialize if we have valid data
        if burst_timings and isinstance(burst_timings, list) and len(burst_timings) > 0:
            print(f"   📋 Found {len(burst_timings)} burst timing records")
            try:
                burst_timing_json = json.dumps(burst_timings)
            except Exception as e:
                print(f"   ⚠️  Failed to serialize burst timings: {e}")
                burst_timing_json = None
        else:
            # Check if this is IW TOPSAR data that requires burst timings
            is_iw_topsar = False
            product_type = ""
            mode = ""
            if isinstance(processor.metadata, dict):
                product_type = processor.metadata.get("product_type", "").upper()
                mode = processor.metadata.get("mode", "").upper()
                is_iw_topsar = "IW" in mode or "IW" in product_type

            # Fallback heuristic: multiple subswaths imply IW-style merged data even if type/mode missing
            if not is_iw_topsar:
                def _has_multi_subswath(candidate: Any) -> bool:
                    return isinstance(candidate, dict) and len(candidate) >= 2

                subs_from_cache = None
                if cached_meta:
                    subs_from_cache = cached_meta.get("subswaths") or cached_meta.get("subswath_metadata")
                subs_from_meta = None
                if isinstance(processor.metadata, dict):
                    subs_from_meta = processor.metadata.get("subswaths")

                if _has_multi_subswath(subs_from_cache) or _has_multi_subswath(subs_from_meta):
                    is_iw_topsar = True

            if is_iw_topsar:
                # SCIENTIFIC FIX: This is a FATAL error, not a warning
                # Missing burst timing for merged IW TOPSAR causes 3× azimuth pixel errors
                # which completely invalidates geocoding accuracy. We cannot continue.
                step_duration = time.time() - step_start
                processor.log_step(
                    step_number,
                    "Terrain Correction",
                    "error",
                    "FATAL: No burst timing data available for IW TOPSAR product",
                    step_duration,
                )
                raise RuntimeError(
                    "SCIENTIFIC MODE FAILURE: No burst timing data available for IW TOPSAR product!\n"
                    "This WILL cause 3× azimuth pixel errors in geocoded output.\n"
                    "Burst timing is REQUIRED for accurate merged IW terrain correction.\n"
                    "Check that deburst stage populated burst_timing_records correctly.\n"
                    "This error cannot be bypassed as it would produce scientifically invalid results."
                )
            else:
                print("   ℹ️  No burst timing data (not IW TOPSAR or single burst)")

        # TASK 2: Choose bbox source (sar_footprint vs metadata)
        # NOTE: SAR footprint bbox requires subswath geometry in Rust (not yet implemented)
        # Current implementation lacks per-subswath slant_range_time, causing incorrect geocoding
        # In strict_science mode, requesting SAR footprint will fail hard (see _compute_sar_footprint_bbox)
        # In non-strict mode, it falls back to metadata bbox with warning
        processor_options = getattr(processor, 'options', {}) or {}
        use_sar_footprint = processor_options.get('use_sar_footprint_bbox', False)  # Default False until Rust support added
        margin_percent = processor_options.get('bbox_margin_percent', 2.0)  # Default 2% margin
        
        if use_sar_footprint:
            try:
                # Compute bbox from actual SAR image footprint (post-merge, post-multilook)
                footprint_result = _compute_sar_footprint_bbox(
                    processor=processor,
                    working_data=working_data,
                    cached_meta=cached_meta,
                    reader=reader,
                    orbit_times=orbit_times,
                    orbit_positions=orbit_positions,
                    orbit_velocities=orbit_velocities,
                    real_metadata=real_metadata,
                    burst_timing_json=burst_timing_json,
                    margin_percent=margin_percent
                )
                
                sar_bbox = footprint_result
                bbox_source = "sar_footprint"
                logger.info(
                    f"🎯 Using SAR footprint bbox: {sar_bbox} "
                    f"(margin={margin_percent}%)"
                )
                
                # Store provenance metadata
                if hasattr(processor, '_bbox_provenance'):
                    processor._bbox_provenance.update({
                        'bbox_used': sar_bbox,
                        'bbox_source': bbox_source,
                        'bbox_margin_percent': margin_percent,
                    })
                else:
                    processor._bbox_provenance = {
                        'bbox_used': sar_bbox,
                        'bbox_source': bbox_source,
                        'bbox_margin_percent': margin_percent,
                    }
                
            except Exception as e:
                # Fallback to metadata bbox if footprint computation fails
                logger.warning(
                    f"⚠️ SAR footprint bbox computation failed: {e}. "
                    f"Falling back to metadata bbox."
                )
                sar_bbox = [float(metadata_bbox[0]), float(metadata_bbox[1]), 
                           float(metadata_bbox[2]), float(metadata_bbox[3])]
                bbox_source = "metadata_fallback"
                if hasattr(processor, '_bbox_provenance'):
                    processor._bbox_provenance.update({
                        'bbox_used': sar_bbox,
                        'bbox_source': bbox_source,
                        'bbox_error': str(e),
                    })
                else:
                    processor._bbox_provenance = {
                        'bbox_used': sar_bbox,
                        'bbox_source': bbox_source,
                        'bbox_error': str(e),
                    }
        else:
            # Use metadata bbox (original approach)
            sar_bbox = [float(metadata_bbox[0]), float(metadata_bbox[1]), 
                       float(metadata_bbox[2]), float(metadata_bbox[3])]
            bbox_source = "metadata"
            logger.info(f"⚠️ Using metadata bbox: {sar_bbox}")
            if hasattr(processor, '_bbox_provenance'):
                processor._bbox_provenance.update({
                    'bbox_used': sar_bbox,
                    'bbox_source': bbox_source,
                })
            else:
                processor._bbox_provenance = {
                    'bbox_used': sar_bbox,
                    'bbox_source': bbox_source,
                }

        # Refine bbox using actual image dimensions and pixel spacing to avoid DEM overreach
        try:
            range_spacing_m = real_metadata.get("range_pixel_spacing")
            if range_spacing_m is not None and np.isfinite(range_spacing_m):
                range_spacing_m = float(range_spacing_m) * range_ml
            az_spacing_m = real_metadata.get("native_azimuth_pixel_spacing")
            if az_spacing_m is None:
                az_spacing_m = real_metadata.get("azimuth_pixel_spacing")
            if az_spacing_m is not None and np.isfinite(az_spacing_m):
                az_spacing_m = float(az_spacing_m) * azimuth_ml

            refined = None
            if range_spacing_m is not None and az_spacing_m is not None:
                refined = refine_geocoding_bbox(
                    sar_bbox,
                    actual_azimuth_lines,
                    actual_range_samples,
                    range_spacing_m,
                    az_spacing_m,
                )
            if refined is not None:
                refined_bbox, metrics = refined
                shrink_lat_pct = metrics.get("shrink_lat_pct", 0.0) * 100.0
                shrink_lon_pct = metrics.get("shrink_lon_pct", 0.0) * 100.0
                sar_bbox = refined_bbox
                logger.info(
                    f"📐 Refined bbox to match image footprint: {refined_bbox} "
                    f"(shrink lat {shrink_lat_pct:.1f}%, lon {shrink_lon_pct:.1f}%)"
                )
                if hasattr(processor, '_bbox_provenance'):
                    processor._bbox_provenance.update({
                        'bbox_used': sar_bbox,
                        'bbox_source': f"{bbox_source}_refined",
                        'bbox_refine_metrics': metrics,
                    })
        except Exception as bbox_refine_exc:
            logger.debug(f"BBox refinement skipped: {bbox_refine_exc}")

        # VALIDATION: Check SAR image has valid data before terrain correction
        total_pixels = working_data.size
        zeros = np.sum(working_data == 0)
        negatives = np.sum(working_data < 0)
        finite_positive = np.sum((working_data > 0) & np.isfinite(working_data))
        nan_count = np.sum(np.isnan(working_data))
        inf_count = np.sum(np.isinf(working_data))
        
        finite_pos = working_data[(working_data > 0) & np.isfinite(working_data)]
        if len(finite_pos) > 0:
            min_val, max_val, mean_val = np.min(finite_pos), np.max(finite_pos), np.mean(finite_pos)
        else:
            min_val = max_val = mean_val = np.nan
        
        valid_data_pct = (finite_positive / total_pixels * 100) if total_pixels > 0 else 0.0
        
        print(f"   🔍 SAR IMAGE STATISTICS (before terrain correction):")
        print(f"      Shape: {working_data.shape}")
        print(f"      Total pixels: {total_pixels:,}")
        print(f"      Zeros: {zeros:,} ({zeros/total_pixels*100:.1f}%)")
        print(f"      Negatives: {negatives:,} ({negatives/total_pixels*100:.1f}%)")
        print(f"      Finite positive: {finite_positive:,} ({valid_data_pct:.1f}%)")
        print(f"      NaN: {nan_count:,} ({nan_count/total_pixels*100:.1f}%)")
        print(f"      Inf: {inf_count:,} ({inf_count/total_pixels*100:.1f}%)")
        
        # Fail early if no valid data
        if finite_positive == 0:
            step_duration = time.time() - step_start
            processor.log_step(
                10,
                "Terrain Correction",
                "error",
                "No finite positive values in SAR data - cannot geocode",
                step_duration,
            )
            raise RuntimeError(
                "SCIENTIFIC MODE FAILURE: SAR image contains no finite positive values. "
                "Check calibration and upstream processing stages."
            )
        
        # Warn if very few valid pixels
        if valid_data_pct < 10.0:
            print(f"      ⚠️  WARNING: Only {valid_data_pct:.1f}% valid data - geocoding quality may be poor!")
            processor._record_validation(
                "terrain_correction",
                "input_valid_data_percentage",
                value=valid_data_pct,
                expected=">=10%",
                passed=False,
                severity="warning",
                message=f"Low valid data percentage ({valid_data_pct:.1f}%) before terrain correction",
            )
        else:
            print(f"      Finite positive range: [{min_val:.6e}, {max_val:.6e}], mean: {mean_val:.6e}")
        
        # Perform range-Doppler terrain correction via Rust binding
        # CRITICAL: Filter real_metadata to only float values - Rust expects HashMap<String, f64>
        # Complex objects like 'subswaths' will cause "must be real number, not str/dict" errors
        # 
        # NOTE: Subswath metadata is NOT passed through real_metadata. Instead, Rust reads
        # subswath geometry directly from cached_metadata.sub_swaths via the slc_reader parameter.
        # This ensures Python and Rust use the same subswath source from the cached metadata.
        # --- Metadata Bridge: Python dict -> Rust HashMap<String, f64> ---
        # The Rust terrain_correction binding accepts only f64 values.
        # Complex objects (lists, dicts, strings) are passed via the slc_reader's
        # cached metadata instead.  We validate that the critical scalar keys
        # survive the conversion so that silent data loss is impossible.
        REQUIRED_FLOAT_KEYS = {
            "native_range_pixel_spacing",
            "slant_range_time",
            "prf",
        }

        float_metadata = {}
        skipped_keys = []
        for key, value in real_metadata.items():
            try:
                float_metadata[key] = float(value)
            except (TypeError, ValueError):
                skipped_keys.append(key)

        # Validate that every required scalar key made it through
        missing_required = REQUIRED_FLOAT_KEYS - set(float_metadata.keys())
        if missing_required:
            strict_science = bool(getattr(processor, 'strict_science', False))
            msg = (
                f"METADATA CONTRACT VIOLATION: Required float keys missing after "
                f"bridge conversion: {sorted(missing_required)}.\n"
                f"  Available float keys: {sorted(float_metadata.keys())}\n"
                f"  Skipped (non-float) keys: {sorted(skipped_keys)}\n"
                f"  These keys are required by the Rust terrain_correction binding."
            )
            if strict_science:
                raise RuntimeError(
                    f"STRICT SCIENCE MODE FAILURE: {msg}\n\n"
                    "FIX: Ensure SAFE annotation parsing populates these keys "
                    "in real_metadata before terrain correction."
                )
            logger.error(f"⚠️  {msg}")

        if skipped_keys:
            logger.debug(f"Filtered out non-float metadata keys for Rust bridge: {skipped_keys}")
            # Validate that subswaths are available in cached metadata if this is IW data
            if 'subswaths' in skipped_keys and cached_meta:
                subswaths_in_cache = cached_meta.get('subswaths')
                if subswaths_in_cache:
                    logger.debug("Subswath metadata available in cached_metadata (not passed via real_metadata)")
                else:
                    logger.warning(
                        "⚠️  Subswath metadata not found in cached_metadata. "
                        "This may cause incorrect range-Doppler calculations for merged IW data."
                    )
        
        # Determine RTC mode - use processor setting or default to area projection
        rtc_mode = getattr(processor, "rtc_mode", "area")  # Default: scientifically correct
        output_lia = getattr(processor, "output_lia", False)
        output_masks = getattr(processor, "output_masks", False)
        
        # Interpolation method (required by Rust binding; default bilinear)
        interpolation_method = (
            (getattr(processor, "options", {}) or {}).get("interpolation_method")
            or (getattr(processor, "options", {}) or {}).get("interpolation")
            or "bilinear"
        )
        interpolation_method = str(interpolation_method).strip().lower()
        allowed_methods = {"nearest", "bilinear", "bicubic", "sinc", "lanczos"}
        if interpolation_method not in allowed_methods:
            logger.warning(
                f"Unknown interpolation method {interpolation_method!r}; falling back to 'bilinear'. "
                f"Allowed: {sorted(allowed_methods)}"
            )
            interpolation_method = "bilinear"

        print(f"   🎯 RTC mode: {rtc_mode} (integrated with geocoding)")
        print(f"   🔧 Interpolation: {interpolation_method}")
        
        try:
            tc_result = sardine.terrain_correction(
                working_data,
                sar_bbox,
                orbit_times,
                orbit_positions,
                orbit_velocities,
                str(dem_cache_dir),
                float(processor.target_resolution),
                float_metadata,  # Use filtered dict with only float values
                reader,
                interpolation_method,  # interpolation_method (required)
                burst_timing_json,
                rtc_mode,     # RTC mode: "area" (default), "cosine", or "none"
                output_lia,   # Whether to output LIA array
                output_masks, # Whether to output shadow/layover masks
            )
        except TypeError as type_exc:
            # TypeError indicates Python-Rust binding mismatch (wrong signature, invalid metadata format, etc.)
            # In strict science mode, this is a hard failure because skipping geocoding
            # produces invalid output (data remains in SAR geometry, not geocoded)
            
            strict_science = bool(getattr(processor, 'strict_science', False))
            if strict_science:
                raise RuntimeError(
                    "STRICT SCIENCE MODE FAILURE: Terrain correction binding error.\n"
                    f"TypeError from Rust binding: {type_exc}\n\n"
                    "This likely indicates:\n"
                    "  1. Metadata format mismatch (float_metadata contains invalid types), OR\n"
                    "  2. Required parameter missing/wrong type in terrain_correction call, OR\n"
                    "  3. Rust binding signature changed but Python code not updated\n\n"
                    "Skipping terrain correction would produce SAR-geometry output without geocoding,\n"
                    "which is scientifically invalid for most applications.\n\n"
                    "To proceed:\n"
                    "  1. Check real_metadata normalization (should contain only float values)\n"
                    "  2. Verify all required parameters match Rust binding signature\n"
                    "  3. Set processor.strict_science = False to allow graceful skip (not recommended)\n"
                ) from type_exc
            
            # Non-strict mode: log prominent warning and skip gracefully
            logger.error(
                "⚠️  TERRAIN CORRECTION BINDING ERROR - SKIPPING GEOCODING\n"
                f"   TypeError: {type_exc}\n"
                "   ⚠️  OUTPUT WILL REMAIN IN SAR GEOMETRY (NOT GEOCODED)\n"
                "   This is likely due to metadata format mismatch or binding signature changes.\n"
                "   Enable strict_science mode to prevent this fallback."
            )
            step_duration = time.time() - step_start
            processor.log_step(
                step_number,
                "Terrain Correction",
                "skipped",
                f"Binding error (see logs): {type_exc}",
                step_duration,
            )
            return

        if isinstance(tc_result, dict):
            geocoded = tc_result.get("data")
            geo_transform = tc_result.get("geo_transform")
            projection = tc_result.get("projection")
            
            # Handle RTC outputs if present
            lia_array = tc_result.get("local_incidence_angle")
            shadow_mask = tc_result.get("shadow_mask")
            layover_mask = tc_result.get("layover_mask")
            rtc_mode_used = tc_result.get("rtc_mode", "unknown")
            
            if lia_array is not None:
                processor._lia_array = lia_array
                print(f"   📐 Local incidence angle array stored ({lia_array.shape})")
            if shadow_mask is not None:
                processor._shadow_mask = shadow_mask
                print(f"   🌑 Shadow mask stored ({shadow_mask.shape})")
            if layover_mask is not None:
                processor._layover_mask = layover_mask
                print(f"   📊 Layover mask stored ({layover_mask.shape})")
            print(f"   🎯 RTC mode applied: {rtc_mode_used}")
        else:
            geocoded = tc_result
            geo_transform = None
            projection = None

        if not isinstance(geocoded, np.ndarray):
            raise RuntimeError("SCIENTIFIC MODE FAILURE: Terrain correction did not return numpy array")

        # CRITICAL FIX: Ensure geotransform is north-up, east-positive orientation
        # If the geotransform comes out flipped, flip the data to match the expected orientation.
        # This prevents mirrored/flipped output that would be scientifically incorrect.
        if geo_transform is not None:
            try:
                geocoded, geo_transform, orientation_msg = _apply_north_up_orientation(
                    geocoded, geo_transform
                )
                if orientation_msg:
                    print(f"   {orientation_msg}")
            except Exception as orient_exc:
                logger.warning(f"Could not validate/fix geotransform orientation: {orient_exc}")

        # CRITICAL FIX: Validate coverage to prevent silent 0% valid pixels failure
        # This addresses the critical bug where terrain correction could "succeed"
        # while producing all-NaN output. Account for explicit layover/shadow masks
        # when present so expected invalid areas do not trigger false failures.
        valid_mask = np.isfinite(geocoded) & (geocoded > 0)

        # Areas blocked by layover/shadow masks are excluded from the coverage denominator
        available_mask = np.ones_like(valid_mask, dtype=bool)
        masked_out = {}
        for mask_name, mask_array in (("shadow", shadow_mask), ("layover", layover_mask)):
            if mask_array is None:
                continue
            try:
                mask_np = np.asarray(mask_array)
                if mask_np.shape == geocoded.shape:
                    blocked = mask_np > 0
                    blocked_count = int(np.count_nonzero(blocked))
                    if blocked_count > 0:
                        available_mask &= ~blocked
                        masked_out[mask_name] = blocked_count
            except Exception:
                pass

        effective_total_pixels = int(np.count_nonzero(available_mask))
        coverage_denominator = effective_total_pixels if effective_total_pixels > 0 else geocoded.size

        valid_pixel_count = int(np.count_nonzero(valid_mask & available_mask))
        valid_percentage = (
            (valid_pixel_count / coverage_denominator) * 100.0
            if coverage_denominator > 0
            else 0.0
        )
        if masked_out:
            masked_pct = sum(masked_out.values()) / float(geocoded.size) * 100.0
            logger.info(
                f"🌓 Masked out shadow/layover pixels: {masked_out} "
                f"({masked_pct:.2f}% of grid)"
            )
        
        # Use module-level coverage thresholds for validation
        if valid_percentage < CRITICAL_COVERAGE_THRESHOLD:
            step_duration = time.time() - step_start
            from sardine.processors.backscatter.processor import STEP_NUMBERS
            processor.log_step(
                STEP_NUMBERS["Terrain Correction"],
                "Terrain Correction",
                "error",
                f"COVERAGE FAILURE: Only {valid_percentage:.2f}% valid pixels ({valid_pixel_count}/{coverage_denominator})",
                step_duration,
            )
            raise RuntimeError(
                f"SCIENTIFIC MODE FAILURE: Terrain correction produced only {valid_percentage:.2f}% "
                f"valid pixels ({valid_pixel_count}/{coverage_denominator}). This indicates a fundamental error in "
                f"the geocoding process (timing mismatch, DEM bbox error, or invalid orbit data). "
                f"Check native_range_pixel_spacing vs range_pixel_spacing, burst_timings keys, and orbit vector count."
            )
        elif valid_percentage < WARNING_COVERAGE_THRESHOLD:
            logger.warning(
                f"⚠️  Low coverage in terrain correction: {valid_percentage:.1f}% valid pixels "
                f"({valid_pixel_count}/{coverage_denominator}). This may indicate partial geocoding issues."
            )

        if geo_transform is not None:
            processor.geo_transform = geo_transform
            
            # CRITICAL FIX: Compute product_bbox from actual geocoded image geotransform
            # This ensures we use the correct bbox from the geocoded image, not metadata
            rows, cols = geocoded.shape
            
            # Get CRS/EPSG code for proper bbox computation and validation
            output_epsg = getattr(processor, 'output_epsg', None)
            if output_epsg is None and projection is not None:
                # Try to extract EPSG from projection string (e.g., "EPSG:32632")
                import re
                match = re.search(r'EPSG[:\s]+(\d+)', str(projection), re.IGNORECASE)
                if match:
                    output_epsg = int(match.group(1))
            
            computed_bbox = compute_bbox_from_transform(geo_transform, rows, cols, epsg=output_epsg)
            
            if computed_bbox is not None:
                is_valid, error = _validate_bbox(computed_bbox, "geocoded_image_geotransform", epsg=output_epsg)
                if is_valid:
                    processor.product_bbox = computed_bbox
                    logger.info(
                        f"✅ Updated product_bbox from geocoded image: "
                        f"[{computed_bbox[0]:.6f}, {computed_bbox[1]:.6f}, "
                        f"{computed_bbox[2]:.6f}, {computed_bbox[3]:.6f}]"
                    )
                else:
                    logger.warning(
                        f"⚠️  Computed bbox from geotransform failed validation: {error}. "
                        f"Keeping existing product_bbox."
                    )
            else:
                logger.warning(
                    "⚠️  Could not compute bbox from geotransform. "
                    "Keeping existing product_bbox."
                )
        
        if projection is not None:
            processor.coordinate_system = projection

        geocoded_array = np.asarray(geocoded, dtype=np.float32)
        
        # Apply aggressive cropping to remove empty borders
        # This improves valid pixel percentage by removing large empty regions
        if geo_transform is not None and processor.product_bbox is not None:
            crop_result = crop_geocoded_output(geocoded_array, geo_transform, processor.product_bbox)
            if crop_result is not None:
                logger.info(
                    f"✂️  Cropped geocoded output: {geocoded_array.shape} → {crop_result['data'].shape} "
                    f"(valid pixels: {crop_result['valid_pct']:.1f}%)"
                )
                geocoded_array = crop_result['data']
                processor.geo_transform = crop_result['transform']
                processor.product_bbox = crop_result['bbox']
                logger.info(
                    f"✅ Updated geotransform and bbox after cropping"
                )
        
        # Diagnostic logging: After terrain correction
        processor._log_diagnostic_statistics(
            "After Terrain Correction",
            geocoded_array,
            context="geocoded"
        )
        
        # TASK 3: Log provenance metadata for bbox source and coverage metrics
        if hasattr(processor, '_bbox_provenance'):
            # Compute grid/input pixel ratio for validation
            # working_data is the pre-geocode SAR array; use it to size the input grid
            input_pixels = working_data.size if working_data is not None else 0
            grid_pixels = geocoded_array.size
            grid_input_ratio = grid_pixels / input_pixels if input_pixels > 0 else 0
            effective_grid_pixels = coverage_denominator
            
            # Update provenance with actual grid and coverage metrics
            processor._bbox_provenance.update({
                'input_pixels': input_pixels,
                'expected_grid_pixels': grid_pixels,
                'effective_grid_pixels': effective_grid_pixels,
                'actual_valid_pixels': valid_pixel_count,
                'coverage_percent': valid_percentage,
                'grid_input_ratio': grid_input_ratio,
            })
            
            # Log comprehensive bbox provenance
            bbox_source = processor._bbox_provenance.get('bbox_source', 'unknown')
            metadata_bbox = processor._bbox_provenance.get('bbox_metadata', None)
            used_bbox = processor._bbox_provenance.get('bbox_used', None)
            
            logger.info(
                f"📊 BBox Provenance - Source: {bbox_source}, "
                f"Metadata: {metadata_bbox}, Used: {used_bbox}"
            )
            logger.info(
                f"📊 Grid Metrics - Input: {input_pixels:,} px, "
                f"Grid: {grid_pixels:,} px, Valid: {valid_pixel_count:,} px ({valid_percentage:.2f}%), "
                f"Grid/Input: {grid_input_ratio:.3f}"
            )
            
            # TASK 5: Add guardrails for low coverage
            if valid_percentage < 80.0:
                logger.warning(
                    f"⚠️  COVERAGE WARNING: Valid pixel percentage ({valid_percentage:.2f}%) is below 80%. "
                    f"This may indicate bbox or geocoding issues."
                )
                processor._record_validation(
                    "terrain_correction",
                    "coverage_percentage",
                    value=valid_percentage,
                    expected=f">={CRITICAL_COVERAGE_THRESHOLD}% (warn <80%)",
                    # Warn but do not fail the validation gate unless below the critical threshold.
                    passed=valid_percentage >= CRITICAL_COVERAGE_THRESHOLD,
                    severity="warning",
                    message=(
                        f"Coverage {valid_percentage:.2f}% below recommended 80% threshold; "
                        f"critical threshold is {CRITICAL_COVERAGE_THRESHOLD}%"
                    ),
                )
            
            # Check if scientific-strict mode is enabled
            scientific_strict = processor_options.get('scientific_strict', False)
            if scientific_strict and valid_percentage < 75.0:
                step_duration = time.time() - step_start
                from sardine.processors.backscatter.processor import STEP_NUMBERS
                processor.log_step(
                    STEP_NUMBERS["Terrain Correction"],
                    "Terrain Correction",
                    "error",
                    f"SCIENTIFIC-STRICT MODE: Coverage {valid_percentage:.2f}% < 75% threshold",
                    step_duration,
                )
                raise RuntimeError(
                    f"SCIENTIFIC-STRICT MODE FAILURE: Coverage {valid_percentage:.2f}% is below the "
                    f"75% critical threshold required in scientific-strict mode. "
                    f"BBox source: {bbox_source}, Grid/Input ratio: {grid_input_ratio:.3f}"
                )
        
        processor._working_data = geocoded_array

        step_duration = time.time() - step_start
        from sardine.processors.backscatter.processor import STEP_NUMBERS
        processor.log_step(
            STEP_NUMBERS["Terrain Correction"],
            "Terrain Correction",
            "success",
            f"Geocoded: {geocoded.shape}, resolution ~{processor.target_resolution}m",
            step_duration,
        )
        processor._record_stage_timing("Terrain Correction", step_duration)

    except Exception as exc:
        step_duration = time.time() - step_start
        from sardine.processors.backscatter.processor import STEP_NUMBERS
        processor.log_step(STEP_NUMBERS["Terrain Correction"], "Terrain Correction", "error", f"Terrain correction failed: {exc}", step_duration)
        raise RuntimeError(f"SCIENTIFIC MODE FAILURE: Terrain correction failed. Error: {exc}")
