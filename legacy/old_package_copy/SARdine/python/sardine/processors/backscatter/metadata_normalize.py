"""
Metadata normalization for Python → Rust PyO3 bridge compatibility.

This module ensures that metadata keys from Python processing stages
match exactly what the Rust terrain correction code expects.

Key Mappings:
    Python Key               → Rust Expected Key
    ─────────────────────────────────────────────
    subswath_geometry        → subswaths
    subswath_metadata        → subswaths
    burst_timing_records     → burst_timings
    total_azimuth_lines      → number_of_lines
    num_lines                → number_of_lines
"""

from typing import Dict, Any, List, Tuple, Optional
from datetime import datetime
import logging

logger = logging.getLogger(__name__)

# ============================================================================
# Key Mapping Contract: Python → Rust
# ============================================================================

KEY_MAPPINGS: Dict[str, str] = {
    # Subswath geometry
    "subswath_geometry": "subswaths",
    "subswath_metadata": "subswaths",
    "iw_subswaths": "subswaths",
    
    # Burst timing
    "burst_timing_records": "burst_timings",
    "burst_records": "burst_timings",
    
    # Line counts
    "total_azimuth_lines": "number_of_lines",
    "num_lines": "number_of_lines",
    "azimuth_lines": "number_of_lines",
    
    # Sample counts
    "total_range_samples": "number_of_samples",
    "num_samples": "number_of_samples",
    "range_samples": "number_of_samples",
}

# Required keys for Rust terrain correction
REQUIRED_KEYS_TERRAIN_CORRECTION: List[str] = [
    "subswaths",
    "burst_timings",
    "number_of_lines",
    "native_range_pixel_spacing",
    "native_azimuth_pixel_spacing",
    "range_pixel_spacing",
    "azimuth_pixel_spacing",
    "slant_range_time",
    "azimuth_time_interval",
]

# Minimum orbit vectors required for accurate terrain correction
# 
# RATIONALE: 
# - Python uses 10 vectors as minimum for scientific processing (ensures adequate temporal coverage)
# - Rust uses 3 vectors minimum for interpolation (cubic spline requires 3+ points)
# - Precise orbit files (POEORB) typically have 9000+ vectors (10s intervals over 24h)
# - SAFE annotation orbit has ~17 vectors (10s intervals around acquisition)
# - 10 vectors provides ~100s of coverage, sufficient for most scenes
# - For very long scenes or edge cases, more vectors are better
MINIMUM_ORBIT_VECTORS = 10

# IW mode native pixel spacings (known good values for Sentinel-1)
IW_NATIVE_RANGE_SPACING_M = 2.329562    # meters, slant range
IW_NATIVE_AZIMUTH_SPACING_M = 13.968    # meters, approximate


class MissingMetadataError(ValueError):
    """Raised when required metadata is missing and strict mode is enabled."""
    pass


# ============================================================================
# Main Normalization Function
# ============================================================================

def normalize_metadata_for_rust(
    metadata: Dict[str, Any],
    is_multilooked: bool = False,
    range_looks: int = 1,
    azimuth_looks: int = 1,
    log_changes: bool = True,
    strict: bool = True,
) -> Dict[str, Any]:
    """
    Normalize metadata keys and values for Rust PyO3 bridge.
    
    This function ensures that:
    1. Key names match Rust expectations (see KEY_MAPPINGS)
    2. Native pixel spacing values are present
    3. Orbit times are in ISO8601 format
    4. Nested structures are flattened appropriately
    
    By default (`strict=True`), this function will raise MissingMetadataError
    if required values are missing, ensuring upstream parsing problems are not
    silently ignored. When `strict=False`, IW defaults are used as fallbacks
    for missing or invalid spacing values to allow robust processing.
    
    Args:
        metadata: Raw metadata dictionary from Python processing
        is_multilooked: Whether data has been multilooked
        range_looks: Number of range looks applied (default 1)
        azimuth_looks: Number of azimuth looks applied (default 1)
        log_changes: Whether to log key renames and additions
        
    Returns:
        Normalized metadata dictionary with correct keys for Rust
        
    Raises:
        MissingMetadataError: If required values are missing
        
    Example:
        >>> raw = {"subswath_geometry": {...}, "total_azimuth_lines": 1500}
        >>> normalized = normalize_metadata_for_rust(raw)
        >>> assert "subswaths" in normalized
        >>> assert "number_of_lines" in normalized
    """
    if metadata is None:
        logger.warning("normalize_metadata_for_rust received None, returning empty dict")
        return {}
    
    normalized = dict(metadata)
    changes_made = []
    
    # Step 1: Rename keys to match Rust expectations
    for old_key, new_key in KEY_MAPPINGS.items():
        if old_key in normalized and new_key not in normalized:
            normalized[new_key] = normalized.pop(old_key)
            changes_made.append(f"renamed '{old_key}' → '{new_key}'")
    
    # Step 1.5: Deserialize burst_timing_json if present
    if "burst_timing_json" in normalized and "burst_timings" not in normalized:
        import json
        try:
            parsed = json.loads(normalized.pop("burst_timing_json"))
            normalized["burst_timings"] = parsed
            changes_made.append(f"deserialized burst_timing_json → burst_timings ({len(parsed)} records)")
        except (json.JSONDecodeError, TypeError) as e:
            logger.warning(f"Failed to deserialize burst_timing_json: {e}")
    
    # Step 2: Ensure native spacing keys exist
    spacing_changes = _ensure_native_spacing(
        normalized, is_multilooked, range_looks, azimuth_looks, strict
    )
    changes_made.extend(spacing_changes)
    
    # Step 3: Normalize orbit time formats to ISO8601
    time_changes = _normalize_orbit_times(normalized)
    changes_made.extend(time_changes)
    
    # Step 4: Flatten nested subswath structures
    flatten_changes = _flatten_subswaths(normalized)
    changes_made.extend(flatten_changes)
    
    # Step 5: Ensure number_of_lines from subswath data if missing
    lines_changes = _ensure_number_of_lines(normalized)
    changes_made.extend(lines_changes)
    
    if log_changes and changes_made:
        logger.info(f"Metadata normalization: {', '.join(changes_made)}")
    
    return normalized


# ============================================================================
# Validation Functions
# ============================================================================

def validate_metadata_completeness(
    metadata: Dict[str, Any],
    required_keys: Optional[List[str]] = None,
) -> Tuple[bool, List[str]]:
    """
    Validate that all required metadata keys are present.
    
    Args:
        metadata: Metadata dictionary to validate
        required_keys: List of required keys (defaults to REQUIRED_KEYS_TERRAIN_CORRECTION)
    
    Returns:
        Tuple of (is_valid, list_of_missing_keys)
        
    Example:
        >>> is_valid, missing = validate_metadata_completeness({"range_pixel_spacing": 2.33})
        >>> print(f"Valid: {is_valid}, Missing: {missing}")
    """
    if required_keys is None:
        required_keys = REQUIRED_KEYS_TERRAIN_CORRECTION
    
    missing = []
    for key in required_keys:
        if key not in metadata or metadata[key] is None:
            missing.append(key)
    
    return (len(missing) == 0, missing)


def validate_orbit_vector_count(
    orbit_times: Optional[List[Any]] = None,
    orbit_positions: Optional[List[Any]] = None,
    metadata: Optional[Dict[str, Any]] = None,
    minimum_vectors: int = MINIMUM_ORBIT_VECTORS,
) -> Tuple[bool, int, str]:
    """
    Validate that sufficient orbit vectors are available.
    
    Per project requirements, we need at least 10 orbit vectors for
    accurate terrain correction.
    
    Args:
        orbit_times: List of orbit times (optional)
        orbit_positions: List of orbit positions (optional)
        metadata: Metadata dict that may contain orbit_state_vectors (optional)
        minimum_vectors: Minimum required vectors (default 10)
        
    Returns:
        Tuple of (is_valid, vector_count, message)
    """
    vector_count = 0
    
    # Try to get count from provided lists
    if orbit_times is not None:
        vector_count = len(orbit_times)
    elif orbit_positions is not None:
        vector_count = len(orbit_positions)
    elif metadata is not None:
        # Try to extract from metadata
        if "orbit_state_vectors" in metadata:
            vectors = metadata["orbit_state_vectors"]
            if isinstance(vectors, list):
                vector_count = len(vectors)
        elif "orbit_vectors_count" in metadata:
            vector_count = int(metadata["orbit_vectors_count"])
    
    if vector_count < minimum_vectors:
        msg = f"Insufficient orbit data: {vector_count} vectors, need >= {minimum_vectors}"
        return (False, vector_count, msg)
    
    # Cross-validate if both lists provided
    if orbit_times is not None and orbit_positions is not None:
        if len(orbit_times) != len(orbit_positions):
            msg = f"Orbit data mismatch: {len(orbit_times)} times but {len(orbit_positions)} positions"
            return (False, vector_count, msg)
    
    return (True, vector_count, f"Orbit data OK: {vector_count} vectors")


# ============================================================================
# Helper Functions
# ============================================================================

def _ensure_native_spacing(
    metadata: Dict[str, Any],
    is_multilooked: bool,
    range_looks: int,
    azimuth_looks: int,
    strict: bool,
) -> List[str]:
    """
    Ensure native_*_pixel_spacing keys exist.
    
    Native spacing refers to the original single-look spacing before any
    multilooking. This is critical for Range-Doppler geocoding.
    
    Returns list of changes made.
    
    Raises:
        MissingMetadataError: If spacing data is unavailable
    """
    changes = []
    
    # Native range spacing
    if "native_range_pixel_spacing" not in metadata:
        if "range_pixel_spacing" in metadata:
            try:
                current_spacing = float(metadata["range_pixel_spacing"])
                if is_multilooked and range_looks > 1:
                    # Current spacing is multilooked, compute native
                    native = current_spacing / range_looks
                    metadata["native_range_pixel_spacing"] = native
                    changes.append(f"computed native_range_pixel_spacing={native:.4f}")
                else:
                    # Assume current spacing IS native
                    metadata["native_range_pixel_spacing"] = current_spacing
                    changes.append(f"set native_range_pixel_spacing={current_spacing:.4f}")
            except (ValueError, TypeError) as e:
                if strict:
                    raise MissingMetadataError(
                        f"Cannot parse range_pixel_spacing: {e}. "
                        "This indicates a parsing failure upstream. "
                        "Fix the metadata extraction."
                    )
                else:
                    metadata["native_range_pixel_spacing"] = IW_NATIVE_RANGE_SPACING_M
                    changes.append(
                        f"fallback native_range_pixel_spacing={IW_NATIVE_RANGE_SPACING_M:.4f} (IW default)"
                    )
        else:
            if strict:
                raise MissingMetadataError(
                    "Missing range_pixel_spacing in metadata. "
                    "This is required for geocoding. Check SAFE manifest parsing."
                )
            else:
                metadata["native_range_pixel_spacing"] = IW_NATIVE_RANGE_SPACING_M
                changes.append(
                    f"fallback native_range_pixel_spacing={IW_NATIVE_RANGE_SPACING_M:.4f} (IW default)"
                )
    
    # Native azimuth spacing
    if "native_azimuth_pixel_spacing" not in metadata:
        if "azimuth_pixel_spacing" in metadata:
            try:
                current_spacing = float(metadata["azimuth_pixel_spacing"])
                if is_multilooked and azimuth_looks > 1:
                    native = current_spacing / azimuth_looks
                    metadata["native_azimuth_pixel_spacing"] = native
                    changes.append(f"computed native_azimuth_pixel_spacing={native:.4f}")
                else:
                    metadata["native_azimuth_pixel_spacing"] = current_spacing
                    changes.append(f"set native_azimuth_pixel_spacing={current_spacing:.4f}")
            except (ValueError, TypeError) as e:
                if strict:
                    raise MissingMetadataError(
                        f"Cannot parse azimuth_pixel_spacing: {e}. "
                        "This indicates a parsing failure upstream."
                    )
                else:
                    metadata["native_azimuth_pixel_spacing"] = IW_NATIVE_AZIMUTH_SPACING_M
                    changes.append(
                        f"fallback native_azimuth_pixel_spacing={IW_NATIVE_AZIMUTH_SPACING_M:.4f} (IW default)"
                    )
        else:
            if strict:
                raise MissingMetadataError(
                    "Missing azimuth_pixel_spacing in metadata. "
                    "This is required for geocoding. Check SAFE manifest parsing."
                )
            else:
                metadata["native_azimuth_pixel_spacing"] = IW_NATIVE_AZIMUTH_SPACING_M
                changes.append(
                    f"fallback native_azimuth_pixel_spacing={IW_NATIVE_AZIMUTH_SPACING_M:.4f} (IW default)"
                )
    
    return changes


def _normalize_orbit_times(metadata: Dict[str, Any]) -> List[str]:
    """
    Convert orbit times to ISO8601 strings if they are Unix timestamps.
    
    Special case: keep `orbit_ref_epoch_utc` numeric for Rust geocoder while
    adding an `orbit_ref_epoch_iso` companion string for logging/QA.
    
    Returns list of changes made.
    """
    changes = []
    time_keys = [
        "orbit_ref_epoch_utc",
        "orbit_time",
        "state_vector_time",
        "first_line_time",
        "last_line_time",
        "azimuth_start_time",
        "azimuth_stop_time",
    ]
    
    for key in time_keys:
        if key in metadata:
            val = metadata[key]
            if isinstance(val, (int, float)):
                # Keep orbit_ref_epoch_utc numeric for Rust; add ISO copy instead of replacing
                if key == "orbit_ref_epoch_utc":
                    try:
                        iso_str = datetime.utcfromtimestamp(val).strftime(
                            "%Y-%m-%dT%H:%M:%S.%fZ"
                        )
                        metadata["orbit_ref_epoch_iso"] = iso_str
                        changes.append("added orbit_ref_epoch_iso from numeric epoch")
                    except (ValueError, OSError) as e:
                        logger.warning(f"Failed to format orbit_ref_epoch_utc={val} to ISO8601: {e}")
                    # Do not overwrite the numeric value
                    continue

                # For other time fields, continue converting to ISO string
                try:
                    iso_str = datetime.utcfromtimestamp(val).strftime(
                        "%Y-%m-%dT%H:%M:%S.%fZ"
                    )
                    metadata[key] = iso_str
                    changes.append(f"converted {key} to ISO8601")
                except (ValueError, OSError) as e:
                    logger.warning(f"Failed to convert {key}={val} to ISO8601: {e}")
            # If already string, assume correct format
    
    # Also handle orbit_state_vectors list if present
    if "orbit_state_vectors" in metadata:
        vectors = metadata["orbit_state_vectors"]
        converted_count = 0
        if isinstance(vectors, list):
            for vec in vectors:
                if isinstance(vec, dict) and "time" in vec:
                    t = vec["time"]
                    if isinstance(t, (int, float)):
                        try:
                            vec["time"] = datetime.utcfromtimestamp(t).strftime(
                                "%Y-%m-%dT%H:%M:%S.%fZ"
                            )
                            converted_count += 1
                        except (ValueError, OSError):
                            pass
        if converted_count > 0:
            changes.append(f"converted {converted_count} orbit vector times to ISO8601")
    
    return changes


def _flatten_subswaths(metadata: Dict[str, Any]) -> List[str]:
    """
    Flatten nested subswath structures.
    
    Handles cases where subswaths are double-nested like:
        {"subswaths": {"subswath_metadata": {"iw1": {...}}}}
    
    Should become:
        {"subswaths": {"iw1": {...}}}
        
    Returns list of changes made.
    """
    changes = []
    
    if "subswaths" in metadata:
        subswaths = metadata["subswaths"]
        
        if isinstance(subswaths, dict):
            # Handle double-nested case: {"subswath_metadata": {...}}
            if "subswath_metadata" in subswaths and len(subswaths) == 1:
                metadata["subswaths"] = subswaths["subswath_metadata"]
                changes.append("flattened nested subswath_metadata")
            
            # Handle other potential nesting patterns
            elif "subswath_geometry" in subswaths and len(subswaths) == 1:
                metadata["subswaths"] = subswaths["subswath_geometry"]
                changes.append("flattened nested subswath_geometry")
    
    return changes


def _ensure_number_of_lines(metadata: Dict[str, Any]) -> List[str]:
    """
    Ensure number_of_lines is present, extracting from subswaths if needed.
    
    Returns list of changes made.
    """
    changes = []
    
    if "number_of_lines" not in metadata or metadata["number_of_lines"] is None:
        # Try to extract from subswaths
        if "subswaths" in metadata and isinstance(metadata["subswaths"], dict):
            subswaths = metadata["subswaths"]
            # Look for lines in any subswath
            for sw_key, sw_data in subswaths.items():
                if isinstance(sw_data, dict):
                    for line_key in ["number_of_lines", "lines", "num_lines", "azimuth_lines"]:
                        if line_key in sw_data:
                            try:
                                metadata["number_of_lines"] = int(sw_data[line_key])
                                changes.append(f"extracted number_of_lines from subswaths[{sw_key}]")
                                return changes
                            except (ValueError, TypeError):
                                pass
        
        # Try image dimensions
        for dim_key in ["image_lines", "azimuth_size", "height"]:
            if dim_key in metadata:
                try:
                    metadata["number_of_lines"] = int(metadata[dim_key])
                    changes.append(f"set number_of_lines from {dim_key}")
                    return changes
                except (ValueError, TypeError):
                    pass
    
    return changes


# ============================================================================
# Convenience Functions
# ============================================================================

def prepare_metadata_for_terrain_correction(
    processor_metadata: Dict[str, Any],
    geocoding_cache: Optional[Dict[str, Any]] = None,
    is_multilooked: bool = False,
    range_looks: int = 1,
    azimuth_looks: int = 1,
) -> Tuple[Dict[str, Any], bool, List[str]]:
    """
    High-level function to prepare metadata for Rust terrain correction.
    
    Combines metadata from processor and geocoding cache, normalizes keys,
    and validates completeness.
    
    Args:
        processor_metadata: Main metadata from BackscatterProcessor
        geocoding_cache: Cached geometry data from IW split stage
        is_multilooked: Whether data has been multilooked
        range_looks: Number of range looks
        azimuth_looks: Number of azimuth looks
        
    Returns:
        Tuple of (prepared_metadata, is_valid, missing_keys)
    """
    # Start with processor metadata
    combined = dict(processor_metadata) if processor_metadata else {}
    
    # Merge in geocoding cache
    if geocoding_cache:
        for key, value in geocoding_cache.items():
            if key not in combined or combined[key] is None:
                combined[key] = value
    
    # Normalize
    normalized = normalize_metadata_for_rust(
        combined,
        is_multilooked=is_multilooked,
        range_looks=range_looks,
        azimuth_looks=azimuth_looks,
    )
    
    # Validate
    is_valid, missing = validate_metadata_completeness(normalized)
    
    return (normalized, is_valid, missing)


# ============================================================================
# Physics Validation Functions
# ============================================================================

def validate_timing_consistency(
    metadata: Dict[str, Any],
    tolerance: float = 0.02,  # 2% relative tolerance
) -> Tuple[bool, List[str]]:
    """
    Validate timing parameter consistency (PRF vs azimuth_time_interval).
    
    For raw SLC data, PRF ≈ 1/azimuth_time_interval.
    For TOPS merged data, the relationship is more complex due to burst gaps.
    
    Args:
        metadata: Metadata dictionary containing timing parameters
        tolerance: Relative tolerance for validation (default 2%)
        
    Returns:
        Tuple of (is_consistent, list_of_warnings)
        
    Example:
        >>> meta = {"prf": 1685.8, "azimuth_time_interval": 0.000593}
        >>> is_ok, warnings = validate_timing_consistency(meta)
    """
    warnings = []
    
    prf = metadata.get("prf")
    ati = metadata.get("azimuth_time_interval")
    
    if prf is None or ati is None:
        # Can't validate without both parameters
        if prf is None:
            warnings.append("PRF not available for timing validation")
        if ati is None:
            warnings.append("azimuth_time_interval not available for timing validation")
        return (True, warnings)  # Not a failure, just can't validate
    
    try:
        prf = float(prf)
        ati = float(ati)
    except (TypeError, ValueError) as e:
        warnings.append(f"Cannot parse timing parameters: {e}")
        return (False, warnings)
    
    if prf <= 0:
        warnings.append(f"Invalid PRF value: {prf} Hz (must be positive)")
        return (False, warnings)
    
    if ati <= 0:
        warnings.append(f"Invalid azimuth_time_interval: {ati} s (must be positive)")
        return (False, warnings)
    
    # For raw SLC: expected_ati ≈ 1/prf
    expected_ati_raw = 1.0 / prf
    ratio = ati / expected_ati_raw
    
    # For IW TOPS merged data, ratio is typically ~2-4x due to burst gaps
    # This is expected and not an error
    if ratio < 0.5:
        # ATI is much smaller than expected - unusual
        warnings.append(
            f"Unusual timing ratio: azimuth_time_interval ({ati:.6f}s) < 1/PRF ({expected_ati_raw:.6f}s). "
            f"Ratio = {ratio:.3f}x"
        )
    elif 0.5 <= ratio <= 1.0 + tolerance:
        # Close to raw SLC timing - this is fine
        pass
    elif ratio > 4.0:
        # Ratio is very large - might indicate an issue
        warnings.append(
            f"Large timing ratio: azimuth_time_interval/PRF = {ratio:.2f}x. "
            f"This is unusual even for TOPS merged data (expected 1-4x)."
        )
    # Ratio between 1 and 4 is normal for TOPS merged data
    
    return (len(warnings) == 0, warnings)


def validate_doppler_centroid(
    metadata: Dict[str, Any],
) -> Tuple[bool, List[str]]:
    """
    Validate Doppler centroid metadata for Range-Doppler geocoding.
    
    Checks that Doppler polynomial coefficients are present and reasonable.
    
    Args:
        metadata: Metadata dictionary
        
    Returns:
        Tuple of (is_valid, list_of_issues)
    """
    issues = []
    
    # Check for Doppler centroid polynomial
    dc_poly = metadata.get("doppler_centroid_coefficients")
    dc_poly_alt = metadata.get("dc_estimate_poly")
    
    # Use dc_poly if it's not None (even if empty list), else try dc_poly_alt
    poly = dc_poly if dc_poly is not None else dc_poly_alt
    
    if poly is None:
        # Doppler centroid is optional for some processing modes
        # Return valid but note the absence
        issues.append("Doppler centroid polynomial not found (optional for some modes)")
        return (True, issues)
    
    if not isinstance(poly, (list, tuple)):
        issues.append(f"Doppler centroid polynomial has invalid type: {type(poly)}")
        return (False, issues)
    
    if len(poly) < 1:
        issues.append("Doppler centroid polynomial is empty")
        return (False, issues)
    
    # Check polynomial coefficients are numeric
    try:
        coeffs = [float(c) for c in poly]
    except (TypeError, ValueError) as e:
        issues.append(f"Doppler centroid polynomial contains non-numeric values: {e}")
        return (False, issues)
    
    # Check for reasonable DC value (first coefficient is DC at t=0)
    # Typical range for IW: -50 to +50 Hz
    dc_at_zero = coeffs[0]
    if abs(dc_at_zero) > 500:
        issues.append(
            f"Doppler centroid at t=0 is unusually large: {dc_at_zero:.1f} Hz "
            f"(typical range: ±50 Hz for zero-Doppler steering)"
        )
    
    # Validate polynomial degree (typically 3-5 coefficients)
    if len(coeffs) > 10:
        issues.append(
            f"Unusual Doppler centroid polynomial degree: {len(coeffs)} coefficients"
        )
    
    return (len(issues) == 0, issues)


def validate_geometry_metadata(
    metadata: Dict[str, Any],
) -> Tuple[bool, List[str]]:
    """
    Validate geometry metadata for geocoding.
    
    Checks pixel spacing, slant range time, and other geometric parameters.
    
    Args:
        metadata: Metadata dictionary
        
    Returns:
        Tuple of (is_valid, list_of_issues)
    """
    issues = []
    
    # Required geometry parameters
    required_params = [
        ("range_pixel_spacing", 1.0, 100.0, "meters"),
        ("azimuth_pixel_spacing", 1.0, 100.0, "meters"),
        ("slant_range_time", 0.001, 0.01, "seconds"),
    ]
    
    for param, min_val, max_val, units in required_params:
        val = metadata.get(param)
        if val is None:
            issues.append(f"Missing required geometry parameter: {param}")
            continue
        
        try:
            val = float(val)
        except (TypeError, ValueError):
            issues.append(f"Cannot parse {param}: {val}")
            continue
        
        if val <= 0:
            issues.append(f"Invalid {param}: {val} {units} (must be positive)")
        elif val < min_val or val > max_val:
            issues.append(
                f"Unusual {param}: {val} {units} (expected {min_val}-{max_val} {units})"
            )
    
    # Validate native spacing if present
    native_range = metadata.get("native_range_pixel_spacing")
    current_range = metadata.get("range_pixel_spacing")
    
    if native_range is not None and current_range is not None:
        try:
            native_range = float(native_range)
            current_range = float(current_range)
            
            if native_range > current_range * 1.1:  # Allow 10% tolerance
                issues.append(
                    f"Native range spacing ({native_range:.2f}m) > current spacing ({current_range:.2f}m). "
                    f"This is invalid - native should be <= current."
                )
        except (TypeError, ValueError):
            pass
    
    return (len(issues) == 0, issues)


def log_metadata_summary(
    metadata: Dict[str, Any],
    stage: str = "unknown",
    log_fn=None,
) -> None:
    """
    Log a summary of metadata state for debugging.
    
    Args:
        metadata: Metadata dictionary to summarize
        stage: Pipeline stage name for context
        log_fn: Optional logging function (defaults to logger.info)
    """
    if log_fn is None:
        log_fn = logger.info
    
    log_fn(f"=== Metadata Summary for {stage} ===")
    
    # Key categories
    categories = {
        "Timing": ["azimuth_time_interval", "prf", "first_line_time", "slant_range_time"],
        "Geometry": ["range_pixel_spacing", "azimuth_pixel_spacing", 
                     "native_range_pixel_spacing", "native_azimuth_pixel_spacing"],
        "Subswaths": ["subswaths", "burst_timings", "number_of_lines"],
        "Orbit": ["orbit_state_vectors", "orbit_ref_epoch_utc"],
        "Calibration": ["calibration_type", "calibration_applied"],
    }
    
    for category, keys in categories.items():
        present = [k for k in keys if k in metadata and metadata[k] is not None]
        missing = [k for k in keys if k not in metadata or metadata[k] is None]
        status = "✅" if len(missing) == 0 else "⚠️" if len(present) > 0 else "❌"
        log_fn(f"  {status} {category}: {len(present)}/{len(keys)} present")
        if missing:
            logger.debug(f"    Missing: {missing}")
    
    # Log key values for debugging
    key_values = [
        ("range_pixel_spacing", "m"),
        ("native_range_pixel_spacing", "m"),
        ("azimuth_time_interval", "s"),
        ("prf", "Hz"),
        ("slant_range_time", "s"),
    ]
    
    for key, unit in key_values:
        val = metadata.get(key)
        if val is not None:
            try:
                logger.debug(f"  • {key}: {float(val):.6f} {unit}")
            except (TypeError, ValueError):
                logger.debug(f"  • {key}: {val} (non-numeric)")
