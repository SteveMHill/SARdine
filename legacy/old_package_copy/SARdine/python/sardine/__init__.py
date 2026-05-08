"""
SARdine: High-Performance SAR Processing Library

A modern SAR data processing library for Sentinel-1 data, implemented in Rust 
with Python bindings. Provides complete processing pipeline from SLC to 
analysis-ready backscatter products.

IMPORTANT: Environment variables (SARDINE_SERDE_ONLY, SARDINE_REQUIRE_SUBSWATHS, 
SARDINE_ORBIT_CACHE, etc.) must be set BEFORE importing this module to take effect.
"""

__version__ = "0.2.1"
__author__ = "SARdine Contributors"

import os
import warnings
from pathlib import Path

# Track initial env var state to warn about post-import changes
_INITIAL_ENV_STATE = {
    'SARDINE_SERDE_ONLY': os.environ.get('SARDINE_SERDE_ONLY'),
    'SARDINE_REQUIRE_SUBSWATHS': os.environ.get('SARDINE_REQUIRE_SUBSWATHS'),
    'SARDINE_ORBIT_CACHE': os.environ.get('SARDINE_ORBIT_CACHE'),
    'SARDINE_STRICT': os.environ.get('SARDINE_STRICT'),
}
_MODULE_INITIALIZED = True


def check_env_var_timing():
    """
    Check if environment variables were changed after module import.
    
    Call this function to verify that critical env vars haven't been
    modified after sardine was imported (which would have no effect).
    
    Returns:
        List of warnings about env vars that changed post-import
    """
    changed = []
    for key, initial_value in _INITIAL_ENV_STATE.items():
        current_value = os.environ.get(key)
        if current_value != initial_value:
            changed.append(
                f"WARNING: {key} was changed after sardine import "
                f"('{initial_value}' → '{current_value}'). "
                f"This change has NO EFFECT. Set env vars before importing sardine."
            )
            warnings.warn(changed[-1], RuntimeWarning, stacklevel=2)
    return changed

# Import core module but NOT with star import - we'll control the public API
from . import _core

# =============================================================================
# PUBLIC API - These are the officially supported functions
# =============================================================================

# --- Classes ---
from ._core import SlcReader, BackgroundIoPool

# Try to import validation gateway (may not exist in all builds)
try:
    from ._core import ValidationGateway
except ImportError:
    ValidationGateway = None

# Try to import calibration job (may not exist in all builds)  
try:
    from ._core import CalibrationJob
except ImportError:
    CalibrationJob = None

# --- SLC Reading & Metadata ---
# Primary: create_cached_slc_reader() -> reader.get_cached_metadata()
from ._core import get_product_info

# --- Orbit ---
from ._core import apply_precise_orbit_file, load_orbit_file

# --- Deburst ---
from ._core import deburst_topsar_cached, read_slc_data_for_subswath_only

# --- Calibration ---
from ._core import (
    radiometric_calibration_with_denoising_cached,
    radiometric_calibration,  # Basic variant
    prepare_calibration_job_cached,
    run_calibration_job,
)

# --- Merge ---
from ._core import merge_subswaths_cached, topsar_merge_cached

# --- Multilook ---
from ._core import apply_multilooking, estimate_num_looks

# --- Terrain Flattening & Correction ---
from ._core import (
    terrain_correction,
    load_dem_for_bbox,
    get_dem_pixel_spacing,
)

# --- Speckle Filtering ---
from ._core import apply_speckle_filter

# --- Masking ---
from ._core import apply_masking

# --- dB Conversion ---
from ._core import convert_to_db_real

# --- STEP-2 Diagnostics (Jan 2026) ---
from ._core import (
    get_step2_diagnostics_config,
    enable_step2_diagnostics,
    disable_step2_diagnostics,
    read_step2_diagnostics_summary,
)

# --- Export ---
from ._core import (
    export_geotiff,
    export_cog_with_stac,
    generate_metadata,
    export_metadata_json,
)

# --- Quality Assessment ---
from ._core import perform_quality_assessment

# --- Data Download/Search ---
from ._core import download_sentinel1_products, search_sentinel1_products

# --- Utilities (internal, but commonly used) ---
from ._core import latlon_to_ecef

# --- Terrain Correction Footprint ---
from ._core import compute_sar_footprint_bbox

# =============================================================================
# INTERNAL APIs - Available but not part of public API contract
# Use at your own risk, may change without notice
# =============================================================================
# These are still importable via sardine._core.function_name() if needed:
# - extract_subswath_complex_data
# - extract_calibration_vectors  
# - radiometric_calibration_with_denoising (non-cached, use _cached version)
# - radiometric_calibration_direct_luts
# - prepare_calibration_job_cached, run_calibration_job
# - merge_subswaths (use merge_subswaths_cached)
# - db_to_linear_inplace_py, linear_to_db_inplace_py
# - export_db_parallel_py
# - dem_matching_refinement, apply_geocoding_offset
# - test_srtm_download, test_dem_reading
# - get_product_info_cached
# - export_metadata_xml
# - process_all_subswaths_batch

# =============================================================================
# PUBLIC API LIST - Controls 'from sardine import *'
# =============================================================================
__all__ = [
    # Version info
    "__version__",
    "__author__",
    
    # Classes
    "SlcReader",
    "BackgroundIoPool", 
    "ValidationGateway",
    "CalibrationJob",
    
    # Convenience functions
    "create_cached_slc_reader",
    "setup_default_orbit_cache",
    "check_env_var_timing",
    "topsar_merge",  # Python wrapper
    
    # SLC Reading & Metadata
    "get_product_info",
    
    # Orbit
    "apply_precise_orbit_file",
    "load_orbit_file",
    
    # Deburst  
    "deburst_topsar_cached",
    "read_slc_data_for_subswath_only",
    
    # Calibration
    "radiometric_calibration_with_denoising_cached",
    "radiometric_calibration",
    
    # Merge
    "merge_subswaths_cached",
    "topsar_merge_cached",
    
    # Multilook
    "apply_multilooking",
    "estimate_num_looks",
    
    # Terrain
    "terrain_correction",
    "load_dem_for_bbox",
    "get_dem_pixel_spacing",
    "compute_sar_footprint_bbox",
    
    # Speckle
    "apply_speckle_filter",
    
    # Masking
    "apply_masking",
    
    # dB Conversion
    "convert_to_db_real",
    
    # STEP-2 Diagnostics (Jan 2026)
    "get_step2_diagnostics_config",
    "enable_step2_diagnostics",
    "disable_step2_diagnostics",
    "read_step2_diagnostics_summary",
    
    # Export
    "export_geotiff",
    "export_cog_with_stac",
    "generate_metadata",
    "export_metadata_json",
    
    # Quality
    "perform_quality_assessment",
    
    # Download
    "download_sentinel1_products",
    "search_sentinel1_products",
    
    # Utilities
    "latlon_to_ecef",
]

def setup_default_orbit_cache(output_dir=None):
    """
    Set up default orbit cache for cached reader architecture.
    
    Creates a default orbit cache directory and sets the SARDINE_ORBIT_CACHE
    environment variable if not already configured. This ensures cached reader
    operations work seamlessly without manual environment setup.
    
    Args:
        output_dir: Base directory for orbit cache. If None, uses current working directory.
        
    Returns:
        str: Path to the orbit cache directory
        
    Scientific Justification:
        - Does not compromise orbit file accuracy (still requires real orbit data)
        - Maintains traceability of orbit data sources
        - Enables cached reader architecture without manual configuration
        - Follows graceful degradation principle without accuracy loss
    """
    if 'SARDINE_ORBIT_CACHE' not in os.environ:
        if output_dir is None:
            cache_dir = str(Path.cwd() / "orbit_cache")
        else:
            cache_dir = str(Path(output_dir) / "orbit_cache")
        
        # Create the cache directory
        Path(cache_dir).mkdir(parents=True, exist_ok=True)
        
        # Set environment variable for Rust code
        os.environ['SARDINE_ORBIT_CACHE'] = cache_dir
        
        return cache_dir
    else:
        return os.environ['SARDINE_ORBIT_CACHE']

# PERFORMANCE OPTIMIZATION: Convenience wrapper for cached SLC reader
def create_cached_slc_reader(slc_path: str):
    """
    Create a high-performance cached SLC reader for optimal metadata access.
    
    PERFORMANCE BENEFIT: ~93% faster metadata extraction compared to traditional method.
    
    Args:
        slc_path: Path to Sentinel-1 SLC product (ZIP or SAFE format)
        
    Returns:
        SlcReader: Cached reader instance for optimal performance
        
    Example:
        >>> reader = sardine.create_cached_slc_reader("S1A_IW_SLC_*.zip")
        >>> metadata = reader.get_cached_metadata()  # ~93% faster
        >>> print(f"Mission: {metadata['mission']}")
    """
    return SlcReader.new_with_full_cache(slc_path)

# Backward-compatible Python wrapper to preserve CLI API
def topsar_merge(input_path: str,
				 polarization: str,
				 calibration_type: str,
				 overlap_method: str,
				 output_grid: str = "auto",
				 annotation_metadata: dict | None = None):
	"""
	Backward-compatible wrapper matching the CLI call signature.
	It orchestrates per-subswath extraction + calibration, then calls the
	Rust core `topsar_merge(subswath_data, polarization, zip_path, annotation_metadata)`.

	Parameters match the existing CLI usage to avoid breaking changes.
	Currently `overlap_method` and `output_grid` are accepted but not used.
	"""
	import numpy as np

	pol = polarization.upper()
	subswaths = ["IW1", "IW2", "IW3"]

	# 1) Extract complex SLC and calibrate per subswath
	subswath_data = {}
	for sw in subswaths:
		try:
			ext = _core.extract_subswath_complex_data(str(input_path), sw, pol)
			if isinstance(ext, dict) and ext.get("status") != "success":
				continue
			slc_complex = ext["data"]  # numpy complex64 array

			cal = radiometric_calibration(
				str(input_path),
				sw,
				pol,
				calibration_type,
				slc_complex,
			)
			if isinstance(cal, dict) and cal.get("status") != "success":
				continue
			calib_arr = cal["data"]  # numpy float32 array (sigma0/beta0/gamma0)
			subswath_data[sw] = calib_arr
		except Exception:
			# Skip subswath on any failure to stay compatible with CLI behavior
			continue

	if len(subswath_data) < 2:
		raise ValueError(f"TOPSAR merge requires at least 2 subswaths, found {len(subswath_data)}")

	# 2) Call the Rust core merge with real metadata (zip_path)
	core_res = _core.topsar_merge(subswath_data, pol, str(input_path), annotation_metadata)

	# 3) Repackage to match CLI expectations (compat keys)
	#    - CLI expects 'intensity_data' instead of core 'data'
	#    - CLI expects nested 'metadata' with 'num_swaths' and 'overlap_count'
	out = {}
	out["intensity_data"] = core_res.get("data")
	out["valid_mask"] = core_res.get("valid_mask")
	out["quality_mask"] = core_res.get("quality_mask")

	meta = {
		"num_swaths": core_res.get("subswaths_processed"),
		"overlap_count": core_res.get("overlap_regions"),
		"valid_pixels": core_res.get("valid_pixels"),
		"processing_time": core_res.get("processing_time"),
		"algorithm": core_res.get("algorithm"),
		"overall_quality": core_res.get("overall_quality"),
		"radiometric_consistency": core_res.get("radiometric_consistency"),
	}
	out["metadata"] = meta

	# Pass-through some grid info if present
	for k in ("range_samples", "azimuth_samples", "range_pixel_spacing", "azimuth_pixel_spacing"):
		if k in core_res:
			out[k] = core_res[k]

	return out

# Feature availability flags
_TERRAIN_CORRECTION_AVAILABLE = True
_TERRAIN_FLATTENING_AVAILABLE = True
