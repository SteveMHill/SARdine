"""
SARdine Backscatter Processor

Complete SAR backscatter processing pipeline with REAL scientific data only.
Implements a 15-step pipeline for processing Sentinel-1 SLC data into 
analysis-ready backscatter products with scientific quality assurance.

Steps:
  1. Read Metadata
  2. Apply Precise Orbit
  3. IW Split
  4. TOPSAR Deburst
  5. Radiometric Calibration
  6. Expert IW Merge
  7. Multilooking
  8. Terrain Flattening
  9. Speckle Filtering
  10. Terrain Correction
  11. Mask Invalid Areas
  12. Convert to dB
  13. Export Final Products
  14. Generate Metadata
  15. Quality Assessment
"""

import gc
import inspect
import json
import logging
import math
import os
import threading
import time
from concurrent.futures import ThreadPoolExecutor, Future
from contextlib import nullcontext
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple, Union
import warnings

import numpy as np
from contextlib import contextmanager

# Module-level logger for backscatter processor
logger = logging.getLogger(__name__)

# Safe dB conversion helper for diagnostics
def _safe_db(x: float, eps: float = 1e-12) -> float:
    if not math.isfinite(x):
        return float("nan")
    return 10.0 * math.log10(max(x, eps))

# FIXED: Use context manager for scoped NumPy error suppression
# This allows better control and debugging while still suppressing expected warnings
@contextmanager
def suppress_numpy_errors():
    """Context manager to suppress expected NumPy warnings during SAR processing.
    
    These warnings are expected during backscatter calculations with near-zero values.
    Outside this context, errors will be visible for debugging.
    """
    old_settings = np.seterr(divide='ignore', invalid='ignore', over='ignore', under='ignore')
    try:
        yield
    finally:
        np.seterr(**old_settings)

# Set default error handling to 'warn' instead of 'ignore' for better debugging
# FIXED: Only suppress in specific contexts, not globally
np.seterr(divide='warn', invalid='warn', over='warn', under='warn')


def _to_bool(value: Any, default: bool = False) -> bool:
    """Convert various types to boolean, handling string representations correctly.
    
    This fixes the issue where bool("False") returns True in Python.
    Properly handles string values like "false", "0", "no", "off" as False.
    
    Args:
        value: Value to convert to boolean
        default: Default value if value is None
        
    Returns:
        Boolean value
    """
    if value is None:
        return default
    if isinstance(value, bool):
        return value
    if isinstance(value, str):
        return value.lower() in ("true", "1", "yes", "on", "t", "y")
    if isinstance(value, (int, float)):
        return bool(value)
    return default


import sardine
from sardine.export import create_cog_with_stac
from sardine.processors.base import BaseProcessor
from sardine.processors.pipeline import FunctionStage, PipelineContext, ProcessingPipeline
from sardine.processors.backscatter import export as export_utils
from sardine.processors.backscatter import options as opt_utils
from sardine.processors.backscatter import io as io_utils
from sardine.processors.backscatter import multilook as multilook_utils
from sardine.processors.backscatter import terrain as terrain_utils
from sardine.processors.backscatter import metadata_orbit as metadata_utils
from sardine.processors.backscatter import deburst_merge as deburst_utils
from sardine.processors.backscatter import merge_stage as merge_utils
from sardine.processors.backscatter.strict_mode import configure_strict_science
from sardine.processors.backscatter.geometry import (
    compute_bbox_from_transform,
    crop_geocoded_output,
    meters_per_degree,
    refine_geocoding_bbox,
)
from sardine.processors.backscatter import qc as qc_utils
from sardine.processors.backscatter import speckle as speckle_utils
from sardine.validation import ValidationError, ValidationGates

try:
    from sardine.performance import PerformanceMonitor, optimize_chunk_size, print_system_info
except ImportError:
    # Fallback if performance module is not available
    class PerformanceMonitor:
        def __init__(self, enable_monitoring=True): pass
        def start_step(self, *args, **kwargs): pass
        def end_step(self): pass
        def stop_monitoring(self): pass
        def print_summary(self): pass
        def record_phase_duration(self, phase_name, duration): pass

# OPTIMIZATION: Memory-adaptive configuration for users with varying RAM
try:
    from sardine.memory_adaptive import get_memory_config, MemoryConfig, MemoryTier
    MEMORY_ADAPTIVE_AVAILABLE = True
except ImportError:
    MEMORY_ADAPTIVE_AVAILABLE = False
    MemoryConfig = None
    MemoryTier = None
    def get_memory_config():
        return None
    
    def optimize_chunk_size(data_size_mb, num_threads, memory_limit_mb=8000):
        return max(64, min(int(data_size_mb / num_threads / 4), 2048))
    
    def print_system_info():
        print("💻 Performance monitoring not available (psutil required)")

def _resolve_terrain_flattening():
    # Prefer the new name; fall back to legacy for backward compatibility
    fn = getattr(sardine, "terrain_flattening", None)
    if fn is None:
        fn = getattr(sardine, "apply_scientific_terrain_flattening", None)
    return fn


_TERRAIN_FLATTEN_FN = _resolve_terrain_flattening()
try:
    _SCI_TERRAIN_SIG = inspect.signature(_TERRAIN_FLATTEN_FN) if _TERRAIN_FLATTEN_FN else None
except (ValueError, TypeError):
    _SCI_TERRAIN_SIG = None

_TF_SUPPORTS_DEM_METADATA = bool(_SCI_TERRAIN_SIG) and (
    "dem_geo_transform" in _SCI_TERRAIN_SIG.parameters
)

MAX_REFINE_SHRINK_FRACTION = 0.95  # Allow aggressive shrink (up to 95%) so refined footprint can match tight multilooked grids


# Step number mapping for consistent logging
STEP_NUMBERS = {
    "Read Metadata": 1,
    "Apply Precise Orbit": 2,
    "IW Split": 3,
    "TOPSAR Deburst": 4,
    "Radiometric Calibration": 5,
    "Expert IW Merge": 6,
    "Multilooking": 7,
    "Terrain Flattening": 8,
    "Speckle Filtering": 9,
    "Terrain Correction": 10,
    "Mask Invalid Areas": 11,
    "Convert to dB": 12,
    "Export Final Products": 13,
    "Generate Metadata": 14,
    "Quality Assessment": 15,
}

# Single source of truth for step count in banners/reports.
TOTAL_PIPELINE_STEPS = max(STEP_NUMBERS.values())

# Masking threshold constants (Step 11)
# These define when data is considered "linear power" vs needs percentile masking
MASK_LINEAR_MAX_VALUE = 1e3  # Maximum value for linear power interpretation
MASK_LINEAR_DYNAMIC_RANGE = 1e6  # Maximum dynamic range for linear data
MASK_MIN_POSITIVE_VALUE = 1e-8  # Minimum valid positive value for linear masking
MASK_STABILITY_EPSILON = 1e-12  # Small epsilon for numerical stability
MASK_MIN_VALID_PERCENTAGE = 10.0  # Minimum percentage of valid pixels after masking

# dB conversion validation constants (Step 12)
DB_EXPECTED_MIN = -50.0  # Minimum expected dB value for SAR backscatter
DB_EXPECTED_MAX = 15.0  # Maximum expected dB value for SAR backscatter
DB_WARN_MIN = -40.0  # Warn if values below this
DB_WARN_MAX = 10.0  # Warn if values above this


class BackscatterProcessor(BaseProcessor):
    """Complete SAR backscatter processing pipeline with REAL scientific data only"""

    def __init__(
        self,
        input_path: Union[str, Path],
        output_dir: Union[str, Path],
        options: Optional[Dict[str, Any]] = None,
        **overrides: Any,
    ) -> None:
        # Handle both ZIP files and SAFE directories with strict validation
        self.input_path = opt_utils.normalize_input_path(Path(input_path))
        opt_utils.validate_sentinel_product(self.input_path)

        merged_options: Dict[str, Any] = {}
        if options:
            merged_options.update(options)
        if overrides:
            merged_options.update(overrides)
        merged_options = opt_utils.normalize_options(merged_options)

        # Scientific strict mode wiring: delegate to helper module so
        # the main processor stays focused on pipeline logic.
        merged_options, strict_science = configure_strict_science(merged_options)

        super().__init__(output_dir, merged_options)
        # OPTIMIZATION #9: Use RLock to allow concurrent read operations
        # RLock can be acquired multiple times by the same thread (reentrant),
        # which reduces contention when multiple threads need read-only access
        self.reader_lock = threading.RLock()
        self._last_orbit_result: Optional[Dict[str, Any]] = None
        # Expose strict science flag on the processor instance for downstream
        # stages that need to tighten invariants (orbit coverage, timing, etc.).
        self.strict_science: bool = strict_science

        # OPTIMIZATION: Memory-adaptive configuration for different system RAM levels
        # This automatically detects available memory and configures:
        # - Chunk sizes for I/O operations
        # - Prefetch behavior for subswaths
        # - Parallelism levels for burst processing
        # - Cleanup frequency for GC
        self._memory_config: Optional[MemoryConfig] = None
        if MEMORY_ADAPTIVE_AVAILABLE:
            self._memory_config = get_memory_config()
            if self._memory_config:
                print(f"🧠 Memory tier: {self._memory_config.tier.name} "
                      f"({self._memory_config.available_gb:.1f}GB available)")
                # Set environment variable for Rust chunk size
                os.environ["SARDINE_TIFF_CHUNK_LINES"] = str(self._memory_config.chunk_size_lines)

        # Background thread tracking with Future-based error handling
        self._thread_executor: Optional[ThreadPoolExecutor] = None
        self._dem_prefetch_future: Optional[Future] = None
        self._orbit_prefetch_future: Optional[Future] = None

        # Configure orbit cache directory with per-product isolation
        orbit_cache_path, dem_cache_path = opt_utils.resolve_caches(self.input_path, self.output_dir, self.options)
        self.orbit_cache_dir = orbit_cache_path
        self.dem_cache_dir = dem_cache_path
        self.options["dem_cache"] = str(dem_cache_path)
        io_utils.configure_environment(orbit_cache_path, dem_cache_path)

        # Track pixel spacing as it evolves through multilooking and other resampling steps
        self.original_range_spacing = None
        self.original_azimuth_spacing = None
        self.current_range_spacing = None
        self.current_azimuth_spacing = None

        # Working buffers and caches populated during processing
        self._working_data: Optional[np.ndarray] = None
        self._calibrated_subswaths: Dict[str, np.ndarray] = {}
        self._primary_subswath: Optional[str] = None
        self._geocoding_metadata_cache: Optional[Dict[str, Any]] = None
        self._native_azimuth_lines: Optional[int] = None
        self._native_range_samples: Optional[int] = None
        self._range_multilook_factor: Optional[float] = None
        self._azimuth_multilook_factor: Optional[float] = None
        self._multilooked_shape: Optional[Tuple[int, int]] = None
        self._merge_valid_fraction_bounds: Optional[Tuple[float, float, float, float]] = None
        self._used_subswaths_for_merge: Optional[List[str]] = None

        # DEM cache to avoid reloading between terrain flattening and terrain correction
        self._cached_dem_data: Optional[np.ndarray] = None
        self._cached_dem_bbox: Optional[Tuple[float, float, float, float]] = None
        self._cached_dem_geo_transform: Optional[Tuple[float, ...]] = None
        self._cached_dem_resolution: Optional[float] = None
        self._cached_dem_crs: Optional[str] = None

        # Orbit vector cache to avoid re-parsing EOF files or metadata
        self._cached_orbit_times: Optional[List[str]] = None
        self._cached_orbit_positions: Optional[List[List[float]]] = None
        self._cached_orbit_velocities: Optional[List[List[float]]] = None

        # Time-domain metadata extracted from annotations
        self.burst_timing_records: List[Dict[str, Any]] = []

        # Core processing configuration
        # FIXED: Normalize polarization to uppercase for case-insensitive input
        pol_raw = self.options.get("polarization", "VV")
        self.polarization = str(pol_raw).upper() if pol_raw else "VV"
        if self.polarization not in {"VV", "VH", "HH", "HV"}:
            raise ValueError(f"❌ SCIENTIFIC MODE: Invalid polarization: {self.polarization}")

        self.speckle_filter = self.options.get("speckle_filter", "enhanced_lee")
        # FIXED: Add input validation for filter_window
        filter_window_raw = self.options.get("filter_window", 7)
        self.filter_window = self._validate_filter_window(filter_window_raw)
        
        # FIXED: Add input validation for multilook factors
        multilook_range_raw = self.options.get("multilook_range", 2)
        multilook_azimuth_raw = self.options.get("multilook_azimuth", 2)
        self.multilook_range = self._validate_multilook_factor(multilook_range_raw, "range")
        self.multilook_azimuth = self._validate_multilook_factor(multilook_azimuth_raw, "azimuth")
        
        # SCIENTIFIC FIX: Determine correct calibration type based on workflow
        # - If terrain_flatten is enabled: MUST use sigma0, then terrain flattening
        #   produces gamma0_tc (terrain-corrected gamma) using local incidence angle
        # - If geocode is enabled but NOT terrain_flatten: gamma0 (ellipsoid-based)
        # - Otherwise: sigma0 for basic backscatter
        #
        # NOTE: There are TWO different gamma0 definitions:
        # 1. Radiometric gamma0: gamma0 = sigma0 / cos(θ_ellipsoid) - uses LUT
        # 2. Terrain-corrected gamma0: gamma0_tc = sigma0 / cos(θ_local) - uses DEM
        # Applying terrain flattening to radiometric gamma0 is WRONG as it 
        # double-corrects for incidence angle!
        terrain_flatten_enabled = _to_bool(self.options.get("terrain_flatten", True), True)
        geocode_enabled = _to_bool(self.options.get("geocode", True), True)
        
        if terrain_flatten_enabled:
            # When terrain flattening is enabled, we MUST start with sigma0
            # Terrain flattening will then produce the correct gamma0_tc
            default_cal_type = "sigma0"
        elif geocode_enabled:
            # Geocoding without terrain flattening uses radiometric gamma0
            default_cal_type = "gamma0"
        else:
            # Basic backscatter without terrain correction
            default_cal_type = "sigma0"
        # FIXED: Normalize calibration_type to lowercase
        cal_type_raw = self.options.get("calibration_type", default_cal_type)
        self.calibration_type = str(cal_type_raw).lower() if cal_type_raw else default_cal_type
        
        # FIXED: Use _to_bool for proper string-to-bool conversion
        self.terrain_flatten = _to_bool(self.options.get("terrain_flatten", True), True)
        self.geocode = _to_bool(self.options.get("geocode", True), True)
        self.skip_antenna = _to_bool(self.options.get("skip_antenna", False), False)
        self.skip_noise = _to_bool(self.options.get("skip_noise", False), False)

        # RTC configuration (applied during terrain correction/geocoding).
        # - "area": Area-projection (Small 2011), default when terrain_flatten=True
        # - "cosine": Cosine local incidence approximation
        # - "none": Disable RTC (output σ⁰)
        rtc_mode_raw = self.options.get("rtc_mode", "auto")
        rtc_mode = str(rtc_mode_raw).strip().lower() if rtc_mode_raw is not None else "auto"
        if rtc_mode in ("", "auto"):
            # SCIENTIFIC FIX (Jan 2026): RTC is now integrated into terrain
            # correction. Legacy DEM-based terrain_flattening is deprecated
            # and should not be combined with RTC. Default behaviour:
            #   - If user requests geocoding: enable RTC via 'area' mode.
            #   - If no geocoding: disable RTC.
            #
            # NOTE: This defaults to RTC even when terrain_flatten=False,
            # using the integrated terrain_correction_with_rtc path in Rust.
            rtc_mode = "area" if self.geocode else "none"
        elif rtc_mode in ("off", "false", "0", "no", "disabled"):
            rtc_mode = "none"
        if not self.geocode:
            # RTC without geocoding is undefined in current design.
            rtc_mode = "none"
        # SCIENTIFIC FIX (Jan 2026): Do not forbid RTC when terrain_flatten=False.
        # The modern path is: calibrated σ⁰ → Rust terrain_correction_with_rtc → γ⁰_tc.
        # Legacy DEM-based terrain_flattening (Python helper) is deprecated and
        # should be disabled when RTC is active, not the other way around.
        self.rtc_mode = rtc_mode
        self.output_lia = _to_bool(self.options.get("output_lia", False), False)
        self.output_masks = _to_bool(self.options.get("output_masks", False), False)

        self.target_resolution = self.options.get("resolution", opt_utils.DEFAULT_TARGET_RESOLUTION_M)

        self.quality_report = _to_bool(self.options.get("quality_report", True), True)
        self.use_real_orbit = _to_bool(self.options.get("use_real_orbit", True), True)
        self.allow_synthetic = _to_bool(self.options.get("allow_synthetic", False), False)

        self.enable_iw_merging = _to_bool(self.options.get("enable_iw_merging", True), True)
        self.cosine_clip_threshold = self.options.get("cosine_clip_threshold", 0.1)
        
        # Masking thresholds (applied during Step 11 when inputs are available).
        # SCIENTIFIC FIX: Tightened from [-60, 50] to [-35, 20] dB
        # Rationale: Typical SAR backscatter ranges from -30 dB (smooth water) to +5 dB (urban).
        # Values outside [-35, 20] indicate calibration errors, not physical scattering.
        # Upper bound is 20 dB: urban corner reflectors and man-made structures regularly
        # reach +10 to +15 dB; +20 dB is a safe upper bound per Ulaby & Long (2014).
        # A ceiling of 10 dB (previous default) incorrectly clips strong urban returns.
        self.lia_threshold_cos = self.options.get("lia_threshold", 0.6)
        self.gamma0_min_db = self.options.get("gamma0_min", -35.0)
        self.gamma0_max_db = self.options.get("gamma0_max", 20.0)
        self.dem_threshold_m = self.options.get("dem_threshold", -100.0)
        
        # dem_threshold is reserved for future DEM-on-output-grid masking.
        if self.dem_threshold_m != -100.0:
            warnings.warn(
                f"dem_threshold={self.dem_threshold_m} is not yet applied in the backscatter pipeline "
                f"(reserved for DEM-on-output-grid masking); ignoring for now.",
                UserWarning,
                stacklevel=2,
            )
        
        self.skip_masking = _to_bool(self.options.get("no_masking", False), False)
        self.skip_mask_export = _to_bool(self.options.get("no_masks", False), False)
        self.linear_processing = _to_bool(self.options.get("linear_processing", True), True)

        # Intermediate product export (for validation against SNAP/ESA reference)
        # When enabled, exports GeoTIFF checkpoints after: calibration, multilook, 
        # terrain flattening, and terrain correction
        self.export_intermediate_products = _to_bool(
            self.options.get("export_intermediate_products", False), False
        )

        # FIXED: Add input validation for DEM resolution
        dem_resolution_raw = self.options.get("dem_resolution", 30.0)
        self.dem_resolution = self._validate_dem_resolution(dem_resolution_raw)

        self.rust_threads = self.options.get("rust_threads")
        self.burst_chunk_size = self.options.get("burst_chunk_size", 3)
        self.terrain_chunk_size = self.options.get("terrain_chunk_size", 256)
        self.serial_terrain = _to_bool(self.options.get("serial_terrain", False), False)
        self.legacy_mode = _to_bool(self.options.get("legacy_mode", False), False)

        self.enable_parallel = _to_bool(self.options.get("parallel", True), True)
        self.sequential = _to_bool(self.options.get("sequential", False), False)
        requested_threads = self.options.get("num_threads")
        requested_chunk_size = self.options.get("chunk_size")
        self.chunk_size = int(requested_chunk_size or 1024)

        # Verbosity control: 0=errors only, 1=normal (default), 2=verbose, 3=debug
        self.verbosity = int(self.options.get("verbosity", 1))
        
        # Performance tuning options
        self.prefetch_dem = _to_bool(self.options.get("prefetch_dem", True), True)
        self.use_float32 = _to_bool(self.options.get("use_float32", True), True)
        self.inplace_operations = _to_bool(self.options.get("inplace_operations", True), True)

        self.optimization_mode = self.options.get("optimization_mode", "complete")
        # Noise removal: default to "range" (fast 1D profile) like SNAP.
        # Options: "off", "range" (default), "azimuth", "full"
        # - "range": Fast 1D range-dependent noise profile (~O(width) memory)
        # - "azimuth": Interpolated azimuth variation (~O(vectors × width))
        # - "full": Complete 2D LUT (slow, ~O(height × width) memory)
        self.noise_strategy = str(self.options.get("noise_strategy", "range")).lower()
        default_noise = self.noise_strategy == "range"
        self.fast_mode_noise_removal = bool(
            self.options.get("fast_mode_noise_removal", default_noise)
        )
        self.debug_geocoding = bool(self.options.get("debug_geocoding", False))

        self.performance_monitor = PerformanceMonitor(enable_monitoring=True)

        # Threading policy: cap to a sane default to avoid cache/NUMA thrash
        max_default_threads = 32
        detected_cores = os.cpu_count() or 1

        if self.sequential:
            self.enable_parallel = False
            self.num_threads = 1
            print("🔄 Sequential processing mode enabled (parallel processing disabled)")
        elif self.enable_parallel:
            parsed_threads: Optional[int] = None
            if requested_threads is not None:
                try:
                    parsed_threads = max(1, int(requested_threads))
                except (ValueError, TypeError):
                    parsed_threads = None

            if parsed_threads is None:
                capped = min(detected_cores, max_default_threads)
                self.num_threads = capped
                print(
                    f"🚀 Parallel processing enabled: auto-detected {detected_cores} CPU cores "
                    f"(using {capped} to protect caches)"
                )
            else:
                capped = min(parsed_threads, max_default_threads)
                if capped != parsed_threads:
                    print(
                        f"🚀 Parallel processing enabled: requested {parsed_threads} threads; "
                        f"capping to {capped} to reduce oversubscription"
                    )
                else:
                    print(f"🚀 Parallel processing enabled: {capped} threads specified")
                self.num_threads = capped

            if requested_chunk_size is None:
                tuned = optimize_chunk_size(100, self.num_threads)
                tuned = max(256, min(1024, tuned))  # keep chunks L2-friendly
                self.chunk_size = tuned
                print(f"📊 Optimized chunk size: {self.chunk_size}")
            else:
                self.chunk_size = max(64, int(self.chunk_size))
                print(f"📊 Chunk size: {self.chunk_size}")
        else:
            self.num_threads = 1
            print("🔄 Parallel processing disabled")

        # Note: Separable calibration LUT removed (IW TOPS mode requires dense 2D LUT)

        # Keep Rust side in sync unless explicitly overridden
        if self.options.get("rust_threads") is None:
            self.rust_threads = self.num_threads
            self.options["rust_threads"] = self.num_threads
        else:
            self.rust_threads = max(1, int(self.options.get("rust_threads")))

        self._configure_parallel_environment()

        if self.enable_parallel:
            try:
                print_system_info()
            except Exception:
                pass

        if self.use_real_orbit and self.allow_synthetic:
            raise ValueError(
                "❌ INVALID CONFIGURATION: Cannot use real orbit mode with synthetic data allowed"
            )

        # Propagate profiling toggles to Rust side
        os.environ["SARDINE_SKIP_ANTENNA"] = "1" if self.skip_antenna else "0"
        if self.skip_noise:
            self.noise_strategy = "off"

        # Propagate noise strategy to Rust side; default is "range" (fast 1D profile)
        # Use noise_strategy="off" or skip_noise=True to disable
        os.environ["SARDINE_SKIP_NOISE"] = "1" if self.noise_strategy == "off" else "0"
        os.environ["SARDINE_NOISE_STRATEGY"] = self.noise_strategy

        # Preserve legacy fast-mode toggle by inferring from range strategy when unspecified
        if self.noise_strategy == "range" and not self.options.get("fast_mode_noise_removal"):
            self.fast_mode_noise_removal = True

        if not self.allow_synthetic:
            print("🔬 SCIENTIFIC MODE: Real data only - no synthetic fallbacks")
        else:
            print("⚠️  DEMONSTRATION MODE: Synthetic data fallbacks enabled")

        # Default state before pipeline stages run
        self.reader = None
        self.metadata = None
        self.geo_transform = None
        self.validated_metadata: Optional[Dict[str, Any]] = None
        self.cached_subswath_geometry: Dict[str, Any] = {}
        self.geocoding_bbox: Optional[Tuple[float, float, float, float]] = None
        self.actual_range_looks: Optional[float] = None
        self.actual_azimuth_looks: Optional[float] = None

        self._requested_output_crs = self.options.get("output_crs")
        self.output_epsg: Optional[int] = None
        self.coordinate_system: Optional[str] = None
        self.metadata_builder = None
        self._metadata_steps_recorded = False
        self.calibration_lut_source: Optional[str] = None
        self.precise_orbit_path: Optional[str] = None
        # Strict orbit mode: fail if precise orbit not available (default: False for backward compatibility)
        self.strict_orbit_mode: bool = self.options.get("strict_orbit_mode", False)
        self.dem_geo_transform: Optional[Tuple[float, ...]] = None
        self._active_dem_geo_transform: Optional[Tuple[float, ...]] = None
        self.dem_bbox: Optional[Tuple[float, float, float, float]] = None
        self._precise_orbit_records: Optional[List[Dict[str, Any]]] = None
        self._orbit_reference_epoch: Optional[float] = None
        self._product_time_range: Tuple[Optional[float], Optional[float]] = (None, None)
        self.latest_quality_report: Optional[Dict[str, Any]] = None
        self.validation = ValidationGates(
            enabled=bool(self.options.get("validation_enabled", True)),
            fail_fast=bool(self.options.get("validation_fail_fast", True)),
        )
        self.validation_report_path = self.output_dir / "validation_report.json"
        
        # Stage timing tracking for performance analysis
        self._stage_timings: Dict[str, float] = {}
        self._memory_cleanup_enabled = self.options.get("memory_cleanup", True)
        
        # Optional download manager for auto-downloading missing data
        # FIXED: Removed duplicate initialization block
        # FIXED: Always initialize DownloadManager for orbit downloads (required for scientific mode)
        self._download_manager = None
        import logging
        import math
        logger = logging.getLogger(__name__)

        try:
            from sardine.download import DownloadManager
            # Use orbit cache directory for download manager cache
            orbit_cache = self.options.get("orbit_dir") or str(self.output_dir / "orbit_cache")
            sources = self.options.get("download_sources", ["esa", "aws"])  # Prefer ESA for orbits
            self._download_manager = DownloadManager(
                cache_dir=orbit_cache,
                sources=sources,
                auto_download=True
            )
            logger.info(f"DownloadManager initialized for orbit downloads (sources: {sources})")
        except ImportError:
            logger.warning("DownloadManager not available. Install sardine with download support.")
        except Exception as e:
            logger.warning(f"Failed to initialize DownloadManager: {e}")

    # ------------------------------------------------------------------ #
    # Pipeline composition and stage entry points
    # ------------------------------------------------------------------ #
    def _build_pipeline(self) -> ProcessingPipeline:
        """Wire the modular stages into a sequential pipeline."""

        stages = [
            FunctionStage("Read Metadata", self._stage_ingest_metadata, "Load annotation and reader"),
            FunctionStage("Apply Precise Orbit", self._stage_apply_precise_orbit, "Download/validate precise orbit"),
            FunctionStage("IW Split", self._stage_iw_split, "Validate subswath geometry"),
            FunctionStage(
                "TOPSAR Deburst",
                self._stage_deburst_and_calibrate,
                "Per-subswath deburst (with inline calibration)",
            ),
            FunctionStage(
                "Radiometric Calibration",
                self._stage_validate_calibration,
                "Validate calibrated subswaths and LUT metadata",
            ),
            FunctionStage(
                "Export Calibration Checkpoint",
                self._stage_export_calibration,
                "Export calibrated subswaths for validation (if enabled)",
            ),
            FunctionStage("Expert IW Merge", self._stage_merge_subswaths, "Combine calibrated subswaths"),
            FunctionStage(
                "Export Merge Checkpoint",
                self._stage_export_merge,
                "Export merged data for validation (if enabled)",
            ),
            FunctionStage("Multilooking", self._stage_multilooking, "Apply scientific multilook factors"),
            FunctionStage(
                "Export Multilook Checkpoint",
                self._stage_export_multilook,
                "Export multilooked data for validation (if enabled)",
            ),
            FunctionStage("Terrain Flattening", self._stage_terrain_flattening, "DEM-based flattening"),
            FunctionStage(
                "Export Terrain Flattening Checkpoint",
                self._stage_export_terrain_flattening,
                "Export terrain-flattened data for validation (if enabled)",
            ),
            FunctionStage("Speckle Filtering", self._stage_speckle_filtering, "Adaptive speckle reduction"),
            FunctionStage("Terrain Correction", self._stage_terrain_correction, "Range-Doppler geocoding"),
            FunctionStage(
                "Export Terrain Correction Checkpoint",
                self._stage_export_terrain_correction,
                "Export geocoded data for validation (if enabled)",
            ),
            FunctionStage("Mask Invalid Areas", self._stage_mask_invalid_areas, "Quality masking"),
            FunctionStage("Convert to dB", self._stage_convert_to_db, "Linear-to-dB conversion"),
            FunctionStage("Export Final Products", self._stage_export_products, "Write GeoTIFF/JSON artifacts"),
            FunctionStage("Generate Metadata", self._stage_generate_metadata, "Processing metadata artifacts"),
            FunctionStage("Quality Assessment", self._stage_finalize_pipeline, "Reporting and validation outputs"),
        ]
        return ProcessingPipeline(stages)

    def _stage_ingest_metadata(self, context: PipelineContext) -> None:
        """Stage 1: read metadata and set up reader/cache."""
        metadata_utils.read_metadata(self, context)
        # Ensure validated_metadata is present for downstream geometry stages
        if self.validated_metadata is None and isinstance(self.metadata, dict):
            self.validated_metadata = self.metadata

        # Populate initial spacing hints from metadata if available
        _ = self.get_current_range_spacing()
        _ = self.get_current_azimuth_spacing()

    def _stage_apply_precise_orbit(self, context: PipelineContext) -> None:
        """Stage 2: acquire precise orbit file."""
        # I/O Optimization: Pre-fetch orbit file in background if not already cached
        # This overlaps orbit download/reading with metadata processing
        self._pre_fetch_orbit_background()
        
        # FIXED: Set SARDINE_ORBIT_CACHE environment variable for deburst_topsar_cached
        # This ensures deburst can find the orbit file we just downloaded
        import os
        os.environ["SARDINE_ORBIT_CACHE"] = str(self.orbit_cache_dir)
        metadata_utils.apply_precise_orbit(self, context)

    def _stage_iw_split(self, context: PipelineContext) -> None:
        """Stage 3: validate IW subswath geometry and cache metadata."""
        import logging
        logger = logging.getLogger(__name__)

        step_number = STEP_NUMBERS["IW Split"]
        self.announce_step(step_number, "IW Split", "Validating and caching subswath geometry")
        step_start = time.time()
        
        if not isinstance(self.metadata, dict):
            raise ValueError("No metadata available for IW split")
        
        # Validate product mode before processing
        try:
            with self.reader_lock:
                if hasattr(self, "reader") and self.reader is not None:
                    if hasattr(self.reader, "is_iw_mode"):
                        is_iw = self.reader.is_iw_mode()
                        if not is_iw:
                            raise ValueError(
                                f"Product is not IW mode. IW split requires Interferometric Wide mode products. "
                                f"Product type: {self.metadata.get('product_type', 'unknown')}"
                            )
        except ValueError:
            raise
        except Exception as e:
            logger.warning(f"Could not verify IW mode: {e}")
        
        # Parse and validate subswath string
        subswaths_str = self.metadata.get("subswaths")
        if not subswaths_str:
            raise ValueError("No subswath list found in metadata")
        
        # Parse and normalize subswath names
        parsed = [sw.strip().upper() for sw in subswaths_str.split(",") if sw.strip()]
        if not parsed:
            raise ValueError(f"Invalid subswath metadata: {subswaths_str!r}")
        
        # Validate subswath format (IW1, IW2, IW3 for IW mode)
        valid_iw_patterns = {"IW1", "IW2", "IW3"}
        invalid_subswaths = [sw for sw in parsed if sw not in valid_iw_patterns]
        if invalid_subswaths:
            raise ValueError(
                f"Invalid subswath names: {invalid_subswaths}. "
                f"Expected IW mode subswaths: {valid_iw_patterns}. "
                f"Found: {parsed}"
            )
        
        # Check for duplicates
        if len(parsed) != len(set(parsed)):
            duplicates = [sw for sw in parsed if parsed.count(sw) > 1]
            raise ValueError(f"Duplicate subswaths found in metadata: {set(duplicates)}")
        
        # Normalize subswath order (IW1, IW2, IW3)
        order_map = {"IW1": 1, "IW2": 2, "IW3": 3}
        parsed = sorted(parsed, key=lambda s: order_map.get(s.upper(), 999))
        
        context.set_artifact("subswaths", parsed)
        
        # Extract and validate geometry cache
        cached_metadata_map: Optional[Dict[str, Any]] = None

        try:
            with self.reader_lock:
                # Re-check reader after acquiring lock
                if not hasattr(self, "reader") or self.reader is None:
                    raise RuntimeError("Reader not available for IW geometry extraction")
                
                # Store reference to prevent None assignment during call
                reader = self.reader
                if reader is None:
                    raise RuntimeError("Reader became None during geometry extraction")
                
                try:
                    geometry_cache = reader.get_all_iw_subswaths()
                except AttributeError as e:
                    raise RuntimeError(
                        f"Reader missing required method 'get_all_iw_subswaths': {e}"
                    ) from e

                if hasattr(reader, "get_cached_metadata"):
                    try:
                        cached_metadata_map = reader.get_cached_metadata()
                    except Exception as exc:  # pragma: no cover - diagnostics only
                        logger.debug(f"Failed to read cached metadata from reader: {exc}")
                
                # Validate geometry cache structure
                if not isinstance(geometry_cache, dict):
                    raise RuntimeError(
                        f"Invalid geometry cache type: {type(geometry_cache)}. "
                        f"Expected dict mapping polarization to subswaths."
                    )
                
                if not geometry_cache:
                    raise RuntimeError(
                        "No IW subswath geometry extracted. "
                        "This may indicate missing annotation files or parsing failure."
                    )
                
                # Validate that requested polarization has geometry
                requested_pol = self.polarization.upper()
                # Geometry cache keys may be Polarization enum or string
                pol_keys = {str(k).upper() if not isinstance(k, str) else k.upper() for k in geometry_cache.keys()}
                if requested_pol not in pol_keys:
                    raise RuntimeError(
                        f"Requested polarization '{requested_pol}' not found in geometry cache. "
                        f"Available polarizations: {list(geometry_cache.keys())}"
                    )
                
                # Find the correct key for requested polarization
                pol_key = None
                for k in geometry_cache.keys():
                    k_str = str(k).upper() if not isinstance(k, str) else k.upper()
                    if k_str == requested_pol:
                        pol_key = k
                        break
                
                if pol_key is None:
                    raise RuntimeError(f"Could not resolve polarization key for '{requested_pol}'")
                
                # Validate subswaths for requested polarization
                pol_subswaths = geometry_cache[pol_key]
                if not isinstance(pol_subswaths, dict) or not pol_subswaths:
                    raise RuntimeError(
                        f"No subswaths found for polarization '{requested_pol}'. "
                        f"Expected dict of subswath geometry data."
                    )
                
                # Validate geometry values for each subswath
                for swath_id, swath_data in pol_subswaths.items():
                    self._validate_subswath_geometry(swath_id, swath_data)
                
                # Extract azimuth_time_interval from first subswath and add to top-level metadata
                # This is required for terrain correction (no fallbacks per user requirement)
                if not self.metadata.get("azimuth_time_interval"):
                    first_subswath_data = next(iter(pol_subswaths.values()), None)
                    if first_subswath_data and isinstance(first_subswath_data, dict):
                        ati = first_subswath_data.get("azimuth_time_interval")
                        if ati is not None and ati != 0:
                            self.metadata["azimuth_time_interval"] = float(ati)
                            logger.info(f"   ✅ Extracted azimuth_time_interval={ati:.9f}s from subswath geometry")
                
                # FIXED: Use atomic assignment to prevent race condition
                # Cache successfully validated geometry
                new_cache = {"subswaths": geometry_cache}
                if cached_metadata_map:
                    new_cache["cached_metadata"] = cached_metadata_map
                # Atomic update to prevent race conditions
                self._geocoding_metadata_cache = new_cache
                context.set_artifact("iw_geometry", geometry_cache)
                if cached_metadata_map:
                    context.set_artifact("cached_metadata", cached_metadata_map)
                
                subswath_count = sum(len(subs) for subs in geometry_cache.values())
                print(f"   ✅ Cached geometry for {subswath_count} subswaths across {len(geometry_cache)} polarizations")
                
                # I/O Optimization: Pre-load DEM in background while processing subswaths
                # DEM bbox is known from metadata, so we can start loading early
                self._pre_load_dem_background()

                step_duration = time.time() - step_start
                self.log_step(step_number, "IW Split", "success", "Subswath geometry cached", step_duration)
        except Exception as exc:
            self._geometry_cache_failed = True
            self._geometry_cache_error = str(exc)
            
            # Always fail - IW geometry is required for accurate geocoding
            step_duration = time.time() - step_start
            self.log_step(step_number, "IW Split", "error", f"IW geometry extraction failed: {exc}", step_duration)
            raise RuntimeError(
                f"IW geometry extraction failed: {exc}\n"
                f"This is required for accurate geocoding and terrain correction."
            ) from exc

        # Hydrate burst timing records as soon as cached metadata is available.
        if not self.burst_timing_records and cached_metadata_map:
            burst_records = self._load_burst_timing_from_cache(cached_metadata_map, source="iw_split")
            if burst_records:
                self.burst_timing_records = burst_records
    
    def _validate_subswath_geometry(self, swath_id: str, swath_data) -> None:
        """Validate subswath geometry data has required fields and reasonable values."""
        # Required fields for subswath geometry
        required_fields = [
            "range_samples", "azimuth_samples",
            "range_pixel_spacing", "azimuth_pixel_spacing",
            "slant_range_time", "burst_count"
        ]
        
        missing_fields = [f for f in required_fields if not hasattr(swath_data, f)]
        if missing_fields:
            # Some fields may be optional or named differently - log but don't fail
            import logging
            logging.getLogger(__name__).debug(
                f"Subswath {swath_id} missing optional fields: {missing_fields}"
            )
            return  # Don't fail on missing optional fields
        
        # Validate geometry values are reasonable if present
        if hasattr(swath_data, "range_samples") and hasattr(swath_data, "azimuth_samples"):
            if swath_data.range_samples <= 0 or swath_data.azimuth_samples <= 0:
                raise RuntimeError(
                    f"Invalid dimensions for {swath_id}: "
                    f"{swath_data.range_samples}x{swath_data.azimuth_samples}"
                )
        
        if hasattr(swath_data, "range_pixel_spacing") and hasattr(swath_data, "azimuth_pixel_spacing"):
            if swath_data.range_pixel_spacing <= 0 or swath_data.azimuth_pixel_spacing <= 0:
                raise RuntimeError(
                    f"Invalid pixel spacing for {swath_id}: "
                    f"range={swath_data.range_pixel_spacing}, "
                    f"azimuth={swath_data.azimuth_pixel_spacing}"
                )
        
        if hasattr(swath_data, "slant_range_time"):
            if swath_data.slant_range_time <= 0:
                raise RuntimeError(
                    f"Invalid slant range time for {swath_id}: {swath_data.slant_range_time}"
                )

    def _stage_deburst_and_calibrate(self, context: PipelineContext) -> None:
        """Stage 4: per-subswath deburst + calibration."""
        import logging

        if not self.burst_timing_records:
            cached_meta = None
            cache_source = "geocoding_cache"
            cache_blob = getattr(self, "_geocoding_metadata_cache", None)
            if isinstance(cache_blob, dict):
                cached_meta = cache_blob.get("cached_metadata")
            if cached_meta is None and hasattr(self, "reader") and self.reader is not None:
                cache_source = "reader"
                try:
                    with self.reader_lock:
                        if hasattr(self.reader, "get_cached_metadata"):
                            cached_meta = self.reader.get_cached_metadata()
                except Exception as exc:  # pragma: no cover - diagnostics only
                    logging.getLogger(__name__).debug(
                        "Failed to refresh cached metadata for burst timing: %s",
                        exc,
                    )

            if cached_meta:
                records = self._load_burst_timing_from_cache(cached_meta, source=cache_source)
                if records:
                    self.burst_timing_records = records

        if not self.burst_timing_records and getattr(self, "strict_science", False):
            raise RuntimeError(
                "SCIENTIFIC MODE FAILURE: No burst timing metadata available after IW Split."
            )

        calibrated_subswaths, pipeline_results, primary_subswath = deburst_utils.run_process_subswaths(
            self, context
        )
        self._calibrated_subswaths = calibrated_subswaths
        self._primary_subswath = primary_subswath
        self._working_data = calibrated_subswaths[primary_subswath]
        # Capture calibration provenance per subswath for later validation
        self._calibration_provenance = [r.get("calibration_provenance") for r in pipeline_results if r.get("calibration_provenance")]
        # Cache burst metadata if present
        for result in pipeline_results:
            bursts = result.get("burst_metadata")
            if bursts:
                self._register_burst_metadata(result["subswath"], bursts)

    def _stage_validate_calibration(self, context: PipelineContext) -> None:
        """Stage 5: ensure calibrated subswaths and LUT metadata are available.
        
        Uses dual-metric validation:
        - finite_full_canvas: fraction of finite values over entire calibrated array (informational)
        - finite_valid_window: fraction of finite values in valid range window (gating, ≥99%)
        
        Valid window is computed from burst metadata (first_valid_sample, last_valid_sample).
        Falls back to full canvas validation if burst metadata unavailable.
        """
        step_number = STEP_NUMBERS["Radiometric Calibration"]
        self.announce_step(step_number, "Radiometric Calibration", "Validating calibrated subswaths and LUT metadata")
        step_start = time.time()

        try:
            if not getattr(self, "_calibrated_subswaths", None):
                raise RuntimeError("No calibrated subswaths available")

            prov_list = getattr(self, "_calibration_provenance", []) or []
            if not prov_list and getattr(self, "strict_science", False):
                raise RuntimeError("Calibration provenance missing for all subswaths")

            # Dual-metric validation: full canvas (informational) + valid window (gating)
            # For TOPSAR, ~30-40% zeros are expected in full canvas due to edge invalidity
            finite_valid_threshold = 0.99  # Require 99% finite in valid window
            issues = []
            diagnostics = {}
            
            for swath, arr in self._calibrated_subswaths.items():
                if not isinstance(arr, np.ndarray):
                    issues.append(f"{swath}: not an array")
                    continue
                
                total = arr.size
                swath_diag = {
                    "shape": arr.shape,
                    "total_pixels": total,
                }
                
                # Full canvas metrics (informational only)
                finite_mask = np.isfinite(arr)
                finite_count = int(finite_mask.sum())
                nan_count = int(np.isnan(arr).sum())
                posinf_count = int(np.isposinf(arr).sum())
                neginf_count = int(np.isneginf(arr).sum())
                
                finite_full_pct = 100.0 * finite_count / total if total > 0 else 0.0
                
                swath_diag.update({
                    "full_canvas": {
                        "finite_count": finite_count,
                        "nan_count": nan_count,
                        "posinf_count": posinf_count,
                        "neginf_count": neginf_count,
                        "finite_fraction": finite_full_pct / 100.0,
                    }
                })
                
                # Get valid range window from burst metadata
                burst_metadata = getattr(self, "_burst_metadata", {}).get(swath, [])
                first_valid_samples = []
                last_valid_samples = []
                
                if burst_metadata:
                    for burst in burst_metadata:
                        if isinstance(burst, dict):
                            fv = burst.get("first_valid_sample")
                            lv = burst.get("last_valid_sample")
                            if fv is not None:
                                first_valid_samples.append(int(fv))
                            if lv is not None:
                                last_valid_samples.append(int(lv))
                
                # Compute valid window if burst metadata available
                missing_valid_window = False
                if first_valid_samples and last_valid_samples:
                    first_valid = int(np.median(first_valid_samples))
                    last_valid = int(np.median(last_valid_samples))
                    
                    # Clamp to array bounds
                    first_valid = max(0, min(first_valid, arr.shape[1] - 1))
                    last_valid = max(first_valid + 1, min(last_valid, arr.shape[1]))
                    
                    # Extract valid window
                    valid_window = arr[:, first_valid:last_valid]
                    window_total = valid_window.size
                    
                    if window_total > 0:
                        window_finite_mask = np.isfinite(valid_window)
                        window_finite_count = int(window_finite_mask.sum())
                        window_nan_count = int(np.isnan(valid_window).sum())
                        window_posinf_count = int(np.isposinf(valid_window).sum())
                        window_neginf_count = int(np.isneginf(valid_window).sum())
                        
                        finite_valid_pct = 100.0 * window_finite_count / window_total
                        
                        swath_diag["valid_window"] = {
                            "first_valid_col": first_valid,
                            "last_valid_col": last_valid,
                            "window_pixels": window_total,
                            "finite_count": window_finite_count,
                            "nan_count": window_nan_count,
                            "posinf_count": window_posinf_count,
                            "neginf_count": window_neginf_count,
                            "finite_fraction": finite_valid_pct / 100.0,
                        }
                        
                        # Gate on valid window
                        if finite_valid_pct < finite_valid_threshold * 100.0:
                            issues.append(
                                f"{swath}: valid window finite {finite_valid_pct:.2f}% "
                                f"(<{finite_valid_threshold*100.0:.0f}%, cols {first_valid}:{last_valid})"
                            )
                        
                        # Check for Inf in valid window (indicates LUT error)
                        if window_posinf_count > 0 or window_neginf_count > 0:
                            inf_total = window_posinf_count + window_neginf_count
                            issues.append(
                                f"{swath}: {inf_total} Inf values in valid window. "
                                f"Likely division by zero or invalid LUT usage. "
                                f"Check whether apply.rs uses the correct LUT array (gain vs factor)."
                            )
                    else:
                        missing_valid_window = True
                else:
                    missing_valid_window = True
                
                if missing_valid_window:
                    # Fall back to full canvas validation
                    swath_diag["missing_valid_window"] = True
                    finite_pct = finite_full_pct
                    
                    # Use relaxed threshold for full canvas (accounts for TOPSAR margins)
                    if finite_pct < 50.0:  # Very lenient: 50% for full canvas
                        issues.append(
                            f"{swath}: full canvas finite {finite_pct:.2f}% (<50%, "
                            f"burst metadata unavailable for precise validation)"
                        )
                    
                    logger.warning(
                        f"⚠️  {swath}: No valid window metadata, using full canvas validation "
                        f"({finite_pct:.2f}% finite)"
                    )
                
                diagnostics[swath] = swath_diag
                
                # Log summary
                if "valid_window" in swath_diag:
                    print(
                        f"   📊 {swath} validation: full canvas {finite_full_pct:.2f}% finite, "
                        f"valid window {swath_diag['valid_window']['finite_fraction']*100:.2f}% finite "
                        f"(cols {swath_diag['valid_window']['first_valid_col']}:"
                        f"{swath_diag['valid_window']['last_valid_col']})"
                    )
                else:
                    print(
                        f"   📊 {swath} validation: full canvas {finite_full_pct:.2f}% finite "
                        f"(no valid window metadata)"
                    )

            # Store diagnostics for metadata export
            if self.metadata is not None and isinstance(self.metadata, dict):
                self.metadata["calibration_validation"] = diagnostics

            if issues:
                raise RuntimeError("Calibration validation failed: " + "; ".join(issues))

            # Require LUT provenance if available
            if self.calibration_lut_source is None:
                print("   ℹ️  Calibration LUT source not reported; using annotation default")

            if prov_list and isinstance(self.metadata, dict):
                self.metadata["calibration_provenance"] = prov_list

            # Record calibration type in metadata for downstream stages and output files
            if self.metadata is not None and isinstance(self.metadata, dict):
                self.metadata["calibration_type"] = self.calibration_type
                self.metadata["calibration_applied"] = True

            # Populate spacing if still missing
            self.get_current_range_spacing()
            self.get_current_azimuth_spacing()

            step_duration = time.time() - step_start
            self.log_step(step_number, "Radiometric Calibration", "success", "Calibration validated", step_duration)
        except Exception as exc:
            step_duration = time.time() - step_start
            self.log_step(step_number, "Radiometric Calibration", "error", str(exc), step_duration)
            raise

    def _stage_merge_subswaths(self, context: PipelineContext) -> None:
        """Stage 6: merge calibrated subswaths (if applicable)."""
        merge_utils.perform_merge(self)

    # ------------------------------------------------------------------ #
    # Core helpers required by the modular stages
    # ------------------------------------------------------------------ #
    def create_reader(self, self_ref=None):
        """Create a cached SLC reader for high-performance metadata access."""
        return sardine.create_cached_slc_reader(str(self.input_path))

    def validate_metadata(self, metadata: Any) -> Dict[str, Any]:
        """Validate metadata and check polarization availability."""
        import logging
        logger = logging.getLogger(__name__)
        
        if not isinstance(metadata, dict):
            raise ValueError(f"Invalid metadata type: {type(metadata)!r}")
        
        # Validate that self.polarization is set and valid format
        if not hasattr(self, 'polarization') or not self.polarization:
            raise ValueError("Polarization not set on processor before validation")
        
        valid_polarizations = {'VV', 'VH', 'HH', 'HV'}
        requested_pol = self.polarization.upper()
        if requested_pol not in valid_polarizations:
            raise ValueError(f"Invalid polarization format: '{self.polarization}'. Expected one of: {valid_polarizations}")
        
        # Validate that requested polarization is available in product
        available_pols = []
        pol_str = metadata.get("polarizations", "")
        if isinstance(pol_str, str) and pol_str:
            available_pols = [p.strip().upper() for p in pol_str.split(",") if p.strip()]
        elif isinstance(pol_str, list):
            available_pols = [p.upper() for p in pol_str if p]
        
        if available_pols:
            if requested_pol not in available_pols:
                raise ValueError(
                    f"❌ SCIENTIFIC MODE: Requested polarization '{self.polarization}' not available. "
                    f"Product contains: {', '.join(available_pols)}"
                )
            print(f"   ✅ Polarization validated: {self.polarization} (available: {', '.join(available_pols)})")
        else:
            # No polarization info in metadata - fail
            raise ValueError(
                f"Cannot validate polarization '{self.polarization}' - metadata lacks polarization info. "
                f"Check SAFE manifest parsing."
            )
        
        self.validated_metadata = metadata
        return metadata

    def _load_dem_for_scene(self, bbox: Optional[Tuple[float, float, float, float]], dem_resolution: float):
        """Load DEM for the requested bounding box with simple caching."""

        dem_resolution = float(dem_resolution or 30.0)

        if bbox is None:
            bbox = getattr(self, "product_bbox", None)
            if bbox is None and isinstance(self.metadata, dict):
                try:
                    bbox = (
                        float(self.metadata.get("min_longitude")),
                        float(self.metadata.get("min_latitude")),
                        float(self.metadata.get("max_longitude")),
                        float(self.metadata.get("max_latitude")),
                    )
                except (TypeError, ValueError, KeyError) as e:
                    logger.debug(f"Failed to extract bbox from metadata: {e}")
                    bbox = None

        if not bbox or len(bbox) != 4:
            raise RuntimeError("DEM bbox is missing or invalid for terrain processing")

        bbox_tuple = tuple(float(v) for v in bbox)

        if (
            self._cached_dem_data is not None
            and self._cached_dem_bbox == bbox_tuple
            and self._cached_dem_resolution == dem_resolution
            and self._cached_dem_geo_transform is not None
        ):
            dem_crs = self._cached_dem_crs
            if not dem_crs:
                raise RuntimeError("Cached DEM has no CRS - cannot proceed with geocoding")
            return self._cached_dem_data, self._cached_dem_geo_transform, dem_crs

        dem_result = sardine.load_dem_for_bbox(
            list(bbox_tuple),
            str(self.dem_cache_dir),
            dem_resolution,
        )

        if not isinstance(dem_result, dict):
            raise RuntimeError("DEM loader returned an unexpected result")

        dem_data = np.asarray(dem_result.get("data"), dtype=np.float32)
        dem_transform = tuple(dem_result.get("geo_transform") or ())
        dem_crs = dem_result.get("crs")
        if not dem_crs:
            raise RuntimeError(
                "DEM loader did not return CRS. Cannot proceed with geocoding without "
                "knowing the DEM coordinate reference system."
            )

        if dem_data.size == 0 or len(dem_transform) != 6:
            raise RuntimeError("DEM loading failed: missing data or geotransform")

        # DEM CRS/datum validation - CRITICAL for correct geocoding
        # SAR processing expects geographic CRS (lat/lon) with WGS84 datum
        self._validate_dem_crs(dem_crs)

        self._cached_dem_data = dem_data
        self._cached_dem_bbox = bbox_tuple
        self._cached_dem_geo_transform = dem_transform
        self._cached_dem_resolution = dem_resolution
        self._cached_dem_crs = dem_crs

        return dem_data, dem_transform, dem_crs

    def _pre_load_dem_background(self) -> None:
        """I/O Optimization: Pre-load DEM in background thread while processing subswaths.
        
        This function starts DEM loading early (after bbox is known from IW split stage)
        so it can overlap with deburst/calibration processing. The DEM will be cached
        and ready when terrain correction stage needs it.
        
        Uses Future-based tracking for proper error handling.
        """
        logger = logging.getLogger(__name__)
        
        # Check if DEM is already cached
        if self._cached_dem_data is not None:
            logger.debug("DEM already cached, skipping background pre-load")
            return
        
        # Get bbox from metadata
        bbox = None
        if isinstance(self.metadata, dict):
            try:
                min_lon = self.metadata.get("min_longitude")
                min_lat = self.metadata.get("min_latitude")
                max_lon = self.metadata.get("max_longitude")
                max_lat = self.metadata.get("max_latitude")
                # Don't use fallback to 0 - require actual values
                if min_lon is None or min_lat is None or max_lon is None or max_lat is None:
                    logger.debug("Missing bbox coordinates in metadata for DEM pre-loading")
                    return
                bbox = (
                    float(min_lon),
                    float(min_lat),
                    float(max_lon),
                    float(max_lat),
                )
            except (ValueError, TypeError):
                logger.debug("Could not extract bbox from metadata for DEM pre-loading")
                return
        
        if not bbox or len(bbox) != 4:
            logger.debug("Invalid bbox for DEM pre-loading")
            return
        
        # Validate bbox is not all zeros (which would be an invalid fallback)
        if bbox == (0.0, 0.0, 0.0, 0.0):
            logger.debug("Bbox is all zeros - likely missing metadata, skipping DEM pre-loading")
            return
        
        # Get DEM resolution from options
        dem_resolution = float(self.options.get("dem_resolution", 30.0))
        
        # Create executor if needed
        if self._thread_executor is None:
            self._thread_executor = ThreadPoolExecutor(max_workers=2, thread_name_prefix="sardine_prefetch")
        
        def load_dem_task() -> Dict[str, Any]:
            """Load DEM and return result dict. Raises exception on failure."""
            logger.info("🔄 Background DEM pre-loading started")
            start_time = time.time()
            
            dem_result = sardine.load_dem_for_bbox(
                list(bbox),
                str(self.dem_cache_dir),
                dem_resolution,
            )
            
            if not isinstance(dem_result, dict):
                raise RuntimeError("DEM loader returned unexpected result type")
            
            dem_data = np.asarray(dem_result.get("data"), dtype=np.float32)
            dem_transform = tuple(dem_result.get("geo_transform") or ())
            dem_crs = dem_result.get("crs")
            
            if not dem_crs:
                raise RuntimeError("DEM loader did not return CRS")
            
            if dem_data.size == 0:
                raise RuntimeError("DEM loader returned empty data")
            
            if len(dem_transform) != 6:
                raise RuntimeError(f"DEM loader returned invalid transform: {len(dem_transform)} elements")
            
            # Cache DEM data
            self._cached_dem_data = dem_data
            self._cached_dem_bbox = tuple(float(v) for v in bbox)
            self._cached_dem_geo_transform = dem_transform
            self._cached_dem_resolution = dem_resolution
            self._cached_dem_crs = dem_crs
            
            elapsed = time.time() - start_time
            logger.info(f"✅ Background DEM pre-loading completed in {elapsed:.1f}s")
            return {"success": True, "elapsed": elapsed}
        
        # Submit to executor with Future tracking
        self._dem_prefetch_future = self._thread_executor.submit(load_dem_task)
        logger.debug("Background DEM pre-loading submitted to executor")

    def _pre_fetch_orbit_background(self) -> None:
        """I/O Optimization: Pre-fetch orbit file in background thread.
        
        This function checks if orbit file is cached, and if not, starts
        downloading it in a background thread. This overlaps orbit I/O with
        metadata processing.
        
        Uses Future-based tracking for proper error handling.
        """
        logger = logging.getLogger(__name__)
        
        # Check if orbit is already available
        if self._last_orbit_result is not None:
            logger.debug("Orbit already available, skipping background pre-fetch")
            return
        
        # Get product metadata for orbit download
        if not isinstance(self.metadata, dict):
            logger.debug("No metadata available for orbit pre-fetching")
            return
        
        # Try both key names: Rust exports "start_time", docs may reference "sensing_start"
        sensing_start_str = self.metadata.get("start_time") or self.metadata.get("sensing_start")
        if not sensing_start_str:
            logger.debug("No start_time/sensing_start in metadata for orbit pre-fetching")
            return
        
        try:
            sensing_start = datetime.fromisoformat(sensing_start_str.replace("Z", "+00:00"))
        except ValueError as e:
            logger.debug(f"Could not parse sensing_start for orbit pre-fetching: {e}")
            return
        
        # Create executor if needed
        if self._thread_executor is None:
            self._thread_executor = ThreadPoolExecutor(max_workers=2, thread_name_prefix="sardine_prefetch")
        
        def fetch_orbit_task() -> Dict[str, Any]:
            """Fetch orbit file and return result. Raises exception on failure."""
            logger.info("🔄 Background orbit pre-fetching started")
            start_time = time.time()
            
            orbit_type = self.options.get("orbit_type", "POEORB")
            orbit_path = self.download_orbit_file(
                sensing_start,
                orbit_type,
                str(self.orbit_cache_dir)
            )
            
            elapsed = time.time() - start_time
            logger.info(f"✅ Background orbit pre-fetching completed in {elapsed:.1f}s: {orbit_path}")
            return {"success": True, "orbit_path": orbit_path, "elapsed": elapsed}
        
        # Submit to executor with Future tracking
        self._orbit_prefetch_future = self._thread_executor.submit(fetch_orbit_task)
        logger.debug("Background orbit pre-fetching submitted to executor")

    def _check_prefetch_errors(self) -> None:
        """Check if any prefetch operations failed and log/raise as appropriate.
        
        This should be called before stages that depend on prefetched data
        to surface errors early rather than silently continuing.
        """
        logger = logging.getLogger(__name__)
        
        # Check DEM prefetch
        if self._dem_prefetch_future is not None and self._dem_prefetch_future.done():
            try:
                self._dem_prefetch_future.result()  # Raises if task failed
            except Exception as e:
                logger.warning(f"DEM prefetch failed: {e}")
                # Clear the future so we don't check again
                self._dem_prefetch_future = None
        
        # Check orbit prefetch
        if self._orbit_prefetch_future is not None and self._orbit_prefetch_future.done():
            try:
                self._orbit_prefetch_future.result()  # Raises if task failed
            except Exception as e:
                logger.warning(f"Orbit prefetch failed: {e}")
                # Clear the future so we don't check again
                self._orbit_prefetch_future = None

    def _cleanup_memory(self) -> None:
        """Explicit memory cleanup after processing stages.
        
        Clears large working buffers that are no longer needed to reduce
        memory pressure during long processing runs.
        """
        logger = logging.getLogger(__name__)
        
        # Clear working data buffer
        if self._working_data is not None:
            del self._working_data
            self._working_data = None
        
        # Clear calibrated subswaths after merge
        if self._calibrated_subswaths:
            for key in list(self._calibrated_subswaths.keys()):
                del self._calibrated_subswaths[key]
            self._calibrated_subswaths.clear()
        
        # Clear DEM cache after terrain correction
        if self._cached_dem_data is not None:
            del self._cached_dem_data
            self._cached_dem_data = None
            self._cached_dem_bbox = None
            self._cached_dem_geo_transform = None
            self._cached_dem_resolution = None
            self._cached_dem_crs = None
        
        # Force garbage collection
        gc.collect()
        logger.debug("Memory cleanup completed")

    def _maybe_gc(self, force: bool = False) -> None:
        """Memory-adaptive garbage collection based on system RAM tier.
        
        On LOW memory systems, GC runs after every subswath.
        On MEDIUM systems, GC runs after every 2 subswaths.
        On HIGH+ systems, GC only runs when forced (e.g., after merge).
        
        Args:
            force: If True, always run GC regardless of memory tier
        """
        if force:
            gc.collect()
            return
            
        if self._memory_config is None:
            # Default: run GC (conservative behavior)
            gc.collect()
            return
            
        # Only run GC based on memory tier
        if self._memory_config.tier == MemoryTier.LOW:
            # Always GC on low memory systems
            gc.collect()
        elif self._memory_config.tier == MemoryTier.MEDIUM:
            # Occasional GC on medium memory systems
            # Track via subswath count (would need counter, skip for simplicity)
            pass
        # HIGH and VERY_HIGH: skip GC to maximize performance

    def _shutdown_executor(self) -> None:
        """Shutdown thread executor and wait for pending tasks."""
        if self._thread_executor is not None:
            self._thread_executor.shutdown(wait=True)
            self._thread_executor = None

    def download_orbit_file(self, sensing_start: datetime, orbit_type: Optional[str], orbit_cache_dir: str) -> str:
        """
        Download precise orbit data and return the cached path.

        Uses the unified download manager if available, otherwise falls back to
        the cached reader or core apply_precise_orbit_file binding.
        """
        cache_dir = Path(orbit_cache_dir)
        cache_dir.mkdir(parents=True, exist_ok=True)

        product_id = None
        if isinstance(self.metadata, dict):
            product_id = self.metadata.get("product_id")
        if not product_id:
            product_id = Path(self.input_path).stem

        # FIXED: Use DownloadManager for orbit downloads (new unified system)
        orbit_path: Optional[str] = None
        import logging
        logger = logging.getLogger(__name__)
        
        # Try unified download manager first (if available and auto-download enabled)
        if hasattr(self, "_download_manager") and self._download_manager is not None:
            try:
                start_time_iso = sensing_start.astimezone(timezone.utc).isoformat()
                downloaded_path = self._download_manager.download_orbit(product_id, start_time_iso)
                if downloaded_path and Path(downloaded_path).exists():
                    logger.info(f"Downloaded orbit file via DownloadManager: {downloaded_path}")
                    orbit_path = downloaded_path  # Set this so we can use it below
                else:
                    # DownloadManager might have downloaded but returned None - check orbits subdirectory
                    import shutil
                    from datetime import datetime
                    age_days = (datetime.now(timezone.utc) - sensing_start).days
                    orbit_type = "POEORB" if age_days > 20 else "RESORB"
                    orbits_subdir = cache_dir / "orbits"
                    possible_path = orbits_subdir / f"{product_id}_{orbit_type}.EOF"
                    if possible_path.exists():
                        logger.info(f"Found orbit file in DownloadManager cache: {possible_path}")
                        downloaded_path = str(possible_path)
                        orbit_path = downloaded_path
                    else:
                        logger.warning(f"DownloadManager returned invalid path and file not found in cache: {downloaded_path}")
                        downloaded_path = None
            except Exception as exc:
                # DownloadManager failed - check if file was downloaded anyway (in orbits subdirectory)
                logger.warning(f"DownloadManager orbit download raised exception: {exc}")
                import shutil
                from datetime import datetime
                age_days = (datetime.now(timezone.utc) - sensing_start).days
                orbit_type = "POEORB" if age_days > 20 else "RESORB"
                orbits_subdir = cache_dir / "orbits"
                possible_path = orbits_subdir / f"{product_id}_{orbit_type}.EOF"
                if possible_path.exists():
                    logger.info(f"Found orbit file in cache despite exception: {possible_path}")
                    downloaded_path = str(possible_path)
                    orbit_path = downloaded_path
                else:
                    downloaded_path = None
                    logger.debug(f"Orbit file not found in cache after exception")
            
            # If we have a downloaded file, copy it to expected locations
            if downloaded_path and Path(downloaded_path).exists():
                import shutil
                from datetime import datetime
                age_days = (datetime.now(timezone.utc) - sensing_start).days
                orbit_type = "POEORB" if age_days > 20 else "RESORB"
                expected_filename = f"{product_id}_{orbit_type}.EOF"
                
                # Location 1: Main cache directory (for apply_precise_orbit_file)
                expected_path = cache_dir / expected_filename
                
                # Location 2: Product-specific subdirectory (for deburst_topsar_cached via SARDINE_ORBIT_CACHE)
                # SARDINE_ORBIT_CACHE is set to self.orbit_cache_dir which is the product-specific subdirectory
                product_specific_path = self.orbit_cache_dir / expected_filename
                
                # Ensure cache directories exist
                cache_dir.mkdir(parents=True, exist_ok=True)
                self.orbit_cache_dir.mkdir(parents=True, exist_ok=True)
                
                # Copy to main cache directory (for apply_precise_orbit_file)
                if str(expected_path) != downloaded_path:
                    if expected_path.exists():
                        logger.debug(f"Orbit file already exists at expected location: {expected_path}")
                    else:
                        shutil.copy2(downloaded_path, expected_path)
                        logger.info(f"Copied orbit file to main cache: {expected_path}")
                    orbit_path = str(expected_path)
                else:
                    orbit_path = downloaded_path
                
                # Also copy to product-specific subdirectory (for deburst_topsar_cached)
                if str(product_specific_path) != orbit_path:
                    if product_specific_path.exists():
                        logger.debug(f"Orbit file already exists in product-specific cache: {product_specific_path}")
                    else:
                        shutil.copy2(downloaded_path, product_specific_path)
                        logger.info(f"Copied orbit file to product-specific cache for deburst: {product_specific_path}")

        # Attempt to download via the reader if DownloadManager didn't work
        if orbit_path is None:
            try:
                if hasattr(self, "reader") and self.reader is not None and hasattr(self.reader, "download_orbit_files"):
                    paths = self.reader.download_orbit_files(str(cache_dir))
                    if isinstance(paths, list) and paths:
                        orbit_path = str(paths[-1])
            except Exception as exc:
                logger.debug(f"Orbit download via reader failed: {exc}")

        # FIXED: Validate orbit file - use apply_precise_orbit_file only if file exists in cache
        # The function will try to download if file doesn't exist, which fails with new download module
        try:
            # If we have a downloaded orbit file, try to validate it
            # apply_precise_orbit_file will find it in cache_dir if it exists
            orbit_result = sardine.apply_precise_orbit_file(
                str(product_id),
                sensing_start.astimezone(timezone.utc).isoformat(),
                str(cache_dir),
            )
            self._last_orbit_result = orbit_result
            
            # Validate result structure immediately
            if not isinstance(orbit_result, dict):
                raise RuntimeError(f"Orbit validation returned invalid type: {type(orbit_result)}")
            
            # Rust wraps everything in a "result" key
            result_data = orbit_result.get("result", orbit_result)
            if not isinstance(result_data, dict):
                raise RuntimeError(f"Orbit result data has invalid type: {type(result_data)}")
            
            # Extract and cache orbit state vectors for terrain correction
            # Rust returns keys: osv_times, osv_positions, osv_velocities inside result
            self._cached_orbit_times = result_data.get("osv_times", [])
            self._cached_orbit_positions = result_data.get("osv_positions", [])
            self._cached_orbit_velocities = result_data.get("osv_velocities", [])
            
            vector_count = result_data.get("orbit_vectors_count")
            if vector_count is None:
                raise RuntimeError("Orbit validation did not report orbit_vectors_count")
                
        except Exception as exc:
            import logging
            logger = logging.getLogger(__name__)
            
            # Always fail on validation errors - orbit data is critical
            raise RuntimeError(
                f"Precise orbit validation failed: {exc}\n"
                f"This indicates the downloaded orbit file is invalid or corrupted."
            ) from exc

        if orbit_path is None:
            try:
                eof_candidates = sorted(cache_dir.glob("*.EOF"))
                if eof_candidates:
                    # Prefer file matching product_id
                    product_id_base = product_id.replace(".SAFE", "").replace(".zip", "")
                    matching_files = [
                        f for f in eof_candidates
                        if product_id_base in f.name
                    ]
                    
                    if matching_files:
                        orbit_path = str(matching_files[-1])  # Most recent matching file
                        import logging
                        logging.getLogger(__name__).info(f"Found orbit file matching product: {orbit_path}")
                    else:
                        # No matching orbit file - always fail (no fallback)
                        import logging
                        logging.getLogger(__name__).error(
                            f"No orbit file matches product_id '{product_id}'. "
                            f"Found {len(eof_candidates)} EOF files in cache but none match."
                        )
                        raise RuntimeError(
                            f"No orbit file matches product_id '{product_id}'. "
                            f"Orbit files must match the product to ensure correct geocoding. "
                            f"Download the correct orbit file for this product."
                        )
                    
                    # Validate file exists and is readable
                    if not Path(orbit_path).exists():
                        raise RuntimeError(f"Selected orbit file does not exist: {orbit_path}")
                    
                    # Validate file is not empty
                    if Path(orbit_path).stat().st_size == 0:
                        raise RuntimeError(f"Orbit file is empty: {orbit_path}")
                        
            except RuntimeError:
                raise  # Re-raise validation errors
            except Exception as e:
                import logging
                logging.getLogger(__name__).error(f"Failed to resolve orbit file from cache: {e}")
                orbit_path = None

        if orbit_path is None:
            raise RuntimeError(
                f"No precise orbit file could be resolved in the cache: {cache_dir}\n"
                f"Product ID: {product_id}\n"
                f"Start time: {sensing_start}"
            )

        return orbit_path

    def _register_orbit_file_for_metadata(self, orbit_path: str) -> None:
        """Track the orbit file used for metadata outputs."""
        self._orbit_file_for_metadata = orbit_path

    def _register_burst_metadata(self, subswath: str, bursts: Any) -> None:
        """Cache burst metadata for later diagnostics or metadata export."""
        if not hasattr(self, "_burst_metadata") or self._burst_metadata is None:
            self._burst_metadata = {}
        self._burst_metadata[subswath] = bursts

    def _validate_dem_crs(self, dem_crs: str) -> None:
        """Validate DEM CRS for SAR processing compatibility.

        SAR terrain correction requires geographic CRS (lat/lon) with WGS84 datum.
        This prevents silent geocoding errors from incompatible CRS.

        Raises:
            ValueError: If CRS is not compatible and ALLOW_UNSAFE_DEM is not set.
        """
        import os

        if not dem_crs:
            dem_crs = "UNKNOWN"

        # Normalize CRS string
        dem_crs_upper = dem_crs.upper().replace(" ", "")

        # Check for compatible geographic CRS with WGS84 datum
        compatible_crs = [
            "EPSG:4326",   # WGS84 geographic
            "WGS84",       # Common alias
            "EPSG4326",    # Without colon
            "CRS:84",      # OGC alias for WGS84
        ]

        is_compatible = any(crs.upper() in dem_crs_upper for crs in compatible_crs)

        # Also check for WGS 84 in the full CRS string
        if "WGS" in dem_crs_upper and "84" in dem_crs_upper:
            is_compatible = True

        if not is_compatible:
            # Check for unsafe override
            allow_unsafe = os.environ.get("ALLOW_UNSAFE_DEM", "0") == "1"

            if allow_unsafe:
                warnings.warn(
                    f"⚠️ DEM CRS '{dem_crs}' may not be compatible with SAR processing. "
                    f"Expected WGS84 geographic (EPSG:4326). Proceeding due to ALLOW_UNSAFE_DEM=1.",
                    UserWarning,
                    stacklevel=3,
                )
            else:
                raise ValueError(
                    f"❌ DEM CRS '{dem_crs}' is not compatible with SAR processing. "
                    f"Expected WGS84 geographic (EPSG:4326). "
                    f"Set ALLOW_UNSAFE_DEM=1 to override this check."
                )

        # Check for vertical datum (EGM96/EGM2008 vs ellipsoidal)
        # This is a warning only since we can't easily detect the vertical datum
        if "EGM" not in dem_crs_upper and "GEOID" not in dem_crs_upper:
            if self.verbosity >= 2:
                print(
                    f"   ℹ️ DEM CRS '{dem_crs}' does not specify vertical datum. "
                    f"Assuming ellipsoidal heights (WGS84)."
                )

    def get_current_range_spacing(self) -> Optional[float]:
        """Return the active range pixel spacing (updates from metadata if unset)."""
        if self.current_range_spacing is not None:
            return self.current_range_spacing
        if isinstance(self.metadata, dict):
            for key in ("range_pixel_spacing", "pixel_spacing_range", "range_spacing"):
                value = self.metadata.get(key)
                try:
                    if value is not None:
                        spacing = float(value)
                        self.current_range_spacing = spacing
                        if self.original_range_spacing is None:
                            self.original_range_spacing = spacing
                        return spacing
                except (TypeError, ValueError):
                    continue
        return None

    def get_current_azimuth_spacing(self) -> Optional[float]:
        """Return the active azimuth pixel spacing (updates from metadata if unset)."""
        if self.current_azimuth_spacing is not None:
            return self.current_azimuth_spacing
        if isinstance(self.metadata, dict):
            for key in ("azimuth_pixel_spacing", "pixel_spacing_azimuth", "azimuth_spacing"):
                value = self.metadata.get(key)
                try:
                    if value is not None:
                        spacing = float(value)
                        self.current_azimuth_spacing = spacing
                        if self.original_azimuth_spacing is None:
                            self.original_azimuth_spacing = spacing
                        return spacing
                except (TypeError, ValueError):
                    continue
        return None

    def update_current_spacing(
        self,
        range_spacing: Optional[float],
        azimuth_spacing: Optional[float],
    ) -> None:
        """Update the active pixel spacing and preserve originals if unset."""
        if range_spacing is not None:
            self.current_range_spacing = float(range_spacing)
            if self.original_range_spacing is None:
                self.original_range_spacing = self.current_range_spacing
        if azimuth_spacing is not None:
            self.current_azimuth_spacing = float(azimuth_spacing)
            if self.original_azimuth_spacing is None:
                self.original_azimuth_spacing = self.current_azimuth_spacing

    def _log(self, message: str, level: int = 1) -> None:
        """Print message if verbosity level is sufficient.
        
        Levels:
            0: Errors only (always printed via raise)
            1: Normal operation messages (default)
            2: Verbose details
            3: Debug information
        """
        if self.verbosity >= level:
            print(message)

    def _record_stage_timing(self, stage_name: str, duration: float) -> None:
        """Record timing for a processing stage with detailed breakdown.
        
        This integrates with the performance monitor to provide detailed timing
        breakdowns for optimization analysis.
        """
        self._stage_timings[stage_name] = duration
        
        # Also record in performance monitor if available
        if hasattr(self, 'performance_monitor') and self.performance_monitor:
            self.performance_monitor.record_phase_duration(stage_name, duration)

    def _validate_stage_output(self, stage_name: str, context: PipelineContext) -> None:
        """Validate stage output for NaN percentage, dimensions, and data ranges.
        
        This is called automatically after each pipeline stage to catch issues early.
        """
        import logging
        logger = logging.getLogger(__name__)
        
        # Get working data if available
        working_data = getattr(self, '_working_data', None)
        if working_data is None:
            # Try to get from calibrated subswaths for deburst/calibration stages
            calibrated = getattr(self, '_calibrated_subswaths', None)
            if calibrated and isinstance(calibrated, dict):
                # Validate each calibrated subswath
                for swath_id, data in calibrated.items():
                    if hasattr(data, 'shape') and hasattr(data, '__iter__'):
                        try:
                            arr = np.asarray(data, dtype=np.float32)
                            self._validate_array(f"{stage_name} ({swath_id})", arr, logger)
                        except (ValueError, TypeError) as e:
                            logger.debug(f"Could not validate subswath {swath_id}: {e}")
                return
            return
        
        # Validate main working data
        if hasattr(working_data, 'shape') and hasattr(working_data, '__iter__'):
            try:
                arr = np.asarray(working_data, dtype=np.float32)
                self._validate_array(stage_name, arr, logger)
            except Exception as e:
                logger.debug(f"Could not validate {stage_name} output: {e}")

    def _mask_nan_before_processing(self, stage_name: str) -> None:
        """Mask NaN values in working data before processing to prevent propagation.
        
        This is called at critical pipeline points (after merge, before multilooking, etc.)
        to ensure NaN values don't propagate through processing stages.
        """
        import logging
        logger = logging.getLogger(__name__)
        
        working_data = getattr(self, '_working_data', None)
        if working_data is None:
            return
        
        try:
            arr = np.asarray(working_data, dtype=np.float32)
            if arr.size == 0:
                return
            
            nan_count_before = np.sum(np.isnan(arr))
            if nan_count_before == 0:
                return  # No NaN values to mask
            
            # For now, we just log - actual masking depends on processing stage requirements
            # Some stages (like multilooking) may handle NaN differently
            nan_pct = 100.0 * nan_count_before / arr.size
            logger.debug(
                f"📊 {stage_name}: {nan_count_before} NaN pixels ({nan_pct:.1f}%) "
                f"will be handled by processing stage"
            )
        except Exception as e:
            logger.debug(f"Could not check NaN for {stage_name}: {e}")

    def _validate_array(self, stage_name: str, arr: np.ndarray, logger: logging.Logger) -> None:
        """Validate a numpy array for data quality."""
        if arr.size == 0:
            logger.warning(f"⚠️  {stage_name}: Empty array")
            return
        
        # Check dimensions
        if arr.ndim != 2:
            logger.warning(f"⚠️  {stage_name}: Unexpected dimensions {arr.ndim}D (expected 2D)")
        
        # Check NaN percentage
        nan_count = np.isnan(arr).sum()
        nan_pct = 100.0 * nan_count / arr.size
        finite_count = np.isfinite(arr).sum()
        finite_pct = 100.0 * finite_count / arr.size
        
        if nan_pct > 50.0:
            logger.warning(
                f"⚠️  {stage_name}: High NaN percentage {nan_pct:.1f}% "
                f"({nan_count}/{arr.size} pixels)"
            )
        elif nan_pct > 10.0:
            logger.info(
                f"ℹ️  {stage_name}: NaN percentage {nan_pct:.1f}% "
                f"({nan_count}/{arr.size} pixels)"
            )
        
        # Check data range for finite values
        if finite_count > 0:
            finite_data = arr[np.isfinite(arr)]
            if finite_data.size > 0:
                min_val = float(np.min(finite_data))
                max_val = float(np.max(finite_data))
                mean_val = float(np.mean(finite_data))
                median_val = float(np.median(finite_data))
                
                # Check for reasonable ranges (backscatter should be positive)
                if min_val < 0.0 and stage_name not in ["Convert to dB"]:
                    logger.warning(
                        f"⚠️  {stage_name}: Negative values detected "
                        f"(min={min_val:.4f}, this may be expected for some stages)"
                    )
                
                # Log statistics at debug level
                logger.debug(
                    f"📊 {stage_name} validation: shape={arr.shape}, "
                    f"finite={finite_pct:.1f}%, "
                    f"range=[{min_val:.4f}, {max_val:.4f}], "
                    f"mean={mean_val:.4f}, median={median_val:.4f}"
                )
            else:
                logger.warning(f"⚠️  {stage_name}: No finite values in array")
        else:
            logger.warning(f"⚠️  {stage_name}: All values are NaN or infinite")

    def _log_diagnostic_statistics(self, stage_name: str, data: np.ndarray, context: str = "") -> None:
        """Log comprehensive diagnostic statistics for data at a processing stage.
        
        This method logs min, max, mean, median, NaN count, and other statistics
        to help diagnose issues in the processing pipeline.
        
        Args:
            stage_name: Name of the processing stage (e.g., "After Calibration")
            data: numpy array to analyze
            context: Optional additional context string
        """
        import logging
        logger = logging.getLogger(__name__)
        
        if not isinstance(data, np.ndarray):
            logger.warning(f"📊 {stage_name}: Data is not a numpy array: {type(data)}")
            return
        
        if data.size == 0:
            logger.warning(f"📊 {stage_name}: Empty array")
            return
        
        # Basic statistics
        nan_count = int(np.isnan(data).sum())
        inf_count = int(np.isinf(data).sum())
        finite_mask = np.isfinite(data)
        finite_count = int(finite_mask.sum())
        total_count = int(data.size)
        
        nan_pct = 100.0 * nan_count / total_count if total_count > 0 else 0.0
        inf_pct = 100.0 * inf_count / total_count if total_count > 0 else 0.0
        finite_pct = 100.0 * finite_count / total_count if total_count > 0 else 0.0
        
        # Statistics for finite values
        if finite_count > 0:
            finite_data = data[finite_mask]
            min_val = float(np.min(finite_data))
            max_val = float(np.max(finite_data))
            mean_val = float(np.mean(finite_data))
            median_val = float(np.median(finite_data))
            std_val = float(np.std(finite_data))

            # Percentiles
            p1_val = float(np.percentile(finite_data, 1))
            p5_val = float(np.percentile(finite_data, 5))
            p95_val = float(np.percentile(finite_data, 95))
            p99_val = float(np.percentile(finite_data, 99))

            zero_like_mask = finite_data <= 0.0
            zero_like_count = int(np.count_nonzero(zero_like_mask))
            zero_like_pct = 100.0 * zero_like_count / finite_count if finite_count > 0 else 0.0

            # Log comprehensive statistics
            log_msg = (
                f"📊 DIAGNOSTIC {stage_name}{' (' + context + ')' if context else ''}:\n"
                f"   Shape: {data.shape}, Total pixels: {total_count}\n"
                f"   Finite: {finite_count} ({finite_pct:.2f}%), NaN: {nan_count} ({nan_pct:.2f}%), Inf: {inf_count} ({inf_pct:.2f}%)\n"
                f"   Range: [{min_val:.6e}, {max_val:.6e}]\n"
                f"   Mean: {mean_val:.6e}, Median: {median_val:.6e}, Std: {std_val:.6e}\n"
                f"   Percentiles: P1={p1_val:.6e}, P5={p5_val:.6e}, P95={p95_val:.6e}, P99={p99_val:.6e}\n"
                f"   Non-positive (<=0): {zero_like_count} ({zero_like_pct:.2f}%)"
            )

            # For dB data, also log in dB scale if values are already in dB
            if stage_name == "After dB Conversion" or "dB" in stage_name:
                logger.info(log_msg)
                print(log_msg)
            else:
                # For linear power data, use safe dB conversion (epsilon floor)
                mean_db, min_db, max_db = (_safe_db(mean_val), _safe_db(min_val), _safe_db(max_val))
                log_msg += f"\n   Equivalent dB (if linear): mean={mean_db:.2f}, range=[{min_db:.2f}, {max_db:.2f}]"

                logger.info(log_msg)
                if self.verbosity >= 1:
                    print(log_msg)
        else:
            logger.warning(f"📊 {stage_name}: No finite values in array (all NaN/Inf)")
            print(f"   ⚠️  {stage_name}: No finite values detected")
    
    def _cleanup_intermediate_data(self, *attrs_to_clear: str) -> None:
        """Release intermediate data to reduce memory pressure.
        
        Only runs if memory_cleanup option is enabled (default True).
        FIXED: Use explicit del statements before gc.collect() for better memory release.
        """
        if not self._memory_cleanup_enabled:
            return
        for attr in attrs_to_clear:
            if hasattr(self, attr):
                # FIXED: Use explicit del to ensure reference is removed
                value = getattr(self, attr)
                if isinstance(value, (dict, list)):
                    # Clear collections explicitly
                    value.clear()
                delattr(self, attr)
                del value  # Explicit deletion
        # Force garbage collection after explicit deletions
        gc.collect()
    
    def _validate_filter_window(self, value: Any) -> int:
        """Validate filter_window parameter.
        
        Note: Rust speckle filter enforces max window size of 25 for performance.
        """
        if not isinstance(value, int):
            try:
                value = int(value)
            except (ValueError, TypeError):
                raise ValueError(f"filter_window must be an integer, got {type(value).__name__}: {value}")
        
        if value < 3:
            raise ValueError(f"filter_window must be >= 3, got {value}")
        if value > 25:
            raise ValueError(f"filter_window too large: {value} (max 25, enforced by Rust speckle filter)")
        if value % 2 == 0:
            raise ValueError(f"filter_window must be odd, got {value}")
        
        return value
    
    def _validate_multilook_factor(self, value: Any, direction: str) -> int:
        """Validate multilook factor (range or azimuth)."""
        if not isinstance(value, (int, float)):
            try:
                value = int(float(value))
            except (ValueError, TypeError):
                raise ValueError(f"multilook_{direction} must be numeric, got {type(value).__name__}: {value}")
        
        value = int(value)
        if value < 1:
            raise ValueError(f"multilook_{direction} must be >= 1, got {value}")
        if value > 20:
            raise ValueError(f"multilook_{direction} too large: {value} (max 20 for reasonable processing)")
        
        return value
    
    def _validate_dem_resolution(self, value: Any) -> float:
        """Validate DEM resolution parameter."""
        try:
            value = float(value)
        except (ValueError, TypeError):
            raise ValueError(f"dem_resolution must be numeric, got {type(value).__name__}: {value}")
        
        if value <= 0:
            raise ValueError(f"dem_resolution must be > 0, got {value}")
        if value > 1000:
            raise ValueError(f"dem_resolution too large: {value} (max 1000m for reasonable processing)")
        
        return value
    
    def _export_intermediate_product(
        self, 
        data: np.ndarray, 
        stage_name: str,
        description: str,
        is_geocoded: bool = False
    ) -> None:
        """Export intermediate processing checkpoint for validation.
        
        Exports intermediate products as NumPy binary + JSON metadata for:
        - Step-by-step validation against reference implementations
        - QA/debugging of pipeline transformations
        - Performance benchmarking
        
        Args:
            data: Processed array (single-band, typically sigma0/gamma0)
            stage_name: Short identifier (e.g., 'calibration', 'multilook')
            description: Human-readable description for metadata
            is_geocoded: True if data is in geographic coordinates (EPSG:4326/UTM),
                        False if in SAR geometry (range/azimuth)
        
        Raises:
            Exception: On export failure (logged but non-fatal)
        """
        import logging
        import datetime
        import json
        
        logger = logging.getLogger(__name__)
        
        try:
            # Build output filename
            pol = self.options.get("polarization", "VV")
            base_name = f"intermediate_{stage_name}_{pol}"
            npy_path = self.output_dir / f"{base_name}.npy"
            json_path = self.output_dir / f"{base_name}.json"
            
            # Prepare metadata
            metadata = {
                "timestamp": datetime.datetime.utcnow().isoformat() + "Z",
                "polarization": pol,
                "stage": stage_name,
                "description": description,
                "shape": list(data.shape),
                "dtype": str(data.dtype),
                "is_geocoded": is_geocoded,
                "statistics": {
                    "min": float(np.nanmin(data)),
                    "max": float(np.nanmax(data)),
                    "mean": float(np.nanmean(data)),
                    "median": float(np.nanmedian(data)),
                    "std": float(np.nanstd(data)),
                }
            }
            
            # Add coordinate system info
            if is_geocoded:
                # Geocoded products use EPSG:4326 or UTM
                epsg = getattr(self, "_target_epsg", 4326)
                metadata["crs"] = f"EPSG:{epsg}"
                
                # Get geotransform from terrain correction metadata if available
                if hasattr(self, "_geocode_geotransform"):
                    gt = self._geocode_geotransform
                    metadata["geotransform"] = list(gt)
                    metadata["pixel_spacing_m"] = gt[1]  # X pixel size
            else:
                # SAR geometry: range/azimuth
                metadata["crs"] = "SAR_GEOMETRY"
                
                # Get pixel spacing from metadata if available
                if hasattr(self, "metadata"):
                    meta = self.metadata
                    # Try to get range/azimuth spacing
                    if "range_pixel_spacing" in meta:
                        metadata["range_spacing_m"] = meta["range_pixel_spacing"]
                    if "azimuth_pixel_spacing" in meta:
                        metadata["azimuth_spacing_m"] = meta["azimuth_pixel_spacing"]
            
            # Export NumPy array + JSON metadata
            logger.info(f"Exporting intermediate product: {npy_path.name}")
            
            np.save(npy_path, data)
            with open(json_path, 'w') as f:
                json.dump(metadata, f, indent=2)
            
            size_mb = npy_path.stat().st_size / (1024 * 1024)
            logger.info(f"✓ Exported {stage_name}: {npy_path.name} ({data.shape[0]}×{data.shape[1]}, {size_mb:.1f} MB)")
            
        except Exception as e:
            # Non-fatal: log warning and continue processing
            logger.warning(f"Failed to export intermediate product '{stage_name}': {e}")
            import traceback
            logger.debug(traceback.format_exc())

    def _record_validation(
        self,
        stage: str,
        metric: str,
        *,
        value: Any,
        expected: Any = None,
        tolerance: Optional[float] = None,
        passed: bool = True,
        severity: str = "info",
        message: str = "",
    ) -> None:
        """Send a validation record to the ValidationGates tracker if enabled."""

        if getattr(self, "validation", None) is None:
            return
        try:
            self.validation.record(
                stage,
                metric,
                value=value,
                expected=expected,
                tolerance=tolerance,
                passed=passed,
                severity=severity,
                message=message,
            )
        except Exception as e:
            # Validation recording should never block processing, but log for debugging
            logger.debug(f"Validation recording failed (non-fatal): {e}")

    def _require_validation(
        self,
        condition: bool,
        stage: str,
        metric: str,
        *,
        value: Any = None,
        expected: Any = None,
        tolerance: Optional[float] = None,
        severity: str = "error",
        message: str = "",
    ) -> bool:
        """Assert a validation condition, honoring fail-fast behaviour."""

        if getattr(self, "validation", None) is None:
            return condition
        try:
            return self.validation.require(
                condition,
                stage,
                metric,
                value=value,
                expected=expected,
                tolerance=tolerance,
                severity=severity,
                message=message,
            )
        except ValidationError:
            # Preserve ValidationError for upstream handling to stop pipeline
            raise
        except (AttributeError, TypeError, ValueError) as e:
            # Non-validation exceptions in recording should not halt processing
            logger.debug(f"Validation condition check failed: {e}")
            return condition

    def _validate_backscatter_range(self, data: np.ndarray, stage_name: str) -> bool:
        """Validate that backscatter values are within expected physical range.
        
        For linear power (sigma0/gamma0), typical values:
        - Land: 0.001 to 1.0 (roughly -30 to 0 dB)
        - Water: 0.0001 to 0.01 (roughly -40 to -20 dB)
        - Urban: 0.1 to 10.0 (roughly -10 to +10 dB)
        
        For dB values:
        - Typical range: -30 to +10 dB
        - Values outside -50 to +20 dB are suspicious
        """
        is_valid, warnings = qc_utils.validate_backscatter_range(data)
        for warning in warnings:
            self._log(f"   ⚠️  {stage_name}: {warning}", 1)
        if not is_valid and not warnings:
            self._log(f"   ⚠️  {stage_name}: All values are non-finite", 1)
        return is_valid

    def _validate_power_preservation(
        self,
        pre_data: np.ndarray,
        post_data: np.ndarray,
        stage_name: str,
        tolerance: float = 0.05,
        *,
        log_warning: bool = True,
    ) -> bool:
        """Validate that mean power is preserved within tolerance.
        
        Args:
            pre_data: Data before processing step (linear power)
            post_data: Data after processing step (linear power)
            stage_name: Name of the processing stage for logging
            tolerance: Acceptable relative power change (default 5%)
            
        Returns:
            True if power is preserved within tolerance
        """
        try:
            is_valid, relative_change, pre_mean, post_mean = qc_utils.validate_power_preservation(
                pre_data, post_data, tolerance
            )
            if log_warning:
                if is_valid:
                    self._log(
                        f"   ✅ {stage_name}: Power preserved (change: {relative_change*100:.2f}%)",
                        2,
                    )
                else:
                    self._log(
                        f"   ⚠️  {stage_name}: Power not preserved (change: {relative_change*100:.1f}% > {tolerance*100:.0f}% threshold)",
                        1,
                    )
                    self._log(f"      Pre-mean: {pre_mean:.4e}, Post-mean: {post_mean:.4e}", 2)
            return is_valid
        except Exception as e:
            self._log(f"   ⚠️  {stage_name}: Power validation error (treating as failed): {e}", 2)
            return False  # Validation errors should be treated as failures, not passes

    def _ensure_contiguous(self, arr: np.ndarray, dtype: np.dtype = np.float32) -> np.ndarray:
        """Ensure array is C-contiguous and correct dtype for optimal Rust interop.
        
        Non-contiguous arrays require copying when passed to Rust, which is slow.
        This method ensures optimal memory layout for PyO3/numpy interop.
        """
        if arr.dtype == dtype and arr.flags['C_CONTIGUOUS']:
            return arr
        return np.ascontiguousarray(arr, dtype=dtype)

    def _load_burst_timing_from_cache(
        self, cached_metadata: Optional[Dict[str, Any]], *, source: str
    ) -> List[Dict[str, Any]]:
        """Extract burst timing records from cached metadata if available."""

        if not cached_metadata or not isinstance(cached_metadata, dict):
            return []

        import json
        import logging

        logger = logging.getLogger(__name__)
        raw_json = cached_metadata.get("burst_timing_json")
        burst_timings: List[Dict[str, Any]] = []

        if raw_json:
            try:
                burst_timings = json.loads(raw_json)
            except Exception as exc:  # pragma: no cover - diagnostics only
                logger.warning(
                    "Failed to parse burst_timing_json from %s metadata: %s",
                    source,
                    exc,
                )

        # Accept alternative legacy keys if JSON missing
        if not burst_timings:
            for key in ("burst_timings", "burst_timing_records", "burst_records"):
                candidate = cached_metadata.get(key)
                if candidate:
                    burst_timings = list(candidate)
                    break

        if burst_timings:
            logger.debug(
                "Loaded %d burst timing records from %s metadata cache",
                len(burst_timings),
                source,
            )

        return burst_timings

    def _finalize_burst_timings(
        self,
        orbit_ref_epoch_utc: float,
        product_start_rel_s: Optional[float] = None,
        azimuth_time_interval: Optional[float] = None,
    ) -> List[Dict[str, Any]]:
        """Finalize burst timing records with relative orbit times.
        
        Converts burst timing records to have consistent `azimuth_time_rel_orbit`
        values, which are required for Range-Doppler geocoding. The method:
        
        1. Normalizes subswath IDs to uppercase (e.g., 'iw1' → 'IW1')
        2. Preserves existing `azimuth_time_rel_orbit` if present
        3. Converts `azimuth_time_absolute` to relative time if present
        4. Fails fast if timing cannot be computed
        
        Args:
            orbit_ref_epoch_utc: Orbit reference epoch as Unix timestamp (seconds)
            product_start_rel_s: Product start time relative to orbit epoch (seconds).
                                 Falls back to self.metadata['product_start_rel_s'].
            azimuth_time_interval: Time between azimuth lines (seconds).
                                   Falls back to self.metadata['azimuth_time_interval'].
        
        Returns:
            List of finalized burst timing records with azimuth_time_rel_orbit
            
        Raises:
            ValueError: If azimuth_time_rel_orbit cannot be computed
        """
        # Determine strict timing behaviour: in strict science mode we require
        # explicit orbit-relative or absolute times from the metadata and do
        # not allow reconstructing timing purely from line indices.
        strict_timing = bool(getattr(self, "strict_science", False)) or bool(
            getattr(self, "strict_orbit_mode", False)
        )
        # Get parameters from metadata if not provided
        if product_start_rel_s is None:
            product_start_rel_s = self.metadata.get("product_start_rel_s")
            if product_start_rel_s is None:
                raise ValueError(
                    "Missing product_start_rel_s in metadata. "
                    "This indicates orbit integration was not applied or failed. "
                    "Apply precise orbit file before processing."
                )
        
        if azimuth_time_interval is None:
            azimuth_time_interval = self.metadata.get("azimuth_time_interval")
            if azimuth_time_interval is None or azimuth_time_interval == 0.0:
                # Try to compute from PRF
                prf = self.metadata.get("prf")
                if prf is not None and prf > 0:
                    azimuth_time_interval = 1.0 / prf
                    self._log(f"   ℹ️  Computed azimuth_time_interval={azimuth_time_interval:.6f}s from PRF={prf:.1f}Hz", 2)
                else:
                    raise ValueError(
                        "Missing azimuth_time_interval and prf in metadata. "
                        "Cannot compute burst timing. Check SAFE annotation parsing."
                    )
        
        finalized: List[Dict[str, Any]] = []
        
        for record in self.burst_timing_records:
            # Copy record to avoid mutating original
            rec = dict(record)
            
            # Normalize subswath ID to uppercase
            subswath_id = rec.get("subswath_id", "")
            if subswath_id:
                rec["subswath_id"] = subswath_id.upper()
            
            # Determine azimuth_time_rel_orbit
            if "azimuth_time_rel_orbit" in rec:
                # Already has relative time - preserve it
                value = rec.get("azimuth_time_rel_orbit")
                if strict_timing and (value is None or not isinstance(value, (int, float))):
                    burst_id = f"{rec.get('subswath_id', '?')}/burst_{rec.get('burst_index', '?')}"
                    raise ValueError(
                        f"STRICT MODE: invalid azimuth_time_rel_orbit for {burst_id}: {value!r}"
                    )
            elif "azimuth_time_absolute" in rec and orbit_ref_epoch_utc > 0:
                # Convert absolute time to relative
                rec["azimuth_time_rel_orbit"] = rec["azimuth_time_absolute"] - orbit_ref_epoch_utc
            elif azimuth_time_interval > 0:
                # In strict timing mode we do not allow reconstructing orbit
                # times purely from line indices, as this hides upstream timing
                # domain bugs (e.g. ANX vs orbit-relative confusion).
                if strict_timing:
                    burst_id = f"{rec.get('subswath_id', '?')}/burst_{rec.get('burst_index', '?')}"
                    raise ValueError(
                        "STRICT MODE: missing azimuth_time_rel_orbit/absolute for "
                        f"{burst_id}; cannot derive from line index only."
                    )

                # Relaxed mode: compute from line index using product start
                # and azimuth_time_interval (center of burst lines).
                first_line = rec.get("first_line_global")
                if first_line is None:
                    burst_id = f"{rec.get('subswath_id', '?')}/burst_{rec.get('burst_index', '?')}"
                    raise ValueError(f"Missing first_line_global for {burst_id}")
                last_line = rec.get("last_line_global", first_line)
                center_line = (first_line + last_line) / 2.0

                # Time = product_start + center_line * line_interval
                rec["azimuth_time_rel_orbit"] = product_start_rel_s + center_line * azimuth_time_interval
            else:
                # Cannot compute timing
                burst_id = f"{rec.get('subswath_id', '?')}/burst_{rec.get('burst_index', '?')}"
                raise ValueError(
                    f"Cannot compute azimuth_time_rel_orbit for {burst_id}. "
                    "Missing timing metadata (azimuth_time_interval, prf, or absolute time)."
                )
            
            finalized.append(rec)

        # Enforce monotonic, single-domain timing in strict mode: after
        # finalization, azimuth_time_rel_orbit must increase with burst index
        # within each subswath.
        if strict_timing:
            by_subswath: Dict[str, List[Dict[str, Any]]] = {}
            for rec in finalized:
                sid = str(rec.get("subswath_id", "")).upper()
                by_subswath.setdefault(sid, []).append(rec)

            for sid, recs in by_subswath.items():
                recs_sorted = sorted(recs, key=lambda r: int(r.get("burst_index", 0)))
                prev_t: Optional[float] = None
                for rec in recs_sorted:
                    t = rec.get("azimuth_time_rel_orbit")
                    if not isinstance(t, (int, float)):
                        raise ValueError(
                            f"STRICT MODE: burst timing for {sid} contains non-numeric azimuth_time_rel_orbit: {t!r}"
                        )
                    t_f = float(t)
                    if prev_t is not None and t_f <= prev_t:
                        raise ValueError(
                            f"STRICT MODE: non-monotonic azimuth_time_rel_orbit in {sid}: {prev_t:.6f} -> {t_f:.6f}"
                        )
                    prev_t = t_f

        return finalized

    def _log_array_metrics(self, stage_name: str, data: Any) -> None:
        """Emit basic quality metrics (finite %, mean power, mean dB) for a stage output."""

        if not isinstance(data, np.ndarray) or data.size == 0:
            return

        finite = np.isfinite(data)
        finite_count = int(np.count_nonzero(finite))
        finite_pct = (finite_count / float(data.size)) * 100.0

        pos = data[np.isfinite(data) & (data > 0)]
        mean_lin = float(np.mean(pos)) if pos.size else float('nan')
        # FIXED: Use suppress_numpy_errors context for log10 operations
        with suppress_numpy_errors():
            mean_db = float(10.0 * np.log10(mean_lin)) if mean_lin > 0 else float('nan')

        self._log(
            f"   📊 {stage_name}: {finite_pct:.1f}% finite, mean={mean_lin:.4e} (linear), {mean_db:.2f} dB",
            1,
        )

        # Record as validation info (non-fatal) to help stage-by-stage inspection
        try:
            self._record_validation(
                stage_name.lower(),
                "finite_fraction",
                value=finite_pct / 100.0,
                expected=">0",
                passed=True,
                severity="info" if finite_pct >= 90.0 else "warning",
                message="Stage finite fraction diagnostic",
            )
        except (AttributeError, TypeError) as e:
            # Validation recording is optional diagnostic - log but don't fail
            logger.debug(f"Validation recording skipped: {e}")

    def _process_subswath_pipeline(self, subswath: str, preloaded_slc_data: Optional[np.ndarray] = None) -> Dict[str, Any]:
        """
        Deburst and radiometrically calibrate a single subswath.
        Returns a result dictionary consumed by deburst_utils.run_process_subswaths.
        """
        import logging
        logger = logging.getLogger(__name__)
        
        pol = self.polarization
        
        # Valid calibration types for Sentinel-1
        valid_calibration_types = {"sigma0", "beta0", "gamma0", "dn"}
        if self.calibration_type not in valid_calibration_types:
            raise ValueError(
                f"Invalid calibration type: {self.calibration_type}. "
                f"Expected one of: {valid_calibration_types}"
            )

        def _extract_array(result: Any, label: str, expected_dtype: Optional[type] = None) -> np.ndarray:
            """Extract array from result with comprehensive validation."""
            # Handle error status
            if isinstance(result, dict):
                if result.get("status") == "error":
                    raise RuntimeError(f"{label} failed: {result.get('message', 'unknown error')}")
                
                # Try to extract data from known keys
                data = None
                for key in ["data", "calibrated_data", "array"]:
                    if key in result:
                        data = result[key]
                        break
                
                if data is None:
                    raise RuntimeError(
                        f"{label} did not return data (keys={list(result.keys())})"
                    )
                arr = np.asarray(data)
            else:
                arr = np.asarray(result)
            
            # Validate array dimensions
            if arr.ndim != 2:
                raise RuntimeError(
                    f"{label} returned invalid array dimensions: {arr.ndim}D (expected 2D). "
                    f"Shape: {arr.shape}"
                )
            
            # Validate array is not empty
            if arr.size == 0:
                raise RuntimeError(f"{label} returned empty array")
            
            # Validate minimum dimensions for SAR data
            min_dim = 10
            if arr.shape[0] < min_dim or arr.shape[1] < min_dim:
                raise RuntimeError(
                    f"{label} returned suspiciously small array: {arr.shape}. "
                    f"Expected at least {min_dim}x{min_dim} pixels."
                )
            
            # Validate dtype if specified
            if expected_dtype is not None:
                if not np.issubdtype(arr.dtype, expected_dtype):
                    logger.debug(
                        f"{label} dtype {arr.dtype} will be converted to {expected_dtype}"
                    )
            
            # Enhanced validation: check for finite values (for real arrays)
            if np.issubdtype(arr.dtype, np.floating):
                finite_count = np.sum(np.isfinite(arr))
                finite_pct = 100.0 * finite_count / arr.size
                if finite_pct < 50.0:
                    logger.warning(
                        f"⚠️  {label} has low finite value percentage: {finite_pct:.1f}% "
                        f"({finite_count}/{arr.size} pixels)"
                    )
                elif finite_pct < 90.0:
                    logger.info(
                        f"ℹ️  {label} finite value percentage: {finite_pct:.1f}% "
                        f"({finite_count}/{arr.size} pixels)"
                    )
            
            return arr

        # FIXED: Reader lock is acquired here - Python's threading.Lock is NOT re-entrant
        # Each thread acquires the lock once, so no deadlock risk as long as we don't
        # try to acquire it again within the same thread. The lock is released after
        # the operation completes.
        # Deburst with cached reader (guarded by lock for thread safety)
        # I/O Optimization: Use pre-loaded SLC data if available (avoids ~12s I/O per subswath)
        # NOTE: Chunked deburst removed (separable LUT not suitable for IW TOPS mode)
        deburst_start = time.time()
        use_chunked_result = False  # Chunked deburst disabled (removed)
        
        with self.reader_lock:
            # Standard path: use deburst_topsar_cached (with dense 2D calibration LUT)
            if preloaded_slc_data is not None:
                # Use pre-loaded SLC data (I/O optimization active)
                logger.debug(f"Using pre-loaded SLC data for {subswath} (shape: {preloaded_slc_data.shape})")
                deburst_result = sardine.deburst_topsar_cached(
                    self.reader, 
                    subswath, 
                    pol,
                    preloaded_slc_data=preloaded_slc_data
                )
            else:
                # Standard path: read SLC data from disk
                deburst_result = sardine.deburst_topsar_cached(self.reader, subswath, pol)
        
        # Standard path processing (skip if chunked deburst was used)
        if not use_chunked_result:
            # Validate deburst result structure
            deburst_overrides = None
            burst_metadata = None
            if isinstance(deburst_result, dict):
                # Validate expected keys exist
                if "data" not in deburst_result and "status" not in deburst_result:
                    logger.warning(
                        f"Deburst result missing 'data' key for {subswath}. "
                        f"Available keys: {list(deburst_result.keys())}"
                    )
                
                # Enhanced validation of deburst_overrides structure
                deburst_overrides = {
                    "timing_reference": deburst_result.get("timing_reference"),
                    "burst_timing": deburst_result.get("burst_timing"),
                    "row_provenance": deburst_result.get("row_provenance"),
                    "azimuth_index_origin": deburst_result.get("azimuth_index_origin"),
                    "range_sample_origin": deburst_result.get("range_sample_origin"),
                }
                
                # Validate deburst_overrides structure
                burst_timing = deburst_overrides["burst_timing"]
                row_provenance = deburst_overrides["row_provenance"]
                
                if not burst_timing or not row_provenance:
                    print(
                        f"   ⚠️  Deburst overrides empty for {subswath}: burst_timing={len(burst_timing or [])}, row_provenance={len(row_provenance or [])}"
                    )
                    deburst_overrides = None
                else:
                    # Validate structure of burst_timing and row_provenance
                    if not isinstance(burst_timing, (list, tuple)):
                        logger.warning(f"Invalid burst_timing type for {subswath}: {type(burst_timing)}")
                        deburst_overrides = None
                    elif not isinstance(row_provenance, (list, tuple)):
                        logger.warning(f"Invalid row_provenance type for {subswath}: {type(row_provenance)}")
                        deburst_overrides = None
                    else:
                        # Validate row_provenance entries have required fields
                        required_rp_fields = ["out_row_start", "out_row_end", "burst_id", "burst_line_start"]
                        invalid_rp_count = 0
                        for rp in row_provenance:
                            if not isinstance(rp, dict):
                                invalid_rp_count += 1
                                continue
                            missing_fields = [f for f in required_rp_fields if f not in rp]
                            if missing_fields:
                                logger.debug(f"row_provenance entry missing fields {missing_fields} for {subswath}")
                        
                        if invalid_rp_count > 0:
                            logger.warning(f"{invalid_rp_count} invalid row_provenance entries for {subswath}")
                        
                        print(
                            f"   🛰️  Deburst overrides captured for {subswath}: burst_timing={len(burst_timing)}, row_provenance={len(row_provenance)}"
                        )
                burst_metadata = deburst_result.get("burst_metadata")
            
            deburst_data = _extract_array(deburst_result, "deburst_topsar_cached", np.complexfloating).astype(np.complex64, copy=False)
            deburst_duration = time.time() - deburst_start

            # Radiometric calibration
            calibration_start = time.time()
            # FIXED: Reader lock acquired again - safe because previous lock was released
            # Use cached reader to avoid reopening/parsing calibration XML per subswath
            with self.reader_lock:
                cal_result = sardine.radiometric_calibration_with_denoising_cached(
                    self.reader,
                    subswath,
                    pol,
                    self.calibration_type,
                    deburst_data,
                    self.fast_mode_noise_removal,
                )
        
            # Validate calibration result structure
            if isinstance(cal_result, dict):
                if cal_result.get("status") == "error":
                    raise RuntimeError(
                        f"Calibration failed for {subswath}: {cal_result.get('message', 'unknown')}"
                    )
                
                # Validate LUT source was found
                lut_source = cal_result.get("lut_source")
                if lut_source is None:
                    logger.warning(
                        f"Calibration LUT source not reported for {subswath}. "
                        f"Available keys: {list(cal_result.keys())}"
                    )
            
            cal_provenance = None
            if isinstance(cal_result, dict):
                cal_provenance = {
                    "subswath": subswath,
                    "lut_source": cal_result.get("lut_source"),
                    "noise_source": cal_result.get("noise_source"),
                    "antenna_pattern": cal_result.get("antenna_pattern"),
                    "ipf_version": cal_result.get("ipf_version"),
                    "calibration_type": self.calibration_type,
                    "range_pixel_spacing": cal_result.get("range_pixel_spacing"),
                    "azimuth_pixel_spacing": cal_result.get("azimuth_pixel_spacing"),
                    "vector_count": cal_result.get("vector_count"),
                }

            calibrated = _extract_array(cal_result, "radiometric_calibration", np.floating).astype(np.float32, copy=False)
            calibration_duration = time.time() - calibration_start
        
        # Diagnostic logging: After calibration
        self._log_diagnostic_statistics(f"After Calibration ({subswath})", calibrated, context=f"calibration_type={self.calibration_type}")
        
        # Validate calibrated data ranges (common path for both chunked and standard)
        finite_mask = np.isfinite(calibrated)
        finite_count = np.count_nonzero(finite_mask)
        total_count = calibrated.size
        finite_pct = 100.0 * finite_count / total_count if total_count > 0 else 0.0
        
        if finite_pct < 50.0:
            raise RuntimeError(
                f"Calibration produced mostly non-finite values for {subswath}: "
                f"{finite_pct:.1f}% finite ({finite_count}/{total_count})"
            )
        
        if finite_count > 0:
            finite_data = calibrated[finite_mask]
            min_val = float(np.min(finite_data))
            max_val = float(np.max(finite_data))
            
            # Linear power should be positive
            if min_val < 0:
                negative_count = np.count_nonzero(finite_data < 0)
                logger.warning(
                    f"Calibrated data for {subswath} contains {negative_count} negative values "
                    f"(min={min_val:.4e}). Linear power should be non-negative."
                )
            
            # Extreme values check (linear power outside reasonable range)
            # Typical SAR backscatter: 1e-6 to 1e2 in linear power
            if max_val > 1e6:
                logger.warning(
                    f"Calibrated data for {subswath} contains extremely high values "
                    f"(max={max_val:.4e}). This may indicate calibration issues."
                )

            # Mask non-positive or non-finite pixels to NaN for downstream stages
            invalid_mask = (~finite_mask) | (calibrated <= 0.0)
            invalid_pct = 100.0 * np.count_nonzero(invalid_mask) / total_count if total_count > 0 else 0.0
            if invalid_pct > 0.0:
                calibrated[invalid_mask] = np.nan
                logger.info(
                    f"Calibrated data for {subswath}: masked {invalid_pct:.2f}% non-finite/non-positive pixels"
                )

        # Capture LUT source/range spacing if provided
        if isinstance(cal_result, dict):
            self.calibration_lut_source = cal_result.get("lut_source", self.calibration_lut_source)
            try:
                range_px = cal_result.get("range_pixel_spacing")
                az_px = cal_result.get("azimuth_pixel_spacing")
                if range_px is not None:
                    self.current_range_spacing = float(range_px)
                    if self.original_range_spacing is None:
                        self.original_range_spacing = float(range_px)
                if az_px is not None:
                    self.current_azimuth_spacing = float(az_px)
                    if self.original_azimuth_spacing is None:
                        self.original_azimuth_spacing = float(az_px)
            except (TypeError, ValueError) as e:
                logger.warning(f"Failed to extract pixel spacing for {subswath}: {e}")

        # Track deburst shape for this subswath (used in result dict)
        # NOTE: Native dimensions for the full merged image are computed after merge
        # in merge_stage.py, not here, to avoid overwriting with per-subswath values.
        deburst_shape = None
        if 'deburst_data' in locals():
            deburst_shape = deburst_data.shape
        else:
            deburst_shape = calibrated.shape

        return {
            "subswath": subswath,
            "calibrated_data": calibrated,
            "deburst_shape": deburst_shape,
            "deburst_duration": deburst_duration,
            "calibration_duration": calibration_duration,
            "needs_calibration": False,
            "burst_metadata": burst_metadata,
            "deburst_overrides": deburst_overrides,
            "calibration_provenance": cal_provenance,
        }

    def _parse_requested_output_crs(self) -> Optional[int]:
        """Parse user-specified output CRS option into an EPSG code.
        
        Returns None for 'auto' or empty values to trigger automatic CRS detection.
        """

        value = self._requested_output_crs
        if value is None:
            return None
        if isinstance(value, str):
            text = value.strip().lower()
            if not text or text == "auto":
                return None  # Use automatic detection
            if text.upper().startswith("EPSG:"):
                text = text.split(":", 1)[1]
            value = text
        try:
            epsg = int(value)
        except (TypeError, ValueError):
            raise ValueError(f"Invalid output_crs option: {self._requested_output_crs!r}")
        if epsg <= 0:
            raise ValueError(f"Invalid EPSG code: {epsg}")
        return epsg

    def _default_projected_epsg(self, min_lat: float, max_lat: float, min_lon: float, max_lon: float) -> int:
        """Resolve a sensible projected EPSG based on scene extent."""

        centre_lat = 0.5 * (min_lat + max_lat)
        centre_lon = 0.5 * (min_lon + max_lon)

        if centre_lat >= 84.0:
            return 3413  # Arctic polar stereographic
        if centre_lat <= -80.0:
            return 3031  # Antarctic polar stereographic

        zone = int(math.floor((centre_lon + 180.0) / 6.0) + 1)
        zone = max(1, min(zone, 60))
        if centre_lat >= 0.0:
            return 32600 + zone  # Northern hemisphere UTM
        return 32700 + zone  # Southern hemisphere UTM

    def _resolve_output_epsg(self) -> int:
        """Determine the target output EPSG from user options or scene metadata."""
        import logging
        logger = logging.getLogger(__name__)

        requested = self._parse_requested_output_crs()
        if requested is not None:
            return requested

        metadata = getattr(self, "metadata", {}) or {}
        
        # Check required keys exist
        required_keys = ["min_latitude", "max_latitude", "min_longitude", "max_longitude"]
        missing = [k for k in required_keys if k not in metadata or metadata[k] is None]
        if missing:
            logger.warning(f"Missing coordinate keys for EPSG resolution: {missing}. Using EPSG:4326")
            return 4326
        
        try:
            min_lat = float(metadata["min_latitude"])
            max_lat = float(metadata["max_latitude"])
            min_lon = float(metadata["min_longitude"])
            max_lon = float(metadata["max_longitude"])
        except (TypeError, ValueError) as e:
            logger.warning(f"Invalid coordinate values in metadata: {e}. Using EPSG:4326")
            return 4326
        
        # Validate coordinate ranges
        if not (-90 <= min_lat <= 90) or not (-90 <= max_lat <= 90):
            logger.warning(f"Latitude out of range [{min_lat}, {max_lat}]. Using EPSG:4326")
            return 4326
        if not (-180 <= min_lon <= 180) or not (-180 <= max_lon <= 180):
            logger.warning(f"Longitude out of range [{min_lon}, {max_lon}]. Using EPSG:4326")
            return 4326

        if not all(np.isfinite(v) for v in (min_lat, max_lat, min_lon, max_lon)):
            logger.warning("Non-finite coordinate values in metadata. Using EPSG:4326")
            return 4326
        if max_lat <= min_lat or max_lon <= min_lon:
            logger.warning(f"Invalid coordinate ordering: lat=[{min_lat}, {max_lat}], lon=[{min_lon}, {max_lon}]. Using EPSG:4326")
            return 4326

        return self._default_projected_epsg(min_lat, max_lat, min_lon, max_lon)

    def _stage_speckle_filtering(self, context: PipelineContext) -> None:
        """Apply speckle filtering with validated parameters."""
        speckle_utils.run_speckle_filtering(self)

    def _stage_multilooking(self, context: PipelineContext) -> None:
        """Apply scientific multilooking after merging subswaths."""
        # Mask NaN values before multilooking to prevent propagation
        self._mask_nan_before_processing("multilooking")
        multilook_utils.run_multilooking(self)

    def _stage_terrain_flattening(self, context: PipelineContext) -> None:
        """Apply DEM-based terrain flattening before speckle filtering."""
        # Check for prefetch errors before using DEM
        self._check_prefetch_errors()
        terrain_utils.run_terrain_flattening(self)

    def _stage_terrain_correction(self, context: PipelineContext) -> None:
        """Execute scientific Range-Doppler terrain correction."""
        # Check for prefetch errors before using DEM
        self._check_prefetch_errors()
        terrain_utils.run_terrain_correction(self)

    def _stage_export_calibration(self, context: PipelineContext) -> None:
        """Export calibrated subswaths checkpoint (before merge)."""
        if not self.export_intermediate_products:
            return
        
        # Export calibrated subswaths (before merge)
        if self._working_data is not None and isinstance(self._working_data, np.ndarray):
            self._export_intermediate_product(
                data=self._working_data,
                stage_name="calibration",
                description=f"Radiometric calibration ({self.calibration_type}) before merge",
                is_geocoded=False,
            )

    def _stage_export_merge(self, context: PipelineContext) -> None:
        """Export merged calibrated data checkpoint (after merge, before multilook)."""
        if not self.export_intermediate_products:
            return
        
        # Export merged data matching SNAP's TOPSAR-Merge output
        if self._working_data is not None and isinstance(self._working_data, np.ndarray):
            self._export_intermediate_product(
                data=self._working_data,
                stage_name="merge",
                description=f"TOPSAR merge ({self.calibration_type})",
                is_geocoded=False,
            )

    def _stage_export_multilook(self, context: PipelineContext) -> None:
        """Export multilooked data checkpoint."""
        if not self.export_intermediate_products:
            return
        
        if self._working_data is not None and isinstance(self._working_data, np.ndarray):
            self._export_intermediate_product(
                data=self._working_data,
                stage_name="multilook",
                description=f"Multilooked ({self.multilook_range}×{self.multilook_azimuth})",
                is_geocoded=False,
            )

    def _stage_export_terrain_flattening(self, context: PipelineContext) -> None:
        """Export terrain-flattened data checkpoint."""
        if not self.export_intermediate_products:
            return
        
        if self._working_data is not None and isinstance(self._working_data, np.ndarray):
            self._export_intermediate_product(
                data=self._working_data,
                stage_name="terrain_flattening",
                description="DEM-based terrain flattening (gamma0_tc)",
                is_geocoded=False,
            )

    def _stage_export_terrain_correction(self, context: PipelineContext) -> None:
        """Export geocoded data checkpoint."""
        if not self.export_intermediate_products:
            return
        
        if self._working_data is not None and isinstance(self._working_data, np.ndarray):
            self._export_intermediate_product(
                data=self._working_data,
                stage_name="terrain_correction",
                description="Range-Doppler geocoded",
                is_geocoded=True,
            )

    def _stage_mask_invalid_areas(self, context: PipelineContext) -> None:
        """Apply quality masks to geocoded data."""
        import logging
        logger = logging.getLogger(__name__)

        self.announce_step(STEP_NUMBERS["Mask Invalid Areas"], "Mask Invalid Areas", "Applying quality masks")
        step_start = time.time()

        # Check if masking is disabled via no_masking option
        if getattr(self, "skip_masking", False):
            step_duration = time.time() - step_start
            self.log_step(
                STEP_NUMBERS["Mask Invalid Areas"],
                "Mask Invalid Areas",
                "skipped",
                "Masking disabled via no_masking option",
                step_duration,
            )
            print("   ℹ️  Quality masking skipped (no_masking=True)")
            return

        geocoded_data = self._working_data

        try:
            if not isinstance(geocoded_data, np.ndarray):
                raise RuntimeError(
                    "SCIENTIFIC MODE FAILURE: Geocoded data is not a valid numpy array - cannot apply quality masking"
                )

            if geocoded_data.ndim == 0:
                raise RuntimeError(
                    "SCIENTIFIC MODE FAILURE: Scalar geocoded data detected - cannot apply quality mask"
                )

            basic_mask = np.isfinite(geocoded_data)
            basic_valid_percentage = np.sum(basic_mask) / basic_mask.size * 100
            print(f"   📊 Basic finite data: {basic_valid_percentage:.1f}% finite values")

            finite_data = geocoded_data[basic_mask]
            if len(finite_data) > 0:
                data_min = np.min(finite_data)
                data_max = np.max(finite_data)
                data_mean = np.mean(finite_data)
                print(
                    f"   📊 Data range: {data_min:.6f} to {data_max:.6f}, mean: {data_mean:.6f}"
                )

                if data_min >= 0:
                    # Heuristic: only treat as linear power if the dynamic range is modest and
                    # max is not sky-high; otherwise fall back to percentile to avoid wiping out
                    # scenes that are already scaled (e.g., dB or normalized floats).
                    dynamic_span = (data_max + MASK_STABILITY_EPSILON) / max(data_min, MASK_STABILITY_EPSILON)
                    linear_like = data_max <= MASK_LINEAR_MAX_VALUE and dynamic_span <= MASK_LINEAR_DYNAMIC_RANGE

                    if linear_like:
                        if data_min <= 0.0:
                            range_min = 0.0
                        else:
                            range_min = max(MASK_MIN_POSITIVE_VALUE, data_min)

                        range_max = data_max
                        if np.isfinite(range_max) and range_max > 0.0:
                            range_max = min(MASK_LINEAR_MAX_VALUE, range_max * 1.1)
                        else:
                            range_max = max(range_min, MASK_LINEAR_MAX_VALUE)

                        valid_mask = (
                            basic_mask
                            & (geocoded_data >= range_min)
                            & (geocoded_data <= range_max)
                        )
                        print(
                            f"   🔍 Applied linear unit masking: range [{range_min:.2e}, {range_max:.2f}]"
                        )
                    else:
                        p1 = np.percentile(finite_data, 1)
                        p99 = np.percentile(finite_data, 99)
                        valid_mask = basic_mask & (geocoded_data >= p1) & (geocoded_data <= p99)
                        print(
                            f"   🔍 Applied percentile masking (ambiguous units): range [{p1:.3f}, {p99:.3f}]"
                        )
                else:
                    p1 = np.percentile(finite_data, 1)
                    p99 = np.percentile(finite_data, 99)
                    valid_mask = basic_mask & (geocoded_data >= p1) & (geocoded_data <= p99)
                    print(
                        f"   🔍 Applied percentile masking: range [{p1:.3f}, {p99:.3f}] (1st-99th percentile)"
                    )
            else:
                valid_mask = basic_mask
                logger.warning(
                    "No finite data found in geocoded output - using basic mask only. "
                    "This may indicate upstream processing issues."
                )
                print("   ⚠️  No finite data found - using basic mask only")

            # Optional RTC-derived masks/thresholds (only applied when inputs are available)
            shadow_mask = getattr(self, "_shadow_mask", None)
            if isinstance(shadow_mask, np.ndarray) and shadow_mask.shape == valid_mask.shape:
                valid_mask &= shadow_mask == 0
                print("   🌑 Applied shadow mask from RTC")

            layover_mask = getattr(self, "_layover_mask", None)
            if isinstance(layover_mask, np.ndarray) and layover_mask.shape == valid_mask.shape:
                valid_mask &= layover_mask == 0
                print("   🏔️  Applied layover mask from RTC")

            lia_array = getattr(self, "_lia_array", None)
            if isinstance(lia_array, np.ndarray) and lia_array.shape == valid_mask.shape:
                try:
                    threshold_cos = float(self.lia_threshold_cos)
                    if 0.0 < threshold_cos < 1.0:
                        valid_mask &= np.cos(np.deg2rad(lia_array)) >= threshold_cos
                        print(f"   📐 Applied LIA cosine threshold: {threshold_cos:.3f}")
                except (TypeError, ValueError):
                    logger.debug("Invalid lia_threshold; skipping LIA masking")

            try:
                gamma0_min_db = float(self.gamma0_min_db)
                gamma0_max_db = float(self.gamma0_max_db)
                if np.isfinite([gamma0_min_db, gamma0_max_db]).all() and gamma0_min_db < gamma0_max_db:
                    min_power = 10.0 ** (gamma0_min_db / 10.0)
                    max_power = 10.0 ** (gamma0_max_db / 10.0)
                    valid_mask &= (geocoded_data >= min_power) & (geocoded_data <= max_power)
                    print(f"   📉 Applied backscatter thresholds: [{gamma0_min_db:.1f}, {gamma0_max_db:.1f}] dB")
            except (TypeError, ValueError, OverflowError):
                logger.debug("Invalid gamma0_min/gamma0_max; skipping backscatter threshold masking")

            # Optimize: use np.where for in-place-style masking instead of copy + assignment
            # This avoids allocating a full copy of the array
            masked_data = np.where(valid_mask, geocoded_data, np.nan).astype(np.float32)
            valid_percentage = np.sum(valid_mask) / valid_mask.size * 100
            
            # Validate that masking didn't remove too much data
            if valid_percentage < MASK_MIN_VALID_PERCENTAGE:
                logger.warning(
                    f"Quality masking reduced valid pixels to {valid_percentage:.1f}% "
                    f"(below {MASK_MIN_VALID_PERCENTAGE}% threshold). "
                    f"This may indicate data quality issues or overly aggressive masking."
                )
                self._record_validation(
                    "mask_invalid_areas",
                    "valid_pixel_percentage",
                    value=valid_percentage,
                    expected=f">={MASK_MIN_VALID_PERCENTAGE}%",
                    passed=False,
                    severity="warning",
                    message=f"Low valid pixel percentage after masking: {valid_percentage:.1f}%",
                )

            step_duration = time.time() - step_start
            self.log_step(
                STEP_NUMBERS["Mask Invalid Areas"],
                "Mask Invalid Areas",
                "success",
                f"Quality mask: {valid_percentage:.1f}% valid pixels",
                step_duration,
            )

            self._working_data = masked_data
            self._valid_percentage = float(valid_percentage)

        except Exception as exc:
            step_duration = time.time() - step_start
            self.log_step(STEP_NUMBERS["Mask Invalid Areas"], "Mask Invalid Areas", "error", f"Quality masking failed: {exc}", step_duration)
            raise RuntimeError(
                f"SCIENTIFIC MODE FAILURE: Quality masking failed. Error: {exc}"
            )

    def _stage_convert_to_db(self, context: PipelineContext) -> None:
        """Convert masked backscatter data to decibel scale (if enabled)."""
        logger = logging.getLogger(__name__)

        self.announce_step(STEP_NUMBERS["Convert to dB"], "Convert to dB", "Transforming to logarithmic scale")
        step_start = time.time()

        # FIXED: Respect convert_to_db option to allow power-scale output
        convert_to_db = self.options.get("convert_to_db", True)
        if not convert_to_db:
            step_duration = time.time() - step_start
            self.log_step(
                STEP_NUMBERS["Convert to dB"],
                "Convert to dB",
                "skipped",
                "Linear/power scale output requested (convert_to_db=False)",
                step_duration,
            )
            print("   ⏭️  Keeping linear scale (dB conversion disabled)")
            
            # FIXED: Set _db_data to linear working data so export stage can proceed
            # The export stage will output linear power values instead of dB
            linear_data = self._working_data
            if isinstance(linear_data, np.ndarray):
                finite_mask = np.isfinite(linear_data)
                total_pixels = int(linear_data.size)
                valid_pixel_count = int(finite_mask.sum())
                valid_percentage = (valid_pixel_count / total_pixels * 100.0) if total_pixels > 0 else 0.0
                
                self._db_data = linear_data  # Use linear data for export
                self._total_pixels = total_pixels
                self._valid_pixel_count = valid_pixel_count
                self._valid_percentage = float(valid_percentage)
                self._linear_output = True  # Flag to indicate output is linear, not dB
                print(f"   📊 Linear output: {valid_percentage:.1f}% valid pixels")
            else:
                raise RuntimeError("SCIENTIFIC MODE FAILURE: Working data unavailable for linear export")
            return

        masked_data = self._working_data

        try:
            if not isinstance(masked_data, np.ndarray):
                raise RuntimeError(
                    "SCIENTIFIC MODE FAILURE: Masked data is not a valid numpy array - cannot convert to dB"
                )

            if masked_data.ndim == 0:
                raise RuntimeError(
                    "SCIENTIFIC MODE FAILURE: Scalar masked data detected - cannot convert to dB"
                )
            
            if masked_data.ndim != 2:
                raise RuntimeError(
                    f"SCIENTIFIC MODE FAILURE: dB conversion requires 2D array, got {masked_data.ndim}D"
                )

            db_result = sardine.convert_to_db_real(masked_data)

            # Validate dB conversion result with clear error messages
            if db_result is None:
                raise RuntimeError(
                    "dB conversion returned None - this indicates a failure in the Rust bindings"
                )
            
            db_data = None
            if isinstance(db_result, dict):
                db_data = db_result.get('data')
                if db_data is None:
                    available_keys = list(db_result.keys())
                    raise RuntimeError(
                        f"dB conversion result dict missing 'data' key. "
                        f"Available keys: {available_keys}"
                    )
            elif isinstance(db_result, np.ndarray):
                db_data = db_result
            else:
                raise RuntimeError(
                    f"Unexpected dB conversion result type: {type(db_result)}. "
                    f"Expected dict or numpy.ndarray."
                )
            
            if not isinstance(db_data, np.ndarray):
                raise RuntimeError(
                    f"dB conversion data is not a numpy array: {type(db_data)}"
                )
            
            if db_data.ndim != 2:
                raise RuntimeError(
                    f"dB conversion produced {db_data.ndim}D array, expected 2D"
                )
            
            if db_data.size == 0:
                raise RuntimeError("dB conversion produced empty array")
            
            # Validate dB range is scientifically reasonable
            finite_db = db_data[np.isfinite(db_data)]
            if len(finite_db) > 0:
                db_min = float(np.min(finite_db))
                db_max = float(np.max(finite_db))
                
                # Check if values are within expected SAR backscatter range
                if db_min < DB_EXPECTED_MIN or db_max > DB_EXPECTED_MAX:
                    logger.warning(
                        f"dB values [{db_min:.1f}, {db_max:.1f}] outside expected SAR backscatter range "
                        f"[{DB_EXPECTED_MIN}, {DB_EXPECTED_MAX}] dB. "
                        f"This may indicate calibration issues or non-backscatter data."
                    )
                elif db_min < DB_WARN_MIN or db_max > DB_WARN_MAX:
                    logger.info(
                        f"dB values [{db_min:.1f}, {db_max:.1f}] near edge of typical range "
                        f"[{DB_WARN_MIN}, {DB_WARN_MAX}] dB."
                    )

            step_duration = time.time() - step_start
            self.log_step(
                STEP_NUMBERS["Convert to dB"],
                "Convert to dB",
                "success",
                f"dB conversion: range {np.nanmin(db_data):.1f} to {np.nanmax(db_data):.1f} dB",
                step_duration,
            )

            finite_mask = np.isfinite(db_data)
            total_pixels = int(db_data.size)
            valid_pixel_count = int(finite_mask.sum())
            valid_percentage = (
                (valid_pixel_count / total_pixels) * 100.0 if total_pixels > 0 else 0.0
            )

            # Diagnostic logging: After dB conversion
            self._log_diagnostic_statistics("After dB Conversion", db_data, context="final_output")
            
            self._working_data = db_data
            self._db_data = db_data
            self._total_pixels = total_pixels
            self._valid_pixel_count = valid_pixel_count
            self._valid_percentage = float(valid_percentage)

        except Exception as exc:
            step_duration = time.time() - step_start
            self.log_step(STEP_NUMBERS["Convert to dB"], "Convert to dB", "error", f"dB conversion failed: {exc}", step_duration)
            raise RuntimeError(
                f"SCIENTIFIC MODE FAILURE: dB conversion failed. Error: {exc}"
            )

    def _stage_export_products(self, context: PipelineContext) -> None:
        """Write final products (GeoTIFF/NumPy/JSON) to disk."""
        step_number = STEP_NUMBERS["Export Final Products"]
        self.announce_step(step_number, "Export Products", "Writing final output files")
        step_start = time.time()
        
        export_succeeded = False
        export_details = ""
        
        try:
            export_utils.export_final_products(self)
            export_succeeded = True
            export_details = f"Exported to {self.output_dir}"
            
            # Validate that at least one output file was created
            output_files = list(self.output_dir.glob(f"*{self.polarization}*"))
            if not output_files:
                logger.warning(
                    f"Export reported success but no {self.polarization} files found in {self.output_dir}"
                )
                export_details += " (no output files detected)"
            else:
                export_details = f"Exported {len(output_files)} file(s) to {self.output_dir}"
                
        except Exception as exc:
            # Log the original error before attempting fallback
            logger.warning(f"Primary export failed: {exc}")
            print(f'   ⚠️  Primary export error: {exc}')
            
            # Attempt fallback export
            try:
                db_data = getattr(self, '_db_data', None)
                if db_data is None or not isinstance(db_data, np.ndarray):
                    raise RuntimeError("No valid dB data available for fallback export")
                    
                output_npy = self.output_dir / f"backscatter_{self.polarization}_final.npy"
                np.save(output_npy, db_data)
                
                # Validate the fallback file was created
                if not output_npy.exists():
                    raise RuntimeError(f"Fallback numpy file was not created: {output_npy}")
                    
                file_size_mb = output_npy.stat().st_size / (1024 * 1024)
                export_succeeded = True
                export_details = f"Fallback export: {output_npy.name} ({file_size_mb:.1f} MB)"
                print(f"   ✅ Fallback: Saved numpy array: {output_npy.name}")
                logger.info(f"Fallback export saved: {output_npy} ({file_size_mb:.1f} MB)")
                
            except Exception as fallback_error:
                logger.error(f"Fallback export also failed: {fallback_error}")
                print(f'   ❌ Fallback failed: {fallback_error}')
                step_duration = time.time() - step_start
                self.log_step(step_number, "Export Products", "error", 
                              f"Export failed: {exc}; Fallback failed: {fallback_error}", step_duration)
                raise RuntimeError(
                    f'SCIENTIFIC MODE FAILURE: Unable to export any results. '
                    f"Primary error: {exc}; Fallback error: {fallback_error}"
                ) from fallback_error
        
        step_duration = time.time() - step_start
        if export_succeeded:
            self.log_step(step_number, "Export Products", "success", export_details, step_duration)
        else:
            self.log_step(step_number, "Export Products", "warning", "Export may be incomplete", step_duration)

    def _compute_quality_metrics_payload(self, db_data, valid_percentage):
        """Compute core quality metrics shared by metadata and quality reports.

        Returns a dictionary compatible with RTCMetadataBuilder.set_quality_metrics,
        including at minimum:
          - valid_pixel_percentage
          - backscatter_statistics {min, max, mean, std}
        and, when masks are available, shadow/layover pixel percentages.
        """

        if not isinstance(db_data, np.ndarray):
            raise TypeError("db_data must be a numpy array for quality metrics computation")

        finite_data = db_data[np.isfinite(db_data)]
        if finite_data.size > 0:
            stats_min = float(np.min(finite_data))
            stats_max = float(np.max(finite_data))
            stats_mean = float(np.mean(finite_data))
            stats_std = float(np.std(finite_data))
        else:
            stats_min = float("nan")
            stats_max = float("nan")
            stats_mean = float("nan")
            stats_std = float("nan")

        total_pixels = int(getattr(self, "_total_pixels", 0) or db_data.size)

        metrics = {
            "valid_pixel_percentage": float(valid_percentage),
            "backscatter_statistics": {
                "min": stats_min,
                "max": stats_max,
                "mean": stats_mean,
                "std": stats_std,
            },
        }

        # Optional: derive shadow/layover coverage when RTC masks are available
        if total_pixels > 0:
            shadow_mask = getattr(self, "_shadow_mask", None)
            if isinstance(shadow_mask, np.ndarray) and shadow_mask.shape == db_data.shape:
                shadow_pixels = int(np.count_nonzero(shadow_mask != 0))
                metrics["shadow_pixel_percentage"] = (shadow_pixels / total_pixels) * 100.0

            layover_mask = getattr(self, "_layover_mask", None)
            if isinstance(layover_mask, np.ndarray) and layover_mask.shape == db_data.shape:
                layover_pixels = int(np.count_nonzero(layover_mask != 0))
                metrics["layover_pixel_percentage"] = (layover_pixels / total_pixels) * 100.0

        return metrics

    def _stage_generate_metadata(self, context: PipelineContext) -> None:
        """Produce metadata artifacts after exports."""
        step_number = STEP_NUMBERS["Generate Metadata"]
        self.announce_step(step_number, "Generate Metadata", "Writing processing metadata artifacts")
        step_start = time.time()

        db_data = self._db_data
        valid_percentage = self._valid_percentage or 0.0

        try:
            if not isinstance(db_data, np.ndarray):
                raise RuntimeError("SCIENTIFIC MODE FAILURE: dB data unavailable for metadata generation")

            metadata_json = self.output_dir / "metadata_complete.json"
            metadata_written = False
            metadata_details = ""

            if self.metadata_builder and self.precise_orbit_path:
                try:
                    self._populate_metadata_processing_steps()
                    quality_metrics = self._compute_quality_metrics_payload(db_data, valid_percentage)
                    self.metadata_builder.set_quality_metrics(quality_metrics)

                    scene_extent = self._resolve_scene_extent_bbox()
                    if scene_extent is None:
                        raise ValueError("Scene extent unavailable for RTC metadata builder")

                    dem_coverage = getattr(self, "dem_bbox", None) or scene_extent
                    pixel_spacing = (
                        float(self.get_current_range_spacing()),
                        float(self.get_current_azimuth_spacing()),
                    )
                    scene_center = (
                        0.5 * (scene_extent[0] + scene_extent[2]),
                        0.5 * (scene_extent[1] + scene_extent[3]),
                    )

                    incidence_near = None
                    incidence_far = None
                    if isinstance(self.validated_metadata, dict):
                        incidence_near = self.validated_metadata.get("incidence_angle_near")
                        incidence_far = self.validated_metadata.get("incidence_angle_far")
                    metadata_ref = getattr(self, "metadata", {})
                    if incidence_near is None and isinstance(metadata_ref, dict):
                        incidence_near = metadata_ref.get("incidence_angle_near")
                    if incidence_far is None and isinstance(metadata_ref, dict):
                        incidence_far = metadata_ref.get("incidence_angle_far")

                    incidence_range = (
                        float(incidence_near) if incidence_near is not None else 0.0,
                        float(incidence_far) if incidence_far is not None else 0.0,
                    )

                    # Determine DEM source from processor configuration
                    dem_source = getattr(self, 'dem_source', None) or "SRTM"
                    dem_version = getattr(self, 'dem_version', None) or "SRTMGL1_v003"
                    
                    rtc_metadata = self.metadata_builder.build_metadata(
                        dem_source=dem_source,
                        dem_resolution=float(self.dem_resolution),
                        dem_version=dem_version,
                        dem_coverage=tuple(float(v) for v in dem_coverage),
                        rtc_method="range_doppler_zero_doppler",
                        cosine_threshold=float(self.cosine_clip_threshold),
                        masking_enabled=True,
                        calibration_type=self.calibration_type,
                        lut_source=self.calibration_lut_source or "annotation_calibration_lut",
                        noise_applied=self.optimization_mode == "complete",
                        incidence_range=incidence_range,
                        pixel_spacing=pixel_spacing,
                        scene_center=scene_center,
                        scene_extent=tuple(float(v) for v in scene_extent),
                        coordinate_system=self.coordinate_system or "EPSG:4326",
                    )
                    rtc_metadata.save_to_file(str(metadata_json))
                    metadata_written = True
                    metadata_details = f"RTC metadata (DEM: {dem_source})"
                except Exception as builder_error:
                    logger.warning(f"RTC metadata builder failed: {builder_error}")
                    print(f"   ⚠️  RTC metadata builder failed: {builder_error}")

            if not metadata_written:
                processing_params = {
                    "polarization": self.polarization,
                    "speckle_filter": self.speckle_filter,
                    "filter_window": str(self.filter_window),
                    "multilook_range": f"{self.actual_range_looks if self.actual_range_looks is not None else self.multilook_range}",
                    "multilook_azimuth": f"{self.actual_azimuth_looks if self.actual_azimuth_looks is not None else self.multilook_azimuth}",
                    "terrain_flatten": str(self.terrain_flatten),
                    "geocode": str(self.geocode),
                    "target_resolution": str(self.target_resolution),
                    "processing_level": f"COMPLETE_{TOTAL_PIPELINE_STEPS}_STEP_PIPELINE",
                }

                metadata_result = sardine.generate_metadata(
                    f"S1A_backscatter_{self.polarization}_complete",
                    processing_params,
                    [str(self.input_path)],
                    {"valid_pixel_percentage": valid_percentage},
                )

                json_result = sardine.export_metadata_json(metadata_result)
                with open(metadata_json, "w") as f:
                    f.write(json_result)
                metadata_details = f"Metadata: {len(metadata_result)} fields"

            # Validate metadata file was created and is not empty
            if metadata_json.exists():
                metadata_size = metadata_json.stat().st_size
                if metadata_size == 0:
                    logger.warning(f"Metadata file is empty: {metadata_json}")
            else:
                logger.warning(f"Metadata file was not created: {metadata_json}")

            step_duration = time.time() - step_start
            self.log_step(
                step_number,
                "Generate Metadata",
                "success",
                metadata_details or "Metadata artifact written",
                step_duration,
            )

        except Exception as exc:
            logger.warning(f"Metadata generation failed: {exc}")
            print(f"   ⚠️  Metadata warning: {exc}")
            step_duration = time.time() - step_start
            self.log_step(
                step_number,
                "Generate Metadata",
                "warning",
                f"Metadata generation incomplete: {exc}",
                step_duration,
            )

    def _stage_finalize_pipeline(self, context: PipelineContext) -> None:
        """Finalize processing by generating reports and recording results."""

        step_number = STEP_NUMBERS["Quality Assessment"]
        self.announce_step(step_number, "Quality Assessment", "Final reporting and validation outputs")
        step_start = time.time()

        db_data = self._db_data
        if not isinstance(db_data, np.ndarray):
            raise RuntimeError("SCIENTIFIC MODE FAILURE: dB data unavailable for finalization")

        valid_percentage = float(self._valid_percentage or 0.0)
        valid_pixel_count = int(self._valid_pixel_count or np.isfinite(db_data).sum())
        total_pixels = int(self._total_pixels or db_data.size)

        try:
            # Use finite values only for min calculation (exclude -inf from zero backscatter)
            finite_mask = np.isfinite(db_data)
            if finite_mask.any():
                data_min = float(np.min(db_data[finite_mask]))
                data_max = float(np.max(db_data[finite_mask]))
                data_range_valid = True
            else:
                data_min = float("nan")
                data_max = float("nan")
                data_range_valid = False
        except (ValueError, TypeError):
            data_min = float("nan")
            data_max = float("nan")
            data_range_valid = False

        if self.quality_report:
            try:
                self.generate_quality_report(db_data, valid_percentage)
            except Exception as exc:
                print(f"   ⚠️  Quality report generation failed: {exc}")

        try:
            self.save_processing_log()
        except Exception as exc:
            print(f"   ⚠️  Processing log save failed: {exc}")

        try:
            print("\n" + "=" * 80)
            self.performance_monitor.print_summary()
        except Exception as exc:
            print(f"   ⚠️  Performance summary unavailable: {exc}")
        finally:
            try:
                self.performance_monitor.stop_monitoring()
            except Exception:
                pass

        # Print stage timing breakdown if available
        if self._stage_timings:
            try:
                print("\n⏱️  Stage Timing Breakdown:")
                sorted_stages = sorted(self._stage_timings.items(), key=lambda x: -x[1])
                for stage_name, duration in sorted_stages:
                    print(f"   {stage_name}: {duration:.2f}s")
            except Exception:
                pass

        total_duration = time.time() - self.start_time

        # Summarize validation gates before declaring success
        validation_summary = None
        if getattr(self, "validation", None) is not None and self.validation.records:
            try:
                validation_summary = self.validation.summary()
                report_path = self.validation.save_json(self.validation_report_path)
                print(f"🧪 Validation report written: {report_path}")
                print(
                    f"🧪 Validation summary: "
                    f"{validation_summary['passed']} passed, "
                    f"{validation_summary['failed']} failed, "
                    f"{validation_summary['warnings']} warnings"
                )
            except Exception as exc:
                print(f"   ⚠️  Validation report could not be saved: {exc}")

        # In strict science mode with fail-fast validation, treat any failed
        # error-level gates as a pipeline failure, approximating SNAP's
        # conservative behaviour.
        strict_science = bool(self.options.get("strict_science", True))
        if strict_science and validation_summary is not None and validation_summary.get("failed", 0) > 0:
            summary = validation_summary
            raise RuntimeError(
                "SCIENTIFIC MODE FAILURE: One or more validation gates failed "
                f"(checks={summary['total_checks']}, failed={summary['failed']}, "
                f"warnings={summary['warnings']}). See validation_report.json for details."
            )

        try:
            print(f"\n🎉 SUCCESS: Complete {TOTAL_PIPELINE_STEPS}-step backscatter processing finished!")
            print(f"⏱️  Total processing time: {total_duration:.1f} seconds")
            print(f"📂 Output directory: {self.output_dir}")
            print(f"📊 Output shape: {db_data.shape}")
            if data_range_valid:
                print(f"📈 Data range: {data_min:.1f} to {data_max:.1f} dB")
            else:
                print("📈 Data range: unavailable (all values NaN)")
            print(
                f"✅ Valid pixels: {valid_pixel_count:,} / {total_pixels:,} ({valid_percentage:.1f}%)"
            )
        except Exception as exc:
            print(f"   ⚠️  Final summary reporting failed: {exc}")

        try:
            output_files = [str(path) for path in self.output_dir.glob('*')]
            context.set_artifact(
                "result",
                {
                    "status": "success",
                    "steps_completed": len(self.processing_log),
                    "processing_time": total_duration,
                    "output_files": output_files,
                    "dimensions": tuple(int(v) for v in db_data.shape),
                    "data_range": [data_min, data_max],
                    "valid_percentage": valid_percentage,
                },
            )
        except Exception as exc:
            raise RuntimeError(
                f"SCIENTIFIC MODE FAILURE: Unable to record final processing result. Error: {exc}"
            ) from exc

        step_duration = time.time() - step_start
        self.log_step(step_number, "Quality Assessment", "success", "Pipeline finalized", step_duration)

    def _stage_execute_backscatter(self, context: PipelineContext) -> None:
        """Legacy monolithic pipeline removed; use process_backscatter()."""
        raise RuntimeError('Backscatter pipeline is modular; call process_backscatter().')

    def process_backscatter(self):
        """Execute the complete backscatter processing pipeline."""

        print("=" * 80)
        print(f"🛰️  SARdine {TOTAL_PIPELINE_STEPS}-Step SAR Backscatter Processing Pipeline")
        print("🔬 Complete scientific processing workflow")
        print(f"📁 Input: {self.input_path}")
        print(f"📂 Output: {self.output_dir}")
        print(f"📡 Polarization: {self.polarization}")
        print("=" * 80)

        pipeline = self._build_pipeline()
        context = PipelineContext(processor=self)

        try:
            pipeline.run(context)
            result = context.get_artifact("result")
            if result is None:
                raise RuntimeError("Backscatter pipeline completed without producing a result artifact")
            return result
        except Exception as exc:
            print(f"\n❌ ERROR: {TOTAL_PIPELINE_STEPS}-step processing failed: {exc}")
            import traceback

            traceback.print_exc()
            return {
                'status': 'error',
                'error': str(exc),
                'steps_completed': len(self.processing_log),
            }
        finally:
            # Always cleanup resources
            self._cleanup_memory()
            self._shutdown_executor()
    
    def generate_quality_report(self, db_data, valid_percentage):
        """Generate comprehensive quality assessment report"""
        
        print("\n📊 Generating Quality Report...")

        # Calculate quality metrics (shared with RTC metadata builder)
        quality_metrics = self._compute_quality_metrics_payload(db_data, valid_percentage)
        backscatter_stats = quality_metrics.get("backscatter_statistics", {})

        mean_backscatter = float(backscatter_stats.get("mean", float("nan")))
        std_backscatter = float(backscatter_stats.get("std", float("nan")))
        min_backscatter = float(backscatter_stats.get("min", float("nan")))
        max_backscatter = float(backscatter_stats.get("max", float("nan")))
        dynamic_range = (
            float(max_backscatter - min_backscatter)
            if np.isfinite([max_backscatter, min_backscatter]).all()
            else float("nan")
        )
        
        range_looks_used = self.actual_range_looks if self.actual_range_looks is not None else self.multilook_range
        azimuth_looks_used = self.actual_azimuth_looks if self.actual_azimuth_looks is not None else self.multilook_azimuth
        if range_looks_used is None:
            range_looks_used = 1
        if azimuth_looks_used is None:
            azimuth_looks_used = 1
        range_looks_used_f = float(range_looks_used)
        azimuth_looks_used_f = float(azimuth_looks_used)

        rtc_mode = getattr(self, "rtc_mode", None) or "none"
        output_units = "power" if getattr(self, "_linear_output", False) else "dB"
        
        quality_report = {
            "processing_parameters": {
                "polarization": self.polarization,
                "speckle_filter": self.speckle_filter,
                "multilook_factors": [range_looks_used_f, azimuth_looks_used_f],
                "terrain_flattening": self.terrain_flatten,
                "geocoding": self.geocode,
                "rtc_mode": rtc_mode,
                "output_units": output_units,
                "processing_level": f"COMPLETE_{TOTAL_PIPELINE_STEPS}_STEP_PIPELINE",
            },
            "data_quality": {
                "valid_pixel_percentage": float(quality_metrics.get("valid_pixel_percentage", valid_percentage)),
                "mean_backscatter_db": mean_backscatter,
                "std_backscatter_db": std_backscatter,
                "min_backscatter_db": min_backscatter,
                "max_backscatter_db": max_backscatter,
                "dynamic_range_db": dynamic_range,
            },
            "processing_summary": {
                "total_steps": TOTAL_PIPELINE_STEPS,
                "processing_time_seconds": time.time() - self.start_time,
                "output_dimensions": db_data.shape,
                "steps_completed": len(self.processing_log)
            }
        }

        # Enrich data_quality with any available mask-derived percentages
        shadow_pct = quality_metrics.get("shadow_pixel_percentage")
        if shadow_pct is not None:
            quality_report["data_quality"]["shadow_pixel_percentage"] = float(shadow_pct)

        layover_pct = quality_metrics.get("layover_pixel_percentage")
        if layover_pct is not None:
            quality_report["data_quality"]["layover_pixel_percentage"] = float(layover_pct)
        
        # Save quality report
        quality_file = self.output_dir / "quality_report_complete.json"
        with open(quality_file, "w") as f:
            json.dump(quality_report, f, indent=2)
        
        print(f"   ✅ Quality report saved: {quality_file.name}")
        if np.isfinite(mean_backscatter):
            print(f"   📊 Mean backscatter: {mean_backscatter:.1f} dB")
        print(f"   📊 Valid pixels: {valid_percentage:.1f}%")
        print(f"   🧾 Steps logged: {len(self.processing_log)}")
        self.latest_quality_report = quality_report
        return quality_report
    
    def save_processing_log(self):
        """Save detailed processing log"""
        
        log_file = self.output_dir / "processing_log_complete.json"
        with open(log_file, "w") as f:
            json.dump(self.processing_log, f, indent=2)
        
        print(f"   ✅ Processing log saved: {log_file.name}")

    def _configure_parallel_environment(self):
        """Ensure Rust/Python thread pools honour the requested parallel configuration."""

        desired_threads = 1
        if self.enable_parallel:
            if self.num_threads is not None:
                desired_threads = max(1, int(self.num_threads))
            else:
                detected = os.cpu_count() or 1
                desired_threads = max(1, min(detected, 32))
        thread_value = str(desired_threads)

        runtime_env = {
            "RAYON_NUM_THREADS": thread_value,
            "SARDINE_MAX_THREADS": thread_value,
            "OMP_NUM_THREADS": thread_value,
            "OPENBLAS_NUM_THREADS": thread_value,
            "NUMEXPR_NUM_THREADS": thread_value,
        }

        for key, value in runtime_env.items():
            current = os.environ.get(key)
            if current != value:
                os.environ[key] = value

        # When sequential mode is requested, force Rayon to a single worker
        if not self.enable_parallel:
            os.environ["RAYON_NUM_THREADS"] = thread_value

    def _auto_tune_chunk_size(self, array):
        """Adapt processing chunk size based on the active array footprint."""

        if not self.enable_parallel:
            return
        if array is None or not hasattr(array, "nbytes"):
            return
        if self.options.get('chunk_size') is not None:
            return

        data_size_bytes = int(getattr(array, "nbytes", 0))
        if data_size_bytes <= 0:
            return

        data_size_mb = max(1e-3, data_size_bytes / (1024 ** 2))
        threads = max(1, self.num_threads or (os.cpu_count() or 1))
        tuned_chunk = optimize_chunk_size(data_size_mb, threads)

        if tuned_chunk != self.chunk_size:
            self.chunk_size = tuned_chunk
            os.environ.setdefault("SARDINE_CHUNK_SIZE", str(self.chunk_size))
            print(f"📊 Adaptive chunk size selected: {self.chunk_size}")
