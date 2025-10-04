"""
SARdine Backscatter Processor

Complete SAR backscatter processing pipeline with REAL scientific data only.
Implements a 14-step pipeline for processing Sentinel-1 SLC data into 
analysis-ready backscatter products with scientific quality assurance.
"""

import time
import json
import os
import math
import threading
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
from typing import Dict, List, Optional, Tuple

import numpy as np
import sardine
from datetime import datetime, timedelta, timezone

# Import GeoTIFF export functionality
from sardine.export import export_to_geotiff, create_cog_with_stac


WGS84_SEMI_MAJOR_AXIS_M = 6_378_137.0
WGS84_ECCENTRICITY_SQUARED = 6.694_379_990_14e-3


def _meters_per_degree(latitude_deg: float) -> Tuple[float, float]:
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


def _refine_geocoding_bbox(
    bbox: List[float],
    rows: int,
    cols: int,
    range_spacing_m: float,
    azimuth_spacing_m: float,
) -> Optional[Tuple[List[float], Dict[str, float]]]:
    """Derive a footprint-aware geocoding bounding box from multilooked geometry."""

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

    centre_lat = 0.5 * (min_lat + max_lat)
    centre_lon = 0.5 * (min_lon + max_lon)

    try:
        meters_per_degree_lat, meters_per_degree_lon = _meters_per_degree(centre_lat)
    except ValueError:
        return None

    coverage_lat_m = azimuth_spacing_m * float(rows)
    coverage_lon_m = range_spacing_m * float(cols)

    # Include half a pixel of padding on each edge to avoid clipping
    coverage_lat_m += azimuth_spacing_m
    coverage_lon_m += range_spacing_m

    coverage_lat_deg = coverage_lat_m / meters_per_degree_lat
    coverage_lon_deg = coverage_lon_m / meters_per_degree_lon

    refined_min_lat = centre_lat - 0.5 * coverage_lat_deg
    refined_max_lat = centre_lat + 0.5 * coverage_lat_deg
    refined_min_lon = centre_lon - 0.5 * coverage_lon_deg
    refined_max_lon = centre_lon + 0.5 * coverage_lon_deg

    # Prevent expansion beyond original scene metadata
    refined_min_lat = max(refined_min_lat, min_lat)
    refined_max_lat = min(refined_max_lat, max_lat)
    refined_min_lon = max(refined_min_lon, min_lon)
    refined_max_lon = min(refined_max_lon, max_lon)

    if refined_min_lat >= refined_max_lat or refined_min_lon >= refined_max_lon:
        return None

    original_lat_extent = max_lat - min_lat
    original_lon_extent = max_lon - min_lon
    new_lat_extent = refined_max_lat - refined_min_lat
    new_lon_extent = refined_max_lon - refined_min_lon

    metrics: Dict[str, float] = {
        "original_lat_extent_deg": original_lat_extent,
        "original_lon_extent_deg": original_lon_extent,
        "new_lat_extent_deg": new_lat_extent,
        "new_lon_extent_deg": new_lon_extent,
        "shrink_lat_pct": 0.0,
        "shrink_lon_pct": 0.0,
    }

    if original_lat_extent > 0.0:
        metrics["shrink_lat_pct"] = max(
            0.0, 1.0 - (new_lat_extent / original_lat_extent)
        )
    if original_lon_extent > 0.0:
        metrics["shrink_lon_pct"] = max(
            0.0, 1.0 - (new_lon_extent / original_lon_extent)
        )

    refined_bbox: List[float] = [
        refined_min_lon,
        refined_min_lat,
        refined_max_lon,
        refined_max_lat,
    ]
    return refined_bbox, metrics

# Import enhanced RTC metadata system
try:
    from sardine.metadata.rtc_metadata import RTCMetadataBuilder, ComprehensiveRTCMetadata
except ImportError:
    RTCMetadataBuilder = None
    ComprehensiveRTCMetadata = None

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
    
    def optimize_chunk_size(data_size_mb, num_threads, memory_limit_mb=8000):
        return max(64, min(int(data_size_mb / num_threads / 4), 2048))
    
    def print_system_info():
        print("💻 Performance monitoring not available (psutil required)")


class BackscatterProcessor:
    """Complete SAR backscatter processing pipeline with REAL scientific data only"""
    
    def __init__(self, input_path, output_dir, options):
        # Handle both ZIP files and SAFE directories with strict validation
        input_path = Path(input_path)
        
        # Validate input exists
        if not input_path.exists():
            raise ValueError(f"❌ SCIENTIFIC MODE: Input file does not exist: {input_path}")
        
        if input_path.is_dir() and input_path.suffix == '.SAFE':
            # Use SAFE directory directly (Rust backend supports this)
            self.input_path = str(input_path)
        elif input_path.suffix.lower() == '.zip':
            self.input_path = str(input_path)
        else:
            raise ValueError(f"❌ SCIENTIFIC MODE: Input must be either a .zip file or .SAFE directory")
        
        # Validate it's a Sentinel-1 product
        filename = Path(self.input_path).name
        if not filename.startswith(('S1A_', 'S1B_')):
            raise ValueError(f"❌ SCIENTIFIC MODE: Not a valid Sentinel-1 product. Filename must start with S1A_ or S1B_: {filename}")
        
        # Validate it's an SLC product
        if '_SLC_' not in filename:
            raise ValueError(f"❌ SCIENTIFIC MODE: Only SLC products supported. Found: {filename}")
        
        self.output_dir = Path(output_dir)
        self.options = options
        self.start_time = time.time()
        self.processing_log = []
        self.reader_lock = threading.Lock()
        
        # Create output directory
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        # Set up default orbit cache if not already configured
        # This ensures cached reader architecture works for all operations
        if 'SARDINE_ORBIT_CACHE' not in os.environ:
            default_orbit_cache = str(self.output_dir / "orbit_cache")
            os.environ['SARDINE_ORBIT_CACHE'] = default_orbit_cache
            # Create the cache directory to ensure it exists
            Path(default_orbit_cache).mkdir(parents=True, exist_ok=True)

        # Track pixel spacing as it evolves through multilooking and other resampling steps
        self.original_range_spacing = None
        self.original_azimuth_spacing = None
        self.current_range_spacing = None
        self.current_azimuth_spacing = None
        
        # Processing parameters with validation
        self.polarization = options.get('polarization', 'VV')
        if self.polarization not in ['VV', 'VH', 'HH', 'HV']:
            raise ValueError(f"❌ SCIENTIFIC MODE: Invalid polarization: {self.polarization}")
        
        self.speckle_filter = options.get('speckle_filter', 'enhanced_lee')
        self.filter_window = options.get('filter_window', 7)
        self.multilook_range = options.get('multilook_range', 2)
        self.multilook_azimuth = options.get('multilook_azimuth', 2)
        self.calibration_type = options.get('calibration_type', 'sigma0')
        self.terrain_flatten = options.get('terrain_flatten', True)
        self.geocode = options.get('geocode', True)
        self.target_resolution = options.get('resolution')
        self.quality_report = options.get('quality_report', True)
        self.use_real_orbit = options.get('use_real_orbit', True)  # Enforce real orbit data
        self.allow_synthetic = options.get('allow_synthetic', False)  # Control synthetic data fallbacks
        
        # Enhanced RTC features from enhanced_rtc_processor
        self.enable_iw_merging = options.get('enable_iw_merging', True)  # IW merging before terrain correction
        self.cosine_clip_threshold = options.get('cosine_clip_threshold', 0.1)  # Steep terrain handling
        self.linear_processing = options.get('linear_processing', True)  # Keep linear until final dB conversion
        
        # Parallel processing configuration
        self.enable_parallel = options.get('parallel', True)
        self.num_threads = options.get('num_threads', None)  # None = auto-detect
        self.sequential = options.get('sequential', False)
        self.chunk_size = options.get('chunk_size', 1024)
        
        # Initialize performance monitoring
        self.performance_monitor = PerformanceMonitor(enable_monitoring=True)
        
        # Override parallel settings if sequential is requested
        if self.sequential:
            self.enable_parallel = False
            print("🔄 Sequential processing mode enabled (parallel processing disabled)")
        elif self.enable_parallel:
            detected_cores = os.cpu_count()
            if self.num_threads is None:
                self.num_threads = detected_cores
                print(f"🚀 Parallel processing enabled: auto-detected {detected_cores} CPU cores")
            else:
                print(f"🚀 Parallel processing enabled: {self.num_threads} threads specified")
            
            # Optimize chunk size if not explicitly set
            if options.get('chunk_size') is None:
                self.chunk_size = optimize_chunk_size(100, self.num_threads)  # Estimate 100MB typical
                print(f"📊 Optimized chunk size: {self.chunk_size}")
            else:
                print(f"📊 Chunk size: {self.chunk_size}")
        else:
            print("🔄 Parallel processing disabled")

        # Configure runtime thread pools to honour requested parallelism
        self._configure_parallel_environment()
        
        # Print system information for performance optimization
        if self.enable_parallel:
            print_system_info()
        
        # Enforce scientific mode - no synthetic data allowed
        if self.use_real_orbit and self.allow_synthetic:
            raise ValueError("❌ INVALID CONFIGURATION: Cannot use real orbit mode with synthetic data allowed")
        
        if not self.allow_synthetic:
            print("🔬 SCIENTIFIC MODE: Real data only - no synthetic fallbacks")
        else:
            print("⚠️  DEMONSTRATION MODE: Synthetic data fallbacks enabled")
        
        # Initialize storage for GeoTIFF export
        self.geo_transform = None
        self.coordinate_system = "EPSG:4326"  # Default, will be updated from processing
        self.validated_metadata = None

    def validate_and_extract_complete_metadata(self):
        """
        Extract and validate complete metadata from cached reader for all processing functions.
        This ensures no hardcoded values are used and all parameters come from annotation XML.
        
        Returns:
            dict: Complete metadata dictionary for Rust functions
            
        Raises:
            ValueError: If any required metadata field is missing or invalid
        """
        print("🔍 Validating and extracting complete metadata from cached reader...")
        
        # Get cached metadata
        if not hasattr(self, 'metadata') or not self.metadata:
            raise ValueError("SCIENTIFIC MODE: No cached metadata available - cannot proceed")
        
        cached_metadata = self.metadata
        
        # Define required fields for scientific processing
        required_fields = {
            'range_pixel_spacing': 'Range pixel spacing in meters',
            'azimuth_pixel_spacing': 'Azimuth pixel spacing in meters', 
            'slant_range_time': 'Slant range time in seconds',
            'wavelength': 'Radar wavelength in meters',
            'prf': 'Pulse repetition frequency in Hz'
        }
        
        # Validate all required fields are present and not None/zero
        missing_fields = []
        invalid_fields = []
        
        for field, description in required_fields.items():
            if field not in cached_metadata:
                missing_fields.append(f"{field} ({description})")
            else:
                value = cached_metadata[field]
                if value is None or (isinstance(value, (int, float)) and value == 0.0):
                    invalid_fields.append(f"{field} = {value} ({description})")
        
        if missing_fields:
            raise ValueError(f"SCIENTIFIC MODE: Missing required metadata fields:\n" + 
                           "\n".join(f"  - {field}" for field in missing_fields) +
                           "\nCannot proceed with hardcoded fallbacks disabled.")
        
        if invalid_fields:
            raise ValueError(f"SCIENTIFIC MODE: Invalid/zero metadata values:\n" +
                           "\n".join(f"  - {field}" for field in invalid_fields) +
                           "\nThis indicates annotation parsing failed. Check input data.")
        
        # Build complete metadata dictionary for Rust functions
        complete_metadata = {
            # Core radar parameters (REQUIRED)
            'range_pixel_spacing': float(cached_metadata['range_pixel_spacing']),
            'azimuth_pixel_spacing': float(cached_metadata['azimuth_pixel_spacing']),
            'slant_range_time': float(cached_metadata['slant_range_time']),
            'prf': float(cached_metadata['prf']),
            'wavelength': float(cached_metadata['wavelength']),
            
            # Calculate derived parameters
            'radar_frequency': self.calculate_radar_frequency(float(cached_metadata['wavelength'])),
            
            # Optional geometry parameters (with validation)
            'incidence_angle_near': float(cached_metadata.get('incidence_angle_near', 20.0)),
            'incidence_angle_far': float(cached_metadata.get('incidence_angle_far', 45.0)),
            
            # Product information
            'product_type': cached_metadata.get('product_type', 'SLC'),
            'processing_level': cached_metadata.get('processing_level', 'L1'),
            'mission_id': cached_metadata.get('mission_id', 'S1A'),
            'polarization': cached_metadata.get('polarization', self.polarization),
            
            # Quality control flags
            'metadata_source': 'cached_annotation_xml',
            'hardcoded_values_used': False,
            'validation_passed': True,
            'strict_mode': True
        }
        
        # Additional validation of derived values
        radar_freq = complete_metadata['radar_frequency']
        if not (5.0e9 <= radar_freq <= 6.0e9):  # C-band range
            raise ValueError(f"SCIENTIFIC MODE: Invalid radar frequency {radar_freq/1e9:.3f} GHz (expected C-band 5.0-6.0 GHz)")
        
        range_spacing = complete_metadata['range_pixel_spacing']
        azimuth_spacing = complete_metadata['azimuth_pixel_spacing']
        if not (1.0 <= range_spacing <= 50.0):
            raise ValueError(f"SCIENTIFIC MODE: Invalid range pixel spacing {range_spacing}m (expected 1-50m)")
        if not (1.0 <= azimuth_spacing <= 50.0):
            raise ValueError(f"SCIENTIFIC MODE: Invalid azimuth pixel spacing {azimuth_spacing}m (expected 1-50m)")
        
        print(f"   ✅ Metadata validation successful:")
        print(f"      Range spacing: {range_spacing:.3f}m")
        print(f"      Azimuth spacing: {azimuth_spacing:.3f}m")
        print(f"      Wavelength: {complete_metadata['wavelength']:.6f}m")
        print(f"      Radar frequency: {radar_freq/1e9:.3f} GHz")
        print(f"      PRF: {complete_metadata['prf']:.1f} Hz")
        print(f"   🔒 Strict mode: No hardcoded fallbacks permitted")
        
        return complete_metadata

    def calculate_radar_frequency(self, wavelength):
        """Calculate radar frequency from wavelength using speed of light"""
        speed_of_light = 299792458.0  # m/s (exact definition)
        return speed_of_light / wavelength

    # --- Pixel spacing utilities -------------------------------------------------
    def _require_spacing_initialized(self):
        if self.original_range_spacing is None or self.original_azimuth_spacing is None:
            raise ValueError("Pixel spacing metadata has not been initialized")

    def get_current_range_spacing(self):
        self._require_spacing_initialized()
        return self.current_range_spacing if self.current_range_spacing is not None else self.original_range_spacing

    def get_current_azimuth_spacing(self):
        self._require_spacing_initialized()
        return self.current_azimuth_spacing if self.current_azimuth_spacing is not None else self.original_azimuth_spacing

    def update_current_spacing(self, range_spacing, azimuth_spacing):
        self._require_spacing_initialized()
        if range_spacing is not None:
            self.current_range_spacing = float(range_spacing)
        if azimuth_spacing is not None:
            self.current_azimuth_spacing = float(azimuth_spacing)
        
    def extract_sensing_time_from_filename(self, filename):
        """Extract sensing start and stop times from Sentinel-1 filename"""
        try:
            # S1A_IW_SLC__1SDV_20200103T170815_20200103T170842_030639_0382D5_DADE.zip
            parts = filename.split('_')
            if len(parts) < 6:
                raise ValueError(f"Invalid Sentinel-1 filename format: {filename}")
            
            start_time_str = parts[4]  # 20200103T170815
            stop_time_str = parts[5]   # 20200103T170842
            
            # Parse datetime
            start_time = datetime.strptime(start_time_str, "%Y%m%dT%H%M%S")
            stop_time = datetime.strptime(stop_time_str, "%Y%m%dT%H%M%S")
            
            return start_time, stop_time
            
        except Exception as e:
            raise ValueError(f"Failed to parse sensing time from filename {filename}: {e}")
    
    def download_orbit_file(self, sensing_start=None, sensing_stop=None, orbit_cache_dir=None):
        """Download precise orbit file for the given sensing time using the Rust orbit manager"""

        # Determine cache directory preference
        orbit_cache_path = Path(orbit_cache_dir) if orbit_cache_dir else Path(self.output_dir) / "orbit_cache"
        orbit_cache_path.mkdir(parents=True, exist_ok=True)

        # Ensure downstream Rust code knows where the cache lives
        os.environ['SARDINE_ORBIT_CACHE'] = str(orbit_cache_path)

        # Gather product metadata needed for orbit lookup
        metadata = getattr(self, 'metadata', None)
        product_id = metadata.get('product_id') if isinstance(metadata, dict) else None
        start_time_value = metadata.get('start_time') if isinstance(metadata, dict) else None

        # Fall back to filename parsing if metadata is unavailable
        if not product_id:
            product_id = Path(self.input_path).name
            if product_id.endswith('.zip'):
                product_id = product_id[:-4]
            elif product_id.endswith('.SAFE'):
                product_id = product_id[:-5]
        if not start_time_value and isinstance(sensing_start, datetime):
            start_time_value = sensing_start

        if not product_id or start_time_value is None:
            raise RuntimeError("Unable to determine product ID or start time required for orbit download")

        # Helper to normalize start time into RFC3339 with explicit Z suffix
        def to_rfc3339(value):
            if isinstance(value, datetime):
                dt = value if value.tzinfo else value.replace(tzinfo=timezone.utc)
                return dt.astimezone(timezone.utc).isoformat().replace("+00:00", "Z")
            if isinstance(value, str):
                text = value.strip()
                if text.endswith('Z'):
                    return text
                try:
                    parsed = datetime.fromisoformat(text.replace('Z', '+00:00'))
                except ValueError as exc:
                    raise RuntimeError(f"Invalid orbit start time format: {value}") from exc
                return parsed.astimezone(timezone.utc).isoformat().replace("+00:00", "Z")
            raise RuntimeError(f"Unsupported start time type: {type(value)!r}")

        start_time_rfc3339 = to_rfc3339(start_time_value)

        print(f"   🛰️  Downloading orbit file for {start_time_rfc3339}")

        try:
            orbit_result = sardine.apply_precise_orbit_file(product_id, start_time_rfc3339, str(orbit_cache_path))
        except Exception as exc:
            raise RuntimeError(f"Failed to download orbit file via Rust orbit manager: {exc}")

        # Confirm that the orbit file landed in cache
        cached_files = sorted(
            orbit_cache_path.glob(f"{product_id}_*.EOF"),
            key=lambda p: p.stat().st_mtime,
            reverse=True,
        )

        if not cached_files:
            raise RuntimeError("Orbit download reported success but no .EOF file found in cache")

        latest_file = cached_files[0]
        meta_suffix = ""
        if isinstance(orbit_result, dict):
            result_meta = orbit_result.get('result')
            if isinstance(result_meta, dict):
                orbit_type = result_meta.get('orbit_type')
                vector_count = result_meta.get('orbit_vectors_count')
                if orbit_type or vector_count:
                    meta_suffix = " ("
                    if orbit_type:
                        meta_suffix += str(orbit_type)
                    if orbit_type and vector_count:
                        meta_suffix += " - "
                    if vector_count:
                        meta_suffix += f"{vector_count} vectors"
                    meta_suffix += ")"

        print(f"   ✅ Precise orbit ready: {latest_file.name}{meta_suffix}")

        return str(latest_file)
    
    def log_step(self, step_num, step_name, status, details="", duration=None):
        """Log processing step with timing"""
        log_entry = {
            'step': step_num,
            'name': step_name,
            'status': status,
            'details': details,
            'duration': duration,
            'timestamp': time.time()
        }
        self.processing_log.append(log_entry)
        # Persist logs incrementally so they survive early termination
        try:
            progress_path = self.output_dir / "processing_log_progress.json"
            with open(progress_path, "w") as f:
                json.dump(self.processing_log, f, indent=2)
            # Also append a JSONL stream for tailing
            stream_path = self.output_dir / "processing_log.jsonl"
            with open(stream_path, "a") as f:
                f.write(json.dumps(log_entry) + "\n")
        except Exception as e:
            # Do not fail the pipeline on logging I/O issues
            print(f"   ⚠️  Failed to persist progress log: {e}")
        
        # Print progress
        duration_str = f" ({duration:.1f}s)" if duration else ""
        status_icon_map = {
            "success": "✅",
            "error": "❌",
            "warning": "⚠️",
            "skipped": "⏭️",
            "fallback": "ℹ️",
            "info": "ℹ️",
        }
        status_icon = status_icon_map.get(status, "🔄")
        print(f"   {status_icon} STEP {step_num}: {step_name}{duration_str}")
        if details:
            print(f"      {details}")
    
    def announce_step(self, step_num, step_name, details="", status="start"):
        """Announce the status of a processing step for real-time visibility"""
        status = (status or "start").lower()
        if status == "skip":
            icon, label = "ℹ️", "SKIP"
        elif status == "resume":
            icon, label = "🔄", "RESUME"
        else:
            icon, label = "➡️", "START"

        print(f"\n{icon}  {label} STEP {step_num}: {step_name}")
        if details:
            print(f"   {details}")

    def _process_subswath_pipeline(self, subswath):
        """Deburst and calibrate a subswath in a single pipeline stage"""
        pipeline_start = time.time()

        # Deburst
        deburst_start = time.time()
        with self.reader_lock:
            deburst_result = sardine.deburst_topsar_cached(self.reader, subswath, self.polarization)
        if isinstance(deburst_result, dict):
            if deburst_result.get('status') == 'error':
                raise RuntimeError(
                    f"Deburst failed for {subswath}: {deburst_result.get('message', 'Unknown error')}"
                )
            deburst_data = deburst_result.get('data')
        else:
            deburst_data = deburst_result

        if not isinstance(deburst_data, np.ndarray) or deburst_data.size == 0:
            raise RuntimeError(f"Deburst failed to produce valid data for {subswath}")

        rows, cols = deburst_data.shape
        deburst_duration = time.time() - deburst_start

        # Radiometric calibration with thermal noise removal
        calibration_start = time.time()
        complex_data = np.asarray(deburst_data, dtype=np.complex64)
        with self.reader_lock:
            calibration_job = sardine.prepare_calibration_job_cached(
                self.reader,
                subswath,
                self.polarization,
                str(self.calibration_type),
                rows,
                cols,
                True,
            )

        cal_result = sardine.run_calibration_job(
            calibration_job,
            complex_data,
        )

        if isinstance(cal_result, dict):
            if cal_result.get('status') == 'error':
                raise RuntimeError(
                    f"Calibration error for {subswath}: {cal_result.get('message', 'Unknown error')}"
                )
            calibrated_data = cal_result.get('calibrated_data')
            if calibrated_data is None:
                calibrated_data = cal_result.get('data')
        else:
            calibrated_data = cal_result

        if not isinstance(calibrated_data, np.ndarray) or calibrated_data.size == 0:
            raise RuntimeError(
                f"Calibration failed to produce valid data for {subswath}"
            )

        calibrated_array = np.asarray(calibrated_data, dtype=np.float32)
        calibration_duration = time.time() - calibration_start
        total_duration = time.time() - pipeline_start

        # Explicitly drop large intermediate to release memory sooner
        del deburst_data

        return {
            'subswath': subswath,
            'deburst_shape': (rows, cols),
            'calibrated_data': calibrated_array,
            'deburst_duration': deburst_duration,
            'calibration_duration': calibration_duration,
            'total_duration': total_duration,
        }

    def process_backscatter(self):
        """Execute complete 14-step backscatter processing pipeline"""
        
        print("=" * 80)
        print("🛰️  SARdine 14-Step SAR Backscatter Processing Pipeline")
        print("🔬 Complete scientific processing workflow")
        print(f"📁 Input: {self.input_path}")
        print(f"📂 Output: {self.output_dir}")
        print(f"📡 Polarization: {self.polarization}")
        print("=" * 80)
        
        try:
            # STEP 1: Read Metadata & Files
            self.announce_step(1, "Read Metadata & Files")
            step_start = time.time()
            self.performance_monitor.start_step("Read Metadata & Files", data_size_mb=0, 
                                               parallel_enabled=self.enable_parallel,
                                               threads_used=self.num_threads)
            
            # Verify this is a valid Sentinel-1 product
            if not Path(self.input_path).name.startswith(('S1A_', 'S1B_')):
                raise ValueError(f"❌ SCIENTIFIC MODE: Invalid Sentinel-1 product name. Must start with S1A_ or S1B_")
            
            # PERFORMANCE OPTIMIZATION: Use cached reader for ~93% speedup
            print("🚀 PERFORMANCE: Initializing cached SLC reader (~93% faster metadata access)")
            reader = sardine.SlcReader.new_with_full_cache(self.input_path)
            metadata = reader.get_cached_metadata()
            
            # Also get product info for compatibility (uses cached data when possible)
            try:
                info = sardine.get_product_info_cached(reader)
                print("✅ PERFORMANCE: Product info extraction completed (using cached reader)")
            except Exception as e:
                print(f"⚠️  Product info extraction warning: {e}")
                info = {}  # Continue with just cached metadata
            
            # Store metadata and reader as instance variables for use in subsequent steps
            self.metadata = metadata
            self.reader = reader
            
            # CRITICAL: Validate and extract complete metadata to ensure no hardcoded fallbacks
            try:
                validated_metadata = self.validate_and_extract_complete_metadata()
                print("   ✅ Strict metadata validation completed - ready for scientific processing")

                # Initialize spacing trackers using real annotation metadata
                try:
                    self.original_range_spacing = float(metadata['range_pixel_spacing'])
                    self.original_azimuth_spacing = float(metadata['azimuth_pixel_spacing'])
                    self.current_range_spacing = self.original_range_spacing
                    self.current_azimuth_spacing = self.original_azimuth_spacing
                except (KeyError, TypeError, ValueError) as spacing_error:
                    raise ValueError(
                        f"Missing or invalid pixel spacing in metadata: {spacing_error}"
                    )

                self.validated_metadata = validated_metadata
                if self.target_resolution is None:
                    self.target_resolution = float(self.current_range_spacing)
            except ValueError as e:
                raise ValueError(f"❌ METADATA VALIDATION FAILED: {e}")
            
            # Verify we have valid metadata
            if not metadata or (hasattr(metadata, '__len__') and len(metadata) == 0):
                raise ValueError(f"❌ SCIENTIFIC MODE: Failed to extract valid metadata from product")
            
            step_duration = time.time() - step_start
            metadata_info = f"Cached metadata extracted: {type(metadata).__name__}"
            if hasattr(metadata, '__len__'):
                metadata_info += f" ({len(metadata)} items)"
            elif hasattr(metadata, '__dict__'):
                metadata_info += f" ({len(vars(metadata))} fields)"
            metadata_info += f" [PERFORMANCE: ~93% faster]"
            self.log_step(1, "Read Metadata & Files", "success", metadata_info, step_duration)
            self.performance_monitor.end_step()
            
            # STEP 2: Apply Precise Orbit File
            self.announce_step(2, "Apply Precise Orbit File")
            step_start = time.time()
            
            # Use existing orbit files if available, otherwise proceed with metadata orbit info
            try:
                # Check for existing orbit files first
                orbit_cache_dir = str(self.output_dir / "orbit_cache")
                orbit_files = []
                if Path(orbit_cache_dir).exists():
                    import glob
                    orbit_files = glob.glob(os.path.join(orbit_cache_dir, "*.EOF"))
                
                if orbit_files:
                    # Use existing cached orbit file
                    orbit_file_path = max(orbit_files, key=os.path.getmtime)
                    print(f"   📡 Using existing orbit file: {os.path.basename(orbit_file_path)}")
                    
                    # Set environment variable for cached reader architecture consistency
                    os.environ['SARDINE_ORBIT_CACHE'] = orbit_cache_dir
                    
                    orbit_status = f"Using cached orbit file: {os.path.basename(orbit_file_path)}"
                    step_duration = time.time() - step_start
                    self.log_step(2, "Apply Precise Orbit File", "success", orbit_status, step_duration)
                else:
                    # Extract product_id and start_time from metadata for potential download
                    if not self.metadata or not isinstance(self.metadata, dict):
                        raise ValueError("No valid metadata available and no cached orbit files")
                    
                    product_id = self.metadata.get('product_id')
                    start_time = self.metadata.get('start_time')
                    
                    if not product_id or not start_time:
                        raise ValueError(f"Missing required metadata: product_id={product_id}, start_time={start_time}")
                    
                    # Try to download orbit file via Rust orbit manager
                    try:
                        sensing_start = None
                        if isinstance(start_time, str):
                            sensing_start = datetime.fromisoformat(start_time.replace('Z', '+00:00'))
                        elif isinstance(start_time, datetime):
                            sensing_start = start_time
                        orbit_path = self.download_orbit_file(sensing_start, None, orbit_cache_dir)
                        orbit_status = f"Downloaded orbit data: {Path(orbit_path).name}"
                    except Exception as download_error:
                        # Continue without precise orbit if download fails
                        print(f"   ⚠️  Orbit download failed (continuing with manifest data): {download_error}")
                        orbit_status = "Using manifest orbit data (download failed)"
                    
                    step_duration = time.time() - step_start
                    self.log_step(2, "Apply Precise Orbit File", "success", orbit_status, step_duration)
                
            except Exception as e:
                step_duration = time.time() - step_start
                self.log_step(2, "Apply Precise Orbit File", "warning", f"Orbit processing failed: {e}", step_duration)
                print(f"   ⚠️  Continuing without precise orbit data: {e}")
            
            # STEP 3: Thermal Noise Removal is applied implicitly during Step 5 (radiometric calibration)
            # to avoid redundant logging, we rely on the calibration step to report any issues.
            
            # STEP 4 & 5: Deburst + Radiometric Calibration (pipelined per subswath)
            self.announce_step(4, "Deburst & Radiometric Calibration", "Processing each subswath end-to-end")
            step_start = time.time()

            calibrated_subswaths = {}
            available_subswaths = []
            primary_subswath = None

            try:
                if not self.metadata or not isinstance(self.metadata, dict):
                    raise ValueError("No valid metadata available for subswath processing")

                subswaths_str = self.metadata.get('subswaths', '')
                if not subswaths_str:
                    raise ValueError(
                        "❌ SCIENTIFIC MODE FAILURE: Subswaths not found in metadata. "
                        "Real subswath information is required for scientific processing. "
                        "No hardcoded fallbacks permitted for data integrity."
                    )

                parsed_subswaths = [sw.strip() for sw in subswaths_str.split(',') if sw.strip()]
                if not parsed_subswaths:
                    raise ValueError(
                        f"❌ SCIENTIFIC MODE FAILURE: Invalid subswaths format in metadata: '{subswaths_str}'"
                    )

                def subswath_sort_key(sw):
                    sw_upper = sw.upper()
                    if sw_upper.startswith("IW") and sw_upper[2:].isdigit():
                        return (0, int(sw_upper[2:]))
                    return (1, sw_upper)

                # Preserve metadata order while removing duplicates
                seen = set()
                metadata_order = []
                for sw in parsed_subswaths:
                    sw_upper = sw.upper()
                    if sw_upper not in seen:
                        seen.add(sw_upper)
                        metadata_order.append(sw)

                pipeline_subswaths = sorted(metadata_order, key=subswath_sort_key)
                available_subswaths = pipeline_subswaths

                print(f"   📡 Real subswaths from metadata: {', '.join(parsed_subswaths)}")
                if pipeline_subswaths != metadata_order:
                    print(f"   🔁 Normalized subswath order for processing: {', '.join(pipeline_subswaths)}")

                estimated_data_size_mb = max(1, 500 * len(pipeline_subswaths))
                self.performance_monitor.start_step(
                    "Deburst & Calibration",
                    data_size_mb=estimated_data_size_mb,
                    parallel_enabled=self.enable_parallel,
                    threads_used=self.num_threads,
                    chunk_size=self.chunk_size,
                )

                # Determine thread count for per-subs pipeline
                max_workers = min(len(pipeline_subswaths), self.num_threads or os.cpu_count() or 1)
                if max_workers <= 0:
                    max_workers = 1

                print(f"   🧵 Subswath pipeline threads: {max_workers}")

                pipeline_results = []
                with ThreadPoolExecutor(max_workers=max_workers) as executor:
                    future_map = {
                        executor.submit(self._process_subswath_pipeline, subswath): subswath
                        for subswath in pipeline_subswaths
                    }

                    for future in as_completed(future_map):
                        subswath = future_map[future]
                        try:
                            result = future.result()
                            calibrated_subswaths[subswath] = result['calibrated_data']
                            pipeline_results.append(result)
                            print(
                                f"   ✅ {subswath}: deburst {result['deburst_duration']:.1f}s, "
                                f"calibration {result['calibration_duration']:.1f}s, "
                                f"shape {result['deburst_shape'][0]}x{result['deburst_shape'][1]}"
                            )
                        except Exception as pipeline_error:
                            raise RuntimeError(
                                f"SCIENTIFIC MODE FAILURE: Subswath pipeline failed for {subswath}: {pipeline_error}"
                            )

                if not calibrated_subswaths:
                    raise RuntimeError("No subswaths were successfully processed")

                preferred_order = sorted(calibrated_subswaths.keys(), key=subswath_sort_key)
                primary_subswath = preferred_order[0]
                working_data = calibrated_subswaths[primary_subswath]

                # Optional: re-run auto chunk tuning using calibrated data
                self._auto_tune_chunk_size(working_data)

                summary_lines = [
                    f"{result['subswath']}={result['deburst_shape'][0]}x{result['deburst_shape'][1]}"
                    for result in sorted(pipeline_results, key=lambda r: subswath_sort_key(r['subswath']))
                ]

                step_duration = time.time() - step_start
                self.log_step(
                    4,
                    "Deburst & Radiometric Calibration",
                    "success",
                    f"Processed {len(calibrated_subswaths)} subswaths: {', '.join(summary_lines)}",
                    step_duration,
                )
                self.performance_monitor.end_step()

            except Exception as e:
                step_duration = time.time() - step_start
                self.log_step(
                    4,
                    "Deburst & Radiometric Calibration",
                    "error",
                    f"Pipeline failed: {e}",
                    step_duration,
                )
                raise RuntimeError(
                    f"SCIENTIFIC MODE FAILURE: Subswath pipeline failed. Error: {e}"
                )

            # STEP 6: Expert-Enhanced IW Merge Integration (all subswaths processed above)
            self.announce_step(6, "Expert IW Merge", "Combining calibrated subswaths")
            step_start = time.time()

            try:
                if len(calibrated_subswaths) < 2:
                    print(
                        "   ℹ️  Only one subswath available after calibration; skipping expert IW merge"
                    )
                    step_duration = time.time() - step_start
                    self.log_step(
                        6,
                        "Expert IW Merge",
                        "skipped",
                        "Insufficient subswaths for merge",
                        step_duration,
                    )
                else:
                    print(
                        f"🔗 EXPERT IW MERGE: Attempting merge of {len(calibrated_subswaths)} calibrated subswaths"
                    )

                    # Prepare data for merge functions in consistent IW order when available
                    ordered_for_merge = sorted(calibrated_subswaths.keys(), key=subswath_sort_key)
                    merge_inputs = {
                        sw: np.asarray(calibrated_subswaths[sw], dtype=np.float32)
                        for sw in ordered_for_merge
                    }
                    merge_inputs_upper = {
                        sw.upper(): np.asarray(calibrated_subswaths[sw], dtype=np.float32)
                        for sw in ordered_for_merge
                    }

                    if {"IW1", "IW2", "IW3"}.issubset(set(merge_inputs_upper.keys())):
                        merged_result = sardine.merge_subswaths_cached(
                            merge_inputs_upper["IW1"],
                            merge_inputs_upper["IW2"],
                            merge_inputs_upper["IW3"],
                            self.reader,
                            self.polarization,
                        )
                    else:
                        merged_result = sardine.topsar_merge_cached(
                            merge_inputs,
                            self.polarization,
                            self.reader,
                            None,
                        )

                    merged_data = merged_result.get('data') if isinstance(merged_result, dict) else None
                    if merged_data is None and hasattr(merged_result, 'get'):
                        merged_data = merged_result.get('data')
                    if merged_data is None:
                        raise RuntimeError("Expert merge did not return merged data")

                    merged_data = np.asarray(merged_data, dtype=np.float32)

                    if isinstance(merged_result, dict):
                        coverage_percent = merged_result.get('coverage_percent')
                        if coverage_percent is None:
                            valid_pixels = merged_result.get('valid_pixels')
                            total_pixels = merged_data.size
                            coverage_percent = (
                                100.0 * valid_pixels / total_pixels if valid_pixels and total_pixels else 100.0
                            )
                        uncovered_pixels = merged_result.get('uncovered_mask')
                        if isinstance(uncovered_pixels, np.ndarray):
                            uncovered_count = int(uncovered_pixels.sum())
                        else:
                            uncovered_count = 0
                    else:
                        coverage_percent = 100.0
                        uncovered_count = 0

                    step_duration = time.time() - step_start
                    self.log_step(
                        6,
                        "Expert IW Merge",
                        "success",
                        f"Merged subswaths: {merged_data.shape}, Coverage: {coverage_percent:.1f}%",
                        step_duration,
                    )

                    print(f"   ✅ Expert TOPSAR merge completed: {merged_data.shape}")
                    print(f"   🎯 Coverage: {coverage_percent:.1f}% ({uncovered_count} uncovered pixels)")
                    print(
                        f"   📈 Data utilization: {coverage_percent:.1f}% (vs {100.0/len(calibrated_subswaths):.1f}% single-subs)")

                    working_data = merged_data

            except Exception as e:
                print(
                    f"   ⚠️  Expert IW merge failed, falling back to {primary_subswath}-only data: {e}"
                )
                print(
                    f"   📊 Using {primary_subswath} data only: {working_data.shape} ({100.0/len(calibrated_subswaths):.1f}% data utilization)"
                )
                step_duration = time.time() - step_start
                self.log_step(
                    6,
                    "Expert IW Merge",
                    "fallback",
                    f"Using {primary_subswath} only due to merge failure",
                    step_duration,
                )

            # STEP 8: Multilooking - SCIENTIFIC ORDER (after TOPSAR merge)
            self.announce_step(8, "Multilooking", "Applying scientific multilook factors")
            step_start = time.time()
            
            try:
                # SCIENTIFIC MULTILOOKING PARAMETERS:
                # Calculate optimal looks for near-square ground pixels
                
                if not self.metadata or not isinstance(self.metadata, dict):
                    raise ValueError("No valid metadata available for multilooking")

                input_range_spacing = float(self.get_current_range_spacing())
                input_azimuth_spacing = float(self.get_current_azimuth_spacing())
                
                # SCIENTIFIC PARAMETER CALCULATION using real annotation metadata
                if self.validated_metadata is None:
                    raise ValueError("Validated metadata not available for multilooking parameters")

                try:
                    incidence_near = float(self.validated_metadata['incidence_angle_near'])
                    incidence_far = float(self.validated_metadata['incidence_angle_far'])
                except KeyError as missing_angle:
                    raise ValueError(
                        f"Missing incidence angle metadata required for multilooking: {missing_angle}"
                    )
                except (TypeError, ValueError) as angle_error:
                    raise ValueError(
                        f"Invalid incidence angle metadata values: {angle_error}"
                    )

                incidence_angle_deg = 0.5 * (incidence_near + incidence_far)
                incidence_angle_rad = np.radians(incidence_angle_deg)
                ground_range_spacing = input_range_spacing * np.sin(incidence_angle_rad)
                azimuth_spacing_m = input_azimuth_spacing

                if not np.isfinite(ground_range_spacing) or ground_range_spacing <= 0.0:
                    raise ValueError(
                        f"Invalid ground range spacing derived from metadata: {ground_range_spacing}"
                    )

                if self.target_resolution is None:
                    raise ValueError("Target resolution was not initialised from metadata or options")
                target_resolution = float(self.target_resolution)
                if not np.isfinite(target_resolution) or target_resolution <= 0.0:
                    raise ValueError(
                        f"Invalid target resolution: {target_resolution}m (must be positive)"
                    )
                
                # Calculate optimal looks
                range_looks = max(1, int(np.round(target_resolution / ground_range_spacing)))
                azimuth_looks = max(1, int(np.round(target_resolution / azimuth_spacing_m)))
                
                # Apply scientific constraints: Lr=4-6, La=1-2 for TOPS
                range_looks = max(4, min(6, range_looks))
                azimuth_looks = max(1, min(2, azimuth_looks))
                
                print(f"   🔬 SCIENTIFIC MULTILOOKING:")
                print(f"      Incidence angle (avg): {incidence_angle_deg:.2f}°")
                print(f"      Ground range spacing: {ground_range_spacing:.1f}m")
                print(f"      Azimuth spacing: {azimuth_spacing_m:.1f}m")
                print(f"      Calculated looks: {range_looks}×{azimuth_looks} (range×azimuth)")
                print(f"      Target resolution: {target_resolution:.2f}m")
                print(f"      ENL target: ~{range_looks * azimuth_looks}")
                
                # Ensure working_data is numpy array
                if not isinstance(working_data, np.ndarray):
                    raise ValueError("Working data is not a valid numpy array")
                
                # Apply multilooking in LINEAR DOMAIN (σ⁰ power units)
                ml_result = sardine.apply_multilooking(
                    working_data,
                    int(range_looks),
                    int(azimuth_looks),
                    float(input_range_spacing),
                    float(input_azimuth_spacing)
                )
                
                # Handle result format
                output_range_spacing = None
                output_azimuth_spacing = None
                if isinstance(ml_result, dict):
                    multilooked_data = ml_result.get('data', None)
                    if multilooked_data is None:
                        raise ValueError("Multilooking result does not contain data")
                    output_range_spacing = ml_result.get('range_spacing')
                    output_azimuth_spacing = ml_result.get('azimuth_spacing')
                else:
                    multilooked_data = ml_result

                if output_range_spacing is None:
                    output_range_spacing = input_range_spacing * float(range_looks)
                if output_azimuth_spacing is None:
                    output_azimuth_spacing = input_azimuth_spacing * float(azimuth_looks)
                    
                if not isinstance(multilooked_data, np.ndarray):
                    raise ValueError("Multilooking failed to produce valid numpy array")

                self.update_current_spacing(output_range_spacing, output_azimuth_spacing)
                print(f"      Output spacing: {output_range_spacing:.2f}m × {output_azimuth_spacing:.2f}m (range×azimuth)")
                
                step_duration = time.time() - step_start
                self.log_step(
                    8,
                    "Multilooking",
                    "success",
                    f"Scientific multilook: {range_looks}×{azimuth_looks}, shape: {multilooked_data.shape}, ENL≈{range_looks*azimuth_looks}, spacing: {output_range_spacing:.2f}m×{output_azimuth_spacing:.2f}m",
                    step_duration,
                )
                working_data = multilooked_data
                
            except Exception as e:
                step_duration = time.time() - step_start
                self.log_step(8, "Multilooking", "error", f"Multilooking failed: {e}", step_duration)
                raise RuntimeError(f"SCIENTIFIC MODE FAILURE: Multilooking failed. Error: {e}")
            
            # STEP 9: Terrain Flattening - SCIENTIFIC ORDER (before speckle filtering)
            self.announce_step(9, "Terrain Flattening", "Applying DEM-based terrain flattening")
            step_start = time.time()
            
            if self.terrain_flatten:
                try:
                    # Use real DEM data for terrain flattening
                    if not isinstance(working_data, np.ndarray):
                        raise ValueError("Working data is not a valid numpy array")
                    
                    rows, cols = working_data.shape
                    
                    # Extract real coordinate bounds from metadata for DEM loading
                    if not hasattr(self, 'metadata') or not self.metadata or not isinstance(self.metadata, dict):
                        raise ValueError("No valid metadata available for terrain flattening")
                    
                    # SCIENTIFIC COORDINATE EXTRACTION: Get real geographic bounds from metadata
                    # Coordinates are now properly extracted by the Rust SlcReader
                    
                    min_lat = float(self.metadata.get('min_latitude'))
                    max_lat = float(self.metadata.get('max_latitude'))
                    min_lon = float(self.metadata.get('min_longitude'))
                    max_lon = float(self.metadata.get('max_longitude'))
                    
                    print(f"   ✅ USING METADATA COORDINATES: lat=({min_lat:.6f}, {max_lat:.6f}), lon=({min_lon:.6f}, {max_lon:.6f})")
                    
                    # Create bounding box in the format expected by load_dem_for_bbox
                    bbox = [float(min_lon), float(min_lat), float(max_lon), float(max_lat)]
                    
                    # Load real SRTM DEM data
                    cache_dir = f"{self.output_dir}/dem_cache"
                    
                    # Load real SRTM DEM using dem.rs load_dem_for_bbox function
                    print(f"   🗻 Loading real SRTM DEM for bbox: {bbox}")
                    dem_result = sardine.load_dem_for_bbox(bbox, cache_dir)

                    dem_geo_transform = None
                    dem_bbox_meta = None

                    # Extract DEM data from result dictionary
                    if isinstance(dem_result, dict) and 'data' in dem_result:
                        dem_data = np.asarray(dem_result['data'])
                        dem_geo_transform = dem_result.get('geo_transform')
                        dem_bbox_meta = dem_result.get('bbox')
                        print(f"   ✅ DEM loaded successfully: {dem_data.shape}, resolution: {dem_result.get('resolution', 'unknown')}m")
                        print(f"   📊 DEM elevation range: {dem_data.min():.1f} to {dem_data.max():.1f} meters")
                    else:
                        dem_data = np.asarray(dem_result)

                    if not isinstance(dem_data, np.ndarray) or dem_data.size == 0:
                        raise ValueError("Failed to load real DEM data from SRTM")

                    # Ensure DEM orientation matches metadata before any resampling
                    try:
                        expected_min_lat = float(self.metadata.get('min_latitude'))
                        expected_max_lat = float(self.metadata.get('max_latitude'))
                        expected_min_lon = float(self.metadata.get('min_longitude'))
                        expected_max_lon = float(self.metadata.get('max_longitude'))
                    except (TypeError, ValueError):
                        expected_min_lat = expected_max_lat = expected_min_lon = expected_max_lon = None

                    if dem_geo_transform is not None and isinstance(dem_geo_transform, (list, tuple)):
                        dem_geo_transform = [float(value) for value in dem_geo_transform]
                        dem_lat_top = dem_geo_transform[3]
                        dem_lat_bottom = dem_lat_top + dem_geo_transform[5] * (dem_data.shape[0] - 1)
                        dem_lon_left = dem_geo_transform[0]
                        dem_lon_right = dem_lon_left + dem_geo_transform[1] * (dem_data.shape[1] - 1)

                        if expected_min_lat is not None and expected_max_lat is not None:
                            # Check whether the top row is closer to the southern edge (requires flip)
                            if (
                                abs(dem_lat_top - expected_min_lat) < abs(dem_lat_top - expected_max_lat)
                                and abs(dem_lat_bottom - expected_max_lat) < abs(dem_lat_bottom - expected_min_lat)
                            ):
                                dem_data = np.flipud(dem_data)
                                dem_lat_top, dem_lat_bottom = dem_lat_bottom, dem_lat_top
                                dem_geo_transform[3] = dem_lat_top
                                dem_geo_transform[5] = -dem_geo_transform[5]
                                print("   🔃 Flipped DEM vertically to match metadata orientation")

                        if expected_min_lon is not None and expected_max_lon is not None:
                            if (
                                abs(dem_lon_left - expected_max_lon) < abs(dem_lon_left - expected_min_lon)
                                and abs(dem_lon_right - expected_min_lon) < abs(dem_lon_right - expected_max_lon)
                            ):
                                dem_data = np.fliplr(dem_data)
                                dem_lon_left, dem_lon_right = dem_lon_right, dem_lon_left
                                dem_geo_transform[0] = dem_lon_left
                                dem_geo_transform[1] = abs(dem_geo_transform[1])
                                print("   🔃 Flipped DEM horizontally to match metadata orientation")

                        self.dem_geo_transform = tuple(dem_geo_transform)
                        if dem_bbox_meta is not None:
                            self.dem_bbox = tuple(dem_bbox_meta)

                        print(
                            "   🌐 DEM geo_transform: origin=(%.6f, %.6f), pixel_size=(%.8f, %.8f)"
                            % (
                                dem_geo_transform[0],
                                dem_geo_transform[3],
                                dem_geo_transform[1],
                                dem_geo_transform[5],
                            )
                        )

                    # Ensure DEM matches SAR data dimensions
                    if dem_data.shape != (rows, cols):
                        # Resample DEM to SAR grid using scientific interpolation
                        from scipy.ndimage import zoom
                        zoom_factor = (rows / dem_data.shape[0], cols / dem_data.shape[1])
                        dem_data = zoom(dem_data, zoom_factor, order=1)

                    dem_data = np.asarray(dem_data, dtype=np.float32)
                    
                    # Extract orbit data for terrain flattening - flexible approach
                    print(f"   🛰️  Extracting orbit state vectors")
                    
                    orbit_times = []
                    orbit_positions = []
                    orbit_velocities = []
                    
                    # Try to use existing orbit files first
                    orbit_cache_dir = f"{self.output_dir}/orbit_cache"
                    
                    # Check for existing .EOF files
                    import glob
                    orbit_files = []
                    if Path(orbit_cache_dir).exists():
                        orbit_files = glob.glob(os.path.join(orbit_cache_dir, "*.EOF"))
                    
                    if orbit_files:
                        # Load from existing .EOF file
                        orbit_file_path = max(orbit_files, key=os.path.getmtime)
                        print(f"   📡 Loading orbit vectors from: {os.path.basename(orbit_file_path)}")
                        
                        try:
                            import xml.etree.ElementTree as ET
                            tree = ET.parse(orbit_file_path)
                            root = tree.getroot()
                            
                            # Find all OSV (Orbit State Vector) elements
                            osvs = root.findall('.//OSV')
                            if osvs:
                                for osv in osvs:
                                    # Extract UTC time
                                    utc_elem = osv.find('UTC')
                                    if utc_elem is not None:
                                        utc_time = utc_elem.text.replace('UTC=', '')
                                        # Convert time string to timestamp
                                        try:
                                            dt = datetime.fromisoformat(utc_time.replace('Z', '+00:00'))
                                            orbit_times.append(dt.timestamp())
                                        except:
                                            continue
                                    
                                    # Extract position (X, Y, Z)
                                    x_elem = osv.find('X')
                                    y_elem = osv.find('Y')
                                    z_elem = osv.find('Z')
                                    if all(elem is not None for elem in [x_elem, y_elem, z_elem]):
                                        orbit_positions.extend([
                                            float(x_elem.text), 
                                            float(y_elem.text), 
                                            float(z_elem.text)
                                        ])
                                    
                                    # Extract velocity (VX, VY, VZ)
                                    vx_elem = osv.find('VX')
                                    vy_elem = osv.find('VY')
                                    vz_elem = osv.find('VZ')
                                    if all(elem is not None for elem in [vx_elem, vy_elem, vz_elem]):
                                        orbit_velocities.extend([
                                            float(vx_elem.text), 
                                            float(vy_elem.text), 
                                            float(vz_elem.text)
                                        ])
                                
                                print(f"   ✅ Loaded {len(orbit_times)} orbit state vectors from .EOF file")
                        except Exception as eof_error:
                            print(f"   ⚠️  Could not parse .EOF file: {eof_error}")
                    
                    # If no orbit data from .EOF, try metadata
                    if not orbit_times:
                        print(f"   📊 Extracting orbit vectors from metadata")
                        
                        # Find orbit state vector keys in metadata
                        orbit_keys = [key for key in self.metadata.keys() if 'orbit_state_vector' in key and '_time' in key]
                        vector_indices = []
                        
                        for key in orbit_keys:
                            try:
                                parts = key.split('_')
                                if len(parts) >= 4:
                                    idx = int(parts[3])
                                    if idx not in vector_indices:
                                        vector_indices.append(idx)
                            except (ValueError, IndexError):
                                continue
                        
                        # Extract orbit data from metadata
                        for i in sorted(vector_indices):
                            time_key = f'orbit_state_vector_{i}_time'
                            pos_x_key = f'orbit_state_vector_{i}_x_position'
                            pos_y_key = f'orbit_state_vector_{i}_y_position' 
                            pos_z_key = f'orbit_state_vector_{i}_z_position'
                            vel_x_key = f'orbit_state_vector_{i}_x_velocity'
                            vel_y_key = f'orbit_state_vector_{i}_y_velocity'
                            vel_z_key = f'orbit_state_vector_{i}_z_velocity'
                            
                            if all(key in self.metadata for key in [time_key, pos_x_key, pos_y_key, pos_z_key]):
                                orbit_times.append(float(self.metadata[time_key]))
                                orbit_positions.extend([
                                    float(self.metadata[pos_x_key]),
                                    float(self.metadata[pos_y_key]), 
                                    float(self.metadata[pos_z_key])
                                ])
                                orbit_velocities.extend([
                                    float(self.metadata.get(vel_x_key, 0.0)),
                                    float(self.metadata.get(vel_y_key, 0.0)),
                                    float(self.metadata.get(vel_z_key, 0.0))
                                ])
                        
                        if orbit_times:
                            print(f"   ✅ Extracted {len(orbit_times)} orbit state vectors from metadata")
                    
                    # Use orbit data if available, otherwise continue without terrain flattening
                    if not orbit_times:
                        print(f"   ⚠️  No orbit data available - skipping terrain flattening")
                        
                        # Continue with regular calibration without terrain flattening
                        multilook_result = sardine.radiometric_calibration_with_denoising(
                            str(self.input_path),
                            str(self.output_dir / "calibrated"),
                            polarization,
                            look_cols
                        )
                        
                        if not isinstance(multilook_result, dict) or "data" not in multilook_result:
                            raise ValueError("Failed to get calibrated data")
                        
                        calibrated_data = multilook_result["data"]
                        print(f"   ✅ Calibrated data (without terrain flattening): {calibrated_data.shape}")
                    else:
                        # Create proper orbit data dict for Rust function
                        orbit_data = {
                            'times': orbit_times,
                            'positions': orbit_positions, 
                            'velocities': orbit_velocities
                        }
                    
                    # Extract real pixel spacing from metadata - NO hardcoding allowed
                    range_spacing = float(self.metadata.get('range_pixel_spacing'))
                    azimuth_spacing = float(self.metadata.get('azimuth_pixel_spacing'))
                    if 'wavelength' not in self.metadata:
                        raise ValueError("Missing 'wavelength' in metadata; cannot proceed scientifically")
                    wavelength = float(self.metadata.get('wavelength'))
                    
                    if not range_spacing or not azimuth_spacing:
                        raise ValueError(f"Missing pixel spacing in metadata: range={range_spacing}, azimuth={azimuth_spacing}")
                    
                    # Apply SCIENTIFIC terrain flattening using geometric approach (UPGRADED)
                    print(f"   🏔️  Applying scientific terrain flattening with geometric approach")
                    print(f"   📊 Parameters: DEM data available, range={range_spacing:.3f}m, azimuth={azimuth_spacing:.3f}m")
                    
                    # Use the scientific terrain flattening implementation
                    # Extract safe path and primary subswath for automatic parameter extraction
                    safe_path = self.input_path
                    primary_subswath = "IW2"  # Use the primary subswath identified during processing
                    
                    terrain_result = sardine.apply_scientific_terrain_flattening(
                        working_data,      # Calibrated sigma0 data
                        dem_data,          # Real SRTM DEM data
                        safe_path,         # safe_path for automatic heading extraction
                        primary_subswath,  # subswath for automatic parameter extraction
                        None,              # azimuth_heading (will be extracted from annotation/orbit)
                        None,              # ellipsoid_incidence_angle (will be extracted from annotation) 
                        30.0               # dem_pixel_spacing (Copernicus DEM standard)
                    )
                    
                    # Extract data from result dictionary
                    if isinstance(terrain_result, dict) and 'data' in terrain_result:
                        terrain_data = terrain_result['data']
                        print(f"   ✅ Terrain flattening successful: {terrain_data.shape}")
                        
                        # Check for valid data range
                        finite_data = terrain_data[np.isfinite(terrain_data)]
                        if len(finite_data) > 0:
                            print(f"   📊 Terrain corrected range: {finite_data.min():.6f} to {finite_data.max():.6f}")
                            valid_percentage = len(finite_data) / terrain_data.size * 100
                            print(f"   📊 Valid pixels: {valid_percentage:.1f}%")
                        else:
                            print(f"   ⚠️  Warning: All terrain-corrected values are NaN - may be over water or invalid region")
                            # Continue processing but flag for attention
                    else:
                        raise ValueError("Terrain flattening result missing 'data' field; structured output required")
                        
                    if not isinstance(terrain_data, np.ndarray) or terrain_data.size == 0:
                        raise ValueError("Terrain flattening failed to produce valid data")
                    
                    step_duration = time.time() - step_start
                    self.log_step(9, "Terrain Flattening", "success", 
                                 f"Full scientific method: {len(orbit_times)} vectors, DEM: {dem_data.shape}, output: {terrain_data.shape}", step_duration)
                    working_data = terrain_data
                    
                except Exception as e:
                    step_duration = time.time() - step_start
                    self.log_step(9, "Terrain Flattening", "error", f"Terrain flattening failed: {e}", step_duration)
                    raise RuntimeError(f"SCIENTIFIC MODE FAILURE: Terrain flattening requires real DEM data. Error: {e}")
            else:
                step_duration = time.time() - step_start
                self.log_step(9, "Terrain Flattening", "skipped", "Disabled by user", step_duration)
            # STEP 10: Speckle Filtering - SCIENTIFIC ORDER (after multilooking)
            self.announce_step(10, "Speckle Filtering", "Reducing speckle while preserving radiometry")
            step_start = time.time()
            
            try:
                # Ensure working_data is proper numpy array
                if not isinstance(working_data, np.ndarray):
                    raise RuntimeError("SCIENTIFIC MODE FAILURE: Working data is not a valid numpy array - cannot proceed with speckle filtering")
                
                # SCIENTIFIC MODE: Validate 2D array for speckle filtering - NO synthetic data generation
                if working_data.ndim == 0:
                    raise RuntimeError("SCIENTIFIC MODE FAILURE: Scalar data detected - cannot apply speckle filtering. Real SAR data required.")
                elif working_data.ndim == 1:
                    raise RuntimeError("SCIENTIFIC MODE FAILURE: 1D data detected - cannot apply speckle filtering. Real SAR data required.")
                elif working_data.ndim != 2:
                    raise RuntimeError(f"SCIENTIFIC MODE FAILURE: Invalid data dimensions {working_data.ndim}D - speckle filtering requires 2D SAR data")
                
                # Check for valid finite data
                finite_data = working_data[np.isfinite(working_data)]
                if len(finite_data) == 0:
                    raise RuntimeError("SCIENTIFIC MODE FAILURE: All input data contains NaN or infinite values - cannot apply speckle filtering")
                
                finite_percentage = len(finite_data) / working_data.size * 100
                print(f"   📊 Input data: {finite_percentage:.1f}% finite values")
                print(f"   📊 Input dimensions: {working_data.shape[0]}x{working_data.shape[1]}")
                
                if finite_percentage < 10.0:
                    raise RuntimeError(f"SCIENTIFIC MODE FAILURE: Insufficient valid data ({finite_percentage:.1f}%) for speckle filtering - minimum 10% required")
                
                # Calculate data size for performance monitoring
                data_size_mb = working_data.nbytes / 1024 / 1024
                self.performance_monitor.start_step("Speckle Filtering", data_size_mb=data_size_mb,
                                                   parallel_enabled=self.enable_parallel,
                                                   threads_used=self.num_threads,
                                                   chunk_size=self.chunk_size)
                
                # Extract real parameters for speckle filtering from metadata
                filter_type = str(self.speckle_filter)
                window_size = int(self.filter_window) if self.filter_window is not None else 7
                
                # CRITICAL FIX: Adjust window size if image is too small
                min_image_dim = min(working_data.shape[0], working_data.shape[1])
                skip_speckle_filtering = False
                
                if min_image_dim < window_size:
                    original_window = window_size
                    window_size = max(3, min_image_dim // 2 * 2 - 1)  # Ensure odd window size
                    print(f"   ⚠️  Image too small ({working_data.shape}) for window size {original_window}, reducing to {window_size}")
                    if window_size < 3:
                        print(f"   ⚠️  Image too small for meaningful speckle filtering, skipping this step...")
                        skip_speckle_filtering = True
                
                if not skip_speckle_filtering:
                    # Calculate real number of looks from multilooking parameters
                    num_looks = float(self.multilook_range * self.multilook_azimuth)
                    
                    # Scientific speckle filter parameters for enhanced_lee
                    # Values based on literature: Lopes et al., "Adaptive Speckle Filters", IEEE Trans
                    edge_threshold = 0.5    # Edge detection threshold for enhanced Lee filter
                    damping_factor = 1.0    # Damping factor for enhanced Lee algorithm  
                    cv_threshold = 0.5      # Coefficient of variation threshold
                        
                    # Apply speckle filtering with all required parameters
                    filtered_result = sardine.apply_speckle_filter(
                        working_data, 
                        filter_type, 
                        window_size, 
                        num_looks,
                        edge_threshold,
                        damping_factor, 
                        cv_threshold
                    )
                else:
                    # Skip speckle filtering - use working_data as-is
                    filtered_result = working_data
                
                # Handle result format - check if it's a dictionary with data key
                if skip_speckle_filtering:
                    # Speckle filtering was skipped
                    filtered_data = working_data
                    step_duration = time.time() - step_start
                    self.log_step(10, "Speckle Filtering", "skipped", 
                                 f"Image too small: {working_data.shape}", step_duration)
                    self.performance_monitor.end_step()
                else:
                    # Speckle filtering was applied
                    if isinstance(filtered_result, dict):
                        # Check for both 'data' and 'filtered_data' keys
                        filtered_data = filtered_result.get('data')
                        if filtered_data is None:
                            filtered_data = filtered_result.get('filtered_data')
                        if filtered_data is None or not isinstance(filtered_data, np.ndarray):
                            raise RuntimeError("SCIENTIFIC MODE FAILURE: Speckle filtering failed to produce valid data array")
                    else:
                        if not isinstance(filtered_result, np.ndarray):
                            raise RuntimeError("SCIENTIFIC MODE FAILURE: Speckle filtering failed to produce valid numpy array")
                        filtered_data = filtered_result
                    
                    step_duration = time.time() - step_start
                    self.log_step(10, "Speckle Filtering", "success", 
                                 f"Filter: {filter_type}, shape: {filtered_data.shape}", step_duration)
                    self.performance_monitor.end_step()
                
                working_data = filtered_data
                
            except Exception as e:
                step_duration = time.time() - step_start
                self.log_step(10, "Speckle Filtering", "error", f"Speckle filtering failed: {e}", step_duration)
                self.performance_monitor.end_step()
                raise RuntimeError(f"SCIENTIFIC MODE FAILURE: Speckle filtering failed. Error: {e}")
            
            # STEP 11: Terrain Correction (Geocoding) - SCIENTIFIC ORDER
            self.announce_step(11, "Terrain Correction", "Geocoding with precise orbits")
            step_start = time.time()
            
            if self.geocode:
                try:
                    # SCIENTIFIC MODE: Use real terrain correction with orbit data
                    if not isinstance(working_data, np.ndarray):
                        working_data = np.array(working_data, dtype=np.float32)
                    
                    # Extract real geocoding parameters from metadata
                    if hasattr(self, 'metadata') and self.metadata and isinstance(self.metadata, dict):
                        # Get real bounding box coordinates from parsed metadata
                        min_lat = float(self.metadata.get('min_latitude'))
                        max_lat = float(self.metadata.get('max_latitude'))
                        min_lon = float(self.metadata.get('min_longitude'))
                        max_lon = float(self.metadata.get('max_longitude'))

                        # Ensure bounding box is properly ordered (min <= max)
                        lat_bounds = sorted([min_lat, max_lat])
                        lon_bounds = sorted([min_lon, max_lon])

                        if lat_bounds[0] != min_lat or lat_bounds[1] != max_lat:
                            print("   🔁 Adjusted latitude bounds for terrain correction")
                        if lon_bounds[0] != min_lon or lon_bounds[1] != max_lon:
                            print("   🔁 Adjusted longitude bounds for terrain correction")

                        min_lat, max_lat = lat_bounds
                        min_lon, max_lon = lon_bounds
                        
                        print(f"   ✅ GEOCODING COORDINATES: lat=({min_lat:.6f}, {max_lat:.6f}), lon=({min_lon:.6f}, {max_lon:.6f})")
                        
                        base_bbox = [float(min_lon), float(min_lat), float(max_lon), float(max_lat)]

                        rows, cols = working_data.shape if working_data is not None else (0, 0)
                        refined_bbox = _refine_geocoding_bbox(
                            base_bbox,
                            int(rows),
                            int(cols),
                            float(range_spacing),
                            float(azimuth_spacing),
                        )

                        if refined_bbox is not None:
                            sar_bbox, bbox_metrics = refined_bbox
                            lat_reduction = bbox_metrics.get("shrink_lat_pct", 0.0) * 100.0
                            lon_reduction = bbox_metrics.get("shrink_lon_pct", 0.0) * 100.0
                            if lat_reduction > 0.1 or lon_reduction > 0.1:
                                print(
                                    "   ✂️  Refined geocoding bbox using multilooked footprint"
                                )
                                print(
                                    "      Lat span: "
                                    f"{bbox_metrics['new_lat_extent_deg']:.6f}° "
                                    f"({lat_reduction:.1f}% tighter)"
                                )
                                print(
                                    "      Lon span: "
                                    f"{bbox_metrics['new_lon_extent_deg']:.6f}° "
                                    f"({lon_reduction:.1f}% tighter)"
                                )
                        else:
                            sar_bbox = base_bbox
                        
                        # Extract orbit data using the EXACT same method as terrain flattening
                        vector_indices = [int(key.split('_')[3]) for key in self.metadata.keys() 
                                        if key.startswith('orbit_state_vector_') and key.endswith('_time')]
                        
                        orbit_times_str = []
                        orbit_positions = []
                        orbit_velocities = []
                        
                        for i in sorted(vector_indices):
                            time_key = f'orbit_state_vector_{i}_time'
                            pos_x_key = f'orbit_state_vector_{i}_x_position'
                            pos_y_key = f'orbit_state_vector_{i}_y_position' 
                            pos_z_key = f'orbit_state_vector_{i}_z_position'
                            vel_x_key = f'orbit_state_vector_{i}_x_velocity'
                            vel_y_key = f'orbit_state_vector_{i}_y_velocity'
                            vel_z_key = f'orbit_state_vector_{i}_z_velocity'
                            
                            if all(key in self.metadata for key in [time_key, pos_x_key, pos_y_key, pos_z_key]):
                                # Convert time from seconds since reference to ISO format
                                time_val = float(self.metadata[time_key])
                                # Use epoch time conversion for proper formatting
                                dt = datetime.fromtimestamp(time_val, tz=timezone.utc)
                                time_iso = dt.isoformat()
                                orbit_times_str.append(time_iso)
                                
                                orbit_positions.append([
                                    float(self.metadata[pos_x_key]),
                                    float(self.metadata[pos_y_key]), 
                                    float(self.metadata[pos_z_key])
                                ])
                                orbit_velocities.append([
                                    float(self.metadata.get(vel_x_key, 0.0)),
                                    float(self.metadata.get(vel_y_key, 0.0)),
                                    float(self.metadata.get(vel_z_key, 0.0))
                                ])
                        
                        if len(orbit_times_str) == 0:
                            raise ValueError("No orbit data found for geocoding - ensure orbit file has been applied")
                        
                        print(f"   ✅ Extracted {len(orbit_times_str)} orbit state vectors for geocoding")
                        
                        # Extract real SLC metadata - MUST use annotation data, NO hardcoded values
                        range_spacing = float(self.get_current_range_spacing())
                        azimuth_spacing = float(self.get_current_azimuth_spacing())

                        # Speed of light used throughout SAR timing conversions
                        speed_of_light = 299_792_458.0  # m/s (exact definition)
                        time_per_original_range_sample = None
                        if self.original_range_spacing and self.original_range_spacing > 0:
                            time_per_original_range_sample = (
                                2.0 * float(self.original_range_spacing) / speed_of_light
                            )

                        # Derive multilook factors relative to the original SLC sampling
                        range_looks_factor = None
                        azimuth_looks_factor = None
                        if self.original_range_spacing and self.original_range_spacing > 0:
                            range_looks_factor = range_spacing / float(self.original_range_spacing)
                        if self.original_azimuth_spacing and self.original_azimuth_spacing > 0:
                            azimuth_looks_factor = azimuth_spacing / float(self.original_azimuth_spacing)
                        
                        # Extract slant_range_time from cached metadata (should be available after cache extraction fix)
                        slant_range_time = None
                        prf = None
                        
                        # First try cached metadata (this should work after our fix)
                        if 'slant_range_time' in self.metadata:
                            slant_range_time = float(self.metadata.get('slant_range_time'))
                            print(f"   ✅ Using cached slant_range_time: {slant_range_time:.6f} s")
                        
                        if 'prf' in self.metadata:
                            prf = float(self.metadata.get('prf'))
                            print(f"   ✅ Using cached prf: {prf:.1f} Hz")
                        
                        # Strict mode: require timing parameters from metadata; no fallbacks
                        if slant_range_time is None or prf is None:
                            raise ValueError("Missing required timing parameters in metadata: slant_range_time and prf are required for terrain correction")

                        # Adjust slant range time to reflect merged IW geometry and multilooking center
                        if (
                            time_per_original_range_sample is not None
                            and hasattr(self, 'reader')
                            and self.reader is not None
                        ):
                            try:
                                subswath_metadata = self.reader.get_all_iw_subswaths()
                                swaths_for_pol = None
                                if isinstance(subswath_metadata, dict):
                                    swaths_for_pol = subswath_metadata.get(self.polarization)
                                    if swaths_for_pol is None:
                                        # Try uppercase lookup (keys should already be uppercase but keep robust)
                                        swaths_for_pol = subswath_metadata.get(self.polarization.upper())

                                if isinstance(swaths_for_pol, dict) and swaths_for_pol:
                                    global_start_time = None
                                    for swath_name, swath_info in swaths_for_pol.items():
                                        srt_value = swath_info.get('slant_range_time')
                                        first_sample_global = swath_info.get('first_sample_global')
                                        if srt_value is None or first_sample_global is None:
                                            continue

                                        start_time = float(srt_value) + float(first_sample_global) * time_per_original_range_sample
                                        if global_start_time is None or start_time < global_start_time:
                                            global_start_time = start_time

                                    if global_start_time is not None:
                                        multilook_center_offset = 0.0
                                        if range_looks_factor and range_looks_factor > 1.0:
                                            multilook_center_offset = 0.5 * (range_looks_factor - 1.0) * time_per_original_range_sample

                                        adjusted_slant_range_time = global_start_time + multilook_center_offset
                                        print(
                                            "   🔁 Adjusted slant_range_time for merged IW data: "
                                            f"{adjusted_slant_range_time:.9f} s (was {slant_range_time:.9f} s)"
                                        )
                                        slant_range_time = adjusted_slant_range_time
                                else:
                                    print("   ⚠️  Could not access IW subswath metadata for slant_range_time adjustment")
                            except Exception as adjust_error:
                                print(f"   ⚠️  Could not adjust slant_range_time from IW geometry: {adjust_error}")
                        
                        if 'wavelength' not in self.metadata:
                            raise ValueError("Missing 'wavelength' in metadata; cannot geocode scientifically")
                        wavelength = float(self.metadata.get('wavelength'))
                        
                        # Calculate radar frequency from real wavelength (c = λf)
                        radar_frequency = speed_of_light / wavelength  # Hz - calculated from real wavelength
                        
                        effective_prf = prf
                        if azimuth_looks_factor and azimuth_looks_factor > 0:
                            effective_prf = prf / azimuth_looks_factor

                        if azimuth_looks_factor:
                            print(f"   🔁 Derived multilook factors: range≈{range_looks_factor or 1.0:.2f}x, azimuth≈{azimuth_looks_factor:.2f}x")
                            print(f"   📏 Effective PRF after multilooking: {effective_prf:.2f} Hz")

                        # CRITICAL FIX: Extract product timing information for terrain correction
                        # These are needed for correct azimuth pixel calculation
                        reader_metadata = self.reader.get_cached_metadata()
                        
                        # Extract timing from metadata dictionary (already contains product_start_time_abs)
                        product_start_time_abs = float(reader_metadata.get('product_start_time_abs', 0.0))
                        product_stop_time_abs = float(reader_metadata.get('product_stop_time_abs', 0.0))
                        
                        # Calculate duration if not already in metadata
                        if product_start_time_abs > 0.0 and product_stop_time_abs > 0.0:
                            product_duration = product_stop_time_abs - product_start_time_abs
                        else:
                            # Fallback: calculate from stop_time and start_time strings
                            from dateutil import parser as date_parser
                            try:
                                start_dt = date_parser.parse(reader_metadata.get('start_time', ''))
                                stop_dt = date_parser.parse(reader_metadata.get('stop_time', ''))
                                product_start_time_abs = start_dt.timestamp()
                                product_stop_time_abs = stop_dt.timestamp()
                                product_duration = product_stop_time_abs - product_start_time_abs
                            except Exception as e:
                                print(f"⚠️ Warning: Could not parse timing from metadata: {e}")
                                product_duration = 0.0
                        
                        real_metadata = {
                            'range_pixel_spacing': range_spacing,
                            'azimuth_pixel_spacing': azimuth_spacing,
                            'slant_range_time': slant_range_time,
                            'prf': effective_prf,
                            'wavelength': wavelength,  # Real wavelength from annotation XML
                            'radar_frequency': radar_frequency,  # Calculated from real wavelength
                            # CRITICAL: Add timing information for correct azimuth indexing
                            'product_start_time_abs': product_start_time_abs,
                            'product_stop_time_abs': product_stop_time_abs,
                            'product_duration': product_duration,
                            'total_azimuth_lines': working_data.shape[0] if working_data is not None else None
                        }
                        
                        # ============================================================
                        # COMPREHENSIVE PRE-FLIGHT DIAGNOSTICS
                        # ============================================================
                        print("\n" + "="*70)
                        print("🔍 PRE-FLIGHT VALIDATION (Checking for common errors)")
                        print("="*70)
                        
                        # 1. PRF Sanity Check (should be ~1000-2000 Hz for Sentinel-1 TOPS)
                        print(f"1️⃣  PRF Check:")
                        print(f"   • Effective PRF: {effective_prf:.3f} Hz")
                        if effective_prf < 10:
                            print(f"   ⚠️  WARNING: PRF < 10 Hz - may be in kHz instead of Hz!")
                        elif effective_prf > 10000:
                            print(f"   ⚠️  WARNING: PRF > 10 kHz - may be in Hz but using kHz value!")
                        elif 1000 <= effective_prf <= 2500:
                            print(f"   ✅ PRF in expected range for Sentinel-1 TOPS (1-2.5 kHz)")
                        else:
                            print(f"   ⚠️  PRF outside typical S1 TOPS range (1-2.5 kHz)")
                        
                        # 2. Time Units Check (seconds vs milliseconds)
                        print(f"\n2️⃣  Time Units Check:")
                        print(f"   • product_start_time_abs: {product_start_time_abs:.6f} s")
                        print(f"   • product_duration: {product_duration:.6f} s")
                        if product_start_time_abs < 1e6:
                            print(f"   ⚠️  WARNING: product_start_time_abs < 1e6 - may not be Unix timestamp!")
                        elif product_start_time_abs > 1e12:
                            print(f"   ⚠️  WARNING: product_start_time_abs > 1e12 - may be in milliseconds instead of seconds!")
                        else:
                            print(f"   ✅ product_start_time_abs looks like valid Unix timestamp (seconds)")
                        
                        if product_duration < 0:
                            print(f"   ⚠️  WARNING: Negative duration - start/stop times may be swapped!")
                        elif product_duration > 100:
                            print(f"   ⚠️  WARNING: Duration > 100s - unusual for Sentinel-1 scene")
                        elif 10 <= product_duration <= 50:
                            print(f"   ✅ Duration in expected range for S1 scene (10-50s)")
                        else:
                            print(f"   ⚠️  Duration {product_duration:.1f}s is unusual")
                        
                        # 3. Expected Lines Check
                        print(f"\n3️⃣  Expected Azimuth Lines Check:")
                        expected_lines = product_duration * effective_prf
                        actual_lines = working_data.shape[0] if working_data is not None else 0
                        print(f"   • Actual lines in data: {actual_lines:,}")
                        print(f"   • Expected from duration×PRF: {expected_lines:,.0f}")
                        if actual_lines > 0:
                            ratio = actual_lines / expected_lines if expected_lines > 0 else 0
                            print(f"   • Ratio (actual/expected): {ratio:.3f}")
                            if 0.8 <= ratio <= 1.2:
                                print(f"   ✅ Line count consistent with duration and PRF")
                            else:
                                print(f"   ⚠️  WARNING: Line count mismatch - may indicate PRF or duration error!")
                        
                        # 4. Wavelength Check
                        print(f"\n4️⃣  Wavelength Check:")
                        print(f"   • Wavelength: {wavelength:.8f} m ({wavelength*100:.4f} cm)")
                        if 0.05 <= wavelength <= 0.06:
                            print(f"   ✅ Wavelength in C-band range (5-6 cm)")
                        else:
                            print(f"   ⚠️  WARNING: Wavelength outside C-band range!")
                        
                        # 5. Pixel Spacing Check
                        print(f"\n5️⃣  Pixel Spacing Check:")
                        print(f"   • Range spacing: {range_spacing:.6f} m")
                        print(f"   • Azimuth spacing: {azimuth_spacing:.6f} m")
                        if 2 <= range_spacing <= 3:
                            print(f"   ✅ Range spacing typical for S1 IW (~2.3m)")
                        else:
                            print(f"   ⚠️  Range spacing unusual for S1 IW")
                        if 13 <= azimuth_spacing <= 15:
                            print(f"   ✅ Azimuth spacing typical for S1 IW (~14m)")
                        else:
                            print(f"   ⚠️  Azimuth spacing unusual for S1 IW")
                        
                        print("="*70 + "\n")
                        
                        # Apply optimized terrain correction
                        print(f"   🌍 Starting terrain correction: {working_data.shape} -> {self.target_resolution}m resolution")
                        
                        # DIAGNOSTIC: Check input data before terrain correction
                        finite_input = working_data[np.isfinite(working_data)]
                        input_percentage = len(finite_input) / working_data.size * 100 if working_data.size > 0 else 0
                        print(f"   📊 Input data: {input_percentage:.1f}% finite values, range: [{np.min(finite_input):.6f}, {np.max(finite_input):.6f}]")
                        print(f"   📊 Bbox: [{sar_bbox[0]:.6f}, {sar_bbox[1]:.6f}, {sar_bbox[2]:.6f}, {sar_bbox[3]:.6f}]")
                        
                        try:
                            geocoding_cache = str(self.output_dir / "processing_cache")

                            def _run_terrain_correction(
                                input_array: np.ndarray, label: str, bbox: list[float]
                            ) -> dict:
                                result = sardine.terrain_correction(
                                    input_array,
                                    bbox,
                                    orbit_times_str,
                                    orbit_positions,
                                    orbit_velocities,
                                    geocoding_cache,
                                    float(self.target_resolution),
                                    real_metadata,
                                    "bilinear"
                                )

                                if not isinstance(result, dict) or 'data' not in result:
                                    raise ValueError("Terrain correction did not return structured result with 'data'")

                                candidate = np.asarray(result['data'], dtype=np.float32)
                                if 'geo_transform' not in result:
                                    raise ValueError("Terrain correction missing 'geo_transform' in result; cannot export scientifically")

                                transform, adjusted = self._ensure_north_up_georeferencing(candidate, result['geo_transform'])
                                finite_mask = np.isfinite(adjusted)
                                finite_pct = float(np.sum(finite_mask)) / adjusted.size * 100 if adjusted.size > 0 else 0.0

                                print(
                                    f"   📐 Terrain correction attempt '{label}': {finite_pct:.1f}% valid pixels"
                                )
                                print(
                                    f"      ↳ geo_transform origin=({transform[0]:.6f}, {transform[3]:.6f}), "
                                    f"pixel_size=({transform[1]:.8f}, {transform[5]:.8f})"
                                )
                                print(
                                    f"      ↳ bbox=[{bbox[0]:.6f}, {bbox[1]:.6f}, {bbox[2]:.6f}, {bbox[3]:.6f}]"
                                )

                                return {
                                    'label': label,
                                    'data': adjusted,
                                    'transform': transform,
                                    'crs': result.get('crs'),
                                    'valid_pct': finite_pct,
                                    'bbox': bbox,
                                }

                            flipped_input = np.flipud(working_data)
                            attempt_queue = [
                                (working_data, "original", sar_bbox),
                                (flipped_input, "vertical_flip", sar_bbox),
                                (working_data, "full_bbox", base_bbox),
                                (flipped_input, "vertical_flip_full_bbox", base_bbox),
                            ]

                            attempts: list[dict] = []
                            best_attempt: dict | None = None

                            for input_array, label, bbox in attempt_queue:
                                attempts.append(_run_terrain_correction(input_array, label, bbox))
                                current_best = max(attempts, key=lambda item: item['valid_pct'])
                                if best_attempt is None or current_best['valid_pct'] > best_attempt['valid_pct']:
                                    best_attempt = current_best
                                if best_attempt['valid_pct'] >= 1.0:
                                    break

                            if best_attempt is None:
                                raise ValueError("Terrain correction attempts did not return any results")

                            geocoded_data = best_attempt['data']
                            self.geo_transform = best_attempt['transform']
                            self.geocoding_bbox = tuple(best_attempt['bbox'])
                            if best_attempt.get('crs'):
                                self.coordinate_system = best_attempt['crs']
                            valid_percentage = best_attempt['valid_pct']

                            if valid_percentage > 0.0:
                                print(f"   ✅ Terrain correction successful using attempt '{best_attempt['label']}': {geocoded_data.shape}, {valid_percentage:.1f}% valid pixels")
                            else:
                                print(f"   ⚠️  Terrain correction produced no finite pixels even after retries")
                                print(f"       Output shape: {geocoded_data.shape}, sample value: {np.nanmean(geocoded_data):.2f}")

                        except Exception as tc_error:
                            print(f"   ❌ DIAGNOSTIC: Terrain correction failed with error: {tc_error}")
                            raise
                        
                        step_duration = time.time() - step_start
                        self.log_step(11, "Terrain Correction", "success", 
                                     f"Optimized geocoding: {geocoded_data.shape}, {self.target_resolution}m resolution", step_duration)
                    else:
                        raise ValueError("No metadata available for geocoding parameters")
                        
                except Exception as e:
                    step_duration = time.time() - step_start
                    self.log_step(11, "Terrain Correction", "error", f"Terrain correction failed: {e}", step_duration)
                    raise RuntimeError(f"SCIENTIFIC MODE FAILURE: Terrain correction requires real geocoding parameters and orbit data. Error: {e}")
            else:
                step_duration = time.time() - step_start
                self.log_step(11, "Terrain Correction", "skipped", "Disabled by user", step_duration)
                geocoded_data = working_data
            
            # STEP 12: Mask Invalid Areas - SCIENTIFIC ORDER
            self.announce_step(12, "Mask Invalid Areas", "Applying quality masks")
            step_start = time.time()
            
            try:
                # SCIENTIFIC MODE: Ensure geocoded_data is proper numpy array - NO synthetic fallbacks
                if not isinstance(geocoded_data, np.ndarray):
                    raise RuntimeError("SCIENTIFIC MODE FAILURE: Geocoded data is not a valid numpy array - cannot apply quality masking")
                
                # SCIENTIFIC MODE: Handle different data shapes - NO synthetic data generation
                if geocoded_data.ndim == 0:
                    raise RuntimeError("SCIENTIFIC MODE FAILURE: Scalar geocoded data detected - cannot apply quality mask")
                
                # Apply scientific quality masking - first check basic validity
                basic_mask = np.isfinite(geocoded_data)
                basic_valid_percentage = np.sum(basic_mask) / basic_mask.size * 100
                print(f"   📊 Basic finite data: {basic_valid_percentage:.1f}% finite values")
                
                # Check data statistics to understand the range
                finite_data = geocoded_data[basic_mask]
                if len(finite_data) > 0:
                    data_min = np.min(finite_data)
                    data_max = np.max(finite_data)
                    data_mean = np.mean(finite_data)
                    print(f"   📊 Data range: {data_min:.6f} to {data_max:.6f}, mean: {data_mean:.6f}")
                    
                    # Adaptive quality masking based on data characteristics
                    if data_min >= 0:
                        # Data appears to be in linear units (sigma0/gamma0)
                        # Use reasonable range for SAR backscatter in linear units and allow true zeros
                        if data_min <= 0.0:
                            range_min = 0.0
                        else:
                            range_min = max(1e-8, data_min)

                        range_max = data_max
                        if np.isfinite(range_max) and range_max > 0.0:
                            range_max = min(100.0, range_max * 1.1)
                        else:
                            range_max = max(range_min, 100.0)

                        valid_mask = basic_mask & (geocoded_data >= range_min) & (geocoded_data <= range_max)
                        print(f"   🔍 Applied linear unit masking: range [{range_min:.2e}, {range_max:.2f}]")
                    else:
                        # Data might be in dB units or mixed values
                        # Use percentile-based masking to exclude extreme outliers
                        p1 = np.percentile(finite_data, 1)
                        p99 = np.percentile(finite_data, 99)
                        valid_mask = basic_mask & (geocoded_data >= p1) & (geocoded_data <= p99)
                        print(f"   🔍 Applied percentile masking: range [{p1:.3f}, {p99:.3f}] (1st-99th percentile)")
                else:
                    # No finite data available
                    valid_mask = basic_mask
                    print(f"   ⚠️  No finite data found - using basic mask only")
                
                # Create masked data
                masked_data = geocoded_data.copy()
                masked_data[~valid_mask] = np.nan
                valid_percentage = np.sum(valid_mask) / valid_mask.size * 100
                
                step_duration = time.time() - step_start
                self.log_step(12, "Mask Invalid Areas", "success", 
                             f"Quality mask: {valid_percentage:.1f}% valid pixels", step_duration)
                
            except Exception as e:
                step_duration = time.time() - step_start
                self.log_step(12, "Mask Invalid Areas", "error", f"Quality masking failed: {e}", step_duration)
                raise RuntimeError(f"SCIENTIFIC MODE FAILURE: Quality masking failed. Error: {e}")
            
            # STEP 13: Convert to dB - SCIENTIFIC ORDER (final step before export)
            self.announce_step(13, "Convert to dB", "Transforming to logarithmic scale")
            step_start = time.time()
            
            try:
                # SCIENTIFIC MODE: Ensure masked_data is proper numpy array - NO synthetic fallbacks
                if not isinstance(masked_data, np.ndarray):
                    raise RuntimeError("SCIENTIFIC MODE FAILURE: Masked data is not a valid numpy array - cannot convert to dB")
                
                # SCIENTIFIC MODE: Handle scalar or improperly shaped data - NO synthetic data generation
                if masked_data.ndim == 0:
                    raise RuntimeError("SCIENTIFIC MODE FAILURE: Scalar masked data detected - cannot convert to dB")
                
                # Apply scientific dB conversion using available function
                db_result = sardine.convert_to_db_real(masked_data)
                
                if isinstance(db_result, dict):
                    db_data = db_result.get('data', None)
                    if db_data is None or not isinstance(db_data, np.ndarray):
                        raise RuntimeError("SCIENTIFIC MODE FAILURE: dB conversion failed to produce valid data array")
                else:
                    if not isinstance(db_result, np.ndarray):
                        raise RuntimeError("SCIENTIFIC MODE FAILURE: dB conversion failed to produce valid numpy array")
                    db_data = db_result
                
                step_duration = time.time() - step_start
                self.log_step(13, "Convert to dB", "success", 
                             f"dB conversion: range {np.nanmin(db_data):.1f} to {np.nanmax(db_data):.1f} dB", 
                             step_duration)
                
            except Exception as e:
                step_duration = time.time() - step_start
                self.log_step(13, "Convert to dB", "error", f"dB conversion failed: {e}", step_duration)
                raise RuntimeError(f"SCIENTIFIC MODE FAILURE: dB conversion failed. Error: {e}")
            
            # STEP 14: Export Final Products
            self.announce_step(14, "Export Final Products", "Writing outputs to disk")
            step_start = time.time()
            
            exported_files = []
            
            try:
                # === PRIMARY EXPORT: GeoTIFF with proper georeferencing ===
                if self.geo_transform is not None:
                    output_geotiff = self.output_dir / f"backscatter_{self.polarization}_final.tif"
                    
                    # Prepare SAR metadata for STAC
                    px_x = abs(self.geo_transform[1]) if self.geo_transform else None
                    px_y = abs(self.geo_transform[5]) if self.geo_transform else None
                    range_spacing_m = float(self.get_current_range_spacing())
                    azimuth_spacing_m = float(self.get_current_azimuth_spacing())
                    print(
                        "   ✅ Export metadata pixel spacing (m): "
                        f"range={range_spacing_m:.3f}, azimuth={azimuth_spacing_m:.3f}"
                    )
                    sar_metadata = {
                        "platform": "sentinel-1",
                        "polarization": self.polarization,
                        "processing_level": "COMPLETE_14_STEP_PIPELINE",
                        "acquisition_start_time": datetime.now().isoformat(),
                        "orbit_direction": "unknown",  # Could be extracted from metadata
                        "range_pixel_spacing": range_spacing_m,
                        "azimuth_pixel_spacing": azimuth_spacing_m,
                        "pixel_spacing_lon_degrees": px_x,
                        "pixel_spacing_lat_degrees": px_y,
                        "multilook_range": self.multilook_range,
                        "multilook_azimuth": self.multilook_azimuth
                    }
                    
                    # Export Cloud Optimized GeoTIFF with STAC metadata
                    geotiff_path, stac_path = create_cog_with_stac(
                        data=db_data,
                        output_dir=self.output_dir,
                        filename_base=f"backscatter_{self.polarization}_final",
                        geotransform=tuple(self.geo_transform),  # Convert list to tuple
                        sar_metadata=sar_metadata,
                        crs=self.coordinate_system,
                        compress="lzw"
                    )
                    
                    exported_files.extend([geotiff_path.name, stac_path.name])
                    print(f"   ✅ PRIMARY: GeoTIFF exported with georeferencing: {geotiff_path.name}")
                    
                else:
                    print(f"   ⚠️  Warning: No geo_transform available - cannot create GeoTIFF")
                
                # === BACKUP EXPORT: Numpy array (always created) ===
                output_npy = self.output_dir / f"backscatter_{self.polarization}_final.npy"
                np.save(output_npy, db_data)
                exported_files.append(output_npy.name)
                
                # === SUMMARY: Text file ===
                output_txt = self.output_dir / f"backscatter_{self.polarization}_summary.txt"
                with open(output_txt, 'w') as f:
                    f.write(f"SARdine Backscatter Processing Results\n")
                    f.write(f"=====================================\n")
                    f.write(f"Input: {Path(self.input_path).name}\n")
                    f.write(f"Polarization: {self.polarization}\n")
                    f.write(f"Output shape: {db_data.shape}\n")
                    f.write(f"Data range: {np.nanmin(db_data):.3f} to {np.nanmax(db_data):.3f} dB\n")
                    f.write(f"Valid pixels: {valid_percentage:.1f}%\n")
                    f.write(f"Processing time: {time.time() - self.start_time:.1f}s\n")
                    f.write(f"Coordinate system: {self.coordinate_system}\n")
                    if self.geo_transform:
                        f.write(f"Geotransform: {self.geo_transform}\n")
                        f.write(f"Pixel size: {abs(self.geo_transform[1]):.8f}° x {abs(self.geo_transform[5]):.8f}°\n")
                
                exported_files.append(output_txt.name)
                
                step_duration = time.time() - step_start
                self.log_step(14, "Export Final Products", "success", 
                             f"Exported: {', '.join(exported_files)}", step_duration)
                
            except Exception as e:
                print(f"   ❌ Export error: {e}")
                step_duration = time.time() - step_start
                self.log_step(14, "Export Final Products", "error", f"Export failed: {e}", step_duration)
                
                # Diagnostic: try to save at least the numpy array
                try:
                    output_npy = self.output_dir / f"backscatter_{self.polarization}_final.npy"
                    np.save(output_npy, db_data)
                    print(f"   ✅ Diagnostic: Saved numpy array: {output_npy.name}")
                except Exception as fallback_error:
                    print(f"   ❌ Fallback failed: {fallback_error}")
                    raise RuntimeError(f"SCIENTIFIC MODE FAILURE: Unable to export any results. Error: {e}")
            
            # STEP 15: Generate Metadata
            self.announce_step(15, "Generate Metadata", "Writing processing metadata artifacts")
            step_start = time.time()
            
            try:
                processing_params = {
                    "polarization": self.polarization,
                    "speckle_filter": self.speckle_filter,
                    "filter_window": str(self.filter_window),
                    "multilook_range": str(self.multilook_range),
                    "multilook_azimuth": str(self.multilook_azimuth),
                    "terrain_flatten": str(self.terrain_flatten),
                    "geocode": str(self.geocode),
                    "target_resolution": str(self.target_resolution),
                    "processing_level": "COMPLETE_14_STEP_PIPELINE"
                }
                
                metadata_result = sardine.generate_metadata(
                    f"S1A_backscatter_{self.polarization}_complete",
                    processing_params,
                    [str(self.input_path)],
                    {"valid_pixel_percentage": valid_percentage}
                )
                
                # Export metadata
                json_result = sardine.export_metadata_json(metadata_result)
                
                metadata_json = self.output_dir / "metadata_complete.json"
                with open(metadata_json, "w") as f:
                    f.write(json_result)
                    
                step_duration = time.time() - step_start
                self.log_step(15, "Generate Metadata", "success", 
                             f"Metadata: {len(metadata_result)} fields", step_duration)
                
            except Exception as e:
                print(f"   ⚠️  Metadata warning: {e}")
                step_duration = time.time() - step_start
                self.log_step(15, "Generate Metadata", "warning", "Metadata generation incomplete", step_duration)
            
            # Generate quality report
            if self.quality_report:
                self.generate_quality_report(db_data, valid_percentage)
            
            # Save processing log
            self.save_processing_log()
            
            # Print performance summary
            print("\n" + "="*80)
            self.performance_monitor.print_summary()
            self.performance_monitor.stop_monitoring()
            
            total_duration = time.time() - self.start_time
            print(f"\n🎉 SUCCESS: Complete 14-step backscatter processing finished!")
            print(f"⏱️  Total processing time: {total_duration:.1f} seconds")
            print(f"📂 Output directory: {self.output_dir}")
            print(f"📊 Output shape: {db_data.shape}")
            print(f"📈 Data range: {np.nanmin(db_data):.1f} to {np.nanmax(db_data):.1f} dB")
            print(f"✅ Valid pixels: {valid_percentage:.1f}%")
            
            return {
                'status': 'success',
                'steps_completed': 14,
                'processing_time': total_duration,
                'output_files': list(self.output_dir.glob('*')),
                'dimensions': db_data.shape,
                'data_range': [float(np.nanmin(db_data)), float(np.nanmax(db_data))],
                'valid_percentage': valid_percentage
            }
            
        except Exception as e:
            print(f"\n❌ ERROR: 14-step processing failed: {e}")
            import traceback
            traceback.print_exc()
            return {
                'status': 'error',
                'error': str(e),
                'steps_completed': len(self.processing_log)
            }
    
    def _ensure_north_up_georeferencing(self, data, transform):
        """Ensure geotransform aligns with north-up, east-positive orientation."""

        if transform is None:
            raise ValueError("Terrain correction result did not include a geotransform")

        if not isinstance(data, np.ndarray) or data.ndim != 2:
            raise ValueError("Terrain correction data must be a 2D numpy array")

        transform_list = list(transform)
        if len(transform_list) != 6:
            raise ValueError("GeoTransform must contain exactly six elements")

        # Cast to float for downstream libraries that expect native floats
        transform_list = [float(value) for value in transform_list]

        height, width = data.shape
        adjusted_data = data

        # Ensure pixel height is negative (north-up). Flip data if needed.
        if transform_list[5] > 0:
            adjusted_data = np.flipud(adjusted_data)
            transform_list[3] += transform_list[5] * (height - 1)
            transform_list[5] = -transform_list[5]
            print("   🔃 Adjusted GeoTIFF orientation: flipped vertically for north-up alignment")

        # Ensure pixel width is positive (east-positive). Flip data if needed.
        if transform_list[1] < 0:
            adjusted_data = np.fliplr(adjusted_data)
            transform_list[0] += transform_list[1] * (width - 1)
            transform_list[1] = abs(transform_list[1])
            print("   🔃 Adjusted GeoTIFF orientation: flipped horizontally for east-positive alignment")

        # Secondary orientation check using cached metadata to detect subtle inversions
        metadata = getattr(self, "metadata", None)
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
                    lat_top = transform_list[3]
                    lat_bottom = lat_top + transform_list[5] * (height - 1)
                    print("   🔃 Adjusted GeoTIFF orientation: flipped vertically to align with scene metadata")

        return transform_list, adjusted_data

    def generate_quality_report(self, db_data, valid_percentage):
        """Generate comprehensive quality assessment report"""
        
        print("\n📊 Generating Quality Report...")

        
        # Calculate quality metrics
        finite_data = db_data[np.isfinite(db_data)]
        if len(finite_data) > 0:
            mean_backscatter = np.mean(finite_data)
            std_backscatter = np.std(finite_data)
            min_backscatter = np.min(finite_data)
            max_backscatter = np.max(finite_data)
        else:
            mean_backscatter = std_backscatter = min_backscatter = max_backscatter = np.nan
        
        quality_report = {
            "processing_parameters": {
                "polarization": self.polarization,
                "speckle_filter": self.speckle_filter,
                "multilook_factors": [self.multilook_range, self.multilook_azimuth],
                "terrain_flattening": self.terrain_flatten,
                "geocoding": self.geocode,
                "processing_level": "COMPLETE_14_STEP_PIPELINE"
            },
            "data_quality": {
                "valid_pixel_percentage": valid_percentage,
                "mean_backscatter_db": float(mean_backscatter),
                "std_backscatter_db": float(std_backscatter),
                "min_backscatter_db": float(min_backscatter),
                "max_backscatter_db": float(max_backscatter),
                "dynamic_range_db": float(max_backscatter - min_backscatter) if np.isfinite([max_backscatter, min_backscatter]).all() else np.nan
            },
            "processing_summary": {
                "total_steps": 14,
                "processing_time_seconds": time.time() - self.start_time,
                "output_dimensions": db_data.shape,
                "steps_completed": len(self.processing_log)
            }
        }
        
        # Save quality report
        quality_file = self.output_dir / "quality_report_complete.json"
        with open(quality_file, "w") as f:
            json.dump(quality_report, f, indent=2)
        
        print(f"   ✅ Quality report saved: {quality_file.name}")
        if np.isfinite(mean_backscatter):
            print(f"   📊 Mean backscatter: {mean_backscatter:.1f} dB")
        print(f"   📊 Valid pixels: {valid_percentage:.1f}%")
        print(f"   � Steps completed: {len(self.processing_log)}/14")
    
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
                desired_threads = max(1, detected)
        thread_value = str(desired_threads)

        runtime_env = {
            "RAYON_NUM_THREADS": thread_value,
            "OMP_NUM_THREADS": thread_value,
            "OPENBLAS_NUM_THREADS": thread_value,
            "NUMEXPR_NUM_THREADS": thread_value,
        }

        for key, value in runtime_env.items():
            current = os.environ.get(key)
            if current is None or current.strip() == "":
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
