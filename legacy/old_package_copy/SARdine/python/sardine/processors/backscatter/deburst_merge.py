"""
Deburst and radiometric calibration stage helpers for backscatter processing.
"""

import logging
import os
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from typing import Any, Dict, List, Tuple

import numpy as np
import sardine

# Valid IW subswath names
VALID_IW_SUBSWATHS = {"IW1", "IW2", "IW3"}


def subswath_sort_key(sw: str):
    """Sort key for IW subswaths - ensures IW1, IW2, IW3 ordering."""
    sw_upper = sw.upper()
    if sw_upper.startswith("IW") and sw_upper[2:].isdigit():
        return (0, int(sw_upper[2:]))
    return (1, sw_upper)


def normalize_subswath_order(subswaths: List[str]) -> List[str]:
    """Normalize subswath order to IW1, IW2, IW3."""
    return sorted(subswaths, key=subswath_sort_key)


def validate_subswath_names(subswaths: List[str], context: str = "") -> List[str]:
    """Validate and normalize subswath names.
    
    Args:
        subswaths: List of subswath names to validate
        context: Context string for error messages
        
    Returns:
        List of normalized, validated subswath names
        
    Raises:
        ValueError: If subswaths are invalid
    """
    logger = logging.getLogger(__name__)
    
    # Normalize to uppercase
    normalized = [sw.strip().upper() for sw in subswaths if sw.strip()]
    
    if not normalized:
        raise ValueError(f"Empty subswath list{' in ' + context if context else ''}")
    
    # Check for invalid subswath names
    invalid = [sw for sw in normalized if sw not in VALID_IW_SUBSWATHS]
    if invalid:
        raise ValueError(
            f"Invalid subswath names: {invalid}. "
            f"Expected IW mode subswaths: {VALID_IW_SUBSWATHS}"
        )
    
    # Check for duplicates
    if len(normalized) != len(set(normalized)):
        duplicates = [sw for sw in normalized if normalized.count(sw) > 1]
        logger.warning(f"Duplicate subswaths found: {set(duplicates)}")
        normalized = list(dict.fromkeys(normalized))  # Remove duplicates, preserve order
    
    return normalize_subswath_order(normalized)


def _populate_burst_timing_records(processor, deburst_overrides: Dict[str, Any]) -> None:
    """Populate processor.burst_timing_records from deburst overrides.
    
    The deburst stage emits timing information that terrain correction needs for
    accurate azimuth-to-time mapping. This function converts the deburst timing
    format to the BurstTiming format expected by the Rust terrain correction.
    
    Note: This is a SECONDARY path. The primary path is via cached_metadata.burst_timing_json
    which contains annotation-derived timing. This function provides a backup using
    deburst-derived timing which may be more accurate for actual data extent.
    
    Args:
        processor: BackscatterProcessor instance
        deburst_overrides: Dictionary mapping subswath IDs to deburst override dicts
    """
    logger = logging.getLogger(__name__)
    
    if not deburst_overrides:
        return
    
    timing_records = []
    
    for subswath_id, override in deburst_overrides.items():
        burst_timing = override.get("burst_timing")
        if not burst_timing or not isinstance(burst_timing, (list, tuple)):
            continue
        
        timing_reference = override.get("timing_reference", 0.0)
        azimuth_index_origin = override.get("azimuth_index_origin", 0)
        
        for bt in burst_timing:
            if not isinstance(bt, dict):
                continue
            
            burst_id = bt.get("burst_id", 0)
            t_start_rel = bt.get("t_start_rel", 0.0)
            t_end_rel = bt.get("t_end_rel", 0.0)
            line_count = bt.get("line_count_emitted", 0)
            prf_hz = bt.get("prf_hz", 0.0)
            
            # Convert to BurstTiming format expected by terrain correction
            # Estimate first/last line from timing info
            dt = bt.get("dt", 1.0 / prf_hz if prf_hz > 0 else 0.0)
            
            timing_records.append({
                "subswath_id": str(subswath_id).upper(),
                "burst_index": int(burst_id),
                "azimuth_time_rel_orbit": float(t_start_rel) if timing_reference else float(t_start_rel),
                "first_line_global": int(azimuth_index_origin),  # Approximate - actual comes from annotation
                "last_line_global": int(azimuth_index_origin + line_count),
                "first_valid_sample": None,
                "last_valid_sample": None,
            })
    
    if timing_records:
        processor.burst_timing_records = timing_records
        logger.info(f"Populated {len(timing_records)} burst_timing_records from deburst overrides")


def run_process_subswaths(processor, context):
    """
    Run deburst and calibration over each IW subswath.
    Returns (calibrated_subswaths, pipeline_results, primary_subswath).
    
    I/O Optimization: Prefetches next subswath's metadata while processing current subswath
    to overlap I/O with computation.
    """
    logger = logging.getLogger(__name__)
    processor.announce_step(4, "TOPSAR Deburst", "Deburst + radiometric calibration per subswath")
    step_start = time.time()
    
    # I/O Optimization: Prefetch cache for next subswath's metadata
    # This allows overlapping I/O (reading next subswath) with computation (processing current)
    prefetch_cache = {}
    prefetch_futures = {}

    calibrated_subswaths: Dict[str, np.ndarray] = {}
    deburst_overrides: Dict[str, Any] = {}
    available_subswaths: List[str] = []

    try:
        if not processor.metadata or not isinstance(processor.metadata, dict):
            raise ValueError("No valid metadata available for subswath processing")

        subswaths_str = processor.metadata.get("subswaths", "")
        if not subswaths_str:
            raise ValueError(
                "❌ SCIENTIFIC MODE FAILURE: Subswaths not found in metadata. "
                "Real subswath information is required for scientific processing. "
                "No hardcoded fallbacks permitted for data integrity."
            )

        parsed_subswaths = [sw.strip() for sw in subswaths_str.split(",") if sw.strip()]
        if not parsed_subswaths:
            raise ValueError(
                f"❌ SCIENTIFIC MODE FAILURE: Invalid subswaths format in metadata: '{subswaths_str}'"
            )

        # Validate and normalize subswaths using shared utility
        try:
            pipeline_subswaths = validate_subswath_names(parsed_subswaths, "metadata")
        except ValueError as e:
            raise ValueError(f"❌ SCIENTIFIC MODE FAILURE: {e}") from e
        
        available_subswaths = pipeline_subswaths

        print(f"   📡 Real subswaths from metadata: {', '.join(parsed_subswaths)}")
        if pipeline_subswaths != parsed_subswaths:
            print(f"   🔁 Normalized subswath order for processing: {', '.join(pipeline_subswaths)}")

        # Apply user-requested subswath filter if specified
        requested_subswaths = processor.options.get("subswaths") if hasattr(processor, "options") else None
        if requested_subswaths:
            # Normalize user request (accept list or comma-separated string)
            if isinstance(requested_subswaths, str):
                requested_subswaths = [s.strip().upper() for s in requested_subswaths.split(",") if s.strip()]
            else:
                requested_subswaths = [s.strip().upper() for s in requested_subswaths if s]
            
            # Filter to only requested subswaths that are available
            filtered = [sw for sw in pipeline_subswaths if sw.upper() in requested_subswaths]
            if not filtered:
                raise ValueError(
                    f"❌ Requested subswaths {requested_subswaths} not available. "
                    f"Available: {pipeline_subswaths}"
                )
            if len(filtered) < len(requested_subswaths):
                missing = [sw for sw in requested_subswaths if sw not in [f.upper() for f in filtered]]
                print(f"   ⚠️  Some requested subswaths not available: {missing}")
            
            pipeline_subswaths = filtered
            print(f"   🎯 Filtered to user-requested subswaths: {', '.join(pipeline_subswaths)}")

        estimated_data_size_mb = max(1, 500 * len(pipeline_subswaths))
        processor.performance_monitor.start_step(
            "Deburst & Calibration",
            data_size_mb=estimated_data_size_mb,
            parallel_enabled=processor.enable_parallel,
            threads_used=processor.num_threads,
            chunk_size=processor.chunk_size,
        )

        # Allow overriding subswath concurrency via env to reduce contention.
        env_subswath_threads = os.environ.get("SARDINE_SUBSWATH_THREADS")
        if env_subswath_threads:
            try:
                requested = max(1, int(env_subswath_threads))
            except Exception:
                requested = None
        else:
            requested = None

        max_workers = min(len(pipeline_subswaths), processor.num_threads or os.cpu_count() or 1)
        if requested is not None:
            max_workers = max(1, min(requested, len(pipeline_subswaths)))
        if max_workers <= 0:
            max_workers = 1

        print(f"   🧵 Subswath pipeline threads: {max_workers}")

        # Track geometry cache status for later validation
        if not hasattr(processor, "_geocoding_metadata_cache"):
            processor._geocoding_metadata_cache = None
        processor._geometry_cache_failed = False
        processor._geometry_cache_error = None
        
        if processor.geocode and hasattr(processor, "reader") and processor.reader is not None:
            try:
                print("   📊 Checking for pre-cached geocoding metadata...")
                
                # OPTIMIZATION: Reuse geometry from IW Split stage if available (~40s savings)
                cached_geometry = getattr(processor, "_geocoding_metadata_cache", None)
                if cached_geometry and isinstance(cached_geometry, dict) and "subswaths" in cached_geometry:
                    subswath_data = cached_geometry["subswaths"]
                    cached_meta = cached_geometry.get("cached_metadata")
                    cached_names = list(subswath_data.keys()) if isinstance(subswath_data, dict) else []
                    missing = [sw for sw in pipeline_subswaths if sw not in cached_names]
                    if missing:
                        print(f"   ⚠️  Cached geocoding metadata missing subswaths: {missing}; re-extracting from reader")
                        subswath_data = None
                    else:
                        print("   ⚡ Reusing cached geometry from IW Split stage (skipping redundant extraction)")
                if subswath_data is None:
                    # Fallback: extract geometry from reader
                    print("   📊 Extracting geocoding metadata from reader...")
                    with processor.reader_lock:
                        subswath_data = processor.reader.get_all_iw_subswaths()
                        cached_meta = (
                            processor.reader.get_cached_metadata()
                            if hasattr(processor.reader, "get_cached_metadata")
                            else None
                        )
                    
                    # Validate we got actual data
                    if not subswath_data:
                        raise ValueError("get_all_iw_subswaths() returned empty data")
                    
                    # Use correct key name that matches Rust expectations
                    processor._geocoding_metadata_cache = {
                        "subswaths": subswath_data,
                        "cached_metadata": cached_meta,
                    }
                    
                # Enumerate available geocoding metadata keys for diagnostics
                pol_keys = list(subswath_data.keys()) if isinstance(subswath_data, dict) else []
                pol_used = None
                pol_subswaths = None
                current_pol = getattr(processor, "polarization", None) or "VV"

                if isinstance(subswath_data, dict):
                    # Case 1: dict already keyed by IW* (no polarization nesting)
                    if all(isinstance(k, str) and k.upper().startswith("IW") for k in pol_keys):
                        pol_subswaths = subswath_data
                        pol_used = "<direct IW keys>"
                    else:
                        # Case 2: dict keyed by polarization -> subswath dict
                        for k in pol_keys:
                            k_upper = str(k).upper()
                            if k_upper == current_pol or current_pol in k_upper:
                                pol_used = k
                                pol_subswaths = subswath_data[k]
                                break
                        # Fallback: single entry
                        if pol_subswaths is None and len(pol_keys) == 1:
                            pol_used = pol_keys[0]
                            pol_subswaths = subswath_data[pol_used]

                swath_keys = list(pol_subswaths.keys()) if isinstance(pol_subswaths, dict) else []
                swath_count = len(swath_keys)
                total_entries = len(pol_keys) if pol_keys else swath_count
                print(
                    f"   ✅ Geocoding metadata extracted successfully: {swath_count} subswaths for {current_pol} "
                    f"(entries={total_entries}, pol_keys={pol_keys}, used_key={pol_used}, swaths={swath_keys})"
                )

                # Invariant: geocoding metadata must cover all pipeline subswaths for the selected polarization
                if not swath_keys:
                    raise RuntimeError(
                        f"Geocoding metadata missing subswaths for polarization {current_pol}. "
                        f"Available entries: {pol_keys}"
                    )

                missing_geo = [sw for sw in pipeline_subswaths if sw not in swath_keys]
                if missing_geo:
                    raise RuntimeError(
                        f"Geocoding metadata incomplete for polarization {current_pol}: missing {missing_geo}. "
                        f"Available subswaths: {swath_keys}; entries: {pol_keys}"
                    )
                
            except Exception as exc:
                # Track the failure for later stages to check
                processor._geometry_cache_failed = True
                processor._geometry_cache_error = str(exc)
                processor._geocoding_metadata_cache = None
                
                # Log with full context for debugging
                import traceback
                print(f"   ⚠️  Warning: Could not pre-extract geocoding metadata: {exc}")
                if processor.verbosity >= 2:
                    print(f"      Traceback: {traceback.format_exc()}")
                
                # Record validation warning
                processor._record_validation(
                    "deburst",
                    "geocoding_metadata_cache",
                    value=False,
                    expected=True,
                    passed=False,
                    severity="warning",
                    message=f"Geocoding metadata extraction failed: {exc}",
                )

        pipeline_results = []
        pipeline_errors: Dict[str, str] = {}  # Collect all errors
        
        # I/O Optimization Phase 2: Pre-read all SLC data in parallel
        # This overlaps I/O (reading all subswaths) before computation (processing)
        # Expected speedup: 44% for 2 subswaths, 59% for 3 subswaths
        slc_cache: Dict[str, np.ndarray] = {}
        enable_slc_prefetching = getattr(processor, 'enable_slc_prefetching', True)
        
        if enable_slc_prefetching and len(pipeline_subswaths) > 1:
            logger.info(f"🚀 I/O Optimization: Pre-reading {len(pipeline_subswaths)} subswaths in parallel")
            prefetch_start = time.time()
            
            def _pre_read_slc_data(subswath: str) -> tuple[str, Any]:
                """Pre-read SLC data for a subswath in background thread."""
                try:
                    with processor.reader_lock:
                        result = sardine.read_slc_data_for_subswath_only(
                            processor.reader,
                            subswath,
                            processor.polarization
                        )
                        
                        if isinstance(result, dict) and result.get("status") == "success":
                            slc_data = result.get("data")
                            if slc_data is not None:
                                return (subswath, np.asarray(slc_data, dtype=np.complex64))
                            else:
                                logger.warning(f"Pre-read succeeded but no data for {subswath}")
                                return (subswath, None)
                        else:
                            error_msg = result.get("message", "Unknown error") if isinstance(result, dict) else "Invalid result"
                            logger.warning(f"Pre-read failed for {subswath}: {error_msg}")
                            return (subswath, None)
                except Exception as e:
                    logger.warning(f"Pre-read exception for {subswath}: {e}")
                    return (subswath, None)
            
            # Pre-read all subswaths in parallel
            with ThreadPoolExecutor(max_workers=len(pipeline_subswaths)) as prefetch_executor:
                prefetch_futures = {
                    prefetch_executor.submit(_pre_read_slc_data, subswath): subswath
                    for subswath in pipeline_subswaths
                }
                
                for future in as_completed(prefetch_futures):
                    subswath = prefetch_futures[future]
                    try:
                        sw, slc_data = future.result(timeout=60.0)  # 60s timeout per subswath
                        if slc_data is not None:
                            slc_cache[sw] = slc_data
                            logger.info(f"✅ Pre-read SLC for {sw}: shape {slc_data.shape}")
                        else:
                            logger.warning(f"⚠️  Pre-read returned None for {sw}, will read synchronously")
                    except Exception as e:
                        logger.warning(f"⚠️  Pre-read failed for {subswath}: {e}, will read synchronously")
            
            prefetch_elapsed = time.time() - prefetch_start
            logger.info(f"📊 Pre-reading completed in {prefetch_elapsed:.1f}s: {len(slc_cache)}/{len(pipeline_subswaths)} subswaths cached")
        else:
            if not enable_slc_prefetching:
                logger.info("I/O optimization disabled (enable_slc_prefetching=False)")
            else:
                logger.info(f"Skipping pre-reading (only {len(pipeline_subswaths)} subswath)")
        
        # I/O Optimization: Prefetch next subswath's metadata while processing current
        # This overlaps I/O (reading next subswath) with computation (processing current)
        # The prefetch happens in a separate thread and caches results for later use
        def _prefetch_subswath_metadata(subswath: str) -> Dict[str, Any]:
            """Prefetch metadata for a subswath in background thread."""
            try:
                with processor.reader_lock:
                    # Pre-read annotation files to warm cache
                    if hasattr(processor.reader, "find_all_annotation_files"):
                        try:
                            processor.reader.find_all_annotation_files()
                        except Exception:
                            pass  # Ignore prefetch errors, will retry during actual processing
                return {"subswath": subswath, "status": "prefetched"}
            except Exception as e:
                logger.debug(f"Prefetch failed for {subswath}: {e}")
                return {"subswath": subswath, "status": "failed", "error": str(e)}
        
        # Start prefetching metadata for all subswaths in background
        prefetch_executor = ThreadPoolExecutor(max_workers=min(3, len(pipeline_subswaths)))
        prefetch_futures = {
            prefetch_executor.submit(_prefetch_subswath_metadata, subswath): subswath
            for subswath in pipeline_subswaths
        }
        
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            future_map = {
                executor.submit(
                    processor._process_subswath_pipeline,
                    subswath,
                    slc_cache.get(subswath)  # Use pre-read SLC data if available
                ): subswath
                for subswath in pipeline_subswaths
            }

            for future in as_completed(future_map):
                subswath = future_map[future]
                try:
                    result = future.result()

                    if result.get("needs_calibration", False):
                        pipeline_errors[subswath] = (
                            f"SCIENTIFIC MODE FAILURE: {subswath} data is not calibrated! "
                            f"Phase {result.get('phase', 'unknown')} returned uncalibrated power. "
                            f"This breaks all downstream processing. Set optimization_mode='complete'."
                        )
                        continue

                    calibrated_subswaths[subswath] = result["calibrated_data"]
                    pipeline_results.append(result)
                    if result.get("deburst_overrides"):
                        deburst_overrides[subswath] = result["deburst_overrides"]

                    if "fused_operations" in result:
                        ops = "+".join(result["fused_operations"])
                        print(
                            f"   ✅ {subswath}: {ops} fused {result['deburst_duration']:.1f}s, "
                            f"shape {result['deburst_shape'][0]}x{result['deburst_shape'][1]}"
                        )
                    else:
                        print(
                            f"   ✅ {subswath}: deburst {result['deburst_duration']:.1f}s, "
                            f"calibration {result['calibration_duration']:.1f}s, "
                            f"shape {result['deburst_shape'][0]}x{result['deburst_shape'][1]}"
                        )
                except Exception as pipeline_error:
                    # FIXED: Preserve full error context including stack trace
                    import traceback
                    error_info = {
                        "error": str(pipeline_error),
                        "type": type(pipeline_error).__name__,
                        "traceback": traceback.format_exc(),
                    }
                    pipeline_errors[subswath] = error_info
                    logger.error(
                        f"Subswath {subswath} failed: {pipeline_error}",
                        exc_info=True  # Include full traceback in logs
                    )
        
        # I/O Optimization: Wait for prefetch futures to complete (cleanup)
        # Prefetch is best-effort, so we don't fail if it errors
        for future in prefetch_futures:
            try:
                future.result(timeout=1.0)  # Don't wait long, prefetch is best-effort
            except Exception:
                pass  # Ignore prefetch errors, they're not critical
        
        # FIXED: Properly shutdown prefetch executor to prevent resource leak
        prefetch_executor.shutdown(wait=False)
        
        # Report all errors after all futures complete
        if pipeline_errors:
            # FIXED: Handle both string and dict error formats
            error_lines = []
            for sw, err in sorted(pipeline_errors.items()):
                if isinstance(err, dict):
                    error_lines.append(f"  - {sw}: {err['type']}: {err['error']}")
                    if processor.verbosity >= 2:
                        error_lines.append(f"    Traceback:\n{err['traceback']}")
                else:
                    error_lines.append(f"  - {sw}: {err}")
            error_summary = "\n".join(error_lines)
            raise RuntimeError(
                f"SCIENTIFIC MODE FAILURE: {len(pipeline_errors)} subswath(s) failed:\n{error_summary}"
            )

        if not calibrated_subswaths:
            raise RuntimeError("No subswaths were successfully processed")

        preferred_order = sorted(calibrated_subswaths.keys(), key=subswath_sort_key)
        primary_subswath = preferred_order[0]
        working_data = calibrated_subswaths[primary_subswath]

        processor._auto_tune_chunk_size(working_data)

        summary_lines = [
            f"{result['subswath']}={result['deburst_shape'][0]}x{result['deburst_shape'][1]}"
            for result in sorted(pipeline_results, key=lambda r: subswath_sort_key(r["subswath"]))
        ]

        step_duration = time.time() - step_start
        processor.log_step(
            4,
            "Deburst & Radiometric Calibration",
            "success",
            f"Processed {len(calibrated_subswaths)} subswaths: {', '.join(summary_lines)}",
            step_duration,
        )
        processor.performance_monitor.end_step()

        if deburst_overrides:
            processor._deburst_overrides = {k.upper(): v for k, v in deburst_overrides.items()}
        else:
            processor._deburst_overrides = {}

        for result in pipeline_results:
            bursts = result.get("burst_metadata")
            if bursts:
                processor._register_burst_metadata(result["subswath"], bursts)

        # CRITICAL FIX: Populate burst_timing_records from deburst overrides for terrain correction
        # The deburst stage emits burst_timing which terrain correction needs for azimuth mapping.
        # Previously this was stored in _deburst_overrides but never propagated to burst_timing_records.
        if processor._deburst_overrides:
            _populate_burst_timing_records(processor, processor._deburst_overrides)

        return calibrated_subswaths, pipeline_results, primary_subswath

    except Exception as exc:
        step_duration = time.time() - step_start
        processor.log_step(
            4,
            "Deburst & Radiometric Calibration",
            "error",
            f"Pipeline failed: {exc}",
            step_duration,
        )
        raise RuntimeError(f"SCIENTIFIC MODE FAILURE: Subswath pipeline failed. Error: {exc}")
