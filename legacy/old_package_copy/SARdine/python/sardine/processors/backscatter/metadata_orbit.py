"""
Metadata and precise orbit acquisition helpers for backscatter.
"""

import logging
import os
import time
from datetime import datetime, timezone, timedelta
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

# ESA requirement: minimum 10 vectors for interpolation
# See: ESA Sentinel-1 Product Specification S1-RS-MDA-52-7441
MINIMUM_ORBIT_VECTORS = 10


def read_metadata(processor, context) -> None:
    """Stage 1: read metadata and reader handle."""
    import logging
    logger = logging.getLogger(__name__)
    
    processor.announce_step(1, "Read Metadata & Files", "Parsing Sentinel-1 annotation and SLC info")
    step_start = time.time()

    try:
        # Create reader with validation
        reader = processor.create_reader(self_ref=processor)
        if reader is None:
            raise RuntimeError("Failed to create SLC reader - reader is None")
        if not hasattr(reader, "get_cached_metadata") and not hasattr(reader, "get_metadata"):
            raise RuntimeError("Reader missing required metadata methods")
        
        # Safe metadata extraction with proper fallback handling
        metadata = None
        cache_error = None
        if hasattr(reader, "get_cached_metadata"):
            try:
                metadata = reader.get_cached_metadata()
            except Exception as e:
                cache_error = e
                logger.warning(f"Cache metadata failed: {e}, falling back to get_metadata()")
        
        if metadata is None:
            try:
                metadata = reader.get_metadata()
            except Exception as fallback_err:
                if cache_error:
                    raise RuntimeError(
                        f"Both cached and fallback metadata extraction failed: "
                        f"cache={cache_error}, fallback={fallback_err}"
                    ) from cache_error
                raise
        
        # Validate structure before assignment
        if not isinstance(metadata, dict):
            raise ValueError(f"Metadata must be dict, got {type(metadata)}")
        
        # Check for critical metadata fields
        required_fields = ["product_id", "start_time"]
        missing_fields = [f for f in required_fields if f not in metadata or metadata[f] is None]
        if missing_fields:
            raise ValueError(
                f"Missing required metadata fields: {missing_fields}. "
                f"Check SAFE manifest parsing."
            )
        
        # Validate content before assignment
        validated = processor.validate_metadata(metadata)
        
        # Only assign after successful validation
        processor.metadata = metadata
        processor.reader = reader
        processor.validated_metadata = validated

        processor.output_epsg = processor._resolve_output_epsg()
        processor.coordinate_system = f"EPSG:{processor.output_epsg}"

        step_duration = time.time() - step_start
        
        # Safe metadata info string construction
        try:
            if isinstance(metadata, dict):
                metadata_info = f"Cached metadata extracted: dict ({len(metadata)} keys)"
            elif hasattr(metadata, "__len__"):
                try:
                    metadata_info = f"Cached metadata extracted: {type(metadata).__name__} ({len(metadata)} items)"
                except (TypeError, AttributeError):
                    metadata_info = f"Cached metadata extracted: {type(metadata).__name__}"
            else:
                metadata_info = f"Cached metadata extracted: {type(metadata).__name__}"
        except Exception as e:
            metadata_info = f"Cached metadata extracted: {type(metadata).__name__} (info extraction failed: {e})"
        
        # Conditional performance indicator
        if cache_error is None and hasattr(reader, "get_cached_metadata"):
            metadata_info += " [PERFORMANCE: cached mode]"
        else:
            metadata_info += " [PERFORMANCE: uncached mode]"
            
        processor.log_step(1, "Read Metadata & Files", "success", metadata_info, step_duration)
        processor.performance_monitor.end_step()

        context.set_artifact("metadata", metadata)
        context.set_artifact("reader", reader)
    except Exception as exc:
        step_duration = time.time() - step_start
        input_path = getattr(processor, "input_path", "unknown")
        reader_type = type(reader).__name__ if 'reader' in dir() and reader is not None else "not created"
        has_cached = hasattr(reader, 'get_cached_metadata') if 'reader' in dir() and reader is not None else False
        error_msg = (
            f"Metadata read failed for {input_path}: {exc}\n"
            f"Reader type: {reader_type}\n"
            f"Has cached method: {has_cached}"
        )
        processor.log_step(1, "Read Metadata & Files", "error", error_msg, step_duration)
        raise


def apply_precise_orbit(processor, context) -> None:
    """Stage 2: download and validate precise orbit.
    
    Precise orbit files are ALWAYS required for accurate SAR processing:
    - Improves satellite state vectors from meter to cm-level accuracy
    - Required for proper burst timing in deburst stage
    - Required for accurate azimuth/range calculations
    - Required for Range-Doppler geocoding (when enabled)
    
    This stage should NOT be skipped even with --no-geocode.
    """
    processor.announce_step(2, "Apply Precise Orbit File")
    step_start = time.time()

    try:
        if not processor.metadata or not isinstance(processor.metadata, dict):
            raise ValueError("No cached metadata available for orbit acquisition")

        product_id = processor.metadata.get("product_id")
        start_time = processor.metadata.get("start_time")
        if not product_id or not start_time:
            raise ValueError(
                "Missing required metadata for orbit retrieval: "
                f"product_id={product_id}, start_time={start_time}"
            )

        # Parse start_time with comprehensive validation
        if isinstance(start_time, datetime):
            sensing_start = start_time
            if sensing_start.tzinfo is None:
                sensing_start = sensing_start.replace(tzinfo=timezone.utc)
            else:
                sensing_start = sensing_start.astimezone(timezone.utc)
        elif isinstance(start_time, str):
            # Normalize ISO format
            normalized_time = start_time.replace("Z", "+00:00")
            try:
                sensing_start = datetime.fromisoformat(normalized_time)
            except ValueError as e:
                raise ValueError(
                    f"Invalid start_time format: {start_time!r}. "
                    f"Expected ISO 8601 format (e.g., '2023-01-01T12:00:00Z'). "
                    f"Error: {e}"
                ) from e
            
            if sensing_start.tzinfo is None:
                sensing_start = sensing_start.replace(tzinfo=timezone.utc)
            else:
                sensing_start = sensing_start.astimezone(timezone.utc)
        else:
            raise ValueError(
                f"Unsupported start_time type for orbit retrieval: {type(start_time)!r}. "
                f"Expected datetime or ISO 8601 string."
            )
        
        # Validate datetime is reasonable (not too far in past/future)
        now = datetime.now(timezone.utc)
        age_days = (now - sensing_start).total_seconds() / 86400
        if age_days > 365 * 10:  # More than 10 years old
            logger = logging.getLogger(__name__)
            logger.warning(
                f"Start time is very old ({age_days:.0f} days ago): {sensing_start}. "
                f"Orbit files may not be available for very old products."
            )
        if age_days < -1:  # More than 1 day in future
            raise ValueError(
                f"Start time is in the future: {sensing_start}. "
                f"Orbit files are not available for future acquisitions."
            )

        # Use user-provided orbit directory or fall back to output_dir/orbit_cache
        orbit_cache_dir = processor.options.get("orbit_dir") or str(processor.output_dir / "orbit_cache")
        orbit_path = processor.download_orbit_file(sensing_start, None, orbit_cache_dir)
        processor.precise_orbit_path = orbit_path
        processor._register_orbit_file_for_metadata(orbit_path)

        # Validate orbit result structure
        if not hasattr(processor, "_last_orbit_result") or processor._last_orbit_result is None:
            raise RuntimeError(
                "Orbit validation result not available. "
                "This indicates download_orbit_file() failed silently or validation was skipped."
            )
        
        orbit_meta = processor._last_orbit_result
        if not isinstance(orbit_meta, dict):
            raise RuntimeError(f"Invalid orbit result type: {type(orbit_meta)}")
        
        # Extract result with proper validation
        result_meta = orbit_meta.get("result")
        if result_meta is None:
            # Try direct access if "result" key doesn't exist
            result_meta = orbit_meta
        elif not isinstance(result_meta, dict):
            raise RuntimeError(f"Orbit result has invalid 'result' field type: {type(result_meta)}")
        
        # Extract vector count with validation
        vector_count = result_meta.get("orbit_vectors_count")
        if vector_count is None:
            raise RuntimeError(
                "Orbit validation did not report orbit_vectors_count. "
                f"Available keys: {list(result_meta.keys())}"
            )
        
        # Ensure vector_count is an integer
        try:
            vector_count = int(vector_count)
        except (TypeError, ValueError) as e:
            raise RuntimeError(
                f"Invalid orbit_vectors_count type: {type(vector_count)}, value: {vector_count}"
            ) from e
        
        if vector_count < MINIMUM_ORBIT_VECTORS:
            raise RuntimeError(
                f"Precise orbit file contains only {vector_count} state vectors "
                f"(minimum required: {MINIMUM_ORBIT_VECTORS} for scientific processing)"
            )
        
        # Validate orbit coverage if time range is available
        orbit_start_str = result_meta.get("start_time") or result_meta.get("validity_start")
        orbit_stop_str = result_meta.get("stop_time") or result_meta.get("validity_stop")
        if orbit_start_str and orbit_stop_str:
            try:
                orbit_start_dt = datetime.fromisoformat(str(orbit_start_str).replace("Z", "+00:00"))
                orbit_stop_dt = datetime.fromisoformat(str(orbit_stop_str).replace("Z", "+00:00"))

                logger = logging.getLogger(__name__)

                # Derive product time range; stop_time may not be present in all metadata
                product_start_dt = sensing_start
                product_stop_raw = processor.metadata.get("stop_time")
                if isinstance(product_stop_raw, datetime):
                    product_stop_dt = product_stop_raw.astimezone(timezone.utc)
                elif isinstance(product_stop_raw, str):
                    try:
                        product_stop_dt = datetime.fromisoformat(str(product_stop_raw).replace("Z", "+00:00"))
                        if product_stop_dt.tzinfo is None:
                            product_stop_dt = product_stop_dt.replace(tzinfo=timezone.utc)
                        else:
                            product_stop_dt = product_stop_dt.astimezone(timezone.utc)
                    except Exception:
                        product_stop_dt = product_start_dt
                else:
                    # Fallback: treat product as effectively instantaneous
                    product_stop_dt = product_start_dt

                # Require orbit to cover full sensing window with a safety margin
                margin_seconds = float(processor.options.get("orbit_time_margin_s", 60.0))
                margin = margin_seconds
                covered = (
                    orbit_start_dt <= (product_start_dt - timedelta(seconds=margin))
                    and orbit_stop_dt >= (product_stop_dt + timedelta(seconds=margin))
                )

                # Always log coverage diagnostics
                logger.info(
                    "Orbit coverage check: product=[%s, %s], orbit=[%s, %s], margin=%.1fs, covered=%s",
                    product_start_dt,
                    product_stop_dt,
                    orbit_start_dt,
                    orbit_stop_dt,
                    margin_seconds,
                    covered,
                )

                # In strict orbit or strict science mode, treat insufficient coverage as a hard error
                strict_orbit = bool(getattr(processor, "strict_orbit_mode", False)) or bool(
                    getattr(processor, "strict_science", False)
                )
                if strict_orbit:
                    processor._require_validation(
                        covered,
                        "orbit",
                        "orbit_coverage",
                        value={
                            "product_start": product_start_dt.isoformat(),
                            "product_stop": product_stop_dt.isoformat(),
                            "orbit_start": orbit_start_dt.isoformat(),
                            "orbit_stop": orbit_stop_dt.isoformat(),
                            "margin_s": margin_seconds,
                        },
                        expected="orbit coverage must fully span sensing window with margin",
                        severity="error",
                        message=(
                            "Orbit file does not fully cover product sensing window with required margin. "
                            "This can cause invalid burst timing and TOPS deburst stripes."
                        ),
                    )
                else:
                    if not covered:
                        logger.warning(
                            "Orbit file may not fully cover product time range. "
                            "Product=[%s, %s], Orbit=[%s, %s]",
                            product_start_dt,
                            product_stop_dt,
                            orbit_start_dt,
                            orbit_stop_dt,
                        )

            except (ValueError, TypeError) as e:
                logger = logging.getLogger(__name__)
                logger.debug(f"Could not validate orbit coverage: {e}")

        processor._require_validation(
            vector_count >= MINIMUM_ORBIT_VECTORS,
            "orbit",
            "orbit_vectors_count",
            value=vector_count,
            expected=f">={MINIMUM_ORBIT_VECTORS}",
            message="Insufficient orbit vectors for precise processing",
        )
        processor._record_validation(
            "orbit",
            "orbit_file",
            value=orbit_path,
            expected="cached",
            passed=Path(orbit_path).exists(),
            severity="info" if Path(orbit_path).exists() else "warning",
            message="Orbit file resolved" if Path(orbit_path).exists() else "Orbit file should exist in cache",
        )

        step_duration = time.time() - step_start
        orbit_status = f"Precise orbit ready: {Path(orbit_path).name} ({vector_count} vectors)"
        processor.log_step(2, "Apply Precise Orbit File", "success", orbit_status, step_duration)
    except Exception as exc:
        step_duration = time.time() - step_start
        processor.log_step(2, "Apply Precise Orbit File", "error", f"Orbit processing failed: {exc}", step_duration)
        raise RuntimeError(f"SCIENCE MODE FAILURE: precise orbit acquisition failed: {exc}")
