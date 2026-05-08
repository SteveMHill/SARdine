"""
Speckle filtering stage helper for the backscatter pipeline.
"""

import logging
import time
import numpy as np
import sardine


# Speckle filter parameters - configurable via processor options
DEFAULT_EDGE_THRESHOLD = 0.5
DEFAULT_DAMPING_FACTOR = 1.0
DEFAULT_CV_THRESHOLD = 0.5

# ENL (Equivalent Number of Looks) valid range for speckle filter
ENL_MIN = 0.1
ENL_MAX = 100.0

# Minimum image dimensions for meaningful speckle filtering
MIN_IMAGE_DIMENSION = 16
MIN_WINDOW_SIZE = 3

# Minimum percentage of finite pixels required
MIN_FINITE_PERCENTAGE = 10.0

logger = logging.getLogger(__name__)


def run_speckle_filtering(processor) -> None:
    """Apply speckle filtering with validated parameters."""
    processor.announce_step(9, "Speckle Filtering", "Reducing speckle while preserving radiometry")
    step_start = time.time()

    working_data = getattr(processor, "_working_data", None)

    try:
        if not isinstance(working_data, np.ndarray):
            raise RuntimeError(
                "SCIENTIFIC MODE FAILURE: Working data is not a valid numpy array - cannot proceed with speckle filtering"
            )

        if working_data.ndim == 0:
            raise RuntimeError(
                "SCIENTIFIC MODE FAILURE: Scalar data detected - cannot apply speckle filtering. Real SAR data required."
            )
        if working_data.ndim == 1:
            raise RuntimeError(
                "SCIENTIFIC MODE FAILURE: 1D data detected - cannot apply speckle filtering. Real SAR data required."
            )
        if working_data.ndim != 2:
            raise RuntimeError(
                f"SCIENTIFIC MODE FAILURE: Invalid data dimensions {working_data.ndim}D - speckle filtering requires 2D SAR data"
            )

        finite_mask = np.isfinite(working_data)
        if not np.any(finite_mask):
            raise RuntimeError(
                "SCIENTIFIC MODE FAILURE: All input data contains NaN or infinite values - cannot apply speckle filtering"
            )

        finite_percentage = float(np.sum(finite_mask)) / float(finite_mask.size) * 100.0
        if finite_percentage < MIN_FINITE_PERCENTAGE:
            raise RuntimeError(
                f"SCIENTIFIC MODE FAILURE: Insufficient valid data ({finite_percentage:.1f}%) for speckle filtering - minimum {MIN_FINITE_PERCENTAGE}% required"
            )

        data_size_mb = max(1.0, working_data.nbytes / (1024 * 1024))
        processor.performance_monitor.start_step(
            "Speckle Filtering",
            data_size_mb=data_size_mb,
            parallel_enabled=processor.enable_parallel,
            threads_used=processor.num_threads,
            chunk_size=processor.chunk_size,
        )

        filter_type = str(processor.speckle_filter)
        window_size = int(processor.filter_window) if processor.filter_window is not None else 7

        rows, cols = working_data.shape
        if rows < MIN_IMAGE_DIMENSION or cols < MIN_IMAGE_DIMENSION:
            logger.warning(
                f"Image too small ({rows}x{cols}) for meaningful speckle filtering "
                f"(minimum {MIN_IMAGE_DIMENSION}x{MIN_IMAGE_DIMENSION}). Skipping this step."
            )
            print(f"   ⚠️  Image too small ({rows}x{cols}) for meaningful speckle filtering, skipping this step...")
            step_duration = time.time() - step_start
            processor.log_step(
                9,
                "Speckle Filtering",
                "skipped",
                f"Image too small: {working_data.shape}",
                step_duration,
            )
            processor.performance_monitor.end_step()
            processor._working_data = np.asarray(working_data, dtype=np.float32)
            return

        min_image_dim = min(rows, cols)
        if min_image_dim < window_size:
            original_window = window_size
            window_size = max(MIN_WINDOW_SIZE, (min_image_dim // 2) * 2 - 1)
            logger.warning(
                f"Image size ({rows}x{cols}) requires reduced window size: {original_window} -> {window_size}"
            )
            print(
                f"   ⚠️  Image too small ({working_data.shape}) for window size {original_window}, reducing to {window_size}"
            )
            if window_size < MIN_WINDOW_SIZE:
                logger.warning(
                    f"Cannot reduce window size below {MIN_WINDOW_SIZE}. Skipping speckle filtering."
                )
                print("   ⚠️  Image too small for meaningful speckle filtering, skipping this step...")
                step_duration = time.time() - step_start
                processor.log_step(
                    9,
                    "Speckle Filtering",
                    "skipped",
                    f"Image too small: {working_data.shape}",
                    step_duration,
                )
                processor.performance_monitor.end_step()
                processor._working_data = np.asarray(working_data, dtype=np.float32)
                return

        if processor.actual_range_looks is None or processor.actual_azimuth_looks is None:
            raise RuntimeError(
                "SCIENTIFIC MODE FAILURE: Multilook factors were not initialised before speckle filtering"
            )

        num_looks = float(processor.actual_range_looks) * float(processor.actual_azimuth_looks)

        # SCIENTIFIC DECISION: Skip speckle filtering for very high ENL
        # At ENL > 100, speckle is already well-suppressed (σ/μ < 10%), and additional
        # filtering may degrade spatial resolution without meaningful benefit.
        # Reference: Lee et al. (1994), Oliver & Quegan (2004)
        if num_looks > ENL_MAX:
            logger.info(
                f"ENL {num_looks:.1f} exceeds {ENL_MAX} - speckle already well-suppressed. "
                f"Skipping speckle filtering to preserve spatial resolution."
            )
            print(
                f"   ℹ️  ENL {num_looks:.1f} > {ENL_MAX}: Speckle already suppressed (σ/μ < 10%), "
                f"skipping filter to preserve resolution"
            )
            step_duration = time.time() - step_start
            processor.log_step(
                9,
                "Speckle Filtering",
                "skipped",
                f"ENL {num_looks:.1f} exceeds {ENL_MAX} - speckle already well-suppressed",
                step_duration,
            )
            processor._working_data = np.asarray(working_data, dtype=np.float32)
            return
            
        if num_looks < ENL_MIN:
            logger.warning(
                f"ENL {num_looks:.2f} below minimum {ENL_MIN}. Raising to {ENL_MIN}. "
                f"This indicates insufficient multilooking for accurate speckle characterization."
            )
            print(
                f"   ⚠️  Speckle filter ENL {num_looks:.2f} raised to {ENL_MIN} minimum"
            )
            num_looks = ENL_MIN

        # Get filter parameters from options or use defaults
        options = getattr(processor, "options", {}) or {}
        edge_threshold = float(options.get("edge_threshold", DEFAULT_EDGE_THRESHOLD))
        damping_factor = float(options.get("damping_factor", DEFAULT_DAMPING_FACTOR))
        cv_threshold = float(options.get("cv_threshold", DEFAULT_CV_THRESHOLD))

        working_data = processor._ensure_contiguous(working_data)

        filtered_result = sardine.apply_speckle_filter(
            working_data,
            filter_type,
            window_size,
            num_looks,
            edge_threshold,
            damping_factor,
            cv_threshold,
        )

        # Validate speckle filter result with clear error messages
        if filtered_result is None:
            raise RuntimeError(
                "Speckle filter returned None - this indicates a failure in the Rust bindings"
            )
        
        filtered_data = None
        if isinstance(filtered_result, dict):
            # Try primary key first
            filtered_data = filtered_result.get("data")
            if filtered_data is None:
                # Try alternative key
                filtered_data = filtered_result.get("filtered_data")
            if filtered_data is None:
                available_keys = list(filtered_result.keys())
                raise RuntimeError(
                    f"Speckle filter result dict missing data. "
                    f"Expected 'data' or 'filtered_data' key. "
                    f"Available keys: {available_keys}"
                )
        elif isinstance(filtered_result, np.ndarray):
            filtered_data = filtered_result
        else:
            raise RuntimeError(
                f"Unexpected speckle filter result type: {type(filtered_result)}. "
                f"Expected dict or numpy.ndarray."
            )

        if not isinstance(filtered_data, np.ndarray):
            raise RuntimeError(
                f"Speckle filter data is not a numpy array: {type(filtered_data)}"
            )
        
        if filtered_data.ndim != 2:
            raise RuntimeError(
                f"Speckle filter produced {filtered_data.ndim}D array, expected 2D"
            )
        
        if filtered_data.size == 0:
            raise RuntimeError("Speckle filter produced empty array")

        # Validate output finite pixel count
        finite_after = np.isfinite(filtered_data)
        valid_after = (
            float(np.sum(finite_after)) / float(finite_after.size) * 100.0 if finite_after.size > 0 else 0.0
        )
        
        # Warn if significant data loss occurred
        if valid_after < finite_percentage - 5.0:
            logger.warning(
                f"Speckle filtering caused significant data loss: "
                f"{finite_percentage:.1f}% -> {valid_after:.1f}% finite pixels"
            )
        
        print(
            f"   ✅ Speckle filtering applied ({filter_type}): {valid_after:.1f}% finite pixels remaining"
        )

        step_duration = time.time() - step_start
        processor.log_step(
            9,
            "Speckle Filtering",
            "success",
            f"Filter: {filter_type}, shape: {filtered_data.shape}",
            step_duration,
        )
        processor.performance_monitor.end_step()

        processor._working_data = np.asarray(filtered_data, dtype=np.float32)

    except Exception as exc:
        step_duration = time.time() - step_start
        processor.log_step(9, "Speckle Filtering", "error", f"Speckle filtering failed: {exc}", step_duration)
        try:
            processor.performance_monitor.end_step()
        except Exception:
            pass
        raise RuntimeError(
            f"SCIENTIFIC MODE FAILURE: Speckle filtering failed. Error: {exc}"
        )
