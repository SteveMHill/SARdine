"""
Multilooking stage helper for the backscatter pipeline.
"""

import logging
import numpy as np
import sardine
import time


# Typical look factor ranges for Sentinel-1 TOPS mode
TYPICAL_RANGE_LOOKS = (4, 6)
TYPICAL_AZIMUTH_LOOKS = (1, 2)
# Maximum valid look factor to prevent unreasonable downsampling
MAX_VALID_LOOKS = 20
# Power preservation tolerance (5%)
POWER_PRESERVATION_TOLERANCE = 0.05
# Valid scale range for power rescaling
POWER_SCALE_MIN = 0.05
POWER_SCALE_MAX = 20.0

logger = logging.getLogger(__name__)


def run_multilooking(processor) -> None:
    """Apply scientific multilooking after merging subswaths."""
    processor.announce_step(7, "Multilooking", "Applying scientific multilook factors")
    step_start = time.time()

    working_data = getattr(processor, "_working_data", None)

    try:
        if not isinstance(working_data, np.ndarray):
            raise ValueError("Working data is not a valid numpy array")

        if not processor.metadata or not isinstance(processor.metadata, dict):
            raise ValueError("No valid metadata available for multilooking")

        input_range_spacing = float(processor.get_current_range_spacing())
        input_azimuth_spacing = float(processor.get_current_azimuth_spacing())

        if processor.validated_metadata is None:
            raise ValueError("Validated metadata not available for multilooking parameters")

        # Extract incidence angles with graceful fallback
        incidence_near = None
        incidence_far = None
        
        # Try validated_metadata first
        try:
            incidence_near = float(processor.validated_metadata.get("incidence_angle_near"))
            incidence_far = float(processor.validated_metadata.get("incidence_angle_far"))
        except (TypeError, ValueError):
            pass
        
        # Fallback to raw metadata if not in validated
        if incidence_near is None or incidence_far is None:
            try:
                if isinstance(processor.metadata, dict):
                    incidence_near = incidence_near or float(processor.metadata.get("incidence_angle_near", 0))
                    incidence_far = incidence_far or float(processor.metadata.get("incidence_angle_far", 0))
            except (TypeError, ValueError):
                pass
        
        # SCIENTIFIC FALLBACK: Use typical Sentinel-1 swath-based incidence angles
        # IW mode has well-defined incidence angle ranges per subswath:
        # IW1: 29.1°-35.5°, IW2: 34.1°-41.0°, IW3: 39.3°-46.0°
        # For merged IW, the full range is approximately 29°-46° (mid ~37.5°)
        # SM mode: typically 20°-45° depending on beam
        if not incidence_near or incidence_near <= 0 or not incidence_far or incidence_far <= 0:
            # --- strict_science gate: reject fallback angles ---
            strict_science = bool(getattr(processor, 'strict_science', False))
            if strict_science:
                raise RuntimeError(
                    "STRICT SCIENCE MODE FAILURE: Incidence angles missing from metadata.\n"
                    "incidence_angle_near and incidence_angle_far are required for accurate\n"
                    "multilook factor calculation. Mode-based fallbacks would introduce\n"
                    "5-15% errors in ground-range spacing and therefore in the output\n"
                    "pixel size.\n\n"
                    "FIX: Check SAFE annotation parsing; incidence angles should be in\n"
                    "the <product> element of the annotation XML."
                )

            # Try to get mode from metadata to pick appropriate defaults
            mode = ""
            if isinstance(processor.metadata, dict):
                mode = str(processor.metadata.get("mode", "")).upper()
            
            if "IW" in mode:
                fallback_near = 29.1  # IW1 near edge
                fallback_far = 46.0   # IW3 far edge
            elif "EW" in mode:
                fallback_near = 18.9
                fallback_far = 47.0
            elif "SM" in mode:
                fallback_near = 20.0
                fallback_far = 45.0
            else:
                fallback_near = 29.0
                fallback_far = 46.0
            
            logger.warning(
                f"⚠️  SCIENTIFIC WARNING: incidence_angle_near/far not found in metadata. "
                f"Using mode-based fallback values: near={fallback_near}°, far={fallback_far}° (mode={mode or 'unknown'}). "
                f"This may affect multilook factor calculation by 5-15%%. "
                f"Enable strict_science mode to prevent this fallback."
            )
            
            incidence_near = incidence_near if (incidence_near and incidence_near > 0) else fallback_near
            incidence_far = incidence_far if (incidence_far and incidence_far > 0) else fallback_far

        incidence_angle_deg = 0.5 * (incidence_near + incidence_far)
        incidence_angle_rad = np.radians(incidence_angle_deg)
        ground_range_spacing = input_range_spacing / np.sin(incidence_angle_rad)
        azimuth_spacing_m = input_azimuth_spacing

        if not np.isfinite(ground_range_spacing) or ground_range_spacing <= 0.0:
            raise ValueError(
                f"Invalid ground range spacing derived from metadata: {ground_range_spacing}"
            )

        if processor.target_resolution is None:
            raise ValueError("Target resolution was not initialised from metadata or options")
        target_resolution = float(processor.target_resolution)
        if not np.isfinite(target_resolution) or target_resolution <= 0.0:
            raise ValueError(f"Invalid target resolution: {target_resolution}m (must be positive)")

        range_looks_optimal = max(1, int(np.round(target_resolution / ground_range_spacing)))
        azimuth_looks_optimal = max(1, int(np.round(target_resolution / azimuth_spacing_m)))

        range_looks = range_looks_optimal
        azimuth_looks = azimuth_looks_optimal

        if range_looks > MAX_VALID_LOOKS:
            print(
                f"      ⚠️  Range looks {range_looks} capped to {MAX_VALID_LOOKS} to satisfy multilooking limits"
            )
            range_looks = MAX_VALID_LOOKS
        if azimuth_looks > MAX_VALID_LOOKS:
            print(
                f"      ⚠️  Azimuth looks {azimuth_looks} capped to {MAX_VALID_LOOKS} to satisfy multilooking limits"
            )
            azimuth_looks = MAX_VALID_LOOKS

        warnings_ml = []
        if range_looks < TYPICAL_RANGE_LOOKS[0] or range_looks > TYPICAL_RANGE_LOOKS[1]:
            warnings_ml.append(
                f"Range looks {range_looks} outside typical TOPS range {TYPICAL_RANGE_LOOKS} "
                f"(driven by {target_resolution}m target resolution)"
            )
        if azimuth_looks < TYPICAL_AZIMUTH_LOOKS[0] or azimuth_looks > TYPICAL_AZIMUTH_LOOKS[1]:
            warnings_ml.append(
                f"Azimuth looks {azimuth_looks} outside typical TOPS range {TYPICAL_AZIMUTH_LOOKS} "
                f"(driven by {target_resolution}m target resolution)"
            )

        print(f"   🔬 SCIENTIFIC MULTILOOKING:")
        print(f"      Incidence angle (avg): {incidence_angle_deg:.2f}°")
        print(f"      Slant spacing (input): {input_range_spacing:.2f}m × {azimuth_spacing_m:.2f}m (range×az)")
        print(f"      Ground spacing (input): {ground_range_spacing:.2f}m (range), {azimuth_spacing_m:.2f}m (az)")
        print(f"      Calculated looks: {range_looks}×{azimuth_looks} (range×azimuth)")
        print(f"      Target resolution: {target_resolution:.2f}m")
        for warning in warnings_ml:
            print(f"      ⚠️  {warning}")
        print(f"      ENL target: ~{range_looks * azimuth_looks}")

        # CRITICAL: Store native spacing BEFORE multilooking for terrain correction
        # The Rust terrain correction needs native SLC spacing, not multilooked spacing
        if processor.metadata is not None and isinstance(processor.metadata, dict):
            if "native_range_pixel_spacing" not in processor.metadata:
                processor.metadata["native_range_pixel_spacing"] = input_range_spacing
                logger.debug(f"Stored native_range_pixel_spacing={input_range_spacing:.4f}m before multilooking")
            if "native_azimuth_pixel_spacing" not in processor.metadata:
                processor.metadata["native_azimuth_pixel_spacing"] = input_azimuth_spacing
                logger.debug(f"Stored native_azimuth_pixel_spacing={input_azimuth_spacing:.4f}m before multilooking")

        ml_result = sardine.apply_multilooking(
            working_data,
            int(range_looks),
            int(azimuth_looks),
            float(input_range_spacing),
            float(input_azimuth_spacing),
        )

        # Validate multilook result structure
        output_range_spacing = None
        output_azimuth_spacing = None
        
        if ml_result is None:
            raise ValueError(
                "Multilooking returned None - this indicates a failure in the Rust bindings"
            )
        
        if isinstance(ml_result, dict):
            multilooked_data = ml_result.get("data")
            if multilooked_data is None:
                available_keys = list(ml_result.keys())
                raise ValueError(
                    f"Multilooking result dict missing 'data' key. "
                    f"Available keys: {available_keys}"
                )
            output_range_spacing = ml_result.get("range_spacing")
            output_azimuth_spacing = ml_result.get("azimuth_spacing")
        elif isinstance(ml_result, np.ndarray):
            multilooked_data = ml_result
        else:
            raise ValueError(
                f"Unexpected multilooking result type: {type(ml_result)}. "
                f"Expected dict or numpy.ndarray."
            )

        # Validate output spacing - calculate if missing and warn
        expected_range_spacing = input_range_spacing * float(range_looks)
        expected_azimuth_spacing = input_azimuth_spacing * float(azimuth_looks)
        
        if output_range_spacing is None:
            logger.debug(
                f"Multilook result missing range_spacing; using calculated value: {expected_range_spacing:.2f}m"
            )
            output_range_spacing = expected_range_spacing
        else:
            # Validate returned spacing matches expected
            spacing_diff_range = abs(output_range_spacing - expected_range_spacing)
            if spacing_diff_range > 0.5:  # Allow 0.5m tolerance
                logger.warning(
                    f"Multilook range spacing mismatch: returned {output_range_spacing:.2f}m, "
                    f"expected {expected_range_spacing:.2f}m (diff={spacing_diff_range:.2f}m)"
                )
        
        if output_azimuth_spacing is None:
            logger.debug(
                f"Multilook result missing azimuth_spacing; using calculated value: {expected_azimuth_spacing:.2f}m"
            )
            output_azimuth_spacing = expected_azimuth_spacing
        else:
            spacing_diff_azimuth = abs(output_azimuth_spacing - expected_azimuth_spacing)
            if spacing_diff_azimuth > 0.5:
                logger.warning(
                    f"Multilook azimuth spacing mismatch: returned {output_azimuth_spacing:.2f}m, "
                    f"expected {expected_azimuth_spacing:.2f}m (diff={spacing_diff_azimuth:.2f}m)"
                )

        if not isinstance(multilooked_data, np.ndarray):
            raise ValueError(
                f"Multilooking data is not a numpy array: {type(multilooked_data)}"
            )
        
        if multilooked_data.ndim != 2:
            raise ValueError(
                f"Multilooking produced {multilooked_data.ndim}D array, expected 2D"
            )
        
        if multilooked_data.size == 0:
            raise ValueError("Multilooking produced empty array")

        processor.update_current_spacing(output_range_spacing, output_azimuth_spacing)
        processor._range_multilook_factor = float(range_looks)
        processor._azimuth_multilook_factor = float(azimuth_looks)
        processor._multilooked_shape = tuple(int(v) for v in multilooked_data.shape)
        ground_output_spacing = (output_range_spacing / np.sin(incidence_angle_rad)) if np.sin(incidence_angle_rad) != 0 else float("nan")
        print(
            f"      Output spacing (slant): {output_range_spacing:.2f}m × {output_azimuth_spacing:.2f}m (range×azimuth)"
        )
        print(
            f"      Output spacing (ground equiv.): {ground_output_spacing:.2f}m × {output_azimuth_spacing:.2f}m"
        )
        
        # Diagnostic logging: After multilooking
        processor._log_diagnostic_statistics(
            "After Multilooking",
            multilooked_data,
            context=f"looks={range_looks}×{azimuth_looks}"
        )

        if processor._native_range_samples and float(range_looks) > 0:
            expected_range = processor._native_range_samples / float(range_looks)
            diff_range = abs(expected_range - multilooked_data.shape[1])
            range_aligned = diff_range <= 2.0
            if range_aligned:
                processor._record_validation(
                    "multilook",
                    "range_dimension_alignment",
                    value=multilooked_data.shape[1],
                    expected=f"{expected_range:.1f}±2",
                    passed=True,
                    severity="info",
                    message="Multilooked range dimension within tolerance",
                )
            else:
                processor._require_validation(
                    False,
                    "multilook",
                    "range_dimension_alignment",
                    value=multilooked_data.shape[1],
                    expected=f"{expected_range:.1f}±2",
                    severity="warning",
                    message="Multilooked range dimension deviates from expected",
                )

        if processor._native_azimuth_lines and float(azimuth_looks) > 0:
            expected_az = processor._native_azimuth_lines / float(azimuth_looks)
            diff_az = abs(expected_az - multilooked_data.shape[0])
            az_aligned = diff_az <= 2.0
            if az_aligned:
                processor._record_validation(
                    "multilook",
                    "azimuth_dimension_alignment",
                    value=multilooked_data.shape[0],
                    expected=f"{expected_az:.1f}±2",
                    passed=True,
                    severity="info",
                    message="Multilooked azimuth dimension within tolerance",
                )
            else:
                processor._require_validation(
                    False,
                    "multilook",
                    "azimuth_dimension_alignment",
                    value=multilooked_data.shape[0],
                    expected=f"{expected_az:.1f}±2",
                    severity="warning",
                    message="Multilooked azimuth dimension deviates from expected",
                )

        processor.actual_range_looks = range_looks
        processor.actual_azimuth_looks = azimuth_looks

        step_duration = time.time() - step_start
        processor.log_step(
            7,
            "Multilooking",
            "success",
            "Scientific multilook: "
            f"{range_looks}×{azimuth_looks}, shape: {multilooked_data.shape}, "
            f"ENL≈{range_looks*azimuth_looks}, spacing: {output_range_spacing:.2f}m×{output_azimuth_spacing:.2f}m",
            step_duration,
        )
        processor._record_stage_timing("Multilooking", step_duration)

        # Validate power preservation with tolerance for multilooking
        # Expected behavior: mean power should be preserved within ~5% tolerance
        power_ok = processor._validate_power_preservation(
            working_data, multilooked_data, "Multilooking", 
            tolerance=0.05,  # 5% tolerance for multilooking
            log_warning=False
        )

        if not power_ok:
            # Compute diagnostic statistics to understand the failure
            import logging
            logger = logging.getLogger(__name__)
            try:
                # Use >= 0 to match Rust multilooking behavior (includes zeros in averaging)
                pre_vals = working_data[np.isfinite(working_data) & (working_data >= 0)]
                post_vals = multilooked_data[np.isfinite(multilooked_data) & (multilooked_data >= 0)]
                if pre_vals.size >= 100 and post_vals.size >= 100:
                    pre_mean = float(np.mean(pre_vals))
                    post_mean = float(np.mean(post_vals))
                    power_ratio = post_mean / pre_mean if pre_mean > 0 else 0.0
                    logger.warning(
                        f"Power preservation check failed for multilooking. "
                        f"Pre-multilook mean: {pre_mean:.6g}, Post-multilook mean: {post_mean:.6g}, "
                        f"Ratio: {power_ratio:.4f} (expected ~1.00 ± 0.05)"
                    )
                else:
                    logger.error(
                        "Power preservation check failed for multilooking. "
                        "Insufficient valid pixels for diagnostics."
                    )
            except Exception as diag_err:
                logger.error(
                    f"Power preservation check failed for multilooking. "
                    f"Diagnostic computation failed: {diag_err}"
                )
            
            # Only rescale if explicitly allowed in options
            allow_power_rescale = getattr(processor, "options", {}).get("allow_power_rescale", False)
            if not allow_power_rescale:
                logger.warning(
                    "Power rescaling disabled (default). Set allow_power_rescale=True to enable automatic correction. "
                    "Proceeding with unscaled data - results may have incorrect radiometry."
                )
            else:
                try:
                    # Use >= 0 to match Rust multilooking behavior
                    pre_vals = working_data[np.isfinite(working_data) & (working_data >= 0)]
                    post_vals = multilooked_data[np.isfinite(multilooked_data) & (multilooked_data >= 0)]
                    if pre_vals.size >= 100 and post_vals.size >= 100:
                        pre_mean = float(np.mean(pre_vals))
                        post_mean = float(np.mean(post_vals))
                        if post_mean > 0:
                            scale = pre_mean / post_mean
                            if np.isfinite(scale) and POWER_SCALE_MIN <= scale <= POWER_SCALE_MAX:
                                multilooked_data = np.asarray(multilooked_data, dtype=np.float32) * scale
                                logger.warning(
                                    f"Multilooking: applied power rescale factor {scale:.4f} to enforce preservation. "
                                    f"This indicates the multilooking algorithm is not properly preserving power."
                                )
                                processor._log(
                                    f"   ⚠️  Multilooking: scaled output by {scale:.4f} to enforce power preservation",
                                    1,
                                )
                                power_ok = processor._validate_power_preservation(
                                    working_data,
                                    multilooked_data,
                                    "Multilooking (after rescale)",
                                    tolerance=POWER_PRESERVATION_TOLERANCE,
                                    log_warning=True,
                                )
                            else:
                                logger.error(
                                    f"Multilooking power scale {scale:.3g} outside safe range "
                                    f"[{POWER_SCALE_MIN}, {POWER_SCALE_MAX}]; skipping rescale. "
                                    f"This indicates a serious issue with the multilooking output."
                                )
                                processor._log(
                                    f"   ❌  Multilooking power scale {scale:.3g} outside safe range; skipping rescale",
                                    2,
                                )
                    else:
                        logger.warning(
                            f"Insufficient valid pixels for power rescale calculation: "
                            f"pre={pre_vals.size}, post={post_vals.size} (need ≥100 each)"
                        )
                except Exception as power_err:
                    logger.error(f"Multilooking power rescale failed: {power_err}")
                    processor._log(
                        f"   ❌  Multilooking power rescale skipped due to error: {power_err}",
                        2,
                    )

        processor._working_data = multilooked_data
        processor._log_array_metrics("Multilooking", multilooked_data)

    except Exception as exc:
        step_duration = time.time() - step_start
        processor.log_step(7, "Multilooking", "error", f"Multilooking failed: {exc}", step_duration)
        raise RuntimeError(f"SCIENTIFIC MODE FAILURE: Multilooking failed. Error: {exc}")
