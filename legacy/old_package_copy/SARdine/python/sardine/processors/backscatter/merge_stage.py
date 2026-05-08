"""
Merge stage helper for the backscatter pipeline.
"""

import logging
import time
import inspect
import numpy as np
import sardine

from .deburst_merge import subswath_sort_key, VALID_IW_SUBSWATHS


def perform_merge(processor) -> None:
    """Combine calibrated subswaths after deburst/calibration."""
    logger = logging.getLogger(__name__)
    self = processor
    calibrated_subswaths = getattr(self, "_calibrated_subswaths", {})

    if not calibrated_subswaths:
        raise RuntimeError("No calibrated subswaths available for merge")

    self.announce_step(6, "Expert IW Merge", "Combining calibrated subswaths")
    step_start = time.time()

    try:
        if len(calibrated_subswaths) < 2:
            print("   ℹ️  Only one subswath available after calibration; promoting to working data")
            # FIXED: Promote single subswath data to _working_data so later stages can proceed
            single_sw_name = list(calibrated_subswaths.keys())[0]
            single_sw_data = calibrated_subswaths[single_sw_name]
            self._working_data = np.asarray(single_sw_data, dtype=np.float32)
            self._used_subswaths_for_merge = [single_sw_name]
            logger.info(f"Single subswath {single_sw_name} promoted to working data: {self._working_data.shape}")
            step_duration = time.time() - step_start
            self.log_step(
                6,
                "Expert IW Merge",
                "skipped",
                f"Single subswath {single_sw_name} promoted (no merge needed)",
                step_duration,
            )
            return

        print(f"🔗 EXPERT IW MERGE: Attempting merge of {len(calibrated_subswaths)} calibrated subswaths")

        # Validate and sort subswaths
        ordered_for_merge = sorted(calibrated_subswaths.keys(), key=subswath_sort_key)
        
        # Validate subswath order and completeness
        available_subswaths = {sw.upper() for sw in ordered_for_merge}
        missing = VALID_IW_SUBSWATHS - available_subswaths
        if missing:
            # FIXED: Allow partial scenes with 2 subswaths instead of requiring all 3
            # This enables processing of cropped/partial acquisitions
            logger.info(
                f"Partial coverage: missing subswaths {missing}. "
                f"Available: {available_subswaths}. Proceeding with available subswaths."
            )
        
        # Validate order is correct
        expected_order = ["IW1", "IW2", "IW3"]
        actual_order = [sw.upper() for sw in ordered_for_merge]
        if actual_order != expected_order[:len(actual_order)]:
            logger.warning(
                f"Subswath order may be incorrect: {actual_order}. "
                f"Expected: {expected_order[:len(actual_order)]}"
            )
        
        # Validate and convert subswath data
        merge_inputs = {}
        merge_inputs_upper = {}
        for sw in ordered_for_merge:
            sw_data = calibrated_subswaths.get(sw)
            if sw_data is None:
                raise RuntimeError(f"Missing calibrated data for subswath: {sw}")
            
            # OPTIMIZATION #2: Avoid redundant array copy if already float32 ndarray
            try:
                if isinstance(sw_data, np.ndarray) and sw_data.dtype == np.float32:
                    array = sw_data  # Reuse existing array (zero-copy)
                else:
                    array = np.asarray(sw_data, dtype=np.float32)
            except Exception as e:
                raise RuntimeError(f"Failed to convert {sw} data to array: {e}") from e
            
            # Validate array structure
            if array.ndim != 2:
                raise RuntimeError(
                    f"Invalid dimensions for {sw}: {array.ndim}D. Expected 2D array."
                )
            
            if array.size == 0:
                raise RuntimeError(f"Empty array for subswath: {sw}")
            
            # OPTIMIZATION #5: Single-pass validation - compute finite mask once, reuse for all stats
            finite_mask = np.isfinite(array)
            finite_count = np.sum(finite_mask)
            finite_fraction = finite_count / array.size if array.size > 0 else 0.0
            if finite_count == 0:
                raise RuntimeError(
                    f"Subswath {sw} contains no finite values - cannot merge"
                )
            
            # Enhanced validation: check finite fraction
            if finite_fraction < 0.5:
                logger.warning(
                    f"⚠️  Subswath {sw} has low finite value fraction: {finite_fraction*100:.1f}% "
                    f"({finite_count}/{array.size} pixels). Merge may produce many NaN values."
                )
            elif finite_fraction < 0.9:
                logger.info(
                    f"ℹ️  Subswath {sw} finite value fraction: {finite_fraction*100:.1f}% "
                    f"({finite_count}/{array.size} pixels)"
                )
            
            # Validate data range (should be positive for backscatter)
            # OPTIMIZATION #5 continued: Reuse finite_mask instead of recomputing
            if finite_count > 0:
                finite_data = array[finite_mask]
                min_val = float(np.min(finite_data))
                max_val = float(np.max(finite_data))
                mean_val = float(np.mean(finite_data))
                
                if min_val < 0.0:
                    logger.warning(
                        f"⚠️  Subswath {sw} contains negative values: min={min_val:.4f} "
                        f"(may be expected for some processing stages)"
                    )
                
                logger.debug(
                    f"📊 Subswath {sw} validation: shape={array.shape}, "
                    f"finite={finite_fraction*100:.1f}%, "
                    f"range=[{min_val:.4f}, {max_val:.4f}], mean={mean_val:.4f}"
                )
            
            merge_inputs[sw] = array
            merge_inputs_upper[sw.upper()] = array
        
        # Get deburst overrides with proper validation
        deburst_overrides_raw = getattr(self, "_deburst_overrides", None)
        if deburst_overrides_raw is None:
            print("   ⚠️ Deburst overrides missing; merge will use annotation timing")
            deburst_overrides = None
        elif not isinstance(deburst_overrides_raw, dict):
            logger.warning(
                f"Invalid deburst_overrides type: {type(deburst_overrides_raw)}. "
                f"Expected dict. Using annotation timing instead."
            )
            deburst_overrides = None
        elif not deburst_overrides_raw:
            print("   ⚠️ Deburst overrides empty; merge will use annotation timing")
            deburst_overrides = None
        else:
            # Validate deburst_overrides structure
            try:
                keys = list(deburst_overrides_raw.keys())
                print(f"   🛰️ Deburst overrides available for {sorted(keys)}")
                
                # Validate each override entry has required fields
                for swath_id, override in deburst_overrides_raw.items():
                    if not isinstance(override, dict):
                        logger.warning(f"Invalid override type for {swath_id}: {type(override)}")
                        continue
                    
                    required_fields = ["burst_timing", "row_provenance"]
                    missing_fields = [f for f in required_fields if f not in override or not override[f]]
                    if missing_fields:
                        logger.warning(f"Override for {swath_id} missing fields: {missing_fields}")
                
                deburst_overrides = deburst_overrides_raw
            except Exception as e:
                logger.warning(
                    f"Failed to validate deburst_overrides: {e}. "
                    f"Using annotation timing instead."
                )
                deburst_overrides = None
        
        # Validate reader is available
        if not hasattr(self, "reader") or self.reader is None:
            raise RuntimeError(
                "Reader not available for merge. "
                "This indicates Step 1 (metadata extraction) failed."
            )
        
        # Validate reader has required methods
        if not hasattr(self.reader, "get_cached_metadata"):
            raise RuntimeError(
                "Reader missing get_cached_metadata() method. "
                "Merge requires cached metadata for geometry."
            )

        supports_deburst_overrides = False
        try:
            sig = inspect.signature(sardine.merge_subswaths_cached)
            supports_deburst_overrides = "deburst_overrides" in sig.parameters
        except (TypeError, ValueError):
            supports_deburst_overrides = False

        if {"IW1", "IW2", "IW3"}.issubset(set(merge_inputs_upper.keys())):
            if supports_deburst_overrides:
                merged_result = sardine.merge_subswaths_cached(
                    merge_inputs_upper["IW1"],
                    merge_inputs_upper["IW2"],
                    merge_inputs_upper["IW3"],
                    self.reader,
                    self.polarization,
                    deburst_overrides,
                )
            else:
                merged_result = sardine.merge_subswaths_cached(
                    merge_inputs_upper["IW1"],
                    merge_inputs_upper["IW2"],
                    merge_inputs_upper["IW3"],
                    self.reader,
                    self.polarization,
                )
        else:
            merged_args = [merge_inputs, self.polarization, self.reader]
            if supports_deburst_overrides:
                merged_args.append(None)
                merged_args.append(deburst_overrides)
            merged_result = sardine.topsar_merge_cached(*merged_args)

        # Validate merge result structure
        if not isinstance(merged_result, dict):
            raise RuntimeError(
                f"Invalid merge result type: {type(merged_result)}. "
                f"Expected dict with 'data' key."
            )
        
        # Extract and validate merged data
        merged_data = merged_result.get("data")
        if merged_data is None:
            available_keys = list(merged_result.keys())
            raise RuntimeError(
                f"Merge result missing 'data' key. "
                f"Available keys: {available_keys}"
            )
        
        # Validate data is a valid array
        try:
            merged_array = np.asarray(merged_data, dtype=np.float32)
        except Exception as e:
            raise RuntimeError(
                f"Failed to convert merge data to array: {e}. "
                f"Data type: {type(merged_data)}"
            ) from e
        
        # Validate array dimensions
        if merged_array.ndim != 2:
            raise RuntimeError(
                f"Invalid merge output dimensions: {merged_array.ndim}D. "
                f"Expected 2D array (lines x samples)."
            )
        
        if merged_array.size == 0:
            raise RuntimeError("Merge produced empty array")
        
        # Validate array values are reasonable
        finite_count = np.sum(np.isfinite(merged_array))
        finite_fraction = finite_count / merged_array.size if merged_array.size > 0 else 0.0
        nan_count = np.sum(np.isnan(merged_array))
        nan_fraction = nan_count / merged_array.size if merged_array.size > 0 else 0.0
        
        if finite_fraction < 0.5:
            raise RuntimeError(
                f"Merge output has too many invalid values: "
                f"{finite_fraction*100:.1f}% finite (minimum 50% required). "
                f"NaN count: {nan_count} ({nan_fraction*100:.1f}%)"
            )
        
        # Enhanced validation: check NaN distribution by subswath region
        if finite_fraction < 0.9:
            logger.warning(
                f"⚠️  Merge output has elevated NaN percentage: {nan_fraction*100:.1f}% "
                f"({nan_count}/{merged_array.size} pixels). "
                f"This may indicate segment execution issues."
            )
            
            # Analyze NaN distribution by column regions (approximate subswath positions)
            cols = merged_array.shape[1]
            if cols > 0:
                # Approximate IW1, IW2, IW3 column ranges (these are approximate)
                iw1_end = cols // 3
                iw2_start = iw1_end
                iw2_end = 2 * cols // 3
                iw3_start = iw2_end
                
                iw1_nan = np.sum(np.isnan(merged_array[:, :iw1_end]))
                iw2_nan = np.sum(np.isnan(merged_array[:, iw2_start:iw2_end]))
                iw3_nan = np.sum(np.isnan(merged_array[:, iw3_start:]))
                
                iw1_total = merged_array.shape[0] * iw1_end
                iw2_total = merged_array.shape[0] * (iw2_end - iw2_start)
                iw3_total = merged_array.shape[0] * (cols - iw3_start)
                
                logger.info(
                    f"📊 NaN distribution (approximate): "
                    f"IW1={100.0*iw1_nan/iw1_total:.1f}%, "
                    f"IW2={100.0*iw2_nan/iw2_total:.1f}%, "
                    f"IW3={100.0*iw3_nan/iw3_total:.1f}%"
                )
        
        # Validate merge output quality
        coverage = merged_result.get("coverage_percent")
        if coverage is not None:
            coverage_val = float(coverage)
            if coverage_val < 50.0:
                raise RuntimeError(
                    f"Merge produced insufficient coverage: {coverage_val:.1f}%. "
                    f"Minimum 50% required for scientific processing."
                )
            elif coverage_val < 80.0:
                logger.warning(
                    f"Low merge coverage: {coverage_val:.1f}%. "
                    f"May indicate subswath alignment issues."
                )
        
        # Validate hit count if available
        hit_count = merged_result.get("hit_count")
        if hit_count is not None:
            hit_array = np.asarray(hit_count)
            zero_hit_count = np.sum(hit_array == 0)
            if zero_hit_count > merged_array.size * 0.2:
                logger.warning(
                    f"High number of zero-hit pixels: {zero_hit_count} "
                    f"({zero_hit_count/merged_array.size*100:.1f}%)"
                )

        # Seam QA: check for gaps and amplitude jumps across subswath seams
        seam_report = {}
        try:
            rows, cols = merged_array.shape
            valid_mask = np.isfinite(merged_array) & (merged_array > 0)
            col_valid_frac = np.sum(valid_mask, axis=0) / float(rows)

            low_cols = col_valid_frac < 0.05
            if low_cols.any():
                low_int = low_cols.astype(np.int32)
                transitions = np.diff(np.concatenate(([0], low_int, [0])))
                starts = np.where(transitions == 1)[0]
                ends = np.where(transitions == -1)[0]
                run_lengths = ends - starts
                longest_gap = int(run_lengths.max()) if run_lengths.size else 0
                longest_gap_pct = longest_gap / float(cols) * 100.0 if cols else 0.0
                seam_report["longest_low_valid_run_cols"] = longest_gap
                seam_report["longest_low_valid_run_pct"] = longest_gap_pct
                if longest_gap_pct > 20.0:
                    logger.warning(
                        f"⚠️  Large gap of low-valid columns detected: {longest_gap} cols "
                        f"({longest_gap_pct:.1f}% of width)."
                    )

            # Evaluate amplitude continuity at approximate seam boundaries
            seam_stats = []
            num_swaths = len(ordered_for_merge)
            if num_swaths >= 2:
                step = cols / float(num_swaths)
                window = max(32, int(step * 0.05))

                def _safe_db(arr: np.ndarray) -> np.ndarray:
                    return 10.0 * np.log10(np.maximum(arr, 1e-12))

                for idx in range(num_swaths - 1):
                    boundary = int(round(step * (idx + 1)))
                    left_slice = merged_array[:, max(0, boundary - window):boundary]
                    right_slice = merged_array[:, boundary:min(cols, boundary + window)]

                    left_valid = left_slice[np.isfinite(left_slice) & (left_slice > 0)]
                    right_valid = right_slice[np.isfinite(right_slice) & (right_slice > 0)]

                    stat = {
                        "boundary_col": boundary,
                        "left_count": int(left_valid.size),
                        "right_count": int(right_valid.size),
                    }

                    if left_valid.size and right_valid.size:
                        left_db = float(np.nanmedian(_safe_db(left_valid)))
                        right_db = float(np.nanmedian(_safe_db(right_valid)))
                        jump_db = abs(left_db - right_db)
                        stat["jump_db"] = jump_db
                        seam_stats.append(stat)
                        if jump_db > 3.0:
                            logger.warning(
                                f"⚠️  Amplitude jump across seam {idx+1}: {jump_db:.2f} dB "
                                f"(boundary col {boundary}, window {window} cols)"
                            )
                    else:
                        stat["jump_db"] = None
                        seam_stats.append(stat)
                        logger.warning(
                            f"⚠️  Missing valid samples to assess seam {idx+1} near column {boundary}."
                        )

            seam_report["seam_stats"] = seam_stats
            seam_report["subswaths"] = ordered_for_merge

            if hasattr(self, "_record_validation"):
                for idx, stat in enumerate(seam_stats, start=1):
                    jump = stat.get("jump_db")
                    msg = (
                        f"Seam {idx} jump {jump:.2f} dB" if jump is not None else
                        "Seam jump unavailable (insufficient data)"
                    )
                    self._record_validation(
                        "merge",
                        f"seam_{idx}_jump_db",
                        value=jump if jump is not None else float("nan"),
                        expected="<=3 dB",
                        passed=(jump is not None and jump <= 3.0),
                        severity="warning" if jump is None or jump > 3.0 else "info",
                        message=msg,
                    )
        except Exception as seam_exc:
            logger.debug(f"Seam QA skipped due to error: {seam_exc}")

        if seam_report:
            self._merge_seam_report = seam_report

        self._working_data = merged_array
        
        # CRITICAL FIX: Update metadata dimensions to reflect merged dimensions
        # This ensures downstream stages (terrain correction, etc.) use correct dimensions
        merged_height, merged_width = merged_array.shape
        
        # CRITICAL FIX: Set native dimensions from merged output
        # This is the proper place to set these, after all subswaths are merged.
        # The per-subswath processing should NOT set these as they would overwrite each other.
        self._native_azimuth_lines = merged_height
        self._native_range_samples = merged_width
        logger.info(
            f"📐 Set native dimensions from merged output: "
            f"{merged_height} azimuth lines × {merged_width} range samples"
        )
        
        # Validate merged dimensions are reasonable
        # FIXED: Reduced minimum to support multilooked and cropped data
        min_dim = 100  # Allow small arrays for heavy multilooking or test crops
        if merged_width < min_dim or merged_height < min_dim:
            raise RuntimeError(
                f"Merged array dimensions ({merged_height}x{merged_width}) are too small. "
                f"Expected at least {min_dim}x{min_dim}. Check subswath data."
            )
        if merged_width > 200000 or merged_height > 200000:
            raise RuntimeError(
                f"Merged array dimensions ({merged_height}x{merged_width}) are unexpectedly large. "
                f"Expected at most 200000x200000. Check subswath alignment."
            )
        
        if hasattr(self, 'metadata') and isinstance(self.metadata, dict):
            # Update range_samples and azimuth_samples to merged dimensions
            old_range = self.metadata.get('range_samples', 'unknown')
            old_azimuth = self.metadata.get('azimuth_samples', 'unknown')
            self.metadata['range_samples'] = merged_width
            self.metadata['azimuth_samples'] = merged_height
            logger.info(
                f"📊 Updated metadata dimensions: "
                f"range_samples {old_range} → {merged_width}, "
                f"azimuth_samples {old_azimuth} → {merged_height}"
            )
        
        # Also update cached metadata if available
        if hasattr(self, 'reader') and self.reader:
            try:
                cached_meta = getattr(self.reader, 'get_cached_metadata', None)
                if cached_meta and callable(cached_meta):
                    cached = cached_meta()
                    if isinstance(cached, dict):
                        cached['range_samples'] = merged_width
                        cached['azimuth_samples'] = merged_height
                        logger.debug(f"Updated cached metadata dimensions to {merged_height}x{merged_width}")
            except Exception as e:
                logger.debug(f"Could not update cached metadata dimensions: {e}")
        
        self._used_subswaths_for_merge = ordered_for_merge

        if isinstance(merged_result, dict):
            geo_transform = merged_result.get("geo_transform") or merged_result.get("transform")
            if geo_transform is not None:
                self.geo_transform = geo_transform

        if hasattr(merged_result, "get") and merged_result.get("valid_fraction_bounds") is not None:
            bounds = merged_result.get("valid_fraction_bounds")
            try:
                self._merge_valid_fraction_bounds = tuple(float(v) for v in bounds)
            except Exception:
                pass

        step_duration = time.time() - step_start
        self.log_step(
            6,
            "Expert IW Merge",
            "success",
            f"Merged {len(ordered_for_merge)} subswaths",
            step_duration,
        )
        self._record_stage_timing("Expert IW Merge", step_duration)

    except Exception as exc:
        step_duration = time.time() - step_start
        self.log_step(6, "Expert IW Merge", "error", f"Merge failed: {exc}", step_duration)
        raise
