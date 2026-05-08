#![allow(dead_code)]
//! Merge plan construction: deterministic, purely geometric.
//!
//! The merge plan decomposes the merge operation into a sequence of
//! copy/blend row segments that can be executed efficiently.

use crate::types::{SarError, SarResult, SubSwath};

use super::overlap::OverlapRegion;
use super::assert_subswath_geometry;

/// Weighting strategy for a merge segment
#[derive(Clone, Debug)]
pub enum MergeWeight {
    /// Constant weight applied to every sample in the segment
    Constant(f32),
    /// Weight derived from an overlap region weight matrix
    Overlap {
        overlap_index: usize,
        row: usize,
        col_offset: usize,
        /// If true, use (1 - w) instead of w
        inverse: bool,
    },
}

/// Row segment describing how to copy/blend a portion of a subswath into the output grid.
#[derive(Clone, Debug)]
pub struct MergeRowSegment {
    pub swath_idx: usize,
    pub src_row: usize,
    pub src_col_start: usize,
    pub dst_col_start: usize,
    pub len: usize,
    pub weight: MergeWeight,
}

/// Precomputed merge plan decomposing the merge into copy/blend row segments.
#[derive(Clone, Debug)]
pub struct MergePlan {
    pub rows: usize,
    pub cols: usize,
    pub rows_plan: Vec<Vec<MergeRowSegment>>,
}

/// Build a deterministic merge plan composed of copy/blend row segments.
pub fn build_merge_plan(
    subswaths: &[SubSwath],
    overlap_regions: &[OverlapRegion],
    output_rows: usize,
    output_cols: usize,
) -> SarResult<MergePlan> {
    if output_rows == 0 || output_cols == 0 {
        return Err(SarError::Processing(
            "Output grid dimensions must be positive to build merge plan".to_string(),
        ));
    }

    #[cfg(debug_assertions)]
    for sw in subswaths {
        assert_subswath_geometry(sw);
    }

    let mut rows_plan: Vec<Vec<MergeRowSegment>> = vec![Vec::new(); output_rows];

    // Build swath index lookup
    let mut swath_index: std::collections::HashMap<&str, usize> = std::collections::HashMap::new();
    for (idx, swath) in subswaths.iter().enumerate() {
        swath_index.insert(swath.id.as_str(), idx);
    }

    // Build overlap lookup per swath
    let mut overlaps_by_swath: Vec<Vec<(usize, bool)>> = vec![Vec::new(); subswaths.len()];
    for (idx, overlap) in overlap_regions.iter().enumerate() {
        if let Some(&sw_idx) = swath_index.get(overlap.swath1_id.as_str()) {
            overlaps_by_swath[sw_idx].push((idx, true)); // is_first = true
        }
        if let Some(&sw_idx) = swath_index.get(overlap.swath2_id.as_str()) {
            overlaps_by_swath[sw_idx].push((idx, false)); // is_first = false
        }
    }

    for (swath_idx, swath) in subswaths.iter().enumerate() {
        let src_rows = swath.azimuth_samples;
        let src_cols = swath.range_samples;

        if src_rows == 0 || src_cols == 0 {
            continue;
        }

        let global_row_start = swath.first_line_global;
        let global_col_start = swath.first_sample_global;

        // Respect per-swath valid windows
        let valid_row_start_global = swath
            .valid_first_line
            .unwrap_or(swath.first_line_global)
            .max(swath.first_line_global);
        let valid_row_end_global = swath
            .valid_last_line
            .unwrap_or(swath.last_line_global)
            .min(swath.last_line_global);
        let valid_col_start_global = swath
            .valid_first_sample
            .unwrap_or(swath.first_sample_global)
            .max(swath.first_sample_global);
        let valid_col_end_global = swath
            .valid_last_sample
            .unwrap_or(swath.last_sample_global)
            .min(swath.last_sample_global);

        if valid_row_end_global <= valid_row_start_global
            || valid_col_end_global <= valid_col_start_global
        {
            // DIAGNOSTIC: Log when subswath is skipped due to invalid valid ranges
            log::debug!(
                "Skipping {}: invalid valid ranges (rows: {}..{}, cols: {}..{})",
                swath.id,
                valid_row_start_global,
                valid_row_end_global,
                valid_col_start_global,
                valid_col_end_global
            );
            continue;
        }

        if global_row_start >= output_rows || global_col_start >= output_cols {
            continue;
        }

        // Clamp to actual array extent
        let effective_last_row = swath.last_line_global.min(
            swath
                .first_line_global
                .saturating_add(swath.azimuth_samples),
        );
        let effective_last_col = swath.last_sample_global.min(
            swath
                .first_sample_global
                .saturating_add(swath.range_samples),
        );

        let valid_row_start_local = valid_row_start_global.saturating_sub(swath.first_line_global);
        let valid_row_end_local = valid_row_end_global
            .min(effective_last_row)
            .saturating_sub(swath.first_line_global);

        let valid_col_start_local =
            valid_col_start_global.saturating_sub(swath.first_sample_global);
        let valid_col_end_local = valid_col_end_global
            .min(effective_last_col)
            .saturating_sub(swath.first_sample_global);

        // TARGETED DIAGNOSTIC: Log IW3's valid ranges for gap investigation
        if swath.id == "IW3" {
            log::info!(
                "🔍 IW3 MERGE PLAN: global_row_start={}, valid_row_start_global={}, valid_row_end_global={}, effective_last_row={}",
                global_row_start,
                valid_row_start_global,
                valid_row_end_global,
                effective_last_row
            );
            log::info!(
                "🔍 IW3 MERGE PLAN: valid_row_start_local={}, valid_row_end_local={}, src_rows={}",
                valid_row_start_local,
                valid_row_end_local,
                src_rows
            );
            // Check what dst_row range will be generated
            if valid_row_end_local > valid_row_start_local {
                let min_dst_row = global_row_start + valid_row_start_local;
                let max_dst_row = global_row_start + valid_row_end_local.saturating_sub(1);
                log::info!(
                    "🔍 IW3 MERGE PLAN: Will generate segments for dst_row range [{}, {}] (output_rows={})",
                    min_dst_row,
                    max_dst_row,
                    output_rows
                );
                // Check specifically for gap rows
                let gap_row_start = 12426;
                let gap_row_end = 12470;
                if max_dst_row < gap_row_start {
                    log::warn!(
                        "⚠️  IW3 MERGE PLAN: max_dst_row {} < gap_row_start {} - IW3 will NOT contribute to gap rows!",
                        max_dst_row,
                        gap_row_start
                    );
                } else if min_dst_row > gap_row_end {
                    log::warn!(
                        "⚠️  IW3 MERGE PLAN: min_dst_row {} > gap_row_end {} - IW3 will NOT contribute to gap rows!",
                        min_dst_row,
                        gap_row_end
                    );
                } else {
                    log::info!(
                        "✅ IW3 MERGE PLAN: Will contribute to gap rows [{}, {}]",
                        gap_row_start.max(min_dst_row),
                        gap_row_end.min(max_dst_row)
                    );
                }
            }
        }

        if valid_row_end_local <= valid_row_start_local
            || valid_col_end_local <= valid_col_start_local
        {
            continue;
        }

        let mut col_limit = valid_col_end_local
            .saturating_sub(valid_col_start_local)
            .min(src_cols.saturating_sub(valid_col_start_local));
        if global_col_start + valid_col_start_local + col_limit > output_cols {
            col_limit = output_cols.saturating_sub(global_col_start + valid_col_start_local);
        }
        if col_limit == 0 {
            continue;
        }

        // Track segment creation for IW3 in gap region
        let mut iw3_gap_segments = 0usize;
        for local_row in valid_row_start_local..valid_row_end_local {
            let dst_row = global_row_start + local_row;
            if dst_row >= output_rows {
                if swath.id == "IW3" && dst_row < output_rows + 10 {
                    log::debug!(
                        "🔍 IW3: Skipping local_row {} -> dst_row {} (>= output_rows {})",
                        local_row,
                        dst_row,
                        output_rows
                    );
                }
                continue;
            }

            // Track IW3 segments in gap region
            if swath.id == "IW3" && dst_row >= 12426 && dst_row <= 12470 {
                iw3_gap_segments += 1;
            }

            let mut segments = vec![MergeRowSegment {
                swath_idx,
                src_row: local_row,
                src_col_start: valid_col_start_local,
                dst_col_start: global_col_start + valid_col_start_local,
                len: col_limit,
                weight: MergeWeight::Constant(1.0),
            }];

            // Apply overlap modifications
            for &(overlap_idx, is_first) in &overlaps_by_swath[swath_idx] {
                segments = apply_overlap_to_segments(
                    segments,
                    overlap_idx,
                    dst_row,
                    is_first,
                    overlap_regions,
                    subswaths,
                    output_cols,
                )?;
            }

            segments.retain(|seg| seg.len > 0);
            if segments.is_empty() {
                continue;
            }

            segments.sort_by_key(|seg| seg.dst_col_start);
            rows_plan[dst_row].extend(segments);
        }

        // Report IW3 gap segment count
        if swath.id == "IW3" && iw3_gap_segments > 0 {
            log::info!(
                "✅ IW3 MERGE PLAN: Created {} segments for gap rows [12426, 12470]",
                iw3_gap_segments
            );
        } else if swath.id == "IW3" {
            log::warn!(
                "⚠️  IW3 MERGE PLAN: Created 0 segments for gap rows [12426, 12470] - this explains the gaps!"
            );
        }
    }

    Ok(MergePlan {
        rows: output_rows,
        cols: output_cols,
        rows_plan,
    })
}

/// Split constant-weight segments so overlap regions are represented explicitly.
fn apply_overlap_to_segments(
    segments: Vec<MergeRowSegment>,
    overlap_idx: usize,
    dst_row: usize,
    is_first: bool,
    overlap_regions: &[OverlapRegion],
    subswaths: &[SubSwath],
    output_cols: usize,
) -> SarResult<Vec<MergeRowSegment>> {
    let overlap = overlap_regions.get(overlap_idx).ok_or_else(|| {
        SarError::Processing(format!(
            "Invalid overlap index {} while building merge plan",
            overlap_idx
        ))
    })?;

    if segments.is_empty() {
        return Ok(segments);
    }

    // Check if this row is in the overlap's azimuth range
    if dst_row < overlap.azimuth_start || dst_row >= overlap.azimuth_end {
        return Ok(segments);
    }

    let row_in_overlap = dst_row.checked_sub(overlap.azimuth_start).ok_or_else(|| {
        SarError::Processing(
            "Overlap azimuth start exceeds destination row during merge plan".to_string(),
        )
    })?;

    if row_in_overlap >= overlap.weights.nrows() {
        return Ok(segments);
    }

    let swath_idx = segments.first().map(|s| s.swath_idx).unwrap_or_default();
    #[cfg(debug_assertions)]
    if let Some(sw) = subswaths.get(swath_idx) {
        assert_subswath_geometry(sw);
    }
    let swath_first_sample_global = subswaths
        .get(swath_idx)
        .map(|s| s.first_sample_global)
        .unwrap_or(0);

    let (range_start_local, range_end_local) = if is_first {
        (overlap.swath1_range_start, overlap.swath1_range_end)
    } else {
        (overlap.swath2_range_start, overlap.swath2_range_end)
    };

    let range_start = swath_first_sample_global.saturating_add(range_start_local);
    let range_end = swath_first_sample_global.saturating_add(range_end_local);

    if range_end <= range_start {
        return Ok(segments);
    }

    let mut result = Vec::new();

    for seg in segments {
        if seg.len == 0 {
            continue;
        }

        let seg_start = seg.dst_col_start;
        let seg_end = seg_start.saturating_add(seg.len);

        // No overlap with this segment
        if seg_end <= range_start || seg_start >= range_end {
            result.push(seg);
            continue;
        }

        match seg.weight {
            MergeWeight::Constant(weight) => {
                let mut cursor_dst = seg_start;
                let mut cursor_src = seg.src_col_start;
                let segment_end = seg_end.min(output_cols);

                // Left (non-overlap) portion
                if cursor_dst < range_start {
                    let left_end = range_start.min(segment_end);
                    let left_len = left_end.saturating_sub(cursor_dst);
                    if left_len > 0 {
                        result.push(MergeRowSegment {
                            swath_idx: seg.swath_idx,
                            src_row: seg.src_row,
                            src_col_start: cursor_src,
                            dst_col_start: cursor_dst,
                            len: left_len,
                            weight: MergeWeight::Constant(weight),
                        });
                        cursor_dst += left_len;
                        cursor_src += left_len;
                    }
                }

                // Overlap portion
                let overlap_start = cursor_dst.max(range_start);
                let overlap_end = segment_end.min(range_end);

                if overlap_end > overlap_start {
                    let base_col_offset =
                        overlap_start.checked_sub(range_start).ok_or_else(|| {
                            SarError::Processing(
                                "Overlap range underflow during merge plan".to_string(),
                            )
                        })?;
                    let weight_cols = overlap.weights.ncols();

                    if base_col_offset < weight_cols {
                        let desired_overlap_len = overlap_end - overlap_start;
                        let available = weight_cols - base_col_offset;

                        // FIX for bug B: Single check, fail if insufficient
                        if available < desired_overlap_len {
                            return Err(SarError::Processing(format!(
                                "Overlap weight width {} too short for desired length {} at row {} (swath {})",
                                available,
                                desired_overlap_len,
                                dst_row,
                                subswaths[seg.swath_idx].id
                            )));
                        }

                        let src_offset = overlap_start.saturating_sub(seg_start);
                        result.push(MergeRowSegment {
                            swath_idx: seg.swath_idx,
                            src_row: seg.src_row,
                            src_col_start: seg.src_col_start + src_offset,
                            dst_col_start: overlap_start,
                            len: desired_overlap_len,
                            weight: MergeWeight::Overlap {
                                overlap_index: overlap_idx,
                                row: row_in_overlap,
                                col_offset: base_col_offset,
                                inverse: !is_first,
                            },
                        });

                        cursor_dst = overlap_start + desired_overlap_len;
                        cursor_src = seg.src_col_start + (cursor_dst - seg_start);
                    }
                }

                // Right (non-overlap) portion
                if cursor_dst < segment_end {
                    let remaining_len = segment_end - cursor_dst;
                    if remaining_len > 0 {
                        result.push(MergeRowSegment {
                            swath_idx: seg.swath_idx,
                            src_row: seg.src_row,
                            src_col_start: cursor_src,
                            dst_col_start: cursor_dst,
                            len: remaining_len,
                            weight: MergeWeight::Constant(weight),
                        });
                    }
                }
            }
            _ => {
                // Already has overlap weight, pass through
                result.push(seg);
            }
        }
    }

    Ok(result)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn simple_swath(
        id: &str,
        first_line: usize,
        first_sample: usize,
        rows: usize,
        cols: usize,
    ) -> SubSwath {
        SubSwath {
            id: id.to_string(),
            burst_count: 1,
            lines_per_burst: rows,
            range_samples: cols,
            azimuth_samples: rows,
            first_line_global: first_line,
            last_line_global: first_line + rows,
            first_sample_global: first_sample,
            last_sample_global: first_sample + cols,
            full_range_samples: cols,
            valid_first_line: Some(first_line),
            valid_last_line: Some(first_line + rows),
            valid_first_sample: Some(first_sample),
            valid_last_sample: Some(first_sample + cols),
            range_pixel_spacing: 1.0,
            azimuth_pixel_spacing: 1.0,
            slant_range_time: 0.0,
            burst_duration: rows as f64,
            near_range_m: 0.0,
            prf_hz: Some(1000.0),
            dc_polynomial: Some(vec![0.0]),
            azimuth_time_interval: Some(0.001),
            dc_polynomial_t0: Some(0.0),
            fm_rate_estimates: None,
        }
    }

    #[test]
    fn plan_single_swath() {
        let sw = simple_swath("IW1", 0, 0, 10, 20);
        let plan = build_merge_plan(&[sw], &[], 10, 20).unwrap();

        assert_eq!(plan.rows, 10);
        assert_eq!(plan.cols, 20);

        // Each row should have exactly one segment covering all columns
        for row_plan in &plan.rows_plan {
            assert_eq!(row_plan.len(), 1);
            assert_eq!(row_plan[0].len, 20);
        }
    }

    #[test]
    fn plan_fails_on_zero_dims() {
        let sw = simple_swath("IW1", 0, 0, 10, 20);

        let result = build_merge_plan(&[sw.clone()], &[], 0, 20);
        assert!(result.is_err());

        let result = build_merge_plan(&[sw], &[], 10, 0);
        assert!(result.is_err());
    }
}
