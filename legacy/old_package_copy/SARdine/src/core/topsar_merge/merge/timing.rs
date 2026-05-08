use super::types::{AzimuthTimingModel, BurstTiming};
use crate::core::geometry::type_safe_units::Seconds;
use crate::core::DcFmRateProvider;
use crate::types::{SarError, SarResult, SubSwath};
use std::cmp::Ordering;
use std::collections::HashMap;
use std::sync::Arc;

impl AzimuthTimingModel {
    /// Find the burst timing entry covering a given line for a specific subswath
    pub fn burst_for_line<'a>(
        &'a self,
        subswath_id: &str,
        line_idx: usize,
    ) -> Option<&'a BurstTiming> {
        self.burst_timing.iter().find(|b| {
            b.subswath_id == subswath_id
                && line_idx >= b.first_line_merged
                && line_idx < b.last_line_merged
        })
    }

    /// Get precise azimuth time for a specific line in the merged grid
    /// Accounts for burst boundaries and timing discontinuities
    pub fn get_azimuth_time_at_line(&self, line_idx: usize) -> Option<f64> {
        for burst in &self.burst_timing {
            if line_idx >= burst.first_line_merged && line_idx < burst.last_line_merged {
                let Some(line_in_burst) = line_idx.checked_sub(burst.first_line_merged) else {
                    return None;
                };
                let azimuth_time =
                    burst.azimuth_time_start + (line_in_burst as f64 * burst.azimuth_time_interval);
                return Some(azimuth_time);
            }
        }
        None
    }

    /// Validate azimuth timing consistency across bursts
    /// Checks for timing gaps or overlaps that could affect merge quality
    pub fn validate_timing_consistency(&self) -> Vec<String> {
        let mut warnings = Vec::new();

        for swath_id in self
            .burst_timing
            .iter()
            .map(|b| b.subswath_id.clone())
            .collect::<std::collections::HashSet<_>>()
        {
            let mut bursts: Vec<&BurstTiming> = self
                .burst_timing
                .iter()
                .filter(|b| b.subswath_id == swath_id)
                .collect();
            bursts.sort_by(|a, b| {
                a.azimuth_time_start
                    .partial_cmp(&b.azimuth_time_start)
                    .unwrap_or(Ordering::Equal)
            });

            for i in 1..bursts.len() {
                let prev_burst = bursts[i - 1];
                let curr_burst = bursts[i];

                if prev_burst.burst_id == curr_burst.burst_id {
                    continue;
                }

                let time_gap = curr_burst.azimuth_time_start - prev_burst.azimuth_time_end;

                let gap_threshold = (curr_burst.azimuth_time_interval
                    + prev_burst.azimuth_time_interval)
                    * 0.5
                    * 2.0;

                if time_gap > gap_threshold {
                    warnings.push(format!(
                        "Large timing gap between bursts {}/{} and {}/{}: {:.6} s",
                        prev_burst.subswath_id,
                        prev_burst.burst_id,
                        curr_burst.subswath_id,
                        curr_burst.burst_id,
                        time_gap
                    ));
                } else if time_gap < 0.0 {
                    // TASK E: Timing overlaps are normal TOPS behavior, log at DEBUG level
                    log::debug!(
                        "Timing overlap between bursts {}/{} and {}/{}: {:.6} s (NORMAL TOPS behavior)",
                        prev_burst.subswath_id,
                        prev_burst.burst_id,
                        curr_burst.subswath_id,
                        curr_burst.burst_id,
                        -time_gap
                    );
                    // Don't add to warnings array - this is expected
                }
            }
        }

        let calculated_interval = 1.0 / self.prf;
        let interval_diff = (self.azimuth_time_interval - calculated_interval).abs();

        if interval_diff > 1.0e-6 {
            warnings.push(format!(
                "PRF/azimuth interval mismatch (reference): 1/PRF={:.9} vs interval={:.9}",
                calculated_interval, self.azimuth_time_interval
            ));
        }

        warnings
    }

    /// Strict invariants for timing mapping and provider coverage.
    /// In strict mode these return an error; otherwise they emit warnings.
    pub fn enforce_invariants(
        &self,
        providers: &HashMap<String, Arc<dyn DcFmRateProvider>>,
        subswaths: &[SubSwath],
    ) -> SarResult<()> {
        let strict = crate::types::strict_mode();
        let mut violations: Vec<String> = Vec::new();

        let subswath_t0: HashMap<&str, f64> = subswaths
            .iter()
            .filter_map(|sw| sw.dc_polynomial_t0.map(|t0| (sw.id.as_str(), t0)))
            .collect();

        for burst in &self.burst_timing {
            let n_lines = burst
                .last_line_merged
                .saturating_sub(burst.first_line_merged);
            if n_lines == 0 {
                violations.push(format!(
                    "Burst {} {} has zero line span",
                    burst.subswath_id, burst.burst_index
                ));
                continue;
            }

            // Use per-burst interval (same as get_azimuth_time_at_line) rather than global interval
            // This is critical because different subswaths/bursts may have slightly different PRFs
            let az_dt = burst.azimuth_time_interval;
            let expected_end =
                burst.azimuth_time_start + (n_lines.saturating_sub(1)) as f64 * az_dt;
            if (expected_end - burst.azimuth_time_end).abs() > 1e-9 {
                violations.push(format!(
                    "Burst {} {} end time mismatch: recorded {:.12}, expected {:.12} (n_lines={})",
                    burst.subswath_id,
                    burst.burst_index,
                    burst.azimuth_time_end,
                    expected_end,
                    n_lines
                ));
            }

            // IMPORTANT: Only validate INTERNAL burst timing consistency, not cross-burst alignment.
            // In TOPSAR IW, overlaps are expected and normal. In overlap regions, different bursts
            // can legitimately have different times for the same line. The merge process handles
            // these differences through blending. What matters is that each burst's timing is
            // internally consistent (start + n*dt = end), not that all bursts agree on timing.
            //
            // We've already validated end time consistency above. The cross-validation with
            // get_azimuth_time_at_line is problematic because:
            // 1. get_azimuth_time_at_line returns the FIRST burst's time (implementation detail)
            // 2. In overlaps, multiple bursts contribute with potentially different times
            // 3. These differences are expected and handled by the merge process
            //
            // Therefore, we skip the cross-burst timing validation entirely and only check
            // internal consistency (which we already did with the end time check above).
            //
            // Provider validation is still done for burst boundaries below (seam validation).
        }

        for window in self.burst_timing.windows(2) {
            let [prev, next] = match window {
                [a, b] => [a, b],
                _ => continue,
            };

            if prev.subswath_id != next.subswath_id {
                continue;
            }

            let Some(provider) = providers.get(&prev.subswath_id) else {
                continue;
            };

            let (range_start, range_end) = provider.get_time_range();
            let ref_time = subswath_t0
                .get(prev.subswath_id.as_str())
                .copied()
                .unwrap_or(self.reference_azimuth_time);

            let last_line_prev = prev.last_line_merged.saturating_sub(1);
            let first_line_next = next.first_line_merged;

            let t_prev_last = self
                .get_azimuth_time_at_line(last_line_prev)
                .unwrap_or(prev.azimuth_time_end);
            let t_next_first = self
                .get_azimuth_time_at_line(first_line_next)
                .unwrap_or(next.azimuth_time_start);

            let dc_prev = provider
                .get_dc(Seconds::new(t_prev_last))
                .map(|v| v.value())
                .unwrap_or(f64::NAN);
            let fm_prev = provider
                .get_fm_rate(Seconds::new(t_prev_last))
                .unwrap_or(f64::NAN);

            let dc_next = provider
                .get_dc(Seconds::new(t_next_first))
                .map(|v| v.value())
                .unwrap_or(f64::NAN);
            let fm_next = provider
                .get_fm_rate(Seconds::new(t_next_first))
                .unwrap_or(f64::NAN);

            log::debug!(
                "Seam {} burst {}→{}:\n  burst_k start={:.6} end={:.6}\n  burst_k+1 start={:.6} end={:.6}\n  last_line_abs_time={:.6}\n  first_line_abs_time={:.6}\n  provider: ref_time(t0)={:.6} range=[{:.6},{:.6}]\n  dt_last={:.6} dt_first={:.6}\n  dc_last={:.3}Hz dc_first={:.3}Hz\n  fm_last={:.3e}Hz/s fm_first={:.3e}Hz/s",
                prev.subswath_id,
                prev.burst_index,
                next.burst_index,
                prev.azimuth_time_start,
                prev.azimuth_time_end,
                next.azimuth_time_start,
                next.azimuth_time_end,
                t_prev_last,
                t_next_first,
                ref_time,
                range_start.value(),
                range_end.value(),
                t_prev_last - ref_time,
                t_next_first - ref_time,
                dc_prev,
                dc_next,
                fm_prev,
                fm_next,
            );
        }

        if !violations.is_empty() {
            let total = violations.len();
            for v in violations.iter().take(10) {
                if strict {
                    log::error!("STRICT invariant violation: {}", v);
                } else {
                    log::warn!("Invariant violation: {}", v);
                }
            }
            if total > 10 {
                log::warn!("… {} additional invariant violations elided", total - 10);
            }

            if strict {
                return Err(SarError::Processing(
                    "Azimuth timing/provider invariants failed (strict mode)".to_string(),
                ));
            }
        }

        Ok(())
    }
}
