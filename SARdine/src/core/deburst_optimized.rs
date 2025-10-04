//! Deprecated optimized TOPSAR deburst module.
//! The high-performance path has been merged into `core::deburst`.
//! This shim preserves the previous API surface while delegating to
//! the unified `TopSarDeburstProcessor` implementation.

use crate::core::deburst::{self, DeburstConfig, DeburstResult, TopSarDeburstProcessor};
use crate::types::{SarComplex, SarResult};
use ndarray::Array2;
use std::time::Instant;

/// Row mapping information retained for compatibility. The underlying
/// implementation now derives this from the canonical deburst processor.
#[derive(Clone, Default)]
pub struct RowMapping {
    pub output_to_burst: Vec<(u16, u16)>,
    pub slc_line_map: Vec<i32>,
}

/// Performance metrics structure retained for callers relying on the legacy API.
#[derive(Default, Clone)]
pub struct PerformanceMetrics {
    pub precompute_time_ms: f64,
    pub deburst_time_ms: f64,
    pub normalize_time_ms: f64,
    pub total_time_ms: f64,
    pub pixels_processed: u64,
    pub throughput_mpixels_per_sec: f64,
}

/// Compatibility wrapper around `DeburstResult`.
pub struct OptimizedDeburstResult {
    pub image: Array2<SarComplex>,
    pub row_mapping: RowMapping,
    pub performance_metrics: PerformanceMetrics,
    pub quality_score: f32,
    pub hit_count: Array2<u16>,
    pub power_ratio: f64,
    pub uncovered_pixels: usize,
    pub total_azimuth_lines: usize,
    pub total_range_samples: usize,
    pub azimuth_index_origin: usize,
}

impl OptimizedDeburstResult {
    fn from_unified(
        result: DeburstResult,
        mapping: RowMapping,
        metrics: PerformanceMetrics,
    ) -> Self {
        let DeburstResult {
            image,
            hit_count,
            power_ratio,
            uncovered_pixels,
            blend_quality_score,
            total_azimuth_lines,
            total_range_samples,
            azimuth_index_origin,
        } = result;

        Self {
            image,
            row_mapping: mapping,
            performance_metrics: metrics,
            quality_score: blend_quality_score as f32,
            hit_count,
            power_ratio,
            uncovered_pixels,
            total_azimuth_lines,
            total_range_samples,
            azimuth_index_origin,
        }
    }
}

/// Compatibility wrapper that preserves the legacy optimized processor name.
pub struct OptimizedTopsDeburstProcessor {
    inner: TopSarDeburstProcessor,
    burst_info: Vec<deburst::BurstInfo>,
    config: DeburstConfig,
}

impl OptimizedTopsDeburstProcessor {
    pub fn new(
        burst_info: &[deburst::BurstInfo],
        config: &DeburstConfig,
        satellite_velocity: f64,
    ) -> SarResult<Self> {
        let inner =
            TopSarDeburstProcessor::new(burst_info.to_vec(), config.clone(), satellite_velocity);

        Ok(Self {
            inner,
            burst_info: burst_info.to_vec(),
            config: config.clone(),
        })
    }

    pub fn process_optimized(
        &self,
        slc_data: &Array2<SarComplex>,
    ) -> SarResult<OptimizedDeburstResult> {
        let start_time = Instant::now();
        let result = self.inner.deburst_topsar_enhanced(slc_data)?;
        let total_time = start_time.elapsed();

        let performance_metrics = self.build_performance_metrics(&result, total_time);
        let row_mapping = self.build_row_mapping(result.image.dim().0);

        Ok(OptimizedDeburstResult::from_unified(
            result,
            row_mapping,
            performance_metrics,
        ))
    }

    fn build_performance_metrics(
        &self,
        result: &DeburstResult,
        total_time: std::time::Duration,
    ) -> PerformanceMetrics {
        let pixels_processed = (result.image.dim().0 * result.image.dim().1) as u64;
        let total_time_ms = total_time.as_secs_f64() * 1_000.0;
        let throughput = if total_time_ms > 0.0 {
            pixels_processed as f64 / (total_time_ms / 1_000.0) / 1_000_000.0
        } else {
            0.0
        };

        PerformanceMetrics {
            precompute_time_ms: 0.0,
            deburst_time_ms: total_time_ms,
            normalize_time_ms: 0.0,
            total_time_ms,
            pixels_processed,
            throughput_mpixels_per_sec: throughput,
        }
    }

    fn build_row_mapping(&self, output_height: usize) -> RowMapping {
        let mut output_to_burst = vec![(0u16, 0u16); output_height];
        let mut slc_line_map = vec![0i32; output_height];

        let mut cumulative_offset = 0usize;

        for (idx, burst) in self.burst_info.iter().enumerate() {
            let burst_lines = burst.lines();
            let skip_front = if idx == 0 {
                0
            } else {
                self.config.blend_lines.min(burst_lines)
            };

            for line_in_burst in skip_front..burst_lines {
                let output_line = cumulative_offset + line_in_burst;
                if output_line >= output_height {
                    break;
                }

                output_to_burst[output_line] = (idx as u16, line_in_burst as u16);
                slc_line_map[output_line] = (burst.start_line + line_in_burst) as i32;
            }

            if idx == 0 {
                cumulative_offset += burst_lines;
            } else {
                cumulative_offset += burst_lines.saturating_sub(self.config.blend_lines);
            }
        }

        RowMapping {
            output_to_burst,
            slc_line_map,
        }
    }
}

pub fn deburst_optimized(
    slc_data: &Array2<SarComplex>,
    burst_info: &[deburst::BurstInfo],
    config: &DeburstConfig,
    satellite_velocity: f64,
) -> SarResult<Array2<SarComplex>> {
    let processor = OptimizedTopsDeburstProcessor::new(burst_info, config, satellite_velocity)?;
    let result = processor.process_optimized(slc_data)?;
    Ok(result.image)
}

pub fn deburst_optimized_enhanced(
    slc_data: &Array2<SarComplex>,
    burst_info: &[deburst::BurstInfo],
    config: &DeburstConfig,
    satellite_velocity: f64,
) -> SarResult<OptimizedDeburstResult> {
    let processor = OptimizedTopsDeburstProcessor::new(burst_info, config, satellite_velocity)?;
    processor.process_optimized(slc_data)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::deburst::{BurstInfo, DeburstConfig};
    use crate::types::SarComplex;
    use ndarray::Array2;

    #[test]
    fn test_wrapper_roundtrip() {
        let burst_info = vec![BurstInfo {
            burst_id: 0,
            start_line: 0,
            end_line: 9,
            start_sample: 0,
            end_sample: 9,
            azimuth_time: String::new(),
            sensing_time: String::new(),
            first_valid_sample: vec![0; 10],
            last_valid_sample: vec![9; 10],
            byte_offset: 0,
            azimuth_fm_rate: 0.0,
            azimuth_steering_rate: 0.0,
            slant_range_time: 0.0,
            doppler_centroid: 0.0,
            azimuth_bandwidth: 0.0,
            range_sampling_rate: 0.0,
            range_pixel_spacing: 0.0,
            azimuth_pixel_spacing: 0.0,
            azimuth_time_interval: 0.0002,
            dc_polynomial: vec![0.0],
            fm_polynomial: vec![0.0],
        }];

        let config = DeburstConfig::default();
        let satellite_velocity = 7500.0;
        let slc_data = Array2::<SarComplex>::zeros((10, 10));

        let enhanced =
            deburst_optimized_enhanced(&slc_data, &burst_info, &config, satellite_velocity)
                .expect("wrapper call should succeed");

        assert_eq!(enhanced.image.dim(), (10, 10));
        assert_eq!(enhanced.performance_metrics.pixels_processed, 100);
    }
}
