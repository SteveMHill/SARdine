use crate::core::topsar_merge::{
    BlendingMethod, MergeParameters, MergedSwathData, QualityControl, TopsarMerge,
};
use crate::types::{SarError, SarImage, SarRealImage, SarResult, SubSwath};
use ndarray::Array2;
use num_complex::Complex32;
use std::collections::HashMap;

#[derive(Clone, Debug)]
struct PhaseRampConfig {
    base_phase: f32,
    dphi_dx: f32,
    dphi_dy: f32,
}

#[derive(Clone, Debug)]
struct InterpolationConfig {
    src_len: usize,
    dst_len: usize,
    scale: f32,
    offset: f32,
}

#[derive(Clone, Debug)]
struct BlendConfig {
    overlap_width: usize,
    height: usize,
}

/// Compatibility wrapper preserving the legacy optimized API surface while
/// delegating to the enhanced `TopsarMerge` implementation.
pub struct OptimizedTopsarMerge {
    width: usize,
    height: usize,
    chunk_size: usize,
    phase_config: Option<PhaseRampConfig>,
    interpolation_config: Option<InterpolationConfig>,
    blend_config: Option<BlendConfig>,
}

impl OptimizedTopsarMerge {
    /// Create a new compatibility wrapper.
    pub fn new(width: usize, height: usize, chunk_size: Option<usize>) -> SarResult<Self> {
        let chunk_size = chunk_size.unwrap_or(2048).min(height).max(256);
        log::info!(
            "♻️  Consolidated TOPSAR merge wrapper initialized: {}×{}, chunk_size={}",
            width,
            height,
            chunk_size
        );

        Ok(Self {
            width,
            height,
            chunk_size,
            phase_config: None,
            interpolation_config: None,
            blend_config: None,
        })
    }

    /// Store phase ramp configuration for diagnostics.
    pub fn precompute_phase_ramps(
        &mut self,
        height: usize,
        width: usize,
        base_phase: f32,
        dphi_dx: f32,
        dphi_dy: f32,
    ) -> SarResult<()> {
        self.phase_config = Some(PhaseRampConfig {
            base_phase,
            dphi_dx,
            dphi_dy,
        });

        log::info!(
            "🧮 Cached phase ramp config: target {}×{}, base_phase={:.4}, dphi_dx={:.5}, dphi_dy={:.5}",
            height,
            width,
            base_phase,
            dphi_dx,
            dphi_dy
        );
        Ok(())
    }

    /// Store interpolation configuration for transparency.
    pub fn precompute_interpolation_weights(
        &mut self,
        src_len: usize,
        dst_len: usize,
        scale: f32,
        offset: f32,
    ) -> SarResult<()> {
        self.interpolation_config = Some(InterpolationConfig {
            src_len,
            dst_len,
            scale,
            offset,
        });

        log::info!(
            "⚙️  Cached interpolation config: src_len={}, dst_len={}, scale={:.5}, offset={:.5}",
            src_len,
            dst_len,
            scale,
            offset
        );
        Ok(())
    }

    /// Store blend configuration for diagnostics.
    pub fn precompute_blend_weights(
        &mut self,
        overlap_width: usize,
        height: usize,
    ) -> SarResult<()> {
        self.blend_config = Some(BlendConfig {
            overlap_width,
            height,
        });

        log::info!(
            "🌈 Cached blend weights config: overlap_width={}, height={}",
            overlap_width,
            height
        );
        Ok(())
    }

    /// Legacy real-valued merge entry point.
    pub fn merge_subswaths_optimized(
        &mut self,
        subswath_data: &HashMap<String, SarRealImage>,
        subswaths: &[SubSwath],
    ) -> SarResult<SarRealImage> {
        log::info!(
            "🚀 Delegating optimized TOPSAR merge to unified processor ({}×{}, chunk_size={})",
            self.width,
            self.height,
            self.chunk_size
        );

        if let Some(cfg) = &self.phase_config {
            log::debug!(
                "  • phase ramp: base_phase={:.4}, dphi_dx={:.5}, dphi_dy={:.5}",
                cfg.base_phase,
                cfg.dphi_dx,
                cfg.dphi_dy
            );
        }
        if let Some(cfg) = &self.interpolation_config {
            log::debug!(
                "  • interpolation: src_len={}, dst_len={}, scale={:.5}, offset={:.5}",
                cfg.src_len,
                cfg.dst_len,
                cfg.scale,
                cfg.offset
            );
        }
        if let Some(cfg) = &self.blend_config {
            log::debug!(
                "  • blend weights: overlap_width={}, height={}",
                cfg.overlap_width,
                cfg.height
            );
        }

        let mut merge_params = MergeParameters::default();
        merge_params.chunk_size = self.chunk_size;
        if self.phase_config.is_some() {
            merge_params.preserve_phase = true;
            merge_params.blending_method = BlendingMethod::PhaseCoherent;
        }

        let merger = TopsarMerge::new_with_params(
            subswaths.to_vec(),
            merge_params,
            QualityControl::default(),
        )?;

        let result = merger.merge_subswaths(subswath_data, false, None)?;
        let MergedSwathData {
            merged_intensity, ..
        } = result;
        Ok(merged_intensity)
    }

    /// Legacy complex merge entry point.
    pub fn merge_complex_optimized(
        &mut self,
        subswath_data: &HashMap<String, SarImage>,
        subswaths: &[SubSwath],
    ) -> SarResult<Array2<Complex32>> {
        log::info!(
            "🔄 Delegating optimized complex TOPSAR merge to unified processor (chunk_size={})",
            self.chunk_size
        );

        // Derive intensity inputs from complex data for the unified pipeline.
        let intensity_inputs: HashMap<String, SarRealImage> = subswath_data
            .iter()
            .map(|(id, data)| (id.clone(), data.mapv(|z| z.norm_sqr())))
            .collect();

        let mut merge_params = MergeParameters::default();
        merge_params.chunk_size = self.chunk_size;
        merge_params.preserve_phase = true;
        merge_params.blending_method = BlendingMethod::PhaseCoherent;

        let merger = TopsarMerge::new_with_params(
            subswaths.to_vec(),
            merge_params,
            QualityControl::default(),
        )?;

        let result = merger.merge_subswaths(&intensity_inputs, true, Some(subswath_data))?;
        match result.merged_complex {
            Some(complex) => Ok(complex),
            None => Err(SarError::Processing(
                "Unified TOPSAR merge did not produce complex output".to_string(),
            )),
        }
    }
}

/// Performance benchmarking helper retained for compatibility with the legacy API.
pub struct PerformanceBenchmark {
    stage_times: HashMap<String, f64>,
}

impl PerformanceBenchmark {
    pub fn new() -> Self {
        Self {
            stage_times: HashMap::new(),
        }
    }

    pub fn time_stage<F, R>(&mut self, stage_name: &str, f: F) -> R
    where
        F: FnOnce() -> R,
    {
        let start = std::time::Instant::now();
        let result = f();
        let elapsed = start.elapsed().as_secs_f64();

        self.stage_times.insert(stage_name.to_string(), elapsed);
        log::info!("⏱️ {}: {:.3}s", stage_name, elapsed);

        result
    }

    pub fn report(&self) {
        let total: f64 = self.stage_times.values().sum();

        log::info!("🎯 Performance Report:");
        for (stage, time) in &self.stage_times {
            let percent = if total > 0.0 {
                (time / total) * 100.0
            } else {
                0.0
            };
            log::info!("  {}: {:.3}s ({:.1}%)", stage, time, percent);
        }
        log::info!("  Total: {:.3}s", total);
    }
}
