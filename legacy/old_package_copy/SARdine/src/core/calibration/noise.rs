use ndarray::{Array2, Axis};
use num_complex::Complex;
use rayon::prelude::*;
use std::time::Instant;

use crate::core::calibration::model::{
    AzimuthBracketCache, NoiseLUT, NoiseLutMode, MIN_VALID_POWER,
};
use crate::core::calibration::CalibrationType;
use crate::core::global_clamp_debug::ClampDebug;
use crate::types::{SarError, SarResult};

use super::model::{CalibrationCoefficients, NoiseCoefficients, NoiseVector, ValidSampleRanges};

/// Noise removal strategies. Default is resolved from environment variables:
/// - `SARDINE_SKIP_NOISE=1` → Disabled
/// - `SARDINE_NOISE_STRATEGY`/`SARDINE_NOISE_MODE` in {off,range,knots,full}
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum NoiseStrategy {
    Disabled,
    RangeOnly,
    AzimuthInterpolated,
    Full2DLut,
}

/// Noise scaling mode for ESA-compliant processing.
/// Controls whether noise removal and calibration are fused or separate.
///
/// # ESA Specification (Sentinel-1 Product Specification)
/// The correct formula for calibrated denoised backscatter is:
/// ```text
/// σ⁰ = (DN² - noiseLUT) / calibrationLUT²
///    = (DN² - noiseLUT) × gain    where gain = 1/calibrationLUT²
/// ```
///
/// The noise is subtracted in DN² domain FIRST, then calibration is applied.
///
/// # Environment Variable
/// Set `SARDINE_NOISE_SCALING` to control:
/// - `esa` or `compliant`: ESA-compliant fused denoise+calibrate (default)
/// - `simple` or `legacy`: Separate denoise then calibrate (faster)
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum NoiseScalingMode {
    /// ESA-compliant fused denoise+calibrate:
    /// σ⁰ = (DN² - noiseLUT) × gain
    /// Performs both operations in one pass for efficiency and correctness.
    #[default]
    EsaCompliant,

    /// Separate denoise then calibrate:
    /// 1. denoised = DN² - noise
    /// 2. σ⁰ = denoised × gain
    /// Mathematically equivalent but processed in two passes.
    SimpleLegacy,
}

impl NoiseScalingMode {
    /// Parse noise scaling mode from environment variable
    pub fn from_env() -> Self {
        let raw = std::env::var("SARDINE_NOISE_SCALING").unwrap_or_default();
        match raw.to_ascii_lowercase().as_str() {
            "esa" | "compliant" | "esa_compliant" | "esa-compliant" => Self::EsaCompliant,
            "simple" | "legacy" | "simple_legacy" => Self::SimpleLegacy,
            // Default: ESA-compliant for scientific accuracy
            _ => Self::EsaCompliant,
        }
    }

    /// Human-readable label
    pub fn label(&self) -> &'static str {
        match self {
            Self::EsaCompliant => "ESA-compliant ((DN²-noise)×gain)",
            Self::SimpleLegacy => "separate-legacy (denoise then calibrate)",
        }
    }
}

impl NoiseStrategy {
    pub fn from_env() -> Self {
        let skip = std::env::var("SARDINE_SKIP_NOISE").unwrap_or_default();
        if matches!(skip.as_str(), "1" | "true" | "TRUE") {
            return NoiseStrategy::Disabled;
        }

        let raw = std::env::var("SARDINE_NOISE_STRATEGY")
            .or_else(|_| std::env::var("SARDINE_NOISE_MODE"))
            .unwrap_or_default();

        match raw.to_ascii_lowercase().as_str() {
            "off" | "none" | "disabled" => NoiseStrategy::Disabled,
            "range" | "range_only" | "1d" => NoiseStrategy::RangeOnly,
            "knots" | "interp" | "azimuth" | "azimuth_interp" => NoiseStrategy::AzimuthInterpolated,
            "full" | "full2d" | "lut" => NoiseStrategy::Full2DLut,
            // Default: RangeOnly is fast and scientifically appropriate for Sentinel-1 IW mode.
            // Thermal noise is predominantly range-dependent in TOPSAR acquisitions.
            // Use SARDINE_NOISE_STRATEGY=off to disable, or =full for 2D variation.
            _ => NoiseStrategy::RangeOnly,
        }
    }

    pub fn label(&self) -> &'static str {
        match self {
            NoiseStrategy::Disabled => "disabled",
            NoiseStrategy::RangeOnly => "range-only",
            NoiseStrategy::AzimuthInterpolated => "azimuth-interpolated",
            NoiseStrategy::Full2DLut => "full-2d-lut",
        }
    }
}

/// Pre-compute noise LUT by forwarding to the legacy implementation.
pub fn precompute_noise_lut(
    coeffs: &mut NoiseCoefficients,
    image_dims: (usize, usize),
) -> SarResult<()> {
    let strategy = NoiseStrategy::from_env();
    precompute_noise_lut_with_strategy(coeffs, image_dims, strategy)
}

/// Pre-compute noise representation using a specific strategy.
pub fn precompute_noise_lut_with_strategy(
    coeffs: &mut NoiseCoefficients,
    image_dims: (usize, usize),
    strategy: NoiseStrategy,
) -> SarResult<()> {
    let t_start = Instant::now();
    let (height, width) = image_dims;

    if coeffs.vectors.is_empty() {
        return Err(SarError::Processing(
            "Noise vectors missing; cannot build noise model".to_string(),
        ));
    }

    log::info!(
        "🔇 Noise strategy: {} ({}×{})",
        strategy.label(),
        height,
        width
    );

    if matches!(strategy, NoiseStrategy::Disabled) {
        coeffs.lut = Some(NoiseLUT::disabled(image_dims));
        return Ok(());
    }

    let total = height
        .checked_mul(width)
        .ok_or_else(|| SarError::Processing("Noise LUT dimension overflow".to_string()))?;

    let lut_bytes = total.saturating_mul(std::mem::size_of::<f32>());
    let lut_mb = lut_bytes as f64 / (1024.0 * 1024.0);
    if lut_mb > 2_048.0 {
        return Err(SarError::Processing(format!(
            "Noise LUT too large: {:.1} MiB ({} px)",
            lut_mb, total
        )));
    }

    // Map output grid to SLC coordinates when a mapper is present
    let slc_lines: Vec<f64> = (0..height)
        .map(|row| {
            if let Some(mapper) = &coeffs.coordinate_mapper {
                mapper.map_coordinates(row as f64, 0.0).0
            } else {
                row as f64
            }
        })
        .collect();

    // Precompute per-vector range profiles (width samples)
    let mut vectors = coeffs.vectors.clone();
    vectors.sort_by(|a, b| {
        a.line
            .partial_cmp(&b.line)
            .unwrap_or(std::cmp::Ordering::Equal)
    });
    let precomputed_rows: Vec<Vec<f32>> = vectors
        .par_iter()
        .map(|v| interpolate_noise_row_for_vector(v, width))
        .collect();

    match strategy {
        NoiseStrategy::Full2DLut => {
            let lines: Vec<f64> = vectors.iter().map(|v| v.line).collect();
            let (slopes, intercepts) = fit_azimuth_ramps(&lines, &precomputed_rows);
            let mut noise_lut = ndarray::Array2::<f32>::zeros((height, width));

            // Precompute azimuth brackets once over sorted vectors
            let mut az_brackets: Vec<(usize, usize, f32)> = Vec::with_capacity(height);
            for &slc_line in slc_lines.iter() {
                let mut lower = 0usize;
                let mut upper = vectors.len() - 1;
                for (idx, v) in vectors.iter().enumerate() {
                    if v.line <= slc_line {
                        lower = idx;
                    }
                    if v.line >= slc_line {
                        upper = idx;
                        break;
                    }
                }
                if lower > upper {
                    lower = upper;
                }

                let line1 = vectors[lower].line as f32;
                let line2 = vectors[upper].line as f32;
                let weight = if (line2 - line1).abs() > f32::EPSILON {
                    ((slc_line as f32 - line1) / (line2 - line1)).clamp(0.0, 1.0)
                } else {
                    0.0
                };

                az_brackets.push((lower, upper, weight));
            }

            // Read SARDINE_NOISE_RAMP_WEIGHT once before the parallel section so that:
            // (a) we avoid a syscall per image row (~25 000+ calls per subswath), and
            // (b) every row uses the same value, giving deterministic results.
            // Override via SARDINE_NOISE_RAMP_WEIGHT (0.0-1.0, default 0.7).
            // - 1.0 = pure azimuth ramp model (smooth, may miss local features)
            // - 0.0 = pure interpolation (captures local variation but noisier)
            // - 0.7 = recommended default for Sentinel-1 IW
            let ramp_weight: f32 = std::env::var("SARDINE_NOISE_RAMP_WEIGHT")
                .ok()
                .and_then(|v| v.parse().ok())
                .map(|w: f32| w.clamp(0.0, 1.0))
                .unwrap_or(0.7);
            let interp_weight = 1.0 - ramp_weight;

            noise_lut
                .axis_iter_mut(Axis(0))
                .into_par_iter()
                .enumerate()
                .for_each(|(row_idx, mut row_view)| {
                    let slc_line = slc_lines[row_idx];
                    let (lower, upper, weight) = az_brackets[row_idx];
                    let lower_row = &precomputed_rows[lower];
                    let upper_row = &precomputed_rows[upper];

                    // ramp_weight and interp_weight are captured from the outer scope.
                    // SCIENTIFIC FIX: Configurable blending between azimuth ramp model and
                    // direct interpolation. Default 0.7/0.3 favors the ramp model which captures
                    // smooth azimuth variation, while interpolation handles local anomalies.
                    // Override via SARDINE_NOISE_RAMP_WEIGHT (0.0-1.0, default 0.7).
                    // - 1.0 = pure ramp model (smooth, may miss local features)
                    // - 0.0 = pure interpolation (noisy but captures local variation)
                    // - 0.7 = default (recommended for most Sentinel-1 IW data)
                    for col_idx in 0..width {
                        let n_lower = *lower_row.get(col_idx).unwrap_or(&0.0);
                        let n_upper = *upper_row.get(col_idx).unwrap_or(&0.0);
                        let interp = (1.0 - weight) * n_lower + weight * n_upper;

                        let ramp = intercepts[col_idx] + slopes[col_idx] * slc_line;
                        let blended = ramp_weight * ramp as f32 + interp_weight * interp;
                        row_view[col_idx] = blended.max(0.0);
                    }
                });

            coeffs.lut = Some(NoiseLUT {
                noise_values: noise_lut,
                range_profile: None,
                precomputed_rows: None,
                azimuth_brackets: None,
                azimuth_axis: (0..height).map(|v| v as f64).collect(),
                range_axis: (0..width).map(|v| v as f64).collect(),
                height,
                width,
                mode: NoiseLutMode::Full2D,
                is_precomputed: true,
            });

            log::info!(
                "🔇 Noise LUT precomputed (full 2D): {}x{} ({:.1} MiB) in {:.2?}",
                height,
                width,
                lut_mb,
                t_start.elapsed()
            );
        }
        NoiseStrategy::RangeOnly => {
            // Average per-range noise across azimuth vectors (fast, low memory)
            let mut range_profile = vec![0.0f32; width];
            for row in &precomputed_rows {
                for (col, val) in row.iter().enumerate().take(width) {
                    range_profile[col] += *val;
                }
            }
            let denom = precomputed_rows.len().max(1) as f32;
            for v in range_profile.iter_mut() {
                *v = (*v / denom).max(0.0);
            }

            coeffs.lut = Some(NoiseLUT {
                noise_values: ndarray::Array2::zeros((0, 0)),
                range_profile: Some(range_profile),
                precomputed_rows: None,
                azimuth_brackets: None,
                azimuth_axis: vec![],
                range_axis: (0..width).map(|v| v as f64).collect(),
                height,
                width,
                mode: NoiseLutMode::RangeOnly,
                is_precomputed: true,
            });

            log::info!(
                "🔇 Noise LUT (range-only) ready: ~{:.1} MiB in {:.2?}",
                (width * std::mem::size_of::<f32>()) as f64 / (1024.0 * 1024.0),
                t_start.elapsed()
            );
        }
        NoiseStrategy::AzimuthInterpolated => {
            // Cache azimuth brackets for each output row to blend per-vector rows at runtime
            let mut az_brackets: Vec<AzimuthBracketCache> = Vec::with_capacity(height);
            for &slc_line in slc_lines.iter() {
                let mut lower = 0usize;
                let mut upper = vectors.len() - 1;
                for (idx, v) in vectors.iter().enumerate() {
                    if v.line <= slc_line {
                        lower = idx;
                    }
                    if v.line >= slc_line {
                        upper = idx;
                        break;
                    }
                }
                if lower > upper {
                    lower = upper;
                }

                let line1 = vectors[lower].line as f32;
                let line2 = vectors[upper].line as f32;
                let weight = if (line2 - line1).abs() > f32::EPSILON {
                    ((slc_line as f32 - line1) / (line2 - line1)).clamp(0.0, 1.0)
                } else {
                    0.0
                };

                az_brackets.push(AzimuthBracketCache {
                    lower,
                    upper,
                    weight: weight.dbg_clamp(0.0, 1.0, "noise_az_weight"),
                });
            }

            coeffs.lut = Some(NoiseLUT {
                noise_values: ndarray::Array2::zeros((0, 0)),
                range_profile: None,
                precomputed_rows: Some(precomputed_rows.clone()),
                azimuth_brackets: Some(az_brackets),
                azimuth_axis: (0..height).map(|v| v as f64).collect(),
                range_axis: (0..width).map(|v| v as f64).collect(),
                height,
                width,
                mode: NoiseLutMode::AzimuthInterpolated,
                is_precomputed: true,
            });

            let approx_mb = (precomputed_rows.len() * width * std::mem::size_of::<f32>()) as f64
                / (1024.0 * 1024.0);
            log::info!(
                "🔇 Noise LUT (azimuth-interpolated) cached rows: {} vectors × {} cols (~{:.1} MiB) in {:.2?}",
                precomputed_rows.len(),
                width,
                approx_mb,
                t_start.elapsed()
            );
        }
        NoiseStrategy::Disabled => {
            coeffs.lut = Some(NoiseLUT::disabled(image_dims));
        }
    }

    Ok(())
}

/// Apply thermal noise removal (power domain) in-place.
pub fn apply_thermal_noise_removal_inplace(
    power_data: &mut Array2<f32>,
    noise_coefficients: &NoiseCoefficients,
) -> SarResult<()> {
    let (h, w) = power_data.dim();

    let lut = match &noise_coefficients.lut {
        Some(l) if l.is_precomputed => l,
        Some(_) => {
            log::warn!("⚠️  Noise LUT not precomputed; skipping noise removal");
            return Ok(());
        }
        None => {
            log::warn!("⚠️  Noise coefficients missing LUT; skipping noise removal");
            return Ok(());
        }
    };

    if lut.height != h || lut.width != w {
        // SCIENTIFIC FIX: Always fail on dimension mismatch - partial noise removal is wrong
        return Err(SarError::Processing(format!(
            "Noise LUT dimension mismatch: LUT is {}x{}, data is {}x{}. \
            Noise removal requires exact dimension match. \
            Check that noise LUT was computed for the correct image dimensions.",
            lut.height, lut.width, h, w
        )));
    }

    if matches!(lut.mode, NoiseLutMode::Disabled) {
        log::info!("🔇 Noise removal disabled (noise strategy=off)");
        return Ok(());
    }

    // Pre-flight statistics to detect pathological LUT/data ratios.
    let total = (h * w).max(1) as f64;
    let mut denoised = 0usize;

    match lut.mode {
        NoiseLutMode::Full2D => {
            let stride = (h * w / 4096).max(1);
            let mut sample_power = 0.0f64;
            let mut sample_noise = 0.0f64;
            let mut would_zero = 0usize;

            for (idx, (&p, &n)) in power_data.iter().zip(lut.noise_values.iter()).enumerate() {
                if idx % stride == 0 {
                    if p.is_finite() {
                        sample_power += p as f64;
                    }
                    if n.is_finite() {
                        sample_noise += n as f64;
                    }
                }
                if !p.is_finite() || p - n <= MIN_VALID_POWER {
                    would_zero += 1;
                }
            }

            let zero_ratio = would_zero as f64 / total;
            let sample_norm = (total / stride as f64).max(1.0);
            let mean_power_sample = sample_power / sample_norm;
            let mean_noise_sample = sample_noise / sample_norm;

            // CONFIGURABLE THRESHOLDS: Allow tuning via environment variables
            // SARDINE_NOISE_ZERO_THRESHOLD: Maximum fraction of pixels that would be zeroed (default: 0.95)
            // SARDINE_NOISE_RATIO_THRESHOLD: Maximum noise/power ratio before skipping (default: 0.9)
            let zero_threshold: f64 = std::env::var("SARDINE_NOISE_ZERO_THRESHOLD")
                .ok()
                .and_then(|v| v.parse().ok())
                .map(|t: f64| t.clamp(0.5, 0.99))
                .unwrap_or(0.95);

            let ratio_threshold: f64 = std::env::var("SARDINE_NOISE_RATIO_THRESHOLD")
                .ok()
                .and_then(|v| v.parse().ok())
                .map(|t: f64| t.clamp(0.5, 0.99))
                .unwrap_or(0.9);

            if zero_ratio > zero_threshold
                || (mean_power_sample > 0.0
                    && mean_noise_sample > ratio_threshold * mean_power_sample)
            {
                log::warn!(
                    "⚠️  Skipping thermal noise removal: unsafe ratio (would zero {:.1}% pixels, threshold {:.0}%) \
                    or noise mean {:.3e} too high vs power {:.3e} (ratio threshold {:.0}%)",
                    zero_ratio * 100.0,
                    zero_threshold * 100.0,
                    mean_noise_sample,
                    mean_power_sample,
                    ratio_threshold * 100.0
                );

                // FALL-2/FALL-3: Emit quality flag instead of silent skip
                use crate::core::quality_flags::{global_quality_flags, QualityFlag};
                let noise_power_ratio = if mean_power_sample > 0.0 {
                    Some(mean_noise_sample / mean_power_sample)
                } else {
                    None
                };
                global_quality_flags().add(QualityFlag::NoiseRemovalSkipped {
                    reason: if zero_ratio > zero_threshold {
                        format!(
                            "Would zero {:.1}% of pixels (threshold {:.0}%)",
                            zero_ratio * 100.0,
                            zero_threshold * 100.0
                        )
                    } else {
                        format!(
                            "Noise/power ratio {:.3} exceeds threshold {:.2}",
                            noise_power_ratio.unwrap_or(0.0),
                            ratio_threshold
                        )
                    },
                    zero_ratio: Some(zero_ratio),
                    noise_power_ratio,
                    zero_threshold,
                    ratio_threshold,
                });

                return Ok(());
            }

            ndarray::Zip::from(power_data.view_mut())
                .and(lut.noise_values.view())
                .for_each(|p, &n| {
                    let val = *p - n;
                    if val.is_finite() && val > MIN_VALID_POWER {
                        *p = val;
                    } else {
                        // Scientific requirement: treat noise-dominated pixels as invalid
                        // and propagate them as NaN rather than silent zeros.
                        *p = f32::NAN;
                        denoised += 1;
                    }
                });

            let denoise_ratio = denoised as f64 / total;
            if denoise_ratio > 0.9 {
                log::warn!(
                    "⚠️  Thermal noise removal zeroed {:.1}% of pixels — check noise LUT units/coverage",
                    denoise_ratio * 100.0
                );
            } else {
                log::info!(
                    "🔇 Thermal noise removal zeroed {:.1}% of pixels (mean noise {:.3e} vs power {:.3e})",
                    denoise_ratio * 100.0,
                    mean_noise_sample,
                    mean_power_sample
                );
            }
        }
        NoiseLutMode::RangeOnly => {
            let profile = match &lut.range_profile {
                Some(p) if p.len() == w => p,
                _ => {
                    log::warn!("⚠️  Range-only noise profile missing/mismatched; skipping");
                    return Ok(());
                }
            };

            for mut row in power_data.axis_iter_mut(Axis(0)) {
                for col in 0..w {
                    let val = row[col] - profile[col];
                    if val.is_finite() && val > MIN_VALID_POWER {
                        row[col] = val;
                    } else {
                        row[col] = f32::NAN;
                        denoised += 1;
                    }
                }
            }

            let denoise_ratio = denoised as f64 / total;
            log::info!(
                "🔇 Thermal noise removal (range-only) zeroed {:.1}% of pixels",
                denoise_ratio * 100.0
            );
        }
        NoiseLutMode::AzimuthInterpolated => {
            let rows = match &lut.precomputed_rows {
                Some(r) if !r.is_empty() => r,
                _ => {
                    log::warn!("⚠️  Noise rows missing for azimuth-interpolated mode; skipping");
                    return Ok(());
                }
            };
            let brackets = match &lut.azimuth_brackets {
                Some(b) if b.len() == h => b,
                _ => {
                    log::warn!("⚠️  Noise azimuth brackets missing/mismatched; skipping");
                    return Ok(());
                }
            };

            // OPTIMIZATION: Parallelize row processing with rayon (Phase 3)
            // Each row's noise subtraction is independent
            use std::sync::atomic::{AtomicUsize, Ordering};
            let denoised_atomic = AtomicUsize::new(0);

            power_data
                .axis_iter_mut(Axis(0))
                .into_par_iter()
                .enumerate()
                .for_each(|(row_idx, mut row)| {
                    let bracket = &brackets[row_idx];
                    let lower_row = &rows[bracket.lower];
                    let upper_row = &rows[bracket.upper];
                    let weight = bracket.weight;
                    let mut local_denoised = 0usize;

                    for col in 0..w {
                        let n_lower = *lower_row.get(col).unwrap_or(&0.0);
                        let n_upper = *upper_row.get(col).unwrap_or(&0.0);
                        let noise = (1.0 - weight) * n_lower + weight * n_upper;

                        let val = row[col] - noise;
                        if val.is_finite() && val > MIN_VALID_POWER {
                            row[col] = val;
                        } else {
                            row[col] = f32::NAN;
                            local_denoised += 1;
                        }
                    }
                    denoised_atomic.fetch_add(local_denoised, Ordering::Relaxed);
                });

            denoised = denoised_atomic.load(Ordering::Relaxed);
            let denoise_ratio = denoised as f64 / total;
            log::info!(
                "🔇 Thermal noise removal (azimuth-interpolated) zeroed {:.1}% of pixels",
                denoise_ratio * 100.0
            );
        }
        NoiseLutMode::Disabled => {
            log::info!("🔇 Noise removal disabled; no-op");
        }
    }

    Ok(())
}

/// Apply thermal noise removal (power domain) returning a new array.
pub fn apply_thermal_noise_removal(
    power_data: &Array2<f32>,
    noise_coefficients: &NoiseCoefficients,
) -> SarResult<Array2<f32>> {
    let mut out = power_data.clone();
    apply_thermal_noise_removal_inplace(&mut out, noise_coefficients)?;
    Ok(out)
}

/// Apply ESA-compliant fused denoise+calibrate in a single pass.
///
/// # ESA Specification (Sentinel-1 Product Specification)
/// The correct formula for calibrated denoised backscatter is:
/// ```text
/// σ⁰ = (DN² - noiseLUT) / calibrationLUT²
///    = (DN² - noiseLUT) × gain    where gain = 1/calibrationLUT²
/// ```
///
/// CRITICAL: The noise is subtracted in DN² domain FIRST, then calibration
/// is applied to the denoised result. This is mathematically equivalent to:
/// ```text
/// σ⁰ = DN²/LUT² - noiseLUT/LUT²
/// ```
///
/// # Fused Implementation
/// This function performs both operations in one pass for efficiency:
/// ```text
/// σ⁰[i] = (power[i] - noise[i]) × gain[i]
/// ```
///
/// # Arguments
/// * `power_data` - Input power data (DN²)
/// * `noise_coefficients` - Thermal noise coefficients with pre-computed LUT
/// * `calibration_coefficients` - Calibration coefficients with pre-computed LUT
/// * `calibration_type` - Which calibration LUT to use
///
/// # Returns
/// * ESA-compliant calibrated denoised backscatter (σ⁰, β⁰, or γ⁰)
pub fn apply_esa_compliant_fused_denoise_calibrate(
    power_data: &Array2<f32>,
    noise_coefficients: &NoiseCoefficients,
    calibration_coefficients: &CalibrationCoefficients,
    calibration_type: CalibrationType,
) -> SarResult<Array2<f32>> {
    let (h, w) = power_data.dim();

    // Get noise LUT
    let noise_lut = match &noise_coefficients.lut {
        Some(l) if l.is_precomputed => l,
        _ => {
            log::warn!("⚠️  Noise LUT not precomputed; falling back to simple subtraction");
            let mut denoised = power_data.clone();
            apply_thermal_noise_removal_inplace(&mut denoised, noise_coefficients)?;
            return crate::core::calibration::apply_calibration_to_denoised(
                &denoised,
                calibration_coefficients,
                calibration_type,
                None,
            );
        }
    };

    // Get calibration LUT
    let cal_lut = match &calibration_coefficients.lut {
        Some(l) if l.is_precomputed => l,
        _ => {
            log::warn!("⚠️  Calibration LUT not precomputed; falling back to simple subtraction");
            let mut denoised = power_data.clone();
            apply_thermal_noise_removal_inplace(&mut denoised, noise_coefficients)?;
            return crate::core::calibration::apply_calibration_to_denoised(
                &denoised,
                calibration_coefficients,
                calibration_type,
                None,
            );
        }
    };

    // Verify dimensions
    if noise_lut.height != h || noise_lut.width != w {
        return Err(crate::types::SarError::Processing(format!(
            "Noise LUT dimension mismatch: LUT is {}x{}, data is {}x{}",
            noise_lut.height, noise_lut.width, h, w
        )));
    }

    let cal_dims = cal_lut.sigma_values.dim();
    if cal_dims.0 != h || cal_dims.1 != w {
        return Err(crate::types::SarError::Processing(format!(
            "Calibration LUT dimension mismatch: LUT is {}x{}, data is {}x{}",
            cal_dims.0, cal_dims.1, h, w
        )));
    }

    if matches!(noise_lut.mode, NoiseLutMode::Disabled) {
        log::info!("🔇 Noise removal disabled; applying calibration only");
        return crate::core::calibration::apply_calibration_to_denoised(
            power_data,
            calibration_coefficients,
            calibration_type,
            None,
        );
    }

    log::info!(
        "🔇 Applying ESA-compliant fused denoise+calibrate ({:?}): σ⁰ = (DN² - noise) × gain",
        calibration_type
    );

    let mut result = Array2::<f32>::zeros((h, w));
    let total = (h * w).max(1) as f64;
    let mut zeroed_count = 0usize;

    // Select calibration gain array based on type (gain = 1/LUT²)
    let cal_gains = match calibration_type {
        CalibrationType::Sigma0 => &cal_lut.sigma_values,
        CalibrationType::Beta0 => &cal_lut.beta_values,
        CalibrationType::Gamma0 => &cal_lut.gamma_values,
        CalibrationType::Dn => &cal_lut.dn_values,
    };

    match noise_lut.mode {
        NoiseLutMode::Full2D => {
            for row in 0..h {
                for col in 0..w {
                    let power = power_data[[row, col]];
                    let noise = noise_lut.noise_values[[row, col]];
                    let cal_gain = cal_gains[[row, col]];

                    // ESA formula: σ⁰ = (DN² - noiseLUT) / calibrationLUT²
                    //            = (DN² - noiseLUT) × gain
                    // Denoise in DN² domain first, then calibrate
                    let denoised_dn2 = power - noise;

                    if denoised_dn2.is_finite() && denoised_dn2 > MIN_VALID_POWER {
                        // Apply calibration to denoised DN²
                        result[[row, col]] = denoised_dn2 * cal_gain;
                    } else {
                        result[[row, col]] = 0.0;
                        zeroed_count += 1;
                    }
                }
            }
        }
        NoiseLutMode::RangeOnly => {
            let profile = match &noise_lut.range_profile {
                Some(p) if p.len() == w => p,
                _ => {
                    log::warn!("⚠️  Range-only noise profile missing/mismatched");
                    let mut denoised = power_data.clone();
                    apply_thermal_noise_removal_inplace(&mut denoised, noise_coefficients)?;
                    return crate::core::calibration::apply_calibration_to_denoised(
                        &denoised,
                        calibration_coefficients,
                        calibration_type,
                        None,
                    );
                }
            };

            for row in 0..h {
                for col in 0..w {
                    let power = power_data[[row, col]];
                    let noise = profile[col];
                    let cal_gain = cal_gains[[row, col]];

                    let denoised_dn2 = power - noise;

                    if denoised_dn2.is_finite() && denoised_dn2 > MIN_VALID_POWER {
                        result[[row, col]] = denoised_dn2 * cal_gain;
                    } else {
                        result[[row, col]] = 0.0;
                        zeroed_count += 1;
                    }
                }
            }
        }
        NoiseLutMode::AzimuthInterpolated => {
            let rows = match &noise_lut.precomputed_rows {
                Some(r) if !r.is_empty() => r,
                _ => {
                    log::warn!("⚠️  Noise rows missing for azimuth-interpolated mode");
                    let mut denoised = power_data.clone();
                    apply_thermal_noise_removal_inplace(&mut denoised, noise_coefficients)?;
                    return crate::core::calibration::apply_calibration_to_denoised(
                        &denoised,
                        calibration_coefficients,
                        calibration_type,
                        None,
                    );
                }
            };
            let brackets = match &noise_lut.azimuth_brackets {
                Some(b) if b.len() == h => b,
                _ => {
                    log::warn!("⚠️  Noise azimuth brackets missing/mismatched");
                    let mut denoised = power_data.clone();
                    apply_thermal_noise_removal_inplace(&mut denoised, noise_coefficients)?;
                    return crate::core::calibration::apply_calibration_to_denoised(
                        &denoised,
                        calibration_coefficients,
                        calibration_type,
                        None,
                    );
                }
            };

            for row in 0..h {
                let bracket = &brackets[row];
                let lower_row = &rows[bracket.lower];
                let upper_row = &rows[bracket.upper];
                let weight = bracket.weight;

                for col in 0..w {
                    let power = power_data[[row, col]];
                    let n_lower = *lower_row.get(col).unwrap_or(&0.0);
                    let n_upper = *upper_row.get(col).unwrap_or(&0.0);
                    let noise = (1.0 - weight) * n_lower + weight * n_upper;
                    let cal_gain = cal_gains[[row, col]];

                    let denoised_dn2 = power - noise;

                    if denoised_dn2.is_finite() && denoised_dn2 > MIN_VALID_POWER {
                        result[[row, col]] = denoised_dn2 * cal_gain;
                    } else {
                        result[[row, col]] = 0.0;
                        zeroed_count += 1;
                    }
                }
            }
        }
        NoiseLutMode::Disabled => {
            return crate::core::calibration::apply_calibration_to_denoised(
                power_data,
                calibration_coefficients,
                calibration_type,
                None,
            );
        }
    }

    let zero_ratio = zeroed_count as f64 / total;
    log::info!(
        "🔇 ESA-compliant fused denoise+calibrate: zeroed {:.1}% of pixels",
        zero_ratio * 100.0
    );

    Ok(result)
}

/// Apply fused thermal noise removal + calibration (power domain).
///
/// Uses ESA-compliant fused denoise+calibrate when `SARDINE_NOISE_SCALING=esa` (default).
/// Falls back to separate denoise then calibrate with `SARDINE_NOISE_SCALING=simple`.
///
/// # ESA Formula
/// σ⁰ = (DN² - noiseLUT) / calibrationLUT²
///    = (DN² - noiseLUT) × gain
///
/// The noise is subtracted in DN² domain, then calibration is applied once.
pub fn apply_fused_noise_calibration(
    power_data: &Array2<f32>,
    noise_coefficients: &NoiseCoefficients,
    calibration_coefficients: &CalibrationCoefficients,
    calibration_type: CalibrationType,
    valid_ranges: Option<&ValidSampleRanges>,
) -> SarResult<Array2<f32>> {
    let scaling_mode = NoiseScalingMode::from_env();
    log::info!("🔇 Noise scaling mode: {}", scaling_mode.label());

    match scaling_mode {
        NoiseScalingMode::EsaCompliant => {
            // ESA-compliant: fused denoise+calibrate in one pass
            // σ⁰ = (DN² - noise) × gain
            apply_esa_compliant_fused_denoise_calibrate(
                power_data,
                noise_coefficients,
                calibration_coefficients,
                calibration_type,
            )
        }
        NoiseScalingMode::SimpleLegacy => {
            // Legacy: separate denoise then calibrate
            let mut denoised = power_data.clone();
            apply_thermal_noise_removal_inplace(&mut denoised, noise_coefficients)?;
            crate::core::calibration::apply_calibration_to_denoised(
                &denoised,
                calibration_coefficients,
                calibration_type,
                valid_ranges,
            )
        }
    }
}

/// Apply ultra-fused SLC→backscatter processing via legacy implementation.
pub fn apply_fused_slc_calibration(
    slc_data: &Array2<Complex<f32>>,
    noise_coefficients: &NoiseCoefficients,
    calibration_coefficients: &CalibrationCoefficients,
    calibration_type: CalibrationType,
    valid_ranges: Option<&ValidSampleRanges>,
) -> SarResult<Array2<f32>> {
    let (h, w) = slc_data.dim();
    let mut power = Array2::<f32>::zeros((h, w));
    ndarray::Zip::from(power.view_mut())
        .and(slc_data.view())
        .for_each(|p, &c| *p = c.norm_sqr());
    apply_fused_noise_calibration(
        &power,
        noise_coefficients,
        calibration_coefficients,
        calibration_type,
        valid_ranges,
    )
}

/// Interpolate a noise row for a given noise vector (helper used in validation paths).
pub fn interpolate_noise_row_for_vector(vector: &NoiseVector, width: usize) -> Vec<f32> {
    if vector.range_pixels.is_empty() || vector.noise_range_lut.is_empty() {
        return vec![0.0; width];
    }
    let xs = &vector.range_pixels;
    let ys = &vector.noise_range_lut;

    let mut out = vec![0.0f32; width];

    for x in 0..width {
        let pos = x as f64;

        if pos <= xs[0] {
            out[x] = ys[0];
            continue;
        }
        // SAFETY FIX: Use safer pattern instead of unwrap() to prevent panic on edge cases
        let last_x = match xs.last() {
            Some(&val) => val,
            None => continue, // Already handled by is_empty check above, but be defensive
        };
        let last_y = match ys.last() {
            Some(&val) => val,
            None => continue,
        };
        if pos >= last_x {
            out[x] = last_y;
            continue;
        }

        let mut lo = 0usize;
        let mut hi = xs.len() - 1;
        while hi - lo > 1 {
            let mid = (lo + hi) / 2;
            if xs[mid] <= pos {
                lo = mid;
            } else {
                hi = mid;
            }
        }

        let x1 = xs[lo];
        let x2 = xs[hi];
        let y1 = ys[lo];
        let y2 = ys[hi];

        if (x2 - x1).abs() < f64::EPSILON {
            out[x] = y1;
        } else {
            let t = ((pos - x1) / (x2 - x1)).clamp(0.0, 1.0) as f32;
            out[x] = y1 + t * (y2 - y1);
        }
    }

    out
}

fn fit_azimuth_ramps(lines: &[f64], rows: &[Vec<f32>]) -> (Vec<f64>, Vec<f64>) {
    let width = rows.first().map(|r| r.len()).unwrap_or(0);
    let mut slopes = vec![0.0f64; width];
    let mut intercepts = vec![0.0f64; width];
    if lines.is_empty() || width == 0 {
        return (slopes, intercepts);
    }

    let n = lines.len() as f64;
    let mean_x = lines.iter().sum::<f64>() / n;

    for col in 0..width {
        let mut sum_y = 0.0;
        for row in rows {
            sum_y += row.get(col).copied().unwrap_or(0.0) as f64;
        }
        let mean_y = sum_y / n;

        let mut num = 0.0;
        let mut den = 0.0;
        for (x, row) in lines.iter().zip(rows.iter()) {
            let y = row.get(col).copied().unwrap_or(0.0) as f64;
            let dx = x - mean_x;
            num += dx * (y - mean_y);
            den += dx * dx;
        }
        let slope = if den.abs() > 1e-9 { num / den } else { 0.0 };
        let intercept = mean_y - slope * mean_x;

        slopes[col] = slope;
        intercepts[col] = intercept;
    }

    (slopes, intercepts)
}
