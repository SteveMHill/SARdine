use ndarray::parallel::prelude::*;
use ndarray::{Array2, Axis};
use num_complex::Complex;

use crate::types::SarResult;

use super::model::{sane_gain, CalibrationCoefficients, ValidSampleRanges};

/// Types of radiometric calibration
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CalibrationType {
    Sigma0, // Radar cross section per unit area
    Beta0,  // Radar brightness
    Gamma0, // Backscatter coefficient
    Dn,     // Digital numbers (uncalibrated)
}

/// Apply calibration to pre-denoised power data with optional valid sample ranges.
pub fn apply_calibration_to_denoised(
    power_data: &Array2<f32>,
    calibration: &CalibrationCoefficients,
    calibration_type: CalibrationType,
    valid_ranges: Option<&ValidSampleRanges>,
) -> SarResult<Array2<f32>> {
    let effective_ranges = valid_ranges.or(calibration.valid_sample_ranges.as_ref());

    // Apply calibration gain per pixel using the precomputed dense 2D LUT.
    // Parallelized row-wise for multi-core scalability.
    if let Some(lut) = &calibration.lut {
        log::debug!("🔍 Using precomputed LUT path, dimensions: {:?}", lut.beta_values.dim());
        let (height, width) = power_data.dim();
        let mut out = Array2::<f32>::zeros((height, width));

        out.axis_iter_mut(Axis(0))
            .into_par_iter()
            .enumerate()
            .for_each(|(row, mut row_view)| {
                let in_range = effective_ranges.map(|ranges| ranges.ranges.get(row));

                for col in 0..width {
                    // Check valid range – samples outside ValidSampleRanges are
                    // genuine zero-fill and should be treated as no-data in
                    // radiometric outputs.
                    if let Some(Some((start, end_valid))) = in_range {
                        if col < *start || col > *end_valid {
                            row_view[col] = f32::NAN;
                            continue;
                        }
                    }

                    let raw_gain = match calibration_type {
                        CalibrationType::Sigma0 => lut.sigma_values[[row, col]],
                        CalibrationType::Beta0 => lut.beta_values[[row, col]],  // Beta0 uses 1/β²
                        CalibrationType::Gamma0 => lut.gamma_values[[row, col]],
                        CalibrationType::Dn => lut.dn_values[[row, col]],
                    };

                    let gain = sane_gain(raw_gain);
                    if !gain.is_finite() {
                        // Invalid gain → mark calibrated pixel as NaN instead of
                        // silently forcing it to zero. This preserves the
                        // radiometric invariant and allows diagnostics to
                        // quantify invalid coverage.
                        row_view[col] = f32::NAN;
                        continue;
                    }

                    let power_val = power_data[[row, col]];
                    row_view[col] = power_val * gain;
                }
            });
        return Ok(out);
    }

    // Fallback: use dynamic lookup via helper; default to unity gain when lookup fails.
    log::warn!("🔍 Using FALLBACK dynamic lookup path (LUT not precomputed)");
    let (height, width) = power_data.dim();
    let mut out = Array2::<f32>::zeros((height, width));
    for row in 0..height {
        for col in 0..width {
            let in_range = effective_ranges
                .map(|ranges| {
                    ranges
                        .ranges
                        .get(row)
                        .map(|(start, end)| col >= *start && col <= *end)
                        .unwrap_or(true)
                })
                .unwrap_or(true);

            if !in_range {
                // Consistent with the LUT path: mark out-of-valid-range pixels as NaN
                // (burst zero-fill edge region) rather than setting them to 0.0.
                out[[row, col]] = f32::NAN;
                continue;
            }

            let raw_gain = match calibration.get_calibration_value(row as i32, col, calibration_type) {
                Ok(g) => g,
                Err(e) => {
                    // SCIENTIFIC FIX: Always fail on calibration lookup errors
                    // Unity gain fallback produces scientifically meaningless results
                    return Err(crate::types::SarError::Processing(format!(
                        "Calibration lookup failed at ({}, {}): {}. \
                        Calibration LUT must be complete for valid radiometric output.",
                        row, col, e
                    )));
                }
            };
            let gain = sane_gain(raw_gain);
            out[[row, col]] = if gain.is_finite() {
                power_data[[row, col]] * gain
            } else {
                f32::NAN
            };
        }
    }

    Ok(out)
}

/// Apply σ⁰ calibration to magnitude-squared data by forwarding to the legacy implementation.
pub fn apply_sigma0_calibration(
    power_data: &Array2<f32>,
    calibration: &CalibrationCoefficients,
) -> SarResult<Array2<f32>> {
    apply_calibration_to_denoised(power_data, calibration, CalibrationType::Sigma0, None)
}

/// Apply β⁰ calibration to magnitude-squared data by forwarding to the legacy implementation.
pub fn apply_beta0_calibration(
    power_data: &Array2<f32>,
    calibration: &CalibrationCoefficients,
) -> SarResult<Array2<f32>> {
    apply_calibration_to_denoised(power_data, calibration, CalibrationType::Beta0, None)
}

/// Apply γ⁰ calibration to magnitude-squared data by forwarding to the legacy implementation.
pub fn apply_gamma0_calibration(
    power_data: &Array2<f32>,
    calibration: &CalibrationCoefficients,
) -> SarResult<Array2<f32>> {
    apply_calibration_to_denoised(power_data, calibration, CalibrationType::Gamma0, None)
}

/// Return DN values (identity calibration).
pub fn apply_dn_calibration(slc_data: &Array2<Complex<f32>>) -> SarResult<Array2<Complex<f32>>> {
    Ok(slc_data.clone())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::calibration::model::{CalibrationCoefficients, CalibrationLUT};

    fn make_calibration(lut: CalibrationLUT) -> CalibrationCoefficients {
        CalibrationCoefficients {
            lut: Some(lut),
            ..CalibrationCoefficients::new()
        }
    }

    #[test]
    fn lut_selection_matches_calibration_type() {
        let mut lut = CalibrationLUT::new((2, 2));
        lut.sigma_values.fill(2.0);
        lut.beta_values.fill(3.0);
        lut.gamma_values.fill(4.0);
        lut.dn_values.fill(5.0);
        lut.is_precomputed = true;

        let calib = make_calibration(lut);
        let power = Array2::from_elem((2, 2), 10.0f32);

        let sigma = apply_calibration_to_denoised(&power, &calib, CalibrationType::Sigma0, None)
            .unwrap();
        let beta = apply_calibration_to_denoised(&power, &calib, CalibrationType::Beta0, None)
            .unwrap();
        let gamma = apply_calibration_to_denoised(&power, &calib, CalibrationType::Gamma0, None)
            .unwrap();
        let dn = apply_calibration_to_denoised(&power, &calib, CalibrationType::Dn, None).unwrap();

        assert_eq!(sigma[[0, 0]], 20.0); // 10 * 2.0
        assert_eq!(beta[[0, 0]], 30.0);  // 10 * 3.0
        assert_eq!(gamma[[0, 0]], 40.0); // 10 * 4.0
        assert_eq!(dn[[0, 0]], 50.0);    // 10 * 5.0
    }

    #[test]
    fn constant_gain_scales_power_without_introducing_zeros() {
        let mut lut = CalibrationLUT::new((4, 4));
        let gain = 1.0e-5_f32;
        lut.sigma_values.fill(gain);
        lut.beta_values.fill(gain);
        lut.gamma_values.fill(gain);
        lut.dn_values.fill(1.0);
        lut.is_precomputed = true;

        let calib = make_calibration(lut);

        let mut power = Array2::from_elem((4, 4), 10.0f32);
        // Include a few exact zeros to ensure they stay zeros after calibration
        power[(0, 0)] = 0.0;
        power[(1, 1)] = 0.0;

        let sigma = apply_calibration_to_denoised(&power, &calib, CalibrationType::Sigma0, None)
            .unwrap();

        for ((r, c), &p) in power.indexed_iter() {
            let s = sigma[(r, c)];
            if p > 0.0 {
                // Positive power must be scaled by the constant gain and must
                // not be forced to zero.
                assert!(s > 0.0);
                let expected = p * gain;
                let diff = (s - expected).abs();
                assert!(diff <= expected * 1.0e-5 + 1.0e-7);
            } else {
                // Zero input power is allowed to remain zero.
                assert_eq!(s, 0.0);
            }
        }
    }

    #[test]
    fn non_finite_gain_produces_nan_output() {
        let mut lut = CalibrationLUT::new((2, 2));
        lut.sigma_values.fill(1.0);
        lut.beta_values.fill(1.0);
        lut.gamma_values.fill(1.0);
        lut.dn_values.fill(1.0);
        // Inject a NaN gain at a single location
        lut.sigma_values[(0, 1)] = f32::NAN;
        lut.is_precomputed = true;

        let calib = make_calibration(lut);
        let power = Array2::from_elem((2, 2), 10.0f32);

        let sigma = apply_calibration_to_denoised(&power, &calib, CalibrationType::Sigma0, None)
            .unwrap();

        assert!(sigma[(0, 1)].is_nan());
        // Neighbouring pixels with finite gains should be finite and positive.
        assert!(sigma[(0, 0)].is_finite() && sigma[(0, 0)] > 0.0);
        assert!(sigma[(1, 0)].is_finite() && sigma[(1, 0)] > 0.0);
        assert!(sigma[(1, 1)].is_finite() && sigma[(1, 1)] > 0.0);
    }
}
