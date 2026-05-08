#![allow(dead_code)]
//! DC/FM provider construction and validation.
//!
//! Builds Doppler Centroid (DC) and FM rate providers for each subswath
//! using burst timing and polynomial metadata.

use crate::core::geometry::dc_fm_provider::OutOfRangePolicy;
use crate::core::geometry::type_safe_units::Seconds;
use crate::core::{DcFmRateProvider, DcPolynomial, FmPolynomial, PolynomialDcFmProvider};
use crate::types::{SarError, SarResult, SubSwath};
use std::collections::HashMap;
use std::sync::Arc;

use super::super::types::AzimuthTimingModel;
use super::super::assert_subswath_geometry;

/// Default satellite velocity when not available from metadata
const DEFAULT_SATELLITE_VELOCITY_MPS: f64 = 7_500.0;

/// Sentinel-1 antenna length in azimuth (meters)
const SENTINEL1_ANTENNA_LENGTH_AZIMUTH_M: f64 = 12.3;

/// Validate DC prerequisites before building providers.
///
/// CRITICAL: Prevents silent alignment failures from zero-DC fallback.
pub fn validate_dc_prerequisites(subswaths: &[SubSwath]) -> SarResult<()> {
    for sw in subswaths {
        #[cfg(debug_assertions)]
        assert_subswath_geometry(sw);

        match &sw.dc_polynomial {
            None => {
                return Err(SarError::Processing(format!(
                    "❌ CRITICAL: Subswath {} lacks DC polynomial!\n\
                     This indicates deburst did NOT use DC-aware processing.\n\
                     Inter-subswath alignment will be INCORRECT.\n\
                     Solution: Re-run deburst with DC extraction enabled.",
                    sw.id
                )));
            }
            Some(coeffs) if coeffs.is_empty() => {
                return Err(SarError::Processing(format!(
                    "❌ CRITICAL: Subswath {} has empty DC polynomial!\n\
                     This indicates DC extraction failed during deburst.\n\
                     Inter-subswath alignment will be INCORRECT.",
                    sw.id
                )));
            }
            Some(coeffs) => {
                log::info!(
                    "✅ Subswath {} has DC polynomial: {} coefficients",
                    sw.id,
                    coeffs.len()
                );
            }
        }

        if sw.dc_polynomial_t0.is_none() {
            return Err(SarError::Processing(format!(
                "❌ CRITICAL: Subswath {} missing dc_polynomial_t0 reference time; DC evaluation would silently fall back to zero.",
                sw.id
            )));
        }

        if sw.azimuth_time_interval.is_none() {
            log::warn!(
                "⚠️ Subswath {} lacks azimuth time interval - will use PRF fallback",
                sw.id
            );
        }
    }

    log::info!("✅ DC prerequisites validated - proceeding with DC-aware merge");
    Ok(())
}

/// Calculate azimuth FM rate for TOPSAR steering correction.
///
/// Requires metadata PRF - no fallback estimation allowed.
///
/// Implements: K_a = 2 * PRF² / (V_s * L_antenna)
pub fn calculate_azimuth_fm_rate(subswath: &SubSwath, satellite_velocity: f64) -> SarResult<f64> {
    #[cfg(debug_assertions)]
    assert_subswath_geometry(subswath);

    // SCIENTIFIC FIX: Always require annotation-derived PRF
    let prf = subswath.prf_hz.ok_or_else(|| {
        SarError::Metadata(format!(
            "PRF missing for subswath {} - required for FM rate calculation. \
            Check product integrity.",
            subswath.id
        ))
    })?;

    let antenna_length_azimuth = SENTINEL1_ANTENNA_LENGTH_AZIMUTH_M;

    // K_a = 2 * PRF² / (V_s * L_antenna)
    let azimuth_fm_rate = 2.0 * prf.powi(2) / (satellite_velocity * antenna_length_azimuth);

    log::debug!(
        "🔬 Azimuth FM rate calculation: PRF={:.1} Hz, V_s={:.1} m/s, L_ant={:.1} m → K_a={:.2e} Hz/s",
        prf, satellite_velocity, antenna_length_azimuth, azimuth_fm_rate
    );

    Ok(azimuth_fm_rate)
}

/// Build a DC/FM provider for a single subswath.
pub fn build_dc_fm_provider_for_swath(
    sw: &SubSwath,
    time_range: (f64, f64),
    azimuth_fm_rate: f64,
) -> SarResult<Arc<dyn DcFmRateProvider>> {
    #[cfg(debug_assertions)]
    assert_subswath_geometry(sw);

    let coefficients = sw.dc_polynomial.clone().ok_or_else(|| {
        SarError::Processing(format!(
            "Subswath {} missing DC polynomial when building provider",
            sw.id
        ))
    })?;

    if coefficients.is_empty() {
        return Err(SarError::Processing(format!(
            "Subswath {} has empty DC polynomial when building provider",
            sw.id
        )));
    }

    let t0 = sw.dc_polynomial_t0.ok_or_else(|| {
        SarError::Processing(format!(
            "Subswath {} missing dc_polynomial_t0 when building provider",
            sw.id
        ))
    })?;

    let reference_time = Seconds::new(t0);

    // SCIENTIFIC FIX (Jan 2026): Require absolute azimuth times here.
    // Silent promotion of relative time ranges to epoch using t0 made it
    // difficult to diagnose inconsistent timing domains. The azimuth_timing
    // model is now responsible for providing absolute seconds since epoch.
    let (start_raw, end_raw) = time_range;
    if start_raw.abs() < 1.0e6 {
        return Err(SarError::Metadata(format!(
            "DC/FM time range for subswath {} appears non-epoch [{:.6}, {:.6}]s; expected absolute seconds since annotation epoch.",
            sw.id, start_raw, end_raw
        )));
    }

    let (start_abs, end_abs) = (start_raw, end_raw);

    // Allow tolerance around burst timing to absorb rounding drift
    const AZ_TIME_EPS: f64 = 2.0e-3;
    let start = (start_abs - AZ_TIME_EPS).max(0.0);
    let end = (end_abs + AZ_TIME_EPS).max(start);
    let time_range = (Seconds::new(start), Seconds::new(end));

    let dc_poly = DcPolynomial::new(coefficients, reference_time, time_range);
    let fm_poly = FmPolynomial::new(vec![azimuth_fm_rate], reference_time, time_range);

    let provider = PolynomialDcFmProvider::with_policy(
        vec![dc_poly],
        vec![fm_poly],
        vec![time_range],
        OutOfRangePolicy::Clamp,
    )?;

    Ok(Arc::new(provider))
}

/// Build DC/FM providers for all subswaths.
pub fn build_dc_fm_provider_map(
    subswaths: &[SubSwath],
    azimuth_timing: &AzimuthTimingModel,
) -> SarResult<HashMap<String, Arc<dyn DcFmRateProvider>>> {
    let mut providers = HashMap::new();

    for sw in subswaths {
        #[cfg(debug_assertions)]
        assert_subswath_geometry(sw);

        let swath_bursts: Vec<_> = azimuth_timing
            .burst_timing
            .iter()
            .filter(|b| b.subswath_id == sw.id)
            .collect();

        if swath_bursts.is_empty() {
            return Err(SarError::Processing(format!(
                "Missing burst timing for subswath {}",
                sw.id
            )));
        }

        let start = swath_bursts
            .iter()
            .map(|b| b.azimuth_time_start)
            .fold(f64::INFINITY, |acc, v| acc.min(v));
        let end = swath_bursts
            .iter()
            .map(|b| b.azimuth_time_end)
            .fold(f64::NEG_INFINITY, |acc, v| acc.max(v));

        if !start.is_finite() || !end.is_finite() || end < start {
            return Err(SarError::Processing(format!(
                "Invalid burst timing range for subswath {} (start={:.6}, end={:.6})",
                sw.id, start, end
            )));
        }

        let fm_rate = calculate_azimuth_fm_rate(sw, DEFAULT_SATELLITE_VELOCITY_MPS)?;

        let provider = build_dc_fm_provider_for_swath(sw, (start, end), fm_rate)?;

        providers.insert(sw.id.clone(), provider);
    }

    Ok(providers)
}

/// Validate provider coverage spans the merged azimuth timeline.
pub fn validate_dc_fm_providers(
    providers: &HashMap<String, Arc<dyn DcFmRateProvider>>,
    azimuth_timing: &AzimuthTimingModel,
    subswaths: &[SubSwath],
) -> SarResult<()> {
    const AZ_TIME_EPS: f64 = 2.0e-3;

    for sw in subswaths {
        #[cfg(debug_assertions)]
        assert_subswath_geometry(sw);

        let provider = providers.get(&sw.id).ok_or_else(|| {
            SarError::Processing(format!(
                "Missing DC/FM provider for subswath {} during coverage validation",
                sw.id
            ))
        })?;

        let (start, end) = provider.get_time_range();
        for burst in azimuth_timing
            .burst_timing
            .iter()
            .filter(|b| b.subswath_id == sw.id)
        {
            let line_count = burst
                .last_line_merged
                .saturating_sub(burst.first_line_merged)
                .saturating_add(1);
            let expected_end = burst.azimuth_time_start
                + (line_count.saturating_sub(1)) as f64 * burst.azimuth_time_interval;

            if start.value() > burst.azimuth_time_start + AZ_TIME_EPS
                || end.value() + AZ_TIME_EPS < expected_end
            {
                return Err(SarError::Processing(format!(
                    "DC/FM provider for {} does not cover burst [{:.6}, {:.6}] s (expected_end {:.6}); provider [{:.6}, {:.6}]",
                    sw.id,
                    burst.azimuth_time_start,
                    burst.azimuth_time_end,
                    expected_end,
                    start.value(),
                    end.value()
                )));
            }
        }
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn test_swath() -> SubSwath {
        SubSwath {
            id: "IW1".to_string(),
            burst_count: 1,
            lines_per_burst: 50,
            range_samples: 100,
            azimuth_samples: 50,
            first_line_global: 0,
            last_line_global: 50,
            first_sample_global: 0,
            last_sample_global: 100,
            full_range_samples: 100,
            valid_first_line: Some(0),
            valid_last_line: Some(50),
            valid_first_sample: Some(0),
            valid_last_sample: Some(100),
            range_pixel_spacing: 2.3,
            azimuth_pixel_spacing: 14.0,
            slant_range_time: 0.005,
            burst_duration: 2.7,
            near_range_m: 800000.0,
            prf_hz: Some(1680.0),
            dc_polynomial: Some(vec![100.0, 0.1, 0.001]),
            azimuth_time_interval: Some(0.000595),
            dc_polynomial_t0: Some(1000.0),
            fm_rate_estimates: None,
        }
    }

    #[test]
    fn validate_dc_prereqs_passes_valid() {
        let sw = test_swath();
        let result = validate_dc_prerequisites(&[sw]);
        assert!(result.is_ok());
    }

    #[test]
    fn validate_dc_prereqs_fails_missing_poly() {
        let mut sw = test_swath();
        sw.dc_polynomial = None;
        let result = validate_dc_prerequisites(&[sw]);
        assert!(result.is_err());
    }

    #[test]
    fn validate_dc_prereqs_fails_missing_t0() {
        let mut sw = test_swath();
        sw.dc_polynomial_t0 = None;
        let result = validate_dc_prerequisites(&[sw]);
        assert!(result.is_err());
    }

    #[test]
    fn fm_rate_uses_metadata_prf() {
        let sw = test_swath();
        let fm = calculate_azimuth_fm_rate(&sw, 7500.0).unwrap();

        // With PRF=1680, V=7500, L=12.3:
        // K_a = 2 * 1680² / (7500 * 12.3) ≈ 61.2 Hz/s
        assert!(
            fm > 50.0 && fm < 70.0,
            "FM rate should be ~61 Hz/s, got {}",
            fm
        );
    }

    #[test]
    fn fm_rate_errors_when_no_prf() {
        let mut sw = test_swath();
        sw.prf_hz = None;

        // Scientific mode requires annotation-derived PRF; missing PRF must error.
        assert!(calculate_azimuth_fm_rate(&sw, 7500.0).is_err());
    }
}
