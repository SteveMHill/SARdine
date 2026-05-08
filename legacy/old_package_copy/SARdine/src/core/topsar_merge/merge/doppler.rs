//! Doppler centroid and FM rate evaluation methods for TOPSAR merge.
//!
//! This module provides methods for evaluating Doppler centroid frequencies
//! and FM (frequency modulation) rates, which are essential for phase-coherent
//! TOPSAR processing and sub-swath merging.

use super::TopsarMerge;
use crate::core::geometry::type_safe_units::Seconds;
use crate::types::{SarError, SarResult, SubSwath};

/// Sentinel-1 antenna length in azimuth (meters)
const SENTINEL1_ANTENNA_LENGTH_AZIMUTH_M: f64 = 12.3;

impl TopsarMerge {
    /// Evaluate Doppler centroid for a subswath at absolute azimuth time using the provider map
    pub(super) fn evaluate_dc_hz(&self, sw: &SubSwath, azimuth_time: f64) -> SarResult<f64> {
        let provider = self.dc_fm_providers.get(&sw.id).ok_or_else(|| {
            SarError::Processing(format!("Missing DC/FM provider for subswath {}", sw.id))
        })?;

        let (start, end) = provider.get_time_range();
        let clamped = Seconds::new(azimuth_time.clamp(start.value(), end.value()));
        let dc = provider.get_dc(clamped)?;
        Ok(dc.value())
    }

    /// Fetch FM rate for a swath with clamping+logging
    pub(super) fn fm_rate_for_swath(&self, swath_idx: usize, azimuth_time: f64) -> SarResult<f64> {
        let swath = self.subswaths.get(swath_idx).ok_or_else(|| {
            SarError::Processing(format!("Invalid swath index {} for FM lookup", swath_idx))
        })?;
        let provider = self.dc_fm_providers.get(&swath.id).ok_or_else(|| {
            SarError::Processing(format!("Missing DC/FM provider for subswath {}", swath.id))
        })?;

        let (start, end) = provider.get_time_range();
        let clamped = Seconds::new(azimuth_time.clamp(start.value(), end.value()));
        provider.get_fm_rate(clamped).map_err(|e| {
            SarError::Processing(format!(
                "FM provider failure for {} at az_time={:.6}s (range [{:.6}, {:.6}]): {}",
                swath.id,
                azimuth_time,
                start.value(),
                end.value(),
                e
            ))
        })
    }

    /// Fetch Doppler centroid for a swath with clamping+logging
    pub(super) fn doppler_centroid_for_swath(
        &self,
        swath_idx: usize,
        azimuth_time: f64,
    ) -> SarResult<f64> {
        let swath = self.subswaths.get(swath_idx).ok_or_else(|| {
            SarError::Processing(format!("Invalid swath index {} for DC lookup", swath_idx))
        })?;
        let provider = self.dc_fm_providers.get(&swath.id).ok_or_else(|| {
            SarError::Processing(format!("Missing DC/FM provider for subswath {}", swath.id))
        })?;

        let (start, end) = provider.get_time_range();
        let clamped = Seconds::new(azimuth_time.clamp(start.value(), end.value()));
        provider.get_dc(clamped).map(|hz| hz.value()).map_err(|e| {
            SarError::Processing(format!(
                "DC provider failure for {} at az_time={:.6}s (range [{:.6}, {:.6}]): {}",
                swath.id,
                azimuth_time,
                start.value(),
                end.value(),
                e
            ))
        })
    }

    /// Compute azimuth FM product (K_a * eta) for a swath
    pub(super) fn azimuth_fm(&self, swath_idx: usize, azimuth_time: f64, line_idx: usize) -> f64 {
        let swath = match self.subswaths.get(swath_idx) {
            Some(sw) => sw,
            None => return 0.0,
        };

        let burst = match self
            .output_grid
            .azimuth_timing
            .burst_for_line(&swath.id, line_idx)
        {
            Some(b) => b,
            None => return 0.0,
        };

        let eta = azimuth_time - burst.sensing_time_center;
        let fm_rate = match self.fm_rate_for_swath(swath_idx, azimuth_time) {
            Ok(v) => v,
            Err(e) => {
                log::error!("{}", e);
                return 0.0;
            }
        };

        fm_rate * eta
    }

    /// Calculate azimuth FM rate for TOPSAR steering correction
    /// Implements: K_a = 2 * PRF² / (V_s * L_antenna)
    /// This is essential for phase-coherent TOPSAR processing
    pub(super) fn calculate_azimuth_fm_rate(
        subswath: &SubSwath,
        satellite_velocity: f64,
    ) -> SarResult<f64> {
        // Calculate PRF from pixel spacing
        let prf = satellite_velocity / subswath.azimuth_pixel_spacing;

        // Typical Sentinel-1 antenna length in azimuth
        let antenna_length_azimuth = SENTINEL1_ANTENNA_LENGTH_AZIMUTH_M; // meters (Sentinel-1 specification)

        // Calculate azimuth FM rate using SAR theory
        // K_a = 2 * PRF² / (V_s * L_antenna)
        let azimuth_fm_rate = 2.0 * prf.powi(2) / (satellite_velocity * antenna_length_azimuth);

        log::debug!("🔬 Azimuth FM rate calculation: PRF={:.1} Hz, V_s={:.1} m/s, L_ant={:.1} m → K_a={:.2e} Hz/s",
                   prf, satellite_velocity, antenna_length_azimuth, azimuth_fm_rate);

        Ok(azimuth_fm_rate)
    }
}
