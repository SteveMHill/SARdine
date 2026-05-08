#![allow(dead_code, unused_variables)]
//! Validation functions for terrain correction processing.
//!
//! This module contains validation logic for:
//! - Output spacing and pixel size calculations
//! - Processing parameters (wavelength, spacing, PRF)
//! - Orbit data quality and coverage
//! - SAR imaging parameters

use super::range_doppler::RangeDopplerParams;
use super::TerrainCorrector;
use crate::types::{OrbitData, SarError, SarResult};

impl TerrainCorrector {
    /// Validate and potentially fix output spacing for geographic coordinate systems.
    ///
    /// For projected CRS (non-4326), output spacing is used directly in meters.
    /// For geographic CRS (4326), validates that pixel size in degrees is correctly
    /// calculated from the target resolution in meters.
    pub fn validate_and_fix_output_spacing(
        &mut self,
        target_resolution_m: f64,
        scene_center_lat: f64,
        _scene_center_lon: f64,
    ) -> SarResult<f64> {
        use super::config::TerrainCorrectionConfig;

        if self.output_crs != 4326 {
            self.output_spacing = target_resolution_m;
            log::info!(
                "🔧 Projected output spacing fixed at {:.3} m for EPSG:{}",
                self.output_spacing,
                self.output_crs
            );
            return Ok(self.output_spacing);
        }

        // Calculate what the pixel size should be based on target resolution
        let expected_pixel_size = TerrainCorrectionConfig::calculate_pixel_size_degrees(
            target_resolution_m,
            scene_center_lat,
        );

        // Calculate what the current output_spacing would produce
        let center_lat_rad = scene_center_lat.to_radians();
        let a = crate::constants::geodetic::WGS84_SEMI_MAJOR_AXIS_M;
        let e2 = crate::constants::geodetic::WGS84_ECCENTRICITY_SQUARED;
        let sin_lat = center_lat_rad.sin();
        let prime_vertical_radius = a / (1.0 - e2 * sin_lat * sin_lat).sqrt();
        let meters_per_degree_lon =
            prime_vertical_radius * center_lat_rad.cos() * std::f64::consts::PI / 180.0;
        let current_pixel_size = self.output_spacing / meters_per_degree_lon;

        // Check for the critical bug (factor of ~73,000 error)
        let ratio = current_pixel_size / expected_pixel_size;

        log::info!("🔍 PIXEL SIZE VALIDATION:");
        log::info!("   🎯 Target resolution: {:.1}m", target_resolution_m);
        log::info!(
            "   📐 Expected pixel size: {:.8}° ({:.6} arcsec)",
            expected_pixel_size,
            expected_pixel_size * 3600.0
        );
        log::info!("   📊 Current output_spacing: {:.6}m", self.output_spacing);
        log::info!(
            "   📏 Resulting pixel size: {:.12}° ({:.8} arcsec)",
            current_pixel_size,
            current_pixel_size * 3600.0
        );
        log::info!("   ⚖️  Ratio (current/expected): {:.2e}", ratio);

        if ratio < 1e-4 || ratio > 1e4 {
            log::error!("🚨 CRITICAL GEOREFERENCING BUG DETECTED!");
            log::error!(
                "   📉 Pixel size error factor: {:.0}x",
                1.0 / ratio.min(1.0 / ratio)
            );
            log::error!(
                "   🔧 Correcting output_spacing from {:.6}m to {:.1}m",
                self.output_spacing,
                target_resolution_m
            );

            // FIX: Set output_spacing to target resolution
            self.output_spacing = target_resolution_m;

            // Recalculate corrected pixel size
            let corrected_pixel_size = self.output_spacing / meters_per_degree_lon;
            log::info!(
                "✅ CORRECTED pixel size: {:.8}° ({:.6} arcsec)",
                corrected_pixel_size,
                corrected_pixel_size * 3600.0
            );

            Ok(corrected_pixel_size)
        } else {
            log::info!("✅ Pixel size calculation is correct");
            Ok(current_pixel_size)
        }
    }

    /// Comprehensive parameter validation using the validation framework.
    ///
    /// Validates processing parameters against known physical constraints:
    /// - Radar frequency must come from annotation metadata (no hardcoded values)
    /// - Wavelength must match frequency
    /// - Pixel spacings must be physically reasonable
    /// - PRF must be within valid range
    /// - Output spacing and CRS must be supported
    pub fn validate_processing_parameters(
        &self,
        params: &RangeDopplerParams,
        metadata: &crate::types::SarMetadata,
    ) -> SarResult<()> {
        use crate::validation::ParameterValidator;

        let validator = ParameterValidator::new();

        // Extract radar frequency from metadata - NO hardcoded values permitted
        let radar_frequency = metadata.radar_frequency.ok_or_else(|| {
            SarError::ParameterError("Radar frequency must be extracted from annotation XML metadata; no hardcoded values permitted".to_string())
        })?;

        // Validate all parameters comprehensively using real annotation data
        validator.validate_all_parameters(
            radar_frequency, // Real frequency from annotation XML (not 5.405e9 hardcode)
            params.wavelength,
            params.range_pixel_spacing,
            params.azimuth_pixel_spacing,
            params.prf,
            "annotation XML",
        )?;

        // Additional terrain correction specific validations
        if self.output_spacing <= 0.0 || self.output_spacing > 1000.0 {
            return Err(SarError::InvalidParameter(format!(
                "Invalid output spacing: {:.3}m. Must be positive and ≤1000m",
                self.output_spacing
            )));
        }

        if self.output_crs != 4326 && !(self.output_crs >= 32601 && self.output_crs <= 32760) {
            return Err(SarError::InvalidParameter(format!(
                "Unsupported output CRS: EPSG:{}. Must be WGS84 (4326) or UTM (32601-32760)",
                self.output_crs
            )));
        }

        log::info!("✅ All processing parameters validated successfully");
        Ok(())
    }

    /// Validate orbit data quality and coverage.
    ///
    /// Checks:
    /// - Minimum number of state vectors for interpolation (≥3)
    /// - Orbit time span coverage (≥60 seconds)
    /// - Orbital velocity magnitude within expected range (7000-8000 m/s for LEO)
    ///
    /// NOTE: Minimum 3 vectors required for cubic spline interpolation (mathematical minimum).
    /// Python layer enforces MINIMUM_ORBIT_VECTORS=10 for scientific processing (ensures adequate
    /// temporal coverage). Precise orbit files typically have 9000+ vectors.
    pub(super) fn validate_orbit_data(&self, orbit_data: &OrbitData) -> SarResult<()> {
        // Enforce scientific minimum: 10 state vectors for adequate temporal coverage.
        // Cubic spline needs 3, but 10 ensures proper interpolation quality across the scene.
        // Precise orbit files typically have 9000+ vectors.
        if orbit_data.state_vectors.len() < 10 {
            return Err(SarError::Processing(format!(
                "Insufficient orbit vectors for terrain correction: {} (minimum 10 required for scientific quality)",
                orbit_data.state_vectors.len()
            )));
        }

        // Check orbit vector spacing (should be ~10 seconds for Sentinel-1)
        // AUDIT FIX: Safe access after length check above
        let last_sv = orbit_data
            .state_vectors
            .last()
            .ok_or_else(|| SarError::Processing("No orbit vectors available".to_string()))?;
        let first_sv = orbit_data
            .state_vectors
            .first()
            .ok_or_else(|| SarError::Processing("No orbit vectors available".to_string()))?;
        let time_span = last_sv.time - first_sv.time;
        let time_span_seconds = time_span.num_seconds() as f64;

        if time_span_seconds < 60.0 {
            return Err(SarError::Processing(format!(
                "Orbit time span too short: {:.1}s (minimum 60s required)",
                time_span_seconds
            )));
        }

        // Validate orbital velocity magnitude (Sentinel-1A: ~7.5 km/s)
        for vector in &orbit_data.state_vectors {
            let velocity_magnitude = (vector.velocity[0].powi(2)
                + vector.velocity[1].powi(2)
                + vector.velocity[2].powi(2))
            .sqrt();

            if velocity_magnitude < 7000.0 || velocity_magnitude > 8000.0 {
                log::warn!(
                    "Unusual orbital velocity: {:.1} m/s (expected 7000-8000 m/s)",
                    velocity_magnitude
                );
            }
        }

        log::info!(
            "✅ Orbit data validation passed: {} vectors, {:.1}s span",
            orbit_data.state_vectors.len(),
            time_span_seconds
        );
        Ok(())
    }

    /// Validate SAR imaging parameters.
    ///
    /// Checks parameter ranges against typical SAR sensor specifications:
    /// - Range pixel spacing: 1-10m (IW mode ~2.3m)
    /// - Azimuth pixel spacing: 5-50m (IW mode ~14m)
    /// - Wavelength: 0.015-0.30m (Ku to L band)
    pub(super) fn validate_sar_parameters(&self, params: &RangeDopplerParams) -> SarResult<()> {
        // Validate range pixel spacing (Sentinel-1 IW: ~2.3m)
        if params.range_pixel_spacing < 1.0 || params.range_pixel_spacing > 10.0 {
            log::warn!(
                "Unusual range pixel spacing: {:.2}m (expected 1-10m)",
                params.range_pixel_spacing
            );
        }

        // Validate azimuth pixel spacing (Sentinel-1 IW: ~14m)
        if params.azimuth_pixel_spacing < 5.0 || params.azimuth_pixel_spacing > 50.0 {
            log::warn!(
                "Unusual azimuth pixel spacing: {:.2}m (expected 5-50m)",
                params.azimuth_pixel_spacing
            );
        }

        // Validate wavelength (typical SAR bands: L=0.24m, S=0.10m, C=0.055m, X=0.031m, Ku=0.019m)
        let min_wavelength = 0.015; // Ku-band lower limit
        let max_wavelength = 0.30; // L-band upper limit
        if params.wavelength < min_wavelength || params.wavelength > max_wavelength {
            return Err(SarError::Processing(format!(
                "Invalid wavelength: {:.3}m (expected {:.3}-{:.3}m)",
                params.wavelength, min_wavelength, max_wavelength
            )));
        }

        log::info!("✅ SAR parameters validation passed");
        Ok(())
    }
}
