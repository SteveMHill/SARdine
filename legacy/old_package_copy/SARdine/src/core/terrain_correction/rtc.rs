#![allow(dead_code, unused_variables)]
//! Radiometric Terrain Correction (RTC) functions
//!
//! This module contains functions for computing radiometric terrain correction factors,
//! which are essential for accurate SAR backscatter analysis over terrain.
//!
//! # RTC Modes
//!
//! Two correction methods are implemented:
//!
//! 1. **Area-Projection (Small 2011)** - DEFAULT
//!    - γ⁰ = σ⁰ × (A_flat / A_slope)
//!    - Physically correct normalization using projected pixel areas
//!    - Accounts for true terrain surface area vs ellipsoidal area
//!    - Recommended for scientific applications
//!
//! 2. **Cosine Local Incidence Angle**
//!    - γ⁰ = σ⁰ / cos(θ_local)
//!    - Simpler approximation, faster computation
//!    - Good for flat-to-moderate terrain
//!    - May be less accurate on steep slopes
//!
//! # Key Functions
//! - `compute_surface_normal` - Computes surface normal vectors from DEM gradients
//! - `look_vector_at` - Computes radar look vector at a ground point
//! - `rtc_gamma0_scale` - Computes gamma0 RTC scale factor (cosine method)
//! - `rtc_area_projection_scale` - Computes gamma0 RTC scale factor (area projection)
//! - `detect_shadow_layover` - Detects shadow and layover areas
//! - `compute_local_incidence_angle` - Computes local incidence angle
//! - `compute_rtc_for_pixel` - Unified RTC computation for geocoding integration
//!
//! # Literature References
//! - Small, D. (2011): "Flattening Gamma: Radiometric Terrain Correction for SAR Imagery"
//!   IEEE Transactions on Geoscience and Remote Sensing, Vol. 49, No. 8, pp. 3081-3093
//! - Zebker & Goldstein (1986): "Topographic mapping from interferometric SAR"
//! - Bamler & Hartl (1998): "Synthetic aperture radar interferometry"
//! - ESA Sentinel-1 Level 1 Product Specification

use crate::constants::geodetic::WGS84_SEMI_MAJOR_AXIS_M;
use crate::core::geometry::type_safe_units::Meters;
use crate::types::{OrbitData, SarError, SarResult, SurfaceNormal};
use ndarray::Array2;

/// Result of RTC computation for a single pixel
///
/// Contains all outputs needed for terrain-corrected backscatter and quality assessment.
#[derive(Debug, Clone, Copy)]
pub struct RtcResult {
    /// RTC scale factor to apply to calibrated σ⁰ to get γ⁰
    /// γ⁰ = σ⁰ × scale
    pub scale: f32,

    /// Local incidence angle in degrees (0-90)
    pub local_incidence_angle_deg: f32,

    /// Quality flag (bitwise OR of RtcQualityFlag values)
    pub quality_flag: u8,

    /// True if pixel is in radar shadow
    pub is_shadow: bool,

    /// True if pixel is in layover zone
    pub is_layover: bool,
}

impl Default for RtcResult {
    fn default() -> Self {
        Self {
            scale: 1.0,
            local_incidence_angle_deg: 0.0,
            quality_flag: 0,
            is_shadow: false,
            is_layover: false,
        }
    }
}

/// Extended Range-Doppler transformation result with geometric information for RTC
///
/// This struct contains the SAR coordinates plus all geometric quantities needed
/// for radiometric terrain correction, computed during the Range-Doppler solve.
#[derive(Debug, Clone)]
pub struct RangeDopplerWithGeometry {
    /// Native range pixel coordinate
    pub range_pixel_native: f64,

    /// Native azimuth pixel coordinate
    pub azimuth_pixel_native: f64,

    /// Look vector from target to satellite (unit vector, ECEF)
    pub look_vector: super::types::Vector3,

    /// Target position in ECEF coordinates
    pub target_ecef: [f64; 3],

    /// Satellite position at zero-Doppler time (ECEF)
    pub sat_position: super::types::Position3D,

    /// Slant range in meters
    pub slant_range: f64,
}

/// RTC mode selection
///
/// Determines which radiometric terrain correction algorithm to use.
/// Default is AreaProjection (Small 2011) for scientific accuracy.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum RtcMode {
    /// Area-Projection method (Small 2011) - DEFAULT
    ///
    /// γ⁰ = σ⁰ × (A_flat / A_slope)
    ///
    /// Computes the ratio of ellipsoidal pixel area to terrain surface area.
    /// This is the physically correct approach that accounts for true terrain geometry.
    ///
    /// Area ratio formula: A_slope/A_flat = 1 / (cos(θ_local) × √(1 + (∂z/∂r)² + (∂z/∂a)²))
    /// Simplified: A_flat/A_slope = cos(θ_local) × √(1 + slope²) for moderate slopes
    ///
    /// Recommended for:
    /// - Scientific analysis requiring accurate γ⁰ values
    /// - Mountainous or complex terrain
    /// - Quantitative backscatter studies
    #[default]
    AreaProjection,

    /// Cosine Local Incidence Angle method
    ///
    /// γ⁰ = σ⁰ / cos(θ_local)
    ///
    /// Simple division by cosine of local incidence angle.
    /// Faster but less physically accurate than area projection.
    ///
    /// Recommended for:
    /// - Quick processing where speed matters more than accuracy
    /// - Relatively flat terrain
    /// - Qualitative analysis
    CosineLocalIncidenceAngle,

    /// SNAP-compatible Sigma-Sin method
    ///
    /// σ⁰_TC = σ⁰ × sin(θ_DEM) / sin(θ_ellipsoid)
    ///
    /// This matches the default behavior of SNAP's Range-Doppler Terrain Correction
    /// operator when using "Use projected local incidence angle from DEM".
    ///
    /// The formula normalizes by the ratio of sines of incidence angles:
    /// - θ_DEM: Local incidence angle computed from DEM surface normal
    /// - θ_ellipsoid: Reference incidence angle on the ellipsoid (from metadata)
    ///
    /// Use this mode when:
    /// - Exact compatibility with SNAP outputs is required
    /// - Comparing results with SNAP-processed products
    /// - Workflows that expect SNAP-style terrain correction
    ///
    /// Note: AreaProjection (Small 2011) is more physically rigorous, but this
    /// mode produces outputs that match SNAP for validation/comparison purposes.
    SnapSigmaSin,
}

impl RtcMode {
    /// Parse RTC mode from environment variable or string
    pub fn from_env() -> Self {
        std::env::var("SARDINE_RTC_MODE")
            .ok()
            .and_then(|s| Self::from_str(&s))
            .unwrap_or_default()
    }

    /// Parse from string (case-insensitive)
    pub fn from_str(s: &str) -> Option<Self> {
        match s.to_lowercase().as_str() {
            "area" | "area-projection" | "area_projection" | "small2011" | "small" => {
                Some(Self::AreaProjection)
            }
            "cosine" | "lia" | "cos" | "cosine-lia" | "cosine_lia" => {
                Some(Self::CosineLocalIncidenceAngle)
            }
            "snap" | "snap-sigma-sin" | "snap_sigma_sin" | "sigma-sin" | "sigma_sin" | "sin" => {
                Some(Self::SnapSigmaSin)
            }
            _ => None,
        }
    }

    /// Get human-readable name
    pub fn name(&self) -> &'static str {
        match self {
            Self::AreaProjection => "Area-Projection (Small 2011)",
            Self::CosineLocalIncidenceAngle => "Cosine Local Incidence Angle",
            Self::SnapSigmaSin => "SNAP-compatible Sigma-Sin",
        }
    }
}

/// Quality flags for RTC pixels
///
/// Bitwise flags indicating quality conditions for each pixel.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[repr(u8)]
pub enum RtcQualityFlag {
    /// Pixel is valid with no quality issues
    Valid = 0,
    /// Scale factor was clamped low (potential layover)
    ClampedLow = 1,
    /// Scale factor was clamped high (steep slope / shadow edge)
    ClampedHigh = 2,
    /// Pixel is in radar shadow
    Shadow = 4,
    /// Pixel is in layover zone
    Layover = 8,
    /// DEM has no data at this pixel
    DemNoData = 16,
}

use super::range_doppler::RangeDopplerParams;
use super::types::Vector3;
use super::TerrainCorrector;

impl TerrainCorrector {
    /// Compute surface normal from DEM gradient at a pixel location
    ///
    /// Uses central differences for gradient estimation with proper handling of
    /// geographic vs projected coordinate systems.
    ///
    /// # Scientific Implementation
    /// 1. Compute partial derivatives dz/dx and dz/dy using central differences
    /// 2. Construct tangent basis vectors in meters
    /// 3. Compute normal via cross product: n = normalize(u × v)
    ///
    /// # Arguments
    /// * `dem_array` - DEM elevation data
    /// * `row` - Pixel row index
    /// * `col` - Pixel column index
    ///
    /// # Returns
    /// * Unit surface normal vector
    pub fn compute_surface_normal(
        &self,
        dem_array: &Array2<f32>,
        row: usize,
        col: usize,
    ) -> SurfaceNormal {
        let (height, width) = dem_array.dim();

        // CRITICAL FIX: Use DEM pixel spacing in meters, handling geographic vs projected CRS
        // Based on expert recommendations for scientifically accurate slope computation
        let (sx, sy) = if self.dem_crs == 4326 {
            // Geographic DEM (degrees) - convert to meters at pixel latitude
            let pixel_lat =
                self.dem_transform.top_left_y + (row as f64) * self.dem_transform.pixel_height;

            // Use WGS84 ellipsoid formulas (already implemented for pixel size calculation)
            let lat_rad = pixel_lat.to_radians();
            let sin_lat = lat_rad.sin();
            let cos_lat = lat_rad.cos();

            // Prime vertical radius of curvature
            let n = WGS84_SEMI_MAJOR_AXIS_M
                / (1.0
                    - crate::constants::geodetic::WGS84_ECCENTRICITY_SQUARED * sin_lat * sin_lat)
                    .sqrt();

            // Convert degrees to meters at this latitude
            let meters_per_deg_lon = n * cos_lat * std::f64::consts::PI / 180.0;
            let meters_per_deg_lat = n
                * (1.0 - crate::constants::geodetic::WGS84_ECCENTRICITY_SQUARED)
                / (1.0
                    - crate::constants::geodetic::WGS84_ECCENTRICITY_SQUARED * sin_lat * sin_lat)
                    .powf(1.5)
                * std::f64::consts::PI
                / 180.0;

            (
                Meters::new(self.dem_transform.pixel_width.abs() * meters_per_deg_lon),
                Meters::new(self.dem_transform.pixel_height.abs() * meters_per_deg_lat),
            )
        } else {
            // Projected DEM (already in meters)
            (
                Meters::new(self.dem_transform.pixel_width.abs()),
                Meters::new(self.dem_transform.pixel_height.abs()),
            )
        };

        let sx_val = sx.value();
        let sy_val = sy.value();

        // Central differences with boundary handling - same algorithm but correct scaling
        let dz_dx = if col == 0 {
            (dem_array[[row, col + 1]] - dem_array[[row, col]]) / sx_val as f32
        } else if col == width - 1 {
            (dem_array[[row, col]] - dem_array[[row, col - 1]]) / sx_val as f32
        } else {
            (dem_array[[row, col + 1]] - dem_array[[row, col - 1]]) / (2.0 * sx_val as f32)
        };

        let dz_dy = if row == 0 {
            (dem_array[[row + 1, col]] - dem_array[[row, col]]) / sy_val as f32
        } else if row == height - 1 {
            (dem_array[[row, col]] - dem_array[[row - 1, col]]) / sy_val as f32
        } else {
            (dem_array[[row + 1, col]] - dem_array[[row - 1, col]]) / (2.0 * sy_val as f32)
        };

        // Tangent basis vectors (meters) - following expert RTC recommendations
        let ux = SurfaceNormal {
            x: sx_val,
            y: 0.0,
            z: (dz_dx as f64) * sx_val,
        };
        let vy = SurfaceNormal {
            x: 0.0,
            y: sy_val,
            z: (dz_dy as f64) * sy_val,
        };

        // Normal = normalize(u × v) - cross product gives true surface normal
        let nx = ux.y * vy.z - ux.z * vy.y;
        let ny = ux.z * vy.x - ux.x * vy.z;
        let nz = ux.x * vy.y - ux.y * vy.x;

        let mut normal = SurfaceNormal {
            x: nx,
            y: ny,
            z: nz,
        };
        normal.normalize();
        normal
    }

    /// Compute true radar look vector at a ground point
    ///
    /// This replaces the placeholder vertical look vector with physically accurate
    /// sensor-to-target vector computation based on expert RTC recommendations.
    ///
    /// # Scientific Implementation
    /// 1. Ground ECEF: Convert (lat, lon, h) to Earth-Centered Earth-Fixed coordinates
    /// 2. Zero-Doppler time: Find azimuth time when target is closest to satellite
    /// 3. Satellite position: Interpolate orbit state at zero-Doppler time
    /// 4. Look vector: Unit vector from target to sensor = normalize(X_s - X_t)
    ///
    /// # Literature References
    /// - Zebker & Goldstein (1986): "Topographic mapping from interferometric SAR"
    /// - Bamler & Hartl (1998): "Synthetic aperture radar interferometry"
    /// - ESA Sentinel-1 Level 1 Product Specification, Section 4.2.1
    ///
    /// # Arguments
    /// * `lat` - Target latitude in degrees
    /// * `lon` - Target longitude in degrees
    /// * `h` - Target height in meters (ellipsoidal)
    /// * `orbit` - Orbit state vectors
    /// * `params` - Range-Doppler parameters
    ///
    /// # Returns
    /// * Unit look vector from target to sensor
    pub fn look_vector_at(
        &self,
        lat: f64,
        lon: f64,
        h: f64,
        orbit: &OrbitData,
        params: &RangeDopplerParams,
    ) -> SarResult<Vector3> {
        // 1) Ground ECEF coordinates
        let xt = self.latlon_to_ecef(lat, lon, h);

        // 2) Zero-Doppler time (use unified solver)
        let t0 = self
            .solve_zero_doppler_default(&xt, orbit, params)
            .ok_or_else(|| {
                SarError::Processing(format!(
                    "Zero-Doppler failed for target ({:.6}, {:.6}, {:.1})",
                    lat, lon, h
                ))
            })?;

        // 3) Satellite position at zero-Doppler time
        let absolute_t0 = t0 + params.product_start_absolute();
        let (xs, _vs) = self.scientific_orbit_interpolation(orbit, absolute_t0)?;

        // 4) Unit look vector l = normalize(X_s - X_t)
        let dx = xs.x - xt[0];
        let dy = xs.y - xt[1];
        let dz = xs.z - xt[2];
        let norm = (dx * dx + dy * dy + dz * dz).sqrt();

        if norm <= 0.0 {
            return Err(SarError::Processing(format!(
                "Degenerate look vector at ({:.6}, {:.6}, {:.1})",
                lat, lon, h
            )));
        }

        Ok(Vector3 {
            x: dx / norm,
            y: dy / norm,
            z: dz / norm,
        })
    }

    /// Apply simple cosine-based gamma naught RTC scaling
    ///
    /// Based on expert RTC recommendations for robust "flattening" appearance.
    /// This is NOT true area normalization but gives stable results over terrain.
    ///
    /// # Scientific Implementation
    /// scale = cos(θ_ref) / max(cos(θ_lia), ε)
    /// I_rtc = I_calibrated * scale
    ///
    /// where:
    /// - θ_ref: reference incidence angle (ellipsoid or scene-average)
    /// - θ_lia: local incidence angle (from surface normal and look vector)
    /// - ε: small value to avoid wild gains in near-shadow areas
    ///
    /// # Literature References
    /// - Small, D. (2011): "Flattening Gamma: Radiometric Terrain Correction"
    /// - ESA Sentinel-1 Level 1 Product Specification, Section 4.4.2
    /// - GAMMA Software documentation: "Radiometric Terrain Correction"
    ///
    /// # Arguments
    /// * `cos_lia` - Cosine of local incidence angle [0, 1]
    /// * `cos_ref` - Cosine of reference incidence angle [0, 1]
    ///
    /// # Returns
    /// * RTC scale factor (clamped to configured bounds for stability)
    pub fn rtc_gamma0_scale(cos_lia: f64, cos_ref: f64) -> f32 {
        // SCIENTIFIC FIX: Use wider bounds to preserve radiometric accuracy.
        // Previous [0.1, 10.0] caused 20 dB bias on slopes >75°.
        // New [0.01, 100.0] = 40 dB range, sufficient for slopes up to ~89°.
        // Pixels hitting bounds should be quality-flagged, not silently clipped.
        let eps = 1e-4; // Prevent division by zero (cos(89.994°) ≈ 1e-4)
        let scale = cos_ref / cos_lia.max(eps);

        // Configurable bounds via environment (default: 0.01 to 100.0 = 40 dB range)
        let mut min_scale: f64 = std::env::var("SARDINE_RTC_MIN_SCALE")
            .ok()
            .and_then(|v| v.parse().ok())
            .unwrap_or(0.01);
        let mut max_scale: f64 = std::env::var("SARDINE_RTC_MAX_SCALE")
            .ok()
            .and_then(|v| v.parse().ok())
            .unwrap_or(100.0);

        // Guard against inverted bounds (e.g. user sets MIN > MAX)
        if min_scale >= max_scale {
            log::error!(
                "RTC scale bounds inverted: min={} >= max={} — swapping to prevent silent corruption",
                min_scale, max_scale
            );
            std::mem::swap(&mut min_scale, &mut max_scale);
        }

        // Track clamping for quality assessment (log first few occurrences)
        // These counters should be reset via reset_rtc_counters() at start of each run
        static CLAMP_LOW_COUNT: std::sync::atomic::AtomicUsize =
            std::sync::atomic::AtomicUsize::new(0);
        static CLAMP_HIGH_COUNT: std::sync::atomic::AtomicUsize =
            std::sync::atomic::AtomicUsize::new(0);

        if scale < min_scale {
            let count = CLAMP_LOW_COUNT.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
            if count < 10 {
                log::warn!(
                    "⚠️  RTC scale {:.4} clamped to {:.4} (cos_lia={:.4}, cos_ref={:.4}) - potential layover",
                    scale, min_scale, cos_lia, cos_ref
                );
            }
        } else if scale > max_scale {
            let count = CLAMP_HIGH_COUNT.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
            if count < 10 {
                log::warn!(
                    "⚠️  RTC scale {:.4} clamped to {:.4} (cos_lia={:.4}, cos_ref={:.4}) - steep slope/shadow edge",
                    scale, max_scale, cos_lia, cos_ref
                );
            }
        }

        scale.clamp(min_scale, max_scale) as f32
    }

    /// Extended RTC scale computation with quality flag output
    ///
    /// Returns both the clamped scale and a quality flag indicating if clamping occurred.
    /// Use this version when building quality masks.
    ///
    /// # Quality Flags
    /// * 0 = valid (no clamping)
    /// * 1 = scale clamped low (potential layover)
    /// * 2 = scale clamped high (steep slope/shadow edge)
    pub fn rtc_gamma0_scale_with_quality(cos_lia: f64, cos_ref: f64) -> (f32, u8) {
        let eps = 1e-4;
        let scale = cos_ref / cos_lia.max(eps);

        let min_scale: f64 = std::env::var("SARDINE_RTC_MIN_SCALE")
            .ok()
            .and_then(|v| v.parse().ok())
            .unwrap_or(0.01);
        let max_scale: f64 = std::env::var("SARDINE_RTC_MAX_SCALE")
            .ok()
            .and_then(|v| v.parse().ok())
            .unwrap_or(100.0);

        let quality_flag = if scale < min_scale {
            1u8 // Clamped low
        } else if scale > max_scale {
            2u8 // Clamped high
        } else {
            0u8 // Valid
        };

        (scale.clamp(min_scale, max_scale) as f32, quality_flag)
    }

    /// SNAP-compatible Sigma-Sin RTC scale computation
    ///
    /// This implements the same algorithm used by SNAP's Range-Doppler Terrain Correction
    /// operator when "Use projected local incidence angle from DEM" is enabled (default).
    ///
    /// # Formula (SNAP RangeDopplerGeocodingOp.java)
    /// ```text
    /// σ⁰_TC = σ⁰ × sin(θ_DEM) / sin(θ_ellipsoid)
    /// ```
    ///
    /// where:
    /// - θ_DEM: Local incidence angle from DEM surface normal (θ_local)
    /// - θ_ellipsoid: Reference incidence angle on the ellipsoid (θ_ref from metadata)
    ///
    /// Since we work with cosines: sin(θ) = √(1 - cos²(θ))
    /// ```text
    /// scale = sin(θ_local) / sin(θ_ref)
    ///       = √(1 - cos²(θ_local)) / √(1 - cos²(θ_ref))
    /// ```
    ///
    /// # Arguments
    /// * `cos_lia` - Cosine of local incidence angle [0, 1]
    /// * `cos_ref` - Cosine of reference incidence angle [0, 1]
    ///
    /// # Returns
    /// * RTC scale factor for SNAP-compatible terrain correction
    ///
    /// # Reference
    /// SNAP S1TBX: s1tbx-op-sar-processing/src/main/java/org/esa/s1tbx/sar/gpf/geometric/RangeDopplerGeocodingOp.java
    pub fn rtc_snap_sigma_sin_scale(cos_lia: f64, cos_ref: f64) -> f32 {
        let eps = 1e-6;

        // sin(θ) = √(1 - cos²(θ))
        let sin_lia = (1.0 - cos_lia * cos_lia).max(0.0).sqrt();
        let sin_ref = (1.0 - cos_ref * cos_ref).max(eps).sqrt();

        // SNAP formula: scale = sin(θ_local) / sin(θ_ref)
        let scale = sin_lia / sin_ref;

        // Apply same bounds as other methods
        let min_scale: f64 = std::env::var("SARDINE_RTC_MIN_SCALE")
            .ok()
            .and_then(|v| v.parse().ok())
            .unwrap_or(0.01);
        let max_scale: f64 = std::env::var("SARDINE_RTC_MAX_SCALE")
            .ok()
            .and_then(|v| v.parse().ok())
            .unwrap_or(100.0);

        scale.clamp(min_scale, max_scale) as f32
    }

    /// SNAP-compatible Sigma-Sin RTC with quality flags
    ///
    /// Returns both the scale factor and quality information.
    /// See `rtc_snap_sigma_sin_scale` for algorithm details.
    ///
    /// # Quality Flags
    /// * 0 = valid (no clamping)
    /// * 1 = scale clamped low (steep backslope)
    /// * 2 = scale clamped high (grazing incidence)
    pub fn rtc_snap_sigma_sin_scale_with_quality(cos_lia: f64, cos_ref: f64) -> (f32, u8) {
        let eps = 1e-6;

        // sin(θ) = √(1 - cos²(θ))
        let sin_lia = (1.0 - cos_lia * cos_lia).max(0.0).sqrt();
        let sin_ref = (1.0 - cos_ref * cos_ref).max(eps).sqrt();

        // SNAP formula: scale = sin(θ_local) / sin(θ_ref)
        let scale = sin_lia / sin_ref;

        let min_scale: f64 = std::env::var("SARDINE_RTC_MIN_SCALE")
            .ok()
            .and_then(|v| v.parse().ok())
            .unwrap_or(0.01);
        let max_scale: f64 = std::env::var("SARDINE_RTC_MAX_SCALE")
            .ok()
            .and_then(|v| v.parse().ok())
            .unwrap_or(100.0);

        let quality_flag = if scale < min_scale {
            1u8 // Clamped low
        } else if scale > max_scale {
            2u8 // Clamped high
        } else {
            0u8 // Valid
        };

        (scale.clamp(min_scale, max_scale) as f32, quality_flag)
    }

    /// Compute RTC scale using Area-Projection method (Small 2011)
    ///
    /// This is the **scientifically correct** method for radiometric terrain correction.
    /// It computes the ratio of flat (ellipsoidal) pixel area to sloped terrain area.
    ///
    /// # Formula
    /// **γ⁰ = σ⁰ × (A_flat / A_slope)**
    ///
    /// where the area ratio is computed from:
    /// - cos(θ_local): cosine of local incidence angle
    /// - slope_magnitude: √(1 + (∂z/∂x)² + (∂z/∂y)²) - the terrain surface normal magnitude factor
    ///
    /// The full derivation (Small 2011):
    /// A_slope / A_flat = 1 / (cos(θ_local) × sec(slope))
    /// => A_flat / A_slope = cos(θ_local) / √(1 + (∂z/∂x)² + (∂z/∂y)²)
    ///
    /// But we need to also account for the projection onto the look direction:
    /// scale = A_flat / A_slope = |n̂ · l̂| × √(1 + slope²) / slope_factor
    ///
    /// For practical implementation, we use:
    /// scale = cos(θ_local) × √(1 + (∂z/∂x)² + (∂z/∂y)²)
    ///
    /// This represents the effective area compression due to terrain slope
    /// relative to the ellipsoidal assumption.
    ///
    /// # Arguments
    /// * `cos_lia` - Cosine of local incidence angle [0, 1]
    /// * `dz_dx` - Terrain slope in x (range) direction (dimensionless)
    /// * `dz_dy` - Terrain slope in y (azimuth) direction (dimensionless)
    /// * `cos_ref` - Cosine of reference incidence angle (for normalization)
    ///
    /// # Returns
    /// * RTC scale factor for area projection method
    ///
    /// # Literature Reference
    /// Small, D. (2011): "Flattening Gamma: Radiometric Terrain Correction for SAR Imagery"
    /// IEEE TGRS, Vol. 49, No. 8, pp. 3081-3093
    /// DOI: 10.1109/TGRS.2011.2120616
    pub fn rtc_area_projection_scale(cos_lia: f64, dz_dx: f64, dz_dy: f64, cos_ref: f64) -> f32 {
        let eps = 1e-6;

        // Terrain surface area factor: √(1 + (∂z/∂x)² + (∂z/∂y)²)
        // This represents how much larger the terrain surface is compared to horizontal
        let slope_magnitude_sq = 1.0 + dz_dx * dz_dx + dz_dy * dz_dy;
        let slope_magnitude = slope_magnitude_sq.sqrt();

        // Area ratio: A_flat / A_slope
        // For ellipsoidal reference, flat area = pixel_size²
        // Terrain area = flat_area × slope_magnitude / cos(θ_local)
        // But we want to normalize by incidence: σ⁰ is already referenced to ellipsoid
        // So the correction factor is: cos(θ_local) × slope_magnitude / cos(θ_ref)
        //
        // Wait - let's derive this carefully:
        // γ⁰ = σ⁰ × (A_flat / A_slope) = σ⁰ × (A_flat / (A_flat × slope_magnitude))
        //    = σ⁰ / slope_magnitude
        // But this doesn't account for look direction...
        //
        // Actually Small 2011 formula:
        // γ⁰ = σ⁰ × (A_ellipsoid / A_terrain) × (sin(θ_local) / sin(θ_ref))
        //
        // Simplified for our case (using cos formulation):
        // scale = (cos_ref / slope_magnitude) for terrain area normalization
        //       × (slope_magnitude × cos_lia) for projection adjustment
        //       = cos_ref × cos_lia
        //
        // But that's just the cosine method! The difference is subtle:
        // Area projection uses the FULL terrain normal including azimuth slope,
        // while cosine-LIA only uses the angle between normal and look vector.
        //
        // The correct formula accounting for both:
        // scale = cos_lia × slope_magnitude / (cos_ref × slope_magnitude_ref)
        // where slope_magnitude_ref ≈ 1 for flat reference
        //
        // So: scale ≈ cos_lia × slope_magnitude
        // And the γ⁰ = σ⁰ × (1 / scale) = σ⁰ / (cos_lia × slope_magnitude)
        //
        // But wait - this INCREASES values on slopes, opposite of what we want!
        // Let me re-derive:
        //
        // σ⁰ = backscatter per unit GROUND area
        // γ⁰ = backscatter per unit SLOPE area
        // A_slope > A_ground for sloped terrain
        // So γ⁰ = σ⁰ × (A_ground / A_slope) = σ⁰ / slope_magnitude (for range-parallel slopes)
        //
        // But incidence angle also matters:
        // γ⁰ = σ⁰ × cos(θ_ref) / cos(θ_local) × (1 / slope_magnitude)
        //
        // This is the area-projection formula:
        // scale = cos(θ_ref) / (cos(θ_local) × slope_magnitude)

        let denominator = (cos_lia.abs() + eps) * slope_magnitude;
        let scale = cos_ref / denominator;

        // Apply the same clamping bounds as cosine method
        let min_scale: f64 = std::env::var("SARDINE_RTC_MIN_SCALE")
            .ok()
            .and_then(|v| v.parse().ok())
            .unwrap_or(0.01);
        let max_scale: f64 = std::env::var("SARDINE_RTC_MAX_SCALE")
            .ok()
            .and_then(|v| v.parse().ok())
            .unwrap_or(100.0);

        scale.clamp(min_scale, max_scale) as f32
    }

    /// Area-Projection RTC with quality flags
    ///
    /// Returns both the scale factor and quality information.
    /// See `rtc_area_projection_scale` for algorithm details.
    ///
    /// # Scientific Basis
    ///
    /// Following Small (2011) "Flattening Gamma: Radiometric Terrain Correction for SAR Imagery",
    /// IEEE TGRS 49(8):3081-3093, the area-projection RTC converts sigma-nought to gamma-nought:
    ///
    /// ```text
    /// γ⁰ = σ⁰ × cos(θ_ref) / (cos(θ_local) × A_norm)
    /// ```
    ///
    /// where:
    /// - θ_ref: reference incidence angle on ellipsoid (varies with range)
    /// - θ_local: local incidence angle (terrain-dependent)  
    /// - A_norm: normalized terrain surface area = sqrt(1 + (∂z/∂x)² + (∂z/∂y)²)
    ///
    /// The scale factor normalizes the backscatter by removing terrain-induced radiometric
    /// variations, producing consistent gamma-nought values independent of local topography.
    ///
    /// # Returns
    /// * (scale, quality_flag, area_ratio)
    /// * quality_flag: 0=valid, 1=clamped_low, 2=clamped_high, 4=shadow, 8=layover
    /// * area_ratio: A_flat / A_slope ratio (for diagnostics)
    pub fn rtc_area_projection_scale_with_quality(
        cos_lia: f64,
        dz_dx: f64,
        dz_dy: f64,
        cos_ref: f64,
    ) -> (f32, u8, f64) {
        let eps = 1e-6;

        // Terrain surface area factor
        let slope_magnitude_sq = 1.0 + dz_dx * dz_dx + dz_dy * dz_dy;
        let slope_magnitude = slope_magnitude_sq.sqrt();

        // Area ratio for diagnostics
        let area_ratio = 1.0 / slope_magnitude;

        // Scale factor
        let denominator = (cos_lia.abs() + eps) * slope_magnitude;
        let scale = cos_ref / denominator;

        let min_scale: f64 = std::env::var("SARDINE_RTC_MIN_SCALE")
            .ok()
            .and_then(|v| v.parse().ok())
            .unwrap_or(0.01);
        let max_scale: f64 = std::env::var("SARDINE_RTC_MAX_SCALE")
            .ok()
            .and_then(|v| v.parse().ok())
            .unwrap_or(100.0);

        // Determine quality flag
        let mut quality_flag = 0u8;

        if cos_lia <= 0.0 {
            quality_flag |= RtcQualityFlag::Shadow as u8;
        }

        if scale < min_scale {
            quality_flag |= RtcQualityFlag::ClampedLow as u8;
        } else if scale > max_scale {
            quality_flag |= RtcQualityFlag::ClampedHigh as u8;
        }

        (
            scale.clamp(min_scale, max_scale) as f32,
            quality_flag,
            area_ratio,
        )
    }

    /// Compute DEM gradients for area projection RTC
    ///
    /// Returns (dz_dx, dz_dy) at the specified pixel location using central differences.
    ///
    /// # Arguments
    /// * `dem_array` - DEM elevation data
    /// * `row` - Pixel row index
    /// * `col` - Pixel column index
    /// * `pixel_spacing` - DEM pixel spacing in meters
    ///
    /// # Returns
    /// * (dz_dx, dz_dy) - Terrain gradients (dimensionless, rise/run)
    pub fn compute_dem_gradients(
        dem_array: &Array2<f32>,
        row: usize,
        col: usize,
        pixel_spacing: f64,
    ) -> (f64, f64) {
        let (height, width) = dem_array.dim();

        // Central differences with boundary handling
        let dz_dx = if col == 0 {
            (dem_array[[row, col + 1]] - dem_array[[row, col]]) as f64 / pixel_spacing
        } else if col == width - 1 {
            (dem_array[[row, col]] - dem_array[[row, col - 1]]) as f64 / pixel_spacing
        } else {
            (dem_array[[row, col + 1]] - dem_array[[row, col - 1]]) as f64 / (2.0 * pixel_spacing)
        };

        let dz_dy = if row == 0 {
            (dem_array[[row + 1, col]] - dem_array[[row, col]]) as f64 / pixel_spacing
        } else if row == height - 1 {
            (dem_array[[row, col]] - dem_array[[row - 1, col]]) as f64 / pixel_spacing
        } else {
            (dem_array[[row + 1, col]] - dem_array[[row - 1, col]]) as f64 / (2.0 * pixel_spacing)
        };

        (dz_dx, dz_dy)
    }

    /// Detect shadow and layover conditions for RTC masking
    ///
    /// Based on expert recommendations for proper RTC quality control.
    ///
    /// Shadow occurs when the radar look vector points away from the surface normal,
    /// meaning the surface is not illuminated by the radar.
    ///
    /// Layover occurs when the terrain slope IN THE RANGE DIRECTION is steeper than
    /// the radar incidence angle, causing multiple ground positions to map to the same range.
    /// For right-looking SAR: layover when dz/dx > cot(θ_i) (terrain facing sensor).
    ///
    /// # Arguments
    /// * `surface_normal` - Unit surface normal vector (ENU: x=East, y=North, z=Up)
    /// * `look_vector` - Unit radar look vector (sensor to target in ENU)
    ///
    /// # Returns
    /// * Tuple of (is_shadow, is_layover, cos_lia)
    pub fn detect_shadow_layover(
        &self,
        surface_normal: &SurfaceNormal,
        look_vector: &Vector3,
    ) -> (bool, bool, f64) {
        // Compute local incidence angle cosine
        let cos_lia = self.compute_local_incidence_angle(surface_normal, look_vector);

        // Shadow: dot(l, n) < 0 ⇒ sensor behind terrain, no direct illumination
        let dot_product = look_vector.x * surface_normal.x
            + look_vector.y * surface_normal.y
            + look_vector.z * surface_normal.z;
        let is_shadow = dot_product < 0.0;

        // SCIENTIFIC FIX: Proper layover detection using range direction slope
        // Layover occurs when terrain slope in range direction exceeds incidence angle.
        // For unit look vector l = (lx, ly, lz), the range direction on ground is:
        //   range_dir = normalize(lx, ly, 0) (horizontal projection of look vector)
        // Terrain slope in range direction: slope_range = -(nx/nz, ny/nz) · range_dir
        // Layover if: slope_range > cot(incidence_angle) = lz / sqrt(lx² + ly²)
        //
        // For slopes facing the sensor (positive slope_range), layover condition is:
        //   slope_range > cot(θ_i)
        // which simplifies to: cos_lia < 0 OR (n̂ · l̂_horizontal) > cos(θ_i)

        let horizontal_mag = (look_vector.x * look_vector.x + look_vector.y * look_vector.y).sqrt();

        let is_layover = if horizontal_mag > 1e-6 && surface_normal.z.abs() > 1e-6 {
            // Range direction unit vector (horizontal)
            let range_x = look_vector.x / horizontal_mag;
            let range_y = look_vector.y / horizontal_mag;

            // Terrain slope in range direction (positive = facing sensor)
            // slope = -(nx/nz * range_x + ny/nz * range_y)
            let slope_range = -(surface_normal.x / surface_normal.z * range_x
                + surface_normal.y / surface_normal.z * range_y);

            // Cotangent of incidence angle = vertical / horizontal component
            let cot_incidence = look_vector.z.abs() / horizontal_mag;

            // Layover: terrain slopes toward sensor steeper than incidence angle
            // Also flag very steep slopes (cos_lia < 0.05 ≈ 87°) as potential layover
            slope_range > cot_incidence || cos_lia < 0.05
        } else {
            // Fallback for edge cases: flag very steep angles
            cos_lia < 0.05
        };

        (is_shadow, is_layover, cos_lia)
    }

    /// Compute layover/shadow mask using SNAP-style elevation angle traversal
    ///
    /// This implements the algorithm used by SNAP's RangeDopplerGeocodingOp to detect
    /// layover and shadow regions by traversing elevation angles along range direction.
    ///
    /// # Algorithm (based on SNAP RangeDopplerGeocodingOp.java)
    ///
    /// For each azimuth line, traverse from near range to far range:
    /// 1. Compute elevation angle at each range position
    /// 2. Track maximum elevation angle seen so far (from near range)
    /// 3. **Layover**: If current elevation angle > previous maximum, mark as layover
    ///    (Multiple terrain points map to same range bin)
    /// 4. **Shadow**: If current elevation angle < global minimum along the line
    ///    (Terrain is not illuminated by radar)
    ///
    /// # Arguments
    /// * `dem_array` - DEM elevation data in geocoded coordinates
    /// * `sensor_position` - Sensor position in ECEF coordinates
    /// * `near_range_m` - Near range distance in meters
    /// * `range_spacing_m` - Range pixel spacing in meters
    /// * `azimuth_line` - Which azimuth line to process
    ///
    /// # Returns
    /// * Tuple of (layover_mask, shadow_mask) as boolean vectors along range
    ///
    /// # Reference
    /// SNAP S1TBX: RangeDopplerGeocodingOp.java, computeLayoverShadow() method
    pub fn compute_layover_shadow_traversal(
        &self,
        dem_array: &Array2<f32>,
        sensor_position: &super::types::Position3D,
        near_range_m: f64,
        range_spacing_m: f64,
        azimuth_line: usize,
    ) -> (Vec<bool>, Vec<bool>) {
        let (height, width) = dem_array.dim();

        if azimuth_line >= height {
            return (vec![false; width], vec![false; width]);
        }

        let mut layover_mask = vec![false; width];
        let mut shadow_mask = vec![false; width];

        // Phase 1: Near-to-far traversal for layover detection
        // Layover occurs when a closer terrain point has a higher elevation angle
        // than the current point (foreshortening becomes extreme)
        let mut max_elev_angle = f64::NEG_INFINITY;

        for col in 0..width {
            let elevation = dem_array[[azimuth_line, col]] as f64;

            // Compute ground position (approximate - use DEM pixel center)
            let (lon, lat) = self.pixel_to_geo(azimuth_line, col);
            let target_ecef = self.latlon_to_ecef(lat, lon, elevation);

            // Compute slant range and elevation angle
            let dx = target_ecef[0] - sensor_position.x;
            let dy = target_ecef[1] - sensor_position.y;
            let dz = target_ecef[2] - sensor_position.z;
            let slant_range = (dx * dx + dy * dy + dz * dz).sqrt();

            // Elevation angle = angle from horizontal plane to target
            // Using simplified formula: sin(elev) ≈ (h_target - h_sensor) / slant_range
            // More accurate: use proper geometry with Earth center
            let horizontal_dist = (dx * dx + dy * dy).sqrt();
            let elev_angle = dz.atan2(horizontal_dist); // radians

            // Layover detection: if elevation angle exceeds previous maximum
            // (terrain rising faster than range increases)
            if elev_angle > max_elev_angle + 1e-6 {
                // This could be layover - mark if slope is significant
                if col > 0 && elev_angle - max_elev_angle > 0.01 {
                    // ~0.6 degrees
                    layover_mask[col] = true;
                }
                max_elev_angle = elev_angle;
            }
        }

        // Phase 2: Far-to-near traversal for shadow detection
        // Shadow occurs when terrain blocks radar view from sensor
        let mut min_elev_angle = f64::INFINITY;

        for col in (0..width).rev() {
            let elevation = dem_array[[azimuth_line, col]] as f64;

            let (lon, lat) = self.pixel_to_geo(azimuth_line, col);
            let target_ecef = self.latlon_to_ecef(lat, lon, elevation);

            let dx = target_ecef[0] - sensor_position.x;
            let dy = target_ecef[1] - sensor_position.y;
            let dz = target_ecef[2] - sensor_position.z;

            let horizontal_dist = (dx * dx + dy * dy).sqrt();
            let elev_angle = dz.atan2(horizontal_dist);

            // Shadow detection: if elevation angle is below the minimum seen from far range
            // (terrain is hidden behind closer terrain)
            if elev_angle < min_elev_angle - 1e-6 {
                if col < width - 1 && min_elev_angle - elev_angle > 0.01 {
                    shadow_mask[col] = true;
                }
            }
            min_elev_angle = min_elev_angle.min(elev_angle);
        }

        (layover_mask, shadow_mask)
    }

    /// Compute full layover/shadow mask for entire DEM using elevation angle traversal
    ///
    /// This is the full SNAP-compatible implementation that processes all azimuth lines.
    ///
    /// # Arguments
    /// * `dem_array` - DEM elevation data
    /// * `sensor_positions` - Sensor positions per azimuth line (or interpolated)
    ///
    /// # Returns
    /// * Tuple of (layover_mask, shadow_mask) as 2D boolean arrays
    pub fn compute_full_layover_shadow_mask(
        &self,
        dem_array: &Array2<f32>,
        sensor_position: &super::types::Position3D, // Simplified: assume constant sensor pos
        near_range_m: f64,
        range_spacing_m: f64,
    ) -> (Array2<bool>, Array2<bool>) {
        let (height, width) = dem_array.dim();

        let mut layover_mask = Array2::from_elem((height, width), false);
        let mut shadow_mask = Array2::from_elem((height, width), false);

        // Process each azimuth line in parallel
        use rayon::prelude::*;

        let results: Vec<(usize, Vec<bool>, Vec<bool>)> = (0..height)
            .into_par_iter()
            .map(|row| {
                let (lo, sh) = self.compute_layover_shadow_traversal(
                    dem_array,
                    sensor_position,
                    near_range_m,
                    range_spacing_m,
                    row,
                );
                (row, lo, sh)
            })
            .collect();

        for (row, lo, sh) in results {
            for col in 0..width {
                layover_mask[[row, col]] = lo[col];
                shadow_mask[[row, col]] = sh[col];
            }
        }

        // Log statistics
        let layover_count = layover_mask.iter().filter(|&&v| v).count();
        let shadow_count = shadow_mask.iter().filter(|&&v| v).count();
        let total = height * width;

        log::info!(
            "🏔️  Layover/shadow mask computed: {:.2}% layover, {:.2}% shadow",
            layover_count as f64 / total as f64 * 100.0,
            shadow_count as f64 / total as f64 * 100.0
        );

        (layover_mask, shadow_mask)
    }

    /// Helper to convert pixel coordinates to geographic coordinates
    fn pixel_to_geo(&self, row: usize, col: usize) -> (f64, f64) {
        let lon =
            self.dem_transform.top_left_x + (col as f64 + 0.5) * self.dem_transform.pixel_width;
        let lat =
            self.dem_transform.top_left_y + (row as f64 + 0.5) * self.dem_transform.pixel_height;
        (lon, lat)
    }

    /// Compute local incidence angle cosine
    ///
    /// The local incidence angle (LIA) is the angle between the radar look vector
    /// and the surface normal. It accounts for terrain slope effects on backscatter.
    ///
    /// # Arguments
    /// * `surface_normal` - Unit surface normal vector
    /// * `radar_look_vector` - Unit radar look vector (sensor to target)
    ///
    /// # Returns
    /// * Cosine of local incidence angle [0, 1]
    pub fn compute_local_incidence_angle(
        &self,
        surface_normal: &SurfaceNormal,
        radar_look_vector: &Vector3,
    ) -> f64 {
        // Local incidence angle is angle between radar look vector and surface normal
        // cos(θ_lia) = |look_vector · surface_normal|
        let dot_product = radar_look_vector.x * surface_normal.x
            + radar_look_vector.y * surface_normal.y
            + radar_look_vector.z * surface_normal.z;
        dot_product.abs()
    }

    /// Compute RTC for a single pixel during geocoding
    ///
    /// This is the unified RTC computation function that should be called from the
    /// geocoding loop. It computes all RTC-related quantities from the geometric
    /// information obtained during Range-Doppler transformation.
    ///
    /// # Scientific Implementation (Small 2011)
    ///
    /// For **AreaProjection** mode (default):
    /// ```text
    /// γ⁰ = σ⁰ × (A_flat / A_slope)
    ///    = σ⁰ × cos(θ_ref) / (cos(θ_local) × M)
    /// ```
    /// where:
    /// - θ_ref: reference ellipsoid incidence angle (from metadata)
    /// - θ_local: local incidence angle (from surface normal and look vector)
    /// - M: slope magnitude = √(1 + (∂z/∂x)² + (∂z/∂y)²)
    ///
    /// For **CosineLocalIncidenceAngle** mode:
    /// ```text
    /// γ⁰ = σ⁰ × cos(θ_ref) / cos(θ_local)
    /// ```
    ///
    /// # Arguments
    /// * `dem_array` - DEM elevation data
    /// * `dem_row` - DEM pixel row
    /// * `dem_col` - DEM pixel column
    /// * `look_vector` - Unit look vector from target to satellite (ECEF)
    /// * `cos_ref` - Cosine of reference incidence angle
    /// * `mode` - RTC mode to use
    ///
    /// # Returns
    /// * `RtcResult` containing scale factor, LIA, and quality flags
    pub fn compute_rtc_for_pixel(
        &self,
        dem_array: &Array2<f32>,
        dem_row: usize,
        dem_col: usize,
        look_vector: &Vector3,
        cos_ref: f64,
        mode: RtcMode,
    ) -> RtcResult {
        // 1) Compute surface normal from DEM gradients
        let surface_normal = self.compute_surface_normal(dem_array, dem_row, dem_col);

        // 2) Detect shadow and layover conditions
        let (is_shadow, is_layover, cos_lia) =
            self.detect_shadow_layover(&surface_normal, look_vector);

        // 3) Compute local incidence angle in degrees
        let lia_deg = cos_lia.clamp(0.0, 1.0).acos().to_degrees() as f32;

        // 4) Handle shadow/layover cases - return NaN scale to indicate invalid pixel
        if is_shadow {
            return RtcResult {
                scale: f32::NAN,
                local_incidence_angle_deg: lia_deg,
                quality_flag: RtcQualityFlag::Shadow as u8,
                is_shadow: true,
                is_layover: false,
            };
        }

        // 5) Compute RTC scale based on mode
        let (scale, quality_flag) = match mode {
            RtcMode::AreaProjection => {
                // Get DEM pixel spacing for gradient computation
                let (height, width) = dem_array.dim();

                // Compute DEM gradients
                let pixel_spacing = self.get_dem_pixel_spacing_meters(dem_row);
                let (dz_dx, dz_dy) =
                    Self::compute_dem_gradients(dem_array, dem_row, dem_col, pixel_spacing);

                // Compute area projection scale with quality
                let (scale, flag, _area_ratio) =
                    Self::rtc_area_projection_scale_with_quality(cos_lia, dz_dx, dz_dy, cos_ref);

                (scale, flag)
            }
            RtcMode::CosineLocalIncidenceAngle => {
                // Simple cosine correction
                Self::rtc_gamma0_scale_with_quality(cos_lia, cos_ref)
            }
            RtcMode::SnapSigmaSin => {
                // SNAP-compatible sigma-sin correction
                // σ⁰_TC = σ⁰ × sin(θ_DEM) / sin(θ_ellipsoid)
                Self::rtc_snap_sigma_sin_scale_with_quality(cos_lia, cos_ref)
            }
        };

        // 6) Combine quality flags
        let mut final_flag = quality_flag;
        if is_layover {
            final_flag |= RtcQualityFlag::Layover as u8;
        }

        RtcResult {
            scale,
            local_incidence_angle_deg: lia_deg,
            quality_flag: final_flag,
            is_shadow,
            is_layover,
        }
    }

    /// Get DEM pixel spacing in meters at a given row
    ///
    /// Handles both geographic (degrees) and projected (meters) DEMs.
    fn get_dem_pixel_spacing_meters(&self, row: usize) -> f64 {
        if self.dem_crs == 4326 {
            // Geographic DEM (degrees) - convert to meters at pixel latitude
            let pixel_lat =
                self.dem_transform.top_left_y + (row as f64) * self.dem_transform.pixel_height;
            let lat_rad = pixel_lat.to_radians();
            let cos_lat = lat_rad.cos();

            // Average of x and y spacing at this latitude
            let n = WGS84_SEMI_MAJOR_AXIS_M
                / (1.0
                    - crate::constants::geodetic::WGS84_ECCENTRICITY_SQUARED
                        * lat_rad.sin().powi(2))
                .sqrt();
            let meters_per_deg = n * cos_lat * std::f64::consts::PI / 180.0;

            self.dem_transform.pixel_width.abs() * meters_per_deg
        } else {
            // Projected DEM (already in meters)
            self.dem_transform.pixel_width.abs()
        }
    }
}

/// Compute look vector from target ECEF to satellite position
///
/// This is a standalone function for use when the TerrainCorrector is not available.
///
/// # Arguments
/// * `target_ecef` - Target position in ECEF coordinates [x, y, z]
/// * `sat_pos` - Satellite position in ECEF coordinates
///
/// # Returns
/// * Unit look vector from target to satellite
pub fn compute_look_vector(target_ecef: &[f64; 3], sat_pos: &super::types::Position3D) -> Vector3 {
    let dx = sat_pos.x - target_ecef[0];
    let dy = sat_pos.y - target_ecef[1];
    let dz = sat_pos.z - target_ecef[2];
    let norm = (dx * dx + dy * dy + dz * dz).sqrt();

    if norm <= 0.0 {
        // Degenerate case - return vertical
        return Vector3 {
            x: 0.0,
            y: 0.0,
            z: 1.0,
        };
    }

    Vector3 {
        x: dx / norm,
        y: dy / norm,
        z: dz / norm,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::types::GeoTransform;

    /// Test surface normal computation on flat terrain
    #[test]
    fn test_surface_normal_flat() {
        // Create a flat DEM
        let dem = Array2::from_elem((10, 10), 100.0f32);
        let transform = GeoTransform {
            top_left_x: 0.0,
            top_left_y: 0.0,
            pixel_width: 10.0,
            pixel_height: -10.0,
            rotation_x: 0.0,
            rotation_y: 0.0,
        };

        let corrector = TerrainCorrector::new(dem.clone(), transform, -9999.0, 32632, 32632, 10.0);

        let normal = corrector.compute_surface_normal(&dem, 5, 5);

        // Flat terrain should have vertical normal (0, 0, 1)
        assert!(normal.z > 0.99, "Flat terrain should have upward normal");
        assert!(normal.x.abs() < 0.01, "Flat terrain should have no x-tilt");
        assert!(normal.y.abs() < 0.01, "Flat terrain should have no y-tilt");
    }

    /// Test local incidence angle computation
    #[test]
    fn test_local_incidence_angle() {
        let dem = Array2::from_elem((10, 10), 100.0f32);
        let transform = GeoTransform {
            top_left_x: 0.0,
            top_left_y: 0.0,
            pixel_width: 10.0,
            pixel_height: -10.0,
            rotation_x: 0.0,
            rotation_y: 0.0,
        };

        let corrector = TerrainCorrector::new(dem, transform, -9999.0, 32632, 32632, 10.0);

        // Vertical surface normal
        let normal = SurfaceNormal {
            x: 0.0,
            y: 0.0,
            z: 1.0,
        };

        // 45-degree look vector
        let look = Vector3 {
            x: 0.7071,
            y: 0.0,
            z: 0.7071,
        };

        let cos_lia = corrector.compute_local_incidence_angle(&normal, &look);
        assert!(
            (cos_lia - 0.7071).abs() < 0.01,
            "45-degree look should give cos(45°) ≈ 0.707"
        );
    }

    /// Test RTC gamma0 scale factor
    #[test]
    fn test_rtc_gamma0_scale() {
        // Equal angles should give scale = 1.0
        let scale = TerrainCorrector::rtc_gamma0_scale(0.5, 0.5);
        assert!(
            (scale - 1.0).abs() < 0.01,
            "Equal cos values should give scale ≈ 1.0"
        );

        // Steeper local angle (smaller cos_lia) should increase scale
        let scale = TerrainCorrector::rtc_gamma0_scale(0.25, 0.5);
        assert!(scale > 1.0, "Steeper angle should increase scale");

        // Scale should be clamped to new wider bounds (default: 0.01 to 100.0)
        // Old test checked <= 10.0, new bounds allow up to 100.0
        let scale = TerrainCorrector::rtc_gamma0_scale(0.001, 0.5);
        assert!(
            scale <= 100.0,
            "Scale should be clamped to max (default 100.0)"
        );
        assert!(
            scale >= 0.01,
            "Scale should be clamped to min (default 0.01)"
        );

        // Test the quality-flagged version
        let (scale, flag) = TerrainCorrector::rtc_gamma0_scale_with_quality(0.0001, 0.5);
        assert!(scale <= 100.0, "Scale should be clamped");
        assert_eq!(flag, 2, "High clamping should set flag to 2");

        let (scale, flag) = TerrainCorrector::rtc_gamma0_scale_with_quality(0.5, 0.5);
        assert_eq!(flag, 0, "No clamping should set flag to 0");
    }

    /// Test incidence angle calculation against known values
    ///
    /// Validates that local incidence angle calculation is correct:
    /// cos(θ_lia) = n⃗ · l⃗
    /// Where n⃗ is surface normal and l⃗ is look vector
    #[test]
    fn test_incidence_angle_calculation_accuracy() {
        let dem = Array2::from_elem((10, 10), 100.0f32);
        let transform = GeoTransform {
            top_left_x: 0.0,
            top_left_y: 0.0,
            pixel_width: 10.0,
            pixel_height: -10.0,
            rotation_x: 0.0,
            rotation_y: 0.0,
        };

        let corrector = TerrainCorrector::new(dem, transform, -9999.0, 32632, 32632, 10.0);

        // Test case 1: Normal perpendicular to look vector (90° incidence)
        let normal_90 = SurfaceNormal {
            x: 0.0,
            y: 0.0,
            z: 1.0, // Vertical normal
        };
        let look_90 = Vector3 {
            x: 1.0,
            y: 0.0,
            z: 0.0, // Horizontal look
        };
        let cos_lia_90 = corrector.compute_local_incidence_angle(&normal_90, &look_90);
        assert!(
            cos_lia_90.abs() < 0.01,
            "90° incidence should give cos ≈ 0: got {:.6}",
            cos_lia_90
        );

        // Test case 2: Normal parallel to look vector (0° incidence)
        let normal_0 = SurfaceNormal {
            x: 0.0,
            y: 0.0,
            z: 1.0,
        };
        let look_0 = Vector3 {
            x: 0.0,
            y: 0.0,
            z: 1.0, // Same direction as normal
        };
        let cos_lia_0 = corrector.compute_local_incidence_angle(&normal_0, &look_0);
        assert!(
            (cos_lia_0 - 1.0).abs() < 0.01,
            "0° incidence should give cos ≈ 1: got {:.6}",
            cos_lia_0
        );

        // Test case 3: 45° incidence
        let normal_45 = SurfaceNormal {
            x: 0.0,
            y: 0.0,
            z: 1.0,
        };
        let look_45 = Vector3 {
            x: 0.7071067811865476, // cos(45°)
            y: 0.0,
            z: 0.7071067811865476, // sin(45°)
        };
        let cos_lia_45 = corrector.compute_local_incidence_angle(&normal_45, &look_45);
        assert!(
            (cos_lia_45 - 0.7071067811865476).abs() < 0.01,
            "45° incidence should give cos ≈ 0.707: got {:.6}",
            cos_lia_45
        );

        // Test case 4: Tilted surface (30° slope)
        let slope_angle = 30.0_f64.to_radians();
        let normal_tilted = SurfaceNormal {
            x: slope_angle.sin(),
            y: 0.0,
            z: slope_angle.cos(),
        };
        let look_tilted = Vector3 {
            x: 0.0,
            y: 0.0,
            z: 1.0, // Vertical look
        };
        let cos_lia_tilted = corrector.compute_local_incidence_angle(&normal_tilted, &look_tilted);
        let expected_cos = slope_angle.cos(); // cos(30°) ≈ 0.866
        assert!(
            (cos_lia_tilted - expected_cos).abs() < 0.01,
            "Tilted surface (30°) should give cos ≈ {:.6}: got {:.6}",
            expected_cos,
            cos_lia_tilted
        );
    }
}
