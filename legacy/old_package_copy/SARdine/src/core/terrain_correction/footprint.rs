//! SAR image footprint computation
//!
//! This module provides functions to compute the geographic bounding box
//! (footprint) of a SAR image by geocoding its corners using Range-Doppler
//! transformation. This is critical for loading the correct DEM coverage
//! before terrain correction.

use crate::core::terrain_correction::range_doppler::RangeDopplerParams;
use crate::core::terrain_correction::types::{Position3D, Velocity3D};
use crate::io::orbit::OrbitReader;
use crate::types::{BoundingBox, OrbitData, SarError, SarResult};
use chrono::DateTime;

/// Epsilon for near-zero comparisons in velocity/vector magnitude calculations
/// Used to detect degenerate cases (parallel vectors, zero magnitude)
const EPSILON: f64 = 1e-6;

// Vector magnitude helper functions (optimized for 3D vectors)
#[inline]
fn vec3_norm(x: f64, y: f64, z: f64) -> f64 {
    (x * x + y * y + z * z).sqrt()
}

/// Convert ECEF coordinates to WGS84 latitude/longitude
///
/// Returns (latitude, longitude) in degrees, matching standard geographic coordinate order.
fn ecef_to_wgs84(x: f64, y: f64, z: f64) -> (f64, f64) {
    let a = crate::constants::geodetic::WGS84_SEMI_MAJOR_AXIS_M;
    let f = 1.0 / 298.257223563; // WGS84 flattening
    let e_sq = 2.0 * f - f * f; // First eccentricity squared

    let p = (x * x + y * y).sqrt();
    let theta = (z * a / (p * (1.0_f64 - e_sq).sqrt())).atan();

    // OPTIMIZED: Use sin_cos() for better precision and performance
    let (sin_theta, cos_theta) = theta.sin_cos();

    let lat = (z + e_sq * a / (1.0 - e_sq).sqrt() * sin_theta * sin_theta * sin_theta)
        / (p - e_sq * a * cos_theta * cos_theta * cos_theta);
    let lat = lat.atan();

    let lon = y.atan2(x);

    // Return (lat, lon) to match standard geographic coordinate order
    (lat.to_degrees(), lon.to_degrees())
}

/// Standalone orbit interpolation (doesn't need TerrainCorrector instance)
///
/// This allows computing bbox before DEM loading, which is critical for
/// ensuring the DEM is loaded for the correct geographic region.
pub(crate) fn standalone_orbit_interpolation(
    orbit_data: &OrbitData,
    time_seconds: f64,
) -> SarResult<(Position3D, Velocity3D)> {
    // Check for empty orbit vectors before accessing
    let first_sv = orbit_data.state_vectors.first().ok_or_else(|| {
        SarError::DataProcessingError(
            "Empty orbit state vectors in standalone_orbit_interpolation".to_string(),
        )
    })?;
    let last_sv = orbit_data.state_vectors.last().ok_or_else(|| {
        SarError::DataProcessingError(
            "Empty orbit state vectors in standalone_orbit_interpolation".to_string(),
        )
    })?;

    // Check temporal coverage with 5-second margin
    let first = crate::types::datetime_to_utc_seconds(first_sv.time);
    let last = crate::types::datetime_to_utc_seconds(last_sv.time);

    // CRITICAL DEBUG: Log orbit coverage and target time
    log::debug!("🛰️  [ORBIT INTERPOLATION DEBUG]");
    log::debug!(
        "   Target time: {:.6}s ({})",
        time_seconds,
        chrono::DateTime::from_timestamp(
            time_seconds as i64,
            ((time_seconds.fract()) * 1e9) as u32
        )
        .map(|dt| dt.to_rfc3339())
        .unwrap_or_else(|| "INVALID".to_string())
    );
    log::debug!(
        "   Orbit coverage: [{:.6}, {:.6}]s ({} state vectors)",
        first,
        last,
        orbit_data.state_vectors.len()
    );
    log::debug!("   Reference epoch: {:?}", orbit_data.reference_time);

    // Log first few state vectors for debugging
    for (i, sv) in orbit_data.state_vectors.iter().take(3).enumerate() {
        let sv_time = crate::types::datetime_to_utc_seconds(sv.time);
        log::debug!(
            "   SV[{}]: time={:.6}s pos=({:.1}, {:.1}, {:.1})",
            i,
            sv_time,
            sv.position[0],
            sv.position[1],
            sv.position[2]
        );
    }
    if orbit_data.state_vectors.len() > 3 {
        log::debug!(
            "   ... ({} more state vectors)",
            orbit_data.state_vectors.len() - 3
        );
    }

    if time_seconds < first - 5.0 || time_seconds > last + 5.0 {
        log::error!(
            "❌ Interpolation time {:.3}s OUTSIDE orbit coverage [{:.3}, {:.3}]!",
            time_seconds,
            first,
            last
        );
        return Err(SarError::DataProcessingError(format!(
            "Interpolation time {:.3}s outside orbit coverage [{:.3}, {:.3}]",
            time_seconds, first, last
        )));
    }

    // Convert time_seconds to DateTime<Utc> using modern chrono API
    let time = DateTime::from_timestamp(time_seconds as i64, (time_seconds.fract() * 1e9) as u32)
        .ok_or_else(|| {
        SarError::DataProcessingError(format!("Invalid timestamp: {}", time_seconds))
    })?;

    // STRICT REQUIREMENT: Must have ≥4 state vectors for cubic spline interpolation
    // No linear fallbacks - fail fast for scientific compliance
    if orbit_data.state_vectors.len() >= 4 {
        let position = OrbitReader::interpolate_position(orbit_data, time).map_err(|e| {
            SarError::DataProcessingError(format!(
                "CUBIC SPLINE ONLY: Orbit position interpolation failed at time {}: {}",
                time, e
            ))
        })?;

        let velocity = OrbitReader::interpolate_velocity(orbit_data, time).map_err(|e| {
            SarError::DataProcessingError(format!(
                "CUBIC SPLINE ONLY: Orbit velocity interpolation failed at time {}: {}",
                time, e
            ))
        })?;

        // DEBUG: Log interpolated position
        log::debug!(
            "   Interpolated position: ({:.1}, {:.1}, {:.1}) m",
            position[0],
            position[1],
            position[2]
        );

        // Validate finiteness of interpolated values
        if !position[0].is_finite()
            || !position[1].is_finite()
            || !position[2].is_finite()
            || !velocity[0].is_finite()
            || !velocity[1].is_finite()
            || !velocity[2].is_finite()
        {
            return Err(SarError::DataProcessingError(format!(
                "Non-finite orbit interpolation result at time {}",
                time
            )));
        }

        Ok((
            Position3D {
                x: position[0],
                y: position[1],
                z: position[2],
            },
            Velocity3D {
                x: velocity[0],
                y: velocity[1],
                z: velocity[2],
            },
        ))
    } else {
        Err(SarError::DataProcessingError(format!(
            "Insufficient orbit state vectors ({} < 4) for cubic spline interpolation",
            orbit_data.state_vectors.len()
        )))
    }
}

/// Compute bounding box from SAR image corners
///
/// This function geocodes the four corners of a SAR image using Range-Doppler
/// transformation to determine the actual geographic footprint. This is critical
/// for loading the correct DEM coverage before terrain correction.
///
/// # Arguments
/// * `sar_height` - Height of the SAR image in pixels (multilooked)
/// * `sar_width` - Width of the SAR image in pixels (multilooked)
/// * `params` - Range-Doppler parameters for the SAR image
/// * `orbit_data` - Precise orbit data for satellite position/velocity
///
/// # Returns
/// Geographic bounding box covering the SAR image footprint
pub(crate) fn compute_bbox_from_sar_image(
    sar_height: usize,
    sar_width: usize,
    params: &RangeDopplerParams,
    orbit_data: &OrbitData,
) -> SarResult<BoundingBox> {
    log::debug!(
        "🔍 [BBOX DIAG] Computing bounding box from SAR image corners ({}x{})",
        sar_width, sar_height
    );
    log::debug!(
        "   Orbit reference epoch: {:.6}s UTC",
        params.orbit_ref_epoch_utc
    );
    log::debug!(
        "   Product start (rel): {:.6}s, duration: {:.3}s",
        params.product_start_rel_s, params.product_duration
    );
    log::debug!(
        "   Burst segments: {}, Subswaths: {}",
        params.burst_segments.len(),
        params.subswaths.len()
    );
    log::debug!(
        "🔍 Computing bounding box from SAR image corners ({}x{})",
        sar_width,
        sar_height
    );
    log::debug!(
        "   Orbit reference epoch: {:.6}s UTC",
        params.orbit_ref_epoch_utc
    );
    log::debug!(
        "   Product start (rel): {:.6}s, duration: {:.3}s",
        params.product_start_rel_s,
        params.product_duration
    );
    log::debug!(
        "   Burst segments: {}, Subswaths: {}",
        params.burst_segments.len(),
        params.subswaths.len()
    );

    // Define four corners: (row, col) in multilooked coordinates
    let corners = [
        (0, 0),                          // Top-left
        (0, sar_width - 1),              // Top-right
        (sar_height - 1, 0),             // Bottom-left
        (sar_height - 1, sar_width - 1), // Bottom-right
    ];

    let mut min_lat = f64::INFINITY;
    let mut max_lat = f64::NEG_INFINITY;
    let mut min_lon = f64::INFINITY;
    let mut max_lon = f64::NEG_INFINITY;
    let mut successful_corners = 0;

    // Convert multilooked coordinates to native (use pre-computed safe values)
    let range_ml = params.range_multilook_safe;
    let azimuth_ml = params.azimuth_multilook_safe;

    for (corner_idx, (row_ml, col_ml)) in corners.iter().enumerate() {
        // Convert to native SAR coordinates
        let row_native = (*row_ml as f64) * azimuth_ml;
        let col_native = (*col_ml as f64) * range_ml;

        log::debug!(
            "🔍 [BBOX DIAG] Corner {}: multilooked=({}, {}), native=({:.1}, {:.1})",
            corner_idx, row_ml, col_ml, row_native, col_native
        );
        log::debug!(
            "🔍 [FOOTPRINT DEBUG] Corner {}: multilooked=({}, {}), native=({:.1}, {:.1})",
            corner_idx,
            row_ml,
            col_ml,
            row_native,
            col_native
        );

        // Convert native azimuth pixel to azimuth time
        let azimuth_time_rel = if !params.burst_segments.is_empty() {
            log::debug!(
                "   Using burst segments ({} available)",
                params.burst_segments.len()
            );
            // Use burst segments for accurate timing
            let line_ml = *row_ml as f64;
            if let Some(segment) = params
                .burst_segments
                .iter()
                .find(|seg| line_ml >= seg.start_line && line_ml < seg.end_line)
            {
                let line_offset = line_ml - segment.start_line;
                let time_offset = line_offset * segment.line_time_interval;
                let azimuth_time_rel = segment.start_time_rel + time_offset;
                log::debug!("   Found burst segment: subswath={}, burst={}, start_line={:.1}, end_line={:.1}",
                           segment.subswath_id, segment.burst_index, segment.start_line, segment.end_line);
                log::debug!("   Segment timing: start_time_rel={:.9}s, end_time_rel={:.9}s, line_offset={:.1}, time_offset={:.9}s",
                           segment.start_time_rel, segment.end_time_rel, line_offset, time_offset);
                log::debug!("   ⚠️  TIMING MISMATCH CHECK: segment.start_time_rel={:.9}s vs product_start_rel_s={:.9}s (diff={:.9}s)",
                           segment.start_time_rel, params.product_start_rel_s,
                           segment.start_time_rel - params.product_start_rel_s);
                azimuth_time_rel
            } else {
                // Fallback: use product start + linear mapping
                log::warn!("   ⚠️  Corner {}: No burst segment found for line {:.1}, using linear fallback", corner_idx, row_ml);
                params.product_start_rel_s + (row_native * params.azimuth_time_interval)
            }
        } else {
            // No burst segments: linear mapping
            log::warn!("   ⚠️  Corner {}: No burst segments available, using linear mapping (may be inaccurate for TOPSAR)", corner_idx);
            params.product_start_rel_s + (row_native * params.azimuth_time_interval)
        };

        log::debug!(
            "   Azimuth time (rel): {:.9}s (product_start={:.9}s, interval={:.9}s)",
            azimuth_time_rel,
            params.product_start_rel_s,
            params.azimuth_time_interval
        );

        // Convert native range pixel to slant range
        // CRITICAL: Use per-subswath slant_range_time for merged IW data
        // For merged IW (IW1+IW2+IW3), each subswath has its own slant_range_time.
        // Using the global slant_range_time (from IW1) for IW2/IW3 corners causes incorrect geocoding.
        let range_pixel_spacing_time = 2.0 * params.range_pixel_spacing / params.speed_of_light;

        // Find which subswath this native range pixel belongs to
        // CRITICAL FIX: For merged IW, using wrong subswath's slant_range_time causes incorrect geocoding
        // Try exact match first, then try nearby pixels (handle rounding errors), then find closest subswath

        // EXPLICIT VALIDATION: For right edge corners (high col_native), we expect IW3
        // Log all subswath ranges to help debug matching issues
        if col_native > 40000.0 {
            log::debug!("🔍 [BBOX VALIDATION] Corner {}: col_native={:.1}, checking subswath match. Available subswaths:", corner_idx, col_native);
            for (id, sw) in &params.subswaths {
                log::error!(
                    "   {}: samples {}..{} (exclusive), slant_range_time={:.9}s",
                    id,
                    sw.first_sample_global,
                    sw.last_sample_global,
                    sw.slant_range_time
                );
                // Explicit check: should col_native match this subswath?
                let sample = col_native as usize;
                let should_match =
                    sample >= sw.first_sample_global && sample < sw.last_sample_global;
                if should_match {
                    log::error!(
                        "      → col_native={:.1} SHOULD match {} ({} <= {} < {})",
                        col_native,
                        id,
                        sw.first_sample_global,
                        sample,
                        sw.last_sample_global
                    );
                }
            }
        }

        let two_way_time = if let Some(subswath) = params.find_subswath_for_sample(col_native) {
            // VALIDATION: For high col_native, verify we matched IW3
            if col_native > 40000.0 && subswath.id != "IW3" {
                log::debug!("⚠️  [BBOX VALIDATION] Corner {}: col_native={:.1} matched {} instead of IW3! This may cause wrong coordinates!",
                           corner_idx, col_native, subswath.id);
            }
            // Local pixel within this subswath (0-based within the subswath)
            let local_pixel = col_native - (subswath.first_sample_global as f64);
            // Two-way time = subswath's slant_range_time + local_pixel * spacing
            let twt = subswath.slant_range_time + (local_pixel * range_pixel_spacing_time);
            log::debug!("✅ [BBOX DIAG] Corner {}: EXACT MATCH - Subswath {} (samples {}..{}), local_pixel={:.1}, slant_range_time={:.9}s, twt={:.9}s",
                       corner_idx, subswath.id, subswath.first_sample_global, subswath.last_sample_global,
                       local_pixel, subswath.slant_range_time, twt);
            log::debug!("   Subswath found: {} (samples {}..{}), local_pixel={:.1}, slant_range_time={:.9}s",
                       subswath.id, subswath.first_sample_global, subswath.last_sample_global,
                       local_pixel, subswath.slant_range_time);
            twt
        } else if let Some(subswath) = params
            .find_subswath_for_sample((col_native - 1.0).max(0.0))
            .or_else(|| params.find_subswath_for_sample(col_native + 1.0))
        {
            // Try adjacent pixels to handle rounding errors at subswath boundaries
            let local_pixel = col_native - (subswath.first_sample_global as f64);
            let twt = subswath.slant_range_time + (local_pixel * range_pixel_spacing_time);
            log::debug!("⚠️  [BBOX DIAG] Corner {}: ADJACENT PIXEL MATCH - Subswath {} (samples {}..{}), local_pixel={:.1}, slant_range_time={:.9}s, twt={:.9}s",
                       corner_idx, subswath.id, subswath.first_sample_global, subswath.last_sample_global,
                       local_pixel, subswath.slant_range_time, twt);
            log::debug!("   Subswath found (adjacent pixel): {} (samples {}..{}), local_pixel={:.1}, slant_range_time={:.9}s",
                       subswath.id, subswath.first_sample_global, subswath.last_sample_global,
                       local_pixel, subswath.slant_range_time);
            twt
        } else if !params.subswaths.is_empty() {
            // Find closest subswath by distance to col_native
            // This handles edge cases where corners are just outside subswath bounds
            let closest_subswath = params.subswaths.values().min_by(|a, b| {
                let dist_a = if col_native < a.first_sample_global as f64 {
                    a.first_sample_global as f64 - col_native
                } else if col_native > a.last_sample_global as f64 {
                    col_native - a.last_sample_global as f64
                } else {
                    0.0
                };
                let dist_b = if col_native < b.first_sample_global as f64 {
                    b.first_sample_global as f64 - col_native
                } else if col_native > b.last_sample_global as f64 {
                    col_native - b.last_sample_global as f64
                } else {
                    0.0
                };
                dist_a
                    .partial_cmp(&dist_b)
                    .unwrap_or(std::cmp::Ordering::Equal)
            });

            if let Some(subswath) = closest_subswath {
                let local_pixel = col_native - (subswath.first_sample_global as f64);
                let twt = subswath.slant_range_time + (local_pixel * range_pixel_spacing_time);

                // CRITICAL WARNING: Corner pixel is outside all subswath bounds
                // Using closest subswath may introduce geocoding errors at image edges
                let distance_to_boundary = if col_native < subswath.first_sample_global as f64 {
                    subswath.first_sample_global as f64 - col_native
                } else {
                    col_native - subswath.last_sample_global as f64
                };

                log::debug!("⚠️  [BBOX DIAG] Corner {}: CLOSEST SUBSWATH FALLBACK - col_native={:.1} not in any subswath, using closest: {} (samples {}..{}), local_pixel={:.1}, slant_range_time={:.9}s, twt={:.9}s",
                          corner_idx, col_native, subswath.id, subswath.first_sample_global, subswath.last_sample_global, local_pixel, subswath.slant_range_time, twt);
                log::warn!("⚠️  [GEOCODING ACCURACY] Corner {} range pixel {:.1} is {:.1} pixels outside subswath bounds",
                          corner_idx, col_native, distance_to_boundary);
                log::warn!("   Using closest subswath {} (samples {}..{}) may introduce geocoding errors of ~{:.1}m at image edges",
                          subswath.id, subswath.first_sample_global, subswath.last_sample_global,
                          distance_to_boundary * params.range_pixel_spacing);
                log::warn!("   This is expected for corner pixels in merged IW data and should not affect interior pixels");
                twt
            } else {
                // This shouldn't happen (subswaths is not empty)
                log::error!(
                    "❌ Corner {}: No subswaths available for range pixel {:.1}",
                    corner_idx,
                    col_native
                );
                return Err(SarError::Processing(format!(
                    "Failed to geocode corner {}: no subswaths available",
                    corner_idx
                )));
            }
        } else {
            // No subswaths (single subswath mode like SM/EW) - use global slant_range_time
            let twt = params.slant_range_time + (col_native * range_pixel_spacing_time);
            log::debug!("   ⚠️  [BBOX DIAG] Corner {}: NO SUBSWATHS MODE - using global slant_range_time={:.9}s, col_native={:.1}, twt={:.9}s",
                      corner_idx, params.slant_range_time, col_native, twt);
            log::debug!(
                "   Corner {}: No subswaths (single subswath mode), using global slant_range_time",
                corner_idx
            );
            twt
        };

        let slant_range = two_way_time * params.speed_of_light / 2.0;
        log::debug!(
            "   [BBOX DIAG] Corner {}: two_way_time={:.9}s, slant_range={:.1}m",
            corner_idx, two_way_time, slant_range
        );
        log::debug!(
            "   Two-way time: {:.9}s, slant_range: {:.1}m",
            two_way_time,
            two_way_time * params.speed_of_light / 2.0
        );

        // Calculate absolute time for orbit interpolation
        // azimuth_time_rel is relative to orbit_ref_epoch_utc
        // The math: absolute_time = orbit_ref_epoch_utc + azimuth_time_rel
        // We use the equivalent: product_start_time_abs + (azimuth_time_rel - product_start_rel_s)
        // where product_start_rel_s = product_start_time_abs - orbit_ref_epoch_utc
        let absolute_time = params.orbit_ref_epoch_utc + azimuth_time_rel;

        // Validate: Check if absolute_time is within scene acquisition window
        let scene_start_abs = params.orbit_ref_epoch_utc + params.product_start_rel_s;
        let scene_end_abs = scene_start_abs + params.product_duration;
        if absolute_time < scene_start_abs - 60.0 || absolute_time > scene_end_abs + 60.0 {
            log::debug!(
                "   ⚠️  WARNING: Absolute time {:.6}s outside scene window [{:.6}, {:.6}]s",
                absolute_time,
                scene_start_abs,
                scene_end_abs
            );
            log::debug!("      This indicates time reference mismatch!");
        }

        log::debug!(
            "   Absolute time: {:.6}s (orbit_epoch={:.6}s + azimuth_time_rel={:.6}s)",
            absolute_time,
            params.orbit_ref_epoch_utc,
            azimuth_time_rel
        );

        // Interpolate orbit to get satellite position
        let (sat_pos, sat_vel) = match standalone_orbit_interpolation(orbit_data, absolute_time) {
            Ok(pos_vel) => {
                log::debug!("   [BBOX DIAG] Corner {}: Satellite position ECEF({:.1}, {:.1}, {:.1})m, velocity({:.1}, {:.1}, {:.1})m/s",
                           corner_idx, pos_vel.0.x, pos_vel.0.y, pos_vel.0.z, pos_vel.1.x, pos_vel.1.y, pos_vel.1.z);
                log::debug!(
                    "   Satellite position: ECEF({:.1}, {:.1}, {:.1})m",
                    pos_vel.0.x,
                    pos_vel.0.y,
                    pos_vel.0.z
                );
                pos_vel
            }
            Err(e) => {
                log::debug!(
                    "   ❌ [BBOX DIAG] Corner {}: Failed to interpolate orbit: {}",
                    corner_idx, e
                );
                log::debug!(
                    "   ❌ Failed to interpolate orbit for corner {}: {}",
                    corner_idx,
                    e
                );
                continue;
            }
        };

        // Solve for ground point using reverse Range-Doppler transformation
        // Given: satellite position, velocity, and slant range at azimuth time
        // Find: ground point P where:
        //   1. |sat_pos - P| = slant_range
        //   2. (sat_pos - P) · sat_vel = 0 (zero-Doppler condition)
        //
        // The zero-Doppler condition means the range vector is perpendicular to velocity.
        // This defines a plane. The ground point is the intersection of:
        //   - Sphere: |sat_pos - P| = slant_range
        //   - Plane: (sat_pos - P) · sat_vel = 0
        //   - Earth surface: |P| = R_earth + elevation

        let sat_ecef = [sat_pos.x, sat_pos.y, sat_pos.z];
        let sat_vel_ecef = [sat_vel.x, sat_vel.y, sat_vel.z];

        // Normalize velocity vector
        let vel_mag = vec3_norm(sat_vel_ecef[0], sat_vel_ecef[1], sat_vel_ecef[2]);
        if vel_mag < EPSILON {
            log::warn!("⚠️  Zero velocity at corner {}, skipping", corner_idx);
            continue;
        }
        let vel_unit = [
            sat_vel_ecef[0] / vel_mag,
            sat_vel_ecef[1] / vel_mag,
            sat_vel_ecef[2] / vel_mag,
        ];

        // =================================================================
        // ISCE2-STYLE TCN RANGE-DOPPLER SOLUTION FOR GROUND POINT
        // =================================================================
        // Reference: ISCE2 topozero.f90, lines 483-523
        //
        // Given: satellite position S, velocity V, slant range R
        // Find: ground point P using TCN (Track-Cross-Nadir) coordinate system
        //
        // TCN Basis:
        //   N (nhat) = unit vector toward Earth center (nadir)
        //   T (that) = velocity projected perpendicular to N, normalized (along-track)
        //   C (chat) = N × T (cross-track, by right-hand rule points LEFT of track)
        //
        // The range vector from satellite to ground:
        //   delta = gamma*N + alpha*T + beta*C
        // where:
        //   gamma = cos(theta) * R  (nadir component)
        //   alpha = (dopfact - gamma * (N·V_unit)) / (V_unit·T)  (along-track, zero for zero-Doppler)
        //   beta = ILRL * sqrt(R² * sin²(theta) - alpha²)  (cross-track)
        //
        // ILRL = -1 for right-looking (Sentinel-1), +1 for left-looking
        // This directly encodes the look side without ambiguous theta selection!

        // Sentinel-1 is always right-looking
        const ILRL: f64 = -1.0;

        let sat_mag = vec3_norm(sat_ecef[0], sat_ecef[1], sat_ecef[2]);

        // =================================================================
        // ISCE2 TCN BASIS: Nadir is LOCAL ELLIPSOID NORMAL
        // =================================================================
        // ISCE2 Ellipsoid.cpp uses the ellipsoid normal at the sub-satellite point,
        // NOT the geocentric direction (which would be -sat_pos/|sat_pos|).
        //
        // From Ellipsoid.cpp line 140-143:
        //   latlon(pos,llh,XYZ_2_LLH);
        //   n[0] = -cos(llh[0]) * cos(llh[1]);
        //   n[1] = -cos(llh[0]) * sin(llh[1]);
        //   n[2] = -sin(llh[0]);
        //
        // This ensures n·v ≈ 0 for a satellite following a realistic orbit.

        // Get satellite geodetic lat/lon
        let (sat_lat, sat_lon) = ecef_to_wgs84(sat_ecef[0], sat_ecef[1], sat_ecef[2]);
        let sat_lat_rad = sat_lat.to_radians();
        let sat_lon_rad = sat_lon.to_radians();

        // Compute local ellipsoid normal (pointing toward Earth = inward)
        let nhat = [
            -(sat_lat_rad.cos() * sat_lon_rad.cos()),
            -(sat_lat_rad.cos() * sat_lon_rad.sin()),
            -sat_lat_rad.sin(),
        ];

        // =================================================================
        // ISCE2 TCN BASIS CONSTRUCTION (Ellipsoid.cpp lines 138-152)
        // =================================================================
        // CRITICAL: ISCE2 computes TCN as:
        //   chat = normalize(nhat × velocity)  -- cross-track
        //   that = normalize(chat × nhat)      -- along-track
        // This is DIFFERENT from: that = vel_perp, chat = nhat × that

        // C (chat) = normalize(N × V) - cross-track direction
        let n_cross_v = [
            nhat[1] * vel_unit[2] - nhat[2] * vel_unit[1],
            nhat[2] * vel_unit[0] - nhat[0] * vel_unit[2],
            nhat[0] * vel_unit[1] - nhat[1] * vel_unit[0],
        ];
        let n_cross_v_mag = vec3_norm(n_cross_v[0], n_cross_v[1], n_cross_v[2]);

        if n_cross_v_mag < EPSILON {
            log::warn!(
                "⚠️  Velocity parallel to nadir at corner {}, skipping",
                corner_idx
            );
            continue;
        }

        let chat = [
            n_cross_v[0] / n_cross_v_mag,
            n_cross_v[1] / n_cross_v_mag,
            n_cross_v[2] / n_cross_v_mag,
        ];

        // T (that) = normalize(C × N) - along-track direction
        let c_cross_n = [
            chat[1] * nhat[2] - chat[2] * nhat[1],
            chat[2] * nhat[0] - chat[0] * nhat[2],
            chat[0] * nhat[1] - chat[1] * nhat[0],
        ];
        let c_cross_n_mag = vec3_norm(c_cross_n[0], c_cross_n[1], c_cross_n[2]);
        let that = [
            c_cross_n[0] / c_cross_n_mag,
            c_cross_n[1] / c_cross_n_mag,
            c_cross_n[2] / c_cross_n_mag,
        ];

        // =================================================================
        // ISCE2 GEOMETRY: Compute theta (off-nadir angle) from range sphere
        // =================================================================
        // Use local spherical approximation with Earth radius at satellite position
        //
        // From triangle: satellite (height H), center of curvature, ground point (height 0)
        //   aa = H + R_curv (satellite radius from center of curvature)
        //   bb = R_curv + h_ground (ground point radius, h_ground ≈ 0 for first iteration)
        //   costheta = 0.5 * ((aa/rng) + (rng/aa) - (bb/aa)*(bb/rng))
        //
        // Simplified for h_ground = 0:
        //   aa = satellite altitude above Earth center = sat_mag
        //   bb = Earth radius = R_earth (approximate)
        //   We solve the law of cosines for the off-nadir angle

        // Compute WGS84 ellipsoid radius at satellite latitude for accurate geometry
        // This improves accuracy compared to using constant equatorial radius
        let a = crate::constants::geodetic::WGS84_SEMI_MAJOR_AXIS_M;
        let e2 = crate::constants::geodetic::WGS84_ECCENTRICITY_SQUARED;
        let sin_lat = sat_lat_rad.sin();
        // Prime vertical radius of curvature at satellite latitude
        let nu = a / (1.0 - e2 * sin_lat * sin_lat).sqrt();
        // Use prime vertical radius as Earth radius approximation (more accurate than constant)
        let earth_radius = nu;

        let aa = sat_mag; // Distance from Earth center to satellite
        let bb = earth_radius; // Distance from Earth center to ground (ellipsoid-aware)
        let r = slant_range;

        // Law of cosines: r² = aa² + bb² - 2*aa*bb*cos(angle_at_center)
        // We want theta = off-nadir angle at satellite
        // Using: aa² = bb² + r² - 2*bb*r*cos(theta_at_ground)
        // And:   cos(theta) = (aa² + r² - bb²) / (2*aa*r)
        // This is the angle at the satellite in the sat-center-ground triangle

        let cos_theta = (aa * aa + r * r - bb * bb) / (2.0 * aa * r);
        log::debug!("   [BBOX DIAG] Corner {}: ISCE2-style params: aa={:.1}m (sat_mag), bb={:.1}m (R_earth at lat={:.2}°), r={:.1}m (slant_range)",
                   corner_idx, aa, bb, sat_lat, r);
        log::debug!(
            "   [BBOX DIAG] Corner {}: cos_theta={:.6}, computed as (aa²+r²-bb²)/(2·aa·r)",
            corner_idx, cos_theta
        );

        // Validate cos_theta is in valid range
        if cos_theta.abs() > 1.0 {
            log::debug!(
                "   ❌ Corner {}: Invalid geometry (cos_theta={:.4})",
                corner_idx,
                cos_theta
            );
            log::debug!(
                "      Slant range {} doesn't intersect Earth from altitude {}m",
                r,
                aa - bb
            );
            continue;
        }

        let sin_theta = (1.0 - cos_theta * cos_theta).sqrt();

        // =================================================================
        // ISCE2 RANGE VECTOR DECOMPOSITION
        // =================================================================
        // For zero-Doppler (dopfact = 0):
        //   alpha = (0 - gamma * (N·V_unit)) / (V_unit·T)
        //
        // Since V_unit is mostly along T (velocity perpendicular to nadir):
        //   V_unit·T ≈ 1, V_unit·N ≈ 0
        //   So alpha ≈ 0 for zero-Doppler

        let gamma = cos_theta * r; // Nadir component of range

        // Zero-Doppler assumption: range vector perpendicular to velocity
        // This means the along-track component alpha ≈ 0
        // (For non-zero Doppler: alpha = (dopfact - gamma*(nhat·vhat)) / (vhat·that))
        let n_dot_v = nhat[0] * vel_unit[0] + nhat[1] * vel_unit[1] + nhat[2] * vel_unit[2];
        let t_dot_v = that[0] * vel_unit[0] + that[1] * vel_unit[1] + that[2] * vel_unit[2];
        let dopfact = 0.0; // Zero-Doppler assumption

        let alpha = if t_dot_v.abs() > EPSILON {
            (dopfact - gamma * n_dot_v) / t_dot_v
        } else {
            0.0 // Degenerate case, assume alpha = 0
        };

        // Cross-track component: beta = -ILRL * sqrt(r² * sin²(theta) - alpha²)
        // ISCE2 topozero.f90 line 505: beta = -ilrl * sqrt(...)
        // ILRL = -1 for right-looking Sentinel-1
        // The negative sign is critical for correct look direction!
        let cross_track_sq = r * r * sin_theta * sin_theta - alpha * alpha;
        if cross_track_sq < 0.0 {
            log::debug!(
                "   ❌ Corner {}: Invalid cross-track (alpha² > r²sin²θ)",
                corner_idx
            );
            continue;
        }
        let beta = -ILRL * cross_track_sq.sqrt(); // NOTE: negative sign per ISCE2!

        log::debug!("   [BBOX DIAG] Corner {}: TCN components: gamma={:.1}m (nadir), alpha={:.1}m (along-track), beta={:.1}m (cross-track, ILRL={})",
                   corner_idx, gamma, alpha, beta, ILRL);

        // =================================================================
        // GROUND POINT COMPUTATION
        // =================================================================
        // Ground point = satellite + range_vector
        // range_vector = gamma*nhat + alpha*that + beta*chat
        // Note: range_vector points FROM satellite TO ground (opposite of ISCE2's delta)
        // So we ADD it to satellite position

        let range_vec = [
            gamma * nhat[0] + alpha * that[0] + beta * chat[0],
            gamma * nhat[1] + alpha * that[1] + beta * chat[1],
            gamma * nhat[2] + alpha * that[2] + beta * chat[2],
        ];

        // Ground = Satellite + range_vector (range_vector points toward ground)
        let ground_ecef = [
            sat_ecef[0] + range_vec[0],
            sat_ecef[1] + range_vec[1],
            sat_ecef[2] + range_vec[2],
        ];

        // Verify computed slant range matches input
        let computed_range = vec3_norm(range_vec[0], range_vec[1], range_vec[2]);
        let range_error = (computed_range - slant_range).abs();
        if range_error > 1.0 {
            // Allow 1m error
            log::warn!(
                "   ⚠️ Corner {}: Range error {:.1}m (computed={:.1}m, expected={:.1}m)",
                corner_idx,
                range_error,
                computed_range,
                slant_range
            );
        }

        log::debug!(
            "   [BBOX DIAG] Corner {}: Ground point ECEF({:.1}, {:.1}, {:.1})m",
            corner_idx, ground_ecef[0], ground_ecef[1], ground_ecef[2]
        );
        log::debug!(
            "   Ground ECEF: ({:.1}, {:.1}, {:.1})m",
            ground_ecef[0],
            ground_ecef[1],
            ground_ecef[2]
        );

        // Convert ECEF to lat/lon
        let (lat, lon) = ecef_to_wgs84(ground_ecef[0], ground_ecef[1], ground_ecef[2]);
        log::debug!(
            "   [BBOX DIAG] Corner {}: Ground point lat={:.6}°, lon={:.6}°",
            corner_idx, lat, lon
        );

        if lat.is_finite() && lon.is_finite() {
            min_lat = min_lat.min(lat);
            max_lat = max_lat.max(lat);
            min_lon = min_lon.min(lon);
            max_lon = max_lon.max(lon);
            successful_corners += 1;
            // CRITICAL DIAGNOSTIC: Always print corner coordinates to identify which corner is wrong
            log::error!(
                "✅ [BBOX DIAG] Corner {}: SAR({}, {}) -> geo(lat={:.6}°, lon={:.6}°)",
                corner_idx,
                row_ml,
                col_ml,
                lat,
                lon
            );
            log::debug!(
                "  ✅ Corner {}: SAR({}, {}) -> geo(lat={:.6}°, lon={:.6}°)",
                corner_idx,
                row_ml,
                col_ml,
                lat,
                lon
            );
        } else {
            log::debug!(
                "  ❌ Corner {}: Invalid lat/lon result (lat={}, lon={})",
                corner_idx,
                lat,
                lon
            );
        }
    }

    if successful_corners == 0 {
        return Err(SarError::Processing(
            "Failed to geocode any SAR image corners - cannot compute output grid extent. \
             This indicates a critical error in orbit interpolation or coordinate transformation."
                .to_string(),
        ));
    }

    if successful_corners < 4 {
        log::warn!(
            "⚠️  Only geocoded {}/4 corners, bbox may be incomplete",
            successful_corners
        );
        log::warn!("   Proceeding with partial bounding box - output may have coverage gaps");
    }

    // Add small margin (1%) to account for DEM elevation variations and approximation errors
    let lat_margin = (max_lat - min_lat).max(0.01) * 0.01;
    let lon_margin = (max_lon - min_lon).max(0.01) * 0.01;

    let computed_bbox = BoundingBox {
        min_lat: min_lat - lat_margin,
        max_lat: max_lat + lat_margin,
        min_lon: min_lon - lon_margin,
        max_lon: max_lon + lon_margin,
    };

    log::debug!("\n=== COMPUTED BBOX FROM SAR CORNERS ===");
    log::debug!(
        "   min_lat={:.6}°, max_lat={:.6}°",
        computed_bbox.min_lat, computed_bbox.max_lat
    );
    log::debug!(
        "   min_lon={:.6}°, max_lon={:.6}°",
        computed_bbox.min_lon, computed_bbox.max_lon
    );
    log::debug!(
        "   Coverage: {:.3}° lat × {:.3}° lon ({}/4 corners)",
        computed_bbox.max_lat - computed_bbox.min_lat,
        computed_bbox.max_lon - computed_bbox.min_lon,
        successful_corners
    );

    log::info!(
        "✅ Computed bbox from SAR image: [{:.6}, {:.6}, {:.6}, {:.6}]",
        computed_bbox.min_lon,
        computed_bbox.min_lat,
        computed_bbox.max_lon,
        computed_bbox.max_lat
    );
    log::debug!(
        "   Coverage: {:.3}° × {:.3}° ({}/4 corners successful)",
        computed_bbox.max_lon - computed_bbox.min_lon,
        computed_bbox.max_lat - computed_bbox.min_lat,
        successful_corners
    );

    Ok(computed_bbox)
}
