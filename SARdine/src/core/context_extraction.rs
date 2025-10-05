//! Context Extraction Utilities
//!
//! Helpers to extract ProcessingContext from existing annotation and metadata structures

// use std::collections::HashMap; // Unused import

use crate::core::processing_context::*;
use crate::io::annotation::AnnotationRoot;
use crate::types::{OrbitData, SarResult};

/// Extract ProcessingContext from annotation and orbit data
pub fn extract_from_annotation(
    annotation: &AnnotationRoot,
    orbit_data: &OrbitData,
) -> SarResult<ProcessingContext> {
    // Extract product info
    let product_info = extract_product_info(annotation)?;
    
    // Extract timing context
    let timing = extract_timing_context(annotation, orbit_data)?;
    
    // Extract radar parameters
    let radar = extract_radar_parameters(annotation)?;
    
    // Extract geometry context
    let geometry = extract_geometry_context(annotation, orbit_data)?;
    
    // Build context
    let mut ctx = ProcessingContext::new(product_info, timing, radar, geometry);
    
    // Extract burst context if TOPS mode
    if let Some(burst_ctx) = extract_burst_context(annotation)? {
        ctx.bursts = Some(burst_ctx);
    }
    
    Ok(ctx)
}

fn extract_product_info(annotation: &AnnotationRoot) -> SarResult<ProductInfo> {
    let ads_header = annotation.ads_header.as_ref()
        .ok_or_else(|| crate::types::SarError::Metadata("Missing ads_header".to_string()))?;
    
    let image_info = annotation.image_annotation.as_ref()
        .and_then(|ia| ia.image_information.as_ref());
    
    // Parse start time - MUST be present and valid
    let start_time = ads_header.start_time.as_ref()
        .and_then(|s| crate::io::annotation::parse_time_robust(s))
        .ok_or_else(|| crate::types::SarError::Metadata(
            "Missing or invalid product start time".to_string()
        ))?;
    
    // Parse stop time - MUST be present and valid
    let stop_time = ads_header.stop_time.as_ref()
        .and_then(|s| crate::io::annotation::parse_time_robust(s))
        .ok_or_else(|| crate::types::SarError::Metadata(
            "Missing or invalid product stop time".to_string()
        ))?;
    
    // Validate time ordering
    if stop_time < start_time {
        return Err(crate::types::SarError::Metadata(
            format!("Product stop time {:?} is before start time {:?}", stop_time, start_time)
        ));
    }
    
    Ok(ProductInfo {
        platform: ads_header.mission_id.clone().unwrap_or_else(|| "UNKNOWN".to_string()),
        acquisition_mode: ads_header.mode.clone().unwrap_or_else(|| "UNKNOWN".to_string()),
        product_type: ads_header.product_type.clone().unwrap_or_else(|| "UNKNOWN".to_string()),
        polarisation: ads_header.polarisation.clone().unwrap_or_else(|| "UNKNOWN".to_string()),
        swath: ads_header.swath.clone(),
        product_id: format!("{}_{}_{}", 
            ads_header.mission_id.as_deref().unwrap_or("UNK"),
            ads_header.mode.as_deref().unwrap_or("UNK"),
            start_time.format("%Y%m%dT%H%M%S")
        ),
        start_time,
        stop_time,
        pass_direction: "UNKNOWN".to_string(), // TODO: Extract from annotation if available
    })
}

fn extract_timing_context(
    annotation: &AnnotationRoot,
    orbit_data: &OrbitData,
) -> SarResult<TimingContext> {
    use crate::types::datetime_to_utc_seconds;
    
    let image_info = annotation.image_annotation.as_ref()
        .and_then(|ia| ia.image_information.as_ref())
        .ok_or_else(|| crate::types::SarError::Metadata("Missing image_information".to_string()))?;
    
    // Get orbit reference epoch - must be finite
    let orbit_ref_epoch = datetime_to_utc_seconds(orbit_data.reference_time);
    if !orbit_ref_epoch.is_finite() {
        return Err(crate::types::SarError::Metadata(
            "Orbit reference epoch is not finite".to_string()
        ));
    }
    
    // Parse product start time
    let product_start_time_abs = {
        let start_opt = image_info.product_first_line_utc_time.clone()
            .or_else(|| annotation.ads_header.as_ref().and_then(|h| h.start_time.clone()));
        let dt = start_opt
            .and_then(|s| crate::io::annotation::parse_time_robust(&s))
            .ok_or_else(|| crate::types::SarError::Metadata("Missing product start time".to_string()))?;
        (dt.timestamp() as f64) + (dt.timestamp_subsec_nanos() as f64) * 1e-9
    };
    
    // Parse product stop time - MUST be present and valid
    let product_stop_time_abs = {
        let stop_opt = image_info.product_last_line_utc_time.clone()
            .or_else(|| annotation.ads_header.as_ref().and_then(|h| h.stop_time.clone()));
        let dt = stop_opt
            .and_then(|s| crate::io::annotation::parse_time_robust(&s))
            .ok_or_else(|| crate::types::SarError::Metadata("Missing product stop time".to_string()))?;
        (dt.timestamp() as f64) + (dt.timestamp_subsec_nanos() as f64) * 1e-9
    };
    
    // Validate product times are finite
    if !product_start_time_abs.is_finite() || !product_stop_time_abs.is_finite() {
        return Err(crate::types::SarError::Metadata(
            "Product times are not finite".to_string()
        ));
    }
    
    // Validate time ordering
    if product_stop_time_abs < product_start_time_abs {
        return Err(crate::types::SarError::Metadata(
            format!("Product stop time {} is before start time {}", 
                product_stop_time_abs, product_start_time_abs)
        ));
    }
    
    let product_duration = product_stop_time_abs - product_start_time_abs;
    
    // Get PRF
    let prf = annotation.get_pulse_repetition_frequency()
        .ok_or_else(|| crate::types::SarError::Metadata("Missing PRF".to_string()))?;
    
    // Get azimuth time interval (CRITICAL - must be from annotation, not computed!)
    let azimuth_time_interval = image_info.azimuth_time_interval
        .ok_or_else(|| crate::types::SarError::Metadata(
            "Missing azimuthTimeInterval in annotation - cannot use 1/PRF fallback for TOPS".to_string()
        ))?;
    
    // Validate azimuth time interval is positive and finite
    if !azimuth_time_interval.is_finite() || azimuth_time_interval <= 0.0 {
        return Err(crate::types::SarError::Metadata(
            format!("Invalid azimuthTimeInterval: {}", azimuth_time_interval)
        ));
    }
    
    // Get range sampling rate from annotation - MUST be present
    let range_sampling_rate = annotation.general_annotation.as_ref()
        .and_then(|ga| ga.product_information.as_ref())
        .map(|pi| pi.range_sampling_rate)
        .ok_or_else(|| crate::types::SarError::Metadata(
            "Missing rangeSamplingRate in productInformation".to_string()
        ))?;
    
    // Validate range sampling rate is positive and finite
    if !range_sampling_rate.is_finite() || range_sampling_rate <= 0.0 {
        return Err(crate::types::SarError::Metadata(
            format!("Invalid rangeSamplingRate: {}", range_sampling_rate)
        ));
    }
    
    // Get slant range time
    let slant_range_time = annotation.get_slant_range_time()
        .ok_or_else(|| crate::types::SarError::Metadata("Missing slant range time".to_string()))?;
    
    // Get total lines
    let total_azimuth_lines = image_info.number_of_lines.map(|v| v as usize);
    
    Ok(TimingContext {
        orbit_ref_epoch,
        product_start_time_abs,
        product_stop_time_abs,
        product_duration,
        azimuth_time_interval,
        range_sampling_rate,
        slant_range_time,
        total_azimuth_lines,
        prf,
    })
}

fn extract_radar_parameters(annotation: &AnnotationRoot) -> SarResult<RadarParameters> {
    let (range_ps, az_ps) = annotation.get_pixel_spacing()
        .ok_or_else(|| crate::types::SarError::Metadata("Missing pixel spacing".to_string()))?;
    
    let carrier_frequency = annotation.get_radar_frequency_hz()
        .ok_or_else(|| crate::types::SarError::Metadata("Missing carrier frequency".to_string()))?;
    
    let speed_of_light = crate::constants::physical::SPEED_OF_LIGHT_M_S;
    let wavelength = speed_of_light / carrier_frequency;
    
    Ok(RadarParameters {
        carrier_frequency,
        wavelength,
        range_pixel_spacing: range_ps,
        azimuth_pixel_spacing: az_ps,
        speed_of_light,
        chirp_bandwidth: None, // TODO: Extract if available
        range_bandwidth: None, // TODO: Extract if available
    })
}

fn extract_geometry_context(
    annotation: &AnnotationRoot,
    orbit_data: &OrbitData,
) -> SarResult<GeometryContext> {
    // Try to get heading from annotation or calculate from orbit
    let azimuth_heading_rad = estimate_heading_from_orbit(orbit_data);
    
    // Try to get incidence angle from geolocation grid
    let ellipsoid_incidence_angle_rad = annotation.image_annotation.as_ref()
        .and_then(|ia| ia.image_information.as_ref())
        .and_then(|ii| ii.incidence_angle_mid_swath)
        .map(|deg| deg.to_radians())
        .unwrap_or(35.0_f64.to_radians()); // Default to typical SAR incidence
    
    let look_direction = "Right".to_string(); // Sentinel-1 is always right-looking
    
    // Try to get scene center from geolocation grid
    let (scene_center_lat, scene_center_lon) = get_scene_center(annotation);
    
    // Extract Doppler centroid
    let doppler_centroid = extract_doppler_centroid(annotation);
    
    Ok(GeometryContext {
        azimuth_heading_rad,
        ellipsoid_incidence_angle_rad,
        look_direction,
        scene_center_lat,
        scene_center_lon,
        doppler_centroid,
    })
}

fn estimate_heading_from_orbit(orbit_data: &OrbitData) -> f64 {
    // Get velocity at mid-orbit
    if orbit_data.state_vectors.is_empty() {
        return 0.0;
    }
    
    let mid_idx = orbit_data.state_vectors.len() / 2;
    let sv = &orbit_data.state_vectors[mid_idx];
    
    // Calculate heading from velocity vector (simple approximation)
    let heading = sv.velocity[1].atan2(sv.velocity[0]);
    heading
}

fn get_scene_center(annotation: &AnnotationRoot) -> (Option<f64>, Option<f64>) {
    // Try to get from geolocation grid
    if let Some(geo_grid) = annotation.geolocation_grid.as_ref()
        .and_then(|gg| gg.geolocation_grid_point_list.as_ref())
        .and_then(|gpl| gpl.geolocation_grid_points.as_ref()) {
        
        if !geo_grid.is_empty() {
            let mid_idx = geo_grid.len() / 2;
            let point = &geo_grid[mid_idx];
            return (Some(point.latitude), Some(point.longitude));
        }
    }
    
    (None, None)
}

fn extract_doppler_centroid(annotation: &AnnotationRoot) -> Option<DopplerCentroidInfo> {
    let est = annotation.general_annotation.as_ref()
        .and_then(|ga| ga.dc_estimate_list.as_ref())
        .and_then(|dl| dl.dc_estimates.as_ref())
        .and_then(|v| v.first().cloned())
        .or_else(|| {
            annotation.doppler_centroid.as_ref()
                .and_then(|dc| dc.dc_estimate_list.as_ref())
                .and_then(|dl| dl.dc_estimates.as_ref())
                .and_then(|v| v.first().cloned())
        })?;
    
    Some(DopplerCentroidInfo {
        t0: est.t0,
        coeffs: est.data_dc_polynomial,
        reference_epoch: "product_start".to_string(),
    })
}

fn extract_burst_context(annotation: &AnnotationRoot) -> SarResult<Option<BurstContext>> {
    // Check if this is TOPS mode with bursts
    let swath_timing = match annotation.swath_timing.as_ref() {
        Some(st) => st,
        None => return Ok(None),
    };
    
    let burst_list = match swath_timing.burst_list.as_ref() {
        Some(bl) => bl,
        None => return Ok(None),
    };
    
    let bursts_data = match burst_list.bursts.as_ref() {
        Some(b) => b,
        None => return Ok(None),
    };
    
    if bursts_data.is_empty() {
        return Ok(None);
    }
    
    // Get azimuth time interval for burst duration computation
    let azimuth_time_interval = annotation.image_annotation.as_ref()
        .and_then(|ia| ia.image_information.as_ref())
        .and_then(|ii| ii.azimuth_time_interval)
        .ok_or_else(|| crate::types::SarError::Metadata(
            "Missing azimuthTimeInterval for burst duration computation".to_string()
        ))?;
    
    // Get lines per burst
    let lines_per_burst = swath_timing.lines_per_burst
        .ok_or_else(|| crate::types::SarError::Metadata(
            "Missing linesPerBurst in SwathTiming".to_string()
        ))? as f64;
    
    // Extract burst metadata
    let mut bursts = Vec::new();
    for (idx, burst) in bursts_data.iter().enumerate() {
        // MUST have azimuth_time from burst metadata - no synthesis allowed
        let sensing_start = burst.azimuth_time.as_ref()
            .and_then(|s| crate::io::annotation::parse_time_robust(s))
            .ok_or_else(|| crate::types::SarError::Metadata(
                format!("Missing or invalid azimuth_time for burst {}", idx)
            ))?;
        
        // Compute sensing_stop from lines_per_burst and azimuth_time_interval
        let burst_duration_s = lines_per_burst * azimuth_time_interval;
        let sensing_stop = sensing_start + chrono::Duration::milliseconds((burst_duration_s * 1000.0) as i64);
        
        bursts.push(BurstMetadata {
            burst_id: idx,
            start_line: 0, // Will be computed during deburst
            end_line: 0,
            sensing_start,
            sensing_stop,
            azimuth_fm_rate: 0.0, // TODO: Extract if available
            azimuth_steering_rate: 0.0, // TODO: Extract if available
            doppler_centroid: 0.0, // TODO: Extract if available
        });
    }
    
    Ok(Some(BurstContext {
        num_bursts: bursts.len(),
        bursts,
        overlap_windows: Vec::new(), // Will be computed during deburst
    }))
}

/// Update context with calibration information
pub fn update_with_calibration(
    ctx: &mut ProcessingContext,
    cal_type: &str,
    cal_constant: f64,
    lut_ids: Vec<String>,
    noise_removal: bool,
) {
    ctx.calibration = CalibrationContext {
        calibration_type: cal_type.to_string(),
        calibration_constant: cal_constant,
        lut_ids,
        noise_removal_applied: noise_removal,
        nesz_stats: None,
        output_units: "linear".to_string(),
    };
}

/// Update context with DEM information
pub fn update_with_dem(
    ctx: &mut ProcessingContext,
    dem_source: String,
    dem_crs: String,
    dem_crs_epsg: u32,
    dem_spacing: f32,
    vertical_datum: String,
    min_elev: f32,
    max_elev: f32,
    mean_elev: f32,
    std_elev: f32,
    bounds: (f64, f64, f64, f64),
) {
    ctx.dem = Some(DemContext {
        dem_source,
        dem_crs,
        dem_crs_epsg,
        dem_spacing_meters: dem_spacing,
        vertical_datum,
        elevation_stats: ElevationStatistics {
            min_m: min_elev,
            max_m: max_elev,
            mean_m: mean_elev,
            std_m: std_elev,
        },
        dem_bounds: DemBounds {
            min_lat: bounds.0,
            max_lat: bounds.1,
            min_lon: bounds.2,
            max_lon: bounds.3,
        },
    });
}

#[cfg(test)]
mod tests {
    use super::*;
    use chrono::{Utc, TimeZone};
    
    #[test]
    fn test_heading_estimation() {
        // Use fixed time instead of Utc::now() for deterministic tests
        let fixed_time = Utc.with_ymd_and_hms(2020, 1, 1, 0, 0, 0).unwrap();
        
        let orbit_data = OrbitData {
            reference_time: fixed_time,
            state_vectors: vec![
                crate::types::StateVector {
                    time: fixed_time,
                    position: [0.0, 0.0, 7000000.0],
                    velocity: [1000.0, 7000.0, 0.0],
                },
            ],
        };
        
        let heading = estimate_heading_from_orbit(&orbit_data);
        assert!(heading.abs() < std::f64::consts::PI);
    }
}
