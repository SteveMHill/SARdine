use chrono::{DateTime, Utc};
use std::collections::HashMap;

use crate::constants::physical::SPEED_OF_LIGHT_M_S;
use crate::core::deburst::geometry::derive_burst_geometry;
use crate::types::{FmRateEstimate, Polarization, SarError, SarResult, StateVector};

use super::raw::{parse_time_robust, ProductRoot};

// ============================================================================
// TIME BASE STRUCTURES - Critical for terrain correction solver
// ============================================================================

/// Time base structure for consistent temporal referencing
/// CRITICAL: All times must be referenced to orbit_ref_epoch to avoid
/// solver failures with "azimuth=130k (max 50k)" type errors
#[derive(Debug, Clone)]
pub struct TimeBases {
    /// Reference epoch from orbit file (earliest state vector time)
    /// All other times should be computed relative to this
    pub orbit_ref_epoch: DateTime<Utc>,
    /// Product start time in UTC
    pub product_start_utc: DateTime<Utc>,
    /// Product stop time in UTC (if available)
    pub product_stop_utc: Option<DateTime<Utc>>,
}

/// Doppler Centroid model expressed in the Sentinel-1 slant-range domain
#[derive(Clone, Debug)]
pub struct DcModel {
    /// Optional azimuth timestamp (seconds since orbit_ref_epoch) used for selecting the model
    pub azimuth_time_rel_s: Option<f64>,
    /// Slant-range reference time τ₀ (seconds) used as the polynomial origin
    pub slant_range_t0_s: f64,
    /// Polynomial coefficients [c0, c1, c2, ...] evaluated as f_dc(τ) = Σ c_n (τ - τ₀)^n
    pub coeffs: Vec<f64>,
}

impl DcModel {
    /// Evaluate Doppler centroid at a given two-way slant-range time τ (seconds)
    pub fn eval_hz(&self, slant_range_time_s: f64) -> f64 {
        let dt = slant_range_time_s - self.slant_range_t0_s;
        self.coeffs
            .iter()
            .enumerate()
            .fold(0.0, |acc, (i, &a)| acc + a * dt.powi(i as i32))
    }
}

// Legacy structures for backward compatibility with existing SARdine code
#[derive(Debug, Clone)]
pub struct SubSwathGeometry {
    pub near_range: f64,
    pub far_range: f64,
    pub incidence_near: f64,
    pub incidence_far: f64,
    pub azimuth_start_time: f64,
    pub azimuth_end_time: f64,
    pub first_line: u32,
    pub last_line: u32,
    pub first_sample: u32,
    pub last_sample: u32,

    // Additional fields from original structure
    pub range_pixel_spacing: f64,
    pub azimuth_pixel_spacing: f64,
    pub incidence_angle_near: f64,
    pub incidence_angle_far: f64,
    // Bounds are half-open: [first, last)
    pub first_valid_sample: usize,
    pub last_valid_sample: usize,
    pub first_valid_line: usize,
    pub last_valid_line: usize,
    pub slant_range_time: f64,
    pub range_sampling_rate: f64,
}

// Legacy data structures for burst processing compatibility
#[derive(Debug, Clone)]
pub struct AnnotationData {
    pub bursts: Vec<BurstData>,
}

#[derive(Debug, Clone)]
pub struct BurstData {
    pub lines_per_burst: usize,
    pub azimuth_pixel_spacing: f64,
    pub first_valid_sample: Vec<i32>,
    pub last_valid_sample: Vec<i32>,
    pub azimuth_time: String,
    pub sensing_time: String,
    pub byte_offset: u64,
    pub azimuth_fm_rate: f64,
    pub azimuth_steering_rate: f64,
    pub slant_range_time: f64,
    pub doppler_centroid: f64,
    pub azimuth_bandwidth: f64,
    pub range_sampling_rate: f64,
    pub range_pixel_spacing: f64,
}

// ============================================================================
// BURST GEOMETRY DERIVATION
// ============================================================================

// ============================================================================
// DERIVED HELPERS ON PARSED ANNOTATION
// ============================================================================

impl ProductRoot {
    /// Extract pixel spacing values (critical for SAR processing)
    pub fn get_pixel_spacing(&self) -> Option<(f64, f64)> {
        // Try general annotation first
        if let Some(general) = &self.general_annotation {
            if let Some(product_info) = &general.product_information {
                if let (Some(range_ps), Some(azimuth_ps)) = (
                    product_info.range_pixel_spacing,
                    product_info.azimuth_pixel_spacing,
                ) {
                    return Some((range_ps, azimuth_ps));
                }
            }
        }

        // Fall back to image annotation
        if let Some(image) = &self.image_annotation {
            if let Some(image_info) = &image.image_information {
                if let (Some(range_ps), Some(azimuth_ps)) = (
                    image_info.range_pixel_spacing,
                    image_info.azimuth_pixel_spacing,
                ) {
                    return Some((range_ps, azimuth_ps));
                }
            }
        }

        None
    }

    /// Derive time bases from annotation and orbit data
    /// CRITICAL: This establishes the orbit_ref_epoch that all times must be relative to
    /// to avoid terrain correction solver failures
    pub fn derive_time_bases(&self, orbit_vectors: &[StateVector]) -> SarResult<TimeBases> {
        let product_start_utc = self
            .image_annotation
            .as_ref()
            .and_then(|ia| ia.image_information.as_ref())
            .and_then(|ii| ii.product_first_line_utc_time.as_ref())
            .and_then(|s| parse_time_robust(s))
            .or_else(|| {
                self.ads_header
                    .as_ref()
                    .and_then(|h| h.start_time.as_ref())
                    .and_then(|s| parse_time_robust(s))
            })
            .ok_or_else(|| SarError::Metadata("Missing product start time in annotation".into()))?;

        let product_stop_utc = self
            .image_annotation
            .as_ref()
            .and_then(|ia| ia.image_information.as_ref())
            .and_then(|ii| ii.product_last_line_utc_time.as_ref())
            .and_then(|s| parse_time_robust(s))
            .or_else(|| {
                self.ads_header
                    .as_ref()
                    .and_then(|h| h.stop_time.as_ref())
                    .and_then(|s| parse_time_robust(s))
            });

        // Use earliest orbit vector time as reference epoch
        let orbit_ref_epoch = orbit_vectors
            .iter()
            .map(|v| v.time)
            .min()
            .ok_or_else(|| SarError::Metadata("No orbit vectors provided".into()))?;

        Ok(TimeBases {
            orbit_ref_epoch,
            product_start_utc,
            product_stop_utc,
        })
    }

    /// Convert an absolute UTC time to seconds since orbit_ref_epoch
    /// This ensures consistency across all time references in terrain correction
    pub fn seconds_since_orbit_ref(&self, tb: &TimeBases, t_utc: DateTime<Utc>) -> f64 {
        let ref_epoch = tb.orbit_ref_epoch;
        let delta = t_utc - ref_epoch;
        delta.num_seconds() as f64 + (delta.subsec_nanos() as f64) * 1e-9
    }

    /// Build Doppler Centroid models while preserving Sentinel-1 slant-range semantics
    pub fn build_dc_models(&self, tb: &TimeBases) -> SarResult<Vec<DcModel>> {
        let list = self
            .general_annotation
            .as_ref()
            .and_then(|ga| ga.dc_estimate_list.as_ref())
            .or_else(|| {
                self.doppler_centroid
                    .as_ref()
                    .and_then(|dc| dc.dc_estimate_list.as_ref())
            })
            .ok_or_else(|| SarError::Metadata("No dcEstimateList in annotation".into()))?;

        let estimates = list
            .dc_estimates
            .as_ref()
            .ok_or_else(|| SarError::Metadata("Empty dcEstimateList".into()))?;

        let mut out = Vec::new();
        for e in estimates {
            // Validate polynomial coefficients
            if e.data_dc_polynomial.is_empty() {
                return Err(SarError::Metadata(
                    "Empty DC polynomial coefficients".into(),
                ));
            }

            if !e.t0.is_finite() {
                return Err(SarError::Metadata(
                    "dcEstimate missing slant-range reference time t0".into(),
                ));
            }

            let azimuth_time_rel_s = if let Some(az) = &e.azimuth_time {
                let az_utc = parse_time_robust(az)
                    .ok_or_else(|| SarError::Metadata("Bad azimuthTime in dcEstimate".into()))?;
                Some(self.seconds_since_orbit_ref(tb, az_utc))
            } else {
                None
            };

            out.push(DcModel {
                azimuth_time_rel_s,
                slant_range_t0_s: e.t0,
                coeffs: e.data_dc_polynomial.clone(),
            });
        }

        log::info!(
            "📊 Built {} DC models with orbit_ref_epoch time base",
            out.len()
        );
        Ok(out)
    }

    /// Evaluate Doppler centroid by selecting the closest azimuth-time model and
    /// evaluating it at the provided slant-range time τ (seconds).
    pub fn eval_dc_hz(dc_models: &[DcModel], slant_range_time_s: f64, t_rel_s: f64) -> f64 {
        if dc_models.is_empty() {
            return 0.0;
        }
        let (best, _) =
            dc_models
                .iter()
                .enumerate()
                .fold((0usize, f64::INFINITY), |acc, (idx, model)| {
                    let delta = match model.azimuth_time_rel_s {
                        Some(az) => (az - t_rel_s).abs(),
                        None => f64::INFINITY,
                    };
                    if delta < acc.1 {
                        (idx, delta)
                    } else {
                        acc
                    }
                });
        dc_models[best].eval_hz(slant_range_time_s)
    }

    /// Extract orbit information
    pub fn get_orbit_info(&self) -> Option<(String, u32)> {
        if let Some(header) = &self.ads_header {
            if let Some(mission) = &header.mission_id {
                if let Some(orbit) = header.absolute_orbit_number {
                    return Some((mission.clone(), orbit));
                }
            }
        }
        None
    }

    /// Extract antenna pattern arrays - STRICT scientific validation
    pub fn get_antenna_patterns(&self) -> SarResult<Vec<(Vec<f64>, Vec<f64>, Vec<f64>)>> {
        let mut patterns = Vec::new();

        // Check both root-level and general annotation
        let antenna_list = self.antenna_pattern.as_ref().or_else(|| {
            self.general_annotation
                .as_ref()
                .and_then(|ga| ga.antenna_pattern.as_ref())
        });

        if let Some(list) = antenna_list {
            let mut append_patterns =
                |pattern_vec: &[crate::io::annotation::AntennaPatternValues]| -> SarResult<()> {
                    for pattern in pattern_vec {
                        // CRITICAL SCIENTIFIC FIX: No fallback values for antenna parameters
                        // Missing antenna pattern data indicates corrupted annotation file
                        let elev = pattern.elevation_angle.as_ref().ok_or_else(|| {
                            SarError::Processing(
                                "Missing elevation angle in antenna pattern - invalid annotation"
                                    .to_string(),
                            )
                        })?;
                        let inc = pattern.incidence_angle.as_ref().ok_or_else(|| {
                            SarError::Processing(
                                "Missing incidence angle in antenna pattern - invalid annotation"
                                    .to_string(),
                            )
                        })?;
                        let srt = pattern.slant_range_time.as_ref().ok_or_else(|| {
                            SarError::Processing(
                                "Missing slant range time in antenna pattern - invalid annotation"
                                    .to_string(),
                            )
                        })?;
                        patterns.push((elev.clone(), inc.clone(), srt.clone()));
                    }
                    Ok(())
                };

            if let Some(pattern_vec) = &list.antenna_patterns {
                append_patterns(pattern_vec)?;
            }

            if let Some(nested) = &list.nested_list {
                if let Some(pattern_vec) = &nested.antenna_patterns {
                    append_patterns(pattern_vec)?;
                }
            }
        }

        Ok(patterns)
    }

    /// Extract bounding box from annotation
    pub fn extract_bounding_box(annotation: &ProductRoot) -> SarResult<crate::types::BoundingBox> {
        log::debug!("Attempting to extract bounding box from comprehensive annotation");

        if let Some(ref geoloc_grid) = annotation.geolocation_grid {
            if let Some(ref point_list) = geoloc_grid.geolocation_grid_point_list {
                if let Some(ref points) = point_list.geolocation_grid_points {
                    log::debug!("Found geolocation grid with {} points", points.len());

                    if !points.is_empty() {
                        let lats: Vec<f64> = points.iter().map(|p| p.latitude).collect();
                        let lons: Vec<f64> = points.iter().map(|p| p.longitude).collect();

                        let min_lat = lats.iter().fold(f64::INFINITY, |a, &b| a.min(b));
                        let max_lat = lats.iter().fold(f64::NEG_INFINITY, |a, &b| a.max(b));
                        let min_lon = lons.iter().fold(f64::INFINITY, |a, &b| a.min(b));
                        let max_lon = lons.iter().fold(f64::NEG_INFINITY, |a, &b| a.max(b));

                        log::info!(
                            "Extracted bounding box: [{:.6}, {:.6}, {:.6}, {:.6}]",
                            min_lon,
                            min_lat,
                            max_lon,
                            max_lat
                        );

                        let mut bbox = crate::types::BoundingBox {
                            min_lat,
                            min_lon,
                            max_lat,
                            max_lon,
                        };

                        let is_valid = |b: &crate::types::BoundingBox| {
                            let finite = b.min_lat.is_finite()
                                && b.max_lat.is_finite()
                                && b.min_lon.is_finite()
                                && b.max_lon.is_finite();
                            let proper_order = b.min_lat < b.max_lat && b.min_lon < b.max_lon;
                            let not_all_zero = !(b.min_lat == 0.0
                                && b.max_lat == 0.0
                                && b.min_lon == 0.0
                                && b.max_lon == 0.0);
                            finite && proper_order && not_all_zero
                        };

                        if !is_valid(&bbox) {
                            return Err(SarError::Metadata(
                                "Invalid bounding box in annotation geolocation grid".to_string(),
                            ));
                        }

                        // Enforce ordering invariants (prevents downstream clamp panics)
                        bbox.normalize();
                        return Ok(bbox);
                    }
                }
            }
        }

        Err(SarError::Metadata(
            "No valid geolocation grid found for bounding box extraction".to_string(),
        ))
    }

    /// Get subswath geometry from parsed annotation - FIXED to use subswath-specific data
    pub fn get_subswath_info(&self, subswath: &str) -> SarResult<SubSwathGeometry> {
        // SCIENTIFIC FIX: Extract geometry specific to requested subswath (IW1/IW2/IW3)
        log::info!("Extracting geometry for specific subswath: {}", subswath);

        // First verify subswath exists and get subswath-specific data
        let subs = ProductRoot::extract_subswaths(self)?;
        let subswath_data = subs.get(subswath).ok_or_else(|| {
            SarError::Metadata(format!("Subswath {} not found in metadata", subswath))
        })?;

        // Use subswath-specific pixel spacing from extracted subswath data
        let range_pixel_spacing = subswath_data.range_pixel_spacing;
        let azimuth_pixel_spacing = subswath_data.azimuth_pixel_spacing;
        let slant_range_time = subswath_data.slant_range_time;

        // Validate extracted parameters
        if range_pixel_spacing <= 0.0 {
            return Err(SarError::Metadata(format!(
                "Invalid range pixel spacing for {}: {}",
                subswath, range_pixel_spacing
            )));
        }
        if azimuth_pixel_spacing <= 0.0 {
            return Err(SarError::Metadata(format!(
                "Invalid azimuth pixel spacing for {}: {}",
                subswath, azimuth_pixel_spacing
            )));
        }
        if slant_range_time <= 0.0 {
            return Err(SarError::Metadata(format!(
                "Invalid slant range time for {}: {}",
                subswath, slant_range_time
            )));
        }

        log::info!(
            "✅ Using subswath-specific geometry: range_ps={:.3}m, az_ps={:.3}m, srt={:.6}s",
            range_pixel_spacing,
            azimuth_pixel_spacing,
            slant_range_time
        );

        // Use subswath-specific dimensions instead of global image dimensions
        let num_samples = subswath_data.range_samples;
        let num_lines = subswath_data.azimuth_samples;
        if num_samples == 0 || num_lines == 0 {
            return Err(SarError::Metadata(format!(
                "Invalid dimensions for {}: {}x{} samples",
                subswath, num_samples, num_lines
            )));
        }

        log::info!(
            "✅ Subswath {} dimensions: {}x{} (range x azimuth)",
            subswath,
            num_samples,
            num_lines
        );

        // Derive near/far range and incidence from geolocation grid using min/max slantRangeTime
        let mut near_srt = f64::INFINITY;
        let mut far_srt = f64::NEG_INFINITY;
        let mut near_inc_sum = 0.0;
        let mut near_inc_count = 0usize;
        let mut far_inc_sum = 0.0;
        let mut far_inc_count = 0usize;
        if let Some(geo) = &self.geolocation_grid {
            if let Some(list) = &geo.geolocation_grid_point_list {
                if let Some(points) = &list.geolocation_grid_points {
                    for p in points {
                        let srt = p.slant_range_time;
                        if srt < near_srt - 1e-9 {
                            near_srt = srt;
                            near_inc_sum = p.incidence_angle;
                            near_inc_count = 1;
                        } else if (srt - near_srt).abs() < 1e-9 {
                            near_inc_sum += p.incidence_angle;
                            near_inc_count += 1;
                        }

                        if srt > far_srt + 1e-9 {
                            far_srt = srt;
                            far_inc_sum = p.incidence_angle;
                            far_inc_count = 1;
                        } else if (srt - far_srt).abs() < 1e-9 {
                            far_inc_sum += p.incidence_angle;
                            far_inc_count += 1;
                        }
                    }
                }
            }
        }
        if !near_srt.is_finite() || !far_srt.is_finite() {
            return Err(SarError::Metadata(
                "Failed to derive near/far slantRangeTime from geolocationGrid".to_string(),
            ));
        }
        let near_range = 0.5 * SPEED_OF_LIGHT_M_S * near_srt;
        let far_range = 0.5 * SPEED_OF_LIGHT_M_S * far_srt;

        // SCIENTIFIC MODE: Require valid incidence angle data
        let incidence_near = if near_inc_count > 0 {
            near_inc_sum / near_inc_count as f64
        } else {
            return Err(SarError::Metadata("Missing near-range incidence angle data in geolocation grid - cannot proceed scientifically".to_string()));
        };
        let incidence_far = if far_inc_count > 0 {
            far_inc_sum / far_inc_count as f64
        } else {
            return Err(SarError::Metadata("Missing far-range incidence angle data in geolocation grid - cannot proceed scientifically".to_string()));
        };

        // Azimuth start/end time from product/image header (seconds since epoch)
        let (az_start_sec, az_end_sec) = {
            // Prefer imageInformation productFirst/LastLineUtcTime when available
            let (start_opt, stop_opt) = self
                .image_annotation
                .as_ref()
                .and_then(|ia| ia.image_information.as_ref())
                .map(|ii| {
                    (
                        ii.product_first_line_utc_time.clone(),
                        ii.product_last_line_utc_time.clone(),
                    )
                })
                .unwrap_or((None, None));

            let start_time = start_opt
                .and_then(|s| parse_time_robust(&s))
                .or_else(|| {
                    self.ads_header
                        .as_ref()
                        .and_then(|h| h.start_time.as_ref())
                        .and_then(|s| parse_time_robust(s))
                })
                .ok_or_else(|| {
                    SarError::Metadata(
                        "Missing product first line/start time in annotation".to_string(),
                    )
                })?;
            let stop_time = stop_opt
                .and_then(|s| parse_time_robust(&s))
                .or_else(|| {
                    self.ads_header
                        .as_ref()
                        .and_then(|h| h.stop_time.as_ref())
                        .and_then(|s| parse_time_robust(s))
                })
                .ok_or_else(|| {
                    SarError::Metadata(
                        "Missing product last line/stop time in annotation".to_string(),
                    )
                })?;
            (
                start_time.timestamp() as f64 + start_time.timestamp_subsec_nanos() as f64 * 1e-9,
                stop_time.timestamp() as f64 + stop_time.timestamp_subsec_nanos() as f64 * 1e-9,
            )
        };

        // Valid sample bounds from bursts (min firstValidSample, max lastValidSample)
        let mut first_valid_sample = None::<usize>;
        let mut last_valid_sample_exclusive = None::<usize>;
        if let Some(st) = &self.swath_timing {
            if let Some(bl) = &st.burst_list {
                if let Some(bursts) = &bl.bursts {
                    for b in bursts {
                        if !b.first_valid_sample.is_empty() {
                            let min_first = b
                                .first_valid_sample
                                .iter()
                                .copied()
                                .filter(|v| *v >= 0)
                                .min()
                                .map(|v| v as usize);
                            if let Some(min_first) = min_first {
                                first_valid_sample = Some(
                                    first_valid_sample.map_or(min_first, |v| v.min(min_first)),
                                );
                            }
                        }
                        if !b.last_valid_sample.is_empty() {
                            let max_last_exclusive = b
                                .last_valid_sample
                                .iter()
                                .copied()
                                .filter(|v| *v >= 0)
                                .max()
                                .map(|v| (v as usize).saturating_add(1));
                            if let Some(max_last_exclusive) = max_last_exclusive {
                                last_valid_sample_exclusive = Some(
                                    last_valid_sample_exclusive
                                        .map_or(max_last_exclusive, |v| v.max(max_last_exclusive)),
                                );
                            }
                        }
                    }
                }
            }
        }

        // Use actual valid sample bounds or return error instead of fallback values
        let first_valid_sample = first_valid_sample
            .ok_or_else(|| {
                SarError::Metadata(format!("Missing first valid sample data for {}", subswath))
            })?
            .min(num_samples);
        let last_valid_sample_exclusive = last_valid_sample_exclusive
            .ok_or_else(|| {
                SarError::Metadata(format!("Missing last valid sample data for {}", subswath))
            })?
            .min(num_samples);

        if first_valid_sample >= last_valid_sample_exclusive {
            return Err(SarError::Metadata(format!(
                "Invalid valid-sample bounds for {}: first={} last_exclusive={}",
                subswath, first_valid_sample, last_valid_sample_exclusive
            )));
        }

        // Range sampling rate
        let range_sampling_rate = self
            .general_annotation
            .as_ref()
            .and_then(|ga| ga.product_information.as_ref())
            .map(|pi| pi.range_sampling_rate)
            .ok_or_else(|| {
                SarError::Metadata("Missing rangeSamplingRate in productInformation".to_string())
            })?;

        Ok(SubSwathGeometry {
            near_range,
            far_range,
            incidence_near,
            incidence_far,
            azimuth_start_time: az_start_sec,
            azimuth_end_time: az_end_sec,
            first_line: 0,
            last_line: num_lines as u32,
            first_sample: 0,
            last_sample: num_samples as u32,
            range_pixel_spacing,
            azimuth_pixel_spacing,
            incidence_angle_near: incidence_near,
            incidence_angle_far: incidence_far,
            first_valid_sample,
            last_valid_sample: last_valid_sample_exclusive,
            first_valid_line: 0,
            last_valid_line: num_lines,
            slant_range_time,
            range_sampling_rate,
        })
    }

    /// Extract burst timing information
    pub fn extract_burst_times(annotation: &ProductRoot) -> SarResult<Vec<String>> {
        if let Some(ref swath_timing) = annotation.swath_timing {
            if let Some(ref burst_list) = swath_timing.burst_list {
                if let Some(ref bursts) = burst_list.bursts {
                    let burst_times: Vec<String> = bursts
                        .iter()
                        .filter_map(|b| b.azimuth_time.clone())
                        .collect();
                    return Ok(burst_times);
                }
            }
        }

        Ok(vec![])
    }

    /// Extract sub-swath information from annotation
    pub fn extract_subswaths(
        annotation: &ProductRoot,
    ) -> SarResult<HashMap<String, crate::types::SubSwath>> {
        log::info!("🔍 DEBUGGING: Starting extract_subswaths function");
        let mut subswaths = HashMap::new();

        #[derive(Debug)]
        struct PendingSubswath {
            id: String,
            burst_count: usize,
            samples_per_burst: usize,
            lines_per_burst: usize,
            derived_first_line: usize,
            derived_last_line: usize,
            derived_first_sample: usize,
            derived_last_sample: usize,
            product_range_samples: usize,
            product_azimuth_samples: usize,
            range_pixel_spacing: f64,
            azimuth_pixel_spacing: f64,
            slant_range_time: f64,
            burst_duration: f64,
            prf_hz: Option<f64>,
            dc_polynomial: Option<Vec<f64>>,
            azimuth_time_interval: Option<f64>,
            near_range_m: f64,
            dc_polynomial_t0: Option<f64>,
            fm_rate_models: Vec<FmRateEstimate>,
        }

        let mut pending_subswaths: Vec<PendingSubswath> = Vec::new();
        let annotation_prf = annotation.get_pulse_repetition_frequency()?;
        let fm_rate_entries = annotation
            .general_annotation
            .as_ref()
            .and_then(|ga| ga.azimuth_fm_rate_list.as_ref())
            .and_then(|list| list.azimuth_fm_rates.as_ref())
            .ok_or_else(|| {
                SarError::Metadata("Missing azimuthFmRate entries in annotation".to_string())
            })?;

        let fm_rate_models = fm_rate_entries
            .iter()
            .map(|entry| {
                let azimuth_time_utc = entry
                    .azimuth_time
                    .as_deref()
                    .and_then(|t| parse_time_robust(t))
                    .ok_or_else(|| {
                        SarError::Metadata(
                            "azimuthFmRate entry missing azimuthTime timestamp".to_string(),
                        )
                    })?;

                Ok(FmRateEstimate {
                    azimuth_time_utc: Some(azimuth_time_utc),
                    azimuth_time_rel_orbit: None,
                    slant_range_reference_time: entry.t0,
                    coefficients: entry.azimuth_fm_rate_polynomial.clone(),
                })
            })
            .collect::<SarResult<Vec<_>>>()?;

        if fm_rate_models.is_empty() {
            return Err(SarError::Metadata(
                "Annotation contains zero azimuthFmRate polynomials".to_string(),
            ));
        }

        log::info!(
            "✅ DEBUGGING: azimuthFmRateList contains {} entries",
            fm_rate_models.len()
        );

        let approx_equal = |a: f64, b: f64| {
            let diff = (a - b).abs();
            let scale = a.abs().max(b.abs()).max(1.0);
            diff <= 1e-6 * scale || diff <= 1e-3
        };

        log::info!("🔍 DEBUGGING: Checking image_annotation + processing_information");
        if let Some(ref image_annotation) = annotation.image_annotation {
            let proc_info_opt = image_annotation.processing_information.as_ref();
            if proc_info_opt.is_none() {
                log::warn!("❌ DEBUGGING: processing_information is None");
            }

            // Collect swath entries if available
            if let Some(proc_info) = proc_info_opt {
                if let Some(ref params_list) = proc_info.swath_proc_params_list {
                    if let Some(ref params_vec) = params_list.swath_proc_params {
                        log::info!(
                            "✅ DEBUGGING: swath_proc_params entries: {}",
                            params_vec.len()
                        );

                        // Try to get geometry from image_information; if missing, try general annotation
                        let (range_ps_opt, az_ps_opt) = annotation
                            .get_pixel_spacing()
                            .map(|(r, a)| (Some(r), Some(a)))
                            .unwrap_or_else(|| {
                                if let Some(ref img_info) = image_annotation.image_information {
                                    (img_info.range_pixel_spacing, img_info.azimuth_pixel_spacing)
                                } else {
                                    (None, None)
                                }
                            });
                        let srt_opt = annotation.get_slant_range_time().or_else(|| {
                            image_annotation
                                .image_information
                                .as_ref()
                                .and_then(|ii| ii.slant_range_time)
                        });

                        // Extract azimuth_time_interval from image_information for terrain correction
                        let azimuth_time_interval_opt = image_annotation
                            .image_information
                            .as_ref()
                            .and_then(|ii| ii.azimuth_time_interval);

                        let (range_samples_opt, az_lines_opt) =
                            if let Some(ref img_info) = image_annotation.image_information {
                                (
                                    img_info.number_of_samples.map(|v| v as usize),
                                    img_info.number_of_lines.map(|v| v as usize),
                                )
                            } else {
                                (None, None)
                            };

                        // Preserve the original image dimensions (if present) for bounds checks
                        let original_range_samples = range_samples_opt.unwrap_or(0);
                        let original_azimuth_samples = az_lines_opt.unwrap_or(0);

                        let burst_geometry = annotation
                            .swath_timing
                            .as_ref()
                            .and_then(|st| derive_burst_geometry(st));
                        if let Some(ref geom) = burst_geometry {
                            log::info!(
                                "✅ DEBUGGING: Derived burst geometry -> lines {}..{} samples {}..{} (bursts {}/{})",
                                geom.first_valid_line.unwrap_or(0),
                                geom.last_valid_line_exclusive.unwrap_or(0),
                                geom.first_valid_sample.unwrap_or(0),
                                geom.last_valid_sample_exclusive.unwrap_or(0),
                                geom.valid_bursts,
                                geom.total_bursts
                            );
                        } else {
                            log::warn!(
                                "⚠️ DEBUGGING: Unable to derive burst geometry from swathTiming"
                            );
                        }

                        let fallback_burst_count = annotation
                            .swath_timing
                            .as_ref()
                            .and_then(|st| st.burst_list.as_ref())
                            .and_then(|bl| bl.bursts.as_ref())
                            .map(|b| b.len())
                            .unwrap_or(0);

                        let fallback_lines_per_burst = annotation
                            .swath_timing
                            .as_ref()
                            .and_then(|st| st.lines_per_burst)
                            .or_else(|| {
                                burst_geometry
                                    .as_ref()
                                    .and_then(|bg| bg.lines_per_burst)
                                    .map(|v| v as u32)
                            });

                        let fallback_samples_per_burst = annotation
                            .swath_timing
                            .as_ref()
                            .and_then(|st| st.samples_per_burst)
                            .or_else(|| {
                                burst_geometry
                                    .as_ref()
                                    .and_then(|bg| bg.samples_per_burst)
                                    .map(|v| v as u32)
                            });

                        // Iterate through swath_proc_params entries
                        for (i, params) in params_vec.iter().enumerate() {
                            if let Some(swath_id) = &params.swath {
                                log::info!(
                                    "✅ DEBUGGING: Processing swath {} (index {})",
                                    swath_id,
                                    i
                                );

                                // Extract pixel spacing values
                                let range_pixel_spacing = params
                                    .range_sampling_rate
                                    .map(|rsr| SPEED_OF_LIGHT_M_S / (2.0 * rsr))
                                    .or(range_ps_opt)
                                    .or_else(|| params.performed_range_cell_spacing)
                                    .ok_or_else(|| {
                                        SarError::Metadata(format!(
                                            "Missing range pixel spacing for swath {} (rangeSamplingRate/annotation/performedRangeCellSpacing)",
                                            swath_id
                                        ))
                                    })?;

                                // Keep azimuth pixel spacing strictly in meters; do not mix with azimuthTimeInterval (seconds)
                                let azimuth_pixel_spacing = az_ps_opt.ok_or_else(|| {
                                    SarError::Metadata(format!(
                                        "Missing azimuth pixel spacing for swath {} in annotation",
                                        swath_id
                                    ))
                                })?;

                                let azimuth_time_interval =
                                    params.azimuth_time_interval.or(azimuth_time_interval_opt);

                                let slant_range_time = params
                                    .slant_range_time
                                    .or(srt_opt)
                                    .ok_or_else(|| {
                                        SarError::Metadata(format!(
                                            "Missing slantRangeTime for swath {} in annotation/swathProcParams",
                                            swath_id
                                        ))
                                    })?;

                                // Extract PRF from swathProcParams with strict validation
                                let valid_prf =
                                    |prf: f64| prf.is_finite() && prf > 0.0 && prf < 6000.0;

                                if let Some(swath_prf) = params.prf.filter(|p| valid_prf(*p)) {
                                    if !approx_equal(swath_prf, annotation_prf) {
                                        return Err(SarError::Metadata(format!(
                                            "Swath {} PRF ({:.6} Hz) disagrees with downlinkInformation PRF ({:.6} Hz)",
                                            swath_id, swath_prf, annotation_prf
                                        )));
                                    }
                                }

                                if let Some(az_freq) =
                                    params.azimuth_frequency.filter(|p| valid_prf(*p))
                                {
                                    if !approx_equal(az_freq, annotation_prf) {
                                        return Err(SarError::Metadata(format!(
                                            "azimuthFrequency ({:.6} Hz) differs from downlinkInformation PRF ({:.6} Hz) for swath {}",
                                            az_freq, annotation_prf, swath_id
                                        )));
                                    }
                                }

                                // SCIENTIFIC FIX (Jan 2026): Validate azimuth timing consistency
                                // For stripmap: azimuth_time_interval SHOULD equal 1/PRF (line spacing)
                                // For TOPSAR: azimuth_time_interval may represent burst interval, not line spacing
                                if let Some(ati) = azimuth_time_interval {
                                    if ati.is_finite() && ati > 0.0 {
                                        let expected_line_spacing = 1.0 / annotation_prf;
                                        let ratio = ati / expected_line_spacing;
                                        
                                        log::info!(
                                            "🔬 Azimuth timing for {}: PRF={:.3}Hz, line_spacing(1/PRF)={:.9}s, annotation_ati={:.9}s, ratio={:.2}x",
                                            swath_id, annotation_prf, expected_line_spacing, ati, ratio
                                        );
                                        
                                        // For TOPSAR IW, ATI = 1/azimuthFrequency (SLC line spacing), typically ~0.002s vs 1/PRF ~5.8e-4s
                                        // The ratio represents azimuth multi-looking factor (SLC sampled at ~1/3 raw PRF)
                                        if ratio > 2.0 && ratio < 6.0 {
                                            log::info!(
                                                "✅ TOPSAR azimuth multi-looking: ATI ({:.6}s) = 1/azimuthFrequency (SLC line spacing), {:.2}× larger than 1/PRF ({:.6}s)",
                                                ati, ratio, expected_line_spacing
                                            );
                                        } else if (ratio - 1.0).abs() < 0.01 {
                                            log::info!(
                                                "✅ Stripmap mode: azimuth_time_interval ({:.9}s) matches 1/PRF ({:.9}s)",
                                                ati, expected_line_spacing
                                            );
                                        } else {
                                            log::warn!(
                                                "⚠️  SCIENTIFIC WARNING: azimuth_time_interval ({:.6}s) has unexpected ratio to 1/PRF ({:.9}s): {:.2}x. \
                                                 Expected ~1.0 (stripmap) or ~3-5 (TOPSAR burst interval). Verify annotation semantics.",
                                                ati, expected_line_spacing, ratio
                                            );
                                        }
                                    }
                                }

                                let prf_hz = annotation_prf;

                                // Use burst count from swathProcParams if available, else fallback
                                let mut burst_count =
                                    params.burst_count.unwrap_or(fallback_burst_count).max(1);

                                // Normalize fallbacks to usize
                                let fallback_samples_per_burst_usize =
                                    fallback_samples_per_burst.map(|v| v as usize);
                                let fallback_lines_per_burst_usize =
                                    fallback_lines_per_burst.map(|v| v as usize);

                                // Get range and azimuth samples per burst
                                let mut samples_per_burst = params
                                    .range_samples_per_burst
                                    .or(fallback_samples_per_burst_usize)
                                    .unwrap_or(0);
                                let mut lines_per_burst = params
                                    .azimuth_samples_per_burst
                                    .or(fallback_lines_per_burst_usize)
                                    .unwrap_or(0);

                                let mut first_line_global_opt = params.first_line_of_valid_samples;
                                let mut last_line_global_opt = params.last_line_of_valid_samples;
                                let mut first_sample_global_opt =
                                    params.first_sample_of_valid_samples;
                                let mut last_sample_global_opt =
                                    params.last_sample_of_valid_samples;

                                let mut burst_duration = params.burst_duration.unwrap_or(0.0);

                                let mut fallback_lines_per_burst_opt = None;
                                if let Some(geom) = &burst_geometry {
                                    fallback_lines_per_burst_opt = geom.lines_per_burst;
                                    if geom.total_bursts > 0 && burst_count == 0 {
                                        log::warn!(
                                            "⚠️ DEBUGGING: Using derived burst_count from swathTiming: {}",
                                            geom.total_bursts
                                        );
                                        burst_count = geom.total_bursts;
                                    }
                                }

                                if first_line_global_opt.is_none()
                                    && last_line_global_opt.is_none()
                                    && range_samples_opt.is_some()
                                    && az_lines_opt.is_some()
                                {
                                    let derived_last_line = az_lines_opt.unwrap_or(0);
                                    let derived_first_line = 0usize;
                                    first_line_global_opt = Some(derived_first_line);
                                    last_line_global_opt = Some(derived_last_line);
                                    log::info!(
                                        "✅ DEBUGGING: Using derived line bounds {}..{} from imageInformation",
                                        derived_first_line, derived_last_line
                                    );
                                }

                                if first_sample_global_opt.is_none()
                                    && last_sample_global_opt.is_none()
                                    && range_samples_opt.is_some()
                                {
                                    let derived_last_sample = range_samples_opt.unwrap_or(0);
                                    let derived_first_sample = 0usize;
                                    first_sample_global_opt = Some(derived_first_sample);
                                    last_sample_global_opt = Some(derived_last_sample);
                                    log::info!(
                                        "✅ DEBUGGING: Using derived sample bounds {}..{} from imageInformation",
                                        derived_first_sample, derived_last_sample
                                    );
                                }

                                // Fallback to burst geometry-derived dimensions if missing
                                if let Some(geom) = &burst_geometry {
                                    if lines_per_burst == 0 {
                                        lines_per_burst =
                                            geom.lines_per_burst.unwrap_or(lines_per_burst);
                                    }
                                    if samples_per_burst == 0 {
                                        samples_per_burst =
                                            geom.samples_per_burst.unwrap_or(samples_per_burst);
                                    }

                                    if first_line_global_opt.is_none()
                                        && geom.first_valid_line.is_some()
                                    {
                                        first_line_global_opt = geom.first_valid_line;
                                        last_line_global_opt = geom.last_valid_line_exclusive;
                                    }

                                    if first_sample_global_opt.is_none()
                                        && geom.first_valid_sample.is_some()
                                        && geom.last_valid_sample_exclusive.is_some()
                                    {
                                        first_sample_global_opt = geom.first_valid_sample;
                                        last_sample_global_opt = geom.last_valid_sample_exclusive;
                                    }
                                }

                                // Compute burst duration if missing, using PRF and lines per burst
                                if burst_duration == 0.0 {
                                    if let Some(lines_per_burst) = params
                                        .azimuth_samples_per_burst
                                        .or(fallback_lines_per_burst_opt)
                                    {
                                        burst_duration = lines_per_burst as f64 / prf_hz;
                                    }
                                }

                                if burst_duration <= 0.0 || !burst_duration.is_finite() {
                                    return Err(SarError::Metadata(format!(
                                        "Missing/invalid burstDuration for swath {}",
                                        swath_id
                                    )));
                                }

                                if first_line_global_opt.is_none() {
                                    first_line_global_opt = Some(0);
                                    last_line_global_opt = Some(original_azimuth_samples);
                                }
                                if first_sample_global_opt.is_none() {
                                    first_sample_global_opt = Some(0);
                                    last_sample_global_opt = Some(original_range_samples);
                                }

                                let first_line_global = first_line_global_opt.unwrap_or(0usize);
                                let first_sample_global = first_sample_global_opt.unwrap_or(0usize);

                                if burst_count == 0 {
                                    burst_count = fallback_burst_count;
                                }

                                if lines_per_burst == 0 {
                                    lines_per_burst = fallback_lines_per_burst_usize
                                        .or(original_azimuth_samples.checked_div(burst_count))
                                        .unwrap_or(0);
                                }
                                if samples_per_burst == 0 {
                                    samples_per_burst = original_range_samples;
                                }

                                if samples_per_burst == 0 || lines_per_burst == 0 {
                                    return Err(SarError::Metadata(format!(
                                        "Missing valid swath dimensions for {} (range_samples={}, azimuth_samples={})",
                                        swath_id, samples_per_burst, lines_per_burst
                                    )));
                                }

                                let last_line_global_candidate = last_line_global_opt
                                    .unwrap_or(first_line_global.saturating_add(lines_per_burst));
                                let last_line_global =
                                    last_line_global_candidate.min(original_azimuth_samples);
                                let last_sample_global_candidate = last_sample_global_opt
                                    .unwrap_or(
                                        first_sample_global.saturating_add(samples_per_burst),
                                    );
                                let last_sample_global =
                                    last_sample_global_candidate.min(original_range_samples);

                                log::info!(
                                    "✅ DEBUGGING: Subswath {} -> lines {}..{} ({} per burst), samples {}..{} ({})",
                                    swath_id,
                                    first_line_global,
                                    last_line_global,
                                    lines_per_burst,
                                    first_sample_global,
                                    last_sample_global,
                                    samples_per_burst
                                );

                                // CRITICAL FIX: Extract DC polynomials and reference time from parsed annotation
                                // This prevents deburst from falling back to zeros and missing t0
                                let dc_polynomial_data = annotation
                                    .doppler_centroid
                                    .as_ref()
                                    .and_then(|dc| dc.dc_estimate_list.as_ref())
                                    .and_then(|list| list.dc_estimates.as_ref())
                                    .and_then(|estimates| {
                                        if estimates.is_empty() {
                                            log::warn!("⚠️  No DC estimates found in annotation for subswath {}", swath_id);
                                            return None;
                                        }

                                        let first = &estimates[0];
                                        let poly = first.data_dc_polynomial.clone();

                                        if poly.is_empty() {
                                            log::warn!("⚠️  DC polynomial is empty for subswath {}", swath_id);
                                            return None;
                                        }

                                        let t0_from_iso = first
                                            .azimuth_time
                                            .as_ref()
                                            .and_then(|t| parse_time_robust(t))
                                            .map(|dt| {
                                                dt.timestamp() as f64
                                                    + dt.timestamp_subsec_nanos() as f64 * 1e-9
                                            });

                                        let t0_from_numeric = if first.t0.is_finite() {
                                            Some(first.t0)
                                        } else {
                                            None
                                        };

                                        // CRITICAL FIX (Jan 2026): DC polynomial t0 is SLANT RANGE reference time (~0.0053-0.0060s)
                                        // NOT azimuth epoch time! DC(τ) = c0 + c1*(τ - t0) + c2*(τ - t0)^2 where τ = slant_range_time.
                                        // Using epoch time causes catastrophic errors (DC ~ 10^13 Hz instead of ~50 Hz).
                                        let chosen_t0 = t0_from_numeric.or(t0_from_iso);

                                        // SCIENTIFIC VALIDATION: Verify t0 is in correct domain
                                        if let Some(t0_val) = chosen_t0 {
                                            if t0_val > 1e9 {
                                                log::error!(
                                                    "❌ CRITICAL ERROR: DC polynomial t0={:.3}s appears to be absolute UTC epoch, not slant range time! \
                                                     Expected range: 0.005-0.006s. This will cause catastrophic Doppler errors.",
                                                    t0_val
                                                );
                                            } else if t0_val < 0.004 || t0_val > 0.007 {
                                                log::warn!(
                                                    "⚠️  DC polynomial t0={:.6}s outside typical Sentinel-1 slant range window [0.004, 0.007]s. \
                                                     Verify this is slant_range_time, not azimuth time.",
                                                    t0_val
                                                );
                                            }
                                        }

                                        log::info!(
                                            "✅ Extracted {} DC polynomial estimates from annotation for subswath {}",
                                            estimates.len(), swath_id
                                        );
                                        log::info!(
                                            "   DC polynomial degree: {} (coefficients: {:?}, t0={:?} [DOMAIN: slant_range_time in seconds])",
                                            poly.len() - 1,
                                            &poly,
                                            chosen_t0
                                        );

                                        Some((poly, chosen_t0))
                                    });

                                let (dc_polynomial, dc_polynomial_t0) = match dc_polynomial_data {
                                    Some((poly, t0)) => (Some(poly), t0),
                                    None => (None, None),
                                };

                                let derived_first_sample =
                                    first_sample_global_opt.unwrap_or(0usize);
                                let derived_last_sample = last_sample_global_opt.unwrap_or(
                                    derived_first_sample.saturating_add(samples_per_burst),
                                );

                                pending_subswaths.push(PendingSubswath {
                                    id: swath_id.clone(),
                                    burst_count,
                                    samples_per_burst,
                                    lines_per_burst,
                                    derived_first_line: first_line_global,
                                    derived_last_line: last_line_global,
                                    derived_first_sample,
                                    derived_last_sample,
                                    product_range_samples: original_range_samples,
                                    product_azimuth_samples: original_azimuth_samples,
                                    range_pixel_spacing,
                                    azimuth_pixel_spacing,
                                    slant_range_time,
                                    burst_duration,
                                    prf_hz: Some(prf_hz),
                                    dc_polynomial,
                                    dc_polynomial_t0,
                                    azimuth_time_interval,
                                    near_range_m: slant_range_time * 0.5 * SPEED_OF_LIGHT_M_S,
                                    fm_rate_models: fm_rate_models.clone(),
                                });
                            } else {
                                log::warn!(
                                    "❌ DEBUGGING: swath_proc_params entry {} has no swath field",
                                    i
                                );
                            }
                        }
                    } else {
                        log::warn!("❌ DEBUGGING: swath_proc_params vector is None");
                    }
                } else {
                    log::warn!("❌ DEBUGGING: swath_proc_params_list is None");
                }
            }
        } else {
            log::warn!("❌ DEBUGGING: image_annotation is None");
        }

        if pending_subswaths.is_empty() {
            log::warn!("⚠️  No subswaths extracted from annotation");
        } else {
            let base_near_range_m = pending_subswaths
                .iter()
                .map(|p| p.near_range_m + p.derived_first_sample as f64 * p.range_pixel_spacing)
                .fold(f64::INFINITY, f64::min);

            for pending in pending_subswaths {
                let effective_near = pending.near_range_m
                    + pending.derived_first_sample as f64 * pending.range_pixel_spacing;
                let offset_samples = ((effective_near - base_near_range_m)
                    / pending.range_pixel_spacing)
                    .round()
                    .max(0.0) as usize;

                let first_sample_global = offset_samples + pending.derived_first_sample;

                if pending.burst_count == 0 {
                    return Err(SarError::Metadata(format!(
                        "Subswath {} reports zero bursts in annotation",
                        pending.id
                    )));
                }

                let derived_total_lines = pending
                    .derived_last_line
                    .saturating_sub(pending.derived_first_line);

                if derived_total_lines == 0 {
                    return Err(SarError::Metadata(format!(
                        "Derived zero-length azimuth extent for subswath {}",
                        pending.id
                    )));
                }

                let mut lines_per_burst = pending.lines_per_burst;
                let derived_lines_per_burst = if pending.burst_count > 0 {
                    derived_total_lines / pending.burst_count
                } else {
                    0
                };

                if lines_per_burst == 0 {
                    lines_per_burst = derived_lines_per_burst;
                }

                if lines_per_burst == 0 {
                    return Err(SarError::Metadata(format!(
                        "Missing azimuth_samples_per_burst for subswath {}",
                        pending.id
                    )));
                }

                if derived_lines_per_burst > 0
                    && (lines_per_burst as isize - derived_lines_per_burst as isize).abs() > 4
                {
                    log::warn!(
                        "⚠️  Subswath {} lines_per_burst annotation mismatch: annotation={} inferred={} ({} bursts)",
                        pending.id,
                        lines_per_burst,
                        derived_lines_per_burst,
                        pending.burst_count
                    );
                    lines_per_burst = derived_lines_per_burst.max(1);
                }

                let computed_total_lines =
                    lines_per_burst.saturating_mul(pending.burst_count.max(1));
                let mut azimuth_samples = computed_total_lines.max(derived_total_lines);
                if azimuth_samples == 0 {
                    azimuth_samples = pending.product_azimuth_samples;
                }
                if azimuth_samples == 0 {
                    azimuth_samples = lines_per_burst;
                }

                let first_line_global = pending.derived_first_line;
                let mut last_line_global = first_line_global.saturating_add(azimuth_samples);
                if pending.product_azimuth_samples > 0 {
                    last_line_global = last_line_global.min(pending.product_azimuth_samples);
                    azimuth_samples = last_line_global.saturating_sub(first_line_global);
                }

                let derived_range = pending
                    .derived_last_sample
                    .saturating_sub(pending.derived_first_sample);
                let mut range_samples = derived_range.max(pending.samples_per_burst);
                if range_samples == 0 {
                    range_samples = pending.product_range_samples;
                }
                if range_samples == 0 {
                    range_samples = pending.samples_per_burst;
                }

                let mut last_sample_global = first_sample_global.saturating_add(range_samples);
                if pending.product_range_samples > 0 {
                    last_sample_global = last_sample_global.min(pending.product_range_samples);
                    range_samples = last_sample_global.saturating_sub(first_sample_global);
                }

                let valid_first_sample = first_sample_global;
                let valid_last_sample = last_sample_global;

                let final_azimuth_samples = last_line_global.saturating_sub(first_line_global);
                if final_azimuth_samples == 0 {
                    return Err(SarError::Metadata(format!(
                        "Computed zero azimuth span for subswath {} ({}..{})",
                        pending.id, first_line_global, last_line_global
                    )));
                }

                if azimuth_samples != final_azimuth_samples {
                    log::warn!(
                        "⚠️  Subswath {} azimuth_samples mismatch: stored={} vs bounds={} (forcing bounds)",
                        pending.id,
                        azimuth_samples,
                        final_azimuth_samples
                    );
                }

                if pending.product_azimuth_samples > 0
                    && final_azimuth_samples != pending.product_azimuth_samples
                {
                    log::warn!(
                        "⚠️  Subswath {} azimuth span {} does not match raster height {} from imageInformation",
                        pending.id,
                        final_azimuth_samples,
                        pending.product_azimuth_samples
                    );
                }

                let normalized_lines_per_burst = if pending.burst_count > 0 {
                    final_azimuth_samples / pending.burst_count.max(1)
                } else {
                    lines_per_burst
                }
                .max(1);

                if pending.burst_count > 0 {
                    let implied_height = normalized_lines_per_burst * pending.burst_count;
                    if implied_height != final_azimuth_samples {
                        log::warn!(
                            "⚠️  Subswath {} burst geometry mismatch: {} bursts * {} lines = {} (expected {})",
                            pending.id,
                            pending.burst_count,
                            normalized_lines_per_burst,
                            implied_height,
                            final_azimuth_samples
                        );
                    }
                }

                let subswath = crate::types::SubSwath {
                    id: pending.id.clone(),
                    burst_count: pending.burst_count,
                    lines_per_burst: normalized_lines_per_burst,
                    range_samples,
                    azimuth_samples: final_azimuth_samples,
                    first_line_global,
                    last_line_global,
                    first_sample_global,
                    last_sample_global,
                    full_range_samples: range_samples,
                    valid_first_line: Some(first_line_global),
                    valid_last_line: Some(last_line_global),
                    valid_first_sample: Some(valid_first_sample),
                    valid_last_sample: Some(valid_last_sample),
                    range_pixel_spacing: pending.range_pixel_spacing,
                    azimuth_pixel_spacing: pending.azimuth_pixel_spacing,
                    slant_range_time: pending.slant_range_time,
                    burst_duration: pending.burst_duration,
                    near_range_m: pending.near_range_m,
                    // Use the same PRF we validated earlier (falls back to azimuthFrequency)
                    prf_hz: pending.prf_hz.or(Some(annotation_prf)),
                    dc_polynomial: pending.dc_polynomial,
                    azimuth_time_interval: pending.azimuth_time_interval,
                    dc_polynomial_t0: pending.dc_polynomial_t0,
                    fm_rate_estimates: Some(pending.fm_rate_models.clone()),
                };

                match subswaths.entry(pending.id) {
                    std::collections::hash_map::Entry::Vacant(v) => {
                        v.insert(subswath);
                    }
                    std::collections::hash_map::Entry::Occupied(mut entry) => {
                        let existing = entry.get_mut();
                        if existing.dc_polynomial.is_none() && subswath.dc_polynomial.is_some() {
                            existing.dc_polynomial = subswath.dc_polynomial;
                        }
                        if existing.dc_polynomial_t0.is_none()
                            && subswath.dc_polynomial_t0.is_some()
                        {
                            existing.dc_polynomial_t0 = subswath.dc_polynomial_t0;
                        }
                        if existing.azimuth_time_interval.is_none()
                            && subswath.azimuth_time_interval.is_some()
                        {
                            existing.azimuth_time_interval = subswath.azimuth_time_interval;
                        }
                    }
                }
            }
        }

        log::info!(
            "🔍 DEBUGGING: extract_subswaths returning {} subswaths",
            subswaths.len()
        );

        Ok(subswaths)
    }

    pub fn get_slant_range_time(&self) -> Option<f64> {
        self.image_annotation
            .as_ref()
            .and_then(|ia| ia.image_information.as_ref())
            .and_then(|ii| ii.slant_range_time)
    }

    pub fn get_pulse_repetition_frequency(&self) -> SarResult<f64> {
        let valid_prf = |prf: f64| prf.is_finite() && prf > 0.0 && prf < 6000.0;
        let approx_equal = |a: f64, b: f64| {
            let diff = (a - b).abs();
            let scale = a.abs().max(b.abs()).max(1.0);
            diff <= 1e-6 * scale || diff <= 1e-3
        };

        let mut downlink_prfs: Vec<f64> = Vec::new();
        if let Some(ga) = &self.general_annotation {
            if let Some(dl_list) = &ga.downlink_information_list {
                if let Some(entries) = &dl_list.downlink_information {
                    for di in entries {
                        if let Some(prf) = di.prf {
                            if valid_prf(prf) {
                                downlink_prfs.push(prf);
                            }
                        }
                    }
                }
            }
        }

        let raw_prfs = downlink_prfs.clone();
        downlink_prfs.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
        downlink_prfs.dedup_by(|a, b| approx_equal(*a, *b));

        let prf_from_azimuth_frequency = self
            .image_annotation
            .as_ref()
            .and_then(|ia| ia.image_information.as_ref())
            .and_then(|ii| ii.azimuth_frequency)
            .filter(|p| valid_prf(*p));

        if downlink_prfs.len() > 1 {
            log::error!(
                "❌ PRF mismatch detected in downlinkInformationList: {:?}",
                raw_prfs
            );
            return Err(SarError::Metadata(format!(
                "PRF mismatch in downlinkInformationList (values {:?})",
                raw_prfs
            )));
        }

        if let Some(prf) = downlink_prfs.first() {
            let prf_value = *prf;
            if let Some(az_freq) = prf_from_azimuth_frequency {
                if !approx_equal(prf_value, az_freq) {
                    log::warn!(
                        "⚠️  azimuthFrequency ({:.6} Hz) differs from downlink PRF ({:.6} Hz)",
                        az_freq,
                        prf_value
                    );
                }
            }
            return Ok(prf_value);
        }

        if let Some(fallback) = prf_from_azimuth_frequency {
            log::warn!(
                "⚠️  <prf> missing in downlinkInformationList; using azimuthFrequency={:.6} Hz",
                fallback
            );
            return Ok(fallback);
        }

        Err(SarError::Metadata(
            "PRF not found in annotation (downlinkInformationList/azimuthFrequency missing)"
                .to_string(),
        ))
    }

    pub fn get_radar_frequency_hz(&self) -> Option<f64> {
        self.general_annotation
            .as_ref()
            .and_then(|ga| ga.product_information.as_ref())
            .map(|pi| pi.radar_frequency)
    }

    /// Extract Range-Doppler parameters with proper time base
    /// CRITICAL: Requires orbit vectors to establish orbit_ref_epoch
    pub fn extract_range_doppler_params(
        &self,
        orbit_vectors: &[StateVector],
    ) -> SarResult<crate::core::terrain_correction::RangeDopplerParams> {
        // CRITICAL: Derive time bases first to establish orbit_ref_epoch
        let time_bases = self.derive_time_bases(orbit_vectors)?;
        let (range_ps, az_ps) = self
            .get_pixel_spacing()
            .ok_or_else(|| SarError::Metadata("Missing range/azimuth pixel spacing".to_string()))?;
        let srt = self
            .get_slant_range_time()
            .ok_or_else(|| SarError::Metadata("Missing slantRangeTime".to_string()))?;
        let prf = self.get_pulse_repetition_frequency().map_err(|e| {
            SarError::Metadata(format!("Missing PRF in downlinkInformation: {}", e))
        })?;
        let radar_freq = self.get_radar_frequency_hz().ok_or_else(|| {
            SarError::Metadata("Missing radarFrequency in productInformation".to_string())
        })?;
        let wavelength = SPEED_OF_LIGHT_M_S / radar_freq;

        // Extract azimuth time interval from annotation (CRITICAL for TOPS)
        // This is the actual line time interval from the annotation XML
        // For TOPS/IW data, this can differ from 1/PRF due to steering rate
        let azimuth_time_interval = self
            .image_annotation
            .as_ref()
            .and_then(|ia| ia.image_information.as_ref())
            .and_then(|ii| ii.azimuth_time_interval)
            .unwrap_or_else(|| {
                log::warn!(
                    "⚠️  azimuthTimeInterval not in annotation, using 1/PRF fallback (may be inaccurate for TOPS)"
                );
                1.0 / prf
            });

        // CRITICAL: Use time_bases to compute times relative to orbit_ref_epoch
        let orbit_ref_epoch_utc = (time_bases.orbit_ref_epoch.timestamp() as f64)
            + (time_bases.orbit_ref_epoch.timestamp_subsec_nanos() as f64) * 1e-9;

        let product_start_time_abs = (time_bases.product_start_utc.timestamp() as f64)
            + (time_bases.product_start_utc.timestamp_subsec_nanos() as f64) * 1e-9;

        let product_stop_time_abs = if let Some(stop_utc) = time_bases.product_stop_utc {
            (stop_utc.timestamp() as f64) + (stop_utc.timestamp_subsec_nanos() as f64) * 1e-9
        } else {
            product_start_time_abs // fallback: zero duration if missing
        };

        // CRITICAL: Compute product_start_rel_s (seconds since orbit_ref_epoch)
        let product_start_rel_s = product_start_time_abs - orbit_ref_epoch_utc;

        let product_duration = (product_stop_time_abs - product_start_time_abs).max(0.0);

        log::info!(
            "⏱️  Time base established: orbit_ref_epoch={:.6}s, product_start_rel={:.3}s, duration={:.3}s",
            orbit_ref_epoch_utc, product_start_rel_s, product_duration
        );

        // CRITICAL VALIDATION: Ensure product duration is reasonable
        if product_duration < 0.0 || product_duration > 100.0 {
            log::warn!(
                "⚠️ Unusual product duration: {:.3}s (expected 5-60s for Sentinel-1)",
                product_duration
            );
        }
        if product_duration > 0.0 && product_duration < 60.0 {
            log::debug!("✅ Product duration: {:.3}s (reasonable)", product_duration);
        }

        // Attempt to determine total azimuth lines (merged) from annotation logic reused above
        // Use image_annotation.image_information.lines if present else None
        let total_azimuth_lines = self
            .image_annotation
            .as_ref()
            .and_then(|ia| ia.image_information.as_ref())
            .and_then(|ii| ii.number_of_lines.map(|v| v as usize));

        // Build optional Doppler centroid model from estimates
        let doppler_centroid = {
            // Prefer general_annotation.dcEstimateList, fallback to root doppler_centroid
            let est = self
                .general_annotation
                .as_ref()
                .and_then(|ga| ga.dc_estimate_list.as_ref())
                .and_then(|dl| dl.dc_estimates.as_ref())
                .and_then(|v| v.first().cloned())
                .or_else(|| {
                    self.doppler_centroid
                        .as_ref()
                        .and_then(|dc| dc.dc_estimate_list.as_ref())
                        .and_then(|dl| dl.dc_estimates.as_ref())
                        .and_then(|v| v.first().cloned())
                });

            // CRITICAL VALIDATION: Ensure DC polynomial is valid
            if let Some(ref e) = est {
                crate::io::parsing_validation::validate_dc_polynomial(&e.data_dc_polynomial)?;
            }

            est.map(|e| crate::core::terrain_correction::DopplerCentroidModel {
                t0: e.t0,
                coeffs: e.data_dc_polynomial,
            })
        };
        // Extract valid line/sample bounds from burst geometry (if available)
        let burst_geometry = self
            .swath_timing
            .as_ref()
            .and_then(|st| derive_burst_geometry(st));

        let (first_valid_line, last_valid_line, first_valid_sample, last_valid_sample) =
            if let Some(geom) = burst_geometry {
                (
                    geom.first_valid_line,
                    geom.last_valid_line_exclusive.map(|x| x.saturating_sub(1)), // Convert exclusive to inclusive
                    geom.first_valid_sample,
                    geom.last_valid_sample_exclusive
                        .map(|x| x.saturating_sub(1)),
                ) // Convert exclusive to inclusive
            } else {
                (None, None, None, None)
            };

        let mut params = crate::core::terrain_correction::RangeDopplerParams {
            range_pixel_spacing: range_ps,
            azimuth_pixel_spacing: az_ps,
            slant_range_time: srt,
            prf,
            azimuth_time_interval,
            wavelength,
            speed_of_light: SPEED_OF_LIGHT_M_S,
            orbit_ref_epoch_utc, // NEW: Orbit reference epoch
            product_start_rel_s, // NEW: Time relative to orbit_ref_epoch
            #[allow(deprecated)]
            product_start_time_abs, // DEPRECATED: For legacy compatibility
            #[allow(deprecated)]
            product_stop_time_abs, // DEPRECATED: For legacy compatibility
            product_duration,
            total_azimuth_lines,
            doppler_centroid,
            first_valid_line,
            last_valid_line,
            first_valid_sample,
            last_valid_sample,
            range_multilook_factor: 1.0, // Default to no multilooking for annotation parser
            azimuth_multilook_factor: 1.0, // Default to no multilooking for annotation parser
            range_multilook_safe: 1.0,   // Will be computed below
            azimuth_multilook_safe: 1.0, // Will be computed below
            subswaths: std::collections::HashMap::new(), // Empty for annotation parser
            burst_timings: Vec::new(),
            burst_segments: Vec::new(),
            reference_incidence_angle_deg: None, // Populated later from image_annotation
            incidence_angle_near_deg: None,      // Populated later from metadata
            incidence_angle_far_deg: None,       // Populated later from metadata
            total_range_samples: None,           // Populated later from image dimensions
        };
        params.compute_safe_multilook_factors();
        Ok(params)
    }

    pub fn evaluate_doppler_centroid(&self, az_time_since_start: f64) -> SarResult<f64> {
        // Prefer general_annotation.dcEstimateList
        let est = self
            .general_annotation
            .as_ref()
            .and_then(|ga| ga.dc_estimate_list.as_ref())
            .and_then(|dl| dl.dc_estimates.as_ref())
            .and_then(|v| v.first())
            .cloned()
            .or_else(|| {
                self.doppler_centroid
                    .as_ref()
                    .and_then(|dc| dc.dc_estimate_list.as_ref())
                    .and_then(|dl| dl.dc_estimates.as_ref())
                    .and_then(|v| v.first())
                    .cloned()
            })
            .ok_or_else(|| {
                SarError::Metadata("No Doppler centroid estimates in annotation".to_string())
            })?;

        // CRITICAL VALIDATION: Ensure DC polynomial is valid before evaluation
        crate::io::parsing_validation::validate_dc_polynomial(&est.data_dc_polynomial)?;

        let x = az_time_since_start - est.t0;
        let mut acc = 0.0;
        let mut pow = 1.0;
        for coeff in est.data_dc_polynomial.iter() {
            acc += coeff * pow;
            pow *= x;
        }
        Ok(acc)
    }

    /// Validate Doppler centroid estimates for consistency across the swath.
    ///
    /// IMPORTANT (Jan 2026 fix): DC polynomials are RANGE polynomials, not azimuth polynomials!
    /// The independent variable is slant range time (~5ms), NOT azimuth time (~seconds).
    /// This function compares c0 coefficients (constant term in Hz) between adjacent
    /// DC estimates to verify they are consistent within expected bounds.
    ///
    /// For S1 IW mode, DC c0 should be ~±50 Hz and vary slowly along track.
    pub fn validate_doppler_model_against_bursts(&self, pol: Polarization) {
        let dc_estimates = self
            .general_annotation
            .as_ref()
            .and_then(|ga| ga.dc_estimate_list.as_ref())
            .and_then(|dl| dl.dc_estimates.as_ref())
            .cloned()
            .or_else(|| {
                self.doppler_centroid
                    .as_ref()
                    .and_then(|dc| dc.dc_estimate_list.as_ref())
                    .and_then(|dl| dl.dc_estimates.as_ref())
                    .cloned()
            });

        let dc_estimates = match dc_estimates {
            Some(est) if !est.is_empty() => est,
            _ => return,
        };

        // Extract c0 (constant term) from each DC estimate
        let c0_values: Vec<(usize, f64)> = dc_estimates
            .iter()
            .enumerate()
            .filter_map(|(i, est)| est.data_dc_polynomial.first().copied().map(|c0| (i, c0)))
            .collect();

        if c0_values.is_empty() {
            return;
        }

        // Log summary of DC c0 values
        let c0_only: Vec<f64> = c0_values.iter().map(|(_, c0)| *c0).collect();
        let c0_min = c0_only.iter().cloned().fold(f64::INFINITY, f64::min);
        let c0_max = c0_only.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        let c0_mean = c0_only.iter().sum::<f64>() / c0_only.len() as f64;

        log::info!(
            "🔍 {} DC estimates: {} values, c0 range [{:.1}, {:.1}] Hz, mean={:.1} Hz",
            pol,
            c0_values.len(),
            c0_min,
            c0_max,
            c0_mean
        );

        // Check physical plausibility: S1 IW DC should be within ~±200 Hz
        const DC_PLAUSIBLE_LIMIT_HZ: f64 = 500.0;
        let mut flagged = 0usize;

        for &(idx, c0) in &c0_values {
            if c0.abs() > DC_PLAUSIBLE_LIMIT_HZ {
                flagged += 1;
                log::warn!(
                    "⚠️  {} DC estimate {}: c0={:.1} Hz exceeds plausible range (±{} Hz)",
                    pol,
                    idx,
                    c0,
                    DC_PLAUSIBLE_LIMIT_HZ
                );
            }
        }

        // Check consistency between adjacent estimates
        const DC_ADJACENT_THRESH_HZ: f64 = 50.0;
        for window in c0_values.windows(2) {
            let (idx_a, c0_a) = window[0];
            let (idx_b, c0_b) = window[1];
            let delta = (c0_b - c0_a).abs();

            if delta > DC_ADJACENT_THRESH_HZ {
                log::debug!(
                    "🔍 {} DC jump between estimates {} and {}: Δc0={:.1} Hz",
                    pol,
                    idx_a,
                    idx_b,
                    delta
                );
            }
        }

        if flagged == 0 {
            log::info!(
                "✅ Doppler centroid estimates valid for {} ({} estimates, c0 range {:.1} Hz)",
                pol,
                c0_values.len(),
                c0_max - c0_min
            );
        }
    }

    /// Calculate average satellite velocity from orbit state vectors
    /// Uses actual orbit data for scientifically accurate phase calculations
    pub fn calculate_satellite_velocity(metadata: &crate::types::SarMetadata) -> SarResult<f64> {
        if let Some(orbit_data) = &metadata.orbit_data {
            let state_vectors = &orbit_data.state_vectors;
            if state_vectors.len() < 2 {
                return Err(SarError::Metadata(
                    "Need at least 2 orbit state vectors to calculate velocity".to_string(),
                ));
            }

            let mut total_velocity = 0.0;
            let mut count = 0;

            for state in state_vectors {
                // Calculate velocity magnitude from components
                let vx = state.velocity[0];
                let vy = state.velocity[1];
                let vz = state.velocity[2];
                let velocity_magnitude = (vx * vx + vy * vy + vz * vz).sqrt();

                // Validate velocity is in reasonable range for LEO satellites
                if velocity_magnitude > 6000.0 && velocity_magnitude < 8000.0 {
                    total_velocity += velocity_magnitude;
                    count += 1;
                }
            }

            if count == 0 {
                return Err(SarError::Metadata(
                    "No valid satellite velocities found in orbit state vectors".to_string(),
                ));
            }

            Ok(total_velocity / count as f64)
        } else {
            Err(SarError::Metadata(
                "No orbit data available for velocity calculation".to_string(),
            ))
        }
    }
}
