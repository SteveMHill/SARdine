//! Validation logic for scene metadata invariants.
//!
//! The [`check`] function returns all violations found, not just the first.
//! An empty result means the metadata is valid.

use chrono::{DateTime, Utc};

use crate::types::*;

// ─── Error Types ─────────────────────────────────────────────────────

/// Individual validation violation.
#[derive(Debug, Clone, thiserror::Error)]
pub enum ValidationError {
    #[error("product_id must not be empty")]
    EmptyProductId,

    #[error("at least one polarization is required")]
    NoPolarizations,

    #[error("start_time ({start}) must be before stop_time ({stop})")]
    TimingInverted {
        start: DateTime<Utc>,
        stop: DateTime<Utc>,
    },

    #[error("radar_frequency_hz must be finite and positive, got {value}")]
    InvalidRadarFrequency { value: f64 },

    #[error("range_sampling_rate_hz must be finite and positive, got {value}")]
    InvalidRangeSamplingRate { value: f64 },

    #[error("bounding_box.{field}: {detail}")]
    InvalidBoundingBox {
        field: &'static str,
        detail: String,
    },

    #[error("{mode} mode requires at least one sub-swath")]
    NoSubSwaths { mode: AcquisitionMode },

    #[error("sub-swath {id}: {detail}")]
    InvalidSubSwath { id: SubSwathId, detail: String },

    #[error("{mode} mode requires at least one burst entry")]
    NoBursts { mode: AcquisitionMode },

    #[error("burst {subswath}/{index}: {detail}")]
    InvalidBurst {
        subswath: SubSwathId,
        index: usize,
        detail: String,
    },

    #[error("orbit: {detail}")]
    InvalidOrbit { detail: String },
}

/// Collection of all validation errors found.
#[derive(Debug, Clone)]
pub struct ValidationErrors(pub Vec<ValidationError>);

impl std::fmt::Display for ValidationErrors {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{} validation error(s):", self.0.len())?;
        for e in &self.0 {
            write!(f, "\n  - {}", e)?;
        }
        Ok(())
    }
}

impl std::error::Error for ValidationErrors {}

// ─── Main Validation Entry Point ─────────────────────────────────────

/// Check all invariants on scene metadata. Returns a list of violations (empty = valid).
pub fn check(m: &SceneMetadata) -> Vec<ValidationError> {
    let mut errors = Vec::new();

    // Product identity
    if m.product_id.trim().is_empty() {
        errors.push(ValidationError::EmptyProductId);
    }

    // Polarizations
    if m.polarizations.is_empty() {
        errors.push(ValidationError::NoPolarizations);
    }

    // Timing
    if m.start_time >= m.stop_time {
        errors.push(ValidationError::TimingInverted {
            start: m.start_time,
            stop: m.stop_time,
        });
    }

    // Radar parameters
    if !m.radar_frequency_hz.is_finite() || m.radar_frequency_hz <= 0.0 {
        errors.push(ValidationError::InvalidRadarFrequency {
            value: m.radar_frequency_hz,
        });
    }
    if !m.range_sampling_rate_hz.is_finite() || m.range_sampling_rate_hz <= 0.0 {
        errors.push(ValidationError::InvalidRangeSamplingRate {
            value: m.range_sampling_rate_hz,
        });
    }

    // Bounding box
    check_bounding_box(&m.bounding_box, &mut errors);

    // Sub-swaths
    if m.acquisition_mode.is_tops() && m.sub_swaths.is_empty() {
        errors.push(ValidationError::NoSubSwaths {
            mode: m.acquisition_mode,
        });
    }
    for sw in &m.sub_swaths {
        check_subswath(sw, m.acquisition_mode, &mut errors);
    }

    // Bursts
    if m.acquisition_mode.is_tops() && m.bursts.is_empty() {
        errors.push(ValidationError::NoBursts {
            mode: m.acquisition_mode,
        });
    }
    for burst in &m.bursts {
        check_burst(burst, m, &mut errors);
    }

    // Cross-check: burst count per subswath must match SubSwathMetadata.burst_count
    for sw in &m.sub_swaths {
        let actual = m.bursts.iter().filter(|b| b.subswath_id == sw.id).count();
        if actual != sw.burst_count {
            errors.push(ValidationError::InvalidSubSwath {
                id: sw.id,
                detail: format!(
                    "burst_count is {} but found {} burst entries",
                    sw.burst_count, actual
                ),
            });
        }
    }

    // Burst monotonicity: within each subswath, bursts must be time-ordered
    for sw in &m.sub_swaths {
        let sw_bursts: Vec<_> = m.bursts.iter().filter(|b| b.subswath_id == sw.id).collect();
        for w in sw_bursts.windows(2) {
            if w[1].azimuth_time_utc <= w[0].azimuth_time_utc {
                errors.push(ValidationError::InvalidBurst {
                    subswath: sw.id,
                    index: w[1].burst_index,
                    detail: format!(
                        "burst time {} is not after previous burst time {}",
                        w[1].azimuth_time_utc, w[0].azimuth_time_utc
                    ),
                });
                break; // one violation per subswath is enough
            }
        }
    }

    // Orbit
    check_orbit(&m.orbit, &mut errors);

    errors
}

// ─── Component Validators ────────────────────────────────────────────

fn check_bounding_box(bb: &BoundingBox, errors: &mut Vec<ValidationError>) {
    let all_finite = bb.min_lat_deg.is_finite()
        && bb.max_lat_deg.is_finite()
        && bb.min_lon_deg.is_finite()
        && bb.max_lon_deg.is_finite();

    if !all_finite {
        errors.push(ValidationError::InvalidBoundingBox {
            field: "coordinates",
            detail: "all coordinates must be finite".into(),
        });
        return;
    }

    if bb.min_lat_deg >= bb.max_lat_deg {
        errors.push(ValidationError::InvalidBoundingBox {
            field: "latitude",
            detail: format!(
                "min_lat ({}) >= max_lat ({})",
                bb.min_lat_deg, bb.max_lat_deg
            ),
        });
    }
    if bb.min_lon_deg >= bb.max_lon_deg {
        errors.push(ValidationError::InvalidBoundingBox {
            field: "longitude",
            detail: format!(
                "min_lon ({}) >= max_lon ({})",
                bb.min_lon_deg, bb.max_lon_deg
            ),
        });
    }
}

fn check_subswath(
    sw: &SubSwathMetadata,
    mode: AcquisitionMode,
    errors: &mut Vec<ValidationError>,
) {
    // Mode consistency
    if mode.is_tops() && !sw.id.matches_mode(mode) {
        errors.push(ValidationError::InvalidSubSwath {
            id: sw.id,
            detail: format!("sub-swath {} is not valid for {} mode", sw.id, mode),
        });
    }

    // Dimension checks
    if sw.burst_count == 0 {
        errors.push(ValidationError::InvalidSubSwath {
            id: sw.id,
            detail: "burst_count must be > 0".into(),
        });
    }
    if sw.lines_per_burst == 0 {
        errors.push(ValidationError::InvalidSubSwath {
            id: sw.id,
            detail: "lines_per_burst must be > 0".into(),
        });
    }
    if sw.range_samples == 0 {
        errors.push(ValidationError::InvalidSubSwath {
            id: sw.id,
            detail: "range_samples must be > 0".into(),
        });
    }
    if sw.azimuth_samples == 0 {
        errors.push(ValidationError::InvalidSubSwath {
            id: sw.id,
            detail: "azimuth_samples must be > 0".into(),
        });
    }

    // Bound ordering
    if sw.first_line >= sw.last_line {
        errors.push(ValidationError::InvalidSubSwath {
            id: sw.id,
            detail: format!(
                "first_line ({}) must be < last_line ({})",
                sw.first_line, sw.last_line
            ),
        });
    }
    if sw.first_sample >= sw.last_sample {
        errors.push(ValidationError::InvalidSubSwath {
            id: sw.id,
            detail: format!(
                "first_sample ({}) must be < last_sample ({})",
                sw.first_sample, sw.last_sample
            ),
        });
    }

    // Physical parameters: all must be finite and positive
    check_positive_finite(sw.range_pixel_spacing_m, "range_pixel_spacing_m", sw.id, errors);
    check_positive_finite(sw.azimuth_pixel_spacing_m, "azimuth_pixel_spacing_m", sw.id, errors);
    check_positive_finite(sw.slant_range_time_s, "slant_range_time_s", sw.id, errors);
    check_positive_finite(sw.azimuth_time_interval_s, "azimuth_time_interval_s", sw.id, errors);
    check_positive_finite(sw.prf_hz, "prf_hz", sw.id, errors);
    check_positive_finite(sw.burst_cycle_time_s, "burst_cycle_time_s", sw.id, errors);
}

fn check_positive_finite(
    value: f64,
    name: &str,
    id: SubSwathId,
    errors: &mut Vec<ValidationError>,
) {
    if !value.is_finite() || value <= 0.0 {
        errors.push(ValidationError::InvalidSubSwath {
            id,
            detail: format!("{} must be finite and positive, got {}", name, value),
        });
    }
}

fn check_burst(burst: &BurstEntry, scene: &SceneMetadata, errors: &mut Vec<ValidationError>) {
    // Burst azimuth time must be within the scene acquisition window.
    // A small tolerance (1 s) allows for timing quantization at scene edges.
    let tolerance = chrono::Duration::seconds(1);
    if burst.azimuth_time_utc < scene.start_time - tolerance
        || burst.azimuth_time_utc > scene.stop_time + tolerance
    {
        errors.push(ValidationError::InvalidBurst {
            subswath: burst.subswath_id,
            index: burst.burst_index,
            detail: format!(
                "azimuth_time_utc {} is outside scene window [{}, {}]",
                burst.azimuth_time_utc, scene.start_time, scene.stop_time
            ),
        });
    }

    // Bound ordering
    if burst.first_line >= burst.last_line {
        errors.push(ValidationError::InvalidBurst {
            subswath: burst.subswath_id,
            index: burst.burst_index,
            detail: format!(
                "first_line ({}) must be < last_line ({})",
                burst.first_line, burst.last_line
            ),
        });
    }
    if burst.first_valid_sample >= burst.last_valid_sample {
        errors.push(ValidationError::InvalidBurst {
            subswath: burst.subswath_id,
            index: burst.burst_index,
            detail: format!(
                "first_valid_sample ({}) must be < last_valid_sample ({})",
                burst.first_valid_sample, burst.last_valid_sample
            ),
        });
    }

    // Cross-reference: burst subswath must exist in the scene sub_swaths
    if !scene.sub_swaths.iter().any(|sw| sw.id == burst.subswath_id) {
        errors.push(ValidationError::InvalidBurst {
            subswath: burst.subswath_id,
            index: burst.burst_index,
            detail: format!(
                "references sub-swath {} which is not in scene.sub_swaths",
                burst.subswath_id
            ),
        });
    }
}

fn check_orbit(orbit: &OrbitData, errors: &mut Vec<ValidationError>) {
    // Minimum vector count
    if orbit.state_vectors.len() < MIN_ORBIT_VECTORS {
        errors.push(ValidationError::InvalidOrbit {
            detail: format!(
                "need at least {} state vectors, got {}",
                MIN_ORBIT_VECTORS,
                orbit.state_vectors.len()
            ),
        });
    }

    // Strict time monotonicity
    for window in orbit.state_vectors.windows(2) {
        if window[0].time >= window[1].time {
            errors.push(ValidationError::InvalidOrbit {
                detail: format!(
                    "state vector times must be strictly increasing: {} >= {}",
                    window[0].time, window[1].time
                ),
            });
            break; // one violation is enough to report
        }
    }

    // All positions and velocities must be finite
    for (i, sv) in orbit.state_vectors.iter().enumerate() {
        let pos_finite = sv.position_m.iter().all(|v| v.is_finite());
        let vel_finite = sv.velocity_m_s.iter().all(|v| v.is_finite());
        if !pos_finite || !vel_finite {
            errors.push(ValidationError::InvalidOrbit {
                detail: format!("state vector {} has non-finite position or velocity", i),
            });
            break;
        }
    }
}

// ─── Tests ───────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use chrono::TimeZone;

    // ── Test Helpers ─────────────────────────────────────────────────

    fn sample_orbit() -> OrbitData {
        let epoch = Utc.with_ymd_and_hms(2020, 10, 5, 17, 6, 0).unwrap();
        let state_vectors: Vec<StateVector> = (0..17u32)
            .map(|i| StateVector {
                time: epoch + chrono::Duration::seconds(i64::from(i) * 10),
                position_m: [
                    4_000_000.0 + f64::from(i) * 1000.0,
                    5_000_000.0,
                    3_000_000.0,
                ],
                velocity_m_s: [100.0, 200.0, 7500.0],
            })
            .collect();
        OrbitData {
            reference_epoch: epoch,
            state_vectors,
        }
    }

    fn sample_subswaths() -> Vec<SubSwathMetadata> {
        vec![SubSwathMetadata {
            id: SubSwathId::IW1,
            burst_count: 9,
            lines_per_burst: 1526,
            range_samples: 21608,
            azimuth_samples: 13743,
            first_line: 0,
            last_line: 13743,
            first_sample: 0,
            last_sample: 21608,
            range_pixel_spacing_m: 2.329562,
            azimuth_pixel_spacing_m: 13.968,
            slant_range_time_s: 0.005387,
            azimuth_time_interval_s: 0.002063,
            prf_hz: 1717.0,
            burst_cycle_time_s: 2.7578,
            dc_estimates: Vec::new(),
            fm_rates: Vec::new(),
        }]
    }

    fn sample_bursts() -> Vec<BurstEntry> {
        let base = Utc.with_ymd_and_hms(2020, 10, 5, 17, 8, 26).unwrap();
        (0..9usize)
            .map(|i| BurstEntry {
                subswath_id: SubSwathId::IW1,
                burst_index: i,
                azimuth_time_utc: base + chrono::Duration::milliseconds((i as i64) * 2758),
                first_line: i * 1526,
                last_line: (i + 1) * 1526,
                first_valid_sample: 100,
                last_valid_sample: 21600,
                slice_index: 0,
            })
            .collect()
    }

    fn valid_metadata() -> SceneMetadata {
        SceneMetadata {
            product_id:
                "S1A_IW_SLC__1SDV_20201005T170824_20201005T170851_034664_04098A_1E66".into(),
            mission: Mission::S1A,
            acquisition_mode: AcquisitionMode::IW,
            polarizations: vec![Polarization::VV, Polarization::VH],
            start_time: Utc.with_ymd_and_hms(2020, 10, 5, 17, 8, 24).unwrap(),
            stop_time: Utc.with_ymd_and_hms(2020, 10, 5, 17, 8, 51).unwrap(),
            radar_frequency_hz: 5.405e9,
            range_sampling_rate_hz: 64.345e6,
            bounding_box: BoundingBox {
                min_lat_deg: 40.0,
                max_lat_deg: 42.0,
                min_lon_deg: 10.0,
                max_lon_deg: 12.5,
            },
            sub_swaths: sample_subswaths(),
            bursts: sample_bursts(),
            orbit: sample_orbit(),
        }
    }

    fn has_error<F: Fn(&ValidationError) -> bool>(errors: &[ValidationError], pred: F) -> bool {
        errors.iter().any(pred)
    }

    // ── Tests ────────────────────────────────────────────────────────

    #[test]
    fn valid_metadata_passes() {
        let errors = check(&valid_metadata());
        assert!(errors.is_empty(), "expected no errors, got: {:?}", errors);
    }

    #[test]
    fn validated_returns_ok_for_valid() {
        let m = valid_metadata();
        assert!(m.validated().is_ok());
    }

    #[test]
    fn empty_product_id_rejected() {
        let mut m = valid_metadata();
        m.product_id = String::new();
        let errors = check(&m);
        assert!(has_error(&errors, |e| matches!(e, ValidationError::EmptyProductId)));
    }

    #[test]
    fn empty_polarizations_rejected() {
        let mut m = valid_metadata();
        m.polarizations.clear();
        let errors = check(&m);
        assert!(has_error(&errors, |e| matches!(e, ValidationError::NoPolarizations)));
    }

    #[test]
    fn inverted_timing_rejected() {
        let mut m = valid_metadata();
        std::mem::swap(&mut m.start_time, &mut m.stop_time);
        let errors = check(&m);
        assert!(has_error(&errors, |e| matches!(e, ValidationError::TimingInverted { .. })));
    }

    #[test]
    fn nan_radar_frequency_rejected() {
        let mut m = valid_metadata();
        m.radar_frequency_hz = f64::NAN;
        let errors = check(&m);
        assert!(has_error(&errors, |e| matches!(e, ValidationError::InvalidRadarFrequency { .. })));
    }

    #[test]
    fn zero_radar_frequency_rejected() {
        let mut m = valid_metadata();
        m.radar_frequency_hz = 0.0;
        let errors = check(&m);
        assert!(has_error(&errors, |e| matches!(e, ValidationError::InvalidRadarFrequency { .. })));
    }

    #[test]
    fn no_subswaths_rejected_for_iw() {
        let mut m = valid_metadata();
        m.sub_swaths.clear();
        m.bursts.clear(); // also clear bursts to avoid cross-reference errors
        let errors = check(&m);
        assert!(has_error(&errors, |e| matches!(e, ValidationError::NoSubSwaths { .. })));
    }

    #[test]
    fn no_bursts_rejected_for_iw() {
        let mut m = valid_metadata();
        m.bursts.clear();
        let errors = check(&m);
        assert!(has_error(&errors, |e| matches!(e, ValidationError::NoBursts { .. })));
    }

    #[test]
    fn burst_outside_scene_window_rejected() {
        let mut m = valid_metadata();
        m.bursts[0].azimuth_time_utc = Utc.with_ymd_and_hms(2000, 1, 1, 0, 0, 0).unwrap();
        let errors = check(&m);
        assert!(has_error(&errors, |e| {
            matches!(e, ValidationError::InvalidBurst { detail, .. } if detail.contains("outside scene window"))
        }));
    }

    #[test]
    fn inverted_burst_line_bounds_rejected() {
        let mut m = valid_metadata();
        m.bursts[0].first_line = 1000;
        m.bursts[0].last_line = 500;
        let errors = check(&m);
        assert!(has_error(&errors, |e| {
            matches!(e, ValidationError::InvalidBurst { detail, .. } if detail.contains("first_line"))
        }));
    }

    #[test]
    fn inverted_subswath_bounds_rejected() {
        let mut m = valid_metadata();
        m.sub_swaths[0].first_line = 20000;
        m.sub_swaths[0].last_line = 100;
        let errors = check(&m);
        assert!(has_error(&errors, |e| {
            matches!(e, ValidationError::InvalidSubSwath { detail, .. } if detail.contains("first_line"))
        }));
    }

    #[test]
    fn zero_pixel_spacing_rejected() {
        let mut m = valid_metadata();
        m.sub_swaths[0].range_pixel_spacing_m = 0.0;
        let errors = check(&m);
        assert!(has_error(&errors, |e| {
            matches!(e, ValidationError::InvalidSubSwath { detail, .. } if detail.contains("range_pixel_spacing_m"))
        }));
    }

    #[test]
    fn insufficient_orbit_vectors_rejected() {
        let mut m = valid_metadata();
        m.orbit.state_vectors.truncate(3);
        let errors = check(&m);
        assert!(has_error(&errors, |e| {
            matches!(e, ValidationError::InvalidOrbit { detail } if detail.contains("at least"))
        }));
    }

    #[test]
    fn nonmonotonic_orbit_times_rejected() {
        let mut m = valid_metadata();
        let dup_time = m.orbit.state_vectors[0].time;
        m.orbit.state_vectors[1].time = dup_time;
        let errors = check(&m);
        assert!(has_error(&errors, |e| {
            matches!(e, ValidationError::InvalidOrbit { detail } if detail.contains("strictly increasing"))
        }));
    }

    #[test]
    fn inverted_bbox_latitude_rejected() {
        let mut m = valid_metadata();
        m.bounding_box.min_lat_deg = 50.0;
        m.bounding_box.max_lat_deg = 40.0;
        let errors = check(&m);
        assert!(has_error(&errors, |e| {
            matches!(e, ValidationError::InvalidBoundingBox { field: "latitude", .. })
        }));
    }

    #[test]
    fn burst_referencing_missing_subswath_rejected() {
        let mut m = valid_metadata();
        m.bursts[0].subswath_id = SubSwathId::IW3; // IW3 not in sub_swaths
        let errors = check(&m);
        assert!(has_error(&errors, |e| {
            matches!(e, ValidationError::InvalidBurst { detail, .. } if detail.contains("not in scene"))
        }));
    }

    #[test]
    fn ew_subswath_in_iw_mode_rejected() {
        let mut m = valid_metadata();
        m.sub_swaths[0].id = SubSwathId::EW1;
        let errors = check(&m);
        assert!(has_error(&errors, |e| {
            matches!(e, ValidationError::InvalidSubSwath { detail, .. } if detail.contains("not valid for IW"))
        }));
    }

    #[test]
    fn wavelength_derived_correctly() {
        let m = valid_metadata();
        let expected = SPEED_OF_LIGHT_M_S / 5.405e9;
        assert!((m.wavelength_m() - expected).abs() < 1e-10);
    }

    #[test]
    fn near_range_derived_correctly() {
        let m = valid_metadata();
        let expected = 0.005387 * SPEED_OF_LIGHT_M_S / 2.0;
        assert!((m.sub_swaths[0].near_range_m() - expected).abs() < 1e-3);
    }

    #[test]
    fn sm_mode_allows_empty_subswaths() {
        let mut m = valid_metadata();
        m.acquisition_mode = AcquisitionMode::SM;
        m.sub_swaths.clear();
        m.bursts.clear();
        let errors = check(&m);
        // Should not have NoSubSwaths or NoBursts errors
        assert!(!has_error(&errors, |e| matches!(
            e,
            ValidationError::NoSubSwaths { .. } | ValidationError::NoBursts { .. }
        )));
    }

    #[test]
    fn burst_count_mismatch_rejected() {
        let mut m = valid_metadata();
        // Remove one burst so count (8) != sub_swath.burst_count (9)
        m.bursts.pop();
        let errors = check(&m);
        assert!(has_error(&errors, |e| {
            matches!(e, ValidationError::InvalidSubSwath { detail, .. } if detail.contains("burst entries"))
        }));
    }

    #[test]
    fn nonmonotonic_burst_times_rejected() {
        let mut m = valid_metadata();
        // Swap burst 0 and 1 times so they are out of order
        let t0 = m.bursts[0].azimuth_time_utc;
        let t1 = m.bursts[1].azimuth_time_utc;
        m.bursts[0].azimuth_time_utc = t1;
        m.bursts[1].azimuth_time_utc = t0;
        let errors = check(&m);
        assert!(has_error(&errors, |e| {
            matches!(e, ValidationError::InvalidBurst { detail, .. } if detail.contains("not after previous"))
        }));
    }
}
