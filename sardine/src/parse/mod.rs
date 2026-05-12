//! Annotation XML parser for Sentinel-1 SAFE products.
//!
//! Public entry point: [`parse_safe_directory`] reads all annotation XMLs
//! from a `.SAFE` directory and returns a validated [`SceneMetadata`].

mod xml;

use std::collections::BTreeMap;
use std::path::Path;

use chrono::{DateTime, NaiveDateTime, Utc};

use crate::types::*;
use crate::validate::ValidationErrors;
use xml::ProductXml;

// ── Error type ──────────────────────────────────────────────────────

#[derive(Debug, thiserror::Error)]
pub enum ParseError {
    #[error("IO error: {0}")]
    Io(#[from] std::io::Error),

    #[error("XML deserialization failed: {0}")]
    Xml(String),

    #[error("missing required field: {0}")]
    MissingField(&'static str),

    #[error("invalid value for {field}: {detail}")]
    InvalidValue {
        field: &'static str,
        detail: String,
    },

    #[error("no annotation XMLs found in {0}")]
    NoAnnotations(String),

    #[error("polynomial has {found} coefficients, expected 3 for {context}")]
    IncompatiblePolynomial {
        context: &'static str,
        found: usize,
    },

    #[error("validation failed: {0}")]
    Validation(#[from] ValidationErrors),
}

// ── Public API ──────────────────────────────────────────────────────

/// Parse a Sentinel-1 `.SAFE` directory into validated [`SceneMetadata`].
///
/// Reads all annotation XMLs in `<safe_path>/annotation/`, groups by sub-swath,
/// and assembles a validated metadata structure.
pub fn parse_safe_directory(safe_path: &Path) -> Result<SceneMetadata, ParseError> {
    // 1. Discover annotation XMLs
    let annotation_dir = safe_path.join("annotation");
    let mut xml_paths: Vec<_> = std::fs::read_dir(&annotation_dir)
        .map_err(|e| {
            ParseError::Io(std::io::Error::new(
                e.kind(),
                format!("{}: {}", annotation_dir.display(), e),
            ))
        })?
        .filter_map(|entry| entry.ok())
        .filter(|entry| {
            entry
                .file_type()
                .map(|t| t.is_file())
                .unwrap_or(false) // SAFETY-OK: failed file_type → skip the entry; not a numeric path
        })
        .filter(|entry| {
            entry
                .file_name()
                .to_string_lossy()
                .ends_with(".xml")
        })
        .map(|entry| entry.path())
        .collect();

    if xml_paths.is_empty() {
        return Err(ParseError::NoAnnotations(
            annotation_dir.display().to_string(),
        ));
    }

    // Sort for deterministic processing order
    xml_paths.sort();

    // 2. Parse all annotation XMLs
    let mut annotations = Vec::with_capacity(xml_paths.len());
    for path in &xml_paths {
        let content = std::fs::read_to_string(path)?;
        let parsed = parse_annotation_xml(&content)?;
        annotations.push(parsed);
    }

    // 3. Extract SceneMetadata
    extract_scene_metadata(safe_path, &annotations)
}

/// Parse a single annotation XML string into the raw XML structure.
fn parse_annotation_xml(xml_content: &str) -> Result<ProductXml, ParseError> {
    // Strip BOM if present
    let content = xml_content
        .strip_prefix('\u{feff}')
        .unwrap_or(xml_content); // SAFETY-OK: strip_prefix-on-no-match returns input unchanged by design (BOM optional)

    quick_xml::de::from_str::<ProductXml>(content)
        .map_err(|e| ParseError::Xml(e.to_string()))
}

// ── Extraction logic ────────────────────────────────────────────────

fn extract_scene_metadata(
    safe_path: &Path,
    annotations: &[ProductXml],
) -> Result<SceneMetadata, ParseError> {
    let first = annotations
        .first()
        .ok_or(ParseError::NoAnnotations("empty annotations".into()))?;

    // Product ID from directory name
    let product_id = safe_path
        .file_name()
        .and_then(|n| n.to_str())
        .map(|n| n.strip_suffix(".SAFE").unwrap_or(n)) // SAFETY-OK: strip_suffix-on-no-match returns input unchanged by design
        .ok_or(ParseError::MissingField("product_id (SAFE dir name)"))?
        .to_string();

    // Mission
    let mission = parse_mission(&first.ads_header.mission_id)?;

    // Acquisition mode
    let mode = parse_mode(&first.ads_header.mode)?;

    // Radar parameters (identical across all annotations in a scene)
    let radar_frequency_hz = first.general_annotation.product_information.radar_frequency;
    let range_sampling_rate_hz = first.general_annotation.product_information.range_sampling_rate;

    // Collect unique polarizations
    let mut polarizations = Vec::new();
    for ann in annotations {
        let pol = parse_polarization(&ann.ads_header.polarisation)?;
        if !polarizations.contains(&pol) {
            polarizations.push(pol);
        }
    }

    // Orbit from first annotation (identical across all)
    let orbit = build_orbit(&first.general_annotation.orbit_list)?;

    // Start/stop times: overall min/max across all annotations
    let mut start_time: Option<DateTime<Utc>> = None;
    let mut stop_time: Option<DateTime<Utc>> = None;
    for ann in annotations {
        let st = parse_s1_timestamp(&ann.ads_header.start_time, "startTime")?;
        let et = parse_s1_timestamp(&ann.ads_header.stop_time, "stopTime")?;
        start_time = Some(match start_time {
            Some(prev) => prev.min(st),
            None => st,
        });
        stop_time = Some(match stop_time {
            Some(prev) => prev.max(et),
            None => et,
        });
    }
    let start_time = start_time.ok_or(ParseError::MissingField("startTime"))?;
    let stop_time = stop_time.ok_or(ParseError::MissingField("stopTime"))?;

    // Group annotations by sub-swath
    let mut by_swath: BTreeMap<String, Vec<&ProductXml>> = BTreeMap::new();
    for ann in annotations {
        by_swath
            .entry(ann.ads_header.swath.clone())
            .or_default()
            .push(ann);
    }

    let mut sub_swaths = Vec::new();
    let mut bursts = Vec::new();
    let mut all_lats = Vec::new();
    let mut all_lons = Vec::new();

    for (_swath_name, swath_anns) in &by_swath {
        let ann = swath_anns[0]; // First annotation per swath (all pols share geometry)
        let swath_id = parse_subswath_id(&ann.ads_header.swath)?;
        let img = &ann.image_annotation.image_information;
        let timing = &ann.swath_timing;

        // PRF from downlink info
        let prf_hz = ann
            .general_annotation
            .downlink_information_list
            .items
            .first()
            .map(|di| di.prf)
            .ok_or(ParseError::MissingField("prf"))?;

        let lines_per_burst = timing.lines_per_burst as usize;
        let num_bursts = timing.burst_list.bursts.len();

        // Parse burst azimuth times for duration calculation
        let burst_utc_times: Vec<DateTime<Utc>> = timing
            .burst_list
            .bursts
            .iter()
            .map(|b| parse_s1_timestamp(&b.azimuth_time, "burst azimuthTime"))
            .collect::<Result<Vec<_>, _>>()?;

        // Per-subswath burst cycle time: median inter-burst interval.
        // ESA annotation XML has no burstDuration element; the cycle time must be
        // derived from consecutive burst azimuth times. Verified on S1A and S1B:
        // deltas are quantized to PRI ticks (~2.7565, 2.7586, 2.7606 s for IW),
        // all within ~4 ms, so the median is a robust estimator.
        let burst_cycle_time_s = if num_bursts >= 2 {
            let mut deltas: Vec<f64> = burst_utc_times
                .windows(2)
                .map(|w| (w[1] - w[0]).num_microseconds().unwrap_or(0) as f64 / 1e6) // SAFETY-OK: chrono microseconds cannot overflow for burst-cycle deltas
                .collect();
            deltas.sort_by(|a, b| a.partial_cmp(b).unwrap());
            deltas[deltas.len() / 2] // median
        } else {
            lines_per_burst as f64 * img.azimuth_time_interval
        };

        // Parse Doppler centroid estimates (top-level `<dopplerCentroid>` block).
        // Empty Vec when absent — not an error here; the InSAR entry point checks.
        let dc_estimates: Vec<DcEstimate> = match &ann.doppler_centroid {
            Some(dc_block) => {
                let mut out = Vec::with_capacity(dc_block.dc_estimate_list.estimates.len());
                for e in &dc_block.dc_estimate_list.estimates {
                    let t = parse_s1_timestamp(&e.azimuth_time, "dcEstimate azimuthTime")?;
                    let poly = &e.data_dc_polynomial;
                    if poly.len() != 3 {
                        return Err(ParseError::IncompatiblePolynomial {
                            context: "dataDcPolynomial",
                            found: poly.len(),
                        });
                    }
                    out.push(DcEstimate {
                        azimuth_time: t,
                        t0_s: e.t0,
                        data_dc_poly: [poly[0], poly[1], poly[2]],
                    });
                }
                out
            }
            None => Vec::new(),
        };

        // Parse azimuth FM rate estimates (inside `<generalAnnotation>`).
        let fm_rates: Vec<AzimuthFmRate> =
            match &ann.general_annotation.azimuth_fm_rate_list {
                Some(fm_list) => {
                    let mut out = Vec::with_capacity(fm_list.rates.len());
                    for r in &fm_list.rates {
                        let t =
                            parse_s1_timestamp(&r.azimuth_time, "azimuthFmRate azimuthTime")?;
                        let poly = &r.polynomial;
                        if poly.len() != 3 {
                            return Err(ParseError::IncompatiblePolynomial {
                                context: "azimuthFmRatePolynomial",
                                found: poly.len(),
                            });
                        }
                        out.push(AzimuthFmRate {
                            azimuth_time: t,
                            t0_s: r.t0,
                            poly: [poly[0], poly[1], poly[2]],
                        });
                    }
                    out
                }
                None => Vec::new(),
            };

        sub_swaths.push(SubSwathMetadata {
            id: swath_id,
            burst_count: num_bursts,
            lines_per_burst,
            range_samples: img.number_of_samples as usize,
            azimuth_samples: img.number_of_lines as usize,
            first_line: 0,
            last_line: img.number_of_lines as usize,
            first_sample: 0,
            last_sample: img.number_of_samples as usize,
            range_pixel_spacing_m: img.range_pixel_spacing,
            azimuth_pixel_spacing_m: img.azimuth_pixel_spacing,
            slant_range_time_s: img.slant_range_time,
            azimuth_time_interval_s: img.azimuth_time_interval,
            prf_hz,
            burst_cycle_time_s,
            dc_estimates,
            fm_rates,
        });

        // Build burst entries
        for (i, burst_xml) in timing.burst_list.bursts.iter().enumerate() {
            let burst_utc = burst_utc_times[i];

            // firstValidSample / lastValidSample are per-line arrays (one i32 per
            // azimuth line). Value -1 marks invalid guard lines at burst edges.
            // Verified on S1A + S1B: all valid (non-negative) entries within a
            // single burst are identical, so any non-negative value is representative.
            let first_valid = match first_nonneg(&burst_xml.first_valid_sample) {
                Some(v) => v,
                None => {
                    return Err(ParseError::InvalidValue {
                        field: "firstValidSample",
                        detail: format!(
                            "burst {} in {} has no valid samples (all -1)",
                            i, swath_id
                        ),
                    });
                }
            };
            let last_valid = match first_nonneg(&burst_xml.last_valid_sample) {
                Some(v) => v + 1, // exclusive end
                None => {
                    return Err(ParseError::InvalidValue {
                        field: "lastValidSample",
                        detail: format!(
                            "burst {} in {} has no valid samples (all -1)",
                            i, swath_id
                        ),
                    });
                }
            };

            bursts.push(BurstEntry {
                subswath_id: swath_id,
                burst_index: i,
                azimuth_time_utc: burst_utc,
                first_line: i * lines_per_burst,
                last_line: (i + 1) * lines_per_burst,
                first_valid_sample: first_valid,
                last_valid_sample: last_valid,
                slice_index: 0,
            });
        }

        // Bounding box from geolocation grid
        for pt in &ann.geolocation_grid.point_list.points {
            all_lats.push(pt.latitude);
            all_lons.push(pt.longitude);
        }
    }

    // Sort for deterministic output
    sub_swaths.sort_by_key(|s| s.id);
    bursts.sort_by(|a, b| {
        a.subswath_id
            .cmp(&b.subswath_id)
            .then(a.burst_index.cmp(&b.burst_index))
    });

    // Bounding box
    let bounding_box = if all_lats.is_empty() {
        return Err(ParseError::MissingField("geolocation grid points"));
    } else {
        let min_lat_deg = all_lats.iter().copied().fold(f64::INFINITY, f64::min);
        let max_lat_deg = all_lats.iter().copied().fold(f64::NEG_INFINITY, f64::max);
        let raw_min_lon = all_lons.iter().copied().fold(f64::INFINITY, f64::min);
        let raw_max_lon = all_lons.iter().copied().fold(f64::NEG_INFINITY, f64::max);

        // If the naive longitude span > 180°, the scene straddles the
        // anti-meridian (e.g. GCPs at 178° and −178° give span = 356°,
        // not a real 356°-wide scene).  Re-wrap negative longitudes by
        // +360° so the range becomes e.g. [178°, 182°] — max_lon_deg > 180
        // is the canonical internal representation for anti-meridian
        // crossing scenes.  Consumers that need standard [−180, 180]
        // coordinates (STAC JSON, tile names) normalise on output.
        let (min_lon_deg, max_lon_deg) = if raw_max_lon - raw_min_lon > 180.0 {
            let rewrapped: Vec<f64> = all_lons
                .iter()
                .map(|&l| if l < 0.0 { l + 360.0 } else { l })
                .collect();
            (
                rewrapped.iter().copied().fold(f64::INFINITY, f64::min),
                rewrapped.iter().copied().fold(f64::NEG_INFINITY, f64::max),
            )
        } else {
            (raw_min_lon, raw_max_lon)
        };

        BoundingBox { min_lat_deg, max_lat_deg, min_lon_deg, max_lon_deg }
    };

    let scene = SceneMetadata {
        product_id,
        mission,
        acquisition_mode: mode,
        polarizations,
        start_time,
        stop_time,
        radar_frequency_hz,
        range_sampling_rate_hz,
        bounding_box,
        sub_swaths,
        bursts,
        orbit,
    };

    scene.validated().map_err(ParseError::Validation)
}

// ── Geolocation grid parsing ────────────────────────────────────────

/// Per-subswath geolocation grid, keyed by `(SubSwathId, Polarization)`.
///
/// Each subswath has a sparse grid (~10 azimuth × 21 range points) with
/// incidence angle, elevation angle, slant range time, geographic coordinates,
/// and height above ellipsoid.
///
/// Follows the same pattern as [`crate::calibration::parse_calibration_noise`]:
/// a standalone function that reads from the SAFE directory and returns structured data.
pub fn parse_geolocation_grids(
    safe_path: &Path,
) -> Result<Vec<(SubSwathId, Vec<GeolocationGridPoint>)>, ParseError> {
    let annotation_dir = safe_path.join("annotation");
    let mut xml_paths: Vec<_> = std::fs::read_dir(&annotation_dir)
        .map_err(|e| {
            ParseError::Io(std::io::Error::new(
                e.kind(),
                format!("{}: {}", annotation_dir.display(), e),
            ))
        })?
        .filter_map(|entry| entry.ok())
        .filter(|entry| {
            entry
                .file_type()
                .map(|t| t.is_file())
                .unwrap_or(false) // SAFETY-OK: failed file_type → skip the entry; not a numeric path
        })
        .filter(|entry| entry.file_name().to_string_lossy().ends_with(".xml"))
        .map(|entry| entry.path())
        .collect();

    xml_paths.sort();

    // Group by subswath; use first polarization per swath (geometry is identical across pols)
    let mut seen_swaths = std::collections::HashSet::new();
    let mut result = Vec::new();

    for path in &xml_paths {
        let content = std::fs::read_to_string(path)?;
        let ann = parse_annotation_xml(&content)?;
        let swath_id = parse_subswath_id(&ann.ads_header.swath)?;

        // Geolocation grid is identical across polarizations for the same swath
        if !seen_swaths.insert(swath_id) {
            continue;
        }

        let mut points = Vec::with_capacity(ann.geolocation_grid.point_list.points.len());
        for pt in &ann.geolocation_grid.point_list.points {
            let azimuth_time_utc =
                parse_s1_timestamp(&pt.azimuth_time, "geolocationGrid azimuthTime")?;
            points.push(GeolocationGridPoint {
                azimuth_time_utc,
                slant_range_time_s: pt.slant_range_time,
                line: pt.line,
                pixel: pt.pixel,
                latitude_deg: pt.latitude,
                longitude_deg: pt.longitude,
                height_m: pt.height,
                incidence_angle_deg: pt.incidence_angle,
                elevation_angle_deg: pt.elevation_angle,
            });
        }

        result.push((swath_id, points));
    }

    result.sort_by_key(|(id, _)| *id);
    Ok(result)
}

// ── Helpers ─────────────────────────────────────────────────────────

fn parse_s1_timestamp(s: &str, field: &'static str) -> Result<DateTime<Utc>, ParseError> {
    let naive = NaiveDateTime::parse_from_str(s.trim(), "%Y-%m-%dT%H:%M:%S%.f").map_err(|e| {
        ParseError::InvalidValue {
            field,
            detail: format!("'{}': {}", s, e),
        }
    })?;
    Ok(naive.and_utc())
}

fn parse_mission(s: &str) -> Result<Mission, ParseError> {
    match s.trim() {
        "S1A" => Ok(Mission::S1A),
        "S1B" => Ok(Mission::S1B),
        other => Err(ParseError::InvalidValue {
            field: "missionId",
            detail: other.to_string(),
        }),
    }
}

fn parse_mode(s: &str) -> Result<AcquisitionMode, ParseError> {
    match s.trim() {
        "IW" => Ok(AcquisitionMode::IW),
        "EW" => Ok(AcquisitionMode::EW),
        "SM" => Ok(AcquisitionMode::SM),
        "WV" => Ok(AcquisitionMode::WV),
        other => Err(ParseError::InvalidValue {
            field: "mode",
            detail: other.to_string(),
        }),
    }
}

fn parse_polarization(s: &str) -> Result<Polarization, ParseError> {
    match s.trim() {
        "VV" => Ok(Polarization::VV),
        "VH" => Ok(Polarization::VH),
        "HV" => Ok(Polarization::HV),
        "HH" => Ok(Polarization::HH),
        other => Err(ParseError::InvalidValue {
            field: "polarisation",
            detail: other.to_string(),
        }),
    }
}

fn parse_subswath_id(s: &str) -> Result<SubSwathId, ParseError> {
    match s.trim() {
        "IW1" => Ok(SubSwathId::IW1),
        "IW2" => Ok(SubSwathId::IW2),
        "IW3" => Ok(SubSwathId::IW3),
        "EW1" => Ok(SubSwathId::EW1),
        "EW2" => Ok(SubSwathId::EW2),
        "EW3" => Ok(SubSwathId::EW3),
        "EW4" => Ok(SubSwathId::EW4),
        "EW5" => Ok(SubSwathId::EW5),
        other => Err(ParseError::InvalidValue {
            field: "swath",
            detail: other.to_string(),
        }),
    }
}

fn build_orbit(
    orbit_list: &xml::OrbitListXml,
) -> Result<OrbitData, ParseError> {
    if orbit_list.orbits.is_empty() {
        return Err(ParseError::MissingField("orbit state vectors"));
    }

    let mut state_vectors: Vec<StateVector> = orbit_list
        .orbits
        .iter()
        .map(|o| {
            let time = parse_s1_timestamp(&o.time, "orbit time")?;
            Ok(StateVector {
                time,
                position_m: [o.position.x, o.position.y, o.position.z],
                velocity_m_s: [o.velocity.x, o.velocity.y, o.velocity.z],
            })
        })
        .collect::<Result<Vec<_>, ParseError>>()?;

    // Sort by time
    state_vectors.sort_by_key(|sv| sv.time);

    let reference_epoch = state_vectors[0].time;

    Ok(OrbitData {
        reference_epoch,
        state_vectors,
    })
}

/// Representative valid-sample value from a per-line burst array.
///
/// Within a single burst, all non-negative entries in `firstValidSample` and
/// `lastValidSample` are identical (verified across S1A/S1B IW products).
/// The -1 entries mark invalid guard lines at burst edges (~35–54 of ~1500 lines).
fn first_nonneg(samples: &[i32]) -> Option<usize> {
    samples.iter().find(|&&v| v >= 0).map(|&v| v as usize)
}

// ── Tests ───────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::PathBuf;

    /// Path to S1A test scene (may not exist in CI).
    fn s1a_safe() -> PathBuf {
        PathBuf::from("/home/datacube/dev/SARdine/data/SLC/S1A_IW_SLC__1SDV_20201005T170824_20201005T170851_034664_04098A_1E66.SAFE")
    }

    /// Path to S1B test scene (may not exist in CI).
    fn s1b_safe() -> PathBuf {
        PathBuf::from("/home/datacube/dev/SARdine/data/SLC/S1B_IW_SLC__1SDV_20190123T053348_20190123T053415_014617_01B3D4_E833.SAFE")
    }

    fn have_test_data() -> bool {
        s1a_safe().join("annotation").is_dir()
    }

    // ── Unit tests ──────────────────────────────────────────────────

    #[test]
    fn test_parse_s1_timestamp() {
        let dt = parse_s1_timestamp("2020-10-05T17:08:24.055395", "test").unwrap();
        assert_eq!(dt.format("%Y-%m-%d %H:%M:%S").to_string(), "2020-10-05 17:08:24");
    }

    #[test]
    fn test_parse_mission() {
        assert_eq!(parse_mission("S1A").unwrap(), Mission::S1A);
        assert_eq!(parse_mission("S1B").unwrap(), Mission::S1B);
        assert!(parse_mission("S2A").is_err());
    }

    #[test]
    fn test_parse_mode() {
        assert_eq!(parse_mode("IW").unwrap(), AcquisitionMode::IW);
        assert_eq!(parse_mode("EW").unwrap(), AcquisitionMode::EW);
        assert!(parse_mode("XX").is_err());
    }

    #[test]
    fn test_parse_polarization() {
        assert_eq!(parse_polarization("VV").unwrap(), Polarization::VV);
        assert_eq!(parse_polarization("VH").unwrap(), Polarization::VH);
        assert!(parse_polarization("XY").is_err());
    }

    #[test]
    fn test_parse_subswath_id() {
        assert_eq!(parse_subswath_id("IW1").unwrap(), SubSwathId::IW1);
        assert_eq!(parse_subswath_id("EW5").unwrap(), SubSwathId::EW5);
        assert!(parse_subswath_id("IW4").is_err());
    }

    #[test]
    fn test_first_nonneg() {
        assert_eq!(first_nonneg(&[-1, -1, -1, 1184, 1184, -1]), Some(1184));
        assert_eq!(first_nonneg(&[-1, -1, -1]), None);
        assert_eq!(first_nonneg(&[0, 5, 10]), Some(0));
        assert_eq!(first_nonneg(&[]), None);
    }

    // ── Integration tests with real data ────────────────────────────

    #[test]
    fn test_parse_single_annotation_xml() {
        if !have_test_data() {
            eprintln!("Skipping: test data not available");
            return;
        }

        // Read IW1-VV annotation
        let ann_dir = s1a_safe().join("annotation");
        let xml_path = std::fs::read_dir(&ann_dir)
            .unwrap()
            .filter_map(|e| e.ok())
            .find(|e| {
                let n = e.file_name().to_string_lossy().to_string();
                n.contains("iw1") && n.contains("vv") && n.ends_with(".xml")
            })
            .expect("IW1-VV annotation not found")
            .path();

        let content = std::fs::read_to_string(&xml_path).unwrap();
        let parsed = parse_annotation_xml(&content).unwrap();

        // Verify header
        assert_eq!(parsed.ads_header.mission_id, "S1A");
        assert_eq!(parsed.ads_header.mode, "IW");
        assert_eq!(parsed.ads_header.swath, "IW1");
        assert_eq!(parsed.ads_header.polarisation, "VV");

        // Verify radar params
        let pi = &parsed.general_annotation.product_information;
        assert!((pi.radar_frequency - 5.405e9).abs() < 1e6);
        assert!(pi.range_sampling_rate > 6e7);

        // Verify orbit
        assert_eq!(parsed.general_annotation.orbit_list.orbits.len(), 17);

        // Verify image dimensions
        let img = &parsed.image_annotation.image_information;
        assert_eq!(img.number_of_samples, 22111);
        assert_eq!(img.number_of_lines, 13482);
        assert!(img.slant_range_time > 0.005);
        assert!(img.azimuth_time_interval > 0.002);

        // Verify bursts
        let timing = &parsed.swath_timing;
        assert_eq!(timing.lines_per_burst, 1498);
        assert_eq!(timing.burst_list.bursts.len(), 9);

        // First burst valid samples
        let b0 = &timing.burst_list.bursts[0];
        assert_eq!(b0.first_valid_sample.len(), 1498);
        assert!(b0.first_valid_sample.iter().any(|&v| v >= 0));
        assert!(b0.azimuth_anx_time > 800.0);

        // Geolocation grid
        assert!(parsed.geolocation_grid.point_list.points.len() >= 100);
    }

    #[test]
    fn test_parse_s1a_full_scene() {
        if !have_test_data() {
            eprintln!("Skipping: test data not available");
            return;
        }

        let meta = parse_safe_directory(&s1a_safe()).expect("S1A parse failed");

        // Mission & mode
        assert_eq!(meta.mission, Mission::S1A);
        assert_eq!(meta.acquisition_mode, AcquisitionMode::IW);

        // Polarizations
        assert_eq!(meta.polarizations.len(), 2);
        assert!(meta.polarizations.contains(&Polarization::VV));
        assert!(meta.polarizations.contains(&Polarization::VH));

        // Sub-swaths
        assert_eq!(meta.sub_swaths.len(), 3);
        assert_eq!(meta.sub_swaths[0].id, SubSwathId::IW1);
        assert_eq!(meta.sub_swaths[1].id, SubSwathId::IW2);
        assert_eq!(meta.sub_swaths[2].id, SubSwathId::IW3);

        // IW1 specifics
        let iw1 = &meta.sub_swaths[0];
        assert_eq!(iw1.burst_count, 9);
        assert_eq!(iw1.lines_per_burst, 1498);
        assert_eq!(iw1.range_samples, 22111);
        assert_eq!(iw1.azimuth_samples, 13482);
        assert!(iw1.prf_hz > 1700.0 && iw1.prf_hz < 1720.0);
        assert!(iw1.burst_cycle_time_s > 2.7 && iw1.burst_cycle_time_s < 2.8);
        assert!(iw1.slant_range_time_s > 0.005);
        assert!(iw1.range_pixel_spacing_m > 2.0 && iw1.range_pixel_spacing_m < 3.0);

        // Orbit
        assert_eq!(meta.orbit.state_vectors.len(), 17);

        // Bursts: 9 per sub-swath × 3 sub-swaths
        let iw1_bursts: Vec<_> = meta
            .bursts
            .iter()
            .filter(|b| b.subswath_id == SubSwathId::IW1)
            .collect();
        assert_eq!(iw1_bursts.len(), 9);

        // Burst times should be within scene time window
        assert!(iw1_bursts[0].azimuth_time_utc >= meta.start_time);
        assert!(iw1_bursts[8].azimuth_time_utc <= meta.stop_time);
        // Bursts should be monotonically increasing in time
        for w in iw1_bursts.windows(2) {
            assert!(w[1].azimuth_time_utc > w[0].azimuth_time_utc);
        }

        // Valid samples
        assert!(iw1_bursts[0].first_valid_sample > 0);
        assert!(iw1_bursts[0].last_valid_sample > iw1_bursts[0].first_valid_sample);

        // Radar frequency
        assert!((meta.radar_frequency_hz - 5.405e9).abs() < 1e6);

        // Bounding box sanity
        assert!(meta.bounding_box.max_lat_deg > meta.bounding_box.min_lat_deg);
        assert!(meta.bounding_box.max_lon_deg > meta.bounding_box.min_lon_deg);

        // Product ID
        assert!(meta.product_id.starts_with("S1A_IW_SLC"));
    }

    #[test]
    fn test_parse_s1b_full_scene() {
        let path = s1b_safe();
        if !path.join("annotation").is_dir() {
            eprintln!("Skipping: S1B test data not available");
            return;
        }

        let meta = parse_safe_directory(&path).expect("S1B parse failed");

        assert_eq!(meta.mission, Mission::S1B);
        assert_eq!(meta.acquisition_mode, AcquisitionMode::IW);
        assert_eq!(meta.polarizations.len(), 2);
        assert_eq!(meta.sub_swaths.len(), 3);
        assert!(meta.product_id.starts_with("S1B_IW_SLC"));
    }

    /// Verify that valid sample values are constant within each burst.
    ///
    /// This invariant underpins `first_nonneg()`: since all valid (≥0) entries
    /// in firstValidSample / lastValidSample are identical within a burst,
    /// picking any one of them is representative.
    #[test]
    fn test_valid_sample_constant_within_burst() {
        if !have_test_data() {
            eprintln!("Skipping: test data not available");
            return;
        }

        for safe_path in [s1a_safe(), s1b_safe()] {
            if !safe_path.join("annotation").is_dir() {
                continue;
            }
            let ann_dir = safe_path.join("annotation");
            for entry in std::fs::read_dir(&ann_dir).unwrap() {
                let path = entry.unwrap().path();
                if !path.extension().map_or(false, |e| e == "xml") {
                    continue;
                }
                let content = std::fs::read_to_string(&path).unwrap();
                let parsed = parse_annotation_xml(&content).unwrap();
                let swath = &parsed.ads_header.swath;

                for (bi, burst) in parsed.swath_timing.burst_list.bursts.iter().enumerate() {
                    let valid_fvs: Vec<i32> =
                        burst.first_valid_sample.iter().copied().filter(|&v| v >= 0).collect();
                    let valid_lvs: Vec<i32> =
                        burst.last_valid_sample.iter().copied().filter(|&v| v >= 0).collect();

                    if let (Some(&first_f), Some(&first_l)) = (valid_fvs.first(), valid_lvs.first())
                    {
                        assert!(
                            valid_fvs.iter().all(|&v| v == first_f),
                            "{} burst {}: firstValidSample not constant (min={}, max={})",
                            swath,
                            bi,
                            valid_fvs.iter().min().unwrap(),
                            valid_fvs.iter().max().unwrap(),
                        );
                        assert!(
                            valid_lvs.iter().all(|&v| v == first_l),
                            "{} burst {}: lastValidSample not constant (min={}, max={})",
                            swath,
                            bi,
                            valid_lvs.iter().min().unwrap(),
                            valid_lvs.iter().max().unwrap(),
                        );
                    }
                }
            }
        }
    }

    /// Verify burst duration is within expected IW cycle time range.
    ///
    /// For IW mode, the per-subswath burst cycle time should be ~2.758 s
    /// (quantized to PRI ticks). All inter-burst deltas should be within
    /// ~10 ms of each other.
    #[test]
    fn test_burst_duration_iw_range() {
        if !have_test_data() {
            eprintln!("Skipping: test data not available");
            return;
        }

        for safe_path in [s1a_safe(), s1b_safe()] {
            if !safe_path.join("annotation").is_dir() {
                continue;
            }
            let meta = parse_safe_directory(&safe_path).unwrap();
            for sw in &meta.sub_swaths {
                assert!(
                    sw.burst_cycle_time_s > 2.72 && sw.burst_cycle_time_s < 2.80,
                    "{}: burst_cycle_time_s={} outside expected IW range [2.72, 2.80]",
                    sw.id,
                    sw.burst_cycle_time_s,
                );
            }
        }
    }

    // ── Geolocation grid tests ──────────────────────────────────────

    #[test]
    fn test_parse_geolocation_grids_s1a() {
        if !have_test_data() {
            eprintln!("Skipping: test data not available");
            return;
        }

        let grids = parse_geolocation_grids(&s1a_safe()).unwrap();
        // IW mode: expect 3 subswaths
        assert_eq!(grids.len(), 3, "expected 3 subswath grids");

        let ids: Vec<_> = grids.iter().map(|(id, _)| *id).collect();
        assert_eq!(ids, vec![SubSwathId::IW1, SubSwathId::IW2, SubSwathId::IW3]);

        for (id, points) in &grids {
            // 10 azimuth × 21 range = 210 points
            assert_eq!(
                points.len(),
                210,
                "{}: expected 210 grid points, got {}",
                id,
                points.len()
            );

            // Check basic value ranges
            for pt in points {
                assert!(
                    pt.incidence_angle_deg > 25.0 && pt.incidence_angle_deg < 50.0,
                    "{}: incidence angle {:.2} outside [25, 50] deg",
                    id,
                    pt.incidence_angle_deg
                );
                assert!(
                    pt.elevation_angle_deg > 20.0 && pt.elevation_angle_deg < 50.0,
                    "{}: elevation angle {:.2} outside [20, 50] deg",
                    id,
                    pt.elevation_angle_deg
                );
                assert!(
                    pt.slant_range_time_s > 4e-3 && pt.slant_range_time_s < 8e-3,
                    "{}: slant range time {} outside [4ms, 8ms]",
                    id,
                    pt.slant_range_time_s
                );
                assert!(
                    pt.latitude_deg > 40.0 && pt.latitude_deg < 60.0,
                    "{}: latitude {:.2} outside expected range",
                    id,
                    pt.latitude_deg
                );
            }
        }
    }

    #[test]
    fn test_geolocation_grid_incidence_increases_in_range() {
        if !have_test_data() {
            eprintln!("Skipping: test data not available");
            return;
        }

        let grids = parse_geolocation_grids(&s1a_safe()).unwrap();
        for (id, points) in &grids {
            // Group by line, check incidence angle increases with pixel
            let first_line = points[0].line;
            let row: Vec<_> = points.iter().filter(|p| p.line == first_line).collect();
            for w in row.windows(2) {
                assert!(
                    w[1].incidence_angle_deg > w[0].incidence_angle_deg,
                    "{}: incidence should increase in range: pixel {} ({:.2}) >= pixel {} ({:.2})",
                    id,
                    w[1].pixel,
                    w[1].incidence_angle_deg,
                    w[0].pixel,
                    w[0].incidence_angle_deg,
                );
            }
        }
    }

    #[test]
    fn test_geolocation_grid_iw_incidence_ranges() {
        if !have_test_data() {
            eprintln!("Skipping: test data not available");
            return;
        }

        let grids = parse_geolocation_grids(&s1a_safe()).unwrap();
        // IW1 should be near-range (~30°), IW3 far-range (~46°)
        let angles: Vec<(SubSwathId, f64, f64)> = grids
            .iter()
            .map(|(id, pts)| {
                let min = pts.iter().map(|p| p.incidence_angle_deg).fold(f64::INFINITY, f64::min);
                let max = pts.iter().map(|p| p.incidence_angle_deg).fold(f64::NEG_INFINITY, f64::max);
                (*id, min, max)
            })
            .collect();

        // IW1 near < IW2 near < IW3 near
        assert!(angles[0].1 < angles[1].1, "IW1 min inc should < IW2 min inc");
        assert!(angles[1].1 < angles[2].1, "IW2 min inc should < IW3 min inc");

        // Approximate expected ranges
        assert!(angles[0].1 > 29.0 && angles[0].1 < 32.0, "IW1 min inc: {:.2}", angles[0].1);
        assert!(angles[2].2 > 44.0 && angles[2].2 < 47.0, "IW3 max inc: {:.2}", angles[2].2);
    }

    /// Verify that the Doppler centroid and azimuth FM rate blocks are parsed
    /// correctly from a real S1B IW1 VV annotation XML.
    ///
    /// Checks:
    /// - Exactly 10 estimates are parsed (matches `count="10"` in XML).
    /// - Reference time `t0` is within the expected near-range slant-range window
    ///   (~5.35 ms for S1B IW1).
    /// - First estimate polynomial a0 ≈ −18.7 Hz (data DC centroid at near range).
    /// - FM rate polynomial a0 ≈ −2318 Hz/s.
    /// - `evaluate()` produces physically plausible values (|f_dc| < 3000 Hz).
    #[test]
    fn test_parse_dc_and_fm_rate_s1b_iw1() {
        let path = s1b_safe();
        if !path.join("annotation").is_dir() {
            eprintln!("Skipping: S1B test data not available");
            return;
        }

        let meta = parse_safe_directory(&path).expect("S1B parse failed");

        let iw1 = meta
            .sub_swaths
            .iter()
            .find(|s| s.id == SubSwathId::IW1)
            .expect("IW1 not found");

        // DC estimates: should have 10 entries for IW1 SLC
        assert_eq!(
            iw1.dc_estimates.len(),
            10,
            "IW1 dc_estimates count: expected 10, got {}",
            iw1.dc_estimates.len()
        );

        // FM rate estimates: should also have 10 entries
        assert_eq!(
            iw1.fm_rates.len(),
            10,
            "IW1 fm_rates count: expected 10, got {}",
            iw1.fm_rates.len()
        );

        // First DC estimate: t0 should be near-range slant-range time (≈5.35 ms for S1B IW1)
        let dc0 = &iw1.dc_estimates[0];
        assert!(
            dc0.t0_s > 5.0e-3 && dc0.t0_s < 6.0e-3,
            "dc0.t0_s = {:.6e} — expected in [5ms, 6ms]",
            dc0.t0_s
        );

        // First DC estimate polynomial a0 ≈ −18.7 Hz (data Doppler centroid at t0)
        assert!(
            (dc0.data_dc_poly[0] + 18.7).abs() < 5.0,
            "dc0 a0 = {:.3} Hz — expected ≈ −18.7 Hz",
            dc0.data_dc_poly[0]
        );

        // evaluate() at t0 should equal a0 (polynomial zero-offset)
        let f_at_t0 = dc0.evaluate(dc0.t0_s);
        assert!(
            (f_at_t0 - dc0.data_dc_poly[0]).abs() < 1e-6,
            "evaluate(t0) = {f_at_t0:.6} ≠ a0 = {:.6}",
            dc0.data_dc_poly[0]
        );

        // Doppler centroid at any pixel should be physically plausible (|f_dc| < 3000 Hz)
        // Check near-range and far-range of IW1
        let near_tau = iw1.slant_range_time_s;
        let far_tau = iw1.slant_range_time_s
            + (iw1.range_samples as f64) * iw1.range_pixel_spacing_m
                / crate::types::SPEED_OF_LIGHT_M_S
                * 2.0;
        let f_near = dc0.evaluate(near_tau);
        let f_far = dc0.evaluate(far_tau);
        assert!(
            f_near.abs() < 3000.0,
            "f_dc at near range = {f_near:.1} Hz — exceeds 3000 Hz"
        );
        assert!(
            f_far.abs() < 3000.0,
            "f_dc at far range = {f_far:.1} Hz — exceeds 3000 Hz"
        );

        // First FM rate estimate: a0 ≈ −2318 Hz/s
        let fm0 = &iw1.fm_rates[0];
        assert!(
            (fm0.poly[0] + 2318.0).abs() < 50.0,
            "fm0 a0 = {:.1} Hz/s — expected ≈ −2318 Hz/s",
            fm0.poly[0]
        );

        // FM rate evaluate() at t0 should equal a0
        let ka_at_t0 = fm0.evaluate(fm0.t0_s);
        assert!(
            (ka_at_t0 - fm0.poly[0]).abs() < 1e-6,
            "fm rate evaluate(t0) = {ka_at_t0:.6} ≠ a0 = {:.6}",
            fm0.poly[0]
        );

        // DC estimate azimuth times should be strictly increasing
        for w in iw1.dc_estimates.windows(2) {
            assert!(
                w[1].azimuth_time > w[0].azimuth_time,
                "dc_estimates not monotonically increasing in time"
            );
        }
    }

    // ── Anti-meridian bounding-box tests ──────────────────────────

    /// The bounding-box logic should re-wrap GCP longitudes when the naive
    /// min/max span exceeds 180°, producing max_lon > 180 rather than a
    /// spurious 356°-wide box.
    #[test]
    fn test_bbox_antimeridian_rewrap() {
        // Simulate GCP lons straddling ±180°.
        // Naive min = −179, max = 179, span = 358° → triggers re-wrap.
        // Re-wrapped: [179, 179, 181, 181] → min = 179, max = 181.
        let all_lons = vec![179.0_f64, 179.5, -179.5, -179.0];

        let raw_min = all_lons.iter().copied().fold(f64::INFINITY, f64::min);
        let raw_max = all_lons.iter().copied().fold(f64::NEG_INFINITY, f64::max);
        let (min_lon, max_lon) = if raw_max - raw_min > 180.0 {
            let rewrapped: Vec<f64> =
                all_lons.iter().map(|&l| if l < 0.0 { l + 360.0 } else { l }).collect();
            (
                rewrapped.iter().copied().fold(f64::INFINITY, f64::min),
                rewrapped.iter().copied().fold(f64::NEG_INFINITY, f64::max),
            )
        } else {
            (raw_min, raw_max)
        };

        // True span should be ≈ 2°, not 358°
        assert!(max_lon - min_lon < 180.0, "span should be < 180°, got {}", max_lon - min_lon);
        assert!(max_lon > 180.0, "anti-meridian crossing: max_lon should be > 180, got {max_lon}");
        assert!((min_lon - 179.0).abs() < 1e-10, "expected min_lon ≈ 179, got {min_lon}");
        assert!((max_lon - 181.0).abs() < 1e-10, "expected max_lon ≈ 181, got {max_lon}");
    }

    /// Normal scenes (not crossing ±180°) should be unaffected.
    #[test]
    fn test_bbox_normal_scene_unaffected() {
        let all_lons = vec![7.5_f64, 8.0, 8.5, 9.0];

        let raw_min = all_lons.iter().copied().fold(f64::INFINITY, f64::min);
        let raw_max = all_lons.iter().copied().fold(f64::NEG_INFINITY, f64::max);
        let (min_lon, max_lon) = if raw_max - raw_min > 180.0 {
            let rewrapped: Vec<f64> =
                all_lons.iter().map(|&l| if l < 0.0 { l + 360.0 } else { l }).collect();
            (
                rewrapped.iter().copied().fold(f64::INFINITY, f64::min),
                rewrapped.iter().copied().fold(f64::NEG_INFINITY, f64::max),
            )
        } else {
            (raw_min, raw_max)
        };

        assert!((min_lon - 7.5).abs() < 1e-10);
        assert!((max_lon - 9.0).abs() < 1e-10);
    }
}
