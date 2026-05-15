//! Calibration and noise annotation parsing for Sentinel-1 SAFE products.
//!
//! Provides [`parse_calibration_noise`] to read all calibration and noise XMLs
//! from a `.SAFE` directory and return structured [`CalibrationNoiseData`].
//!
//! # SAR-specific notes
//!
//! - Calibration LUT values (`sigma_nought`, `beta_nought`, `gamma`, `dn`) are in
//!   **linear units** (NOT dB). Typical sigma0 values are 280–340.
//! - The LUTs already incorporate `absoluteCalibrationConstant`. When calibrating
//!   via LUTs the correct formula is simply `σ⁰ = |DN|² / K²` where K comes from
//!   the `sigmaNought` LUT. Do **not** apply `absoluteCalibrationConstant` again.
//! - `absoluteCalibrationConstant` is stored for traceability/validation only.
//!   It equals 1.0 for S1A but can differ (e.g. 1.393 for S1B). See ESA
//!   Sentinel-1 Product Specification: the constant is provided "for completeness"
//!   and is already built into the per-pixel LUT values.
//! - Noise LUT values (`noise_range_lut`, `noise_azimuth_lut`) are in **linear power**.
//! - Pixel indices are sampled mostly every 40 range samples, with a shorter
//!   final step (e.g. 30 or 7) to reach the last pixel. The grid is constant
//!   across azimuth vectors within a swath. Interpolation to the full image grid
//!   is a processing step, not done here.
//! - Noise line indices can be negative (extrapolation beyond image start).
//! - One calibration file and one noise file exist per (swath, polarization).

use std::path::Path;

use chrono::{DateTime, NaiveDateTime, Utc};
use serde::Deserialize;

use crate::types::{Polarization, SubSwathId};

// ═══════════════════════════════════════════════════════════════════════
// Public types
// ═══════════════════════════════════════════════════════════════════════

/// A single calibration vector: LUT values at sampled pixel positions for one azimuth line.
#[derive(Debug, Clone)]
pub struct CalibrationVector {
    /// UTC time of this azimuth line.
    pub azimuth_time: DateTime<Utc>,
    /// Azimuth line index in the measurement image.
    pub line: i32,
    /// Range pixel indices where LUT values are sampled.
    pub pixels: Vec<u32>,
    /// σ⁰ (sigma nought) calibration LUT, linear units.
    pub sigma_nought: Vec<f32>,
    /// β⁰ (beta nought) calibration LUT, linear units.
    pub beta_nought: Vec<f32>,
    /// γ⁰ (gamma) calibration LUT, linear units.
    pub gamma: Vec<f32>,
    /// DN calibration LUT, linear units.
    pub dn: Vec<f32>,
}

/// Calibration annotation for one sub-swath and polarization.
#[derive(Debug, Clone)]
pub struct SwathCalibration {
    pub subswath_id: SubSwathId,
    pub polarization: Polarization,
    /// Absolute calibration constant from the XML header.
    ///
    /// **For LUT-based calibration this value should NOT be applied separately.**
    /// The per-pixel LUT values (`sigma_nought`, `beta_nought`, `gamma`, `dn`)
    /// already incorporate this constant. It is stored here for traceability
    /// and validation only.  Equals 1.0 for S1A; 1.393 for S1B (IPF 2.x).
    pub absolute_calibration_constant: f64,
    /// Calibration vectors sorted by azimuth line.
    pub vectors: Vec<CalibrationVector>,
}

/// A single noise range vector: noise LUT at sampled pixel positions for one azimuth line.
#[derive(Debug, Clone)]
pub struct NoiseRangeVector {
    /// UTC time of this azimuth line.
    pub azimuth_time: DateTime<Utc>,
    /// Azimuth line index (can be negative for extrapolation beyond image start).
    pub line: i32,
    /// Range pixel indices where noise LUT values are sampled.
    pub pixels: Vec<u32>,
    /// Noise range LUT values, linear power.
    pub noise_range_lut: Vec<f32>,
}

/// Azimuth noise vector: noise variation along azimuth for one sub-swath.
#[derive(Debug, Clone)]
pub struct NoiseAzimuthVector {
    pub swath: SubSwathId,
    pub first_azimuth_line: i32,
    pub first_range_sample: i32,
    pub last_azimuth_line: i32,
    pub last_range_sample: i32,
    /// Azimuth line indices where noise values are sampled.
    pub lines: Vec<i32>,
    /// Noise azimuth LUT values, linear units.
    pub noise_azimuth_lut: Vec<f32>,
}

/// Noise annotation for one sub-swath and polarization.
#[derive(Debug, Clone)]
pub struct SwathNoise {
    pub subswath_id: SubSwathId,
    pub polarization: Polarization,
    /// Noise range vectors sorted by azimuth line.
    pub range_vectors: Vec<NoiseRangeVector>,
    /// Noise azimuth vectors (typically one per swath in modern IPF).
    pub azimuth_vectors: Vec<NoiseAzimuthVector>,
}

/// Combined calibration and noise data for an entire SAFE product.
#[derive(Debug, Clone)]
pub struct CalibrationNoiseData {
    /// Calibration annotations, one per (swath, polarization).
    pub calibrations: Vec<SwathCalibration>,
    /// Noise annotations, one per (swath, polarization).
    pub noises: Vec<SwathNoise>,
}

// ═══════════════════════════════════════════════════════════════════════
// Error type
// ═══════════════════════════════════════════════════════════════════════

#[derive(Debug, thiserror::Error)]
pub enum CalibrationParseError {
    #[error("IO error: {0}")]
    Io(#[from] std::io::Error),

    #[error("calibration XML deserialization failed: {0}")]
    Xml(String),

    #[error("invalid value in calibration/noise XML for {field}: {detail}")]
    InvalidValue {
        field: &'static str,
        detail: String,
    },

    #[error("no calibration/noise XMLs found in {0}")]
    NoFiles(String),

    #[error("array length mismatch in {field}: expected {expected}, got {actual}")]
    ArrayLengthMismatch {
        field: &'static str,
        expected: usize,
        actual: usize,
    },
}

// ═══════════════════════════════════════════════════════════════════════
// Serde XML structs (private)
// ═══════════════════════════════════════════════════════════════════════

// ── Calibration XML ──

#[derive(Debug, Deserialize)]
struct CalibrationXml {
    #[serde(rename = "adsHeader")]
    ads_header: CalAdsHeaderXml,
    #[serde(rename = "calibrationInformation")]
    calibration_information: CalibrationInformationXml,
    #[serde(rename = "calibrationVectorList")]
    calibration_vector_list: CalibrationVectorListXml,
}

#[derive(Debug, Deserialize)]
struct CalAdsHeaderXml {
    polarisation: String,
    swath: String,
}

#[derive(Debug, Deserialize)]
struct CalibrationInformationXml {
    #[serde(rename = "absoluteCalibrationConstant")]
    absolute_calibration_constant: String,
}

#[derive(Debug, Deserialize)]
struct CalibrationVectorListXml {
    #[serde(rename = "calibrationVector", default)]
    vectors: Vec<CalibrationVectorXml>,
}

#[derive(Debug, Deserialize)]
struct CalibrationVectorXml {
    #[serde(rename = "azimuthTime")]
    azimuth_time: String,
    line: String,
    pixel: String,
    #[serde(rename = "sigmaNought")]
    sigma_nought: String,
    #[serde(rename = "betaNought")]
    beta_nought: String,
    gamma: String,
    dn: String,
}

// ── Noise XML ──

#[derive(Debug, Deserialize)]
struct NoiseXml {
    #[serde(rename = "adsHeader")]
    ads_header: CalAdsHeaderXml,
    #[serde(rename = "noiseRangeVectorList")]
    noise_range_vector_list: NoiseRangeVectorListXml,
    #[serde(rename = "noiseAzimuthVectorList")]
    noise_azimuth_vector_list: NoiseAzimuthVectorListXml,
}

#[derive(Debug, Deserialize)]
struct NoiseRangeVectorListXml {
    #[serde(rename = "noiseRangeVector", default)]
    vectors: Vec<NoiseRangeVectorXml>,
}

#[derive(Debug, Deserialize)]
struct NoiseRangeVectorXml {
    #[serde(rename = "azimuthTime")]
    azimuth_time: String,
    line: String,
    pixel: String,
    #[serde(rename = "noiseRangeLut")]
    noise_range_lut: String,
}

#[derive(Debug, Deserialize)]
struct NoiseAzimuthVectorListXml {
    #[serde(rename = "noiseAzimuthVector", default)]
    vectors: Vec<NoiseAzimuthVectorXml>,
}

#[derive(Debug, Deserialize)]
struct NoiseAzimuthVectorXml {
    swath: String,
    #[serde(rename = "firstAzimuthLine")]
    first_azimuth_line: String,
    #[serde(rename = "firstRangeSample")]
    first_range_sample: String,
    #[serde(rename = "lastAzimuthLine")]
    last_azimuth_line: String,
    #[serde(rename = "lastRangeSample")]
    last_range_sample: String,
    line: String,
    #[serde(rename = "noiseAzimuthLut")]
    noise_azimuth_lut: String,
}

// ═══════════════════════════════════════════════════════════════════════
// Public API
// ═══════════════════════════════════════════════════════════════════════

/// Parse all calibration and noise annotation XMLs from a `.SAFE` directory.
///
/// Reads files matching `calibration-*.xml` and `noise-*.xml` from
/// `<safe_path>/annotation/calibration/`.
pub fn parse_calibration_noise(
    safe_path: &Path,
) -> Result<CalibrationNoiseData, CalibrationParseError> {
    let cal_dir = safe_path.join("annotation").join("calibration");

    let entries: Vec<_> = std::fs::read_dir(&cal_dir)
        .map_err(|e| {
            CalibrationParseError::Io(std::io::Error::new(
                e.kind(),
                format!("{}: {}", cal_dir.display(), e),
            ))
        })?
        .filter_map(|entry| entry.ok())
        .filter(|entry| {
            entry
                .file_type()
                .map(|t| t.is_file())
                .unwrap_or(false) // SAFETY-OK: failed file_type → skip the entry; not a numeric path
        })
        .collect();

    let mut cal_paths = Vec::new();
    let mut noise_paths = Vec::new();

    for entry in &entries {
        let name = entry.file_name().to_string_lossy().to_string();
        if name.starts_with("calibration-") && name.ends_with(".xml") {
            cal_paths.push(entry.path());
        } else if name.starts_with("noise-") && name.ends_with(".xml") {
            noise_paths.push(entry.path());
        }
    }

    if cal_paths.is_empty() && noise_paths.is_empty() {
        return Err(CalibrationParseError::NoFiles(
            cal_dir.display().to_string(),
        ));
    }

    // Sort for deterministic order
    cal_paths.sort();
    noise_paths.sort();

    // Parse calibration files
    let mut calibrations = Vec::with_capacity(cal_paths.len());
    for path in &cal_paths {
        let content = std::fs::read_to_string(path)?;
        let cal = parse_calibration_xml(&content)?;
        calibrations.push(cal);
    }

    // Parse noise files
    let mut noises = Vec::with_capacity(noise_paths.len());
    for path in &noise_paths {
        let content = std::fs::read_to_string(path)?;
        let noise = parse_noise_xml(&content)?;
        noises.push(noise);
    }

    // Sort by (swath, polarization) for deterministic output
    calibrations.sort_by(|a, b| {
        a.subswath_id
            .cmp(&b.subswath_id)
            .then_with(|| a.polarization.cmp(&b.polarization))
    });
    noises.sort_by(|a, b| {
        a.subswath_id
            .cmp(&b.subswath_id)
            .then_with(|| a.polarization.cmp(&b.polarization))
    });

    Ok(CalibrationNoiseData {
        calibrations,
        noises,
    })
}

// ═══════════════════════════════════════════════════════════════════════
// Internal parsing
// ═══════════════════════════════════════════════════════════════════════

fn parse_calibration_xml(xml_content: &str) -> Result<SwathCalibration, CalibrationParseError> {
    let content = xml_content
        .strip_prefix('\u{feff}')
        .unwrap_or(xml_content); // SAFETY-OK: strip_prefix-on-no-match returns input unchanged by design (BOM optional)

    let raw: CalibrationXml =
        quick_xml::de::from_str(content).map_err(|e| CalibrationParseError::Xml(e.to_string()))?;

    let subswath_id = parse_subswath_id(&raw.ads_header.swath)?;
    let polarization = parse_polarization(&raw.ads_header.polarisation)?;
    let absolute_calibration_constant = parse_f64(
        &raw.calibration_information.absolute_calibration_constant,
        "absoluteCalibrationConstant",
    )?;

    let mut vectors = Vec::with_capacity(raw.calibration_vector_list.vectors.len());
    for v in &raw.calibration_vector_list.vectors {
        let azimuth_time = parse_timestamp(&v.azimuth_time, "calibrationVector/azimuthTime")?;
        let line = parse_i32(&v.line, "calibrationVector/line")?;
        let pixels = parse_space_u32(&v.pixel, "calibrationVector/pixel")?;
        let sigma_nought =
            parse_space_f32(&v.sigma_nought, "calibrationVector/sigmaNought")?;
        let beta_nought =
            parse_space_f32(&v.beta_nought, "calibrationVector/betaNought")?;
        let gamma = parse_space_f32(&v.gamma, "calibrationVector/gamma")?;
        let dn = parse_space_f32(&v.dn, "calibrationVector/dn")?;

        // Verify all LUT arrays have same length as pixel array
        let n = pixels.len();
        check_array_len("sigmaNought", n, sigma_nought.len())?;
        check_array_len("betaNought", n, beta_nought.len())?;
        check_array_len("gamma", n, gamma.len())?;
        check_array_len("dn", n, dn.len())?;

        vectors.push(CalibrationVector {
            azimuth_time,
            line,
            pixels,
            sigma_nought,
            beta_nought,
            gamma,
            dn,
        });
    }

    // Sort by line for consistent ordering
    vectors.sort_by_key(|v| v.line);

    // Validate pixel arrays: binary_search requires strictly ascending indices.
    for v in &vectors {
        if v.pixels.windows(2).any(|w| w[0] >= w[1]) {
            return Err(CalibrationParseError::InvalidValue {
                field: "pixel",
                detail: "pixel indices must be strictly ascending (required for interpolation)".into(),
            });
        }
    }

    // Validate LUT values: non-finite values would propagate silently into
    // the calibrated output as NaN/Inf and produce corrupt rasters.
    for v in &vectors {
        let lut_ok = v.sigma_nought.iter().all(|x| x.is_finite())
            && v.beta_nought.iter().all(|x| x.is_finite())
            && v.gamma.iter().all(|x| x.is_finite())
            && v.dn.iter().all(|x| x.is_finite());
        if !lut_ok {
            return Err(CalibrationParseError::InvalidValue {
                field: "calibrationVector",
                detail: format!(
                    "non-finite LUT value in vector at line {}",
                    v.line
                ),
            });
        }
    }

    Ok(SwathCalibration {
        subswath_id,
        polarization,
        absolute_calibration_constant,
        vectors,
    })
}

fn parse_noise_xml(xml_content: &str) -> Result<SwathNoise, CalibrationParseError> {
    let content = xml_content
        .strip_prefix('\u{feff}')
        .unwrap_or(xml_content); // SAFETY-OK: strip_prefix-on-no-match returns input unchanged by design (BOM optional)

    let raw: NoiseXml =
        quick_xml::de::from_str(content).map_err(|e| CalibrationParseError::Xml(e.to_string()))?;

    let subswath_id = parse_subswath_id(&raw.ads_header.swath)?;
    let polarization = parse_polarization(&raw.ads_header.polarisation)?;

    // Range vectors
    let mut range_vectors = Vec::with_capacity(raw.noise_range_vector_list.vectors.len());
    for v in &raw.noise_range_vector_list.vectors {
        let azimuth_time = parse_timestamp(&v.azimuth_time, "noiseRangeVector/azimuthTime")?;
        let line = parse_i32(&v.line, "noiseRangeVector/line")?;
        let pixels = parse_space_u32(&v.pixel, "noiseRangeVector/pixel")?;
        let noise_range_lut =
            parse_space_f32(&v.noise_range_lut, "noiseRangeVector/noiseRangeLut")?;

        check_array_len("noiseRangeLut", pixels.len(), noise_range_lut.len())?;

        range_vectors.push(NoiseRangeVector {
            azimuth_time,
            line,
            pixels,
            noise_range_lut,
        });
    }
    range_vectors.sort_by_key(|v| v.line);

    // Azimuth vectors
    let mut azimuth_vectors = Vec::with_capacity(raw.noise_azimuth_vector_list.vectors.len());
    for v in &raw.noise_azimuth_vector_list.vectors {
        let swath = parse_subswath_id(&v.swath)?;
        let first_azimuth_line = parse_i32(&v.first_azimuth_line, "firstAzimuthLine")?;
        let first_range_sample = parse_i32(&v.first_range_sample, "firstRangeSample")?;
        let last_azimuth_line = parse_i32(&v.last_azimuth_line, "lastAzimuthLine")?;
        let last_range_sample = parse_i32(&v.last_range_sample, "lastRangeSample")?;
        let lines = parse_space_i32(&v.line, "noiseAzimuthVector/line")?;
        let noise_azimuth_lut =
            parse_space_f32(&v.noise_azimuth_lut, "noiseAzimuthVector/noiseAzimuthLut")?;

        check_array_len("noiseAzimuthLut", lines.len(), noise_azimuth_lut.len())?;

        azimuth_vectors.push(NoiseAzimuthVector {
            swath,
            first_azimuth_line,
            first_range_sample,
            last_azimuth_line,
            last_range_sample,
            lines,
            noise_azimuth_lut,
        });
    }

    Ok(SwathNoise {
        subswath_id,
        polarization,
        range_vectors,
        azimuth_vectors,
    })
}

// ═══════════════════════════════════════════════════════════════════════
// Helpers
// ═══════════════════════════════════════════════════════════════════════

fn parse_timestamp(
    s: &str,
    field: &'static str,
) -> Result<DateTime<Utc>, CalibrationParseError> {
    let naive =
        NaiveDateTime::parse_from_str(s.trim(), "%Y-%m-%dT%H:%M:%S%.f").map_err(|e| {
            CalibrationParseError::InvalidValue {
                field,
                detail: format!("'{}': {}", s, e),
            }
        })?;
    Ok(naive.and_utc())
}

fn parse_f64(s: &str, field: &'static str) -> Result<f64, CalibrationParseError> {
    s.trim().parse::<f64>().map_err(|e| {
        CalibrationParseError::InvalidValue {
            field,
            detail: format!("'{}': {}", s, e),
        }
    })
}

fn parse_i32(s: &str, field: &'static str) -> Result<i32, CalibrationParseError> {
    s.trim().parse::<i32>().map_err(|e| {
        CalibrationParseError::InvalidValue {
            field,
            detail: format!("'{}': {}", s, e),
        }
    })
}

fn parse_space_u32(s: &str, field: &'static str) -> Result<Vec<u32>, CalibrationParseError> {
    s.split_whitespace()
        .map(|t| {
            t.parse::<u32>().map_err(|e| {
                CalibrationParseError::InvalidValue {
                    field,
                    detail: format!("'{}': {}", t, e),
                }
            })
        })
        .collect()
}

fn parse_space_i32(s: &str, field: &'static str) -> Result<Vec<i32>, CalibrationParseError> {
    s.split_whitespace()
        .map(|t| {
            t.parse::<i32>().map_err(|e| {
                CalibrationParseError::InvalidValue {
                    field,
                    detail: format!("'{}': {}", t, e),
                }
            })
        })
        .collect()
}

fn parse_space_f32(s: &str, field: &'static str) -> Result<Vec<f32>, CalibrationParseError> {
    s.split_whitespace()
        .map(|t| {
            t.parse::<f32>().map_err(|e| {
                CalibrationParseError::InvalidValue {
                    field,
                    detail: format!("'{}': {}", t, e),
                }
            })
        })
        .collect()
}

fn parse_polarization(s: &str) -> Result<Polarization, CalibrationParseError> {
    match s.trim() {
        "VV" => Ok(Polarization::VV),
        "VH" => Ok(Polarization::VH),
        "HV" => Ok(Polarization::HV),
        "HH" => Ok(Polarization::HH),
        other => Err(CalibrationParseError::InvalidValue {
            field: "polarisation",
            detail: other.to_string(),
        }),
    }
}

fn parse_subswath_id(s: &str) -> Result<SubSwathId, CalibrationParseError> {
    match s.trim() {
        "IW1" => Ok(SubSwathId::IW1),
        "IW2" => Ok(SubSwathId::IW2),
        "IW3" => Ok(SubSwathId::IW3),
        "EW1" => Ok(SubSwathId::EW1),
        "EW2" => Ok(SubSwathId::EW2),
        "EW3" => Ok(SubSwathId::EW3),
        "EW4" => Ok(SubSwathId::EW4),
        "EW5" => Ok(SubSwathId::EW5),
        other => Err(CalibrationParseError::InvalidValue {
            field: "swath",
            detail: other.to_string(),
        }),
    }
}

fn check_array_len(
    field: &'static str,
    expected: usize,
    actual: usize,
) -> Result<(), CalibrationParseError> {
    if expected != actual {
        Err(CalibrationParseError::ArrayLengthMismatch {
            field,
            expected,
            actual,
        })
    } else {
        Ok(())
    }
}

// ═══════════════════════════════════════════════════════════════════════
// Tests
// ═══════════════════════════════════════════════════════════════════════

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::PathBuf;

    fn s1a_safe() -> PathBuf {
        PathBuf::from("/home/datacube/dev/SARdine/data/SLC/S1A_IW_SLC__1SDV_20201005T170824_20201005T170851_034664_04098A_1E66.SAFE")
    }

    fn s1b_safe() -> PathBuf {
        PathBuf::from("/home/datacube/dev/SARdine/data/SLC/S1B_IW_SLC__1SDV_20190123T053348_20190123T053415_014617_01B3D4_E833.SAFE")
    }

    fn have_test_data() -> bool {
        s1a_safe()
            .join("annotation")
            .join("calibration")
            .is_dir()
    }

    // ── Unit tests: space-separated parsing ──

    #[test]
    fn test_parse_space_u32() {
        let vals = parse_space_u32("0 40 80 120 160", "test").unwrap();
        assert_eq!(vals, vec![0, 40, 80, 120, 160]);
    }

    #[test]
    fn test_parse_space_f32() {
        let vals = parse_space_f32("3.339e+02 2.370e+02 1.0", "test").unwrap();
        assert_eq!(vals.len(), 3);
        assert!((vals[0] - 333.9).abs() < 0.1);
        assert!((vals[1] - 237.0).abs() < 0.1);
        assert!((vals[2] - 1.0).abs() < 1e-6);
    }

    #[test]
    fn test_parse_space_i32() {
        let vals = parse_space_i32("-1498 0 1498 2996", "test").unwrap();
        assert_eq!(vals, vec![-1498, 0, 1498, 2996]);
    }

    #[test]
    fn test_parse_space_u32_empty() {
        let vals = parse_space_u32("", "test").unwrap();
        assert!(vals.is_empty());
    }

    // ── Integration: single calibration XML ──

    #[test]
    fn test_parse_single_calibration_xml() {
        if !have_test_data() {
            eprintln!("Skipping: test data not available");
            return;
        }

        let cal_dir = s1a_safe().join("annotation").join("calibration");
        let xml_path = std::fs::read_dir(&cal_dir)
            .unwrap()
            .filter_map(|e| e.ok())
            .find(|e| {
                let n = e.file_name().to_string_lossy().to_string();
                n.starts_with("calibration-") && n.contains("iw1") && n.contains("vv")
            })
            .expect("calibration IW1-VV not found")
            .path();

        let content = std::fs::read_to_string(&xml_path).unwrap();
        let cal = parse_calibration_xml(&content).unwrap();

        assert_eq!(cal.subswath_id, SubSwathId::IW1);
        assert_eq!(cal.polarization, Polarization::VV);
        assert!((cal.absolute_calibration_constant - 1.0).abs() < 1e-10);

        // S1A IW1 has 30 calibration vectors (1 per second over ~25s)
        assert!(!cal.vectors.is_empty());
        assert!(cal.vectors.len() >= 20);

        // Each vector should have 554 pixel samples
        let v0 = &cal.vectors[0];
        assert_eq!(v0.pixels.len(), 554);
        assert_eq!(v0.sigma_nought.len(), 554);
        assert_eq!(v0.beta_nought.len(), 554);
        assert_eq!(v0.gamma.len(), 554);
        assert_eq!(v0.dn.len(), 554);

        // Pixel sampling: starts at 0, sampled every 40
        assert_eq!(v0.pixels[0], 0);
        assert_eq!(v0.pixels[1], 40);

        // Values are in linear domain (~50-5000 for sigma, ~237 for beta)
        assert!(v0.sigma_nought[0] > 100.0 && v0.sigma_nought[0] < 1000.0);
        assert!((v0.beta_nought[0] - 237.0).abs() < 1.0);
        assert!(v0.gamma[0] > 100.0 && v0.gamma[0] < 1000.0);
        assert!((v0.dn[0] - 237.0).abs() < 1.0);

        // Lines are sorted and increasing
        for w in cal.vectors.windows(2) {
            assert!(w[1].line > w[0].line);
        }
    }

    // ── Integration: single noise XML ──

    #[test]
    fn test_parse_single_noise_xml() {
        if !have_test_data() {
            eprintln!("Skipping: test data not available");
            return;
        }

        let cal_dir = s1a_safe().join("annotation").join("calibration");
        let xml_path = std::fs::read_dir(&cal_dir)
            .unwrap()
            .filter_map(|e| e.ok())
            .find(|e| {
                let n = e.file_name().to_string_lossy().to_string();
                n.starts_with("noise-") && n.contains("iw1") && n.contains("vv")
            })
            .expect("noise IW1-VV not found")
            .path();

        let content = std::fs::read_to_string(&xml_path).unwrap();
        let noise = parse_noise_xml(&content).unwrap();

        assert_eq!(noise.subswath_id, SubSwathId::IW1);
        assert_eq!(noise.polarization, Polarization::VV);

        // Range vectors
        assert!(!noise.range_vectors.is_empty());
        assert!(noise.range_vectors.len() >= 8);

        let rv0 = &noise.range_vectors[0];
        assert_eq!(rv0.pixels.len(), rv0.noise_range_lut.len());
        // First range vector line can be negative (extrapolation)
        assert!(rv0.line < 0 || rv0.line == 0);
        // Noise values are positive (linear power)
        assert!(rv0.noise_range_lut.iter().all(|&v| v > 0.0));

        // Azimuth vectors
        assert!(!noise.azimuth_vectors.is_empty());
        let av0 = &noise.azimuth_vectors[0];
        assert_eq!(av0.swath, SubSwathId::IW1);
        assert!(av0.first_azimuth_line >= 0);
        assert!(av0.last_azimuth_line > av0.first_azimuth_line);
        assert_eq!(av0.lines.len(), av0.noise_azimuth_lut.len());
        // Azimuth noise is a multiplicative correction factor (~1.0)
        assert!(av0.noise_azimuth_lut.iter().all(|&v| v > 0.0 && v < 10.0));
    }

    // ── Integration: full SAFE calibration/noise ──

    #[test]
    fn test_parse_s1a_calibration_noise() {
        if !have_test_data() {
            eprintln!("Skipping: test data not available");
            return;
        }

        let data = parse_calibration_noise(&s1a_safe()).unwrap();

        // S1A dual-pol IW: 3 swaths × 2 polarizations = 6 calibration + 6 noise files
        assert_eq!(data.calibrations.len(), 6);
        assert_eq!(data.noises.len(), 6);

        // Verify all swath/pol combinations are present
        let cal_keys: Vec<_> = data
            .calibrations
            .iter()
            .map(|c| (c.subswath_id, c.polarization))
            .collect();
        assert!(cal_keys.contains(&(SubSwathId::IW1, Polarization::VH)));
        assert!(cal_keys.contains(&(SubSwathId::IW1, Polarization::VV)));
        assert!(cal_keys.contains(&(SubSwathId::IW2, Polarization::VH)));
        assert!(cal_keys.contains(&(SubSwathId::IW2, Polarization::VV)));
        assert!(cal_keys.contains(&(SubSwathId::IW3, Polarization::VH)));
        assert!(cal_keys.contains(&(SubSwathId::IW3, Polarization::VV)));

        let noise_keys: Vec<_> = data
            .noises
            .iter()
            .map(|n| (n.subswath_id, n.polarization))
            .collect();
        assert!(noise_keys.contains(&(SubSwathId::IW1, Polarization::VV)));
        assert!(noise_keys.contains(&(SubSwathId::IW3, Polarization::VH)));
    }

    #[test]
    fn test_parse_s1b_calibration_noise() {
        let path = s1b_safe();
        if !path.join("annotation").join("calibration").is_dir() {
            eprintln!("Skipping: S1B test data not available");
            return;
        }

        let data = parse_calibration_noise(&path).unwrap();

        assert_eq!(data.calibrations.len(), 6);
        assert_eq!(data.noises.len(), 6);

        // S1B should have the same modern noise format
        for noise in &data.noises {
            assert!(!noise.range_vectors.is_empty());
            assert!(!noise.azimuth_vectors.is_empty());
        }
    }

    /// Verify that calibration LUT values are consistently linear (not dB).
    ///
    /// ESA Sentinel-1 calibration vectors contain linear values in the range ~50-5000.
    /// If they were dB-encoded, they would be in the range ~17-37 dB.
    /// This test ensures we don't misinterpret the domain.
    #[test]
    fn test_calibration_values_are_linear() {
        if !have_test_data() {
            eprintln!("Skipping: test data not available");
            return;
        }

        let data = parse_calibration_noise(&s1a_safe()).unwrap();

        for cal in &data.calibrations {
            for vec in &cal.vectors {
                // sigma_nought: linear values should be > 50
                let sigma_min = vec
                    .sigma_nought
                    .iter()
                    .copied()
                    .fold(f32::INFINITY, f32::min);
                let sigma_max = vec
                    .sigma_nought
                    .iter()
                    .copied()
                    .fold(f32::NEG_INFINITY, f32::max);
                assert!(
                    sigma_min > 50.0,
                    "sigma0 too small for linear: {} (swath {:?})",
                    sigma_min,
                    cal.subswath_id
                );
                assert!(
                    sigma_max < 5000.0,
                    "sigma0 too large: {} (swath {:?})",
                    sigma_max,
                    cal.subswath_id
                );

                // gamma: should be similar range to sigma
                let gamma_min = vec
                    .gamma
                    .iter()
                    .copied()
                    .fold(f32::INFINITY, f32::min);
                assert!(
                    gamma_min > 50.0,
                    "gamma too small for linear: {} (swath {:?})",
                    gamma_min,
                    cal.subswath_id
                );
            }
        }
    }
}
