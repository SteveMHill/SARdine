// Comprehensive Sentinel-1 Annotation XML Parser
// Built with serde-xml-rs for robust space-separated value handling
// Covers all major annotation structures needed for SAR processing

use crate::constants::physical::SPEED_OF_LIGHT_M_S;
use crate::types::{SarError, SarResult, StateVector};
use chrono::{DateTime, NaiveDateTime, Utc};
use quick_xml::events::Event;
use quick_xml::Reader;
use regex::Regex;
use serde::{Deserialize, Deserializer};
use std::collections::HashMap; // For simple XML parsing

// ============================================================================
// TIME PARSING HELPERS - Handle different timestamp formats
// ============================================================================

/// Robust time parsing function that handles multiple timestamp formats
/// Common formats in Sentinel-1 XML:
/// - "2020-12-28T21:59:42+00:00" (with timezone, no microseconds)
/// - "2020-12-28T21:59:42.123456+00:00" (with timezone and microseconds)
/// - "2020-12-28T21:59:42Z" (UTC format)
/// - "2020-12-28T21:59:42.123456Z" (UTC with microseconds)
pub fn parse_time_robust(time_str: &str) -> Option<DateTime<Utc>> {
    // Try multiple common formats
    let formats = [
        "%Y-%m-%dT%H:%M:%S%.f%z", // With timezone and optional microseconds
        "%Y-%m-%dT%H:%M:%S%z",    // With timezone, no microseconds
        "%Y-%m-%dT%H:%M:%S%.fZ",  // UTC with optional microseconds
        "%Y-%m-%dT%H:%M:%SZ",     // UTC, no microseconds
        "%Y-%m-%dT%H:%M:%S%.f",   // No timezone, optional microseconds (assume UTC)
        "%Y-%m-%dT%H:%M:%S",      // No timezone, no microseconds (assume UTC)
    ];

    for format in &formats {
        if let Ok(dt) = DateTime::parse_from_str(time_str, format) {
            return Some(dt.with_timezone(&Utc));
        }
        // For formats without timezone, parse as naive and assume UTC
        if format.ends_with('f') || format.ends_with('S') {
            if let Ok(ndt) = NaiveDateTime::parse_from_str(time_str, format) {
                return Some(DateTime::from_naive_utc_and_offset(ndt, Utc));
            }
        }
    }

    // Fallback: try the basic chrono parser
    if let Ok(dt) = time_str.parse::<DateTime<Utc>>() {
        return Some(dt);
    }

    None
}

// ============================================================================
// XML PREPROCESSING HELPERS - Handle BOM and namespaces
// ============================================================================

fn strip_bom(s: &str) -> &str {
    const BOM: &str = "\u{FEFF}";
    s.strip_prefix(BOM).unwrap_or(s)
}

fn strip_xml_namespaces(xml: &str) -> String {
    // remove xmlns declarations
    let re_xmlns = Regex::new(r#"xmlns(:\w+)?="[^"]*""#).unwrap();
    let no_xmlns = re_xmlns.replace_all(xml, "");
    // drop element prefixes like <s1:tag> or </s1:tag>  ->  <tag> / </tag>
    let re_prefix = Regex::new(r#"<(/?)(\w+):"#).unwrap();
    re_prefix.replace_all(&no_xmlns, "<$1").to_string()
}

// ============================================================================
// XML EXTRACTION HELPERS - Handle different XML structures
// ============================================================================

/// Extract the content of a nested XML element (for SAFE format orbit data)
fn extract_element_content(xml: &str, tag: &str) -> Option<String> {
    let pattern = format!(r"(?s)<{}[^>]*>(.*?)</{}>", tag, tag);
    if let Ok(re) = Regex::new(&pattern) {
        if let Some(caps) = re.captures(xml) {
            return Some(caps[1].to_string());
        }
    }
    None
}

/// Extract a simple value from XML content
fn extract_value(xml: &str, tag: &str) -> Option<String> {
    let pattern = format!(r"<{}[^>]*>(.*?)</{}>", tag, tag);
    if let Ok(re) = Regex::new(&pattern) {
        if let Some(caps) = re.captures(xml) {
            return Some(caps[1].trim().to_string());
        }
    }
    None
}

/// Extract a float value from XML content
fn extract_float(xml: &str, tag: &str) -> Option<f64> {
    extract_value(xml, tag)?.parse().ok()
}

/// Extract nth occurrence of a float value from XML content
fn extract_float_nth(xml: &str, tag: &str, n: usize) -> Option<f64> {
    let pattern = format!(r"<{}[^>]*>(.*?)</{}>", tag, tag);
    if let Ok(re) = Regex::new(&pattern) {
        let matches: Vec<_> = re.captures_iter(xml).collect();
        if matches.len() >= n {
            return matches[n - 1][1].trim().parse().ok();
        }
    }
    None
}

/// Extract orbit state vectors from anywhere in XML content
/// Handles both SAFE (<orbit><position/><velocity/></orbit>) and ZIP (<stateVector>) formats
fn extract_orbits_anywhere(xml: &str) -> Vec<StateVector> {
    #[derive(Default, Debug)]
    struct PartialVector {
        time: Option<String>,
        position: [Option<f64>; 3],
        velocity: [Option<f64>; 3],
    }

    impl PartialVector {
        fn push_component(target: &mut [Option<f64>; 3], idx: usize, text: &str) {
            if let Ok(value) = text.trim().parse::<f64>() {
                target[idx] = Some(value);
            }
        }

        fn finalize(self) -> Option<StateVector> {
            let time = self.time.and_then(|t| parse_time_robust(&t));
            let position = match (self.position[0], self.position[1], self.position[2]) {
                (Some(x), Some(y), Some(z)) => [x, y, z],
                _ => return None,
            };
            let velocity = match (self.velocity[0], self.velocity[1], self.velocity[2]) {
                (Some(x), Some(y), Some(z)) => [x, y, z],
                _ => return None,
            };
            time.map(|time| StateVector {
                time,
                position,
                velocity,
            })
        }
    }

    #[derive(Copy, Clone, Debug, Eq, PartialEq)]
    enum Container {
        Orbit,
        StateVector,
    }

    #[derive(Copy, Clone, Debug, Eq, PartialEq)]
    enum Section {
        None,
        Position,
        Velocity,
    }

    let mut reader = Reader::from_str(xml);
    reader.trim_text(true);

    let mut buf = Vec::new();
    let mut container_stack: Vec<Container> = Vec::new();
    let mut current = Vec::new();
    let mut section = Section::None;

    #[derive(Copy, Clone, Debug, Eq, PartialEq)]
    enum ValueTarget {
        None,
        Time,
        PosX,
        PosY,
        PosZ,
        VelX,
        VelY,
        VelZ,
    }

    let mut target = ValueTarget::None;
    let mut results = Vec::new();

    loop {
        match reader.read_event_into(&mut buf) {
            Ok(Event::Start(e)) => {
                let name = e.local_name();
                match name.as_ref() {
                    b"orbit" => {
                        container_stack.push(Container::Orbit);
                        current.push(PartialVector::default());
                    }
                    b"stateVector" => {
                        container_stack.push(Container::StateVector);
                        current.push(PartialVector::default());
                    }
                    b"position" => {
                        section = Section::Position;
                    }
                    b"velocity" => {
                        section = Section::Velocity;
                    }
                    b"time" | b"utc" | b"utcTime" | b"utcDateTime" | b"azimuthTime" => {
                        target = ValueTarget::Time;
                    }
                    b"x" => {
                        target = match section {
                            Section::Position => ValueTarget::PosX,
                            Section::Velocity => ValueTarget::VelX,
                            Section::None => ValueTarget::PosX,
                        };
                    }
                    b"y" => {
                        target = match section {
                            Section::Position => ValueTarget::PosY,
                            Section::Velocity => ValueTarget::VelY,
                            Section::None => ValueTarget::PosY,
                        };
                    }
                    b"z" => {
                        target = match section {
                            Section::Position => ValueTarget::PosZ,
                            Section::Velocity => ValueTarget::VelZ,
                            Section::None => ValueTarget::PosZ,
                        };
                    }
                    b"vx" => target = ValueTarget::VelX,
                    b"vy" => target = ValueTarget::VelY,
                    b"vz" => target = ValueTarget::VelZ,
                    _ => {
                        target = ValueTarget::None;
                    }
                }
            }
            Ok(Event::Text(e)) => {
                if let Some(partial) = current.last_mut() {
                    let text = e
                        .unescape()
                        .map(|cow| cow.into_owned())
                        .unwrap_or_else(|_| String::from_utf8_lossy(e.as_ref()).into_owned());
                    match target {
                        ValueTarget::Time => {
                            partial.time = Some(text.trim().to_string());
                        }
                        ValueTarget::PosX => {
                            PartialVector::push_component(&mut partial.position, 0, &text);
                        }
                        ValueTarget::PosY => {
                            PartialVector::push_component(&mut partial.position, 1, &text);
                        }
                        ValueTarget::PosZ => {
                            PartialVector::push_component(&mut partial.position, 2, &text);
                        }
                        ValueTarget::VelX => {
                            PartialVector::push_component(&mut partial.velocity, 0, &text);
                        }
                        ValueTarget::VelY => {
                            PartialVector::push_component(&mut partial.velocity, 1, &text);
                        }
                        ValueTarget::VelZ => {
                            PartialVector::push_component(&mut partial.velocity, 2, &text);
                        }
                        ValueTarget::None => {}
                    }
                }
            }
            Ok(Event::End(e)) => {
                let name = e.local_name();
                match name.as_ref() {
                    b"position" | b"velocity" => {
                        section = Section::None;
                    }
                    b"orbit" | b"stateVector" => {
                        if let Some(container) = container_stack.pop() {
                            if let Some(partial) = current.pop() {
                                if let Some(vector) = partial.finalize() {
                                    results.push(vector);
                                } else {
                                    log::debug!(
                                        "🔍 ORBIT ENRICHMENT: Skipping incomplete {:?} entry",
                                        container
                                    );
                                }
                            } else {
                                log::warn!("⚠️ Missing state vector data during XML parsing - skipping incomplete entry");
                            }
                        }
                    }
                    _ => {}
                }
                target = ValueTarget::None;
            }
            Ok(Event::Eof) => break,
            Ok(_) => {}
            Err(e) => {
                log::warn!("🔍 ORBIT ENRICHMENT: quick-xml streaming error: {}", e);
                break;
            }
        }
        buf.clear();
    }

    results
}

// regex-based swath extraction removed (strict serde parsing only)

// ============================================================================
// CUSTOM DESERIALIZERS - The key to handling space-separated values
// ============================================================================

// Custom deserializer for space-separated f64 values
fn deserialize_space_separated_f64<'de, D>(deserializer: D) -> Result<Vec<f64>, D::Error>
where
    D: Deserializer<'de>,
{
    let s = String::deserialize(deserializer)?;
    let trimmed = s.trim();

    if trimmed.is_empty() {
        return Ok(Vec::new());
    }

    // Split by whitespace and parse each value
    let values: Result<Vec<f64>, _> = trimmed
        .split_whitespace()
        .filter(|part| !part.is_empty())
        .map(|part| {
            // Handle scientific notation like 5.336427094578738e-03
            part.parse::<f64>()
        })
        .collect();

    values.map_err(|e| {
        serde::de::Error::custom(format!("Failed to parse space-separated f64 values: {}", e))
    })
}

// Custom deserializer for single f64 that might contain spaces (takes first value)
fn deserialize_flexible_f64<'de, D>(deserializer: D) -> Result<f64, D::Error>
where
    D: Deserializer<'de>,
{
    let s = String::deserialize(deserializer)?;
    let trimmed = s.trim();

    // If it contains spaces, take the first number
    if trimmed.contains(' ') {
        let first = trimmed.split_whitespace().next().unwrap_or("");
        first.parse::<f64>().map_err(serde::de::Error::custom)
    } else {
        trimmed.parse::<f64>().map_err(serde::de::Error::custom)
    }
}

// Custom deserializer for optional flexible f64
fn deserialize_optional_flexible_f64<'de, D>(deserializer: D) -> Result<Option<f64>, D::Error>
where
    D: Deserializer<'de>,
{
    let opt = Option::<String>::deserialize(deserializer)?;
    match opt {
        Some(s) if !s.trim().is_empty() => {
            let trimmed = s.trim();
            // Treat common non-numeric placeholders as missing
            let lower = trimmed.to_ascii_lowercase();
            if matches!(lower.as_str(), "n/a" | "na" | "none" | "null") {
                return Ok(None);
            }
            // If it contains spaces, take the first token and try parsing
            let token = if trimmed.contains(' ') {
                trimmed.split_whitespace().next().unwrap_or("")
            } else {
                trimmed
            };
            match token.parse::<f64>() {
                Ok(v) => Ok(Some(v)),
                Err(_) => Ok(None),
            }
        }
        _ => Ok(None),
    }
}

// Custom deserializer for flexible u32 (handles optional values)
#[allow(dead_code)]
fn deserialize_flexible_u32<'de, D>(deserializer: D) -> Result<u32, D::Error>
where
    D: Deserializer<'de>,
{
    let s = String::deserialize(deserializer)?;
    s.trim().parse::<u32>().map_err(serde::de::Error::custom)
}

// Custom deserializer for optional space-separated f64 values
fn deserialize_optional_space_separated_f64<'de, D>(
    deserializer: D,
) -> Result<Option<Vec<f64>>, D::Error>
where
    D: Deserializer<'de>,
{
    let opt = Option::<String>::deserialize(deserializer)?;
    match opt {
        Some(s) if !s.trim().is_empty() => {
            let trimmed = s.trim();
            let lower = trimmed.to_ascii_lowercase();
            if matches!(lower.as_str(), "n/a" | "na" | "none" | "null") {
                return Ok(None);
            }
            let mut out = Vec::new();
            for part in trimmed.split_whitespace().filter(|p| !p.is_empty()) {
                if let Ok(v) = part.parse::<f64>() {
                    out.push(v);
                } else {
                    // On any non-numeric token, treat the entire field as missing
                    return Ok(None);
                }
            }
            Ok(Some(out))
        }
        _ => Ok(None),
    }
}

// Custom deserializer for optional u32
fn deserialize_optional_u32<'de, D>(deserializer: D) -> Result<Option<u32>, D::Error>
where
    D: Deserializer<'de>,
{
    let opt = Option::<String>::deserialize(deserializer)?;
    match opt {
        Some(s) if !s.trim().is_empty() => Ok(Some(
            s.trim().parse::<u32>().map_err(serde::de::Error::custom)?,
        )),
        _ => Ok(None),
    }
}

// Custom deserializer for optional u64 values (for large byte offsets)
fn deserialize_optional_u64<'de, D>(deserializer: D) -> Result<Option<u64>, D::Error>
where
    D: Deserializer<'de>,
{
    let opt = Option::<String>::deserialize(deserializer)?;
    match opt {
        Some(s) if !s.trim().is_empty() => Ok(Some(
            s.trim().parse::<u64>().map_err(serde::de::Error::custom)?,
        )),
        _ => Ok(None),
    }
}

// Custom deserializer for space-separated i32 values
fn deserialize_space_separated_i32<'de, D>(deserializer: D) -> Result<Vec<i32>, D::Error>
where
    D: Deserializer<'de>,
{
    let s = String::deserialize(deserializer)?;
    let values: Result<Vec<i32>, _> = s
        .split_whitespace()
        .filter(|s| !s.is_empty())
        .map(|s| s.parse::<i32>())
        .collect();

    values.map_err(serde::de::Error::custom)
}

// ============================================================================
// ROOT STRUCTURES
// ============================================================================

#[derive(Debug, Deserialize)]
pub struct ProductRoot {
    #[serde(rename = "adsHeader")]
    pub ads_header: Option<AdsHeader>,

    #[serde(rename = "qualityInformation")]
    pub quality_information: Option<QualityInformation>,

    #[serde(rename = "generalAnnotation")]
    pub general_annotation: Option<GeneralAnnotation>,

    #[serde(rename = "imageAnnotation")]
    pub image_annotation: Option<ImageAnnotation>,

    #[serde(rename = "dopplerCentroid")]
    pub doppler_centroid: Option<DopplerCentroid>,

    #[serde(rename = "antennaPattern")]
    pub antenna_pattern: Option<AntennaPatternList>,

    #[serde(rename = "swathTiming")]
    pub swath_timing: Option<SwathTiming>,

    #[serde(rename = "geolocationGrid")]
    pub geolocation_grid: Option<GeolocationGrid>,

    #[serde(rename = "coordinateConversion")]
    pub coordinate_conversion: Option<CoordinateConversion>,

    // swathMerging element is not required for our processing; ignore it

    // Additional fields for comprehensive metadata
    pub orbit_list: Option<Vec<StateVector>>,
}

// Legacy alias for compatibility with existing codebase
pub type AnnotationRoot = ProductRoot;

// ============================================================================
// HEADER INFORMATION
// ============================================================================

#[derive(Debug, Deserialize, Clone)]
pub struct AdsHeader {
    #[serde(rename = "missionId")]
    pub mission_id: Option<String>,

    #[serde(rename = "productType")]
    pub product_type: Option<String>,

    #[serde(rename = "polarisation")]
    pub polarisation: Option<String>,

    #[serde(rename = "mode")]
    pub mode: Option<String>,

    #[serde(rename = "swath")]
    pub swath: Option<String>,

    #[serde(rename = "startTime")]
    pub start_time: Option<String>,

    #[serde(rename = "stopTime")]
    pub stop_time: Option<String>,

    #[serde(rename = "absoluteOrbitNumber")]
    #[serde(deserialize_with = "deserialize_optional_u32")]
    pub absolute_orbit_number: Option<u32>,

    #[serde(rename = "missionDataTakeId")]
    #[serde(deserialize_with = "deserialize_optional_u32")]
    pub mission_data_take_id: Option<u32>,

    #[serde(rename = "imageNumber")]
    #[serde(deserialize_with = "deserialize_optional_u32")]
    pub image_number: Option<u32>,
}

// ============================================================================
// QUALITY INFORMATION
// ============================================================================

#[derive(Debug, Deserialize)]
pub struct QualityInformation {
    #[serde(rename = "productQualityIndex")]
    #[serde(deserialize_with = "deserialize_optional_flexible_f64")]
    pub product_quality_index: Option<f64>,

    #[serde(rename = "qualityDataList")]
    pub quality_data_list: Option<QualityDataList>,
}

#[derive(Debug, Deserialize)]
pub struct QualityDataList {
    #[serde(rename = "@count")]
    #[serde(deserialize_with = "deserialize_optional_u32")]
    pub count: Option<u32>,

    #[serde(rename = "qualityData")]
    pub quality_data: Option<Vec<QualityData>>,
}

#[derive(Debug, Deserialize)]
pub struct QualityData {
    #[serde(rename = "azimuthTime")]
    pub azimuth_time: Option<String>,

    #[serde(rename = "downlinkQuality")]
    pub downlink_quality: Option<DownlinkQuality>,

    #[serde(rename = "rawDataAnalysisQuality")]
    pub raw_data_analysis_quality: Option<RawDataAnalysisQuality>,

    #[serde(rename = "imageQuality")]
    pub image_quality: Option<ImageQuality>,
}

#[derive(Debug, Deserialize)]
pub struct DownlinkQuality {
    #[serde(rename = "iInputDataMean")]
    #[serde(deserialize_with = "deserialize_optional_flexible_f64")]
    pub i_input_data_mean: Option<f64>,

    #[serde(rename = "qInputDataMean")]
    #[serde(deserialize_with = "deserialize_optional_flexible_f64")]
    pub q_input_data_mean: Option<f64>,

    #[serde(rename = "iInputDataStdDev")]
    #[serde(deserialize_with = "deserialize_optional_flexible_f64")]
    pub i_input_data_std_dev: Option<f64>,

    #[serde(rename = "qInputDataStdDev")]
    #[serde(deserialize_with = "deserialize_optional_flexible_f64")]
    pub q_input_data_std_dev: Option<f64>,
}

#[derive(Debug, Deserialize)]
pub struct RawDataAnalysisQuality {
    #[serde(rename = "iBias")]
    #[serde(deserialize_with = "deserialize_optional_flexible_f64")]
    pub i_bias: Option<f64>,

    #[serde(rename = "qBias")]
    #[serde(deserialize_with = "deserialize_optional_flexible_f64")]
    pub q_bias: Option<f64>,

    #[serde(rename = "iqGainImbalance")]
    #[serde(deserialize_with = "deserialize_optional_flexible_f64")]
    pub iq_gain_imbalance: Option<f64>,

    #[serde(rename = "iqQuadratureDeparture")]
    #[serde(deserialize_with = "deserialize_optional_flexible_f64")]
    pub iq_quadrature_departure: Option<f64>,
}

#[derive(Debug, Deserialize)]
pub struct ImageQuality {
    #[serde(rename = "imageStatistics")]
    pub image_statistics: Option<ImageStatistics>,
}

#[derive(Debug, Deserialize)]
pub struct ImageStatistics {
    #[serde(rename = "outputDataMean")]
    pub output_data_mean: Option<ComplexValue>,

    #[serde(rename = "outputDataStdDev")]
    pub output_data_std_dev: Option<ComplexValue>,
}

#[derive(Debug, Deserialize, Clone)]
pub struct ComplexValue {
    #[serde(rename = "re")]
    #[serde(deserialize_with = "deserialize_flexible_f64")]
    pub re: f64,

    #[serde(rename = "im")]
    #[serde(deserialize_with = "deserialize_flexible_f64")]
    pub im: f64,
}

// Legacy alias for ComplexValue
pub type ComplexNumber = ComplexValue;

// ============================================================================
// GENERAL ANNOTATION - Core metadata
// ============================================================================

#[derive(Debug, Deserialize)]
pub struct GeneralAnnotation {
    #[serde(rename = "productInformation")]
    pub product_information: Option<ProductInformation>,

    #[serde(rename = "downlinkInformationList")]
    pub downlink_information_list: Option<DownlinkInformationList>,

    #[serde(rename = "antennaPattern")]
    pub antenna_pattern: Option<AntennaPatternList>,

    #[serde(rename = "azimuthFmRateList")]
    pub azimuth_fm_rate_list: Option<AzimuthFmRateList>,

    #[serde(rename = "dcEstimateList")]
    pub dc_estimate_list: Option<DcEstimateList>,
}

#[derive(Debug, Deserialize)]
pub struct ProductInformation {
    #[serde(rename = "pass")]
    pub pass: Option<String>, // Ascending/Descending

    #[serde(rename = "timelinessCategory")]
    pub timeliness_category: Option<String>,

    #[serde(rename = "platformHeading")]
    #[serde(deserialize_with = "deserialize_flexible_f64")]
    pub platform_heading: f64,

    #[serde(rename = "projection")]
    pub projection: Option<String>,

    #[serde(rename = "rangeSamplingRate")]
    #[serde(deserialize_with = "deserialize_flexible_f64")]
    pub range_sampling_rate: f64,

    #[serde(rename = "radarFrequency")]
    #[serde(deserialize_with = "deserialize_flexible_f64")]
    pub radar_frequency: f64,

    #[serde(rename = "azimuthSteeringRate")]
    #[serde(deserialize_with = "deserialize_flexible_f64")]
    pub azimuth_steering_rate: f64,

    // Critical for pixel spacing calculations
    #[serde(rename = "rangePixelSpacing")]
    #[serde(deserialize_with = "deserialize_optional_flexible_f64")]
    #[serde(default)]
    pub range_pixel_spacing: Option<f64>,

    #[serde(rename = "azimuthPixelSpacing")]
    #[serde(deserialize_with = "deserialize_optional_flexible_f64")]
    #[serde(default)]
    pub azimuth_pixel_spacing: Option<f64>,
}

#[derive(Debug, Deserialize)]
pub struct DownlinkInformationList {
    #[serde(rename = "@count")]
    #[serde(deserialize_with = "deserialize_optional_u32")]
    pub count: Option<u32>,

    #[serde(rename = "downlinkInformation")]
    pub downlink_information: Option<Vec<DownlinkInformation>>,
}

#[derive(Debug, Deserialize)]
pub struct DownlinkInformation {
    #[serde(rename = "swath")]
    pub swath: Option<String>,

    #[serde(rename = "azimuthTime")]
    pub azimuth_time: Option<String>,

    #[serde(rename = "firstLineSensingTime")]
    pub first_line_sensing_time: Option<String>,

    #[serde(rename = "lastLineSensingTime")]
    pub last_line_sensing_time: Option<String>,

    #[serde(rename = "prf")]
    #[serde(deserialize_with = "deserialize_flexible_f64")]
    pub prf: f64, // Pulse Repetition Frequency
}

// ============================================================================
// ANTENNA PATTERN - Space-separated values critical here!
// ============================================================================

#[derive(Debug, Deserialize, Clone)]
pub struct AntennaPatternList {
    #[serde(rename = "antennaPatternListHeader")]
    pub header: Option<AntennaPatternListHeader>,

    #[serde(rename = "antennaPattern")]
    pub antenna_patterns: Option<Vec<AntennaPatternValues>>,
}

// Legacy alias for compatibility
pub type AntennaPattern = AntennaPatternList;

#[derive(Debug, Deserialize, Clone)]
pub struct AntennaPatternListHeader {
    #[serde(rename = "count")]
    pub count: Option<u32>,
}

#[derive(Debug, Deserialize, Clone)]
pub struct AntennaPatternValues {
    #[serde(rename = "elevationAngle")]
    #[serde(deserialize_with = "deserialize_optional_space_separated_f64")]
    pub elevation_angle: Option<Vec<f64>>, // CRITICAL: Space-separated values for 3 IW subswaths

    #[serde(rename = "incidenceAngle")]
    #[serde(deserialize_with = "deserialize_optional_space_separated_f64")]
    pub incidence_angle: Option<Vec<f64>>, // CRITICAL: Space-separated values

    #[serde(rename = "slantRangeTime")]
    #[serde(deserialize_with = "deserialize_optional_space_separated_f64")]
    pub slant_range_time: Option<Vec<f64>>, // CRITICAL: Space-separated values

    #[serde(rename = "elevationPattern")]
    #[serde(deserialize_with = "deserialize_optional_space_separated_f64")]
    pub elevation_pattern: Option<Vec<f64>>,

    #[serde(rename = "terrainHeight")]
    #[serde(deserialize_with = "deserialize_optional_space_separated_f64")]
    pub terrain_height: Option<Vec<f64>>,

    #[serde(rename = "roll")]
    #[serde(deserialize_with = "deserialize_flexible_f64")]
    pub roll: f64,
}

// ============================================================================
// AZIMUTH FM RATE - Polynomial coefficients (space-separated)
// ============================================================================

#[derive(Debug, Deserialize)]
pub struct AzimuthFmRateList {
    #[serde(rename = "azimuthFmRate")]
    pub azimuth_fm_rates: Option<Vec<AzimuthFmRate>>,
}

#[derive(Debug, Deserialize)]
pub struct AzimuthFmRate {
    #[serde(rename = "azimuthTime")]
    pub azimuth_time: Option<String>,

    #[serde(rename = "t0")]
    #[serde(deserialize_with = "deserialize_flexible_f64")]
    pub t0: f64,

    #[serde(rename = "azimuthFmRatePolynomial")]
    #[serde(deserialize_with = "deserialize_space_separated_f64")]
    pub azimuth_fm_rate_polynomial: Vec<f64>, // CRITICAL: Polynomial coefficients
}

// ============================================================================
// DOPPLER CENTROID
// ============================================================================

#[derive(Debug, Deserialize)]
pub struct DopplerCentroid {
    #[serde(rename = "dcEstimateList")]
    pub dc_estimate_list: Option<DcEstimateList>,
}

#[derive(Debug, Deserialize)]
pub struct DcEstimateList {
    #[serde(rename = "@count")]
    #[serde(deserialize_with = "deserialize_optional_u32")]
    pub count: Option<u32>,

    #[serde(rename = "dcEstimate")]
    pub dc_estimates: Option<Vec<DcEstimate>>,
}

#[derive(Debug, Deserialize, Clone)]
pub struct DcEstimate {
    #[serde(rename = "azimuthTime")]
    pub azimuth_time: Option<String>,

    #[serde(rename = "t0")]
    #[serde(deserialize_with = "deserialize_flexible_f64")]
    pub t0: f64,

    #[serde(rename = "dataDcPolynomial")]
    #[serde(deserialize_with = "deserialize_space_separated_f64")]
    pub data_dc_polynomial: Vec<f64>, // CRITICAL: Polynomial coefficients

    #[serde(rename = "fineDceAzimuthStartTime")]
    pub fine_dce_azimuth_start_time: Option<String>,

    #[serde(rename = "fineDceAzimuthStopTime")]
    pub fine_dce_azimuth_stop_time: Option<String>,
}

// ============================================================================
// SWATH MERGE
// ============================================================================

#[derive(Debug, Deserialize)]
pub struct SwathMerging {
    #[serde(rename = "swathMergeList")]
    pub swath_merge_list: Option<SwathMergeList>,
}

#[derive(Debug, Deserialize)]
pub struct SwathMergeList {
    #[serde(rename = "@count")]
    #[serde(deserialize_with = "deserialize_optional_u32")]
    pub count: Option<u32>,

    #[serde(rename = "swathMerge")]
    pub swath_merges: Option<Vec<SwathMerge>>,
}

#[derive(Debug, Deserialize)]
pub struct SwathMerge {
    #[serde(rename = "swath")]
    pub swath: Option<String>,

    #[serde(rename = "azimuthTime")]
    pub azimuth_time: Option<String>,

    #[serde(rename = "firstAzimuthLine")]
    #[serde(deserialize_with = "deserialize_optional_u32")]
    pub first_azimuth_line: Option<u32>,

    #[serde(rename = "lastAzimuthLine")]
    #[serde(deserialize_with = "deserialize_optional_u32")]
    pub last_azimuth_line: Option<u32>,
}

// ============================================================================
// COORDINATE CONVERSION - Coordinate transformation information
// ============================================================================

#[derive(Debug, Deserialize)]
pub struct CoordinateConversion {
    #[serde(rename = "coordinateConversionList")]
    pub coordinate_conversion_list: Option<CoordinateConversionList>,
}

#[derive(Debug, Deserialize)]
pub struct CoordinateConversionList {
    #[serde(rename = "@count")]
    #[serde(deserialize_with = "deserialize_optional_u32")]
    pub count: Option<u32>,

    #[serde(rename = "coordinateConversion")]
    pub coordinate_conversions: Option<Vec<CoordinateConversionItem>>,
}

#[derive(Debug, Deserialize)]
pub struct CoordinateConversionItem {
    // Add fields as needed when we encounter actual coordinate conversion data
    // For now, this is just a placeholder since count=0 in our test data
}

// ============================================================================
// IMAGE ANNOTATION - Image-specific metadata
// ============================================================================

#[derive(Debug, Deserialize)]
pub struct ImageAnnotation {
    #[serde(rename = "imageInformation")]
    pub image_information: Option<ImageInformation>,

    #[serde(rename = "processingInformation")]
    pub processing_information: Option<ProcessingInformation>,
}

#[derive(Debug, Deserialize, Default)]
#[serde(default)]
pub struct ImageInformation {
    #[serde(rename = "slantRangeTime")]
    #[serde(deserialize_with = "deserialize_optional_flexible_f64")]
    #[serde(default)]
    pub slant_range_time: Option<f64>,

    #[serde(rename = "rangePixelSpacing")]
    #[serde(deserialize_with = "deserialize_optional_flexible_f64")]
    #[serde(default)]
    pub range_pixel_spacing: Option<f64>,

    #[serde(rename = "azimuthPixelSpacing")]
    #[serde(deserialize_with = "deserialize_optional_flexible_f64")]
    #[serde(default)]
    pub azimuth_pixel_spacing: Option<f64>,

    #[serde(rename = "numberOfSamples")]
    #[serde(deserialize_with = "deserialize_optional_u32")]
    pub number_of_samples: Option<u32>,

    #[serde(rename = "numberOfLines")]
    #[serde(deserialize_with = "deserialize_optional_u32")]
    pub number_of_lines: Option<u32>,

    #[serde(rename = "azimuthTimeInterval")]
    #[serde(deserialize_with = "deserialize_optional_flexible_f64")]
    #[serde(default)]
    pub azimuth_time_interval: Option<f64>,

    #[serde(rename = "rangeTimeInterval")]
    #[serde(deserialize_with = "deserialize_optional_flexible_f64")]
    #[serde(default)]
    pub range_time_interval: Option<f64>,

    #[serde(rename = "azimuthFrequency")]
    #[serde(deserialize_with = "deserialize_optional_flexible_f64")]
    #[serde(default)]
    pub azimuth_frequency: Option<f64>,

    // Additional fields from the original parser
    #[serde(rename = "productFirstLineUtcTime")]
    pub product_first_line_utc_time: Option<String>,

    #[serde(rename = "productLastLineUtcTime")]
    pub product_last_line_utc_time: Option<String>,

    #[serde(rename = "ascendingNodeTime")]
    pub ascending_node_time: Option<String>,

    #[serde(rename = "anchorTime")]
    pub anchor_time: Option<String>,

    #[serde(rename = "productComposition")]
    pub product_composition: Option<String>,

    #[serde(rename = "sliceNumber")]
    #[serde(default)]
    pub slice_number: Option<String>, // Changed to String to be more permissive

    #[serde(rename = "sliceList")]
    pub slice_list: Option<SliceList>,

    #[serde(rename = "pixelValue")]
    pub pixel_value: Option<String>,

    #[serde(rename = "outputPixels")]
    pub output_pixels: Option<String>,

    #[serde(rename = "zeroDopMinusAcqTime")]
    #[serde(deserialize_with = "deserialize_optional_flexible_f64")]
    pub zero_dop_minus_acq_time: Option<f64>,

    #[serde(rename = "incidenceAngleMidSwath")]
    #[serde(deserialize_with = "deserialize_optional_flexible_f64")]
    pub incidence_angle_mid_swath: Option<f64>,

    #[serde(rename = "imageStatistics")]
    pub image_statistics: Option<ImageStatistics>,
}

#[derive(Debug, Deserialize)]
pub struct ProcessingInformation {
    #[serde(rename = "swathProcParamsList")]
    pub swath_proc_params_list: Option<SwathProcParamsList>,
}

#[derive(Debug, Deserialize)]
pub struct SwathProcParamsList {
    #[serde(rename = "swathProcParams")]
    pub swath_proc_params: Option<Vec<SwathProcParams>>,
}

#[derive(Debug, Deserialize)]
pub struct SwathProcParams {
    #[serde(rename = "swath")]
    pub swath: Option<String>,

    #[serde(rename = "rangeProcessing")]
    pub range_processing: Option<RangeProcessing>,

    #[serde(rename = "azimuthProcessing")]
    pub azimuth_processing: Option<AzimuthProcessing>,
}

#[derive(Debug, Deserialize)]
pub struct RangeProcessing {
    #[serde(rename = "numberOfLooks")]
    #[serde(deserialize_with = "deserialize_optional_u32")]
    pub number_of_looks: Option<u32>,

    #[serde(rename = "lookBandwidth")]
    #[serde(deserialize_with = "deserialize_flexible_f64")]
    pub look_bandwidth: f64,

    #[serde(rename = "processingBandwidth")]
    #[serde(deserialize_with = "deserialize_flexible_f64")]
    pub processing_bandwidth: f64,
}

#[derive(Debug, Deserialize)]
pub struct AzimuthProcessing {
    #[serde(rename = "numberOfLooks")]
    #[serde(deserialize_with = "deserialize_optional_u32")]
    pub number_of_looks: Option<u32>,

    #[serde(rename = "lookBandwidth")]
    #[serde(deserialize_with = "deserialize_flexible_f64")]
    pub look_bandwidth: f64,

    #[serde(rename = "processingBandwidth")]
    #[serde(deserialize_with = "deserialize_flexible_f64")]
    pub processing_bandwidth: f64,
}

// ============================================================================
// SWATH TIMING
// ============================================================================

#[derive(Debug, Deserialize)]
pub struct SwathTiming {
    #[serde(rename = "linesPerBurst")]
    #[serde(deserialize_with = "deserialize_optional_u32")]
    pub lines_per_burst: Option<u32>,

    #[serde(rename = "samplesPerBurst")]
    #[serde(deserialize_with = "deserialize_optional_u32")]
    pub samples_per_burst: Option<u32>,

    #[serde(rename = "burstList")]
    pub burst_list: Option<BurstList>,
}

#[derive(Debug, Deserialize)]
pub struct BurstList {
    #[serde(rename = "@count")]
    #[serde(deserialize_with = "deserialize_optional_u32")]
    pub count: Option<u32>,

    #[serde(rename = "burst")]
    pub bursts: Option<Vec<Burst>>,
}

#[derive(Debug, Deserialize)]
pub struct Burst {
    #[serde(rename = "azimuthTime")]
    pub azimuth_time: Option<String>,

    #[serde(rename = "azimuthAnxTime")]
    #[serde(deserialize_with = "deserialize_flexible_f64")]
    pub azimuth_anx_time: f64,

    #[serde(rename = "sensingTime")]
    pub sensing_time: Option<String>,

    #[serde(rename = "byteOffset")]
    #[serde(deserialize_with = "deserialize_optional_u64")]
    pub byte_offset: Option<u64>,

    #[serde(rename = "firstValidSample")]
    #[serde(deserialize_with = "deserialize_space_separated_i32")]
    pub first_valid_sample: Vec<i32>, // CRITICAL: Space-separated values

    #[serde(rename = "lastValidSample")]
    #[serde(deserialize_with = "deserialize_space_separated_i32")]
    pub last_valid_sample: Vec<i32>, // CRITICAL: Space-separated values
}

#[derive(Debug, Clone)]
struct DerivedBurstGeometry {
    total_bursts: usize,
    valid_bursts: usize,
    lines_per_burst: Option<usize>,
    samples_per_burst: Option<usize>,
    total_lines: usize,
    first_valid_line: Option<usize>,
    last_valid_line_exclusive: Option<usize>,
    first_valid_sample: Option<usize>,
    last_valid_sample_exclusive: Option<usize>,
}

fn derive_burst_geometry(swath_timing: &SwathTiming) -> Option<DerivedBurstGeometry> {
    let bursts = swath_timing
        .burst_list
        .as_ref()
        .and_then(|list| list.bursts.as_ref())?;

    if bursts.is_empty() {
        return None;
    }

    let mut line_cursor = 0usize;
    let mut first_valid_line = None;
    let mut last_valid_line = None;
    let mut first_valid_sample = None;
    let mut last_valid_sample = None;
    let mut valid_bursts = 0usize;

    for burst in bursts {
        let lines_in_burst = burst
            .first_valid_sample
            .len()
            .max(burst.last_valid_sample.len());
        let mut burst_has_valid = false;

        for line_idx in 0..lines_in_burst {
            let first_sample = burst
                .first_valid_sample
                .get(line_idx)
                .copied()
                .unwrap_or(-1);
            let last_sample = burst.last_valid_sample.get(line_idx).copied().unwrap_or(-1);

            if first_sample >= 0 && last_sample >= first_sample {
                let global_line = line_cursor + line_idx;
                first_valid_line =
                    Some(first_valid_line.map_or(global_line, |v: usize| v.min(global_line)));
                last_valid_line = Some(
                    last_valid_line.map_or(global_line + 1, |v: usize| v.max(global_line + 1)),
                );

                let first_usize = first_sample as usize;
                let last_usize = (last_sample as usize).saturating_add(1);
                first_valid_sample =
                    Some(first_valid_sample.map_or(first_usize, |v: usize| v.min(first_usize)));
                last_valid_sample =
                    Some(last_valid_sample.map_or(last_usize, |v: usize| v.max(last_usize)));
                burst_has_valid = true;
            }
        }

        if burst_has_valid {
            valid_bursts += 1;
        }

        line_cursor += lines_in_burst;
    }

    let lines_per_burst = swath_timing
        .lines_per_burst
        .map(|v| v as usize)
        .or_else(|| {
            let len = bursts[0]
                .first_valid_sample
                .len()
                .max(bursts[0].last_valid_sample.len());
            if len > 0 {
                Some(len)
            } else {
                None
            }
        });

    let samples_per_burst = swath_timing.samples_per_burst.map(|v| v as usize);

    Some(DerivedBurstGeometry {
        total_bursts: bursts.len(),
        valid_bursts,
        lines_per_burst,
        samples_per_burst,
        total_lines: line_cursor,
        first_valid_line,
        last_valid_line_exclusive: last_valid_line,
        first_valid_sample,
        last_valid_sample_exclusive: last_valid_sample,
    })
}

// ============================================================================
// GEOLOCATION GRID
// ============================================================================

#[derive(Debug, Deserialize)]
pub struct GeolocationGrid {
    #[serde(rename = "geolocationGridPointList")]
    pub geolocation_grid_point_list: Option<GeolocationGridPointList>,
}

#[derive(Debug, Deserialize)]
pub struct GeolocationGridPointList {
    #[serde(rename = "@count")]
    #[serde(deserialize_with = "deserialize_optional_u32")]
    pub count: Option<u32>,

    #[serde(rename = "geolocationGridPoint")]
    pub geolocation_grid_points: Option<Vec<GeolocationGridPoint>>,
}

#[derive(Debug, Deserialize)]
pub struct GeolocationGridPoint {
    #[serde(rename = "azimuthTime")]
    pub azimuth_time: Option<String>,

    #[serde(rename = "slantRangeTime")]
    #[serde(deserialize_with = "deserialize_flexible_f64")]
    pub slant_range_time: f64,

    #[serde(rename = "line")]
    #[serde(deserialize_with = "deserialize_flexible_f64")]
    pub line: f64,

    #[serde(rename = "pixel")]
    #[serde(deserialize_with = "deserialize_flexible_f64")]
    pub pixel: f64,

    #[serde(rename = "latitude")]
    #[serde(deserialize_with = "deserialize_flexible_f64")]
    pub latitude: f64,

    #[serde(rename = "longitude")]
    #[serde(deserialize_with = "deserialize_flexible_f64")]
    pub longitude: f64,

    #[serde(rename = "height")]
    #[serde(deserialize_with = "deserialize_flexible_f64")]
    pub height: f64,

    #[serde(rename = "incidenceAngle")]
    #[serde(deserialize_with = "deserialize_flexible_f64")]
    pub incidence_angle: f64,

    #[serde(rename = "elevationAngle")]
    #[serde(deserialize_with = "deserialize_flexible_f64")]
    pub elevation_angle: f64,
}

// ============================================================================
// LEGACY COMPATIBILITY STRUCTURES
// ============================================================================

#[derive(Debug, Deserialize)]
pub struct SliceList {
    #[serde(rename = "@count")]
    #[serde(deserialize_with = "deserialize_optional_u32")]
    pub count: Option<u32>,

    #[serde(rename = "slice")]
    pub slices: Option<Vec<Slice>>,
}

#[derive(Debug, Deserialize)]
pub struct Slice {
    #[serde(rename = "sliceNumber")]
    #[serde(default)]
    pub slice_number: Option<String>, // Changed to String to be more permissive

    #[serde(rename = "sensingStartTime")]
    pub sensing_start_time: Option<String>,

    #[serde(rename = "sensingStopTime")]
    pub sensing_stop_time: Option<String>,
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
    pub first_valid_sample: usize,
    pub last_valid_sample: usize,
    pub first_valid_line: usize,
    pub last_valid_line: usize,
    pub slant_range_time: f64,
    pub range_sampling_rate: f64,
}

// Flexible XML element for catch-all fields
#[derive(Debug, Deserialize, Clone)]
pub struct FlexibleXmlElement {
    #[serde(rename = "re")]
    pub re: Option<f64>,
    #[serde(rename = "im")]
    pub im: Option<f64>,
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
// MAIN PARSER FUNCTION
// ============================================================================

/// Parse Sentinel-1 annotation XML using SIMPLE regex approach
/// NO MORE COMPLEX DESERIALIZERS - just extract what we need!
pub fn parse_annotation_xml(xml_content: &str) -> Result<ProductRoot, Box<dyn std::error::Error>> {
    log::debug!("🔍 ANNOTATION PARSING: Starting annotation XML parsing");

    // 1) normalize BOM/newlines
    let mut s = strip_bom(xml_content).replace("\r\n", "\n");
    // 2) strip namespaces so our plain tag names match
    s = strip_xml_namespaces(&s);

    // Env toggle: force serde-only and fail fast (no regex fallback)
    let serde_only = std::env::var("SARDINE_SERDE_ONLY")
        .map(|v| v == "1" || v.eq_ignore_ascii_case("true"))
        .unwrap_or(false);
    let serde_require_subswaths = std::env::var("SARDINE_REQUIRE_SUBSWATHS")
        .map(|v| v == "1" || v.eq_ignore_ascii_case("true"))
        .unwrap_or(serde_only); // default to strict when serde_only enabled

    // 3) try serde parsing
    match quick_xml::de::from_str::<ProductRoot>(&s) {
        Ok(mut ok) => {
            log::debug!(
                "🔍 ANNOTATION PARSING: Successfully parsed with serde - checking orbit data"
            );
            // Optional strict validation when requested
            if serde_require_subswaths {
                let has_params = ok
                    .image_annotation
                    .as_ref()
                    .and_then(|ia| ia.processing_information.as_ref())
                    .and_then(|pi| pi.swath_proc_params_list.as_ref())
                    .and_then(|pl| pl.swath_proc_params.as_ref())
                    .map(|v| !v.is_empty())
                    .unwrap_or(false);
                if !has_params {
                    return Err("Serde parse succeeded but missing processingInformation/swathProcParamsList (strict mode)".into());
                }
            }
            // No regex injection; trust serde-only parsing for swathProcParamsList

            let need_orbits = ok.orbit_list.as_ref().map(|v| v.is_empty()).unwrap_or(true);
            if need_orbits {
                let orbits = extract_orbits_anywhere(&s);
                if !orbits.is_empty() {
                    log::info!(
                        "🔍 ORBIT ENRICHMENT: Filling orbit_list with {} vectors",
                        orbits.len()
                    );
                    ok.orbit_list = Some(orbits);
                } else {
                    log::warn!("🔍 ORBIT ENRICHMENT: No SAFE/ZIP orbits found in XML");
                }
            } else if let Some(ref orbit_list) = ok.orbit_list {
                log::debug!(
                    "🔍 ANNOTATION PARSING: Serde found {} orbit vectors",
                    orbit_list.len()
                );
            }
            Ok(ok)
        }
        Err(e) => {
            // Fail fast: regex-based parser removed by policy
            return Err(format!("quick_xml error: {e}").into());
        }
    }
}

// Legacy alias for the main parsing function
pub fn parse_annotation(xml_content: &str) -> SarResult<ProductRoot> {
    // Use the unified serde-based parser
    parse_annotation_xml(xml_content).map_err(|e| SarError::XmlParsing(e.to_string()))
}

// ============================================================================
// UTILITY FUNCTIONS FOR COMPATIBILITY WITH EXISTING CODE
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
            if let Some(pattern_vec) = &list.antenna_patterns {
                for pattern in pattern_vec {
                    // CRITICAL SCIENTIFIC FIX: No fallback values for antenna parameters
                    // Missing antenna pattern data indicates corrupted annotation file
                    let elev = pattern.elevation_angle.as_ref().ok_or_else(|| {
                        SarError::Processing("Missing elevation angle in antenna pattern - invalid annotation".to_string())
                    })?;
                    let inc = pattern.incidence_angle.as_ref().ok_or_else(|| {
                        SarError::Processing("Missing incidence angle in antenna pattern - invalid annotation".to_string())
                    })?;
                    let srt = pattern.slant_range_time.as_ref().ok_or_else(|| {
                        SarError::Processing("Missing slant range time in antenna pattern - invalid annotation".to_string())
                    })?;
                    patterns.push((elev.clone(), inc.clone(), srt.clone()));
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
        let num_samples = subswath_data.range_samples as usize;
        let num_lines = subswath_data.azimuth_samples as usize;
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
        let near_range = 0.5 * SPEED_OF_LIGHT_M_S as f64 * near_srt;
        let far_range = 0.5 * SPEED_OF_LIGHT_M_S as f64 * far_srt;

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
        let mut last_valid_sample = None::<usize>;
        if let Some(st) = &self.swath_timing {
            if let Some(bl) = &st.burst_list {
                if let Some(bursts) = &bl.bursts {
                    for b in bursts {
                        if !b.first_valid_sample.is_empty() {
                            let min_first = b
                                .first_valid_sample
                                .iter()
                                .cloned()
                                .min()
                                .unwrap_or(0)
                                .max(0) as usize;
                            first_valid_sample =
                                Some(first_valid_sample.map_or(min_first, |v| v.min(min_first)));
                        }
                        if !b.last_valid_sample.is_empty() {
                            let max_last = b
                                .last_valid_sample
                                .iter()
                                .cloned()
                                .max()
                                .unwrap_or(0)
                                .max(0) as usize;
                            last_valid_sample =
                                Some(last_valid_sample.map_or(max_last, |v| v.max(max_last)));
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
        let last_valid_sample = last_valid_sample
            .ok_or_else(|| {
                SarError::Metadata(format!("Missing last valid sample data for {}", subswath))
            })?
            .min(num_samples);

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
            last_valid_sample,
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

                        let (range_samples_opt, az_lines_opt) =
                            if let Some(ref img_info) = image_annotation.image_information {
                                (
                                    img_info.number_of_samples.map(|v| v as usize),
                                    img_info.number_of_lines.map(|v| v as usize),
                                )
                            } else {
                                (None, None)
                            };

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
                            .and_then(|st| st.lines_per_burst.map(|v| v as usize));

                        let fallback_samples_per_burst = annotation
                            .swath_timing
                            .as_ref()
                            .and_then(|st| st.samples_per_burst.map(|v| v as usize));

                        let prf_opt = annotation.get_pulse_repetition_frequency();

                        for (i, params) in params_vec.iter().enumerate() {
                            log::info!("🔍 DEBUGGING: Processing swath_proc_params entry {}", i);
                            if let Some(ref swath_id) = params.swath {
                                log::info!("✅ DEBUGGING: Found swath_id: {}", swath_id);

                                // STRICT MODE: No fallbacks allowed for critical parameters
                                let range_pixel_spacing = range_ps_opt.ok_or_else(|| {
                                    SarError::Metadata(format!(
                                        "Missing range_pixel_spacing for swath {}",
                                        swath_id
                                    ))
                                })?;
                                let azimuth_pixel_spacing = az_ps_opt.ok_or_else(|| {
                                    SarError::Metadata(format!(
                                        "Missing azimuth_pixel_spacing for swath {}",
                                        swath_id
                                    ))
                                })?;
                                let slant_range_time = srt_opt.ok_or_else(|| {
                                    SarError::Metadata(format!(
                                        "Missing slant_range_time for swath {}",
                                        swath_id
                                    ))
                                })?;
                                let mut range_samples = range_samples_opt.ok_or_else(|| {
                                    SarError::Metadata(format!(
                                        "Missing range_samples for swath {}",
                                        swath_id
                                    ))
                                })?;
                                let mut azimuth_samples = az_lines_opt.ok_or_else(|| {
                                    SarError::Metadata(format!(
                                        "Missing azimuth_samples for swath {}",
                                        swath_id
                                    ))
                                })?;

                                let original_range_samples = range_samples;
                                let original_azimuth_samples = azimuth_samples;

                                let mut first_line_global = 0usize;
                                let mut last_line_global = original_azimuth_samples;
                                let mut first_sample_global = 0usize;
                                let mut last_sample_global = original_range_samples;
                                let mut burst_count = fallback_burst_count;
                                let mut burst_duration = 0.0f64;

                                if let Some(ref geom) = burst_geometry {
                                    burst_count = if geom.valid_bursts > 0 {
                                        geom.valid_bursts
                                    } else {
                                        geom.total_bursts
                                    };

                                    let derived_first_line = geom
                                        .first_valid_line
                                        .unwrap_or(0)
                                        .min(original_azimuth_samples);
                                    let derived_last_line = geom
                                        .last_valid_line_exclusive
                                        .or_else(|| Some(geom.total_lines))
                                        .unwrap_or(original_azimuth_samples)
                                        .min(original_azimuth_samples);

                                    first_line_global = derived_first_line;
                                    last_line_global = if derived_last_line > derived_first_line {
                                        derived_last_line
                                    } else {
                                        original_azimuth_samples
                                    };

                                    let mut extent_lines =
                                        last_line_global.saturating_sub(first_line_global);
                                    if extent_lines == 0 || extent_lines > original_azimuth_samples
                                    {
                                        first_line_global = 0;
                                        last_line_global = original_azimuth_samples;
                                        extent_lines = original_azimuth_samples;
                                    }
                                    azimuth_samples = extent_lines;

                                    let derived_first_sample = geom
                                        .first_valid_sample
                                        .unwrap_or(0)
                                        .min(original_range_samples);
                                    let derived_last_sample = geom
                                        .last_valid_sample_exclusive
                                        .or(geom.samples_per_burst)
                                        .or(fallback_samples_per_burst)
                                        .unwrap_or(original_range_samples)
                                        .min(original_range_samples);

                                    first_sample_global = derived_first_sample;
                                    last_sample_global =
                                        if derived_last_sample > derived_first_sample {
                                            derived_last_sample
                                        } else {
                                            original_range_samples
                                        };

                                    let mut extent_samples =
                                        last_sample_global.saturating_sub(first_sample_global);
                                    if extent_samples == 0
                                        || extent_samples > original_range_samples
                                    {
                                        first_sample_global = 0;
                                        last_sample_global = original_range_samples;
                                        extent_samples = original_range_samples;
                                    }
                                    range_samples = extent_samples;

                                    if let (Some(lines_per_burst), Some(prf)) =
                                        (geom.lines_per_burst, prf_opt)
                                    {
                                        if prf > 0.0 {
                                            burst_duration = lines_per_burst as f64 / prf;
                                        }
                                    }
                                } else if let (Some(lines_per_burst), Some(prf)) =
                                    (fallback_lines_per_burst, prf_opt)
                                {
                                    if prf > 0.0 {
                                        burst_duration = lines_per_burst as f64 / prf;
                                    }
                                }

                                if burst_count == 0 {
                                    burst_count = fallback_burst_count;
                                }

                                if azimuth_samples == 0 {
                                    azimuth_samples = original_azimuth_samples;
                                }
                                if range_samples == 0 {
                                    range_samples = original_range_samples;
                                }

                                last_line_global = (first_line_global + azimuth_samples)
                                    .min(original_azimuth_samples);
                                last_sample_global = (first_sample_global + range_samples)
                                    .min(original_range_samples);

                                log::info!(
                                    "✅ DEBUGGING: Subswath {} -> lines {}..{} ({}), samples {}..{} ({})",
                                    swath_id,
                                    first_line_global,
                                    last_line_global,
                                    azimuth_samples,
                                    first_sample_global,
                                    last_sample_global,
                                    range_samples
                                );

                                let subswath = crate::types::SubSwath {
                                    id: swath_id.clone(),
                                    burst_count,
                                    range_samples,
                                    azimuth_samples,
                                    first_line_global,
                                    last_line_global,
                                    first_sample_global,
                                    last_sample_global,
                                    range_pixel_spacing,
                                    azimuth_pixel_spacing,
                                    slant_range_time,
                                    burst_duration,
                                    prf_hz: annotation
                                        .general_annotation
                                        .as_ref()
                                        .and_then(|ga| ga.downlink_information_list.as_ref())
                                        .and_then(|dl| dl.downlink_information.as_ref())
                                        .and_then(|vec| vec.first())
                                        .map(|di| di.prf),
                                    dc_polynomial: None, // DC polynomial extracted during deburst
                                    azimuth_time_interval: None, // Timing extracted during deburst
                                };
                                subswaths.insert(swath_id.clone(), subswath);
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

    pub fn get_pulse_repetition_frequency(&self) -> Option<f64> {
        self.general_annotation
            .as_ref()
            .and_then(|ga| ga.downlink_information_list.as_ref())
            .and_then(|dl| dl.downlink_information.as_ref())
            .and_then(|vec| vec.first())
            .map(|di| di.prf)
    }

    pub fn get_radar_frequency_hz(&self) -> Option<f64> {
        self.general_annotation
            .as_ref()
            .and_then(|ga| ga.product_information.as_ref())
            .map(|pi| pi.radar_frequency)
    }

    pub fn extract_range_doppler_params(
        &self,
    ) -> SarResult<crate::core::terrain_correction::RangeDopplerParams> {
        let (range_ps, az_ps) = self
            .get_pixel_spacing()
            .ok_or_else(|| SarError::Metadata("Missing range/azimuth pixel spacing".to_string()))?;
        let srt = self
            .get_slant_range_time()
            .ok_or_else(|| SarError::Metadata("Missing slantRangeTime".to_string()))?;
        let prf = self
            .get_pulse_repetition_frequency()
            .ok_or_else(|| SarError::Metadata("Missing PRF in downlinkInformation".to_string()))?;
        let radar_freq = self.get_radar_frequency_hz().ok_or_else(|| {
            SarError::Metadata("Missing radarFrequency in productInformation".to_string())
        })?;
        let wavelength = SPEED_OF_LIGHT_M_S as f64 / radar_freq;

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

        // Parse product start time (grid epoch)
        let product_start_time_abs = {
            // Prefer imageInformation productFirstLineUtcTime, fallback to ads_header.start_time
            let start_opt = self
                .image_annotation
                .as_ref()
                .and_then(|ia| ia.image_information.as_ref())
                .and_then(|ii| ii.product_first_line_utc_time.clone())
                .or_else(|| self.ads_header.as_ref().and_then(|h| h.start_time.clone()));
            let dt = start_opt
                .and_then(|s| parse_time_robust(&s))
                .ok_or_else(|| {
                    SarError::Metadata(
                        "Missing product start time for Doppler centroid time base".to_string(),
                    )
                })?;
            let time_seconds = (dt.timestamp() as f64) + (dt.timestamp_subsec_nanos() as f64) * 1e-9;
            
            // CRITICAL VALIDATION: Ensure product start time is reasonable (2000-2100)
            const Y2000: f64 = 946684800.0;  // 2000-01-01 00:00:00 UTC
            const Y2100: f64 = 4102444800.0; // 2100-01-01 00:00:00 UTC
            if time_seconds < Y2000 || time_seconds > Y2100 {
                return Err(SarError::Metadata(format!(
                    "Product start time {:.3} ({}) outside reasonable range [2000-2100]",
                    time_seconds, dt
                )));
            }
            
            log::debug!("✅ Product start time: {} ({:.3}s since Unix epoch)", dt, time_seconds);
            time_seconds
        };

        // Parse product stop time (last line)
        let product_stop_time_abs = {
            let stop_opt = self
                .image_annotation
                .as_ref()
                .and_then(|ia| ia.image_information.as_ref())
                .and_then(|ii| ii.product_last_line_utc_time.clone())
                .or_else(|| self.ads_header.as_ref().and_then(|h| h.stop_time.clone()));
            let dt = stop_opt.and_then(|s| parse_time_robust(&s));
            dt.map(|d| (d.timestamp() as f64) + (d.timestamp_subsec_nanos() as f64) * 1e-9)
                .unwrap_or(product_start_time_abs) // fallback: zero duration if missing
        };

        let product_duration = (product_stop_time_abs - product_start_time_abs).max(0.0);

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
        let burst_geometry = self.swath_timing
            .as_ref()
            .and_then(|st| crate::io::annotation::derive_burst_geometry(st));
        
        let (first_valid_line, last_valid_line, first_valid_sample, last_valid_sample) = 
            if let Some(geom) = burst_geometry {
                (geom.first_valid_line, 
                 geom.last_valid_line_exclusive.map(|x| x.saturating_sub(1)), // Convert exclusive to inclusive
                 geom.first_valid_sample,
                 geom.last_valid_sample_exclusive.map(|x| x.saturating_sub(1))) // Convert exclusive to inclusive
            } else {
                (None, None, None, None)
            };
        
        Ok(crate::core::terrain_correction::RangeDopplerParams {
            range_pixel_spacing: range_ps,
            azimuth_pixel_spacing: az_ps,
            slant_range_time: srt,
            prf,
            azimuth_time_interval,
            wavelength,
            speed_of_light: SPEED_OF_LIGHT_M_S as f64,
            product_start_time_abs,
            product_stop_time_abs,
            product_duration,
            total_azimuth_lines,
            doppler_centroid,
            first_valid_line,
            last_valid_line,
            first_valid_sample,
            last_valid_sample,
        })
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
