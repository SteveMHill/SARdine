#![allow(dead_code)]
// Comprehensive Sentinel-1 Annotation XML Parser
// Built with serde-xml-rs for robust space-separated value handling
// Covers all major annotation structures needed for SAR processing

use crate::types::{SarError, SarResult, StateVector};
use chrono::{DateTime, NaiveDateTime, Utc};
use quick_xml::events::Event;
use quick_xml::Reader;
use regex::Regex;
use serde::{Deserialize, Deserializer};

// ============================================================================
// TIME PARSING HELPERS - Handle different timestamp formats
// ============================================================================

/// Robust time parsing function that handles multiple timestamp formats
/// Common formats in Sentinel-1 XML:
/// - "2020-12-28T21:59:42+00:00" (with timezone, no microseconds)
/// - "2020-12-28T21:59:42.123456+00:00" (with timezone and microseconds)
/// - "2020-12-28T21:59:42Z" (UTC format)
/// - "2020-12-28T21:59:42.123456Z" (UTC with microseconds)
///
/// FIXED: Deterministic parsing without brittle ends_with checks
pub fn parse_time_robust(time_str: &str) -> Option<DateTime<Utc>> {
    use chrono::TimeZone;

    // With offset (Z or ±hh:mm), with/without fractional seconds
    const WITH_TZ: &[&str] = &[
        "%Y-%m-%dT%H:%M:%S%.f%:z",
        "%Y-%m-%dT%H:%M:%S%:z",
        "%Y-%m-%dT%H:%M:%S%.fZ",
        "%Y-%m-%dT%H:%M:%SZ",
    ];
    for fmt in WITH_TZ {
        if let Ok(dt) = DateTime::parse_from_str(time_str, fmt) {
            return Some(dt.with_timezone(&Utc));
        }
    }

    // Naive → assume UTC (deterministic conversion)
    const NAIVE: &[&str] = &["%Y-%m-%dT%H:%M:%S%.f", "%Y-%m-%dT%H:%M:%S"];
    for fmt in NAIVE {
        if let Ok(ndt) = NaiveDateTime::parse_from_str(time_str, fmt) {
            return Some(Utc.from_utc_datetime(&ndt));
        }
    }

    // Last resort: try direct parse
    time_str.parse::<DateTime<Utc>>().ok()
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

// Custom deserializer for optional usize values
fn deserialize_optional_usize<'de, D>(deserializer: D) -> Result<Option<usize>, D::Error>
where
    D: Deserializer<'de>,
{
    let opt = Option::<String>::deserialize(deserializer)?;
    match opt {
        Some(s) if !s.trim().is_empty() => Ok(Some(
            s.trim()
                .parse::<usize>()
                .map_err(serde::de::Error::custom)?,
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
// DERIVED HELPERS ON PARSED ANNOTATION
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
    #[serde(default)]
    #[serde(deserialize_with = "deserialize_optional_u32")]
    pub absolute_orbit_number: Option<u32>,

    #[serde(rename = "missionDataTakeId")]
    #[serde(default)]
    #[serde(deserialize_with = "deserialize_optional_u32")]
    pub mission_data_take_id: Option<u32>,

    #[serde(rename = "imageNumber")]
    #[serde(default)]
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
    #[serde(default)]
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
    // BURST GEOMETRY DERIVATION
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
    #[serde(deserialize_with = "deserialize_optional_flexible_f64")]
    #[serde(default)]
    pub platform_heading: Option<f64>,

    #[serde(rename = "projection")]
    pub projection: Option<String>,

    #[serde(rename = "rangeSamplingRate")]
    #[serde(deserialize_with = "deserialize_flexible_f64")]
    pub range_sampling_rate: f64,

    #[serde(rename = "radarFrequency")]
    #[serde(deserialize_with = "deserialize_flexible_f64")]
    pub radar_frequency: f64,

    #[serde(rename = "azimuthSteeringRate")]
    #[serde(deserialize_with = "deserialize_optional_flexible_f64")]
    #[serde(default)]
    pub azimuth_steering_rate: Option<f64>,

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
    #[serde(default)]
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
    #[serde(deserialize_with = "deserialize_optional_flexible_f64")]
    #[serde(default)]
    pub prf: Option<f64>, // Pulse Repetition Frequency
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

    #[serde(rename = "antennaPatternList")]
    pub nested_list: Option<AntennaPatternNestedList>,
}

// Legacy alias for compatibility
pub type AntennaPattern = AntennaPatternList;

#[derive(Debug, Deserialize, Clone)]
pub struct AntennaPatternNestedList {
    #[serde(rename = "@count")]
    #[serde(default)]
    #[serde(deserialize_with = "deserialize_optional_u32")]
    pub count: Option<u32>,

    #[serde(rename = "antennaPattern")]
    pub antenna_patterns: Option<Vec<AntennaPatternValues>>,
}

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
    #[serde(default)]
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
    #[serde(default)]
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
    #[serde(default)]
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

    #[serde(rename = "rangeSamplingRate")]
    #[serde(deserialize_with = "deserialize_optional_flexible_f64")]
    #[serde(default)]
    pub range_sampling_rate: Option<f64>,

    #[serde(rename = "performedRangeCellSpacing")]
    #[serde(deserialize_with = "deserialize_optional_flexible_f64")]
    #[serde(default)]
    pub performed_range_cell_spacing: Option<f64>,

    #[serde(rename = "azimuthTimeInterval")]
    #[serde(deserialize_with = "deserialize_optional_flexible_f64")]
    #[serde(default)]
    pub azimuth_time_interval: Option<f64>,

    #[serde(rename = "slantRangeTime")]
    #[serde(deserialize_with = "deserialize_optional_flexible_f64")]
    #[serde(default)]
    pub slant_range_time: Option<f64>,

    #[serde(rename = "prf")]
    #[serde(deserialize_with = "deserialize_optional_flexible_f64")]
    #[serde(default)]
    pub prf: Option<f64>,

    #[serde(rename = "azimuthFrequency")]
    #[serde(deserialize_with = "deserialize_optional_flexible_f64")]
    #[serde(default)]
    pub azimuth_frequency: Option<f64>,

    #[serde(rename = "burstCount")]
    #[serde(deserialize_with = "deserialize_optional_usize")]
    #[serde(default)]
    pub burst_count: Option<usize>,

    #[serde(rename = "rangeSamplesPerBurst")]
    #[serde(deserialize_with = "deserialize_optional_usize")]
    #[serde(default)]
    pub range_samples_per_burst: Option<usize>,

    #[serde(rename = "azimuthSamplesPerBurst")]
    #[serde(deserialize_with = "deserialize_optional_usize")]
    #[serde(default)]
    pub azimuth_samples_per_burst: Option<usize>,

    #[serde(rename = "firstLineOfValidSamples")]
    #[serde(deserialize_with = "deserialize_optional_usize")]
    #[serde(default)]
    pub first_line_of_valid_samples: Option<usize>,

    #[serde(rename = "lastLineOfValidSamples")]
    #[serde(deserialize_with = "deserialize_optional_usize")]
    #[serde(default)]
    pub last_line_of_valid_samples: Option<usize>,

    #[serde(rename = "firstSampleOfValidSamples")]
    #[serde(deserialize_with = "deserialize_optional_usize")]
    #[serde(default)]
    pub first_sample_of_valid_samples: Option<usize>,

    #[serde(rename = "lastSampleOfValidSamples")]
    #[serde(deserialize_with = "deserialize_optional_usize")]
    #[serde(default)]
    pub last_sample_of_valid_samples: Option<usize>,

    #[serde(rename = "burstDuration")]
    #[serde(deserialize_with = "deserialize_optional_flexible_f64")]
    #[serde(default)]
    pub burst_duration: Option<f64>,
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
    #[serde(default)]
    #[serde(deserialize_with = "deserialize_optional_u32")]
    pub lines_per_burst: Option<u32>,

    #[serde(rename = "samplesPerBurst")]
    #[serde(default)]
    #[serde(deserialize_with = "deserialize_optional_u32")]
    pub samples_per_burst: Option<u32>,

    #[serde(rename = "burstList")]
    pub burst_list: Option<BurstList>,
}

#[derive(Debug, Deserialize)]
pub struct BurstList {
    #[serde(rename = "@count")]
    #[serde(default)]
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
    #[serde(default)]
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
    #[serde(default)]
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

#[derive(Debug, Deserialize, Clone)]
pub struct FlexibleXmlElement {
    #[serde(rename = "re")]
    pub re: Option<f64>,
    #[serde(rename = "im")]
    pub im: Option<f64>,
}

// MAIN PARSER FUNCTION
// ============================================================================
/// Parse Sentinel-1 annotation XML using SIMPLE regex approach
/// NO MORE COMPLEX DESERIALIZERS - just extract what we need!
///
/// This function returns Result<ProductRoot, Box<dyn Error>> for backwards compatibility.
/// For new code, prefer parse_annotation() which returns SarResult.
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

/// Canonical annotation parsing function with consistent error handling
///
/// Returns SarResult for consistency with the rest of the codebase.
/// This is the function to use for all new code.
pub fn parse_annotation(xml_content: &str) -> SarResult<ProductRoot> {
    // Call the legacy function and convert error type
    parse_annotation_xml(xml_content).map_err(|e| SarError::XmlParsing(e.to_string()))
}

// ============================================================================
// UTILITY FUNCTIONS FOR COMPATIBILITY WITH EXISTING CODE
// ============================================================================
