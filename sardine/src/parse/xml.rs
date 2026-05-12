//! Minimal serde structs for Sentinel-1 annotation XML.
//!
//! Only the fields required to populate [`SceneMetadata`] are included.
//! Unknown elements are silently ignored by quick-xml.

use serde::{Deserialize, Deserializer};

// ── Custom deserializers ─────────────────────────────────────────────

/// Deserialize space-separated i32 values (burst valid sample arrays).
pub(crate) fn de_space_i32<'de, D>(deserializer: D) -> Result<Vec<i32>, D::Error>
where
    D: Deserializer<'de>,
{
    let s = String::deserialize(deserializer)?;
    s.split_whitespace()
        .map(|t| t.parse::<i32>().map_err(serde::de::Error::custom))
        .collect()
}

/// Deserialize f64 from XML text content, trimming whitespace.
///
/// quick-xml may pass text with leading/trailing whitespace; `str::parse::<f64>`
/// would reject that, so we trim first.
pub(crate) fn de_f64<'de, D>(deserializer: D) -> Result<f64, D::Error>
where
    D: Deserializer<'de>,
{
    let s = String::deserialize(deserializer)?;
    s.trim()
        .parse::<f64>()
        .map_err(serde::de::Error::custom)
}

/// Deserialize u32 from XML text content, trimming whitespace.
pub(crate) fn de_u32<'de, D>(deserializer: D) -> Result<u32, D::Error>
where
    D: Deserializer<'de>,
{
    let s = String::deserialize(deserializer)?;
    s.trim()
        .parse::<u32>()
        .map_err(serde::de::Error::custom)
}

/// Deserialize space-separated f64 values (e.g. polynomial coefficients).
///
/// Used for `<dataDcPolynomial count="3">a b c</dataDcPolynomial>` where the
/// `count` attribute is ignored and the text content is parsed as floats.
pub(crate) fn de_space_f64<'de, D>(deserializer: D) -> Result<Vec<f64>, D::Error>
where
    D: Deserializer<'de>,
{
    let s = String::deserialize(deserializer)?;
    s.split_whitespace()
        .map(|t| t.parse::<f64>().map_err(serde::de::Error::custom))
        .collect()
}

// ── Root ─────────────────────────────────────────────────────────────

#[derive(Debug, Deserialize)]
pub(crate) struct ProductXml {
    #[serde(rename = "adsHeader")]
    pub ads_header: AdsHeaderXml,

    #[serde(rename = "generalAnnotation")]
    pub general_annotation: GeneralAnnotationXml,

    #[serde(rename = "imageAnnotation")]
    pub image_annotation: ImageAnnotationXml,

    /// Top-level `<dopplerCentroid>` block.  Present in all IW/EW SLC products
    /// but absent in SM/WV and GRD, so `Option` + `default` for safe fallback.
    #[serde(rename = "dopplerCentroid", default)]
    pub doppler_centroid: Option<DopplerCentroidXml>,

    #[serde(rename = "swathTiming")]
    pub swath_timing: SwathTimingXml,

    #[serde(rename = "geolocationGrid")]
    pub geolocation_grid: GeolocationGridXml,
}

// ── ADS Header ───────────────────────────────────────────────────────

#[derive(Debug, Deserialize)]
pub(crate) struct AdsHeaderXml {
    #[serde(rename = "missionId")]
    pub mission_id: String,

    pub polarisation: String,
    pub mode: String,
    pub swath: String,

    #[serde(rename = "startTime")]
    pub start_time: String,

    #[serde(rename = "stopTime")]
    pub stop_time: String,
}

// ── General Annotation ───────────────────────────────────────────────

#[derive(Debug, Deserialize)]
pub(crate) struct GeneralAnnotationXml {
    #[serde(rename = "productInformation")]
    pub product_information: ProductInformationXml,

    #[serde(rename = "downlinkInformationList")]
    pub downlink_information_list: DownlinkInfoListXml,

    #[serde(rename = "orbitList")]
    pub orbit_list: OrbitListXml,

    /// Azimuth FM rate polynomial list, nested inside `<generalAnnotation>`.
    /// Absent in some product types, so `Option` + `default`.
    #[serde(rename = "azimuthFmRateList", default)]
    pub azimuth_fm_rate_list: Option<AzimuthFmRateListXml>,
}

#[derive(Debug, Deserialize)]
pub(crate) struct ProductInformationXml {
    #[serde(rename = "radarFrequency")]
    #[serde(deserialize_with = "de_f64")]
    pub radar_frequency: f64,

    #[serde(rename = "rangeSamplingRate")]
    #[serde(deserialize_with = "de_f64")]
    pub range_sampling_rate: f64,
}

#[derive(Debug, Deserialize)]
pub(crate) struct DownlinkInfoListXml {
    #[serde(rename = "downlinkInformation", default)]
    pub items: Vec<DownlinkInfoXml>,
}

#[derive(Debug, Deserialize)]
pub(crate) struct DownlinkInfoXml {
    #[serde(deserialize_with = "de_f64")]
    pub prf: f64,
}

// ── Orbit ────────────────────────────────────────────────────────────

#[derive(Debug, Deserialize)]
pub(crate) struct OrbitListXml {
    #[serde(rename = "orbit", default)]
    pub orbits: Vec<OrbitXml>,
}

#[derive(Debug, Deserialize)]
pub(crate) struct OrbitXml {
    pub time: String,
    pub position: Vec3Xml,
    pub velocity: Vec3Xml,
}

#[derive(Debug, Deserialize)]
pub(crate) struct Vec3Xml {
    #[serde(deserialize_with = "de_f64")]
    pub x: f64,
    #[serde(deserialize_with = "de_f64")]
    pub y: f64,
    #[serde(deserialize_with = "de_f64")]
    pub z: f64,
}

// ── Image Annotation ─────────────────────────────────────────────────

#[derive(Debug, Deserialize)]
pub(crate) struct ImageAnnotationXml {
    #[serde(rename = "imageInformation")]
    pub image_information: ImageInformationXml,
}

#[derive(Debug, Deserialize)]
#[allow(dead_code)]
pub(crate) struct ImageInformationXml {
    #[serde(rename = "productFirstLineUtcTime")]
    pub product_first_line_utc_time: String,

    #[serde(rename = "productLastLineUtcTime")]
    pub product_last_line_utc_time: String,

    #[serde(rename = "ascendingNodeTime")]
    pub ascending_node_time: String,

    #[serde(rename = "slantRangeTime")]
    #[serde(deserialize_with = "de_f64")]
    pub slant_range_time: f64,

    #[serde(rename = "rangePixelSpacing")]
    #[serde(deserialize_with = "de_f64")]
    pub range_pixel_spacing: f64,

    #[serde(rename = "azimuthPixelSpacing")]
    #[serde(deserialize_with = "de_f64")]
    pub azimuth_pixel_spacing: f64,

    #[serde(rename = "azimuthTimeInterval")]
    #[serde(deserialize_with = "de_f64")]
    pub azimuth_time_interval: f64,

    #[serde(rename = "azimuthFrequency")]
    #[serde(deserialize_with = "de_f64")]
    pub azimuth_frequency: f64,

    #[serde(rename = "numberOfSamples")]
    #[serde(deserialize_with = "de_u32")]
    pub number_of_samples: u32,

    #[serde(rename = "numberOfLines")]
    #[serde(deserialize_with = "de_u32")]
    pub number_of_lines: u32,
}

// ── Swath Timing ─────────────────────────────────────────────────────

#[derive(Debug, Deserialize)]
#[allow(dead_code)]
pub(crate) struct SwathTimingXml {
    #[serde(rename = "linesPerBurst")]
    #[serde(deserialize_with = "de_u32")]
    pub lines_per_burst: u32,

    #[serde(rename = "samplesPerBurst")]
    #[serde(deserialize_with = "de_u32")]
    pub samples_per_burst: u32,

    #[serde(rename = "burstList")]
    pub burst_list: BurstListXml,
}

#[derive(Debug, Deserialize)]
pub(crate) struct BurstListXml {
    #[serde(rename = "burst", default)]
    pub bursts: Vec<BurstXml>,
}

#[derive(Debug, Deserialize)]
#[allow(dead_code)]
pub(crate) struct BurstXml {
    #[serde(rename = "azimuthTime")]
    pub azimuth_time: String,

    #[serde(rename = "azimuthAnxTime")]
    #[serde(deserialize_with = "de_f64")]
    pub azimuth_anx_time: f64,

    #[serde(rename = "firstValidSample")]
    #[serde(deserialize_with = "de_space_i32")]
    pub first_valid_sample: Vec<i32>,

    #[serde(rename = "lastValidSample")]
    #[serde(deserialize_with = "de_space_i32")]
    pub last_valid_sample: Vec<i32>,
}

// ── Geolocation Grid ─────────────────────────────────────────────────

#[derive(Debug, Deserialize)]
pub(crate) struct GeolocationGridXml {
    #[serde(rename = "geolocationGridPointList")]
    pub point_list: GeoGridPointListXml,
}

#[derive(Debug, Deserialize)]
pub(crate) struct GeoGridPointListXml {
    #[serde(rename = "geolocationGridPoint", default)]
    pub points: Vec<GeoGridPointXml>,
}

#[derive(Debug, Deserialize)]
pub(crate) struct GeoGridPointXml {
    #[serde(rename = "azimuthTime")]
    pub azimuth_time: String,
    #[serde(rename = "slantRangeTime", deserialize_with = "de_f64")]
    pub slant_range_time: f64,
    #[serde(deserialize_with = "de_u32")]
    pub line: u32,
    #[serde(deserialize_with = "de_u32")]
    pub pixel: u32,
    #[serde(deserialize_with = "de_f64")]
    pub latitude: f64,
    #[serde(deserialize_with = "de_f64")]
    pub longitude: f64,
    #[serde(deserialize_with = "de_f64")]
    pub height: f64,
    #[serde(rename = "incidenceAngle", deserialize_with = "de_f64")]
    pub incidence_angle: f64,
    #[serde(rename = "elevationAngle", deserialize_with = "de_f64")]
    pub elevation_angle: f64,
}

// ── Doppler Centroid ─────────────────────────────────────────────────

/// Top-level `<dopplerCentroid>` block (direct child of `<product>`).
#[derive(Debug, Deserialize)]
pub(crate) struct DopplerCentroidXml {
    #[serde(rename = "dcEstimateList")]
    pub dc_estimate_list: DcEstimateListXml,
}

#[derive(Debug, Deserialize)]
pub(crate) struct DcEstimateListXml {
    #[serde(rename = "dcEstimate", default)]
    pub estimates: Vec<DcEstimateXml>,
}

/// A single `<dcEstimate>` entry.
///
/// The `<dataDcPolynomial count="3">` element contains three space-separated
/// f64 coefficients [a0, a1, a2].  The `count` attribute is ignored by the
/// deserializer; `de_space_f64` reads only the text content.
#[derive(Debug, Deserialize)]
pub(crate) struct DcEstimateXml {
    #[serde(rename = "azimuthTime")]
    pub azimuth_time: String,

    /// Reference slant-range time (seconds) — polynomial zero-point in range.
    #[serde(deserialize_with = "de_f64")]
    pub t0: f64,

    /// Data-derived Doppler centroid polynomial coefficients [a0, a1, a2].
    #[serde(rename = "dataDcPolynomial", deserialize_with = "de_space_f64")]
    pub data_dc_polynomial: Vec<f64>,
}

// ── Azimuth FM Rate ───────────────────────────────────────────────────

/// The `<azimuthFmRateList>` block nested inside `<generalAnnotation>`.
#[derive(Debug, Deserialize)]
pub(crate) struct AzimuthFmRateListXml {
    #[serde(rename = "azimuthFmRate", default)]
    pub rates: Vec<AzimuthFmRateXml>,
}

/// A single `<azimuthFmRate>` entry.
#[derive(Debug, Deserialize)]
pub(crate) struct AzimuthFmRateXml {
    #[serde(rename = "azimuthTime")]
    pub azimuth_time: String,

    /// Reference slant-range time (seconds).
    #[serde(deserialize_with = "de_f64")]
    pub t0: f64,

    /// FM rate polynomial coefficients [a0, a1, a2].
    #[serde(rename = "azimuthFmRatePolynomial", deserialize_with = "de_space_f64")]
    pub polynomial: Vec<f64>,
}
