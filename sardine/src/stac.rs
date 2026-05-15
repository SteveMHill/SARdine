//! STAC Item JSON generation for sardine processed outputs.
//!
//! Generates a [STAC 1.0.0](https://stacspec.org/) Feature item alongside the
//! main GeoTIFF output.  The item embeds:
//!
//! - SAR extension (`sar:*`) fields from the Sentinel-1 annotation.
//! - Projection extension (`proj:epsg`) when a CRS is known.
//! - Backscatter, LIA, and mask assets (whichever were written).
//! - WGS84 `bbox` and GeoJSON `geometry` polygon derived from the output
//!   raster extent (TC mode) or the input SLC footprint (GRD mode).
//!
//! # Sidecar path convention
//!
//! The STAC item is written to `<output>.stac.json` (same stem, new suffix).
//! Example: `scene_VV.tif` → `scene_VV.stac.json`.
//!
//! # Orbit source serialisation
//!
//! `sardine:orbit_source` is `"poeorb"` or `"annotation"` — same strings as
//! the provenance JSON so they can be joined without remapping.

use std::path::Path;

use serde::Serialize;
use thiserror::Error;

use crate::provenance::Provenance;

// ─── Error type ───────────────────────────────────────────────────────────────

#[derive(Debug, Error)]
pub enum StacError {
    #[error("serialising STAC item to JSON: {0}")]
    Serialize(#[from] serde_json::Error),

    #[error("writing STAC item file `{path}`: {source}")]
    Io {
        path: String,
        #[source]
        source: std::io::Error,
    },

    #[error("projecting output raster corners to WGS84: {0}")]
    Projection(String),
}

// ─── STAC structs ─────────────────────────────────────────────────────────────

#[derive(Debug, Serialize)]
struct StacItem {
    #[serde(rename = "type")]
    item_type: &'static str,
    stac_version: &'static str,
    stac_extensions: Vec<&'static str>,
    id: String,
    geometry: GeoJsonPolygon,
    /// [west, south, east, north] in WGS84 degrees.
    bbox: [f64; 4],
    properties: StacProperties,
    assets: serde_json::Map<String, serde_json::Value>,
    links: Vec<()>,
}

#[derive(Debug, Serialize)]
struct GeoJsonPolygon {
    #[serde(rename = "type")]
    geometry_type: &'static str,
    /// Outer ring only; first and last vertex are identical (GeoJSON spec).
    coordinates: Vec<Vec<[f64; 2]>>,
}

/// STAC Item `properties` block.
///
/// All `sar:*` fields follow the
/// [SAR extension v1.0.0 spec](https://github.com/stac-extensions/sar).
/// All `proj:*` fields follow the
/// [Projection extension v1.1.0 spec](https://github.com/stac-extensions/projection).
#[derive(Debug, Serialize)]
struct StacProperties {
    /// Primary timestamp (RFC 3339): scene start time.
    datetime: String,
    start_datetime: String,
    end_datetime: String,

    /// Platform identifier, e.g. `"sentinel-1b"`.
    platform: String,
    /// SAR instrument.
    instruments: Vec<&'static str>,

    #[serde(rename = "sar:instrument_mode")]
    sar_instrument_mode: String,
    #[serde(rename = "sar:frequency_band")]
    sar_frequency_band: &'static str,
    /// Radar centre frequency in GHz, derived from the annotation `radar_frequency_hz`.
    #[serde(rename = "sar:center_frequency")]
    sar_center_frequency: f64,
    /// SAR look direction; `"right"` for all Sentinel-1 acquisitions.
    #[serde(rename = "sar:observation_direction")]
    sar_observation_direction: &'static str,
    /// One or both of `"VV"`, `"VH"`.
    #[serde(rename = "sar:polarizations")]
    sar_polarizations: Vec<String>,
    /// `"RTC"` (terrain-corrected) or `"GRD"` (ground range).
    #[serde(rename = "sar:product_type")]
    sar_product_type: String,
    /// Range multilook factor.
    #[serde(rename = "sar:looks_range")]
    sar_looks_range: usize,
    /// Azimuth multilook factor.
    #[serde(rename = "sar:looks_azimuth")]
    sar_looks_azimuth: usize,
    /// Output pixel spacing in range direction, metres; omitted for TC (use proj:transform).
    #[serde(rename = "sar:pixel_spacing_range", skip_serializing_if = "Option::is_none")]
    sar_pixel_spacing_range: Option<f64>,
    /// Output pixel spacing in azimuth direction, metres; omitted for TC.
    #[serde(rename = "sar:pixel_spacing_azimuth", skip_serializing_if = "Option::is_none")]
    sar_pixel_spacing_azimuth: Option<f64>,

    /// Output CRS EPSG code (projection extension).
    #[serde(rename = "proj:epsg", skip_serializing_if = "Option::is_none")]
    proj_epsg: Option<u32>,

    /// Orbit pass direction per SAT extension v1.0.0: `"ascending"` or `"descending"`.
    /// Omitted for InSAR items (no provenance struct available in that path).
    #[serde(rename = "sat:orbit_state", skip_serializing_if = "Option::is_none")]
    sat_orbit_state: Option<String>,
    /// Absolute orbit counter per SAT extension v1.0.0.
    /// Omitted for InSAR items.
    #[serde(rename = "sat:absolute_orbit", skip_serializing_if = "Option::is_none")]
    sat_absolute_orbit: Option<u32>,

    /// sardine crate version (e.g. `"0.1.0"`).
    #[serde(rename = "sardine:version")]
    sardine_version: String,
    /// Orbit source: `"poeorb"` or `"annotation"`.
    #[serde(rename = "sardine:orbit_source")]
    sardine_orbit_source: String,

    // ── CARD4L NRB / processing extension fields ──────────────────────────

    /// Processing software identifier (CARD4L NRB §3.5, processing extension).
    /// Value: `"sardine v{version}"`.
    #[serde(rename = "processing:software")]
    processing_software: serde_json::Value,
    /// Processing level (CARD4L NRB §3.5, processing extension). Always `"L2"`.
    #[serde(rename = "processing:level")]
    processing_level: &'static str,
    /// DEM source identifier, e.g. `"SRTM-1"` or `"GLO-30"` (CARD4L NRB §3.6).
    #[serde(rename = "sardine:dem_source")]
    sardine_dem_source: String,
    /// Geoid model identifier, e.g. `"EGM96"` or `"none"` (CARD4L NRB §3.6).
    #[serde(rename = "sardine:geoid")]
    sardine_geoid: String,
    /// Whether terrain flattening (σ⁰ → γ⁰, Small 2011) was applied (CARD4L NRB §3.3).
    /// `null` for GRD mode (no terrain correction).
    #[serde(rename = "sardine:terrain_flattening", skip_serializing_if = "Option::is_none")]
    sardine_terrain_flattening: Option<bool>,
    /// Noise floor threshold in dB used for masking sub-noise-floor pixels.
    /// `0.0` means no noise-floor masking was applied (CARD4L NRB §3.4).
    #[serde(rename = "sardine:noise_floor_db")]
    sardine_noise_floor_db: f32,
}

// ─── Public API ───────────────────────────────────────────────────────────────

/// Generate and write a STAC 1.0.0 Feature item from a completed provenance
/// record.
///
/// The `path` argument is the destination file (typically `<output>.stac.json`
/// obtained via [`crate::run::sidecar_path`]).
///
/// # Errors
///
/// Returns [`StacError::Projection`] when a UTM output raster's corners cannot
/// be projected back to WGS84 (this should never happen for valid S-1 scenes
/// over land).  Returns [`StacError::Io`] on write failure.
pub fn write_stac_item(prov: &Provenance, path: &Path) -> Result<(), StacError> {
    let bbox = compute_wgs84_bbox(prov)?;
    let [west, south, east, north] = bbox;

    // Rectangular GeoJSON polygon from bbox corners.
    let ring = vec![
        [west, south],
        [east, south],
        [east, north],
        [west, north],
        [west, south], // close ring
    ];

    let item_id = format!("{}_{}", prov.input.product_id, prov.input.polarization);

    // Platform string: lower-case, e.g. "sentinel-1b".
    let platform = prov.input.mission.to_lowercase();

    // Product type: follow SAR STAC extension conventions.
    // mode strings are set by build_provenance / build_grd_provenance in run_provenance.rs.
    let sar_product_type = match prov.processing.mode.as_deref() {
        Some("rtc") => "RTC",
        Some("nrb") => "NRB",
        _ => "GRD", // covers "grd" and any unknown mode
    }
    .to_owned();

    let properties = StacProperties {
        datetime: prov.input.scene_start_utc.clone(),
        start_datetime: prov.input.scene_start_utc.clone(),
        end_datetime: prov.input.scene_stop_utc.clone(),
        platform,
        instruments: vec!["c-sar"],
        sar_instrument_mode: prov.input.acquisition_mode.to_uppercase(),
        sar_frequency_band: "C",
        sar_center_frequency: prov.input.radar_frequency_hz / 1e9,
        sar_observation_direction: "right",
        sar_polarizations: vec![prov.input.polarization.clone()],
        sar_product_type,
        sar_looks_range: prov.output.range_looks.unwrap_or(1),
        sar_looks_azimuth: prov.output.azimuth_looks.unwrap_or(1),
        sar_pixel_spacing_range: prov.output.range_pixel_spacing_m,
        sar_pixel_spacing_azimuth: prov.output.azimuth_pixel_spacing_m,
        proj_epsg: prov.output.crs_epsg,
        sat_orbit_state: Some(prov.input.orbit_pass_direction.clone()),
        sat_absolute_orbit: Some(prov.input.absolute_orbit_number),
        sardine_version: prov.sardine.version.clone(),
        sardine_orbit_source: match prov.orbit.source {
            crate::provenance::OrbitSource::Poeorb => "poeorb".to_owned(),
            crate::provenance::OrbitSource::Annotation => "annotation".to_owned(),
        },
        processing_software: serde_json::json!({ "sardine": prov.sardine.version.clone() }),
        processing_level: "L2",
        sardine_dem_source: infer_dem_source(&prov.dem.directory),
        sardine_geoid: infer_geoid_name(&prov.geoid.spec),
        sardine_terrain_flattening: prov.processing.flatten,
        sardine_noise_floor_db: prov.processing.noise_floor_db,
    };

    let mut assets: serde_json::Map<String, serde_json::Value> = serde_json::Map::new();
    assets.insert(
        "backscatter".to_owned(),
        stac_asset(
            &prov.output.raster_path,
            "image/tiff; application=geotiff",
            match prov.output.units.as_str() {
                "dB" => "Backscatter (dB)",
                _ => "Backscatter (linear σ⁰)",
            },
            &["data"],
        ),
    );
    if let Some(p) = &prov.output.lia_sidecar_path {
        assets.insert(
            "lia".to_owned(),
            stac_asset(p, "image/tiff; application=geotiff", "Local incidence angle cosine", &["data"]),
        );
    }
    if let Some(p) = &prov.output.mask_sidecar_path {
        assets.insert(
            "mask".to_owned(),
            stac_asset(p, "image/tiff; application=geotiff", "Quality mask", &["data"]),
        );
    }

    let item = StacItem {
        item_type: "Feature",
        stac_version: "1.0.0",
        stac_extensions: vec![
            "https://stac-extensions.github.io/sar/v1.0.0/schema.json",
            "https://stac-extensions.github.io/projection/v1.1.0/schema.json",
            "https://stac-extensions.github.io/processing/v1.2.0/schema.json",
            "https://stac-extensions.github.io/sat/v1.0.0/schema.json",
        ],
        id: item_id,
        geometry: GeoJsonPolygon {
            geometry_type: "Polygon",
            coordinates: vec![ring],
        },
        bbox,
        properties,
        assets,
        links: vec![],
    };

    let bytes = serde_json::to_vec_pretty(&item)?;
    std::fs::write(path, &bytes).map_err(|source| StacError::Io {
        path: path.display().to_string(),
        source,
    })?;
    Ok(())
}

// ─── InSAR STAC writer ────────────────────────────────────────────────────────

/// Inputs to [`write_insar_stac_item`].  Constructed by `run_insar` from the
/// already-computed per-subswath data — no `Provenance` struct required.
pub struct InsarStacInput<'a> {
    /// Reference scene product ID (e.g. `"S1B_IW_SLC__1SDV_20190123T053348…"`).
    pub ref_product_id: &'a str,
    /// Secondary scene product ID.
    pub sec_product_id: &'a str,
    /// Platform string, lower-case (e.g. `"sentinel-1b"`).
    pub platform: &'a str,
    /// Polarization channel used (e.g. `"VV"`).
    pub polarization: &'a str,
    /// Reference scene start UTC as RFC 3339 string.
    pub ref_start_utc: &'a str,
    /// Reference scene stop UTC as RFC 3339 string.
    pub ref_stop_utc: &'a str,
    /// Sub-swath label (e.g. `"iw1"`).
    pub subswath: &'a str,
    /// Coherence window azimuth looks.
    pub az_looks: usize,
    /// Coherence window range looks.
    pub rg_looks: usize,
    /// Output CRS EPSG code.
    pub crs_epsg: u32,
    /// GDAL geotransform of the geocoded coherence raster.
    pub geotransform: [f64; 6],
    /// Width of the geocoded coherence raster in pixels.
    pub cols: usize,
    /// Height of the geocoded coherence raster in pixels.
    pub rows: usize,
    /// `true` when a POEORB file was supplied for the reference scene.
    pub ref_orbit_is_poeorb: bool,
    /// Acquisition mode of the reference scene (e.g. `"iw"`, `"ew"`).
    pub acquisition_mode: &'a str,
    /// DEM source identifier (e.g. `"srtm1"`).
    pub dem_source: &'a str,
    /// Geoid model spec (e.g. `"auto"`, `"zero"`, or a file path).
    pub geoid_spec: &'a str,
    /// Orbit pass direction of the reference scene (`"ascending"` or `"descending"`).
    pub orbit_pass_direction: &'a str,
    /// Absolute orbit number of the reference scene.
    pub absolute_orbit_number: u32,
    /// Filesystem path of the coherence GeoTIFF (for the `coherence` asset).
    pub coherence_path: &'a str,
    /// Filesystem path of the phase GeoTIFF, if it was written.
    pub phase_path: Option<&'a str>,
}

/// Generate and write a STAC 1.0.0 Feature item for one InSAR subswath output.
///
/// The item ID is `<ref_product_id>_<sec_product_id>_<subswath>_<polarization>`.
/// Assets include `coherence` (always) and `phase` (when `input.phase_path` is
/// `Some`).  The bounding box is derived from the geocoded coherence geotransform.
///
/// The `path` argument is the destination file (e.g.
/// `<out_base>_iw1_coherence.stac.json`).
pub fn write_insar_stac_item(input: &InsarStacInput<'_>, path: &Path) -> Result<(), StacError> {
    // Compute WGS84 bbox from the geocoded coherence geotransform.
    let bbox = compute_wgs84_bbox_from_gt(
        input.geotransform,
        input.cols,
        input.rows,
        input.crs_epsg,
    )?;
    let [west, south, east, north] = bbox;

    let ring = vec![
        [west, south],
        [east, south],
        [east, north],
        [west, north],
        [west, south],
    ];

    let item_id = format!(
        "{}_{}_{}_{}",
        input.ref_product_id,
        input.sec_product_id,
        input.subswath,
        input.polarization,
    );

    let properties = StacProperties {
        datetime: input.ref_start_utc.to_owned(),
        start_datetime: input.ref_start_utc.to_owned(),
        end_datetime: input.ref_stop_utc.to_owned(),
        platform: input.platform.to_owned(),
        instruments: vec!["c-sar"],
        sar_instrument_mode: input.acquisition_mode.to_uppercase(),
        sar_frequency_band: "C",
        // S-1 C-band centre frequency; InSAR path does not carry radar_frequency_hz.
        sar_center_frequency: 5.405,
        sar_observation_direction: "right", // Sentinel-1 is always right-looking
        sar_polarizations: vec![input.polarization.to_owned()],
        sar_product_type: "InSAR".to_owned(),
        sar_looks_range: input.rg_looks,
        sar_looks_azimuth: input.az_looks,
        sar_pixel_spacing_range: None,
        sar_pixel_spacing_azimuth: None,
        proj_epsg: Some(input.crs_epsg),
        sat_orbit_state: Some(input.orbit_pass_direction.to_owned()),
        sat_absolute_orbit: Some(input.absolute_orbit_number),
        sardine_version: env!("CARGO_PKG_VERSION").to_owned(),
        sardine_orbit_source: if input.ref_orbit_is_poeorb {
            "poeorb"
        } else {
            "annotation"
        }
        .to_owned(),
        processing_software: serde_json::json!({ "sardine": env!("CARGO_PKG_VERSION") }),
        processing_level: "L2",
        sardine_dem_source: input.dem_source.to_owned(),
        sardine_geoid: infer_geoid_name(input.geoid_spec),
        sardine_terrain_flattening: None,
        sardine_noise_floor_db: 0.0,
    };

    let mut assets: serde_json::Map<String, serde_json::Value> = serde_json::Map::new();
    assets.insert(
        "coherence".to_owned(),
        stac_asset(
            input.coherence_path,
            "image/tiff; application=geotiff",
            "InSAR coherence magnitude",
            &["data"],
        ),
    );
    if let Some(phase_path) = input.phase_path {
        assets.insert(
            "phase".to_owned(),
            stac_asset(
                phase_path,
                "image/tiff; application=geotiff",
                "Wrapped interferometric phase (radians)",
                &["data"],
            ),
        );
    }

    let item = StacItem {
        item_type: "Feature",
        stac_version: "1.0.0",
        stac_extensions: vec![
            "https://stac-extensions.github.io/sar/v1.0.0/schema.json",
            "https://stac-extensions.github.io/projection/v1.1.0/schema.json",
        ],
        id: item_id,
        geometry: GeoJsonPolygon {
            geometry_type: "Polygon",
            coordinates: vec![ring],
        },
        bbox,
        properties,
        assets,
        links: vec![],
    };

    let bytes = serde_json::to_vec_pretty(&item)?;
    std::fs::write(path, &bytes).map_err(|source| StacError::Io {
        path: path.display().to_string(),
        source,
    })?;
    Ok(())
}

// ─── PolSAR STAC writer ──────────────────────────────────────────────────────

/// Inputs to [`write_polsar_stac_item`].  Constructed by `run_polsar` from the
/// geocoded H-band output — all three bands share the same CRS and dimensions.
pub struct PolsarStacInput<'a> {
    /// Scene product ID (e.g. `"S1B_IW_SLC__1SDV_20190123T053348…"`).
    pub product_id: &'a str,
    /// Platform string, lower-case (e.g. `"sentinel-1b"`).
    pub platform: &'a str,
    /// Co-pol channel label, upper-case (e.g. `"VV"`).
    pub pol_co: &'a str,
    /// Cross-pol channel label, upper-case (e.g. `"VH"`).
    pub pol_cross: &'a str,
    /// Scene start UTC as RFC 3339 string.
    pub scene_start_utc: &'a str,
    /// Scene stop UTC as RFC 3339 string.
    pub scene_stop_utc: &'a str,
    /// Azimuth multilook factor applied to the C2 covariance matrix.
    pub az_looks: usize,
    /// Range multilook factor applied to the C2 covariance matrix.
    pub rg_looks: usize,
    /// Output CRS EPSG code (same for all three bands).
    pub crs_epsg: u32,
    /// GDAL geotransform of the geocoded output (same for all three bands).
    pub geotransform: [f64; 6],
    /// Raster width in pixels.
    pub cols: usize,
    /// Raster height in pixels.
    pub rows: usize,
    /// `true` when a POEORB file was supplied.
    pub orbit_is_poeorb: bool,
    /// Acquisition mode, lower-case (e.g. `"iw"`).
    pub acquisition_mode: &'a str,
    /// DEM source identifier (e.g. `"srtm1"`).
    pub dem_source: &'a str,
    /// Geoid model spec (e.g. `"auto"`, `"zero"`, or a file path).
    pub geoid_spec: &'a str,
    /// Orbit pass direction (`"ascending"` or `"descending"`).
    pub orbit_pass_direction: &'a str,
    /// Absolute orbit number.
    pub absolute_orbit_number: u32,
    /// Filesystem path of the entropy (H) GeoTIFF.
    pub h_path: &'a str,
    /// Filesystem path of the anisotropy (A) GeoTIFF.
    pub a_path: &'a str,
    /// Filesystem path of the mean-alpha (α) GeoTIFF.
    pub alpha_path: &'a str,
}

/// Generate and write a STAC 1.0.0 Feature item for a PolSAR H/A/Alpha output.
///
/// The item ID is `<product_id>_polsar_<pol_co>+<pol_cross>`.
/// Assets include `entropy` (H), `anisotropy` (A), and `mean_alpha` (α).
/// The bounding box is derived from the geocoded geotransform.
///
/// The `path` argument is the destination file (e.g. `<output>.polsar.stac.json`).
pub fn write_polsar_stac_item(
    input: &PolsarStacInput<'_>,
    path: &Path,
) -> Result<(), StacError> {
    let bbox = compute_wgs84_bbox_from_gt(
        input.geotransform,
        input.cols,
        input.rows,
        input.crs_epsg,
    )?;
    let [west, south, east, north] = bbox;

    let ring = vec![
        [west, south],
        [east, south],
        [east, north],
        [west, north],
        [west, south],
    ];

    let item_id = format!(
        "{}_polsar_{}+{}",
        input.product_id,
        input.pol_co,
        input.pol_cross,
    );

    let properties = StacProperties {
        datetime: input.scene_start_utc.to_owned(),
        start_datetime: input.scene_start_utc.to_owned(),
        end_datetime: input.scene_stop_utc.to_owned(),
        platform: input.platform.to_owned(),
        instruments: vec!["c-sar"],
        sar_instrument_mode: input.acquisition_mode.to_uppercase(),
        sar_frequency_band: "C",
        sar_center_frequency: 5.405,
        sar_observation_direction: "right",
        sar_polarizations: vec![
            input.pol_co.to_owned(),
            input.pol_cross.to_owned(),
        ],
        sar_product_type: "PolSAR".to_owned(),
        sar_looks_range: input.rg_looks,
        sar_looks_azimuth: input.az_looks,
        sar_pixel_spacing_range: None,
        sar_pixel_spacing_azimuth: None,
        proj_epsg: Some(input.crs_epsg),
        sat_orbit_state: Some(input.orbit_pass_direction.to_owned()),
        sat_absolute_orbit: Some(input.absolute_orbit_number),
        sardine_version: env!("CARGO_PKG_VERSION").to_owned(),
        sardine_orbit_source: if input.orbit_is_poeorb {
            "poeorb"
        } else {
            "annotation"
        }
        .to_owned(),
        processing_software: serde_json::json!({ "sardine": env!("CARGO_PKG_VERSION") }),
        processing_level: "L2",
        sardine_dem_source: infer_dem_source(input.dem_source),
        sardine_geoid: infer_geoid_name(input.geoid_spec),
        sardine_terrain_flattening: Some(false),
        sardine_noise_floor_db: 0.0,
    };

    let mut assets: serde_json::Map<String, serde_json::Value> = serde_json::Map::new();
    assets.insert(
        "entropy".to_owned(),
        stac_asset(
            input.h_path,
            "image/tiff; application=geotiff",
            "Polarimetric entropy H ∈ [0, 1]",
            &["data"],
        ),
    );
    assets.insert(
        "anisotropy".to_owned(),
        stac_asset(
            input.a_path,
            "image/tiff; application=geotiff",
            "Polarimetric anisotropy A ∈ [0, 1]",
            &["data"],
        ),
    );
    assets.insert(
        "mean_alpha".to_owned(),
        stac_asset(
            input.alpha_path,
            "image/tiff; application=geotiff",
            "Mean scattering angle alpha ∈ [0°, 90°]",
            &["data"],
        ),
    );

    let item = StacItem {
        item_type: "Feature",
        stac_version: "1.0.0",
        stac_extensions: vec![
            "https://stac-extensions.github.io/sar/v1.0.0/schema.json",
            "https://stac-extensions.github.io/projection/v1.1.0/schema.json",
            "https://stac-extensions.github.io/processing/v1.2.0/schema.json",
            "https://stac-extensions.github.io/sat/v1.0.0/schema.json",
        ],
        id: item_id,
        geometry: GeoJsonPolygon {
            geometry_type: "Polygon",
            coordinates: vec![ring],
        },
        bbox,
        properties,
        assets,
        links: vec![],
    };

    let bytes = serde_json::to_vec_pretty(&item)?;
    std::fs::write(path, &bytes).map_err(|source| StacError::Io {
        path: path.display().to_string(),
        source,
    })?;
    Ok(())
}

// ─── Helpers ──────────────────────────────────────────────────────────────────

/// Infer a human-readable DEM source identifier from the DEM directory path.
///
/// Checks the directory path for well-known substrings; falls back to the
/// directory's base name when none match.
fn infer_dem_source(dem_dir: &str) -> String {
    let lower = dem_dir.to_lowercase();
    if lower.contains("srtm1") || lower.contains("srtm_1") {
        return "SRTM-1".to_owned();
    }
    if lower.contains("srtm3") || lower.contains("srtm_3") {
        return "SRTM-3".to_owned();
    }
    if lower.contains("glo30") || lower.contains("cop30") || lower.contains("copernicus30") {
        return "GLO-30".to_owned();
    }
    if lower.contains("glo90") || lower.contains("cop90") {
        return "GLO-90".to_owned();
    }
    // Fall back to the last path component.
    std::path::Path::new(dem_dir)
        .file_name()
        .and_then(|n| n.to_str())
        .unwrap_or(dem_dir)
        .to_owned()
}

/// Infer a human-readable geoid model name from the `--geoid` spec string.
///
/// `"auto"` and paths containing `"egm96"` map to `"EGM96"`.
/// `"zero"` maps to `"none"` (no geoid correction).
fn infer_geoid_name(spec: &str) -> String {
    let lower = spec.to_lowercase();
    if lower == "auto" || lower.contains("egm96") {
        return "EGM96".to_owned();
    }
    if lower.contains("egm2008") {
        return "EGM2008".to_owned();
    }
    if lower == "zero" {
        return "none".to_owned();
    }
    // Unknown spec — return verbatim so nothing is silently dropped.
    spec.to_owned()
}

/// Build a STAC asset JSON object.
fn stac_asset(href: &str, media_type: &str, title: &str, roles: &[&str]) -> serde_json::Value {
    // Use the filename only (relative href) so the STAC item is portable:
    // an absolute path embeds the processing machine's filesystem layout.
    let relative_href = Path::new(href)
        .file_name()
        .and_then(|s| s.to_str())
        .unwrap_or(href);
    serde_json::json!({
        "href": relative_href,
        "type": media_type,
        "title": title,
        "roles": roles,
    })
}

/// Derive a WGS84 [west, south, east, north] bounding box from a provenance
/// record.
///
/// For TC outputs with a known geotransform and EPSG code, the four raster
/// corners are projected back to WGS84 via [`crate::output_crs`].  For all
/// other cases (GRD outputs, or TC outputs without geotransform) the input
/// scene footprint stored in the provenance is used instead.
fn compute_wgs84_bbox(prov: &Provenance) -> Result<[f64; 4], StacError> {
    if let (Some(gt), Some(epsg)) = (prov.output.geotransform, prov.output.crs_epsg) {
        return compute_wgs84_bbox_from_gt(gt, prov.output.cols, prov.output.rows, epsg);
    }

    // GRD mode or TC without geotransform: fall back to scene footprint.
    let bb = &prov.input.scene_bbox_deg;
    Ok([bb.min_lon, bb.min_lat, bb.max_lon, bb.max_lat])
}

/// Compute WGS84 [west, south, east, north] from a GDAL geotransform + raster
/// dimensions + EPSG code.  Used by both [`compute_wgs84_bbox`] (provenance
/// path) and [`write_insar_stac_item`] (InSAR direct path).
fn compute_wgs84_bbox_from_gt(
    gt: [f64; 6],
    cols: usize,
    rows: usize,
    epsg: u32,
) -> Result<[f64; 4], StacError> {
    let x0 = gt[0];
    let dx = gt[1];
    let y0 = gt[3];
    let dy = gt[5];

    let x1 = x0 + cols as f64 * dx;
    let y1 = y0 + rows as f64 * dy;

    if epsg == 4326 {
        let west = x0.min(x1);
        let east = x0.max(x1);
        let south = y0.min(y1);
        let north = y0.max(y1);
        return Ok([west, south, east, north]);
    }

    // UTM or other projected CRS: project all four corners to WGS84.
    let crs = crate::output_crs::OutputCrs::from_epsg(epsg)
        .map_err(|e| StacError::Projection(e.to_string()))?;
    let projector = crs
        .projector()
        .map_err(|e| StacError::Projection(e.to_string()))?;

    let corners = [(x0, y0), (x1, y0), (x0, y1), (x1, y1)];
    let mut lons = [0f64; 4];
    let mut lats = [0f64; 4];
    for (i, &(x, y)) in corners.iter().enumerate() {
        let (lon, lat) = projector
            .xy_to_lonlat(x, y)
            .map_err(|e| StacError::Projection(e.to_string()))?;
        lons[i] = lon;
        lats[i] = lat;
    }

    let west = lons.iter().cloned().fold(f64::INFINITY, f64::min);
    let east = lons.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    let south = lats.iter().cloned().fold(f64::INFINITY, f64::min);
    let north = lats.iter().cloned().fold(f64::NEG_INFINITY, f64::max);

    // An anti-meridian crossing yields west > east in a simple min/max bbox.
    // Warn loudly; the caller may need to split or clamp the bbox for
    // downstream STAC consumers that do not support anti-meridian geometries.
    if west > east {
        tracing::warn!(
            "STAC bbox anti-meridian crossing detected (west={:.4}°, east={:.4}°); \
             bbox [west, south, east, north] = [{:.4}, {:.4}, {:.4}, {:.4}] may be \
             invalid for consumers that do not handle anti-meridian splits",
            west, east, west, south, east, north
        );
    }

    Ok([west, south, east, north])
}

// ─── Tests ────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use crate::provenance::{
        BoundingBoxJson, DemInfo, GeoidInfo, InputInfo, OrbitInfo, OrbitSource, OutputInfo,
        ProcessingInfo, Provenance, SardineInfo, StatsInfo, WarningsInfo, SCHEMA_VERSION,
    };

    fn sample_prov_wgs84() -> Provenance {
        Provenance {
            schema_version: SCHEMA_VERSION,
            generated_utc: "2026-04-29T12:00:00Z".to_owned(),
            sardine: SardineInfo {
                package: "sardine".to_owned(),
                version: "0.1.0".to_owned(),
            },
            input: InputInfo {
                safe_path: "/data/X.SAFE".to_owned(),
                product_id: "S1B_IW_SLC__1SDV_20190123T053348_20190123T053415_014617_01B3D4_E833".to_owned(),
                mission: "S1B".to_owned(),
                acquisition_mode: "IW".to_owned(),
                orbit_pass_direction: "descending".to_owned(),
                absolute_orbit_number: 14617,
                polarization: "VV".to_owned(),
                radar_frequency_hz: 5_405_000_450.0,
                scene_start_utc: "2019-01-23T05:33:48Z".to_owned(),
                scene_stop_utc: "2019-01-23T05:34:15Z".to_owned(),
                scene_bbox_deg: BoundingBoxJson {
                    min_lat: 47.0, max_lat: 48.5,
                    min_lon: 7.0,  max_lon: 9.5,
                },
            },
            orbit: OrbitInfo {
                source: OrbitSource::Poeorb,
                file_path: None,
                reference_epoch_utc: "2019-01-22T22:59:42Z".to_owned(),
                state_vector_count: 17,
            },
            dem: DemInfo { directory: "/dem".to_owned(), tile_count: 4 },
            geoid: GeoidInfo { spec: "zero".to_owned() },
            processing: ProcessingInfo {
                mode: Some("rtc".to_owned()),
                pixel_spacing_deg: Some(0.0001),
                target_spacing_m: None,
                flatten: Some(true),
                compute_lia: Some(false),
                noise_floor_db: 0.0,
                threads: 0,
                speckle_filter: None,
                speckle_window: None,
                speckle_enl: None,
                speckle_damping: None,
                speckle_order: None,
                resampling_kernel: None,
            },
            output: OutputInfo {
                raster_path: "/out/scene_VV.tif".to_owned(),
                lia_sidecar_path: None,
                mask_sidecar_path: None,
                cols: 1000,
                rows: 800,
                // WGS84: origin = top-left (7.0 lon, 48.5 lat), dx=0.001, dy=-0.001
                geotransform: Some([7.0, 0.001, 0.0, 48.5, 0.0, -0.001]),
                crs_epsg: Some(4326),
                range_pixel_spacing_m: None,
                azimuth_pixel_spacing_m: None,
                range_looks: Some(4),
                azimuth_looks: None,
                gcps_count: None,
                units: "dB".to_owned(),
                nodata: "NaN".to_owned(),
            },
            stats: StatsInfo {
                total_pixel_count: 800_000,
                valid_pixel_count: 600_000,
                dem_missing_count: Some(0),
                non_converged_count: Some(0),
                flat_masked_count: Some(0),
                noise_masked_count: Some(0),
                db_converted_count: Some(600_000),
                db_min: Some(-25.0),
                db_max: Some(-4.0),
                db_mean: Some(-14.5),
                db_median: Some(-15.0),
            },
            warnings: WarningsInfo {
                cal_lut_extrapolation_gap_px: 0,
                noise_lut_extrapolation_gap_px: 0,
            },
            quality_flags: vec![],
        }
    }

    #[test]
    fn wgs84_bbox_from_geotransform() {
        let prov = sample_prov_wgs84();
        // origin (7.0, 48.5), dx=0.001, dy=-0.001, 1000 cols, 800 rows
        // → x1 = 7.0 + 1000*0.001 = 8.0, y1 = 48.5 + 800*(-0.001) = 47.7
        let bbox = compute_wgs84_bbox(&prov).unwrap();
        assert!((bbox[0] - 7.0).abs() < 1e-9, "west={}", bbox[0]);
        assert!((bbox[1] - 47.7).abs() < 1e-6, "south={}", bbox[1]);
        assert!((bbox[2] - 8.0).abs() < 1e-9, "east={}", bbox[2]);
        assert!((bbox[3] - 48.5).abs() < 1e-9, "north={}", bbox[3]);
    }

    #[test]
    fn utm_bbox_projects_to_wgs84() {
        let mut prov = sample_prov_wgs84();
        // UTM32N (EPSG:32632): Munich area ≈ (690000 E, 5335000 N)
        // Use a small 100×100 m tile so corners are close.
        prov.output.geotransform = Some([690000.0, 10.0, 0.0, 5336000.0, 0.0, -10.0]);
        prov.output.crs_epsg = Some(32632);
        prov.output.cols = 10;
        prov.output.rows = 10;

        let bbox = compute_wgs84_bbox(&prov).unwrap();
        // Munich is around lon 11.5°, lat 48.1°; a 100×100 m tile is a tiny sliver.
        assert!(bbox[0] > 11.0 && bbox[0] < 12.0, "west={}", bbox[0]);
        assert!(bbox[2] > 11.0 && bbox[2] < 12.0, "east={}", bbox[2]);
        assert!(bbox[1] > 48.0 && bbox[1] < 49.0, "south={}", bbox[1]);
        assert!(bbox[3] > 48.0 && bbox[3] < 49.0, "north={}", bbox[3]);
        assert!(bbox[2] > bbox[0], "east > west");
        assert!(bbox[3] > bbox[1], "north > south");
    }

    #[test]
    fn grd_mode_falls_back_to_scene_footprint() {
        let mut prov = sample_prov_wgs84();
        prov.output.geotransform = None;
        prov.output.crs_epsg = None;
        prov.input.scene_bbox_deg = BoundingBoxJson {
            min_lat: 47.1, max_lat: 48.4,
            min_lon: 7.2,  max_lon: 9.3,
        };

        let bbox = compute_wgs84_bbox(&prov).unwrap();
        assert!((bbox[0] - 7.2).abs() < 1e-9);
        assert!((bbox[1] - 47.1).abs() < 1e-9);
        assert!((bbox[2] - 9.3).abs() < 1e-9);
        assert!((bbox[3] - 48.4).abs() < 1e-9);
    }

    #[test]
    fn write_stac_item_creates_valid_json() {
        let prov = sample_prov_wgs84();
        let dir = tempfile::tempdir().expect("tempdir");
        let path = dir.path().join("scene_VV.stac.json");

        write_stac_item(&prov, &path).expect("write");

        let bytes = std::fs::read(&path).expect("read back");
        let v: serde_json::Value = serde_json::from_slice(&bytes).expect("parse");

        assert_eq!(v["type"], "Feature");
        assert_eq!(v["stac_version"], "1.0.0");
        assert!(v["id"].as_str().unwrap().contains("S1B"));
        assert_eq!(v["geometry"]["type"], "Polygon");

        // bbox sanity
        let bbox = v["bbox"].as_array().unwrap();
        assert_eq!(bbox.len(), 4);

        // properties
        assert_eq!(v["properties"]["platform"], "s1b");
        assert_eq!(v["properties"]["sar:instrument_mode"], "IW");
        assert_eq!(v["properties"]["sar:product_type"], "RTC");
        assert_eq!(v["properties"]["proj:epsg"], 4326);
        assert_eq!(v["properties"]["sardine:orbit_source"], "poeorb");
        assert_eq!(v["properties"]["sar:looks_range"], 4);

        // assets
        assert!(v["assets"]["backscatter"]["href"].is_string());
        assert!(v["assets"].get("lia").is_none(), "no lia sidecar expected");
        assert!(v["assets"].get("mask").is_none(), "no mask sidecar expected");
    }

    #[test]
    fn write_stac_item_includes_lia_and_mask_assets() {
        let mut prov = sample_prov_wgs84();
        prov.output.lia_sidecar_path = Some("/out/scene_VV.lia.tif".to_owned());
        prov.output.mask_sidecar_path = Some("/out/scene_VV.mask.tif".to_owned());

        let dir = tempfile::tempdir().expect("tempdir");
        let path = dir.path().join("scene_VV.stac.json");
        write_stac_item(&prov, &path).expect("write");

        let bytes = std::fs::read(&path).expect("read back");
        let v: serde_json::Value = serde_json::from_slice(&bytes).expect("parse");
        assert!(v["assets"]["lia"]["href"].is_string());
        assert!(v["assets"]["mask"]["href"].is_string());
    }

    #[test]
    fn stac_geometry_ring_is_closed() {
        let prov = sample_prov_wgs84();
        let dir = tempfile::tempdir().expect("tempdir");
        let path = dir.path().join("x.stac.json");
        write_stac_item(&prov, &path).expect("write");

        let bytes = std::fs::read(&path).expect("read back");
        let v: serde_json::Value = serde_json::from_slice(&bytes).expect("parse");
        let ring = v["geometry"]["coordinates"][0].as_array().unwrap();
        assert_eq!(ring.len(), 5, "ring must have 5 points (closed)");
        assert_eq!(ring[0], ring[4], "first and last vertex must match (closed ring)");
    }

    // ── InSAR STAC item tests ─────────────────────────────────────────────────

    fn sample_insar_input() -> InsarStacInput<'static> {
        // Minimal WGS84 geotransform: origin=(10.0,55.0), dx=0.0001, dy=-0.0001
        InsarStacInput {
            ref_product_id: "S1B_IW_SLC__1SDV_20190123T053348_RFF",
            sec_product_id: "S1B_IW_SLC__1SDV_20190204T053348_SEC",
            platform: "s1b",
            polarization: "VV",
            ref_start_utc: "2019-01-23T05:33:48Z",
            ref_stop_utc: "2019-01-23T05:34:15Z",
            subswath: "iw1",
            az_looks: 1,
            rg_looks: 4,
            crs_epsg: 4326,
            geotransform: [10.0, 0.0001, 0.0, 55.0, 0.0, -0.0001],
            cols: 100,
            rows: 50,
            ref_orbit_is_poeorb: true,
            acquisition_mode: "iw",
            dem_source: "srtm1",
            geoid_spec: "auto",
            orbit_pass_direction: "descending",
            absolute_orbit_number: 14617,
            coherence_path: "out_iw1_coh.tif",
            phase_path: None,
        }
    }

    #[test]
    fn insar_stac_coh_only_valid_json() {
        let dir = tempfile::tempdir().expect("tempdir");
        let path = dir.path().join("insar.stac.json");
        let input = sample_insar_input();
        write_insar_stac_item(&input, &path).expect("write");

        let bytes = std::fs::read(&path).expect("read back");
        let v: serde_json::Value = serde_json::from_slice(&bytes).expect("parse");

        assert_eq!(v["type"], "Feature");
        assert_eq!(v["stac_version"], "1.0.0");
        // coherence asset present
        assert!(v["assets"]["coherence"]["href"].is_string());
        // phase asset absent
        assert!(v["assets"]["phase"].is_null());
        // product type
        assert_eq!(v["properties"]["sar:product_type"], "InSAR");
        // projection EPSG
        assert_eq!(v["properties"]["proj:epsg"], 4326);
        // orbit source
        assert_eq!(v["properties"]["sardine:orbit_source"], "poeorb");
        // looks
        assert_eq!(v["properties"]["sar:looks_azimuth"], 1);
        assert_eq!(v["properties"]["sar:looks_range"], 4);
    }

    #[test]
    fn insar_stac_with_phase_asset() {
        let dir = tempfile::tempdir().expect("tempdir");
        let path = dir.path().join("insar_phase.stac.json");
        let mut input = sample_insar_input();
        input.phase_path = Some("out_iw1_phase.tif");
        write_insar_stac_item(&input, &path).expect("write");

        let bytes = std::fs::read(&path).expect("read back");
        let v: serde_json::Value = serde_json::from_slice(&bytes).expect("parse");

        assert!(v["assets"]["coherence"]["href"].is_string());
        assert_eq!(v["assets"]["phase"]["href"], "out_iw1_phase.tif");
    }

    #[test]
    fn insar_stac_annotation_orbit_source() {
        let dir = tempfile::tempdir().expect("tempdir");
        let path = dir.path().join("insar_ann.stac.json");
        let mut input = sample_insar_input();
        input.ref_orbit_is_poeorb = false;
        write_insar_stac_item(&input, &path).expect("write");

        let bytes = std::fs::read(&path).expect("read back");
        let v: serde_json::Value = serde_json::from_slice(&bytes).expect("parse");
        assert_eq!(v["properties"]["sardine:orbit_source"], "annotation");
    }

    #[test]
    fn insar_stac_geometry_ring_is_closed() {
        let dir = tempfile::tempdir().expect("tempdir");
        let path = dir.path().join("insar_ring.stac.json");
        let input = sample_insar_input();
        write_insar_stac_item(&input, &path).expect("write");

        let bytes = std::fs::read(&path).expect("read back");
        let v: serde_json::Value = serde_json::from_slice(&bytes).expect("parse");
        let ring = v["geometry"]["coordinates"][0].as_array().unwrap();
        assert_eq!(ring.len(), 5, "ring must have 5 points (closed)");
        assert_eq!(ring[0], ring[4], "first and last vertex must match");
    }
}
