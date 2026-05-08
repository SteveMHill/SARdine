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
    sar_instrument_mode: &'static str,
    #[serde(rename = "sar:frequency_band")]
    sar_frequency_band: &'static str,
    /// One or both of `"VV"`, `"VH"`.
    #[serde(rename = "sar:polarizations")]
    sar_polarizations: Vec<String>,
    /// `"RTC"` (terrain-corrected) or `"GRD"` (ground range).
    #[serde(rename = "sar:product_type")]
    sar_product_type: String,
    /// Range multilook factor; omitted when equal to 1.
    #[serde(rename = "sar:looks_range", skip_serializing_if = "Option::is_none")]
    sar_looks_range: Option<usize>,
    /// Azimuth multilook factor; omitted when equal to 1.
    #[serde(rename = "sar:looks_azimuth", skip_serializing_if = "Option::is_none")]
    sar_looks_azimuth: Option<usize>,

    /// Output CRS EPSG code (projection extension).
    #[serde(rename = "proj:epsg", skip_serializing_if = "Option::is_none")]
    proj_epsg: Option<u32>,

    /// sardine crate version (e.g. `"0.1.0"`).
    #[serde(rename = "sardine:version")]
    sardine_version: String,
    /// Orbit source: `"poeorb"` or `"annotation"`.
    #[serde(rename = "sardine:orbit_source")]
    sardine_orbit_source: String,
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

    // Product type: γ⁰ RTC, σ⁰ RTC, or GRD.
    let sar_product_type = match prov.processing.mode.as_deref() {
        Some("tc") => match prov.processing.flatten {
            Some(true) => "RTC-Gamma0",
            _ => "RTC",
        },
        _ => "GRD",
    }
    .to_owned();

    let properties = StacProperties {
        datetime: prov.input.scene_start_utc.clone(),
        start_datetime: prov.input.scene_start_utc.clone(),
        end_datetime: prov.input.scene_stop_utc.clone(),
        platform,
        instruments: vec!["c-sar"],
        sar_instrument_mode: "IW",
        sar_frequency_band: "C",
        sar_polarizations: vec![prov.input.polarization.clone()],
        sar_product_type,
        sar_looks_range: prov.output.range_looks,
        sar_looks_azimuth: prov.output.azimuth_looks,
        proj_epsg: prov.output.crs_epsg,
        sardine_version: prov.sardine.version.clone(),
        sardine_orbit_source: match prov.orbit.source {
            crate::provenance::OrbitSource::Poeorb => "poeorb".to_owned(),
            crate::provenance::OrbitSource::Annotation => "annotation".to_owned(),
        },
    };

    let mut assets: serde_json::Map<String, serde_json::Value> = serde_json::Map::new();
    assets.insert(
        "backscatter".to_owned(),
        stac_asset(
            &prov.output.raster_path,
            "image/tiff; application=geotiff",
            if prov.processing.mode.as_deref() == Some("grd") {
                "Backscatter (linear σ⁰)"
            } else {
                "Backscatter (dB)"
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

/// Build a STAC asset JSON object.
fn stac_asset(href: &str, media_type: &str, title: &str, roles: &[&str]) -> serde_json::Value {
    serde_json::json!({
        "href": href,
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
        let cols = prov.output.cols as f64;
        let rows = prov.output.rows as f64;

        // gt = [origin_x, dx, 0, origin_y, 0, dy] — standard GDAL convention.
        // dy is negative for north-up rasters.
        let x0 = gt[0];
        let dx = gt[1];
        let y0 = gt[3];
        let dy = gt[5];

        let x1 = x0 + cols * dx;
        let y1 = y0 + rows * dy;

        if epsg == 4326 {
            // Geotransform is already in lon/lat degrees — no projection needed.
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
        return Ok([west, south, east, north]);
    }

    // GRD mode or TC without geotransform: fall back to scene footprint.
    let bb = &prov.input.scene_bbox_deg;
    Ok([bb.min_lon, bb.min_lat, bb.max_lon, bb.max_lat])
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
                package: "sardine-scene".to_owned(),
                version: "0.1.0".to_owned(),
            },
            input: InputInfo {
                safe_path: "/data/X.SAFE".to_owned(),
                product_id: "S1B_IW_SLC__1SDV_20190123T053348_20190123T053415_014617_01B3D4_E833".to_owned(),
                mission: "S1B".to_owned(),
                acquisition_mode: "IW".to_owned(),
                polarization: "VV".to_owned(),
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
                mode: Some("tc".to_owned()),
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
        assert_eq!(v["properties"]["sar:product_type"], "RTC-Gamma0");
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
}
