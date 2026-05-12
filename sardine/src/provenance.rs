//! Provenance sidecar for processed scenes.
//!
//! Writes a small JSON document next to the main GeoTIFF output that records
//! the inputs, parameters, and outcome counts of a single `sardine process`
//! invocation.  The schema is versioned (`schema_version`) so downstream
//! tooling can reject incompatible documents.
//!
//! # Faithful, not exhaustive
//!
//! Only fields whose values are unambiguously available at the point of
//! writing are recorded.  In particular:
//!
//! * **DEM tile paths and checksums** are *not* recorded — the
//!   [`crate::dem::DemMosaic`] currently does not expose tile paths through
//!   its public API, and synthesising them from the input directory would
//!   be a fabrication (some tiles in the directory may not be loaded, e.g.
//!   wrong extension).  The DEM directory and tile *count* are recorded
//!   instead.  See `dem.directory` and `dem.tile_count`.
//! * **Git SHA** is *not* recorded.  Adding it requires a `build.rs` script
//!   that shells out to `git`, which would couple the build to a working
//!   `git` binary and a non-shallow clone.  The crate version
//!   (`CARGO_PKG_VERSION`) is recorded instead.
//! * **Geoid model identity** is recorded as the user-supplied `--geoid`
//!   spec string, not as a hash of the loaded grid.  Hashing would require
//!   exposing the grid's backing buffer; the spec is sufficient to
//!   reproduce the load.
//!
//! Every omission is intentional and follows the project's "small explicit
//! failure over large plausible-but-wrong implementation" rule.

use std::path::Path;

use serde::Serialize;
use thiserror::Error;

/// Current provenance schema version.
///
/// Increment when adding a required field, removing a field, or changing
/// the meaning of an existing field.  Adding a new optional field does
/// not require a version bump.
pub const SCHEMA_VERSION: u32 = 1;

/// Errors writing a provenance sidecar.
#[derive(Debug, Error)]
pub enum ProvenanceError {
    #[error("serialising provenance to JSON: {0}")]
    Serialize(#[from] serde_json::Error),

    #[error("writing provenance file `{path}`: {source}")]
    Io {
        path: String,
        #[source]
        source: std::io::Error,
    },
}

/// Type of quality flag raised during processing.
///
/// Each variant corresponds to a specific condition that downstream consumers
/// may wish to filter on.  New variants may be added without a schema version
/// bump because consumers should ignore unknown variants.
#[derive(Debug, Clone, Copy, Serialize, PartialEq, Eq)]
#[serde(rename_all = "SCREAMING_SNAKE_CASE")]
pub enum QualityFlagType {
    /// Annotation orbit was used; geometric accuracy is metre-level.
    /// Supply a POEORB `.EOF` file via `--orbit` for cm-level accuracy.
    OrbitAnnotation,
    /// Calibration LUT did not fully cover the range extent; the last valid
    /// sample was extrapolated beyond the swath edge.
    CalibrationLutExtrapolation,
    /// Noise LUT did not fully cover the range extent; extrapolated at the
    /// swath edge.
    NoiseLutExtrapolation,
    /// More than 1% of output pixels had no DEM elevation value.
    DemPartialCoverage,
    /// One or more output pixels' zero-Doppler Newton iteration did not
    /// converge within the configured iteration limit.
    SolverNonConvergence,
}

/// Severity of a quality flag.
#[derive(Debug, Clone, Copy, Serialize, PartialEq, Eq)]
#[serde(rename_all = "SCREAMING_SNAKE_CASE")]
pub enum QualityFlagSeverity {
    Info,
    Warning,
    Error,
}

/// A single quality flag entry embedded in the provenance document.
#[derive(Debug, Clone, Serialize)]
pub struct QualityFlagEntry {
    pub flag_type: QualityFlagType,
    pub severity: QualityFlagSeverity,
    /// Human-readable description of the condition, including relevant
    /// numeric context (counts, percentages, pixel gaps).
    pub message: String,
}

/// Source of the orbit state vectors used during processing.
#[derive(Debug, Clone, Copy, Serialize, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
pub enum OrbitSource {
    /// Precise orbit (`POEORB` `.EOF` file) supplied via `--orbit`.
    Poeorb,
    /// Annotation orbit (lower accuracy) used because no `--orbit` was
    /// supplied and `SARDINE_ALLOW_ANNOTATION_ORBIT=1` was set.
    Annotation,
}

/// Top-level provenance document written next to the main output GeoTIFF.
///
/// The struct mirrors the JSON shape directly; nested structs name the
/// JSON object keys.  All field names are stable across schema versions
/// unless the version is incremented.
#[derive(Debug, Clone, Serialize)]
pub struct Provenance {
    /// Schema version of this document.  See [`SCHEMA_VERSION`].
    pub schema_version: u32,
    /// UTC timestamp when this document was written, RFC 3339.
    pub generated_utc: String,

    pub sardine: SardineInfo,
    pub input: InputInfo,
    pub orbit: OrbitInfo,
    pub dem: DemInfo,
    pub geoid: GeoidInfo,
    pub processing: ProcessingInfo,
    pub output: OutputInfo,
    pub stats: StatsInfo,
    pub warnings: WarningsInfo,
    /// Quality flags raised during processing.  An empty array means no
    /// conditions were detected.  The array is always present in the JSON
    /// (never omitted); consumers should treat an absent key as equivalent
    /// to an empty array for forward-compatibility with old documents.
    pub quality_flags: Vec<QualityFlagEntry>,
}

#[derive(Debug, Clone, Serialize)]
pub struct SardineInfo {
    /// Cargo package name (e.g. `"sardine"`).
    pub package: String,
    /// Cargo package version string (e.g. `"0.1.0"`).
    pub version: String,
}

#[derive(Debug, Clone, Serialize)]
pub struct InputInfo {
    /// Absolute path to the input `.SAFE` directory.
    pub safe_path: String,
    /// SAFE product identifier (the `.SAFE` directory basename without the
    /// `.SAFE` suffix).
    pub product_id: String,
    /// Mission ID (e.g. `"S1A"`, `"S1B"`).
    pub mission: String,
    /// Acquisition mode (e.g. `"IW"`).
    pub acquisition_mode: String,
    /// Polarization channel processed (e.g. `"VV"`).
    pub polarization: String,
    /// Scene first-line UTC time (RFC 3339).
    pub scene_start_utc: String,
    /// Scene last-line UTC time (RFC 3339).
    pub scene_stop_utc: String,
    /// Geographic bounding box of the scene footprint (degrees, WGS84).
    pub scene_bbox_deg: BoundingBoxJson,
}

#[derive(Debug, Clone, Serialize)]
pub struct BoundingBoxJson {
    pub min_lat: f64,
    pub max_lat: f64,
    pub min_lon: f64,
    pub max_lon: f64,
}

#[derive(Debug, Clone, Serialize)]
pub struct OrbitInfo {
    /// Source of the orbit state vectors.
    pub source: OrbitSource,
    /// Absolute path to the POEORB `.EOF` file when
    /// `source == OrbitSource::Poeorb`; `null` otherwise.
    pub file_path: Option<String>,
    /// Reference epoch of the orbit data (earliest state vector time, RFC 3339).
    pub reference_epoch_utc: String,
    /// Number of state vectors in the orbit data used for processing.
    pub state_vector_count: usize,
}

#[derive(Debug, Clone, Serialize)]
pub struct DemInfo {
    /// Absolute path to the DEM directory passed via `--dem`.
    pub directory: String,
    /// Number of DEM tiles successfully loaded from the directory.
    ///
    /// Tile paths and checksums are intentionally not recorded; see the
    /// module-level documentation.
    pub tile_count: usize,
}

#[derive(Debug, Clone, Serialize)]
pub struct GeoidInfo {
    /// Verbatim `--geoid` spec string supplied on the CLI
    /// (`"auto"`, `"zero"`, or a path).
    pub spec: String,
}

#[derive(Debug, Clone, Serialize)]
pub struct ProcessingInfo {
    /// Pipeline mode: `"tc"` for terrain-corrected (geocoded) output,
    /// `"grd"` for ground-range radar-geometry output.
    ///
    /// Optional only for backward compatibility with v1 documents written
    /// before the GRD subcommand existed; new writers always populate it.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub mode: Option<String>,
    /// Output pixel spacing (degrees) on both axes.  Populated for `mode=tc`;
    /// `null` for `mode=grd` (use [`OutputInfo::range_pixel_spacing_m`] and
    /// [`OutputInfo::azimuth_pixel_spacing_m`] instead).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub pixel_spacing_deg: Option<f64>,
    /// Target ground-range pixel spacing in metres.  Populated for `mode=grd`;
    /// `null` for `mode=tc`.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub target_spacing_m: Option<f64>,
    /// Whether Small (2011) terrain flattening was applied.
    /// `null` for `mode=grd` (flattening requires a DEM and geocoding).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub flatten: Option<bool>,
    /// Whether per-pixel local incidence angle and quality mask were computed.
    /// `null` for `mode=grd` (LIA requires a DEM and geocoding).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub compute_lia: Option<bool>,
    /// Noise floor in dB used for `--noise-floor-db` masking.  A value of
    /// `0.0` means "no noise-floor masking" (only non-positive σ⁰ are
    /// removed during the dB conversion).
    pub noise_floor_db: f32,
    /// Number of Rayon threads requested via `--threads`.  `0` means
    /// "all logical CPUs".
    pub threads: usize,
    /// Speckle filter applied to the linear-power output prior to
    /// dB conversion (`mode=tc`) or directly (`mode=grd`).  Variants
    /// are `"none"`, `"boxcar"`, `"lee"`, `"frost"`, `"gamma_map"`.
    /// Field is optional for backward compatibility with documents
    /// written before the speckle stage existed.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub speckle_filter: Option<String>,
    /// Window size of the speckle filter (odd integer).  `null` when
    /// `speckle_filter == "none"`.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub speckle_window: Option<usize>,
    /// Effective number of looks fed to Lee / Gamma-MAP.  `null` for
    /// boxcar, Frost, or `none`.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub speckle_enl: Option<f32>,
    /// Frost damping factor.  `null` for everything else.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub speckle_damping: Option<f32>,
    /// Where the speckle filter was applied relative to terrain correction.
    /// `"before_tc"` = slant-range geometry; `"after_tc"` = map geometry (default).
    /// `null` when `speckle_filter` is `null` (no speckle applied) or for `mode=grd`.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub speckle_order: Option<String>,
}

#[derive(Debug, Clone, Serialize)]
pub struct OutputInfo {
    /// Absolute path of the main output raster (dB GeoTIFF for `mode=tc`,
    /// linear-scale ground-range TIFF for `mode=grd`).
    pub raster_path: String,
    /// Path of the LIA sidecar (cos LIA, Float32) if `--write-lia` was set.
    pub lia_sidecar_path: Option<String>,
    /// Path of the quality mask sidecar (uint8) if `--write-mask` was set.
    pub mask_sidecar_path: Option<String>,
    /// Output raster width in pixels.
    pub cols: usize,
    /// Output raster height in pixels.
    pub rows: usize,
    /// GDAL geotransform of the output raster, populated for `mode=tc`:
    /// `[origin_lon, dx, 0, origin_lat, 0, -dy]`.  `null` for `mode=grd`
    /// (radar geometry has no map projection).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub geotransform: Option<[f64; 6]>,
    /// EPSG code of the output CRS, populated for `mode=tc`.  `null` for
    /// `mode=grd` (radar geometry, not a map projection).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub crs_epsg: Option<u32>,
    /// Ground-range pixel spacing in metres, populated for `mode=grd`.
    /// `null` for `mode=tc`.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub range_pixel_spacing_m: Option<f64>,
    /// Azimuth pixel spacing in metres, populated for `mode=grd`.
    /// `null` for `mode=tc`.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub azimuth_pixel_spacing_m: Option<f64>,
    /// Range multilook factor applied during ground-range projection,
    /// populated for `mode=grd`.  `null` for `mode=tc`.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub range_looks: Option<usize>,
    /// Azimuth multilook factor applied during ground-range projection,
    /// populated for `mode=grd`.  `null` for `mode=tc`.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub azimuth_looks: Option<usize>,
    /// Units of the main output raster (`"dB"` for `mode=tc`,
    /// `"linear"` for `mode=grd`).
    pub units: String,
    /// Sentinel for invalid pixels in the main raster (always `"NaN"`).
    pub nodata: String,
}

#[derive(Debug, Clone, Serialize)]
pub struct StatsInfo {
    /// Total output pixels (`rows * cols`).
    pub total_pixel_count: usize,
    /// Pixels with a valid σ⁰/γ⁰ value in the output raster.
    pub valid_pixel_count: usize,
    /// Pixels rejected because the DEM had no value at the geocoded location.
    /// `null` for `mode=grd` (no DEM lookup).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub dem_missing_count: Option<usize>,
    /// Pixels rejected because the zero-Doppler iteration did not converge.
    /// `null` for `mode=grd` (no geocoding).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub non_converged_count: Option<usize>,
    /// Pixels masked by the terrain-flattening guard (foreshortening or
    /// shadow/layover).  Always `0` when `processing.flatten == false`;
    /// `null` for `mode=grd`.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub flat_masked_count: Option<usize>,
    /// Pixels masked by the dB-conversion noise floor.  `null` for `mode=grd`
    /// (no dB conversion is performed; the GRD raster is in linear σ⁰).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub noise_masked_count: Option<usize>,
    /// Pixels successfully converted to dB.  `null` for `mode=grd`.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub db_converted_count: Option<usize>,
    /// Minimum dB value among all converted (finite) output pixels.
    /// `null` for `mode=grd` (output is linear) and when no pixels were
    /// converted.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub db_min: Option<f32>,
    /// Maximum dB value among all converted (finite) output pixels.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub db_max: Option<f32>,
    /// Arithmetic mean of all converted (finite) output pixels in dB.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub db_mean: Option<f32>,
    /// Approximate median of converted pixels in dB, estimated by
    /// sorting a stride-sampled subset of at most 100 000 finite values.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub db_median: Option<f32>,
}

#[derive(Debug, Clone, Serialize)]
pub struct WarningsInfo {
    /// Maximum extrapolation gap (range pixels) of the calibration LUT
    /// observed during merging.  `0` means the LUT covered every pixel
    /// without extrapolation.
    pub cal_lut_extrapolation_gap_px: u32,
    /// Maximum extrapolation gap (range pixels) of the noise LUT observed
    /// during merging.
    pub noise_lut_extrapolation_gap_px: u32,
}

impl Provenance {
    /// Construct a [`SardineInfo`] populated from the compile-time Cargo
    /// package metadata of this crate.
    pub fn sardine_info() -> SardineInfo {
        SardineInfo {
            package: env!("CARGO_PKG_NAME").to_owned(),
            version: env!("CARGO_PKG_VERSION").to_owned(),
        }
    }

    /// Serialise this provenance document to a pretty-printed JSON file at
    /// `path`.  The file is overwritten if it already exists.
    pub fn write_json(&self, path: &Path) -> Result<(), ProvenanceError> {
        let bytes = serde_json::to_vec_pretty(self)?;
        std::fs::write(path, bytes).map_err(|source| ProvenanceError::Io {
            path: path.display().to_string(),
            source,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn sample_provenance() -> Provenance {
        Provenance {
            schema_version: SCHEMA_VERSION,
            generated_utc: "2026-04-25T12:34:56Z".to_owned(),
            sardine: SardineInfo {
                package: "sardine".to_owned(),
                version: "0.1.0".to_owned(),
            },
            input: InputInfo {
                safe_path: "/data/SLC/X.SAFE".to_owned(),
                product_id: "X".to_owned(),
                mission: "S1B".to_owned(),
                acquisition_mode: "IW".to_owned(),
                polarization: "VV".to_owned(),
                scene_start_utc: "2019-01-23T05:33:48Z".to_owned(),
                scene_stop_utc: "2019-01-23T05:34:15Z".to_owned(),
                scene_bbox_deg: BoundingBoxJson {
                    min_lat: 47.0,
                    max_lat: 48.5,
                    min_lon: 7.0,
                    max_lon: 9.5,
                },
            },
            orbit: OrbitInfo {
                source: OrbitSource::Poeorb,
                file_path: Some("/data/orbits/X.EOF".to_owned()),
                reference_epoch_utc: "2019-01-22T22:59:42Z".to_owned(),
                state_vector_count: 17,
            },
            dem: DemInfo {
                directory: "/data/dem/srtm1".to_owned(),
                tile_count: 4,
            },
            geoid: GeoidInfo {
                spec: "auto".to_owned(),
            },
            processing: ProcessingInfo {
                mode: Some("tc".to_owned()),
                pixel_spacing_deg: Some(0.0001),
                target_spacing_m: None,
                flatten: Some(true),
                compute_lia: Some(true),
                noise_floor_db: -17.0,
                threads: 0,
                speckle_filter: Some("refined_lee".to_owned()),
                speckle_window: Some(7),
                speckle_enl: Some(1.0),
                speckle_damping: None,
                speckle_order: None,
            },
            output: OutputInfo {
                raster_path: "/out/scene.tif".to_owned(),
                lia_sidecar_path: Some("/out/scene.lia.tif".to_owned()),
                mask_sidecar_path: Some("/out/scene.mask.tif".to_owned()),
                cols: 1234,
                rows: 567,
                geotransform: Some([7.0, 0.0001, 0.0, 48.5, 0.0, -0.0001]),
                crs_epsg: Some(4326),
                range_pixel_spacing_m: None,
                azimuth_pixel_spacing_m: None,
                range_looks: None,
                azimuth_looks: None,
                units: "dB".to_owned(),
                nodata: "NaN".to_owned(),
            },
            stats: StatsInfo {
                total_pixel_count: 1234 * 567,
                valid_pixel_count: 600_000,
                dem_missing_count: Some(100),
                non_converged_count: Some(0),
                flat_masked_count: Some(250),
                noise_masked_count: Some(1234),
                db_converted_count: Some(598_000),
                db_min: Some(-25.3),
                db_max: Some(-4.1),
                db_mean: Some(-14.7),
                db_median: Some(-15.2),
            },
            warnings: WarningsInfo {
                cal_lut_extrapolation_gap_px: 0,
                noise_lut_extrapolation_gap_px: 0,
            },
            quality_flags: vec![],
        }
    }

    /// Pin every top-level JSON key.  If the schema evolves, this test
    /// must be updated *and* `SCHEMA_VERSION` must be bumped.
    #[test]
    fn provenance_json_top_level_keys_are_pinned() {
        let p = sample_provenance();
        let v: serde_json::Value =
            serde_json::from_slice(&serde_json::to_vec(&p).expect("serialize")).expect("parse");
        let obj = v.as_object().expect("top-level object");
        let mut keys: Vec<&str> = obj.keys().map(String::as_str).collect();
        keys.sort();
        assert_eq!(
            keys,
            vec![
                "dem",
                "generated_utc",
                "geoid",
                "input",
                "orbit",
                "output",
                "processing",
                "quality_flags",
                "sardine",
                "schema_version",
                "stats",
                "warnings",
            ]
        );
    }

    #[test]
    fn provenance_schema_version_is_serialized_as_integer() {
        let p = sample_provenance();
        let v: serde_json::Value = serde_json::to_value(&p).expect("serialize");
        let n = v
            .get("schema_version")
            .and_then(serde_json::Value::as_u64)
            .expect("schema_version is a u64");
        assert_eq!(n, SCHEMA_VERSION as u64);
    }

    #[test]
    fn orbit_source_serialises_as_snake_case() {
        let v = serde_json::to_value(OrbitSource::Poeorb).expect("serialize");
        assert_eq!(v, serde_json::Value::String("poeorb".to_owned()));
        let v = serde_json::to_value(OrbitSource::Annotation).expect("serialize");
        assert_eq!(v, serde_json::Value::String("annotation".to_owned()));
    }

    #[test]
    fn write_json_round_trips_through_disk() {
        let dir = tempfile::tempdir().expect("tempdir");
        let path = dir.path().join("p.json");
        let p = sample_provenance();
        p.write_json(&path).expect("write");

        let bytes = std::fs::read(&path).expect("read back");
        let v: serde_json::Value = serde_json::from_slice(&bytes).expect("parse round-trip");

        // Spot-check a deeply nested field to confirm the full tree round-tripped.
        assert_eq!(
            v.get("output")
                .and_then(|o| o.get("geotransform"))
                .and_then(serde_json::Value::as_array)
                .map(|a| a.len()),
            Some(6)
        );
        assert_eq!(
            v.get("orbit")
                .and_then(|o| o.get("file_path"))
                .and_then(serde_json::Value::as_str),
            Some("/data/orbits/X.EOF")
        );
        assert_eq!(
            v.get("processing")
                .and_then(|o| o.get("flatten"))
                .and_then(serde_json::Value::as_bool),
            Some(true)
        );
    }

    #[test]
    fn optional_paths_serialise_as_null_when_absent() {
        let mut p = sample_provenance();
        p.orbit.source = OrbitSource::Annotation;
        p.orbit.file_path = None;
        p.output.lia_sidecar_path = None;
        p.output.mask_sidecar_path = None;

        let v: serde_json::Value = serde_json::to_value(&p).expect("serialize");
        assert!(v
            .get("orbit")
            .and_then(|o| o.get("file_path"))
            .expect("file_path present")
            .is_null());
        assert!(v
            .get("output")
            .and_then(|o| o.get("lia_sidecar_path"))
            .expect("lia_sidecar_path present")
            .is_null());
        assert!(v
            .get("output")
            .and_then(|o| o.get("mask_sidecar_path"))
            .expect("mask_sidecar_path present")
            .is_null());
    }

    #[test]
    fn sardine_info_reads_from_cargo_env() {
        let info = Provenance::sardine_info();
        assert_eq!(info.package, "sardine");
        assert!(!info.version.is_empty());
        // Versions are SemVer-ish; first character should be a digit.
        assert!(
            info.version.chars().next().map(|c| c.is_ascii_digit()) == Some(true),
            "unexpected version string: {}",
            info.version
        );
    }

    /// In `mode=grd` the geocoding-specific fields (`geotransform`,
    /// `crs_epsg`, `pixel_spacing_deg`, `flatten`, `compute_lia`,
    /// `dem_missing_count`, `non_converged_count`, `flat_masked_count`,
    /// `noise_masked_count`, `db_converted_count`) MUST be absent from the
    /// JSON, and the GRD-specific fields (`target_spacing_m`,
    /// `range_pixel_spacing_m`, `azimuth_pixel_spacing_m`, `range_looks`,
    /// `azimuth_looks`) MUST be present.  This is the safety guarantee
    /// that downstream tooling cannot mistake a radar-geometry GRD product
    /// for a geocoded TC product on the basis of provenance alone.
    #[test]
    fn grd_mode_provenance_omits_tc_fields_and_includes_grd_fields() {
        let mut p = sample_provenance();
        p.processing = ProcessingInfo {
            mode: Some("grd".to_owned()),
            pixel_spacing_deg: None,
            target_spacing_m: Some(10.0),
            flatten: None,
            compute_lia: None,
            noise_floor_db: 0.0,
            threads: 0,
            speckle_filter: None,
            speckle_window: None,
            speckle_enl: None,
            speckle_damping: None,
            speckle_order: None,
        };
        p.output = OutputInfo {
            raster_path: "/out/scene.grd.tif".to_owned(),
            lia_sidecar_path: None,
            mask_sidecar_path: None,
            cols: 25_000,
            rows: 1_500,
            geotransform: None,
            crs_epsg: None,
            range_pixel_spacing_m: Some(10.0),
            azimuth_pixel_spacing_m: Some(13.95),
            range_looks: Some(1),
            azimuth_looks: Some(1),
            units: "linear".to_owned(),
            nodata: "NaN".to_owned(),
        };
        p.stats = StatsInfo {
            total_pixel_count: 25_000 * 1_500,
            valid_pixel_count: 25_000 * 1_500,
            dem_missing_count: None,
            non_converged_count: None,
            flat_masked_count: None,
            noise_masked_count: None,
            db_converted_count: None,
            db_min: None,
            db_max: None,
            db_mean: None,
            db_median: None,
        };

        let v: serde_json::Value = serde_json::to_value(&p).expect("serialize");
        let processing = v.get("processing").and_then(|x| x.as_object()).expect("processing");
        let output = v.get("output").and_then(|x| x.as_object()).expect("output");
        let stats = v.get("stats").and_then(|x| x.as_object()).expect("stats");

        // TC fields MUST be absent (skip_serializing_if removes them).
        for key in ["pixel_spacing_deg", "flatten", "compute_lia"] {
            assert!(
                !processing.contains_key(key),
                "processing must NOT contain `{key}` in GRD mode"
            );
        }
        for key in ["geotransform", "crs_epsg"] {
            assert!(
                !output.contains_key(key),
                "output must NOT contain `{key}` in GRD mode"
            );
        }
        for key in [
            "dem_missing_count",
            "non_converged_count",
            "flat_masked_count",
            "noise_masked_count",
            "db_converted_count",
        ] {
            assert!(
                !stats.contains_key(key),
                "stats must NOT contain `{key}` in GRD mode"
            );
        }

        // GRD fields MUST be present.
        assert_eq!(
            processing.get("mode").and_then(serde_json::Value::as_str),
            Some("grd")
        );
        assert_eq!(
            processing
                .get("target_spacing_m")
                .and_then(serde_json::Value::as_f64),
            Some(10.0)
        );
        assert_eq!(
            output
                .get("range_pixel_spacing_m")
                .and_then(serde_json::Value::as_f64),
            Some(10.0)
        );
        assert_eq!(
            output
                .get("azimuth_pixel_spacing_m")
                .and_then(serde_json::Value::as_f64),
            Some(13.95)
        );
        assert_eq!(
            output.get("range_looks").and_then(serde_json::Value::as_u64),
            Some(1)
        );
        assert_eq!(
            output.get("azimuth_looks").and_then(serde_json::Value::as_u64),
            Some(1)
        );
        assert_eq!(
            output.get("units").and_then(serde_json::Value::as_str),
            Some("linear")
        );
    }
}
