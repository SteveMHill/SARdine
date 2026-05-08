//! Library entry points for the `process` and `grd` pipelines.
//!
//! This module centralises the orchestration that used to live inside
//! the `sardine` binary so that:
//!
//! * the CLI binary stays a thin clap-only shim, and
//! * external embedders (Python bindings via PyO3, integration tests,
//!   downstream Rust crates) can run the same pipeline by constructing
//!   plain Rust [`ProcessOptions`] / [`GrdOptions`] structs without
//!   depending on `clap`.
//!
//! ## Domain notes
//!
//! These functions are direct extractions of the CLI subcommands; no
//! algorithmic change is intended.  All field semantics and default values
//! match the CLI 1:1.  Progress and warning messages are emitted via
//! `tracing::info!` / `tracing::warn!`; callers must initialise a
//! [`tracing_subscriber`] to see output (the CLI binary does this automatically).

use std::path::{Path, PathBuf};
use std::time::Instant;

use anyhow::{anyhow, bail, Context, Result};
use rayon::prelude::*;

use crate::output_crs::OutputCrs;
use crate::types::Polarization;

/// Emit a single `[timing]` line on stderr.  Centralised so the format stays
/// consistent across stages and so future migration to a structured
/// timing record (sidecar JSON, tracing span) only has to touch one
/// place.  All times are wall-clock seconds since the supplied `Instant`.
fn report_timing(stage: &str, started: Instant) {
    tracing::info!(
        "[timing] {} = {:.3} s",
        stage,
        started.elapsed().as_secs_f64(),
    );
}

// ─────────────────────────────────────────────────────────────────────────────
// Plain options (no clap)
// ─────────────────────────────────────────────────────────────────────────────

/// Sub-swath and burst selection for IW Split.
///
/// If `subswaths` is empty, all three IW subswaths (IW1, IW2, IW3) are
/// processed.  `burst_range`, when set, is applied as a 0-based inclusive
/// range `[start, end]` to the sorted burst list of each selected subswath.
#[derive(Debug, Clone, Default)]
pub struct IwSelection {
    /// Subswaths to process.  Empty = all three.
    pub subswaths: Vec<crate::types::SubSwathId>,
    /// Burst range within each selected subswath (0-based, inclusive).
    /// `None` = all bursts.
    pub burst_range: Option<(usize, usize)>,
}

/// Output mode for the `sardine process` pipeline.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum OutputMode {
    /// Ground-range radar geometry (no geocoding, no CRS).
    Grd,
    /// Terrain-corrected and terrain-flattened (σ⁰ → γ⁰), geocoded.  Default.
    Rtc,
    /// RTC with mandatory LIA, quality mask, and STAC sidecars (NRB profile).
    Nrb,
}

impl Default for OutputMode {
    fn default() -> Self {
        OutputMode::Rtc
    }
}

/// Where in the pipeline to apply the speckle filter.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SpeckleOrder {
    /// Apply speckle in slant-range geometry, before terrain correction.
    BeforeTc,
    /// Apply speckle in map geometry, after terrain correction.  Default.
    AfterTc,
}

impl Default for SpeckleOrder {
    fn default() -> Self {
        SpeckleOrder::AfterTc
    }
}

/// Plain (non-clap) input set for [`run_process`].
///
/// Every field mirrors a `sardine process` CLI flag of the same name.
/// Defaults match the CLI defaults — see [`ProcessOptions::new`].
#[derive(Debug, Clone)]
pub struct ProcessOptions {
    pub safe: PathBuf,
    pub dem: PathBuf,
    pub output: PathBuf,
    pub orbit: Option<PathBuf>,
    pub polarization: String,
    pub no_flatten: bool,
    pub noise_floor_db: f32,
    pub pixel_spacing_deg: f64,
    pub pixel_spacing_m: f64,
    pub crs: String,
    pub write_mask: bool,
    pub write_lia: bool,
    pub no_provenance: bool,
    /// Convert the output GeoTIFF to Cloud-Optimised GeoTIFF in-place after
    /// writing.  Requires `gdal_translate` (GDAL ≥ 3.1) in PATH.
    pub cog: bool,
    pub threads: usize,
    pub speckle: String,
    pub speckle_window: usize,
    pub enl: f32,
    pub frost_damping: f32,
    /// Geoid spec — REQUIRED, no default.  Accepted: `auto`, `zero`, or
    /// a path to an EGM96 grid file (`.bin` / `.gtx` / `.GRD`).
    pub geoid: String,
    /// Range multilook factor applied *before* terrain correction, on the
    /// merged slant-range σ⁰ image.  1 = no-op (default).
    pub multilook_range: usize,
    /// Azimuth multilook factor applied *before* terrain correction, on the
    /// merged slant-range σ⁰ image.  1 = no-op (default).
    pub multilook_azimuth: usize,
    /// Additional SAFE slice paths for multi-slice (assembled) processing.
    ///
    /// When non-empty, all paths `[safe] ++ extra_safe_paths` are assembled
    /// via [`crate::slice_assembly::assemble_slices`] before processing.
    /// Paths must be ordered ascending in time (earliest slice first).
    pub extra_safe_paths: Vec<PathBuf>,
    /// Sub-swath and burst selection (IW Split).  Default: all three IWs, all bursts.
    pub iw_selection: IwSelection,
    /// Output mode: `Grd`, `Rtc` (default), or `Nrb`.
    pub mode: OutputMode,
    /// Where to apply the speckle filter: `BeforeTc` (slant-range) or `AfterTc` (default).
    pub speckle_order: SpeckleOrder,
}

impl ProcessOptions {
    /// Construct with the same defaults the CLI uses.  `safe`, `dem`,
    /// `output`, and `geoid` have no CLI default and must be supplied.
    pub fn new(safe: PathBuf, dem: PathBuf, output: PathBuf, geoid: String) -> Self {
        Self {
            safe,
            dem,
            output,
            orbit: None,
            polarization: "VV".to_owned(),
            no_flatten: false,
            noise_floor_db: 0.0,
            pixel_spacing_deg: 0.0001,
            pixel_spacing_m: 10.0,
            crs: "wgs84".to_owned(),
            write_mask: false,
            write_lia: false,
            no_provenance: false,
            cog: false,
            threads: 0,
            speckle: "none".to_owned(),
            speckle_window: 7,
            enl: 1.0,
            frost_damping: 1.0,
            geoid,
            multilook_range: 1,
            multilook_azimuth: 1,
            extra_safe_paths: vec![],
            iw_selection: IwSelection::default(),
            mode: OutputMode::default(),
            speckle_order: SpeckleOrder::default(),
        }
    }

    /// When `mode == OutputMode::Grd`, the effective ground-range target spacing.
    /// Uses `pixel_spacing_m` (the metric spacing field shared with TC mode).
    pub fn target_spacing_m_for_grd(&self) -> f64 {
        self.pixel_spacing_m
    }
}

/// Plain (non-clap) input set for [`run_grd`].
///
/// Mirrors the `sardine grd` CLI flags one-to-one.
#[derive(Debug, Clone)]
pub struct GrdOptions {
    pub safe: PathBuf,
    pub output: PathBuf,
    pub orbit: Option<PathBuf>,
    pub polarization: String,
    pub target_spacing_m: f64,
    pub no_provenance: bool,
    /// Convert the output GeoTIFF to Cloud-Optimised GeoTIFF in-place after
    /// writing.  Requires `gdal_translate` (GDAL ≥ 3.1) in PATH.
    pub cog: bool,
    pub threads: usize,
    pub speckle: String,
    pub speckle_window: usize,
    pub enl: f32,
    pub frost_damping: f32,
    /// Sub-swath and burst selection (IW Split).  Default: all three IWs, all bursts.
    pub iw_selection: IwSelection,
}

impl GrdOptions {
    /// Construct with CLI defaults.  `safe` and `output` are required.
    pub fn new(safe: PathBuf, output: PathBuf) -> Self {
        Self {
            safe,
            output,
            orbit: None,
            polarization: "VV".to_owned(),
            target_spacing_m: 10.0,
            no_provenance: false,
            cog: false,
            threads: 0,
            speckle: "none".to_owned(),
            speckle_window: 7,
            enl: 1.0,
            frost_damping: 1.0,
            iw_selection: IwSelection::default(),
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Helpers (extracted verbatim from the CLI)
// ─────────────────────────────────────────────────────────────────────────────

/// Compute summary statistics from a dB-converted raster (NaN = masked).
///
/// Returns `(min, max, mean, median)`.  The median is approximated by
/// stride-sampling at most 100 000 finite values, sorting that sample, and
/// taking the middle element.  For the large rasters produced by this
/// pipeline (~400 M valid pixels) this gives an estimate accurate to within
/// ~0.1 dB.  Returns `(NaN, NaN, NaN, NaN)` when no finite pixels exist.
fn compute_db_stats(data: &[f32]) -> (f32, f32, f32, f32) {
    let mut min = f32::INFINITY;
    let mut max = f32::NEG_INFINITY;
    let mut sum = 0.0_f64;
    let mut n = 0usize;
    for &v in data {
        if v.is_finite() {
            if v < min {
                min = v;
            }
            if v > max {
                max = v;
            }
            sum += v as f64;
            n += 1;
        }
    }
    if n == 0 {
        return (f32::NAN, f32::NAN, f32::NAN, f32::NAN);
    }
    let mean = (sum / n as f64) as f32;

    const MAX_SAMPLE: usize = 100_000;
    let stride = (n / MAX_SAMPLE).max(1);
    let mut sample: Vec<f32> = data
        .iter()
        .filter(|v| v.is_finite())
        .step_by(stride)
        .copied()
        .collect();
    sample.sort_unstable_by(|a, b| a.partial_cmp(b).expect("filtered to finite")); // SAFETY-OK: all values in `sample` are finite (guaranteed by the filter above)
    let median = sample[sample.len() / 2];

    (min, max, mean, median)
}

/// Build the quality flags list for one processing run.
///
/// Pass `dem_missing_count = None` and `non_converged_count = None` for GRD
/// mode, which has no geocoding.
fn collect_quality_flags(
    orbit_source: crate::provenance::OrbitSource,
    cal_lut_extrapolation_gap_px: u32,
    noise_lut_extrapolation_gap_px: u32,
    total_pixel_count: usize,
    dem_missing_count: Option<usize>,
    non_converged_count: Option<usize>,
) -> Vec<crate::provenance::QualityFlagEntry> {
    use crate::provenance::{
        OrbitSource, QualityFlagEntry, QualityFlagSeverity, QualityFlagType,
    };
    let mut flags: Vec<QualityFlagEntry> = Vec::new();

    if orbit_source == OrbitSource::Annotation {
        flags.push(QualityFlagEntry {
            flag_type: QualityFlagType::OrbitAnnotation,
            severity: QualityFlagSeverity::Warning,
            message: "Annotation orbit used; geometric accuracy is metre-level. \
                     Supply a POEORB .EOF file via --orbit for cm-level accuracy."
                .to_owned(),
        });
    }
    if cal_lut_extrapolation_gap_px > 0 {
        flags.push(QualityFlagEntry {
            flag_type: QualityFlagType::CalibrationLutExtrapolation,
            severity: QualityFlagSeverity::Warning,
            message: format!(
                "Calibration LUT extrapolated by up to {} range pixels at the swath edge.",
                cal_lut_extrapolation_gap_px
            ),
        });
    }
    if noise_lut_extrapolation_gap_px > 0 {
        flags.push(QualityFlagEntry {
            flag_type: QualityFlagType::NoiseLutExtrapolation,
            severity: QualityFlagSeverity::Warning,
            message: format!(
                "Noise LUT extrapolated by up to {} range pixels at the swath edge.",
                noise_lut_extrapolation_gap_px
            ),
        });
    }
    if let Some(dem_missing) = dem_missing_count {
        let pct = dem_missing as f64 / total_pixel_count.max(1) as f64 * 100.0;
        if pct > 1.0 {
            flags.push(QualityFlagEntry {
                flag_type: QualityFlagType::DemPartialCoverage,
                severity: QualityFlagSeverity::Warning,
                message: format!(
                    "{} pixels ({:.1}% of total) had no DEM coverage and are set to NaN.",
                    dem_missing, pct
                ),
            });
        }
    }
    if let Some(nc) = non_converged_count {
        if nc > 0 {
            flags.push(QualityFlagEntry {
                flag_type: QualityFlagType::SolverNonConvergence,
                severity: QualityFlagSeverity::Warning,
                message: format!(
                    "{} pixels' zero-Doppler iteration did not converge and are set to NaN.",
                    nc
                ),
            });
        }
    }
    flags
}

/// Build the sidecar output path by stripping the trailing extension of
/// `output` and appending `suffix`.  Example: `out.tif` + `.mask.tif` →
/// `out.mask.tif`; `out` (no extension) + `.lia.tif` → `out.lia.tif`.
pub fn sidecar_path(output: &Path, suffix: &str) -> Result<String> {
    let mut base = output.to_path_buf();
    if base.extension().is_some() {
        base.set_extension("");
    }
    let mut s = base
        .to_str()
        .ok_or_else(|| anyhow!("output path contains non-UTF-8 characters"))?
        .to_owned();
    s.push_str(suffix);
    Ok(s)
}

/// Resolve the `--crs <SPEC>` argument into an [`OutputCrs`].
pub fn resolve_crs(
    spec: &str,
    scene_centre_lon_lat: Option<(f64, f64)>,
) -> Result<OutputCrs> {
    let trimmed = spec.trim();
    if trimmed.is_empty() {
        bail!(
            "--crs SPEC must not be empty.  Use 'wgs84' (default), 'auto', \
             'EPSG:4326', or 'EPSG:326NN' (UTM)."
        );
    }
    if trimmed.eq_ignore_ascii_case("auto") {
        let (lon, lat) = scene_centre_lon_lat.ok_or_else(|| {
            anyhow!("--crs auto requires a scene centre but none was provided")
        })?;
        return OutputCrs::auto_for_lon_lat(lon, lat).map_err(|e| {
            anyhow!("--crs auto: cannot pick UTM zone for ({lon:.4}, {lat:.4}): {e}")
        });
    }
    if trimmed.eq_ignore_ascii_case("wgs84") {
        return Ok(OutputCrs::Wgs84LatLon);
    }
    OutputCrs::from_spec(trimmed).map_err(|e| anyhow!("--crs {trimmed:?}: {e}"))
}

/// Resolve the `--speckle` knobs into a [`SpeckleFilter`] enum, or
/// `Ok(None)` for `none`.
pub fn resolve_speckle(
    kind: &str,
    window: usize,
    enl: f32,
    damping: f32,
) -> Result<Option<crate::speckle::SpeckleFilter>> {
    use crate::speckle::SpeckleFilter;
    match kind.trim().to_ascii_lowercase().as_str() {
        "none" => Ok(None),
        "boxcar" => Ok(Some(SpeckleFilter::Boxcar { window })),
        "lee" => Ok(Some(SpeckleFilter::Lee { window, enl })),
        "frost" => Ok(Some(SpeckleFilter::Frost { window, damping })),
        "gamma_map" | "gamma-map" | "gammamap" => {
            Ok(Some(SpeckleFilter::GammaMap { window, enl }))
        }
        "refined_lee" | "refined-lee" | "refinedlee" => {
            Ok(Some(SpeckleFilter::RefinedLee { enl }))
        }
        other => bail!(
            "--speckle: unknown filter '{other}'.  Expected one of: \
             none, boxcar, lee, frost, gamma_map, refined_lee."
        ),
    }
}

/// Apply a speckle filter in-place on a linear-power raster.
pub fn apply_speckle_step(
    data: &mut Vec<f32>,
    cols: usize,
    rows: usize,
    choice: Option<crate::speckle::SpeckleFilter>,
) -> Result<()> {
    let Some(filter) = choice else {
        return Ok(());
    };
    tracing::info!("applying speckle filter {filter:?} …");
    let out = crate::speckle::apply_speckle_filter(data, cols, rows, filter)
        .with_context(|| format!("speckle filter {filter:?}"))?;
    *data = out;
    Ok(())
}

/// Project the in-effect speckle configuration onto the four optional
/// provenance fields.
pub fn speckle_provenance_fields(
    choice: Option<crate::speckle::SpeckleFilter>,
) -> (Option<String>, Option<usize>, Option<f32>, Option<f32>) {
    use crate::speckle::SpeckleFilter;
    match choice {
        None => (None, None, None, None),
        Some(SpeckleFilter::Boxcar { window }) => {
            (Some("boxcar".to_owned()), Some(window), None, None)
        }
        Some(SpeckleFilter::Lee { window, enl }) => {
            (Some("lee".to_owned()), Some(window), Some(enl), None)
        }
        Some(SpeckleFilter::Frost { window, damping }) => (
            Some("frost".to_owned()),
            Some(window),
            None,
            Some(damping),
        ),
        Some(SpeckleFilter::GammaMap { window, enl }) => {
            (Some("gamma_map".to_owned()), Some(window), Some(enl), None)
        }
        Some(SpeckleFilter::RefinedLee { enl }) => {
            (Some("refined_lee".to_owned()), Some(7), Some(enl), None)
        }
    }
}

/// Resolve the `--geoid <SPEC>` argument into a [`GeoidModel`].
pub fn resolve_geoid(spec: &str) -> Result<crate::geoid::GeoidModel> {
    use crate::geoid::{Egm96Grid, GeoidModel};

    let trimmed = spec.trim();
    if trimmed.is_empty() {
        bail!(
            "--geoid SPEC must not be empty.  Use 'auto' to fetch EGM96 \
             (requires --features geoid-fetch), 'egm2008' to fetch EGM2008, \
             'zero' to disable geoid correction explicitly, \
             or a path to an EGM96 binary/ASCII grid file."
        );
    }

    match trimmed.to_ascii_lowercase().as_str() {
        "zero" => {
        tracing::warn!(
            "geoid disabled (--geoid zero). \
             Expect up to ±80 m systematic vertical bias and a horizontal \
             geolocation error of similar magnitude.  Use 'auto' or a \
             grid path for production."
        );
            Ok(GeoidModel::Zero)
        }
        "auto" => {
            #[cfg(feature = "geoid-fetch")]
            {
                use crate::geoid_fetch::fetch_egm96;
                tracing::info!(
                    "locating or downloading EGM96 grid \
                     ($SARDINE_GEOID_DIR or ~/.sardine/geoid/) …"
                );
                let grid = fetch_egm96()
                    .with_context(|| "fetching EGM96 grid for --geoid auto")?;
                Ok(GeoidModel::Egm96(grid))
            }
            #[cfg(not(feature = "geoid-fetch"))]
            {
                bail!(
                    "--geoid auto requires the 'geoid-fetch' Cargo feature.\n\
                     Rebuild with: cargo build --release --features geoid-fetch\n\
                     Or pass an explicit grid path: --geoid /path/to/egm96.bin"
                );
            }
        }
        "egm2008" => {
            #[cfg(feature = "geoid-fetch")]
            {
                use crate::geoid_fetch::fetch_egm2008;
                tracing::info!(
                    "locating or downloading EGM2008 grid \
                     ($SARDINE_GEOID_DIR or ~/.sardine/geoid/) …"
                );
                let grid = fetch_egm2008()
                    .with_context(|| "fetching EGM2008 grid for --geoid egm2008")?;
                Ok(GeoidModel::Egm2008(grid))
            }
            #[cfg(not(feature = "geoid-fetch"))]
            {
                bail!(
                    "--geoid egm2008 requires the 'geoid-fetch' Cargo feature.\n\
                     Rebuild with: cargo build --release --features geoid-fetch\n\
                     Or pass an explicit grid path: --geoid /path/to/egm2008.bin"
                );
            }
        }
        _ => {
            let path = std::path::Path::new(trimmed);
            if !path.exists() {
                bail!(
                    "--geoid: '{}' is neither 'auto', 'zero', nor an existing file path",
                    trimmed
                );
            }
            let ext = path
                .extension()
                .and_then(|e| e.to_str())
                .map(|s| s.to_ascii_lowercase());
            let grid = match ext.as_deref() {
                Some("bin") => Egm96Grid::load_binary(path).with_context(|| {
                    format!("loading EGM96 binary cache: {}", path.display())
                })?,
                Some("gtx") | Some("grd") => Egm96Grid::from_ww15mgh(path).with_context(|| {
                    format!("loading NGA WW15MGH grid: {}", path.display())
                })?,
                _ => bail!(
                    "--geoid path '{}' has unrecognised extension; expected \
                     .bin (compact cache) or .gtx / .GRD (NGA WW15MGH)",
                    path.display()
                ),
            };
            Ok(GeoidModel::Egm96(grid))
        }
    }
}

/// Find a measurement TIFF in `dir` whose name contains both `iw` and `pol`.
fn find_measurement_tiff(dir: &Path, iw: &str, pol: &str) -> Result<PathBuf> {
    let entries = std::fs::read_dir(dir)
        .with_context(|| format!("reading measurement dir: {}", dir.display()))?;
    for entry in entries {
        let entry = entry?;
        let name = entry.file_name().to_string_lossy().to_lowercase();
        if name.contains(iw) && name.contains(pol) && name.ends_with(".tiff") {
            return Ok(entry.path());
        }
    }
    bail!(
        "No measurement TIFF found for iw={} pol={} in {}",
        iw,
        pol,
        dir.display()
    )
}

// ─────────────────────────────────────────────────────────────────────────────
// Shared scene preparation
// ─────────────────────────────────────────────────────────────────────────────

/// Output of [`prepare_merged_scene`].
pub struct PreparedScene {
    pub scene: crate::types::SceneMetadata,
    pub orbit_source: crate::provenance::OrbitSource,
    pub merged: crate::merge_subswaths::MergedSigma0,
    pub grids: Vec<(
        crate::types::SubSwathId,
        Vec<crate::types::GeolocationGridPoint>,
    )>,
    pub pol: Polarization,
}

/// Run the steps shared by `process` and `grd`: configure Rayon, parse SAFE
/// metadata, apply the precise (or annotation) orbit, parse calibration/noise
/// LUTs, parse geolocation grids, deburst + calibrate every IW subswath in
/// parallel, and merge the subswaths into one slant-range σ⁰ image.
///
/// `iw_selection` controls which sub-swaths are processed and optionally
/// limits the burst range within each.  Pass `&IwSelection::default()` to
/// process all three IW subswaths with all bursts (the original behaviour).
pub fn prepare_merged_scene(
    safe: &Path,
    orbit_path: Option<&Path>,
    polarization: &str,
    threads: usize,
    iw_selection: &IwSelection,
) -> Result<PreparedScene> {
    use crate::apply_calibration::apply_calibration;
    use crate::calibration::parse_calibration_noise;
    use crate::deburst::deburst_subswath;
    use crate::merge_subswaths::{merge_single_subswath_owned, merge_subswaths, SwathInput};
    use crate::orbit::{apply_precise_orbit, parse_eof_file};
    use crate::parse::{parse_geolocation_grids, parse_safe_directory};
    use crate::slc_reader::SlcReader;
    use crate::types::SubSwathId;

    if threads > 0 {
        rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build_global()
            .with_context(|| format!("setting Rayon thread count to {threads}"))?;
        tracing::info!("using {threads} thread(s)");
    } else {
        tracing::info!("using all available CPU threads (Rayon default)");
    }

    let pol = match polarization.to_uppercase().as_str() {
        "VV" => Polarization::VV,
        "VH" => Polarization::VH,
        other => bail!("Unsupported polarization: {}. Use VV or VH.", other),
    };

    tracing::info!("parsing SAFE metadata …");
    let t_parse = Instant::now();
    let scene = parse_safe_directory(safe)
        .with_context(|| format!("parsing SAFE: {}", safe.display()))?;

    let (scene, orbit_source) = if let Some(orbit_path) = orbit_path {
        let orbit = parse_eof_file(orbit_path)
            .with_context(|| format!("parsing orbit file: {}", orbit_path.display()))?;
        let scene = apply_precise_orbit(scene, &orbit)
            .with_context(|| "applying precise orbit")?;
        (scene, crate::provenance::OrbitSource::Poeorb)
    } else {
        let allow = std::env::var("SARDINE_ALLOW_ANNOTATION_ORBIT")
            .map(|v| v == "1")
            .unwrap_or(false); // SAFETY-OK: env var parse; fallback to deny is the safe direction
        if !allow {
            bail!(
                "No --orbit file provided and SARDINE_ALLOW_ANNOTATION_ORBIT != 1.\n\
                 Provide a POEORB file with --orbit, or set SARDINE_ALLOW_ANNOTATION_ORBIT=1 \
                 to use the annotation orbit (metre-level accuracy)."
            );
        }
        tracing::warn!(
            "using annotation orbit (metre-level accuracy). \
             Provide a POEORB file for cm-level accuracy."
        );
        (scene, crate::provenance::OrbitSource::Annotation)
    };

    tracing::info!("parsing calibration/noise LUTs …");
    let cal_noise = parse_calibration_noise(safe)
        .with_context(|| format!("parsing calibration/noise: {}", safe.display()))?;

    tracing::info!("parsing geolocation grids …");
    let grids = parse_geolocation_grids(safe)
        .with_context(|| format!("parsing geolocation grids: {}", safe.display()))?;
    report_timing("parse_safe_meta", t_parse);

    let all_iw_ids = [SubSwathId::IW1, SubSwathId::IW2, SubSwathId::IW3];
    let iw_ids: Vec<SubSwathId> = if iw_selection.subswaths.is_empty() {
        all_iw_ids.to_vec()
    } else {
        let mut ids = iw_selection.subswaths.clone();
        ids.sort();
        ids
    };

    // Warn if selected subswaths are non-contiguous (e.g. IW1+IW3 skips IW2).
    if iw_ids.len() >= 2 {
        let selected_indices: Vec<usize> = iw_ids
            .iter()
            .filter_map(|id| all_iw_ids.iter().position(|x| x == id))
            .collect();
        if selected_indices.windows(2).any(|w| w[1] - w[0] > 1) {
            let names: Vec<String> = iw_ids.iter().map(|id| format!("{id:?}")).collect();
            tracing::warn!(
                "selected subswaths {} are non-contiguous \
                 — output will contain a NaN gap between them",
                names.join(", ")
            );
        }
    }

    let measurement_dir = safe.join("measurement");

    struct SwathData {
        swath: crate::types::SubSwathMetadata,
        sigma0: crate::apply_calibration::Sigma0Array,
        azimuth_start_time: chrono::DateTime<chrono::Utc>,
    }

    let t_swaths = Instant::now();
    let mut swath_data: Vec<SwathData> = iw_ids
        .par_iter()
        .map(|&iw_id| -> Result<SwathData> {
            let iw_str = format!("{}", iw_id).to_lowercase();
            let pol_str = format!("{}", pol).to_lowercase();
            let tiff_path = find_measurement_tiff(&measurement_dir, &iw_str, &pol_str)?;

            let swath_meta = scene
                .sub_swaths
                .iter()
                .find(|s| s.id == iw_id)
                .with_context(|| format!("subswath {:?} not found in scene metadata", iw_id))?
                .clone();

            let bursts: Vec<_> = {
                let mut b: Vec<_> = scene
                    .bursts
                    .iter()
                    .filter(|b| b.subswath_id == iw_id)
                    .cloned()
                    .collect();
                b.sort_by_key(|b| b.burst_index);
                // Apply burst range selection if requested.
                if let Some((start, end)) = iw_selection.burst_range {
                    let end_excl = (end + 1).min(b.len());
                    let start = start.min(b.len());
                    b = b[start..end_excl].to_vec();
                }
                b
            };
            if bursts.is_empty() {
                return Err(anyhow::anyhow!(
                    "no bursts found for subswath {:?} in scene metadata",
                    iw_id
                ));
            }
            let azimuth_start_time = bursts[0].azimuth_time_utc;
            let tiff_line_origin = bursts[0].first_line;

            let cal = cal_noise
                .calibrations
                .iter()
                .find(|c| c.subswath_id == iw_id && c.polarization == pol)
                .with_context(|| {
                    format!("calibration LUT not found for {:?} {:?}", iw_id, pol)
                })?;
            let noise = cal_noise
                .noises
                .iter()
                .find(|n| n.subswath_id == iw_id && n.polarization == pol)
                .with_context(|| format!("noise LUT not found for {:?} {:?}", iw_id, pol))?;

            tracing::info!("debursting {:?} {} …", iw_id, pol_str.to_uppercase());
            let mut reader = SlcReader::open(&tiff_path)
                .with_context(|| format!("opening SLC TIFF: {}", tiff_path.display()))?;

            let deburst = deburst_subswath(&mut reader, &swath_meta, &bursts)
                .with_context(|| format!("debursting {:?}", iw_id))?;

            tracing::info!("calibrating {:?} {} …", iw_id, pol_str.to_uppercase());
            let sigma0 = apply_calibration(&deburst, cal, noise, tiff_line_origin)
                .with_context(|| format!("calibration for {:?}", iw_id))?;

            if sigma0.cal_lut_extrapolation_gap_px > 0
                || sigma0.noise_lut_extrapolation_gap_px > 0
            {
                tracing::warn!(
                    "{:?} {} LUT extrapolation gap: cal={} px, noise={} px",
                    iw_id,
                    pol_str.to_uppercase(),
                    sigma0.cal_lut_extrapolation_gap_px,
                    sigma0.noise_lut_extrapolation_gap_px,
                );
            }

            Ok(SwathData {
                swath: swath_meta,
                sigma0,
                azimuth_start_time,
            })
        })
        .collect::<Result<Vec<_>>>()?;
    report_timing("read_cal_deburst_subswaths", t_swaths);

    tracing::info!("merging subswaths …");
    let t_merge = Instant::now();
    let merged = if swath_data.len() == 1 {
        let sd = swath_data.remove(0);
        merge_single_subswath_owned(sd.sigma0, &sd.swath, sd.azimuth_start_time)
            .with_context(|| "merging subswaths")?
    } else {
        let swath_inputs: Vec<SwathInput<'_>> = swath_data
            .iter()
            .map(|sd| SwathInput {
                sigma0: &sd.sigma0,
                swath: &sd.swath,
                azimuth_start_time: sd.azimuth_start_time,
            })
            .collect();
        merge_subswaths(&swath_inputs)
            .with_context(|| "merging subswaths")?
    };
    report_timing("merge_subswaths", t_merge);

    if merged.cal_lut_extrapolation_gap_px > 0 || merged.noise_lut_extrapolation_gap_px > 0 {
        tracing::warn!(
            "merged LUT extrapolation gap: cal={} px, noise={} px",
            merged.cal_lut_extrapolation_gap_px,
            merged.noise_lut_extrapolation_gap_px,
        );
    }

    Ok(PreparedScene {
        scene,
        orbit_source,
        merged,
        grids,
        pol,
    })
}

// ─────────────────────────────────────────────────────────────────────────────
// Multi-slice (assembled) scene preparation
// ─────────────────────────────────────────────────────────────────────────────

/// Stitch per-slice [`Sigma0Array`]s into one, applying midpoint seam trimming
/// in the calibrated power domain.
///
/// # Seam trimming
///
/// At boundary between slice k and slice k+1 with overlap `ov` lines:
/// - slice k loses its last `⌊ov/2⌋` lines
/// - slice k+1 loses its first `⌈ov/2⌉` lines
///
/// This is the same midpoint selection as [`crate::deburst::deburst_subswath`]
/// but applied to already-calibrated σ⁰ values instead of CInt16 samples.
/// For incoherent power (σ⁰ = |DN|²/K²) this is equivalent because
/// calibration is a per-pixel multiplicative operation (no spatial mixing).
///
/// # Errors
///
/// - If `slices` is empty.
/// - If any two slices have different `samples` or `valid_sample_offset`.
/// - If a seam overlap is so large that a slice would contribute zero or
///   fewer lines after trimming.
fn stitch_sigma0(
    slices: &[crate::apply_calibration::Sigma0Array],
    seam_overlaps: &[usize],
) -> Result<crate::apply_calibration::Sigma0Array> {
    if slices.is_empty() {
        bail!("stitch_sigma0: empty slice list");
    }
    if seam_overlaps.len() != slices.len() - 1 {
        bail!(
            "stitch_sigma0: {} slices but {} seam overlaps (expected {})",
            slices.len(),
            seam_overlaps.len(),
            slices.len() - 1
        );
    }

    let samples = slices[0].samples;
    let vso = slices[0].valid_sample_offset;
    for (k, s) in slices.iter().enumerate().skip(1) {
        if s.samples != samples {
            bail!(
                "stitch_sigma0: slice {k} has {} samples != slice 0 {} samples; \
                 all slices must have the same valid range window",
                s.samples, samples
            );
        }
        if s.valid_sample_offset != vso {
            bail!(
                "stitch_sigma0: slice {k} valid_sample_offset {} != slice 0 {}; \
                 cannot concatenate",
                s.valid_sample_offset, vso
            );
        }
    }

    if slices.len() == 1 {
        let s = &slices[0];
        return Ok(crate::apply_calibration::Sigma0Array {
            data: s.data.clone(),
            lines: s.lines,
            samples: s.samples,
            valid_sample_offset: s.valid_sample_offset,
            cal_lut_extrapolation_gap_px: s.cal_lut_extrapolation_gap_px,
            noise_lut_extrapolation_gap_px: s.noise_lut_extrapolation_gap_px,
            nesz: s.nesz.clone(),
        });
    }

    let n = slices.len();
    let mut trim_start = vec![0usize; n];
    let mut trim_end = vec![0usize; n];
    for k in 0..n - 1 {
        let ov = seam_overlaps[k];
        trim_end[k] = ov / 2;
        trim_start[k + 1] = (ov + 1) / 2;
    }

    let total_lines: usize = slices
        .iter()
        .enumerate()
        .map(|(k, s)| s.lines.saturating_sub(trim_start[k] + trim_end[k]))
        .sum();
    let mut data = Vec::with_capacity(total_lines * samples);
    let mut nesz = Vec::with_capacity(total_lines * samples);
    let mut max_cal_gap: u32 = 0;
    let mut max_nr_gap: u32 = 0;

    for (k, s) in slices.iter().enumerate() {
        let start_line = trim_start[k];
        let end_line = s.lines.saturating_sub(trim_end[k]);
        if end_line <= start_line {
            bail!(
                "stitch_sigma0: slice {k} has only {} lines after trimming \
                 ({} start + {} end); seam overlap too large",
                end_line as i64 - start_line as i64,
                trim_start[k],
                trim_end[k]
            );
        }
        data.extend_from_slice(&s.data[start_line * samples..end_line * samples]);
        nesz.extend_from_slice(&s.nesz[start_line * samples..end_line * samples]);
        max_cal_gap = max_cal_gap.max(s.cal_lut_extrapolation_gap_px);
        max_nr_gap = max_nr_gap.max(s.noise_lut_extrapolation_gap_px);
    }

    Ok(crate::apply_calibration::Sigma0Array {
        data,
        lines: total_lines,
        samples,
        valid_sample_offset: vso,
        cal_lut_extrapolation_gap_px: max_cal_gap,
        noise_lut_extrapolation_gap_px: max_nr_gap,
        nesz,
    })
}

/// Merge geolocation grids from all slices of an assembled scene.
///
/// For each subswath, collects grid points from every slice and sorts them
/// by azimuth time.  The `line` field is not used by terrain correction so
/// no offset correction is applied.
fn merge_geolocation_grids(
    all_grids: Vec<Vec<(crate::types::SubSwathId, Vec<crate::types::GeolocationGridPoint>)>>,
) -> Vec<(crate::types::SubSwathId, Vec<crate::types::GeolocationGridPoint>)> {
    use crate::types::SubSwathId;
    let ids = [SubSwathId::IW1, SubSwathId::IW2, SubSwathId::IW3];
    ids.iter()
        .filter_map(|&id| {
            let mut pts: Vec<crate::types::GeolocationGridPoint> = all_grids
                .iter()
                .flat_map(|slice_grids| {
                    slice_grids
                        .iter()
                        .find(|(sw, _)| *sw == id)
                        .map(|(_, pts)| pts.as_slice())
                        .unwrap_or(&[]) // SAFETY-OK: absent subswath in a slice contributes no grid points; empty slice is correct design intent
                })
                .cloned()
                .collect();
            if pts.is_empty() {
                return None;
            }
            pts.sort_by_key(|p| p.azimuth_time_utc);
            Some((id, pts))
        })
        .collect()
}

/// Variant of [`prepare_merged_scene`] that accepts a multi-slice
/// [`AssembledScene`] produced by [`crate::slice_assembly::assemble_slices`].
///
/// # Calibration correctness
///
/// Calibration LUT line indices in S-1 annotation files are TIFF-local
/// (relative to that slice's measurement TIFF).  This function processes each
/// slice independently — deburst + calibrate — then stitches the resulting
/// calibrated σ⁰ arrays using midpoint seam trimming at each inter-slice
/// boundary.  This avoids any TIFF-line offset translation in
/// [`crate::apply_calibration::apply_calibration`].
///
/// # Seam handling
///
/// At each inter-slice boundary the burst overlap is computed from azimuth
/// time differences via [`crate::deburst::compute_overlap_lines`] and midpoint
/// selection is applied in the calibrated power domain, which is equivalent
/// to applying it in the raw CInt16 domain for incoherent σ⁰ processing.
pub fn prepare_merged_scene_assembled(
    assembled: &crate::slice_assembly::AssembledScene,
    orbit_path: Option<&Path>,
    polarization: &str,
    threads: usize,
    iw_selection: &IwSelection,
) -> Result<PreparedScene> {
    use crate::apply_calibration::apply_calibration;
    use crate::calibration::parse_calibration_noise;
    use crate::deburst::{compute_overlap_lines, deburst_subswath};
    use crate::merge_subswaths::{merge_single_subswath_owned, merge_subswaths, SwathInput};
    use crate::orbit::{apply_precise_orbit, parse_eof_file};
    use crate::parse::parse_geolocation_grids;
    use crate::slc_reader::SlcReader;
    use crate::types::SubSwathId;

    let n_slices = assembled.safe_paths.len();

    if threads > 0 {
        rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build_global()
            .with_context(|| format!("setting Rayon thread count to {threads}"))?;
        tracing::info!("using {threads} thread(s)");
    } else {
        tracing::info!("using all available CPU threads (Rayon default)");
    }

    let pol = match polarization.to_uppercase().as_str() {
        "VV" => Polarization::VV,
        "VH" => Polarization::VH,
        other => bail!("Unsupported polarization: {}. Use VV or VH.", other),
    };

    let scene = assembled.scene.clone();
    let (scene, orbit_source) = if let Some(orbit_path) = orbit_path {
        let orbit = parse_eof_file(orbit_path)
            .with_context(|| format!("parsing orbit file: {}", orbit_path.display()))?;
        let scene = apply_precise_orbit(scene, &orbit)
            .with_context(|| "applying precise orbit")?;
        (scene, crate::provenance::OrbitSource::Poeorb)
    } else {
        let allow = std::env::var("SARDINE_ALLOW_ANNOTATION_ORBIT")
            .map(|v| v == "1")
            .unwrap_or(false); // SAFETY-OK: env var parse; fallback to deny is the safe direction
        if !allow {
            bail!(
                "No --orbit file provided and SARDINE_ALLOW_ANNOTATION_ORBIT != 1.\n\
                 Provide a POEORB file with --orbit, or set SARDINE_ALLOW_ANNOTATION_ORBIT=1 \
                 to use the annotation orbit (metre-level accuracy)."
            );
        }
        tracing::warn!(
            "using annotation orbit (metre-level accuracy). \
             Provide a POEORB file for cm-level accuracy."
        );
        (scene, crate::provenance::OrbitSource::Annotation)
    };

    // Geolocation grids: parse from every slice and merge.
    tracing::info!("parsing geolocation grids from {n_slices} slices …");
    let all_grids: Vec<_> = assembled
        .safe_paths
        .iter()
        .enumerate()
        .map(|(k, p)| {
            parse_geolocation_grids(p)
                .with_context(|| format!("geolocation grids slice {k}: {}", p.display()))
        })
        .collect::<Result<Vec<_>>>()?;
    let grids = merge_geolocation_grids(all_grids);

    let all_iw_ids = [SubSwathId::IW1, SubSwathId::IW2, SubSwathId::IW3];
    let iw_ids: Vec<SubSwathId> = if iw_selection.subswaths.is_empty() {
        all_iw_ids.to_vec()
    } else {
        let mut ids = iw_selection.subswaths.clone();
        ids.sort();
        ids
    };

    // Warn if selected subswaths are non-contiguous.
    if iw_ids.len() >= 2 {
        let selected_indices: Vec<usize> = iw_ids
            .iter()
            .filter_map(|id| all_iw_ids.iter().position(|x| x == id))
            .collect();
        if selected_indices.windows(2).any(|w| w[1] - w[0] > 1) {
            let names: Vec<String> = iw_ids.iter().map(|id| format!("{id:?}")).collect();
            tracing::warn!(
                "selected subswaths {} are non-contiguous \
                 — output will contain a NaN gap between them",
                names.join(", ")
            );
        }
    }

    struct SwathData {
        swath: crate::types::SubSwathMetadata,
        sigma0: crate::apply_calibration::Sigma0Array,
        azimuth_start_time: chrono::DateTime<chrono::Utc>,
    }

    let t_swaths = Instant::now();
    let mut swath_data: Vec<SwathData> = iw_ids
        .par_iter()
        .map(|&iw_id| -> Result<SwathData> {
            let iw_str = format!("{iw_id}").to_lowercase();
            let pol_str = format!("{pol}").to_lowercase();

            // Assembled subswath metadata (covers all slices).
            let assembled_sw = scene
                .sub_swaths
                .iter()
                .find(|s| s.id == iw_id)
                .with_context(|| {
                    format!("subswath {iw_id:?} not found in assembled scene")
                })?
                .clone();

            let ati = assembled_sw.azimuth_time_interval_s;
            let lps = assembled_sw.lines_per_burst;

            // Per-slice deburst + calibrate.
            let mut slice_sigma0: Vec<crate::apply_calibration::Sigma0Array> =
                Vec::with_capacity(n_slices);
            let mut last_bursts_per_slice: Vec<crate::types::BurstEntry> =
                Vec::with_capacity(n_slices);
            let mut first_bursts_per_slice: Vec<crate::types::BurstEntry> =
                Vec::with_capacity(n_slices);

            for slice_idx in 0..n_slices {
                let safe = &assembled.safe_paths[slice_idx];

                // Bursts belonging to (iw_id, slice_idx) in the assembled scene.
                let global_bursts: Vec<crate::types::BurstEntry> = scene
                    .bursts
                    .iter()
                    .filter(|b| b.subswath_id == iw_id && b.slice_index == slice_idx)
                    .cloned()
                    .collect();
                if global_bursts.is_empty() {
                    bail!("no bursts for {iw_id:?} in slice {slice_idx}");
                }

                // The first burst's first_line encodes the logical offset added by
                // merge_bursts() in slice_assembly.  Subtracting it gives the
                // TIFF-local first_line (= burst_index_within_slice * lines_per_burst).
                let slice_line_offset = global_bursts[0].first_line;

                // Build TIFF-local burst list.
                let local_bursts: Vec<crate::types::BurstEntry> = global_bursts
                    .iter()
                    .map(|b| crate::types::BurstEntry {
                        first_line: b.first_line - slice_line_offset,
                        last_line: b.last_line - slice_line_offset,
                        burst_index: b.burst_index - global_bursts[0].burst_index,
                        slice_index: 0, // SAFETY-OK: local burst list is always single-slice context
                        ..*b
                    })
                    .collect();

                let n_k = local_bursts.len();
                let sw_k = crate::types::SubSwathMetadata {
                    id: iw_id,
                    burst_count: n_k,
                    azimuth_samples: n_k * lps,
                    first_line: 0,
                    last_line: n_k * lps,
                    ..assembled_sw.clone()
                };

                // Calibration LUTs for this slice.
                let cal_noise = parse_calibration_noise(safe).with_context(|| {
                    format!("calibration slice {slice_idx}: {}", safe.display())
                })?;
                let cal = cal_noise
                    .calibrations
                    .iter()
                    .find(|c| c.subswath_id == iw_id && c.polarization == pol)
                    .with_context(|| {
                        format!(
                            "no calibration LUT for {iw_id:?} {pol:?} in slice {slice_idx}"
                        )
                    })?;
                let noise = cal_noise
                    .noises
                    .iter()
                    .find(|n| n.subswath_id == iw_id && n.polarization == pol)
                    .with_context(|| {
                        format!(
                            "no noise LUT for {iw_id:?} {pol:?} in slice {slice_idx}"
                        )
                    })?;

                tracing::info!(
                    "debursting {iw_id:?} {} slice {slice_idx} …",
                    pol_str.to_uppercase()
                );
                let mdir = safe.join("measurement");
                let tiff = find_measurement_tiff(&mdir, &iw_str, &pol_str)?;
                let mut reader = SlcReader::open(&tiff)
                    .with_context(|| format!("opening SLC TIFF: {}", tiff.display()))?;

                let tiff_line_origin = local_bursts[0].first_line; // = 0 for all slices

                let deburst = deburst_subswath(&mut reader, &sw_k, &local_bursts)
                    .with_context(|| {
                        format!("debursting {iw_id:?} slice {slice_idx}")
                    })?;

                tracing::info!(
                    "calibrating {iw_id:?} {} slice {slice_idx} …",
                    pol_str.to_uppercase()
                );
                let sigma0 = apply_calibration(&deburst, cal, noise, tiff_line_origin)
                    .with_context(|| {
                        format!("calibration {iw_id:?} slice {slice_idx}")
                    })?;

                first_bursts_per_slice.push(global_bursts[0].clone());
                last_bursts_per_slice.push(
                    global_bursts
                        .last() // SAFETY-OK: global_bursts.is_empty() checked above
                        .expect("non-empty burst list")
                        .clone(),
                );
                slice_sigma0.push(sigma0);
            }

            // Compute seam overlaps between consecutive slices.
            let mut seam_overlaps: Vec<usize> = Vec::with_capacity(n_slices.saturating_sub(1));
            for k in 0..n_slices.saturating_sub(1) {
                let ov = compute_overlap_lines(
                    &last_bursts_per_slice[k],
                    &first_bursts_per_slice[k + 1],
                    ati,
                    lps,
                )
                .with_context(|| {
                    format!("seam overlap {iw_id:?} slice {k}→{}", k + 1)
                })?;
                seam_overlaps.push(ov);
            }

            // Stitch calibrated slices (midpoint selection in power domain).
            let stitched = stitch_sigma0(&slice_sigma0, &seam_overlaps)
                .with_context(|| format!("stitching {iw_id:?}"))?;

            let azimuth_start_time = first_bursts_per_slice[0].azimuth_time_utc;

            if stitched.cal_lut_extrapolation_gap_px > 0
                || stitched.noise_lut_extrapolation_gap_px > 0
            {
                tracing::warn!(
                    "{iw_id:?} {} LUT extrapolation gap: cal={} px, noise={} px",
                    pol_str.to_uppercase(),
                    stitched.cal_lut_extrapolation_gap_px,
                    stitched.noise_lut_extrapolation_gap_px,
                );
            }

            Ok(SwathData {
                swath: assembled_sw,
                sigma0: stitched,
                azimuth_start_time,
            })
        })
        .collect::<Result<Vec<_>>>()?;
    report_timing("assembled_read_cal_deburst", t_swaths);

    tracing::info!("merging assembled subswaths …");
    let t_merge = Instant::now();
    let merged = if swath_data.len() == 1 {
        let sd = swath_data.remove(0);
        merge_single_subswath_owned(sd.sigma0, &sd.swath, sd.azimuth_start_time)
            .with_context(|| "merging assembled subswaths")?
    } else {
        let swath_inputs: Vec<SwathInput<'_>> = swath_data
            .iter()
            .map(|sd| SwathInput {
                sigma0: &sd.sigma0,
                swath: &sd.swath,
                azimuth_start_time: sd.azimuth_start_time,
            })
            .collect();
        merge_subswaths(&swath_inputs)
            .with_context(|| "merging assembled subswaths")?
    };
    report_timing("assembled_merge_subswaths", t_merge);

    if merged.cal_lut_extrapolation_gap_px > 0 || merged.noise_lut_extrapolation_gap_px > 0 {
        tracing::warn!(
            "assembled merged LUT extrapolation gap: cal={} px, noise={} px",
            merged.cal_lut_extrapolation_gap_px,
            merged.noise_lut_extrapolation_gap_px,
        );
    }

    Ok(PreparedScene {
        scene,
        orbit_source,
        merged,
        grids,
        pol,
    })
}

// ─────────────────────────────────────────────────────────────────────────────
// run_process
// ─────────────────────────────────────────────────────────────────────────────

// ─────────────────────────────────────────────────────────────────────────────
// Pre-TC multilook
// ─────────────────────────────────────────────────────────────────────────────

/// Apply a boxcar power-domain multilook to a merged σ⁰ image in slant-range
/// geometry.
///
/// This implements the **Multilook** step of the agreed processing chain:
///
/// ```text
/// Deburst → Calibration → Merge → Multilook → TF → TC → Speckle → dB
/// ```
///
/// Output pixel `(ol, os)` is the average of the `azimuth_looks × range_looks`
/// input pixels in block `[ol*A, (ol+1)*A) × [os*R, (os+1)*R)`.  NaN pixels
/// (seam gaps, swath edges) are excluded from the average.  If *all* pixels in
/// a block are NaN the output is NaN.
///
/// # Coordinate update
///
/// Output pixels are centre-aligned to their input blocks:
///
/// - Range: `new_near_t = old_near_t + (R−1)/2 × old_range_spacing_s`
/// - Azimuth: `new_az_start = old_az_start + (A−1)/2 × ati_s`
///
/// The downstream `RadarGeometry::from_scene_and_merged` reads
/// `near_slant_range_time_s`, `range_pixel_spacing_m`, `azimuth_start_time`,
/// `lines`, and `samples` from the returned struct, so this update is the
/// only place that needs to change for TC to operate correctly on the
/// multilooked image.
///
/// # Arguments
///
/// * `merged`                    — consumed input.
/// * `range_looks`               — range pixels averaged per output pixel (≥ 1).
/// * `azimuth_looks`             — azimuth lines averaged per output pixel (≥ 1).
/// * `azimuth_time_interval_s`   — ATI of the *input* image in seconds / line.
///   Only used when `azimuth_looks > 1`; pass 0.0 when `azimuth_looks == 1`.
pub(crate) fn apply_multilook(
    merged: crate::merge_subswaths::MergedSigma0,
    range_looks: usize,
    azimuth_looks: usize,
    azimuth_time_interval_s: f64,
) -> Result<crate::merge_subswaths::MergedSigma0> {
    use chrono::Duration;

    if range_looks == 1 && azimuth_looks == 1 {
        return Ok(merged);
    }
    if range_looks == 0 || azimuth_looks == 0 {
        bail!(
            "multilook factors must be ≥ 1 (range={range_looks}, azimuth={azimuth_looks})"
        );
    }

    let in_lines = merged.lines;
    let in_samples = merged.samples;
    let out_lines = in_lines / azimuth_looks;
    let out_samples = in_samples / range_looks;

    if out_lines == 0 || out_samples == 0 {
        bail!(
            "multilook ({range_looks}×{azimuth_looks}) is larger than \
             image dimensions ({in_samples} samples × {in_lines} lines)"
        );
    }

    // ── Update range geometry ────────────────────────────────────────────────
    // Centre-align: new pixel 0 corresponds to the centre of original pixels
    // [0, R), i.e. at offset (R−1)/2 from old pixel 0.
    let new_range_spacing_s = merged.range_pixel_spacing_s * range_looks as f64;
    let new_range_spacing_m = merged.range_pixel_spacing_m * range_looks as f64;
    let new_near_slant_range_time_s = merged.near_slant_range_time_s
        + (range_looks - 1) as f64 / 2.0 * merged.range_pixel_spacing_s;

    // ── Update azimuth geometry ──────────────────────────────────────────────
    let new_azimuth_start_time = if azimuth_looks > 1 {
        if !azimuth_time_interval_s.is_finite() || azimuth_time_interval_s <= 0.0 {
            bail!(
                "azimuth_time_interval_s must be > 0 when azimuth_looks > 1 \
                 (got {azimuth_time_interval_s})"
            );
        }
        // Centre of first azimuth block: offset (A−1)/2 lines forward.
        let half_block_us =
            ((azimuth_looks - 1) as f64 / 2.0 * azimuth_time_interval_s * 1e6) as i64;
        merged.azimuth_start_time + Duration::microseconds(half_block_us)
    } else {
        merged.azimuth_start_time
    };

    // ── Boxcar average (NaN-aware) ───────────────────────────────────────────
    let mut out_data = vec![0.0f32; out_lines * out_samples];
    let mut out_nesz = vec![0.0f32; out_lines * out_samples];

    for ol in 0..out_lines {
        for os in 0..out_samples {
            let mut sum = 0.0f64;
            let mut sum_nesz = 0.0f64;
            let mut count = 0u32;
            for da in 0..azimuth_looks {
                for dr in 0..range_looks {
                    let il = ol * azimuth_looks + da;
                    let is_ = os * range_looks + dr;
                    let v = merged.data[il * in_samples + is_];
                    if !v.is_nan() {
                        sum += v as f64;
                        sum_nesz += merged.nesz[il * in_samples + is_] as f64;
                        count += 1;
                    }
                }
            }
            out_data[ol * out_samples + os] = if count == 0 {
                f32::NAN
            } else {
                (sum / count as f64) as f32
            };
            out_nesz[ol * out_samples + os] = if count == 0 {
                f32::NAN
            } else {
                (sum_nesz / count as f64) as f32
            };
        }
    }

    Ok(crate::merge_subswaths::MergedSigma0 {
        data: out_data,
        lines: out_lines,
        samples: out_samples,
        near_slant_range_time_s: new_near_slant_range_time_s,
        range_pixel_spacing_s: new_range_spacing_s,
        range_pixel_spacing_m: new_range_spacing_m,
        azimuth_start_time: new_azimuth_start_time,
        // Propagate worst-case LUT gap flags unchanged — the multilook does
        // not affect LUT geometry.
        cal_lut_extrapolation_gap_px: merged.cal_lut_extrapolation_gap_px,
        noise_lut_extrapolation_gap_px: merged.noise_lut_extrapolation_gap_px,
        nesz: out_nesz,
    })
}

/// Run the full backscatter pipeline (deburst → calibration → merge →
/// optional multilook → terrain-correction → optional speckle → dB → GeoTIFF + sidecars).
///
/// This is the library equivalent of the `sardine process` CLI subcommand;
/// the byte output is identical for equivalent inputs.
pub fn run_process(opts: &ProcessOptions) -> Result<()> {
    use crate::dem::DemMosaic;
    use crate::export::{to_db_inplace, write_geotiff_with_crs};
    use crate::terrain_correction::{terrain_correction, TerrainCorrectionConfig};

    let t_total = Instant::now();
    let prepared = if !opts.extra_safe_paths.is_empty() {
        use crate::slice_assembly::assemble_slices;
        let mut all: Vec<&Path> = vec![opts.safe.as_path()];
        all.extend(opts.extra_safe_paths.iter().map(|p| p.as_path()));
        tracing::info!(
            "assembling {} slices …",
            all.len()
        );
        let assembled = assemble_slices(&all)
            .with_context(|| "assembling SAFE slices")?;
        prepare_merged_scene_assembled(
            &assembled,
            opts.orbit.as_deref(),
            &opts.polarization,
            opts.threads,
            &opts.iw_selection,
        )?
    } else {
        prepare_merged_scene(
            &opts.safe,
            opts.orbit.as_deref(),
            &opts.polarization,
            opts.threads,
            &opts.iw_selection,
        )?
    };
    let PreparedScene {
        scene,
        orbit_source,
        mut merged,
        grids,
        pol,
    } = prepared;

    // ── Pre-TC multilook ─────────────────────────────────────────────────────
    if opts.multilook_range > 1 || opts.multilook_azimuth > 1 {
        let ati = scene
            .sub_swaths
            .first()
            .ok_or_else(|| anyhow!("scene has no subswaths"))?
            .azimuth_time_interval_s;
        tracing::info!(
            "multilook {}×{} (range×azimuth) on merged σ⁰ ({} × {} → {} × {}) …",
            opts.multilook_range,
            opts.multilook_azimuth,
            merged.samples,
            merged.lines,
            merged.samples / opts.multilook_range,
            merged.lines / opts.multilook_azimuth,
        );
        let t_ml = Instant::now();
        merged = apply_multilook(merged, opts.multilook_range, opts.multilook_azimuth, ati)
            .with_context(|| "multilook")?;
        report_timing("multilook", t_ml);
    }

    // ── GRD mode: ground-range output (no DEM, no geocoding) ─────────────────
    if opts.mode == OutputMode::Grd {
        run_grd_from_prepared(
            &opts.safe,
            opts.orbit.as_deref(),
            &opts.output,
            opts.target_spacing_m_for_grd(),
            opts.cog,
            opts.no_provenance,
            &opts.speckle,
            opts.speckle_window,
            opts.enl,
            opts.frost_damping,
            opts.threads,
            &scene,
            orbit_source,
            &merged,
            &grids,
            pol,
        )?;
        tracing::info!("done.");
        report_timing("total", t_total);
        return Ok(());
    }

    // ── Pre-TC speckle (slant-range, before terrain correction) ──────────────
    let speckle_choice = resolve_speckle(
        &opts.speckle,
        opts.speckle_window,
        opts.enl,
        opts.frost_damping,
    )?;
    if opts.speckle_order == SpeckleOrder::BeforeTc {
        let t_speckle = Instant::now();
        apply_speckle_step(&mut merged.data, merged.samples, merged.lines, speckle_choice)?;
        report_timing("speckle_before_tc", t_speckle);
    }

    tracing::info!("loading DEM mosaic from {} …", opts.dem.display());
    let t_dem = Instant::now();
    let dem = DemMosaic::load_directory(&opts.dem)
        .with_context(|| format!("loading DEM from: {}", opts.dem.display()))?;
    report_timing("dem_load", t_dem);

    let bb = scene.bounding_box;
    tracing::info!(
        "scene footprint lat=[{:.4}, {:.4}], lon=[{:.4}, {:.4}]; \
         DEM mosaic has {} tile(s)",
        bb.min_lat_deg, bb.max_lat_deg, bb.min_lon_deg, bb.max_lon_deg,
        dem.tile_count(),
    );
    dem.covers_bbox(
        bb.min_lat_deg,
        bb.max_lat_deg,
        bb.min_lon_deg,
        bb.max_lon_deg,
        0.05,
    )
    .with_context(|| {
        format!(
            "DEM mosaic at {} does not cover the scene footprint",
            opts.dem.display()
        )
    })?;

    let geoid = resolve_geoid(&opts.geoid)?;
    let bb = &scene.bounding_box;
    let scene_centre = (
        0.5 * (bb.min_lon_deg + bb.max_lon_deg),
        0.5 * (bb.min_lat_deg + bb.max_lat_deg),
    );
    let output_crs = resolve_crs(&opts.crs, Some(scene_centre))?;
    let spacing = if output_crs.is_metric() {
        opts.pixel_spacing_m
    } else {
        opts.pixel_spacing_deg
    };
    tracing::info!(
        "terrain correction (CRS=EPSG:{}, spacing={} {}) …",
        output_crs.epsg(),
        spacing,
        if output_crs.is_metric() { "m" } else { "deg" },
    );
    let mut tc_config = TerrainCorrectionConfig::new(geoid);
    tc_config.crs = output_crs;
    tc_config.pixel_spacing_deg = spacing;
    tc_config.flatten = !opts.no_flatten;
    // NRB mode forces LIA and mask computation.
    let effective_write_lia = opts.write_lia || opts.mode == OutputMode::Nrb;
    let effective_write_mask = opts.write_mask || opts.mode == OutputMode::Nrb;
    let _effective_cog = opts.cog || opts.mode == OutputMode::Nrb;
    tc_config.compute_lia = effective_write_mask || effective_write_lia;
    let t_tc = Instant::now();
    let geocoded = terrain_correction(&merged, &scene, &dem, &grids, &tc_config)
        .with_context(|| "terrain correction")?;
    report_timing("terrain_correction", t_tc);

    // Newton zero-Doppler iteration histogram.  Slot `i` (1≤i≤size−1) =
    // pixels that converged in exactly `i` iterations; final slot is the
    // saturating-clamped overflow bin.  This is the data input for any
    // future tightening of `max_newton_iterations`.
    {
        let hist = &geocoded.newton_iter_histogram;
        let total: u64 = hist.iter().sum();
        if total > 0 {
            let mut parts = Vec::with_capacity(hist.len());
            for (i, n) in hist.iter().enumerate().skip(1) {
                if *n > 0 {
                    parts.push(format!("{i}={n}"));
                }
            }
            tracing::info!(
                "[timing] newton_iter_hist (total={}): {}",
                total,
                parts.join(" "),
            );
        }
    }

    tracing::info!(
        "geocoded {}×{} pixels — valid={}, dem_missing={}, \
         not_converged={}, flat_masked={}",
        geocoded.cols,
        geocoded.rows,
        geocoded.valid_pixel_count,
        geocoded.dem_missing_count,
        geocoded.non_converged_count,
        geocoded.flat_masked_count,
    );

    let lia_sidecar_path: Option<String> = if effective_write_lia {
        let lia_path = sidecar_path(&opts.output, ".lia.tif")?;
        tracing::info!("writing LIA sidecar → {}", lia_path);
        crate::export::write_geotiff_with_crs(
            &lia_path,
            &geocoded.cos_lia,
            geocoded.cols,
            geocoded.rows,
            geocoded.geotransform,
            &geocoded.crs,
        )
        .with_context(|| format!("writing LIA GeoTIFF: {lia_path}"))?;
        Some(lia_path)
    } else {
        None
    };
    let mask_sidecar_path: Option<String> = if effective_write_mask {
        let mask_path = sidecar_path(&opts.output, ".mask.tif")?;
        tracing::info!("writing quality mask sidecar → {}", mask_path);
        crate::export::write_geotiff_u8_with_crs(
            &mask_path,
            &geocoded.mask,
            geocoded.cols,
            geocoded.rows,
            geocoded.geotransform,
            &geocoded.crs,
        )
        .with_context(|| format!("writing mask GeoTIFF: {mask_path}"))?;
        Some(mask_path)
    } else {
        None
    };

    let mut geocoded = geocoded;

    // ── Post-TC speckle (map geometry, after terrain correction, default) ─────
    if opts.speckle_order == SpeckleOrder::AfterTc {
        let t_speckle = Instant::now();
        apply_speckle_step(&mut geocoded.data, geocoded.cols, geocoded.rows, speckle_choice)?;
        report_timing("speckle", t_speckle);
    }

    let noise_floor_linear = if opts.noise_floor_db <= 0.0 {
        0.0_f32
    } else {
        10.0_f32.powf(opts.noise_floor_db / 10.0)
    };
    let t_db = Instant::now();
    let (converted, masked) = to_db_inplace(&mut geocoded.data, noise_floor_linear)
        .with_context(|| "dB conversion")?;
    report_timing("db_conversion", t_db);
    tracing::info!(
        "dB conversion — converted={}, noise_masked={}",
        converted, masked
    );

    let (db_min, db_max, db_mean, db_median) = compute_db_stats(&geocoded.data);
    let total_pixel_count = geocoded.rows * geocoded.cols;
    let quality_flags = collect_quality_flags(
        orbit_source,
        merged.cal_lut_extrapolation_gap_px,
        merged.noise_lut_extrapolation_gap_px,
        total_pixel_count,
        Some(geocoded.dem_missing_count),
        Some(geocoded.non_converged_count),
    );
    if !quality_flags.is_empty() {
        tracing::warn!("{} quality flag(s) raised:", quality_flags.len());
        for f in &quality_flags {
            tracing::warn!("  [{:?}] {:?}: {}", f.severity, f.flag_type, f.message);
        }
    }

    let output_str = opts
        .output
        .to_str()
        .ok_or_else(|| anyhow!("output path contains non-UTF-8 characters"))?;
    tracing::info!("writing GeoTIFF → {}", output_str);
    let t_export = Instant::now();
    write_geotiff_with_crs(
        output_str,
        &geocoded.data,
        geocoded.cols,
        geocoded.rows,
        geocoded.geotransform,
        &geocoded.crs,
    )
    .with_context(|| format!("writing GeoTIFF: {}", output_str))?;
    report_timing("export_geotiff", t_export);

    if opts.cog {
        tracing::info!("converting to COG via gdal_translate ...");
        crate::cog::convert_to_cog(&opts.output)
            .with_context(|| format!("COG conversion of {output_str}"))?;
        tracing::info!("COG conversion done.");
    }

    if !opts.no_provenance {
        let prov = build_provenance(
            opts,
            &scene,
            orbit_source,
            dem.tile_count(),
            &merged,
            &geocoded,
            output_str,
            lia_sidecar_path.as_deref(),
            mask_sidecar_path.as_deref(),
            converted,
            masked,
            (db_min, db_max, db_mean, db_median),
            quality_flags,
            pol,
            opts.multilook_range,
            opts.multilook_azimuth,
            opts.speckle_order,
            opts.mode,
        )?;
        let prov_path_str = sidecar_path(&opts.output, ".provenance.json")?;
        tracing::info!("writing provenance sidecar → {}", prov_path_str);
        let prov_path = std::path::PathBuf::from(&prov_path_str);
        prov.write_json(&prov_path)
            .with_context(|| format!("writing provenance JSON: {prov_path_str}"))?;

        let stac_path_str = sidecar_path(&opts.output, ".stac.json")?;
        tracing::info!("writing STAC sidecar → {}", stac_path_str);
        crate::stac::write_stac_item(&prov, &std::path::PathBuf::from(&stac_path_str))
            .with_context(|| format!("writing STAC JSON: {stac_path_str}"))?;
    }

    tracing::info!("done.");
    report_timing("total", t_total);
    Ok(())
}

#[allow(clippy::too_many_arguments)]
fn build_provenance(
    opts: &ProcessOptions,
    scene: &crate::types::SceneMetadata,
    orbit_source: crate::provenance::OrbitSource,
    dem_tile_count: usize,
    merged: &crate::merge_subswaths::MergedSigma0,
    geocoded: &crate::terrain_correction::GeocodedImage,
    output_path: &str,
    lia_sidecar_path: Option<&str>,
    mask_sidecar_path: Option<&str>,
    db_converted_count: usize,
    noise_masked_count: usize,
    db_stats: (f32, f32, f32, f32),
    quality_flags: Vec<crate::provenance::QualityFlagEntry>,
    pol: Polarization,
    multilook_range: usize,
    multilook_azimuth: usize,
    speckle_order: SpeckleOrder,
    output_mode: OutputMode,
) -> Result<crate::provenance::Provenance> {
    use crate::provenance::{
        BoundingBoxJson, DemInfo, GeoidInfo, InputInfo, OrbitInfo, OutputInfo, ProcessingInfo,
        Provenance, StatsInfo, WarningsInfo, SCHEMA_VERSION,
    };

    let bb = scene.bounding_box;
    let (db_min, db_max, db_mean, db_median) = db_stats;

    let (speckle_filter_field, speckle_window_field, speckle_enl_field, speckle_damping_field) =
        speckle_provenance_fields(resolve_speckle(
            &opts.speckle,
            opts.speckle_window,
            opts.enl,
            opts.frost_damping,
        )?);

    let safe_path_str = opts
        .safe
        .to_str()
        .ok_or_else(|| anyhow!("--safe path contains non-UTF-8 characters"))?
        .to_owned();
    let dem_dir_str = opts
        .dem
        .to_str()
        .ok_or_else(|| anyhow!("--dem path contains non-UTF-8 characters"))?
        .to_owned();
    let orbit_file_path = match orbit_source {
        crate::provenance::OrbitSource::Poeorb => Some(
            opts.orbit
                .as_ref()
                .expect("orbit_source=Poeorb implies opts.orbit is Some")
                .to_str()
                .ok_or_else(|| anyhow!("--orbit path contains non-UTF-8 characters"))?
                .to_owned(),
        ),
        crate::provenance::OrbitSource::Annotation => None,
    };

    Ok(Provenance {
        schema_version: SCHEMA_VERSION,
        generated_utc: chrono::Utc::now().to_rfc3339(),
        sardine: Provenance::sardine_info(),
        input: InputInfo {
            safe_path: safe_path_str,
            product_id: scene.product_id.clone(),
            mission: scene.mission.to_string(),
            acquisition_mode: scene.acquisition_mode.to_string(),
            polarization: pol.to_string(),
            scene_start_utc: scene.start_time.to_rfc3339(),
            scene_stop_utc: scene.stop_time.to_rfc3339(),
            scene_bbox_deg: BoundingBoxJson {
                min_lat: bb.min_lat_deg,
                max_lat: bb.max_lat_deg,
                min_lon: bb.min_lon_deg,
                max_lon: bb.max_lon_deg,
            },
        },
        orbit: OrbitInfo {
            source: orbit_source,
            file_path: orbit_file_path,
            reference_epoch_utc: scene.orbit.reference_epoch.to_rfc3339(),
            state_vector_count: scene.orbit.state_vectors.len(),
        },
        dem: DemInfo {
            directory: dem_dir_str,
            tile_count: dem_tile_count,
        },
        geoid: GeoidInfo {
            spec: opts.geoid.clone(),
        },
        processing: ProcessingInfo {
            mode: Some(match output_mode {
                OutputMode::Rtc => "rtc",
                OutputMode::Nrb => "nrb",
                OutputMode::Grd => "grd",
            }.to_owned()),
            pixel_spacing_deg: Some(if geocoded.crs.is_metric() {
                opts.pixel_spacing_m
            } else {
                opts.pixel_spacing_deg
            }),
            target_spacing_m: None,
            flatten: Some(!opts.no_flatten),
            compute_lia: Some(opts.write_mask || opts.write_lia),
            noise_floor_db: opts.noise_floor_db,
            threads: opts.threads,
            speckle_filter: speckle_filter_field,
            speckle_window: speckle_window_field,
            speckle_enl: speckle_enl_field,
            speckle_damping: speckle_damping_field,
            speckle_order: Some(match speckle_order {
                SpeckleOrder::BeforeTc => "before_tc",
                SpeckleOrder::AfterTc => "after_tc",
            }.to_owned()),
        },
        output: OutputInfo {
            raster_path: output_path.to_owned(),
            lia_sidecar_path: lia_sidecar_path.map(str::to_owned),
            mask_sidecar_path: mask_sidecar_path.map(str::to_owned),
            cols: geocoded.cols,
            rows: geocoded.rows,
            geotransform: Some(geocoded.geotransform),
            crs_epsg: Some(geocoded.crs.epsg()),
            range_pixel_spacing_m: None,
            azimuth_pixel_spacing_m: None,
            range_looks: if multilook_range > 1 { Some(multilook_range) } else { None },
            azimuth_looks: if multilook_azimuth > 1 { Some(multilook_azimuth) } else { None },
            units: "dB".to_owned(),
            nodata: "NaN".to_owned(),
        },
        stats: StatsInfo {
            total_pixel_count: geocoded.rows * geocoded.cols,
            valid_pixel_count: geocoded.valid_pixel_count,
            dem_missing_count: Some(geocoded.dem_missing_count),
            non_converged_count: Some(geocoded.non_converged_count),
            flat_masked_count: Some(geocoded.flat_masked_count),
            noise_masked_count: Some(noise_masked_count),
            db_converted_count: Some(db_converted_count),
            db_min: db_min.is_finite().then_some(db_min),
            db_max: db_max.is_finite().then_some(db_max),
            db_mean: db_mean.is_finite().then_some(db_mean),
            db_median: db_median.is_finite().then_some(db_median),
        },
        warnings: WarningsInfo {
            cal_lut_extrapolation_gap_px: merged.cal_lut_extrapolation_gap_px,
            noise_lut_extrapolation_gap_px: merged.noise_lut_extrapolation_gap_px,
        },
        quality_flags,
    })
}

// ─────────────────────────────────────────────────────────────────────────────
// run_grd
// ─────────────────────────────────────────────────────────────────────────────

/// Run the ground-range pipeline (deburst → calibration → merge →
/// ground-range projection → optional speckle → TIFF, no CRS).
///
/// This is the library equivalent of the `sardine grd` CLI subcommand.
pub fn run_grd(opts: &GrdOptions) -> Result<()> {
    let t_total = Instant::now();
    let prepared = prepare_merged_scene(
        &opts.safe,
        opts.orbit.as_deref(),
        &opts.polarization,
        opts.threads,
        &opts.iw_selection,
    )?;
    let PreparedScene {
        scene,
        orbit_source,
        merged,
        grids,
        pol,
    } = prepared;

    run_grd_from_prepared(
        &opts.safe,
        opts.orbit.as_deref(),
        &opts.output,
        opts.target_spacing_m,
        opts.cog,
        opts.no_provenance,
        &opts.speckle,
        opts.speckle_window,
        opts.enl,
        opts.frost_damping,
        opts.threads,
        &scene,
        orbit_source,
        &merged,
        &grids,
        pol,
    )?;
    tracing::info!("done.");
    report_timing("total", t_total);
    Ok(())
}

// ─────────────────────────────────────────────────────────────────────────────
// Shared GRD output kernel
// ─────────────────────────────────────────────────────────────────────────────

/// Project `merged` to ground range, apply optional speckle, write TIFF, and
/// optionally write provenance + STAC sidecars.
///
/// Used by both `run_grd` and the `--mode grd` branch of `run_process`.
#[allow(clippy::too_many_arguments)]
fn run_grd_from_prepared(
    safe_path: &std::path::Path,
    orbit_path: Option<&std::path::Path>,
    output_path: &std::path::Path,
    target_spacing_m: f64,
    cog: bool,
    no_provenance: bool,
    speckle: &str,
    speckle_window: usize,
    enl: f32,
    frost_damping: f32,
    threads: usize,
    scene: &crate::types::SceneMetadata,
    orbit_source: crate::provenance::OrbitSource,
    merged: &crate::merge_subswaths::MergedSigma0,
    grids: &[(
        crate::types::SubSwathId,
        Vec<crate::types::GeolocationGridPoint>,
    )],
    pol: Polarization,
) -> Result<()> {
    use crate::export::write_tiff_no_crs;
    use crate::ground_range::to_ground_range;

    let azimuth_pixel_spacing_m = scene
        .sub_swaths
        .first()
        .ok_or_else(|| anyhow!("scene has no subswath metadata"))?
        .azimuth_pixel_spacing_m;

    tracing::info!(
        "projecting to ground range (target spacing = {:.3} m, az px = {:.3} m) …",
        target_spacing_m, azimuth_pixel_spacing_m,
    );
    let t_grd = Instant::now();
    let grd = to_ground_range(merged, grids, target_spacing_m, azimuth_pixel_spacing_m)
        .with_context(|| "ground-range projection")?;
    report_timing("ground_range_projection", t_grd);

    tracing::info!(
        "ground-range image: {} lines × {} samples (range looks = {}, azimuth looks = {})",
        grd.lines, grd.samples, grd.range_looks, grd.azimuth_looks,
    );

    let mut grd = grd;
    let speckle_choice = resolve_speckle(speckle, speckle_window, enl, frost_damping)?;
    let t_speckle = Instant::now();
    apply_speckle_step(&mut grd.data, grd.samples, grd.lines, speckle_choice)?;
    report_timing("speckle", t_speckle);

    let output_str = output_path
        .to_str()
        .ok_or_else(|| anyhow!("output path contains non-UTF-8 characters"))?;
    tracing::info!("writing ground-range TIFF (no CRS) → {}", output_str);
    let t_export = Instant::now();
    write_tiff_no_crs(output_str, &grd.data, grd.samples, grd.lines)
        .with_context(|| format!("writing ground-range TIFF: {}", output_str))?;
    report_timing("export_tiff", t_export);

    if cog {
        tracing::info!("converting to COG via gdal_translate ...");
        crate::cog::convert_to_cog(output_path)
            .with_context(|| format!("COG conversion of {output_str}"))?;
        tracing::info!("COG conversion done.");
    }

    if !no_provenance {
        let quality_flags = collect_quality_flags(
            orbit_source,
            merged.cal_lut_extrapolation_gap_px,
            merged.noise_lut_extrapolation_gap_px,
            grd.lines * grd.samples,
            None,
            None,
        );
        if !quality_flags.is_empty() {
            tracing::warn!("{} quality flag(s) raised:", quality_flags.len());
            for f in &quality_flags {
                tracing::warn!("  [{:?}] {:?}: {}", f.severity, f.flag_type, f.message);
            }
        }
        let safe_str = safe_path
            .to_str()
            .ok_or_else(|| anyhow!("--safe path contains non-UTF-8 characters"))?;
        let orbit_str = match orbit_source {
            crate::provenance::OrbitSource::Poeorb => Some(
                orbit_path
                    .expect("orbit_source=Poeorb implies orbit_path is Some") // SAFETY-OK: Poeorb variant is only set when opts.orbit is Some (see prepare_merged_scene)
                    .to_str()
                    .ok_or_else(|| anyhow!("--orbit path contains non-UTF-8 characters"))?
                    .to_owned(),
            ),
            crate::provenance::OrbitSource::Annotation => None,
        };
        let prov = build_grd_provenance(
            safe_str,
            orbit_str.as_deref(),
            scene,
            orbit_source,
            merged,
            &grd,
            output_str,
            pol,
            quality_flags,
            target_spacing_m,
            speckle,
            speckle_window,
            enl,
            frost_damping,
            threads,
        )?;
        let prov_path_str = sidecar_path(output_path, ".provenance.json")?;
        tracing::info!("writing provenance sidecar → {}", prov_path_str);
        prov.write_json(&std::path::PathBuf::from(&prov_path_str))
            .with_context(|| format!("writing provenance JSON: {prov_path_str}"))?;
        let stac_path_str = sidecar_path(output_path, ".stac.json")?;
        tracing::info!("writing STAC sidecar → {}", stac_path_str);
        crate::stac::write_stac_item(&prov, &std::path::PathBuf::from(&stac_path_str))
            .with_context(|| format!("writing STAC JSON: {stac_path_str}"))?;
    }

    Ok(())
}

#[allow(clippy::too_many_arguments)]
fn build_grd_provenance(
    safe_path_str: &str,
    orbit_file_path: Option<&str>,
    scene: &crate::types::SceneMetadata,
    orbit_source: crate::provenance::OrbitSource,
    merged: &crate::merge_subswaths::MergedSigma0,
    grd: &crate::ground_range::GroundRangeImage,
    output_path: &str,
    pol: Polarization,
    quality_flags: Vec<crate::provenance::QualityFlagEntry>,
    target_spacing_m: f64,
    speckle: &str,
    speckle_window: usize,
    enl: f32,
    frost_damping: f32,
    threads: usize,
) -> Result<crate::provenance::Provenance> {
    use crate::provenance::{
        BoundingBoxJson, DemInfo, GeoidInfo, InputInfo, OrbitInfo, OutputInfo, ProcessingInfo,
        Provenance, StatsInfo, WarningsInfo, SCHEMA_VERSION,
    };

    let bb = scene.bounding_box;

    let (speckle_filter_field, speckle_window_field, speckle_enl_field, speckle_damping_field) =
        speckle_provenance_fields(resolve_speckle(speckle, speckle_window, enl, frost_damping)?);

    let total = grd.lines.checked_mul(grd.samples).ok_or_else(|| {
        anyhow!(
            "ground-range image dimensions overflow usize: {} × {}",
            grd.lines,
            grd.samples
        )
    })?;
    let valid = grd.data.iter().filter(|x| x.is_finite()).count();

    Ok(Provenance {
        schema_version: SCHEMA_VERSION,
        generated_utc: chrono::Utc::now().to_rfc3339(),
        sardine: Provenance::sardine_info(),
        input: InputInfo {
            safe_path: safe_path_str.to_owned(),
            product_id: scene.product_id.clone(),
            mission: scene.mission.to_string(),
            acquisition_mode: scene.acquisition_mode.to_string(),
            polarization: pol.to_string(),
            scene_start_utc: scene.start_time.to_rfc3339(),
            scene_stop_utc: scene.stop_time.to_rfc3339(),
            scene_bbox_deg: BoundingBoxJson {
                min_lat: bb.min_lat_deg,
                max_lat: bb.max_lat_deg,
                min_lon: bb.min_lon_deg,
                max_lon: bb.max_lon_deg,
            },
        },
        orbit: OrbitInfo {
            source: orbit_source,
            file_path: orbit_file_path.map(str::to_owned),
            reference_epoch_utc: scene.orbit.reference_epoch.to_rfc3339(),
            state_vector_count: scene.orbit.state_vectors.len(),
        },
        dem: DemInfo {
            directory: String::new(),
            tile_count: 0,
        },
        geoid: GeoidInfo {
            spec: String::new(),
        },
        processing: ProcessingInfo {
            mode: Some("grd".to_owned()),
            pixel_spacing_deg: None,
            target_spacing_m: Some(target_spacing_m),
            flatten: None,
            compute_lia: None,
            noise_floor_db: 0.0,
            threads,
            speckle_filter: speckle_filter_field,
            speckle_window: speckle_window_field,
            speckle_enl: speckle_enl_field,
            speckle_damping: speckle_damping_field,
            speckle_order: None,
        },
        output: OutputInfo {
            raster_path: output_path.to_owned(),
            lia_sidecar_path: None,
            mask_sidecar_path: None,
            cols: grd.samples,
            rows: grd.lines,
            geotransform: None,
            crs_epsg: None,
            range_pixel_spacing_m: Some(grd.range_pixel_spacing_m),
            azimuth_pixel_spacing_m: Some(grd.azimuth_pixel_spacing_m),
            range_looks: Some(grd.range_looks),
            azimuth_looks: Some(grd.azimuth_looks),
            units: "linear".to_owned(),
            nodata: "NaN".to_owned(),
        },
        stats: StatsInfo {
            total_pixel_count: total,
            valid_pixel_count: valid,
            dem_missing_count: None,
            non_converged_count: None,
            flat_masked_count: None,
            noise_masked_count: None,
            db_converted_count: None,
            db_min: None,
            db_max: None,
            db_mean: None,
            db_median: None,
        },
        warnings: WarningsInfo {
            cal_lut_extrapolation_gap_px: merged.cal_lut_extrapolation_gap_px,
            noise_lut_extrapolation_gap_px: merged.noise_lut_extrapolation_gap_px,
        },
        quality_flags,
    })
}

// ─────────────────────────────────────────────────────────────────────────────
// Dual-polarization driver
// ─────────────────────────────────────────────────────────────────────────────
//
// `parse_polarizations`, `derive_pol_output_path`, `run_process_multi`, and
// `run_grd_multi` together implement single-invocation dual-pol processing.
//
// Design (deliberately minimal — see AGENTS.md "Do not make broad refactors"):
//
// * The single-pol entries (`run_process`, `run_grd`) are unchanged; they
//   remain the strict, byte-stable single-output contract.
// * The `_multi` entries split the polarization spec, then call the
//   single-pol entry once per pol with a derived output path
//   (`out.tif` + `VV` → `out_VV.tif`).  Geometry (orbit interp, DEM lookup,
//   Newton solver, Small 2011 flattening factor) is recomputed per pol.
//   That is wasteful (~2× wall time for dual-pol) but correct, and avoids
//   a refactor that would split `terrain_correction` into a geometry-pass
//   plus a resample-pass.
// * No silent fallbacks: an unsupported pol token is a typed error.

/// Split a polarization spec into a deduplicated, validated list of
/// polarization tokens.
///
/// Accepted spellings (case-insensitive):
///
/// * Single pol: `"VV"`, `"VH"`.
/// * Dual pol  : `"VV+VH"`, `"VV,VH"`, `"vh+vv"`, `"dual"`, `"both"`.
///
/// `dual` / `both` always expand to `["VV", "VH"]`.  Any other spelling
/// preserves the input order so the caller can choose the output order.
/// Duplicate tokens (e.g. `"VV+VV"`) are silently deduplicated.
///
/// # Errors
///
/// Returns an error for empty input, unknown tokens, or more than two
/// distinct pols (S-1 IW SLC carries at most two).
pub fn parse_polarizations(spec: &str) -> Result<Vec<String>> {
    let upper = spec.to_uppercase();
    let raw_tokens: Vec<String> = match upper.as_str() {
        "DUAL" | "BOTH" => vec!["VV".to_owned(), "VH".to_owned()],
        _ => upper
            .split(|c: char| c == ',' || c == '+')
            .map(|s| s.trim().to_owned())
            .collect(),
    };

    let mut out: Vec<String> = Vec::new();
    for tok in raw_tokens {
        if tok.is_empty() {
            bail!("polarization spec contains an empty token: {:?}", spec);
        }
        if tok != "VV" && tok != "VH" {
            bail!(
                "Unsupported polarization token: {:?} in spec {:?}. \
                 Use VV, VH, VV+VH, dual, or both.",
                tok,
                spec
            );
        }
        if !out.iter().any(|x| *x == tok) {
            out.push(tok);
        }
    }
    if out.is_empty() {
        bail!("polarization spec is empty: {:?}", spec);
    }
    if out.len() > 2 {
        bail!(
            "Sentinel-1 IW SLC carries at most 2 polarizations; spec {:?} \
             yielded {} distinct tokens",
            spec,
            out.len()
        );
    }
    Ok(out)
}

/// Insert a polarization tag before the extension of an output path.
///
/// `out.tif` + `"VV"` → `out_VV.tif`, `out` + `"VH"` → `out_VH`,
/// `dir/sub/x.tiff` + `"VV"` → `dir/sub/x_VV.tiff`.
///
/// # Errors
///
/// Errors if the path has no filename or contains non-UTF-8 components.
pub fn derive_pol_output_path(output: &Path, pol: &str) -> Result<PathBuf> {
    let stem = output
        .file_stem()
        .ok_or_else(|| anyhow!("output path has no filename: {}", output.display()))?
        .to_str()
        .ok_or_else(|| {
            anyhow!(
                "output path filename is not valid UTF-8: {}",
                output.display()
            )
        })?;
    let new_name = match output.extension().and_then(|e| e.to_str()) {
        Some(ext) => format!("{stem}_{pol}.{ext}"),
        None => format!("{stem}_{pol}"),
    };
    Ok(match output.parent() {
        Some(p) if !p.as_os_str().is_empty() => p.join(new_name),
        _ => PathBuf::from(new_name),
    })
}

/// Parse an IW sub-swath selection string into an [`IwSelection`].
///
/// `iw_spec` is a comma-separated list of sub-swath names or empty for
/// "all three":
///
/// ```text
/// ""          → all three sub-swaths, all bursts (default)
/// "IW1"       → IW1 only, all bursts
/// "IW2,IW3"   → IW2 + IW3, all bursts
/// ```
///
/// `burst_range` is an optional `"START-END"` string (0-based, inclusive)
/// that limits the burst range within each selected sub-swath.
///
/// # Errors
///
/// Returns an error for unknown sub-swath names or malformed burst-range strings.
pub fn parse_iw_selection(iw_spec: &str, burst_range: Option<&str>) -> Result<IwSelection> {
    use crate::types::SubSwathId;

    let subswaths: Vec<SubSwathId> = if iw_spec.trim().is_empty() {
        vec![]
    } else {
        iw_spec
            .split(',')
            .map(|s| match s.trim().to_uppercase().as_str() {
                "IW1" => Ok(SubSwathId::IW1),
                "IW2" => Ok(SubSwathId::IW2),
                "IW3" => Ok(SubSwathId::IW3),
                other => Err(anyhow!(
                    "unknown sub-swath '{}'; expected IW1, IW2, or IW3",
                    other
                )),
            })
            .collect::<Result<Vec<_>>>()?
    };

    let burst_range = match burst_range {
        None => None,
        Some(s) => {
            let parts: Vec<&str> = s.splitn(2, '-').collect();
            if parts.len() != 2 {
                bail!(
                    "burst_range must be START-END (e.g. '0-5'), got '{}'",
                    s
                );
            }
            let start: usize = parts[0]
                .parse()
                .with_context(|| {
                    format!("parsing burst range start '{}': not a number", parts[0])
                })?;
            let end: usize = parts[1]
                .parse()
                .with_context(|| {
                    format!("parsing burst range end '{}': not a number", parts[1])
                })?;
            if end < start {
                bail!(
                    "burst_range end ({}) must be >= start ({})",
                    end,
                    start
                );
            }
            Some((start, end))
        }
    };

    Ok(IwSelection {
        subswaths,
        burst_range,
    })
}

/// Parse an output-mode string into [`OutputMode`].
///
/// Accepted values (case-insensitive): `"rtc"` (default), `"grd"`, `"nrb"`.
///
/// # Errors
///
/// Returns an error for unrecognised values.
pub fn parse_output_mode(s: &str) -> Result<OutputMode> {
    match s.to_ascii_lowercase().as_str() {
        "rtc" => Ok(OutputMode::Rtc),
        "grd" => Ok(OutputMode::Grd),
        "nrb" => Ok(OutputMode::Nrb),
        other => bail!(
            "unknown mode '{}'; expected 'rtc', 'grd', or 'nrb'",
            other
        ),
    }
}

/// Parse a speckle-order string into [`SpeckleOrder`].
///
/// Accepted values (case-insensitive): `"after"` (default), `"before"`.
///
/// # Errors
///
/// Returns an error for unrecognised values.
pub fn parse_speckle_order(s: &str) -> Result<SpeckleOrder> {
    match s.to_ascii_lowercase().as_str() {
        "after" => Ok(SpeckleOrder::AfterTc),
        "before" => Ok(SpeckleOrder::BeforeTc),
        other => bail!(
            "unknown speckle_order '{}'; expected 'before' or 'after'",
            other
        ),
    }
}

/// Multi-polarization wrapper around [`run_process`].
///
/// Accepts the same `ProcessOptions`, but reads `opts.polarization` as a
/// spec parsed by [`parse_polarizations`].  Single-pol input is forwarded
/// to `run_process` unchanged (output path preserved exactly).  Dual-pol
/// input runs the full pipeline once per pol and writes to a derived path
/// per pol via [`derive_pol_output_path`].
///
/// Each per-pol invocation writes its own `.provenance.json`,
/// `.lia.tif`, and `.mask.tif` sidecars next to its own raster.
pub fn run_process_multi(opts: &ProcessOptions) -> Result<()> {
    let pols = parse_polarizations(&opts.polarization)?;
    if pols.len() == 1 {
        // Single-pol fast path: keep the user's output path exactly as given.
        let mut sub = opts.clone();
        sub.polarization = pols.into_iter().next().expect("len()==1");
        return run_process(&sub);
    }
    tracing::info!(
        "dual-polarization run — pols = {}; geometry will be recomputed per pol",
        pols.join("+")
    );
    for pol in &pols {
        let mut sub = opts.clone();
        sub.polarization = pol.clone();
        sub.output = derive_pol_output_path(&opts.output, pol)?;
        tracing::info!(
            "--- polarization {} → {} ---",
            pol,
            sub.output.display()
        );
        run_process(&sub)
            .with_context(|| format!("processing polarization {pol}"))?;
    }
    Ok(())
}

/// Multi-polarization wrapper around [`run_grd`].  Same contract as
/// [`run_process_multi`] but for the ground-range pipeline.
pub fn run_grd_multi(opts: &GrdOptions) -> Result<()> {
    let pols = parse_polarizations(&opts.polarization)?;
    if pols.len() == 1 {
        let mut sub = opts.clone();
        sub.polarization = pols.into_iter().next().expect("len()==1");
        return run_grd(&sub);
    }
    tracing::info!(
        "dual-polarization GRD run — pols = {}; per-pol pipeline",
        pols.join("+")
    );
    for pol in &pols {
        let mut sub = opts.clone();
        sub.polarization = pol.clone();
        sub.output = derive_pol_output_path(&opts.output, pol)?;
        tracing::info!(
            "--- polarization {} → {} ---",
            pol,
            sub.output.display()
        );
        run_grd(&sub)
            .with_context(|| format!("processing polarization {pol}"))?;
    }
    Ok(())
}

#[cfg(test)]
mod multipol_tests {
    use super::*;
    use crate::types::SubSwathId;

    #[test]
    fn parse_polarizations_single_vv() {
        assert_eq!(parse_polarizations("VV").unwrap(), vec!["VV"]);
        assert_eq!(parse_polarizations("vv").unwrap(), vec!["VV"]);
    }

    #[test]
    fn parse_polarizations_single_vh() {
        assert_eq!(parse_polarizations("VH").unwrap(), vec!["VH"]);
    }

    #[test]
    fn parse_polarizations_dual_plus_form() {
        assert_eq!(parse_polarizations("VV+VH").unwrap(), vec!["VV", "VH"]);
    }

    #[test]
    fn parse_polarizations_dual_comma_form() {
        assert_eq!(parse_polarizations("vv,vh").unwrap(), vec!["VV", "VH"]);
    }

    #[test]
    fn parse_polarizations_preserves_order() {
        assert_eq!(parse_polarizations("VH+VV").unwrap(), vec!["VH", "VV"]);
    }

    #[test]
    fn parse_polarizations_dual_alias() {
        assert_eq!(parse_polarizations("dual").unwrap(), vec!["VV", "VH"]);
        assert_eq!(parse_polarizations("Both").unwrap(), vec!["VV", "VH"]);
    }

    #[test]
    fn parse_polarizations_dedupe() {
        assert_eq!(parse_polarizations("VV+VV").unwrap(), vec!["VV"]);
    }

    #[test]
    fn parse_polarizations_rejects_unknown() {
        let err = parse_polarizations("HH").unwrap_err().to_string();
        assert!(err.contains("Unsupported polarization token"), "got: {err}");
    }

    #[test]
    fn parse_polarizations_rejects_empty() {
        assert!(parse_polarizations("").is_err());
        assert!(parse_polarizations("VV+").is_err());
        assert!(parse_polarizations("+VV").is_err());
    }

    #[test]
    fn parse_polarizations_rejects_three_distinct() {
        let err = parse_polarizations("VV+VH+HH").unwrap_err().to_string();
        assert!(
            err.contains("Unsupported") || err.contains("at most 2"),
            "got: {err}"
        );
    }

    #[test]
    fn derive_pol_output_path_with_extension() {
        let p = derive_pol_output_path(Path::new("out.tif"), "VV").unwrap();
        assert_eq!(p, PathBuf::from("out_VV.tif"));
    }

    #[test]
    fn derive_pol_output_path_with_dir() {
        let p = derive_pol_output_path(Path::new("/tmp/x/out.tiff"), "VH").unwrap();
        assert_eq!(p, PathBuf::from("/tmp/x/out_VH.tiff"));
    }

    #[test]
    fn derive_pol_output_path_no_extension() {
        let p = derive_pol_output_path(Path::new("out"), "VV").unwrap();
        assert_eq!(p, PathBuf::from("out_VV"));
    }

    #[test]
    fn derive_pol_output_path_multidot_keeps_last_extension() {
        let p = derive_pol_output_path(Path::new("a.b.tif"), "VH").unwrap();
        assert_eq!(p, PathBuf::from("a.b_VH.tif"));
    }

    // --- parse_iw_selection ---

    #[test]
    fn parse_iw_selection_empty_means_all() {
        let sel = parse_iw_selection("", None).unwrap();
        assert!(sel.subswaths.is_empty());
        assert!(sel.burst_range.is_none());
    }

    #[test]
    fn parse_iw_selection_single_subswath() {
        let sel = parse_iw_selection("IW1", None).unwrap();
        assert_eq!(sel.subswaths, vec![SubSwathId::IW1]);
    }

    #[test]
    fn parse_iw_selection_multiple_subswaths() {
        let sel = parse_iw_selection("IW1,IW3", None).unwrap();
        assert_eq!(sel.subswaths, vec![SubSwathId::IW1, SubSwathId::IW3]);
    }

    #[test]
    fn parse_iw_selection_with_burst_range() {
        let sel = parse_iw_selection("IW2", Some("2-5")).unwrap();
        assert_eq!(sel.burst_range, Some((2, 5)));
    }

    #[test]
    fn parse_iw_selection_rejects_invalid_subswath() {
        let err = parse_iw_selection("IW4", None).unwrap_err().to_string();
        assert!(!err.is_empty());
    }

    #[test]
    fn parse_iw_selection_rejects_bad_burst_range_order() {
        let err = parse_iw_selection("IW1", Some("5-2")).unwrap_err().to_string();
        assert!(!err.is_empty());
    }

    // --- parse_output_mode ---

    #[test]
    fn parse_output_mode_rtc() {
        assert!(matches!(parse_output_mode("rtc").unwrap(), OutputMode::Rtc));
        assert!(matches!(parse_output_mode("RTC").unwrap(), OutputMode::Rtc));
    }

    #[test]
    fn parse_output_mode_grd() {
        assert!(matches!(parse_output_mode("grd").unwrap(), OutputMode::Grd));
    }

    #[test]
    fn parse_output_mode_nrb() {
        assert!(matches!(parse_output_mode("nrb").unwrap(), OutputMode::Nrb));
    }

    #[test]
    fn parse_output_mode_rejects_unknown() {
        let err = parse_output_mode("sigma0").unwrap_err().to_string();
        assert!(!err.is_empty());
    }

    // --- parse_speckle_order ---

    #[test]
    fn parse_speckle_order_after() {
        assert!(matches!(
            parse_speckle_order("after").unwrap(),
            SpeckleOrder::AfterTc
        ));
        assert!(matches!(
            parse_speckle_order("After").unwrap(),
            SpeckleOrder::AfterTc
        ));
    }

    #[test]
    fn parse_speckle_order_before() {
        assert!(matches!(
            parse_speckle_order("before").unwrap(),
            SpeckleOrder::BeforeTc
        ));
    }

    #[test]
    fn parse_speckle_order_rejects_unknown() {
        let err = parse_speckle_order("during").unwrap_err().to_string();
        assert!(!err.is_empty());
    }

    #[test]
    fn run_process_multi_reports_unknown_pol_before_io() {
        // No filesystem reads are expected because polarization parsing
        // fails up-front.  Paths can be bogus.
        let opts = ProcessOptions {
            safe: PathBuf::from("/no/such/safe.SAFE"),
            dem: PathBuf::from("/no/such/dem"),
            output: PathBuf::from("/tmp/no_such_out.tif"),
            orbit: None,
            polarization: "HH".to_owned(),
            no_flatten: false,
            noise_floor_db: 0.0,
            pixel_spacing_deg: 0.0001,
            pixel_spacing_m: 10.0,
            crs: "wgs84".to_owned(),
            write_mask: false,
            write_lia: false,
            no_provenance: true,
            cog: false,
            threads: 1,
            speckle: "none".to_owned(),
            speckle_window: 7,
            enl: 1.0,
            frost_damping: 1.0,
            geoid: "zero".to_owned(),
            multilook_range: 1,
            multilook_azimuth: 1,
            extra_safe_paths: vec![],
            iw_selection: crate::run::IwSelection::default(),
            mode: crate::run::OutputMode::default(),
            speckle_order: crate::run::SpeckleOrder::default(),
        };
        let err = run_process_multi(&opts).unwrap_err().to_string();
        assert!(
            err.contains("Unsupported polarization token"),
            "expected polarization parse error, got: {err}"
        );
    }
}

#[cfg(test)]
mod multilook_tests {
    use super::apply_multilook;
    use crate::merge_subswaths::MergedSigma0;

    fn dummy_merged(lines: usize, samples: usize, fill: f32) -> MergedSigma0 {
        MergedSigma0 {
            data: vec![fill; lines * samples],
            lines,
            samples,
            near_slant_range_time_s: 0.005_300,
            range_pixel_spacing_s: 2.329562_f64 / 299_792_458.0_f64 * 2.0,
            range_pixel_spacing_m: 2.329562,
            cal_lut_extrapolation_gap_px: 0,
            noise_lut_extrapolation_gap_px: 0,
            nesz: vec![0.01_f32; lines * samples],
            azimuth_start_time: chrono::DateTime::parse_from_rfc3339("2019-01-23T05:33:48Z")
                .expect("static literal") // SAFETY-OK: constant literal cannot be invalid
                .with_timezone(&chrono::Utc),
        }
    }

    /// R=1, A=1 must be a bit-identical pass-through (no allocation).
    #[test]
    fn multilook_passthrough_r1_a1() {
        let m = dummy_merged(4, 8, 2.0);
        let expected_t0 = m.near_slant_range_time_s;
        let out = apply_multilook(m, 1, 1, 0.0).unwrap(); // SAFETY-OK: r=1,a=1 never uses ati
        assert_eq!(out.lines, 4);
        assert_eq!(out.samples, 8);
        assert_eq!(out.data[0], 2.0);
        assert_eq!(out.near_slant_range_time_s, expected_t0);
    }

    /// R=2, A=1: verify averaging of uniform-value pairs.
    #[test]
    fn multilook_range_2_averages_uniform() {
        // Input row: [a, a, b, b]  → [a, b]
        let mut m = dummy_merged(1, 4, 0.0);
        m.data = vec![1.0, 1.0, 5.0, 5.0];
        let out = apply_multilook(m, 2, 1, 0.0).unwrap(); // SAFETY-OK: a=1, ati unused
        assert_eq!(out.lines, 1);
        assert_eq!(out.samples, 2);
        assert!((out.data[0] - 1.0_f32).abs() < 1e-6);
        assert!((out.data[1] - 5.0_f32).abs() < 1e-6);
    }

    /// R=2, A=1: verify averaging with distinct values.
    #[test]
    fn multilook_range_2_averages_distinct() {
        // Row 0: [1, 3, 5, 7] → [2, 6]
        // Row 1: [2, 4, 6, 8] → [3, 7]
        let mut m = dummy_merged(2, 4, 0.0);
        m.data = vec![1.0, 3.0, 5.0, 7.0, 2.0, 4.0, 6.0, 8.0];
        let out = apply_multilook(m, 2, 1, 0.0).unwrap(); // SAFETY-OK: a=1, ati unused
        assert_eq!(out.lines, 2);
        assert_eq!(out.samples, 2);
        assert!((out.data[0] - 2.0_f32).abs() < 1e-6);
        assert!((out.data[1] - 6.0_f32).abs() < 1e-6);
        assert!((out.data[2] - 3.0_f32).abs() < 1e-6);
        assert!((out.data[3] - 7.0_f32).abs() < 1e-6);
    }

    /// NaN pixels (seam gaps, swath edges) must be excluded from the average.
    #[test]
    fn multilook_nan_excluded_from_average() {
        // [NaN, 4.0] → R=2 → average of valid only → 4.0
        let mut m = dummy_merged(1, 2, 0.0);
        m.data = vec![f32::NAN, 4.0];
        m.nesz = vec![f32::NAN, 0.1];
        let out = apply_multilook(m, 2, 1, 0.0).unwrap(); // SAFETY-OK: a=1, ati unused
        assert!((out.data[0] - 4.0_f32).abs() < 1e-6);
        assert!((out.nesz[0] - 0.1_f32).abs() < 1e-6);
    }

    /// A block where *all* pixels are NaN must produce a NaN output pixel.
    #[test]
    fn multilook_all_nan_block_gives_nan_output() {
        let mut m = dummy_merged(1, 2, 0.0);
        m.data = vec![f32::NAN, f32::NAN];
        let out = apply_multilook(m, 2, 1, 0.0).unwrap(); // SAFETY-OK: a=1, ati unused
        assert!(out.data[0].is_nan());
    }

    /// R=4: verify the geometry fields are updated correctly.
    #[test]
    fn multilook_range_4_updates_geometry() {
        let m = dummy_merged(2, 8, 1.0);
        let old_t0 = m.near_slant_range_time_s;
        let old_dt = m.range_pixel_spacing_s;
        let old_dm = m.range_pixel_spacing_m;
        let out = apply_multilook(m, 4, 1, 0.0).unwrap(); // SAFETY-OK: a=1, ati unused
        assert_eq!(out.samples, 2);
        assert_eq!(out.lines, 2);
        // Spacing must scale by R.
        assert!((out.range_pixel_spacing_s - 4.0 * old_dt).abs() < 1e-25);
        assert!((out.range_pixel_spacing_m - 4.0 * old_dm).abs() < 1e-12);
        // Centre of first block: old_t0 + (4−1)/2 × old_dt = old_t0 + 1.5 × old_dt
        let expected_t0 = old_t0 + 1.5 * old_dt;
        assert!(
            (out.near_slant_range_time_s - expected_t0).abs() < 1e-25,
            "expected {expected_t0}, got {}",
            out.near_slant_range_time_s
        );
    }

    /// A=2: verify azimuth start time is shifted by half a block.
    #[test]
    fn multilook_azimuth_2_updates_start_time() {
        let m = dummy_merged(4, 2, 1.0);
        let old_az = m.azimuth_start_time;
        let ati_s = 2.055e-3; // representative S-1 IW ATI
        let out = apply_multilook(m, 1, 2, ati_s).unwrap();
        assert_eq!(out.lines, 2);
        assert_eq!(out.samples, 2);
        // Expected shift: (2−1)/2 × ati = 0.5 × ati_s = 1027.5 µs
        let shift_us = (0.5 * ati_s * 1e6) as i64;
        let expected_az = old_az + chrono::Duration::microseconds(shift_us);
        assert_eq!(out.azimuth_start_time, expected_az);
    }

    /// R=0 must be rejected with an explicit error, never a panic.
    #[test]
    fn multilook_zero_range_is_error() {
        let m = dummy_merged(2, 4, 1.0);
        let err = apply_multilook(m, 0, 1, 0.0).unwrap_err().to_string();
        assert!(err.contains("multilook factors must be"), "got: {err}");
    }

    /// Factors larger than the image dimensions must be rejected.
    #[test]
    fn multilook_oversized_factor_is_error() {
        let m = dummy_merged(2, 4, 1.0);
        let err = apply_multilook(m, 8, 1, 0.0).unwrap_err().to_string();
        assert!(
            err.contains("image dimensions") || err.contains("larger than"),
            "got: {err}"
        );
    }
}
