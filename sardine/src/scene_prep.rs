//! Shared scene preparation: deburst, calibrate, merge subswaths, and
//! optional multilook.  Used by both the `process` and `grd` pipelines.

use std::path::{Path, PathBuf};
use std::time::Instant;

use anyhow::{bail, Context, Result};
use rayon::prelude::*;

use crate::pipeline_options::IwSelection;
use crate::types::Polarization;

// ─────────────────────────────────────────────────────────────────────────────
// Output type
// ─────────────────────────────────────────────────────────────────────────────

/// Output of [`prepare_merged_scene`] and [`prepare_merged_scene_assembled`].
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

// ─────────────────────────────────────────────────────────────────────────────
// Private helpers
// ─────────────────────────────────────────────────────────────────────────────

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

// ─────────────────────────────────────────────────────────────────────────────
// Orbit resolution helper
// ─────────────────────────────────────────────────────────────────────────────

/// Resolve the orbit for `scene`, applying it and returning the updated
/// `SceneMetadata` together with an [`OrbitSource`] tag for provenance.
///
/// Resolution order:
///
/// 1. **`orbit_override` is `Some(path)`** — parse and apply that file.
///    No network access.
///
/// 2. **`orbit-fetch` feature is enabled (the default)** — locate the orbit
///    cache directory (`$SARDINE_ORBIT_DIR` or `$HOME/.sardine/orbits/`),
///    then call [`crate::orbit_fetch::fetch_poeorb`].  The first successful
///    download is cached on disk so subsequent runs are instant.
///
/// 3. **`SARDINE_ALLOW_ANNOTATION_ORBIT=1` env var** — use the orbit vectors
///    already baked into the SAFE annotation XML (metre-level accuracy).
///    This is intentionally not the default: annotation orbits are not
///    suitable for precision terrain correction.
///
/// Any other combination returns an explicit error — there is no silent
/// fallback to annotation orbit.
///
/// `label` is a human-readable scene name used in error messages (e.g.
/// `"reference"` or `"secondary"`).
pub(crate) fn resolve_orbit(
    scene: crate::types::SceneMetadata,
    orbit_override: Option<&Path>,
    label: &str,
) -> Result<(crate::types::SceneMetadata, crate::provenance::OrbitSource)> {
    use crate::orbit::{apply_precise_orbit, parse_eof_file};

    // ── Path 1: explicit file ─────────────────────────────────────────────
    if let Some(orbit_path) = orbit_override {
        let orbit = parse_eof_file(orbit_path)
            .with_context(|| format!("parsing {label} orbit: {}", orbit_path.display()))?;
        let scene = apply_precise_orbit(scene, &orbit)
            .with_context(|| format!("applying {label} precise orbit"))?;
        return Ok((scene, crate::provenance::OrbitSource::Poeorb));
    }

    // ── Path 2: auto-download (orbit-fetch feature) ───────────────────────
    #[cfg(feature = "orbit-fetch")]
    {
        // Opt-out: if SARDINE_ALLOW_ANNOTATION_ORBIT=1 is set, skip the
        // download and fall through to Path 3 (annotation orbit).
        let skip_fetch = std::env::var("SARDINE_ALLOW_ANNOTATION_ORBIT")
            .map(|v| v == "1")
            .unwrap_or(false); // SAFETY-OK: missing/non-"1" → do not skip (safe direction)

        if !skip_fetch {
        let cache_dir = orbit_cache_dir()
            .with_context(|| format!("resolving {label} orbit cache directory"))?;

        tracing::info!(
            "{label}: auto-downloading POEORB for {} (sensing_start={}) …",
            scene.product_id,
            scene.start_time.format("%Y-%m-%dT%H:%M:%SZ"),
        );

        let eof_path = crate::orbit_fetch::fetch_poeorb(
            &scene.product_id,
            scene.start_time,
            &cache_dir,
        )
        .with_context(|| {
            format!(
                "auto-downloading POEORB for {label} scene '{}'. \
                 Supply --{}orbit <path> to use a local file, or set \
                 SARDINE_ALLOW_ANNOTATION_ORBIT=1 to accept annotation-orbit \
                 accuracy (metre-level).",
                scene.product_id,
                if label == "reference" || label == "secondary" {
                    format!("{label}-")
                } else {
                    String::new()
                }
            )
        })?;

        tracing::info!("{label}: POEORB cached at {}", eof_path.display());

        let orbit = parse_eof_file(&eof_path)
            .with_context(|| format!("parsing downloaded {label} orbit: {}", eof_path.display()))?;
        let scene = apply_precise_orbit(scene, &orbit)
            .with_context(|| format!("applying {label} precise orbit"))?;
        return Ok((scene, crate::provenance::OrbitSource::Poeorb));
        } // if !skip_fetch

        // skip_fetch=true: fall through to annotation orbit.
        tracing::warn!(
            "{label}: SARDINE_ALLOW_ANNOTATION_ORBIT=1 — skipping POEORB download, \
             using annotation orbit (metre-level accuracy)."
        );
        return Ok((scene, crate::provenance::OrbitSource::Annotation));
    }

    // ── Path 3: annotation orbit (opt-in via env var, no orbit-fetch) ─────
    #[cfg(not(feature = "orbit-fetch"))]
    {
        let allow = std::env::var("SARDINE_ALLOW_ANNOTATION_ORBIT")
            .map(|v| v == "1")
            .unwrap_or(false); // SAFETY-OK: env var parse; missing/invalid → deny (safe direction)
        if !allow {
            bail!(
                "No --{lbl}orbit file provided and the `orbit-fetch` feature is disabled.\n\
                 Options:\n  \
                 1. Supply a POEORB file with --{lbl}orbit <path>\n  \
                 2. Set SARDINE_ALLOW_ANNOTATION_ORBIT=1 to use the annotation orbit \
                 (metre-level accuracy — not recommended for geocoded products)",
                lbl = if label == "reference" || label == "secondary" {
                    format!("{label}-")
                } else {
                    String::new()
                }
            );
        }
        tracing::warn!(
            "{label}: using annotation orbit (metre-level accuracy). \
             Provide a POEORB file or enable the `orbit-fetch` feature for cm-level accuracy."
        );
        return Ok((scene, crate::provenance::OrbitSource::Annotation));
    }
}

/// Resolve the orbit cache directory.
///
/// Checks `$SARDINE_ORBIT_DIR` first; falls back to `$HOME/.sardine/orbits/`.
fn orbit_cache_dir() -> Result<PathBuf> {
    if let Ok(dir) = std::env::var("SARDINE_ORBIT_DIR") {
        return Ok(PathBuf::from(dir));
    }
    let home = std::env::var("HOME").with_context(|| {
        "HOME environment variable not set and SARDINE_ORBIT_DIR not set; \
         cannot locate orbit cache directory"
    })?;
    Ok(PathBuf::from(home).join(".sardine").join("orbits"))
}

/// Resolve the DEM for a scene bounding box, ensuring all required tiles are
/// present and returning the directory to pass to
/// [`crate::dem::DemMosaic::load_directory`].
///
/// Resolution order:
///
/// 1. **`dem_override` is `Some(dir)`** — use that directory directly.
///    No network access; the caller is responsible for tile coverage.
///
/// 2. **`dem-fetch` feature is enabled (the default)** — locate the DEM
///    cache directory, then call the appropriate fetcher based on `dem_source`:
///    - `"srtm1"` (default): [`crate::dem_fetch::fetch_srtm1_tiles`], caches
///      into `$SARDINE_DEM_DIR/` (or `$HOME/.sardine/dem/`).
///    - `"glo30"`: [`crate::dem_fetch::fetch_glo30_tiles`], caches into
///      `$SARDINE_DEM_DIR/glo30/` (or `$HOME/.sardine/dem/glo30/`).
///    Already-cached tiles are reused without re-downloading.
///
/// Any other combination returns an explicit error.
///
/// `bb` is the scene bounding box used to compute which tiles are required.
pub(crate) fn resolve_dem(
    dem_override: Option<&Path>,
    dem_source: &str,
    bb: &crate::types::BoundingBox,
) -> Result<PathBuf> {
    // ── Path 1: explicit directory ────────────────────────────────────────
    if let Some(dir) = dem_override {
        // For SRTM-1, verify that all required tiles are actually present
        // in the supplied directory.  This catches the common mistake of
        // pointing --dem at a directory that was prepared for a different
        // scene area (the tiles would be present but for the wrong bbox,
        // and terrain-correction would silently produce all-nodata output).
        #[cfg(feature = "dem-fetch")]
        if dem_source == "srtm1" {
            check_srtm1_tiles_present(dir, bb)
                .with_context(|| {
                    format!(
                        "DEM directory '{}' does not fully cover this scene (--dem-source srtm1). \
                         Pass --dem pointing to a directory with the correct tiles, or \
                         omit --dem to auto-download.",
                        dir.display()
                    )
                })?;
        }
        return Ok(dir.to_path_buf());
    }

    // ── Path 2: auto-download (dem-fetch feature) ─────────────────────────
    #[cfg(feature = "dem-fetch")]
    {
        let cache_dir = dem_cache_dir()
            .with_context(|| "resolving DEM cache directory")?;

        tracing::info!(
            "auto-downloading {} DEM tiles for bbox \
             [{:.3}°N, {:.3}°N] × [{:.3}°E, {:.3}°E] …",
            dem_source,
            bb.min_lat_deg, bb.max_lat_deg,
            bb.min_lon_deg, bb.max_lon_deg,
        );

        let dir = crate::dem_fetch::fetch_dem_tiles_for_source(
            dem_source,
            bb.min_lat_deg,
            bb.max_lat_deg,
            bb.min_lon_deg,
            bb.max_lon_deg,
            &cache_dir,
        )
        .with_context(|| {
            format!(
                "auto-downloading {} DEM tiles for bbox \
                 [{:.3}°N, {:.3}°N] × [{:.3}°E, {:.3}°E]. \
                 Supply --dem <dir> to use a local tile directory.",
                dem_source,
                bb.min_lat_deg, bb.max_lat_deg,
                bb.min_lon_deg, bb.max_lon_deg,
            )
        })?;

        // Post-condition: for SRTM-1, every required tile must have been
        // downloaded (or was already cached).  This fires if the download
        // loop silently skipped a tile due to an unexpected server state.
        if dem_source == "srtm1" {
            check_srtm1_tiles_present(&dir, bb)
                .with_context(|| "SRTM-1 tile completeness check failed after auto-download")?;
        }

        return Ok(dir);
    }

    // ── Path 3: no dem-fetch feature, no override — explicit error ────────
    #[cfg(not(feature = "dem-fetch"))]
    bail!(
        "No --dem directory provided and the `dem-fetch` feature is disabled.\n\
         Supply a directory of DEM tiles with --dem <dir>, or \
         rebuild with `--features dem-fetch` to enable automatic tile download."
    );
}

/// Resolve the DEM cache directory.
///
/// Checks `$SARDINE_DEM_DIR` first; falls back to `$HOME/.sardine/dem/`.
fn dem_cache_dir() -> Result<PathBuf> {
    if let Ok(dir) = std::env::var("SARDINE_DEM_DIR") {
        return Ok(PathBuf::from(dir));
    }
    let home = std::env::var("HOME").with_context(|| {
        "HOME environment variable not set and SARDINE_DEM_DIR not set; \
         cannot locate DEM cache directory"
    })?;
    Ok(PathBuf::from(home).join(".sardine").join("dem"))
}

/// Verify that every SRTM-1 tile required to cover `bb` is present in `dir`.
///
/// SRTM-1 tiles are 1°×1° HGT files named `N{lat}E{lon}.hgt` /
/// `S{lat}W{lon}.hgt`.  Because SRTM-1 provides sea-level tiles for ocean
/// areas, a missing tile is never legitimate — it means either the DEM
/// directory is incomplete or points at the wrong scene area.
///
/// Call this after [`resolve_dem`] returns a directory for `"srtm1"` source.
/// For `"glo30"` do not call this: ocean tiles are legitimately absent (HTTP
/// 404 from the server).
///
/// # Errors
///
/// Returns a human-readable error listing every missing tile filename.
#[cfg(feature = "dem-fetch")]
pub fn check_srtm1_tiles_present(
    dir: &Path,
    bb: &crate::types::BoundingBox,
) -> Result<()> {
    use crate::dem_fetch::{tile_name, tiles_for_bbox};

    let required = tiles_for_bbox(bb.min_lat_deg, bb.max_lat_deg, bb.min_lon_deg, bb.max_lon_deg);
    let missing: Vec<String> = required
        .iter()
        .filter_map(|(lat, lon)| {
            let name = tile_name(*lat, *lon);
            if dir.join(format!("{}.hgt", name)).exists() {
                None
            } else {
                Some(name)
            }
        })
        .collect();

    if !missing.is_empty() {
        bail!(
            "DEM directory '{}' is missing {} SRTM-1 tile(s) required for this scene's \
             bounding box [{:.3}°–{:.3}°N, {:.3}°–{:.3}°E]: {}.\n\
             Either download the missing tiles or point --dem at a directory with \
             complete coverage for this scene.",
            dir.display(),
            missing.len(),
            bb.min_lat_deg, bb.max_lat_deg,
            bb.min_lon_deg, bb.max_lon_deg,
            missing.join(", ")
        );
    }
    Ok(())
}

// ─────────────────────────────────────────────────────────────────────────────
// Public entry points
// ─────────────────────────────────────────────────────────────────────────────

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
        "HH" => Polarization::HH,
        "HV" => Polarization::HV,
        other => bail!(
            "Unsupported polarization: {}. Use VV, VH, HH, or HV.",
            other
        ),
    };

    tracing::info!("parsing SAFE metadata …");
    let t_parse = Instant::now();
    let scene = parse_safe_directory(safe)
        .with_context(|| format!("parsing SAFE: {}", safe.display()))?;

    let (scene, orbit_source) = resolve_orbit(scene, orbit_path, "scene")?;

    tracing::info!("parsing calibration/noise LUTs …");
    let cal_noise = parse_calibration_noise(safe)
        .with_context(|| format!("parsing calibration/noise: {}", safe.display()))?;

    tracing::info!("parsing geolocation grids …");
    let grids = parse_geolocation_grids(safe)
        .with_context(|| format!("parsing geolocation grids: {}", safe.display()))?;
    crate::run::report_timing("parse_safe_meta", t_parse);

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
    crate::run::report_timing("read_cal_deburst_subswaths", t_swaths);

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
    crate::run::report_timing("merge_subswaths", t_merge);

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
        "HH" => Polarization::HH,
        "HV" => Polarization::HV,
        other => bail!(
            "Unsupported polarization: {}. Use VV, VH, HH, or HV.",
            other
        ),
    };

    let scene = assembled.scene.clone();
    let (scene, orbit_source) = resolve_orbit(scene, orbit_path, "scene")?;

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
    crate::run::report_timing("assembled_read_cal_deburst", t_swaths);

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
    crate::run::report_timing("assembled_merge_subswaths", t_merge);

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
