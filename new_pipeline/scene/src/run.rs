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
//! ## Module layout
//!
//! The implementation is split across four sub-modules that are re-exported here:
//!
//! * [`crate::pipeline_options`] — option structs, resolver helpers, and CLI parsers.
//! * [`crate::scene_prep`] — deburst / calibrate / merge / multilook helpers.
//! * [`crate::run_provenance`] — quality flags and provenance builders.
//! * [`crate::multi_pol`] — dual-polarization wrappers.
//!
//! ## Domain notes
//!
//! These functions are direct extractions of the CLI subcommands; no
//! algorithmic change is intended.  All field semantics and default values
//! match the CLI 1:1.  Progress and warning messages are emitted via
//! `tracing::info!` / `tracing::warn!`; callers must initialise a
//! [`tracing_subscriber`] to see output (the CLI binary does this automatically).

#[allow(unused_imports)]
use std::path::{Path, PathBuf};
use std::time::Instant;

use anyhow::{anyhow, Context, Result};

use crate::types::Polarization;

// Re-export the public API of the sub-modules so that callers continue to
// access everything via `crate::run::*` without breaking changes.
pub use crate::pipeline_options::{
    derive_pol_output_path, parse_iw_selection, parse_output_mode, parse_polarizations,
    parse_speckle_order, resolve_crs, resolve_geoid, resolve_speckle, sidecar_path, GrdOptions,
    IwSelection, OutputMode, ProcessOptions, ResamplingKernel, SpeckleOrder,
};
pub use crate::run_provenance::speckle_provenance_fields;
pub(crate) use crate::run_provenance::collect_quality_flags;
pub use crate::scene_prep::{prepare_merged_scene, prepare_merged_scene_assembled, PreparedScene};
pub(crate) use crate::scene_prep::apply_multilook;
pub use crate::multi_pol::{run_grd_multi, run_process_multi};

/// Emit a single `[timing]` line on stderr.  Centralised so the format stays
/// consistent across stages and so future migration to a structured
/// timing record (sidecar JSON, tracing span) only has to touch one
/// place.  All times are wall-clock seconds since the supplied `Instant`.
pub(crate) fn report_timing(stage: &str, started: Instant) {
    tracing::info!(
        "[timing] {} = {:.3} s",
        stage,
        started.elapsed().as_secs_f64(),
    );
}

// ─────────────────────────────────────────────────────────────────────────────
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

// ─────────────────────────────────────────────────────────────────────────────
// Pipeline helpers (remain in run.rs)
// ─────────────────────────────────────────────────────────────────────────────

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

/// Run the full backscatter pipeline for each entry in `batch` sequentially,
/// returning a `Vec` of per-scene results.
///
/// Scenes are processed one after the other in the order they appear in the
/// slice.  Each scene is independent: a failure in scene `i` does not prevent
/// scenes `i+1…` from being attempted.  The caller can zip the results with
/// the input slice to identify which scenes succeeded.
///
/// # Why sequential?
///
/// `run_process` itself runs multi-threaded (rayon) within each scene.
/// Running two scenes simultaneously would create two rayon thread pools
/// competing for the same physical cores, degrading both.  Sequential batch
/// processing keeps core utilisation clean and makes the progress log readable
/// (one scene's log stream at a time).
pub fn run_process_batch(batch: &[ProcessOptions]) -> Vec<Result<()>> {
    batch
        .iter()
        .enumerate()
        .map(|(i, opts)| {
            tracing::info!(
                "batch scene {}/{}: {}",
                i + 1,
                batch.len(),
                opts.output.display(),
            );
            crate::multi_pol::run_process_multi(opts)
        })
        .collect()
}

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
    tc_config.resampling = opts.resampling;
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
    let t_export = Instant::now();
    if opts.cog {
        tracing::info!("writing COG → {}", output_str);
        crate::export::write_cog_with_crs(
            output_str,
            &geocoded.data,
            geocoded.cols,
            geocoded.rows,
            geocoded.geotransform,
            &geocoded.crs,
            512,
        )
        .with_context(|| format!("writing COG: {}", output_str))?;
    } else {
        tracing::info!("writing GeoTIFF → {}", output_str);
        write_geotiff_with_crs(
            output_str,
            &geocoded.data,
            geocoded.cols,
            geocoded.rows,
            geocoded.geotransform,
            &geocoded.crs,
        )
        .with_context(|| format!("writing GeoTIFF: {}", output_str))?;
    }
    report_timing("export_geotiff", t_export);

    if !opts.no_provenance {
        let prov = crate::run_provenance::build_provenance(
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
        let prov = crate::run_provenance::build_grd_provenance(
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
            resampling: crate::run::ResamplingKernel::default(),
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
