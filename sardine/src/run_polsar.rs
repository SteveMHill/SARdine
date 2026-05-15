//! Pipeline entry point for H/A/Alpha dual-pol decomposition.
//!
//! Implements `run_polsar(opts: &ProcessOptions) -> Result<()>`.
//!
//! Pipeline stages:
//!
//! 1. Parse SAFE metadata + orbit.
//! 2. Validate dual-pol availability (VV+VH or HH+HV required).
//! 3. For each IW subswath: deburst both pols, complex-calibrate, build C2.
//! 4. Merge per-subswath C2 arrays into one slant-range image.
//! 5. Multilook the C2 image.
//! 6. Compute H, A, α (mean alpha) for every pixel.
//! 7. Geocode each band with terrain correction and write GeoTIFF / COG.
//! 8. Write provenance JSON sidecar.

use std::path::{Path, PathBuf};
use std::time::Instant;

use anyhow::{bail, Context, Result};
use rayon::prelude::*;

use crate::apply_calibration::apply_calibration_complex;
use crate::deburst::deburst_subswath;
use crate::dem::DemMosaic;
use crate::export::write_geotiff_with_crs;
use crate::merge_subswaths::{build_c2_from_complex, merge_c2_subswaths, C2SubswathData, C2SwathInput};
use crate::merge_subswaths::MergedSigma0;
use crate::pipeline_options::ProcessOptions;
use crate::polsar::{compute_haa_image, multilook_c2};
use crate::run::report_timing;
use crate::scene_prep::{resolve_dem, resolve_orbit};
use crate::slc_reader::SlcReader;
use crate::terrain_correction::{terrain_correction, TerrainCorrectionConfig};
use crate::types::{Polarization, SubSwathId};

// ─────────────────────────────────────────────────────────────────────────────

/// Build a polsar output file path by inserting a band tag before the extension.
///
/// `out.tif` + `"H"` → `out_H.tif`, `out_dir/scene` + `"alpha"` → `out_dir/scene_alpha`.
fn band_path(output: &Path, tag: &str) -> Result<String> {
    let stem = output
        .file_stem()
        .ok_or_else(|| anyhow::anyhow!("output path has no filename: {}", output.display()))?
        .to_str()
        .ok_or_else(|| anyhow::anyhow!("output path is not valid UTF-8: {}", output.display()))?;
    let name = match output.extension().and_then(|e| e.to_str()) {
        Some(ext) => format!("{stem}_{tag}.{ext}"),
        None => format!("{stem}_{tag}"),
    };
    let full: PathBuf = match output.parent() {
        Some(p) if !p.as_os_str().is_empty() => p.join(name),
        _ => PathBuf::from(name),
    };
    full.to_str()
        .ok_or_else(|| anyhow::anyhow!("constructed path is not valid UTF-8"))
        .map(|s| s.to_owned())
}

/// Find a measurement TIFF under `measurement_dir` whose lower-case name
/// contains both `iw` and `pol`.
fn find_tiff(measurement_dir: &Path, iw: &str, pol: &str) -> Result<PathBuf> {
    for entry in std::fs::read_dir(measurement_dir)
        .with_context(|| format!("reading measurement dir: {}", measurement_dir.display()))?
    {
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
        measurement_dir.display()
    )
}

// ─────────────────────────────────────────────────────────────────────────────
// Public entry point
// ─────────────────────────────────────────────────────────────────────────────

/// Run the H/A/Alpha dual-pol decomposition pipeline.
///
/// # Errors
///
/// Returns an error if:
/// - The SAFE product does not contain two co-pol/cross-pol channels.
/// - DEM or geoid resolution fails.
/// - Any I/O or parse step fails.
pub fn run_polsar(opts: &ProcessOptions) -> Result<()> {
    use crate::calibration::parse_calibration_noise;
    use crate::parse::{parse_geolocation_grids, parse_safe_directory};
    use crate::run::resolve_crs;
    use crate::pipeline_options::{resolve_geoid, sidecar_path};

    if opts.threads > 0 {
        rayon::ThreadPoolBuilder::new()
            .num_threads(opts.threads)
            .build_global()
            .with_context(|| format!("setting Rayon thread count to {}", opts.threads))?;
        tracing::info!("using {} thread(s)", opts.threads);
    }

    let t_total = Instant::now();

    // ── Parse SAFE ────────────────────────────────────────────────────────────
    tracing::info!("parsing SAFE metadata …");
    let t_parse = Instant::now();
    let scene_raw = parse_safe_directory(&opts.safe)
        .with_context(|| format!("parsing SAFE: {}", opts.safe.display()))?;

    // ── Validate dual-pol requirement ─────────────────────────────────────────
    //
    // H/A/Alpha decomposition requires two complex channels that form the
    // scattering vector k = [S_co, S_cross]ᵀ.
    // For S-1 IW this is VV+VH (descending/ascending) or HH+HV.
    let pols: Vec<Polarization> = {
        let mut p = scene_raw.polarizations.clone();
        p.sort();
        p
    };

    let (pol_co, pol_cross) = match pols.as_slice() {
        [Polarization::VH, Polarization::VV] | [Polarization::VV, Polarization::VH] => {
            (Polarization::VV, Polarization::VH)
        }
        [Polarization::HH, Polarization::HV] | [Polarization::HV, Polarization::HH] => {
            (Polarization::HH, Polarization::HV)
        }
        other => {
            bail!(
                "H/A/Alpha decomposition requires dual-pol VV+VH or HH+HV, but the \
                 product contains: {:?}. A minimum of two complementary channels is needed.",
                other
            );
        }
    };

    tracing::info!("dual-pol channels: {} (co-pol) + {} (cross-pol)", pol_co, pol_cross);

    let (scene, orbit_source) = resolve_orbit(scene_raw, opts.orbit.as_deref(), "scene")?;

    let cal_noise = parse_calibration_noise(&opts.safe)
        .with_context(|| format!("parsing calibration/noise: {}", opts.safe.display()))?;

    let grids = parse_geolocation_grids(&opts.safe)
        .with_context(|| format!("parsing geolocation grids: {}", opts.safe.display()))?;
    report_timing("parse_safe_meta", t_parse);

    // ── Subswath loop: deburst + complex calibrate + build C2 ────────────────
    let iw_ids = [SubSwathId::IW1, SubSwathId::IW2, SubSwathId::IW3];
    let measurement_dir = opts.safe.join("measurement");

    struct IwC2Data {
        c2: C2SubswathData,
        swath_meta: crate::types::SubSwathMetadata,
        azimuth_start_time: chrono::DateTime<chrono::Utc>,
    }

    let t_swaths = Instant::now();
    let c2_per_iw: Vec<IwC2Data> = iw_ids
        .par_iter()
        .map(|&iw_id| -> Result<IwC2Data> {
            let iw_str = format!("{}", iw_id).to_lowercase();
            let co_str = format!("{}", pol_co).to_lowercase();
            let cx_str = format!("{}", pol_cross).to_lowercase();

            let tiff_co = find_tiff(&measurement_dir, &iw_str, &co_str)?;
            let tiff_cx = find_tiff(&measurement_dir, &iw_str, &cx_str)?;

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
                b
            };
            if bursts.is_empty() {
                bail!("no bursts for subswath {:?}", iw_id);
            }
            let azimuth_start_time = bursts[0].azimuth_time_utc;
            let tiff_line_origin = bursts[0].first_line;

            // Calibration LUTs for both pols.
            let cal_co = cal_noise
                .calibrations
                .iter()
                .find(|c| c.subswath_id == iw_id && c.polarization == pol_co)
                .with_context(|| format!("cal LUT for {:?} {} missing", iw_id, pol_co))?;
            let noise_co = cal_noise
                .noises
                .iter()
                .find(|n| n.subswath_id == iw_id && n.polarization == pol_co)
                .with_context(|| format!("noise LUT for {:?} {} missing", iw_id, pol_co))?;
            let cal_cx = cal_noise
                .calibrations
                .iter()
                .find(|c| c.subswath_id == iw_id && c.polarization == pol_cross)
                .with_context(|| format!("cal LUT for {:?} {} missing", iw_id, pol_cross))?;
            let noise_cx = cal_noise
                .noises
                .iter()
                .find(|n| n.subswath_id == iw_id && n.polarization == pol_cross)
                .with_context(|| format!("noise LUT for {:?} {} missing", iw_id, pol_cross))?;

            // Deburst + complex calibrate co-pol.
            tracing::info!("debursting {:?} {} …", iw_id, pol_co);
            let mut reader_co = SlcReader::open(&tiff_co)
                .with_context(|| format!("opening co-pol TIFF: {}", tiff_co.display()))?;
            let deburst_co = deburst_subswath(&mut reader_co, &swath_meta, &bursts)
                .with_context(|| format!("debursting {:?} {}", iw_id, pol_co))?;
            tracing::info!("calibrating {:?} {} (complex) …", iw_id, pol_co);
            let complex_co = apply_calibration_complex(&deburst_co, cal_co, noise_co, tiff_line_origin)
                .with_context(|| format!("complex calibration {:?} {}", iw_id, pol_co))?;

            // Deburst + complex calibrate cross-pol.
            tracing::info!("debursting {:?} {} …", iw_id, pol_cross);
            let mut reader_cx = SlcReader::open(&tiff_cx)
                .with_context(|| format!("opening cross-pol TIFF: {}", tiff_cx.display()))?;
            let deburst_cx = deburst_subswath(&mut reader_cx, &swath_meta, &bursts)
                .with_context(|| format!("debursting {:?} {}", iw_id, pol_cross))?;
            tracing::info!("calibrating {:?} {} (complex) …", iw_id, pol_cross);
            let complex_cx = apply_calibration_complex(&deburst_cx, cal_cx, noise_cx, tiff_line_origin)
                .with_context(|| format!("complex calibration {:?} {}", iw_id, pol_cross))?;

            // Build C2 covariance at full SLC resolution.
            tracing::info!("building C2 for {:?} …", iw_id);
            let c2 = build_c2_from_complex(&complex_co, &complex_cx);

            Ok(IwC2Data { c2, swath_meta, azimuth_start_time })
        })
        .collect::<Result<Vec<_>>>()?;
    report_timing("deburst_cal_c2_subswaths", t_swaths);

    // ── Merge C2 subswaths ────────────────────────────────────────────────────
    tracing::info!("merging C2 subswaths …");
    let t_merge = Instant::now();
    let c2_inputs: Vec<C2SwathInput<'_>> = c2_per_iw
        .iter()
        .map(|d| C2SwathInput {
            c2: &d.c2,
            swath: &d.swath_meta,
            azimuth_start_time: d.azimuth_start_time,
        })
        .collect();
    let mut merged_c2 = merge_c2_subswaths(&c2_inputs)
        .with_context(|| "merging C2 subswaths")?;
    report_timing("merge_c2", t_merge);

    // ── Multilook C2 ─────────────────────────────────────────────────────────
    if opts.multilook_range > 1 || opts.multilook_azimuth > 1 {
        tracing::info!(
            "multilooking C2 {}×{} (range×azimuth) ({} × {} → {} × {}) …",
            opts.multilook_range,
            opts.multilook_azimuth,
            merged_c2.samples,
            merged_c2.lines,
            merged_c2.samples / opts.multilook_range,
            merged_c2.lines / opts.multilook_azimuth,
        );
        let t_ml = Instant::now();
        merged_c2 = multilook_c2(&merged_c2, opts.multilook_range, opts.multilook_azimuth);
        report_timing("multilook_c2", t_ml);
    }

    // ── Compute H, A, alpha ───────────────────────────────────────────────────
    tracing::info!("computing H/A/Alpha ({} × {} pixels) …", merged_c2.lines, merged_c2.samples);
    let t_haa = Instant::now();
    let (h_buf, a_buf, alpha_buf) = compute_haa_image(&merged_c2);
    report_timing("haa", t_haa);

    // ── Set up terrain correction ─────────────────────────────────────────────
    tracing::info!("resolving DEM …");
    let t_dem = Instant::now();
    let dem_dir = resolve_dem(opts.dem.as_deref(), &opts.dem_source, &scene.bounding_box)
        .with_context(|| "resolving DEM")?;
    let dem = DemMosaic::load_directory(&dem_dir)
        .with_context(|| format!("loading DEM from: {}", dem_dir.display()))?;
    report_timing("dem_load", t_dem);

    let dem_tile_count = dem.tile_count();
    let dem_dir_str = dem_dir.to_str().unwrap_or("").to_owned();

    let bb = &scene.bounding_box;
    dem.covers_bbox(bb.min_lat_deg, bb.max_lat_deg, bb.min_lon_deg, bb.max_lon_deg, 0.05)
        .with_context(|| format!("DEM tiles in {} do not cover the scene footprint", dem_dir.display()))?;

    let geoid = resolve_geoid(&opts.geoid)?;
    let scene_centre = (
        0.5 * (bb.min_lon_deg + bb.max_lon_deg),
        0.5 * (bb.min_lat_deg + bb.max_lat_deg),
    );
    let output_crs = resolve_crs(&opts.crs, Some(scene_centre))?;
    let spacing = if output_crs.is_metric() { opts.pixel_spacing_m } else { opts.pixel_spacing_deg };

    let mut tc_config = TerrainCorrectionConfig::new(geoid);
    tc_config.crs = output_crs;
    tc_config.pixel_spacing_deg = spacing;
    tc_config.flatten = false; // H/A/Alpha are scattering mechanism parameters, not backscatter intensity
    tc_config.resampling = opts.resampling;
    tc_config.compute_lia = false;

    // ── Geocode each band ─────────────────────────────────────────────────────
    let n = merged_c2.lines * merged_c2.samples;

    // Collect (output_path, geocoded) for each band so we can build provenance
    // and STAC sidecars after all three bands have been written.
    let mut band_results: Vec<(String, crate::terrain_correction::GeocodedImage)> =
        Vec::with_capacity(3);

    for (tag, buf) in [("H", &h_buf), ("A", &a_buf), ("alpha", &alpha_buf)] {
        tracing::info!("geocoding band {} …", tag);
        let t_tc = Instant::now();

        // Wrap the flat buffer as a MergedSigma0 so it can be passed to
        // terrain_correction via the RadarImage trait (already implemented for MergedSigma0).
        let band_image = MergedSigma0 {
            data: buf.clone(),
            nesz: vec![0.0f32; n],
            lines: merged_c2.lines,
            samples: merged_c2.samples,
            near_slant_range_time_s: merged_c2.near_slant_range_time_s,
            range_pixel_spacing_s: merged_c2.range_pixel_spacing_s,
            range_pixel_spacing_m: merged_c2.range_pixel_spacing_m,
            cal_lut_extrapolation_gap_px: 0,
            noise_lut_extrapolation_gap_px: 0,
            azimuth_start_time: merged_c2.azimuth_start_time,
            azimuth_time_interval_s: merged_c2.azimuth_time_interval_s,
        };

        let geocoded = terrain_correction(&band_image, &scene, &dem, &grids, &tc_config)
            .with_context(|| format!("terrain correction for band {tag}"))?;
        report_timing(&format!("terrain_correction_{tag}"), t_tc);

        tracing::info!(
            "geocoded {} band: {}×{} pixels, valid={}, dem_missing={}, not_converged={}",
            tag,
            geocoded.cols,
            geocoded.rows,
            geocoded.valid_pixel_count,
            geocoded.dem_missing_count,
            geocoded.non_converged_count,
        );

        let out_path = band_path(&opts.output, tag)?;
        if opts.cog {
            tracing::info!("writing COG band {} → {}", tag, out_path);
            crate::export::write_cog_with_crs(
                &out_path,
                &geocoded.data,
                geocoded.cols,
                geocoded.rows,
                geocoded.geotransform,
                &geocoded.crs,
                512,
            )
            .with_context(|| format!("writing COG for band {tag}: {out_path}"))?;
        } else {
            tracing::info!("writing GeoTIFF band {} → {}", tag, out_path);
            write_geotiff_with_crs(
                &out_path,
                &geocoded.data,
                geocoded.cols,
                geocoded.rows,
                geocoded.geotransform,
                &geocoded.crs,
            )
            .with_context(|| format!("writing GeoTIFF for band {tag}: {out_path}"))?;
        }

        band_results.push((out_path, geocoded));
    }

    // band_results[0] = H, [1] = A, [2] = alpha (insertion order matches loop above).
    let (h_path, geocoded_h) = &band_results[0];
    let (a_path, _) = &band_results[1];
    let (alpha_path, _) = &band_results[2];

    // ── Provenance JSON ───────────────────────────────────────────────────────
    if !opts.no_provenance {
        let quality_flags = crate::run_provenance::collect_quality_flags(
            orbit_source,
            0u32, // polsar complex calibration path does not track per-pixel LUT gaps
            0u32,
            merged_c2.lines * merged_c2.samples,
            Some(geocoded_h.dem_missing_count),
            Some(geocoded_h.non_converged_count),
        );
        let prov = crate::run_provenance::build_polsar_provenance(
            opts,
            &scene,
            orbit_source,
            &dem_dir_str,
            dem_tile_count,
            geocoded_h,
            quality_flags,
            pol_co,
            pol_cross,
            opts.multilook_azimuth,
            opts.multilook_range,
        )?;
        let prov_path = sidecar_path(&opts.output, ".polsar.provenance.json")?;
        tracing::info!("writing provenance → {}", prov_path);
        prov.write_json(std::path::Path::new(&prov_path))
            .with_context(|| format!("writing provenance to {prov_path}"))?;

        // ── STAC sidecar ──────────────────────────────────────────────────────
        let stac_path = sidecar_path(&opts.output, ".polsar.stac.json")?;
        tracing::info!("writing STAC sidecar → {}", stac_path);
        crate::stac::write_polsar_stac_item(
            &crate::stac::PolsarStacInput {
                product_id: &scene.product_id,
                platform: &scene.mission.to_string().to_lowercase(),
                pol_co: &pol_co.to_string().to_uppercase(),
                pol_cross: &pol_cross.to_string().to_uppercase(),
                scene_start_utc: &scene.start_time.format("%Y-%m-%dT%H:%M:%SZ").to_string(),
                scene_stop_utc: &scene.stop_time.format("%Y-%m-%dT%H:%M:%SZ").to_string(),
                az_looks: opts.multilook_azimuth,
                rg_looks: opts.multilook_range,
                crs_epsg: geocoded_h.crs.epsg(),
                geotransform: geocoded_h.geotransform,
                cols: geocoded_h.cols,
                rows: geocoded_h.rows,
                orbit_is_poeorb: opts.orbit.is_some(),
                acquisition_mode: &scene.acquisition_mode.to_string().to_lowercase(),
                dem_source: &opts.dem_source,
                geoid_spec: &opts.geoid,
                orbit_pass_direction: &scene.orbit_pass_direction,
                absolute_orbit_number: scene.absolute_orbit_number,
                h_path,
                a_path,
                alpha_path,
            },
            std::path::Path::new(&stac_path),
        )
        .with_context(|| format!("writing STAC sidecar to {stac_path}"))?;
    }

    tracing::info!("H/A/Alpha decomposition complete.");
    report_timing("total", t_total);
    Ok(())
}
