//! Provenance helpers: quality flags, speckle provenance fields, and the
//! `build_provenance` / `build_grd_provenance` constructors.

use anyhow::{anyhow, Result};

use crate::pipeline_options::{OutputMode, SpeckleOrder};
use crate::types::Polarization;

// ─────────────────────────────────────────────────────────────────────────────
// Shared helpers
// ─────────────────────────────────────────────────────────────────────────────

/// Build the quality flags list for one processing run.
///
/// Pass `dem_missing_count = None` and `non_converged_count = None` for GRD
/// mode, which has no geocoding.
pub(crate) fn collect_quality_flags(
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

// ─────────────────────────────────────────────────────────────────────────────
// Provenance builders
// ─────────────────────────────────────────────────────────────────────────────

#[allow(clippy::too_many_arguments)]
pub(crate) fn build_provenance(
    opts: &crate::pipeline_options::ProcessOptions,
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
    use crate::pipeline_options::resolve_speckle;
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
        .as_deref()
        .map(|p| p.to_str().ok_or_else(|| anyhow!("--dem path contains non-UTF-8 characters")))
        .transpose()?
        .unwrap_or("auto") // SAFETY-OK: provenance display only — None means auto-downloaded
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
                // Normalise anti-meridian crossing scenes (max_lon > 180) back
                // to standard [-180, 180] for JSON / RFC 7946 compatibility.
                max_lon: if bb.max_lon_deg > 180.0 { bb.max_lon_deg - 360.0 } else { bb.max_lon_deg },
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

#[allow(clippy::too_many_arguments)]
pub(crate) fn build_grd_provenance(
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
    use crate::pipeline_options::resolve_speckle;
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
                max_lon: if bb.max_lon_deg > 180.0 { bb.max_lon_deg - 360.0 } else { bb.max_lon_deg },
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
