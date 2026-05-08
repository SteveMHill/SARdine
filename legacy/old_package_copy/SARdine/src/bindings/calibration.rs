#![allow(dead_code)]
//! Calibration bindings for Python
//!
//! Contains PyCalibrationJob class and related calibration functions.

use ndarray::Array2;
use pyo3::prelude::*;
use std::sync::Arc;
use std::time::Instant;

use crate::types::{Polarization, SarError, SarResult};

// ============================================================================
// Types
// ============================================================================

pub(crate) struct CalibrationJob {
    pub subswath: String,
    pub polarization: String,
    pub calibration_type: crate::core::calibration::CalibrationType,
    pub enable_noise_removal: bool,
    pub calibration_coeffs: Arc<crate::core::calibration::CalibrationCoefficients>,
    pub noise_xml: Option<String>,
    pub range_sample_origin: usize,
}

pub(crate) struct CalibrationRunResult {
    pub calibrated: Array2<f32>,
    pub lines: usize,
    pub samples: usize,
    pub valid_pixels: usize,
    pub total_pixels: usize,
    pub min_value: f32,
    pub max_value: f32,
    pub mean_value: f32,
    pub noise_applied: bool,
}

// ============================================================================
// PyCalibrationJob class
// ============================================================================

#[pyclass(name = "CalibrationJob", module = "sardine._core")]
#[derive(Clone)]
pub struct PyCalibrationJob {
    pub(crate) inner: Arc<CalibrationJob>,
}

#[pymethods]
impl PyCalibrationJob {
    #[getter]
    fn subswath(&self) -> String {
        self.inner.subswath.clone()
    }

    #[getter]
    fn polarization(&self) -> String {
        self.inner.polarization.clone()
    }

    #[getter]
    pub fn calibration_type(&self) -> String {
        match self.inner.calibration_type {
            crate::core::calibration::CalibrationType::Sigma0 => "sigma0".to_string(),
            crate::core::calibration::CalibrationType::Beta0 => "beta0".to_string(),
            crate::core::calibration::CalibrationType::Gamma0 => "gamma0".to_string(),
            crate::core::calibration::CalibrationType::Dn => "dn".to_string(),
        }
    }
}

// ============================================================================
// Helper functions
// ============================================================================

pub(crate) fn derive_incidence_angles_from_annotation(
    reader: &crate::io::SlcReader,
    pol: Polarization,
) -> Option<(f32, f32)> {
    let annotations = reader.get_all_cached_annotations(pol).ok()?;

    let mut min_angle = f64::INFINITY;
    let mut max_angle = f64::NEG_INFINITY;

    for annotation in annotations.iter() {
        for (_, incidence, _) in &annotation.get_antenna_patterns().ok()? {
            for &val in incidence.iter() {
                min_angle = min_angle.min(val);
                max_angle = max_angle.max(val);
            }
        }

        if let Some(mid) = annotation
            .image_annotation
            .as_ref()
            .and_then(|img| img.image_information.as_ref())
            .and_then(|info| info.incidence_angle_mid_swath)
        {
            min_angle = min_angle.min(mid);
            max_angle = max_angle.max(mid);
        }
    }

    if min_angle.is_finite() && max_angle.is_finite() && max_angle > min_angle {
        Some((min_angle as f32, max_angle as f32))
    } else {
        None
    }
}

pub(crate) fn configure_calibration_coefficients(
    coeffs: &mut crate::core::calibration::CalibrationCoefficients,
    reader: &mut crate::io::SlcReader,
    pol: Polarization,
    subswath: &str,
    image_dims: (usize, usize),
) -> SarResult<()> {
    use crate::core::calibration::model::EllipsoidIncidenceModel;
    use std::convert::TryFrom;

    let subswath_key = subswath.to_string();
    let burst_start_line = reader
        .get_cached_metadata()
        .ok()
        .and_then(|metadata| metadata.sub_swaths.get(&subswath_key).cloned())
        .and_then(|swath| i32::try_from(swath.first_line_global).ok())
        .ok_or_else(|| {
            SarError::Processing(format!(
                "Subswath {} not found in cached metadata — cannot build coordinate mapper for calibration",
                subswath_key
            ))
        })?;
    let image_start_line = burst_start_line;

    let mapper =
        coeffs.create_auto_coordinate_mapper(burst_start_line, image_start_line, image_dims.1)?;
    coeffs.set_coordinate_mapper(mapper)?;

    if let Some((near_deg, far_deg)) = derive_incidence_angles_from_annotation(reader, pol) {
        let model = EllipsoidIncidenceModel::new(near_deg, far_deg, image_dims.1);
        coeffs.set_incidence_model(Box::new(model))?;
        log::info!(
            "📐 Incidence model configured from annotation: near={:.2}°, far={:.2}°",
            near_deg,
            far_deg
        );
    } else {
        log::warn!(
            "⚠️  Unable to derive incidence angle range from annotation; continuing without incidence model"
        );
    }

    // Antenna pattern correction
    match reader.find_all_annotation_files() {
        Ok(map) => {
            if let Some(files) = map.get(&pol) {
                let subswath_lc = subswath.to_ascii_lowercase();
                let chosen_path = files
                    .iter()
                    .find(|p| p.to_ascii_lowercase().contains(&subswath_lc))
                    .or_else(|| files.first())
                    .cloned();

                if let Some(path) = chosen_path {
                    match reader.read_file_as_string(&path) {
                        Ok(xml) => {
                            if let Err(e) = coeffs.parse_and_set_antenna_pattern(&xml, None) {
                                log::warn!(
                                    "⚠️  Failed to parse antenna pattern from {}: {}",
                                    path,
                                    e
                                );
                            } else if let Err(e) = coeffs.precompute_antenna_pattern_lut(image_dims)
                            {
                                log::warn!("⚠️  Failed to precompute antenna pattern LUT: {}", e);
                            }
                        }
                        Err(e) => log::warn!("⚠️  Could not read annotation XML: {}", e),
                    }
                }
            }
        }
        Err(e) => log::warn!("⚠️  Failed to list annotation files: {}", e),
    }

    Ok(())
}

pub(crate) fn prepare_calibration_coefficients(
    reader: &mut crate::io::SlcReader,
    pol: Polarization,
    subswath: &str,
    image_dims: (usize, usize),
) -> SarResult<crate::core::calibration::CalibrationCoefficients> {
    use crate::core::calibration::CalibrationCoefficients;

    let timer_prepare = Instant::now();
    log::info!(
        "🚀 Step B: preparing calibration coefficients for {} ({:?})",
        subswath,
        pol
    );

    let mut coeffs: CalibrationCoefficients = reader.get_cached_calibration(pol)?.clone();
    log::info!(
        "⌁ Step B: cloned cached calibration vectors ({} vectors)",
        coeffs.vectors.len()
    );

    configure_calibration_coefficients(&mut coeffs, reader, pol, subswath, image_dims)?;

    // Always use dense 2D calibration LUT
    coeffs.precompute_lut(image_dims)?;

    log::info!(
        "⏱️ Step B: calibration coefficient prep finished in {:.2?}",
        timer_prepare.elapsed()
    );
    Ok(coeffs)
}

// REMOVED: prepare_calibration_coefficients_separable (Jan 2026 - dead code, separable LUT removed)

pub(crate) fn build_calibration_job_from_reader(
    reader: &mut crate::io::SlcReader,
    pol: Polarization,
    subswath: &str,
    calibration_type: crate::core::calibration::CalibrationType,
    image_dims: (usize, usize),
    enable_noise_removal: bool,
    range_sample_origin: usize,
) -> SarResult<CalibrationJob> {
    let coeffs = prepare_calibration_coefficients(reader, pol, subswath, image_dims)?;

    let noise_xml = if enable_noise_removal {
        let noise_files = reader.find_noise_files()?;
        let noise_file = noise_files.get(&pol).ok_or_else(|| {
            SarError::Processing(format!("No noise file found for polarization {:?}", pol))
        })?;
        Some(reader.read_file_as_string(noise_file)?)
    } else {
        None
    };

    Ok(CalibrationJob {
        subswath: subswath.to_string(),
        polarization: pol.to_string(),
        calibration_type,
        enable_noise_removal,
        calibration_coeffs: Arc::new(coeffs),
        noise_xml,
        range_sample_origin,
    })
}

// REMOVED: build_calibration_job_from_reader_separable (Jan 2026 - dead code, separable LUT removed)

// NOTE: run_calibration_job_impl and the pyfunctions remain in lib.rs for now
// due to their extensive dependencies and would be moved in a subsequent refactor
