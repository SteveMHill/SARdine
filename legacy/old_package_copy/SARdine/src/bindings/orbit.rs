//! Orbit bindings for Python
//!
//! Contains functions for applying precise orbit files to SAR products.

use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::PyDict;

use crate::io::download::manager::{DownloadConfig, DownloadManager};
use crate::io::orbit::{OrbitReader, OrbitType};
use chrono::{DateTime, Utc};

/// Apply precise orbit file to a SAR product
///
/// This function loads and validates orbit state vectors from a cached orbit file.
/// If the orbit file is not in the cache, it will be automatically downloaded.
///
/// # Arguments
/// * `product_id` - The Sentinel-1 product identifier
/// * `start_time` - The acquisition start time in RFC3339 format
/// * `cache_dir` - Directory containing cached orbit files
///
/// # Returns
/// A dictionary containing:
/// - `status`: "success" or error information
/// - `orbit_vectors_count`: Number of state vectors loaded
/// - `reference_time`: Orbit reference time
/// - `osv_times`: List of state vector times
/// - `osv_positions`: List of [x, y, z] position vectors
/// - `osv_velocities`: List of [vx, vy, vz] velocity vectors
#[pyfunction]
pub fn apply_precise_orbit_file(
    product_id: String,
    start_time: String,
    cache_dir: String,
) -> PyResult<std::collections::HashMap<String, PyObject>> {
    if product_id.is_empty() || start_time.is_empty() {
        return Err(PyValueError::new_err(
            "Product ID and start time cannot be empty",
        ));
    }

    // Parse start time
    let start_dt = DateTime::parse_from_rfc3339(&start_time)
        .map_err(|e| PyValueError::new_err(format!("Invalid start time format: {}", e)))?
        .with_timezone(&Utc);

    // Check cache first before attempting download
    let cache_path = std::path::Path::new(&cache_dir);
    std::fs::create_dir_all(cache_path)
        .map_err(|e| PyValueError::new_err(format!("Failed to create cache directory: {}", e)))?;

    let orbit_type = OrbitReader::determine_orbit_type(start_dt);
    let mut chosen_orbit_type = orbit_type;
    let orbit_filename = format!("{}_{}.EOF", product_id, orbit_type);
    let orbit_path = cache_path.join(&orbit_filename);

    // Also check orbits subdirectory (where DownloadManager places files)
    let orbits_subdir_path = cache_path.join("orbits").join(&orbit_filename);

    // Check if orbit file exists in cache, if not try to download
    let orbit_data = if orbit_path.exists() {
        // File exists in cache root - read and validate it directly
        log::info!("Using cached orbit file: {}", orbit_path.display());
        OrbitReader::read_orbit_file(&orbit_path).map_err(|e| {
            PyValueError::new_err(format!("Failed to read cached orbit file: {}", e))
        })?
    } else if orbits_subdir_path.exists() {
        // File exists in orbits subdirectory
        log::info!(
            "Using cached orbit file from orbits subdir: {}",
            orbits_subdir_path.display()
        );
        OrbitReader::read_orbit_file(&orbits_subdir_path).map_err(|e| {
            PyValueError::new_err(format!("Failed to read cached orbit file: {}", e))
        })?
    } else {
        // File doesn't exist - check for RESORB fallback if we were looking for POEORB
        let resorb_filename = format!("{}_{}.EOF", product_id, OrbitType::RESORB);
        let resorb_path = cache_path.join(&resorb_filename);
        let resorb_subdir_path = cache_path.join("orbits").join(&resorb_filename);

        if resorb_path.exists() {
            log::info!(
                "Using cached RESORB file as fallback: {}",
                resorb_path.display()
            );
            chosen_orbit_type = OrbitType::RESORB;
            OrbitReader::read_orbit_file(&resorb_path).map_err(|e| {
                PyValueError::new_err(format!("Failed to read cached RESORB file: {}", e))
            })?
        } else if resorb_subdir_path.exists() {
            log::info!(
                "Using cached RESORB file from orbits subdir: {}",
                resorb_subdir_path.display()
            );
            chosen_orbit_type = OrbitType::RESORB;
            OrbitReader::read_orbit_file(&resorb_subdir_path).map_err(|e| {
                PyValueError::new_err(format!("Failed to read cached RESORB file: {}", e))
            })?
        } else {
            // No orbit file found - attempt automatic download using DownloadManager
            log::info!("Orbit file not found in cache, attempting automatic download...");

            // Create DownloadManager with cache pointing to the orbit cache directory
            let mut config = DownloadConfig::default();
            config.cache_dir = cache_path.to_path_buf();

            let mut download_manager = DownloadManager::new(config).map_err(|e| {
                PyValueError::new_err(format!("Failed to initialize download manager: {}", e))
            })?;

            // Download orbit file
            match download_manager.download_orbit(&product_id, start_dt, None) {
                Ok(result) => {
                    log::info!(
                        "Successfully downloaded orbit file to: {}",
                        result.path.display()
                    );
                    let path_str = result.path.to_string_lossy();
                    if path_str.contains("RESORB") {
                        chosen_orbit_type = OrbitType::RESORB;
                    }
                    OrbitReader::read_orbit_file(&result.path).map_err(|e| {
                        PyValueError::new_err(format!(
                            "Failed to read downloaded orbit file: {}",
                            e
                        ))
                    })?
                }
                Err(e) => {
                    return Err(PyValueError::new_err(format!(
                        "Failed to download orbit file: {}. \
                        Please check network connectivity or manually download the orbit file to: {}",
                        e, orbit_path.display()
                    )));
                }
            }
        }
    };

    // Validate orbit quality
    if orbit_data.state_vectors.len() < 10 {
        return Err(PyValueError::new_err(
            "Insufficient orbit state vectors for scientific processing",
        ));
    }

    Python::with_gil(|py| {
        let py_result = PyDict::new(py);
        py_result.set_item("status", "success")?;
        py_result.set_item("orbit_vectors_count", orbit_data.state_vectors.len())?;
        py_result.set_item("reference_time", orbit_data.reference_time.to_rfc3339())?;
        let orbit_type_str = match chosen_orbit_type {
            OrbitType::POEORB => "POEORB",
            OrbitType::RESORB => "RESORB",
        };

        py_result.set_item("orbit_type", orbit_type_str)?; // Report actual orbit type used

        // Extract orbit state vectors for terrain correction
        // These are needed by Range-Doppler geocoding algorithm
        let osv_times: Vec<String> = orbit_data
            .state_vectors
            .iter()
            .map(|sv| sv.time.to_rfc3339())
            .collect();
        let osv_positions: Vec<[f64; 3]> = orbit_data
            .state_vectors
            .iter()
            .map(|sv| sv.position)
            .collect();
        let osv_velocities: Vec<[f64; 3]> = orbit_data
            .state_vectors
            .iter()
            .map(|sv| sv.velocity)
            .collect();

        py_result.set_item("osv_times", osv_times)?;
        py_result.set_item("osv_positions", osv_positions)?;
        py_result.set_item("osv_velocities", osv_velocities)?;

        let mut result_map = std::collections::HashMap::new();
        result_map.insert("result".to_string(), py_result.into());
        Ok(result_map)
    })
}

/// Load orbit file directly from disk
///
/// # Arguments
/// * `orbit_file_path` - Path to the orbit EOF file
///
/// # Returns
/// A dictionary containing parsed orbit state vectors
#[pyfunction]
pub fn load_orbit_file(py: Python, orbit_file_path: String) -> PyResult<PyObject> {
    let orbit_path = std::path::Path::new(&orbit_file_path);

    if !orbit_path.exists() {
        return Err(PyValueError::new_err(format!(
            "Orbit file not found: {}",
            orbit_file_path
        )));
    }

    let orbit_data = OrbitReader::read_orbit_file(orbit_path)
        .map_err(|e| PyValueError::new_err(format!("Failed to read orbit file: {}", e)))?;

    let result = PyDict::new(py);
    result.set_item("status", "success")?;
    result.set_item("orbit_vectors_count", orbit_data.state_vectors.len())?;
    result.set_item("reference_time", orbit_data.reference_time.to_rfc3339())?;

    // Extract state vectors
    let osv_times: Vec<String> = orbit_data
        .state_vectors
        .iter()
        .map(|sv| sv.time.to_rfc3339())
        .collect();
    let osv_positions: Vec<[f64; 3]> = orbit_data
        .state_vectors
        .iter()
        .map(|sv| sv.position)
        .collect();
    let osv_velocities: Vec<[f64; 3]> = orbit_data
        .state_vectors
        .iter()
        .map(|sv| sv.velocity)
        .collect();

    result.set_item("osv_times", osv_times)?;
    result.set_item("osv_positions", osv_positions)?;
    result.set_item("osv_velocities", osv_velocities)?;

    Ok(result.into())
}
