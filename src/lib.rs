//! SARdine: A Fast, Modular Sentinel-1 Backscatter Processor
//! 
//! This library provides a modern, open-source alternative to ESA SNAP and GAMMA
//! for processing Sentinel-1 SLC data into calibrated, terrain-corrected backscatter products.

use pyo3::prelude::*;

pub mod types;
pub mod io;
pub mod core;

// Re-export main types and functions for easier access
pub use types::{
    SarImage, SarRealImage, SarMetadata, SarProduct, SarError, SarResult,
    Polarization, AcquisitionMode, CoordinateSystem
};

pub use io::{SlcReader, OrbitReader, DemReader};

/// Python module definition
#[pymodule]
fn _core(_py: Python, m: &PyModule) -> PyResult<()> {
    // Add Python bindings here
    m.add_class::<PySlcReader>()?;
    m.add_class::<PyPolarization>()?;
    m.add_class::<PyMetadata>()?;
    m.add_class::<PyOrbitData>()?;
    m.add_class::<PyStateVector>()?;
    m.add_class::<PySubSwath>()?;
    Ok(())
}

/// Python wrapper for SlcReader
#[pyclass(name = "SlcReader")]
struct PySlcReader {
    inner: SlcReader,
}

#[pymethods]
impl PySlcReader {
    #[new]
    fn new(zip_path: String) -> PyResult<Self> {
        let reader = SlcReader::new(&zip_path)
            .map_err(|e| PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(format!("{}", e)))?;
        
        Ok(PySlcReader { inner: reader })
    }
    
    fn list_files(&mut self) -> PyResult<Vec<String>> {
        self.inner.list_files()
            .map_err(|e| PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(format!("{}", e)))
    }
    
    fn get_metadata(&mut self, polarization: &str) -> PyResult<PyMetadata> {
        let pol = match polarization.to_uppercase().as_str() {
            "VV" => Polarization::VV,
            "VH" => Polarization::VH,
            "HV" => Polarization::HV,
            "HH" => Polarization::HH,
            _ => return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(
                format!("Invalid polarization: {}", polarization)
            )),
        };
        
        let metadata = self.inner.read_annotation(pol)
            .map_err(|e| PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(format!("{}", e)))?;
        
        Ok(PyMetadata { inner: metadata })
    }
    
    fn check_orbit_status(&mut self, cache_dir: Option<String>) -> PyResult<String> {
        let cache_path = match cache_dir {
            Some(dir) => Some(std::path::PathBuf::from(dir)),
            None => None,
        };
        
        let status = self.inner.check_orbit_status(cache_path.as_deref())
            .map_err(|e| PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(format!("{}", e)))?;
        
        Ok(format!("{:?}", status))
    }
    
    fn download_orbit_files(&mut self, cache_dir: Option<String>) -> PyResult<Vec<String>> {
        let cache_path = match cache_dir {
            Some(dir) => Some(std::path::PathBuf::from(dir)),
            None => None,
        };
        
        let result = self.inner.download_orbit_files(cache_path.as_deref())
            .map_err(|e| PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(format!("{}", e)))?;
        
        Ok(result.into_iter().map(|p| p.to_string_lossy().to_string()).collect())
    }
    
    fn get_orbit_data(&mut self, cache_dir: Option<String>) -> PyResult<PyOrbitData> {
        let cache_path = match cache_dir {
            Some(dir) => Some(std::path::PathBuf::from(dir)),
            None => None,
        };
        
        let orbit_data = self.inner.get_orbit_data(cache_path.as_deref())
            .map_err(|e| PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(format!("{}", e)))?;
        
        Ok(PyOrbitData { inner: orbit_data })
    }
    
    fn get_satellite_position_at_pixel(&mut self, polarization: &str, azimuth_line: usize, range_sample: usize, cache_dir: Option<String>) -> PyResult<(f64, f64, f64, f64, f64, f64)> {
        let pol = match polarization.to_uppercase().as_str() {
            "VV" => Polarization::VV,
            "VH" => Polarization::VH,
            "HV" => Polarization::HV,
            "HH" => Polarization::HH,
            _ => return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(
                format!("Invalid polarization: {}", polarization)
            )),
        };
        
        let cache_path = match cache_dir {
            Some(dir) => Some(std::path::PathBuf::from(dir)),
            None => None,
        };
        
        let (position, velocity) = self.inner.get_satellite_position_at_pixel(pol, azimuth_line, range_sample, cache_path.as_deref())
            .map_err(|e| PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(format!("{}", e)))?;
        
        Ok((position[0], position[1], position[2], velocity[0], velocity[1], velocity[2]))
    }
    
    fn extract_iw_subswaths(&mut self, polarization: &str) -> PyResult<std::collections::HashMap<String, PySubSwath>> {
        let pol = match polarization.to_uppercase().as_str() {
            "VV" => Polarization::VV,
            "VH" => Polarization::VH,
            "HV" => Polarization::HV,
            "HH" => Polarization::HH,
            _ => return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(
                format!("Invalid polarization: {}", polarization)
            )),
        };
        
        let subswaths = self.inner.extract_iw_subswaths(pol)
            .map_err(|e| PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(format!("{}", e)))?;
        
        Ok(subswaths.into_iter().map(|(k, v)| (k, PySubSwath { inner: v })).collect())
    }
    
    fn get_all_iw_subswaths(&mut self) -> PyResult<std::collections::HashMap<String, std::collections::HashMap<String, PySubSwath>>> {
        let all_subswaths = self.inner.get_all_iw_subswaths()
            .map_err(|e| PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(format!("{}", e)))?;
        
        let mut result = std::collections::HashMap::new();
        
        for (pol, subswaths) in all_subswaths {
            let pol_str = format!("{}", pol);
            let py_subswaths = subswaths.into_iter()
                .map(|(k, v)| (k, PySubSwath { inner: v }))
                .collect();
            result.insert(pol_str, py_subswaths);
        }
        
        Ok(result)
    }
    
    fn is_iw_mode(&mut self) -> PyResult<bool> {
        self.inner.is_iw_mode()
            .map_err(|e| PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(format!("{}", e)))
    }
    
    fn find_annotation_files(&mut self) -> PyResult<std::collections::HashMap<String, String>> {
        let annotations = self.inner.find_annotation_files()
            .map_err(|e| PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(format!("{}", e)))?;
        
        Ok(annotations.into_iter().map(|(k, v)| (format!("{}", k), v)).collect())
    }
    
    fn deburst_slc(&mut self, polarization: &str) -> PyResult<(Vec<Vec<(f64, f64)>>, (usize, usize))> {
        let pol = match polarization.to_uppercase().as_str() {
            "VV" => Polarization::VV,
            "VH" => Polarization::VH,
            "HV" => Polarization::HV,
            "HH" => Polarization::HH,
            _ => return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(
                format!("Invalid polarization: {}", polarization)
            )),
        };
        
        let deburst_data = self.inner.deburst_slc(pol)
            .map_err(|e| PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(format!("{}", e)))?;
        
        // Convert complex data to Python-compatible format
        let (rows, cols) = deburst_data.dim();
        let mut py_data = Vec::with_capacity(rows);
        
        for i in 0..rows {
            let mut row = Vec::with_capacity(cols);
            for j in 0..cols {
                let complex_val = deburst_data[[i, j]];
                row.push((complex_val.re as f64, complex_val.im as f64));
            }
            py_data.push(row);
        }
        
        Ok((py_data, (rows, cols)))
    }
    
    fn deburst_all_polarizations(&mut self) -> PyResult<std::collections::HashMap<String, (Vec<Vec<(f64, f64)>>, (usize, usize))>> {
        let deburst_results = self.inner.deburst_all_polarizations()
            .map_err(|e| PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(format!("{}", e)))?;
        
        let mut py_results = std::collections::HashMap::new();
        
        for (pol, deburst_data) in deburst_results {
            let pol_str = format!("{}", pol);
            let (rows, cols) = deburst_data.dim();
            let mut py_data = Vec::with_capacity(rows);
            
            for i in 0..rows {
                let mut row = Vec::with_capacity(cols);
                for j in 0..cols {
                    let complex_val = deburst_data[[i, j]];
                    row.push((complex_val.re as f64, complex_val.im as f64));
                }
                py_data.push(row);
            }
            
            py_results.insert(pol_str, (py_data, (rows, cols)));
        }
        
        Ok(py_results)
    }
    
    fn calibrate_slc(&mut self, polarization: &str, calibration_type: &str) -> PyResult<(Vec<Vec<f64>>, (usize, usize))> {
        use crate::core::calibrate::{CalibrationProcessor, CalibrationType};
        
        // Parse polarization
        let pol = match polarization.to_uppercase().as_str() {
            "VV" => Polarization::VV,
            "VH" => Polarization::VH,
            "HV" => Polarization::HV,
            "HH" => Polarization::HH,
            _ => return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(
                format!("Invalid polarization: {}", polarization)
            )),
        };
        
        // Parse calibration type
        let cal_type = match calibration_type.to_lowercase().as_str() {
            "sigma0" => CalibrationType::Sigma0,
            "beta0" => CalibrationType::Beta0,
            "gamma0" => CalibrationType::Gamma0,
            "dn" => CalibrationType::Dn,
            _ => return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(
                format!("Invalid calibration type: {}", calibration_type)
            )),
        };
        
        // Read SLC data
        let slc_data = self.inner.read_slc_data(pol)
            .map_err(|e| PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(format!("{}", e)))?;
        
        // Read calibration data
        let cal_coeffs = self.inner.read_calibration_data(pol)
            .map_err(|e| PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(format!("{}", e)))?;
        
        // Apply calibration
        let processor = CalibrationProcessor::new(cal_coeffs, cal_type);
        let calibrated_data = processor.calibrate(&slc_data)
            .map_err(|e| PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(format!("{}", e)))?;
        
        // Convert to Python format
        let (rows, cols) = calibrated_data.dim();
        let mut py_data = Vec::with_capacity(rows);
        
        for i in 0..rows {
            let mut row = Vec::with_capacity(cols);
            for j in 0..cols {
                row.push(calibrated_data[[i, j]] as f64);
            }
            py_data.push(row);
        }
        
        Ok((py_data, (rows, cols)))
    }
    
    fn get_calibration_info(&mut self, polarization: &str) -> PyResult<std::collections::HashMap<String, String>> {
        // Parse polarization
        let pol = match polarization.to_uppercase().as_str() {
            "VV" => Polarization::VV,
            "VH" => Polarization::VH,
            "HV" => Polarization::HV,
            "HH" => Polarization::HH,
            _ => return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(
                format!("Invalid polarization: {}", polarization)
            )),
        };
        
        // Read calibration data
        let cal_coeffs = self.inner.read_calibration_data(pol)
            .map_err(|e| PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(format!("{}", e)))?;
        
        let mut info = std::collections::HashMap::new();
        info.insert("swath".to_string(), cal_coeffs.swath);
        info.insert("polarization".to_string(), cal_coeffs.polarization);
        info.insert("product_first_line_utc_time".to_string(), cal_coeffs.product_first_line_utc_time);
        info.insert("product_last_line_utc_time".to_string(), cal_coeffs.product_last_line_utc_time);
        info.insert("num_vectors".to_string(), cal_coeffs.vectors.len().to_string());
        
        Ok(info)
    }

    // ...existing code...
}

/// Python wrapper for Polarization enum
#[pyclass(name = "Polarization")]
#[derive(Clone)]
struct PyPolarization {
    inner: Polarization,
}

#[pymethods]
impl PyPolarization {
    #[new]
    fn new(pol_str: String) -> PyResult<Self> {
        let pol = match pol_str.to_uppercase().as_str() {
            "VV" => Polarization::VV,
            "VH" => Polarization::VH,
            "HV" => Polarization::HV,
            "HH" => Polarization::HH,
            _ => return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(
                format!("Invalid polarization: {}", pol_str)
            )),
        };
        
        Ok(PyPolarization { inner: pol })
    }
    
    fn __str__(&self) -> String {
        format!("{}", self.inner)
    }
    
    fn __repr__(&self) -> String {
        format!("Polarization('{}')", self.inner)
    }
}

/// Python wrapper for OrbitData
#[pyclass(name = "OrbitData")]
struct PyOrbitData {
    inner: types::OrbitData,
}

#[pymethods]
impl PyOrbitData {
    #[getter]
    fn state_vectors(&self) -> Vec<PyStateVector> {
        self.inner.state_vectors.iter()
            .map(|sv| PyStateVector { inner: sv.clone() })
            .collect()
    }
    
    #[getter]
    fn reference_time(&self) -> String {
        self.inner.reference_time.to_rfc3339()
    }
    
    fn __len__(&self) -> usize {
        self.inner.state_vectors.len()
    }
    
    fn __str__(&self) -> String {
        format!("OrbitData({} state vectors)", self.inner.state_vectors.len())
    }
}

/// Python wrapper for StateVector
#[pyclass(name = "StateVector")]
struct PyStateVector {
    inner: types::StateVector,
}

#[pymethods]
impl PyStateVector {
    #[getter]
    fn time(&self) -> String {
        self.inner.time.to_rfc3339()
    }
    
    #[getter]
    fn position(&self) -> (f64, f64, f64) {
        (self.inner.position[0], self.inner.position[1], self.inner.position[2])
    }
    
    #[getter]
    fn velocity(&self) -> (f64, f64, f64) {
        (self.inner.velocity[0], self.inner.velocity[1], self.inner.velocity[2])
    }
    
    fn __str__(&self) -> String {
        format!("StateVector(time={}, position=({:.2}, {:.2}, {:.2}))", 
                self.inner.time.to_rfc3339(),
                self.inner.position[0], self.inner.position[1], self.inner.position[2])
    }
}

/// Python wrapper for SarMetadata
#[pyclass(name = "Metadata")]
struct PyMetadata {
    inner: types::SarMetadata,
}

#[pymethods]
impl PyMetadata {
    #[getter]
    fn product_id(&self) -> String {
        self.inner.product_id.clone()
    }
    
    #[getter]
    fn mission(&self) -> String {
        self.inner.mission.clone()
    }
    
    #[getter]
    fn platform(&self) -> String {
        self.inner.platform.clone()
    }
    
    #[getter]
    fn start_time(&self) -> String {
        self.inner.start_time.to_rfc3339()
    }
    
    #[getter]
    fn stop_time(&self) -> String {
        self.inner.stop_time.to_rfc3339()
    }
    
    #[getter]
    fn acquisition_mode(&self) -> String {
        format!("{:?}", self.inner.acquisition_mode)
    }
    
    #[getter]
    fn polarizations(&self) -> Vec<String> {
        self.inner.polarizations.iter()
            .map(|p| format!("{}", p))
            .collect()
    }
    
    #[getter]
    fn pixel_spacing(&self) -> (f64, f64) {
        self.inner.pixel_spacing
    }
    
    #[getter]
    fn bounding_box(&self) -> (f64, f64, f64, f64) {
        let bbox = &self.inner.bounding_box;
        (bbox.min_lon, bbox.min_lat, bbox.max_lon, bbox.max_lat)
    }
    
    #[getter]
    fn sub_swaths(&self) -> std::collections::HashMap<String, PySubSwath> {
        self.inner.sub_swaths.iter()
            .map(|(k, v)| (k.clone(), PySubSwath { inner: v.clone() }))
            .collect()
    }
    
    fn __str__(&self) -> String {
        format!(
            "SARMetadata(product_id='{}', mission='{}', mode={:?}, polarizations={:?})",
            self.inner.product_id,
            self.inner.mission,
            self.inner.acquisition_mode,
            self.inner.polarizations.iter().map(|p| format!("{}", p)).collect::<Vec<_>>()
        )
    }
}

/// Python wrapper for SubSwath
#[pyclass(name = "SubSwath")]
struct PySubSwath {
    inner: types::SubSwath,
}

#[pymethods]
impl PySubSwath {
    #[getter]
    fn id(&self) -> String {
        self.inner.id.clone()
    }
    
    #[getter]
    fn burst_count(&self) -> usize {
        self.inner.burst_count
    }
    
    #[getter]
    fn range_samples(&self) -> usize {
        self.inner.range_samples
    }
    
    #[getter]
    fn azimuth_samples(&self) -> usize {
        self.inner.azimuth_samples
    }
    
    #[getter]
    fn range_pixel_spacing(&self) -> f64 {
        self.inner.range_pixel_spacing
    }
    
    #[getter]
    fn azimuth_pixel_spacing(&self) -> f64 {
        self.inner.azimuth_pixel_spacing
    }
    
    #[getter]
    fn slant_range_time(&self) -> f64 {
        self.inner.slant_range_time
    }
    
    #[getter]
    fn burst_duration(&self) -> f64 {
        self.inner.burst_duration
    }
    
    fn __str__(&self) -> String {
        format!(
            "SubSwath(id='{}', bursts={}, range_samples={}, azimuth_samples={})",
            self.inner.id, self.inner.burst_count, self.inner.range_samples, self.inner.azimuth_samples
        )
    }
}
