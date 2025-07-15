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
