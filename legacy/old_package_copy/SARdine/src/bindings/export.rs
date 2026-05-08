//! Export bindings for Python
//!
//! Contains GeoTIFF export and metadata generation functions.

use numpy::ToPyArray;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::{PyDict, PyTuple};

/// Step 13: Export GeoTIFF with proper georeferencing
#[pyfunction]
pub fn export_geotiff(
    py: Python,
    data: numpy::PyReadonlyArray2<f32>,
    output_path: String,
    geo_transform: Vec<f64>, // [x_origin, pixel_width, 0, y_origin, 0, -pixel_height]
    crs_epsg: i32,
    metadata: Option<std::collections::HashMap<String, String>>,
) -> PyResult<PyObject> {
    // Validate geo_transform
    if geo_transform.len() != 6 {
        return Err(PyValueError::new_err(
            "geo_transform must have exactly 6 elements",
        ));
    }

    // Call the real Python exporter (rasterio-based) to write a COG/GeoTIFF
    let result = PyDict::new(py);

    // Prepare inputs
    let np_array = data.as_array().to_pyarray(py);
    let gt_tuple = PyTuple::new(py, &geo_transform);
    let crs_str = format!("EPSG:{}", crs_epsg);

    let py_export = pyo3::types::PyModule::import(py, "sardine.export").map_err(|e| {
        PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(format!(
            "Failed to import sardine.export module: {:?}",
            e
        ))
    })?;

    let kwargs = PyDict::new(py);
    kwargs.set_item("data", np_array)?;
    kwargs.set_item("output_path", &output_path)?;
    kwargs.set_item("geotransform", gt_tuple)?;
    kwargs.set_item("crs", crs_str)?;

    // Optional metadata
    if let Some(meta) = metadata {
        let meta_dict = PyDict::new(py);
        for (k, v) in meta {
            meta_dict.set_item(k, v)?;
        }
        kwargs.set_item("metadata", meta_dict)?;
    }

    match py_export
        .getattr("export_to_geotiff")
        .and_then(|f| f.call((), Some(kwargs)))
    {
        Ok(_py_path) => {
            // Build a small result dict mirroring the call
            result.set_item("status", "success")?;
            result.set_item("message", "GeoTIFF exported via rasterio")?;
            result.set_item("output_path", &output_path)?;
            result.set_item("crs_epsg", crs_epsg)?;
            result.set_item("geo_transform", &geo_transform)?;
        }
        Err(e) => {
            return Err(PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(format!(
                "GeoTIFF export failed: {:?}",
                e
            )));
        }
    }

    Ok(result.into())
}

/// Export a Cloud-Optimized GeoTIFF (COG) and STAC metadata via Python exporter
#[pyfunction]
pub fn export_cog_with_stac(
    py: Python,
    data: numpy::PyReadonlyArray2<f32>,
    output_dir: String,
    filename_base: String,
    geo_transform: Vec<f64>,
    stac_metadata: std::collections::HashMap<String, PyObject>,
    crs_epsg: i32,
    compress: Option<String>,
) -> PyResult<PyObject> {
    if geo_transform.len() != 6 {
        return Err(PyValueError::new_err(
            "geo_transform must have exactly 6 elements",
        ));
    }

    let py_export = pyo3::types::PyModule::import(py, "sardine.export").map_err(|e| {
        PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(format!(
            "Failed to import sardine.export: {:?}",
            e
        ))
    })?;

    let np_array = data.as_array().to_pyarray(py);
    let gt_tuple = PyTuple::new(py, &geo_transform);
    let crs_str = format!("EPSG:{}", crs_epsg);

    // SCIENTIFIC REQUIREMENT: Compression method must be explicitly specified
    let comp = match compress {
        Some(method) => method,
        None => return Err(PyValueError::new_err(
            "compress parameter is required; no fallback values permitted for scientific accuracy",
        )),
    };

    // Build kwargs for create_cog_with_stac
    let kwargs = PyDict::new(py);
    kwargs.set_item("data", np_array)?;
    kwargs.set_item("output_dir", &output_dir)?;
    kwargs.set_item("filename_base", &filename_base)?;
    kwargs.set_item("geotransform", gt_tuple)?;
    kwargs.set_item("sar_metadata", stac_metadata)?;
    kwargs.set_item("crs", crs_str)?;
    kwargs.set_item("compress", comp)?;

    let func = py_export.getattr("create_cog_with_stac")?;
    let ret = func.call((), Some(kwargs)).map_err(|e| {
        PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(format!(
            "create_cog_with_stac failed: {:?}",
            e
        ))
    })?;

    Ok(ret.into())
}

/// Step 14: Export metadata as JSON
#[pyfunction]
pub fn export_metadata_json(
    metadata: std::collections::HashMap<String, String>,
) -> PyResult<String> {
    use serde_json;

    serde_json::to_string_pretty(&metadata)
        .map_err(|e| PyValueError::new_err(format!("JSON serialization failed: {}", e)))
}

/// Step 14: Export metadata as XML
#[pyfunction]
pub fn export_metadata_xml(
    metadata: std::collections::HashMap<String, String>,
) -> PyResult<String> {
    let mut xml = String::from(
        r#"<?xml version="1.0" encoding="UTF-8"?>
<sar_processing_metadata>
"#,
    );

    for (key, value) in metadata {
        // Escape XML special characters
        let escaped_value = value
            .replace("&", "&amp;")
            .replace("<", "&lt;")
            .replace(">", "&gt;")
            .replace("\"", "&quot;")
            .replace("'", "&apos;");

        xml.push_str(&format!("  <{}>{}</{}>\n", key, escaped_value, key));
    }

    xml.push_str("</sar_processing_metadata>");

    Ok(xml)
}
