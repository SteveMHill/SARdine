//! Python bindings for download module

use crate::io::download::{
    products::{ProductMetadata, SearchParams},
    DownloadConfig, DownloadManager,
};
use crate::types::BoundingBox;
use chrono::{DateTime, Utc};
use pyo3::prelude::*;
use std::collections::HashMap;
use std::path::PathBuf;

/// Python wrapper for DownloadManager
#[pyclass(name = "DownloadManager", module = "sardine._core")]
pub struct PyDownloadManager {
    inner: DownloadManager,
}

#[pymethods]
impl PyDownloadManager {
    #[new]
    fn new(
        cache_dir: Option<String>,
        sources: Option<Vec<String>>,
        auto_download: Option<bool>,
    ) -> PyResult<Self> {
        let mut config = DownloadConfig::default();

        if let Some(dir) = cache_dir {
            config.cache_dir = PathBuf::from(dir);
        }

        if let Some(srcs) = sources {
            config.enabled_sources = srcs;
        }

        if let Some(auto) = auto_download {
            config.auto_download = auto;
        }

        let manager = DownloadManager::new(config).map_err(|e| {
            PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(format!(
                "Failed to create download manager: {}",
                e
            ))
        })?;

        Ok(Self { inner: manager })
    }

    /// Download a Sentinel-1 product
    fn download_product(&mut self, product_id: &str) -> PyResult<String> {
        let result = self.inner.download_product(product_id, None).map_err(|e| {
            PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(format!("Download failed: {}", e))
        })?;
        Ok(result.path.to_string_lossy().to_string())
    }

    /// Download orbit file for a product
    fn download_orbit(&mut self, product_id: &str, start_time: &str) -> PyResult<String> {
        let dt = DateTime::parse_from_rfc3339(start_time)
            .map_err(|e| {
                PyErr::new::<pyo3::exceptions::PyValueError, _>(format!("Invalid datetime: {}", e))
            })?
            .with_timezone(&Utc);

        let result = self
            .inner
            .download_orbit(product_id, dt, None)
            .map_err(|e| {
                PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(format!(
                    "Orbit download failed: {}",
                    e
                ))
            })?;
        Ok(result.path.to_string_lossy().to_string())
    }

    /// Download DEM tiles for a bounding box
    fn download_dem(
        &mut self,
        min_lon: f64,
        min_lat: f64,
        max_lon: f64,
        max_lat: f64,
    ) -> PyResult<Vec<String>> {
        let bbox = BoundingBox {
            min_lon,
            min_lat,
            max_lon,
            max_lat,
        };

        let results = self.inner.download_dem(&bbox, None).map_err(|e| {
            PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(format!("DEM download failed: {}", e))
        })?;

        Ok(results
            .into_iter()
            .map(|r| r.path.to_string_lossy().to_string())
            .collect())
    }

    /// Check if auto-download is enabled
    fn auto_download(&self) -> bool {
        self.inner.auto_download()
    }

    /// Search for Sentinel-1 products
    fn search_products(
        &self,
        start_date: &str,
        end_date: &str,
        product_type: Option<String>,
        platform: Option<String>,
        max_results: Option<usize>,
    ) -> PyResult<Vec<PyProductMetadata>> {
        let start_dt = DateTime::parse_from_rfc3339(start_date)
            .map_err(|e| {
                PyErr::new::<pyo3::exceptions::PyValueError, _>(format!(
                    "Invalid start_date: {}",
                    e
                ))
            })?
            .with_timezone(&Utc);

        let end_dt = DateTime::parse_from_rfc3339(end_date)
            .map_err(|e| {
                PyErr::new::<pyo3::exceptions::PyValueError, _>(format!("Invalid end_date: {}", e))
            })?
            .with_timezone(&Utc);

        let params = SearchParams {
            product_type,
            platform,
            start_date: start_dt,
            end_date: end_dt,
            aoi_wkt: None,
            max_results: max_results.unwrap_or(100),
            acquisition_mode: None,
            polarization: None,
            orbit_direction: None,
            relative_orbit: None,
        };

        let results = self.inner.search_products(&params).map_err(|e| {
            PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(format!("Search failed: {}", e))
        })?;

        Ok(results
            .into_iter()
            .map(|m| PyProductMetadata::from(m))
            .collect())
    }

    /// Download multiple products in batch
    fn download_products_batch(
        &mut self,
        product_ids: Vec<String>,
    ) -> PyResult<HashMap<String, String>> {
        let results = self
            .inner
            .download_products_batch(&product_ids)
            .map_err(|e| {
                PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(format!(
                    "Batch download failed: {}",
                    e
                ))
            })?;

        let mut result_map = HashMap::new();
        for (product_id, download_result) in results {
            result_map.insert(
                product_id,
                download_result.path.to_string_lossy().to_string(),
            );
        }

        Ok(result_map)
    }
}

/// Python wrapper for ProductMetadata
#[pyclass(name = "ProductMetadata")]
pub struct PyProductMetadata {
    #[pyo3(get)]
    pub id: String,
    #[pyo3(get)]
    pub title: String,
    #[pyo3(get)]
    pub product_type: String,
    #[pyo3(get)]
    pub size_mb: f64,
    #[pyo3(get)]
    pub download_url: Option<String>,
    #[pyo3(get)]
    pub checksum_md5: Option<String>,
}

impl From<ProductMetadata> for PyProductMetadata {
    fn from(m: ProductMetadata) -> Self {
        Self {
            id: m.id,
            title: m.title,
            product_type: m.product_type,
            size_mb: m.size_mb,
            download_url: m.download_url,
            checksum_md5: m.checksum_md5,
        }
    }
}

/// Initialize download module for Python
pub fn init_module(m: &PyModule) -> PyResult<()> {
    m.add_class::<PyDownloadManager>()?;
    m.add_class::<PyProductMetadata>()?;
    Ok(())
}
