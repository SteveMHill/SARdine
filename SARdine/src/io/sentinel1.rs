use crate::types::{SarError, SarResult};
use chrono::{DateTime, NaiveDateTime, Utc};
use log::{debug, info, warn};
use reqwest::blocking::Client;
use serde::{Deserialize, Serialize};
use std::path::{Path, PathBuf};
use std::time::Duration;

/// Sentinel-1 data provider configuration
#[derive(Debug, Clone)]
pub enum DataProvider {
    /// ESA Copernicus Open Access Hub
    CopernicusHub {
        username: String,
        password: String,
        base_url: String,
    },
    /// Alaska Satellite Facility DAAC
    ASF {
        username: Option<String>,
        password: Option<String>,
        base_url: String,
    },
    /// Generic provider with custom configuration
    Custom {
        name: String,
        base_url: String,
        auth_method: AuthMethod,
    },
}

/// Authentication methods for different providers
#[derive(Debug, Clone)]
pub enum AuthMethod {
    Basic { username: String, password: String },
    Bearer { token: String },
    None,
}

/// Sentinel-1 product metadata for search and download
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Sentinel1Product {
    pub id: String,
    pub title: String,
    pub platform_name: String,     // S1A or S1B
    pub product_type: String,      // SLC, GRD, etc.
    pub acquisition_mode: String,  // IW, EW, SM, WV
    pub polarization: Vec<String>, // VV, VH, HH, HV
    pub start_time: DateTime<Utc>,
    pub end_time: DateTime<Utc>,
    pub orbit_number: u32,
    pub relative_orbit_number: u32,
    pub orbit_direction: String, // ASCENDING or DESCENDING
    pub footprint: String,       // WKT geometry
    pub size_mb: f64,
    pub download_url: String,
    pub checksum_md5: Option<String>,
}

/// Search parameters for Sentinel-1 products
#[derive(Debug, Clone)]
pub struct SearchParams {
    pub platform: Option<String>,         // S1A, S1B, or both
    pub product_type: Option<String>,     // SLC, GRD, OCN
    pub acquisition_mode: Option<String>, // IW, EW, SM, WV
    pub polarization: Option<String>,     // VV, VH, HH, HV
    pub start_date: DateTime<Utc>,
    pub end_date: DateTime<Utc>,
    pub aoi_wkt: Option<String>,         // Area of interest in WKT format
    pub orbit_direction: Option<String>, // ASCENDING, DESCENDING
    pub relative_orbit: Option<u32>,
    pub max_results: usize,
}

/// Download progress callback
pub type ProgressCallback = Box<dyn Fn(u64, u64) -> () + Send>;

/// Sentinel-1 downloader with support for multiple providers
pub struct Sentinel1Downloader {
    providers: Vec<DataProvider>,
    client: Client,
    output_dir: PathBuf,
    verify_checksums: bool,
}

impl Sentinel1Downloader {
    /// Create a new Sentinel-1 downloader
    pub fn new(output_dir: &Path) -> SarResult<Self> {
        let client = Client::builder()
            .timeout(Duration::from_secs(3600)) // 1 hour for large files
            .user_agent("SARdine/0.2.1 (Sentinel-1 SAR Processing)")
            .build()
            .map_err(|e| SarError::Processing(format!("Failed to create HTTP client: {}", e)))?;

        Ok(Self {
            providers: Vec::new(),
            client,
            output_dir: output_dir.to_path_buf(),
            verify_checksums: true,
        })
    }

    /// Add a data provider
    pub fn add_provider(&mut self, provider: DataProvider) {
        self.providers.push(provider);
    }

    /// Add ESA Copernicus Hub provider
    pub fn add_copernicus_hub(&mut self, username: String, password: String) {
        self.add_provider(DataProvider::CopernicusHub {
            username,
            password,
            base_url: "https://apihub.copernicus.eu".to_string(),
        });
    }

    /// Add Alaska Satellite Facility provider
    pub fn add_asf_provider(&mut self, username: Option<String>, password: Option<String>) {
        self.add_provider(DataProvider::ASF {
            username,
            password,
            base_url: "https://api.daac.asf.alaska.edu".to_string(),
        });
    }

    /// Search for Sentinel-1 products
    pub fn search(&self, params: &SearchParams) -> SarResult<Vec<Sentinel1Product>> {
        let mut all_products = Vec::new();
        let mut search_errors = Vec::new();

        for provider in &self.providers {
            match self.search_provider(provider, params) {
                Ok(mut products) => {
                    info!(
                        "Found {} products from provider {:?}",
                        products.len(),
                        provider
                    );
                    all_products.append(&mut products);
                }
                Err(e) => {
                    warn!("Search failed for provider {:?}: {}", provider, e);
                    search_errors.push(e);
                }
            }
        }

        if all_products.is_empty() && !search_errors.is_empty() {
            return Err(SarError::Processing(format!(
                "All provider searches failed. Last error: {}",
                search_errors.last().unwrap()
            )));
        }

        // Remove duplicates based on product title
        all_products.sort_by(|a, b| a.title.cmp(&b.title));
        all_products.dedup_by(|a, b| a.title == b.title);

        // Limit results
        all_products.truncate(params.max_results);

        Ok(all_products)
    }

    /// Download a Sentinel-1 product
    pub fn download_product(
        &self,
        product: &Sentinel1Product,
        progress_callback: Option<ProgressCallback>,
    ) -> SarResult<PathBuf> {
        info!("Starting download of product: {}", product.title);
        info!("Product size: {:.1} MB", product.size_mb);

        let output_path = self.output_dir.join(format!("{}.zip", product.title));

        // Check if file already exists and is complete
        if self.is_download_complete(&output_path, product)? {
            info!(
                "Product already downloaded and verified: {}",
                output_path.display()
            );
            return Ok(output_path);
        }

        // Try downloading from each provider
        for provider in &self.providers {
            match self.download_from_provider(
                provider,
                product,
                &output_path,
                progress_callback.as_ref(),
            ) {
                Ok(path) => {
                    info!("Successfully downloaded from provider {:?}", provider);

                    if self.verify_checksums {
                        self.verify_download(&path, product)?;
                    }

                    return Ok(path);
                }
                Err(e) => {
                    warn!("Download failed from provider {:?}: {}", provider, e);
                    // Remove partial file
                    if output_path.exists() {
                        let _ = std::fs::remove_file(&output_path);
                    }
                }
            }
        }

        Err(SarError::Processing(format!(
            "Failed to download product {} from all providers",
            product.title
        )))
    }

    /// Download multiple products with parallel downloads
    pub fn download_products(
        &self,
        products: &[Sentinel1Product],
        max_parallel: usize,
        progress_callback: Option<ProgressCallback>,
    ) -> SarResult<Vec<PathBuf>> {
        use rayon::prelude::*;
        use std::sync::Arc;

        let callback = Arc::new(progress_callback);
        let downloader = Arc::new(self);

        let results: Vec<SarResult<PathBuf>> = products
            .par_iter()
            .with_max_len(max_parallel)
            .map(|product| downloader.download_product(product, None))
            .collect();

        let mut success_paths = Vec::new();
        let mut errors = Vec::new();

        for result in results {
            match result {
                Ok(path) => success_paths.push(path),
                Err(e) => errors.push(e),
            }
        }

        if !errors.is_empty() {
            warn!("Some downloads failed: {:?}", errors);
        }

        if success_paths.is_empty() {
            return Err(SarError::Processing("All downloads failed".to_string()));
        }

        Ok(success_paths)
    }

    /// Search products from a specific provider
    fn search_provider(
        &self,
        provider: &DataProvider,
        params: &SearchParams,
    ) -> SarResult<Vec<Sentinel1Product>> {
        match provider {
            DataProvider::CopernicusHub {
                username,
                password,
                base_url,
            } => self.search_copernicus_hub(base_url, username, password, params),
            DataProvider::ASF {
                username,
                password,
                base_url,
            } => self.search_asf(base_url, username.as_deref(), password.as_deref(), params),
            DataProvider::Custom {
                base_url,
                auth_method,
                ..
            } => self.search_custom_provider(base_url, auth_method, params),
        }
    }

    /// Search Copernicus Hub using OpenSearch API
    fn search_copernicus_hub(
        &self,
        base_url: &str,
        username: &str,
        password: &str,
        params: &SearchParams,
    ) -> SarResult<Vec<Sentinel1Product>> {
        let search_url = format!("{}/dhus/search", base_url);

        // Build OpenSearch query
        let mut query_parts = vec!["platformname:Sentinel-1".to_string()];

        if let Some(platform) = &params.platform {
            query_parts.push(format!("platformserialidentifier:{}", platform));
        }

        if let Some(product_type) = &params.product_type {
            query_parts.push(format!("producttype:{}", product_type));
        }

        if let Some(mode) = &params.acquisition_mode {
            query_parts.push(format!("sensoroperationalmode:{}", mode));
        }

        if let Some(polarization) = &params.polarization {
            query_parts.push(format!("polarisationmode:{}", polarization));
        }

        if let Some(orbit_dir) = &params.orbit_direction {
            query_parts.push(format!("orbitdirection:{}", orbit_dir));
        }

        if let Some(rel_orbit) = params.relative_orbit {
            query_parts.push(format!("relativeorbitnumber:{}", rel_orbit));
        }

        // Add time range
        let start_time = params.start_date.format("%Y-%m-%dT%H:%M:%S.%3fZ");
        let end_time = params.end_date.format("%Y-%m-%dT%H:%M:%S.%3fZ");
        query_parts.push(format!("beginposition:[{} TO {}]", start_time, end_time));

        // Add AOI if provided
        if let Some(aoi) = &params.aoi_wkt {
            query_parts.push(format!("footprint:\"Intersects({})\"", aoi));
        }

        let query = query_parts.join(" AND ");

        let request = self
            .client
            .get(&search_url)
            .basic_auth(username, Some(password))
            .query(&[
                ("q", query.as_str()),
                ("format", "json"),
                ("rows", &params.max_results.to_string()),
                ("orderby", "beginposition desc"),
            ]);

        debug!("Performing Copernicus Hub search");

        let response = request
            .send()
            .map_err(|e| SarError::Processing(format!("Copernicus Hub request failed: {}", e)))?;

        if !response.status().is_success() {
            return Err(SarError::Processing(format!(
                "Copernicus Hub search failed with status: {}",
                response.status()
            )));
        }

        let response_text = response
            .text()
            .map_err(|e| SarError::Processing(format!("Failed to read response: {}", e)))?;

        // Parse Copernicus Hub response (simplified)
        self.parse_copernicus_response(&response_text)
    }

    /// Search Alaska Satellite Facility
    fn search_asf(
        &self,
        base_url: &str,
        username: Option<&str>,
        password: Option<&str>,
        params: &SearchParams,
    ) -> SarResult<Vec<Sentinel1Product>> {
        let search_url = format!("{}/services/search/param", base_url);

        let mut query_params = vec![
            ("platform", "SENTINEL-1A,SENTINEL-1B".to_string()),
            ("output", "JSON".to_string()),
            ("maxResults", params.max_results.to_string()),
        ];

        if let Some(product_type) = &params.product_type {
            query_params.push(("processingLevel", product_type.clone()));
        }

        if let Some(mode) = &params.acquisition_mode {
            query_params.push(("beamMode", mode.clone()));
        }

        if let Some(polarization) = &params.polarization {
            query_params.push(("polarization", polarization.clone()));
        }

        if let Some(orbit_dir) = &params.orbit_direction {
            query_params.push(("flightDirection", orbit_dir.clone()));
        }

        // Add time range
        let start_time = params.start_date.format("%Y-%m-%dT%H:%M:%SZ");
        let end_time = params.end_date.format("%Y-%m-%dT%H:%M:%SZ");
        query_params.push(("start", start_time.to_string()));
        query_params.push(("end", end_time.to_string()));

        let mut request = self.client.get(&search_url);

        // Add authentication if provided
        if let (Some(user), Some(pass)) = (username, password) {
            request = request.basic_auth(user, Some(pass));
        }

        let response = request
            .form(&query_params)
            .send()
            .map_err(|e| SarError::Processing(format!("ASF search request failed: {}", e)))?;

        if !response.status().is_success() {
            return Err(SarError::Processing(format!(
                "ASF search failed with status: {}",
                response.status()
            )));
        }

        let response_text = response
            .text()
            .map_err(|e| SarError::Processing(format!("Failed to read ASF response: {}", e)))?;

        // Parse ASF response
        self.parse_asf_response(&response_text)
    }

    /// Search custom provider (placeholder)
    fn search_custom_provider(
        &self,
        _base_url: &str,
        _auth_method: &AuthMethod,
        _params: &SearchParams,
    ) -> SarResult<Vec<Sentinel1Product>> {
        // Implement custom provider search logic here
        Ok(Vec::new())
    }

    /// Download from a specific provider
    fn download_from_provider(
        &self,
        provider: &DataProvider,
        product: &Sentinel1Product,
        output_path: &Path,
        progress_callback: Option<&ProgressCallback>,
    ) -> SarResult<PathBuf> {
        match provider {
            DataProvider::CopernicusHub {
                username, password, ..
            } => self.download_from_copernicus_hub(
                &product.download_url,
                username,
                password,
                output_path,
                progress_callback,
            ),
            DataProvider::ASF {
                username, password, ..
            } => self.download_from_asf(
                &product.download_url,
                username.as_deref(),
                password.as_deref(),
                output_path,
                progress_callback,
            ),
            DataProvider::Custom { auth_method, .. } => self.download_from_custom_provider(
                &product.download_url,
                auth_method,
                output_path,
                progress_callback,
            ),
        }
    }

    /// Download from Copernicus Hub
    fn download_from_copernicus_hub(
        &self,
        url: &str,
        username: &str,
        password: &str,
        output_path: &Path,
        progress_callback: Option<&ProgressCallback>,
    ) -> SarResult<PathBuf> {
        self.download_with_auth(
            url,
            Some((username, password)),
            output_path,
            progress_callback,
        )
    }

    /// Download from ASF
    fn download_from_asf(
        &self,
        url: &str,
        username: Option<&str>,
        password: Option<&str>,
        output_path: &Path,
        progress_callback: Option<&ProgressCallback>,
    ) -> SarResult<PathBuf> {
        let auth = username.zip(password);
        self.download_with_auth(url, auth, output_path, progress_callback)
    }

    /// Download from custom provider
    fn download_from_custom_provider(
        &self,
        url: &str,
        auth_method: &AuthMethod,
        output_path: &Path,
        progress_callback: Option<&ProgressCallback>,
    ) -> SarResult<PathBuf> {
        match auth_method {
            AuthMethod::Basic { username, password } => self.download_with_auth(
                url,
                Some((username, password)),
                output_path,
                progress_callback,
            ),
            AuthMethod::Bearer { token: _ } => {
                // TODO: Implement bearer token authentication
                self.download_with_auth(url, None, output_path, progress_callback)
            }
            AuthMethod::None => self.download_with_auth(url, None, output_path, progress_callback),
        }
    }

    /// Download file with optional authentication and progress reporting
    fn download_with_auth(
        &self,
        url: &str,
        auth: Option<(&str, &str)>,
        output_path: &Path,
        progress_callback: Option<&ProgressCallback>,
    ) -> SarResult<PathBuf> {
        debug!("Downloading from URL: {}", url);

        // Create output directory
        if let Some(parent) = output_path.parent() {
            std::fs::create_dir_all(parent)
                .map_err(|e| SarError::Processing(format!("Failed to create directory: {}", e)))?;
        }

        let mut request = self.client.get(url);

        // Add authentication if provided
        if let Some((username, password)) = auth {
            request = request.basic_auth(username, Some(password));
        }

        // Send request and get response
        let response = request
            .send()
            .map_err(|e| SarError::Processing(format!("Download request failed: {}", e)))?;

        if !response.status().is_success() {
            return Err(SarError::Processing(format!(
                "Download failed with status: {} for URL: {}",
                response.status(),
                url
            )));
        }

        // Get content length for progress reporting
        let total_size = response.content_length().unwrap_or(0);

        // Download with progress tracking
        let mut downloaded = 0u64;
        let mut output_file = std::fs::File::create(output_path)
            .map_err(|e| SarError::Processing(format!("Failed to create output file: {}", e)))?;

        use std::io::{Read, Write};
        let mut reader = response;
        let mut buffer = vec![0u8; 8192];

        loop {
            let bytes_read = reader
                .read(&mut buffer)
                .map_err(|e| SarError::Processing(format!("Failed to read chunk: {}", e)))?;

            if bytes_read == 0 {
                break;
            }

            output_file
                .write_all(&buffer[..bytes_read])
                .map_err(|e| SarError::Processing(format!("Failed to write chunk: {}", e)))?;

            downloaded += bytes_read as u64;

            // Report progress
            if let Some(callback) = progress_callback {
                callback(downloaded, total_size);
            }
        }

        info!(
            "Successfully downloaded {} bytes to {}",
            downloaded,
            output_path.display()
        );
        Ok(output_path.to_path_buf())
    }

    /// Check if download is complete and valid
    fn is_download_complete(&self, path: &Path, product: &Sentinel1Product) -> SarResult<bool> {
        if !path.exists() {
            return Ok(false);
        }

        let metadata = std::fs::metadata(path)
            .map_err(|e| SarError::Processing(format!("Failed to read file metadata: {}", e)))?;

        let file_size_mb = metadata.len() as f64 / (1024.0 * 1024.0);

        // Check if file size matches expected size (with 1% tolerance)
        let size_diff = (file_size_mb - product.size_mb).abs();
        let tolerance = product.size_mb * 0.01; // 1% tolerance

        if size_diff > tolerance {
            warn!(
                "File size mismatch: expected {:.1} MB, found {:.1} MB",
                product.size_mb, file_size_mb
            );
            return Ok(false);
        }

        Ok(true)
    }

    /// Verify download using checksum
    fn verify_download(&self, path: &Path, product: &Sentinel1Product) -> SarResult<()> {
        if let Some(expected_md5) = &product.checksum_md5 {
            debug!("Verifying MD5 checksum for {}", path.display());

            let computed_md5 = self.compute_md5_checksum(path)?;

            if computed_md5.to_lowercase() != expected_md5.to_lowercase() {
                return Err(SarError::Processing(format!(
                    "Checksum verification failed. Expected: {}, Got: {}",
                    expected_md5, computed_md5
                )));
            }

            info!("Checksum verification passed for {}", path.display());
        }

        Ok(())
    }

    /// Compute MD5 checksum of a file
    fn compute_md5_checksum(&self, path: &Path) -> SarResult<String> {
        use std::io::Read;

        let mut file = std::fs::File::open(path).map_err(|e| {
            SarError::Processing(format!("Failed to open file for checksum: {}", e))
        })?;

        let mut hasher = md5::Context::new();
        let mut buffer = [0; 8192];

        loop {
            let bytes_read = file.read(&mut buffer).map_err(|e| {
                SarError::Processing(format!("Failed to read file for checksum: {}", e))
            })?;

            if bytes_read == 0 {
                break;
            }

            hasher.consume(&buffer[..bytes_read]);
        }

        let digest = hasher.compute();
        Ok(format!("{:x}", digest))
    }

    /// Parse Copernicus Hub JSON response
    fn parse_copernicus_response(&self, response: &str) -> SarResult<Vec<Sentinel1Product>> {
        // This is a simplified parser - in practice, you'd use proper JSON parsing
        // with serde_json and define the exact response structure from Copernicus Hub
        debug!("Parsing Copernicus Hub response (simplified)");

        // TODO: Implement proper JSON parsing for Copernicus Hub API response
        // The actual structure depends on the API version and response format

        Ok(Vec::new()) // Placeholder
    }

    /// Parse ASF JSON response
    fn parse_asf_response(&self, response: &str) -> SarResult<Vec<Sentinel1Product>> {
        // This is a simplified parser - in practice, you'd use proper JSON parsing
        debug!("Parsing ASF response (simplified)");

        // TODO: Implement proper JSON parsing for ASF API response

        Ok(Vec::new()) // Placeholder
    }
}

/// Convenience functions for common use cases

/// Download Sentinel-1 products by product names
pub fn download_by_product_names(
    product_names: &[String],
    output_dir: &Path,
    copernicus_credentials: Option<(String, String)>,
    asf_credentials: Option<(String, String)>,
) -> SarResult<Vec<PathBuf>> {
    let mut downloader = Sentinel1Downloader::new(output_dir)?;

    // Add providers
    if let Some((username, password)) = copernicus_credentials {
        downloader.add_copernicus_hub(username, password);
    }

    if let Some((username, password)) = asf_credentials {
        downloader.add_asf_provider(Some(username), Some(password));
    } else {
        // Add ASF without credentials for public data
        downloader.add_asf_provider(None, None);
    }

    // Create search parameters for each product
    let mut downloaded_paths = Vec::new();

    for product_name in product_names {
        // Parse product name to extract search parameters
        let search_params = parse_product_name_to_search_params(product_name)?;

        // Search for the specific product
        let products = downloader.search(&search_params)?;

        // Find exact match
        let product = products
            .into_iter()
            .find(|p| p.title == *product_name)
            .ok_or_else(|| SarError::Processing(format!("Product not found: {}", product_name)))?;

        // Download the product
        let path = downloader.download_product(&product, None)?;
        downloaded_paths.push(path);
    }

    Ok(downloaded_paths)
}

/// Parse Sentinel-1 product name to extract search parameters
fn parse_product_name_to_search_params(product_name: &str) -> SarResult<SearchParams> {
    // Example: S1A_IW_SLC__1SDV_20201230T165244_20201230T165311_035918_0434F0_6788

    let parts: Vec<&str> = product_name.split('_').collect();
    if parts.len() < 7 {
        return Err(SarError::Processing(format!(
            "Invalid product name format: {}",
            product_name
        )));
    }

    let platform = parts[0].to_string(); // S1A or S1B
    let mode = parts[1].to_string(); // IW, EW, etc.
    let product_type = parts[2].to_string(); // SLC, GRD, etc.

    // Parse start time from product name
    let start_time_str = parts[5];
    let start_time = NaiveDateTime::parse_from_str(start_time_str, "%Y%m%dT%H%M%S")
        .map_err(|e| SarError::Processing(format!("Failed to parse start time: {}", e)))?
        .and_utc();

    // Use a narrow time window around the product
    let end_time = start_time + chrono::Duration::minutes(30);

    Ok(SearchParams {
        platform: Some(platform),
        product_type: Some(product_type),
        acquisition_mode: Some(mode),
        polarization: None,
        start_date: start_time,
        end_date: end_time,
        aoi_wkt: None,
        orbit_direction: None,
        relative_orbit: None,
        max_results: 10,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use tempfile::TempDir;

    #[test]
    fn test_downloader_creation() {
        let temp_dir = TempDir::new().unwrap();
        let downloader = Sentinel1Downloader::new(temp_dir.path());
        assert!(downloader.is_ok());
    }

    #[test]
    fn test_product_name_parsing() {
        let product_name = "S1A_IW_SLC__1SDV_20201230T165244_20201230T165311_035918_0434F0_6788";
        let params = parse_product_name_to_search_params(product_name).unwrap();

        assert_eq!(params.platform, Some("S1A".to_string()));
        assert_eq!(params.acquisition_mode, Some("IW".to_string()));
        assert_eq!(params.product_type, Some("SLC".to_string()));
    }

    #[test]
    fn test_add_providers() {
        let temp_dir = TempDir::new().unwrap();
        let mut downloader = Sentinel1Downloader::new(temp_dir.path()).unwrap();

        downloader.add_copernicus_hub("user".to_string(), "pass".to_string());
        downloader.add_asf_provider(Some("user".to_string()), Some("pass".to_string()));

        assert_eq!(downloader.providers.len(), 2);
    }
}
