//! Sentinel-1 product downloader

use crate::types::{SarError, SarResult};
use chrono::{DateTime, Utc};
use log::{debug, info, warn};
use std::fs;
use std::sync::Arc;
use std::time::SystemTime;

use super::cache::{CacheEntry, CacheManager};
use super::manager::{DownloadConfig, DownloadResult};
use super::progress::ProgressCallback;
use super::providers::ProviderRegistry;

// Constants
const MIN_PRODUCT_ID_LENGTH: usize = 50;

/// Sanitize product ID to prevent path traversal attacks
fn sanitize_product_id(product_id: &str) -> String {
    product_id
        .split('/')
        .last()
        .unwrap_or(product_id)
        .replace("..", "")
        .replace("/", "_")
        .replace("\\", "_")
        .replace("\0", "") // Remove null bytes
        .chars()
        .take(200) // Limit length
        .collect()
}

/// Extract date from product ID (format: S1A_IW_SLC__1SDV_20201230T165244_...)
/// FIXED: Added input validation
fn extract_date_from_product_id(product_id: &str) -> Option<DateTime<Utc>> {
    // Validate minimum length
    if product_id.len() < MIN_PRODUCT_ID_LENGTH {
        return None;
    }

    // Product ID format: S1A_IW_SLC__1SDV_20201230T165244_20201230T165311_...
    // Extract the first timestamp
    let parts: Vec<&str> = product_id.split('_').collect();
    for part in parts {
        if part.len() == 15 && part.chars().nth(8) == Some('T') {
            // Validate format: YYYYMMDDTHHMMSS (all digits except T)
            if part.chars().take(8).all(|c| c.is_ascii_digit())
                && part.chars().skip(9).all(|c| c.is_ascii_digit())
            {
                // Format: 20201230T165244
                if let Ok(dt) = DateTime::parse_from_str(&format!("{}Z", part), "%Y%m%dT%H%M%SZ") {
                    return Some(dt.with_timezone(&Utc));
                }
            }
        }
    }
    None
}

/// Extract product type from product ID
fn extract_product_type_from_id(product_id: &str) -> Option<String> {
    if product_id.contains("_SLC_") {
        Some("SLC".to_string())
    } else if product_id.contains("_GRD_") {
        Some("GRD".to_string())
    } else {
        None
    }
}

/// Extract platform from product ID
fn extract_platform_from_id(product_id: &str) -> Option<String> {
    if product_id.starts_with("S1A") {
        Some("S1A".to_string())
    } else if product_id.starts_with("S1B") {
        Some("S1B".to_string())
    } else {
        None
    }
}

/// Product metadata for search and download
#[derive(Debug, Clone)]
pub struct ProductMetadata {
    pub id: String,
    pub title: String,
    pub product_type: String, // SLC, GRD, etc.
    pub size_mb: f64,
    pub download_url: Option<String>,
    pub checksum_md5: Option<String>,
}

/// Search parameters for products
#[derive(Debug, Clone)]
pub struct SearchParams {
    pub product_type: Option<String>, // SLC, GRD
    pub platform: Option<String>,     // S1A, S1B
    pub start_date: DateTime<Utc>,
    pub end_date: DateTime<Utc>,
    pub aoi_wkt: Option<String>, // Area of interest in WKT format
    pub max_results: usize,
    pub acquisition_mode: Option<String>, // IW, EW, SM, WV
    pub polarization: Option<String>,     // VV, VH, HH, HV
    pub orbit_direction: Option<String>,  // ASCENDING, DESCENDING
    pub relative_orbit: Option<u32>,
}

/// Product downloader
pub struct ProductDownloader {
    cache: Arc<CacheManager>,
    providers: Arc<ProviderRegistry>,
    config: DownloadConfig,
}

impl ProductDownloader {
    pub fn new(
        cache: &Arc<CacheManager>,
        providers: &Arc<ProviderRegistry>,
        config: &DownloadConfig,
    ) -> SarResult<Self> {
        Ok(Self {
            cache: cache.clone(),
            providers: providers.clone(),
            config: config.clone(),
        })
    }

    /// Download a product by ID or URL
    pub fn download(
        &self,
        product_id_or_url: &str,
        progress: Option<ProgressCallback>,
    ) -> SarResult<DownloadResult> {
        info!("Downloading product: {}", product_id_or_url);

        // Check cache first
        if let Some(cached) = self.cache.get_product(product_id_or_url) {
            if cached.path.exists() {
                let metadata = fs::metadata(&cached.path).map_err(SarError::Io)?;
                info!("Using cached product: {}", cached.path.display());
                return Ok(DownloadResult {
                    path: cached.path,
                    size: metadata.len(),
                    source: "cache".to_string(),
                });
            }
        }

        // Determine if input is a URL or product ID
        let mut expected_md5: Option<String> = None;
        let download_urls: Vec<String> = if product_id_or_url.starts_with("http://")
            || product_id_or_url.starts_with("https://")
        {
            // Direct URL provided
            vec![product_id_or_url.to_string()]
        } else {
            // Try to resolve product ID to URL(s) using providers
            debug!("Resolving product ID to URL: {}", product_id_or_url);
            let mut resolved_urls = Vec::new();

            for provider_name in self.providers.enabled() {
                if let Some(provider) = self.providers.get_provider(provider_name) {
                    match provider.resolve_product_url(product_id_or_url) {
                        Ok(Some(url)) => {
                            let url_clone = url.clone();
                            resolved_urls.push(url);
                            debug!("Resolved URL via {}: {}", provider.name(), url_clone);
                        }
                        Ok(None) => {
                            // Provider doesn't support direct resolution, try search as fallback
                            debug!(
                                "Provider {} doesn't support direct resolution, trying search",
                                provider.name()
                            );

                            if provider.supports_search() {
                                if let Some(date_str) =
                                    extract_date_from_product_id(product_id_or_url)
                                {
                                    let search_params = SearchParams {
                                        product_type: extract_product_type_from_id(
                                            product_id_or_url,
                                        ),
                                        platform: extract_platform_from_id(product_id_or_url),
                                        start_date: date_str - chrono::Duration::days(1),
                                        end_date: date_str + chrono::Duration::days(1),
                                        aoi_wkt: None,
                                        max_results: 10,
                                        acquisition_mode: None,
                                        polarization: None,
                                        orbit_direction: None,
                                        relative_orbit: None,
                                    };

                                    if let Ok(search_results) =
                                        provider.search_products(&search_params)
                                    {
                                        if let Some(matching_product) =
                                            search_results.iter().find(|p| {
                                                p.id == product_id_or_url
                                                    || p.id.contains(product_id_or_url)
                                                    || p.title == product_id_or_url
                                                    || p.title.contains(product_id_or_url)
                                            })
                                        {
                                            if let Some(ref download_url) =
                                                matching_product.download_url
                                            {
                                                debug!("Found product via search, using download URL: {}", download_url);
                                                resolved_urls.push(download_url.clone());
                                                if expected_md5.is_none() {
                                                    expected_md5 =
                                                        matching_product.checksum_md5.clone();
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                            continue;
                        }
                        Err(e) => {
                            debug!("Provider {} failed to resolve URL: {}", provider.name(), e);
                            continue;
                        }
                    }
                }
            }

            if resolved_urls.is_empty() {
                return Err(SarError::Processing(format!(
                    "Could not resolve product ID to download URL: {}. Try:\n  1. Providing a direct download URL\n  2. Using ASF search API to find available products\n  3. Ensuring ASF credentials are configured (SARDINE_ASF_USERNAME, SARDINE_ASF_PASSWORD)",
                    sanitize_product_id(product_id_or_url) // Sanitize in error message too
                )));
            }

            // Log that we're trying multiple URLs (helpful for AWS)
            if resolved_urls.len() > 1 {
                debug!(
                    "Trying {} URLs for product: {}",
                    resolved_urls.len(),
                    product_id_or_url
                );
            }

            resolved_urls
        };

        // Try downloading from each URL
        // FIXED: Sanitize product ID to prevent path traversal
        let sanitized_id = sanitize_product_id(product_id_or_url);
        let output_path = self
            .cache
            .products_dir()
            .join(format!("{}.zip", sanitized_id));

        // Try each resolved URL
        for download_url in &download_urls {
            debug!("Trying download URL: {}", download_url);

            // Use ASF provider for all downloads
            if let Some(provider) = self.providers.get_provider("asf") {
                match provider.download_product(download_url, &output_path, progress.as_ref()) {
                    Ok(path) => {
                        // Verify checksum if configured and available
                        if self.config.verify_checksums {
                            if let Some(ref md5) = expected_md5 {
                                match super::utils::validate_checksum(&path, md5) {
                                    Ok(true) => { /* ok */ }
                                    Ok(false) => {
                                        warn!("Checksum mismatch for {}. Expected MD5 {}. Retrying...", product_id_or_url, md5);
                                        let _ = std::fs::remove_file(&path);
                                        continue;
                                    }
                                    Err(e) => {
                                        warn!("Checksum validation error: {}. Retrying...", e);
                                        let _ = std::fs::remove_file(&path);
                                        continue;
                                    }
                                }
                            }
                        }

                        // Add to cache
                        let metadata = fs::metadata(&path).map_err(SarError::Io)?;
                        let cache_entry = CacheEntry {
                            path: path.clone(),
                            size: metadata.len(),
                            modified: SystemTime::now(),
                            checksum: None,
                        };
                        self.cache
                            .add_product(product_id_or_url.to_string(), cache_entry);

                        info!("Successfully downloaded product from ASF");
                        return Ok(DownloadResult {
                            path,
                            size: metadata.len(),
                            source: "asf".to_string(),
                        });
                    }
                    Err(e) => {
                        warn!("Download failed from ASF: {}", e);
                    }
                }
            }

            // If provider download failed, try direct download as fallback
            debug!("Trying direct download as fallback for: {}", download_url);
            match super::utils::download_from_url(
                download_url,
                Some(&output_path),
                self.config.timeout_seconds,
            ) {
                Ok(_) => {
                    // Verify checksum if configured and available
                    if self.config.verify_checksums {
                        if let Some(ref md5) = expected_md5 {
                            match super::utils::validate_checksum(&output_path, md5) {
                                Ok(true) => { /* ok */ }
                                Ok(false) => {
                                    warn!("Checksum mismatch for {} via direct download. Expected MD5 {}. Trying next URL...", product_id_or_url, md5);
                                    let _ = std::fs::remove_file(&output_path);
                                    continue;
                                }
                                Err(e) => {
                                    warn!("Checksum validation error: {}. Trying next URL...", e);
                                    let _ = std::fs::remove_file(&output_path);
                                    continue;
                                }
                            }
                        }
                    }
                    let metadata = fs::metadata(&output_path).map_err(SarError::Io)?;
                    let cache_entry = CacheEntry {
                        path: output_path.clone(),
                        size: metadata.len(),
                        modified: SystemTime::now(),
                        checksum: None,
                    };
                    self.cache
                        .add_product(product_id_or_url.to_string(), cache_entry);
                    info!("Successfully downloaded product via direct download");
                    return Ok(DownloadResult {
                        path: output_path,
                        size: metadata.len(),
                        source: "direct".to_string(),
                    });
                }
                Err(e) => {
                    warn!("Direct download also failed: {}", e);
                }
            }

            warn!("Download failed for URL: {}", download_url);
        }

        Err(SarError::Processing(format!(
            "Failed to download product {} from all available providers",
            sanitize_product_id(product_id_or_url) // Sanitize in error message
        )))
    }

    /// Search for products
    pub fn search(&self, params: &SearchParams) -> SarResult<Vec<ProductMetadata>> {
        // Validate search parameters
        if params.start_date > params.end_date {
            return Err(SarError::Processing(format!(
                "Invalid date range: start_date ({}) must be before end_date ({})",
                params.start_date, params.end_date
            )));
        }

        if params.max_results == 0 || params.max_results > 10000 {
            return Err(SarError::Processing(format!(
                "Invalid max_results: must be between 1 and 10000, got {}",
                params.max_results
            )));
        }

        info!("Searching for products with params: {:?}", params);

        let mut all_results = Vec::new();

        // Try each provider that supports search
        for provider_name in self.providers.enabled() {
            if let Some(provider) = self.providers.get_provider(provider_name) {
                if provider.supports_search() && provider.is_available() {
                    debug!("Searching with provider: {}", provider.name());
                    match provider.search_products(params) {
                        Ok(mut results) => {
                            debug!(
                                "Provider {} returned {} results",
                                provider.name(),
                                results.len()
                            );
                            all_results.append(&mut results);

                            // If we have enough results, stop searching
                            if all_results.len() >= params.max_results {
                                all_results.truncate(params.max_results);
                                break;
                            }
                        }
                        Err(e) => {
                            warn!("Search failed with provider {}: {}", provider.name(), e);
                            // Continue to next provider
                        }
                    }
                }
            }
        }

        if all_results.is_empty() {
            warn!("No products found matching search criteria");
        } else {
            info!(
                "Found {} products matching search criteria",
                all_results.len()
            );
        }

        Ok(all_results)
    }
}
