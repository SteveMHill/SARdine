//! ASF DAAC provider

use super::super::utils::create_http_client;
use super::{DownloadProvider, ProductMetadata, ProviderConfig, SearchParams};
use crate::types::{SarError, SarResult};
use log::debug;
use reqwest::blocking::Client;
use serde_json::Value;

/// ASF DAAC provider
pub struct ASFProvider {
    username: Option<String>,
    password: Option<String>,
    token: Option<String>,
    base_url: String,
    client: Client,
}

impl ASFProvider {
    pub fn new(config: &ProviderConfig) -> SarResult<Self> {
        // ASF has a 30-second timeout on searches (as of May 2025)
        // Use a shorter timeout to handle this gracefully
        let client = create_http_client(35)?; // Slightly longer than API timeout to detect it

        Ok(Self {
            username: config.username.clone(),
            password: config.password.clone(),
            token: config.token.clone(),
            base_url: config
                .base_url
                .clone()
                .unwrap_or_else(|| "https://api.daac.asf.alaska.edu".to_string()),
            client,
        })
    }

    /// Resolve product ID to ASF download URL
    /// ASF uses multiple URL patterns:
    /// 1. https://datapool.asf.alaska.edu/{product_type}/{satellite}/{product_id}.zip (most common)
    /// 2. https://datapool.asf.alaska.edu/{product_type}/{product_id}.zip (fallback)
    pub(crate) fn resolve_product_url_impl(&self, product_id: &str) -> SarResult<Option<String>> {
        // ASF URL format: https://datapool.asf.alaska.edu/{product_type}/{satellite}/{product_id}.zip
        // Or: https://datapool.asf.alaska.edu/{product_type}/{product_id}.zip
        // Try to determine product type from product ID
        let product_type = if product_id.contains("_SLC_") {
            "SLC"
        } else if product_id.contains("_GRD_") {
            "GRD"
        } else {
            "SLC" // Default
        };

        // Determine satellite from product ID (S1A -> SA, S1B -> SB)
        let satellite = if product_id.starts_with("S1A") {
            "SA"
        } else if product_id.starts_with("S1B") {
            "SB"
        } else {
            // Try without satellite prefix
            let url = format!(
                "https://datapool.asf.alaska.edu/{}/{}.zip",
                product_type, product_id
            );
            debug!("Resolved ASF URL for {}: {}", product_id, url);
            return Ok(Some(url));
        };

        // Try with satellite prefix first (most common format)
        let url = format!(
            "https://datapool.asf.alaska.edu/{}/{}/{}.zip",
            product_type, satellite, product_id
        );
        debug!("Resolved ASF URL for {}: {}", product_id, url);
        Ok(Some(url))
    }
}

impl DownloadProvider for ASFProvider {
    fn name(&self) -> &str {
        "asf"
    }

    fn supports_search(&self) -> bool {
        true
    }

    fn requires_auth(&self) -> bool {
        false // ASF supports anonymous access
    }

    fn is_available(&self) -> bool {
        true // Always available (supports anonymous)
    }

    fn resolve_product_url(&self, product_id: &str) -> SarResult<Option<String>> {
        // Delegate to the implementation method
        self.resolve_product_url_impl(product_id)
    }

    fn download_product(
        &self,
        url: &str,
        output_path: &std::path::Path,
        _progress: Option<&super::super::progress::ProgressCallback>,
    ) -> SarResult<std::path::PathBuf> {
        // ASF datapool requires Earthdata (URS) authentication for many products.
        // Implement a pragmatic URS login flow with cookie persistence:
        // 1) Try direct GET (may work for public files)
        // 2) If 401/redirect to URS, perform GET on URS login URL with Basic auth
        // 3) Retry original URL with cookies now set

        // Helper to stream response into a .part file and atomically rename
        fn stream_to_file(
            mut resp: reqwest::blocking::Response,
            output_path: &std::path::Path,
        ) -> SarResult<std::path::PathBuf> {
            use std::io::{Read, Write};
            let tmp = output_path
                .parent()
                .map(|p| {
                    p.join(format!(
                        "{}.part",
                        output_path
                            .file_name()
                            .unwrap_or_default()
                            .to_string_lossy()
                    ))
                })
                .unwrap_or_else(|| output_path.with_extension("part"));
            if let Some(parent) = output_path.parent() {
                std::fs::create_dir_all(parent).map_err(SarError::Io)?;
            }
            let mut file = std::fs::File::create(&tmp)
                .map_err(|e| SarError::Processing(format!("Failed to create temp file: {}", e)))?;

            let mut buf = [0u8; 8192];
            loop {
                let n = resp
                    .read(&mut buf)
                    .map_err(|e| SarError::Processing(format!("Failed to read chunk: {}", e)))?;
                if n == 0 {
                    break;
                }
                file.write_all(&buf[..n])
                    .map_err(|e| SarError::Processing(format!("Failed to write chunk: {}", e)))?;
            }
            file.flush().map_err(SarError::Io)?;
            file.sync_all().map_err(SarError::Io)?;
            let _ = std::fs::remove_file(output_path).ok();
            std::fs::rename(&tmp, output_path).map_err(SarError::Io)?;
            Ok(output_path.to_path_buf())
        }

        // First attempt: direct GET, include Basic auth and Bearer token if provided
        let mut req = self.client.get(url);
        if let (Some(user), Some(pass)) = (&self.username, &self.password) {
            req = req.basic_auth(user, Some(pass));
        }
        if let Some(token) = &self.token {
            req = req.bearer_auth(token);
        }
        let mut response = req
            .send()
            .map_err(|e| SarError::Processing(format!("ASF download failed: {}", e)))?;
        if response.status().is_success() {
            return stream_to_file(response, output_path);
        }

        // If we have credentials, try multiple authentication strategies
        if let (Some(user), Some(pass)) = (&self.username, &self.password) {
            debug!("Attempting Earthdata authentication for user: {}", user);

            // Strategy 1: If we have a Bearer token, use it first
            if let Some(token) = &self.token {
                debug!("Trying Bearer token authentication");
                let resp = self
                    .client
                    .get(url)
                    .bearer_auth(token)
                    .send()
                    .map_err(|e| SarError::Processing(format!("ASF Bearer auth failed: {}", e)))?;

                if resp.status().is_success() {
                    return stream_to_file(resp, output_path);
                }
                debug!("Bearer token didn't work, status: {}", resp.status());
            }

            // Strategy 2: Try to get a fresh Bearer token from Earthdata API
            debug!("Trying to acquire Earthdata Bearer token via API");
            let token_resp = self
                .client
                .get("https://urs.earthdata.nasa.gov/api/users/tokens")
                .basic_auth(user, Some(pass))
                .send();

            if let Ok(resp) = token_resp {
                if resp.status().is_success() {
                    if let Ok(text) = resp.text() {
                        // Parse the token array and use the first valid token
                        if let Ok(tokens) = serde_json::from_str::<Vec<serde_json::Value>>(&text) {
                            for token_obj in &tokens {
                                if let Some(access_token) =
                                    token_obj.get("access_token").and_then(|t| t.as_str())
                                {
                                    debug!("Got Earthdata Bearer token, trying download");
                                    let retry = self
                                        .client
                                        .get(url)
                                        .bearer_auth(access_token)
                                        .send()
                                        .map_err(|e| {
                                            SarError::Processing(format!(
                                                "ASF download with token failed: {}",
                                                e
                                            ))
                                        })?;

                                    if retry.status().is_success() {
                                        return stream_to_file(retry, output_path);
                                    }
                                    debug!(
                                        "Token from API didn't work, status: {}",
                                        retry.status()
                                    );
                                    break;
                                }
                            }
                        }
                    }
                }
            }

            // Strategy 3: Create new token if none exist
            debug!("Trying to create new Earthdata token");
            let create_resp = self
                .client
                .post("https://urs.earthdata.nasa.gov/api/users/tokens")
                .basic_auth(user, Some(pass))
                .send();

            if let Ok(resp) = create_resp {
                if resp.status().is_success() {
                    if let Ok(text) = resp.text() {
                        if let Ok(token_obj) = serde_json::from_str::<serde_json::Value>(&text) {
                            if let Some(access_token) =
                                token_obj.get("access_token").and_then(|t| t.as_str())
                            {
                                debug!("Created new Earthdata token, trying download");
                                let retry = self
                                    .client
                                    .get(url)
                                    .bearer_auth(access_token)
                                    .send()
                                    .map_err(|e| {
                                    SarError::Processing(format!(
                                        "ASF download with new token failed: {}",
                                        e
                                    ))
                                })?;

                                if retry.status().is_success() {
                                    return stream_to_file(retry, output_path);
                                }
                            }
                        }
                    }
                }
            }

            // Strategy 4: OAuth cookie flow - warm up session then retry
            debug!("Trying OAuth cookie flow");
            let _ = self
                .client
                .get("https://urs.earthdata.nasa.gov/")
                .basic_auth(user, Some(pass))
                .send();

            // Then hit the OAuth authorize endpoint to establish cookies
            let _ = self.client
                .get(format!(
                    "https://urs.earthdata.nasa.gov/oauth/authorize?client_id=BO_n7nTIlMljdvU6kRRB3g&response_type=code&redirect_uri=https://sentinel1.asf.alaska.edu/login&state=/SLC/SA/{}.zip&app_type=401",
                    output_path.file_stem().unwrap_or_default().to_string_lossy()
                ))
                .basic_auth(user, Some(pass))
                .send();

            // Final attempt with cookies
            response = self
                .client
                .get(url)
                .send()
                .map_err(|e| SarError::Processing(format!("ASF final retry failed: {}", e)))?;

            if response.status().is_success() {
                return stream_to_file(response, output_path);
            }
            debug!(
                "All ASF auth strategies failed, final status: {}",
                response.status()
            );
        }

        // If we reached here, authentication failed or file unavailable
        Err(SarError::Processing(format!(
            "ASF download failed with status: {}",
            response.status()
        )))
    }

    fn search_products(&self, params: &SearchParams) -> SarResult<Vec<ProductMetadata>> {
        let search_url = format!("{}/services/search/param", self.base_url);

        // ASF has a 30-second timeout, so limit max_results to avoid timeouts
        let max_results = std::cmp::min(params.max_results, 1000);

        // ASF supports output formats: JSON, GEOJSON, CSV, KML, METADATA
        // Using JSON for compatibility, but GEOJSON is also available
        let mut query_params = vec![
            ("platform", "SENTINEL-1A,SENTINEL-1B".to_string()),
            ("output", "JSON".to_string()), // Can also use "GEOJSON" for GeoJSON output
            ("maxResults", max_results.to_string()),
        ];

        if max_results < params.max_results {
            log::warn!(
                "ASF search limited to {} results (from {}) to avoid 30-second timeout",
                max_results,
                params.max_results
            );
        }

        if let Some(product_type) = &params.product_type {
            query_params.push(("processingLevel", product_type.clone()));
        }

        if let Some(platform) = &params.platform {
            // Convert S1A/S1B to ASF format
            let platform_str = if platform == "S1A" {
                "SENTINEL-1A"
            } else if platform == "S1B" {
                "SENTINEL-1B"
            } else {
                "SENTINEL-1A,SENTINEL-1B"
            };
            query_params.push(("platform", platform_str.to_string()));
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

        // ASF API uses GET with query parameters, not POST form
        let mut request = self.client.get(&search_url);

        // Add authentication if provided
        if let (Some(user), Some(pass)) = (&self.username, &self.password) {
            request = request.basic_auth(user, Some(pass));
        }

        // Add query parameters
        for (key, value) in &query_params {
            request = request.query(&[(key, value)]);
        }

        // Add proper user agent (ASF may require this)
        request = request.header("User-Agent", "SARdine/1.0 (Scientific SAR Processing)");

        let response = request
            .send()
            .map_err(|e| {
                // Check for timeout errors
                if e.is_timeout() {
                    SarError::Processing(format!(
                        "ASF search timed out (30-second limit). Try reducing max_results or date range. Error: {}",
                        e
                    ))
                } else {
                    SarError::Processing(format!("ASF search request failed: {}", e))
                }
            })?;

        if !response.status().is_success() {
            let status = response.status();
            let error_msg = if status == 403 {
                "ASF search returned 403 Forbidden. This may be due to rate limiting or API access restrictions. Consider using the 'asf_search' Python package as an alternative."
            } else if status == 408 || status == 504 {
                "ASF search timed out (30-second limit). Try reducing max_results or date range."
            } else {
                "ASF search failed"
            };

            return Err(SarError::Processing(format!(
                "{} with status: {}",
                error_msg, status
            )));
        }

        let text = response
            .text()
            .map_err(|e| SarError::Processing(format!("Failed to read ASF response: {}", e)))?;

        let json: Value = serde_json::from_str(&text)
            .map_err(|e| SarError::Processing(format!("Failed to parse ASF JSON: {}", e)))?;

        // Parse search results
        // ASF API returns results in nested array format: [[{...}, {...}]]
        let mut products = Vec::new();

        // Try nested array format first (most common)
        if let Some(outer_array) = json.as_array() {
            for inner_array in outer_array {
                if let Some(results) = inner_array.as_array() {
                    for result in results {
                        if let Some(product) = self.parse_asf_result(result) {
                            products.push(product);
                        }
                    }
                }
            }
        }

        // Fallback: try "results" key (if API changes)
        if products.is_empty() {
            if let Some(results) = json.get("results").and_then(|r| r.as_array()) {
                for result in results {
                    if let Some(product) = self.parse_asf_result(result) {
                        products.push(product);
                    }
                }
            }
        }

        debug!("ASF search returned {} products", products.len());
        Ok(products)
    }
}

impl ASFProvider {
    /// Parse an ASF search result
    fn parse_asf_result(&self, result: &Value) -> Option<ProductMetadata> {
        // ASF uses "sceneId" (camelCase) not "sceneID"
        let scene_id = result
            .get("sceneId")
            .and_then(|s| s.as_str())
            .or_else(|| result.get("sceneID").and_then(|s| s.as_str()))
            .or_else(|| result.get("granuleName").and_then(|s| s.as_str()))
            .or_else(|| result.get("productName").and_then(|s| s.as_str()))
            .map(|s| s.to_string())?;

        // Extract download URL (ASF provides this directly)
        let download_url = result
            .get("downloadUrl")
            .and_then(|u| u.as_str())
            .map(|s| s.to_string())
            .or_else(|| {
                // Fallback: construct URL from scene ID
                let product_type = result
                    .get("processingLevel")
                    .and_then(|t| t.as_str())
                    .unwrap_or("SLC");
                // ASF datapool structure: /SLC/SA/{product_id}.zip or /SLC/{product_id}.zip
                Some(format!(
                    "https://datapool.asf.alaska.edu/{}/{}.zip",
                    product_type, scene_id
                ))
            });

        // Extract size - ASF uses "sizeMB" (camelCase) or "fileSize"
        let size_mb = result
            .get("sizeMB")
            .and_then(|s| s.as_f64())
            .or_else(|| {
                result
                    .get("fileSize")
                    .and_then(|s| s.as_u64())
                    .map(|bytes| bytes as f64 / (1024.0 * 1024.0))
            })
            .unwrap_or(0.0);

        // Extract product type
        let product_type = result
            .get("processingLevel")
            .and_then(|t| t.as_str())
            .map(|s| s.to_string())
            .unwrap_or_else(|| "UNKNOWN".to_string());

        Some(ProductMetadata {
            id: scene_id.clone(),
            title: scene_id,
            product_type,
            size_mb,
            download_url,
            checksum_md5: result
                .get("md5sum")
                .and_then(|m| m.as_str())
                .map(|s| s.to_string()),
        })
    }
}
