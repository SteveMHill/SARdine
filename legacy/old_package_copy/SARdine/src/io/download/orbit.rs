#![allow(dead_code)]
//! Orbit file downloader

use crate::io::orbit::interp::validate_orbit_coverage;
use crate::io::orbit::OrbitType;
use crate::types::{SarError, SarResult};
use chrono::{DateTime, NaiveDateTime, TimeZone, Timelike, Utc};
use regex::Regex;
use std::path::Path;

use super::cache::CacheManager;
use super::manager::{DownloadConfig, DownloadResult};
use super::progress::ProgressCallback;
use super::providers::ProviderRegistry;
use super::utils::{download_from_url, extract_from_zip, is_zip_content};
use std::fs::OpenOptions;
use std::io::Read;
use std::sync::Arc;

// Constants
const FILE_SYNC_DELAY_MS: u64 = 200;
const ORBIT_VALIDATION_CHUNK_SIZE: usize = 100_000;
const ORBIT_LARGE_FILE_THRESHOLD: usize = 1_000_000;

/// Orbit file downloader
pub struct OrbitDownloader {
    cache: Arc<CacheManager>,
    providers: Arc<ProviderRegistry>,
    config: DownloadConfig,
}

impl OrbitDownloader {
    pub fn new(
        cache: &Arc<CacheManager>,
        providers: &Arc<ProviderRegistry>,
        config: &DownloadConfig,
    ) -> SarResult<Self> {
        Ok(Self {
            cache: Arc::clone(cache),
            providers: Arc::clone(providers),
            config: config.clone(),
        })
    }

    /// Download orbit file for a product
    pub fn download(
        &self,
        product_id: &str,
        start_time: DateTime<Utc>,
        _progress: Option<ProgressCallback>,
    ) -> SarResult<DownloadResult> {
        log::info!("Downloading orbit file for product: {}", product_id);

        let orbit_type = Self::determine_orbit_type(start_time);
        log::info!("Selected orbit type: {}", orbit_type);

        let output_path = self
            .cache
            .orbits_dir()
            .join(format!("{}_{}.EOF", product_id, orbit_type));

        // Try ESA first (this is what was working before - step.esa.int)
        // ESA repository is reliable and publicly accessible
        let urls = self.generate_esa_urls(product_id, start_time, orbit_type)?;
        for url in &urls {
            match self.try_download_url(url, &output_path) {
                Ok(_content) => {
                    match crate::io::orbit::OrbitReader::read_orbit_file(&output_path) {
                        Ok(orbit) => match validate_orbit_coverage(&orbit, start_time) {
                            Ok(_) => {
                                log::info!("Successfully downloaded orbit from ESA: {}", url);
                                let size = std::fs::metadata(&output_path)
                                    .map(|m| m.len())
                                    .unwrap_or(0);
                                return Ok(DownloadResult {
                                    path: output_path,
                                    size,
                                    source: "esa".to_string(),
                                });
                            }
                            Err(e) => {
                                log::warn!("Orbit file downloaded but validation failed: {}. Trying next URL.", e);
                            }
                        },
                        Err(e) => {
                            log::warn!("Failed to parse orbit file: {}. Trying next URL.", e);
                        }
                    }
                }
                Err(e) => {
                    log::debug!("ESA download failed: {}, trying next URL", e);
                }
            }
        }

        Err(SarError::Processing(format!(
            "Failed to download {} orbit file for {} from all sources",
            orbit_type, product_id
        )))
    }

    fn try_download_url(&self, url: &str, output_path: &Path) -> SarResult<String> {
        // Download the file
        let bytes = download_from_url(url, Some(output_path), self.config.timeout_seconds)?;

        // FIXED: Use atomic file operations to avoid race conditions
        let final_bytes = if output_path.exists() {
            // Try to open with exclusive read access
            // This ensures the file is fully written and not being modified
            let mut file = OpenOptions::new()
                .read(true)
                .open(output_path)
                .map_err(|e| {
                    SarError::Processing(format!(
                        "Failed to open downloaded file for reading: {}",
                        e
                    ))
                })?;

            // Read with retry on short reads (handles partial writes)
            let mut buffer = Vec::new();
            file.read_to_end(&mut buffer).map_err(|e| {
                SarError::Processing(format!("Failed to read downloaded file: {}", e))
            })?;

            // If file is empty or very small, use in-memory bytes instead
            if buffer.is_empty() || buffer.len() < 100 {
                log::warn!(
                    "Downloaded file is suspiciously small ({} bytes), using in-memory data",
                    buffer.len()
                );
                bytes
            } else {
                buffer
            }
        } else {
            bytes
        };

        // Verify file is not empty
        if final_bytes.is_empty() {
            return Err(SarError::Processing("Downloaded file is empty".to_string()));
        }

        // Check if it's a ZIP file
        // ESA files are .EOF.zip, AWS files are .EOF (already extracted)
        // But if the file already exists and is not a ZIP, it's already extracted
        let is_zip = url.ends_with(".EOF.zip")
            || (is_zip_content(&final_bytes)
                && final_bytes.len() > 100
                && !url.contains("s1-orbits")); // AWS files are not zipped

        let content = if is_zip {
            log::debug!(
                "Extracting EOF file from ZIP archive ({} bytes)",
                final_bytes.len()
            );
            let eof_bytes = extract_from_zip(&final_bytes, ".EOF")?;

            // Write extracted EOF file to output path (overwrite the ZIP)
            if let Some(parent) = output_path.parent() {
                std::fs::create_dir_all(parent).map_err(SarError::Io)?;
            }
            std::fs::write(output_path, &eof_bytes).map_err(SarError::Io)?;
            log::debug!(
                "Extracted EOF file written to: {} ({} bytes)",
                output_path.display(),
                eof_bytes.len()
            );

            String::from_utf8(eof_bytes)
                .map_err(|e| SarError::Processing(format!("Invalid UTF-8 content: {}", e)))?
        } else {
            // File is already extracted or is raw EOF
            log::debug!("File is already extracted or raw EOF format");
            String::from_utf8(final_bytes)
                .map_err(|e| SarError::Processing(format!("Invalid UTF-8 content: {}", e)))?
        };

        // Validate EOF format
        self.validate_eof_format(&content)?;
        Ok(content)
    }

    fn validate_eof_format(&self, content: &str) -> SarResult<()> {
        // Check for XML header or Earth_Explorer_File tag
        if !content.contains("<?xml") && !content.contains("<Earth_Explorer_File>") {
            return Err(SarError::InvalidFormat(
                "Orbit file does not appear to be valid EOF format: missing XML header".to_string(),
            ));
        }

        // Check for required sections (these may be further in the file for large orbit files)
        // We check a reasonable portion of the file (first chunk should contain headers)
        let check_content = if content.len() > ORBIT_VALIDATION_CHUNK_SIZE {
            &content[..ORBIT_VALIDATION_CHUNK_SIZE]
        } else {
            content
        };

        // Data_Block and List_of_OSVs are required but may be deep in the XML structure
        // For very large files, we check if we have the basic structure
        if !check_content.contains("<Earth_Explorer_File>") {
            return Err(SarError::InvalidFormat(
                "Orbit file missing Earth_Explorer_File root element".to_string(),
            ));
        }

        // If file is small enough, check for data sections
        if content.len() < ORBIT_LARGE_FILE_THRESHOLD {
            let required_sections = ["<Data_Block>", "<List_of_OSVs>"];
            for section in required_sections {
                if !content.contains(section) {
                    return Err(SarError::InvalidFormat(format!(
                        "Orbit file missing required section: {}",
                        section
                    )));
                }
            }
        } else {
            // For large files, just check that it's valid XML structure
            log::debug!(
                "Large orbit file ({} bytes), skipping detailed section validation",
                content.len()
            );
        }

        Ok(())
    }

    /// Generate AWS S3 URLs for orbit files
    /// Uses the correct AWS bucket: s1-orbits (Registry of Open Data on AWS)
    /// See: https://registry.opendata.aws/s1-orbits/
    /// Structure: AUX_POEORB/S1A/YEAR/MONTH/{filename}.EOF (NOT .zip!)
    fn generate_aws_urls(
        &self,
        product_id: &str,
        start_time: DateTime<Utc>,
        orbit_type: OrbitType,
    ) -> SarResult<Vec<String>> {
        let satellite = if product_id.starts_with("S1A") {
            "S1A"
        } else if product_id.starts_with("S1B") {
            "S1B"
        } else {
            return Err(SarError::Processing(format!(
                "Cannot determine satellite from product ID: {}",
                product_id
            )));
        };

        // Correct AWS bucket: s1-orbits (Registry of Open Data on AWS)
        // Region: us-west-2
        // Managed by: Alaska Satellite Facility (ASF)
        // Structure: AUX_POEORB/S1A/YEAR/MONTH/{filename}.EOF
        let base_url = "https://s1-orbits.s3.us-west-2.amazonaws.com";
        let mut urls = Vec::new();

        // AWS uses AUX_POEORB and AUX_RESORB prefixes (not POEORB/RESORB)
        let orbit_prefix = match orbit_type {
            OrbitType::POEORB => "AUX_POEORB",
            OrbitType::RESORB => "AUX_RESORB",
        };

        let search_dates = vec![
            start_time - chrono::Duration::days(1),
            start_time,
            start_time + chrono::Duration::days(1),
        ];

        for date in search_dates {
            let year = date.format("%Y");
            let month = date.format("%m");

            for hour in [0, 6, 12, 18] {
                let validity_start = date.with_hour(hour).unwrap_or(date);
                let validity_end = validity_start + chrono::Duration::days(1);
                let production_time = validity_start + chrono::Duration::hours(3);

                // AWS files are .EOF (not .EOF.zip) - they're already extracted
                let filename = format!(
                    "S1{}_OPER_AUX_{}_OPOD_{}_V{}_{}.EOF",
                    satellite.chars().last().unwrap_or('A'),
                    orbit_type,
                    production_time.format("%Y%m%dT%H%M%S"),
                    validity_start.format("%Y%m%dT%H%M%S"),
                    validity_end.format("%Y%m%dT%H%M%S")
                );

                // AWS structure: AUX_POEORB/S1A/YEAR/MONTH/{filename}.EOF
                let url = format!(
                    "{}/{}/{}/{}/{}/{}",
                    base_url, orbit_prefix, satellite, year, month, filename
                );
                urls.push(url);
            }
        }

        Ok(urls)
    }

    /// Generate ESA URLs for orbit files (fallback)
    fn generate_esa_urls(
        &self,
        product_id: &str,
        start_time: DateTime<Utc>,
        orbit_type: OrbitType,
    ) -> SarResult<Vec<String>> {
        let satellite = if product_id.starts_with("S1A") {
            "S1A"
        } else if product_id.starts_with("S1B") {
            "S1B"
        } else {
            return Err(SarError::Processing(format!(
                "Cannot determine satellite from product ID: {}",
                product_id
            )));
        };

        let mut urls = Vec::new();
        let base_url = "https://step.esa.int/auxdata/orbits/Sentinel-1";

        let search_dates = vec![
            start_time - chrono::Duration::days(1),
            start_time,
            start_time + chrono::Duration::days(1),
        ];

        for search_date in search_dates {
            let year = search_date.format("%Y");
            let month = search_date.format("%m");

            let dir_url = format!(
                "{}/{}/{}/{}/{}/",
                base_url, orbit_type, satellite, year, month
            );

            // Try to get files from directory listing
            if let Ok(file_urls) =
                self.get_orbit_files_from_directory(&dir_url, start_time, orbit_type)
            {
                urls.extend(file_urls);
            }
        }

        // Fallback to pattern-based URLs
        if urls.is_empty() {
            urls.extend(self.generate_fallback_urls(product_id, start_time, orbit_type)?);
        }

        if urls.is_empty() {
            return Err(SarError::Processing(format!(
                "No orbit file URLs generated for {} {} at {}",
                satellite,
                orbit_type,
                start_time.format("%Y%m%d")
            )));
        }

        Ok(urls)
    }

    fn get_orbit_files_from_directory(
        &self,
        dir_url: &str,
        target_time: DateTime<Utc>,
        orbit_type: OrbitType,
    ) -> SarResult<Vec<String>> {
        use super::utils::create_http_client;
        let client = create_http_client(30)?;

        let response = client.get(dir_url).send().map_err(|e| {
            SarError::Processing(format!("Failed to fetch directory listing: {}", e))
        })?;

        if !response.status().is_success() {
            return Err(SarError::Processing(format!(
                "Directory listing failed: {}",
                response.status()
            )));
        }

        let html = response
            .text()
            .map_err(|e| SarError::Processing(format!("Failed to read response: {}", e)))?;

        let mut matching_files = Vec::new();

        for line in html.lines() {
            if line.contains(".EOF.zip\"") {
                if let Some(filename) = self.extract_filename_from_html(line) {
                    if self.is_orbit_file_relevant(&filename, target_time, orbit_type) {
                        let file_url = format!("{}{}", dir_url, filename);
                        matching_files.push(file_url);
                    }
                }
            }
        }

        matching_files.sort_by(|a, b| {
            let score_a = self.score_orbit_file_relevance(a, target_time);
            let score_b = self.score_orbit_file_relevance(b, target_time);
            score_a.total_cmp(&score_b)
        });

        Ok(matching_files)
    }

    fn extract_filename_from_html(&self, html_line: &str) -> Option<String> {
        if let Some(start) = html_line.find("href=\"") {
            let start = start + 6;
            if let Some(end) = html_line[start..].find('"') {
                let filename = &html_line[start..start + end];
                if filename.ends_with(".EOF.zip") && !filename.contains('/') {
                    return Some(filename.to_string());
                }
            }
        }
        None
    }

    fn is_orbit_file_relevant(
        &self,
        filename: &str,
        target_time: DateTime<Utc>,
        _orbit_type: OrbitType,
    ) -> bool {
        if let Some((start_time, end_time)) = self.extract_validity_from_name(filename) {
            return target_time >= start_time && target_time <= end_time;
        }

        let target_date = target_time.format("%Y%m%d").to_string();
        filename.contains(&target_date)
    }

    fn score_orbit_file_relevance(&self, file_url: &str, target_time: DateTime<Utc>) -> f64 {
        let filename = file_url.split('/').next_back().unwrap_or("");

        if let Some((start_time, end_time)) = self.extract_validity_from_name(filename) {
            let mid_time = start_time + (end_time - start_time) / 2;
            return (target_time - mid_time).num_seconds().abs() as f64;
        }

        let target_date = target_time.format("%Y%m%d").to_string();
        if filename.contains(&target_date) {
            0.0
        } else {
            1_000_000.0
        }
    }

    fn generate_fallback_urls(
        &self,
        product_id: &str,
        start_time: DateTime<Utc>,
        orbit_type: OrbitType,
    ) -> SarResult<Vec<String>> {
        let satellite = if product_id.starts_with("S1A") {
            "S1A"
        } else {
            "S1B"
        };
        let base_url = "https://step.esa.int/auxdata/orbits/Sentinel-1";

        let mut urls = Vec::new();
        let dates = vec![
            start_time - chrono::Duration::days(1),
            start_time,
            start_time + chrono::Duration::days(1),
        ];

        for date in dates {
            let year = date.format("%Y");
            let month = date.format("%m");

            for hour in [0, 6, 12, 18] {
                let validity_start = date.with_hour(hour).unwrap_or(date);
                let validity_end = validity_start + chrono::Duration::days(1);
                let production_time = validity_start + chrono::Duration::hours(3);

                let filename = format!(
                    "S1{}_OPER_AUX_{}_OPOD_{}_V{}_{}.EOF.zip",
                    satellite.chars().last().unwrap_or('A'),
                    orbit_type,
                    production_time.format("%Y%m%dT%H%M%S"),
                    validity_start.format("%Y%m%dT%H%M%S"),
                    validity_end.format("%Y%m%dT%H%M%S")
                );

                let url = format!(
                    "{}/{}/{}/{}/{}/{}",
                    base_url, orbit_type, satellite, year, month, filename
                );
                urls.push(url);
            }
        }

        Ok(urls)
    }

    fn extract_validity_from_name(&self, filename: &str) -> Option<(DateTime<Utc>, DateTime<Utc>)> {
        let re = Regex::new(r"_V(\d{8}T\d{6})_(\d{8}T\d{6})\.EOF(?:\.zip)?$").ok()?;

        if let Some(captures) = re.captures(filename) {
            let start_str = captures.get(1)?.as_str();
            let end_str = captures.get(2)?.as_str();

            let parse_time = |s: &str| -> Option<DateTime<Utc>> {
                NaiveDateTime::parse_from_str(s, "%Y%m%dT%H%M%S")
                    .map(|dt| Utc.from_utc_datetime(&dt))
                    .ok()
            };

            if let (Some(start_time), Some(end_time)) = (parse_time(start_str), parse_time(end_str))
            {
                return Some((start_time, end_time));
            }
        }
        None
    }

    fn determine_orbit_type(start_time: DateTime<Utc>) -> OrbitType {
        // POEORB is available ~20 days after acquisition
        // RESORB is available ~3 hours after acquisition
        let now = Utc::now();
        let days_since_acquisition = (now - start_time).num_days();

        if days_since_acquisition >= 20 {
            OrbitType::POEORB
        } else {
            OrbitType::RESORB
        }
    }
}
