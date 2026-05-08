//! Shared utilities for download operations

use crate::types::{SarError, SarResult};
use reqwest::blocking::Client;
use std::io::{Read, Write};
use std::path::{Path, PathBuf};
use std::thread;
use std::time::{Duration, Instant};

/// Network retry configuration
pub const MAX_RETRIES: u32 = 3;
pub const INITIAL_TIMEOUT_SECS: u64 = 60;
pub const DOWNLOAD_BUFFER_SIZE: usize = 8192;
pub const OVERALL_TIMEOUT_MULTIPLIER: u64 = 2; // Allow 2x timeout for slow connections

/// Create a configured HTTP client with retry support
pub fn create_http_client(timeout_secs: u64) -> SarResult<Client> {
    let mut builder = reqwest::blocking::Client::builder()
        .timeout(Duration::from_secs(timeout_secs))
        .cookie_store(true)
        .user_agent("SARdine/1.0 (Scientific SAR Processing)");

    // Add proxy support from environment variables
    if let Ok(proxy_url) = std::env::var("HTTP_PROXY") {
        if let Ok(proxy) = reqwest::Proxy::http(&proxy_url) {
            builder = builder.proxy(proxy);
            log::info!("Using HTTP proxy: {}", proxy_url);
        }
    }
    if let Ok(proxy_url) = std::env::var("HTTPS_PROXY") {
        if let Ok(proxy) = reqwest::Proxy::https(&proxy_url) {
            builder = builder.proxy(proxy);
            log::info!("Using HTTPS proxy: {}", proxy_url);
        }
    }

    builder
        .build()
        .map_err(|e| SarError::Processing(format!("HTTP client creation failed: {}", e)))
}

/// Download content from a URL with improved network resilience and progress reporting
pub fn download_from_url(
    url: &str,
    output_path: Option<&Path>,
    timeout_secs: u64,
) -> SarResult<Vec<u8>> {
    download_from_url_with_progress(url, output_path, timeout_secs, None)
}

/// Download content from a URL with progress reporting and resume capability
pub fn download_from_url_with_progress(
    url: &str,
    output_path: Option<&Path>,
    timeout_secs: u64,
    progress: Option<&super::progress::ProgressCallback>,
) -> SarResult<Vec<u8>> {
    log::debug!("Downloading from: {}", url);

    let client = create_http_client(timeout_secs)?;
    let mut last_error: Option<String> = None;
    let overall_timeout = Duration::from_secs(timeout_secs * OVERALL_TIMEOUT_MULTIPLIER);
    let start_time = Instant::now();

    // Check for existing partial file if output_path is specified
    let mut existing_size = 0u64;
    let mut file_handle: Option<std::fs::File> = None;
    let mut final_path: Option<PathBuf> = None;
    let mut temp_path: Option<PathBuf> = None;

    if let Some(path) = output_path {
        final_path = Some(path.to_path_buf());
        // Use a sidecar ".part" file for atomic writes
        let tmp = path
            .parent()
            .map(|p| {
                p.join(format!(
                    "{}.part",
                    path.file_name().unwrap_or_default().to_string_lossy()
                ))
            })
            .unwrap_or_else(|| PathBuf::from(format!("{}.part", path.display())));
        temp_path = Some(tmp.clone());

        // Prefer resuming from .part; if final exists and looks complete, caller may already validate/use it
        let resume_target = if tmp.exists() {
            Some(tmp.as_path())
        } else if path.exists() {
            Some(path)
        } else {
            None
        };
        if let Some(resume_path) = resume_target {
            if let Ok(metadata) = std::fs::metadata(resume_path) {
                existing_size = metadata.len();
                if existing_size > 0 {
                    log::info!(
                        "Resuming download from {} bytes ({}).",
                        existing_size,
                        resume_path.display()
                    );
                    file_handle = Some(
                        std::fs::OpenOptions::new()
                            .create(true)
                            .append(true)
                            .open(resume_path)
                            .map_err(SarError::Io)?,
                    );
                }
            }
        }
    }

    for attempt in 0..MAX_RETRIES {
        if attempt > 0 {
            let delay = Duration::from_secs(2_u64.pow(attempt)); // Exponential backoff
            log::info!(
                "Retrying download (attempt {}/{}) after {}s delay...",
                attempt + 1,
                MAX_RETRIES,
                delay.as_secs()
            );
            thread::sleep(delay);
        }

        // Keep existing partial state across retries to truly resume

        let mut request = client.get(url);

        // Add Range header for resume if we have existing data
        if existing_size > 0 {
            request = request.header("Range", format!("bytes={}-", existing_size));
        }

        match request.send() {
            Ok(response) => {
                let status = response.status();

                // Handle partial content (206) for resume
                let is_partial = status == reqwest::StatusCode::PARTIAL_CONTENT;
                // 416 Range Not Satisfiable means the file is already complete
                let is_range_not_satisfiable = status == reqwest::StatusCode::RANGE_NOT_SATISFIABLE;
                let is_success = status.is_success() || is_partial || is_range_not_satisfiable;

                if is_success {
                    // Handle 416 Range Not Satisfiable - file is already complete
                    if is_range_not_satisfiable {
                        log::info!("File is already complete (416 Range Not Satisfiable)");
                        if let Some(path) = output_path {
                            if path.exists() {
                                let bytes = std::fs::read(path).map_err(SarError::Io)?;
                                if !bytes.is_empty() {
                                    log::info!(
                                        "Using existing complete file: {} ({} bytes)",
                                        path.display(),
                                        bytes.len()
                                    );
                                    return Ok(bytes);
                                }
                            }
                        }
                        // If no file exists, treat as error
                        return Err(SarError::Processing(
                            "416 Range Not Satisfiable but no file exists".to_string(),
                        ));
                    }

                    let total_size = if is_partial {
                        // For partial content, extract total size from Content-Range header
                        response
                            .headers()
                            .get("Content-Range")
                            .and_then(|h| h.to_str().ok())
                            .and_then(|s| s.split('/').last().and_then(|s| s.parse::<u64>().ok()))
                            .unwrap_or_else(|| {
                                existing_size + response.content_length().unwrap_or(0)
                            })
                    } else {
                        response.content_length().unwrap_or(0)
                    };

                    let download_start = Instant::now();
                    let mut downloaded = existing_size;

                    let mut reader = response;
                    let mut buffer = vec![0u8; DOWNLOAD_BUFFER_SIZE];
                    let mut all_bytes = Vec::new();
                    let use_file = file_handle.is_some();

                    loop {
                        let bytes_read = reader.read(&mut buffer).map_err(|e| {
                            SarError::Processing(format!("Failed to read chunk: {}", e))
                        })?;

                        if bytes_read == 0 {
                            break;
                        }

                        if use_file {
                            if let Some(ref mut file) = file_handle {
                                file.write_all(&buffer[..bytes_read]).map_err(|e| {
                                    SarError::Processing(format!(
                                        "Failed to write chunk to {}: {}",
                                        output_path
                                            .map(|p| p.display().to_string())
                                            .unwrap_or_else(|| "file".to_string()),
                                        e
                                    ))
                                })?;
                            }
                        } else {
                            all_bytes.extend_from_slice(&buffer[..bytes_read]);
                        }

                        // FIXED: Use saturating_add to prevent integer overflow
                        downloaded = downloaded.saturating_add(bytes_read as u64);

                        // Check overall timeout
                        if start_time.elapsed() > overall_timeout {
                            return Err(SarError::Processing(format!(
                                "Download timeout exceeded ({}s) for URL: {}",
                                overall_timeout.as_secs(),
                                url
                            )));
                        }

                        // Report progress
                        if let Some(callback) = progress {
                            let elapsed = download_start.elapsed().as_secs_f64();
                            let bytes_per_sec = if elapsed > 0.0 {
                                downloaded as f64 / elapsed
                            } else {
                                0.0
                            };
                            callback(downloaded, total_size, bytes_per_sec);
                        }
                    }

                    if !use_file {
                        // Write all bytes atomically if an output path was provided
                        if let Some(path) = &final_path {
                            if let Some(parent) = path.parent() {
                                std::fs::create_dir_all(parent).map_err(SarError::Io)?;
                            }
                            let tmp = temp_path
                                .as_ref()
                                .cloned()
                                .unwrap_or_else(|| path.with_extension("part"));
                            std::fs::write(&tmp, &all_bytes).map_err(SarError::Io)?;
                            // Atomic rename to final
                            let _ = std::fs::remove_file(path).ok();
                            std::fs::rename(&tmp, path).map_err(SarError::Io)?;
                            log::info!("File saved to: {}", path.display());
                        }
                        return Ok(all_bytes);
                    } else {
                        // File was appended to temp; finalize if complete
                        if let (Some(path), Some(tmp)) = (&final_path, &temp_path) {
                            // Flush and sync the file to ensure it's written
                            if let Some(ref mut file) = file_handle {
                                file.flush().map_err(SarError::Io)?;
                                file.sync_all().map_err(SarError::Io)?;
                            }

                            // Verify file size matches expected total
                            if total_size > 0 && downloaded < total_size {
                                log::warn!(
                                    "Download incomplete: {} / {} bytes",
                                    downloaded,
                                    total_size
                                );
                                if attempt < MAX_RETRIES - 1 {
                                    last_error = Some(format!(
                                        "Incomplete download: {}/{} bytes",
                                        downloaded, total_size
                                    ));
                                    continue;
                                }
                            }

                            // If 416 indicated complete file present, prefer final
                            if status == reqwest::StatusCode::RANGE_NOT_SATISFIABLE && path.exists()
                            {
                                let bytes = std::fs::read(path).map_err(SarError::Io)?;
                                if !bytes.is_empty() {
                                    log::info!(
                                        "Using existing complete file: {} ({} bytes)",
                                        path.display(),
                                        bytes.len()
                                    );
                                    return Ok(bytes);
                                }
                            }

                            // Move temp to final atomically
                            if let Some(parent) = path.parent() {
                                std::fs::create_dir_all(parent).map_err(SarError::Io)?;
                            }
                            let _ = std::fs::remove_file(path).ok();
                            std::fs::rename(tmp, path).map_err(SarError::Io)?;
                            let bytes = std::fs::read(path).map_err(SarError::Io)?;
                            log::info!(
                                "Download completed, file saved to: {} ({} bytes)",
                                path.display(),
                                bytes.len()
                            );
                            return Ok(bytes);
                        }
                    }
                } else if (status.is_server_error() || status == 429) && attempt < MAX_RETRIES - 1 {
                    // Retry on 5xx errors and 429 (Too Many Requests)
                    log::warn!(
                        "Server error {} on attempt {}, will retry",
                        status,
                        attempt + 1
                    );
                    last_error = Some(format!("HTTP {}", status));
                    continue;
                } else {
                    return Err(SarError::Processing(format!(
                        "HTTP request failed with status: {} for URL: {} (attempt {}/{})",
                        status,
                        url,
                        attempt + 1,
                        MAX_RETRIES
                    )));
                }
            }
            Err(e) => {
                if (e.is_timeout() || e.is_connect() || e.is_request()) && attempt < MAX_RETRIES - 1
                {
                    log::warn!(
                        "Network error on attempt {}, will retry: {}",
                        attempt + 1,
                        e
                    );
                    last_error = Some(format!("Network error: {}", e));
                    continue;
                } else {
                    return Err(SarError::Processing(format!(
                        "HTTP request failed after {} attempts: {}. URL: {}",
                        attempt + 1,
                        e,
                        sanitize_url_for_error(url)
                    )));
                }
            }
        }
    }

    Err(SarError::Processing(format!(
        "Failed to download after {} attempts. Last error: {:?}. URL: {}",
        MAX_RETRIES,
        last_error,
        sanitize_url_for_error(url)
    )))
}

/// Sanitize URL for error messages to prevent information leakage
pub fn sanitize_url_for_error(url: &str) -> String {
    // Remove query parameters that might contain tokens
    if let Some(pos) = url.find('?') {
        format!("{}?[REDACTED]", &url[..pos])
    } else {
        url.to_string()
    }
}

/// Check if content is a ZIP file by examining magic bytes
pub fn is_zip_content(bytes: &[u8]) -> bool {
    bytes.len() >= 4 && bytes[0..4] == [0x50, 0x4B, 0x03, 0x04]
}

/// Extract content from ZIP file
pub fn extract_from_zip(zip_bytes: &[u8], extension: &str) -> SarResult<Vec<u8>> {
    use std::io::Cursor;
    use std::io::Read;
    use zip::ZipArchive;

    let cursor = Cursor::new(zip_bytes);
    let mut archive = ZipArchive::new(cursor)
        .map_err(|e| SarError::Processing(format!("Failed to read ZIP archive: {}", e)))?;

    for i in 0..archive.len() {
        let mut file = archive
            .by_index(i)
            .map_err(|e| SarError::Processing(format!("Failed to read ZIP entry {}: {}", i, e)))?;

        if file.name().ends_with(extension) {
            log::debug!("Found {} file in ZIP: {}", extension, file.name());

            let mut contents = Vec::new();
            file.read_to_end(&mut contents)
                .map_err(|e| SarError::Processing(format!("Failed to read file: {}", e)))?;

            return Ok(contents);
        }
    }

    Err(SarError::Processing(format!(
        "No {} file found in ZIP archive",
        extension
    )))
}

/// Validate file checksum (MD5)
pub fn validate_checksum(file_path: &Path, expected_md5: &str) -> SarResult<bool> {
    use std::fs::File;
    use std::io::Read;

    let mut file = File::open(file_path).map_err(SarError::Io)?;
    let mut buffer = Vec::new();
    file.read_to_end(&mut buffer).map_err(SarError::Io)?;

    let digest = md5::compute(&buffer);
    let actual_md5 = format!("{:x}", digest);

    Ok(actual_md5 == expected_md5.to_lowercase())
}
