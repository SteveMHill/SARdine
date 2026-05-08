use std::path::Path;

use crate::types::{SarError, SarResult};

/// Network timeout configuration (used by legacy wrapper)
const INITIAL_TIMEOUT_SECS: u64 = 60;

/// Download content from a URL (legacy orbit helper)
///
/// NOTE: Kept for backward compatibility; delegates to the unified
/// `crate::io::download::utils::download_from_url` for resilience.
pub fn download_from_url(url: &str, output_path: Option<&Path>) -> SarResult<String> {
    let bytes = crate::io::download::utils::download_from_url(
        url,
        output_path,
        INITIAL_TIMEOUT_SECS,
    )?;

    let content = if url.ends_with(".EOF.zip") || is_zip_content(&bytes) {
        log::debug!("Processing ZIP file from: {}", url);
        let eof_content = extract_eof_from_zip(&bytes)?;
        validate_eof_format(&eof_content)?;
        eof_content
    } else {
        let content = String::from_utf8(bytes)
            .map_err(|e| SarError::Processing(format!("Invalid UTF-8 content: {}", e)))?;
        validate_eof_format(&content)?;
        content
    };

    if let Some(path) = output_path {
        if let Some(parent) = path.parent() {
            std::fs::create_dir_all(parent).map_err(SarError::Io)?;
        }
        std::fs::write(path, &content).map_err(SarError::Io)?;
        log::info!("Orbit file saved to: {}", path.display());
    }

    Ok(content)
}

/// Process successful HTTP response
fn process_response(
    response: reqwest::blocking::Response,
    url: &str,
    output_path: Option<&Path>,
) -> SarResult<String> {
    let bytes = response
        .bytes()
        .map_err(|e| SarError::Processing(format!("Failed to read response bytes: {}", e)))?;

    let content = if url.ends_with(".EOF.zip") || is_zip_content(&bytes) {
        log::debug!("Processing ZIP file from: {}", url);
        let eof_content = extract_eof_from_zip(&bytes)?;
        
        // Validate EOF format
        validate_eof_format(&eof_content)?;
        eof_content
    } else {
        let content = String::from_utf8(bytes.to_vec())
            .map_err(|e| SarError::Processing(format!("Invalid UTF-8 content: {}", e)))?;
        
        // Validate EOF format
        validate_eof_format(&content)?;
        content
    };

    if let Some(path) = output_path {
        if let Some(parent) = path.parent() {
            std::fs::create_dir_all(parent).map_err(SarError::Io)?;
        }
        std::fs::write(path, &content).map_err(SarError::Io)?;
        log::info!("Orbit file saved to: {}", path.display());
    }

    Ok(content)
}

/// Validate EOF file format
fn validate_eof_format(content: &str) -> SarResult<()> {
    // Check for EOF header markers
    if !content.contains("<?xml") && !content.contains("<Earth_Explorer_File>") {
        return Err(SarError::InvalidFormat(
            "Orbit file does not appear to be valid EOF format: missing XML header".to_string(),
        ));
    }

    // Check for required sections
    let required_sections = ["<Data_Block>", "<List_of_OSVs>"];
    for section in required_sections {
        if !content.contains(section) {
            return Err(SarError::InvalidFormat(format!(
                "Orbit file missing required section: {}",
                section
            )));
        }
    }

    Ok(())
}

/// Check if content is a ZIP file by examining magic bytes
pub fn is_zip_content(bytes: &[u8]) -> bool {
    bytes.len() >= 4 && bytes[0..4] == [0x50, 0x4B, 0x03, 0x04]
}

/// Extract EOF content from ZIP file
pub fn extract_eof_from_zip(zip_bytes: &[u8]) -> SarResult<String> {
    use std::io::Cursor;
    use zip::ZipArchive;

    let cursor = Cursor::new(zip_bytes);
    let mut archive = ZipArchive::new(cursor)
        .map_err(|e| SarError::Processing(format!("Failed to read ZIP archive: {}", e)))?;

    for i in 0..archive.len() {
        let mut file = archive
            .by_index(i)
            .map_err(|e| SarError::Processing(format!("Failed to read ZIP entry {}: {}", i, e)))?;

        if file.name().ends_with(".EOF") {
            log::debug!("Found EOF file in ZIP: {}", file.name());

            use std::io::Read;
            let mut contents = String::new();
            file.read_to_string(&mut contents)
                .map_err(|e| SarError::Processing(format!("Failed to read EOF file: {}", e)))?;

            return Ok(contents);
        }
    }

    Err(SarError::Processing(
        "No .EOF file found in ZIP archive".to_string(),
    ))
}
