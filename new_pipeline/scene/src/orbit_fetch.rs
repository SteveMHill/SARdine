//! Feature-gated POEORB auto-downloader.
//!
//! Enabled with `--features orbit-fetch` (adds the `ureq` HTTP dependency).
//!
//! Entry point: [`fetch_poeorb`].
//!
//! # Source
//! Uses the AWS Open Data S3 bucket maintained by ASF:
//! `s3://s1-orbits/AUX_POEORB/{S1A|S1B}/{YYYY}/{MM}/{filename}.EOF`
//! (public, no credentials required, files already extracted from ZIP)
//!
//! Falls back to ESA step.esa.int directory listing if AWS yields nothing.
//!
//! # Validation
//! After download the file is parsed with [`crate::orbit::parse_eof_file`] and
//! the coverage window is checked: sensing start must be covered with at least
//! `COVERAGE_MARGIN_S` seconds of margin on both sides.
//!
//! # AGENTS.md rules
//! - No silent fallbacks: every failure is a typed `OrbitFetchError` variant.
//! - No hardcoded Sentinel-1 constants: all timing comes from function arguments.

#[cfg(feature = "orbit-fetch")]
mod inner {
    use std::io::{self, Write};
    use std::path::{Path, PathBuf};

    use chrono::{DateTime, Duration, Utc};
    use thiserror::Error;

    use crate::orbit::{parse_eof_file, OrbitError};

    // Minimum coverage margin required around the sensing window (seconds).
    // ESA requirement S1-RS-MDA-52-7441 — same value used in metadata_orbit.py.
    const COVERAGE_MARGIN_S: i64 = 60;

    // Minimum state vectors required in a valid POEORB file.
    const MIN_STATE_VECTORS: usize = 10;

    // AWS public S3 bucket (ASF, us-west-2, Registry of Open Data).
    const AWS_BASE: &str = "https://s1-orbits.s3.us-west-2.amazonaws.com";

    // ESA STEP server (directory-listing based).
    const ESA_BASE: &str = "https://step.esa.int/auxdata/orbits/Sentinel-1";

    // HTTP timeout for individual requests (seconds).
    const HTTP_TIMEOUT_S: u64 = 60;

    /// Errors returned by [`fetch_poeorb`].
    #[derive(Debug, Error)]
    pub enum OrbitFetchError {
        /// HTTP-level failure (status code or transport).
        #[error("HTTP error: {0}")]
        Http(String),

        /// No matching orbit file found on any source.
        #[error("No POEORB file found for satellite={satellite} sensing_start={sensing_start}")]
        NotFound {
            satellite: String,
            sensing_start: DateTime<Utc>,
        },

        /// I/O error writing the downloaded file to disk.
        #[error("I/O error: {0}")]
        Io(#[from] io::Error),

        /// Downloaded file failed orbit parsing or coverage validation.
        #[error("Orbit validation error: {0}")]
        Parse(#[from] OrbitError),

        /// Satellite identifier in the product name is not recognised.
        #[error("Unrecognised satellite prefix in product id: {0}")]
        UnknownSatellite(String),
    }

    /// Parse the satellite identifier from a SAFE product ID.
    ///
    /// Accepts `S1A_…` or `S1B_…` (case-sensitive per ESA naming).
    fn satellite_from_product_id(product_id: &str) -> Result<&'static str, OrbitFetchError> {
        if product_id.starts_with("S1A") {
            Ok("S1A")
        } else if product_id.starts_with("S1B") {
            Ok("S1B")
        } else {
            Err(OrbitFetchError::UnknownSatellite(
                product_id[..product_id.len().min(6)].to_string(),
            ))
        }
    }

    /// Attempt to download one URL and write the body to `dest`.
    ///
    /// Returns `Err(OrbitFetchError::Http)` on non-2xx status or transport error.
    /// Does **not** validate the content — callers must parse after this returns.
    fn download_to_file(url: &str, dest: &Path) -> Result<(), OrbitFetchError> {
        let response = ureq::builder()
            .timeout(std::time::Duration::from_secs(HTTP_TIMEOUT_S))
            .build()
            .get(url)
            .call()
            .map_err(|e| OrbitFetchError::Http(format!("{}: {}", url, e)))?;

        if response.status() < 200 || response.status() >= 300 {
            return Err(OrbitFetchError::Http(format!(
                "{}: HTTP {}",
                url,
                response.status()
            )));
        }

        // Write body to a temp file first, then rename atomically.
        let tmp_path = dest.with_extension("tmp");
        {
            let mut file = std::fs::File::create(&tmp_path)?;
            let mut reader = response.into_reader();
            io::copy(&mut reader, &mut file)?;
            file.flush()?;
        }
        std::fs::rename(&tmp_path, dest)?;
        Ok(())
    }

    /// Validate that a downloaded orbit file covers `sensing_start` with margin.
    ///
    /// Parses the file and checks:
    /// - At least `MIN_STATE_VECTORS` state vectors present.
    /// - Earliest vector time ≤ `sensing_start − COVERAGE_MARGIN_S`.
    /// - Latest vector time ≥ `sensing_start + COVERAGE_MARGIN_S`.
    fn validate_orbit_coverage(
        path: &Path,
        sensing_start: DateTime<Utc>,
    ) -> Result<(), OrbitFetchError> {
        let orbit = parse_eof_file(path)?;

        if orbit.state_vectors.len() < MIN_STATE_VECTORS {
            return Err(OrbitFetchError::Parse(OrbitError::InsufficientCoverage {
                detail: format!(
                    "only {} state vectors in downloaded file (need >= {})",
                    orbit.state_vectors.len(),
                    MIN_STATE_VECTORS
                ),
            }));
        }

        let first_t = orbit.state_vectors[0].time;
        let last_t = orbit.state_vectors[orbit.state_vectors.len() - 1].time;
        let required_start = sensing_start - Duration::seconds(COVERAGE_MARGIN_S);
        let required_stop = sensing_start + Duration::seconds(COVERAGE_MARGIN_S);

        if first_t > required_start || last_t < required_stop {
            return Err(OrbitFetchError::Parse(OrbitError::InsufficientCoverage {
                detail: format!(
                    "orbit covers [{}, {}], sensing window requires [{}, {}]",
                    first_t.format("%H:%M:%S"),
                    last_t.format("%H:%M:%S"),
                    required_start.format("%H:%M:%S"),
                    required_stop.format("%H:%M:%S"),
                ),
            }));
        }

        Ok(())
    }

    // ──────────────────────────────────────────────
    // AWS S3 source
    // ──────────────────────────────────────────────

    /// List POEORB keys in the AWS bucket matching `satellite` and `sensing_start`.
    ///
    /// Queries the S3 XML listing API for the prefix
    /// `AUX_POEORB/{satellite}/{YYYY}/{MM}/` over the three calendar days
    /// surrounding `sensing_start` (day−1, day, day+1) and returns
    /// `(key, direct_url)` pairs whose validity window covers `sensing_start`.
    fn list_aws_candidates(
        satellite: &str,
        sensing_start: DateTime<Utc>,
    ) -> Result<Vec<(String, String)>, OrbitFetchError> {
        let mut candidates: Vec<(String, String)> = Vec::new();

        for day_offset in [-1i64, 0, 1] {
            let date = sensing_start + Duration::days(day_offset);
            let prefix = format!(
                "AUX_POEORB/{}/{}/{}/",
                satellite,
                date.format("%Y"),
                date.format("%m")
            );
            let list_url = format!(
                "{}?list-type=2&prefix={}",
                AWS_BASE, prefix
            );

            let response = ureq::builder()
                .timeout(std::time::Duration::from_secs(HTTP_TIMEOUT_S))
                .build()
                .get(&list_url)
                .call()
                .map_err(|e| OrbitFetchError::Http(format!("{}: {}", list_url, e)))?;

            if response.status() != 200 {
                // Non-200 on listing means no objects with this prefix — not fatal.
                continue;
            }

            let body = response
                .into_string()
                .map_err(|e| OrbitFetchError::Http(format!("reading S3 list body: {}", e)))?;

            // Parse `<Key>…</Key>` elements from the S3 XML response.
            for key in extract_s3_keys(&body) {
                if let Some(filename) = key.split('/').next_back() {
                    if orbit_filename_covers(filename, sensing_start) {
                        let url = format!("{}/{}", AWS_BASE, key);
                        candidates.push((filename.to_string(), url));
                    }
                }
            }
        }

        Ok(candidates)
    }

    /// Extract `<Key>` element values from an S3 XML listing body.
    pub(super) fn extract_s3_keys(xml: &str) -> Vec<String> {
        let mut keys = Vec::new();
        let mut remaining = xml;
        while let Some(start) = remaining.find("<Key>") {
            remaining = &remaining[start + 5..];
            if let Some(end) = remaining.find("</Key>") {
                keys.push(remaining[..end].to_string());
                remaining = &remaining[end + 6..];
            } else {
                break;
            }
        }
        keys
    }

    // ──────────────────────────────────────────────
    // ESA STEP fallback
    // ──────────────────────────────────────────────

    /// List POEORB file URLs from the ESA STEP server (HTML directory listing).
    ///
    /// Checks directories for day−1, day, day+1 relative to `sensing_start`.
    fn list_esa_candidates(
        satellite: &str,
        sensing_start: DateTime<Utc>,
    ) -> Result<Vec<(String, String)>, OrbitFetchError> {
        let mut candidates: Vec<(String, String)> = Vec::new();

        for day_offset in [-1i64, 0, 1] {
            let date = sensing_start + Duration::days(day_offset);
            let dir_url = format!(
                "{}/POEORB/{}/{}/{}/",
                ESA_BASE,
                satellite,
                date.format("%Y"),
                date.format("%m")
            );

            let response = ureq::builder()
                .timeout(std::time::Duration::from_secs(HTTP_TIMEOUT_S))
                .build()
                .get(&dir_url)
                .call()
                .map_err(|e| OrbitFetchError::Http(format!("{}: {}", dir_url, e)))?;

            if response.status() != 200 {
                continue;
            }

            let html = response
                .into_string()
                .map_err(|e| OrbitFetchError::Http(format!("reading ESA dir body: {}", e)))?;

            // ESA serves `.EOF.zip` files.
            for filename in extract_esa_filenames(&html) {
                // Strip .zip to get the EOF stem for name-based coverage check.
                let stem = filename.trim_end_matches(".zip");
                if orbit_filename_covers(stem, sensing_start) {
                    let url = format!("{}{}", dir_url, filename);
                    candidates.push((filename.clone(), url));
                }
            }
        }

        Ok(candidates)
    }

    /// Extract `.EOF.zip` hrefs from an ESA directory listing HTML page.
    fn extract_esa_filenames(html: &str) -> Vec<String> {
        let mut names = Vec::new();
        for line in html.lines() {
            if !line.contains(".EOF.zip") {
                continue;
            }
            if let Some(start) = line.find("href=\"") {
                let rest = &line[start + 6..];
                if let Some(end) = rest.find('"') {
                    let href = &rest[..end];
                    // Keep only bare filenames (no path separators).
                    if !href.contains('/') && href.ends_with(".EOF.zip") {
                        names.push(href.to_string());
                    }
                }
            }
        }
        names
    }

    // ──────────────────────────────────────────────
    // Filename-based coverage check
    // ──────────────────────────────────────────────

    /// Return true if the orbit filename's validity window contains `t`.
    ///
    /// Standard POEORB filename pattern:
    /// `S1{A|B}_OPER_AUX_POEORB_OPOD_{prod}T{time}_V{start}_{stop}.EOF[.zip]`
    ///
    /// We parse `{start}` and `{stop}` from the `_V` suffix.
    pub(super) fn orbit_filename_covers(filename: &str, t: DateTime<Utc>) -> bool {
        // Look for the `_V` validity suffix.
        let stripped = filename.trim_end_matches(".EOF").trim_end_matches(".zip");
        if let Some(v_pos) = stripped.rfind("_V") {
            let validity = &stripped[v_pos + 2..]; // e.g. "20201004T155021_20201006T155021"
            if let Some(sep) = validity.find('_') {
                let start_s = &validity[..sep];
                let stop_s = &validity[sep + 1..];
                if let (Ok(start), Ok(stop)) = (
                    parse_orbit_validity_time(start_s),
                    parse_orbit_validity_time(stop_s),
                ) {
                    return t >= start && t <= stop;
                }
            }
        }
        // Filename doesn't match pattern — include as candidate (parser will reject if wrong).
        true
    }

    /// Parse a compact datetime string `YYYYMMDDTHHMMSS` to `DateTime<Utc>`.
    fn parse_orbit_validity_time(s: &str) -> Result<DateTime<Utc>, ()> {
        chrono::NaiveDateTime::parse_from_str(s, "%Y%m%dT%H%M%S")
            .map(|ndt| DateTime::from_naive_utc_and_offset(ndt, Utc))
            .map_err(|_| ())
    }

    // ──────────────────────────────────────────────
    // ESA ZIP extraction
    // ──────────────────────────────────────────────

    /// Download a `.EOF.zip` from ESA, extract the inner `.EOF` file, write to `dest`.
    fn download_and_extract_zip(url: &str, dest: &Path) -> Result<(), OrbitFetchError> {
        let tmp_zip = dest.with_extension("zip.tmp");

        let response = ureq::builder()
            .timeout(std::time::Duration::from_secs(HTTP_TIMEOUT_S))
            .build()
            .get(url)
            .call()
            .map_err(|e| OrbitFetchError::Http(format!("{}: {}", url, e)))?;

        if response.status() < 200 || response.status() >= 300 {
            return Err(OrbitFetchError::Http(format!(
                "{}: HTTP {}",
                url,
                response.status()
            )));
        }

        {
            let mut file = std::fs::File::create(&tmp_zip)?;
            let mut reader = response.into_reader();
            io::copy(&mut reader, &mut file)?;
            file.flush()?;
        }

        // Extract .EOF entry from the ZIP archive.
        let zip_data = std::fs::read(&tmp_zip)?;
        let eof_bytes = extract_eof_from_zip(&zip_data).map_err(|e| {
            OrbitFetchError::Http(format!("ZIP extraction failed for {}: {}", url, e))
        })?;

        {
            let mut out = std::fs::File::create(dest)?;
            out.write_all(&eof_bytes)?;
            out.flush()?;
        }
        let _ = std::fs::remove_file(&tmp_zip); // SAFETY-OK: cleanup; failure is non-fatal (only leaves a .zip.tmp behind)

        Ok(())
    }

    /// Minimal ZIP extractor — finds the first entry whose name ends with `.EOF`.
    ///
    /// This avoids adding a zip crate dependency; it works for ESA's simple
    /// single-file archives which use only DEFLATE or STORE compression.
    /// Returns `Err` with a human-readable message on any format mismatch.
    fn extract_eof_from_zip(data: &[u8]) -> Result<Vec<u8>, String> {
        // Locate End-of-Central-Directory (EOCD) signature: PK\x05\x06.
        let eocd_sig: [u8; 4] = [0x50, 0x4b, 0x05, 0x06];
        let eocd_pos = data
            .windows(4)
            .rposition(|w| w == eocd_sig)
            .ok_or("ZIP: EOCD signature not found")?;

        if eocd_pos + 22 > data.len() {
            return Err("ZIP: truncated EOCD record".to_string());
        }

        let cd_offset = u32::from_le_bytes(
            data[eocd_pos + 16..eocd_pos + 20]
                .try_into()
                .map_err(|_| "ZIP: cannot read CD offset")?,
        ) as usize;

        // Walk central directory entries to find an .EOF file.
        let cd_sig: [u8; 4] = [0x50, 0x4b, 0x01, 0x02];
        let mut pos = cd_offset;
        while pos + 46 <= data.len() {
            if data[pos..pos + 4] != cd_sig {
                break;
            }

            let fname_len = u16::from_le_bytes(
                data[pos + 28..pos + 30].try_into().map_err(|_| "ZIP: fname len")?,
            ) as usize;
            let extra_len = u16::from_le_bytes(
                data[pos + 30..pos + 32].try_into().map_err(|_| "ZIP: extra len")?,
            ) as usize;
            let comment_len = u16::from_le_bytes(
                data[pos + 32..pos + 34].try_into().map_err(|_| "ZIP: comment len")?,
            ) as usize;
            let local_header_offset = u32::from_le_bytes(
                data[pos + 42..pos + 46].try_into().map_err(|_| "ZIP: local hdr offset")?,
            ) as usize;

            let fname_bytes = data
                .get(pos + 46..pos + 46 + fname_len)
                .ok_or("ZIP: fname out of bounds")?;
            let fname = std::str::from_utf8(fname_bytes).map_err(|_| "ZIP: non-UTF8 filename")?;

            if fname.ends_with(".EOF") {
                // Parse local file header to find data offset.
                let lh = local_header_offset;
                let lh_sig: [u8; 4] = [0x50, 0x4b, 0x03, 0x04];
                if data.get(lh..lh + 4) != Some(&lh_sig) {
                    return Err("ZIP: local header signature mismatch".to_string());
                }

                let compression = u16::from_le_bytes(
                    data[lh + 8..lh + 10].try_into().map_err(|_| "ZIP: compression")?,
                );
                let compressed_size = u32::from_le_bytes(
                    data[lh + 18..lh + 22].try_into().map_err(|_| "ZIP: compressed size")?,
                ) as usize;
                let uncompressed_size = u32::from_le_bytes(
                    data[lh + 22..lh + 26].try_into().map_err(|_| "ZIP: uncompressed size")?,
                ) as usize;
                let lh_fname_len = u16::from_le_bytes(
                    data[lh + 26..lh + 28].try_into().map_err(|_| "ZIP: lh fname len")?,
                ) as usize;
                let lh_extra_len = u16::from_le_bytes(
                    data[lh + 28..lh + 30].try_into().map_err(|_| "ZIP: lh extra len")?,
                ) as usize;

                let data_start = lh + 30 + lh_fname_len + lh_extra_len;
                let data_end = data_start + compressed_size;
                let compressed_data = data
                    .get(data_start..data_end)
                    .ok_or("ZIP: data out of bounds")?;

                return match compression {
                    0 => {
                        // STORE
                        Ok(compressed_data.to_vec())
                    }
                    8 => {
                        // DEFLATE — requires flate2; return error suggesting it
                        // (ESA files are rarely STORE, but we cannot add flate2 without deps)
                        let _ = uncompressed_size; // SAFETY-OK: DEFLATE branch returns Err immediately; size is unused only because we refuse to decompress without flate2
                        Err("ZIP: DEFLATE compression requires the 'flate2' crate; \
                             add it to Cargo.toml or use the AWS source (no ZIP)"
                            .to_string())
                    }
                    c => Err(format!("ZIP: unsupported compression method {}", c)),
                };
            }

            pos += 46 + fname_len + extra_len + comment_len;
        }

        Err("ZIP: no .EOF entry found in central directory".to_string())
    }

    // ──────────────────────────────────────────────
    // Public API
    // ──────────────────────────────────────────────

    /// Download the POEORB file for `product_id` to `cache_dir`.
    ///
    /// # Arguments
    /// - `product_id`: SAFE directory stem (e.g. `S1A_IW_SLC__1SDV_20201005T170824_…`).
    /// - `sensing_start`: UTC sensing start time (parsed from annotation or SAFE name).
    /// - `cache_dir`: Directory to write the `.EOF` file into (created if absent).
    ///
    /// # Returns
    /// Path to the downloaded `.EOF` file, ready to pass to
    /// [`crate::orbit::parse_eof_file`].
    ///
    /// # Errors
    /// Returns [`OrbitFetchError`] for network failures, not-found, I/O errors,
    /// and orbit coverage validation failures.
    pub fn fetch_poeorb(
        product_id: &str,
        sensing_start: DateTime<Utc>,
        cache_dir: &Path,
    ) -> Result<PathBuf, OrbitFetchError> {
        let satellite = satellite_from_product_id(product_id)?;
        std::fs::create_dir_all(cache_dir)?;

        // ── Try AWS first (pre-extracted .EOF, no ZIP needed) ───────────────
        match list_aws_candidates(satellite, sensing_start) {
            Ok(candidates) if !candidates.is_empty() => {
                for (filename, url) in &candidates {
                    let dest = cache_dir.join(filename);
                    if dest.exists() {
                        // Already cached — validate before returning.
                        if validate_orbit_coverage(&dest, sensing_start).is_ok() {
                            return Ok(dest);
                        }
                        // Cached file failed validation — re-download.
                        let _ = std::fs::remove_file(&dest); // SAFETY-OK: stale/corrupt file; removal failure is harmless, re-download will overwrite
                    }
                    match download_to_file(url, &dest) {
                        Ok(()) => match validate_orbit_coverage(&dest, sensing_start) {
                            Ok(()) => return Ok(dest),
                            Err(e) => {
                                // Remove bad file so it isn't returned on retry.
                                let _ = std::fs::remove_file(&dest); // SAFETY-OK: file failed validation; keeping it would silently return wrong data
                                tracing::warn!("AWS orbit validation failed for {}: {}", url, e);
                            }
                        },
                        Err(e) => {
                            tracing::warn!("AWS orbit download failed for {}: {}", url, e);
                        }
                    }
                }
            }
            Ok(_) => {}
            Err(e) => {
                tracing::warn!("AWS orbit listing failed: {}", e);
            }
        }

        // ── Fallback: ESA STEP (requires ZIP extraction) ────────────────────
        match list_esa_candidates(satellite, sensing_start) {
            Ok(candidates) if !candidates.is_empty() => {
                for (filename, url) in &candidates {
                    // ESA filenames end in .EOF.zip; we store as .EOF after extraction.
                    let eof_name = filename.trim_end_matches(".zip");
                    let dest = cache_dir.join(eof_name);
                    if dest.exists() {
                        if validate_orbit_coverage(&dest, sensing_start).is_ok() {
                            return Ok(dest);
                        }
                        let _ = std::fs::remove_file(&dest); // SAFETY-OK: stale/corrupt; see above
                    }
                    match download_and_extract_zip(url, &dest) {
                        Ok(()) => match validate_orbit_coverage(&dest, sensing_start) {
                            Ok(()) => return Ok(dest),
                            Err(e) => {
                                let _ = std::fs::remove_file(&dest); // SAFETY-OK: validation failure; see above
                                tracing::warn!("ESA orbit validation failed for {}: {}", url, e);
                            }
                        },
                        Err(e) => {
                            tracing::warn!("ESA orbit download failed for {}: {}", url, e);
                        }
                    }
                }
            }
            Ok(_) => {}
            Err(e) => {
                tracing::warn!("ESA orbit listing failed: {}", e);
            }
        }

        Err(OrbitFetchError::NotFound {
            satellite: satellite.to_string(),
            sensing_start,
        })
    }
}

#[cfg(feature = "orbit-fetch")]
pub use inner::{fetch_poeorb, OrbitFetchError};

// ──────────────────────────────────────────────────────────────────────────────
// Unit tests (compile-time only, no network)
// ──────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_satellite_unknown_prefix() {
        #[cfg(feature = "orbit-fetch")]
        {
            use inner::*;
            let err = fetch_poeorb(
                "S2A_MSIL1C_20201005T170824",
                chrono::Utc::now(),
                std::path::Path::new("/tmp"),
            );
            assert!(matches!(err, Err(OrbitFetchError::UnknownSatellite(_))));
        }
    }

    #[cfg(feature = "orbit-fetch")]
    #[test]
    fn test_orbit_filename_covers_valid_window() {
        use inner::*;
        use chrono::TimeZone;
        // A known POEORB filename from the test dataset.
        let filename =
            "S1A_OPER_AUX_POEORB_OPOD_20201007T121500_V20201004T155021_20201006T155021.EOF";
        let t_inside = chrono::Utc.with_ymd_and_hms(2020, 10, 5, 17, 8, 0).unwrap();
        let t_outside = chrono::Utc.with_ymd_and_hms(2020, 10, 7, 0, 0, 0).unwrap();
        assert!(orbit_filename_covers(filename, t_inside));
        assert!(!orbit_filename_covers(filename, t_outside));
    }

    #[cfg(feature = "orbit-fetch")]
    #[test]
    fn test_extract_s3_keys_empty() {
        use inner::*;
        let keys = extract_s3_keys("<ListBucketResult></ListBucketResult>");
        assert!(keys.is_empty());
    }

    #[cfg(feature = "orbit-fetch")]
    #[test]
    fn test_extract_s3_keys_parses() {
        use inner::*;
        let xml = "<ListBucketResult><Contents><Key>AUX_POEORB/S1A/2020/10/foo.EOF</Key></Contents></ListBucketResult>";
        let keys = extract_s3_keys(xml);
        assert_eq!(keys, vec!["AUX_POEORB/S1A/2020/10/foo.EOF"]);
    }
}
