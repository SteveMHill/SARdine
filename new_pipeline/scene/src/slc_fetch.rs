//! Feature-gated Sentinel-1 IW SLC SAFE downloader.
//!
//! Enabled with `--features slc-fetch` (adds the `ureq` HTTP dependency).
//!
//! Entry point: [`fetch_slc`].
//!
//! # Source
//! Downloads from ASF DAAC datapool:
//! `https://datapool.asf.alaska.edu/SLC/{SA|SB}/{product_id}.zip`
//!
//! # Authentication
//! Requires an Earthdata Bearer token.  Obtain one at:
//! <https://urs.earthdata.nasa.gov/user_tokens>
//!
//! Pass via the `EARTHDATA_TOKEN` environment variable, or supply directly to
//! [`fetch_slc`].  Helper [`token_from_env`] reads the env var and returns a
//! typed error if it is absent.
//!
//! # Output
//! Extracts the `.zip` to `{output_dir}/{product_id}.SAFE/`.
//! If the SAFE directory already exists the download is skipped (idempotent).
//! If only the `.zip` already exists (e.g. a previous interrupted run) the
//! download is also skipped and only extraction is re-attempted.
//!
//! # Extraction
//! Uses the system `unzip` binary (universally available on Linux).
//! Flags used: `-q` (quiet), `-n` (never overwrite existing files),
//! `-d {output_dir}` (destination directory).
//!
//! # AGENTS.md rules
//! - No silent fallbacks: every failure is a typed [`SlcFetchError`] variant.
//! - No hardcoded Sentinel-1 constants.

#[cfg(feature = "slc-fetch")]
mod inner {
    use std::io::{self, Read, Write};
    use std::path::{Path, PathBuf};
    use std::process::Command;

    use thiserror::Error;

    // ASF DAAC datapool base URL (public, Earthdata auth required).
    const ASF_DATAPOOL: &str = "https://datapool.asf.alaska.edu";

    // HTTP timeout for downloads in seconds.  SLC products are 3–8 GB, so a
    // generous value is needed; 7200 s = 2 h is still a hard cap.
    const HTTP_TIMEOUT_S: u64 = 7200;

    // Streaming I/O buffer size.
    const BUF_SIZE: usize = 512 * 1024; // 512 KB

    // ─────────────────────────────────────────────────────────────────────────
    // Error type
    // ─────────────────────────────────────────────────────────────────────────

    /// Errors returned by [`fetch_slc`] and [`token_from_env`].
    #[derive(Debug, Error)]
    pub enum SlcFetchError {
        /// HTTP transport or unexpected status code.
        #[error("HTTP error fetching {url}: {detail}")]
        Http { url: String, detail: String },

        /// 401 Unauthorized — token is missing, expired, or invalid.
        #[error(
            "Authentication failed (HTTP 401) for {url}.\n\
             Ensure EARTHDATA_TOKEN is a valid Earthdata Bearer token.\n\
             Obtain one at https://urs.earthdata.nasa.gov/user_tokens"
        )]
        Unauthorized { url: String },

        /// 404 — product ID is not in the ASF catalogue.
        #[error("Product not found in ASF catalogue: {product_id}")]
        NotFound { product_id: String },

        /// OS-level I/O failure.
        #[error("I/O error: {0}")]
        Io(#[from] io::Error),

        /// Product ID prefix is not S1A or S1B.
        #[error("Unrecognised satellite prefix in product ID: {0}")]
        UnknownSatellite(String),

        /// `unzip` process failed or was not found.
        #[error("ZIP extraction failed: {0}")]
        ExtractionFailed(String),

        /// Extraction succeeded but the expected SAFE directory was not found.
        #[error(
            "SAFE directory not found after extraction.\n\
             Expected: {safe_dir}\n\
             Check that the ZIP contains a top-level .SAFE directory."
        )]
        SafeNotFound { safe_dir: PathBuf },

        /// `EARTHDATA_TOKEN` environment variable is not set.
        #[error(
            "EARTHDATA_TOKEN environment variable is not set.\n\
             Obtain a token at https://urs.earthdata.nasa.gov/user_tokens\n\
             and export it: export EARTHDATA_TOKEN=<token>"
        )]
        TokenNotSet,
    }

    // ─────────────────────────────────────────────────────────────────────────
    // Internal helpers
    // ─────────────────────────────────────────────────────────────────────────

    /// Map the S1A/S1B prefix to the ASF datapool two-letter satellite code.
    fn satellite_code(product_id: &str) -> Result<&'static str, SlcFetchError> {
        if product_id.starts_with("S1A") {
            Ok("SA")
        } else if product_id.starts_with("S1B") {
            Ok("SB")
        } else {
            Err(SlcFetchError::UnknownSatellite(
                product_id[..product_id.len().min(6)].to_string(),
            ))
        }
    }

    /// Derive the ASF product type segment from the product ID.
    ///
    /// Returns `"SLC"` unless the ID contains `"_GRD"`.
    fn product_type(product_id: &str) -> &'static str {
        if product_id.contains("_GRD") {
            "GRD"
        } else {
            "SLC"
        }
    }

    /// Build the ASF datapool download URL for a product ID.
    ///
    /// Format: `https://datapool.asf.alaska.edu/{TYPE}/{SAT}/{PRODUCT_ID}.zip`
    fn build_url(product_id: &str) -> Result<String, SlcFetchError> {
        let sat = satellite_code(product_id)?;
        let ptype = product_type(product_id);
        Ok(format!(
            "{}/{}/{}/{}.zip",
            ASF_DATAPOOL, ptype, sat, product_id
        ))
    }

    // Maximum number of redirects to follow manually.
    const MAX_REDIRECTS: usize = 10;

    /// Resolve a `Location` header value against the current request URL.
    fn resolve_redirect_url(current: &str, location: &str) -> String {
        if location.starts_with("http://") || location.starts_with("https://") {
            location.to_string()
        } else if location.starts_with('/') {
            let scheme_end = current.find("://").map(|i| i + 3).unwrap_or(0); // SAFETY-OK: no "://" means scheme_end=0; produces a usable URL that will succeed or fail explicitly on the next hop
            let host_end = current[scheme_end..]
                .find('/')
                .map(|i| scheme_end + i)
                .unwrap_or(current.len()); // SAFETY-OK: host-only URL means host_end=full len; appending abs-path location produces a valid absolute URL
            format!("{}{}", &current[..host_end], location)
        } else {
            let base = current.rfind('/').map(|i| i + 1).unwrap_or(current.len()); // SAFETY-OK: no '/' in URL (host-only) means base=len; appending relative location produces a valid URL that fails explicitly if wrong
            format!("{}{}", &current[..base], location)
        }
    }

    /// Perform a GET request, manually following redirects.
    ///
    /// `ureq` strips `Authorization` on cross-host redirects (security default).
    /// ASF datapool uses a 3-hop chain:
    ///   1. `datapool.asf.alaska.edu` → 307 → `sentinel1.asf.alaska.edu`
    ///   2. `sentinel1.asf.alaska.edu` + Bearer → 303 → pre-signed S3/CloudFront URL
    ///   3. S3/CloudFront pre-signed URL → data (Bearer must NOT be sent here)
    ///
    /// RFC 7231 §6.4.4: 303 See Other = fresh GET at new URI without credentials.
    /// S3 pre-signed URLs return 400 if an Authorization header accompanies the
    /// query-parameter signature.  Auth is therefore stripped on 303 redirects.
    fn get_with_auth(
        agent: &ureq::Agent,
        start_url: &str,
        token: &str,
        product_id: &str,
    ) -> Result<ureq::Response, SlcFetchError> {
        let auth = format!("Bearer {}", token);
        let mut url = start_url.to_string();
        let mut send_auth = true;
        for _ in 0..MAX_REDIRECTS {
            let req = agent.get(&url).set("User-Agent", "sardine-scene/0.1");
            let req = if send_auth { req.set("Authorization", &auth) } else { req };
            let resp = req.call().map_err(|e| match e {
                ureq::Error::Status(401, _) => SlcFetchError::Unauthorized {
                    url: url.clone(),
                },
                ureq::Error::Status(404, _) => SlcFetchError::NotFound {
                    product_id: product_id.to_string(),
                },
                other => SlcFetchError::Http {
                    url: url.clone(),
                    detail: other.to_string(),
                },
            })?;
            let status = resp.status();
            match status {
                301 | 302 | 307 | 308 => {
                    let location =
                        resp.header("Location").ok_or_else(|| SlcFetchError::Http {
                            url: url.clone(),
                            detail: format!("redirect {} without Location header", status),
                        })?;
                    url = resolve_redirect_url(&url, location);
                }
                303 => {
                    // 303 See Other — strip auth for the next hop (S3 pre-signed URL).
                    send_auth = false;
                    let location =
                        resp.header("Location").ok_or_else(|| SlcFetchError::Http {
                            url: url.clone(),
                            detail: "redirect 303 without Location header".to_string(),
                        })?;
                    url = resolve_redirect_url(&url, location);
                }
                _ => return Ok(resp),
            }
        }
        Err(SlcFetchError::Http {
            url: start_url.to_string(),
            detail: format!("exceeded {} redirects", MAX_REDIRECTS),
        })
    }

    /// Stream a URL to a file, authenticating with a Bearer token.
    ///
    /// Downloads to `{dest}.part`, then atomically renames to `dest`.
    /// Any existing `.part` file is removed before the download begins.
    fn download_to_file(
        url: &str,
        dest: &Path,
        product_id: &str,
        token: &str,
    ) -> Result<(), SlcFetchError> {
        let part = {
            let mut p = dest.as_os_str().to_owned();
            p.push(".part");
            PathBuf::from(p)
        };

        // Remove any stale partial file from a previous interrupted run.
        if part.exists() {
            std::fs::remove_file(&part)?;
        }

        let agent = ureq::builder()
            .timeout(std::time::Duration::from_secs(HTTP_TIMEOUT_S))
            // Disable ureq's automatic redirect following so we handle it
            // manually via get_with_auth (which re-attaches the auth header).
            .redirects(0)
            .build();

        let response = get_with_auth(&agent, url, token, product_id)?;

        let mut out = std::fs::File::create(&part)?;
        let mut reader = response.into_reader();
        let mut buf = vec![0u8; BUF_SIZE];
        loop {
            let n = reader.read(&mut buf).map_err(SlcFetchError::Io)?;
            if n == 0 {
                break;
            }
            out.write_all(&buf[..n])?;
        }
        out.flush()?;
        drop(out);

        std::fs::rename(&part, dest)?;
        Ok(())
    }

    /// Extract a ZIP archive to `output_dir` using the system `unzip` binary.
    ///
    /// Uses `-q` (quiet) and `-n` (never overwrite) flags.
    fn extract_zip(zip_path: &Path, output_dir: &Path) -> Result<(), SlcFetchError> {
        let status = Command::new("unzip")
            .args(["-q", "-n"])
            .arg(zip_path)
            .arg("-d")
            .arg(output_dir)
            .status()
            .map_err(|e| {
                SlcFetchError::ExtractionFailed(format!(
                    "failed to run unzip (is it installed?): {}",
                    e
                ))
            })?;

        if !status.success() {
            return Err(SlcFetchError::ExtractionFailed(format!(
                "unzip exited with status {}",
                status
            )));
        }
        Ok(())
    }

    // ─────────────────────────────────────────────────────────────────────────
    // Public API
    // ─────────────────────────────────────────────────────────────────────────

    /// Read the Earthdata Bearer token from the `EARTHDATA_TOKEN` environment variable.
    ///
    /// Returns [`SlcFetchError::TokenNotSet`] if the variable is absent or empty.
    pub fn token_from_env() -> Result<String, SlcFetchError> {
        let val = std::env::var("EARTHDATA_TOKEN")
            .map_err(|_| SlcFetchError::TokenNotSet)?;
        if val.trim().is_empty() {
            return Err(SlcFetchError::TokenNotSet);
        }
        Ok(val)
    }

    /// Download and extract a Sentinel-1 IW SLC SAFE product from ASF datapool.
    ///
    /// # Arguments
    /// - `product_id`: the bare product ID, e.g.
    ///   `S1B_IW_SLC__1SDV_20190123T053348_20190123T053415_014617_01B3D4_E833`
    ///   (without `.SAFE` or `.zip` suffix).
    /// - `output_dir`: directory to place the extracted `.SAFE` into.
    ///   Created if it does not exist.
    /// - `token`: Earthdata Bearer token.  Use [`token_from_env`] to read from
    ///   the `EARTHDATA_TOKEN` environment variable.
    ///
    /// # Returns
    /// Path to the extracted `{product_id}.SAFE` directory.
    ///
    /// # Idempotency
    /// - If `{output_dir}/{product_id}.SAFE` already exists, returns immediately.
    /// - If `{output_dir}/{product_id}.zip` already exists (partial previous run),
    ///   skips the download and retries extraction.
    pub fn fetch_slc(
        product_id: &str,
        output_dir: &Path,
        token: &str,
    ) -> Result<PathBuf, SlcFetchError> {
        // Validate token before any network I/O so the error message is actionable.
        if token.trim().is_empty() {
            return Err(SlcFetchError::TokenNotSet);
        }

        let safe_dir = output_dir.join(format!("{}.SAFE", product_id));

        // Fast path: SAFE already exists.
        if safe_dir.exists() {
            tracing::info!(
                "{} already exists, skipping download",
                safe_dir.display()
            );
            return Ok(safe_dir);
        }

        std::fs::create_dir_all(output_dir)?;

        let zip_path = output_dir.join(format!("{}.zip", product_id));

        if zip_path.exists() {
            tracing::info!(
                "ZIP already present at {}, skipping download",
                zip_path.display()
            );
        } else {
            let url = build_url(product_id)?;
            tracing::info!("downloading {} …", url);
            download_to_file(&url, &zip_path, product_id, token)?;
            tracing::info!("download complete → {}", zip_path.display());
        }

        tracing::info!(
            "extracting {} → {} …",
            zip_path.display(),
            output_dir.display()
        );
        extract_zip(&zip_path, output_dir)?;

        if !safe_dir.exists() {
            return Err(SlcFetchError::SafeNotFound { safe_dir });
        }

        tracing::info!("SAFE ready at {}", safe_dir.display());
        Ok(safe_dir)
    }

    // ─────────────────────────────────────────────────────────────────────────
    // Tests
    // ─────────────────────────────────────────────────────────────────────────

    #[cfg(test)]
    mod tests {
        use super::*;

        #[test]
        fn test_build_url_s1b_slc() {
            let url = build_url(
                "S1B_IW_SLC__1SDV_20190123T053348_20190123T053415_014617_01B3D4_E833",
            )
            .unwrap();
            assert_eq!(
                url,
                "https://datapool.asf.alaska.edu/SLC/SB/\
                 S1B_IW_SLC__1SDV_20190123T053348_20190123T053415_014617_01B3D4_E833.zip"
            );
        }

        #[test]
        fn test_build_url_s1a_slc() {
            let url = build_url(
                "S1A_IW_SLC__1SDV_20201005T170824_20201005T170851_034664_04098A_1E66",
            )
            .unwrap();
            assert_eq!(
                url,
                "https://datapool.asf.alaska.edu/SLC/SA/\
                 S1A_IW_SLC__1SDV_20201005T170824_20201005T170851_034664_04098A_1E66.zip"
            );
        }

        #[test]
        fn test_build_url_unknown_satellite() {
            let err = build_url("S2A_MSIL1C_20190123").unwrap_err();
            assert!(
                matches!(err, SlcFetchError::UnknownSatellite(_)),
                "expected UnknownSatellite, got {:?}",
                err
            );
        }

        #[test]
        fn test_build_url_grd() {
            let url = build_url(
                "S1A_IW_GRDH_1SDV_20201005T170824_20201005T170849_034664_04098A_1234",
            )
            .unwrap();
            assert!(url.contains("/GRD/"), "expected GRD in URL: {}", url);
        }

        #[test]
        fn test_token_from_env_absent() {
            // Only safe to run if the env var is not accidentally set.
            if std::env::var("EARTHDATA_TOKEN").is_ok() {
                return; // Skip: env var set in this environment.
            }
            let err = token_from_env().unwrap_err();
            assert!(
                matches!(err, SlcFetchError::TokenNotSet),
                "expected TokenNotSet, got {:?}",
                err
            );
        }

        #[test]
        fn test_safe_dir_already_exists_is_idempotent() {
            let dir = tempfile::tempdir().unwrap();
            let product_id = "S1B_IW_SLC__1SDV_20190123T053348_20190123T053415_014617_01B3D4_E833";
            let safe = dir.path().join(format!("{}.SAFE", product_id));
            std::fs::create_dir(&safe).unwrap();

            // Should return without error even though token is bogus.
            let result = fetch_slc(product_id, dir.path(), "fake-token");
            assert!(result.is_ok(), "{:?}", result);
            assert_eq!(result.unwrap(), safe);
        }
    }
}

#[cfg(feature = "slc-fetch")]
pub use inner::{fetch_slc, token_from_env, SlcFetchError};
