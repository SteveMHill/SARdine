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
    use std::os::unix::fs::FileExt;
    use std::path::{Path, PathBuf};
    use std::process::Command;
    use std::sync::{
        atomic::{AtomicU64, Ordering},
        Arc,
    };

    use thiserror::Error;
    // ASF Search API endpoint (public, no auth required).
    const ASF_SEARCH_URL: &str = "https://api.daac.asf.alaska.edu/services/search/param";

    // Earthdata Login token endpoint — idempotent: returns existing token or
    // creates new one.  Using /token instead would fail with 403 if the user
    // already holds two tokens.
    const EDL_TOKEN_URL: &str =
        "https://urs.earthdata.nasa.gov/api/users/find_or_create_token";

    // HTTP timeout for downloads in seconds.  SLC products are 3–8 GB, so a
    // generous value is needed; 7200 s = 2 h is still a hard cap.
    const HTTP_TIMEOUT_S: u64 = 7200;

    // HTTP timeout for search and token requests.
    const SEARCH_TIMEOUT_S: u64 = 35;

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

        /// Authentication failed when acquiring a bearer token (username/password path).
        #[error("Token acquisition failed: {detail}")]
        AuthFailed { detail: String },

        /// MD5 checksum of the downloaded ZIP does not match the catalogue value.
        #[error(
            "MD5 mismatch for {product_id}: expected {expected}, got {actual}.\n\
             The download may be corrupt.  Delete the .zip and retry."
        )]
        ChecksumMismatch { product_id: String, expected: String, actual: String },

        /// ASF search API returned an HTTP error.
        #[error("ASF search failed for {url}: {detail}")]
        SearchFailed { url: String, detail: String },

        /// ASF search response could not be parsed.
        #[error("Failed to parse ASF search response: {detail}")]
        SearchParseError { detail: String },

        /// ASF search API returned HTTP 429 — rate limit exceeded.
        #[error(
            "ASF search rate-limited (HTTP 429) for {url}.\n\
             Retry after {retry_after_secs:?} seconds."
        )]
        RateLimited { url: String, retry_after_secs: Option<u64> },

        /// No MD5 checksum available in catalogue; cannot verify download.
        #[error("No MD5 checksum in ASF catalogue for {product_id}")]
        NoChecksum { product_id: String },
    }

    // ─────────────────────────────────────────────────────────────────────────
    // Credential types
    // ─────────────────────────────────────────────────────────────────────────

    /// Earthdata authentication credentials.
    #[derive(Debug)]
    pub enum SlcCredentials {
        /// Pre-obtained bearer token (from `EARTHDATA_TOKEN` env var or `--token` flag).
        BearerToken(String),
        /// Username + password — a bearer token is acquired on first use via
        /// [`acquire_token`].
        UsernamePassword { username: String, password: String },
    }

    // ─────────────────────────────────────────────────────────────────────────
    // Search parameter / result types
    // ─────────────────────────────────────────────────────────────────────────

    /// Parameters for [`search_slc`].
    pub struct SlcSearchParams {
        /// Acquisition start (UTC, inclusive).
        pub start: chrono::DateTime<chrono::Utc>,
        /// Acquisition end (UTC, inclusive).
        pub end: chrono::DateTime<chrono::Utc>,
        /// Optional bounding box `[min_lon, min_lat, max_lon, max_lat]` (WGS-84).
        pub bbox: Option<[f64; 4]>,
        /// Beam mode filter, e.g. `"IW"`, `"EW"`, `"SM"`.
        pub beam_mode: Option<String>,
        /// Polarization filter, e.g. `"VV+VH"`, `"VV"`.
        pub polarization: Option<String>,
        /// Flight direction filter: `"ASCENDING"` or `"DESCENDING"`.
        pub orbit_direction: Option<String>,
        /// Relative orbit (path) number filter.
        pub relative_orbit: Option<u32>,
        /// Maximum number of results (default 100; hard-capped at 250 to stay
        /// within the ASF 30-second query timeout).
        pub max_results: usize,
    }

    impl Default for SlcSearchParams {
        fn default() -> Self {
            Self {
                start: chrono::DateTime::UNIX_EPOCH,
                end: chrono::Utc::now(),
                bbox: None,
                beam_mode: None,
                polarization: None,
                orbit_direction: None,
                relative_orbit: None,
                max_results: 100,
            }
        }
    }

    /// One result record returned by [`search_slc`].
    ///
    /// Field names map directly to the ASF GeoJSON `properties` keys
    /// (as documented in ASFProduct.py from the asf_search package).
    #[derive(Debug)]
    pub struct SlcSearchResult {
        /// Scene name / granule identifier (e.g.
        /// `S1B_IW_SLC__1SDV_20190123T053348_20190123T053415_014617_01B3D4_E833`).
        pub scene_name: String,
        /// Direct download URL for the product ZIP.
        pub url: String,
        /// File size in bytes.
        pub bytes: u64,
        /// MD5 hex digest of the product ZIP, or `None` if not in the catalogue.
        pub md5sum: Option<String>,
        /// Sensing start time (UTC).
        pub start_time: chrono::DateTime<chrono::Utc>,
        /// Sensing stop time (UTC).
        pub stop_time: chrono::DateTime<chrono::Utc>,
        /// Platform string (e.g. `"Sentinel-1B"`).
        pub platform: String,
        /// Beam mode type (e.g. `"IW"`).
        pub beam_mode: String,
        /// Polarization (e.g. `"VV+VH"`).
        pub polarization: String,
        /// Flight direction (`"ASCENDING"` or `"DESCENDING"`).
        pub orbit_direction: String,
        /// Relative orbit (path) number.
        pub relative_orbit: Option<u32>,
        /// ESA frame number.
        pub frame_number: Option<u32>,
    }

    /// Progress callback type for [`fetch_slc_result`] and [`download_to_file_parallel`].
    ///
    /// Arguments: `(bytes_downloaded, total_bytes_or_none)`.
    /// Called after every 512 KB chunk from any worker thread.
    ///
    /// The callback must be `Send + Sync` because parallel downloads invoke it
    /// from multiple threads concurrently.  Wrap a closure in [`Arc::new`].
    pub type ProgressFn = Arc<dyn Fn(u64, Option<u64>) + Send + Sync + 'static>;

    // ─────────────────────────────────────────────────────────────────────────
    // Auth helpers
    // ─────────────────────────────────────────────────────────────────────────

    /// Minimal base64 encoder for constructing HTTP Basic auth headers.
    /// We do not add a `base64` crate dependency; the single use-site is here.
    fn base64_encode(input: &[u8]) -> String {
        const CHARS: &[u8; 64] =
            b"ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";
        let mut out = String::with_capacity((input.len() + 2) / 3 * 4);
        let mut i = 0;
        while i + 2 < input.len() {
            let b = ((input[i] as u32) << 16)
                | ((input[i + 1] as u32) << 8)
                | (input[i + 2] as u32);
            out.push(CHARS[((b >> 18) & 0x3f) as usize] as char);
            out.push(CHARS[((b >> 12) & 0x3f) as usize] as char);
            out.push(CHARS[((b >> 6) & 0x3f) as usize] as char);
            out.push(CHARS[(b & 0x3f) as usize] as char);
            i += 3;
        }
        match input.len() - i {
            1 => {
                let b = (input[i] as u32) << 16;
                out.push(CHARS[((b >> 18) & 0x3f) as usize] as char);
                out.push(CHARS[((b >> 12) & 0x3f) as usize] as char);
                out.push('=');
                out.push('=');
            }
            2 => {
                let b = ((input[i] as u32) << 16) | ((input[i + 1] as u32) << 8);
                out.push(CHARS[((b >> 18) & 0x3f) as usize] as char);
                out.push(CHARS[((b >> 12) & 0x3f) as usize] as char);
                out.push(CHARS[((b >> 6) & 0x3f) as usize] as char);
                out.push('=');
            }
            _ => {} // 0 — no padding needed
        }
        out
    }

    /// Acquire an Earthdata bearer token using username + password.
    ///
    /// Uses `POST /api/users/find_or_create_token` which is idempotent:
    /// it returns the user's existing token or creates a new one if none exists.
    /// This avoids the HTTP 403 `max_token_limit` error returned by
    /// `POST /api/users/token` when the user already holds two tokens.
    pub fn acquire_token(
        username: &str,
        password: &str,
    ) -> Result<String, SlcFetchError> {
        let agent = ureq::builder()
            .timeout(std::time::Duration::from_secs(SEARCH_TIMEOUT_S))
            .build();

        // ureq 2.x does not expose a `.auth()` convenience method.
        // Build the Basic auth header manually (RFC 7617).
        let credentials = base64_encode(format!("{}:{}", username, password).as_bytes());
        let auth_header = format!("Basic {}", credentials);

        let resp = agent
            .post(EDL_TOKEN_URL)
            .set("User-Agent", "sardine-scene/0.1")
            .set("Authorization", &auth_header)
            .call()
            .map_err(|e| match e {
                ureq::Error::Status(401, _) => SlcFetchError::AuthFailed {
                    detail: "invalid Earthdata credentials (HTTP 401)".to_string(),
                },
                other => SlcFetchError::AuthFailed { detail: other.to_string() },
            })?;

        let body: serde_json::Value = resp.into_json().map_err(|e| {
            SlcFetchError::AuthFailed {
                detail: format!("failed to parse token response: {}", e),
            }
        })?;

        body["access_token"]
            .as_str()
            .filter(|s| !s.is_empty())
            .map(|s| s.to_string())
            .ok_or_else(|| SlcFetchError::AuthFailed {
                detail: "token response missing 'access_token' field".to_string(),
            })
    }

    /// Read Earthdata credentials from the environment.
    ///
    /// Tries `EARTHDATA_TOKEN` first (→ [`SlcCredentials::BearerToken`]).
    /// Falls back to `EARTHDATA_USERNAME` + `EARTHDATA_PASSWORD`
    /// (→ [`SlcCredentials::UsernamePassword`]).
    /// Returns [`SlcFetchError::TokenNotSet`] if neither is available.
    pub fn credentials_from_env() -> Result<SlcCredentials, SlcFetchError> {
        if let Ok(tok) = std::env::var("EARTHDATA_TOKEN") {
            if !tok.trim().is_empty() {
                return Ok(SlcCredentials::BearerToken(tok));
            }
        }
        let user = std::env::var("EARTHDATA_USERNAME").ok();
        let pass = std::env::var("EARTHDATA_PASSWORD").ok();
        match (user, pass) {
            (Some(u), Some(p)) if !u.trim().is_empty() && !p.trim().is_empty() => {
                Ok(SlcCredentials::UsernamePassword {
                    username: u,
                    password: p,
                })
            }
            _ => Err(SlcFetchError::TokenNotSet),
        }
    }

    /// Resolve credentials to a bearer token string.
    fn resolve_token(creds: &SlcCredentials) -> Result<String, SlcFetchError> {
        match creds {
            SlcCredentials::BearerToken(t) => Ok(t.clone()),
            SlcCredentials::UsernamePassword { username, password } => {
                acquire_token(username, password)
            }
        }
    }

    // ─────────────────────────────────────────────────────────────────────────
    // Search
    // ─────────────────────────────────────────────────────────────────────────

    /// Search the ASF DAAC catalogue for Sentinel-1 IW SLC products.
    ///
    /// The ASF search endpoint is publicly accessible — no authentication
    /// is required.  Auth is only needed for the download step.
    ///
    /// # Errors
    /// - [`SlcFetchError::RateLimited`] on HTTP 429 (> 250 requests/minute).
    /// - [`SlcFetchError::SearchFailed`] on other HTTP errors.
    /// - [`SlcFetchError::SearchParseError`] if the response cannot be parsed.
    pub fn search_slc(
        params: &SlcSearchParams,
    ) -> Result<Vec<SlcSearchResult>, SlcFetchError> {
        let url = build_search_url(params);
        let agent = ureq::builder()
            .timeout(std::time::Duration::from_secs(SEARCH_TIMEOUT_S))
            .build();

        let resp = agent
            .get(&url)
            .set("User-Agent", "sardine-scene/0.1")
            .call()
            .map_err(|e| match e {
                ureq::Error::Status(429, ref r) => SlcFetchError::RateLimited {
                    url: url.clone(),
                    retry_after_secs: r
                        .header("Retry-After")
                        .and_then(|v| v.parse::<u64>().ok()), // SAFETY-OK: None means no Retry-After header; caller handles None gracefully
                },
                other => SlcFetchError::SearchFailed {
                    url: url.clone(),
                    detail: other.to_string(),
                },
            })?;

        let body = resp.into_string().map_err(|e| SlcFetchError::SearchFailed {
            url: url.clone(),
            detail: e.to_string(),
        })?;

        parse_geojson_results(&body)
    }

    /// Build the ASF search API URL for the given parameters.
    fn build_search_url(params: &SlcSearchParams) -> String {
        use std::fmt::Write as _;

        let mut q = format!(
            "{}?dataset=SENTINEL-1&processingLevel=SLC&output=geojson\
             &start={}&end={}&maxResults={}",
            ASF_SEARCH_URL,
            params.start.format("%Y-%m-%dT%H:%M:%SZ"),
            params.end.format("%Y-%m-%dT%H:%M:%SZ"),
            params.max_results.min(250), // ASF hard cap to stay under 30s timeout
        );

        if let Some([min_lon, min_lat, max_lon, max_lat]) = params.bbox {
            // Closed WKT polygon: first and last vertex must be identical.
            let wkt = format!(
                "POLYGON(({min_lon} {min_lat},{max_lon} {min_lat},\
                 {max_lon} {max_lat},{min_lon} {max_lat},{min_lon} {min_lat}))"
            );
            // URL-encode spaces and parentheses (ASF docs).
            let encoded = wkt
                .replace('(', "%28")
                .replace(')', "%29")
                .replace(' ', "+");
            let _ = write!(q, "&intersectsWith={}", encoded); // SAFETY-OK: write! to String is infallible
        }
        if let Some(bm) = &params.beam_mode {
            let _ = write!(q, "&beamMode={}", bm); // SAFETY-OK: write! to String is infallible
        }
        if let Some(pol) = &params.polarization {
            let _ = write!(q, "&polarization={}", pol.replace('+', "%2B")); // SAFETY-OK: write! to String is infallible
        }
        if let Some(dir) = &params.orbit_direction {
            let _ = write!(q, "&flightDirection={}", dir); // SAFETY-OK: write! to String is infallible
        }
        if let Some(ro) = params.relative_orbit {
            let _ = write!(q, "&relativeOrbit={}", ro); // SAFETY-OK: write! to String is infallible
        }
        q
    }

    /// Parse an ASF GeoJSON FeatureCollection response body into a list of
    /// [`SlcSearchResult`] structs.
    ///
    /// An empty `features` array is a valid result (no matches) — not an error.
    fn parse_geojson_results(body: &str) -> Result<Vec<SlcSearchResult>, SlcFetchError> {
        let v: serde_json::Value =
            serde_json::from_str(body).map_err(|e| SlcFetchError::SearchParseError {
                detail: format!("invalid JSON: {}", e),
            })?;

        let features = v["features"].as_array().ok_or_else(|| {
            SlcFetchError::SearchParseError {
                detail: "missing 'features' array in GeoJSON response".to_string(),
            }
        })?;

        let mut results = Vec::with_capacity(features.len());
        for (i, feat) in features.iter().enumerate() {
            let p = &feat["properties"];
            results.push(parse_one_feature(p, i)?);
        }
        Ok(results)
    }

    /// Parse the `properties` object of one GeoJSON feature.
    fn parse_one_feature(
        p: &serde_json::Value,
        idx: usize,
    ) -> Result<SlcSearchResult, SlcFetchError> {
        let ctx = |field: &str| SlcFetchError::SearchParseError {
            detail: format!("feature[{}]: missing or null '{}' field", idx, field),
        };
        let ctx_parse = |field: &str, detail: String| SlcFetchError::SearchParseError {
            detail: format!("feature[{}]: cannot parse '{}': {}", idx, field, detail),
        };

        let scene_name = p["sceneName"]
            .as_str()
            .ok_or_else(|| ctx("sceneName"))?
            .to_string();

        let url = p["url"]
            .as_str()
            .ok_or_else(|| ctx("url"))?
            .to_string();

        let bytes = p["bytes"]
            .as_f64()
            .ok_or_else(|| ctx("bytes"))? as u64;

        // md5sum is optional — not all ASF products carry it.
        let md5sum = p["md5sum"].as_str().map(|s| s.to_lowercase());

        let start_time = parse_asf_datetime(
            p["startTime"].as_str().ok_or_else(|| ctx("startTime"))?,
        )
        .map_err(|e| ctx_parse("startTime", e))?;

        let stop_time = parse_asf_datetime(
            p["stopTime"].as_str().ok_or_else(|| ctx("stopTime"))?,
        )
        .map_err(|e| ctx_parse("stopTime", e))?;

        let platform = p["platform"]
            .as_str()
            .ok_or_else(|| ctx("platform"))?
            .to_string();

        let beam_mode = p["beamModeType"]
            .as_str()
            .ok_or_else(|| ctx("beamModeType"))?
            .to_string();

        let polarization = p["polarization"]
            .as_str()
            .ok_or_else(|| ctx("polarization"))?
            .to_string();

        let orbit_direction = p["flightDirection"]
            .as_str()
            .ok_or_else(|| ctx("flightDirection"))?
            .to_string();

        let relative_orbit = p["pathNumber"].as_u64().map(|v| v as u32);
        let frame_number = p["frameNumber"].as_u64().map(|v| v as u32);

        Ok(SlcSearchResult {
            scene_name,
            url,
            bytes,
            md5sum,
            start_time,
            stop_time,
            platform,
            beam_mode,
            polarization,
            orbit_direction,
            relative_orbit,
            frame_number,
        })
    }

    /// Parse an ASF datetime string to a UTC [`chrono::DateTime`].
    ///
    /// ASF historically returned bare naive datetimes (`"2019-01-23T05:33:48.000000"`,
    /// no timezone suffix).  As of 2025 the API now returns RFC3339 strings with a
    /// `Z` suffix and no fractional seconds (`"2019-01-23T05:33:48Z"`).  Both
    /// formats are accepted; all values are treated as UTC.
    fn parse_asf_datetime(s: &str) -> Result<chrono::DateTime<chrono::Utc>, String> {
        // Try RFC3339 first (handles 'Z' and '+00:00' suffixes).
        if let Ok(dt) = chrono::DateTime::parse_from_rfc3339(s) {
            return Ok(dt.with_timezone(&chrono::Utc));
        }
        // Fall back to the legacy bare-naive format (no timezone suffix).
        chrono::NaiveDateTime::parse_from_str(s, "%Y-%m-%dT%H:%M:%S%.f")
            .map(|ndt| ndt.and_utc())
            .map_err(|e| format!("{e} (input: {s:?})"))
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
    ///
    /// The `Range` header is carried through all redirects (including 303) —
    /// S3 supports Range on pre-signed URLs.
    fn get_with_auth(
        agent: &ureq::Agent,
        start_url: &str,
        token: &str,
        product_id: &str,
        range_start: Option<u64>,
    ) -> Result<(String, ureq::Response), SlcFetchError> {
        let auth = format!("Bearer {}", token);
        let range_header = range_start.map(|o| format!("bytes={}-", o));
        let mut url = start_url.to_string();
        let mut send_auth = true;
        for _ in 0..MAX_REDIRECTS {
            let req = agent.get(&url).set("User-Agent", "sardine-scene/0.1");
            let req = if send_auth { req.set("Authorization", &auth) } else { req };
            let req = if let Some(ref rh) = range_header {
                req.set("Range", rh)
            } else {
                req
            };
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
                _ => return Ok((url, resp)),
            }
        }
        Err(SlcFetchError::Http {
            url: start_url.to_string(),
            detail: format!("exceeded {} redirects", MAX_REDIRECTS),
        })
    }

    /// Verify the MD5 checksum of a file against an expected hex digest.
    ///
    /// Reads the file in 512 KB chunks to avoid loading it entirely into memory.
    /// Returns [`SlcFetchError::ChecksumMismatch`] if the digests differ.
    fn verify_md5(
        path: &Path,
        expected_hex: &str,
        product_id: &str,
    ) -> Result<(), SlcFetchError> {
        use std::io::Read;
        let mut ctx = md5::Context::new();
        let mut f = std::fs::File::open(path)?;
        let mut buf = vec![0u8; BUF_SIZE];
        loop {
            let n = f.read(&mut buf).map_err(SlcFetchError::Io)?;
            if n == 0 {
                break;
            }
            ctx.consume(&buf[..n]);
        }
        let actual = format!("{:x}", ctx.compute());
        let expected = expected_hex.to_lowercase();
        if actual != expected {
            return Err(SlcFetchError::ChecksumMismatch {
                product_id: product_id.to_string(),
                expected,
                actual,
            });
        }
        Ok(())
    }

    /// Stream a URL to a file, authenticating with a Bearer token.
    ///
    /// Supports HTTP Range resume: if `{dest}.part` already exists, the download
    /// resumes from the existing byte offset.  If the server does not honour the
    /// `Range` header (returns 200 instead of 206), the partial file is discarded
    /// and the download restarts from the beginning.
    ///
    /// After each 512 KB chunk, `progress(bytes_done, bytes_total)` is called
    /// if a callback is provided.  `bytes_total` is `None` when the server omits
    /// `Content-Length`.
    ///
    /// On success, atomically renames `{dest}.part` to `dest`.
    fn download_to_file(
        url: &str,
        dest: &Path,
        product_id: &str,
        token: &str,
        progress: Option<&ProgressFn>,
    ) -> Result<(), SlcFetchError> {
        let part = {
            let mut p = dest.as_os_str().to_owned();
            p.push(".part");
            PathBuf::from(p)
        };

        // Check for an existing partial file to resume from.
        let offset: u64 = if part.exists() {
            std::fs::metadata(&part).map_err(SlcFetchError::Io)?.len()
        } else {
            0
        };

        let agent = ureq::builder()
            .timeout(std::time::Duration::from_secs(HTTP_TIMEOUT_S))
            .redirects(0)
            .build();

        let (response, actual_offset) = if offset > 0 {
            let (_, resp) = get_with_auth(&agent, url, token, product_id, Some(offset))?;
            match resp.status() {
                206 => (resp, offset),
                200 => {
                    // Server does not support Range — restart from the beginning.
                    tracing::warn!(
                        "server returned 200 for Range request on {}; \
                         discarding partial file and restarting",
                        url
                    );
                    std::fs::remove_file(&part)?;
                    let (_, resp2) = get_with_auth(&agent, url, token, product_id, None)?;
                    (resp2, 0)
                }
                s => {
                    return Err(SlcFetchError::Http {
                        url: url.to_string(),
                        detail: format!("unexpected status {} on Range request", s),
                    })
                }
            }
        } else {
            let (_, resp) = get_with_auth(&agent, url, token, product_id, None)?;
            (resp, 0)
        };

        // Parse Content-Length to compute total size if available.
        let bytes_total: Option<u64> = response
            .header("Content-Length")
            .and_then(|v| v.parse::<u64>().ok()) // SAFETY-OK: None means header absent or non-numeric; caller handles None progress gracefully
            .map(|cl| actual_offset + cl);

        // Open the part file: append if resuming, create fresh otherwise.
        let mut out = if actual_offset > 0 {
            std::fs::OpenOptions::new()
                .append(true)
                .open(&part)
                .map_err(SlcFetchError::Io)?
        } else {
            std::fs::File::create(&part).map_err(SlcFetchError::Io)?
        };

        let mut reader = response.into_reader();
        let mut buf = vec![0u8; BUF_SIZE];
        let mut bytes_done: u64 = actual_offset;
        loop {
            let n = reader.read(&mut buf).map_err(SlcFetchError::Io)?;
            if n == 0 {
                break;
            }
            out.write_all(&buf[..n])?;
            bytes_done += n as u64;
            if let Some(cb) = progress {
                cb(bytes_done, bytes_total);
            }
        }
        out.flush()?;
        drop(out);

        std::fs::rename(&part, dest)?;
        Ok(())
    }

    /// Number of parallel range-download threads for SLC ZIP files.
    ///
    /// Reads `SARDINE_SLC_CONCURRENCY` from the environment; falls back to 8.
    /// Values below 1 are treated as 1 (sequential delegation to single-stream).
    fn slc_concurrency() -> usize {
        std::env::var("SARDINE_SLC_CONCURRENCY")
            .ok()
            .and_then(|v| v.parse::<usize>().ok())
            .filter(|&n| n >= 1)
            .unwrap_or(8) // SAFETY-OK: env var parse fallback; 8 is a safe default concurrency, not a physical constant
    }

    /// Download a URL to a file using N parallel byte-range threads.
    ///
    /// # Protocol
    ///
    /// 1. Resolves the final S3 pre-signed URL and `Content-Length` by following
    ///    the ASF 3-hop redirect chain via [`get_with_auth`] (probe request).
    /// 2. Falls back to single-stream [`download_to_file`] if:
    ///    - a `.part` file exists (resume in progress), or
    ///    - the server did not return `Content-Length`, or
    ///    - `SARDINE_SLC_CONCURRENCY=1` or the file is tiny (< 1 MB).
    /// 3. Pre-allocates the destination file via `set_len`, then spawns
    ///    `N` threads with `std::thread::scope`.  Each thread fetches one byte
    ///    range directly from the S3 URL (no `Authorization` header needed —
    ///    S3 pre-signed URLs embed credentials in query parameters).
    ///    Writes are issued via `write_at` (POSIX `pwrite`), which is safe for
    ///    concurrent non-overlapping ranges on the same file descriptor.
    /// 4. On success, atomically renames `{dest}.tmp` to `dest`.
    ///    On failure, removes `{dest}.tmp` and returns the first error.
    fn download_to_file_parallel(
        url: &str,
        dest: &Path,
        product_id: &str,
        token: &str,
        progress: Option<&ProgressFn>,
    ) -> Result<(), SlcFetchError> {
        // If a `.part` file exists a single-stream resume is already in
        // progress — delegate and preserve resume semantics.
        let part = {
            let mut p = dest.as_os_str().to_owned();
            p.push(".part");
            PathBuf::from(p)
        };
        if part.exists() {
            tracing::info!(
                "partial file found at {}, resuming single-stream download",
                part.display()
            );
            return download_to_file(url, dest, product_id, token, progress);
        }

        let n_threads = slc_concurrency();

        // Probe the redirect chain to obtain the final S3 URL + Content-Length.
        let agent = ureq::builder()
            .timeout(std::time::Duration::from_secs(HTTP_TIMEOUT_S))
            .redirects(0)
            .build();
        let (final_url, probe_resp) =
            get_with_auth(&agent, url, token, product_id, None)?;

        let total_bytes = probe_resp
            .header("Content-Length")
            .and_then(|v| v.parse::<u64>().ok()); // SAFETY-OK: None means header absent; fallback to single-stream is the explicit next branch
        drop(probe_resp); // close probe connection before spawning threads

        let total_bytes = match total_bytes {
            Some(len) if n_threads > 1 && len >= 1024 * 1024 => len,
            Some(_) => {
                // Single-thread mode or file too small for parallel benefit.
                tracing::debug!("parallel download skipped (single-thread or small file), using single-stream");
                return download_to_file(url, dest, product_id, token, progress);
            }
            None => {
                tracing::info!(
                    "no Content-Length from server for {}, falling back to single-stream download",
                    product_id
                );
                return download_to_file(url, dest, product_id, token, progress);
            }
        };

        tracing::info!(
            "downloading {} ({:.1} MB) with {} threads …",
            product_id,
            total_bytes as f64 / 1e6,
            n_threads,
        );

        // Pre-allocate the destination in a .tmp file (not .part — we write all
        // ranges before making the file visible, so partial state is never
        // observed as a complete file).
        let tmp = {
            let mut p = dest.as_os_str().to_owned();
            p.push(".tmp");
            PathBuf::from(p)
        };
        {
            let f = std::fs::File::create(&tmp)?;
            f.set_len(total_bytes)?;
        }

        // Re-open in write mode and share across threads.
        // `write_at` (POSIX pwrite) is safe for non-overlapping ranges on the
        // same fd: it does not modify the file position cursor.
        let file = Arc::new(
            std::fs::OpenOptions::new()
                .write(true)
                .open(&tmp)
                .map_err(SlcFetchError::Io)?,
        );

        let chunk_size = (total_bytes + n_threads as u64 - 1) / n_threads as u64;
        let ranges: Vec<(u64, u64)> = (0..n_threads as u64)
            .map(|i| {
                let start = i * chunk_size;
                let end = ((i + 1) * chunk_size - 1).min(total_bytes - 1);
                (start, end)
            })
            .collect();

        let bytes_done = Arc::new(AtomicU64::new(0u64));

        let thread_result: Result<(), SlcFetchError> = std::thread::scope(|scope| {
            let mut handles = Vec::with_capacity(n_threads);

            for (range_start, range_end) in &ranges {
                let range_start = *range_start;
                let range_end = *range_end;
                let file = Arc::clone(&file);
                let bytes_done = Arc::clone(&bytes_done);
                let progress_arc = progress.cloned();
                let final_url = final_url.clone();

                handles.push(scope.spawn(move || -> Result<(), SlcFetchError> {
                    // S3 pre-signed URLs encode credentials in query params.
                    // Do NOT send an Authorization header — S3 returns 400 if
                    // an Authorization header accompanies a query-param signature.
                    let thread_agent = ureq::builder()
                        .timeout(std::time::Duration::from_secs(HTTP_TIMEOUT_S))
                        .build();
                    let range_header = format!("bytes={}-{}", range_start, range_end);
                    let resp = thread_agent
                        .get(&final_url)
                        .set("User-Agent", "sardine-scene/0.1")
                        .set("Range", &range_header)
                        .call()
                        .map_err(|e| SlcFetchError::Http {
                            url: final_url.clone(),
                            detail: format!("range {}-{} request failed: {}", range_start, range_end, e),
                        })?;

                    if resp.status() != 206 {
                        return Err(SlcFetchError::Http {
                            url: final_url,
                            detail: format!(
                                "expected HTTP 206 for range {}-{}, got {}",
                                range_start, range_end, resp.status()
                            ),
                        });
                    }

                    let mut reader = resp.into_reader();
                    let mut buf = vec![0u8; BUF_SIZE];
                    let mut offset = range_start;
                    loop {
                        let n = reader.read(&mut buf).map_err(SlcFetchError::Io)?;
                        if n == 0 {
                            break;
                        }
                        file.write_at(&buf[..n], offset).map_err(SlcFetchError::Io)?;
                        offset += n as u64;
                        let done = bytes_done.fetch_add(n as u64, Ordering::Relaxed) + n as u64;
                        if let Some(ref cb) = progress_arc {
                            cb(done, Some(total_bytes));
                        }
                    }
                    Ok(())
                }));
            }

            for handle in handles {
                handle
                    .join()
                    .map_err(|_| {
                        SlcFetchError::Io(io::Error::new(
                            io::ErrorKind::Other,
                            "SLC download worker thread panicked",
                        ))
                    })??;
            }
            Ok(())
        });

        // Drop the Arc<File> before renaming so the fd is closed.
        drop(file);

        if let Err(e) = thread_result {
            if tmp.exists() {
                let _ = std::fs::remove_file(&tmp); // SAFETY-OK: cleanup of failed parallel download; removal failure leaves only a harmless .tmp file
            }
            return Err(e);
        }

        std::fs::rename(&tmp, dest)?;
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

    /// Download and extract a Sentinel-1 SLC product described by a
    /// [`SlcSearchResult`] from the ASF catalogue.
    ///
    /// # Arguments
    /// - `result`: a record from [`search_slc`] containing the download URL and MD5.
    /// - `output_dir`: directory to place the extracted `.SAFE` into.
    /// - `creds`: Earthdata credentials.  A bearer token is acquired on first use
    ///   if [`SlcCredentials::UsernamePassword`] is provided.
    /// - `progress`: optional callback invoked after every 512 KB chunk with
    ///   `(bytes_downloaded, total_bytes_or_none)`.
    /// - `verify_checksum`: if `true`, verifies the MD5 of the downloaded ZIP
    ///   against `result.md5sum`.  Returns [`SlcFetchError::NoChecksum`] if the
    ///   catalogue provides no MD5.
    ///
    /// # Idempotency
    /// If `{output_dir}/{scene_name}.SAFE` already exists, returns immediately.
    pub fn fetch_slc_result(
        result: &SlcSearchResult,
        output_dir: &Path,
        creds: &SlcCredentials,
        progress: Option<ProgressFn>,
        verify_checksum: bool,
    ) -> Result<PathBuf, SlcFetchError> {
        let safe_dir = output_dir.join(format!("{}.SAFE", result.scene_name));

        if safe_dir.exists() {
            tracing::info!(
                "{} already exists, skipping download",
                safe_dir.display()
            );
            return Ok(safe_dir);
        }

        std::fs::create_dir_all(output_dir)?;

        let zip_path = output_dir.join(format!("{}.zip", result.scene_name));

        if zip_path.exists() {
            tracing::info!(
                "ZIP already present at {}, skipping download",
                zip_path.display()
            );
        } else {
            let token = resolve_token(creds)?;
            tracing::info!("downloading {} …", result.url);
            download_to_file_parallel(
                &result.url,
                &zip_path,
                &result.scene_name,
                &token,
                progress.as_ref(),
            )?;
            tracing::info!("download complete → {}", zip_path.display());
        }

        if verify_checksum {
            match &result.md5sum {
                Some(expected) => {
                    tracing::info!("verifying MD5 …");
                    verify_md5(&zip_path, expected, &result.scene_name)?;
                    tracing::info!("MD5 OK");
                }
                None => {
                    return Err(SlcFetchError::NoChecksum {
                        product_id: result.scene_name.clone(),
                    });
                }
            }
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
    /// This is the backward-compatible entry point.  It locates the product in
    /// the ASF catalogue via `granule_list` (no time-window fuzz, exact match),
    /// then delegates to [`fetch_slc_result`].
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
    /// If `{output_dir}/{product_id}.SAFE` already exists, returns immediately
    /// (no network I/O).
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

        // Fast path: SAFE already exists — no network I/O needed.
        if safe_dir.exists() {
            tracing::info!(
                "{} already exists, skipping download",
                safe_dir.display()
            );
            return Ok(safe_dir);
        }

        // Look up the product in the ASF catalogue to obtain the download URL
        // and MD5.  Using granule_list gives an exact match with no time-window
        // fuzz and no risk of returning multiple products.
        let lookup_url = format!(
            "{}?granule_list={}&output=geojson",
            ASF_SEARCH_URL, product_id
        );
        let agent = ureq::builder()
            .timeout(std::time::Duration::from_secs(SEARCH_TIMEOUT_S))
            .build();
        let resp = agent
            .get(&lookup_url)
            .set("User-Agent", "sardine-scene/0.1")
            .call()
            .map_err(|e| SlcFetchError::SearchFailed {
                url: lookup_url.clone(),
                detail: e.to_string(),
            })?;
        let body = resp.into_string().map_err(|e| SlcFetchError::SearchFailed {
            url: lookup_url.clone(),
            detail: e.to_string(),
        })?;
        let results = parse_geojson_results(&body)?;

        if results.is_empty() {
            return Err(SlcFetchError::NotFound {
                product_id: product_id.to_string(),
            });
        }

        std::fs::create_dir_all(output_dir)?;

        fetch_slc_result(
            &results[0],
            output_dir,
            &SlcCredentials::BearerToken(token.to_string()),
            None,
            true,
        )
    }

    // ─────────────────────────────────────────────────────────────────────────
    // Tests
    // ─────────────────────────────────────────────────────────────────────────

    #[cfg(test)]
    mod tests {
        use super::*;
        use chrono::{Datelike, Timelike};

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
            // Fast path triggers before any network I/O.
            let result = fetch_slc(product_id, dir.path(), "fake-token");
            assert!(result.is_ok(), "{:?}", result);
            assert_eq!(result.unwrap(), safe);
        }

        // ── search URL building ───────────────────────────────────────────────

        #[test]
        fn test_search_url_required_params() {
            use chrono::TimeZone;
            let params = SlcSearchParams {
                start: chrono::Utc.with_ymd_and_hms(2019, 1, 1, 0, 0, 0).unwrap(),
                end: chrono::Utc.with_ymd_and_hms(2019, 2, 1, 0, 0, 0).unwrap(),
                ..Default::default()
            };
            let url = build_search_url(&params);
            assert!(url.contains("dataset=SENTINEL-1"), "url={}", url);
            assert!(url.contains("processingLevel=SLC"), "url={}", url);
            assert!(url.contains("output=geojson"), "url={}", url);
            assert!(url.contains("start=2019-01-01T00%3A00%3A00Z") || url.contains("start=2019-01-01T00:00:00Z"), "url={}", url);
        }

        #[test]
        fn test_search_url_bbox_encoded() {
            use chrono::TimeZone;
            let params = SlcSearchParams {
                start: chrono::Utc.with_ymd_and_hms(2019, 1, 1, 0, 0, 0).unwrap(),
                end: chrono::Utc.with_ymd_and_hms(2019, 2, 1, 0, 0, 0).unwrap(),
                bbox: Some([7.0, 47.0, 9.0, 48.0]),
                ..Default::default()
            };
            let url = build_search_url(&params);
            assert!(url.contains("intersectsWith="), "url={}", url);
            assert!(url.contains("POLYGON"), "url={}", url);
            // WKT parentheses must be encoded
            assert!(!url.contains("POLYGON(("), "unencoded parens in url={}", url);
        }

        #[test]
        fn test_search_url_max_results_capped_at_250() {
            use chrono::TimeZone;
            let params = SlcSearchParams {
                start: chrono::Utc.with_ymd_and_hms(2019, 1, 1, 0, 0, 0).unwrap(),
                end: chrono::Utc.with_ymd_and_hms(2019, 2, 1, 0, 0, 0).unwrap(),
                max_results: 9999,
                ..Default::default()
            };
            let url = build_search_url(&params);
            assert!(url.contains("maxResults=250"), "url={}", url);
        }

        // ── GeoJSON response parsing ──────────────────────────────────────────

        #[test]
        fn test_parse_geojson_empty_features() {
            let body = r#"{"type":"FeatureCollection","features":[]}"#;
            let results = parse_geojson_results(body).unwrap();
            assert!(results.is_empty());
        }

        #[test]
        fn test_parse_geojson_one_result() {
            let body = r#"{
                "type": "FeatureCollection",
                "features": [{
                    "type": "Feature",
                    "geometry": null,
                    "properties": {
                        "sceneName": "S1B_IW_SLC__1SDV_20190123T053348_20190123T053415_014617_01B3D4_E833",
                        "url": "https://datapool.asf.alaska.edu/SLC/SB/S1B_IW_SLC__1SDV.zip",
                        "bytes": 4294967296.0,
                        "md5sum": "ABCDEF1234567890abcdef1234567890",
                        "startTime": "2019-01-23T05:33:48.000000",
                        "stopTime": "2019-01-23T05:34:15.000000",
                        "platform": "Sentinel-1B",
                        "beamModeType": "IW",
                        "polarization": "VV+VH",
                        "flightDirection": "DESCENDING",
                        "pathNumber": 15,
                        "frameNumber": 83
                    }
                }]
            }"#;
            let results = parse_geojson_results(body).unwrap();
            assert_eq!(results.len(), 1);
            let r = &results[0];
            assert_eq!(
                r.scene_name,
                "S1B_IW_SLC__1SDV_20190123T053348_20190123T053415_014617_01B3D4_E833"
            );
            assert_eq!(r.bytes, 4_294_967_296);
            // md5sum is lowercased on ingest
            assert_eq!(r.md5sum.as_deref(), Some("abcdef1234567890abcdef1234567890"));
            assert_eq!(r.platform, "Sentinel-1B");
            assert_eq!(r.beam_mode, "IW");
            assert_eq!(r.polarization, "VV+VH");
            assert_eq!(r.orbit_direction, "DESCENDING");
            assert_eq!(r.relative_orbit, Some(15));
            assert_eq!(r.frame_number, Some(83));
        }

        #[test]
        fn test_parse_geojson_null_md5() {
            // Null md5sum is allowed — stored as None.
            let body = r#"{
                "type": "FeatureCollection",
                "features": [{
                    "type": "Feature",
                    "geometry": null,
                    "properties": {
                        "sceneName": "S1A_IW_SLC__1SDV_20201005T170824_20201005T170851_034664_04098A_1E66",
                        "url": "https://datapool.asf.alaska.edu/SLC/SA/test.zip",
                        "bytes": 3000000000.0,
                        "md5sum": null,
                        "startTime": "2020-10-05T17:08:24.000000",
                        "stopTime": "2020-10-05T17:08:51.000000",
                        "platform": "Sentinel-1A",
                        "beamModeType": "IW",
                        "polarization": "VV+VH",
                        "flightDirection": "ASCENDING",
                        "pathNumber": null,
                        "frameNumber": null
                    }
                }]
            }"#;
            let results = parse_geojson_results(body).unwrap();
            assert_eq!(results.len(), 1);
            assert!(results[0].md5sum.is_none());
            assert!(results[0].relative_orbit.is_none());
        }

        #[test]
        fn test_parse_geojson_missing_required_field_errors() {
            // Missing 'url' field → SearchParseError.
            let body = r#"{
                "type": "FeatureCollection",
                "features": [{
                    "type": "Feature",
                    "geometry": null,
                    "properties": {
                        "sceneName": "S1A_IW_SLC__1SDV",
                        "bytes": 1.0,
                        "md5sum": null,
                        "startTime": "2020-10-05T17:08:24.000000",
                        "stopTime": "2020-10-05T17:08:51.000000",
                        "platform": "Sentinel-1A",
                        "beamModeType": "IW",
                        "polarization": "VV",
                        "flightDirection": "ASCENDING"
                    }
                }]
            }"#;
            let err = parse_geojson_results(body).unwrap_err();
            assert!(
                matches!(err, SlcFetchError::SearchParseError { .. }),
                "expected SearchParseError, got {:?}",
                err
            );
        }

        // ── credential resolution ─────────────────────────────────────────────

        #[test]
        fn test_credentials_from_env_token() {
            // Temporarily set EARTHDATA_TOKEN in this test's environment.
            // Use a unique key to avoid colliding with a real token.
            std::env::set_var("EARTHDATA_TOKEN", "test-bearer-xyz");
            let creds = credentials_from_env().unwrap();
            std::env::remove_var("EARTHDATA_TOKEN");
            assert!(matches!(creds, SlcCredentials::BearerToken(t) if t == "test-bearer-xyz"));
        }

        #[test]
        fn test_credentials_from_env_user_pass() {
            std::env::remove_var("EARTHDATA_TOKEN");
            std::env::set_var("EARTHDATA_USERNAME", "testuser");
            std::env::set_var("EARTHDATA_PASSWORD", "testpass");
            let creds = credentials_from_env().unwrap();
            std::env::remove_var("EARTHDATA_USERNAME");
            std::env::remove_var("EARTHDATA_PASSWORD");
            assert!(matches!(
                creds,
                SlcCredentials::UsernamePassword { username, password }
                if username == "testuser" && password == "testpass"
            ));
        }

        #[test]
        fn test_credentials_from_env_neither_set() {
            std::env::remove_var("EARTHDATA_TOKEN");
            std::env::remove_var("EARTHDATA_USERNAME");
            std::env::remove_var("EARTHDATA_PASSWORD");
            let err = credentials_from_env().unwrap_err();
            assert!(
                matches!(err, SlcFetchError::TokenNotSet),
                "expected TokenNotSet, got {:?}",
                err
            );
        }

        // ── MD5 verification ─────────────────────────────────────────────────

        #[test]
        fn test_verify_md5_match() {
            let dir = tempfile::tempdir().unwrap();
            let path = dir.path().join("test.bin");
            let data = b"hello sardine";
            std::fs::write(&path, data).unwrap();
            // Known MD5 of "hello sardine"
            let expected = format!("{:x}", md5::compute(data));
            verify_md5(&path, &expected, "test").unwrap();
        }

        #[test]
        fn test_verify_md5_mismatch() {
            let dir = tempfile::tempdir().unwrap();
            let path = dir.path().join("test.bin");
            std::fs::write(&path, b"hello sardine").unwrap();
            let wrong = "00000000000000000000000000000000";
            let err = verify_md5(&path, wrong, "test-product").unwrap_err();
            assert!(
                matches!(err, SlcFetchError::ChecksumMismatch { ref product_id, .. }
                    if product_id == "test-product"),
                "expected ChecksumMismatch, got {:?}",
                err
            );
        }

        // ── datetime parsing ─────────────────────────────────────────────────

        #[test]
        fn test_parse_asf_datetime_with_fractional() {
            let dt = parse_asf_datetime("2019-01-23T05:33:48.000000").unwrap();
            assert_eq!(dt.year(), 2019);
            assert_eq!(dt.month(), 1);
            assert_eq!(dt.day(), 23);
            assert_eq!(dt.hour(), 5);
            assert_eq!(dt.minute(), 33);
            assert_eq!(dt.second(), 48);
        }

        #[test]
        fn test_parse_asf_datetime_invalid_errors() {
            let err = parse_asf_datetime("not-a-date");
            assert!(err.is_err(), "expected parse error");
        }
    }
}

#[cfg(feature = "slc-fetch")]
pub use inner::{
    acquire_token, credentials_from_env, fetch_slc, fetch_slc_result, search_slc,
    token_from_env, ProgressFn, SlcCredentials, SlcFetchError, SlcSearchParams,
    SlcSearchResult,
};
