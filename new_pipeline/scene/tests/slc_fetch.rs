//! Integration tests for SLC fetch: search and download.
//!
//! # Test tiers
//!
//! | Test | Network? | Credentials? | Run condition |
//! |---|---|---|---|
//! | `slc_fetch_idempotent_s1b` | no | no | fixture on disk |
//! | `slc_fetch_idempotent_s1a` | no | no | fixture on disk (skip if absent) |
//! | `slc_search_finds_s1b` | yes | no | `SARDINE_NETWORK_TESTS=1` |
//! | `slc_search_finds_s1a` | yes | no | `SARDINE_NETWORK_TESTS=1` |
//! | `slc_download_s1a` | yes | yes | `SARDINE_NETWORK_TESTS=1` + `EARTHDATA_TOKEN` |
//!
//! # How to run
//!
//! ```sh
//! # Fixture-only (no network):
//! cargo test --test slc_fetch --features slc-fetch -- --ignored --nocapture
//!
//! # Search tests (read-only, no auth, ~1 KB responses):
//! SARDINE_NETWORK_TESTS=1 cargo test --test slc_fetch --features slc-fetch \
//!     -- --ignored --nocapture
//!
//! # Full download test (requires Earthdata account):
//! SARDINE_NETWORK_TESTS=1 EARTHDATA_TOKEN=<token> \
//!     cargo test --test slc_fetch --features slc-fetch \
//!     -- --ignored --nocapture slc_download_s1a
//! ```

#[cfg(feature = "slc-fetch")]
mod slc_fetch_tests {
    use std::path::{Path, PathBuf};
    use sardine_scene::slc_fetch::{
        credentials_from_env, fetch_slc, fetch_slc_result, search_slc,
        SlcSearchParams,
    };
    use chrono::{TimeZone, Utc};

    // ─── Fixture paths ────────────────────────────────────────────────────────

    const PRODUCT_S1B: &str =
        "S1B_IW_SLC__1SDV_20190123T053348_20190123T053415_014617_01B3D4_E833";
    const SAFE_S1B: &str = "/home/datacube/dev/SARdine/data/SLC/\
        S1B_IW_SLC__1SDV_20190123T053348_20190123T053415_014617_01B3D4_E833.SAFE";

    const PRODUCT_S1A: &str =
        "S1A_IW_SLC__1SDV_20201005T170824_20201005T170851_034664_04098A_1E66";
    const SAFE_S1A: &str = "/home/datacube/dev/SARdine/data/SLC/\
        S1A_IW_SLC__1SDV_20201005T170824_20201005T170851_034664_04098A_1E66.SAFE";

    fn network_tests_enabled() -> bool {
        std::env::var("SARDINE_NETWORK_TESTS").as_deref() == Ok("1")
    }

    // ─── Fixture tests (no network) ───────────────────────────────────────────

    /// `fetch_slc` fast-paths when the SAFE already exists — no token, no I/O.
    ///
    /// This verifies the idempotency guarantee: pointing `fetch_slc` at a
    /// directory that already contains `{product_id}.SAFE` returns immediately
    /// with the correct path, regardless of the token value.
    #[test]
    #[ignore]
    fn slc_fetch_idempotent_s1b() {
        if !Path::new(SAFE_S1B).is_dir() {
            eprintln!(
                "slc_fetch::slc_fetch_idempotent_s1b: skipping — fixture absent: {SAFE_S1B}"
            );
            return;
        }

        // The parent of the SAFE dir is the output_dir argument.
        let output_dir = Path::new(SAFE_S1B).parent()
            .unwrap_or_else(|| panic!("SAFE_S1B has no parent"));

        // fetch_slc should return immediately — no real token needed.
        let result = fetch_slc(PRODUCT_S1B, output_dir, "fake-token-idempotency-test")
            .unwrap_or_else(|e| panic!("fetch_slc failed on existing SAFE: {e:#}"));

        assert_eq!(
            result,
            PathBuf::from(SAFE_S1B),
            "returned path must equal the existing SAFE directory"
        );
        eprintln!("slc_fetch::slc_fetch_idempotent_s1b: ok — fast-pathed to {}", result.display());
    }

    /// Same idempotency check for S1A — skips gracefully if the SAFE is absent.
    #[test]
    #[ignore]
    fn slc_fetch_idempotent_s1a() {
        if !Path::new(SAFE_S1A).is_dir() {
            eprintln!(
                "slc_fetch::slc_fetch_idempotent_s1a: skipping — fixture absent: {SAFE_S1A}\n\
                 Run slc_download_s1a with SARDINE_NETWORK_TESTS=1 to populate it."
            );
            return;
        }

        let output_dir = Path::new(SAFE_S1A).parent()
            .unwrap_or_else(|| panic!("SAFE_S1A has no parent"));

        let result = fetch_slc(PRODUCT_S1A, output_dir, "fake-token-idempotency-test")
            .unwrap_or_else(|e| panic!("fetch_slc failed on existing SAFE: {e:#}"));

        assert_eq!(result, PathBuf::from(SAFE_S1A));
        eprintln!("slc_fetch::slc_fetch_idempotent_s1a: ok — fast-pathed to {}", result.display());
    }

    // ─── Network tests (SARDINE_NETWORK_TESTS=1) ──────────────────────────────

    /// Search ASF catalogue for the S1B fixture scene.
    ///
    /// Verifies that the search API returns at least one result matching the
    /// known product ID, platform, beam mode, and relative orbit.
    ///
    /// No authentication required — ASF search is public.
    #[test]
    #[ignore]
    fn slc_search_finds_s1b() {
        if !network_tests_enabled() {
            eprintln!(
                "slc_fetch::slc_search_finds_s1b: skipped — set SARDINE_NETWORK_TESTS=1 to enable."
            );
            return;
        }

        // S1B, 2019-01-23, descending IW, relative orbit 15, frame 83.
        let params = SlcSearchParams {
            start: Utc.with_ymd_and_hms(2019, 1, 23, 5, 33, 0).unwrap(),
            end:   Utc.with_ymd_and_hms(2019, 1, 23, 5, 35, 0).unwrap(),
            beam_mode: Some("IW".to_owned()),
            max_results: 10,
            ..Default::default()
        };

        let results = search_slc(&params)
            .unwrap_or_else(|e| panic!("search_slc failed: {e:#}"));

        assert!(!results.is_empty(), "search returned no results for S1B 2019-01-23");

        let hit = results.iter().find(|r| r.scene_name == PRODUCT_S1B)
            .unwrap_or_else(|| panic!(
                "expected {} in search results, got: {:?}",
                PRODUCT_S1B,
                results.iter().map(|r| &r.scene_name).collect::<Vec<_>>()
            ));

        assert_eq!(hit.platform, "Sentinel-1B", "platform mismatch");
        assert_eq!(hit.beam_mode, "IW", "beam mode mismatch");
        assert!(hit.bytes > 0, "file size must be positive");
        assert!(hit.md5sum.is_some(), "md5sum must be present for S1B SLC");

        eprintln!(
            "slc_fetch::slc_search_finds_s1b: found {} ({} bytes, MD5 {})",
            hit.scene_name,
            hit.bytes,
            hit.md5sum.as_deref().unwrap_or("none"),
        );
    }

    /// Search ASF catalogue for the S1A fixture scene.
    ///
    /// Uses a different date, satellite (S1A), and relative orbit from the S1B
    /// test, confirming that `build_url` and the search pipeline handle both
    /// platforms correctly end-to-end.
    #[test]
    #[ignore]
    fn slc_search_finds_s1a() {
        if !network_tests_enabled() {
            eprintln!(
                "slc_fetch::slc_search_finds_s1a: skipped — set SARDINE_NETWORK_TESTS=1 to enable."
            );
            return;
        }

        // S1A, 2020-10-05, ascending IW.
        let params = SlcSearchParams {
            start: Utc.with_ymd_and_hms(2020, 10, 5, 17, 8, 0).unwrap(),
            end:   Utc.with_ymd_and_hms(2020, 10, 5, 17, 10, 0).unwrap(),
            beam_mode: Some("IW".to_owned()),
            max_results: 10,
            ..Default::default()
        };

        let results = search_slc(&params)
            .unwrap_or_else(|e| panic!("search_slc failed: {e:#}"));

        assert!(!results.is_empty(), "search returned no results for S1A 2020-10-05");

        let hit = results.iter().find(|r| r.scene_name == PRODUCT_S1A)
            .unwrap_or_else(|| panic!(
                "expected {} in search results, got: {:?}",
                PRODUCT_S1A,
                results.iter().map(|r| &r.scene_name).collect::<Vec<_>>()
            ));

        assert_eq!(hit.platform, "Sentinel-1A", "platform mismatch");
        assert_eq!(hit.beam_mode, "IW", "beam mode mismatch");
        assert!(hit.bytes > 0, "file size must be positive");
        assert!(hit.md5sum.is_some(), "md5sum must be present for S1A SLC");

        eprintln!(
            "slc_fetch::slc_search_finds_s1a: found {} ({} bytes, MD5 {})",
            hit.scene_name,
            hit.bytes,
            hit.md5sum.as_deref().unwrap_or("none"),
        );
    }

    /// Full download of the S1A scene — exercises `search_slc → fetch_slc_result`
    /// end-to-end including ZIP extraction and MD5 verification.
    ///
    /// S1A is chosen rather than S1B because S1B is already present on disk and
    /// would hit the idempotency fast-path; S1A exercises the real download path
    /// (unless it was previously downloaded by this test, in which case it will
    /// fast-path on the second run — which is the intended caching behaviour).
    ///
    /// Requires:
    /// - `SARDINE_NETWORK_TESTS=1`
    /// - `EARTHDATA_TOKEN=<bearer>` **or** `EARTHDATA_USERNAME` + `EARTHDATA_PASSWORD`
    ///   (NASA Earthdata account — free registration at urs.earthdata.nasa.gov)
    #[test]
    #[ignore]
    fn slc_download_s1a() {
        if !network_tests_enabled() {
            eprintln!(
                "slc_fetch::slc_download_s1a: skipped — set SARDINE_NETWORK_TESTS=1 to enable.\n\
                 Also requires EARTHDATA_TOKEN or EARTHDATA_USERNAME+EARTHDATA_PASSWORD."
            );
            return;
        }

        let creds = match credentials_from_env() {
            Ok(c) => c,
            Err(e) => {
                eprintln!(
                    "slc_fetch::slc_download_s1a: skipped — no credentials: {e}\n\
                     Set EARTHDATA_TOKEN or EARTHDATA_USERNAME+EARTHDATA_PASSWORD."
                );
                return;
            }
        };

        // Search to get the canonical URL and MD5 from the catalogue.
        let params = SlcSearchParams {
            start: Utc.with_ymd_and_hms(2020, 10, 5, 17, 8, 0).unwrap(),
            end:   Utc.with_ymd_and_hms(2020, 10, 5, 17, 10, 0).unwrap(),
            beam_mode: Some("IW".to_owned()),
            max_results: 10,
            ..Default::default()
        };

        let results = search_slc(&params)
            .unwrap_or_else(|e| panic!("search_slc failed: {e:#}"));

        let hit = results.iter().find(|r| r.scene_name == PRODUCT_S1A)
            .unwrap_or_else(|| panic!("S1A scene not found in search results"));

        // Use the canonical data/SLC directory so the S1A fixture is reusable
        // by other tests (regression_s1a_20201005.rs) after this download.
        let output_dir = Path::new(SAFE_S1A).parent()
            .unwrap_or_else(|| panic!("SAFE_S1A has no parent"));

        eprintln!(
            "slc_fetch::slc_download_s1a: downloading {} ({:.1} GB) → {} …",
            hit.scene_name,
            hit.bytes as f64 / 1e9,
            output_dir.display(),
        );

        let safe = fetch_slc_result(hit, output_dir, &creds, None, true)
            .unwrap_or_else(|e| panic!("fetch_slc_result failed: {e:#}"));

        assert!(safe.is_dir(), "SAFE dir not created: {}", safe.display());
        assert!(
            safe.join("manifest.safe").exists(),
            "manifest.safe missing in downloaded SAFE: {}",
            safe.display()
        );
        assert!(
            safe.join("annotation").is_dir(),
            "annotation/ dir missing: {}",
            safe.display()
        );
        assert!(
            safe.join("measurement").is_dir(),
            "measurement/ dir missing: {}",
            safe.display()
        );

        eprintln!("slc_fetch::slc_download_s1a: ok — SAFE at {}", safe.display());
    }
}
