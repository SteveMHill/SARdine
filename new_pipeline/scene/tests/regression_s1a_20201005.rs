//! End-to-end pipeline smoke test: S1A IW SLC 2020-10-05, VV polarization.
//!
//! **Scene**: S1A IW SLC 2020-10-05T170824, orbit 034664, datatake 04098A.
//! **Reference**: None on-disk — this is a smoke test only.  A radiometric
//!   comparison against MPC Catalyst Earth RTC can be run separately via
//!   `scripts/compare_mpc_s1a.py`.
//!
//! # Pass criteria
//!
//! | Check                      | Threshold         |
//! |----------------------------|-------------------|
//! | Pipeline completes         | no panic / error  |
//! | Output file exists         | path is present   |
//! | Fraction finite pixels     | ≥ 50%             |
//! | dB range of valid pixels   | [−40, +5] dB      |
//!
//! These thresholds catch gross failures: missing orbit, wrong DEM extent,
//! inverted LUT sign, NaN propagation bugs, etc.  They do not assert
//! radiometric accuracy (use `compare_mpc_s1a.py` for that).
//!
//! # Opt-in
//!
//! ```sh
//! cargo test --release --features geoid-fetch \
//!     --test regression_s1a_20201005 -- --ignored --nocapture
//! ```
//!
//! The test skips cleanly when any input fixture is absent.
//! On this machine the S1A SAFE has not been downloaded; the test will
//! always skip until the data is present.

#![cfg(feature = "geoid-fetch")]

use std::path::{Path, PathBuf};

use sardine_scene::run::ProcessOptions;

// ─── Fixture paths ────────────────────────────────────────────────────────────

const SAFE: &str = "/home/datacube/dev/SARdine/data/SLC/\
    S1A_IW_SLC__1SDV_20201005T170824_20201005T170851_034664_04098A_1E66.SAFE";

/// POEORB file — follows the same orbit_cache convention as the S1B fixture.
const EOF: &str = "/home/datacube/dev/SARdine/data/SLC/\
    S1A_IW_SLC__1SDV_20201005T170824_20201005T170851_034664_04098A_1E66.SAFE/\
    orbit_cache/\
    S1A_IW_SLC__1SDV_20201005T170824_20201005T170851_034664_04098A_1E66_POEORB.EOF";

const DEM: &str = "/home/datacube/dev/SARdine/data/dem/srtm1";

// ─── Helpers ──────────────────────────────────────────────────────────────────

fn fixtures_present() -> bool {
    Path::new(SAFE).is_dir() && Path::new(DEM).is_dir()
}

fn orbit_present() -> bool {
    Path::new(EOF).is_file()
}

/// Read the first N finite pixel values (dB) from a Float32 strip TIFF for
/// range/sanity checking.  Returns the fraction of finite pixels and the
/// (min, max) of those pixels after dB conversion.
///
/// Skips true NaN and ±Inf.  All finite values are expected to be linear
/// power ≥ 0; negatives are treated as invalid.
fn check_output_tiff(path: &Path) -> (f64, f32, f32) {
    use std::io::BufReader;
    use tiff::decoder::{Decoder, DecodingResult, Limits};

    let file = std::fs::File::open(path)
        .unwrap_or_else(|e| panic!("cannot open {}: {e}", path.display()));
    let mut dec = Decoder::new(BufReader::new(file))
        .unwrap_or_else(|e| panic!("TIFF decode init: {e}"))
        .with_limits(Limits::unlimited());

    let (w, h) = dec.dimensions()
        .unwrap_or_else(|e| panic!("TIFF dimensions: {e}"));
    let total = w as usize * h as usize;

    let img = dec.read_image()
        .unwrap_or_else(|e| panic!("read_image: {e}"));
    let pixels = match img {
        DecodingResult::F32(v) => v,
        other => panic!("expected F32 output, got {other:?}"),
    };
    assert_eq!(pixels.len(), total);

    let mut n_finite = 0usize;
    let mut db_min = f32::INFINITY;
    let mut db_max = f32::NEG_INFINITY;

    for &v in &pixels {
        if v.is_finite() && v > 0.0 {
            n_finite += 1;
            let db = 10.0 * v.log10();
            if db < db_min { db_min = db; }
            if db > db_max { db_max = db; }
        }
    }

    let frac = n_finite as f64 / total as f64;
    (frac, db_min, db_max)
}

// ─── Test ─────────────────────────────────────────────────────────────────────

#[test]
#[ignore]
fn s1a_20201005_vv_pipeline_smoke() {
    if !fixtures_present() {
        eprintln!(
            "regression_s1a_20201005: skipping — S1A SAFE or DEM not present.\n\
             Required:\n  {SAFE}\n  {DEM}\n\
             Download the S1A SLC product to run this test."
        );
        return;
    }

    let output = std::env::temp_dir().join("sardine_s1a_20201005_smoke.tiff");

    let mut opts = ProcessOptions::new(
        PathBuf::from(SAFE),
        PathBuf::from(DEM),
        output.clone(),
        "auto".to_owned(), // EGM96 via geoid-fetch
    );

    if orbit_present() {
        opts.orbit = Some(PathBuf::from(EOF));
        eprintln!("regression_s1a_20201005: using POEORB");
    } else {
        // Allow annotation orbit; sufficient for a smoke test.
        std::env::set_var("SARDINE_ALLOW_ANNOTATION_ORBIT", "1");
        eprintln!(
            "regression_s1a_20201005: POEORB not present at {EOF}\n\
             Falling back to annotation orbit (SARDINE_ALLOW_ANNOTATION_ORBIT=1).\n\
             For precise geolocation, download the POEORB file."
        );
    }

    opts.polarization = "VV".to_owned();
    opts.no_provenance = true;
    opts.crs = "EPSG:32632".to_owned();
    opts.pixel_spacing_m = 10.0;
    opts.speckle = "refined-lee".to_owned();
    opts.enl = 1.0;

    eprintln!("regression_s1a_20201005: running full pipeline …");
    sardine_scene::run::run_process(&opts)
        .unwrap_or_else(|e| panic!("run_process failed: {e:#}"));

    assert!(output.exists(), "output TIFF not created: {}", output.display());

    let (frac_finite, db_min, db_max) = check_output_tiff(&output);

    eprintln!(
        "regression_s1a_20201005: frac_finite={:.3}  dB range=[{db_min:.2}, {db_max:.2}]",
        frac_finite
    );

    assert!(
        frac_finite >= 0.50,
        "only {:.1}% of pixels are finite — pipeline may have failed silently",
        100.0 * frac_finite
    );

    assert!(
        db_min >= -40.0,
        "minimum dB = {db_min:.2} is below −40 dB — likely a calibration error or NaN propagation"
    );

    assert!(
        db_max <= 5.0,
        "maximum dB = {db_max:.2} exceeds +5 dB — likely a calibration error or invalid pixel"
    );
}
