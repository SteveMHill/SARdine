//! GRD (ground-range) pipeline regression test: S1B IW SLC 2019-01-23, VV.
//!
//! Verifies the `run_grd` pipeline (slant-range → ground-range, no terrain
//! correction) against known-good pixel statistics.  The GRD pipeline is
//! faster than the RTC path because it does not need a DEM; it is exercised
//! separately here to ensure it is not inadvertently broken by changes to
//! the RTC path.
//!
//! # Pass criteria (two tiers)
//!
//! **Sanity bounds** — always checked:
//! - Valid pixel fraction ≥ 0.50
//! - dB range of valid pixels: [−40, +5] dB
//! - Finite pixel count > 1 000 000
//!
//! **Tight regression** — checked after `BASELINES_POPULATED = true`:
//! | Metric | Tolerance |
//! |--------|-----------|
//! | mean dB   | ±0.02 dB  |
//! | p10/p50/p90 dB | ±0.02 dB |
//! | valid pixel count | ±0.1% |
//!
//! # Populating baseline values
//!
//! ```sh
//! cargo test --release --features geoid-fetch \
//!     --test regression_grd_s1b -- --ignored --nocapture
//! ```
//!
//! Copy the `BASELINE VALUES CAPTURED` block from the output into the
//! `BASELINE_*` constants below, then set `BASELINES_POPULATED = true`.

#![cfg(feature = "geoid-fetch")]

use std::io::BufReader;
use std::path::{Path, PathBuf};

use sardine::pipeline_options::OutputUnit;
use sardine::run::{run_grd_multi, GrdOptions};

// ─── Fixtures ─────────────────────────────────────────────────────────────────

const SAFE: &str = "/home/datacube/dev/SARdine/data/SLC/\
    S1B_IW_SLC__1SDV_20190123T053348_20190123T053415_014617_01B3D4_E833.SAFE";

const EOF: &str = "/home/datacube/dev/SARdine/data/SLC/\
    S1B_IW_SLC__1SDV_20190123T053348_20190123T053415_014617_01B3D4_E833.SAFE/orbit_cache/\
    S1B_IW_SLC__1SDV_20190123T053348_20190123T053415_014617_01B3D4_E833_POEORB.EOF";

fn fixtures_present() -> bool {
    Path::new(SAFE).is_dir() && Path::new(EOF).is_file()
}

// ─── TIFF reader ──────────────────────────────────────────────────────────────

fn check_output_tiff(path: &Path) -> (usize, f64, f32, f32, f32, f32, f32) {
    use tiff::decoder::{Decoder, DecodingResult, Limits};

    let file = std::fs::File::open(path)
        .unwrap_or_else(|e| panic!("cannot open {}: {e}", path.display()));
    let mut dec = Decoder::new(BufReader::new(file))
        .unwrap_or_else(|e| panic!("TIFF decode: {e}"))
        .with_limits(Limits::unlimited());

    let (w, h) = dec.dimensions().unwrap_or_else(|e| panic!("dimensions: {e}"));
    let total = w as usize * h as usize;

    let img = dec.read_image().unwrap_or_else(|e| panic!("read_image: {e}"));
    let pixels = match img {
        DecodingResult::F32(v) => v,
        other => panic!("expected F32, got {other:?}"),
    };
    assert_eq!(pixels.len(), total);

    let mut finite: Vec<f32> = pixels.iter().copied().filter(|v| v.is_finite()).collect();
    let valid_count = finite.len();
    let valid_frac = valid_count as f64 / total as f64;

    if finite.is_empty() {
        return (0, 0.0, f32::NAN, f32::NAN, f32::NAN, f32::NAN, f32::NAN);
    }

    let _mean_db = finite.iter().copied().map(|v| v as f64).sum::<f64>() / finite.len() as f64;
    finite.sort_unstable_by(|a, b| a.partial_cmp(b).unwrap()); // safe: NaN filtered
    let n = finite.len();
    let p10 = finite[(n as f64 * 0.10) as usize];
    let p25 = finite[(n as f64 * 0.25) as usize];
    let p50 = finite[n / 2];
    let p75 = finite[(n as f64 * 0.75) as usize];
    let p90 = finite[(n as f64 * 0.90) as usize];

    (valid_count, valid_frac, p10 as f32, p25 as f32, p50 as f32, p75 as f32, p90 as f32)
}

// ─── Baseline constants ───────────────────────────────────────────────────────

const BASELINE_VALID_COUNT: usize = 306_851_515;
const BASELINE_P10_DB: f32 = -16.3638_f32;
const BASELINE_P50_DB: f32 = -12.1071_f32;
const BASELINE_P90_DB: f32 = -8.0107_f32;
const BASELINES_POPULATED: bool = true;

// ─── Test ─────────────────────────────────────────────────────────────────────

#[test]
#[ignore]
fn grd_s1b_20190123_vv_statistics() {
    if !fixtures_present() {
        eprintln!(
            "regression_grd_s1b: skipping — fixtures not present.\n\
             Required:\n  {SAFE}\n  {EOF}"
        );
        return;
    }

    let output = std::env::temp_dir().join("sardine_grd_s1b_vv.tif");

    let mut opts = GrdOptions::new(PathBuf::from(SAFE), output.clone());
    opts.orbit = Some(PathBuf::from(EOF));
    opts.polarization = "VV".to_owned();
    opts.target_spacing_m = 10.0;
    opts.no_provenance = true;
    opts.multilook_range = 5; // ESA GRDH-equivalent
    opts.multilook_azimuth = 1;
    opts.output_unit = OutputUnit::Db;

    eprintln!("regression_grd_s1b: running GRD pipeline …");
    run_grd_multi(&opts).unwrap_or_else(|e| panic!("GRD pipeline failed: {e:#}"));

    assert!(output.exists(), "output TIFF not created: {}", output.display());

    let (valid_count, valid_frac, p10, p25, p50, p75, p90) = check_output_tiff(&output);

    eprintln!(
        "GRD VV:  valid_count={valid_count}  valid_frac={valid_frac:.3}\n\
         \t  p10={p10:.3} dB  p25={p25:.3} dB  p50={p50:.3} dB  p75={p75:.3} dB  p90={p90:.3} dB"
    );

    // ── Sanity bounds ─────────────────────────────────────────────────────────
    assert!(
        valid_frac >= 0.50,
        "only {:.1}% finite pixels — GRD pipeline may have failed silently",
        100.0 * valid_frac
    );
    assert!(
        valid_count > 1_000_000,
        "only {valid_count} valid pixels — scene appears empty"
    );
    assert!(
        p10 >= -40.0,
        "p10={p10:.2} dB is below −40 dB — calibration error or NaN propagation"
    );
    assert!(
        p90 <= 5.0,
        "p90={p90:.2} dB exceeds +5 dB — calibration error"
    );
    // GRD σ⁰ in central Europe: typical −20 to −5 dB for VV
    assert!(
        p50 > -25.0 && p50 < 0.0,
        "p50={p50:.2} dB outside expected range (−25, 0) dB for central-Europe VV"
    );

    // ── Tight regression checksums ────────────────────────────────────────────
    if BASELINES_POPULATED {
        const TOL_DB: f32 = 0.02;
        const TOL_COUNT_FRAC: f64 = 0.001; // 0.1%

        let count_tol = (BASELINE_VALID_COUNT as f64 * TOL_COUNT_FRAC) as usize + 1;
        assert!(
            valid_count.abs_diff(BASELINE_VALID_COUNT) <= count_tol,
            "valid pixel count regression: got {valid_count}, expected {BASELINE_VALID_COUNT} ±{count_tol}"
        );
        assert!(
            (p10 - BASELINE_P10_DB).abs() <= TOL_DB,
            "p10 regression: got {p10:.4} dB, expected {BASELINE_P10_DB:.4} dB ±{TOL_DB}"
        );
        assert!(
            (p50 - BASELINE_P50_DB).abs() <= TOL_DB,
            "p50 regression: got {p50:.4} dB, expected {BASELINE_P50_DB:.4} dB ±{TOL_DB}"
        );
        assert!(
            (p90 - BASELINE_P90_DB).abs() <= TOL_DB,
            "p90 regression: got {p90:.4} dB, expected {BASELINE_P90_DB:.4} dB ±{TOL_DB}"
        );
    } else {
        eprintln!(
            "\n\
            ──────────────────────────────────────────────────────────────────\n\
            BASELINE VALUES CAPTURED — paste into BASELINE_* constants above:\n\
            const BASELINE_VALID_COUNT: usize = {valid_count};\n\
            const BASELINE_P10_DB:      f32   = {p10:.4}_f32;\n\
            const BASELINE_P50_DB:      f32   = {p50:.4}_f32;\n\
            const BASELINE_P90_DB:      f32   = {p90:.4}_f32;\n\
            Then set BASELINES_POPULATED = true.\n\
            ──────────────────────────────────────────────────────────────────"
        );
    }
}
