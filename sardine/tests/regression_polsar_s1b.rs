//! PolSAR H/A/α regression test: S1B IW SLC 2019-01-23, VV+VH.
//!
//! Verifies the Cloude–Pottier dual-pol decomposition pipeline end-to-end
//! against known-good per-pixel statistics derived from the first run of
//! the code at commit 77b9393 (all audit fixes applied).
//!
//! # Pass criteria (two tiers)
//!
//! **Sanity bounds** — always checked:
//! - H valid fraction ≥ 0.50
//! - H ∈ [0, 1], A ∈ [0, 1], alpha ∈ [10°, 80°]
//!
//! **Tight regression** — checked after `BASELINES_POPULATED = true`:
//! | Output | Metric        | Tolerance |
//! |--------|---------------|-----------|
//! | H      | mean, p50     | ±0.005    |
//! | A      | mean          | ±0.005    |
//! | alpha  | mean, p50     | ±0.2°     |
//!
//! # Populating baseline values
//!
//! ```sh
//! cargo test --release --features geoid-fetch \
//!     --test regression_polsar_s1b -- --ignored --nocapture
//! ```
//!
//! Copy the `BASELINE VALUES CAPTURED` block from the output into
//! the `BASELINE_*` constants below, then set `BASELINES_POPULATED = true`.
//!
//! # Opt-in
//!
//! Same command as above.  The test skips cleanly when fixtures are absent.

#![cfg(feature = "geoid-fetch")]

use std::io::BufReader;
use std::path::{Path, PathBuf};

use sardine::pipeline_options::OutputMode;
use sardine::run::ProcessOptions;

// ─── Fixtures ─────────────────────────────────────────────────────────────────

const SAFE: &str = "/home/datacube/dev/SARdine/data/SLC/\
    S1B_IW_SLC__1SDV_20190123T053348_20190123T053415_014617_01B3D4_E833.SAFE";

const EOF: &str = "/home/datacube/dev/SARdine/data/SLC/\
    S1B_IW_SLC__1SDV_20190123T053348_20190123T053415_014617_01B3D4_E833.SAFE/orbit_cache/\
    S1B_IW_SLC__1SDV_20190123T053348_20190123T053415_014617_01B3D4_E833_POEORB.EOF";

const DEM: &str = "/home/datacube/dev/SARdine/data/dem/srtm1";

fn fixtures_present() -> bool {
    Path::new(SAFE).is_dir()
        && Path::new(EOF).is_file()
        && Path::new(DEM).is_dir()
}

// ─── TIFF reader ──────────────────────────────────────────────────────────────

fn read_f32_tiff(path: &Path) -> Vec<f32> {
    use tiff::decoder::{Decoder, DecodingResult, Limits};

    let file = std::fs::File::open(path)
        .unwrap_or_else(|e| panic!("cannot open {}: {e}", path.display()));
    let mut dec = Decoder::new(BufReader::new(file))
        .unwrap_or_else(|e| panic!("TIFF decode: {e}"))
        .with_limits(Limits::unlimited());
    let img = dec
        .read_image()
        .unwrap_or_else(|e| panic!("read_image {}: {e}", path.display()));
    match img {
        DecodingResult::F32(v) => v,
        other => panic!("expected F32 TIFF, got {other:?}: {}", path.display()),
    }
}

// ─── Statistics ───────────────────────────────────────────────────────────────

#[derive(Debug)]
struct Stats {
    valid_frac: f64,
    mean: f32,
    p25: f32,
    p50: f32,
    p75: f32,
}

fn compute_stats(pixels: &[f32]) -> Stats {
    let mut finite: Vec<f32> = pixels.iter().copied().filter(|v| v.is_finite()).collect();
    let valid_frac = finite.len() as f64 / pixels.len().max(1) as f64;
    if finite.is_empty() {
        return Stats {
            valid_frac: 0.0,
            mean: f32::NAN,
            p25: f32::NAN,
            p50: f32::NAN,
            p75: f32::NAN,
        };
    }
    let mean = finite.iter().copied().map(|v| v as f64).sum::<f64>() / finite.len() as f64;
    finite.sort_unstable_by(|a, b| a.partial_cmp(b).unwrap()); // safe: NaN filtered
    let n = finite.len();
    Stats {
        valid_frac,
        mean: mean as f32,
        p25: finite[n / 4],
        p50: finite[n / 2],
        p75: finite[3 * n / 4],
    }
}

// ─── Baseline constants ───────────────────────────────────────────────────────
//
// Fill these in from the `BASELINE VALUES CAPTURED` block printed by the first
// run, then set `BASELINES_POPULATED = true`.

const BASELINE_H_MEAN: f32 = 0.28111_f32;
const BASELINE_H_P50: f32 = 0.25659_f32;
const BASELINE_A_MEAN: f32 = 0.87000_f32;
const BASELINE_ALPHA_MEAN: f32 = 22.1663_f32;
const BASELINE_ALPHA_P50: f32 = 20.7998_f32;
const BASELINES_POPULATED: bool = true; // set true after filling

// ─── Test ─────────────────────────────────────────────────────────────────────

#[test]
#[ignore]
fn polsar_s1b_20190123_haalpha_statistics() {
    if !fixtures_present() {
        eprintln!(
            "regression_polsar_s1b: skipping — fixtures not present.\n\
             Required:\n  {SAFE}\n  {EOF}\n  {DEM}"
        );
        return;
    }

    let tmp = std::env::temp_dir();
    let output = tmp.join("sardine_polsar_s1b.tif");

    let h_path = tmp.join("sardine_polsar_s1b_H.tif");
    let a_path = tmp.join("sardine_polsar_s1b_A.tif");
    let alpha_path = tmp.join("sardine_polsar_s1b_alpha.tif");

    // Remove any stale outputs from previous runs.
    for p in [&h_path, &a_path, &alpha_path] {
        let _ = std::fs::remove_file(p);
    }

    let mut opts = ProcessOptions::new(
        PathBuf::from(SAFE),
        Some(PathBuf::from(DEM)),
        output,
        "auto".to_owned(),
    );
    opts.orbit = Some(PathBuf::from(EOF));
    opts.polarization = "VV+VH".to_owned();
    opts.mode = OutputMode::Polsar;
    opts.no_provenance = true;
    opts.pixel_spacing_m = 10.0;
    // 3×1 range×azimuth multilook for C2 — matches CLI default
    opts.multilook_range = 3;
    opts.multilook_azimuth = 1;

    eprintln!("regression_polsar_s1b: running PolSAR pipeline …");
    sardine::run::run_process_multi(&opts)
        .unwrap_or_else(|e| panic!("PolSAR pipeline failed: {e:#}"));

    for p in [&h_path, &a_path, &alpha_path] {
        assert!(p.exists(), "output not created: {}", p.display());
    }

    let h_stats = compute_stats(&read_f32_tiff(&h_path));
    let a_stats = compute_stats(&read_f32_tiff(&a_path));
    let alpha_stats = compute_stats(&read_f32_tiff(&alpha_path));

    eprintln!(
        "H:     valid={:.3}  mean={:.5}  p25={:.5}  p50={:.5}  p75={:.5}",
        h_stats.valid_frac, h_stats.mean, h_stats.p25, h_stats.p50, h_stats.p75
    );
    eprintln!(
        "A:     valid={:.3}  mean={:.5}  p25={:.5}  p50={:.5}  p75={:.5}",
        a_stats.valid_frac, a_stats.mean, a_stats.p25, a_stats.p50, a_stats.p75
    );
    eprintln!(
        "alpha: valid={:.3}  mean={:.4}  p25={:.4}  p50={:.4}  p75={:.4}",
        alpha_stats.valid_frac,
        alpha_stats.mean,
        alpha_stats.p25,
        alpha_stats.p50,
        alpha_stats.p75
    );

    // ── Sanity bounds ─────────────────────────────────────────────────────────
    assert!(
        h_stats.valid_frac >= 0.50,
        "H: only {:.1}% valid pixels — pipeline likely failed",
        100.0 * h_stats.valid_frac
    );
    assert!(
        h_stats.mean >= 0.0 && h_stats.mean <= 1.0,
        "H: mean={:.5} outside [0, 1]",
        h_stats.mean
    );
    assert!(
        h_stats.p50 >= 0.0 && h_stats.p50 <= 1.0,
        "H: median={:.5} outside [0, 1]",
        h_stats.p50
    );
    assert!(
        a_stats.mean >= 0.0 && a_stats.mean <= 1.0,
        "A: mean={:.5} outside [0, 1]",
        a_stats.mean
    );
    assert!(
        alpha_stats.mean > 10.0 && alpha_stats.mean < 80.0,
        "alpha: mean={:.2}° unrealistic for central-Europe land scene",
        alpha_stats.mean
    );

    // ── Tight regression checksums ────────────────────────────────────────────
    if BASELINES_POPULATED {
        const TOL_HA: f32 = 0.005;
        const TOL_DEG: f32 = 0.2;

        assert!(
            (h_stats.mean - BASELINE_H_MEAN).abs() <= TOL_HA,
            "H mean regression: got {:.5}, expected {:.5} ±{TOL_HA}",
            h_stats.mean,
            BASELINE_H_MEAN
        );
        assert!(
            (h_stats.p50 - BASELINE_H_P50).abs() <= TOL_HA,
            "H p50 regression: got {:.5}, expected {:.5} ±{TOL_HA}",
            h_stats.p50,
            BASELINE_H_P50
        );
        assert!(
            (a_stats.mean - BASELINE_A_MEAN).abs() <= TOL_HA,
            "A mean regression: got {:.5}, expected {:.5} ±{TOL_HA}",
            a_stats.mean,
            BASELINE_A_MEAN
        );
        assert!(
            (alpha_stats.mean - BASELINE_ALPHA_MEAN).abs() <= TOL_DEG,
            "alpha mean regression: got {:.4}°, expected {:.4}° ±{TOL_DEG}°",
            alpha_stats.mean,
            BASELINE_ALPHA_MEAN
        );
        assert!(
            (alpha_stats.p50 - BASELINE_ALPHA_P50).abs() <= TOL_DEG,
            "alpha p50 regression: got {:.4}°, expected {:.4}° ±{TOL_DEG}°",
            alpha_stats.p50,
            BASELINE_ALPHA_P50
        );
    } else {
        eprintln!(
            "\n\
            ─────────────────────────────────────────────────────────────────\n\
            BASELINE VALUES CAPTURED — paste into BASELINE_* constants above:\n\
            const BASELINE_H_MEAN:     f32 = {:.5}_f32;\n\
            const BASELINE_H_P50:      f32 = {:.5}_f32;\n\
            const BASELINE_A_MEAN:     f32 = {:.5}_f32;\n\
            const BASELINE_ALPHA_MEAN: f32 = {:.4}_f32;\n\
            const BASELINE_ALPHA_P50:  f32 = {:.4}_f32;\n\
            Then set BASELINES_POPULATED = true.\n\
            ─────────────────────────────────────────────────────────────────",
            h_stats.mean,
            h_stats.p50,
            a_stats.mean,
            alpha_stats.mean,
            alpha_stats.p50
        );
    }
}
