//! `sardine-validate` — end-to-end radiometric regression harness.
//!
//! Reads a scene config (TOML), runs `sardine process`, runs the existing
//! `scripts/multilook_compare.py` script as a subprocess, parses the
//! median bias / std / valid-pixel count from its stdout, and asserts
//! that those metrics are within configured thresholds.  Exits 0 on
//! pass, non-zero on any threshold violation or pipeline error.
//!
//! # Why a Rust orchestrator and not a pure-Rust validator?
//!
//! The radiometric comparison requires GeoTIFF reprojection from the
//! ASF reference's UTM grid to whatever grid the sardine output uses.
//! The proven `multilook_compare.py` script already does this via
//! `rasterio` (GDAL).  Re-implementing reprojection in Rust would be
//! multi-week work for no scientific benefit — the comparison logic
//! is well-tested.  The orchestrator's job is purely:
//!
//!   1. Drive sardine + the comparison subprocess.
//!   2. Parse a small, stable subset of the comparison output.
//!   3. Assert thresholds hard.
//!
//! No silent fallbacks, no partial results: any subprocess failure or
//! parse failure is a hard error.

use std::collections::BTreeMap;
use std::path::{Path, PathBuf};
use std::process::Command;

use clap::Parser;
use serde::Deserialize;

#[derive(Debug, thiserror::Error)]
enum ValidationError {
    #[error("config file {path:?}: {source}")]
    Config { path: PathBuf, #[source] source: anyhow::Error },

    #[error("sardine subprocess failed (exit {exit_code:?}): {stderr_tail}")]
    SardineFailed { exit_code: Option<i32>, stderr_tail: String },

    #[error("comparison subprocess failed (exit {exit_code:?}): {stderr_tail}")]
    ComparisonFailed { exit_code: Option<i32>, stderr_tail: String },

    #[error("could not parse {field} from comparison output (looked for line starting {prefix:?})")]
    ParseFailed { field: &'static str, prefix: &'static str },

    #[error(
        "{scene}: {metric} = {actual:+.4} {unit}, threshold = ±{threshold:.4} {unit} (FAIL)"
    )]
    ThresholdViolation {
        scene: String,
        metric: &'static str,
        actual: f64,
        threshold: f64,
        unit: &'static str,
    },

    #[error("io error: {0}")]
    Io(#[from] std::io::Error),
}

// ── CLI ────────────────────────────────────────────────────────────

/// SARdine end-to-end validation harness.
///
/// Runs the configured scenes through `sardine process` and asserts that
/// the multi-looked radiometric output is within the configured
/// tolerances of the ASF RTC10 reference.
#[derive(Debug, Parser)]
#[command(name = "sardine-validate", version, about, long_about = None)]
struct Cli {
    /// Directory containing scene config TOMLs.  Each `*.toml` file is
    /// loaded, parsed and run.  No silent skipping: any unparseable file
    /// is a hard failure.
    #[arg(long, env = "SARDINE_VALIDATION_SCENES_DIR")]
    scenes_dir: PathBuf,

    /// Path to the `sardine` release binary.  Required — we never
    /// guess based on PATH because the validator must pin the exact
    /// build under test.
    #[arg(long, env = "SARDINE_BIN")]
    sardine_bin: PathBuf,

    /// Path to a Python interpreter that has `numpy` and `rasterio`
    /// importable.  Defaults to `python3`.
    #[arg(long, env = "SARDINE_VALIDATION_PYTHON", default_value = "python3")]
    python: String,

    /// Path to the `multilook_compare.py` script.  Defaults to the
    /// repo-relative path used by the bundled scripts.
    #[arg(
        long,
        env = "SARDINE_MULTILOOK_COMPARE",
        default_value = "scripts/multilook_compare.py"
    )]
    multilook_compare: PathBuf,

    /// Directory to write sardine output rasters into.  One subdirectory
    /// per scene name.  Defaults to a system temp dir.
    #[arg(long, env = "SARDINE_VALIDATION_OUTDIR")]
    outdir: Option<PathBuf>,

    /// Run only the scene with this name (matches the TOML `name` field).
    /// If omitted, runs all scenes in `scenes_dir`.
    #[arg(long)]
    only: Option<String>,
}

// ── Scene config schema ────────────────────────────────────────────

#[derive(Debug, Deserialize)]
struct SceneConfig {
    /// Human-readable scene name.  Must match the TOML file stem.
    name: String,
    /// Path to the `.SAFE` directory.  Absolute or relative to CWD.
    safe: PathBuf,
    /// Path to the POEORB / RESORB orbit file.
    orbit: PathBuf,
    /// Path to the DEM directory or single GeoTIFF.
    dem: PathBuf,
    /// Path to the ASF RTC10 VV (or matching polarisation) reference TIFF.
    asf_reference: PathBuf,
    /// Polarisation to process (`VV`, `VH`, or `dual`).
    polarization: String,
    /// `sardine --crs` argument: `wgs84`, `auto`, or `EPSG:NNNNN`.
    crs: String,
    /// Pixel spacing.  Exactly one of `pixel_spacing_m` /
    /// `pixel_spacing_deg` must be set.
    pixel_spacing_m: Option<f64>,
    pixel_spacing_deg: Option<f64>,
    /// `sardine --geoid` argument: `zero`, `auto`, or a path.
    geoid: String,
    /// Pass-through extra `sardine process` arguments.  Append-only.
    #[serde(default)]
    extra_sardine_args: Vec<String>,
    /// Tolerance bounds for radiometric agreement.
    thresholds: Thresholds,
}

#[derive(Debug, Deserialize)]
struct Thresholds {
    /// Maximum absolute median bias (sardine − ASF), in dB.
    median_bias_db_max: f64,
    /// Maximum multi-looked std deviation (10×10 multilook, ~110 m), in dB.
    /// This residual is dominated by inter-processor spatial registration
    /// differences (SRTM1 vs Copernicus GLO-30) and is not a calibration
    /// metric, but a sudden large jump indicates a regression.
    std_dev_db_max: f64,
    /// Minimum number of jointly-valid multilook samples required.  A
    /// drop below this threshold means too much of the scene was masked
    /// out; sets a floor on the statistical power of the comparison.
    min_joint_samples: u64,
}

// ── Comparison output parsing ──────────────────────────────────────

/// The handful of fields we actually consume from `multilook_compare.py`'s
/// stdout.  Adding new fields is additive; removing one is a breaking
/// change to the validation harness.  See `parse_comparison_output`.
#[derive(Debug)]
struct ComparisonStats {
    median_bias_db: f64,
    std_dev_db: f64,
    n_samples: u64,
}

fn parse_comparison_output(stdout: &str) -> Result<ComparisonStats, ValidationError> {
    fn find_after_prefix<'a>(stdout: &'a str, prefix: &str) -> Option<&'a str> {
        stdout
            .lines()
            .map(str::trim_start)
            .find_map(|line| line.strip_prefix(prefix))
            .map(str::trim)
    }

    fn parse_first_signed_float(s: &str) -> Option<f64> {
        // Examples we must handle:
        //   "+0.0283 dB"
        //   "-10.326 dB"
        //   "2.9572 dB  (single-look was 5.37 dB)"
        let token = s.split_whitespace().next()?;
        token.trim_end_matches(',').parse::<f64>().ok()
    }

    let median = find_after_prefix(stdout, "Median bias:")
        .and_then(parse_first_signed_float)
        .ok_or(ValidationError::ParseFailed {
            field: "median_bias_db",
            prefix: "Median bias:",
        })?;

    let std_dev = find_after_prefix(stdout, "Std dev:")
        .and_then(parse_first_signed_float)
        .ok_or(ValidationError::ParseFailed {
            field: "std_dev_db",
            prefix: "Std dev:",
        })?;

    // "Multi-looked (10×10, ~110 m) comparison: n=4,098,547"
    let line = stdout
        .lines()
        .map(str::trim_start)
        .find(|l| l.starts_with("Multi-looked"))
        .ok_or(ValidationError::ParseFailed {
            field: "n_samples",
            prefix: "Multi-looked",
        })?;
    let n_token = line
        .rsplit("n=")
        .next()
        .ok_or(ValidationError::ParseFailed {
            field: "n_samples",
            prefix: "n=",
        })?;
    let n_samples: u64 = n_token
        .trim()
        .replace(',', "")
        .parse()
        .map_err(|_| ValidationError::ParseFailed {
            field: "n_samples",
            prefix: "n=",
        })?;

    Ok(ComparisonStats {
        median_bias_db: median,
        std_dev_db: std_dev,
        n_samples,
    })
}

// ── Scene loading ──────────────────────────────────────────────────

fn load_scenes(dir: &Path) -> Result<BTreeMap<String, (PathBuf, SceneConfig)>, ValidationError> {
    let mut scenes = BTreeMap::new();
    for entry in std::fs::read_dir(dir)? {
        let entry = entry?;
        let path = entry.path();
        if path.extension().and_then(|s| s.to_str()) != Some("toml") {
            continue;
        }
        let text = std::fs::read_to_string(&path)?;
        let cfg: SceneConfig =
            toml::from_str(&text).map_err(|e| ValidationError::Config {
                path: path.clone(),
                source: anyhow::anyhow!(e),
            })?;
        let stem = path
            .file_stem()
            .and_then(|s| s.to_str())
            .ok_or_else(|| ValidationError::Config {
                path: path.clone(),
                source: anyhow::anyhow!("file stem is not valid UTF-8"),
            })?;
        if cfg.name != stem {
            return Err(ValidationError::Config {
                path: path.clone(),
                source: anyhow::anyhow!(
                    "scene name {:?} does not match file stem {:?}",
                    cfg.name,
                    stem
                ),
            });
        }
        if cfg.pixel_spacing_m.is_some() == cfg.pixel_spacing_deg.is_some() {
            return Err(ValidationError::Config {
                path: path.clone(),
                source: anyhow::anyhow!(
                    "exactly one of pixel_spacing_m / pixel_spacing_deg must be set"
                ),
            });
        }
        scenes.insert(cfg.name.clone(), (path, cfg));
    }
    Ok(scenes)
}

// ── Subprocess drivers ─────────────────────────────────────────────

fn run_sardine(
    cli: &Cli,
    scene: &SceneConfig,
    output_tif: &Path,
) -> Result<(), ValidationError> {
    let mut cmd = Command::new(&cli.sardine_bin);
    cmd.arg("process")
        .arg("--safe")
        .arg(&scene.safe)
        .arg("--orbit")
        .arg(&scene.orbit)
        .arg("--dem")
        .arg(&scene.dem)
        .arg("--output")
        .arg(output_tif)
        .arg("--crs")
        .arg(&scene.crs)
        .arg("--polarization")
        .arg(&scene.polarization)
        .arg("--geoid")
        .arg(&scene.geoid);
    if let Some(m) = scene.pixel_spacing_m {
        cmd.arg("--pixel-spacing-m").arg(format!("{m}"));
    }
    if let Some(d) = scene.pixel_spacing_deg {
        cmd.arg("--pixel-spacing-deg").arg(format!("{d}"));
    }
    for extra in &scene.extra_sardine_args {
        cmd.arg(extra);
    }

    eprintln!("validate: running {cmd:?}");
    let output = cmd.output()?;
    if !output.status.success() {
        let stderr = String::from_utf8_lossy(&output.stderr);
        let tail: String = stderr.lines().rev().take(40).collect::<Vec<_>>().join("\n");
        return Err(ValidationError::SardineFailed {
            exit_code: output.status.code(),
            stderr_tail: tail,
        });
    }
    Ok(())
}

fn run_comparison(
    cli: &Cli,
    scene: &SceneConfig,
    sardine_tif: &Path,
) -> Result<ComparisonStats, ValidationError> {
    let mut cmd = Command::new(&cli.python);
    cmd.arg(&cli.multilook_compare)
        .env("SARDINE", sardine_tif)
        .env("ASF_VV", &scene.asf_reference);
    eprintln!("validate: running {cmd:?}");
    let output = cmd.output()?;
    if !output.status.success() {
        let stderr = String::from_utf8_lossy(&output.stderr);
        let tail: String = stderr.lines().rev().take(40).collect::<Vec<_>>().join("\n");
        return Err(ValidationError::ComparisonFailed {
            exit_code: output.status.code(),
            stderr_tail: tail,
        });
    }
    let stdout = String::from_utf8_lossy(&output.stdout);
    parse_comparison_output(&stdout)
}

// ── Threshold checking ─────────────────────────────────────────────

fn check_thresholds(
    scene: &SceneConfig,
    stats: &ComparisonStats,
) -> Result<(), ValidationError> {
    if stats.median_bias_db.abs() > scene.thresholds.median_bias_db_max {
        return Err(ValidationError::ThresholdViolation {
            scene: scene.name.clone(),
            metric: "|median bias|",
            actual: stats.median_bias_db.abs(),
            threshold: scene.thresholds.median_bias_db_max,
            unit: "dB",
        });
    }
    if stats.std_dev_db > scene.thresholds.std_dev_db_max {
        return Err(ValidationError::ThresholdViolation {
            scene: scene.name.clone(),
            metric: "std dev",
            actual: stats.std_dev_db,
            threshold: scene.thresholds.std_dev_db_max,
            unit: "dB",
        });
    }
    if stats.n_samples < scene.thresholds.min_joint_samples {
        return Err(ValidationError::ThresholdViolation {
            scene: scene.name.clone(),
            metric: "joint samples",
            actual: stats.n_samples as f64,
            threshold: scene.thresholds.min_joint_samples as f64,
            unit: "samples",
        });
    }
    Ok(())
}

// ── Driver ─────────────────────────────────────────────────────────

fn main() -> std::process::ExitCode {
    let cli = Cli::parse();

    let scenes = match load_scenes(&cli.scenes_dir) {
        Ok(s) => s,
        Err(e) => {
            eprintln!("validate: {e}");
            return std::process::ExitCode::from(2);
        }
    };
    if scenes.is_empty() {
        eprintln!("validate: no *.toml scene configs found in {:?}", cli.scenes_dir);
        return std::process::ExitCode::from(2);
    }

    let outdir = cli.outdir.clone().unwrap_or_else(std::env::temp_dir);
    if let Err(e) = std::fs::create_dir_all(&outdir) {
        eprintln!("validate: cannot create outdir {outdir:?}: {e}");
        return std::process::ExitCode::from(2);
    }

    let mut failures: Vec<ValidationError> = Vec::new();
    for (name, (path, scene)) in &scenes {
        if let Some(only) = &cli.only {
            if only != name {
                continue;
            }
        }
        eprintln!("\nvalidate: ── {name} ── ({path:?})");
        let output_tif = outdir.join(format!("sardine_validate_{name}.tif"));
        if let Err(e) = run_sardine(&cli, scene, &output_tif) {
            eprintln!("validate: {name}: SARDINE FAILED: {e}");
            failures.push(e);
            continue;
        }
        let stats = match run_comparison(&cli, scene, &output_tif) {
            Ok(s) => s,
            Err(e) => {
                eprintln!("validate: {name}: COMPARISON FAILED: {e}");
                failures.push(e);
                continue;
            }
        };
        eprintln!(
            "validate: {name}: median_bias = {:+.4} dB, std = {:.4} dB, n = {}",
            stats.median_bias_db, stats.std_dev_db, stats.n_samples,
        );
        if let Err(e) = check_thresholds(scene, &stats) {
            eprintln!("validate: {name}: THRESHOLD VIOLATION: {e}");
            failures.push(e);
        } else {
            eprintln!("validate: {name}: PASS");
        }
    }

    if failures.is_empty() {
        eprintln!("\nvalidate: all scenes PASS ({} run)", scenes.len());
        std::process::ExitCode::SUCCESS
    } else {
        eprintln!("\nvalidate: {} scene(s) FAILED:", failures.len());
        for f in &failures {
            eprintln!("  - {f}");
        }
        std::process::ExitCode::from(1)
    }
}

// ── Tests ──────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    /// Real `multilook_compare.py` output from the post-fix Munich UTM run.
    /// Pinning a sample of the actual output keeps the parser honest if
    /// the Python script ever changes its formatting.
    const REAL_OUTPUT: &str = "\
Reading sardine ...
Reading + reprojecting ASF ...
Block-averaging 10×10 ...

Multi-looked (10×10, ~110 m) comparison: n=4,098,547
  Sardine median: -10.326 dB
  ASF    median:  -10.277 dB
  Mean bias:      +0.0586 dB
  Median bias:    +0.0283 dB
  Std dev:        2.9572 dB  (single-look was 5.37 dB)
  RMSE:           2.9577 dB
";

    #[test]
    fn parse_real_output_extracts_three_metrics() {
        let s = parse_comparison_output(REAL_OUTPUT).expect("parse must succeed");
        assert!((s.median_bias_db - 0.0283).abs() < 1e-9, "median: {}", s.median_bias_db);
        assert!((s.std_dev_db - 2.9572).abs() < 1e-9, "std: {}", s.std_dev_db);
        assert_eq!(s.n_samples, 4_098_547);
    }

    #[test]
    fn parse_missing_median_is_hard_error() {
        let bad = "\
Multi-looked (10×10, ~110 m) comparison: n=1
  Std dev:        1.0 dB
";
        let err = parse_comparison_output(bad).unwrap_err();
        assert!(matches!(err, ValidationError::ParseFailed { field: "median_bias_db", .. }));
    }

    #[test]
    fn parse_missing_n_samples_is_hard_error() {
        let bad = "\
  Median bias:    +0.0 dB
  Std dev:        1.0 dB
";
        let err = parse_comparison_output(bad).unwrap_err();
        assert!(matches!(err, ValidationError::ParseFailed { field: "n_samples", .. }));
    }

    fn munich_thresholds() -> Thresholds {
        Thresholds {
            median_bias_db_max: 0.1,
            std_dev_db_max: 3.5,
            min_joint_samples: 1_000_000,
        }
    }

    fn dummy_scene() -> SceneConfig {
        SceneConfig {
            name: "munich".into(),
            safe: "/tmp".into(),
            orbit: "/tmp".into(),
            dem: "/tmp".into(),
            asf_reference: "/tmp".into(),
            polarization: "VV".into(),
            crs: "auto".into(),
            pixel_spacing_m: Some(10.0),
            pixel_spacing_deg: None,
            geoid: "zero".into(),
            extra_sardine_args: vec![],
            thresholds: munich_thresholds(),
        }
    }

    #[test]
    fn threshold_passes_for_munich_post_fix_numbers() {
        let scene = dummy_scene();
        let stats = ComparisonStats {
            median_bias_db: 0.0283,
            std_dev_db: 2.9572,
            n_samples: 4_098_547,
        };
        check_thresholds(&scene, &stats).expect("post-fix numbers must pass");
    }

    #[test]
    fn threshold_fails_on_excess_bias() {
        let scene = dummy_scene();
        let stats = ComparisonStats {
            median_bias_db: 0.5, // 5x the threshold
            std_dev_db: 2.0,
            n_samples: 4_000_000,
        };
        let err = check_thresholds(&scene, &stats).unwrap_err();
        assert!(matches!(
            err,
            ValidationError::ThresholdViolation { metric: "|median bias|", .. }
        ));
    }

    #[test]
    fn threshold_fails_on_low_sample_count() {
        let scene = dummy_scene();
        let stats = ComparisonStats {
            median_bias_db: 0.0,
            std_dev_db: 2.0,
            n_samples: 100, // far below the floor
        };
        let err = check_thresholds(&scene, &stats).unwrap_err();
        assert!(matches!(
            err,
            ValidationError::ThresholdViolation { metric: "joint samples", .. }
        ));
    }
}
