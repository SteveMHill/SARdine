//! `sardine` — Sentinel-1 backscatter processing CLI.
//!
//! # Subcommands
//! - `process`       — Full pipeline: deburst → calibration → merge → terrain-correction → dB → GeoTIFF
//! - `fetch-orbit`   — Download POEORB for a SAFE product (requires `--features orbit-fetch`)
//! - `download-slc`  — Download a Sentinel-1 IW SLC SAFE from ASF datapool (requires `--features slc-fetch`)
//! - `inspect`       — Print scene metadata (bounding box, polarizations, burst count, timing)
//!
//! All scientific orchestration lives in [`sardine_scene::run`]; this binary
//! is a clap-only shim that maps `*Args` → `*Options` → `run_*`.

use std::path::PathBuf;
use std::process::ExitCode;

use anyhow::{bail, Context, Result};
use clap::{Parser, Subcommand};
use tracing_subscriber::EnvFilter;

use sardine_scene::parse::parse_safe_directory;
use sardine_scene::run::{
    GrdOptions, IwSelection, ProcessOptions,
};

// ─────────────────────────────────────────────────────────────────────────────
// CLI definition
// ─────────────────────────────────────────────────────────────────────────────

#[derive(Parser)]
#[command(
    name = "sardine",
    about = "Sentinel-1 SAR backscatter processing",
    version
)]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {
    /// Run the full backscatter pipeline on a Sentinel-1 IW SLC SAFE product.
    Process(ProcessArgs),

    /// Download a POEORB precise orbit file for a SAFE product.
    ///
    /// Requires the crate to be built with `--features orbit-fetch`.
    FetchOrbit(FetchOrbitArgs),

    /// Download a Sentinel-1 IW SLC SAFE product from ASF DAAC datapool.
    ///
    /// Requires the crate to be built with `--features slc-fetch`.
    /// Authentication: set the `EARTHDATA_TOKEN` environment variable to a
    /// valid Earthdata Bearer token (https://urs.earthdata.nasa.gov/user_tokens).
    DownloadSlc(DownloadSlcArgs),

    /// Print scene metadata for a SAFE product.
    Inspect(InspectArgs),

    /// Produce a ground-range, multi-looked σ⁰ raster in radar geometry
    /// (no DEM, no geocoding, no terrain flattening).  Output is a
    /// non-georeferenced TIFF — radar geometry cannot be honestly
    /// expressed in EPSG:4326.
    Grd(GrdArgs),
}

// ─────────────────────────────────────────────────────────────────────────────
// Process subcommand
// ─────────────────────────────────────────────────────────────────────────────

#[derive(Parser)]
struct ProcessArgs {
    /// Path(s) to Sentinel-1 IW SLC `.SAFE` directories.
    ///
    /// For a single scene, pass one path.  For slice-assembled (multi-frame)
    /// processing, pass multiple paths in ascending temporal order:
    ///
    ///   --safe /path/to/slice1.SAFE --safe /path/to/slice2.SAFE
    ///
    /// Slices must be from the same orbit pass and be temporally consecutive.
    #[arg(long, value_name = "PATH", num_args = 1..)]
    safe: Vec<PathBuf>,

    /// Directory containing SRTM1 DEM GeoTIFF tiles.
    #[arg(long, value_name = "DIR")]
    dem: PathBuf,

    /// Output path for the dB GeoTIFF (written as Float32 GeoTIFF).
    #[arg(long, value_name = "PATH")]
    output: PathBuf,

    /// Path to a POEORB `.EOF` file.
    ///
    /// If omitted, the annotation orbit is used only when the environment
    /// variable `SARDINE_ALLOW_ANNOTATION_ORBIT=1` is set.  Without that
    /// variable, the command exits with an error.
    #[arg(long, value_name = "FILE")]
    orbit: Option<PathBuf>,

    /// Polarization channel(s) to process.  Accepted spellings
    /// (case-insensitive):
    ///
    ///   * Single pol: `VV` (default), `VH`.
    ///   * Dual pol  : `VV+VH`, `VV,VH`, `dual`, `both`.
    ///
    /// For dual pol the pipeline is run once per polarization and the
    /// output path gets a `_VV` / `_VH` tag inserted before the
    /// extension (e.g. `out.tif` → `out_VV.tif` and `out_VH.tif`).
    #[arg(long, value_name = "POL", default_value = "VV")]
    polarization: String,

    /// Disable terrain flattening (Small 2011).  Flattening is on by default.
    #[arg(long)]
    no_flatten: bool,

    /// Noise floor in dB.  Pixels at or below this level are masked to NaN.
    ///
    /// Values ≤ 0.0 dB mean no masking (only truly zero/negative σ⁰ are removed).
    /// Set a negative value to apply thermal noise masking, e.g. `--noise-floor-db -17`
    /// to mask everything below −17 dB (a typical S1 NESZ floor).
    #[arg(long, value_name = "DB", default_value_t = 0.0_f32)]
    noise_floor_db: f32,

    /// Output pixel spacing in degrees (applies to both lat and lon).
    /// Used only when `--crs` resolves to a geographic CRS (the default,
    /// `wgs84`).  For projected CRSs (UTM), use `--pixel-spacing-m` instead.
    #[arg(long, value_name = "DEG", default_value_t = 0.0001_f64)]
    pixel_spacing_deg: f64,

    /// Output pixel spacing in metres.  Used only when `--crs` resolves
    /// to a projected CRS (UTM).  Default 10 m matches the native S-1
    /// IW resolution and the ASF RTC-10 product.
    #[arg(long, value_name = "M", default_value_t = 10.0_f64)]
    pixel_spacing_m: f64,

    /// Output CRS for the geocoded raster.  One of:
    ///   * `wgs84` (default) — EPSG:4326, geographic lat/lon, spacing in degrees.
    ///   * `auto` — pick the UTM zone that contains the scene centre, spacing in metres.
    ///   * `EPSG:NNNNN` — explicit EPSG code (currently only EPSG:4326 and
    ///     UTM zones EPSG:326XX / 327XX are supported).
    #[arg(long, value_name = "SPEC", default_value = "wgs84")]
    crs: String,

    /// Also write a per-pixel quality mask GeoTIFF (uint8) to `<output>.mask.tif`.
    #[arg(long)]
    write_mask: bool,

    /// Also write the local incidence angle cosine to `<output>.lia.tif`
    /// (Float32, NaN where invalid).
    #[arg(long)]
    write_lia: bool,

    /// Suppress the `<output>.provenance.json` sidecar.
    #[arg(long)]
    no_provenance: bool,

    /// Write a Cloud-Optimised GeoTIFF instead of a stripped TIFF.
    ///
    /// Requires `gdal_translate` (GDAL ≥ 3.1) in PATH.
    /// The output file is converted in-place after writing.
    #[arg(long)]
    cog: bool,

    /// Number of CPU threads to use for parallel stages.  0 = Rayon default.
    #[arg(long, value_name = "N", default_value_t = 0_usize)]
    threads: usize,

    /// Speckle filter applied to the linear-power buffer **after**
    /// terrain correction and **before** dB conversion.
    ///
    /// One of: `none`, `boxcar`, `lee`, `frost`, `gamma_map`, `refined_lee`.
    #[arg(long, value_name = "KIND", default_value = "none")]
    speckle: String,

    /// Speckle filter window size (odd integer ≥ 3).
    #[arg(long, value_name = "N", default_value_t = 7_usize)]
    speckle_window: usize,

    /// Effective number of looks fed to Lee / Gamma-MAP.
    #[arg(long, value_name = "L", default_value_t = 1.0_f32)]
    enl: f32,

    /// Frost damping factor `K` (typical 1.0–2.0).
    #[arg(long, value_name = "K", default_value_t = 1.0_f32)]
    frost_damping: f32,

    /// Geoid model used to convert orthometric DEM heights to ellipsoidal
    /// heights for the WGS84 ECEF geocoding solver.  REQUIRED — no default.
    ///
    /// Accepted values: `auto`, `zero`, or a path to an EGM96 grid file
    /// (`.bin` / `.gtx` / `.GRD`).
    #[arg(long, value_name = "SPEC")]
    geoid: String,

    /// Range multilook factor: average N consecutive range samples before
    /// terrain correction.  1 = no-op (default, preserves existing behaviour).
    /// SNAP S1-IW RTC typically uses 4.
    #[arg(long, value_name = "N", default_value_t = 1_usize)]
    multilook_range: usize,

    /// Azimuth multilook factor: average N consecutive azimuth lines before
    /// terrain correction.  1 = no-op (default).
    /// SNAP S1-IW RTC typically uses 1.
    #[arg(long, value_name = "N", default_value_t = 1_usize)]
    multilook_azimuth: usize,

    /// Sub-swaths to process (comma-separated).  Accepted values: `IW1`, `IW2`, `IW3`.
    /// Default: all three.  Example: `--iw IW2` or `--iw IW1,IW2`.
    #[arg(long, value_name = "LIST", default_value = "")]
    iw: String,

    /// Burst range within each selected sub-swath (0-based, inclusive).
    /// Format: `START-END`.  Example: `--burst-range 0-3` selects bursts 0, 1, 2, 3.
    /// Default: all bursts.
    #[arg(long, value_name = "START-END")]
    burst_range: Option<String>,

    /// Where to apply the speckle filter relative to terrain correction.
    /// `after` (default): map geometry, after TC.
    /// `before`: slant-range geometry, before TC.
    #[arg(long, value_name = "ORDER", default_value = "after")]
    speckle_order: String,

    /// Output mode.
    ///   `rtc` (default): terrain-corrected and terrain-flattened σ⁰→γ⁰, geocoded GeoTIFF.
    ///   `grd`: ground-range radar geometry (no DEM, no geocoding, no CRS).
    ///   `nrb`: RTC with mandatory LIA, quality mask, and STAC sidecars (NRB profile).
    #[arg(long, value_name = "MODE", default_value = "rtc")]
    mode: String,
}

// ─────────────────────────────────────────────────────────────────────────────
// Shared CLI helpers
// ─────────────────────────────────────────────────────────────────────────────

/// Thin CLI wrapper around [`sardine_scene::run::parse_iw_selection`].
/// Adds the `--iw` / `--burst-range` flag name to error messages.
fn parse_iw_selection(iw: &str, burst_range: Option<&str>) -> Result<IwSelection> {
    sardine_scene::run::parse_iw_selection(iw, burst_range)
        .with_context(|| format!("--iw / --burst-range: invalid selection (iw={iw:?})"))
}

impl ProcessArgs {
    fn into_options(mut self) -> Result<ProcessOptions> {
        let extra_safe_paths = if self.safe.len() > 1 {
            self.safe.split_off(1)
        } else {
            vec![]
        };
        let primary_safe = self.safe.into_iter().next().expect("clap num_args=1.. guarantees at least one --safe path"); // SAFETY-OK: clap enforces num_args=1.., so at least one element always present

        let iw_selection = parse_iw_selection(&self.iw, self.burst_range.as_deref())?;

        let speckle_order = sardine_scene::run::parse_speckle_order(&self.speckle_order)
            .with_context(|| "--speckle-order")?;

        let mode = sardine_scene::run::parse_output_mode(&self.mode)
            .with_context(|| "--mode")?;

        Ok(ProcessOptions {
            safe: primary_safe,
            extra_safe_paths,
            dem: self.dem,
            output: self.output,
            orbit: self.orbit,
            polarization: self.polarization,
            no_flatten: self.no_flatten,
            noise_floor_db: self.noise_floor_db,
            pixel_spacing_deg: self.pixel_spacing_deg,
            pixel_spacing_m: self.pixel_spacing_m,
            crs: self.crs,
            write_mask: self.write_mask,
            write_lia: self.write_lia,
            no_provenance: self.no_provenance,
            cog: self.cog,
            threads: self.threads,
            speckle: self.speckle,
            speckle_window: self.speckle_window,
            enl: self.enl,
            frost_damping: self.frost_damping,
            geoid: self.geoid,
            multilook_range: self.multilook_range,
            multilook_azimuth: self.multilook_azimuth,
            iw_selection,
            mode,
            speckle_order,
        })
    }
}

fn cmd_process(args: ProcessArgs) -> Result<()> {
    sardine_scene::run::run_process_multi(&args.into_options()?)
}

// ─────────────────────────────────────────────────────────────────────────────
// fetch-orbit subcommand
// ─────────────────────────────────────────────────────────────────────────────

#[derive(Parser)]
struct FetchOrbitArgs {
    /// Path to the Sentinel-1 `.SAFE` directory (used to derive product ID and sensing time).
    #[arg(long, value_name = "PATH")]
    safe: PathBuf,

    /// Directory to write the downloaded `.EOF` file into (created if absent).
    #[arg(long, value_name = "DIR")]
    cache_dir: PathBuf,
}

fn cmd_fetch_orbit(args: FetchOrbitArgs) -> Result<()> {
    #[cfg(feature = "orbit-fetch")]
    {
        use sardine_scene::orbit_fetch::fetch_poeorb;

        let safe_name = args
            .safe
            .file_name()
            .with_context(|| "SAFE path has no filename")?
            .to_string_lossy();
        let product_id = safe_name.trim_end_matches(".SAFE");

        tracing::info!("parsing SAFE metadata to extract sensing start …");
        let scene = parse_safe_directory(&args.safe)
            .with_context(|| format!("parsing SAFE: {}", args.safe.display()))?;
        let sensing_start = scene.start_time;

        tracing::info!(
            "fetching POEORB for {} (sensing_start={}) …",
            product_id,
            sensing_start.format("%Y-%m-%dT%H:%M:%SZ")
        );

        let orbit_path = fetch_poeorb(product_id, sensing_start, &args.cache_dir)
            .with_context(|| "downloading POEORB")?;

        println!("{}", orbit_path.display());
        Ok(())
    }

    #[cfg(not(feature = "orbit-fetch"))]
    {
        let _ = args; // SAFETY-OK: args intentionally unused in non-feature build; no numeric fallback
        bail!(
            "The 'fetch-orbit' subcommand requires the 'orbit-fetch' feature.\n\
             Rebuild with: cargo build --features orbit-fetch"
        );
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// download-slc subcommand
// ─────────────────────────────────────────────────────────────────────────────

#[derive(Parser)]
struct DownloadSlcArgs {
    /// Sentinel-1 product ID (without `.SAFE` or `.zip` suffix).
    #[arg(long, value_name = "PRODUCT_ID")]
    product_id: String,

    /// Directory to place the extracted `.SAFE` into (created if absent).
    #[arg(long, value_name = "DIR")]
    output_dir: PathBuf,

    /// Earthdata Bearer token.  Defaults to the `EARTHDATA_TOKEN` environment variable.
    #[arg(long, value_name = "TOKEN")]
    token: Option<String>,
}

fn cmd_download_slc(args: DownloadSlcArgs) -> Result<()> {
    #[cfg(feature = "slc-fetch")]
    {
        use sardine_scene::slc_fetch::{fetch_slc, token_from_env};

        let token = match args.token {
            Some(t) => t,
            None => token_from_env().with_context(|| {
                "No --token provided and EARTHDATA_TOKEN is not set."
            })?,
        };

        let safe_path = fetch_slc(&args.product_id, &args.output_dir, &token)
            .with_context(|| format!("downloading SLC: {}", args.product_id))?;

        println!("{}", safe_path.display());
        Ok(())
    }

    #[cfg(not(feature = "slc-fetch"))]
    {
        let _ = args; // SAFETY-OK: args intentionally unused in non-feature build; no numeric fallback
        bail!(
            "The 'download-slc' subcommand requires the 'slc-fetch' feature.\n\
             Rebuild with: cargo build --features slc-fetch"
        );
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// inspect subcommand
// ─────────────────────────────────────────────────────────────────────────────

#[derive(Parser)]
struct InspectArgs {
    /// Path to the Sentinel-1 `.SAFE` directory.
    #[arg(value_name = "SAFE")]
    safe: PathBuf,
}

fn cmd_inspect(args: InspectArgs) -> Result<()> {
    let scene = parse_safe_directory(&args.safe)
        .with_context(|| format!("parsing SAFE: {}", args.safe.display()))?;

    let bb = &scene.bounding_box;
    let pols: Vec<String> = scene
        .polarizations
        .iter()
        .map(|p| format!("{}", p))
        .collect();

    println!("Product ID   : {}", scene.product_id);
    println!("Mission      : {}", scene.mission);
    println!("Mode         : {}", scene.acquisition_mode);
    println!("Polarizations: {}", pols.join(", "));
    println!(
        "Sensing start: {}",
        scene.start_time.format("%Y-%m-%dT%H:%M:%S.%6fZ")
    );
    println!(
        "Sensing stop : {}",
        scene.stop_time.format("%Y-%m-%dT%H:%M:%S.%6fZ")
    );
    println!(
        "Duration     : {:.3} s",
        (scene.stop_time - scene.start_time).num_milliseconds() as f64 / 1000.0
    );
    println!(
        "Bounding box : lon=[{:.4}, {:.4}]  lat=[{:.4}, {:.4}]",
        bb.min_lon_deg, bb.max_lon_deg, bb.min_lat_deg, bb.max_lat_deg
    );
    println!("Subswaths    : {}", scene.sub_swaths.len());
    println!("Bursts       : {}", scene.bursts.len());
    println!("Orbit SVs    : {}", scene.orbit.state_vectors.len());
    println!(
        "Radar freq   : {:.6} GHz",
        scene.radar_frequency_hz / 1e9
    );

    Ok(())
}

// ─────────────────────────────────────────────────────────────────────────────
// Grd subcommand
// ─────────────────────────────────────────────────────────────────────────────

#[derive(Parser)]
struct GrdArgs {
    /// Path to the Sentinel-1 IW SLC `.SAFE` directory.
    #[arg(long, value_name = "PATH")]
    safe: PathBuf,

    /// Output path for the ground-range σ⁰ TIFF (Float32, no CRS).
    #[arg(long, value_name = "PATH")]
    output: PathBuf,

    /// Path to a POEORB `.EOF` file.  Same orbit-source policy as
    /// `process`: required unless `SARDINE_ALLOW_ANNOTATION_ORBIT=1`.
    #[arg(long, value_name = "FILE")]
    orbit: Option<PathBuf>,

    /// Polarization channel(s) to process.  Same spec as
    /// `process --polarization` (single `VV` / `VH`, or dual
    /// `VV+VH` / `VV,VH` / `dual` / `both`).  Dual-pol writes one
    /// per-pol output with `_VV` / `_VH` inserted before the extension.
    #[arg(long, value_name = "POL", default_value = "VV")]
    polarization: String,

    /// Target ground-range pixel spacing in metres.
    #[arg(long, value_name = "M", default_value_t = 10.0_f64)]
    target_spacing_m: f64,

    /// Suppress the `<output>.provenance.json` sidecar.
    #[arg(long)]
    no_provenance: bool,

    /// Write a Cloud-Optimised GeoTIFF instead of a stripped TIFF.
    ///
    /// Requires `gdal_translate` (GDAL ≥ 3.1) in PATH.
    /// The output file is converted in-place after writing.
    #[arg(long)]
    cog: bool,

    /// Number of CPU threads to use for parallel stages.  0 = Rayon default.
    #[arg(long, value_name = "N", default_value_t = 0_usize)]
    threads: usize,

    /// Speckle filter applied to the linear-power buffer **after**
    /// ground-range projection and **before** writing the TIFF.
    #[arg(long, value_name = "KIND", default_value = "none")]
    speckle: String,

    /// Speckle filter window size (odd integer ≥ 3).
    #[arg(long, value_name = "N", default_value_t = 7_usize)]
    speckle_window: usize,

    /// Effective number of looks fed to Lee / Gamma-MAP.
    #[arg(long, value_name = "L", default_value_t = 1.0_f32)]
    enl: f32,

    /// Frost damping factor `K` (typical 1.0–2.0).
    #[arg(long, value_name = "K", default_value_t = 1.0_f32)]
    frost_damping: f32,

    /// Sub-swaths to process (comma-separated).  Accepted values: `IW1`, `IW2`, `IW3`.
    /// Default: all three.  Example: `--iw IW2` or `--iw IW1,IW2`.
    #[arg(long, value_name = "LIST", default_value = "")]
    iw: String,

    /// Burst range within each selected sub-swath (0-based, inclusive).
    /// Format: `START-END`.  Example: `--burst-range 0-3` selects bursts 0, 1, 2, 3.
    /// Default: all bursts.
    #[arg(long, value_name = "START-END")]
    burst_range: Option<String>,
}

impl GrdArgs {
    fn into_options(self) -> Result<GrdOptions> {
        let iw_selection = parse_iw_selection(&self.iw, self.burst_range.as_deref())?;
        Ok(GrdOptions {
            safe: self.safe,
            output: self.output,
            orbit: self.orbit,
            polarization: self.polarization,
            target_spacing_m: self.target_spacing_m,
            no_provenance: self.no_provenance,
            cog: self.cog,
            threads: self.threads,
            speckle: self.speckle,
            speckle_window: self.speckle_window,
            enl: self.enl,
            frost_damping: self.frost_damping,
            iw_selection,
        })
    }
}

fn cmd_grd(args: GrdArgs) -> Result<()> {
    sardine_scene::run::run_grd_multi(&args.into_options()?)
}

// ─────────────────────────────────────────────────────────────────────────────
// main
// ─────────────────────────────────────────────────────────────────────────────

fn main() -> ExitCode {
    // Initialise structured logging.  Set RUST_LOG to control verbosity, e.g.:
    //   RUST_LOG=sardine_scene=debug sardine process …
    // Defaults to INFO when RUST_LOG is not set.
    tracing_subscriber::fmt()
        .with_writer(std::io::stderr)
        .with_env_filter(
            EnvFilter::try_from_default_env()
                .unwrap_or_else(|_| EnvFilter::new("info")), // SAFETY-OK: "info" literal is a valid filter
        )
        .init();

    let cli = Cli::parse();
    let result = match cli.command {
        Commands::Process(args) => cmd_process(args),
        Commands::FetchOrbit(args) => cmd_fetch_orbit(args),
        Commands::DownloadSlc(args) => cmd_download_slc(args),
        Commands::Inspect(args) => cmd_inspect(args),
        Commands::Grd(args) => cmd_grd(args),
    };
    match result {
        Ok(()) => ExitCode::SUCCESS,
        Err(e) => {
            tracing::error!("{:#}", e);
            ExitCode::FAILURE
        }
    }
}
