//! `sardine` — Sentinel-1 backscatter processing CLI.
//!
//! # Subcommands
//! - `process`       — Full pipeline: deburst → calibration → merge → terrain-correction → dB → GeoTIFF
//! - `insar`         — InSAR coherence (and optionally phase) for an SLC pair (no geocoding)
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
    GrdOptions, InsarOptions, IwSelection, ProcessOptions,
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

    /// Produce coherence (and optionally wrapped phase) for an InSAR pair.
    ///
    /// Processes all three IW sub-swaths.  Output files are written with
    /// per-subswath suffixes (`_IW1_coherence.tif`, etc.) appended to the
    /// `--output` basename.  No geocoding is performed; output is in
    /// slant-range radar geometry (non-georeferenced TIFF).
    Insar(InsarArgs),

    /// Process multiple scenes from a JSON batch file.
    ///
    /// Each scene is processed sequentially.  A failure in one scene is
    /// recorded and reported at the end; it does not abort later scenes.
    ///
    /// # Batch file format
    ///
    /// A JSON array of objects.  Each object has the same fields as the
    /// `process` subcommand flags (snake_case, types as below).  Fields
    /// with defaults may be omitted.
    ///
    /// Required fields per entry:
    ///   `safe` (string or array of strings), `dem`, `output`, `geoid`
    ///
    /// Example:
    /// ```json
    /// [
    ///   {
    ///     "safe": "/data/S1A.SAFE",
    ///     "dem": "/data/dem",
    ///     "output": "/out/scene1.tif",
    ///     "geoid": "auto"
    ///   }
    /// ]
    /// ```
    Batch(BatchArgs),
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

    /// Directory containing SRTM-1 DEM tiles.
    ///
    /// **Omit for automatic download (default behaviour).**  When not
    /// supplied, the pipeline downloads tiles from the source selected by
    /// `--dem-source` and caches them under `$SARDINE_DEM_DIR`
    /// (or `$HOME/.sardine/dem/`) for future runs.  An explicit path always
    /// takes precedence over the download cache.
    #[arg(long, value_name = "DIR")]
    dem: Option<PathBuf>,

    /// DEM source for automatic tile download.  Ignored when `--dem` is set.
    ///
    /// Accepted values:
    ///   - `srtm1` (default): SRTM 1-arc-second from AWS elevation-tiles-prod.
    ///     EGM96 height reference.  `.hgt` tiles, 3601 × 3601.
    ///   - `glo30`: Copernicus GLO-30 from AWS Open Data.
    ///     EGM2008 height reference.  Cloud-Optimised GeoTIFF, 3600 × 3600.
    ///     Better quality at high latitudes and in mountainous areas.
    #[arg(long, value_name = "SOURCE", default_value = "srtm1")]
    dem_source: String,

    /// Output path for the dB GeoTIFF (written as Float32 GeoTIFF).
    #[arg(long, value_name = "PATH")]
    output: PathBuf,

    /// Path to a POEORB `.EOF` orbit file.
    ///
    /// **Omit for automatic download (default behaviour).**  When not
    /// supplied, the pipeline fetches the matching POEORB from the ASF
    /// AWS Open Data bucket and caches it under `$SARDINE_ORBIT_DIR`
    /// (or `$HOME/.sardine/orbits/`) for future runs.  An explicit path
    /// always takes precedence over the download cache.
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

    /// Resampling kernel for the σ⁰ image interpolation during terrain correction.
    ///   `bilinear` (default): 4-tap bilinear — low-pass, no overshoot.
    ///   `bicubic`           : 16-tap Keys kernel (α=−0.5) — sharper but
    ///                         can ring near bright point targets.
    #[arg(long, value_name = "KERNEL", default_value = "bilinear")]
    resampling: String,
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

/// Parse `--resampling bilinear|bicubic`.
fn parse_resampling(s: &str) -> Result<sardine_scene::run::ResamplingKernel> {
    use sardine_scene::run::ResamplingKernel;
    match s.trim().to_ascii_lowercase().as_str() {
        "bilinear" => Ok(ResamplingKernel::Bilinear),
        "bicubic"  => Ok(ResamplingKernel::Bicubic),
        other      => anyhow::bail!(
            "unknown resampling kernel '{other}'; valid values: bilinear, bicubic"
        ),
    }
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

        let resampling = parse_resampling(&self.resampling)
            .with_context(|| "--resampling")?;

        Ok(ProcessOptions {
            safe: primary_safe,
            extra_safe_paths,
            dem: self.dem,
            dem_source: self.dem_source.clone(),
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
            resampling,
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

    /// Path to a POEORB `.EOF` orbit file.
    ///
    /// **Omit for automatic download (default behaviour).**  Same
    /// caching policy as `process --orbit`.
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
// insar subcommand
// ─────────────────────────────────────────────────────────────────────────────

#[derive(Parser)]
struct InsarArgs {
    /// Path to the reference (master) Sentinel-1 IW SLC `.SAFE` directory.
    #[arg(long, value_name = "PATH")]
    reference: PathBuf,

    /// Path to the secondary (slave) Sentinel-1 IW SLC `.SAFE` directory.
    #[arg(long, value_name = "PATH")]
    secondary: PathBuf,

    /// Output basename.  Per-subswath suffixes are appended automatically:
    ///   `<output>_iw1_coherence.tif`, `<output>_iw2_coherence.tif`, etc.
    ///
    /// If `--output-phase` is set, `<output>_iw1_phase.tif` etc. are also written.
    #[arg(long, value_name = "PATH")]
    output: PathBuf,

    /// Directory containing DEM tiles for geocoding (SRTM-1 or compatible).
    ///
    /// **Omit for automatic download (default behaviour).**  When not
    /// supplied, the pipeline downloads tiles from the source selected by
    /// `--dem-source` and caches them under `$SARDINE_DEM_DIR`
    /// (or `$HOME/.sardine/dem/`).
    #[arg(long, value_name = "DIR")]
    dem: Option<PathBuf>,

    /// DEM source for automatic tile download.  Ignored when `--dem` is set.
    ///
    /// Accepted values: `srtm1` (default), `glo30`.  See `sardine process
    /// --help` for full descriptions.
    #[arg(long, value_name = "SOURCE", default_value = "srtm1")]
    dem_source: String,

    /// Path to a POEORB `.EOF` orbit file for the reference scene.
    ///
    /// **Omit for automatic download (default behaviour).**  When not
    /// supplied, the pipeline fetches the matching POEORB from the ASF
    /// AWS Open Data bucket and caches it under `$SARDINE_ORBIT_DIR`
    /// (or `$HOME/.sardine/orbits/`).  An explicit path always takes
    /// precedence over the download cache.
    #[arg(long, value_name = "FILE")]
    reference_orbit: Option<PathBuf>,

    /// Path to a POEORB `.EOF` orbit file for the secondary scene.
    ///
    /// **Omit for automatic download (default behaviour).**  Same
    /// caching policy as `--reference-orbit`.
    #[arg(long, value_name = "FILE")]
    secondary_orbit: Option<PathBuf>,

    /// Polarization channel: `VV` (default) or `VH`.
    #[arg(long, value_name = "POL", default_value = "VV")]
    polarization: String,

    /// Azimuth window size (lines) for coherence estimation.
    #[arg(long, value_name = "N", default_value_t = sardine_scene::insar::interferogram::DEFAULT_COH_AZ_LOOKS)]
    az_looks: usize,

    /// Range window size (samples) for coherence estimation.
    #[arg(long, value_name = "N", default_value_t = sardine_scene::insar::interferogram::DEFAULT_COH_RG_LOOKS)]
    rg_looks: usize,

    /// Also write the wrapped interferometric phase as `<output>_iw{n}_phase.tif`.
    ///
    /// **Note**: Phase correctness requires validated flat-earth deramping.
    /// For coherence-only products, omit this flag.
    #[arg(long)]
    output_phase: bool,

    /// Geoid model for converting DEM orthometric heights to WGS84 ellipsoidal
    /// heights.  Use `auto` (download EGM96 on first run), `zero` (no
    /// correction — ocean scenes only), or a path to an EGM96 `.bin` grid.
    ///
    /// Omitting the geoid correction introduces a geolocation error equal to
    /// the local geoid undulation (typically 30–50 m at mid-latitudes).
    #[arg(long, value_name = "SPEC", default_value = "auto")]
    geoid: String,

    /// Output pixel spacing in degrees (WGS84 lat/lon).  Approximately 10 m
    /// ≈ 0.0001°; 20 m ≈ 0.0002°.
    #[arg(long, value_name = "DEG", default_value_t = 0.0001_f64)]
    pixel_spacing_deg: f64,

    /// Number of processing threads.  0 = use all available CPU cores (default).
    #[arg(long, value_name = "N", default_value_t = 0usize)]
    threads: usize,
}

fn cmd_insar(args: InsarArgs) -> Result<()> {
    let opts = InsarOptions {
        reference: args.reference,
        secondary: args.secondary,
        output: args.output,
        dem: args.dem,
        dem_source: args.dem_source,
        reference_orbit: args.reference_orbit,
        secondary_orbit: args.secondary_orbit,
        polarization: args.polarization,
        az_looks: args.az_looks,
        rg_looks: args.rg_looks,
        output_phase: args.output_phase,
        geoid: args.geoid,
        pixel_spacing_deg: args.pixel_spacing_deg,
        threads: args.threads,
    };
    sardine_scene::run::run_insar(&opts)
}

// ─────────────────────────────────────────────────────────────────────────────
// batch subcommand
// ─────────────────────────────────────────────────────────────────────────────

#[derive(Parser)]
struct BatchArgs {
    /// Path to the JSON batch file (array of process-spec objects).
    #[arg(value_name = "FILE")]
    file: PathBuf,
}

/// Mirror of `ProcessArgs` fields, owned as plain strings / numbers for
/// JSON deserialisation.  All fields with CLI defaults are `Option<T>` so
/// they can be omitted from the JSON file and the defaults will be applied.
///
/// Fields that have no CLI default (`safe`, `dem`, `output`, `geoid`) are
/// non-optional.  `safe` may be a JSON string (single SAFE path) or a JSON
/// array of strings (multiple slices for slice-assembly).
#[derive(serde::Deserialize)]
struct BatchEntry {
    /// Primary (or only) SAFE path.  May be a string or array of strings
    /// (first element = primary, remainder = extra slices, ascending time order).
    safe: OneOrMany,
    /// Directory containing SRTM-1 DEM tiles.  Omit for automatic download.
    #[serde(default)]
    dem: Option<PathBuf>,
    /// DEM source for automatic download: `"srtm1"` (default) or `"glo30"`.
    #[serde(default = "default_dem_source")]
    dem_source: String,
    /// Output path for the dB GeoTIFF.
    output: PathBuf,
    /// Geoid: `"auto"`, `"zero"`, or path to EGM96 `.bin`/`.gtx`/`.GRD` file.
    geoid: String,
    // --- optional fields (same defaults as the CLI) ---
    #[serde(default)]
    orbit: Option<PathBuf>,
    #[serde(default = "default_polarization")]
    polarization: String,
    #[serde(default)]
    no_flatten: bool,
    #[serde(default)]
    noise_floor_db: f32,
    #[serde(default = "default_pixel_spacing_deg")]
    pixel_spacing_deg: f64,
    #[serde(default = "default_pixel_spacing_m")]
    pixel_spacing_m: f64,
    #[serde(default = "default_crs")]
    crs: String,
    #[serde(default)]
    write_mask: bool,
    #[serde(default)]
    write_lia: bool,
    #[serde(default)]
    no_provenance: bool,
    #[serde(default)]
    cog: bool,
    #[serde(default)]
    threads: usize,
    #[serde(default = "default_speckle")]
    speckle: String,
    #[serde(default = "default_speckle_window")]
    speckle_window: usize,
    #[serde(default = "default_enl")]
    enl: f32,
    #[serde(default = "default_frost_damping")]
    frost_damping: f32,
    #[serde(default = "default_multilook")]
    multilook_range: usize,
    #[serde(default = "default_multilook")]
    multilook_azimuth: usize,
    #[serde(default = "default_iw")]
    iw: String,
    #[serde(default)]
    burst_range: Option<String>,
    #[serde(default = "default_mode")]
    mode: String,
    #[serde(default = "default_speckle_order")]
    speckle_order: String,
    #[serde(default = "default_resampling")]
    resampling: String,
}

/// Accepts either a JSON string `"path"` or a JSON array `["path1", "path2"]`.
#[derive(serde::Deserialize)]
#[serde(untagged)]
enum OneOrMany {
    One(PathBuf),
    Many(Vec<PathBuf>),
}

impl OneOrMany {
    fn into_parts(self) -> (PathBuf, Vec<PathBuf>) {
        match self {
            OneOrMany::One(p) => (p, vec![]),
            OneOrMany::Many(mut v) => {
                if v.is_empty() {
                    // Caller will get a descriptive error from the pipeline,
                    // not a silent success.
                    (PathBuf::new(), vec![])
                } else {
                    let primary = v.remove(0);
                    (primary, v)
                }
            }
        }
    }
}

// Default-value functions required by `#[serde(default = "...")]`.
fn default_polarization() -> String { "VV".to_owned() }
fn default_pixel_spacing_deg() -> f64 { 0.0001 }
fn default_pixel_spacing_m() -> f64 { 10.0 }
fn default_crs() -> String { "wgs84".to_owned() }
fn default_dem_source() -> String { "srtm1".to_owned() }
fn default_speckle() -> String { "none".to_owned() }
fn default_speckle_window() -> usize { 7 }
fn default_enl() -> f32 { 1.0 }
fn default_frost_damping() -> f32 { 1.0 }
fn default_multilook() -> usize { 1 }
fn default_iw() -> String { String::new() }
fn default_mode() -> String { "rtc".to_owned() }
fn default_speckle_order() -> String { "after-tc".to_owned() }
fn default_resampling() -> String { "bilinear".to_owned() }

impl BatchEntry {
    fn into_options(self) -> Result<ProcessOptions> {
        let (primary_safe, extra_safe_paths) = self.safe.into_parts();
        if primary_safe.as_os_str().is_empty() {
            anyhow::bail!("batch entry has an empty 'safe' array — at least one SAFE path is required");
        }

        let iw_selection = parse_iw_selection(&self.iw, self.burst_range.as_deref())?;
        let speckle_order = sardine_scene::run::parse_speckle_order(&self.speckle_order)
            .with_context(|| "speckle_order")?;
        let mode = sardine_scene::run::parse_output_mode(&self.mode)
            .with_context(|| "mode")?;
        let resampling = parse_resampling(&self.resampling)
            .with_context(|| "resampling")?;

        Ok(ProcessOptions {
            safe: primary_safe,
            extra_safe_paths,
            dem: self.dem,
            dem_source: self.dem_source,
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
            resampling,
        })
    }
}

fn cmd_batch(args: BatchArgs) -> Result<()> {
    let json_text = std::fs::read_to_string(&args.file)
        .with_context(|| format!("reading batch file: {}", args.file.display()))?;

    let entries: Vec<BatchEntry> = serde_json::from_str(&json_text)
        .with_context(|| format!("parsing batch file: {}", args.file.display()))?;

    if entries.is_empty() {
        anyhow::bail!("batch file contains an empty array — nothing to do");
    }

    let opts: Vec<ProcessOptions> = entries
        .into_iter()
        .enumerate()
        .map(|(i, e)| {
            e.into_options()
                .with_context(|| format!("batch entry {} (output: {}): invalid options", i, ""))
        })
        .collect::<Result<Vec<_>>>()
        .with_context(|| "one or more batch entries have invalid options")?;

    let results = sardine_scene::run::run_process_batch(&opts);

    let mut failed = 0usize;
    for (i, result) in results.into_iter().enumerate() {
        match result {
            Ok(()) => tracing::info!("batch scene {} OK", i + 1),
            Err(e) => {
                tracing::error!("batch scene {} FAILED: {:#}", i + 1, e);
                failed += 1;
            }
        }
    }

    if failed > 0 {
        anyhow::bail!("{} of {} batch scene(s) failed", failed, opts.len());
    }
    Ok(())
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
        Commands::Insar(args) => cmd_insar(args),
        Commands::Batch(args) => cmd_batch(args),
    };
    match result {
        Ok(()) => ExitCode::SUCCESS,
        Err(e) => {
            tracing::error!("{:#}", e);
            ExitCode::FAILURE
        }
    }
}
