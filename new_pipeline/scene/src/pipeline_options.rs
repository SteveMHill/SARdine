//! Plain (non-clap) option structs and resolver helpers for the `process` and
//! `grd` pipelines.
//!
//! All option structs mirror their corresponding `sardine` CLI flags 1:1.
//! Parser and resolver helpers are pure functions with no I/O side-effects
//! (except `resolve_geoid`, which may fetch or load a grid file).

use std::path::{Path, PathBuf};

use anyhow::{anyhow, bail, Context, Result};

use crate::output_crs::OutputCrs;

// ─────────────────────────────────────────────────────────────────────────────
// Option structs
// ─────────────────────────────────────────────────────────────────────────────

/// Sub-swath and burst selection for IW Split.
///
/// If `subswaths` is empty, all three IW subswaths (IW1, IW2, IW3) are
/// processed.  `burst_range`, when set, is applied as a 0-based inclusive
/// range `[start, end]` to the sorted burst list of each selected subswath.
#[derive(Debug, Clone, Default)]
pub struct IwSelection {
    /// Subswaths to process.  Empty = all three.
    pub subswaths: Vec<crate::types::SubSwathId>,
    /// Burst range within each selected subswath (0-based, inclusive).
    /// `None` = all bursts.
    pub burst_range: Option<(usize, usize)>,
}

/// Output mode for the `sardine process` pipeline.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum OutputMode {
    /// Ground-range radar geometry (no geocoding, no CRS).
    Grd,
    /// Terrain-corrected and terrain-flattened (σ⁰ → γ⁰), geocoded.  Default.
    Rtc,
    /// RTC with mandatory LIA, quality mask, and STAC sidecars (NRB profile).
    Nrb,
}

impl Default for OutputMode {
    fn default() -> Self {
        OutputMode::Rtc
    }
}

/// Where in the pipeline to apply the speckle filter.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SpeckleOrder {
    /// Apply speckle in slant-range geometry, before terrain correction.
    BeforeTc,
    /// Apply speckle in map geometry, after terrain correction.  Default.
    AfterTc,
}

impl Default for SpeckleOrder {
    fn default() -> Self {
        SpeckleOrder::AfterTc
    }
}

/// Resampling kernel used when interpolating the merged σ⁰ image during
/// terrain correction (step 6 of the backward Range-Doppler algorithm).
///
/// The default ([`ResamplingKernel::Bilinear`]) is appropriate for most
/// Sentinel-1 backscatter products.  [`ResamplingKernel::Bicubic`] can
/// slightly reduce aliasing at the cost of mild overshoot (Gibbs ringing)
/// near bright point targets and specular returns — use with care on scenes
/// with strong urban or ship targets.
///
/// CLI flag: `--resampling bilinear` (default) or `--resampling bicubic`.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum ResamplingKernel {
    /// 4-tap bilinear interpolation.  Low-pass, no overshoot.  Default.
    #[default]
    Bilinear,
    /// 16-tap bicubic Keys kernel (parameter α = −0.5, also known as
    /// Catmull-Rom / cubic Hermite spline).  Sharper than bilinear but
    /// can introduce ringing artefacts near high-contrast edges.
    Bicubic,
}

/// Plain (non-clap) input set for [`crate::run::run_process`].
///
/// Every field mirrors a `sardine process` CLI flag of the same name.
/// Defaults match the CLI defaults — see [`ProcessOptions::new`].
#[derive(Debug, Clone)]
pub struct ProcessOptions {
    pub safe: PathBuf,
    pub dem: PathBuf,
    pub output: PathBuf,
    pub orbit: Option<PathBuf>,
    pub polarization: String,
    pub no_flatten: bool,
    pub noise_floor_db: f32,
    pub pixel_spacing_deg: f64,
    pub pixel_spacing_m: f64,
    pub crs: String,
    pub write_mask: bool,
    pub write_lia: bool,
    pub no_provenance: bool,
    /// Convert the output GeoTIFF to Cloud-Optimised GeoTIFF in-place after
    /// writing.  Requires `gdal_translate` (GDAL ≥ 3.1) in PATH.
    pub cog: bool,
    pub threads: usize,
    pub speckle: String,
    pub speckle_window: usize,
    pub enl: f32,
    pub frost_damping: f32,
    /// Geoid spec — REQUIRED, no default.  Accepted: `auto`, `zero`, or
    /// a path to an EGM96 grid file (`.bin` / `.gtx` / `.GRD`).
    pub geoid: String,
    /// Range multilook factor applied *before* terrain correction, on the
    /// merged slant-range σ⁰ image.  1 = no-op (default).
    pub multilook_range: usize,
    /// Azimuth multilook factor applied *before* terrain correction, on the
    /// merged slant-range σ⁰ image.  1 = no-op (default).
    pub multilook_azimuth: usize,
    /// Additional SAFE slice paths for multi-slice (assembled) processing.
    ///
    /// When non-empty, all paths `[safe] ++ extra_safe_paths` are assembled
    /// via [`crate::slice_assembly::assemble_slices`] before processing.
    /// Paths must be ordered ascending in time (earliest slice first).
    pub extra_safe_paths: Vec<PathBuf>,
    /// Sub-swath and burst selection (IW Split).  Default: all three IWs, all bursts.
    pub iw_selection: IwSelection,
    /// Output mode: `Grd`, `Rtc` (default), or `Nrb`.
    pub mode: OutputMode,
    /// Where to apply the speckle filter: `BeforeTc` (slant-range) or `AfterTc` (default).
    pub speckle_order: SpeckleOrder,
    /// Resampling kernel for the σ⁰ image interpolation in terrain correction.
    /// Default: [`ResamplingKernel::Bilinear`].  Set to
    /// [`ResamplingKernel::Bicubic`] for sharper output at the risk of
    /// mild ringing near bright point targets.
    pub resampling: ResamplingKernel,
}

impl ProcessOptions {
    /// Construct with the same defaults the CLI uses.  `safe`, `dem`,
    /// `output`, and `geoid` have no CLI default and must be supplied.
    pub fn new(safe: PathBuf, dem: PathBuf, output: PathBuf, geoid: String) -> Self {
        Self {
            safe,
            dem,
            output,
            orbit: None,
            polarization: "VV".to_owned(),
            no_flatten: false,
            noise_floor_db: 0.0,
            pixel_spacing_deg: 0.0001,
            pixel_spacing_m: 10.0,
            crs: "wgs84".to_owned(),
            write_mask: false,
            write_lia: false,
            no_provenance: false,
            cog: false,
            threads: 0,
            speckle: "none".to_owned(),
            speckle_window: 7,
            enl: 1.0,
            frost_damping: 1.0,
            geoid,
            multilook_range: 1,
            multilook_azimuth: 1,
            extra_safe_paths: vec![],
            iw_selection: IwSelection::default(),
            mode: OutputMode::default(),
            speckle_order: SpeckleOrder::default(),
            resampling: ResamplingKernel::default(),
        }
    }

    /// When `mode == OutputMode::Grd`, the effective ground-range target spacing.
    /// Uses `pixel_spacing_m` (the metric spacing field shared with TC mode).
    pub fn target_spacing_m_for_grd(&self) -> f64 {
        self.pixel_spacing_m
    }
}

/// Plain (non-clap) input set for [`crate::run::run_grd`].
///
/// Mirrors the `sardine grd` CLI flags one-to-one.
#[derive(Debug, Clone)]
pub struct GrdOptions {
    pub safe: PathBuf,
    pub output: PathBuf,
    pub orbit: Option<PathBuf>,
    pub polarization: String,
    pub target_spacing_m: f64,
    pub no_provenance: bool,
    /// Convert the output GeoTIFF to Cloud-Optimised GeoTIFF in-place after
    /// writing.  Requires `gdal_translate` (GDAL ≥ 3.1) in PATH.
    pub cog: bool,
    pub threads: usize,
    pub speckle: String,
    pub speckle_window: usize,
    pub enl: f32,
    pub frost_damping: f32,
    /// Sub-swath and burst selection (IW Split).  Default: all three IWs, all bursts.
    pub iw_selection: IwSelection,
}

impl GrdOptions {
    /// Construct with CLI defaults.  `safe` and `output` are required.
    pub fn new(safe: PathBuf, output: PathBuf) -> Self {
        Self {
            safe,
            output,
            orbit: None,
            polarization: "VV".to_owned(),
            target_spacing_m: 10.0,
            no_provenance: false,
            cog: false,
            threads: 0,
            speckle: "none".to_owned(),
            speckle_window: 7,
            enl: 1.0,
            frost_damping: 1.0,
            iw_selection: IwSelection::default(),
        }
    }
}

/// Plain (non-clap) input set for [`crate::run::run_insar`].
///
/// Mirrors the `sardine insar` CLI flags one-to-one.
#[derive(Debug, Clone)]
pub struct InsarOptions {
    /// Path to the reference (master) Sentinel-1 IW SLC `.SAFE` directory.
    pub reference: PathBuf,
    /// Path to the secondary (slave) Sentinel-1 IW SLC `.SAFE` directory.
    pub secondary: PathBuf,
    /// Output basename.  Per-subswath suffixes are appended:
    /// `<output>_iw1_coherence.tif`, `<output>_iw2_coherence.tif`, etc.
    /// If `output_phase = true`, `<output>_iw1_phase.tif` is also written.
    pub output: PathBuf,
    /// Directory containing DEM tiles for terrain correction and geocoding.
    pub dem: PathBuf,
    /// Path to a POEORB `.EOF` file for the reference scene.
    pub reference_orbit: Option<PathBuf>,
    /// Path to a POEORB `.EOF` file for the secondary scene.
    pub secondary_orbit: Option<PathBuf>,
    /// Polarization channel: `"VV"` or `"VH"`.
    pub polarization: String,
    /// Azimuth multi-look factor (coherence window height in SLC lines).
    pub az_looks: usize,
    /// Range multi-look factor (coherence window width in SLC samples).
    pub rg_looks: usize,
    /// When true, also write the wrapped interferometric phase GeoTIFF.
    pub output_phase: bool,
    /// Geoid model spec: `"auto"`, `"zero"`, or a path to an EGM96 grid.
    pub geoid: String,
    /// Output pixel spacing in degrees (WGS84 lat/lon grid).
    pub pixel_spacing_deg: f64,
    /// Number of Rayon threads.  0 = use all available cores.
    pub threads: usize,
}

impl InsarOptions {
    /// Construct with CLI defaults.  `reference`, `secondary`, `output`, and `dem` are required.
    pub fn new(reference: PathBuf, secondary: PathBuf, output: PathBuf, dem: PathBuf) -> Self {
        Self {
            reference,
            secondary,
            output,
            dem,
            reference_orbit: None,
            secondary_orbit: None,
            polarization: "VV".to_owned(),
            az_looks: crate::insar::interferogram::DEFAULT_COH_AZ_LOOKS,
            rg_looks: crate::insar::interferogram::DEFAULT_COH_RG_LOOKS,
            output_phase: false,
            geoid: "auto".to_owned(),
            pixel_spacing_deg: 0.0001,
            threads: 0,
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Resolver helpers
// ─────────────────────────────────────────────────────────────────────────────

/// Build the sidecar output path by stripping the trailing extension of
/// `output` and appending `suffix`.  Example: `out.tif` + `.mask.tif` →
/// `out.mask.tif`; `out` (no extension) + `.lia.tif` → `out.lia.tif`.
pub fn sidecar_path(output: &Path, suffix: &str) -> Result<String> {
    let mut base = output.to_path_buf();
    if base.extension().is_some() {
        base.set_extension("");
    }
    let mut s = base
        .to_str()
        .ok_or_else(|| anyhow!("output path contains non-UTF-8 characters"))?
        .to_owned();
    s.push_str(suffix);
    Ok(s)
}

/// Resolve the `--crs <SPEC>` argument into an [`OutputCrs`].
pub fn resolve_crs(
    spec: &str,
    scene_centre_lon_lat: Option<(f64, f64)>,
) -> Result<OutputCrs> {
    let trimmed = spec.trim();
    if trimmed.is_empty() {
        bail!(
            "--crs SPEC must not be empty.  Use 'wgs84' (default), 'auto', \
             'EPSG:4326', 'EPSG:326NN' (UTM), 'laea' / 'EPSG:3035' \
             (ETRS89 LAEA Europe), or 'webmercator' / 'EPSG:3857'."
        );
    }
    if trimmed.eq_ignore_ascii_case("auto") {
        let (lon, lat) = scene_centre_lon_lat.ok_or_else(|| {
            anyhow!("--crs auto requires a scene centre but none was provided")
        })?;
        return OutputCrs::auto_for_lon_lat(lon, lat).map_err(|e| {
            anyhow!("--crs auto: cannot pick UTM zone for ({lon:.4}, {lat:.4}): {e}")
        });
    }
    if trimmed.eq_ignore_ascii_case("wgs84") {
        return Ok(OutputCrs::Wgs84LatLon);
    }
    OutputCrs::from_spec(trimmed).map_err(|e| anyhow!("--crs {trimmed:?}: {e}"))
}

/// Resolve the `--speckle` knobs into a [`SpeckleFilter`] enum, or
/// `Ok(None)` for `none`.
pub fn resolve_speckle(
    kind: &str,
    window: usize,
    enl: f32,
    damping: f32,
) -> Result<Option<crate::speckle::SpeckleFilter>> {
    use crate::speckle::SpeckleFilter;
    match kind.trim().to_ascii_lowercase().as_str() {
        "none" => Ok(None),
        "boxcar" => Ok(Some(SpeckleFilter::Boxcar { window })),
        "lee" => Ok(Some(SpeckleFilter::Lee { window, enl })),
        "frost" => Ok(Some(SpeckleFilter::Frost { window, damping })),
        "gamma_map" | "gamma-map" | "gammamap" => {
            Ok(Some(SpeckleFilter::GammaMap { window, enl }))
        }
        "refined_lee" | "refined-lee" | "refinedlee" => {
            Ok(Some(SpeckleFilter::RefinedLee { enl }))
        }
        other => bail!(
            "--speckle: unknown filter '{other}'.  Expected one of: \
             none, boxcar, lee, frost, gamma_map, refined_lee."
        ),
    }
}

/// Resolve the `--geoid <SPEC>` argument into a [`GeoidModel`].
pub fn resolve_geoid(spec: &str) -> Result<crate::geoid::GeoidModel> {
    use crate::geoid::{Egm96Grid, GeoidModel};

    let trimmed = spec.trim();
    if trimmed.is_empty() {
        bail!(
            "--geoid SPEC must not be empty.  Use 'auto' to fetch EGM96 \
             (requires --features geoid-fetch), 'egm2008' to fetch EGM2008, \
             'zero' to disable geoid correction explicitly, \
             or a path to an EGM96 binary/ASCII grid file."
        );
    }

    match trimmed.to_ascii_lowercase().as_str() {
        "zero" => {
        tracing::warn!(
            "geoid disabled (--geoid zero). \
             Expect up to ±80 m systematic vertical bias and a horizontal \
             geolocation error of similar magnitude.  Use 'auto' or a \
             grid path for production."
        );
            Ok(GeoidModel::Zero)
        }
        "auto" => {
            #[cfg(feature = "geoid-fetch")]
            {
                use crate::geoid_fetch::fetch_egm96;
                tracing::info!(
                    "locating or downloading EGM96 grid \
                     ($SARDINE_GEOID_DIR or ~/.sardine/geoid/) …"
                );
                let grid = fetch_egm96()
                    .with_context(|| "fetching EGM96 grid for --geoid auto")?;
                Ok(GeoidModel::Egm96(grid))
            }
            #[cfg(not(feature = "geoid-fetch"))]
            {
                bail!(
                    "--geoid auto requires the 'geoid-fetch' Cargo feature.\n\
                     Rebuild with: cargo build --release --features geoid-fetch\n\
                     Or pass an explicit grid path: --geoid /path/to/egm96.bin"
                );
            }
        }
        "egm2008" => {
            #[cfg(feature = "geoid-fetch")]
            {
                use crate::geoid_fetch::fetch_egm2008;
                tracing::info!(
                    "locating or downloading EGM2008 grid \
                     ($SARDINE_GEOID_DIR or ~/.sardine/geoid/) …"
                );
                let grid = fetch_egm2008()
                    .with_context(|| "fetching EGM2008 grid for --geoid egm2008")?;
                Ok(GeoidModel::Egm2008(grid))
            }
            #[cfg(not(feature = "geoid-fetch"))]
            {
                bail!(
                    "--geoid egm2008 requires the 'geoid-fetch' Cargo feature.\n\
                     Rebuild with: cargo build --release --features geoid-fetch\n\
                     Or pass an explicit grid path: --geoid /path/to/egm2008.bin"
                );
            }
        }
        _ => {
            let path = std::path::Path::new(trimmed);
            if !path.exists() {
                bail!(
                    "--geoid: '{}' is neither 'auto', 'zero', nor an existing file path",
                    trimmed
                );
            }
            let ext = path
                .extension()
                .and_then(|e| e.to_str())
                .map(|s| s.to_ascii_lowercase());
            let grid = match ext.as_deref() {
                Some("bin") => Egm96Grid::load_binary(path).with_context(|| {
                    format!("loading EGM96 binary cache: {}", path.display())
                })?,
                Some("gtx") | Some("grd") => Egm96Grid::from_ww15mgh(path).with_context(|| {
                    format!("loading NGA WW15MGH grid: {}", path.display())
                })?,
                _ => bail!(
                    "--geoid path '{}' has unrecognised extension; expected \
                     .bin (compact cache) or .gtx / .GRD (NGA WW15MGH)",
                    path.display()
                ),
            };
            Ok(GeoidModel::Egm96(grid))
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// CLI-argument parsers
// ─────────────────────────────────────────────────────────────────────────────

/// Split a polarization spec into a deduplicated, validated list of
/// polarization tokens.
///
/// Accepted spellings (case-insensitive):
///
/// * Single pol: `"VV"`, `"VH"`.
/// * Dual pol  : `"VV+VH"`, `"VV,VH"`, `"vh+vv"`, `"dual"`, `"both"`.
///
/// `dual` / `both` always expand to `["VV", "VH"]`.  Any other spelling
/// preserves the input order so the caller can choose the output order.
/// Duplicate tokens (e.g. `"VV+VV"`) are silently deduplicated.
///
/// # Errors
///
/// Returns an error for empty input, unknown tokens, or more than two
/// distinct pols (S-1 IW SLC carries at most two).
pub fn parse_polarizations(spec: &str) -> Result<Vec<String>> {
    let upper = spec.to_uppercase();
    let raw_tokens: Vec<String> = match upper.as_str() {
        "DUAL" | "BOTH" => vec!["VV".to_owned(), "VH".to_owned()],
        _ => upper
            .split(|c: char| c == ',' || c == '+')
            .map(|s| s.trim().to_owned())
            .collect(),
    };

    let mut out: Vec<String> = Vec::new();
    for tok in raw_tokens {
        if tok.is_empty() {
            bail!("polarization spec contains an empty token: {:?}", spec);
        }
        if tok != "VV" && tok != "VH" {
            bail!(
                "Unsupported polarization token: {:?} in spec {:?}. \
                 Use VV, VH, VV+VH, dual, or both.",
                tok,
                spec
            );
        }
        if !out.iter().any(|x| *x == tok) {
            out.push(tok);
        }
    }
    if out.is_empty() {
        bail!("polarization spec is empty: {:?}", spec);
    }
    if out.len() > 2 {
        bail!(
            "Sentinel-1 IW SLC carries at most 2 polarizations; spec {:?} \
             yielded {} distinct tokens",
            spec,
            out.len()
        );
    }
    Ok(out)
}

/// Insert a polarization tag before the extension of an output path.
///
/// `out.tif` + `"VV"` → `out_VV.tif`, `out` + `"VH"` → `out_VH`,
/// `dir/sub/x.tiff` + `"VV"` → `dir/sub/x_VV.tiff`.
///
/// # Errors
///
/// Errors if the path has no filename or contains non-UTF-8 components.
pub fn derive_pol_output_path(output: &Path, pol: &str) -> Result<PathBuf> {
    let stem = output
        .file_stem()
        .ok_or_else(|| anyhow!("output path has no filename: {}", output.display()))?
        .to_str()
        .ok_or_else(|| {
            anyhow!(
                "output path filename is not valid UTF-8: {}",
                output.display()
            )
        })?;
    let new_name = match output.extension().and_then(|e| e.to_str()) {
        Some(ext) => format!("{stem}_{pol}.{ext}"),
        None => format!("{stem}_{pol}"),
    };
    Ok(match output.parent() {
        Some(p) if !p.as_os_str().is_empty() => p.join(new_name),
        _ => PathBuf::from(new_name),
    })
}

/// Parse an IW sub-swath selection string into an [`IwSelection`].
///
/// `iw_spec` is a comma-separated list of sub-swath names or empty for
/// "all three":
///
/// ```text
/// ""          → all three sub-swaths, all bursts (default)
/// "IW1"       → IW1 only, all bursts
/// "IW2,IW3"   → IW2 + IW3, all bursts
/// ```
///
/// `burst_range` is an optional `"START-END"` string (0-based, inclusive)
/// that limits the burst range within each selected sub-swath.
///
/// # Errors
///
/// Returns an error for unknown sub-swath names or malformed burst-range strings.
pub fn parse_iw_selection(iw_spec: &str, burst_range: Option<&str>) -> Result<IwSelection> {
    use crate::types::SubSwathId;

    let subswaths: Vec<SubSwathId> = if iw_spec.trim().is_empty() {
        vec![]
    } else {
        iw_spec
            .split(',')
            .map(|s| match s.trim().to_uppercase().as_str() {
                "IW1" => Ok(SubSwathId::IW1),
                "IW2" => Ok(SubSwathId::IW2),
                "IW3" => Ok(SubSwathId::IW3),
                other => Err(anyhow!(
                    "unknown sub-swath '{}'; expected IW1, IW2, or IW3",
                    other
                )),
            })
            .collect::<Result<Vec<_>>>()?
    };

    let burst_range = match burst_range {
        None => None,
        Some(s) => {
            let parts: Vec<&str> = s.splitn(2, '-').collect();
            if parts.len() != 2 {
                bail!(
                    "burst_range must be START-END (e.g. '0-5'), got '{}'",
                    s
                );
            }
            let start: usize = parts[0]
                .parse()
                .with_context(|| {
                    format!("parsing burst range start '{}': not a number", parts[0])
                })?;
            let end: usize = parts[1]
                .parse()
                .with_context(|| {
                    format!("parsing burst range end '{}': not a number", parts[1])
                })?;
            if end < start {
                bail!(
                    "burst_range end ({}) must be >= start ({})",
                    end,
                    start
                );
            }
            Some((start, end))
        }
    };

    Ok(IwSelection {
        subswaths,
        burst_range,
    })
}

/// Parse an output-mode string into [`OutputMode`].
///
/// Accepted values (case-insensitive): `"rtc"` (default), `"grd"`, `"nrb"`.
///
/// # Errors
///
/// Returns an error for unrecognised values.
pub fn parse_output_mode(s: &str) -> Result<OutputMode> {
    match s.to_ascii_lowercase().as_str() {
        "rtc" => Ok(OutputMode::Rtc),
        "grd" => Ok(OutputMode::Grd),
        "nrb" => Ok(OutputMode::Nrb),
        other => bail!(
            "unknown mode '{}'; expected 'rtc', 'grd', or 'nrb'",
            other
        ),
    }
}

/// Parse a speckle-order string into [`SpeckleOrder`].
///
/// Accepted values (case-insensitive): `"after"` (default), `"before"`.
///
/// # Errors
///
/// Returns an error for unrecognised values.
pub fn parse_speckle_order(s: &str) -> Result<SpeckleOrder> {
    match s.to_ascii_lowercase().as_str() {
        "after" => Ok(SpeckleOrder::AfterTc),
        "before" => Ok(SpeckleOrder::BeforeTc),
        other => bail!(
            "unknown speckle_order '{}'; expected 'before' or 'after'",
            other
        ),
    }
}
