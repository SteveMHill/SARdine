//! Python bindings for SARdine.
//!
//! Exposes two top-level functions that map 1:1 onto the
//! [`sardine::run`] entry points used by the `sardine` CLI:
//!
//! * [`process`] — full backscatter pipeline (equivalent to `sardine process`).
//! * [`grd`]     — ground-range pipeline (equivalent to `sardine grd`).
//!
//! Both functions take only the four required arguments positionally;
//! every other knob is a keyword argument with the same default the CLI
//! uses.  This keeps the Python surface narrow while preserving full
//! parity with the command-line tool.
//!
//! # Errors
//!
//! Any failure inside the Rust pipeline is converted to a Python
//! `RuntimeError` carrying the full chained `anyhow` message
//! (`{:#}`-formatted), so Python callers see the same diagnostic that
//! the CLI prints to stderr after the `error:` prefix.
//!
//! # Threading
//!
//! The pipeline releases the GIL with [`Python::allow_threads`] before
//! invoking Rust, so other Python threads can run in parallel.

use std::path::PathBuf;

use pyo3::exceptions::PyRuntimeError;
use pyo3::prelude::*;

use sardine::run::{
    parse_iw_selection, parse_output_mode, parse_speckle_order,
    run_grd_multi, run_process_multi, GrdOptions, ProcessOptions, ResamplingKernel,
};

/// Convert an `anyhow::Error` into a Python `RuntimeError` with the
/// full chained message.  Matches the CLI's `eprintln!("error: {:#}")`
/// formatting so Python diagnostics are identical to CLI diagnostics.
fn py_err(e: anyhow::Error) -> PyErr {
    PyRuntimeError::new_err(format!("{e:#}"))
}

/// Run the full backscatter pipeline.
///
/// Equivalent to `sardine process …`.  Required positional arguments:
///
/// * ``safe``    — path to the Sentinel-1 IW SLC ``.SAFE`` directory.
/// * ``dem``     — directory containing SRTM1 GeoTIFF tiles.
/// * ``output``  — output GeoTIFF path (Float32, dB).
/// * ``geoid``   — ``"auto"``, ``"zero"``, or a path to an EGM96 grid file.
///
/// All other parameters are keyword-only and default to the CLI defaults.
/// See [`sardine::run::ProcessOptions`] for field semantics.
#[pyfunction]
#[pyo3(signature = (
    safe,
    dem,
    output,
    geoid,
    *,
    orbit = None,
    polarization = "VV".to_owned(),
    no_flatten = false,
    noise_floor_db = 0.0,
    pixel_spacing_deg = 0.0001,
    pixel_spacing_m = 10.0,
    crs = "wgs84".to_owned(),
    write_mask = false,
    write_lia = false,
    no_provenance = false,
    threads = 0,
    speckle = "none".to_owned(),
    speckle_window = 7,
    enl = 1.0,
    frost_damping = 1.0,
    multilook_range = 1_usize,
    multilook_azimuth = 1_usize,
    extra_safe_paths = vec![],
    mode = "rtc".to_owned(),
    speckle_order = "after".to_owned(),
    iw = "".to_owned(),
    burst_range = None,
))]
#[allow(clippy::too_many_arguments)]
fn process(
    py: Python<'_>,
    safe: PathBuf,
    dem: PathBuf,
    output: PathBuf,
    geoid: String,
    orbit: Option<PathBuf>,
    polarization: String,
    no_flatten: bool,
    noise_floor_db: f32,
    pixel_spacing_deg: f64,
    pixel_spacing_m: f64,
    crs: String,
    write_mask: bool,
    write_lia: bool,
    no_provenance: bool,
    threads: usize,
    speckle: String,
    speckle_window: usize,
    enl: f32,
    frost_damping: f32,
    multilook_range: usize,
    multilook_azimuth: usize,
    extra_safe_paths: Vec<PathBuf>,
    mode: String,
    speckle_order: String,
    iw: String,
    burst_range: Option<String>,
) -> PyResult<()> {
    let iw_selection = parse_iw_selection(&iw, burst_range.as_deref()).map_err(py_err)?;
    let mode = parse_output_mode(&mode).map_err(py_err)?;
    let speckle_order = parse_speckle_order(&speckle_order).map_err(py_err)?;
    let opts = ProcessOptions {
        safe,
        dem: Some(dem),
        dem_source: "srtm1".to_string(),
        output,
        orbit,
        polarization,
        no_flatten,
        noise_floor_db,
        pixel_spacing_deg,
        pixel_spacing_m,
        crs,
        write_mask,
        write_lia,
        no_provenance,
        cog: false,
        threads,
        speckle,
        speckle_window,
        enl,
        frost_damping,
        geoid,
        multilook_range,
        multilook_azimuth,
        extra_safe_paths,
        iw_selection,
        mode,
        speckle_order,
        resampling: ResamplingKernel::default(),
    };
    py.allow_threads(|| run_process_multi(&opts)).map_err(py_err)
}

/// Run the ground-range pipeline.
///
/// Equivalent to `sardine grd …`.  Required positional arguments:
///
/// * ``safe``    — path to the Sentinel-1 IW SLC ``.SAFE`` directory.
/// * ``output``  — output TIFF path (Float32, no CRS — radar geometry).
///
/// All other parameters are keyword-only.  See
/// [`sardine::run::GrdOptions`] for field semantics.
#[pyfunction]
#[pyo3(signature = (
    safe,
    output,
    *,
    orbit = None,
    polarization = "VV".to_owned(),
    target_spacing_m = 10.0,
    no_provenance = false,
    threads = 0,
    speckle = "none".to_owned(),
    speckle_window = 7,
    enl = 1.0,
    frost_damping = 1.0,
    iw = "".to_owned(),
    burst_range = None,
))]
#[allow(clippy::too_many_arguments)]
fn grd(
    py: Python<'_>,
    safe: PathBuf,
    output: PathBuf,
    orbit: Option<PathBuf>,
    polarization: String,
    target_spacing_m: f64,
    no_provenance: bool,
    threads: usize,
    speckle: String,
    speckle_window: usize,
    enl: f32,
    frost_damping: f32,
    iw: String,
    burst_range: Option<String>,
) -> PyResult<()> {
    let iw_selection = parse_iw_selection(&iw, burst_range.as_deref()).map_err(py_err)?;
    let opts = GrdOptions {
        safe,
        output,
        orbit,
        polarization,
        target_spacing_m,
        no_provenance,
        cog: false,
        threads,
        speckle,
        speckle_window,
        enl,
        frost_damping,
        iw_selection,
    };
    py.allow_threads(|| run_grd_multi(&opts)).map_err(py_err)
}

// ─────────────────────────────────────────────────────────────────────────────
// Optional fetch helpers (feature-gated)
// ─────────────────────────────────────────────────────────────────────────────

/// Download a POEORB precise orbit `.EOF` file for a Sentinel-1 SAFE
/// product.
///
/// Equivalent to `sardine fetch-orbit …`.  Requires the wheel to have
/// been built with the ``orbit-fetch`` Cargo feature enabled
/// (e.g. ``maturin develop --features orbit-fetch``).  When the feature
/// is not compiled in, calling this function raises ``RuntimeError``.
///
/// Returns the filesystem path of the downloaded `.EOF` as a string.
#[pyfunction]
#[pyo3(signature = (safe, cache_dir))]
fn fetch_orbit(py: Python<'_>, safe: PathBuf, cache_dir: PathBuf) -> PyResult<String> {
    #[cfg(feature = "orbit-fetch")]
    {
        use sardine::orbit_fetch::fetch_poeorb;
        use sardine::parse::parse_safe_directory;

        py.allow_threads(|| -> anyhow::Result<String> {
            let safe_name = safe
                .file_name()
                .ok_or_else(|| anyhow::anyhow!("SAFE path has no filename: {}", safe.display()))?
                .to_string_lossy()
                .into_owned();
            let product_id = safe_name.trim_end_matches(".SAFE");

            let scene = parse_safe_directory(&safe).map_err(|e| {
                anyhow::anyhow!("parsing SAFE {}: {}", safe.display(), e)
            })?;
            let orbit_path = fetch_poeorb(product_id, scene.start_time, &cache_dir)
                .map_err(|e| anyhow::anyhow!("downloading POEORB: {}", e))?;
            orbit_path
                .to_str()
                .map(str::to_owned)
                .ok_or_else(|| anyhow::anyhow!("orbit path is not valid UTF-8"))
        })
        .map_err(py_err)
    }
    #[cfg(not(feature = "orbit-fetch"))]
    {
        let _ = (py, safe, cache_dir); // SAFETY-OK: args intentionally unused without the feature; matches CLI behaviour
        Err(PyRuntimeError::new_err(
            "fetch_orbit requires the 'orbit-fetch' Cargo feature.\n\
             Rebuild with: maturin develop --release --features orbit-fetch",
        ))
    }
}

/// Download and extract a Sentinel-1 IW SLC SAFE product from ASF datapool.
///
/// Equivalent to `sardine download-slc …`.  Requires the ``slc-fetch``
/// Cargo feature.  If ``token`` is ``None``, the Earthdata Bearer token
/// is read from the ``EARTHDATA_TOKEN`` environment variable
/// (https://urs.earthdata.nasa.gov/user_tokens).
///
/// Returns the filesystem path of the extracted ``{product_id}.SAFE``
/// directory as a string.
#[pyfunction]
#[pyo3(signature = (product_id, output_dir, token = None))]
fn download_slc(
    py: Python<'_>,
    product_id: String,
    output_dir: PathBuf,
    token: Option<String>,
) -> PyResult<String> {
    #[cfg(feature = "slc-fetch")]
    {
        use sardine::slc_fetch::{fetch_slc, token_from_env};

        py.allow_threads(|| -> anyhow::Result<String> {
            let token = match token {
                Some(t) => t,
                None => token_from_env().map_err(|e| {
                    anyhow::anyhow!(
                        "no token provided and EARTHDATA_TOKEN is not set: {}",
                        e
                    )
                })?,
            };
            let safe_path = fetch_slc(&product_id, &output_dir, &token)
                .map_err(|e| anyhow::anyhow!("downloading SLC {}: {}", product_id, e))?;
            safe_path
                .to_str()
                .map(str::to_owned)
                .ok_or_else(|| anyhow::anyhow!("SAFE path is not valid UTF-8"))
        })
        .map_err(py_err)
    }
    #[cfg(not(feature = "slc-fetch"))]
    {
        let _ = (py, product_id, output_dir, token); // SAFETY-OK: args intentionally unused without the feature
        Err(PyRuntimeError::new_err(
            "download_slc requires the 'slc-fetch' Cargo feature.\n\
             Rebuild with: maturin develop --release --features slc-fetch",
        ))
    }
}

/// Locate or download the EGM96 geoid grid into the local cache directory
/// (``$SARDINE_GEOID_DIR`` or ``$HOME/.sardine/geoid/``).
///
/// Requires the ``geoid-fetch`` Cargo feature.  Returns ``None`` on
/// success — the function exists so Python callers can pre-warm the
/// cache and surface network failures *before* invoking [`process`].
/// When the feature is compiled in, ``process(geoid="auto")`` will
/// transparently call the same fetch path.
#[pyfunction]
fn fetch_geoid(py: Python<'_>) -> PyResult<()> {
    #[cfg(feature = "geoid-fetch")]
    {
        use sardine::geoid_fetch::fetch_egm96;

        py.allow_threads(|| -> anyhow::Result<()> {
            fetch_egm96()
                .map(|_grid| ())
                .map_err(|e| anyhow::anyhow!("fetching EGM96 grid: {}", e))
        })
        .map_err(py_err)
    }
    #[cfg(not(feature = "geoid-fetch"))]
    {
        let _ = py; // SAFETY-OK: arg intentionally unused without the feature
        Err(PyRuntimeError::new_err(
            "fetch_geoid requires the 'geoid-fetch' Cargo feature.\n\
             Rebuild with: maturin develop --release --features geoid-fetch",
        ))
    }
}

/// Build-time feature flags compiled into this wheel.  Returned as a
/// dict so Python callers can introspect with e.g.
/// ``sardine.features()["geoid_fetch"]``.
#[pyfunction]
fn features(py: Python<'_>) -> PyResult<Py<PyAny>> {
    use pyo3::types::PyDict;
    let dict = PyDict::new_bound(py);
    dict.set_item("geoid_fetch", cfg!(feature = "geoid-fetch"))?;
    dict.set_item("orbit_fetch", cfg!(feature = "orbit-fetch"))?;
    dict.set_item("slc_fetch", cfg!(feature = "slc-fetch"))?;
    Ok(dict.into_any().unbind())
}

/// Module entry point.  The Python import name (`_sardine`) is set in
/// `Cargo.toml` via `[lib].name`; PyO3 derives the `PyInit__sardine`
/// symbol from that.
#[pymodule]
fn _sardine(_py: Python<'_>, m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(process, m)?)?;
    m.add_function(wrap_pyfunction!(grd, m)?)?;
    m.add_function(wrap_pyfunction!(fetch_orbit, m)?)?;
    m.add_function(wrap_pyfunction!(download_slc, m)?)?;
    m.add_function(wrap_pyfunction!(fetch_geoid, m)?)?;
    m.add_function(wrap_pyfunction!(features, m)?)?;
    m.add("__version__", env!("CARGO_PKG_VERSION"))?;
    Ok(())
}
