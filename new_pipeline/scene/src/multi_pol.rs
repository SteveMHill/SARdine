//! Multi-polarization wrappers around the single-pol `run_process` and
//! `run_grd` entry points.

use anyhow::{Context, Result};

use crate::pipeline_options::{derive_pol_output_path, parse_polarizations, GrdOptions,
    ProcessOptions};

/// Multi-polarization wrapper around [`crate::run::run_process`].
///
/// Accepts the same `ProcessOptions`, but reads `opts.polarization` as a
/// spec parsed by [`parse_polarizations`].  Single-pol input is forwarded
/// to `run_process` unchanged (output path preserved exactly).  Dual-pol
/// input runs the full pipeline once per pol and writes to a derived path
/// per pol via [`derive_pol_output_path`].
///
/// Each per-pol invocation writes its own `.provenance.json`,
/// `.lia.tif`, and `.mask.tif` sidecars next to its own raster.
pub fn run_process_multi(opts: &ProcessOptions) -> Result<()> {
    let pols = parse_polarizations(&opts.polarization)?;
    if pols.len() == 1 {
        // Single-pol fast path: keep the user's output path exactly as given.
        let mut sub = opts.clone();
        sub.polarization = pols.into_iter().next().expect("len()==1");
        return crate::run::run_process(&sub);
    }
    tracing::info!(
        "dual-polarization run — pols = {}; geometry will be recomputed per pol",
        pols.join("+")
    );
    for pol in &pols {
        let mut sub = opts.clone();
        sub.polarization = pol.clone();
        sub.output = derive_pol_output_path(&opts.output, pol)?;
        tracing::info!(
            "--- polarization {} → {} ---",
            pol,
            sub.output.display()
        );
        crate::run::run_process(&sub)
            .with_context(|| format!("processing polarization {pol}"))?;
    }
    Ok(())
}

/// Multi-polarization wrapper around [`crate::run::run_grd`].  Same contract as
/// [`run_process_multi`] but for the ground-range pipeline.
pub fn run_grd_multi(opts: &GrdOptions) -> Result<()> {
    let pols = parse_polarizations(&opts.polarization)?;
    if pols.len() == 1 {
        let mut sub = opts.clone();
        sub.polarization = pols.into_iter().next().expect("len()==1");
        return crate::run::run_grd(&sub);
    }
    tracing::info!(
        "dual-polarization GRD run — pols = {}; per-pol pipeline",
        pols.join("+")
    );
    for pol in &pols {
        let mut sub = opts.clone();
        sub.polarization = pol.clone();
        sub.output = derive_pol_output_path(&opts.output, pol)?;
        tracing::info!(
            "--- polarization {} → {} ---",
            pol,
            sub.output.display()
        );
        crate::run::run_grd(&sub)
            .with_context(|| format!("processing polarization {pol}"))?;
    }
    Ok(())
}
