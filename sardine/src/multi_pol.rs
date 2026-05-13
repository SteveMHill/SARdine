//! Multi-polarization wrappers around the single-pol `run_process`, `run_grd`,
//! and `run_insar` entry points.

use anyhow::{Context, Result};

use crate::pipeline_options::{derive_pol_output_path, parse_polarizations, GrdOptions,
    InsarOptions, ProcessOptions};

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

/// Multi-polarization wrapper around [`crate::run::run_insar`].
///
/// Accepts the same [`InsarOptions`], but treats `opts.polarization` as a
/// spec parsed by [`parse_polarizations`].  Single-pol input is forwarded
/// to `run_insar` unchanged (output basename preserved exactly).  Dual-pol
/// input (`"VV+VH"` / `"HH+HV"`) runs the full InSAR pipeline once per
/// polarization; the output basename is suffixed with `_vv` / `_vh` etc.
/// via [`derive_pol_output_path`].
///
/// Each per-pol invocation writes its own coherence, phase, and
/// `.provenance.json` sidecars.
pub fn run_insar_multi(opts: &InsarOptions) -> Result<()> {
    let pols = parse_polarizations(&opts.polarization)?;
    if pols.len() == 1 {
        let mut sub = opts.clone();
        sub.polarization = pols.into_iter().next().expect("len()==1");
        return crate::run::run_insar(&sub);
    }
    tracing::info!(
        "dual-polarization InSAR run — pols = {}; geometry will be recomputed per pol",
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
        crate::run::run_insar(&sub)
            .with_context(|| format!("processing InSAR polarization {pol}"))?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use std::path::PathBuf;
    use crate::pipeline_options::{derive_pol_output_path, parse_polarizations};

    #[test]
    fn parse_polarizations_accepts_dual_insar() {
        let pols = parse_polarizations("VV+VH").unwrap();
        assert_eq!(pols, vec!["VV", "VH"]);
    }

    #[test]
    fn derive_pol_output_path_insar_basename() {
        // No extension on the InSAR output basename.
        let base = PathBuf::from("/tmp/out/scene");
        let vv = derive_pol_output_path(&base, "VV").unwrap();
        let vh = derive_pol_output_path(&base, "VH").unwrap();
        assert_eq!(vv, PathBuf::from("/tmp/out/scene_VV"));
        assert_eq!(vh, PathBuf::from("/tmp/out/scene_VH"));
    }

    #[test]
    fn derive_pol_output_path_insar_tif_extension() {
        let base = PathBuf::from("/tmp/out/scene.tif");
        let hh = derive_pol_output_path(&base, "HH").unwrap();
        assert_eq!(hh, PathBuf::from("/tmp/out/scene_HH.tif"));
    }
}
