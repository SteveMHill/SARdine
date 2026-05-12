//! End-to-end speckle-filter test on the S1B 2019-01-23 SAFE.
//!
//! Verifies that, on a real geocoded γ⁰ buffer (linear power), each of
//! the four speckle filters preserves the contract that the unit tests
//! cannot exercise on synthetic data:
//!
//! 1. Output buffer is the same shape as input.
//! 2. NaN-input pixels remain NaN in the output (nodata preservation).
//! 3. The scene-wide finite-pixel **linear** mean is preserved within
//!    2 % (a bias-free filter cannot shift the local mean by more than
//!    speckle-rms / √N for N independent looks; with ~10⁶ valid pixels
//!    the bound is far below 2 %).
//! 4. Spatial variance over a homogeneous patch (selected from the
//!    interior to avoid edges) is **reduced** by the filter relative to
//!    the unfiltered reference.  Boxcar provides the strongest reduction;
//!    edge-preserving filters less so but still > 0.
//!
//! The test is `#[ignore]`-gated and feature-gated on `geoid-fetch`
//! because (a) it runs the full TC pipeline (~minutes) and (b) it needs
//! the EGM96 grid.  Opt-in:
//!
//! ```sh
//! cargo test --release --features geoid-fetch \
//!     --test speckle_e2e -- --ignored --nocapture
//! ```
//!
//! The test skips cleanly (with an `eprintln!`) if any of the input
//! files (SAFE / POEORB / DEM directory) are not present, so it does
//! not break CI on machines without the data fixture.

#![cfg(feature = "geoid-fetch")]

use std::path::Path;

use sardine_scene::apply_calibration::apply_calibration;
use sardine_scene::calibration::parse_calibration_noise;
use sardine_scene::dem::DemMosaic;
use sardine_scene::deburst::deburst_subswath;
use sardine_scene::geoid::GeoidModel;
use sardine_scene::geoid_fetch::fetch_egm96;
use sardine_scene::merge_subswaths::{merge_subswaths, SwathInput};
use sardine_scene::orbit::{apply_precise_orbit, parse_eof_file};
use sardine_scene::parse::{parse_geolocation_grids, parse_safe_directory};
use sardine_scene::slc_reader::SlcReader;
use sardine_scene::speckle::{apply_speckle_filter, SpeckleFilter};
use sardine_scene::terrain_correction::{terrain_correction, TerrainCorrectionConfig};
use sardine_scene::types::{Polarization, SubSwathId, SubSwathMetadata};

const SAFE: &str = "/home/datacube/dev/SARdine/data/SLC/\
    S1B_IW_SLC__1SDV_20190123T053348_20190123T053415_014617_01B3D4_E833.SAFE";

const TIFFS: [&str; 3] = [
    "/home/datacube/dev/SARdine/data/SLC/\
     S1B_IW_SLC__1SDV_20190123T053348_20190123T053415_014617_01B3D4_E833.SAFE/\
     measurement/s1b-iw1-slc-vv-20190123t053349-20190123t053414-014617-01b3d4-004.tiff",
    "/home/datacube/dev/SARdine/data/SLC/\
     S1B_IW_SLC__1SDV_20190123T053348_20190123T053415_014617_01B3D4_E833.SAFE/\
     measurement/s1b-iw2-slc-vv-20190123t053350-20190123t053415-014617-01b3d4-005.tiff",
    "/home/datacube/dev/SARdine/data/SLC/\
     S1B_IW_SLC__1SDV_20190123T053348_20190123T053415_014617_01B3D4_E833.SAFE/\
     measurement/s1b-iw3-slc-vv-20190123t053348-20190123t053413-014617-01b3d4-006.tiff",
];
const IW_IDS: [SubSwathId; 3] = [SubSwathId::IW1, SubSwathId::IW2, SubSwathId::IW3];

const EOF_FILE: &str = "/home/datacube/dev/SARdine/data/SLC/\
    S1B_IW_SLC__1SDV_20190123T053348_20190123T053415_014617_01B3D4_E833.SAFE/orbit_cache/\
    S1B_IW_SLC__1SDV_20190123T053348_20190123T053415_014617_01B3D4_E833_POEORB.EOF";

const DEM_DIR: &str = "/home/datacube/dev/SARdine/data/dem/srtm1";

/// Skip the test cleanly when the data fixture is absent.  Returning
/// `None` from `prepare_geocoded` triggers an `eprintln!` + early `return`
/// in the caller; this is the simplest skip mechanism that still surfaces
/// the reason in `--nocapture` output.
fn fixture_present() -> bool {
    Path::new(SAFE).is_dir()
        && Path::new(EOF_FILE).is_file()
        && Path::new(DEM_DIR).is_dir()
        && TIFFS.iter().all(|t| Path::new(t).is_file())
}

/// Run the full pipeline up to terrain correction and return the linear
/// γ⁰ geocoded image.  This is the same sequence the `dump_s1b_tc`
/// example uses, minus the dB conversion and TIFF write at the end.
fn build_geocoded() -> sardine_scene::terrain_correction::GeocodedImage {
    let scene = parse_safe_directory(Path::new(SAFE)).expect("parse SAFE");
    let eof = parse_eof_file(Path::new(EOF_FILE)).expect("parse EOF");
    let scene = apply_precise_orbit(scene, &eof).expect("apply POEORB");

    let cal_data = parse_calibration_noise(Path::new(SAFE)).expect("parse cal/noise");
    let geo_grids = parse_geolocation_grids(Path::new(SAFE)).expect("parse geo grids");

    let mut sigma0_arrays = Vec::new();
    let mut swath_metas: Vec<SubSwathMetadata> = Vec::new();
    let mut az_start_times: Vec<chrono::DateTime<chrono::Utc>> = Vec::new();

    for (iw_idx, &iw_id) in IW_IDS.iter().enumerate() {
        let sw = scene
            .sub_swaths
            .iter()
            .find(|s| s.id == iw_id)
            .expect("subswath in scene");

        let mut bursts: Vec<_> = scene
            .bursts
            .iter()
            .filter(|b| b.subswath_id == iw_id)
            .cloned()
            .collect();
        bursts.sort_by_key(|b| b.burst_index);

        let cal = cal_data
            .calibrations
            .iter()
            .find(|c| c.subswath_id == iw_id && c.polarization == Polarization::VV)
            .expect("VV cal");
        let noise = cal_data
            .noises
            .iter()
            .find(|n| n.subswath_id == iw_id && n.polarization == Polarization::VV)
            .expect("VV noise");

        let mut reader = SlcReader::open(TIFFS[iw_idx]).expect("open SLC");
        let deburst = deburst_subswath(&mut reader, sw, &bursts).expect("deburst");
        let sigma0 = apply_calibration(&deburst, cal, noise, 0).expect("calibrate");

        sigma0_arrays.push(sigma0);
        swath_metas.push(sw.clone());
        az_start_times.push(bursts[0].azimuth_time_utc);
    }

    let inputs: Vec<SwathInput<'_>> = sigma0_arrays
        .iter()
        .zip(swath_metas.iter())
        .zip(az_start_times.iter())
        .map(|((s, sw), &t0)| SwathInput {
            sigma0: s,
            swath: sw,
            azimuth_start_time: t0,
        })
        .collect();
    let merged = merge_subswaths(&inputs).expect("merge");
    drop(inputs);
    drop(sigma0_arrays);

    let dem = DemMosaic::load_directory(Path::new(DEM_DIR)).expect("DEM");
    let geoid = GeoidModel::Egm96(fetch_egm96().expect("EGM96"));

    let mut tc_cfg = TerrainCorrectionConfig::new(geoid);
    tc_cfg.pixel_spacing_deg = 0.0001;
    tc_cfg.flatten = true;
    tc_cfg.noise_floor_margin_db = 3.0;

    terrain_correction(&merged, &scene, &dem, &geo_grids, &tc_cfg).expect("terrain correction")
}

/// Mean of the finite (non-NaN) entries of `data`.  Returns NaN if the
/// input has no finite samples — the caller asserts on a real S-1 scene
/// where that cannot happen.
fn finite_mean(data: &[f32]) -> f64 {
    let mut sum = 0.0_f64;
    let mut n = 0_usize;
    for &v in data {
        if v.is_finite() {
            sum += v as f64;
            n += 1;
        }
    }
    if n == 0 { f64::NAN } else { sum / n as f64 }
}

/// Population variance of finite samples in a row-major sub-rectangle.
/// Returns NaN if the patch contains no finite samples.
fn patch_variance(
    data: &[f32],
    cols: usize,
    row0: usize,
    col0: usize,
    h: usize,
    w: usize,
) -> f64 {
    let mut sum = 0.0_f64;
    let mut sq = 0.0_f64;
    let mut n = 0_usize;
    for r in row0..row0 + h {
        for c in col0..col0 + w {
            let v = data[r * cols + c];
            if v.is_finite() {
                sum += v as f64;
                sq += (v as f64) * (v as f64);
                n += 1;
            }
        }
    }
    if n == 0 {
        return f64::NAN;
    }
    let m = sum / n as f64;
    (sq / n as f64 - m * m).max(0.0)
}

/// Locate a `patch × patch` window inside the scene that has no NaN
/// samples and the lowest coefficient of variation among a coarse grid
/// of candidates.  This is our proxy for "homogeneous patch": the
/// filters' variance-reduction property only holds where the underlying
/// scene is locally constant; selecting low-CV windows excludes built-up
/// areas, water/land transitions, and DEM voids.
fn find_homogeneous_patch(
    data: &[f32],
    cols: usize,
    rows: usize,
    patch: usize,
) -> Option<(usize, usize)> {
    let stride = patch; // tile the scene
    let mut best: Option<(usize, usize, f64)> = None;
    for r0 in (0..rows.saturating_sub(patch)).step_by(stride) {
        for c0 in (0..cols.saturating_sub(patch)).step_by(stride) {
            // require a full patch of finite samples
            let mut all_finite = true;
            let mut sum = 0.0_f64;
            let mut sq = 0.0_f64;
            'rows: for r in r0..r0 + patch {
                for c in c0..c0 + patch {
                    let v = data[r * cols + c];
                    if !v.is_finite() {
                        all_finite = false;
                        break 'rows;
                    }
                    sum += v as f64;
                    sq += (v as f64) * (v as f64);
                }
            }
            if !all_finite {
                continue;
            }
            let n = (patch * patch) as f64;
            let m = sum / n;
            if !(m > 0.0) {
                continue;
            }
            let var = (sq / n - m * m).max(0.0);
            let cv = var.sqrt() / m;
            if best.map_or(true, |(_, _, b)| cv < b) {
                best = Some((r0, c0, cv));
            }
        }
    }
    best.map(|(r, c, _)| (r, c))
}

#[test]
#[ignore = "runs the full TC pipeline on the S1B fixture (~minutes); opt in with --ignored"]
fn speckle_filters_preserve_mean_and_reduce_variance_on_real_scene() {
    if !fixture_present() {
        eprintln!(
            "speckle_e2e: skipping — S1B fixture not present at {SAFE}.  \
             This test is expected to skip on machines without the data."
        );
        return;
    }

    eprintln!("speckle_e2e: building geocoded γ⁰ (this takes a while) …");
    let geo = build_geocoded();
    let cols = geo.cols;
    let rows = geo.rows;
    let unfiltered = geo.data;
    let finite_in = unfiltered.iter().filter(|v| v.is_finite()).count();
    assert!(
        finite_in > 100_000,
        "scene has only {finite_in} finite pixels — fixture seems wrong"
    );
    let mean_in = finite_mean(&unfiltered);
    assert!(mean_in.is_finite() && mean_in > 0.0, "mean_in = {mean_in}");

    // Pick a single homogeneous patch once; the same patch is used for
    // every filter so variance ratios are directly comparable.
    let patch = 32_usize;
    let (pr, pc) =
        find_homogeneous_patch(&unfiltered, cols, rows, patch).expect("at least one homogeneous patch");
    let var_in = patch_variance(&unfiltered, cols, pr, pc, patch, patch);
    eprintln!(
        "speckle_e2e: {} × {} scene, finite={}, linear mean={:.4e}, \
         homogeneous patch=({},{})+{}², var={:.4e}",
        rows, cols, finite_in, mean_in, pr, pc, patch, var_in
    );
    assert!(var_in > 0.0, "homogeneous patch had zero variance — pathological");

    // Pre-compute the NaN mask once so we can verify it is preserved
    // pixel-for-pixel by every filter.
    let nan_mask: Vec<bool> = unfiltered.iter().map(|v| v.is_nan()).collect();

    let filters = [
        ("boxcar", SpeckleFilter::Boxcar { window: 7 }),
        ("lee", SpeckleFilter::Lee { window: 7, enl: 1.0 }),
        ("frost", SpeckleFilter::Frost { window: 7, damping: 1.0 }),
        ("gamma_map", SpeckleFilter::GammaMap { window: 7, enl: 4.0 }),
        ("refined_lee", SpeckleFilter::RefinedLee { enl: 1.0 }),
    ];

    for (name, filter) in filters {
        let out = apply_speckle_filter(&unfiltered, cols, rows, filter)
            .unwrap_or_else(|e| panic!("filter {name} failed: {e}"));

        // (1) shape preserved.
        assert_eq!(out.len(), unfiltered.len(), "{name}: length changed");

        // (2) NaN positions preserved.  We allow filters to *introduce*
        // NaN at the edges only if the entire window was NaN; for our
        // scene that does not occur, but the assertion is one-directional
        // ("input NaN ⇒ output NaN") to remain robust to that edge case.
        for (i, &was_nan) in nan_mask.iter().enumerate() {
            if was_nan {
                assert!(
                    out[i].is_nan(),
                    "{name}: NaN at index {i} was filled in by the filter"
                );
            }
        }

        // (3) scene-wide linear mean preserved within 2 %.
        let mean_out = finite_mean(&out);
        let rel_err = ((mean_out - mean_in) / mean_in).abs();
        assert!(
            rel_err < 0.02,
            "{name}: linear mean shifted by {:.3}% \
             (input={mean_in:.4e}, output={mean_out:.4e})",
            100.0 * rel_err
        );

        // (4) homogeneous-patch variance reduced.
        let var_out = patch_variance(&out, cols, pr, pc, patch, patch);
        assert!(
            var_out < var_in,
            "{name}: variance NOT reduced over homogeneous patch \
             (input={var_in:.4e}, output={var_out:.4e})"
        );
        eprintln!(
            "speckle_e2e: {name:>9}  Δmean={:+.3}%  var_ratio={:.3}",
            100.0 * (mean_out - mean_in) / mean_in,
            var_out / var_in
        );
    }
}

/// Verify that Refined Lee 7×7 applied to a geocoded γ⁰ scene does not
/// create row-mean discontinuities that would indicate seam-straddling
/// artefacts.
///
/// After terrain correction the scene is in map geometry; any burst-seam
/// structure from the original TOPS acquisition should already be healed
/// by the TC resampling.  We apply the filter in map space and assert that
/// no row deviates by more than 3 dB from an 11-row sliding neighbourhood
/// mean — a threshold that would catch genuine seam-induced power jumps
/// (≥ 2× power, i.e. ≥ 3 dB) but is loose enough to tolerate normal
/// scene heterogeneity.
///
/// If a seam artefact were present it would appear as a horizontal band
/// of anomalously high or low pixel values, causing the row mean to spike
/// far outside the local neighbourhood trend.
#[test]
#[ignore = "runs the full TC pipeline on the S1B fixture (~minutes); opt in with --ignored"]
fn refined_lee_does_not_create_seam_artifacts() {
    if !fixture_present() {
        eprintln!(
            "refined_lee_seam_check: skipping — S1B fixture not present at {SAFE}"
        );
        return;
    }

    eprintln!("refined_lee_seam_check: building geocoded γ⁰ …");
    let geo = build_geocoded();
    let cols = geo.cols;
    let rows = geo.rows;
    let linear = geo.data;

    eprintln!(
        "refined_lee_seam_check: applying Refined Lee 7×7 to {} × {} scene …",
        rows, cols
    );
    let filtered =
        apply_speckle_filter(&linear, cols, rows, SpeckleFilter::RefinedLee { enl: 1.0 })
            .expect("Refined Lee 7×7 failed");

    // Per-row finite mean in linear power; convert to dB for the deviation
    // check.  Rows where every pixel is NaN (image edges after TC) are
    // excluded from the neighbourhood comparison.
    let row_means_db: Vec<Option<f64>> = (0..rows)
        .map(|r| {
            let m = finite_mean(&filtered[r * cols..(r + 1) * cols]);
            if m > 0.0 { Some(10.0 * m.log10()) } else { None }
        })
        .collect();

    let valid: Vec<(usize, f64)> = row_means_db
        .iter()
        .enumerate()
        .filter_map(|(i, v)| v.map(|db| (i, db)))
        .collect();

    assert!(
        valid.len() > 200,
        "too few valid rows ({}) — fixture may be wrong",
        valid.len()
    );

    // Sliding 11-row neighbourhood mean; check each interior row.
    let half = 5_usize;
    let mut max_dev_db = 0.0_f64;
    let mut max_dev_row = 0_usize;

    for i in half..valid.len().saturating_sub(half) {
        let (row_idx, row_db) = valid[i];
        let nbr_mean = valid[i - half..i + half + 1]
            .iter()
            .map(|(_, v)| *v)
            .sum::<f64>()
            / (2 * half + 1) as f64;
        let dev = (row_db - nbr_mean).abs();
        if dev > max_dev_db {
            max_dev_db = dev;
            max_dev_row = row_idx;
        }
    }

    eprintln!(
        "refined_lee_seam_check: max row-mean deviation = {max_dev_db:.3} dB (row {max_dev_row})"
    );
    assert!(
        max_dev_db < 3.0,
        "Refined Lee 7×7 produces a {max_dev_db:.2} dB row-mean jump at row {max_dev_row} \
         — likely a seam-straddling artefact"
    );
}
