//! InSAR smoke test: single-scene coherence (reference == secondary).
//!
//! Using the same S1B SAFE as both reference and secondary gives a
//! coherence of 1.0 everywhere (zero temporal baseline), which is a
//! tight end-to-end sanity check:
//!
//! * The full chain (parse → co-reg → coherence estimation → terrain
//!   correction → GeoTIFF write) must complete without panic.
//! * Every finite coherence pixel must be in **[0.0, 1.0]**.
//! * Three geocoded coherence GeoTIFFs (`_iw1_coherence.tif` …
//!   `_iw3_coherence.tif`) must exist after the run.
//!
//! # Opt-in
//!
//! ```sh
//! cargo test --release --features geoid-fetch \
//!     --test insar_smoke -- --ignored --nocapture
//! ```
//!
//! The test skips cleanly (no failure) when fixtures are absent.

#![cfg(feature = "geoid-fetch")]

use std::path::Path;

use sardine::run::{run_insar, InsarOptions};

// ─── Fixture paths ────────────────────────────────────────────────────────────

const SAFE: &str = "/home/datacube/dev/SARdine/data/SLC/\
    S1B_IW_SLC__1SDV_20190123T053348_20190123T053415_014617_01B3D4_E833.SAFE";

const EOF: &str = "/home/datacube/dev/SARdine/data/SLC/\
    S1B_IW_SLC__1SDV_20190123T053348_20190123T053415_014617_01B3D4_E833.SAFE/\
    orbit_cache/S1B_IW_SLC__1SDV_20190123T053348_20190123T053415_014617_01B3D4_E833_POEORB.EOF";

const DEM: &str = "/home/datacube/dev/SARdine/data/dem/srtm1";

// ─── Helpers ─────────────────────────────────────────────────────────────────

fn read_coherence_values(path: &str) -> Vec<f32> {
    use std::io::BufReader;
    let f = std::fs::File::open(path).expect("open coherence tif");
    let mut decoder = tiff::decoder::Decoder::new(BufReader::new(f)).expect("tiff decoder");
    match decoder.read_image().expect("read tiff") {
        tiff::decoder::DecodingResult::F32(v) => v,
        other => panic!("expected F32 tiff, got {other:?}"),
    }
}

// ─── Test ─────────────────────────────────────────────────────────────────────

#[test]
#[ignore]
fn insar_smoke_same_scene_coherence_is_one() {
    // Skip gracefully if fixtures are absent.
    for path in [SAFE, EOF, DEM] {
        if !Path::new(path).exists() {
            eprintln!("insar_smoke: skipping — fixture not found: {path}");
            return;
        }
    }

    let out_dir = std::env::temp_dir().join("sardine_insar_smoke");
    std::fs::create_dir_all(&out_dir).unwrap();
    let out_base = out_dir.join("smoke");

    let opts = InsarOptions {
        reference: SAFE.into(),
        secondary: SAFE.into(),
        output: out_base.clone(),
        dem: Some(DEM.into()),
        dem_source: "srtm1".to_owned(),
        reference_orbit: Some(EOF.into()),
        secondary_orbit: Some(EOF.into()),
        polarization: "VV".to_owned(),
        az_looks: sardine::insar::interferogram::DEFAULT_COH_AZ_LOOKS,
        rg_looks: sardine::insar::interferogram::DEFAULT_COH_RG_LOOKS,
        output_phase: false,
        geoid: "auto".to_owned(),
        pixel_spacing_deg: 0.0001,
        pixel_spacing_m: 10.0,
        crs: "wgs84".to_owned(),
        cog: false,
        threads: 0,
    };

    eprintln!("insar_smoke: running pipeline (this may take several minutes) …");
    run_insar(&opts).expect("run_insar failed");

    for iw in ["iw1", "iw2", "iw3"] {
        let coh_path = format!(
            "{}_{}_{}_coherence.tif",
            out_base.display(),
            opts.polarization.to_lowercase(),
            iw
        );
        // Accept either naming convention depending on build variant.
        let alt_path = format!("{}_{}_{}", out_base.display(), iw, "coherence.tif");

        let found_path = [&coh_path, &alt_path]
            .iter()
            .find(|p| Path::new(*p).exists())
            .copied()
            .unwrap_or_else(|| {
                // List what was actually written to help debug.
                let files: Vec<_> = std::fs::read_dir(&out_dir)
                    .unwrap()
                    .filter_map(|e| e.ok())
                    .map(|e| e.file_name().to_string_lossy().into_owned())
                    .collect();
                panic!(
                    "coherence file not found for {iw}. Written files: {files:?}"
                );
            });

        eprintln!("insar_smoke: checking {found_path}");
        let values = read_coherence_values(found_path);
        let finite: Vec<f32> = values.into_iter().filter(|v| v.is_finite()).collect();
        assert!(
            !finite.is_empty(),
            "{iw} coherence has no finite pixels"
        );
        for &v in &finite {
            assert!(
                (0.0..=1.0).contains(&v),
                "{iw} coherence pixel out of [0,1]: {v}"
            );
        }
        // Same-scene pair: median coherence should be close to 1.0.
        let mut sorted = finite.clone();
        sorted.sort_by(|a, b| a.partial_cmp(b).unwrap());
        let median = sorted[sorted.len() / 2];
        eprintln!(
            "insar_smoke: {iw} median coherence = {median:.4} (n_finite={})",
            sorted.len()
        );
        assert!(
            median > 0.8,
            "{iw} median coherence {median:.4} < 0.8 for same-scene pair"
        );
    }

    eprintln!("insar_smoke: PASSED");
}
