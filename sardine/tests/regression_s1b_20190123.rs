//! Radiometric regression test: sardine vs. ASF RTC10-GAMMA reference product.
//!
//! **Scene**: S1B IW SLC 2019-01-23, VV polarization.
//! **Reference**: ASF HyP3 RTC10 GAMMA — σ⁰/γ⁰ linear power, EPSG:32632, 10 m.
//! **Pipeline config**: WGS84 output, 0.0001° spacing, EGM96 geoid, POEORB,
//!   Refined Lee 7×7 (ENL=1), no multilook.
//!
//! **Speckle filter**: Refined Lee 7×7 (ENL=1) — matches ASF RTC10-GAMMA's
//! Enhanced Lee 7×7 filter.  With matching filter, the ENL-mismatch dB offset
//! (~−1 dB unfiltered) should collapse toward 0 dB.
//!
//! # Pass criteria
//!
//! | Metric              | Threshold |
//! |---------------------|-----------|
//! | Linear mean bias    | ≤ ±1.5 dB |
//! | dB median bias      | ≤ ±0.5 dB |
//! | Joint valid pixels  | ≥ 100,000  |
//!
//! The ±0.5 dB dB-domain threshold accommodates DEM differences (SRTM-1 vs
//! GLO-30) and filter-implementation differences (sardine Refined Lee vs
//! ASF Enhanced Lee).  The linear mean threshold is a backstop for gross
//! calibration regressions.
//!
//! # Opt-in
//!
//! ```sh
//! cargo test --release --features geoid-fetch \
//!     --test regression_s1b_20190123 -- --ignored --nocapture
//! ```
//!
//! The test calls `eprintln!` for progress because tracing is not
//! initialised in the integration-test harness.
//!
//! The test skips cleanly (no failure) when the input fixtures are absent.

#![cfg(feature = "geoid-fetch")]

use std::io::BufReader;
use std::path::{Path, PathBuf};

use sardine_scene::run::ProcessOptions;

// ─── Fixture paths ────────────────────────────────────────────────────────────

const SAFE: &str = "/home/datacube/dev/SARdine/data/SLC/\
    S1B_IW_SLC__1SDV_20190123T053348_20190123T053415_014617_01B3D4_E833.SAFE";

/// POEORB file — co-located with the SAFE in the test fixture layout.
const EOF: &str = "/home/datacube/dev/SARdine/data/SLC/\
    S1B_IW_SLC__1SDV_20190123T053348_20190123T053415_014617_01B3D4_E833.SAFE/\
    orbit_cache/S1B_IW_SLC__1SDV_20190123T053348_20190123T053415_014617_01B3D4_E833_POEORB.EOF";

const DEM: &str = "/home/datacube/dev/SARdine/data/dem/srtm1";

const ASF_VV: &str = "/home/datacube/dev/SARdine/data/ASF/\
    S1B_IW_20190123T053348_DVP_RTC10_G_gpufem_6EC9/\
    S1B_IW_20190123T053348_DVP_RTC10_G_gpufem_6EC9_VV.tif";

/// ASF geotransform — origin and pixel size, read from the file header.
/// In GDAL/GeoTIFF convention:
///   x = origin_x + col * px
///   y = origin_y + row * py   (py < 0 for north-up)
#[derive(Clone, Copy, Debug)]
struct Gt {
    origin_x: f64,
    origin_y: f64,
    px: f64,
    py: f64,
}

impl Gt {
    fn pixel_to_xy(&self, col: f64, row: f64) -> (f64, f64) {
        (
            self.origin_x + col * self.px,
            self.origin_y + row * self.py,
        )
    }
}

// ─── Helpers ──────────────────────────────────────────────────────────────────

fn fixtures_present() -> bool {
    Path::new(SAFE).is_dir()
        && Path::new(EOF).is_file()
        && Path::new(DEM).is_dir()
        && Path::new(ASF_VV).is_file()
}

/// Read a single-band Float32 strip TIFF.  Returns (pixels, cols, rows, gt).
///
/// `TAG_MODEL_TIEPOINT` (33922) and `TAG_MODEL_PIXEL_SCALE` (33550) are
/// read to construct the geotransform.  The origin is the upper-left corner
/// of the upper-left pixel (PixelIsArea convention as written by sardine's
/// own export module).
fn read_f32_geotiff(path: &Path) -> (Vec<f32>, usize, usize, Gt) {
    use tiff::decoder::{Decoder, DecodingResult, Limits};
    use tiff::tags::Tag;
    const TAG_MODEL_PIXEL_SCALE: u16 = 33550;
    const TAG_MODEL_TIEPOINT: u16 = 33922;

    let file = std::fs::File::open(path)
        .unwrap_or_else(|e| panic!("cannot open {}: {e}", path.display()));
    // The sardine output can be > 1 GB; lift all decoder size limits.
    let mut dec = Decoder::new(BufReader::new(file))
        .unwrap_or_else(|e| panic!("TIFF decode init {}: {e}", path.display()))
        .with_limits(Limits::unlimited());

    let (w, h) = dec.dimensions()
        .unwrap_or_else(|e| panic!("dimensions {}: {e}", path.display()));
    let cols = w as usize;
    let rows = h as usize;

    let scale = dec.get_tag_f64_vec(Tag::Unknown(TAG_MODEL_PIXEL_SCALE))
        .unwrap_or_else(|e| panic!("ModelPixelScaleTag {}: {e}", path.display()));
    assert!(scale.len() >= 2, "ModelPixelScaleTag must have ≥ 2 values");

    let tp = dec.get_tag_f64_vec(Tag::Unknown(TAG_MODEL_TIEPOINT))
        .unwrap_or_else(|e| panic!("ModelTiepointTag {}: {e}", path.display()));
    assert!(tp.len() >= 6, "ModelTiepointTag must have ≥ 6 values");

    // Tiepoint: (raster_i, raster_j, 0, world_x, world_y, 0)
    // PixelIsArea: upper-left of pixel (raster_i, raster_j) is at (world_x, world_y).
    let origin_x = tp[3] - tp[0] * scale[0];
    let origin_y = tp[4] + tp[1] * scale[1]; // scale[1] is positive; rows go south

    let gt = Gt {
        origin_x,
        origin_y,
        px: scale[0],
        py: -scale[1],
    };

    let img = dec.read_image()
        .unwrap_or_else(|e| panic!("read_image {}: {e}", path.display()));
    let pixels = match img {
        DecodingResult::F32(v) => v,
        other => panic!("expected F32, got {other:?} in {}", path.display()),
    };
    assert_eq!(pixels.len(), rows * cols, "pixel count mismatch");

    (pixels, cols, rows, gt)
}

/// Median of a slice.  Returns NaN for empty input.
fn median(v: &mut [f32]) -> f32 {
    if v.is_empty() {
        return f32::NAN;
    }
    v.sort_unstable_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let mid = v.len() / 2;
    if v.len() % 2 == 0 {
        (v[mid - 1] + v[mid]) / 2.0
    } else {
        v[mid]
    }
}

// ─── Test ─────────────────────────────────────────────────────────────────────

#[test]
#[ignore]
fn s1b_20190123_vv_refined_lee_bias_vs_asf() {
    if !fixtures_present() {
        eprintln!(
            "regression_s1b_20190123: skipping — one or more data fixtures absent.\n\
             Required:\n  {SAFE}\n  {EOF}\n  {DEM}\n  {ASF_VV}"
        );
        return;
    }

    let tmp = tempdir();

    let output = tmp.join("sardine_s1b_vv.tif");

    // ── 1. Run the full pipeline ──────────────────────────────────────────────
    eprintln!("regression: running sardine pipeline …");
    let mut opts = ProcessOptions::new(
        PathBuf::from(SAFE),
        Some(PathBuf::from(DEM)),
        output.clone(),
        "auto".to_owned(), // EGM96 geoid via geoid-fetch
    );
    opts.orbit = Some(PathBuf::from(EOF));
    opts.polarization = "VV".to_owned();
    opts.no_provenance = true;
    opts.crs = "EPSG:4326".to_owned();
    opts.pixel_spacing_deg = 0.0001;
    // Apply Refined Lee 7×7 to match ASF RTC10-GAMMA's Enhanced Lee 7×7 filter.
    // With matching filter ENL ≈ 1 (single-look SLC input) the dB-domain median
    // bias should collapse from ~−1 dB (ENL mismatch) toward ~0 dB.
    opts.speckle = "refined-lee".to_owned();
    opts.enl = 1.0;

    sardine_scene::run::run_process(&opts)
        .unwrap_or_else(|e| panic!("run_process failed: {e:#}"));

    assert!(output.is_file(), "sardine output TIFF must exist after run_process");

    // ── 2. Read sardine output ────────────────────────────────────────────────
    eprintln!("regression: reading sardine output …");
    let (sardine_px, sar_cols, sar_rows, sar_gt) = read_f32_geotiff(&output);
    eprintln!(
        "  sardine: {}×{}, origin=({:.6},{:.6}), spacing=({:.6},{:.6})",
        sar_cols, sar_rows,
        sar_gt.origin_x, sar_gt.origin_y,
        sar_gt.px, sar_gt.py
    );

    // ── 3. Read ASF reference ─────────────────────────────────────────────────
    eprintln!("regression: reading ASF reference …");
    let (asf_px, asf_cols, asf_rows, asf_gt) = read_f32_geotiff(Path::new(ASF_VV));
    eprintln!(
        "  ASF: {}×{}, origin=({:.1},{:.1}), spacing=({:.1},{:.1})",
        asf_cols, asf_rows,
        asf_gt.origin_x, asf_gt.origin_y,
        asf_gt.px, asf_gt.py
    );

    // ── 4. Build coordinate projector UTM32N → WGS84 ─────────────────────────
    // The ASF raster is in EPSG:32632 (UTM zone 32N, WGS84).
    // We need to convert each ASF pixel centre to WGS84 lon/lat, then look
    // up the nearest sardine pixel.
    use proj4rs::Proj;
    use proj4rs::transform::transform;
    let utm32n = Proj::from_proj_string(
        "+proj=utm +zone=32 +datum=WGS84 +units=m +no_defs"
    )
    .expect("UTM32N proj");
    let wgs84 = Proj::from_proj_string(
        "+proj=longlat +datum=WGS84 +no_defs"
    )
    .expect("WGS84 proj");

    // ── 5. Sample every Nth pixel to keep runtime tractable ──────────────────
    // ASF is ~22388 × 28575.  Step = 25 → ~1 M candidate points.
    const STEP: usize = 25;

    // Accumulate linear-domain sums for a filter-invariant calibration check.
    // We compare mean linear power values over the joint valid area.
    //
    // Why linear and not dB-domain?
    //   Sardine outputs unfiltered (ENL ≈ 1) data.  ASF RTC10 applies a 7×7
    //   Enhanced Lee filter (ENL ≈ 30–50).  For single-look speckle:
    //     E[10·log₁₀(x)] ≈ 10·log₁₀(E[x]) − 2.507 dB   (Jensen's inequality)
    //   This causes a systematic ~−1 dB offset in the dB-domain median that
    //   disappears in the linear-domain mean.  See docs/PROGRESS.md §
    //   "filter_bias_check" for the full investigation.
    //   Linear-domain mean bias is constant at ~+0.466 dB across all filter
    //   sizes; dB-domain RMSE between two independent speckle realisations
    //   (sardine SLC vs ASF GRD) is ~5–6 dB and is meaningless as a metric.
    let mut sum_sardine_linear: f64 = 0.0;
    let mut sum_asf_linear: f64 = 0.0;
    let mut diffs_db: Vec<f32> = Vec::with_capacity(
        (asf_rows / STEP + 1) * (asf_cols / STEP + 1)
    );

    let mut n_asf_nodata = 0usize;
    let mut n_outside = 0usize;
    let mut n_sar_nodata = 0usize;

    for row in (0..asf_rows).step_by(STEP) {
        for col in (0..asf_cols).step_by(STEP) {
            let asf_linear = asf_px[row * asf_cols + col];

            // Skip ASF nodata (nodata = 0.0, negative linear values are impossible).
            if asf_linear <= 0.0 {
                n_asf_nodata += 1;
                continue;
            }

            // Pixel centre in UTM32N (PixelIsArea → add 0.5 pixel to corner).
            let (utm_x, utm_y) = asf_gt.pixel_to_xy(col as f64 + 0.5, row as f64 + 0.5);

            // Convert UTM32N → WGS84 via proj4rs (radians out).
            let mut p = (utm_x, utm_y, 0.0);
            if transform(&utm32n, &wgs84, &mut p).is_err() {
                n_outside += 1;
                continue;
            }
            let lon_deg = p.0.to_degrees();
            let lat_deg = p.1.to_degrees();

            // Map WGS84 → sardine pixel (nearest neighbour).
            // sar_gt.py is negative (north-up).
            let sar_col_f = (lon_deg - sar_gt.origin_x) / sar_gt.px;
            let sar_row_f = (lat_deg - sar_gt.origin_y) / sar_gt.py;

            let sar_col = sar_col_f.round() as isize;
            let sar_row = sar_row_f.round() as isize;

            if sar_col < 0
                || sar_row < 0
                || sar_col >= sar_cols as isize
                || sar_row >= sar_rows as isize
            {
                n_outside += 1;
                continue;
            }

            let sar_db = sardine_px[sar_row as usize * sar_cols + sar_col as usize];
            if !sar_db.is_finite() {
                n_sar_nodata += 1;
                continue;
            }

            // Accumulate linear values for the calibration check.
            let sardine_linear = 10.0_f32.powf(sar_db / 10.0);
            sum_sardine_linear += sardine_linear as f64;
            sum_asf_linear += asf_linear as f64;

            // Also collect dB diffs for informational reporting.
            let asf_db = 10.0 * asf_linear.log10();
            diffs_db.push(sar_db - asf_db);
        }
    }

    let n_pairs = diffs_db.len();
    eprintln!(
        "  sample stats: n_pairs={n_pairs}, asf_nodata={n_asf_nodata}, \
         outside={n_outside}, sar_nodata={n_sar_nodata}"
    );

    // ── 6. Assertions ─────────────────────────────────────────────────────────
    assert!(
        n_pairs >= 100_000,
        "too few joint-valid comparison pixels ({n_pairs}) — check extents/CRS"
    );

    // Linear-domain mean bias (filter-invariant baseline check).
    let mean_sardine_linear = (sum_sardine_linear / n_pairs as f64) as f32;
    let mean_asf_linear = (sum_asf_linear / n_pairs as f64) as f32;
    let linear_bias_db = 10.0 * (mean_sardine_linear / mean_asf_linear).log10();

    // dB-domain median bias — with matching Refined Lee 7×7 filter applied,
    // the ENL-mismatch offset (~−1 dB) should be largely removed.  We assert
    // a looser ±0.5 dB threshold to accommodate DEM differences (SRTM-1 vs
    // GLO-30) and filter-implementation differences (sardine Refined Lee vs
    // ASF Enhanced Lee).
    let median_bias_db = median(&mut diffs_db);
    eprintln!(
        "  linear mean bias: {linear_bias_db:+.4} dB  (threshold ±1.5 dB)"
    );
    eprintln!(
        "  dB median bias:   {median_bias_db:+.4} dB  (threshold ±0.5 dB)"
    );

    assert!(
        linear_bias_db.abs() <= 1.5,
        "linear mean bias {linear_bias_db:+.4} dB exceeds ±1.5 dB threshold"
    );
    assert!(
        median_bias_db.abs() <= 0.5,
        "dB median bias {median_bias_db:+.4} dB exceeds ±0.5 dB threshold \
         (Refined Lee 7×7 should remove most ENL mismatch vs ASF Enhanced Lee 7×7)"
    );
}

// ─── Tiny temp-dir helper (no extra dep) ──────────────────────────────────────

fn tempdir() -> PathBuf {
    use std::time::{SystemTime, UNIX_EPOCH};
    let ts = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map(|d| d.as_nanos())
        .unwrap_or(0); // SAFETY-OK: fallback to 0 only on pre-epoch clock — path still unique enough for tests
    let dir = std::env::temp_dir().join(format!("sardine_regression_{ts}"));
    std::fs::create_dir_all(&dir).expect("create temp dir");
    dir
}
