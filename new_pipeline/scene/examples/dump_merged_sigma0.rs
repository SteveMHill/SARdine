//! Run the full pipeline (SLC read → deburst → calibrate → merge → ground range
//! → terrain correction) for IW VV and write outputs as uncompressed float32 TIFFs.
//!
//! Usage:
//!   cargo run --example dump_merged_sigma0 --release [-- <output_path>]
//!
//! Default output path: /tmp/sardine_merged_sigma0.tiff (slant range)
//!                      /tmp/sardine_grd_sigma0.tiff    (ground range, 10 m)
//!                      /tmp/sardine_tc_sigma0.tiff     (terrain-corrected)
//!
//! STRONGLY recommended to run with --release; debug mode is ~30× slower.

use sardine_scene::apply_calibration::apply_calibration;
use sardine_scene::calibration::parse_calibration_noise;
use sardine_scene::dem::DemMosaic;
use sardine_scene::deburst::deburst_subswath;
use sardine_scene::export::{to_db_inplace, write_geotiff};
use sardine_scene::ground_range::to_ground_range;
use sardine_scene::merge_subswaths::{merge_subswaths, SwathInput};
use sardine_scene::orbit::{apply_precise_orbit, parse_eof_file};
use sardine_scene::parse::{parse_geolocation_grids, parse_safe_directory};
use sardine_scene::slc_reader::SlcReader;
use sardine_scene::terrain_correction::{terrain_correction, TerrainCorrectionConfig};
use sardine_scene::types::{Polarization, SubSwathId, SubSwathMetadata};

const SAFE: &str = "/home/datacube/dev/SARdine/data/SLC/\
    S1A_IW_SLC__1SDV_20201005T170824_20201005T170851_034664_04098A_1E66.SAFE";

const TIFFS: [&str; 3] = [
    "/home/datacube/dev/SARdine/data/SLC/\
     S1A_IW_SLC__1SDV_20201005T170824_20201005T170851_034664_04098A_1E66.SAFE/\
     measurement/s1a-iw1-slc-vv-20201005t170824-20201005t170849-034664-04098a-004.tiff",
    "/home/datacube/dev/SARdine/data/SLC/\
     S1A_IW_SLC__1SDV_20201005T170824_20201005T170851_034664_04098A_1E66.SAFE/\
     measurement/s1a-iw2-slc-vv-20201005t170824-20201005t170850-034664-04098a-005.tiff",
    "/home/datacube/dev/SARdine/data/SLC/\
     S1A_IW_SLC__1SDV_20201005T170824_20201005T170851_034664_04098A_1E66.SAFE/\
     measurement/s1a-iw3-slc-vv-20201005t170825-20201005t170851-034664-04098a-006.tiff",
];

const IW_IDS: [SubSwathId; 3] = [SubSwathId::IW1, SubSwathId::IW2, SubSwathId::IW3];

const EOF_FILE: &str = "/home/datacube/dev/SARdine/data/SLC/\
    S1A_IW_SLC__1SDV_20201005T170824_20201005T170851_034664_04098A_1E66.SAFE/orbit_cache/\
    S1A_IW_SLC__1SDV_20201005T170824_20201005T170851_034664_04098A_1E66_POEORB.EOF";

const DEM_DIR: &str = "/home/datacube/dev/SARdine/data/dem/srtm1";

#[cfg(not(feature = "geoid-fetch"))]
fn main() -> Result<(), Box<dyn std::error::Error>> {
    Err("this example requires the `geoid-fetch` Cargo feature. \
         Re-run with: cargo run --release --features geoid-fetch \
         --example dump_merged_sigma0"
        .into())
}

#[cfg(feature = "geoid-fetch")]
fn main() -> Result<(), Box<dyn std::error::Error>> {
    // ── Output path ────────────────────────────────────────────────────────────
    let args: Vec<String> = std::env::args().collect();
    let out_path = args.get(1).map(|s| s.as_str())
        .unwrap_or("/tmp/sardine_merged_sigma0.tiff");

    // ── Step 1: parse metadata ─────────────────────────────────────────────────
    eprintln!("[1/7] Parsing SAFE metadata + precise orbit ...");
    let scene = parse_safe_directory(std::path::Path::new(SAFE))?;
    let scene = {
        let eof_path = std::path::Path::new(EOF_FILE);
        if eof_path.is_file() {
            let eof = parse_eof_file(eof_path)?;
            eprintln!("    using precise orbit: {}", eof_path.display());
            apply_precise_orbit(scene, &eof)?
        } else {
            // Silent fallback to annotation orbit is scientifically unsafe: the
            // annotation orbit is decimetre-to-metre accurate while POEORB is
            // centimetre accurate, and downstream geolocation differs by the
            // same margin.  Require the user to opt in explicitly.
            let allow_annotation = std::env::var("SARDINE_ALLOW_ANNOTATION_ORBIT")
                .map(|v| v == "1" || v.eq_ignore_ascii_case("true"))
                .unwrap_or(false);
            if !allow_annotation {
                return Err(format!(
                    "precise orbit EOF not found at {} and \
                     SARDINE_ALLOW_ANNOTATION_ORBIT is not set.\n\
                     Either (a) provide the POEORB file, or (b) re-run with \
                     SARDINE_ALLOW_ANNOTATION_ORBIT=1 to explicitly accept the \
                     reduced-accuracy annotation orbit.",
                    eof_path.display()
                )
                .into());
            }
            eprintln!(
                "    WARNING: precise EOF not found at {}",
                eof_path.display()
            );
            eprintln!("    SARDINE_ALLOW_ANNOTATION_ORBIT=1 set → using annotation orbit");
            eprintln!("    → geolocation accuracy reduced to ~metre level");
            scene
        }
    };
    let cal_data = parse_calibration_noise(std::path::Path::new(SAFE))?;
    let geo_grids = parse_geolocation_grids(std::path::Path::new(SAFE))?;

    // ── Steps 2–3: per-subswath deburst + calibrate ───────────────────────────
    let mut sigma0_arrays = Vec::new();
    let mut swath_metas: Vec<SubSwathMetadata> = Vec::new();
    let mut az_start_times: Vec<chrono::DateTime<chrono::Utc>> = Vec::new();

    for (iw_idx, &iw_id) in IW_IDS.iter().enumerate() {
        let sw = scene.sub_swaths.iter()
            .find(|s| s.id == iw_id)
            .ok_or_else(|| format!("{iw_id} not found in scene metadata"))?;

        let mut bursts: Vec<_> = scene.bursts.iter()
            .filter(|b| b.subswath_id == iw_id)
            .cloned()
            .collect();
        bursts.sort_by_key(|b| b.burst_index);

        let cal = cal_data.calibrations.iter()
            .find(|c| c.subswath_id == iw_id && c.polarization == Polarization::VV)
            .ok_or_else(|| format!("no VV calibration for {iw_id}"))?;
        let noise = cal_data.noises.iter()
            .find(|n| n.subswath_id == iw_id && n.polarization == Polarization::VV)
            .ok_or_else(|| format!("no VV noise for {iw_id}"))?;

        eprintln!("[2/7] {} debursting {} bursts × {} lines/burst ...",
            iw_id, bursts.len(), sw.lines_per_burst);
        let mut reader = SlcReader::open(TIFFS[iw_idx])?;
        let deburst = deburst_subswath(&mut reader, sw, &bursts)?;
        drop(reader);

        eprintln!("[3/7] {} calibrating ({} lines × {} samples) ...",
            iw_id, deburst.lines, deburst.samples);
        // tiff_line_origin = 0: each subswath has its own TIFF file;
        // burst 0 starts at TIFF row 0.
        let sigma0 = apply_calibration(&deburst, cal, noise, 0)?;
        drop(deburst);  // free ~1–2 GB of raw CInt16 data

        sigma0_arrays.push(sigma0);
        swath_metas.push(sw.clone());
        // Debursted line 0 corresponds to bursts[0] line 0, whose
        // zero-Doppler UTC time is the annotation `azimuthTime` of burst 0.
        az_start_times.push(bursts[0].azimuth_time_utc);
    }

    // ── Step 4: merge ─────────────────────────────────────────────────────────
    eprintln!("[4/7] Merging IW1 + IW2 + IW3 ...");
    let inputs: Vec<SwathInput<'_>> = sigma0_arrays.iter()
        .zip(swath_metas.iter())
        .zip(az_start_times.iter())
        .map(|((s, sw), &t0)| SwathInput {
            sigma0: s,
            swath: sw,
            azimuth_start_time: t0,
        })
        .collect();
    let merged = merge_subswaths(&inputs)?;

    // Free per-subswath arrays before the write (they are now merged into `merged`).
    drop(inputs);
    drop(sigma0_arrays);
    drop(swath_metas);

    // ── Stats ──────────────────────────────────────────────────────────────────
    let total_px = merged.data.len();
    let mut n_valid = 0usize;
    let mut min_val = f32::INFINITY;
    let mut max_val = f32::NEG_INFINITY;
    let mut sum_val = 0.0f64;
    for &v in &merged.data {
        if !v.is_nan() {
            n_valid += 1;
            if v < min_val { min_val = v; }
            if v > max_val { max_val = v; }
            sum_val += v as f64;
        }
    }
    let mean_val = sum_val / n_valid as f64;
    let pct_valid = 100.0 * n_valid as f64 / total_px as f64;

    eprintln!("    {} lines × {} samples  ({:.1}% valid)",
        merged.lines, merged.samples, pct_valid);
    eprintln!("    σ⁰ linear: min={:.4e}  max={:.4e}  mean={:.4e}",
        min_val, max_val, mean_val);
    eprintln!("    near_range_time={:.6} ms  rps={:.4} m",
        merged.near_slant_range_time_s * 1e3, merged.range_pixel_spacing_m);

    // ── Step 5: write slant-range TIFF ─────────────────────────────────────────
    let mb = (merged.data.len() * 4) as f64 / (1 << 20) as f64;
    eprintln!("[5/7] Writing slant-range TIFF {out_path} ({mb:.0} MiB) ...");
    write_f32_tiff(out_path, &merged.data, merged.samples, merged.lines)?;

    // ── Step 6: ground range projection + multilook ────────────────────────────
    let target_spacing_m = 10.0;
    // Azimuth pixel spacing must come from the product metadata (not hardcoded).
    // All IW subswaths share this spacing to within < 0.1 %; use IW1.
    let az_spacing_m = scene
        .sub_swaths
        .first()
        .ok_or("scene has no subswaths")?
        .azimuth_pixel_spacing_m;
    eprintln!(
        "[6/7] Ground range projection → {target_spacing_m} m (input az spacing = {az_spacing_m:.3} m) ..."
    );
    let grd = to_ground_range(&merged, &geo_grids, target_spacing_m, az_spacing_m)?;

    // GRD stats
    let grd_total = grd.data.len();
    let mut grd_valid = 0usize;
    let mut grd_sum = 0.0f64;
    for &v in &grd.data {
        if !v.is_nan() {
            grd_valid += 1;
            grd_sum += v as f64;
        }
    }
    let grd_mean = grd_sum / grd_valid.max(1) as f64;
    let grd_pct = 100.0 * grd_valid as f64 / grd_total as f64;

    eprintln!("    GRD: {} lines × {} samples  ({:.1}% valid)",
        grd.lines, grd.samples, grd_pct);
    eprintln!("    range spacing: {:.2} m, azimuth spacing: {:.2} m, azimuth_looks: {}",
        grd.range_pixel_spacing_m, grd.azimuth_pixel_spacing_m, grd.azimuth_looks);
    eprintln!("    σ⁰ mean={:.4e}", grd_mean);

    let grd_path = out_path.replace("merged_sigma0", "grd_sigma0");
    let grd_path = if grd_path == out_path {
        "/tmp/sardine_grd_sigma0.tiff".to_string()
    } else {
        grd_path
    };
    let grd_mb = (grd.data.len() * 4) as f64 / (1 << 20) as f64;
    eprintln!("    Writing GRD TIFF {grd_path} ({grd_mb:.0} MiB) ...");
    write_f32_tiff(&grd_path, &grd.data, grd.samples, grd.lines)?;

    // ── Step 7: terrain correction (Range-Doppler geocoding) ─────────────────
    eprintln!("[7/7] Terrain correction (Range-Doppler geocoding) ...");
    let dem = DemMosaic::load_directory(std::path::Path::new(DEM_DIR))?;
    eprintln!("    Loaded {} DEM tiles from {}", dem.tile_count(), DEM_DIR);

    // Load EGM96 geoid model.  This example is gated on `geoid-fetch` (see
    // the no-feature stub below) so we never silently fall back to
    // `GeoidModel::Zero`, which would inject a ~40 m geolocation bias on
    // land without any signal in the example's output.  See AGENTS.md
    // "no silent fallbacks".
    use sardine_scene::geoid::GeoidModel;
    use sardine_scene::geoid_fetch::fetch_egm96;
    let geoid = {
        let grid = fetch_egm96().map_err(|e| -> Box<dyn std::error::Error> {
            format!(
                "EGM96 grid fetch failed: {e}. Set $SARDINE_GEOID_DIR or \
                 check your network and rerun. This example refuses to fall \
                 back to GeoidModel::Zero (would inject ~40 m vertical bias)."
            )
            .into()
        })?;
        eprintln!("    EGM96 geoid loaded — applying orthometric→ellipsoidal correction");
        GeoidModel::Egm96(grid)
    };

    let mut tc_cfg = TerrainCorrectionConfig::new(geoid);
    tc_cfg.pixel_spacing_deg = 0.0001;
    tc_cfg.flatten = true; // Small (2011) area-projection: σ⁰ → γ⁰

    let mut tc = terrain_correction(
        &merged,
        &scene,
        &dem,
        &geo_grids,
        &tc_cfg,
    )?;

    let tc_path = out_path.replace("merged_sigma0", "tc_sigma0");
    let tc_path = if tc_path == out_path {
        "/tmp/sardine_tc_sigma0.tiff".to_string()
    } else {
        tc_path
    };
    let tc_mb = (tc.data.len() * 4) as f64 / (1 << 20) as f64;
    eprintln!(
        "    TC: {} rows × {} cols, valid {:.1}%",
        tc.rows,
        tc.cols,
        100.0 * tc.valid_fraction()
    );
    eprintln!(
        "        buckets: valid={} dem_missing={} outside_footprint={} non_converged={} degenerate={} flat_masked={}",
        tc.valid_pixel_count,
        tc.dem_missing_count,
        tc.outside_footprint_count,
        tc.non_converged_count,
        tc.degenerate_geometry_count,
        tc.flat_masked_count,
    );
    eprintln!("    Writing TC TIFF {tc_path} ({tc_mb:.0} MiB) ...");
    write_f32_tiff(&tc_path, &tc.data, tc.cols, tc.rows)?;

    // ── Step 7b: Mask invalid pixels and convert γ⁰ to dB ─────────────────────
    //
    // After terrain flattening, any pixel with value ≤ 0 is a noise-floor
    // artefact (the calibration formula clamps (|DN|²−N) to zero; values at
    // exactly zero represent dark pixels where noise dominated the signal).
    // We mask them to NaN and convert the remainder to dB: 10·log₁₀(γ⁰).
    //
    // noise_floor_linear = 0.0 means mask only the physically impossible values
    // (≤ 0).  Pass a higher threshold (e.g. 10f32.powf(-25.0/10.0) ≈ 0.00316
    // for -25 dB NESZ) if a calibrated noise estimate is available.
    eprintln!("    Applying noise-floor mask and converting to dB ...");
    let (db_valid, db_masked) = to_db_inplace(&mut tc.data, 0.0)
        .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidData, e.to_string()))?;
    eprintln!("        dB conversion: valid={db_valid} noise_masked={db_masked}");

    // ── Step 7c: Write self-contained GeoTIFF with embedded EPSG:4326 CRS ────
    //
    // Uses GeoTIFF extension tags (ModelPixelScale, ModelTiepoint,
    // GeoKeyDirectory) — no .tfw or .prj sidecar files needed.
    let db_path = if tc_path.ends_with(".tiff") {
        tc_path.replace(".tiff", "_db.tiff")
    } else {
        tc_path.replace(".tif", "_db.tif")
    };
    eprintln!("    Writing γ⁰ dB GeoTIFF {db_path} ({tc_mb:.0} MiB, embedded EPSG:4326) ...");
    write_geotiff(&db_path, &tc.data, tc.cols, tc.rows, tc.geotransform)
        .map_err(|e| std::io::Error::new(std::io::ErrorKind::Other, e.to_string()))?;

    eprintln!("Done.");
    eprintln!("Inspect with:");
    eprintln!("  gdalinfo {out_path}");
    eprintln!("  gdalinfo {grd_path}");
    eprintln!("  gdalinfo {db_path}");

    Ok(())
}

// ─── Minimal uncompressed float32 TIFF writer ─────────────────────────────────
//
// Writes a grayscale (single-band) IEEE-754 float32 TIFF compatible with GDAL,
// SNAP, QGIS, ImageJ, etc.  No geo-referencing; pure pixel data with correct
// SampleFormat=3 (IEEE floating point) tag.
//
// File layout (little-endian, classic TIFF):
//   offset 0  : 8-byte TIFF header
//   offset 8  : IFD (162 bytes: 2-byte count + 13×12-byte entries + 4-byte next-IFD)
//   offset 170 : strip offsets array   (n_strips × 4 bytes)
//   offset 170+n_strips×4 : strip byte counts (n_strips × 4 bytes)
//   offset 170+n_strips×8 : XResolution RATIONAL (8 bytes)
//   offset 178+n_strips×8 : YResolution RATIONAL (8 bytes)
//   offset 186+n_strips×8 : image data (height × width × 4 bytes, row-major)
//
// Classic TIFF limit: all strip offsets must fit in u32 (max 4 GB file).
// A full S1 IW scene (~12 k × 67 k × 4 B ≈ 3.2 GiB) stays within this limit.
fn write_f32_tiff(
    path: &str,
    data: &[f32],
    width: usize,
    height: usize,
) -> std::io::Result<()> {
    // Sanity: this function assumes a little-endian host (x86-64).
    // TIFF stores as "II" (little-endian); byte order of f32 must match.
    #[cfg(not(target_endian = "little"))]
    compile_error!("write_f32_tiff assumes a little-endian host");

    use std::io::{BufWriter, Write};

    const ROWS_PER_STRIP: u32 = 64;
    const N_IFD_ENTRIES: u16 = 13;

    let n_strips = (height as u32 + ROWS_PER_STRIP - 1) / ROWS_PER_STRIP;
    let bytes_per_row = (width * 4) as u32;
    let full_strip_bytes = ROWS_PER_STRIP * bytes_per_row;
    let last_strip_rows = height as u32 - (n_strips - 1) * ROWS_PER_STRIP;
    let last_strip_bytes = last_strip_rows * bytes_per_row;

    // ── Compute layout offsets ─────────────────────────────────────────────────
    const IFD_OFFSET: u32 = 8;
    // IFD size = 2 (count) + N_IFD_ENTRIES*12 (entries) + 4 (next_ifd_offset)
    const IFD_SIZE: u32 = 2 + N_IFD_ENTRIES as u32 * 12 + 4;
    let strip_offsets_arr: u32 = IFD_OFFSET + IFD_SIZE;   // = 170
    let strip_bytecounts_arr: u32 = strip_offsets_arr + n_strips * 4;
    let xres_off: u32 = strip_bytecounts_arr + n_strips * 4;
    let yres_off: u32 = xres_off + 8;
    let image_data_off: u32 = yres_off + 8;

    // ── Open file ─────────────────────────────────────────────────────────────
    let file = std::fs::File::create(path)?;
    let mut w = BufWriter::with_capacity(1 << 23, file);  // 8 MiB write buffer

    // ── TIFF header ───────────────────────────────────────────────────────────
    w.write_all(b"II")?;                           // little-endian magic
    w.write_all(&42u16.to_le_bytes())?;            // TIFF magic version
    w.write_all(&IFD_OFFSET.to_le_bytes())?;       // offset to first IFD

    // ── IFD ───────────────────────────────────────────────────────────────────
    w.write_all(&N_IFD_ENTRIES.to_le_bytes())?;

    // Helper: write one 12-byte IFD entry.
    // type codes: SHORT=3, LONG=4, RATIONAL=5.
    // Tags must be in ascending numeric order (TIFF spec requirement).
    fn entry(w: &mut impl Write, tag: u16, typ: u16, count: u32, val: u32) -> std::io::Result<()> {
        w.write_all(&tag.to_le_bytes())?;
        w.write_all(&typ.to_le_bytes())?;
        w.write_all(&count.to_le_bytes())?;
        w.write_all(&val.to_le_bytes())?;
        Ok(())
    }

    entry(&mut w, 256, 4, 1, width as u32)?;                   // ImageWidth
    entry(&mut w, 257, 4, 1, height as u32)?;                  // ImageLength
    entry(&mut w, 258, 3, 1, 32)?;                             // BitsPerSample = 32
    entry(&mut w, 259, 3, 1, 1)?;                              // Compression = None
    entry(&mut w, 262, 3, 1, 1)?;                              // PhotometricInterp = MinIsBlack
    entry(&mut w, 273, 4, n_strips, strip_offsets_arr)?;       // StripOffsets → array
    entry(&mut w, 277, 3, 1, 1)?;                              // SamplesPerPixel = 1
    entry(&mut w, 278, 4, 1, ROWS_PER_STRIP)?;                 // RowsPerStrip
    entry(&mut w, 279, 4, n_strips, strip_bytecounts_arr)?;    // StripByteCounts → array
    entry(&mut w, 282, 5, 1, xres_off)?;                       // XResolution → rational
    entry(&mut w, 283, 5, 1, yres_off)?;                       // YResolution → rational
    entry(&mut w, 296, 3, 1, 1)?;                              // ResolutionUnit = No absolute unit
    entry(&mut w, 339, 3, 1, 3)?;                              // SampleFormat = IEEE float

    w.write_all(&0u32.to_le_bytes())?;                         // next IFD offset = 0

    // ── Strip offsets array ────────────────────────────────────────────────────
    for s in 0..n_strips {
        // Use u64 arithmetic to avoid u32 overflow on large files, then cast.
        let off = image_data_off as u64 + s as u64 * full_strip_bytes as u64;
        w.write_all(&(off as u32).to_le_bytes())?;
    }

    // ── Strip byte counts array ────────────────────────────────────────────────
    for s in 0..n_strips {
        let bc = if s == n_strips - 1 { last_strip_bytes } else { full_strip_bytes };
        w.write_all(&bc.to_le_bytes())?;
    }

    // ── Resolution rationals (72/1 DPI placeholder) ───────────────────────────
    w.write_all(&72u32.to_le_bytes())?;   // X numerator
    w.write_all(&1u32.to_le_bytes())?;    // X denominator
    w.write_all(&72u32.to_le_bytes())?;   // Y numerator
    w.write_all(&1u32.to_le_bytes())?;    // Y denominator

    // ── Image data ────────────────────────────────────────────────────────────
    // Reinterpret the f32 slice as raw bytes.  Safe on little-endian x86:
    // TIFF II and the host both use little-endian byte order for f32.
    let bytes: &[u8] = unsafe {
        std::slice::from_raw_parts(data.as_ptr() as *const u8, data.len() * 4)
    };
    w.write_all(bytes)?;

    Ok(())
}
