//! Run the full pipeline for the S1B scene and compare against the ASF RTC10
//! GAMMA reference product.
//!
//! Output: `/home/datacube/dev/SARdine/sardine_s1b_tc_db.tiff`  (γ⁰ dB, EPSG:4326)
//!
//! This is the file that `scripts/compare_asf.py` expects (set SARDINE env var
//! or edit the script path).
//!
//! Requires POEORB — already present at orbit_cache/:
//!   S1B_IW_SLC__1SDV_20190123T053348…_POEORB.EOF
//!
//! Usage:
//!   cargo run --example dump_s1b_tc --release
//!   cargo run --features geoid-fetch --example dump_s1b_tc --release

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

// Output path matched by scripts/compare_asf.py (SARDINE variable).
const OUT_DB: &str = "/home/datacube/dev/SARdine/sardine_s1b_tc_db.tiff";

#[cfg(not(feature = "geoid-fetch"))]
fn main() -> Result<(), Box<dyn std::error::Error>> {
    Err("this example requires the `geoid-fetch` Cargo feature. \
         Re-run with: cargo run --release --features geoid-fetch \
         --example dump_s1b_tc"
        .into())
}

#[cfg(feature = "geoid-fetch")]
fn main() -> Result<(), Box<dyn std::error::Error>> {
    // ── Step 1: parse metadata + apply POEORB ─────────────────────────────────
    eprintln!("[1/7] Parsing SAFE metadata + precise orbit ...");
    let scene = parse_safe_directory(std::path::Path::new(SAFE))?;
    let eof_path = std::path::Path::new(EOF_FILE);
    if !eof_path.is_file() {
        return Err(format!(
            "POEORB not found: {}\nExpected at orbit_cache/ inside the SAFE.",
            eof_path.display()
        ).into());
    }
    let eof = parse_eof_file(eof_path)?;
    eprintln!("    using precise orbit: {}", eof_path.display());
    let scene = apply_precise_orbit(scene, &eof)?;

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
        let sigma0 = apply_calibration(&deburst, cal, noise, 0)?;
        drop(deburst);

        sigma0_arrays.push(sigma0);
        swath_metas.push(sw.clone());
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
    drop(inputs);
    drop(sigma0_arrays);
    drop(swath_metas);

    let total_px = merged.data.len();
    let n_valid = merged.data.iter().filter(|v| !v.is_nan()).count();
    let mean_val: f64 = merged.data.iter().filter(|v| !v.is_nan()).map(|&v| v as f64).sum::<f64>()
        / n_valid.max(1) as f64;
    eprintln!("    {} lines × {} samples  ({:.1}% valid)",
        merged.lines, merged.samples, 100.0 * n_valid as f64 / total_px as f64);
    eprintln!("    σ⁰ linear mean={:.4e}", mean_val);

    // ── Step 5: write slant-range TIFF ─────────────────────────────────────────
    let slant_path = "/tmp/sardine_s1b_merged_sigma0.tiff";
    let mb = (merged.data.len() * 4) as f64 / (1 << 20) as f64;
    eprintln!("[5/7] Writing slant-range TIFF {slant_path} ({mb:.0} MiB) ...");
    write_f32_tiff(slant_path, &merged.data, merged.samples, merged.lines)?;

    // ── Step 6: ground range projection ───────────────────────────────────────
    let target_spacing_m = 10.0;
    let az_spacing_m = scene
        .sub_swaths
        .first()
        .ok_or("scene has no subswaths")?
        .azimuth_pixel_spacing_m;
    eprintln!("[6/7] Ground range projection → {target_spacing_m} m (az spacing = {az_spacing_m:.3} m) ...");
    let grd = to_ground_range(&merged, &geo_grids, target_spacing_m, az_spacing_m)?;

    let grd_valid = grd.data.iter().filter(|v| !v.is_nan()).count();
    let grd_mean: f64 = grd.data.iter().filter(|v| !v.is_nan()).map(|&v| v as f64).sum::<f64>()
        / grd_valid.max(1) as f64;
    eprintln!("    GRD: {} lines × {} samples  ({:.1}% valid)",
        grd.lines, grd.samples, 100.0 * grd_valid as f64 / grd.data.len() as f64);
    eprintln!("    σ⁰ mean={:.4e}", grd_mean);
    let grd_path = "/tmp/sardine_s1b_grd_sigma0.tiff";
    let grd_mb = (grd.data.len() * 4) as f64 / (1 << 20) as f64;
    eprintln!("    Writing GRD TIFF {grd_path} ({grd_mb:.0} MiB) ...");
    write_f32_tiff(grd_path, &grd.data, grd.samples, grd.lines)?;

    // ── Step 7: terrain correction ────────────────────────────────────────────
    eprintln!("[7/7] Terrain correction (Range-Doppler geocoding) ...");
    let dem = DemMosaic::load_directory(std::path::Path::new(DEM_DIR))?;
    eprintln!("    Loaded {} DEM tiles from {}", dem.tile_count(), DEM_DIR);

    // Load EGM96 geoid model.  Gated on the `geoid-fetch` feature (see
    // no-feature stub `main` above) so we never silently fall back to
    // `GeoidModel::Zero` (would inject ~40 m vertical bias).
    use sardine_scene::geoid::GeoidModel;
    use sardine_scene::geoid_fetch::fetch_egm96;
    let geoid = {
        let grid = fetch_egm96().map_err(|e| -> Box<dyn std::error::Error> {
            format!(
                "EGM96 grid fetch failed: {e}. Set $SARDINE_GEOID_DIR or \
                 check your network and rerun. This example refuses to fall \
                 back to GeoidModel::Zero."
            )
            .into()
        })?;
        eprintln!("    EGM96 geoid loaded");
        GeoidModel::Egm96(grid)
    };

    let mut tc_cfg = TerrainCorrectionConfig::new(geoid);
    tc_cfg.pixel_spacing_deg = 0.0001;
    tc_cfg.flatten = true;
    // Mask pixels within 3 dB of the per-pixel NESZ (typical S-1 noise floor).
    // This removes low-SNR pixels before dB conversion without using a
    // hardcoded global threshold.
    tc_cfg.noise_floor_margin_db = 3.0;

    let mut tc = terrain_correction(&merged, &scene, &dem, &geo_grids, &tc_cfg)?;

    eprintln!(
        "    TC: {} rows × {} cols, valid {:.1}%",
        tc.rows, tc.cols,
        100.0 * tc.valid_fraction()
    );
    eprintln!(
        "        buckets: valid={} dem_missing={} outside_footprint={} non_converged={} degenerate={} flat_masked={} noise_masked={}",
        tc.valid_pixel_count, tc.dem_missing_count, tc.outside_footprint_count,
        tc.non_converged_count, tc.degenerate_geometry_count, tc.flat_masked_count,
        tc.noise_masked_count,
    );

    // Apply noise-floor mask and convert to dB.
    eprintln!("    Converting γ⁰ → dB (masking noise floor) ...");
    let (db_valid, db_masked) = to_db_inplace(&mut tc.data, 0.0)
        .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidData, e.to_string()))?;
    eprintln!("        valid={db_valid} noise_masked={db_masked}");

    // Write GeoTIFF (EPSG:4326, self-contained).
    eprintln!("    Writing γ⁰ dB GeoTIFF {OUT_DB} ...");
    write_geotiff(OUT_DB, &tc.data, tc.cols, tc.rows, tc.geotransform)
        .map_err(|e| std::io::Error::new(std::io::ErrorKind::Other, e.to_string()))?;

    eprintln!("Done.");
    eprintln!();
    eprintln!("Compare with ASF RTC10 reference:");
    eprintln!("  SARDINE={OUT_DB} python3 scripts/compare_asf.py");

    Ok(())
}

// ─── Minimal uncompressed float32 TIFF writer (same as dump_merged_sigma0) ───
fn write_f32_tiff(
    path: &str,
    data: &[f32],
    width: usize,
    height: usize,
) -> std::io::Result<()> {
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

    const IFD_OFFSET: u32 = 8;
    let ifd_size: u32 = 2 + N_IFD_ENTRIES as u32 * 12 + 4;
    let strip_offsets_offset: u32 = IFD_OFFSET + ifd_size;
    let strip_byte_counts_offset: u32 = strip_offsets_offset + n_strips * 4;
    let xres_offset: u32 = strip_byte_counts_offset + n_strips * 4;
    let yres_offset: u32 = xres_offset + 8;
    let image_data_offset: u32 = yres_offset + 8;

    // Build strip offsets.
    let strip_offsets: Vec<u32> = (0..n_strips)
        .map(|i| image_data_offset + i * full_strip_bytes)
        .collect();
    let strip_byte_counts: Vec<u32> = (0..n_strips)
        .map(|i| {
            if i == n_strips - 1 { last_strip_bytes } else { full_strip_bytes }
        })
        .collect();

    let file = std::fs::File::create(path)?;
    let mut w = BufWriter::new(file);

    // TIFF header.
    w.write_all(b"II")?;          // little-endian
    w.write_all(&42u16.to_le_bytes())?;
    w.write_all(&IFD_OFFSET.to_le_bytes())?;

    // IFD.
    w.write_all(&N_IFD_ENTRIES.to_le_bytes())?;
    let w32 = width as u32;
    let h32 = height as u32;

    fn ifd_entry(tag: u16, typ: u16, count: u32, value: u32) -> [u8; 12] {
        let mut e = [0u8; 12];
        e[0..2].copy_from_slice(&tag.to_le_bytes());
        e[2..4].copy_from_slice(&typ.to_le_bytes());
        e[4..8].copy_from_slice(&count.to_le_bytes());
        e[8..12].copy_from_slice(&value.to_le_bytes());
        e
    }

    w.write_all(&ifd_entry(256, 3, 1, w32))?;    // ImageWidth
    w.write_all(&ifd_entry(257, 3, 1, h32))?;    // ImageLength
    w.write_all(&ifd_entry(258, 3, 1, 32))?;     // BitsPerSample
    w.write_all(&ifd_entry(259, 3, 1, 1))?;      // Compression=None
    w.write_all(&ifd_entry(262, 3, 1, 1))?;      // PhotometricInterpretation=BlackIsZero
    w.write_all(&ifd_entry(273, 4, n_strips, strip_offsets_offset))?;  // StripOffsets
    w.write_all(&ifd_entry(277, 3, 1, 1))?;      // SamplesPerPixel
    w.write_all(&ifd_entry(278, 3, 1, ROWS_PER_STRIP))?; // RowsPerStrip
    w.write_all(&ifd_entry(279, 4, n_strips, strip_byte_counts_offset))?; // StripByteCounts
    w.write_all(&ifd_entry(282, 5, 1, xres_offset))?;  // XResolution
    w.write_all(&ifd_entry(283, 5, 1, yres_offset))?;  // YResolution
    w.write_all(&ifd_entry(284, 3, 1, 1))?;      // PlanarConfiguration=Chunky
    w.write_all(&ifd_entry(339, 3, 1, 3))?;      // SampleFormat=IEEE float
    w.write_all(&0u32.to_le_bytes())?;            // NextIFD=0

    // Strip offsets array.
    for &off in &strip_offsets { w.write_all(&off.to_le_bytes())?; }
    // Strip byte counts array.
    for &bc in &strip_byte_counts { w.write_all(&bc.to_le_bytes())?; }
    // XResolution = 1/1, YResolution = 1/1.
    w.write_all(&1u32.to_le_bytes())?; w.write_all(&1u32.to_le_bytes())?;
    w.write_all(&1u32.to_le_bytes())?; w.write_all(&1u32.to_le_bytes())?;

    // Image data.
    for chunk in data.chunks(width) {
        for &v in chunk {
            w.write_all(&v.to_le_bytes())?;
        }
    }
    w.flush()?;
    Ok(())
}
