//! Radiometric output conversion and self-contained GeoTIFF export.
//!
//! Provides two public operations:
//!
//! 1. [`to_db_inplace`] — mask below-noise pixels and convert linear σ⁰/γ⁰ to dB.
//! 2. [`write_geotiff`] — write a Float32 GeoTIFF with *embedded* EPSG:4326 CRS
//!    (GeoKey tags 33550/33922/34735), no sidecar files required.
//!
//! # Design constraints
//!
//! - No hardcoded Sentinel-1 constants: all thresholds are caller-supplied.
//! - `write_geotiff` validates the geotransform; invalid inputs return `ExportError`.
//! - The unsafe byte-reinterpretation block for f32 → u8 carries a `SAFETY-OK`
//!   annotation explaining why it is correct.

use std::io::{BufWriter, Write};

use thiserror::Error;

use crate::output_crs::OutputCrs;

// ── Error ─────────────────────────────────────────────────────────────────────

#[derive(Debug, Error)]
pub enum ExportError {
    #[error("I/O error writing GeoTIFF: {0}")]
    Io(#[from] std::io::Error),

    #[error("noise_floor_linear must be non-negative; got {0}")]
    InvalidNoiseFloor(f32),

    #[error("geotransform pixel width/height must be positive and finite; got {0:?}")]
    InvalidGeotransform([f64; 6]),

    #[error("COG tile size must be a power of two ≥ 64 and ≤ 4096; got {0}")]
    InvalidTileSize(u32),
}

// ── dB conversion ─────────────────────────────────────────────────────────────

/// Mask below-noise pixels and convert linear-scale γ⁰/σ⁰ to dB in-place.
///
/// Applied sequentially to every element of `data`:
///
/// 1. NaN values (flagged upstream by terrain flattening, DEM void, etc.) are
///    kept as NaN; they do not appear in any count.
/// 2. Any value with `v <= noise_floor_linear` is set to NaN and counted as
///    `noise_masked`.  Values ≤ 0 are always masked regardless of the threshold
///    because `log10(x ≤ 0)` is undefined.
/// 3. Remaining values are converted: `10 · log₁₀(v)` and counted as `converted`.
///
/// Pass `noise_floor_linear = 0.0` to mask only physically impossible values
/// (≤ 0) without applying any additional SNR-based noise floor.
///
/// # Returns
/// `(converted_count, noise_masked_count)`
pub fn to_db_inplace(
    data: &mut [f32],
    noise_floor_linear: f32,
) -> Result<(usize, usize), ExportError> {
    if noise_floor_linear < 0.0 {
        return Err(ExportError::InvalidNoiseFloor(noise_floor_linear));
    }

    let mut noise_masked = 0usize;
    let mut converted = 0usize;

    for v in data.iter_mut() {
        if v.is_nan() {
            // already masked upstream — keep as NaN, skip both counters
            continue;
        }
        if *v <= noise_floor_linear {
            *v = f32::NAN;
            noise_masked += 1;
        } else {
            *v = 10.0 * v.log10();
            converted += 1;
        }
    }

    Ok((converted, noise_masked))
}

// ── GeoTIFF writer ────────────────────────────────────────────────────────────

/// Write a single-band Float32 GeoTIFF with embedded EPSG:4326 CRS.
///
/// The output file is a fully self-contained GeoTIFF readable by GDAL, QGIS,
/// SNAP, and ArcGIS without any sidecar files.  Classic TIFF is written when
/// the estimated file size is ≤ 4 GiB; BigTIFF is written automatically for
/// larger outputs (e.g. large UTM scenes at 10 m pixel spacing).  The CRS is
/// encoded via three GeoTIFF extension tags:
///
/// | Tag   | Name                | Value                          |
/// |-------|---------------------|--------------------------------|
/// | 33550 | ModelPixelScaleTag  | `[pixel_w, pixel_h, 0]` (f64) |
/// | 33922 | ModelTiepointTag    | `[0, 0, 0, X₀, Y₀, 0]` (f64) |
/// | 34735 | GeoKeyDirectoryTag  | GTModelType=2, GTRasterType=1, GeographicType=4326 |
///
/// The `geotransform` must be a north-up GDAL affine transform
/// `[origin_lon, pixel_width, 0, origin_lat, 0, -pixel_height]`
/// where `pixel_width > 0` and `pixel_height > 0` (i.e. `gt[5] < 0`).
pub fn write_geotiff(
    path: &str,
    data: &[f32],
    cols: usize,
    rows: usize,
    geotransform: [f64; 6],
) -> Result<(), ExportError> {
    // Backward-compatible wrapper: the historical `write_geotiff` always
    // assumed WGS84 lat/lon.  New callers should use
    // [`write_geotiff_with_crs`] to choose any supported CRS.
    write_geotiff_with_crs(path, data, cols, rows, geotransform, &OutputCrs::Wgs84LatLon)
}

/// Write a single-band Float32 GeoTIFF with embedded CRS chosen by the
/// caller.  See [`write_geotiff`] for the WGS84-only wrapper preserved
/// for backward compatibility.
///
/// `geotransform` is interpreted in the CRS's native units: degrees
/// for `OutputCrs::Wgs84LatLon`, metres for the UTM variants.  Passing
/// a degree-spaced geotransform with a UTM CRS produces a syntactically
/// valid but semantically nonsensical GeoTIFF; the caller is responsible
/// for keeping units consistent.  The terrain-correction pipeline does
/// this consistently because both the grid construction and the CRS
/// flow from the same [`OutputCrs`] value.
pub fn write_geotiff_with_crs(
    path: &str,
    data: &[f32],
    cols: usize,
    rows: usize,
    geotransform: [f64; 6],
    crs: &OutputCrs,
) -> Result<(), ExportError> {
    let pixel_w = geotransform[1];
    let pixel_h = -geotransform[5];
    if !pixel_w.is_finite() || pixel_w <= 0.0 || !pixel_h.is_finite() || pixel_h <= 0.0 {
        return Err(ExportError::InvalidGeotransform(geotransform));
    }

    let geo_keys = crs.geo_key_directory();

    // SAFETY-OK: little-endian host enforced by compile_error in inner; f32
    // backing bytes therefore match the on-disk TIFF II byte order.  The
    // pointer/length arithmetic is sound: data.len() * 4 is the exact byte
    // count of the slice's allocation (f32 is 4 bytes, no padding).
    let bytes: &[u8] = unsafe {
        std::slice::from_raw_parts(data.as_ptr() as *const u8, data.len() * 4)
    };
    if needs_bigtiff(rows, cols, 4, geo_keys.len()) {
        return write_bigtiff_raw_inner(
            path, bytes, cols, rows, geotransform,
            /*bits_per_sample=*/ 32,
            /*sample_format=*/   3, // IEEE float
            /*geo_keys=*/        Some(&geo_keys),
        )
        .map_err(ExportError::Io);
    }
    write_geotiff_raw_inner(
        path, bytes, cols, rows, geotransform,
        /*bits_per_sample=*/ 32,
        /*sample_format=*/   3, // IEEE float
        /*geo_keys=*/        Some(&geo_keys),
    )
    .map_err(ExportError::Io)
}

/// Write a single-band 8-bit unsigned-int GeoTIFF with embedded EPSG:4326 CRS.
///
/// Mirrors [`write_geotiff`] but for `u8` rasters such as the per-pixel
/// quality mask emitted by
/// [`crate::terrain_correction::terrain_correction`] (see
/// [`crate::terrain_correction::mask_code`] for value definitions).
///
/// BigTIFF is written automatically when the estimated file size exceeds 4 GiB
/// (1 byte per pixel, so the threshold is ~65k × 65k — larger than any
/// realistic S-1 IW scene at normal quality-mask resolution).
pub fn write_geotiff_u8(
    path: &str,
    data: &[u8],
    cols: usize,
    rows: usize,
    geotransform: [f64; 6],
) -> Result<(), ExportError> {
    // Backward-compatible WGS84-only wrapper; see [`write_geotiff_u8_with_crs`].
    write_geotiff_u8_with_crs(path, data, cols, rows, geotransform, &OutputCrs::Wgs84LatLon)
}

/// Write a single-band 8-bit unsigned-int GeoTIFF with embedded CRS
/// chosen by the caller.  See [`write_geotiff_with_crs`] for unit
/// conventions.
pub fn write_geotiff_u8_with_crs(
    path: &str,
    data: &[u8],
    cols: usize,
    rows: usize,
    geotransform: [f64; 6],
    crs: &OutputCrs,
) -> Result<(), ExportError> {
    let pixel_w = geotransform[1];
    let pixel_h = -geotransform[5];
    if !pixel_w.is_finite() || pixel_w <= 0.0 || !pixel_h.is_finite() || pixel_h <= 0.0 {
        return Err(ExportError::InvalidGeotransform(geotransform));
    }
    let geo_keys = crs.geo_key_directory();
    if needs_bigtiff(rows, cols, 1, geo_keys.len()) {
        return write_bigtiff_raw_inner(
            path, data, cols, rows, geotransform,
            /*bits_per_sample=*/ 8,
            /*sample_format=*/   1, // unsigned integer
            /*geo_keys=*/        Some(&geo_keys),
        )
        .map_err(ExportError::Io);
    }
    write_geotiff_raw_inner(
        path, data, cols, rows, geotransform,
        /*bits_per_sample=*/ 8,
        /*sample_format=*/   1, // unsigned integer
        /*geo_keys=*/        Some(&geo_keys),
    )
    .map_err(ExportError::Io)
}

/// Write a single-band Float32 plain TIFF with **no CRS** (no GeoTIFF tags).
///
/// Use this for radar-geometry rasters such as the GRD output of
/// [`crate::ground_range::to_ground_range`], where the image is in
/// (azimuth-line × ground-range) coordinates and *not* in a map projection.
/// Writing such a raster with [`write_geotiff`] would silently fabricate an
/// EPSG:4326 CRS that does not describe the data — exactly the silent-
/// fallback failure mode the project rules forbid.
///
/// The pixel spacings are **not** encoded in the TIFF; record them in the
/// provenance sidecar instead so consumers can reason about scale without
/// relying on a fabricated CRS.
///
/// BigTIFF is written automatically when the estimated file size exceeds 4 GiB
/// (slightly conservative estimate because the no-CRS file omits the GeoTIFF
/// tag payload, well below the rounding).
pub fn write_tiff_no_crs(
    path: &str,
    data: &[f32],
    cols: usize,
    rows: usize,
) -> Result<(), ExportError> {
    // SAFETY-OK: little-endian host enforced by compile_error in inner; f32
    // backing bytes therefore match the on-disk TIFF II byte order.  The
    // pointer/length arithmetic is sound: data.len() * 4 is the exact byte
    // count of the slice's allocation (f32 is 4 bytes, no padding).
    let bytes: &[u8] = unsafe {
        std::slice::from_raw_parts(data.as_ptr() as *const u8, data.len() * 4)
    };
    if needs_bigtiff(rows, cols, 4, 0) {
        return write_bigtiff_raw_inner(
            path, bytes, cols, rows,
            /*geotransform=*/ [0.0; 6],
            /*bits_per_sample=*/ 32,
            /*sample_format=*/   3, // IEEE float
            /*geo_keys=*/        None,
        )
        .map_err(ExportError::Io);
    }
    write_geotiff_raw_inner(
        path, bytes, cols, rows,
        /*geotransform=*/ [0.0; 6], // unused on the no-CRS path; placeholder only
        /*bits_per_sample=*/ 32,
        /*sample_format=*/   3, // IEEE float
        /*geo_keys=*/        None,
    )
    .map_err(ExportError::Io)
}

// ── Cloud-Optimised GeoTIFF writer ───────────────────────────────────────────

/// Write a single-band Float32 Cloud-Optimised GeoTIFF (COG) with embedded CRS.
///
/// A COG is a standard GeoTIFF with two structural requirements that enable
/// efficient HTTP range-request access:
///
/// 1. **Tiled layout** — the raster is divided into square tiles (default
///    512 × 512 pixels) so any spatial window can be fetched without reading
///    the whole file.
/// 2. **Overview-before-fullres ordering** — reduced-resolution overviews
///    (halved in both dimensions at each level) are written first; the
///    full-resolution image data follows.  This means low-zoom tile requests
///    can be satisfied by fetching the beginning of the file.
///
/// The file is written in a single sequential pass using a two-phase approach:
/// all byte offsets are computed upfront, then the file is written from
/// beginning to end without seeking.
///
/// # Arguments
/// - `path`          — output file path (created or truncated).
/// - `data`          — full-resolution row-major Float32 raster, length =
///                     `rows × cols`.  `NaN` is used as the nodata sentinel.
/// - `cols`, `rows`  — raster dimensions in pixels.
/// - `geotransform`  — GDAL north-up affine transform
///                     `[x_origin, pixel_w, 0, y_origin, 0, -pixel_h]` in
///                     CRS-native units.
/// - `crs`           — output coordinate reference system; controls the
///                     GeoKey directory embedded in the TIFF.
/// - `tile_size`     — tile edge length in pixels; must be a power of two
///                     between 64 and 4096 inclusive.  Use 512 for the GDAL
///                     default.
///
/// # Overviews
/// Overview levels are generated by 2×2 bilinear box-averaging (i.e. each
/// output pixel is the mean of the four nearest-neighbour source pixels;
/// `NaN` inputs propagate to `NaN` output only when all four are NaN — a
/// 2×2 patch with ≥1 valid pixel produces a finite mean from the valid ones).
/// Levels are added until the overview is smaller than `tile_size` in both
/// dimensions.
///
/// # BigTIFF
/// BigTIFF (version 43) is written automatically when the estimated total file
/// size exceeds 4 GiB.
pub fn write_cog_with_crs(
    path: &str,
    data: &[f32],
    cols: usize,
    rows: usize,
    geotransform: [f64; 6],
    crs: &OutputCrs,
    tile_size: u32,
) -> Result<(), ExportError> {
    let pixel_w = geotransform[1];
    let pixel_h = -geotransform[5];
    if !pixel_w.is_finite() || pixel_w <= 0.0 || !pixel_h.is_finite() || pixel_h <= 0.0 {
        return Err(ExportError::InvalidGeotransform(geotransform));
    }
    if tile_size < 64 || tile_size > 4096 || !tile_size.is_power_of_two() {
        return Err(ExportError::InvalidTileSize(tile_size));
    }

    let geo_keys = crs.geo_key_directory();

    // SAFETY-OK: little-endian host enforced by compile_error in write_cog_raw_inner;
    // f32 backing bytes match on-disk TIFF II byte order.  data.len() * 4 is the
    // exact byte count of the f32 slice (no padding, no alignment gap).
    let bytes: &[u8] = unsafe {
        std::slice::from_raw_parts(data.as_ptr() as *const u8, data.len() * 4)
    };

    drop(bytes); // not needed; tile writer reads from the f32 slice directly
    write_cog_raw_inner(
        path, data, cols, rows, geotransform,
        tile_size, &geo_keys,
    )
    .map_err(ExportError::Io)
}

// ── Overview generation ───────────────────────────────────────────────────────

/// Downsample `src` (dims `src_cols × src_rows`) by 2× in each axis using
/// 2×2 bilinear box-averaging.  NaN inputs are skipped: a 2×2 patch produces
/// a finite mean from whichever of the four pixels are not NaN; a patch where
/// all four are NaN produces NaN.
fn downsample_2x(src: &[f32], src_cols: usize, src_rows: usize) -> (Vec<f32>, usize, usize) {
    let dst_cols = (src_cols + 1) / 2;
    let dst_rows = (src_rows + 1) / 2;
    let mut dst = vec![f32::NAN; dst_cols * dst_rows];
    for dr in 0..dst_rows {
        for dc in 0..dst_cols {
            let sr0 = dr * 2;
            let sc0 = dc * 2;
            let mut sum = 0.0f64;
            let mut count = 0u32;
            for dr2 in 0..2usize {
                for dc2 in 0..2usize {
                    let sr = sr0 + dr2;
                    let sc = sc0 + dc2;
                    if sr < src_rows && sc < src_cols {
                        let v = src[sr * src_cols + sc];
                        if !v.is_nan() {
                            sum += v as f64;
                            count += 1;
                        }
                    }
                }
            }
            if count > 0 {
                dst[dr * dst_cols + dc] = (sum / count as f64) as f32;
            }
            // else: all NaN → dst stays NaN
        }
    }
    (dst, dst_cols, dst_rows)
}

// ── COG inner writer (always BigTIFF-capable) ─────────────────────────────────

/// Write a tiled multi-IFD COG.  Called only from `write_cog_with_crs`.
///
/// # File layout (COG ordering, all overviews before full-res)
///
/// ```text
/// [Header]
/// [Ghost IFD with GDAL_STRUCTURAL_METADATA — optional; omitted here for simplicity]
/// [IFD_ovr_N]  ... overviews in order largest→smallest
/// [IFD_ovr_1]
/// [IFD_full]
/// [GeoTIFF data blocks: pixel scale, tiepoint, GeoKey dir]
/// [Tile data ovr_N]
/// ...
/// [Tile data ovr_1]
/// [Tile data full-res]
/// ```
///
/// We use **BigTIFF unconditionally** when total_bytes > 4 GiB, otherwise
/// classic TIFF.  Because the two formats differ in IFD entry size (12 vs
/// 20 bytes) and offset width (4 vs 8 bytes), we choose the format first,
/// compute all offsets, then emit the file in one pass.
fn write_cog_raw_inner(
    path: &str,
    full_f32: &[f32],    // full-resolution row-major Float32 raster
    cols: usize,
    rows: usize,
    geotransform: [f64; 6],
    tile_size: u32,
    geo_keys: &[u16],
) -> std::io::Result<()> {
    #[cfg(not(target_endian = "little"))]
    compile_error!("write_cog_raw_inner assumes a little-endian host");

    let ts = tile_size as usize;

    // ── Build overview pyramid ────────────────────────────────────────────────
    // Overviews from coarsest (index 0) to finest (index n_ovr-1).
    // The vector stores (f32_data, cols, rows) per level.
    // We stop when the overview is smaller than one tile in both dimensions.
    let mut overviews: Vec<(Vec<f32>, usize, usize)> = Vec::new();
    {
        let (mut ov_data, mut ov_cols, mut ov_rows) =
            downsample_2x(full_f32, cols, rows);
        while ov_cols > ts || ov_rows > ts {
            overviews.push((ov_data.clone(), ov_cols, ov_rows));
            let next = downsample_2x(&ov_data, ov_cols, ov_rows);
            (ov_data, ov_cols, ov_rows) = next;
        }
        // Push the smallest overview (≤ tile_size × tile_size)
        overviews.push((ov_data, ov_cols, ov_rows));
    }
    // Reverse so overviews[0] = coarsest, overviews[last] = finest
    overviews.reverse();

    let n_ovr = overviews.len();

    // ── Decide TIFF format ────────────────────────────────────────────────────
    // Rough image-data size estimate; if > 4 GiB use BigTIFF.
    let image_bytes_estimate: u64 = (rows as u64)
        .saturating_mul(cols as u64)
        .saturating_mul(4);
    let use_bigtiff = image_bytes_estimate > u32::MAX as u64;

    // ── Tile grid helper ──────────────────────────────────────────────────────
    // Returns (n_tiles_x, n_tiles_y, n_tiles_total) for a given image dimension.
    let tile_grid = |c: usize, r: usize| -> (usize, usize, usize) {
        let nx = (c + ts - 1) / ts;
        let ny = (r + ts - 1) / ts;
        (nx, ny, nx * ny)
    };

    // ── Compute IFD sizes and data-block sizes ────────────────────────────────
    //
    // Each IFD has 13 base tags + 3 GeoTIFF tags = 16 tags.
    // Classic TIFF: IFD = 2 + 16×12 + 4 = 198 bytes.
    // BigTIFF:      IFD = 8 + 16×20 + 8 = 336 bytes.
    //
    // Each IFD also has out-of-line data blocks:
    //   - TileOffsets array   (n_tiles × 4 bytes classic, n_tiles × 8 bytes big)
    //   - TileByteCounts array(n_tiles × 4 bytes classic, n_tiles × 8 bytes big)
    //   - GeoTIFF blocks only in the full-res IFD:
    //       ModelPixelScaleTag  24 bytes
    //       ModelTiepointTag    48 bytes
    //       GeoKeyDirectoryTag  geo_keys.len() × 2 bytes
    //
    // IFDs appear in COG order (coarsest overview first → full-res last).
    // After all IFDs come the GeoTIFF auxiliary data, then tile data in the
    // same IFD order (overviews before full-res).

    let ifd_entries: u64 = 17; // 14 base (incl. tile) + 3 GeoTIFF = 17
    let geo_keys_bytes = (geo_keys.len() as u64) * 2;

    // Per-IFD size and tile-array size depend on format.
    let (ifd_fixed_bytes, tile_arr_entry_bytes, header_bytes, next_ifd_ptr_bytes): (u64, u64, u64, u64) =
        if use_bigtiff {
            // BigTIFF: 8 (entry count) + n×20 (entries) + 8 (next IFD ptr)
            (8 + ifd_entries * 20 + 8, 8, 16, 8)
        } else {
            // Classic TIFF: 2 (entry count) + n×12 (entries) + 4 (next IFD ptr)
            (2 + ifd_entries * 12 + 4, 4, 8, 4)
        };
        // next_ifd_ptr_bytes is included in ifd_fixed_bytes; the variable is
        // destructured from the tuple only to document the layout — it is not
        // used independently.
        let _ = next_ifd_ptr_bytes; // SAFETY-OK: dead binding documents tuple layout, no fallback behaviour

    // Compute running file offset, tracking where each IFD and data block go.
    let mut offset: u64 = header_bytes;

    // -- IFD offsets (one per overview level + one for full-res) ----------------
    let mut ifd_offsets: Vec<u64> = Vec::with_capacity(n_ovr + 1);
    // Tile-array offsets (TileOffsets, TileByteCounts) per IFD.
    let mut tile_off_arr_offsets: Vec<u64>   = Vec::with_capacity(n_ovr + 1);
    let mut tile_bc_arr_offsets:  Vec<u64>   = Vec::with_capacity(n_ovr + 1);
    let mut tile_counts:          Vec<usize> = Vec::with_capacity(n_ovr + 1);

    // First lay out all IFDs + their tile arrays (no image data yet).
    // COG order: coarsest overview IFD first, then finest, then full-res.
    let all_levels: Vec<(usize, usize)> = overviews.iter()
        .map(|(_, c, r)| (*c, *r))
        .chain(std::iter::once((cols, rows)))
        .collect();

    for &(c, r) in &all_levels {
        ifd_offsets.push(offset);
        offset += ifd_fixed_bytes;
        let (_, _, n_tiles) = tile_grid(c, r);
        tile_off_arr_offsets.push(offset);
        offset += n_tiles as u64 * tile_arr_entry_bytes;
        tile_bc_arr_offsets.push(offset);
        offset += n_tiles as u64 * tile_arr_entry_bytes;
        tile_counts.push(n_tiles);
    }

    // GeoTIFF auxiliary data blocks (appended after all IFDs, shared by
    // the full-res IFD only).
    let pixel_scale_off: u64 = offset;
    offset += 24;
    let tiepoint_off: u64 = offset;
    offset += 48;
    let geo_key_dir_off: u64 = offset;
    offset += geo_keys_bytes;

    // -- Tile data offsets (same IFD order) ------------------------------------
    // Each tile is exactly tile_size × tile_size × 4 bytes, except edge tiles
    // which are smaller (we will pad to full tile size in the writer below so
    // tile byte counts are uniform and equal to ts × ts × 4).
    let bytes_per_tile: u64 = (ts as u64) * (ts as u64) * 4;

    let mut tile_data_offsets: Vec<Vec<u64>> = Vec::with_capacity(n_ovr + 1);
    for &n in &tile_counts {
        let mut offsets_for_ifd = Vec::with_capacity(n);
        for _ in 0..n {
            offsets_for_ifd.push(offset);
            offset += bytes_per_tile;
        }
        tile_data_offsets.push(offsets_for_ifd);
    }

    // ── Open file and write ───────────────────────────────────────────────────
    let file = std::fs::File::create(path)?;
    let mut w = BufWriter::with_capacity(1 << 23, file);

    // ── Header ────────────────────────────────────────────────────────────────
    if use_bigtiff {
        w.write_all(b"II")?;
        w.write_all(&43u16.to_le_bytes())?;  // BigTIFF version
        w.write_all(&8u16.to_le_bytes())?;   // bytesize of offsets
        w.write_all(&0u16.to_le_bytes())?;   // constant 0
        w.write_all(&(ifd_offsets[0] as u64).to_le_bytes())?;
    } else {
        w.write_all(b"II")?;
        w.write_all(&42u16.to_le_bytes())?;  // classic TIFF version
        w.write_all(&(ifd_offsets[0] as u32).to_le_bytes())?;
    }

    // ── Write IFDs ────────────────────────────────────────────────────────────
    let n_ifds = all_levels.len();
    for (ifd_idx, &(c, r)) in all_levels.iter().enumerate() {
        let is_full_res = ifd_idx == n_ifds - 1;
        let (_, _, n_tiles) = tile_grid(c, r);

        let next_ifd_off: u64 = if ifd_idx + 1 < n_ifds {
            ifd_offsets[ifd_idx + 1]
        } else {
            0 // last IFD
        };

        let tile_off_arr = tile_off_arr_offsets[ifd_idx];
        let tile_bc_arr  = tile_bc_arr_offsets[ifd_idx];

        if use_bigtiff {
            // BigTIFF: entry count as u64
            w.write_all(&(ifd_entries as u64).to_le_bytes())?;

            fn e64(w: &mut impl Write, tag: u16, typ: u16, count: u64, val: u64) -> std::io::Result<()> {
                w.write_all(&tag.to_le_bytes())?;
                w.write_all(&typ.to_le_bytes())?;
                w.write_all(&count.to_le_bytes())?;
                w.write_all(&val.to_le_bytes())?;
                Ok(())
            }

            const RATIONAL_72_1: u64 = (1u64 << 32) | 72u64;

            e64(&mut w, 256,  4, 1,           c as u64)?;           // ImageWidth
            e64(&mut w, 257,  4, 1,           r as u64)?;           // ImageLength
            e64(&mut w, 258,  3, 1,           32)?;                  // BitsPerSample
            e64(&mut w, 259,  3, 1,           1)?;                   // Compression = None
            e64(&mut w, 262,  3, 1,           1)?;                   // PhotometricInterp
            e64(&mut w, 277,  3, 1,           1)?;                   // SamplesPerPixel
            e64(&mut w, 282,  5, 1,           RATIONAL_72_1)?;       // XResolution inline
            e64(&mut w, 283,  5, 1,           RATIONAL_72_1)?;       // YResolution inline
            e64(&mut w, 296,  3, 1,           1)?;                   // ResolutionUnit
            e64(&mut w, 322,  4, 1,           tile_size as u64)?;    // TileWidth
            e64(&mut w, 323,  4, 1,           tile_size as u64)?;    // TileLength
            e64(&mut w, 324, 16, n_tiles as u64, tile_off_arr)?;     // TileOffsets (LONG8)
            e64(&mut w, 325, 16, n_tiles as u64, tile_bc_arr)?;      // TileByteCounts (LONG8)
            e64(&mut w, 339,  3, 1,           3)?;                   // SampleFormat = IEEE float
            if is_full_res {
                e64(&mut w, 33550, 12, 3,              pixel_scale_off)?;
                e64(&mut w, 33922, 12, 6,              tiepoint_off)?;
                e64(&mut w, 34735,  3, geo_keys.len() as u64, geo_key_dir_off)?;
            } else {
                // Placeholder GeoTIFF tags pointing to the full-res GeoTIFF blocks.
                // Overview IFDs share the same spatial reference — GDAL expects them.
                e64(&mut w, 33550, 12, 3, pixel_scale_off)?;
                e64(&mut w, 33922, 12, 6, tiepoint_off)?;
                e64(&mut w, 34735,  3, geo_keys.len() as u64, geo_key_dir_off)?;
            }
            w.write_all(&next_ifd_off.to_le_bytes())?;
        } else {
            // Classic TIFF
            w.write_all(&(ifd_entries as u16).to_le_bytes())?;

            fn e32(w: &mut impl Write, tag: u16, typ: u16, count: u32, val: u32) -> std::io::Result<()> {
                w.write_all(&tag.to_le_bytes())?;
                w.write_all(&typ.to_le_bytes())?;
                w.write_all(&count.to_le_bytes())?;
                w.write_all(&val.to_le_bytes())?;
                Ok(())
            }

            // XRes/YRes RATIONAL need a data block in classic TIFF. Re-use tiepoint_off + 72
            // for the two rational values. But we haven't pre-allocated a slot for this.
            // To keep layout simple: write XRes/YRes as inline LONG (type 4) with value 72
            // (pixels per unit, unit = no absolute unit). This is a common simplification
            // and GDAL ignores XRes/YRes for georeferenced rasters.
            // (The alternative is to allocate a separate 8-byte rational block, which would
            // complicate the offset arithmetic above significantly for a value GDAL ignores.)
            e32(&mut w, 256, 4, 1,           c as u32)?;
            e32(&mut w, 257, 4, 1,           r as u32)?;
            e32(&mut w, 258, 3, 1,           32)?;
            e32(&mut w, 259, 3, 1,           1)?;
            e32(&mut w, 262, 3, 1,           1)?;
            e32(&mut w, 277, 3, 1,           1)?;
            e32(&mut w, 282, 4, 1,           72)?;  // XResolution (simplified LONG)
            e32(&mut w, 283, 4, 1,           72)?;  // YResolution (simplified LONG)
            e32(&mut w, 296, 3, 1,           1)?;
            e32(&mut w, 322, 4, 1,           tile_size)?;
            e32(&mut w, 323, 4, 1,           tile_size)?;
            e32(&mut w, 324, 4, n_tiles as u32, tile_off_arr as u32)?;
            e32(&mut w, 325, 4, n_tiles as u32, tile_bc_arr  as u32)?;
            e32(&mut w, 339, 3, 1,           3)?;
            e32(&mut w, 33550, 12, 3,              pixel_scale_off as u32)?;
            e32(&mut w, 33922, 12, 6,              tiepoint_off    as u32)?;
            e32(&mut w, 34735,  3, geo_keys.len() as u32, geo_key_dir_off as u32)?;
            w.write_all(&(next_ifd_off as u32).to_le_bytes())?;
        }

        // ── Tile offset arrays ────────────────────────────────────────────────
        for &off in &tile_data_offsets[ifd_idx] {
            if use_bigtiff {
                w.write_all(&off.to_le_bytes())?;
            } else {
                w.write_all(&(off as u32).to_le_bytes())?;
            }
        }
        // ── Tile byte-count arrays ────────────────────────────────────────────
        for _ in 0..n_tiles {
            if use_bigtiff {
                w.write_all(&bytes_per_tile.to_le_bytes())?;
            } else {
                w.write_all(&(bytes_per_tile as u32).to_le_bytes())?;
            }
        }
    }

    // ── GeoTIFF auxiliary data (written once, shared by all IFDs) ────────────
    {
        // ModelPixelScaleTag: [pixel_width, pixel_height, 0.0]
        w.write_all(&geotransform[1].to_le_bytes())?;
        w.write_all(&(-geotransform[5]).to_le_bytes())?;
        w.write_all(&0f64.to_le_bytes())?;
        // ModelTiepointTag: [I=0, J=0, K=0, X=x_origin, Y=y_origin, Z=0]
        w.write_all(&0f64.to_le_bytes())?;
        w.write_all(&0f64.to_le_bytes())?;
        w.write_all(&0f64.to_le_bytes())?;
        w.write_all(&geotransform[0].to_le_bytes())?;
        w.write_all(&geotransform[3].to_le_bytes())?;
        w.write_all(&0f64.to_le_bytes())?;
        // GeoKeyDirectoryTag
        for k in geo_keys {
            w.write_all(&k.to_le_bytes())?;
        }
    }

    // ── Tile data (COG order: overviews first, full-res last) ─────────────────
    // Temporary tile buffer (padded to exactly ts × ts × 4 bytes).
    let mut tile_buf = vec![f32::NAN; ts * ts];

    let all_f32_data: Vec<&[f32]> = overviews.iter()
        .map(|(d, _, _)| d.as_slice())
        .chain(std::iter::once(full_f32))
        .collect();

    for (ifd_idx, &(c, r)) in all_levels.iter().enumerate() {
        let src = all_f32_data[ifd_idx];
        let (nx, ny, _) = tile_grid(c, r);

        for ty in 0..ny {
            for tx in 0..nx {
                // Fill tile buffer with NaN (for out-of-raster padding).
                tile_buf.fill(f32::NAN);
                let row_start = ty * ts;
                let col_start = tx * ts;
                let row_end   = (row_start + ts).min(r);
                let col_end   = (col_start + ts).min(c);
                for row in row_start..row_end {
                    let src_row = &src[row * c + col_start .. row * c + col_end];
                    let dst_off = (row - row_start) * ts;
                    tile_buf[dst_off .. dst_off + src_row.len()].copy_from_slice(src_row);
                }
                // Write tile as raw f32 LE bytes.
                // SAFETY-OK: tile_buf is a Vec<f32>; its length is ts*ts;
                // the byte slice is ts*ts*4 bytes which is the correct size;
                // the host is little-endian (enforced by compile_error above).
                let tile_bytes: &[u8] = unsafe {
                    std::slice::from_raw_parts(tile_buf.as_ptr() as *const u8, ts * ts * 4)
                };
                w.write_all(tile_bytes)?;
            }
        }
    }
    drop(tile_buf);
    // Ensure all buffered bytes reach the OS.
    use std::io::Write as _;
    w.flush()?;

    Ok(())
}


/// Returns `true` when the estimated classic-TIFF file size would exceed the
/// 4 GiB (2^32 − 1 byte) limit imposed by 32-bit TIFF offsets, indicating
/// that a BigTIFF should be written instead.
///
/// The size estimate is deliberately conservative: it uses the classic-TIFF
/// header and IFD layout, which is slightly smaller than BigTIFF.  This means
/// BigTIFF is selected slightly earlier than strictly necessary, which is safe
/// (BigTIFF readers handle any size) and simplifies the threshold logic.
///
/// # Arguments
/// - `bytes_per_pixel` allows reuse for `f32` (4) and `u8` (1) writers.
/// - `geo_keys_n_shorts` is the length (in `u16` elements) of the
///   GeoKeyDirectoryTag block, or 0 when no GeoTIFF tags are present.
fn needs_bigtiff(
    rows: usize,
    cols: usize,
    bytes_per_pixel: u64,
    geo_keys_n_shorts: usize,
) -> bool {
    const ROWS_PER_STRIP: u64 = 64;
    let geo_keys_bytes: u64 = (geo_keys_n_shorts as u64) * 2;
    let fixed_prelude_bytes: u64 = if geo_keys_n_shorts > 0 {
        // 8 (TIFF hdr) + 198 (IFD: 2 + 16×12 + 4) + 16 (XYres) + 24 + 48 + geo_keys_bytes
        8 + 198 + 16 + 24 + 48 + geo_keys_bytes
    } else {
        // 8 (TIFF hdr) + 162 (IFD: 2 + 13×12 + 4) + 16 (XYres)
        8 + 162 + 16
    };
    let image_bytes = (rows as u64)
        .saturating_mul(cols as u64)
        .saturating_mul(bytes_per_pixel);
    let n_strips = (rows as u64 + ROWS_PER_STRIP - 1) / ROWS_PER_STRIP;
    let strip_table_bytes = n_strips.saturating_mul(8); // 4 offsets + 4 bytecounts per strip
    let total_bytes = fixed_prelude_bytes
        .saturating_add(strip_table_bytes)
        .saturating_add(image_bytes);
    total_bytes > u32::MAX as u64
}

fn write_geotiff_raw_inner(
    path: &str,
    bytes: &[u8],
    cols: usize,
    rows: usize,
    geotransform: [f64; 6],
    bits_per_sample: u16,
    sample_format: u16,
    geo_keys: Option<&[u16]>,
) -> std::io::Result<()> {
    // All TIFF tag/header data and raster bytes use little-endian byte order.
    // This is only correct on a little-endian host.
    #[cfg(not(target_endian = "little"))]
    compile_error!("write_geotiff assumes a little-endian host");

    debug_assert!(
        bits_per_sample == 8 || bits_per_sample == 32,
        "only 8-bit and 32-bit single-band rasters are supported"
    );
    let bytes_per_pixel = (bits_per_sample / 8) as u32;
    debug_assert_eq!(
        bytes.len() as u64,
        rows as u64 * cols as u64 * bytes_per_pixel as u64,
        "raw byte buffer length must match rows × cols × bytes_per_pixel"
    );

    // ── Layout constants ──────────────────────────────────────────────────────
    //
    // Classic TIFF with either 16 IFD entries (13 base + 3 GeoTIFF tags) when
    // a GeoKey directory is supplied, or 13 entries (base only) when not.
    // The size of the GeoKeyDirectoryTag block varies with the CRS — 16 SHORTs
    // for WGS84 lat/lon (3 keys) or 20 SHORTs for UTM (4 keys) — and shifts
    // the image data offset accordingly.  All offsets are computed from the
    // caller-supplied `geo_keys` slice; the inner writer never assumes a
    // specific CRS.
    //
    // File map (little-endian TIFF, classic, single IFD), GeoTIFF case:
    //
    //   offset 0        : 8-byte TIFF header (II + 42 + IFD_OFFSET)
    //   offset 8        : IFD  (2 + 16×12 + 4 = 198 bytes)
    //   offset 206      : strip offsets array        (n_strips × 4 bytes)
    //   offset 206+S×4  : strip byte-counts array    (n_strips × 4 bytes)
    //   offset 206+S×8  : XResolution RATIONAL        (8 bytes)
    //   offset 214+S×8  : YResolution RATIONAL        (8 bytes)
    //   offset 222+S×8  : ModelPixelScaleTag data     (3 × 8 = 24 bytes)
    //   offset 246+S×8  : ModelTiepointTag data       (6 × 8 = 48 bytes)
    //   offset 294+S×8  : GeoKeyDirectoryTag data     (geo_keys.len() × 2 bytes)
    //   offset 294+S×8 +geo_keys_bytes : image data   (rows × cols × bpp bytes)
    //
    // No-CRS case shrinks the IFD to 13 entries and omits the three GeoTIFF
    // data blocks, so:
    //
    //   offset 8        : IFD  (2 + 13×12 + 4 = 162 bytes)
    //   offset 170      : strip offsets array        (n_strips × 4 bytes)
    //   offset 170+S×4  : strip byte-counts array    (n_strips × 4 bytes)
    //   offset 170+S×8  : XResolution RATIONAL        (8 bytes)
    //   offset 178+S×8  : YResolution RATIONAL        (8 bytes)
    //   offset 186+S×8  : image data                  (rows × cols × bpp bytes)
    //
    // where S = n_strips and bpp = bytes_per_pixel.

    const ROWS_PER_STRIP: u32 = 64;
    let include_geotiff_keys = geo_keys.is_some();
    let n_ifd_entries: u16 = if include_geotiff_keys { 16 } else { 13 };
    let geo_keys_bytes: u32 = match geo_keys {
        Some(k) => (k.len() as u32) * 2,
        None    => 0,
    };

    let n_strips = (rows as u32 + ROWS_PER_STRIP - 1) / ROWS_PER_STRIP;
    let bytes_per_row = (cols as u32) * bytes_per_pixel;
    let full_strip_bytes = ROWS_PER_STRIP * bytes_per_row;
    let last_strip_rows = rows as u32 - (n_strips - 1) * ROWS_PER_STRIP;
    let last_strip_bytes = last_strip_rows * bytes_per_row;

    const IFD_OFFSET: u32 = 8;
    let ifd_size: u32 = 2 + n_ifd_entries as u32 * 12 + 4;

    let strip_offsets_arr: u32   = IFD_OFFSET + ifd_size;
    let strip_bytecounts_arr: u32 = strip_offsets_arr + n_strips * 4;
    let xres_off: u32             = strip_bytecounts_arr + n_strips * 4;
    let yres_off: u32             = xres_off + 8;
    // GeoTIFF extension blocks (only emitted when `include_geotiff_keys`).
    let pixel_scale_off: u32      = yres_off + 8;                    // 3 × f64 = 24 bytes
    let tiepoint_off: u32         = pixel_scale_off + 24;            // 6 × f64 = 48 bytes
    let geo_key_dir_off: u32      = tiepoint_off + 48;               // geo_keys.len() × u16 bytes
    let image_data_off: u32       = if include_geotiff_keys {
        geo_key_dir_off + geo_keys_bytes
    } else {
        yres_off + 8
    };

    // ── Open file ─────────────────────────────────────────────────────────────
    let file = std::fs::File::create(path)?;
    let mut w = BufWriter::with_capacity(1 << 23, file); // 8 MiB write buffer

    // ── TIFF header ───────────────────────────────────────────────────────────
    w.write_all(b"II")?;                           // little-endian magic
    w.write_all(&42u16.to_le_bytes())?;            // TIFF magic
    w.write_all(&IFD_OFFSET.to_le_bytes())?;       // offset to first IFD

    // ── IFD ───────────────────────────────────────────────────────────────────
    // Tags MUST appear in ascending numeric order per the TIFF spec.
    // TIFF type codes: SHORT=3, LONG=4, RATIONAL=5, DOUBLE=12.
    w.write_all(&n_ifd_entries.to_le_bytes())?;

    fn entry(w: &mut impl Write, tag: u16, typ: u16, count: u32, val: u32) -> std::io::Result<()> {
        w.write_all(&tag.to_le_bytes())?;
        w.write_all(&typ.to_le_bytes())?;
        w.write_all(&count.to_le_bytes())?;
        w.write_all(&val.to_le_bytes())?;
        Ok(())
    }

    // ── Base TIFF tags (tags 256–339) ────────────────────────────────────────
    entry(&mut w, 256, 4, 1,        cols as u32)?;             // ImageWidth
    entry(&mut w, 257, 4, 1,        rows as u32)?;             // ImageLength
    entry(&mut w, 258, 3, 1,        bits_per_sample as u32)?; // BitsPerSample (8 or 32)
    entry(&mut w, 259, 3, 1,        1)?;                       // Compression = None
    entry(&mut w, 262, 3, 1,        1)?;                       // PhotometricInterp = MinIsBlack
    entry(&mut w, 273, 4, n_strips, strip_offsets_arr)?;       // StripOffsets → array
    entry(&mut w, 277, 3, 1,        1)?;                       // SamplesPerPixel = 1
    entry(&mut w, 278, 4, 1,        ROWS_PER_STRIP)?;          // RowsPerStrip
    entry(&mut w, 279, 4, n_strips, strip_bytecounts_arr)?;    // StripByteCounts → array
    entry(&mut w, 282, 5, 1,        xres_off)?;                // XResolution → rational
    entry(&mut w, 283, 5, 1,        yres_off)?;                // YResolution → rational
    entry(&mut w, 296, 3, 1,        1)?;                       // ResolutionUnit = No absolute unit
    entry(&mut w, 339, 3, 1,        sample_format as u32)?;    // SampleFormat (1=uint, 3=float)

    // ── GeoTIFF extension tags (33550, 33922, 34735) ─────────────────────────
    // Only emitted when a GeoKey directory was supplied; the no-CRS variant
    // produces a plain TIFF that GIS readers will treat as having no
    // projection, matching the radar-geometry semantics of the GRD writer.
    if let Some(keys) = geo_keys {
        // ModelPixelScaleTag (33550): 3 DOUBLE values stored at pixel_scale_off.
        entry(&mut w, 33550, 12, 3,  pixel_scale_off)?;
        // ModelTiepointTag (33922): 6 DOUBLE values stored at tiepoint_off.
        entry(&mut w, 33922, 12, 6,  tiepoint_off)?;
        // GeoKeyDirectoryTag (34735): N SHORT values stored at geo_key_dir_off.
        entry(&mut w, 34735,  3, keys.len() as u32, geo_key_dir_off)?;
    }

    w.write_all(&0u32.to_le_bytes())?; // next IFD offset = 0 (single IFD)

    // ── Strip offsets ─────────────────────────────────────────────────────────
    for s in 0..n_strips {
        // u64 arithmetic prevents overflow on large files; cast back to u32 is
        // safe because image_data_off + total_image_bytes ≤ 4 GiB (classic TIFF).
        let off = image_data_off as u64 + s as u64 * full_strip_bytes as u64;
        w.write_all(&(off as u32).to_le_bytes())?;
    }

    // ── Strip byte counts ──────────────────────────────────────────────────────
    for s in 0..n_strips {
        let bc = if s == n_strips - 1 { last_strip_bytes } else { full_strip_bytes };
        w.write_all(&bc.to_le_bytes())?;
    }

    // ── Resolution rationals (72/1 DPI placeholder) ───────────────────────────
    w.write_all(&72u32.to_le_bytes())?; // XResolution numerator
    w.write_all(&1u32.to_le_bytes())?;  // XResolution denominator
    w.write_all(&72u32.to_le_bytes())?; // YResolution numerator
    w.write_all(&1u32.to_le_bytes())?;  // YResolution denominator

    // ── ModelPixelScaleTag: [pixel_width, pixel_height, 0.0] ─────────────────
    // Only written when GeoKeys are supplied.  The Tiepoint and GeoKey blocks
    // below are gated on the same `Some(keys)` so the file layout exactly
    // matches the IFD entries emitted above.
    if let Some(keys) = geo_keys {
        // gt[1] = X spacing (positive eastward in CRS native units)
        // -gt[5] = Y spacing magnitude (positive southward)
        let pixel_w = geotransform[1];
        let pixel_h = -geotransform[5];
        w.write_all(&pixel_w.to_le_bytes())?;
        w.write_all(&pixel_h.to_le_bytes())?;
        w.write_all(&0f64.to_le_bytes())?;

        // ── ModelTiepointTag: [I, J, K, X, Y, Z] ─────────────────────────────
        // Tiepoint: pixel (I=0, J=0) corresponds to the upper-left corner of
        // the image, at projected coordinates (X=x_origin, Y=y_origin).
        let x_origin = geotransform[0]; // X of upper-left corner (CRS native units)
        let y_origin = geotransform[3]; // Y of upper-left corner (CRS native units)
        w.write_all(&0f64.to_le_bytes())?;     // I = col 0
        w.write_all(&0f64.to_le_bytes())?;     // J = row 0
        w.write_all(&0f64.to_le_bytes())?;     // K = 0 (unused)
        w.write_all(&x_origin.to_le_bytes())?; // X
        w.write_all(&y_origin.to_le_bytes())?; // Y
        w.write_all(&0f64.to_le_bytes())?;     // Z = 0

        // ── GeoKeyDirectoryTag: caller-supplied SHORT block ─────────────────
        // The block is opaque to this writer; the [`OutputCrs`] caller is
        // responsible for the spec-compliant layout (see
        // [`crate::output_crs::OutputCrs::geo_key_directory`]).
        for k in keys {
            w.write_all(&k.to_le_bytes())?;
        }
    }

    // ── Image data ────────────────────────────────────────────────────────────
    // Caller-supplied raw bytes already match the on-disk byte order (little-
    // endian, enforced by the compile_error above).  For f32 callers, the
    // public `write_geotiff` performs the f32 → bytes reinterpretation; for
    // u8 callers, the bytes are passed straight through.
    w.write_all(bytes)?;

    Ok(())
}

/// BigTIFF writer: same contract as [`write_geotiff_raw_inner`] but emits a
/// BigTIFF (TIFF version 43) with 64-bit file offsets and IFD fields.
///
/// # BigTIFF format differences from classic TIFF
///
/// | Field              | Classic TIFF | BigTIFF  |
/// |--------------------|-------------|----------|
/// | Magic number       | 42          | 43       |
/// | Offset bytesize    | 4 bytes     | 8 bytes  |
/// | Header size        | 8 bytes     | 16 bytes |
/// | IFD entry count    | u16         | u64      |
/// | IFD entry size     | 12 bytes    | 20 bytes |
/// | Next-IFD pointer   | u32         | u64      |
/// | Value/offset field | u32         | u64      |
///
/// Inline value rule: if count × type_size ≤ 8 bytes the value is stored
/// directly in the 8-byte field (left-justified, zero-padded).  Otherwise
/// the field holds a u64 file offset to the data block.
///
/// Consequence relevant here: RATIONAL (2 × u32 = 8 bytes) fits inline in
/// BigTIFF.  XResolution and YResolution therefore need no separate data
/// block, unlike the classic TIFF writer.
///
/// Strip offset and byte-count arrays always use the LONG8 type (code 16,
/// 8-byte unsigned int) and are stored at offset when n_strips > 1, or
/// inline (value = the single strip's actual offset / byte-count) when
/// n_strips == 1.
fn write_bigtiff_raw_inner(
    path: &str,
    bytes: &[u8],
    cols: usize,
    rows: usize,
    geotransform: [f64; 6],
    bits_per_sample: u16,
    sample_format: u16,
    geo_keys: Option<&[u16]>,
) -> std::io::Result<()> {
    // BigTIFF still requires little-endian byte order for the "II" variant.
    #[cfg(not(target_endian = "little"))]
    compile_error!("write_bigtiff assumes a little-endian host");

    debug_assert!(
        bits_per_sample == 8 || bits_per_sample == 32,
        "only 8-bit and 32-bit single-band rasters are supported"
    );
    let bytes_per_pixel = (bits_per_sample / 8) as u64;
    debug_assert_eq!(
        bytes.len() as u64,
        rows as u64 * cols as u64 * bytes_per_pixel,
        "raw byte buffer length must match rows × cols × bytes_per_pixel"
    );

    // ── Layout constants ──────────────────────────────────────────────────────

    const ROWS_PER_STRIP: u64 = 64;
    const BIGTIFF_HEADER_SIZE: u64 = 16;
    const BIGTIFF_ENTRY_SIZE: u64 = 20;

    let include_geotiff_keys = geo_keys.is_some();
    let n_ifd_entries: u64 = if include_geotiff_keys { 16 } else { 13 };
    let geo_keys_bytes: u64 = match geo_keys {
        Some(k) => (k.len() as u64) * 2,
        None    => 0,
    };

    let n_strips = (rows as u64 + ROWS_PER_STRIP - 1) / ROWS_PER_STRIP;
    let bytes_per_row  = cols as u64 * bytes_per_pixel;
    let full_strip_bytes = ROWS_PER_STRIP * bytes_per_row;
    let last_strip_rows  = rows as u64 - (n_strips - 1) * ROWS_PER_STRIP;
    let last_strip_bytes = last_strip_rows * bytes_per_row;

    // IFD: 8-byte count + n_ifd_entries × 20-byte entries + 8-byte next pointer.
    let ifd_size: u64 = 8 + n_ifd_entries * BIGTIFF_ENTRY_SIZE + 8;
    let ifd_end: u64  = BIGTIFF_HEADER_SIZE + ifd_size;

    // Strip arrays are stored at offsets when n_strips > 1.
    // When n_strips == 1, the strip offset and bytecount fit inline in the IFD
    // entry value field (LONG8, 8 bytes each).
    let strips_inline = n_strips == 1;

    let (strip_offsets_arr, strip_bytecounts_arr, after_strips) = if !strips_inline {
        let s_off = ifd_end;
        let s_bc  = ifd_end + n_strips * 8;
        let after = ifd_end + 2 * n_strips * 8;
        (s_off, s_bc, after)
    } else {
        (0u64, 0u64, ifd_end) // 0 = unused sentinel; strips are inline
    };

    // GeoTIFF extension data blocks (only when geo_keys supplied).
    // Note: XResolution and YResolution are RATIONAL (8 bytes) and fit inline;
    //       they do NOT need a separate data block in BigTIFF.
    let (pixel_scale_off, tiepoint_off, geo_key_dir_off, image_data_off) =
        if include_geotiff_keys {
            let ps = after_strips;
            let tp = after_strips + 24;      // after 3 × f64
            let gk = after_strips + 24 + 48; // after 6 × f64
            let im = gk + geo_keys_bytes;
            (ps, tp, gk, im)
        } else {
            (0u64, 0u64, 0u64, after_strips)
        };

    // ── Open file ─────────────────────────────────────────────────────────────
    let file = std::fs::File::create(path)?;
    let mut w = BufWriter::with_capacity(1 << 23, file); // 8 MiB write buffer

    // ── BigTIFF header (16 bytes) ─────────────────────────────────────────────
    w.write_all(b"II")?;                                   // little-endian marker
    w.write_all(&43u16.to_le_bytes())?;                   // BigTIFF version
    w.write_all(&8u16.to_le_bytes())?;                    // bytesize of offsets (always 8)
    w.write_all(&0u16.to_le_bytes())?;                    // constant = 0
    w.write_all(&BIGTIFF_HEADER_SIZE.to_le_bytes())?;     // offset to first IFD

    // ── IFD count (u64) ───────────────────────────────────────────────────────
    w.write_all(&n_ifd_entries.to_le_bytes())?;

    // ── BigTIFF IFD entry helper (20 bytes per entry) ─────────────────────────
    // For scalar or inline values, `val` holds the data directly (left-
    // justified, zero-padded to 8 bytes in little-endian).  For out-of-line
    // arrays, `val` is the file offset to the data block.
    fn entry64(
        w: &mut impl Write,
        tag: u16, typ: u16, count: u64, val: u64,
    ) -> std::io::Result<()> {
        w.write_all(&tag.to_le_bytes())?;
        w.write_all(&typ.to_le_bytes())?;
        w.write_all(&count.to_le_bytes())?;
        w.write_all(&val.to_le_bytes())?;
        Ok(())
    }

    // RATIONAL 72/1 inline encoding (BigTIFF only).
    // RATIONAL = [u32 numerator, u32 denominator] in LE byte order = 8 bytes.
    // In a u64 LE: (denominator as u64) << 32 | (numerator as u64).
    const RATIONAL_72_1: u64 = (1u64 << 32) | 72u64;

    // ── Base TIFF tags (ascending numeric order) ──────────────────────────────
    entry64(&mut w, 256,  4, 1,        cols as u64)?;             // ImageWidth (LONG)
    entry64(&mut w, 257,  4, 1,        rows as u64)?;             // ImageLength (LONG)
    entry64(&mut w, 258,  3, 1,        bits_per_sample as u64)?;  // BitsPerSample (SHORT)
    entry64(&mut w, 259,  3, 1,        1)?;                       // Compression = None
    entry64(&mut w, 262,  3, 1,        1)?;                       // PhotometricInterp = MinIsBlack

    // StripOffsets (LONG8 = type 16).
    // n_strips == 1: inline value is the offset to the image data.
    // n_strips > 1:  value is offset to the strip-offsets array.
    let strip_offsets_ifd_val = if strips_inline { image_data_off } else { strip_offsets_arr };
    entry64(&mut w, 273, 16, n_strips, strip_offsets_ifd_val)?;  // StripOffsets

    entry64(&mut w, 277,  3, 1,        1)?;                       // SamplesPerPixel = 1
    entry64(&mut w, 278,  4, 1,        ROWS_PER_STRIP)?;          // RowsPerStrip (LONG)

    // StripByteCounts (LONG8 = type 16).
    // n_strips == 1: inline value is the strip's byte count.
    // n_strips > 1:  value is offset to the strip-bytecounts array.
    let strip_bytecounts_ifd_val = if strips_inline {
        last_strip_bytes // == full_strip_bytes when n_strips == 1
    } else {
        strip_bytecounts_arr
    };
    entry64(&mut w, 279, 16, n_strips, strip_bytecounts_ifd_val)?; // StripByteCounts

    // XResolution and YResolution: RATIONAL (type 5), count=1, 8 bytes → inline.
    entry64(&mut w, 282,  5, 1, RATIONAL_72_1)?;                   // XResolution
    entry64(&mut w, 283,  5, 1, RATIONAL_72_1)?;                   // YResolution

    entry64(&mut w, 296,  3, 1, 1)?;                               // ResolutionUnit = no abs
    entry64(&mut w, 339,  3, 1, sample_format as u64)?;            // SampleFormat

    // ── GeoTIFF extension tags (33550, 33922, 34735) ─────────────────────────
    if let Some(keys) = geo_keys {
        // ModelPixelScaleTag (33550): 3 DOUBLE values at pixel_scale_off.
        entry64(&mut w, 33550, 12, 3,               pixel_scale_off)?;
        // ModelTiepointTag (33922): 6 DOUBLE values at tiepoint_off.
        entry64(&mut w, 33922, 12, 6,               tiepoint_off)?;
        // GeoKeyDirectoryTag (34735): N SHORT values at geo_key_dir_off.
        entry64(&mut w, 34735,  3, keys.len() as u64, geo_key_dir_off)?;
    }

    // ── Next IFD offset (u64 = 0 → single IFD) ───────────────────────────────
    w.write_all(&0u64.to_le_bytes())?;

    // ── Strip arrays (only when n_strips > 1) ─────────────────────────────────
    if !strips_inline {
        for s in 0..n_strips {
            let off = image_data_off + s * full_strip_bytes;
            w.write_all(&off.to_le_bytes())?;
        }
        for s in 0..n_strips {
            let bc = if s == n_strips - 1 { last_strip_bytes } else { full_strip_bytes };
            w.write_all(&bc.to_le_bytes())?;
        }
    }

    // ── GeoTIFF extension data (only when geo_keys supplied) ─────────────────
    if let Some(keys) = geo_keys {
        // ModelPixelScaleTag: [pixel_width, pixel_height, 0.0]
        w.write_all(&geotransform[1].to_le_bytes())?;   // X spacing
        w.write_all(&(-geotransform[5]).to_le_bytes())?; // Y spacing magnitude
        w.write_all(&0f64.to_le_bytes())?;

        // ModelTiepointTag: [I=0, J=0, K=0, X=x_origin, Y=y_origin, Z=0]
        w.write_all(&0f64.to_le_bytes())?;
        w.write_all(&0f64.to_le_bytes())?;
        w.write_all(&0f64.to_le_bytes())?;
        w.write_all(&geotransform[0].to_le_bytes())?;   // X of upper-left
        w.write_all(&geotransform[3].to_le_bytes())?;   // Y of upper-left
        w.write_all(&0f64.to_le_bytes())?;

        // GeoKeyDirectoryTag: caller-supplied SHORT block.
        for k in keys {
            w.write_all(&k.to_le_bytes())?;
        }
    }

    // ── Image data ────────────────────────────────────────────────────────────
    w.write_all(bytes)?;

    Ok(())
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    // ── to_db_inplace tests ───────────────────────────────────────────────────

    #[test]
    fn test_to_db_known_value() {
        // 10·log10(1.0) = 0 dB exactly.
        // 10·log10(10.0) = 10 dB exactly.
        let mut data = vec![1.0f32, 10.0f32];
        let (conv, masked) = to_db_inplace(&mut data, 0.0).expect("should succeed");
        assert_eq!(conv, 2);
        assert_eq!(masked, 0);
        assert!((data[0] - 0.0f32).abs() < 1e-5, "1.0 linear → 0 dB, got {}", data[0]);
        assert!((data[1] - 10.0f32).abs() < 1e-5, "10.0 linear → 10 dB, got {}", data[1]);
    }

    #[test]
    fn test_to_db_masks_zero_and_negative() {
        let mut data = vec![0.0f32, -0.001f32, 0.5f32];
        let (conv, masked) = to_db_inplace(&mut data, 0.0).expect("should succeed");
        assert_eq!(conv, 1);
        assert_eq!(masked, 2);
        assert!(data[0].is_nan(), "0.0 should be masked to NaN");
        assert!(data[1].is_nan(), "negative should be masked to NaN");
        assert!(data[2].is_finite(), "0.5 should survive");
    }

    #[test]
    fn test_to_db_noise_floor_masks_additional_pixels() {
        // noise_floor = 0.01 (≈ -20 dB) should mask 0.005 but not 0.05.
        let mut data = vec![0.005f32, 0.05f32];
        let (conv, masked) = to_db_inplace(&mut data, 0.01).expect("should succeed");
        assert_eq!(masked, 1, "0.005 is below noise floor 0.01");
        assert_eq!(conv, 1);
        assert!(data[0].is_nan());
        assert!(data[1].is_finite());
    }

    #[test]
    fn test_to_db_preserves_existing_nan() {
        let mut data = vec![f32::NAN, 1.0f32];
        let (conv, masked) = to_db_inplace(&mut data, 0.0).expect("should succeed");
        assert_eq!(conv, 1);
        assert_eq!(masked, 0); // NaN does not count as masked here
        assert!(data[0].is_nan());
    }

    #[test]
    fn test_to_db_rejects_negative_noise_floor() {
        let mut data = vec![1.0f32];
        let result = to_db_inplace(&mut data, -0.001);
        assert!(result.is_err(), "negative noise floor should be rejected");
    }

    // ── write_geotiff tests ───────────────────────────────────────────────────

    #[test]
    fn test_write_geotiff_rejects_invalid_geotransform() {
        let data = vec![1.0f32; 4];
        // Pixel width = 0 is invalid.
        let bad_gt = [8.0, 0.0, 0.0, 51.0, 0.0, -0.0001];
        let r = write_geotiff("/tmp/sardine_test_bad_gt.tiff", &data, 2, 2, bad_gt);
        assert!(r.is_err(), "zero pixel width should be rejected");

        // Positive gt[5] means south-up (not supported).
        let bad_gt2 = [8.0, 0.0001, 0.0, 51.0, 0.0, 0.0001];
        let r2 = write_geotiff("/tmp/sardine_test_bad_gt2.tiff", &data, 2, 2, bad_gt2);
        assert!(r2.is_err(), "positive gt[5] (south-up) should be rejected");
    }

    #[test]
    fn test_write_geotiff_creates_valid_file() {
        use std::io::Read;
        let data: Vec<f32> = vec![1.0, 2.0, 3.0, 4.0];
        let gt = [8.0, 0.0001, 0.0, 52.0, 0.0, -0.0001];
        let path = "/tmp/sardine_test_geotiff.tiff";
        write_geotiff(path, &data, 2, 2, gt).expect("write should succeed");

        let mut buf = Vec::new();
        std::fs::File::open(path).expect("file should exist").read_to_end(&mut buf).unwrap();

        // TIFF magic: "II" (0x49 0x49) + 42 (0x2a 0x00)
        assert_eq!(&buf[0..4], &[0x49, 0x49, 0x2a, 0x00], "TIFF little-endian magic");

        // Check the GeoKey version word (first u16 of geo_key_dir) = 1.
        // For a 2×2 image: n_strips = ceil(2/64) = 1
        // image_data_off = 326 + 1*8 = 334
        // geo_key_dir_off = 294 + 1*8 = 302
        let n_strips: usize = 1;
        let geo_key_dir_off = 294 + n_strips * 8;
        let version = u16::from_le_bytes([buf[geo_key_dir_off], buf[geo_key_dir_off + 1]]);
        assert_eq!(version, 1, "GeoKey directory version should be 1");

        // Check GeographicTypeGeoKey = 4326 (at offset geo_key_dir_off + 28 within SHORT array)
        // Header is 4 shorts (8 bytes), then key 0 is 4 shorts (8 bytes), key 1 is 4 shorts (8 bytes),
        // key 2 value_offset is at position header(8) + key0(8) + key1(8) + 3*2 = 30 bytes into the array.
        let epsg_offset = geo_key_dir_off + 8 + 8 + 8 + 6; // skip header + 2 keys + 3 fields of key 2
        let epsg_val = u16::from_le_bytes([buf[epsg_offset], buf[epsg_offset + 1]]);
        assert_eq!(epsg_val, 4326, "GeographicTypeGeoKey should be 4326 (WGS84)");
    }

    #[test]
    fn test_write_geotiff_large_auto_selects_bigtiff() {
        use std::io::Read;
        // 33k×33k × 4 B ≈ 4.06 GiB — must select BigTIFF.
        assert!(
            needs_bigtiff(33_000, 33_000, 4, 16),
            "33k×33k f32 raster must trigger BigTIFF path"
        );
        // 1×1 must NOT select BigTIFF.
        assert!(
            !needs_bigtiff(1, 1, 4, 16),
            "1×1 f32 raster must NOT trigger BigTIFF"
        );

        // Exercise the BigTIFF writer directly with a small 2×2 GeoTIFF
        // and verify the BigTIFF magic bytes + data round-trip.
        let data: Vec<f32> = vec![1.0, 2.0, 3.0, 4.0];
        let gt = [8.0, 0.0001, 0.0, 52.0, 0.0, -0.0001];
        let path = "/tmp/sardine_test_bigtiff_geotiff.tiff";
        // SAFETY-OK: f32 → u8 reinterpretation on LE host; same invariant as production code.
        let raw: &[u8] = unsafe {
            std::slice::from_raw_parts(data.as_ptr() as *const u8, data.len() * 4)
        };
        let geo_keys = crate::output_crs::OutputCrs::Wgs84LatLon.geo_key_directory();
        write_bigtiff_raw_inner(path, raw, 2, 2, gt, 32, 3, Some(&geo_keys))
            .expect("BigTIFF write should succeed");

        let mut buf = Vec::new();
        std::fs::File::open(path).unwrap().read_to_end(&mut buf).unwrap();

        // BigTIFF magic: II + 43 (u16) + 8 (u16) + 0 (u16) + 16 (u64)
        assert_eq!(&buf[0..2], b"II", "BigTIFF must use II byte-order marker");
        assert_eq!(
            u16::from_le_bytes([buf[2], buf[3]]),
            43,
            "BigTIFF version field must be 43"
        );
        assert_eq!(
            u16::from_le_bytes([buf[4], buf[5]]),
            8,
            "BigTIFF offset-bytesize field must be 8"
        );
        assert_eq!(
            u64::from_le_bytes(buf[8..16].try_into().unwrap()),
            16,
            "BigTIFF first IFD must start at offset 16"
        );

        // Verify f32 data round-trip.  For a 2×2 GeoTIFF BigTIFF (n_strips=1,
        // inline strips), image_data_off = ifd_end + pixel_scale(24) +
        // tiepoint(48) + geo_keys(32 for WGS84 = 16 shorts × 2).
        // IFD: 8(count) + 16×20(entries) + 8(next) = 336 bytes starting at 16
        // ⇒ ifd_end = 352.  after_strips = 352 (inline).  geo_keys_bytes = 32.
        // image_data_off = 352 + 24 + 48 + 32 = 456.
        let img_off: usize = 456;
        for (i, &expected) in data.iter().enumerate() {
            let off = img_off + i * 4;
            let v = f32::from_le_bytes(buf[off..off + 4].try_into().unwrap());
            assert_eq!(v, expected, "f32 sample {i} must round-trip through BigTIFF");
        }
    }

    // ── write_geotiff_u8 tests ────────────────────────────────────────────────

    /// The u8 writer must produce a valid TIFF with `BitsPerSample = 8` and
    /// `SampleFormat = 1` (unsigned int) — the values that GDAL/QGIS use to
    /// distinguish a quality-mask raster from a float backscatter raster.
    #[test]
    fn test_write_geotiff_u8_creates_valid_file_with_correct_tags() {
        use std::io::Read;
        let data: Vec<u8> = vec![0, 1, 5, 8]; // four mask codes
        let gt = [8.0, 0.0001, 0.0, 52.0, 0.0, -0.0001];
        let path = "/tmp/sardine_test_geotiff_u8.tiff";
        write_geotiff_u8(path, &data, 2, 2, gt).expect("u8 write should succeed");

        let mut buf = Vec::new();
        std::fs::File::open(path).expect("file should exist").read_to_end(&mut buf).unwrap();

        // TIFF magic
        assert_eq!(&buf[0..4], &[0x49, 0x49, 0x2a, 0x00], "TIFF little-endian magic");

        // The IFD starts at offset 8.  Its entries are 12 bytes each,
        // beginning at offset 10 (after the 2-byte entry count).  The
        // BitsPerSample tag (258) is the third entry; tag value lives in the
        // last 4 bytes of its 12-byte entry.
        // Entry 0 = ImageWidth (256), entry 1 = ImageLength (257),
        // entry 2 = BitsPerSample (258).  Offset of entry 2 = 10 + 2*12 = 34.
        // Tag value field = entry start + 8 = 42.
        let bits_per_sample = u32::from_le_bytes([buf[42], buf[43], buf[44], buf[45]]);
        assert_eq!(bits_per_sample, 8, "BitsPerSample tag must encode 8");

        // SampleFormat (339) is the 13th entry (index 12, 0-based) in the
        // 16-entry IFD: offset = 10 + 12*12 = 154.  Value field = 154 + 8 = 162.
        let sample_format = u32::from_le_bytes([buf[162], buf[163], buf[164], buf[165]]);
        assert_eq!(sample_format, 1, "SampleFormat tag must encode 1 (unsigned int)");

        // The four image bytes appear contiguously starting at image_data_off
        // = 326 + n_strips*8 = 326 + 8 = 334 for our 2-row image.
        assert_eq!(&buf[334..338], &[0, 1, 5, 8], "u8 raster bytes must round-trip");
    }

    /// The oversize-image check must use the supplied `bytes_per_pixel`, so
    /// a u8 raster four times the area of the f32 size limit is the one
    /// that trips the check, not a 33k × 33k u8 image (~1 GiB, well under
    /// The BigTIFF threshold respects `bytes_per_pixel`: a u8 raster must be
    /// 4× larger (in pixels) than the f32 threshold to trigger the BigTIFF path.
    #[test]
    fn test_write_geotiff_u8_bigtiff_threshold_uses_byte_pixel() {
        // 66k×66k × 1 B ≈ 4.05 GiB → needs BigTIFF.
        assert!(
            needs_bigtiff(66_000, 66_000, 1, 16),
            "66k×66k u8 raster should require BigTIFF"
        );
        // 33k×33k × 1 B ≈ 1 GiB → does NOT need BigTIFF (u8 is 4× smaller than f32).
        assert!(
            !needs_bigtiff(33_000, 33_000, 1, 16),
            "33k×33k u8 raster (~1 GiB) must NOT require BigTIFF"
        );
    }

    // ── write_tiff_no_crs tests ───────────────────────────────────────────────

    /// The no-CRS writer must produce a valid TIFF whose IFD contains
    /// exactly 13 base entries and **no** GeoTIFF extension tags
    /// (33550 ModelPixelScale, 33922 ModelTiepoint, 34735 GeoKeyDirectory).
    /// This is the explicit guarantee that radar-geometry rasters do not
    /// silently advertise a fabricated EPSG:4326 CRS.
    #[test]
    fn test_write_tiff_no_crs_omits_geotiff_tags_and_round_trips_data() {
        use std::io::Read;
        // Two rows × three columns of distinct values for byte-pattern check.
        let data: Vec<f32> = vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0];
        let path = "/tmp/sardine_test_tiff_no_crs.tiff";
        write_tiff_no_crs(path, &data, 3, 2).expect("no-CRS write should succeed");

        let mut buf = Vec::new();
        std::fs::File::open(path).expect("file should exist").read_to_end(&mut buf).unwrap();

        // TIFF II magic.
        assert_eq!(&buf[0..4], &[0x49, 0x49, 0x2a, 0x00], "TIFF little-endian magic");

        // IFD starts at offset 8.  Entry count is u16 → 13 (no GeoTIFF tags).
        let n_entries = u16::from_le_bytes([buf[8], buf[9]]);
        assert_eq!(n_entries, 13, "no-CRS IFD must have exactly 13 entries");

        // Each entry is 12 bytes; first 2 bytes = tag id.  Walk all entries
        // and assert none of them is a GeoTIFF tag.
        for i in 0..(n_entries as usize) {
            let off = 10 + i * 12;
            let tag = u16::from_le_bytes([buf[off], buf[off + 1]]);
            assert!(
                tag != 33550 && tag != 33922 && tag != 34735,
                "no-CRS file unexpectedly contains GeoTIFF tag {tag} at IFD entry {i}"
            );
        }

        // Image data offset for the no-CRS layout:
        //   ifd_size = 2 + 13*12 + 4 = 162
        //   strip_offsets_arr = 8 + 162 = 170
        //   strip_bytecounts_arr = 170 + n_strips*4 = 170 + 4 = 174
        //   xres_off = 174 + 4 = 178
        //   yres_off = 178 + 8 = 186
        //   image_data_off = yres_off + 8 = 194
        // (n_strips == 1 because rows=2 ≤ ROWS_PER_STRIP=64.)
        let img_off = 194usize;
        for (i, &expected) in data.iter().enumerate() {
            let off = img_off + i * 4;
            let v = f32::from_le_bytes([buf[off], buf[off + 1], buf[off + 2], buf[off + 3]]);
            assert_eq!(v, expected, "f32 sample {i} must round-trip");
        }
    }

    // ── write_geotiff_with_crs / OutputCrs integration tests ─────────────────

    /// Verify the WGS84-on-default path is byte-identical between the legacy
    /// [`write_geotiff`] entry point and the new [`write_geotiff_with_crs`]
    /// when called with `OutputCrs::Wgs84LatLon`.  This is the explicit
    /// backward-compat guarantee for callers (e.g. the terrain-correction
    /// pipeline and the e2e tests) that have not been migrated to the
    /// CRS-aware API.
    #[test]
    fn test_write_geotiff_with_crs_wgs84_byte_identical_to_legacy() {
        use crate::output_crs::OutputCrs;
        use std::io::Read;
        let data: Vec<f32> = (0..(7 * 5)).map(|i| i as f32 * 0.5).collect();
        let gt = [8.0, 1e-4, 0.0, 47.0, 0.0, -1e-4];
        let p_legacy = "/tmp/sardine_export_wgs84_legacy.tiff";
        let p_new    = "/tmp/sardine_export_wgs84_new.tiff";
        write_geotiff(p_legacy, &data, 5, 7, gt).expect("legacy WGS84 write");
        write_geotiff_with_crs(p_new, &data, 5, 7, gt, &OutputCrs::Wgs84LatLon)
            .expect("new WGS84 write");

        let mut a = Vec::new();
        let mut b = Vec::new();
        std::fs::File::open(p_legacy).unwrap().read_to_end(&mut a).unwrap();
        std::fs::File::open(p_new).unwrap().read_to_end(&mut b).unwrap();
        assert_eq!(a, b, "WGS84 path must be byte-identical via legacy or _with_crs");
    }

    /// Writing through `write_geotiff_with_crs` with a UTM CRS must produce
    /// a GeoKeyDirectoryTag with 4 keys (header + GTModelTypeGeoKey=Projected,
    /// GTRasterTypeGeoKey=PixelIsArea, ProjectedCSTypeGeoKey=EPSG, and
    /// ProjLinearUnitsGeoKey=Linear_Meter).  This is the minimum that GDAL
    /// requires to recognise a projected CRS.
    #[test]
    fn test_write_geotiff_with_crs_utm_emits_projected_geokeys() {
        use crate::output_crs::OutputCrs;
        use std::io::Read;
        let data: Vec<f32> = vec![0.0; 4 * 4];
        // Plausible UTM 32N geotransform: 10 m pixel, near Munich.
        let gt = [691_600.0, 10.0, 0.0, 5_334_780.0, 0.0, -10.0];
        let path = "/tmp/sardine_export_utm32n.tiff";
        let crs = OutputCrs::UtmNorth { zone: 32 };
        write_geotiff_with_crs(path, &data, 4, 4, gt, &crs).expect("UTM write");

        let mut buf = Vec::new();
        std::fs::File::open(path).unwrap().read_to_end(&mut buf).unwrap();

        // Walk the IFD to find the GeoKeyDirectoryTag (34735) entry.  The
        // tag's count field tells us the number of SHORTs, and its
        // value-offset field points at the SHORT block in the file.
        let n_entries = u16::from_le_bytes([buf[8], buf[9]]) as usize;
        let mut found = None;
        for i in 0..n_entries {
            let off = 10 + i * 12;
            let tag   = u16::from_le_bytes([buf[off],   buf[off + 1]]);
            let count = u32::from_le_bytes([buf[off + 4], buf[off + 5], buf[off + 6], buf[off + 7]]);
            let val   = u32::from_le_bytes([buf[off + 8], buf[off + 9], buf[off + 10], buf[off + 11]]);
            if tag == 34735 {
                found = Some((count, val));
                break;
            }
        }
        let (count, val_off) = found.expect("UTM TIFF must contain GeoKeyDirectoryTag (34735)");
        // 4-key UTM directory = 4 header SHORTs + 4 keys × 4 SHORTs = 20 SHORTs.
        assert_eq!(count, 20, "UTM GeoKey block must contain exactly 20 SHORTs");

        // Decode the 20 SHORTs from the value-offset block.
        let mut keys = Vec::with_capacity(20);
        for i in 0..(count as usize) {
            let off = val_off as usize + i * 2;
            keys.push(u16::from_le_bytes([buf[off], buf[off + 1]]));
        }
        // Header: [version=1, revision_major=1, revision_minor=0, n_keys=4].
        assert_eq!(&keys[0..4], &[1, 1, 0, 4], "GeoKey header for 4 keys");
        // Each key entry is [KeyID, TIFFTagLocation=0, Count=1, Value].
        // Walk the four 4-SHORT key entries and verify each.
        // GTModelTypeGeoKey (1024) = 1 (Projected).
        assert_eq!(&keys[4..8],   &[1024, 0, 1, 1],     "GTModelType must be Projected");
        // GTRasterTypeGeoKey (1025) = 1 (RasterPixelIsArea).
        assert_eq!(&keys[8..12],  &[1025, 0, 1, 1],     "GTRasterType must be PixelIsArea");
        // ProjectedCSTypeGeoKey (3072) = 32632 (UTM 32N WGS84).
        assert_eq!(&keys[12..16], &[3072, 0, 1, 32632], "ProjectedCSType must be EPSG:32632");
        // ProjLinearUnitsGeoKey (3076) = 9001 (Linear_Meter).
        assert_eq!(&keys[16..20], &[3076, 0, 1, 9001],  "ProjLinearUnits must be metre");
    }

    // ── write_cog_with_crs tests ──────────────────────────────────────────────

    #[test]
    fn test_write_cog_rejects_invalid_tile_size() {
        let data = vec![1.0f32; 4];
        let gt = [8.0, 0.0001, 0.0, 52.0, 0.0, -0.0001];
        // Not a power of two.
        let r = write_cog_with_crs("/tmp/sardine_cog_bad_ts.tiff", &data, 2, 2, gt, &OutputCrs::Wgs84LatLon, 100);
        assert!(matches!(r, Err(ExportError::InvalidTileSize(100))));
        // Too small.
        let r2 = write_cog_with_crs("/tmp/sardine_cog_bad_ts2.tiff", &data, 2, 2, gt, &OutputCrs::Wgs84LatLon, 32);
        assert!(matches!(r2, Err(ExportError::InvalidTileSize(32))));
        // Too large.
        let r3 = write_cog_with_crs("/tmp/sardine_cog_bad_ts3.tiff", &data, 2, 2, gt, &OutputCrs::Wgs84LatLon, 8192);
        assert!(matches!(r3, Err(ExportError::InvalidTileSize(8192))));
    }

    #[test]
    fn test_write_cog_rejects_invalid_geotransform() {
        let data = vec![1.0f32; 4];
        let bad_gt = [8.0, 0.0, 0.0, 52.0, 0.0, -0.0001]; // zero pixel_w
        let r = write_cog_with_crs("/tmp/sardine_cog_bad_gt.tiff", &data, 2, 2, bad_gt, &OutputCrs::Wgs84LatLon, 512);
        assert!(matches!(r, Err(ExportError::InvalidGeotransform(_))));
    }

    #[test]
    fn test_write_cog_creates_valid_tiff_header() {
        use std::io::Read;
        // 1024×512 raster, tile_size=512 → 2×1 full-res tiles, one overview level.
        let rows = 512usize;
        let cols = 1024usize;
        let data: Vec<f32> = (0..rows * cols).map(|i| i as f32 * 0.001).collect();
        let gt = [8.0, 0.0001, 0.0, 52.0, 0.0, -0.0001];
        let path = "/tmp/sardine_cog_test.tiff";
        write_cog_with_crs(path, &data, cols, rows, gt, &OutputCrs::Wgs84LatLon, 512)
            .expect("COG write should succeed");

        let mut buf = Vec::new();
        std::fs::File::open(path).unwrap().read_to_end(&mut buf).unwrap();

        // Classic TIFF header: "II" + 42
        assert_eq!(&buf[0..4], &[0x49, 0x49, 0x2a, 0x00], "must be little-endian TIFF");
        // IFD offset is at bytes 4-7; must be non-zero and within the file.
        let ifd_off = u32::from_le_bytes([buf[4], buf[5], buf[6], buf[7]]) as usize;
        assert!(ifd_off > 0 && ifd_off < buf.len(), "IFD offset must be within file");
        // IFD entry count at ifd_off must be 16.
        let n_entries = u16::from_le_bytes([buf[ifd_off], buf[ifd_off + 1]]) as usize;
        assert_eq!(n_entries, 17, "each COG IFD must have 17 entries (14 base + 3 GeoTIFF)");
    }

    #[test]
    fn test_write_cog_tile_data_round_trips() {
        // Write a small COG with known values; read back the full-res tile and
        // verify pixel values survive the round-trip through tile packing.
        use std::io::Read;
        let rows = 64usize;
        let cols = 64usize;
        // Ramp: pixel value = row * cols + col (cast to f32).
        let data: Vec<f32> = (0..rows * cols).map(|i| i as f32).collect();
        let gt = [0.0, 0.001, 0.0, 1.0, 0.0, -0.001];
        let path = "/tmp/sardine_cog_ramp.tiff";
        // tile_size = 64 → exactly one tile, no overview (64 ≤ tile_size).
        write_cog_with_crs(path, &data, cols, rows, gt, &OutputCrs::Wgs84LatLon, 64)
            .expect("COG write should succeed");

        let mut buf = Vec::new();
        std::fs::File::open(path).unwrap().read_to_end(&mut buf).unwrap();

        // With a 64×64 raster and tile_size=64:
        //   - one overview IFD (32×32), one full-res IFD.
        //   - overview tile is one 64×64-padded tile (upper 32×32 filled, rest NaN).
        //   - full-res tile is the whole 64×64 image.
        // Locate TileOffsets tag (324) in the FULL-RES IFD (second IFD, found via
        // next-IFD pointer from the first).
        let ifd0_off = u32::from_le_bytes([buf[4], buf[5], buf[6], buf[7]]) as usize;
        let n0 = u16::from_le_bytes([buf[ifd0_off], buf[ifd0_off + 1]]) as usize;
        // next-IFD pointer is after the last entry of IFD0.
        let ifd1_off = u32::from_le_bytes({
            let p = ifd0_off + 2 + n0 * 12;
            [buf[p], buf[p+1], buf[p+2], buf[p+3]]
        }) as usize;
        assert!(ifd1_off > 0 && ifd1_off < buf.len(), "second IFD must exist");

        // Scan IFD1 for TileOffsets (tag 324).
        let n1 = u16::from_le_bytes([buf[ifd1_off], buf[ifd1_off + 1]]) as usize;
        let mut tile_data_off: Option<usize> = None;
        for i in 0..n1 {
            let e = ifd1_off + 2 + i * 12;
            let tag = u16::from_le_bytes([buf[e], buf[e+1]]);
            if tag == 324 {
                // Value is a u32 offset (count=1 tile).
                let val = u32::from_le_bytes([buf[e+8], buf[e+9], buf[e+10], buf[e+11]]) as usize;
                tile_data_off = Some(val);
                break;
            }
        }
        let off = tile_data_off.expect("TileOffsets tag 324 must be present");
        // The value for count=1 is an offset to a 1-element u32 array,
        // which in turn contains the actual tile image data offset.
        let tile_image_off = u32::from_le_bytes([buf[off], buf[off+1], buf[off+2], buf[off+3]]) as usize;
        // Read 64×64 × 4 bytes from the tile.
        let tile_bytes = &buf[tile_image_off .. tile_image_off + 64*64*4];
        let tile_f32: Vec<f32> = tile_bytes.chunks_exact(4)
            .map(|b| f32::from_le_bytes([b[0], b[1], b[2], b[3]]))
            .collect();
        // First pixel must be 0.0 (row 0, col 0).
        assert!((tile_f32[0] - 0.0f32).abs() < 1e-6, "first pixel must be 0.0, got {}", tile_f32[0]);
        // Pixel at row 1, col 0 must be 64.0 (= 1 * 64 + 0).
        assert!((tile_f32[64] - 64.0f32).abs() < 1e-6, "pixel[64] must be 64.0, got {}", tile_f32[64]);
        // Last pixel must be (64*64 - 1) = 4095.0.
        assert!((tile_f32[64*64 - 1] - 4095.0f32).abs() < 1e-6,
            "last pixel must be 4095.0, got {}", tile_f32[64*64-1]);
    }

    #[test]
    fn test_downsample_2x_basic() {
        // 4×4 ramp, verify that 2×2 average produces the expected 2×2 result.
        let src: Vec<f32> = (0..16).map(|i| i as f32).collect();
        // Row 0: [0, 1, 2, 3]   Row 1: [4, 5, 6, 7]
        // Row 2: [8, 9,10,11]   Row 3: [12,13,14,15]
        //
        // Top-left 2×2 avg:  (0+1+4+5)/4 = 2.5
        // Top-right 2×2 avg: (2+3+6+7)/4 = 4.5
        // Bot-left 2×2 avg:  (8+9+12+13)/4 = 10.5
        // Bot-right 2×2 avg: (10+11+14+15)/4 = 12.5
        let (dst, dst_cols, dst_rows) = downsample_2x(&src, 4, 4);
        assert_eq!((dst_cols, dst_rows), (2, 2));
        assert!((dst[0] - 2.5f32).abs() < 1e-6, "top-left avg = 2.5, got {}", dst[0]);
        assert!((dst[1] - 4.5f32).abs() < 1e-6, "top-right avg = 4.5, got {}", dst[1]);
        assert!((dst[2] - 10.5f32).abs() < 1e-6, "bot-left avg = 10.5, got {}", dst[2]);
        assert!((dst[3] - 12.5f32).abs() < 1e-6, "bot-right avg = 12.5, got {}", dst[3]);
    }

    #[test]
    fn test_downsample_2x_nan_propagation() {
        // If all 4 pixels in a patch are NaN, output is NaN.
        // If only 1 is NaN, the other 3 still produce a finite mean.
        let mut src = vec![f32::NAN; 4];
        src[0] = 1.0;
        let (dst, _, _) = downsample_2x(&src, 2, 2);
        // Only pixel [0,0]=1.0 is valid → mean of 1 = 1.0.
        assert!((dst[0] - 1.0f32).abs() < 1e-6, "partial NaN patch: mean of valid pixels");

        let all_nan = vec![f32::NAN; 4];
        let (dst2, _, _) = downsample_2x(&all_nan, 2, 2);
        assert!(dst2[0].is_nan(), "all-NaN patch → NaN output");
    }
}
