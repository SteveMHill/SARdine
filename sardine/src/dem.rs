//! DEM tile loading and mosaic: SRTM 1" and Copernicus GLO-30.
//!
//! # Supported formats
//!
//! | Format | Extension | Dimensions | Height ref | Notes |
//! |--------|-----------|------------|------------|-------|
//! | SRTM 1" | `.hgt` | 3601 × 3601 (with 1-pixel border overlap) | EGM96 | Big-endian i16 |
//! | Copernicus GLO-30 | `.tif` / `.tiff` | 3600 × 3600 (no overlap) | EGM2008 | GeoTIFF f32 or i16 |
//!
//! [`DemMosaic::load_directory`] auto-detects format by file extension.
//! Use [`DemMosaic::load_srtm_files`] or [`DemMosaic::load_glo30_files`] for
//! explicit control.
//!
//! # Height reference note
//!
//! Both SRTM (EGM96) and GLO-30 (EGM2008) store **orthometric heights** — above
//! the geoid, not the WGS84 ellipsoid.  Range-Doppler geocoding requires
//! ellipsoidal heights.  Apply [`crate::geoid::GeoidModel`] undulation correction
//! **after** reading from this module.  The EGM96/EGM2008 difference is < 4 m
//! globally, so EGM96 correction applied to GLO-30 data is acceptable for most
//! applications.
//!
//! # Void handling
//!
//! SRTM voids (raw value −32768) and GLO-30 nodata pixels are both propagated as
//! `f32::NAN`.  The caller (terrain correction inner loop) skips NaN pixels and
//! counts them as `dem_missing`, which is reported in the output statistics.

use std::cell::Cell;
use std::collections::HashMap;
use std::io::BufReader;
use std::path::Path;

// ── Constants ──────────────────────────────────────────────────────

/// SRTM 1-arc-second tile dimension.  Both rows and columns (includes 1-pixel border overlap).
const SRTM1_DIM: usize = 3601;

/// SRTM void marker value in the raw i16 data.
const SRTM_VOID: i16 = -32768;

/// GeoTIFF tag: ModelPixelScale — pixel size in (scale_x, scale_y, scale_z) world units.
const TAG_MODEL_PIXEL_SCALE: u16 = 33550;

/// GeoTIFF tag: ModelTiepoint — maps a raster (I, J, K) to world (X, Y, Z).
const TAG_MODEL_TIEPOINT: u16 = 33922;

/// GDAL metadata tag: NODATA — ASCII-encoded nodata value.
const TAG_GDAL_NODATA: u16 = 42113;

// ── Error type ─────────────────────────────────────────────────────

#[derive(Debug, thiserror::Error)]
pub enum DemError {
    #[error("IO error reading {path}: {source}")]
    Io {
        path: String,
        source: std::io::Error,
    },

    #[error("tile {path} has wrong size: expected {expected} bytes, got {actual}")]
    BadSize {
        path: String,
        expected: usize,
        actual: usize,
    },

    #[error("TIFF decode error reading {path}: {detail}")]
    TiffDecode {
        path: String,
        detail: String,
    },

    #[error("GLO-30 tile {path} has unsupported pixel type (expected f32 or i16)")]
    UnsupportedPixelType { path: String },

    #[error("cannot determine tile origin for '{name}': filename does not match SRTM (N50E010) or GLO-30 (Copernicus_DSM_COG_10_N50_00_E010_00_DEM) naming conventions")]
    UnrecognisedTileName { name: String },

    #[error("no tiles loaded")]
    NoTiles,

    #[error("query point ({lat:.4}, {lon:.4}) is outside DEM coverage")]
    OutOfBounds { lat: f64, lon: f64 },

    #[error(
        "DEM mosaic does not cover the requested area: \
         requested lat=[{req_min_lat:.4}, {req_max_lat:.4}], \
         lon=[{req_min_lon:.4}, {req_max_lon:.4}]; \
         mosaic lat=[{have_min_lat:.4}, {have_max_lat:.4}], \
         lon=[{have_min_lon:.4}, {have_max_lon:.4}]; \
         missing margin (deg): south={miss_south:.4}, north={miss_north:.4}, \
         west={miss_west:.4}, east={miss_east:.4}"
    )]
    CoverageGap {
        req_min_lat: f64,
        req_max_lat: f64,
        req_min_lon: f64,
        req_max_lon: f64,
        have_min_lat: f64,
        have_max_lat: f64,
        have_min_lon: f64,
        have_max_lon: f64,
        miss_south: f64,
        miss_north: f64,
        miss_west: f64,
        miss_east: f64,
    },
}

// ── SRTM tile ──────────────────────────────────────────────────────

/// A single SRTM 1" tile covering 1° × 1°.
///
/// Height reference: **EGM96** (orthometric).  Apply [`crate::geoid::GeoidModel`]
/// to convert to ellipsoidal height before ECEF conversion.
struct SrtmTile {
    /// SW corner latitude (integer degrees).
    lat_sw: i32,
    /// SW corner longitude (integer degrees).
    lon_sw: i32,
    /// Elevation data, row-major, north-to-south, west-to-east.
    /// Row 0 = north edge (lat_sw + 1°), row 3600 = south edge (lat_sw).
    /// SRTM tiles have a 1-pixel border overlap: both edges are stored,
    /// giving 3601 × 3601 total.
    data: Vec<f32>,
}

impl SrtmTile {
    /// Load a tile from an `.hgt` file.
    ///
    /// The filename must encode the SW corner, e.g. `N50E010.hgt`.
    fn load(path: &Path) -> Result<Self, DemError> {
        let stem = path
            .file_stem()
            .and_then(|s| s.to_str())
            .unwrap_or(""); // SAFETY-OK: empty stem → parse_srtm_tile_name returns UnrecognisedTileName

        let (lat_sw, lon_sw) = parse_srtm_tile_name(stem)?;

        let raw = std::fs::read(path).map_err(|e| DemError::Io {
            path: path.display().to_string(),
            source: e,
        })?;

        let expected = SRTM1_DIM * SRTM1_DIM * 2;
        if raw.len() != expected {
            return Err(DemError::BadSize {
                path: path.display().to_string(),
                expected,
                actual: raw.len(),
            });
        }

        let data: Vec<f32> = raw
            .chunks_exact(2)
            .map(|chunk| {
                let val = i16::from_be_bytes([chunk[0], chunk[1]]);
                if val == SRTM_VOID {
                    // Propagate as NaN; callers skip NaN and count as dem_missing.
                    f32::NAN
                } else {
                    val as f32
                }
            })
            .collect();

        Ok(SrtmTile { lat_sw, lon_sw, data })
    }

    /// Bilinear interpolation at (lat_deg, lon_deg).
    ///
    /// Returns `None` if the query is outside this tile's coverage.
    fn elevation_at(&self, lat_deg: f64, lon_deg: f64) -> Option<f32> {
        let lat_min = self.lat_sw as f64;
        let lon_min = self.lon_sw as f64;
        let frac_lat = lat_deg - lat_min; // 0 = south edge, 1 = north edge
        let frac_lon = lon_deg - lon_min;

        if frac_lat < 0.0 || frac_lat > 1.0 || frac_lon < 0.0 || frac_lon > 1.0 {
            return None;
        }

        let n = SRTM1_DIM - 1; // 3600
        // Row 0 = north edge → invert lat fraction.
        let row_f = (1.0 - frac_lat) * n as f64;
        let col_f = frac_lon * n as f64;
        bilinear_on_grid(&self.data, SRTM1_DIM, row_f, col_f, n)
    }
}

// ── GLO-30 tile ────────────────────────────────────────────────────

/// A single Copernicus GLO-30 tile covering 1° × 1°.
///
/// # Height reference
///
/// GLO-30 uses **EGM2008** (orthometric).  Apply [`crate::geoid::GeoidModel`]
/// to convert to ellipsoidal height.  The EGM2008/EGM96 difference is < 4 m
/// globally, so an EGM96 geoid model applied to GLO-30 data is acceptable.
///
/// # Tile origin
///
/// The top-left (NW) pixel **centre** is at `(origin_lat_deg, origin_lon_deg)`.
/// For a standard 3600 × 3600 tile aligned to a 1° boundary, the origin is
/// approximately one half-pixel inside the NW corner of the 1° cell.
///
/// # Naming conventions
///
/// The loader accepts two filename patterns (case-insensitive lat/lon prefix):
///
/// ```text
/// Copernicus_DSM_COG_10_N50_00_E010_00_DEM.tif   (official ESA/Copernicus)
/// N50E010.tif                                      (flattened, SRTM-style)
/// ```
struct Glo30Tile {
    /// Latitude of the NW pixel centre (degrees).
    origin_lat_deg: f64,
    /// Longitude of the NW pixel centre (degrees).
    origin_lon_deg: f64,
    /// Pixel size in degrees (assumed square; 1 arc-second = 1/3600°).
    pixel_deg: f64,
    /// Number of columns.
    width: usize,
    /// Number of rows.
    height: usize,
    /// Elevation data, row-major, north-to-south, west-to-east.
    /// NaN marks nodata pixels.
    data: Vec<f32>,
}

impl Glo30Tile {
    /// Load a GLO-30 tile from a GeoTIFF file.
    ///
    /// Reads georeferencing from `ModelPixelScaleTag` / `ModelTiepointTag` when
    /// present; otherwise derives the origin from the filename.
    fn load(path: &Path) -> Result<Self, DemError> {
        use tiff::decoder::{Decoder, DecodingResult};
        use tiff::tags::Tag;

        let path_str = path.display().to_string();

        let file = std::fs::File::open(path).map_err(|e| DemError::Io {
            path: path_str.clone(),
            source: e,
        })?;
        let mut decoder = Decoder::new(BufReader::new(file)).map_err(|e| DemError::TiffDecode {
            path: path_str.clone(),
            detail: e.to_string(),
        })?;

        let (width_u32, height_u32) = decoder.dimensions().map_err(|e| DemError::TiffDecode {
            path: path_str.clone(),
            detail: e.to_string(),
        })?;
        let width = width_u32 as usize;
        let height = height_u32 as usize;

        // ── Read GeoTIFF metadata tags (optional; fall back to filename) ──────
        // ModelPixelScaleTag (33550): [scale_x, scale_y, scale_z]
        // ModelTiepointTag (33922): [I, J, K, X, Y, Z] — raster → world mapping
        // I and J are typically 0.0 (pixel corner) or 0.5 (pixel centre).
        // For GLO-30, convention is pixel-corner, so the top-left pixel centre is
        // at (X + 0.5*scale_x, Y - 0.5*scale_y).
        let (origin_lat_deg, origin_lon_deg, pixel_deg) = {
            let scale = decoder.get_tag_f64_vec(Tag::Unknown(TAG_MODEL_PIXEL_SCALE)).ok();
            let tiepoint = decoder.get_tag_f64_vec(Tag::Unknown(TAG_MODEL_TIEPOINT)).ok();

            match (scale, tiepoint) {
                (Some(s), Some(tp)) if s.len() >= 2 && tp.len() >= 6 => {
                    let scale_x = s[0]; // degrees per pixel (longitude)
                    let scale_y = s[1]; // degrees per pixel (latitude, positive)
                    let raster_i = tp[0]; // pixel column of tiepoint
                    let raster_j = tp[1]; // pixel row of tiepoint
                    let world_x = tp[3]; // longitude of tiepoint
                    let world_y = tp[4]; // latitude of tiepoint

                    // Convert tiepoint to the centre of the top-left pixel (row 0, col 0):
                    //   lon_centre = world_x + (0.5 - raster_i) * scale_x
                    //   lat_centre = world_y - (0.5 - raster_j) * scale_y
                    let lon0 = world_x + (0.5 - raster_i) * scale_x;
                    let lat0 = world_y - (0.5 - raster_j) * scale_y;
                    // Use the mean of x and y scale for the pixel size; they should be equal.
                    let px = (scale_x + scale_y) / 2.0;
                    (lat0, lon0, px)
                }
                _ => {
                    // No GeoTIFF tags — derive origin from the filename.
                    let (lat_sw, lon_sw) = parse_glo30_tile_name(path)?;
                    // Standard GLO-30: 3600 pixels over 1°, pixel-centre of row 0 is
                    // one half-pixel below the north edge.
                    let px = 1.0 / width as f64; // ≈ ARC_SECOND_DEG for 3600-wide tiles
                    let lat0 = (lat_sw + 1) as f64 - 0.5 * px;
                    let lon0 = lon_sw as f64 + 0.5 * px;
                    (lat0, lon0, px)
                }
            }
        };

        // ── Read nodata value ─────────────────────────────────────────────────
        let nodata: Option<f32> = decoder
            .get_tag_ascii_string(Tag::Unknown(TAG_GDAL_NODATA))
            .ok()
            .and_then(|s| s.trim().parse::<f32>().ok());

        // ── Read pixel data ───────────────────────────────────────────────────
        let image = decoder.read_image().map_err(|e| DemError::TiffDecode {
            path: path_str.clone(),
            detail: e.to_string(),
        })?;

        let data: Vec<f32> = match image {
            DecodingResult::F32(v) => {
                v.into_iter().map(|x| {
                    if !x.is_finite() { return f32::NAN; }
                    if let Some(nd) = nodata { if x == nd { return f32::NAN; } }
                    x
                }).collect()
            }
            DecodingResult::I16(v) => {
                v.into_iter().map(|x| {
                    if let Some(nd) = nodata { if x == nd as i16 { return f32::NAN; } }
                    x as f32
                }).collect()
            }
            DecodingResult::U16(v) => {
                // Some GLO-30 reprojections use u16
                v.into_iter().map(|x| {
                    if let Some(nd) = nodata { if x == nd as u16 { return f32::NAN; } }
                    x as f32
                }).collect()
            }
            _ => {
                return Err(DemError::UnsupportedPixelType { path: path_str });
            }
        };

        Ok(Glo30Tile {
            origin_lat_deg,
            origin_lon_deg,
            pixel_deg,
            width,
            height,
            data,
        })
    }

    /// Bilinear interpolation at (lat_deg, lon_deg).
    ///
    /// Returns `None` if the query is outside this tile's coverage.
    fn elevation_at(&self, lat_deg: f64, lon_deg: f64) -> Option<f32> {
        // Convert to fractional pixel coordinates (0 = centre of first pixel).
        // Row 0 is north; lat decreases going south.
        let row_f = (self.origin_lat_deg - lat_deg) / self.pixel_deg;
        let col_f = (lon_deg - self.origin_lon_deg) / self.pixel_deg;

        let max_row = (self.height - 1) as f64;
        let max_col = (self.width - 1) as f64;

        // Allow a half-pixel margin beyond the pixel-centre boundary so that
        // queries at the exact tile edge (aligned to the outer grid line) work.
        let margin = 0.5;
        if row_f < -margin || row_f > max_row + margin
            || col_f < -margin || col_f > max_col + margin {
            return None;
        }

        let row_f = row_f.clamp(0.0, max_row);
        let col_f = col_f.clamp(0.0, max_col);
        bilinear_on_grid(&self.data, self.width, row_f, col_f, self.height - 1)
    }
}

// ── Shared bilinear helper ─────────────────────────────────────────

/// Bilinear interpolation on a flat row-major grid.
///
/// `dim` is the stride (number of columns).
/// `max_rc` is the maximum valid row/col index (= dim - 1 for square grids,
/// or height - 1 for rectangular ones).
/// `row_f` and `col_f` are the fractional pixel coordinates (0-based).
///
/// Returns `None` if any of the four surrounding pixels is NaN.
#[inline]
fn bilinear_on_grid(data: &[f32], dim: usize, row_f: f64, col_f: f64, max_rc: usize) -> Option<f32> {
    let r0 = (row_f as usize).min(max_rc.saturating_sub(1));
    let c0 = (col_f as usize).min(dim.saturating_sub(2));
    let r1 = r0 + 1;
    let c1 = c0 + 1;

    let dr = row_f - r0 as f64;
    let dc = col_f - c0 as f64;

    let v00 = data[r0 * dim + c0];
    let v01 = data[r0 * dim + c1];
    let v10 = data[r1 * dim + c0];
    let v11 = data[r1 * dim + c1];

    // Any NaN propagates.
    if v00.is_nan() || v01.is_nan() || v10.is_nan() || v11.is_nan() {
        return None;
    }

    let v = v00 as f64 * (1.0 - dr) * (1.0 - dc)
        + v01 as f64 * (1.0 - dr) * dc
        + v10 as f64 * dr * (1.0 - dc)
        + v11 as f64 * dr * dc;
    Some(v as f32)
}

// ── Tile name parsers ──────────────────────────────────────────────

/// Parse SRTM filename like `N50E010` → `(lat_sw=50, lon_sw=10)`.
fn parse_srtm_tile_name(name: &str) -> Result<(i32, i32), DemError> {
    // Must be at least "N50E010" = 7 chars.
    if name.len() < 7 {
        return Err(DemError::UnrecognisedTileName { name: name.to_string() });
    }
    let bytes = name.as_bytes();
    let lat_sign: i32 = match bytes[0] {
        b'N' | b'n' => 1,
        b'S' | b's' => -1,
        _ => return Err(DemError::UnrecognisedTileName { name: name.to_string() }),
    };
    let lat: i32 = name[1..3].parse().map_err(|_| DemError::UnrecognisedTileName {
        name: name.to_string(),
    })?;
    let lon_sign: i32 = match bytes[3] {
        b'E' | b'e' => 1,
        b'W' | b'w' => -1,
        _ => return Err(DemError::UnrecognisedTileName { name: name.to_string() }),
    };
    let lon: i32 = name[4..7].parse().map_err(|_| DemError::UnrecognisedTileName {
        name: name.to_string(),
    })?;
    Ok((lat * lat_sign, lon * lon_sign))
}

/// Parse a GLO-30 tile path to `(lat_sw, lon_sw)`.
///
/// Accepts two filename conventions:
///
/// 1. **Official Copernicus**: `Copernicus_DSM_COG_10_N50_00_E010_00_DEM.tif`
/// 2. **SRTM-style flat**: `N50E010.tif`
fn parse_glo30_tile_name(path: &Path) -> Result<(i32, i32), DemError> {
    let stem = path
        .file_stem()
        .and_then(|s| s.to_str())
        .unwrap_or(""); // SAFETY-OK: empty stem → UnrecognisedTileName below

    // Try SRTM-style first (short stem like "N50E010").
    if stem.len() == 7 {
        if let Ok(pair) = parse_srtm_tile_name(stem) {
            return Ok(pair);
        }
    }

    // Try Copernicus naming: "Copernicus_DSM_COG_10_N50_00_E010_00_DEM"
    // Split by '_' and look for a token starting with N/S followed by digits,
    // then a token starting with E/W followed by digits.
    let parts: Vec<&str> = stem.split('_').collect();
    let mut lat_val: Option<i32> = None;
    let mut lon_val: Option<i32> = None;
    for (i, &part) in parts.iter().enumerate() {
        if lat_val.is_none() {
            let bytes = part.as_bytes();
            if bytes.len() >= 3 && (bytes[0] == b'N' || bytes[0] == b'S') {
                if let Ok(n) = part[1..].parse::<i32>() {
                    let sign = if bytes[0] == b'N' { 1 } else { -1 };
                    // Skip sub-degree part (the "00" token after lat)
                    // by checking the next token is "00"
                    if matches!(parts.get(i + 1), Some(&"00")) {
                        lat_val = Some(sign * n);
                    }
                }
            }
        } else if lon_val.is_none() {
            let bytes = part.as_bytes();
            if bytes.len() >= 3 && (bytes[0] == b'E' || bytes[0] == b'W') {
                if let Ok(n) = part[1..].parse::<i32>() {
                    let sign = if bytes[0] == b'E' { 1 } else { -1 };
                    if matches!(parts.get(i + 1), Some(&"00")) {
                        lon_val = Some(sign * n);
                    }
                }
            }
        }
        if lat_val.is_some() && lon_val.is_some() {
            break;
        }
    }

    match (lat_val, lon_val) {
        (Some(lat), Some(lon)) => Ok((lat, lon)),
        _ => Err(DemError::UnrecognisedTileName { name: stem.to_string() }),
    }
}

// ── DEM tile enum ──────────────────────────────────────────────────

/// A single DEM tile — either SRTM 1" or Copernicus GLO-30.
enum DemTile {
    Srtm(SrtmTile),
    Glo30(Glo30Tile),
}

impl DemTile {
    fn elevation_at(&self, lat_deg: f64, lon_deg: f64) -> Option<f32> {
        match self {
            DemTile::Srtm(t) => t.elevation_at(lat_deg, lon_deg),
            DemTile::Glo30(t) => t.elevation_at(lat_deg, lon_deg),
        }
    }

    fn lat_sw(&self) -> f64 {
        match self {
            DemTile::Srtm(t) => t.lat_sw as f64,
            // GLO-30 origin is NW pixel centre; lat_sw ≈ origin - 1° + half_pixel
            DemTile::Glo30(t) => t.origin_lat_deg - (t.height as f64 - 0.5) * t.pixel_deg,
        }
    }

    fn lat_ne(&self) -> f64 {
        match self {
            DemTile::Srtm(t) => (t.lat_sw + 1) as f64,
            DemTile::Glo30(t) => t.origin_lat_deg + 0.5 * t.pixel_deg,
        }
    }

    fn lon_sw(&self) -> f64 {
        match self {
            DemTile::Srtm(t) => t.lon_sw as f64,
            DemTile::Glo30(t) => t.origin_lon_deg - 0.5 * t.pixel_deg,
        }
    }

    fn lon_ne(&self) -> f64 {
        match self {
            DemTile::Srtm(t) => (t.lon_sw + 1) as f64,
            DemTile::Glo30(t) => t.origin_lon_deg + (t.width as f64 - 0.5) * t.pixel_deg,
        }
    }
}

// ── DEM Mosaic ─────────────────────────────────────────────────────

/// A DEM mosaic assembled from SRTM 1" and/or Copernicus GLO-30 tiles.
///
/// Provides bilinear-interpolated elevation at any point within the mosaic's
/// combined coverage.
///
/// # Loading
///
/// | Method | What it loads |
/// |--------|---------------|
/// | [`load_directory`](Self::load_directory) | Auto-detects by extension (`.hgt` → SRTM, `.tif`/`.tiff` → GLO-30) |
/// | [`load_srtm_files`](Self::load_srtm_files) | Explicit SRTM `.hgt` files |
/// | [`load_glo30_files`](Self::load_glo30_files) | Explicit GLO-30 `.tif` files |
pub struct DemMosaic {
    tiles: Vec<DemTile>,
    /// Spatial index over the tile set: maps an integer (lat_floor,
    /// lon_floor) 1°×1° cell to the list of tile indices whose bounding
    /// box overlaps that cell.
    ///
    /// SRTM 1" tiles align to integer corners by construction, so each
    /// SRTM tile lands in exactly one bin.  Copernicus GLO-30 tiles centre
    /// their pixels half a pixel inside integer corners, so a single
    /// GLO-30 tile may land in up to four adjacent bins; this is correct
    /// (the search still terminates at the first tile that returns
    /// `Some(_)`) but slightly wastes lookup work — acceptable.
    ///
    /// Replaces the previous O(N_tiles) linear scan in `elevation_at`.
    /// The bin lookup is O(1) amortised; bin contents are typically a
    /// single tile, so per-pixel lookup cost collapses from
    /// `~30× SRTM bbox checks` to `~1× SRTM bbox check` on a typical
    /// Sentinel-1 scene mosaic.
    index: HashMap<(i32, i32), Vec<u32>>,
}

// Per-worker last-hit tile cache.  Terrain correction processes pixels in
// raster order within a row; consecutive pixels almost always hit the same
// tile, so a single-slot cache turns the bin lookup into a bbox check on
// the warm path.  `Cell` is `!Sync` but `thread_local!` storage is
// per-thread by definition, so this is sound across the rayon row
// closures.
thread_local! {
    static LAST_HIT_TILE: Cell<Option<u32>> = Cell::new(None);
}

impl DemMosaic {
    /// Load all DEM tiles from a directory.
    ///
    /// Auto-detects format:
    /// - `.hgt` → SRTM 1-arc-second
    /// - `.tif` / `.tiff` → Copernicus GLO-30 GeoTIFF
    ///
    /// Files with other extensions are silently ignored.
    /// Returns [`DemError::NoTiles`] if no supported files are found.
    pub fn load_directory(dir: &Path) -> Result<Self, DemError> {
        let mut tiles = Vec::new();
        let entries = std::fs::read_dir(dir).map_err(|e| DemError::Io {
            path: dir.display().to_string(),
            source: e,
        })?;

        let mut paths: Vec<_> = entries
            .filter_map(|e| e.ok().map(|e| e.path()))
            .collect();
        // Sort for deterministic ordering (makes tests reproducible).
        paths.sort();

        for path in paths {
            match path.extension().and_then(|s| s.to_str()) {
                Some("hgt") => tiles.push(DemTile::Srtm(SrtmTile::load(&path)?)),
                Some("tif") | Some("tiff") => tiles.push(DemTile::Glo30(Glo30Tile::load(&path)?)),
                _ => {}
            }
        }

        if tiles.is_empty() {
            return Err(DemError::NoTiles);
        }
        let index = Self::build_index(&tiles);
        Ok(DemMosaic { tiles, index })
    }

    /// Load SRTM 1" tiles from explicit `.hgt` file paths.
    pub fn load_srtm_files(paths: &[&Path]) -> Result<Self, DemError> {
        let tiles: Result<Vec<_>, _> = paths
            .iter()
            .map(|p| SrtmTile::load(p).map(DemTile::Srtm))
            .collect();
        let tiles = tiles?;
        if tiles.is_empty() {
            return Err(DemError::NoTiles);
        }
        let index = Self::build_index(&tiles);
        Ok(DemMosaic { tiles, index })
    }

    /// Load Copernicus GLO-30 tiles from explicit GeoTIFF file paths.
    pub fn load_glo30_files(paths: &[&Path]) -> Result<Self, DemError> {
        let tiles: Result<Vec<_>, _> = paths
            .iter()
            .map(|p| Glo30Tile::load(p).map(DemTile::Glo30))
            .collect();
        let tiles = tiles?;
        if tiles.is_empty() {
            return Err(DemError::NoTiles);
        }
        let index = Self::build_index(&tiles);
        Ok(DemMosaic { tiles, index })
    }

    /// Build the spatial index over a set of tiles.  See the
    /// [`DemMosaic::index`] field doc for binning semantics.  A small
    /// epsilon (`1e-9°` ≈ ~0.1 mm) absorbs the GLO-30 "half-pixel
    /// outside the integer corner" pattern so that integer-aligned
    /// SRTM tiles land in a single bin instead of two.
    fn build_index(tiles: &[DemTile]) -> HashMap<(i32, i32), Vec<u32>> {
        const BIN_EPS: f64 = 1e-9;
        let mut idx: HashMap<(i32, i32), Vec<u32>> = HashMap::new();
        for (i, t) in tiles.iter().enumerate() {
            let j_min = t.lat_sw().floor() as i32;
            let j_max = (t.lat_ne() - BIN_EPS).floor() as i32;
            let k_min = t.lon_sw().floor() as i32;
            let k_max = (t.lon_ne() - BIN_EPS).floor() as i32;
            for j in j_min..=j_max {
                for k in k_min..=k_max {
                    idx.entry((j, k)).or_default().push(i as u32);
                }
            }
        }
        idx
    }

    /// Number of tiles in the mosaic.
    pub fn tile_count(&self) -> usize {
        self.tiles.len()
    }

    /// Bilinear-interpolated elevation (metres) at `(lat_deg, lon_deg)`.
    ///
    /// Lookup uses an O(1) integer-cell hash index (see
    /// [`DemMosaic::index`]) plus a per-thread last-hit cache to
    /// exploit raster-order locality in the terrain-correction inner
    /// loop.  For correctness an explicit linear-scan fallback is
    /// available via [`Self::elevation_at_linear`] for diagnostics
    /// and equivalence tests.
    ///
    /// Returns [`DemError::OutOfBounds`] if no tile covers the point.
    /// Returns `Ok(f32::NAN)` if the point is within a tile's coverage but
    /// the pixel is a void/nodata cell.
    pub fn elevation_at(&self, lat_deg: f64, lon_deg: f64) -> Result<f32, DemError> {
        // Normalise longitude to [−180, 180) so that anti-meridian wrapped
        // bounding boxes (max_lon_deg > 180) resolve to standard tile coords.
        // Use rem_euclid to also handle values < −180.
        let lon_norm = ((lon_deg + 180.0).rem_euclid(360.0)) - 180.0;

        // Warm path: try the last tile that satisfied this worker.
        // `tiles[i]` indexing is provably in-bounds because `index` is
        // built once at construction time and `tiles` is never mutated.
        let last = LAST_HIT_TILE.with(|c| c.get());
        if let Some(li) = last {
            if let Some(h) = self.tiles[li as usize].elevation_at(lat_deg, lon_norm) {
                return Ok(h);
            }
        }

        // Cold path: hash bin then per-bin scan.
        let key = (lat_deg.floor() as i32, lon_norm.floor() as i32);
        if let Some(bin) = self.index.get(&key) {
            for &ti in bin {
                if Some(ti) == last {
                    continue;
                }
                if let Some(h) = self.tiles[ti as usize].elevation_at(lat_deg, lon_norm) {
                    LAST_HIT_TILE.with(|c| c.set(Some(ti)));
                    return Ok(h);
                }
            }
        }
        Err(DemError::OutOfBounds { lat: lat_deg, lon: lon_deg })
    }

    /// Linear-scan elevation lookup, retained as a reference
    /// implementation for equivalence testing against the indexed
    /// [`Self::elevation_at`] fast path.  Not intended for production
    /// use; quadratic in the number of tiles per pixel.
    #[doc(hidden)]
    pub fn elevation_at_linear(&self, lat_deg: f64, lon_deg: f64) -> Result<f32, DemError> {
        let lon_norm = ((lon_deg + 180.0).rem_euclid(360.0)) - 180.0;
        for tile in &self.tiles {
            if let Some(h) = tile.elevation_at(lat_deg, lon_norm) {
                return Ok(h);
            }
        }
        Err(DemError::OutOfBounds { lat: lat_deg, lon: lon_deg })
    }

    /// Bounding box of the mosaic: `(min_lat, max_lat, min_lon, max_lon)`.
    pub fn bounds(&self) -> (f64, f64, f64, f64) {
        let mut min_lat = f64::MAX;
        let mut max_lat = f64::MIN;
        let mut min_lon = f64::MAX;
        let mut max_lon = f64::MIN;
        for tile in &self.tiles {
            min_lat = min_lat.min(tile.lat_sw());
            max_lat = max_lat.max(tile.lat_ne());
            min_lon = min_lon.min(tile.lon_sw());
            max_lon = max_lon.max(tile.lon_ne());
        }
        (min_lat, max_lat, min_lon, max_lon)
    }

    /// Verify that this mosaic covers the requested geographic bounding box,
    /// optionally inflated by `margin_deg` in every direction.
    ///
    /// Returns [`DemError::CoverageGap`] if any side of the (inflated) request
    /// box is outside the mosaic's combined extent.  This is intended as a
    /// pre-flight check before terrain correction so that missing tiles fail
    /// loudly with a precise error rather than silently producing NaN-filled
    /// output.
    ///
    /// `margin_deg` should be ≥ the maximum geocoding excursion expected from
    /// terrain (typically a few hundredths of a degree at S-1 incidence
    /// angles).  The check uses tile-corner extents, not per-pixel coverage,
    /// so a passing result does not preclude internal voids — those still
    /// surface as NaN per the `dem_missing` counter.
    ///
    /// Note: the mosaic-bounds computation does not currently handle the
    /// 180° anti-meridian; bounding boxes that straddle it must be split by
    /// the caller before invoking this check.
    pub fn covers_bbox(
        &self,
        min_lat: f64,
        max_lat: f64,
        min_lon: f64,
        max_lon: f64,
        margin_deg: f64,
    ) -> Result<(), DemError> {
        if self.tiles.is_empty() {
            return Err(DemError::NoTiles);
        }
        let (have_min_lat, have_max_lat, have_min_lon, have_max_lon) = self.bounds();
        let req_min_lat = min_lat - margin_deg;
        let req_max_lat = max_lat + margin_deg;
        let req_min_lon = min_lon - margin_deg;
        let req_max_lon = max_lon + margin_deg;

        let miss_south = (have_min_lat - req_min_lat).max(0.0);
        let miss_north = (req_max_lat - have_max_lat).max(0.0);
        let miss_west = (have_min_lon - req_min_lon).max(0.0);
        let miss_east = (req_max_lon - have_max_lon).max(0.0);

        if miss_south > 0.0 || miss_north > 0.0 || miss_west > 0.0 || miss_east > 0.0 {
            return Err(DemError::CoverageGap {
                req_min_lat,
                req_max_lat,
                req_min_lon,
                req_max_lon,
                have_min_lat,
                have_max_lat,
                have_min_lon,
                have_max_lon,
                miss_south,
                miss_north,
                miss_west,
                miss_east,
            });
        }
        Ok(())
    }
}

// ── DemSource trait ────────────────────────────────────────────────

/// Trait for any object that can provide elevation lookups.
///
/// Implement this trait to plug a custom DEM backend into
/// [`crate::terrain_correction::terrain_correction`] without modifying the
/// geocoding code.  [`DemMosaic`] implements this trait; callers that load
/// SRTM or GLO-30 tiles via [`DemMosaic::load_directory`] obtain a
/// `DemSource` automatically.
///
/// # Thread safety
///
/// The terrain correction inner loop is parallelised with Rayon.  All
/// `DemSource` implementations must therefore be `Send + Sync`.
pub trait DemSource: Send + Sync {
    /// Return the bilinear-interpolated elevation (metres) at `(lat_deg, lon_deg)`.
    ///
    /// # Errors
    ///
    /// - [`DemError::OutOfBounds`] — no tile covers the query point.
    /// - Any other `DemError` variant if reading the tile fails.
    ///
    /// # Nodata
    ///
    /// Return `Ok(f32::NAN)` for void/nodata pixels within a tile's coverage.
    /// The caller (terrain correction) skips NaN and counts it as `dem_missing`.
    fn elevation_at(&self, lat_deg: f64, lon_deg: f64) -> Result<f32, DemError>;
}

impl DemSource for DemMosaic {
    fn elevation_at(&self, lat_deg: f64, lon_deg: f64) -> Result<f32, DemError> {
        DemMosaic::elevation_at(self, lat_deg, lon_deg)
    }
}

// ── Tests ──────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parse_tile_names() {
        assert_eq!(parse_srtm_tile_name("N50E010").unwrap(), (50, 10));
        assert_eq!(parse_srtm_tile_name("S01W071").unwrap(), (-1, -71));
        assert_eq!(parse_srtm_tile_name("N49E008").unwrap(), (49, 8));
    }

    #[test]
    fn parse_glo30_names() {
        // Copernicus official naming
        let p = Path::new("Copernicus_DSM_COG_10_N50_00_E010_00_DEM.tif");
        assert_eq!(parse_glo30_tile_name(p).unwrap(), (50, 10));

        let p = Path::new("Copernicus_DSM_COG_10_S01_00_W071_00_DEM.tif");
        assert_eq!(parse_glo30_tile_name(p).unwrap(), (-1, -71));

        // SRTM-style flat name also accepted
        let p = Path::new("N50E010.tif");
        assert_eq!(parse_glo30_tile_name(p).unwrap(), (50, 10));
    }

    #[test]
    fn unrecognised_name_returns_error() {
        assert!(parse_srtm_tile_name("badname").is_err());
        assert!(parse_glo30_tile_name(Path::new("random_file.tif")).is_err());
    }

    // ── GLO-30 synthetic tile unit tests ─────────────────────────────────────

    /// Build a minimal synthetic Glo30Tile for unit tests (no file I/O).
    fn make_synthetic_glo30(lat_sw: i32, lon_sw: i32, fill: f32) -> Glo30Tile {
        // 4×4 pixel synthetic tile at 0.25° spacing (not real GLO-30 resolution,
        // but exercises the same bilinear interpolation path).
        let width = 4usize;
        let height = 4usize;
        let pixel_deg = 0.25;
        // NW pixel centre is one half-pixel inside the north edge.
        let origin_lat_deg = (lat_sw + 1) as f64 - 0.5 * pixel_deg;
        let origin_lon_deg = lon_sw as f64 + 0.5 * pixel_deg;
        Glo30Tile {
            origin_lat_deg,
            origin_lon_deg,
            pixel_deg,
            width,
            height,
            data: vec![fill; width * height],
        }
    }

    #[test]
    fn glo30_flat_tile_returns_fill_value() {
        let tile = make_synthetic_glo30(50, 10, 42.0);
        // Centre of tile
        let h = tile.elevation_at(50.5, 10.5).unwrap();
        assert!((h - 42.0).abs() < 1e-4, "expected 42.0, got {h}");
    }

    #[test]
    fn glo30_outside_tile_returns_none() {
        let tile = make_synthetic_glo30(50, 10, 100.0);
        assert!(tile.elevation_at(49.0, 10.5).is_none()); // south of tile
        assert!(tile.elevation_at(52.0, 10.5).is_none()); // north of tile
        assert!(tile.elevation_at(50.5,  9.0).is_none()); // west of tile
        assert!(tile.elevation_at(50.5, 12.0).is_none()); // east of tile
    }

    #[test]
    fn glo30_nan_pixel_returns_none() {
        let mut tile = make_synthetic_glo30(50, 10, 100.0);
        // Set top-left pixel to NaN — any query that touches it returns None.
        tile.data[0] = f32::NAN;
        // Top-left pixel is row 0, col 0 — query at its centre.
        let h = tile.elevation_at(tile.origin_lat_deg, tile.origin_lon_deg);
        // Bilinear reads 4 neighbours; if any is NaN the result is None.
        assert!(h.is_none(), "NaN pixel should produce None");
    }

    #[test]
    fn covers_bbox_passes_when_request_inside_mosaic() {
        // Two adjacent synthetic GLO-30 tiles cover lat [50, 51], lon [10, 12].
        let mosaic = {
            let tiles = vec![
                DemTile::Glo30(make_synthetic_glo30(50, 10, 0.0)),
                DemTile::Glo30(make_synthetic_glo30(50, 11, 0.0)),
            ];
            let index = DemMosaic::build_index(&tiles);
            DemMosaic { tiles, index }
        };
        // Request fully inside, with a 0.05° margin: must succeed.
        mosaic
            .covers_bbox(50.2, 50.8, 10.2, 11.8, 0.05)
            .expect("request inside mosaic should pass");
    }

    #[test]
    fn covers_bbox_reports_each_missing_side() {
        let mosaic = {
            let tiles = vec![DemTile::Glo30(make_synthetic_glo30(50, 10, 0.0))];
            let index = DemMosaic::build_index(&tiles);
            DemMosaic { tiles, index }
        };
        // Request that pokes 0.1° beyond every edge of the mosaic.
        let res = mosaic.covers_bbox(49.9, 51.1, 9.9, 11.1, 0.0);
        match res {
            Err(DemError::CoverageGap {
                miss_south,
                miss_north,
                miss_west,
                miss_east,
                ..
            }) => {
                assert!((miss_south - 0.1).abs() < 1e-9, "miss_south = {miss_south}");
                assert!((miss_north - 0.1).abs() < 1e-9, "miss_north = {miss_north}");
                assert!((miss_west - 0.1).abs() < 1e-9, "miss_west = {miss_west}");
                assert!((miss_east - 0.1).abs() < 1e-9, "miss_east = {miss_east}");
            }
            other => panic!("expected CoverageGap on all sides, got {other:?}"),
        }
    }

    #[test]
    fn covers_bbox_margin_can_push_request_outside() {
        // Mosaic covers exactly [50, 51] × [10, 11].  Request fits, but a
        // 0.5° margin pushes it outside on every side ⇒ must fail.
        let mosaic = {
            let tiles = vec![DemTile::Glo30(make_synthetic_glo30(50, 10, 0.0))];
            let index = DemMosaic::build_index(&tiles);
            DemMosaic { tiles, index }
        };
        let r = mosaic.covers_bbox(50.1, 50.9, 10.1, 10.9, 0.5);
        assert!(matches!(r, Err(DemError::CoverageGap { .. })));
    }

    #[test]
    fn covers_bbox_empty_mosaic_returns_no_tiles() {
        let mosaic = DemMosaic { tiles: vec![], index: HashMap::new() };
        assert!(matches!(
            mosaic.covers_bbox(0.0, 1.0, 0.0, 1.0, 0.0),
            Err(DemError::NoTiles)
        ));
    }

    #[test]
    fn glo30_bilinear_interpolation_linear_ramp() {
        // Build a tile whose elevation increases linearly from west to east.
        // Elevation at column c = c * 10.0 (for a 4-column tile: 0, 10, 20, 30).
        let width = 4usize;
        let height = 4usize;
        let pixel_deg = 0.25;
        let lat_sw = 50i32;
        let lon_sw = 10i32;
        let origin_lat_deg = (lat_sw + 1) as f64 - 0.5 * pixel_deg;
        let origin_lon_deg = lon_sw as f64 + 0.5 * pixel_deg;
        let data: Vec<f32> = (0..height)
            .flat_map(|_| (0..width).map(|c| c as f32 * 10.0))
            .collect();
        let tile = Glo30Tile { origin_lat_deg, origin_lon_deg, pixel_deg, width, height, data };

        // At the centre of column 0 (lon = 10.125°) the elevation should be 0.
        let h0 = tile.elevation_at(50.5, 10.125).unwrap();
        assert!((h0 - 0.0).abs() < 1e-3, "col 0 centre: {h0}");

        // At the centre of column 1 (lon = 10.375°) the elevation should be 10.
        let h1 = tile.elevation_at(50.5, 10.375).unwrap();
        assert!((h1 - 10.0).abs() < 1e-3, "col 1 centre: {h1}");

        // At the midpoint between col 0 and col 1 (lon = 10.25°) the elevation
        // should interpolate to 5.0.
        let hm = tile.elevation_at(50.5, 10.25).unwrap();
        assert!((hm - 5.0).abs() < 1e-3, "mid col 0–1: {hm}");
    }

    // ── Real data tests ───────────────────────────────────────────────────────

    #[test]
    fn load_real_tile() {
        let path = Path::new("/home/datacube/dev/SARdine/data/dem/srtm1/N50E010.hgt");
        if !path.exists() {
            return; // skip if no data
        }

        let tile = SrtmTile::load(path).unwrap();
        assert_eq!(tile.lat_sw, 50);
        assert_eq!(tile.lon_sw, 10);
        assert_eq!(tile.data.len(), SRTM1_DIM * SRTM1_DIM);

        // Spot check: elevation at centre of tile should be reasonable
        // N50.5, E10.5 is in Thuringia, Germany — elevation ~300-600m
        let h = tile.elevation_at(50.5, 10.5).unwrap();
        assert!(h > 100.0 && h < 1000.0, "unexpected elevation: {}", h);
    }

    #[test]
    fn load_mosaic_directory() {
        let dir = Path::new("/home/datacube/dev/SARdine/data/dem/srtm1");
        if !dir.exists() {
            return;
        }

        let mosaic = DemMosaic::load_directory(dir).unwrap();
        assert_eq!(mosaic.tile_count(), 26);

        // Check bounding box — tiles now cover N47–N51 × E007–E012
        let (min_lat, max_lat, min_lon, max_lon) = mosaic.bounds();
        assert!((min_lat - 47.0).abs() < 0.01);
        assert!((max_lat - 52.0).abs() < 0.01);
        assert!((min_lon - 7.0).abs() < 0.01);
        assert!((max_lon - 13.0).abs() < 0.01);

        // Query elevation at a known point (Kassel, Germany, ~160m)
        let h = mosaic.elevation_at(51.3, 9.5).unwrap();
        assert!(h > 100.0 && h < 400.0, "Kassel elevation: {}", h);
    }

    #[test]
    fn tile_boundary_interpolation() {
        let path = Path::new("/home/datacube/dev/SARdine/data/dem/srtm1/N50E010.hgt");
        if !path.exists() {
            return;
        }

        let tile = SrtmTile::load(path).unwrap();

        // Exactly at SW corner
        let h = tile.elevation_at(50.0, 10.0).unwrap();
        assert!(h.is_finite());

        // Exactly at NE corner
        let h = tile.elevation_at(51.0, 11.0).unwrap();
        assert!(h.is_finite());

        // Outside tile
        assert!(tile.elevation_at(49.9, 10.0).is_none());
        assert!(tile.elevation_at(50.0, 11.1).is_none());
    }

    /// Indexed `elevation_at` must return bit-identical results to the
    /// preserved linear-scan reference implementation, otherwise the
    /// Phase-1 spatial index has changed observable behaviour.  We
    /// sample a deterministic-pseudo-random grid of points across the
    /// real on-disk SRTM mosaic plus a margin of out-of-bounds points
    /// to exercise both `Ok(_)` and `Err(OutOfBounds)` paths.
    #[test]
    fn indexed_lookup_matches_linear_scan_on_real_mosaic() {
        let dir = Path::new("/home/datacube/dev/SARdine/data/dem/srtm1");
        if !dir.exists() {
            return;
        }
        let mosaic = DemMosaic::load_directory(dir).unwrap();
        let (min_lat, max_lat, min_lon, max_lon) = mosaic.bounds();

        // Inflate slightly past the mosaic so a nontrivial fraction of
        // queries lands OOB.
        let lat_lo = min_lat - 0.5;
        let lat_hi = max_lat + 0.5;
        let lon_lo = min_lon - 0.5;
        let lon_hi = max_lon + 0.5;

        // LCG: cheap, deterministic, no extra dep.
        let mut state: u64 = 0x9E3779B97F4A7C15;
        let mut next = || {
            state = state.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
            (state >> 11) as f64 / (1u64 << 53) as f64
        };

        for _ in 0..1024 {
            let lat = lat_lo + next() * (lat_hi - lat_lo);
            let lon = lon_lo + next() * (lon_hi - lon_lo);
            let a = mosaic.elevation_at(lat, lon);
            let b = mosaic.elevation_at_linear(lat, lon);
            match (a, b) {
                (Ok(x), Ok(y)) => {
                    if x.is_nan() {
                        assert!(y.is_nan(), "indexed=NaN but linear={y} at ({lat},{lon})");
                    } else {
                        // Exact equality: same tile resolved, same bilinear weights.
                        assert_eq!(x.to_bits(), y.to_bits(), "mismatch at ({lat},{lon})");
                    }
                }
                (Err(DemError::OutOfBounds { .. }), Err(DemError::OutOfBounds { .. })) => {}
                other => panic!("indexed/linear disagree at ({lat},{lon}): {other:?}"),
            }
        }
    }

    /// The per-thread last-hit cache must produce the same answer on
    /// repeated queries inside a single tile (the raster-order locality
    /// pattern in terrain correction).  Failure here means the cache
    /// is returning a stale tile or short-circuiting incorrectly.
    #[test]
    fn last_hit_cache_returns_correct_value_on_locality_pattern() {
        let dir = Path::new("/home/datacube/dev/SARdine/data/dem/srtm1");
        if !dir.exists() {
            return;
        }
        let mosaic = DemMosaic::load_directory(dir).unwrap();

        // Reset the thread-local cache so this test does not depend on
        // execution order.
        LAST_HIT_TILE.with(|c| c.set(None));

        // Walk a scanline inside one tile, then jump to a different tile,
        // then back, exercising hit-on-cache and cache-update paths.
        let in_a = (50.5, 10.5);
        let in_b = (47.5, 7.5);

        let a0 = mosaic.elevation_at_linear(in_a.0, in_a.1).unwrap();
        let b0 = mosaic.elevation_at_linear(in_b.0, in_b.1).unwrap();

        for _ in 0..5 {
            let a = mosaic.elevation_at(in_a.0, in_a.1).unwrap();
            assert_eq!(a.to_bits(), a0.to_bits());
        }
        for _ in 0..5 {
            let b = mosaic.elevation_at(in_b.0, in_b.1).unwrap();
            assert_eq!(b.to_bits(), b0.to_bits());
        }
        for _ in 0..5 {
            let a = mosaic.elevation_at(in_a.0, in_a.1).unwrap();
            assert_eq!(a.to_bits(), a0.to_bits());
        }
    }
}
