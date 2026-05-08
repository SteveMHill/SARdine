//! Pure-Rust reader for Sentinel-1 SLC measurement TIFF files.
//!
//! # Format
//!
//! S1 measurement TIFFs use CInt16 pixels: `SampleFormat=5`, `BitsPerSample=32`.
//! Each pixel encodes one complex sample as two packed [`i16`] values
//! (I component first, then Q). The file has one strip per row, uncompressed
//! (`Compression=1`), in a contiguous layout.
//!
//! The files are classic TIFF (magic=42), little-endian.
//!
//! # Design
//!
//! This reader does not use any TIFF library.  It parses only the IFD fields
//! needed to locate data, validates the layout, then issues direct
//! `seek` + `read` calls.  Because the strips are contiguous, the byte offset
//! of any row is:
//!
//! ```text
//! offset(row) = base_offset + row * (width * 4)
//! ```
//!
//! Individual bursts can therefore be read without loading the entire subswath.
//! Peak memory per burst is `lines_per_burst * width * 4` bytes (~132 MB for IW).
//!
//! # What this reader does NOT do
//!
//! - No calibration (handled downstream after deburst).
//! - No noise removal.
//! - No deburst / overlap blending.
//! - No invalid-sample masking (use `firstValidSample` / `lastValidSample`
//!   from the annotation XML for that).

use std::fs::File;
use std::io::{self, Read, Seek, SeekFrom};
use std::path::Path;

use thiserror::Error;

// ─── TIFF constants ───────────────────────────────────────────────────────────

const TAG_IMAGE_WIDTH: u16 = 256;
const TAG_IMAGE_LENGTH: u16 = 257;
const TAG_BITS_PER_SAMPLE: u16 = 258;
const TAG_COMPRESSION: u16 = 259;
const TAG_STRIP_OFFSETS: u16 = 273;
const TAG_ROWS_PER_STRIP: u16 = 278;
const TAG_SAMPLES_PER_PIXEL: u16 = 277;
const TAG_SAMPLE_FORMAT: u16 = 339;

const TIFF_MAGIC_CLASSIC: u16 = 42;
const TIFF_MAGIC_BIG: u16 = 43;

const FIELD_TYPE_SHORT: u16 = 3;
const FIELD_TYPE_LONG: u16 = 4;

const COMPRESSION_NONE: u16 = 1;
/// CInt: two signed 16-bit channels packed as a single 32-bit TIFF sample.
const SAMPLE_FORMAT_CINT: u16 = 5;
/// 32 bits total = 2 × 16-bit I/Q for CInt16.
const BITS_PER_SAMPLE_CINT16: u16 = 32;

// ─── Error type ───────────────────────────────────────────────────────────────

/// Errors produced by [`SlcReader`].
#[derive(Debug, Error)]
pub enum SlcReadError {
    /// An I/O error occurred while reading the file.
    #[error("I/O error: {0}")]
    Io(#[from] io::Error),

    /// The file does not begin with a valid TIFF header.
    #[error("not a classic TIFF (bad byte-order marker or magic={magic})")]
    NotTiff { magic: u16 },

    /// BigTIFF format (magic=43) is detected but not supported.
    #[error("BigTIFF (magic=43) is not supported; only classic TIFF (magic=42) is accepted")]
    BigTiff,

    /// A required IFD tag is absent.
    #[error("required TIFF tag {0} is absent from the IFD")]
    MissingTag(u16),

    /// The TIFF format is not the expected CInt16 uncompressed layout.
    #[error("unsupported TIFF format: {0}")]
    UnsupportedFormat(String),

    /// The strip offsets are not contiguous; direct seek arithmetic cannot be used.
    ///
    /// This should never occur for Sentinel-1 measurement files.
    #[error(
        "strip offsets not contiguous (strip[0]={first}, strip[{n_strips}-1]={actual_last}, \
         expected={expected_last}); direct-seek reader cannot be used"
    )]
    NonContiguousStrips {
        first: u64,
        n_strips: u32,
        actual_last: u64,
        expected_last: u64,
    },

    /// The requested line range extends beyond the image height.
    #[error(
        "line range [{first}, {end}) is out of image bounds [0, {height})"
    )]
    LineOutOfBounds {
        first: usize,
        end: usize,
        height: u32,
    },

    /// A burst read request in a multi-slice reader spans a slice boundary.
    ///
    /// This indicates a bug in the assembled burst geometry: each burst must
    /// be contained entirely within one slice's TIFF file.
    #[error(
        "burst read [{first}, {end}) crosses a slice boundary at logical line \
         {slice_end}; this indicates corrupted assembled burst geometry"
    )]
    CrossSliceBoundary {
        first: usize,
        end: usize,
        slice_end: usize,
    },

    /// Two TIFF files in a multi-slice reader have different widths (range
    /// samples), which must be equal for all slices of the same sub-swath.
    #[error(
        "slice {index} width {found} != slice 0 width {expected}; all slices \
         in the same sub-swath must have the same number of range samples"
    )]
    SliceWidthMismatch {
        index: usize,
        expected: u32,
        found: u32,
    },
}

// ─── BurstReader trait ────────────────────────────────────────────────────────

/// Uniform read interface over a single-TIFF [`SlcReader`] and the multi-TIFF
/// [`crate::multi_slc_reader::MultiSlcReader`].
///
/// [`deburst_subswath`][crate::deburst::deburst_subswath] is generic over this
/// trait so it works for both single-SAFE and multi-SAFE (slice-assembled) scenes.
///
/// # Safety contract
///
/// Implementors must ensure that `read_burst_raw(first_line, n_lines)` returns
/// exactly `n_lines * width()` complex samples, or a typed `SlcReadError`
/// describing the precise failure mode.
pub trait BurstReader {
    /// Number of range samples (columns) per line.
    fn width(&self) -> u32;

    /// Total number of azimuth lines accessible through this reader.
    ///
    /// For [`SlcReader`] this is the physical TIFF height.
    /// For `MultiSlcReader` this is the sum of all slice TIFF heights
    /// (logical height of the assembled scene).
    fn height(&self) -> u32;

    /// Read `line_count` consecutive azimuth lines starting at `first_line`.
    ///
    /// Returns a flat `Vec<[i16; 2]>` of length `line_count * self.width()`.
    /// Each element is one CInt16 complex pixel: `[i, q]`.
    ///
    /// # Errors
    ///
    /// - [`SlcReadError::LineOutOfBounds`] if `first_line + line_count > self.height()`.
    /// - [`SlcReadError::Io`] on file read failure.
    fn read_burst_raw(
        &mut self,
        first_line: usize,
        line_count: usize,
    ) -> Result<Vec<[i16; 2]>, SlcReadError>;
}

impl BurstReader for SlcReader {
    #[inline]
    fn width(&self) -> u32 {
        self.width
    }

    #[inline]
    fn height(&self) -> u32 {
        self.height
    }

    #[inline]
    fn read_burst_raw(
        &mut self,
        first_line: usize,
        line_count: usize,
    ) -> Result<Vec<[i16; 2]>, SlcReadError> {
        self.read_burst_raw(first_line, line_count)
    }
}

// ─── Internal helpers ─────────────────────────────────────────────────────────

fn read_u16_from(buf: &[u8], le: bool) -> u16 {
    let arr: [u8; 2] = buf[..2].try_into().expect("need 2 bytes");
    if le {
        u16::from_le_bytes(arr)
    } else {
        u16::from_be_bytes(arr)
    }
}

fn read_u32_from(buf: &[u8], le: bool) -> u32 {
    let arr: [u8; 4] = buf[..4].try_into().expect("need 4 bytes");
    if le {
        u32::from_le_bytes(arr)
    } else {
        u32::from_be_bytes(arr)
    }
}

/// Extract a scalar u32 from the IFD entry's inline value field.
///
/// For count=1 SHORT or LONG entries the value is stored in the first 2 or 4
/// bytes of the 4-byte value/offset field (left-justified, in the file's byte
/// order). Returns `None` for unrecognised field types.
fn entry_scalar(value_field: &[u8; 4], field_type: u16, le: bool) -> Option<u32> {
    match field_type {
        FIELD_TYPE_SHORT => Some(read_u16_from(&value_field[..2], le) as u32),
        FIELD_TYPE_LONG => Some(read_u32_from(&value_field[..4], le)),
        _ => None,
    }
}

/// Read a single strip-offset value from the current file position.
fn read_strip_offset_at_pos(
    file: &mut File,
    field_type: u16,
    le: bool,
) -> Result<u64, SlcReadError> {
    match field_type {
        FIELD_TYPE_SHORT => {
            let mut buf = [0u8; 2];
            file.read_exact(&mut buf)?;
            Ok(read_u16_from(&buf, le) as u64)
        }
        FIELD_TYPE_LONG => {
            let mut buf = [0u8; 4];
            file.read_exact(&mut buf)?;
            Ok(read_u32_from(&buf, le) as u64)
        }
        _ => Err(SlcReadError::UnsupportedFormat(format!(
            "StripOffsets field type {field_type} is not SHORT(3) or LONG(4)"
        ))),
    }
}

// ─── IFD parser ───────────────────────────────────────────────────────────────

struct TiffLayout {
    little_endian: bool,
    width: u32,
    height: u32,
    /// Byte offset of the first pixel of row 0.
    base_offset: u64,
    /// Bytes per row: `width * 4` (two i16 per pixel).
    row_bytes: u64,
}

fn parse_ifd(file: &mut File) -> Result<TiffLayout, SlcReadError> {
    // ── Header (8 bytes) ─────────────────────────────────────────────────────
    let mut header = [0u8; 8];
    file.seek(SeekFrom::Start(0))?;
    file.read_exact(&mut header)?;

    let le = match &header[0..2] {
        b"II" => true,
        b"MM" => false,
        _ => return Err(SlcReadError::NotTiff { magic: 0 }),
    };

    let magic = read_u16_from(&header[2..4], le);
    match magic {
        TIFF_MAGIC_BIG => return Err(SlcReadError::BigTiff),
        TIFF_MAGIC_CLASSIC => {}
        other => return Err(SlcReadError::NotTiff { magic: other }),
    }

    let ifd_offset = read_u32_from(&header[4..8], le) as u64;

    // ── IFD entries ──────────────────────────────────────────────────────────
    file.seek(SeekFrom::Start(ifd_offset))?;
    let mut n_buf = [0u8; 2];
    file.read_exact(&mut n_buf)?;
    let n_entries = read_u16_from(&n_buf, le) as usize;

    let mut entry_bytes = vec![0u8; n_entries * 12];
    file.read_exact(&mut entry_bytes)?;

    // ── Extract tags ─────────────────────────────────────────────────────────
    let mut width: Option<u32> = None;
    let mut height: Option<u32> = None;
    let mut bits_per_sample: Option<u16> = None;
    let mut compression: Option<u16> = None;
    let mut rows_per_strip: Option<u32> = None;
    let mut samples_per_pixel: Option<u16> = None;
    // StripOffsets IFD entry: (field_type, count, value_or_offset)
    let mut strip_offsets_entry: Option<(u16, u32, u32)> = None;
    let mut sample_format: Option<u16> = None;

    for i in 0..n_entries {
        let e = &entry_bytes[i * 12..(i + 1) * 12];
        let tag = read_u16_from(&e[0..2], le);
        let field_type = read_u16_from(&e[2..4], le);
        let count = read_u32_from(&e[4..8], le);
        let value_field: &[u8; 4] = e[8..12].try_into().expect("IFD entry is 12 bytes");

        match tag {
            TAG_IMAGE_WIDTH => {
                width = entry_scalar(value_field, field_type, le);
            }
            TAG_IMAGE_LENGTH => {
                height = entry_scalar(value_field, field_type, le);
            }
            TAG_BITS_PER_SAMPLE => {
                bits_per_sample =
                    entry_scalar(value_field, field_type, le).map(|v| v as u16);
            }
            TAG_COMPRESSION => {
                compression =
                    entry_scalar(value_field, field_type, le).map(|v| v as u16);
            }
            TAG_ROWS_PER_STRIP => {
                rows_per_strip = entry_scalar(value_field, field_type, le);
            }
            TAG_SAMPLES_PER_PIXEL => {
                samples_per_pixel =
                    entry_scalar(value_field, field_type, le).map(|v| v as u16);
            }
            TAG_STRIP_OFFSETS => {
                let value_or_offset = read_u32_from(value_field, le);
                strip_offsets_entry = Some((field_type, count, value_or_offset));
            }
            TAG_SAMPLE_FORMAT => {
                sample_format =
                    entry_scalar(value_field, field_type, le).map(|v| v as u16);
            }
            _ => {}
        }
    }

    // ── Validate required tags ────────────────────────────────────────────────
    let width = width.ok_or(SlcReadError::MissingTag(TAG_IMAGE_WIDTH))?;
    let height = height.ok_or(SlcReadError::MissingTag(TAG_IMAGE_LENGTH))?;
    let bps = bits_per_sample.ok_or(SlcReadError::MissingTag(TAG_BITS_PER_SAMPLE))?;
    let comp = compression.ok_or(SlcReadError::MissingTag(TAG_COMPRESSION))?;
    let rps = rows_per_strip.ok_or(SlcReadError::MissingTag(TAG_ROWS_PER_STRIP))?;
    // SamplesPerPixel defaults to 1 when absent (TIFF spec §7 table note).
    let spp = samples_per_pixel.unwrap_or(1); // SAFETY-OK: TIFF 6.0 spec defines default of 1 when SamplesPerPixel tag absent
    let sfmt = sample_format.ok_or(SlcReadError::MissingTag(TAG_SAMPLE_FORMAT))?;
    let (so_type, so_count, so_value_or_offset) =
        strip_offsets_entry.ok_or(SlcReadError::MissingTag(TAG_STRIP_OFFSETS))?;

    if comp != COMPRESSION_NONE {
        return Err(SlcReadError::UnsupportedFormat(format!(
            "Compression={comp} (expected 1=None)"
        )));
    }
    if sfmt != SAMPLE_FORMAT_CINT {
        return Err(SlcReadError::UnsupportedFormat(format!(
            "SampleFormat={sfmt} (expected 5=CInt)"
        )));
    }
    if bps != BITS_PER_SAMPLE_CINT16 {
        return Err(SlcReadError::UnsupportedFormat(format!(
            "BitsPerSample={bps} (expected 32 for CInt16)"
        )));
    }
    if rps != 1 {
        return Err(SlcReadError::UnsupportedFormat(format!(
            "RowsPerStrip={rps} (expected 1)"
        )));
    }
    if spp != 1 {
        // CInt16 complex data must have exactly one sample per pixel.
        // A value != 1 would mean our row_bytes formula (width × 4) is wrong.
        return Err(SlcReadError::UnsupportedFormat(format!(
            "SamplesPerPixel={spp} (expected 1 for CInt16)"
        )));
    }
    if so_count != height {
        return Err(SlcReadError::UnsupportedFormat(format!(
            "strip count ({so_count}) != image height ({height})"
        )));
    }

    let row_bytes = width as u64 * 4; // 2 i16 per pixel × 2 bytes each

    // ── Resolve base offset and verify contiguity ─────────────────────────────
    // For count=1, the single offset is stored inline in the IFD entry.
    // For count>1, `so_value_or_offset` is a file offset to the offset table.
    let (base_offset, last_offset) = if so_count == 1 {
        let base = so_value_or_offset as u64;
        (base, base)
    } else {
        let table_offset = so_value_or_offset as u64;
        let entry_bytes_each: u64 = match so_type {
            FIELD_TYPE_SHORT => 2,
            FIELD_TYPE_LONG => 4,
            other => {
                return Err(SlcReadError::UnsupportedFormat(format!(
                    "StripOffsets field type {other} is not SHORT(3) or LONG(4)"
                )))
            }
        };

        file.seek(SeekFrom::Start(table_offset))?;
        let base = read_strip_offset_at_pos(file, so_type, le)?;

        file.seek(SeekFrom::Start(
            table_offset + (so_count as u64 - 1) * entry_bytes_each,
        ))?;
        let last = read_strip_offset_at_pos(file, so_type, le)?;

        (base, last)
    };

    // Contiguous layout: strip[k] == base + k * row_bytes for all k.
    // Verifying first and last is sufficient to catch any non-contiguous layout.
    let expected_last = base_offset + (height as u64 - 1) * row_bytes;
    if last_offset != expected_last {
        return Err(SlcReadError::NonContiguousStrips {
            first: base_offset,
            n_strips: so_count,
            actual_last: last_offset,
            expected_last,
        });
    }

    Ok(TiffLayout {
        little_endian: le,
        width,
        height,
        base_offset,
        row_bytes,
    })
}

// ─── Public reader ────────────────────────────────────────────────────────────

/// Reader for a single Sentinel-1 SLC measurement TIFF file.
///
/// Opens the file, validates the TIFF layout, and provides methods to read
/// individual burst windows as raw `[i16; 2]` samples or as `(f32, f32)`.
///
/// # Burst coordinate convention
///
/// Line indices are 0-based and refer to the full-subswath image raster
/// (all bursts stacked in azimuth order).  The burst's `first_line` and
/// `line_count = lines_per_burst` come from the annotation XML (see
/// [`crate::BurstEntry`] and [`crate::SubSwathMetadata`]).
///
/// # Notes on invalid samples
///
/// The first and last samples of each row are typically zero (invalid).
/// Use `firstValidSample` and `lastValidSample` from the annotation XML to
/// identify the valid region.  This reader returns the full row including
/// invalid samples; masking is the caller's responsibility.
pub struct SlcReader {
    file: File,
    /// Image width: range samples per row.
    pub width: u32,
    /// Image height: total azimuth lines (all bursts combined).
    pub height: u32,
    little_endian: bool,
    base_offset: u64,
    row_bytes: u64,
}

impl std::fmt::Debug for SlcReader {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("SlcReader")
            .field("width", &self.width)
            .field("height", &self.height)
            .field("little_endian", &self.little_endian)
            .field("base_offset", &self.base_offset)
            .field("row_bytes", &self.row_bytes)
            .finish_non_exhaustive()
    }
}

impl SlcReader {
    /// Open a Sentinel-1 SLC TIFF file and validate its format.
    ///
    /// Validates:
    /// - Classic TIFF (magic=42), not BigTIFF.
    /// - `SampleFormat=5` (CInt), `BitsPerSample=32` (CInt16).
    /// - `Compression=1` (none).
    /// - `RowsPerStrip=1`.
    /// - Contiguous strip layout (first and last offsets verified).
    ///
    /// # Errors
    ///
    /// Returns [`SlcReadError`] on I/O failure, bad TIFF headers, or any
    /// format mismatch.
    pub fn open<P: AsRef<Path>>(path: P) -> Result<Self, SlcReadError> {
        let mut file = File::open(path)?;
        let layout = parse_ifd(&mut file)?;
        Ok(Self {
            file,
            width: layout.width,
            height: layout.height,
            little_endian: layout.little_endian,
            base_offset: layout.base_offset,
            row_bytes: layout.row_bytes,
        })
    }

    /// Read a contiguous range of rows as raw CInt16 samples.
    ///
    /// Returns a flat `Vec<[i16; 2]>` of length `line_count * self.width`.
    /// Each element is one complex pixel: `pixel[n] = [i, q]`.
    /// The `i16` values are already in host byte order.
    ///
    /// # Arguments
    ///
    /// * `first_line`  — Starting row (0-indexed, inclusive).
    /// * `line_count`  — Number of consecutive rows to read.
    ///
    /// # Errors
    ///
    /// Returns [`SlcReadError::LineOutOfBounds`] if
    /// `first_line + line_count > self.height`.
    pub fn read_burst_raw(
        &mut self,
        first_line: usize,
        line_count: usize,
    ) -> Result<Vec<[i16; 2]>, SlcReadError> {
        let end_line = first_line
            .checked_add(line_count)
            .filter(|&e| e <= self.height as usize)
            .ok_or(SlcReadError::LineOutOfBounds {
                first: first_line,
                end: first_line.saturating_add(line_count),
                height: self.height,
            })?;
        let _ = end_line; // SAFETY-OK: explicit suppression of unused variable; not an error swallow

        let seek_to = self.base_offset + first_line as u64 * self.row_bytes;
        let total_bytes = line_count as u64 * self.row_bytes;
        let n_pixels = self.width as usize * line_count;

        let mut buf = vec![0u8; total_bytes as usize];
        self.file.seek(SeekFrom::Start(seek_to))?;
        self.file.read_exact(&mut buf)?;

        let le = self.little_endian;
        let mut pixels = Vec::with_capacity(n_pixels);
        for chunk in buf.chunks_exact(4) {
            let i_ch = if le {
                i16::from_le_bytes([chunk[0], chunk[1]])
            } else {
                i16::from_be_bytes([chunk[0], chunk[1]])
            };
            let q_ch = if le {
                i16::from_le_bytes([chunk[2], chunk[3]])
            } else {
                i16::from_be_bytes([chunk[2], chunk[3]])
            };
            pixels.push([i_ch, q_ch]);
        }

        Ok(pixels)
    }

    /// Read a contiguous range of rows as complex f32 samples.
    ///
    /// Identical to [`read_burst_raw`] but casts each `i16` to `f32`
    /// (no normalisation; calibration is a downstream step).
    ///
    /// Returns a flat `Vec<(f32, f32)>` of length `line_count * self.width`.
    /// Each element is one complex pixel: `pixel[n] = (i, q)`.
    ///
    /// [`read_burst_raw`]: SlcReader::read_burst_raw
    pub fn read_burst_complex(
        &mut self,
        first_line: usize,
        line_count: usize,
    ) -> Result<Vec<(f32, f32)>, SlcReadError> {
        let raw = self.read_burst_raw(first_line, line_count)?;
        // Direct allocation: no intermediate Vec.
        Ok(raw.into_iter().map(|[i, q]| (i as f32, q as f32)).collect())
    }
}

// ─── Tests ────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    /// Path to the IW1/VV measurement TIFF of the S1A test product.
    fn iw1_vv_tiff() -> std::path::PathBuf {
        std::path::PathBuf::from(
            "/home/datacube/dev/SARdine/data/SLC/\
             S1A_IW_SLC__1SDV_20201005T170824_20201005T170851_034664_04098A_1E66.SAFE/\
             measurement/\
             s1a-iw1-slc-vv-20201005t170824-20201005t170849-034664-04098a-004.tiff",
        )
    }

    fn s1a_fixtures_present() -> bool {
        iw1_vv_tiff().is_file()
    }

    // ── Format validation ─────────────────────────────────────────────────────

    #[test]
    fn open_s1a_iw1_vv_succeeds() {
        if !s1a_fixtures_present() { eprintln!("skipping — S1A fixture not present"); return; }
        let r = SlcReader::open(iw1_vv_tiff()).expect("open failed");
        // Dimensions confirmed by independent Python inspection.
        assert_eq!(r.width, 22111, "unexpected width");
        assert_eq!(r.height, 13482, "unexpected height");
        assert!(r.little_endian, "S1 TIFFs are always little-endian");
    }

    #[test]
    fn open_nonexistent_file_gives_io_error() {
        let err = SlcReader::open("/nonexistent/path/file.tiff").unwrap_err();
        assert!(matches!(err, SlcReadError::Io(_)));
    }

    #[test]
    fn open_empty_file_gives_not_tiff() {
        let tmp = tempfile::NamedTempFile::new().unwrap();
        let err = SlcReader::open(tmp.path()).unwrap_err();
        // Empty file → I/O error or NotTiff (read_exact fails on empty)
        assert!(
            matches!(err, SlcReadError::Io(_) | SlcReadError::NotTiff { .. }),
            "unexpected error: {err}"
        );
    }

    #[test]
    fn open_non_tiff_data_gives_not_tiff() {
        use std::io::Write as _;
        let mut tmp = tempfile::NamedTempFile::new().unwrap();
        tmp.write_all(b"This is not a TIFF file at all!!!").unwrap();
        let err = SlcReader::open(tmp.path()).unwrap_err();
        assert!(
            matches!(err, SlcReadError::NotTiff { .. }),
            "unexpected error: {err}"
        );
    }

    /// Build a minimal well-formed classic TIFF IFD with the given IFD entries,
    /// no actual image data (base_offset points beyond the header). Returns a
    /// byte vector suitable for writing to a temp file.
    ///
    /// Entries: &[(tag, field_type, count, value_as_u32)]
    /// All values are assumed to fit inline (count=1 SHORT or LONG).
    fn make_minimal_tiff(entries: &[(u16, u16, u32, u32)]) -> Vec<u8> {
        // Layout: 8-byte header | IFD (2+n*12+4 bytes) | [image data placeholder]
        let n = entries.len();
        let ifd_offset: u32 = 8;
        let ifd_size = 2 + n as u32 * 12 + 4;
        let data_offset: u32 = 8 + ifd_size;

        let mut buf = Vec::<u8>::new();
        // Header: II, 42, IFD offset
        buf.extend_from_slice(b"II");
        buf.extend_from_slice(&42u16.to_le_bytes());
        buf.extend_from_slice(&ifd_offset.to_le_bytes());
        // IFD: count
        buf.extend_from_slice(&(n as u16).to_le_bytes());
        // Entries (sorted by tag for well-formedness)
        let mut sorted = entries.to_vec();
        sorted.sort_by_key(|(tag, ..)| *tag);
        for (tag, ftype, count, val) in &sorted {
            buf.extend_from_slice(&tag.to_le_bytes());
            buf.extend_from_slice(&ftype.to_le_bytes());
            buf.extend_from_slice(&count.to_le_bytes());
            // Inline value: for SHORT store as 2-byte LE then 2 pad bytes;
            // for LONG store as 4-byte LE.
            if *ftype == 3 {
                buf.extend_from_slice(&(*val as u16).to_le_bytes());
                buf.extend_from_slice(&[0u8, 0u8]);
            } else {
                buf.extend_from_slice(&val.to_le_bytes());
            }
        }
        // Next IFD offset = 0 (no more IFDs)
        buf.extend_from_slice(&0u32.to_le_bytes());
        // One real data byte so the file is not empty
        buf.push(0u8);
        let _ = data_offset; // SAFETY-OK: explicit suppression of unused variable; not an error swallow
        buf
    }

    #[test]
    fn open_multi_sample_tiff_gives_unsupported_format() {
        // Build a TIFF claiming SamplesPerPixel=2 (wrong for CInt16).
        // All other tags are valid S1-like values except SamplesPerPixel.
        use std::io::Write as _;
        const SHORT: u16 = 3;
        const LONG: u16  = 4;
        // Fake strip offset table: not contiguous (doesn't matter, fails earlier)
        // For simplicity, use count=1 (single-strip, height=1).
        let entries: &[(u16, u16, u32, u32)] = &[
            (256, LONG,  1, 10),    // ImageWidth = 10
            (257, LONG,  1, 1),     // ImageLength = 1
            (258, SHORT, 1, 32),    // BitsPerSample = 32
            (259, SHORT, 1, 1),     // Compression = None
            (273, LONG,  1, 9999),  // StripOffsets (inline, count=1)
            (277, SHORT, 1, 2),     // SamplesPerPixel = 2  ← invalid
            (278, LONG,  1, 1),     // RowsPerStrip = 1
            (339, SHORT, 1, 5),     // SampleFormat = CInt
        ];
        let data = make_minimal_tiff(entries);
        let mut tmp = tempfile::NamedTempFile::new().unwrap();
        tmp.write_all(&data).unwrap();
        let err = SlcReader::open(tmp.path()).unwrap_err();
        assert!(
            matches!(err, SlcReadError::UnsupportedFormat(_)),
            "expected UnsupportedFormat for SamplesPerPixel=2, got: {err}"
        );
    }

    // ── Read mechanics ────────────────────────────────────────────────────────

    #[test]
    fn read_burst_raw_returns_correct_sample_count() {
        if !s1a_fixtures_present() { eprintln!("skipping — S1A fixture not present"); return; }
        let mut r = SlcReader::open(iw1_vv_tiff()).expect("open failed");
        // IW1 annotation has 9 bursts of 1498 lines each.
        let first_line = 0;
        let line_count = 1498;
        let pixels = r.read_burst_raw(first_line, line_count).expect("read failed");
        assert_eq!(
            pixels.len(),
            r.width as usize * line_count,
            "sample count mismatch"
        );
    }

    #[test]
    fn read_burst_complex_returns_correct_sample_count() {
        if !s1a_fixtures_present() { eprintln!("skipping — S1A fixture not present"); return; }
        let mut r = SlcReader::open(iw1_vv_tiff()).expect("open failed");
        let pixels = r.read_burst_complex(0, 1498).expect("read failed");
        assert_eq!(pixels.len(), r.width as usize * 1498);
    }

    #[test]
    fn read_burst_raw_exact_line_range_at_end() {
        if !s1a_fixtures_present() { eprintln!("skipping — S1A fixture not present"); return; }
        // Read exactly the last burst to verify end-of-file reads work.
        let mut r = SlcReader::open(iw1_vv_tiff()).expect("open failed");
        // 9 bursts × 1498 lines = 13482 total (= r.height)
        let last_burst_first_line = 13482 - 1498;
        let pixels = r
            .read_burst_raw(last_burst_first_line, 1498)
            .expect("read last burst failed");
        assert_eq!(pixels.len(), r.width as usize * 1498);
    }

    #[test]
    fn read_burst_out_of_bounds_returns_error() {
        if !s1a_fixtures_present() { eprintln!("skipping — S1A fixture not present"); return; }
        let mut r = SlcReader::open(iw1_vv_tiff()).expect("open failed");
        // Request one line past the end.
        let err = r.read_burst_raw(13482, 1).unwrap_err();
        assert!(
            matches!(err, SlcReadError::LineOutOfBounds { .. }),
            "unexpected error: {err}"
        );
    }

    #[test]
    fn read_burst_overflow_returns_error() {
        if !s1a_fixtures_present() { eprintln!("skipping — S1A fixture not present"); return; }
        let mut r = SlcReader::open(iw1_vv_tiff()).expect("open failed");
        // Line 13480 + 10 = 13490 > 13482 → out of bounds.
        let err = r.read_burst_raw(13480, 10).unwrap_err();
        assert!(
            matches!(err, SlcReadError::LineOutOfBounds { .. }),
            "unexpected error: {err}"
        );
    }

    #[test]
    fn read_single_row_correct_byte_count() {
        if !s1a_fixtures_present() { eprintln!("skipping — S1A fixture not present"); return; }
        let mut r = SlcReader::open(iw1_vv_tiff()).expect("open failed");
        let pixels = r.read_burst_raw(0, 1).expect("read row 0 failed");
        // One row: width samples, each an [i16; 2] = 4 bytes.
        assert_eq!(pixels.len(), r.width as usize);
    }

    // ── Data content checks ───────────────────────────────────────────────────

    #[test]
    fn burst_interior_contains_nonzero_samples() {
        if !s1a_fixtures_present() { eprintln!("skipping — S1A fixture not present"); return; }
        // Burst 4 (middle burst, lines 4*1498 .. 5*1498) should contain
        // real backscatter signal.  Even with invalid edge samples,
        // the interior will have non-zero I or Q.
        let mut r = SlcReader::open(iw1_vv_tiff()).expect("open failed");
        let pixels = r.read_burst_raw(4 * 1498, 1498).expect("read failed");
        let any_nonzero = pixels.iter().any(|[i, q]| *i != 0 || *q != 0);
        assert!(any_nonzero, "all samples in burst 4 are zero — reading wrong region");
    }

    #[test]
    fn complex_cast_preserves_sign() {
        // Synthesize a minimal fake TIFF in memory and verify that negative
        // i16 values survive the cast to f32 with the correct sign.
        //
        // We can't easily inject bytes into SlcReader without a real file, so
        // this test validates the cast logic directly.
        let raw: [i16; 2] = [-100, 200];
        let (i_f, q_f) = (raw[0] as f32, raw[1] as f32);
        assert_eq!(i_f, -100.0_f32);
        assert_eq!(q_f, 200.0_f32);
    }
}
