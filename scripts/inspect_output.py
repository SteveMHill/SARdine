#!/usr/bin/env python3
"""
Inspect a sardine output GeoTIFF.

Checks:
  1. File size, dimensions, geotransform
  2. Pixel value distribution (dB histogram, percentiles)
  3. NaN fraction
  4. Known-point spot check — Basel train station (approx 47.547°N, 7.590°E)
     and Freiburg city centre (47.994°N, 7.849°E) — both in S1B footprint
  5. Row-profile to detect burst seams or subswath boundaries

Usage:
    python3 scripts/inspect_output.py /tmp/sardine_s1b_vv_db.tiff
"""

import struct
import math
import sys

TIFF_PATH = sys.argv[1] if len(sys.argv) > 1 else "/tmp/sardine_s1b_vv_db.tiff"

# ── Parse minimal TIFF header ───────────────────────────────────────────────

def read_u16(f, off): f.seek(off); return struct.unpack("<H", f.read(2))[0]
def read_u32(f, off): f.seek(off); return struct.unpack("<I", f.read(4))[0]
def read_f64(f, off): f.seek(off); return struct.unpack("<d", f.read(8))[0]
def read_u16_le(data, i): return struct.unpack_from("<H", data, i)[0]
def read_u32_le(data, i): return struct.unpack_from("<I", data, i)[0]
def read_f64_le(data, i): return struct.unpack_from("<d", data, i)[0]

with open(TIFF_PATH, "rb") as f:
    magic = f.read(4)
    assert magic == b"II\x2a\x00", f"unexpected TIFF magic: {magic!r}"

    ifd_off = read_u32(f, 4)
    f.seek(ifd_off)
    n_entries = struct.unpack("<H", f.read(2))[0]

    tags = {}
    for _ in range(n_entries):
        raw = f.read(12)
        tag  = struct.unpack_from("<H", raw, 0)[0]
        typ  = struct.unpack_from("<H", raw, 2)[0]
        cnt  = struct.unpack_from("<I", raw, 4)[0]
        voff = struct.unpack_from("<I", raw, 8)[0]
        tags[tag] = (typ, cnt, voff, raw[8:12])

    def tag_val_u32(t):
        typ, cnt, voff, raw = tags[t]
        if typ == 3:  # SHORT
            return struct.unpack_from("<H", raw, 0)[0]
        return voff  # LONG stored directly if cnt==1

    width  = tag_val_u32(256)
    height = tag_val_u32(257)
    rps    = tag_val_u32(278)  # rows per strip

    # ModelPixelScaleTag (33550) and ModelTiepointTag (33922)
    _, cnt_ps, off_ps, _ = tags[33550]
    _, cnt_tp, off_tp, _ = tags[33922]

    f.seek(off_ps)
    px_w, px_h = struct.unpack("<dd", f.read(16))[:2]

    f.seek(off_tp)
    tp = struct.unpack("<dddddd", f.read(48))
    origin_lon = tp[3]
    origin_lat = tp[4]

    # Strip offsets (tag 273 = StripOffsets)
    _, cnt_so, off_so, raw_so = tags[273]
    if cnt_so == 1:
        strip_offsets = [struct.unpack_from("<I", raw_so, 0)[0]]
    else:
        f.seek(off_so)
        strip_offsets = [struct.unpack("<I", f.read(4))[0] for _ in range(cnt_so)]

    # image data start = first strip offset
    image_data_off = strip_offsets[0]

# ── Summary ─────────────────────────────────────────────────────────────────

import os
fsize_mb = os.path.getsize(TIFF_PATH) / 1e6

print("=" * 70)
print(f"File    : {TIFF_PATH}")
print(f"Size    : {fsize_mb:.1f} MB")
print(f"Grid    : {width} cols × {height} rows  (rows_per_strip={rps})")
print(f"Pixel   : {px_w*3600:.4f}\"×{px_h*3600:.4f}\"  ({px_w:.6f}°×{px_h:.6f}°)")
print(f"Origin  : lon={origin_lon:.5f}°  lat={origin_lat:.5f}°  (top-left corner)")
lat_south = origin_lat - height * px_h
lon_east  = origin_lon + width  * px_w
print(f"Extent  : lon=[{origin_lon:.4f}, {lon_east:.4f}]  lat=[{lat_south:.4f}, {origin_lat:.4f}]")
print()

# ── Pixel value distribution (sample every 100th row) ────────────────────

def read_row(f, row):
    """Read a full row, handling multi-strip layout."""
    strip_idx = row // rps
    row_in_strip = row % rps
    off = strip_offsets[strip_idx] + row_in_strip * width * 4
    f.seek(off)
    return f.read(width * 4)

print("Sampling pixel values (every 100th row) …")
SAMPLE_STEP = 100
values = []
with open(TIFF_PATH, "rb") as f:
    for row in range(0, height, SAMPLE_STEP):
        raw = read_row(f, row)
        for c in range(0, width, 10):  # every 10th col within sampled row
            v = struct.unpack_from("<f", raw, c * 4)[0]
            values.append(v)

finite = [v for v in values if not math.isnan(v) and math.isfinite(v)]
nan_count = sum(1 for v in values if math.isnan(v))
nan_frac = nan_count / len(values) * 100

print(f"  Sample size : {len(values)} pixels")
print(f"  NaN fraction: {nan_frac:.1f}%  (valid: {len(finite)})")

if finite:
    finite.sort()
    n = len(finite)
    pcts = [1, 5, 10, 25, 50, 75, 90, 95, 99]
    print(f"  dB percentiles:")
    for p in pcts:
        idx = min(int(p/100 * n), n-1)
        print(f"    p{p:02d}: {finite[idx]:+.2f} dB")
    print(f"  min: {finite[0]:+.2f} dB   max: {finite[-1]:+.2f} dB")

    in_range = sum(1 for v in finite if -35.0 < v < 5.0)
    print(f"  Pixels in plausible S1 VV range (-35 to +5 dB): {in_range}/{len(finite)}  ({100*in_range/len(finite):.1f}%)")
print()

# ── Known-point spot checks ──────────────────────────────────────────────────

def sample_point(f, name, lat, lon):
    col = (lon - origin_lon) / px_w
    row = (origin_lat - lat) / px_h
    c = int(round(col)); r = int(round(row))
    if not (0 <= r < height and 0 <= c < width):
        print(f"  {name:30s}: OUTSIDE grid (r={r}, c={c})")
        return
    strip_idx = r // rps
    row_in_strip = r % rps
    offset = strip_offsets[strip_idx] + row_in_strip * width * 4 + c * 4
    f.seek(offset)
    v = struct.unpack("<f", f.read(4))[0]
    label = f"{v:+.2f} dB" if not math.isnan(v) else "NaN"
    print(f"  {name:30s}: r={r:5d} c={c:5d}  σ⁰={label}")

print("Known-point spot checks (σ⁰ VV):")
with open(TIFF_PATH, "rb") as f:
    sample_point(f, "Basel train station",   47.5472,  7.5897)
    sample_point(f, "Freiburg Münster",       47.9944,  7.8526)
    sample_point(f, "Strasbourg centre",      48.5734,  7.7521)
    sample_point(f, "Bern Bundeshaus",        46.9464,  7.4441)
    sample_point(f, "Karlsruhe Schloss",      49.0135,  8.4044)
    sample_point(f, "Stuttgart Hauptbf.",     48.7836,  9.1830)
    sample_point(f, "Zurich HB",              47.3784,  8.5400)
    sample_point(f, "Lake Constance centre",  47.6500,  9.2300)
print()

# ── Row profile: detect burst seams / subswath boundaries ───────────────────

print("Row-mean profile (every 500 rows, centre 20% of cols):")
c_lo = width  * 4 // 10
c_hi = width  * 6 // 10
with open(TIFF_PATH, "rb") as f:
    prev_mean = None
    for row in range(0, height, 500):
        strip_idx = row // rps
        row_in_strip = row % rps
        f.seek(strip_offsets[strip_idx] + row_in_strip * width * 4 + c_lo * 4)
        raw = f.read((c_hi - c_lo) * 4)
        vals = [struct.unpack_from("<f", raw, i*4)[0] for i in range(c_hi - c_lo)]
        finite_v = [v for v in vals if not math.isnan(v)]
        mean = sum(finite_v)/len(finite_v) if finite_v else float("nan")
        lat_r = origin_lat - row * px_h
        jump = ""
        if prev_mean is not None and not math.isnan(mean) and not math.isnan(prev_mean):
            delta = abs(mean - prev_mean)
            if delta > 2.0:
                jump = f"  *** JUMP {delta:+.2f} dB — possible seam ***"
        print(f"  row {row:5d}  lat={lat_r:.4f}°  mean={mean:+.2f} dB{jump}")
        prev_mean = mean

print()
print("Inspection complete.")
