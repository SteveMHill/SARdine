"""
Compare sardine VV dB output against the Microsoft Planetary Computer (MPC)
Sentinel-1 RTC product for the same S1B 2019-01-23 scene.

MPC product: S1B_IW_GRDH_1SDV_20190123T053349_20190123T053414_014617_01B3D4_rtc
  - Processor:   Catalyst Earth
  - DEM:         PlanetDEM
  - Input:       GRD-H (5 range × 1 az looks, ~20 × 22 m resolution)
  - CRS:         EPSG:32632 (UTM 32N)
  - Pixel size:  10 m
  - Radiometry:  γ⁰ linear power, nodata = −32768.0
  - Format:      Float32 COG (~1.75 GB; requires MPC SAS token)

Sardine output:
  - CRS:         EPSG:4326, 0.0001° (~9 m) spacing
  - Radiometry:  γ⁰ dB (10 × log₁₀), NaN = nodata
  - Input:       SLC (single-look complex, terrain-corrected with POEORB + SRTM-1)

Method:
  1. Sign the MPC VV COG URL at runtime via planetary_computer.sign().
  2. Read the MPC raster windowed to the sardine scene's bounding box to avoid
     loading the full 1.75 GB into memory.
  3. Apply nodata mask (MPC nodata = −32768).
  4. Reproject MPC linear γ⁰ → sardine WGS84 grid (bilinear).
  5. Convert MPC linear → dB.
  6. Compute bias, RMSE, histogram on the joint valid pixel set.

Key systematic differences (not bugs):
  - GRD vs SLC input: different multi-looking and speckle pattern.
    Per-pixel RMSE is expected to be larger than in the sardine−ASF comparison
    because the two rasters are *different realisations* of the same scene (not
    co-registered speckle-coherent data).  **Median bias** is the relevant metric.
  - Orbit quality: MPC uses DOWNLINK orbits for this scene; sardine uses POEORB.
    In steep terrain this may cause sub-pixel geolocation offsets, inflating RMSE.
  - DEM: MPC uses PlanetDEM; sardine uses SRTM-1.  Differences affect the
    gamma-flattening correction, particularly on steep slopes.

Usage::

    # With local sardine output already produced:
    SARDINE=/path/to/sardine_s1b_vv.tif python3 scripts/compare_mpc.py

    # Or let the script run sardine first (set SAFE / DEM / ORBIT env vars):
    python3 scripts/compare_mpc.py

Environment variables:
  SARDINE     Path to the sardine γ⁰ dB GeoTIFF to evaluate.
              Default: /home/datacube/dev/SARdine/sardine_s1b_tc_db.tiff
  MPC_VV_TIF  Path to a previously downloaded MPC VV COG (skip download).
              If not set, the script signs and streams the URL.
  SKIP_DOWNLOAD
              If set to "1", expect MPC_VV_TIF to point to a local file.
"""

from __future__ import annotations

import os
import sys
import tempfile
import urllib.request
from pathlib import Path

import numpy as np
import rasterio
import rasterio.transform
from rasterio.crs import CRS
from rasterio.warp import reproject, Resampling
import rasterio.windows

try:
    import planetary_computer
except ImportError:
    print(
        "ERROR: planetary_computer package not installed.\n"
        "       pip install planetary-computer",
        file=sys.stderr,
    )
    sys.exit(1)

# ─── Configuration ────────────────────────────────────────────────────────────

SARDINE = os.environ.get(
    "SARDINE",
    "/home/datacube/dev/SARdine/sardine_s1b_tc_db.tiff",
)

# The MPC COG for this exact acquisition (orbit 14617, datatake 01B3D4).
# Product: S1B_IW_GRDH_1SDV_20190123T053349_20190123T053414 — this matches
# the SLC scene (053348-053415) with ~25 s / ~185 km of spatial overlap.
# Signing this URL at runtime adds a short-lived SAS token — no token is
# stored in this file.
_MPC_VV_UNSIGNED = (
    "https://sentinel1euwestrtc.blob.core.windows.net/sentinel1-grd-rtc/"
    "GRD/2019/1/23/IW/DV/"
    "S1B_IW_GRDH_1SDV_20190123T053349_20190123T053414_014617_01B3D4_39E0/"
    "measurement/iw-vv.rtc.tiff"
)

# Nodata value used by the MPC/Catalyst product.
MPC_NODATA = -32768.0

# Default local cache path for the downloaded COG.
_DEFAULT_CACHE = (
    "/home/datacube/dev/SARdine/data/MPC/S1B_IW_20190123T053349_VV_rtc.tiff"
)

# ─── Helpers ──────────────────────────────────────────────────────────────────


def sign_mpc_url(url: str) -> str:
    """Return a SAS-signed version of an MPC Azure Blob URL."""
    signed = planetary_computer.sign(url)
    # sign() may return the URL unchanged if no subscription key is configured.
    if signed == url:
        print(
            "WARNING: planetary_computer.sign() returned the URL unchanged.\n"
            "         Set the MPC subscription key:\n"
            "           planetary-computer token --subscription-key <KEY>",
            file=sys.stderr,
        )
    return signed


def get_mpc_path() -> tuple[str | None, bool]:
    """Return (path_or_url, is_local).

    Precedence:
      1. MPC_VV_TIF env var (assumed local file or already-signed URL).
      2. Default cache path if the file already exists.
      3. Sign and download from MPC.
    """
    env = os.environ.get("MPC_VV_TIF")
    if env:
        return env, Path(env).exists()

    if Path(_DEFAULT_CACHE).is_file():
        print(f"Using cached MPC VV tiff: {_DEFAULT_CACHE}")
        return _DEFAULT_CACHE, True

    return None, False


def download_mpc_vv(dest: str | Path) -> None:
    """Sign the MPC URL and stream-download to dest."""
    dest = Path(dest)
    dest.parent.mkdir(parents=True, exist_ok=True)

    print("Signing MPC URL …", flush=True)
    signed_url = sign_mpc_url(_MPC_VV_UNSIGNED)

    print(f"Downloading MPC VV COG (~1.75 GB) → {dest}", flush=True)
    print("  (This is a one-time download; future runs use the cached file.)")

    chunk_size = 8 * 1024 * 1024  # 8 MiB
    tmp = str(dest) + ".tmp"
    try:
        req = urllib.request.Request(signed_url)
        with urllib.request.urlopen(req) as response, open(tmp, "wb") as f:
            total = int(response.headers.get("Content-Length", 0))
            downloaded = 0
            while True:
                chunk = response.read(chunk_size)
                if not chunk:
                    break
                f.write(chunk)
                downloaded += len(chunk)
                if total:
                    pct = 100 * downloaded / total
                    mb = downloaded / 1024 / 1024
                    print(
                        f"\r  {mb:,.0f} / {total/1024/1024:,.0f} MB  "
                        f"({pct:.1f}%)",
                        end="",
                        flush=True,
                    )
        print()  # newline after progress
        Path(tmp).rename(dest)
        print(f"  Saved → {dest}")
    except Exception:
        Path(tmp).unlink(missing_ok=True)
        raise


def read_mpc_windowed(
    mpc_path: str,
    sardine_bounds,
    sardine_crs,
) -> tuple[np.ndarray, rasterio.transform.Affine, CRS]:
    """Read the MPC VV raster windowed to the sardine bounding box.

    Projects sardine_bounds into MPC CRS to construct the window, so only
    the relevant ~10% of the 1.75 GB file is read.

    Returns (data_f32, transform, crs) where data has nodata pixels set to NaN.
    """
    from rasterio.warp import transform_bounds

    with rasterio.open(mpc_path) as ds:
        mpc_crs = ds.crs

        # Project sardine bounds → MPC CRS to build the read window.
        left, bottom, right, top = transform_bounds(
            sardine_crs, mpc_crs,
            sardine_bounds.left, sardine_bounds.bottom,
            sardine_bounds.right, sardine_bounds.top,
        )
        # Add 1 % margin on each side.
        dx = (right - left) * 0.01
        dy = (top - bottom) * 0.01
        win = ds.window(
            left - dx, bottom - dy,
            right + dx, top + dy,
        )
        win = win.intersection(rasterio.windows.Window(0, 0, ds.width, ds.height))

        print(
            f"  MPC window: {int(win.height)} rows × {int(win.width)} cols "
            f"(of {ds.height} × {ds.width})",
            flush=True,
        )

        data = ds.read(1, window=win).astype(np.float32)
        win_transform = ds.window_transform(win)

    # Apply nodata mask.
    data[data == MPC_NODATA] = np.nan
    # Also mask any remaining non-positive values (shouldn't exist, but defensive).
    data[data <= 0.0] = np.nan

    return data, win_transform, mpc_crs


# ─── Main ─────────────────────────────────────────────────────────────────────


def main() -> None:
    # ── 1. Read sardine output ────────────────────────────────────────────────
    print(f"Reading sardine output: {SARDINE}", flush=True)
    if not Path(SARDINE).is_file():
        print(f"ERROR: sardine output not found: {SARDINE}", file=sys.stderr)
        print(
            "Set the SARDINE env var to the path of a sardine γ⁰ dB GeoTIFF,\n"
            "or run the pipeline first:\n"
            "  cargo run --release --bin sardine -- process \\\n"
            f"    --safe <SAFE> --dem <DEM> --geoid auto \\\n"
            f"    --output {SARDINE}",
            file=sys.stderr,
        )
        sys.exit(1)

    with rasterio.open(SARDINE) as ds:
        sardine_data = ds.read(1).astype(np.float32)
        sardine_transform = ds.transform
        sardine_crs = ds.crs
        sardine_shape = (ds.height, ds.width)
        sardine_bounds = ds.bounds

    sardine_valid = np.isfinite(sardine_data)
    print(
        f"  {sardine_data.shape[1]} × {sardine_data.shape[0]} pixels, "
        f"valid={sardine_valid.sum():,} ({100 * sardine_valid.mean():.1f}%)"
    )
    print(
        f"  dB range:  {np.nanmin(sardine_data):.2f} .. "
        f"{np.nanmax(sardine_data):.2f}"
    )
    print(f"  median dB: {np.nanmedian(sardine_data):.3f}")

    # ── 2. Obtain MPC VV COG ─────────────────────────────────────────────────
    mpc_path, is_local = get_mpc_path()

    if not is_local or mpc_path is None:
        mpc_path = _DEFAULT_CACHE
        download_mpc_vv(mpc_path)

    print(f"\nReading MPC VV COG (windowed): {mpc_path}", flush=True)
    mpc_raw, mpc_win_transform, mpc_crs = read_mpc_windowed(
        mpc_path, sardine_bounds, sardine_crs
    )

    mpc_valid_frac = np.isfinite(mpc_raw).mean()
    print(
        f"  MPC window shape: {mpc_raw.shape[0]} × {mpc_raw.shape[1]}, "
        f"valid={100 * mpc_valid_frac:.1f}%"
    )
    print(
        f"  linear range: "
        f"{np.nanmin(mpc_raw):.2e} .. {np.nanmax(mpc_raw):.2e}"
    )

    # ── 3. Reproject MPC linear γ⁰ → sardine WGS84 grid ─────────────────────
    print("\nReprojecting MPC → sardine grid (bilinear) …", flush=True)
    mpc_reproj = np.full(sardine_shape, np.nan, dtype=np.float32)
    reproject(
        source=mpc_raw,
        destination=mpc_reproj,
        src_transform=mpc_win_transform,
        src_crs=mpc_crs,
        dst_transform=sardine_transform,
        dst_crs=sardine_crs,
        resampling=Resampling.bilinear,
        src_nodata=np.nan,
        dst_nodata=np.nan,
    )

    # ── 4. Convert MPC linear → dB ───────────────────────────────────────────
    with np.errstate(divide="ignore", invalid="ignore"):
        mpc_db = np.where(mpc_reproj > 0, 10.0 * np.log10(mpc_reproj), np.nan)

    mpc_valid = np.isfinite(mpc_db)
    print(
        f"  MPC (reprojected): valid={mpc_valid.sum():,} "
        f"({100 * mpc_valid.mean():.1f}%)"
    )
    print(
        f"  dB range:  {np.nanmin(mpc_db):.2f} .. {np.nanmax(mpc_db):.2f}"
    )
    print(f"  median dB: {np.nanmedian(mpc_db):.3f}")

    # ── 5. Joint valid mask ───────────────────────────────────────────────────
    joint = sardine_valid & mpc_valid
    n_joint = int(joint.sum())
    print(f"\nJoint valid pixels: {n_joint:,}")
    if n_joint == 0:
        print("ERROR: no overlap — check extents/CRS", file=sys.stderr)
        sys.exit(1)

    s = sardine_data[joint]
    m = mpc_db[joint]
    diff = s - m

    # ── 6. Statistics ─────────────────────────────────────────────────────────
    print()
    print("══════════════════════════════════════════════════════════")
    print("  Radiometric comparison: sardine − MPC RTC (Catalyst)   ")
    print("══════════════════════════════════════════════════════════")
    print()
    print("  NOTE: The MPC product is derived from GRD-H (multi-looked),")
    print("  while sardine processes SLC data.  Both are independent")
    print("  realisations of the same backscatter field with uncorrelated")
    print("  speckle, so per-pixel RMSE will be higher than in the")
    print("  sardine−ASF comparison.  Median bias is the key metric.")
    print()
    print(f"  Mean bias:    {diff.mean():+.4f} dB")
    print(f"  Median bias:  {np.median(diff):+.4f} dB")
    print(f"  Std dev:      {diff.std():.4f} dB")
    print(f"  RMSE:         {np.sqrt((diff ** 2).mean()):.4f} dB")

    print()
    print("── |error| percentiles ─────────────────────────────────")
    for p in [50, 68, 90, 95, 99]:
        print(f"  p{p:2d}: |error| ≤ {np.percentile(np.abs(diff), p):.3f} dB")

    print()
    print("── Bias histogram (sardine − MPC) ──────────────────────")
    bins = np.arange(-5, 5.1, 0.5)
    counts, edges = np.histogram(diff, bins=bins)
    total = counts.sum()
    peak = counts.max()
    for i, c in enumerate(counts):
        bar = "█" * int(40 * c / peak)
        lo, hi = edges[i], edges[i + 1]
        print(f"  [{lo:+.1f},{hi:+.1f})  {100 * c / total:5.1f}%  {bar}")

    print()
    print("── Spatial bias by latitude band (2° bins) ─────────────")
    rows_idx, _ = np.where(joint)
    lats_all = sardine_transform.f + rows_idx * sardine_transform.e
    bin_edges = np.arange(46, 50, 2)
    for i in range(len(bin_edges) - 1):
        mask = (lats_all >= bin_edges[i]) & (lats_all < bin_edges[i + 1])
        if mask.sum() > 0:
            b = diff[mask]
            print(
                f"  lat [{bin_edges[i]:.0f},{bin_edges[i+1]:.0f}): "
                f"n={mask.sum():,}  bias={b.mean():+.3f}  std={b.std():.3f}"
            )

    print()
    print("── Reference comparison (for context) ──────────────────")
    print("  sardine − ASF RTC10 GAMMA:  median bias ≈ +0.017 dB, RMSE ≈ 0.60 dB")
    print(f"  sardine − MPC Catalyst:     median bias = {np.median(diff):+.3f} dB, "
          f"RMSE = {np.sqrt((diff**2).mean()):.2f} dB")


if __name__ == "__main__":
    main()
