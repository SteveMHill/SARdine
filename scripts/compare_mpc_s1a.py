"""
Compare sardine VV dB output against the Microsoft Planetary Computer (MPC)
Sentinel-1 RTC product for the S1A 2020-10-05 scene.

MPC product: S1A_IW_GRDH_1SDV_20201005T170824_20201005T170849_034664_04098A_rtc
  - Processor:   Catalyst Earth
  - DEM:         PlanetDEM
  - Input:       GRD-H (5 range × 1 az looks, ~20 × 22 m resolution)
  - CRS:         EPSG:32632 (UTM 32N)
  - Pixel size:  10 m
  - Radiometry:  γ⁰ linear power, nodata = −32768.0
  - Format:      Float32 COG (~1.75 GB; requires MPC SAS token)

Sardine output (expected):
  - CRS:         EPSG:32632, 10 m spacing
  - Radiometry:  γ⁰ dB (10 × log₁₀), NaN = nodata
  - Input:       SLC (single-look complex, terrain-corrected with POEORB + SRTM-1)
  - Speckle:     Refined Lee 7×7 (ENL=1) recommended to match GRD multilook

Method:
  1. Sign the MPC VV COG URL at runtime via planetary_computer.sign().
  2. Read the MPC raster windowed to the sardine scene bounding box.
  3. Apply nodata mask (MPC nodata = −32768).
  4. Reproject MPC linear γ⁰ → sardine grid (bilinear).
  5. Convert MPC linear → dB.
  6. Compute bias, RMSE, histogram on the joint valid pixel set.

Key systematic differences (not bugs):
  - GRD vs SLC input: MPC input is GRD (multi-looked ~20×22 m); sardine processes
    SLC (single-look).  Per-pixel RMSE will be higher than sardine−ASF because
    the two rasters are *different realisations* (uncorrelated speckle).
    Median bias is the key radiometric metric.
  - DEM: MPC uses PlanetDEM; sardine uses SRTM-1.  Differences matter most on
    steep slopes via terrain-flattening correction.
  - Orbit: sardine uses POEORB; MPC uses GRD-H restituted orbits.

Recommended sardine pipeline command::

    cargo run --release --features geoid-fetch -- process \\
        --safe  data/SLC/S1A_IW_SLC__1SDV_20201005T170824_20201005T170851_034664_04098A_1E66.SAFE \\
        --dem   data/dem/srtm1 \\
        --geoid auto \\
        --crs   EPSG:32632 \\
        --pixel-spacing-m 10 \\
        --speckle refined-lee --enl 1.0 \\
        --polarization VV \\
        --output sardine_s1a_20201005_vv.tiff

Usage::

    # With local sardine output already produced:
    SARDINE=/path/to/sardine_s1a_20201005_vv.tiff python3 scripts/compare_mpc_s1a.py

    # Or set default path:
    python3 scripts/compare_mpc_s1a.py

Environment variables:
  SARDINE     Path to the sardine γ⁰ dB GeoTIFF to evaluate.
              Default: sardine_s1a_20201005_vv.tiff (relative to cwd)
  MPC_VV_TIF  Path to a previously downloaded/cached MPC VV COG.
              If not set, the script signs and streams the URL.
"""

from __future__ import annotations

import os
import sys
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
    "sardine_s1a_20201005_vv.tiff",
)

# Unsigned MPC VV COG URL for S1A orbit 034664, datatake 04098A, granule 2F79.
# Signing this URL at runtime adds a short-lived SAS token — no token stored here.
_MPC_VV_UNSIGNED = (
    "https://sentinel1euwestrtc.blob.core.windows.net/sentinel1-grd-rtc/"
    "GRD/2020/10/5/IW/DV/"
    "S1A_IW_GRDH_1SDV_20201005T170824_20201005T170849_034664_04098A_2F79/"
    "measurement/iw-vv.rtc.tiff"
)

# Nodata value used by the MPC/Catalyst product.
MPC_NODATA = -32768.0

# Default local cache path.
_DEFAULT_CACHE = (
    "/home/datacube/dev/SARdine/data/MPC/S1A_IW_20201005_VV_rtc.tiff"
)

# ─── Helpers ──────────────────────────────────────────────────────────────────


def sign_mpc_url(url: str) -> str:
    """Return a SAS-signed version of an MPC Azure Blob URL."""
    signed = planetary_computer.sign(url)
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
        print()
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
    the relevant portion of the file is read rather than the full ~1.75 GB.

    Returns (data_f32, transform, crs) with nodata pixels set to NaN.
    """
    from rasterio.warp import transform_bounds

    with rasterio.open(mpc_path) as ds:
        mpc_crs = ds.crs

        left, bottom, right, top = transform_bounds(
            sardine_crs, mpc_crs,
            sardine_bounds.left, sardine_bounds.bottom,
            sardine_bounds.right, sardine_bounds.top,
        )
        dx = (right - left) * 0.01
        dy = (top - bottom) * 0.01
        win = ds.window(
            left - dx, bottom - dy,
            right + dx, top + dy,
        )
        win = win.intersection(
            rasterio.windows.Window(0, 0, ds.width, ds.height)
        )

        print(
            f"  MPC window: {int(win.height)} rows × {int(win.width)} cols "
            f"(of {ds.height} × {ds.width})",
            flush=True,
        )

        data = ds.read(1, window=win).astype(np.float32)
        win_transform = ds.window_transform(win)

    data[data == MPC_NODATA] = np.nan
    data[data <= 0.0] = np.nan

    return data, win_transform, mpc_crs


# ─── Main ─────────────────────────────────────────────────────────────────────


def main() -> None:
    # ── 1. Read sardine output ────────────────────────────────────────────────
    print(f"Reading sardine output: {SARDINE}", flush=True)
    if not Path(SARDINE).is_file():
        print(
            f"ERROR: sardine output not found: {SARDINE}\n"
            "\nRun the S1A pipeline first:\n"
            "  cargo run --release --features geoid-fetch -- process \\\n"
            "    --safe  data/SLC/S1A_IW_SLC__1SDV_20201005T170824_20201005T170851"
            "_034664_04098A_1E66.SAFE \\\n"
            "    --dem   data/dem/srtm1 \\\n"
            "    --geoid auto \\\n"
            "    --crs EPSG:32632 --pixel-spacing-m 10 \\\n"
            "    --speckle refined-lee --enl 1.0 \\\n"
            "    --polarization VV \\\n"
            "    --output sardine_s1a_20201005_vv.tiff\n"
            "Then: SARDINE=sardine_s1a_20201005_vv.tiff python3 scripts/compare_mpc_s1a.py",
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

    # ── 3. Reproject MPC linear γ⁰ → sardine grid ────────────────────────────
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
        print(
            "ERROR: no overlapping valid pixels.\n"
            "Check that the sardine output CRS matches the scene extent.\n"
            "Expected: sardine output in EPSG:32632, 10 m, covering the "
            "S1A 2020-10-05T170824 footprint.",
            file=sys.stderr,
        )
        sys.exit(1)

    s = sardine_data[joint]
    m = mpc_db[joint]
    diff = s - m

    # ── 6. Statistics ─────────────────────────────────────────────────────────
    print()
    print("══════════════════════════════════════════════════════════")
    print("  Radiometric comparison: sardine − MPC RTC (Catalyst)   ")
    print("  Scene: S1A 2020-10-05T170824, VV, orbit 034664         ")
    print("══════════════════════════════════════════════════════════")
    print()
    print("  NOTE: MPC input is GRD-H (multi-looked ~20×22 m); sardine")
    print("  processes SLC (single-look).  Both are independent realisations")
    print("  with uncorrelated speckle — per-pixel RMSE is expected to be")
    print("  higher than sardine−ASF.  Median bias is the key metric.")
    print()
    print(f"  n_joint:      {n_joint:,} pixels")
    print(f"  Mean bias:    {diff.mean():+.4f} dB")
    print(f"  Median bias:  {np.median(diff):+.4f} dB")
    print(f"  Std dev:      {diff.std():.4f} dB")
    print(f"  RMSE:         {np.sqrt((diff ** 2).mean()):.4f} dB")

    print()
    print("── |error| percentiles ─────────────────────────────────")
    for p in [50, 68, 90, 95, 99]:
        print(f"  p{p:2d}: |error| ≤ {np.percentile(np.abs(diff), p):.3f} dB")

    print()
    linear_sardine = 10.0 ** (s / 10.0)
    linear_mpc = 10.0 ** (m / 10.0)
    linear_mean_sardine = float(linear_sardine.mean())
    linear_mean_mpc = float(linear_mpc.mean())
    linear_bias_db = (
        10.0 * np.log10(linear_mean_sardine / linear_mean_mpc)
        if linear_mean_mpc > 0 else float("nan")
    )
    print("── Linear-domain (filter-invariant) ──────────────────────")
    print(f"  sardine linear mean:  {linear_mean_sardine:.6e}")
    print(f"  MPC     linear mean:  {linear_mean_mpc:.6e}")
    print(f"  linear mean bias:     {linear_bias_db:+.4f} dB")
    print()
    print("── Interpretation ────────────────────────────────────────")
    median_bias = float(np.median(diff))
    if abs(median_bias) <= 0.5:
        print(f"  PASS  median bias {median_bias:+.3f} dB is within ±0.5 dB")
    elif abs(median_bias) <= 1.5:
        print(f"  WARN  median bias {median_bias:+.3f} dB is within ±1.5 dB")
        print("        (acceptable given GRD vs SLC processing differences)")
    else:
        print(f"  FAIL  median bias {median_bias:+.3f} dB exceeds ±1.5 dB")
        print("        Investigate calibration LUT, terrain-flattening, or CRS.")


if __name__ == "__main__":
    main()
