"""
ASF Reference Validation Test

Compares SARdine RTC/backscatter output against an ASF RTC reference product.

Run only when both env vars are set:
  - SARDINE_TEST_SAFE_PATH: path to Sentinel-1 .SAFE (or .zip)
  - ASF_REF_PATH: path to ASF RTC directory or GeoTIFF

Optional:
  - SARDINE_TEST_ORBIT_CACHE: pre-downloaded orbit dir to avoid network

The test runs the pipeline, locates the SARdine GeoTIFF, reprojects the ASF
product to the same grid, converts to dB as needed, and asserts loose
thresholds on agreement. Intended as a regression/validation guard.
"""

from __future__ import annotations

import json
import os
from pathlib import Path

import numpy as np
import pytest

pytest.importorskip("sardine")
pytest.importorskip("rasterio")

from sardine.processors import BackscatterProcessor
import rasterio
from rasterio.enums import Resampling
from rasterio.warp import reproject


def _find_first_geotiff(path: Path) -> Path | None:
    if path.is_file() and path.suffix.lower() in (".tif", ".tiff"):
        return path
    # Prefer VV polarization if available
    for p in sorted(path.rglob("*_VV.tif")):
        return p
    for p in sorted(path.rglob("*_VV.tiff")):
        return p
    # Fallback to any TIF
    for p in sorted(path.rglob("*.tif")):
        return p
    for p in sorted(path.rglob("*.tiff")):
        return p
    return None


def _find_sardine_output_geotiff(out_dir: Path) -> Path | None:
    candidates = list(out_dir.rglob("*.tif")) or list(out_dir.rglob("*.tiff"))
    if not candidates:
        return None
    filtered = [p for p in candidates if "diff" not in p.name.lower()]
    if filtered:
        candidates = filtered
    ranked = sorted(
        candidates,
        key=lambda p: (
            -int(any(k in p.name.lower() for k in ("merged", "combined", "geocoded", "rtc", "gamma", "sigma"))),
            -p.stat().st_mtime,
        ),
    )
    return ranked[0]


def _reproject_to_match(src: rasterio.DatasetReader, dst_ref: rasterio.DatasetReader) -> np.ndarray:
    src_data = src.read(1, masked=True)
    src_nodata = src.nodata
    if src_nodata is None:
        # If the dataset doesn't declare nodata, fall back to 0.0 (common for many GeoTIFFs)
        # so that masked pixels don't get treated as valid backscatter.
        src_nodata = 0.0
    dst = np.full((dst_ref.height, dst_ref.width), np.nan, dtype=np.float32)
    reproject(
        source=np.asarray(src_data.filled(src_nodata), dtype=np.float32),
        destination=dst,
        src_transform=src.transform,
        src_crs=src.crs,
        dst_transform=dst_ref.transform,
        dst_crs=dst_ref.crs,
        resampling=Resampling.bilinear,
        src_nodata=src_nodata,
        dst_nodata=np.nan,
    )
    return dst


def _reproject_array_to_match(
    src_data: np.ndarray,
    *,
    src_transform: rasterio.Affine,
    src_crs: rasterio.crs.CRS,
    src_nodata: float | None,
    dst_ref: rasterio.DatasetReader,
) -> np.ndarray:
    dst = np.full((dst_ref.height, dst_ref.width), np.nan, dtype=np.float32)

    if src_nodata is None:
        src_nodata = 0.0

    src_arr = np.asarray(src_data, dtype=np.float32)
    if np.ma.isMaskedArray(src_data):
        src_arr = np.asarray(src_data.filled(src_nodata), dtype=np.float32)

    reproject(
        source=src_arr,
        destination=dst,
        src_transform=src_transform,
        src_crs=src_crs,
        dst_transform=dst_ref.transform,
        dst_crs=dst_ref.crs,
        resampling=Resampling.bilinear,
        src_nodata=src_nodata,
        dst_nodata=np.nan,
    )
    return dst


def _safe_to_db(arr: np.ndarray, eps: float = 1e-12) -> np.ndarray:
    arr = np.asarray(arr, dtype=np.float64)
    arr[arr <= 0] = np.nan
    return (10.0 * np.log10(arr + eps)).astype(np.float32)


def _safe_db_to_linear(arr_db: np.ndarray, floor_db: float = -79.0) -> np.ndarray:
    out = np.asarray(arr_db, dtype=np.float64)
    invalid = ~np.isfinite(out) | (out <= float(floor_db))
    out[invalid] = np.nan
    # Power linear = 10^(dB/10)
    return np.power(10.0, out / 10.0).astype(np.float32)


def _dataset_tags_indicate_db(ds: rasterio.DatasetReader) -> bool:
    tags = {**ds.tags(), **ds.tags(1)}
    for key, value in tags.items():
        if not isinstance(value, str):
            continue
        value_lower = value.lower()
        if any(token in value_lower for token in (" db", "db)", "[db]", "dbscale", "decibel")):
            return True
        if key.lower() in {"scale", "data_scale", "value_scale"} and value_lower in {"db", "decibel", "decibels"}:
            return True
    return False


def _looks_like_db_data(arr: np.ndarray) -> bool:
    data = np.asarray(arr, dtype=np.float64)
    data = data[np.isfinite(data)]
    if data.size == 0:
        return False
    neg_fraction = float(np.count_nonzero(data < 0.0)) / data.size
    high_fraction = float(np.count_nonzero(data > 100.0)) / data.size
    if neg_fraction > 0.1 and high_fraction == 0.0:
        return True
    finite_min = float(np.nanmin(data))
    finite_max = float(np.nanmax(data))
    return finite_min <= -70.0 and finite_max <= 40.0


def _to_db(arr: np.ndarray, ds: rasterio.DatasetReader) -> np.ndarray:
    if _dataset_tags_indicate_db(ds) or _looks_like_db_data(arr):
        out = np.asarray(arr, dtype=np.float32)
        floor = -79.0
        invalid = ~np.isfinite(out)
        invalid |= out <= floor
        out[invalid] = np.nan
        return out
    return _safe_to_db(arr)


def _to_linear(arr: np.ndarray, ds: rasterio.DatasetReader) -> np.ndarray:
    """Convert dataset values to linear power.

    We resample/reproject in linear space to avoid bilinear interpolation artifacts in dB.
    """

    if _dataset_tags_indicate_db(ds) or _looks_like_db_data(arr):
        return _safe_db_to_linear(arr)
    out = np.asarray(arr, dtype=np.float32)
    out[~np.isfinite(out) | (out <= 0.0)] = np.nan
    return out


def _diff_stats_db(sard: np.ndarray, ref: np.ndarray) -> dict:
    mask = np.isfinite(sard) & np.isfinite(ref)
    if not np.any(mask):
        return {"count": 0}
    d = sard[mask] - ref[mask]
    return {
        "count": int(mask.sum()),
        "bias_db_mean": float(np.nanmean(d)),
        "bias_db_median": float(np.nanmedian(d)),
        "bias_db_std": float(np.nanstd(d)),
        "p68_within_db": float(np.nanpercentile(np.abs(d), 68)),
        "p95_within_db": float(np.nanpercentile(np.abs(d), 95)),
        "pct_within_1db": float(np.mean(np.abs(d) <= 1.0) * 100.0),
        "pct_within_2db": float(np.mean(np.abs(d) <= 2.0) * 100.0),
        "pct_within_3db": float(np.mean(np.abs(d) <= 3.0) * 100.0),
    }


@pytest.mark.slow
@pytest.mark.integration
def test_compare_asf_reference(tmp_path: Path) -> None:
    safe_input = os.environ.get("SARDINE_TEST_SAFE_PATH")
    asf_ref = os.environ.get("ASF_REF_PATH")
    if not safe_input or not asf_ref:
        pytest.skip("SARDINE_TEST_SAFE_PATH and ASF_REF_PATH must be set to run ASF comparison")

    orbit_cache_override = os.environ.get("SARDINE_TEST_ORBIT_CACHE")
    if orbit_cache_override:
        os.environ["SARDINE_ORBIT_CACHE"] = orbit_cache_override

    # Run SARdine backscatter pipeline (RTC on by default)
    out_dir = tmp_path / "sardine_out"
    opts = {
        "verbose": False,
        "terrain_flatten": True,
        "geocode": True,
        "resolution": 30.0,
        "optimization_mode": "complete",
        "allow_synthetic": False,
        "use_real_orbit": True,
        "sequential": True,
    }
    proc = BackscatterProcessor(safe_input, out_dir, opts)
    proc.process_backscatter()

    sard_tif = _find_sardine_output_geotiff(out_dir)
    assert sard_tif is not None and sard_tif.exists(), "SARdine output GeoTIFF not found"

    asf_path = Path(asf_ref)
    asf_tif = _find_first_geotiff(asf_path)
    assert asf_tif is not None and asf_tif.exists(), "ASF reference GeoTIFF not found"

    # Read and compare in a consistent domain:
    # 1) Convert both datasets to linear power (handling dB products)
    # 2) Reproject ASF to SARdine grid in linear domain
    # 3) Convert both to dB for robust difference stats
    with rasterio.open(sard_tif) as ds_sard, rasterio.open(asf_tif) as ds_asf:
        sard_linear = _to_linear(ds_sard.read(1, masked=True), ds_sard)
        asf_linear_native = _to_linear(ds_asf.read(1, masked=True), ds_asf)

        # Reproject ASF in linear domain to match SARdine output grid
        asf_linear = _reproject_array_to_match(
            asf_linear_native,
            src_transform=ds_asf.transform,
            src_crs=ds_asf.crs,
            src_nodata=ds_asf.nodata,
            dst_ref=ds_sard,
        )

        sard_db = _safe_to_db(sard_linear)
        asf_db = _safe_to_db(asf_linear)

    stats = _diff_stats_db(sard_db, asf_db)

    # Loose acceptance thresholds (scene-dependent; this is a regression guard)
    assert stats["count"] > 1000, "Too few overlapping valid pixels"
    assert stats["p95_within_db"] <= 6.0, f"High 95th percentile error: {stats}"
    assert stats["pct_within_3db"] >= 50.0, f"Low fraction within 3 dB: {stats}"

    # Write summary for debugging
    summary = tmp_path / "asf_compare_summary.json"
    summary.write_text(json.dumps(stats, indent=2))
    print("ASF comparison summary:", stats)

