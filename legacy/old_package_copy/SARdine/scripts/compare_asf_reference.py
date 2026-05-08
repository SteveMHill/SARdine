#!/usr/bin/env python3
"""
Compare SARdine RTC/backscatter output against an ASF Vertex RTC reference.

Workflow
- Run SARdine backscatter pipeline (geocode + terrain flattening ON) for a given SAFE
- Locate SARdine merged GeoTIFF
- Locate ASF RTC GeoTIFF within a provided directory
- Reproject/resample reference to match SARdine grid
- Compute dB-domain stats and emit a JSON summary (+ optional diff GeoTIFF)

Requirements
- Installed SARdine package with CLI available (python -m sardine.cli backscatter)
- rasterio, numpy, tqdm (bundled in the Python deps per repo docs)

Usage
  python scripts/compare_asf_reference.py \
      --safe /path/to/SAFE \
      --asf-ref /path/to/ASF_RTC_dir_or_tif \
      --out /path/to/output_dir \
      [--orbit-cache /path/to/orbit_cache] \
      [--keep-sardine-output] \
      [--save-diff]

Notes
- Set SARDINE_ORBIT_CACHE env var instead of --orbit-cache if preferred
- Script sets strict parsing env flags automatically unless overridden
"""

from __future__ import annotations
import argparse
import json
import os
import sys
import subprocess
from pathlib import Path
from typing import Optional, Tuple

import numpy as np
import rasterio
from rasterio.enums import Resampling
from rasterio.warp import reproject


def find_first_geotiff(path: Path) -> Optional[Path]:
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


def find_sardine_output_geotiff(out_dir: Path) -> Optional[Path]:
    # Look for a plausible merged geotiff name under output dir
    # Try common patterns; else, fallback to newest .tif in tree
    candidates = list(out_dir.rglob("*.tif"))
    if not candidates:
        candidates = list(out_dir.rglob("*.tiff"))
    if not candidates:
        return None
    filtered = [p for p in candidates if "diff" not in p.name.lower()]
    if filtered:
        candidates = filtered
    # Prefer names containing merged/calibrated/geocoded hints
    ranked = sorted(
        candidates,
        key=lambda p: (
            -int(any(k in p.name.lower() for k in ["merged", "combined", "geocoded", "rtc", "gamma", "sigma"])) ,
            -p.stat().st_mtime,
        ),
    )
    return ranked[0]


def run_sardine_pipeline(safe_path: Path, out_dir: Path, orbit_cache: Optional[Path]) -> None:
    env = os.environ.copy()
    env.setdefault("SARDINE_SERDE_ONLY", "1")
    env.setdefault("SARDINE_REQUIRE_SUBSWATHS", "1")
    if orbit_cache is not None:
        env["SARDINE_ORBIT_CACHE"] = str(orbit_cache)

    out_dir.mkdir(parents=True, exist_ok=True)

    # Run CLI backscatter with geocoding & terrain flattening enabled (defaults assumed ON)
    cmd = [
        sys.executable,
        "-m",
        "sardine.cli",
        "backscatter",
        str(safe_path),
        str(out_dir),
        "--sequential",
        "--verbose",
    ]
    print("Running:", " ".join(cmd))
    res = subprocess.run(cmd, env=env, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    print(res.stdout)
    if res.returncode != 0:
        raise SystemExit(f"Sardine pipeline failed with exit code {res.returncode}")


def reproject_to_match(src: rasterio.DatasetReader, dst_ref: rasterio.DatasetReader) -> np.ndarray:
    # Reproject src data to match dst_ref grid
    src_data = src.read(1, masked=True)
    dst = np.full((dst_ref.height, dst_ref.width), np.nan, dtype=np.float32)

    reproject(
        source=src_data,
        destination=dst,
        src_transform=src.transform,
        src_crs=src.crs,
        dst_transform=dst_ref.transform,
        dst_crs=dst_ref.crs,
        resampling=Resampling.bilinear,
        src_nodata=0.0,
        dst_nodata=np.nan,
    )
    return dst


def safe_to_db(arr: np.ndarray, eps: float = 1e-10) -> np.ndarray:
    # Clip to positive and convert power to dB; if amplitude sneaks in, relative diffs still informative
    if isinstance(arr, np.ma.MaskedArray):
        arr = arr.filled(np.nan)
    arr = np.asarray(arr, dtype=np.float64)
    arr[arr <= 0] = np.nan
    return 10.0 * np.log10(arr + eps)


def _dataset_tags_indicate_db(ds: rasterio.DatasetReader) -> bool:
    """Check dataset-level metadata for hints that values are already in dB."""
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


def _looks_like_db_data(arr: np.ma.MaskedArray) -> bool:
    """Heuristic detection of dB-scaled data when explicit metadata is absent."""
    if isinstance(arr, np.ma.MaskedArray):
        data = arr.compressed().astype(np.float64)
    else:
        data = np.asarray(arr, dtype=np.float64)
    data = data[np.isfinite(data)]
    if data.size == 0:
        return False

    neg_fraction = float(np.count_nonzero(data < 0.0)) / data.size
    high_fraction = float(np.count_nonzero(data > 100.0)) / data.size

    # Typical SAR dB backscatter has a majority of negative values and rarely exceeds 100 dB.
    if neg_fraction > 0.1 and high_fraction == 0.0:
        return True

    # Values tightly bounded around a common floor (e.g., -80 dB) likely indicate already-logged data.
    finite_min = float(np.nanmin(data))
    finite_max = float(np.nanmax(data))
    if finite_min <= -70.0 and finite_max <= 40.0:
        return True

    return False


def to_db(arr: np.ma.MaskedArray, ds: rasterio.DatasetReader) -> np.ndarray:
    """Convert input data to decibels, respecting inputs already in dB."""
    if _dataset_tags_indicate_db(ds) or _looks_like_db_data(arr):
        if isinstance(arr, np.ma.MaskedArray):
            db_array = arr.filled(np.nan).astype(np.float32)
        else:
            db_array = np.array(arr, dtype=np.float32)
        # Treat common floor values (e.g., -80 dB) as invalid
        floor = -79.0
        invalid = ~np.isfinite(db_array)
        invalid |= db_array <= floor
        db_array = db_array.astype(np.float32)
        db_array[invalid] = np.nan
        return db_array
    return safe_to_db(arr)


def diff_stats_db(sard: np.ndarray, ref: np.ndarray) -> dict:
    mask = np.isfinite(sard) & np.isfinite(ref)
    if not np.any(mask):
        return {"count": 0}
    d = sard[mask] - ref[mask]
    def pct_within(th):
        return float(np.mean(np.abs(d) <= th) * 100.0)
    return {
        "count": int(mask.sum()),
        "bias_db_mean": float(np.nanmean(d)),
        "bias_db_median": float(np.nanmedian(d)),
        "bias_db_std": float(np.nanstd(d)),
        "p68_within_db": float(np.nanpercentile(np.abs(d), 68)),
        "p95_within_db": float(np.nanpercentile(np.abs(d), 95)),
        "pct_within_1db": pct_within(1.0),
        "pct_within_2db": pct_within(2.0),
        "pct_within_3db": pct_within(3.0),
    }


def main():
    ap = argparse.ArgumentParser(description="Compare SARdine RTC/backscatter against ASF RTC reference")
    ap.add_argument("--safe", required=True, type=Path, help="Path to Sentinel-1 SAFE directory")
    ap.add_argument("--asf-ref", required=True, type=Path, help="Path to ASF RTC directory or GeoTIFF")
    ap.add_argument("--out", required=True, type=Path, help="Output directory for SARdine outputs and comparison")
    ap.add_argument("--orbit-cache", type=Path, default=None, help="Orbit cache directory (overrides SARDINE_ORBIT_CACHE)")
    ap.add_argument("--keep-sardine-output", action="store_true", help="Keep existing SARdine output; skip running pipeline if a GeoTIFF is found")
    ap.add_argument("--save-diff", action="store_true", help="Write GeoTIFF of dB difference (Sardine - ASF)")
    args = ap.parse_args()

    safe = args.safe
    asf_path = args.asf_ref
    out_dir = args.out
    orbit_cache = args.orbit_cache

    # 1) Run SARdine pipeline unless output already present and keep flag set
    out_dir.mkdir(parents=True, exist_ok=True)

    sard_tif = find_sardine_output_geotiff(out_dir)
    if sard_tif is None or not args.keep_sardine_output:
        run_sardine_pipeline(safe, out_dir, orbit_cache)
        sard_tif = find_sardine_output_geotiff(out_dir)

    if sard_tif is None:
        raise SystemExit("Could not locate SARdine output GeoTIFF after pipeline run")

    # 2) Locate ASF RTC GeoTIFF
    asf_tif = find_first_geotiff(asf_path)
    if asf_tif is None:
        raise SystemExit("Could not locate ASF RTC GeoTIFF in provided path")

    # 3) Open datasets
    with rasterio.open(sard_tif) as ds_sard, rasterio.open(asf_tif) as ds_asf:
        # Reproject ASF to match SARdine grid
        asf_on_sard = reproject_to_match(ds_asf, ds_sard)
        sard_data = ds_sard.read(1, masked=True).astype(np.float32)

        # Convert to dB domain
        sard_db = to_db(sard_data, ds_sard)
        asf_db = safe_to_db(asf_on_sard)

        # Compute stats
        stats = diff_stats_db(sard_db, asf_db)
        stats.update({
            "sardine_path": str(sard_tif),
            "asf_path": str(asf_tif),
            "sardine_crs": ds_sard.crs.to_string() if ds_sard.crs else None,
            "asf_crs": ds_asf.crs.to_string() if ds_asf.crs else None,
            "sardine_transform": tuple(ds_sard.transform) if ds_sard.transform else None,
            "asf_transform": tuple(ds_asf.transform) if ds_asf.transform else None,
            "width": ds_sard.width,
            "height": ds_sard.height,
        })

        summary_path = out_dir / "comparison_summary.json"
        with open(summary_path, "w") as f:
            json.dump(stats, f, indent=2)
        print("Wrote:", summary_path)

        # Optional: save dB difference map
        if args.save_diff and stats.get("count", 0) > 0:
            diff_db = sard_db - asf_db
            diff_path = out_dir / "diff_sardine_minus_asf_db.tif"
            profile = ds_sard.profile
            profile.update(dtype=rasterio.float32, count=1, compress="deflate")
            with rasterio.open(diff_path, "w", **profile) as dst:
                dst.write(diff_db.astype(np.float32), 1)
            print("Wrote:", diff_path)

    print("Done.")


if __name__ == "__main__":
    main()
