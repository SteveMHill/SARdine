#!/usr/bin/env python
"""
Minimal test runner for TOPSAR stages (metadata → orbit → IW split → deburst+cal → merge).

It executes the exact early pipeline steps, prints validation stats, and writes two GeoTIFFs
of the merged slant-range backscatter (linear and dB) so you can inspect the output.
"""

import argparse
import json
import time
import os
from pathlib import Path
from typing import Any, Dict, Optional, Tuple

import numpy as np

import sardine
from sardine.geotiff import export_geotiff
from sardine.processors.backscatter import merge_stage
from sardine.processors.backscatter.processor import BackscatterProcessor
from sardine.processors.pipeline import PipelineContext


def _run_stage(label, fn, context):
    start = time.time()
    fn(context)
    return time.time() - start


def _compute_geotransform(processor) -> Optional[Tuple[float, float, float, float, float, float]]:
    """Prefer merge-provided geotransform; fall back to pixel spacing."""
    if getattr(processor, "geo_transform", None):
        try:
            gt = tuple(float(v) for v in processor.geo_transform)
            if len(gt) == 6:
                return gt
        except Exception:
            pass

    rng = processor.get_current_range_spacing()
    azi = processor.get_current_azimuth_spacing()
    if rng and azi:
        # Placeholder transform in slant-range grid units (not geocoded)
        return (0.0, float(rng), 0.0, 0.0, 0.0, -float(azi))
    return None


def _extract_array(result: Any, label: str) -> np.ndarray:
    if isinstance(result, dict):
        data = result.get("data")
        if data is None:
            data = result.get("calibrated_data")
        if data is None:
            data = result.get("power_data")
        if data is None:
            raise RuntimeError(f"{label} did not return data (keys: {list(result.keys())})")
        return np.asarray(data)
    return np.asarray(result)


def _process_subswath_with_fallback(processor: BackscatterProcessor, subswath: str) -> Dict[str, Any]:
    """Deburst + calibrate a subswath using standard path (cached deburst + calibration)."""
    pol = processor.polarization

    # Step A: Deburst (complex) then cached calibration
    with processor.reader_lock:
        deburst_result = sardine.deburst_topsar_cached(processor.reader, subswath, pol)
    deburst_data = _extract_array(deburst_result, "deburst_topsar_cached").astype(np.complex64, copy=False)
    try:
        with processor.reader_lock:
            cal_result = sardine.radiometric_calibration_with_denoising_cached(
                processor.reader,
                subswath,
                pol,
                processor.calibration_type,
                deburst_data,
                processor.fast_mode_noise_removal,
            )
        calibrated = _extract_array(cal_result, "radiometric_calibration").astype(np.float32, copy=False)
        return {"data": calibrated, "deburst_shape": deburst_data.shape, "path": "cached_cal"}
    except Exception as primary_err:
        print(f"   ⚠️  Cached calibration failed on {subswath}: {primary_err}")

    # Step B: Last-resort power only (uncalibrated)
    power = np.abs(deburst_data) ** 2
    print(f"   ⚠️  Using uncalibrated power for {subswath}")
    return {"data": power.astype(np.float32), "deburst_shape": deburst_data.shape, "path": "power_only"}


def main():
    parser = argparse.ArgumentParser(description="Test TOPSAR stages and export merged GeoTIFFs")
    parser.add_argument("input_path", help="Sentinel-1 IW SLC ZIP or SAFE directory")
    parser.add_argument("output_dir", help="Directory for diagnostics and GeoTIFFs")
    parser.add_argument("--polarization", default="VV", help="Polarization (VV/VH/HH/HV)")
    parser.add_argument(
        "--calibration-type",
        default="gamma0",
        choices=["gamma0", "sigma0", "beta0", "dn"],
        help="Radiometric calibration type",
    )
    parser.add_argument(
        "--subswaths",
        help="Comma-separated subswaths to process (e.g., IW1,IW2). Default: all from metadata.",
    )
    parser.add_argument(
        "--noise-removal",
        action="store_true",
        help="Apply thermal noise removal during calibration",
    )
    parser.add_argument(
        "--power-only",
        action="store_true",
        help="Skip calibration entirely and use deburst power for quick inspection.",
    )
    args = parser.parse_args()

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Ensure orbit cache is set for the cached reader path
    sardine.setup_default_orbit_cache(output_dir)

    # Keep only the early stages; disable geocode/terrain/speckle for speed
    options = {
        "polarization": args.polarization.upper(),
        "calibration_type": args.calibration_type,
        "geocode": False,
        "terrain_flatten": False,
        "speckle_filter": "none",
        "multilook_range": 1,
        "multilook_azimuth": 1,
        "enable_parallel": True,
        "num_threads": os.cpu_count() or 4,
        "sequential": False,
        "quality_report": False,
        "memory_cleanup": False,
        "prefetch_dem": False,
        "fast_mode_noise_removal": args.noise_removal,
    }

    processor = BackscatterProcessor(args.input_path, output_dir, options)
    context = PipelineContext(processor=processor)

    timings = {}
    print("\n=== Running TOPSAR test stages ===")
    timings["metadata"] = _run_stage("metadata", processor._stage_ingest_metadata, context)
    timings["orbit"] = _run_stage("orbit", processor._stage_apply_precise_orbit, context)
    timings["iw_split"] = _run_stage("iw_split", processor._stage_iw_split, context)

    # Custom per-subswath processing with fallbacks
    meta_subswaths = []
    if isinstance(processor.metadata, dict):
        subswaths_str = processor.metadata.get("subswaths", "")
        meta_subswaths = [sw.strip() for sw in subswaths_str.split(",") if sw.strip()]
    if args.subswaths:
        subswaths = [sw.strip() for sw in args.subswaths.split(",") if sw.strip()]
    else:
        subswaths = context.get_artifact("subswaths") or meta_subswaths
    if not subswaths:
        raise RuntimeError("No subswaths discovered in metadata")
    print(f"   📡 Subswaths to process: {', '.join(subswaths)}")

    cal_subswaths: Dict[str, np.ndarray] = {}
    for sw in subswaths:
        t0 = time.time()
        if args.power_only:
            # Fast path: deburst only, use power
            with processor.reader_lock:
                deburst_result = sardine.deburst_topsar_cached(processor.reader, sw, processor.polarization)
            deburst_data = _extract_array(deburst_result, "deburst_topsar_cached").astype(np.complex64, copy=False)
            data = (np.abs(deburst_data) ** 2).astype(np.float32)
            result = {"data": data, "path": "power_only"}
        else:
            result = _process_subswath_with_fallback(processor, sw)
        cal_subswaths[sw] = result["data"]
        print(
            f"   ✅ {sw}: {result['data'].shape[0]}x{result['data'].shape[1]} via {result['path']} "
            f"({time.time()-t0:.1f}s)"
        )

    processor._calibrated_subswaths = cal_subswaths
    preferred_order = sorted(cal_subswaths.keys(), key=lambda s: (0, int(s[2:])) if s.upper().startswith("IW") else (1, s))
    processor._primary_subswath = preferred_order[0]
    processor._working_data = cal_subswaths[processor._primary_subswath]

    # Merge
    merge_start = time.time()
    merge_stage.perform_merge(processor)
    timings["merge"] = time.time() - merge_start
    merged = getattr(processor, "_working_data", None)
    if merged is None:
        raise RuntimeError("Merge produced no data; check earlier stage logs.")
    merged = np.asarray(merged)

    finite_mask = np.isfinite(merged)
    finite_pct = 100.0 * finite_mask.sum() / merged.size
    min_val = float(np.nanmin(merged)) if finite_mask.any() else float("nan")
    max_val = float(np.nanmax(merged)) if finite_mask.any() else float("nan")
    mean_val = float(np.nanmean(merged)) if finite_mask.any() else float("nan")

    stats = {
        "finite_pct": finite_pct,
        "min": min_val,
        "max": max_val,
        "mean": mean_val,
        "shape": merged.shape,
        "timings_sec": timings,
    }
    (output_dir / "topsar_stage_stats.json").write_text(json.dumps(stats, indent=2))
    print("\nStage stats:")
    print(json.dumps(stats, indent=2))

    geotransform = _compute_geotransform(processor)
    if geotransform is None:
        raise RuntimeError(
            "GeoTIFF export needs a geotransform or bounds. "
            "Merge did not provide one and pixel spacing is unknown."
        )

    # Export linear power and dB for quick inspection
    lin_tif = output_dir / f"topsar_merged_{args.polarization.lower()}_linear.tif"
    db_tif = output_dir / f"topsar_merged_{args.polarization.lower()}_db.tif"

    export_geotiff(
        merged.astype(np.float32),
        lin_tif,
        geotransform=geotransform,
        crs="EPSG:4326",
        nodata=np.nan,
        description="Merged TOPSAR backscatter (linear power)",
    )

    db_data = 10.0 * np.log10(np.clip(merged, 1e-10, None))
    export_geotiff(
        db_data.astype(np.float32),
        db_tif,
        geotransform=geotransform,
        crs="EPSG:4326",
        nodata=np.nan,
        description="Merged TOPSAR backscatter (dB)",
    )

    print(f"\nGeoTIFFs written:\n- {lin_tif}\n- {db_tif}")


if __name__ == "__main__":
    main()
