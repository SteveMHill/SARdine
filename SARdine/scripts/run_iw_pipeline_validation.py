#!/usr/bin/env python3
"""End-to-end IW pipeline validation harness.

This script drives the critical Sentinel-1 IW processing stages in sequence:

1. Read metadata & cached annotations
2. Apply precise orbit file (POEORB)
3. IW split geometry discovery (via cached slices)
4. Deburst bursts to continuous azimuth lines
5. Radiometric calibration (sigma0/beta0/gamma0)
6. Merge calibrated IW rasters and verify geotransforms

It delegates the heavy lifting to the high-level helpers in ``sardine`` and
``sardine.partial_pipeline`` so the execution mirrors the Python CLI workflow.

Usage
-----

```
python3 scripts/run_iw_pipeline_validation.py \
    --input /path/to/S1?_IW_SLC__*.SAFE \
    --output ./pipeline_output \
    --polarization VV \
    --calibration-type sigma0
```

The script expects that the orbit cache directory either already contains the
relevant precise orbit file or that the environment has network access for the
Rust orbit manager to download it. The cache location defaults to
``<output>/orbit_cache`` and can be overridden with ``--orbit-cache``.
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path
from typing import Dict, Iterable, Mapping, Tuple

import sardine
from sardine.partial_pipeline import export_calibrated_and_merged_geotiffs


def _extract_first(metadata: Mapping[str, str], *keys: Iterable[str]) -> str:
    """Return the first non-empty metadata value among the provided keys."""

    for key in keys:
        value = metadata.get(key)
        if value:
            return value
    raise KeyError(f"None of the metadata keys {keys!r} were present")


def _apply_precise_orbit(
    metadata: Mapping[str, str],
    orbit_cache: Path,
) -> Dict[str, object]:
    """Invoke the Rust precise orbit loader using metadata-derived fields."""

    product_id = _extract_first(
        metadata,
        "product_id",
        "safe_product_id",
        "safe_id",
        "mission_data_take_id",
    )
    start_time = _extract_first(
        metadata,
        "start_time",
        "start_time_iso",
        "start_time_utc",
        "sensing_start",
    )

    orbit_cache.mkdir(parents=True, exist_ok=True)

    print("🛰️  [Step 2/6] Applying precise orbit file…")
    response = sardine.apply_precise_orbit_file(
        product_id,
        start_time,
        str(orbit_cache),
    )
    orbit_info = response.get("result", {}) if isinstance(response, dict) else {}
    print(
        "   • orbit vectors:",
        orbit_info.get("orbit_vectors_count", "unknown"),
    )
    print("   • reference time:", orbit_info.get("reference_time", "unknown"))
    return orbit_info


def main(argv: Iterable[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Run the core IW calibration + merge pipeline with validation",
    )
    parser.add_argument("--input", required=True, help="Path to Sentinel-1 SLC (SAFE or ZIP)")
    parser.add_argument("--output", required=True, help="Directory for intermediate GeoTIFFs")
    parser.add_argument(
        "--polarization",
        default="VV",
        choices=["VV", "VH", "HH", "HV"],
        help="Polarization to process",
    )
    parser.add_argument(
        "--calibration-type",
        default="sigma0",
        choices=["sigma0", "beta0", "gamma0", "dn"],
        help="Radiometric calibration output type",
    )
    parser.add_argument(
        "--orbit-cache",
        help="Explicit orbit cache directory (default: <output>/orbit_cache)",
    )

    args = parser.parse_args(argv)

    input_path = Path(args.input).expanduser().resolve()
    output_dir = Path(args.output).expanduser().resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    orbit_cache = (
        Path(args.orbit_cache).expanduser().resolve()
        if args.orbit_cache
        else output_dir / "orbit_cache"
    )

    if not input_path.exists():
        print(f"❌ Input product not found: {input_path}", file=sys.stderr)
        return 1

    print("📦 [Step 1/6] Initialising cached reader and metadata lookup…")
    reader = sardine.create_cached_slc_reader(str(input_path))
    metadata = reader.get_cached_metadata()

    product_name = metadata.get("product_id", input_path.stem)
    sensing_start = metadata.get("start_time") or metadata.get("start_time_utc")
    print(f"   • product: {product_name}")
    if sensing_start:
        print(f"   • sensing start: {sensing_start}")
    print(f"   • polarization: {args.polarization}")

    orbit_info = _apply_precise_orbit(metadata, orbit_cache)

    print("🪬 [Step 3-6] Running IW deburst, calibration, and merge…")
    start = time.perf_counter()
    summary = export_calibrated_and_merged_geotiffs(
        input_path=input_path,
        output_dir=output_dir,
        polarization=args.polarization,
        calibration_type=args.calibration_type,
        export_subswaths=True,
        export_merged=True,
    )
    elapsed = time.perf_counter() - start

    print(f"   • completed in {elapsed:.1f}s")
    print("   • exported subswaths:")
    for record in summary["subswath_exports"]:
        print(
            f"     - {record['subswath']}: {record['path']} "
            f"({record['rows']}×{record['cols']})",
        )
    if summary.get("merged_raster"):
        merged = summary["merged_raster"]
        print(
            "   • merged raster:",
            merged.get("path"),
            f"({merged.get('rows')}×{merged.get('cols')})",
        )
    print("   • geotransforms (per subswath):")
    for swath, transform in summary["subswath_geotransforms"].items():
        print(f"     - {swath}: {transform}")

    summary_path = output_dir / f"{product_name}_validation_summary.json"
    payload = {
        "metadata": metadata,
        "orbit_info": orbit_info,
        "pipeline_summary": summary,
        "elapsed_seconds": elapsed,
    }
    with summary_path.open("w", encoding="utf-8") as fp:
        json.dump(payload, fp, indent=2)
    print(f"📝 Summary written to {summary_path}")

    print("✅ IW pipeline validation run finished.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
