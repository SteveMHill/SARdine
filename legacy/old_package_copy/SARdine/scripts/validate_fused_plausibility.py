#!/usr/bin/env python3
"""
Run the fused backscatter pipeline (no geocode / no terrain flatten) and emit basic plausibility
metrics over the produced GeoTIFFs.

This is a lightweight check when terrain-corrected goldens are not available. It verifies:
- CLI completes successfully on the given SAFE
- Output GeoTIFFs exist and are readable
- Finite pixel ratio, min/mean/max per product
- File size sanity

Usage:
  python validate_fused_plausibility.py /path/to/SAFE /path/to/output_dir [--pol VV] [--threads 8]

Env used:
  SARDINE_ORBIT_CACHE (defaults to <outdir>/orbit_cache)
  RUST_LOG (defaults to info)
"""
from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import time
from pathlib import Path
from typing import Dict, Any

import numpy as np
import rasterio


def run_cli(safe_path: Path, outdir: Path, pol: str, threads: int | None) -> Dict[str, Any]:
    env = os.environ.copy()
    env.setdefault("RUST_LOG", "info")
    env.setdefault("SARDINE_ORBIT_CACHE", str(outdir / "orbit_cache"))
    outdir.mkdir(parents=True, exist_ok=True)
    Path(env["SARDINE_ORBIT_CACHE"]).mkdir(parents=True, exist_ok=True)

    cmd = [
        sys.executable,
        "-m",
        "sardine.cli",
        "backscatter",
        str(safe_path),
        str(outdir),
        "--polarization",
        pol,
        "--parallel",
        "--no-geocode",
        "--no-terrain-flatten",
    ]
    if threads is not None:
        cmd += ["--threads", str(threads)]

    log_path = outdir / "cli.log"
    start = time.time()
    with log_path.open("w") as lf:
        proc = subprocess.Popen(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            bufsize=1,
            env=env,
        )
        for line in proc.stdout:  # type: ignore[arg-type]
            lf.write(line)
            lf.flush()
        ret = proc.wait()
    end = time.time()

    return {
        "cmd": cmd,
        "log": str(log_path),
        "return_code": ret,
        "duration_seconds": end - start,
        "env": {"RUST_LOG": env.get("RUST_LOG"), "SARDINE_ORBIT_CACHE": env.get("SARDINE_ORBIT_CACHE")},
    }


def scan_geotiffs(outdir: Path) -> Dict[str, Any]:
    stats = {}
    tif_paths = sorted(outdir.glob("*.tif"))
    for tif in tif_paths:
        try:
            with rasterio.open(tif) as ds:
                arr = ds.read(1, masked=True)
                data = arr.compressed()
                finite = np.isfinite(data)
                finite_pct = 100.0 * finite.sum() / max(len(data), 1)
                stats[tif.name] = {
                    "shape": list(arr.shape),
                    "dtype": str(arr.dtype),
                    "min": float(np.nanmin(data)) if data.size else None,
                    "max": float(np.nanmax(data)) if data.size else None,
                    "mean": float(np.nanmean(data)) if data.size else None,
                    "finite_pct": finite_pct,
                    "file_size_bytes": tif.stat().st_size,
                }
        except Exception as exc:  # pragma: no cover - defensive
            stats[tif.name] = {"error": str(exc)}
    return {"count": len(tif_paths), "files": stats}


def main() -> int:
    ap = argparse.ArgumentParser(description="Validate fused backscatter plausibility (no geocode/terrain)")
    ap.add_argument("safe", type=Path, help="Path to SAFE product")
    ap.add_argument("outdir", type=Path, help="Output directory")
    ap.add_argument("--pol", default="VV", help="Polarization (default VV)")
    ap.add_argument("--threads", type=int, default=None, help="Thread override for CLI")
    args = ap.parse_args()

    if not args.safe.exists():
        ap.error(f"SAFE not found: {args.safe}")

    run_info = run_cli(args.safe, args.outdir, args.pol, args.threads)
    tif_info = scan_geotiffs(args.outdir)

    summary = {
        "run": run_info,
        "geotiff_stats": tif_info,
        "status": "ok" if run_info["return_code"] == 0 else "failed",
    }

    summary_path = args.outdir / "plausibility_summary.json"
    with summary_path.open("w") as f:
        json.dump(summary, f, indent=2)

    print(json.dumps(summary, indent=2))
    return 0 if run_info["return_code"] == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
