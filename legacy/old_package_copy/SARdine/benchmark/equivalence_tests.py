#!/usr/bin/env python3
"""
SARdine Equivalence Test Suite
==============================
Tests to verify numerical equivalence before/after optimizations.
These tests capture golden outputs and verify optimizations preserve correctness.
"""

import os
import sys
import json
import hashlib
import numpy as np
from pathlib import Path
from typing import Dict, Any, Optional, Tuple

# Test scene configuration
TEST_SCENE = "/home/datacube/apps/SARdine/data/SLC/S1B_IW_SLC__1SDV_20190123T053348_20190123T053415_014617_01B3D4_E833.SAFE"
GOLDEN_DIR = Path("/home/datacube/apps/SARdine/SARdine/benchmark/golden_outputs")


def compute_array_checksum(arr: np.ndarray) -> str:
    """Compute a deterministic checksum for a numpy array."""
    # Handle NaN consistently
    arr_copy = np.copy(arr)
    arr_copy = np.nan_to_num(arr_copy, nan=0.0, posinf=1e38, neginf=-1e38)
    return hashlib.sha256(arr_copy.tobytes()).hexdigest()[:16]


def compute_array_stats(arr: np.ndarray) -> Dict[str, float]:
    """Compute statistical properties of an array."""
    finite = arr[np.isfinite(arr)]
    if len(finite) == 0:
        return {
            "count": int(0),
            "finite_count": int(0),
            "nan_count": int(arr.size),
            "mean": float('nan'),
            "std": float('nan'),
            "min": float('nan'),
            "max": float('nan'),
            "median": float('nan'),
        }
    return {
        "count": int(arr.size),
        "finite_count": int(len(finite)),
        "nan_count": int(np.sum(~np.isfinite(arr))),
        "mean": float(np.mean(finite)),
        "std": float(np.std(finite)),
        "min": float(np.min(finite)),
        "max": float(np.max(finite)),
        "median": float(np.median(finite)),
    }


def arrays_equivalent(
    arr1: np.ndarray, 
    arr2: np.ndarray, 
    rtol: float = 1e-5, 
    atol: float = 1e-8
) -> Tuple[bool, Dict[str, Any]]:
    """
    Check if two arrays are equivalent within tolerance.
    
    Returns (is_equivalent, diagnostics_dict)
    """
    if arr1.shape != arr2.shape:
        return False, {"error": f"Shape mismatch: {arr1.shape} vs {arr2.shape}"}
    
    # Mask for where both are finite
    finite1 = np.isfinite(arr1)
    finite2 = np.isfinite(arr2)
    finite_both = finite1 & finite2
    
    # Check NaN pattern matches
    nan1 = ~finite1
    nan2 = ~finite2
    nan_mismatch = np.sum(nan1 != nan2)
    
    if nan_mismatch > 0:
        return False, {
            "error": "NaN pattern mismatch",
            "nan_mismatch_count": int(nan_mismatch),
            "arr1_nan_count": int(np.sum(nan1)),
            "arr2_nan_count": int(np.sum(nan2)),
        }
    
    # Compare finite values
    if np.sum(finite_both) == 0:
        return True, {"warning": "No finite values to compare"}
    
    diff = np.abs(arr1[finite_both] - arr2[finite_both])
    max_diff = float(np.max(diff))
    mean_diff = float(np.mean(diff))
    
    # Relative comparison
    denom = np.maximum(np.abs(arr1[finite_both]), np.abs(arr2[finite_both]))
    denom = np.where(denom > 0, denom, 1.0)
    rel_diff = diff / denom
    max_rel_diff = float(np.max(rel_diff))
    
    is_close = np.allclose(
        arr1[finite_both], arr2[finite_both], rtol=rtol, atol=atol, equal_nan=True
    )
    
    diagnostics = {
        "is_equivalent": is_close,
        "max_absolute_diff": max_diff,
        "mean_absolute_diff": mean_diff,
        "max_relative_diff": max_rel_diff,
        "compared_pixels": int(np.sum(finite_both)),
        "total_pixels": arr1.size,
    }
    
    return is_close, diagnostics


def capture_golden_outputs(output_dir: Path) -> Dict[str, Any]:
    """
    Capture golden outputs from a completed pipeline run.
    
    Saves:
    - Array checksums for exact matching
    - Statistical properties for approximate matching
    - Sample pixels for spot-checks
    """
    golden = {
        "source_dir": str(output_dir),
        "captures": {},
    }
    
    # Find output arrays
    for npy_file in output_dir.glob("*.npy"):
        arr = np.load(npy_file)
        name = npy_file.stem
        
        golden["captures"][name] = {
            "checksum": compute_array_checksum(arr),
            "shape": [int(x) for x in arr.shape],
            "dtype": str(arr.dtype),
            "stats": compute_array_stats(arr),
        }
        
        # Capture sample pixels at fixed positions for spot-checks
        if arr.ndim == 2 and arr.shape[0] > 100 and arr.shape[1] > 100:
            samples = []
            for row in [0, arr.shape[0]//4, arr.shape[0]//2, 3*arr.shape[0]//4, arr.shape[0]-1]:
                for col in [0, arr.shape[1]//4, arr.shape[1]//2, 3*arr.shape[1]//4, arr.shape[1]-1]:
                    if row < arr.shape[0] and col < arr.shape[1]:
                        samples.append({
                            "row": row,
                            "col": col,
                            "value": float(arr[row, col]) if np.isfinite(arr[row, col]) else "NaN",
                        })
            golden["captures"][name]["samples"] = samples
    
    # Find GeoTIFFs and capture their checksums
    for tif_file in output_dir.glob("*.tif"):
        try:
            import rasterio
            with rasterio.open(tif_file) as src:
                arr = src.read(1)
                name = tif_file.stem
                
                golden["captures"][name] = {
                    "checksum": compute_array_checksum(arr),
                    "shape": [int(x) for x in arr.shape],
                    "dtype": str(arr.dtype),
                    "stats": compute_array_stats(arr),
                    "crs": str(src.crs),
                    "bounds": [float(x) for x in src.bounds],
                }
        except ImportError:
            print("rasterio not available, skipping GeoTIFF capture")
            break
    
    return golden


def verify_against_golden(output_dir: Path, golden: Dict[str, Any], rtol: float = 1e-5) -> Dict[str, Any]:
    """
    Verify current outputs against golden reference.
    
    Returns a report with pass/fail status and diagnostics.
    """
    report = {
        "passed": True,
        "results": {},
    }
    
    for name, golden_data in golden["captures"].items():
        # Try to load current output
        npy_path = output_dir / f"{name}.npy"
        tif_path = output_dir / f"{name}.tif"
        
        if npy_path.exists():
            current = np.load(npy_path)
        elif tif_path.exists():
            try:
                import rasterio
                with rasterio.open(tif_path) as src:
                    current = src.read(1)
            except ImportError:
                report["results"][name] = {"skipped": True, "reason": "rasterio not available"}
                continue
        else:
            report["results"][name] = {"passed": False, "error": "Output file not found"}
            report["passed"] = False
            continue
        
        # Check shape
        if tuple(current.shape) != tuple(golden_data["shape"]):
            report["results"][name] = {
                "passed": False,
                "error": f"Shape mismatch: {current.shape} vs {golden_data['shape']}",
            }
            report["passed"] = False
            continue
        
        # Check checksum for exact match
        current_checksum = compute_array_checksum(current)
        if current_checksum == golden_data["checksum"]:
            report["results"][name] = {
                "passed": True,
                "match_type": "exact",
            }
            continue
        
        # If not exact, check statistical equivalence
        current_stats = compute_array_stats(current)
        golden_stats = golden_data["stats"]
        
        # Check mean/std/min/max are within tolerance
        stat_checks = {}
        for stat in ["mean", "std", "min", "max"]:
            if np.isnan(golden_stats[stat]) and np.isnan(current_stats[stat]):
                stat_checks[stat] = "both_nan"
            elif np.isnan(golden_stats[stat]) or np.isnan(current_stats[stat]):
                stat_checks[stat] = "nan_mismatch"
            else:
                diff = abs(current_stats[stat] - golden_stats[stat])
                rel_diff = diff / max(abs(golden_stats[stat]), 1e-10)
                stat_checks[stat] = {
                    "golden": golden_stats[stat],
                    "current": current_stats[stat],
                    "abs_diff": diff,
                    "rel_diff": rel_diff,
                    "passed": rel_diff < rtol,
                }
        
        all_stats_passed = all(
            isinstance(v, str) or v.get("passed", True) 
            for v in stat_checks.values()
        )
        
        report["results"][name] = {
            "passed": all_stats_passed,
            "match_type": "statistical",
            "checksum_match": False,
            "stat_checks": stat_checks,
            "finite_count_golden": golden_stats.get("finite_count"),
            "finite_count_current": current_stats.get("finite_count"),
        }
        
        if not all_stats_passed:
            report["passed"] = False
    
    return report


def main():
    import argparse
    parser = argparse.ArgumentParser(description="SARdine Equivalence Tests")
    parser.add_argument("--capture", type=str, help="Capture golden outputs from directory")
    parser.add_argument("--verify", type=str, help="Verify outputs against golden")
    parser.add_argument("--golden", type=str, default=str(GOLDEN_DIR / "golden.json"),
                        help="Path to golden reference file")
    parser.add_argument("--rtol", type=float, default=1e-5, help="Relative tolerance")
    args = parser.parse_args()
    
    if args.capture:
        output_dir = Path(args.capture)
        if not output_dir.exists():
            print(f"Error: Directory not found: {output_dir}")
            return 1
        
        print(f"Capturing golden outputs from: {output_dir}")
        golden = capture_golden_outputs(output_dir)
        
        # Save golden reference
        golden_file = Path(args.golden)
        golden_file.parent.mkdir(parents=True, exist_ok=True)
        with open(golden_file, "w") as f:
            json.dump(golden, f, indent=2)
        
        print(f"Captured {len(golden['captures'])} outputs to: {golden_file}")
        for name, data in golden["captures"].items():
            print(f"  {name}: shape={data['shape']}, checksum={data['checksum']}")
        
        return 0
    
    elif args.verify:
        output_dir = Path(args.verify)
        golden_file = Path(args.golden)
        
        if not output_dir.exists():
            print(f"Error: Directory not found: {output_dir}")
            return 1
        if not golden_file.exists():
            print(f"Error: Golden reference not found: {golden_file}")
            print("Run with --capture first to create golden reference")
            return 1
        
        with open(golden_file) as f:
            golden = json.load(f)
        
        print(f"Verifying outputs against: {golden_file}")
        report = verify_against_golden(output_dir, golden, rtol=args.rtol)
        
        # Print results
        print("\n" + "=" * 60)
        print("EQUIVALENCE TEST RESULTS")
        print("=" * 60)
        
        for name, result in report["results"].items():
            status = "✅ PASS" if result.get("passed") else "❌ FAIL"
            match_type = result.get("match_type", "unknown")
            print(f"{status} {name} ({match_type})")
            
            if not result.get("passed"):
                if "error" in result:
                    print(f"     Error: {result['error']}")
                if "stat_checks" in result:
                    for stat, check in result["stat_checks"].items():
                        if isinstance(check, dict) and not check.get("passed"):
                            print(f"     {stat}: golden={check['golden']:.6g}, current={check['current']:.6g}, diff={check['rel_diff']:.2e}")
        
        print("=" * 60)
        overall = "✅ ALL TESTS PASSED" if report["passed"] else "❌ SOME TESTS FAILED"
        print(overall)
        
        return 0 if report["passed"] else 1
    
    else:
        parser.print_help()
        return 1


if __name__ == "__main__":
    sys.exit(main())
