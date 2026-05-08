#!/usr/bin/env python3
"""
SARdine Baseline Benchmark Script
=================================
Establishes reproducible performance baseline for the backscatter pipeline.

Measures:
- Total wall-clock time
- Per-stage breakdown (from pipeline step announcements)
- Peak memory usage (via /proc/self/status on Linux)
- I/O time vs. compute time where distinguishable
"""

import os
import sys
import time
import json
import subprocess
import re
import gc
import resource
from pathlib import Path
from datetime import datetime

# Fixed benchmark parameters
BENCHMARK_CONFIG = {
    "scene": "/home/datacube/apps/SARdine/data/SLC/S1B_IW_SLC__1SDV_20190123T053348_20190123T053415_014617_01B3D4_E833.SAFE",
    "polarization": "VV",
    "multilook": (2, 2),
    "threads": 8,  # Fixed thread count for reproducibility
    "geocode": False,  # Disable for faster benchmark
    "terrain_flatten": False,  # Disable for faster benchmark
    "resolution": 10.0,
}


def get_machine_info() -> dict:
    """Collect machine information for benchmark documentation."""
    info = {}
    
    # CPU info
    try:
        with open("/proc/cpuinfo", "r") as f:
            cpuinfo = f.read()
        model = re.search(r"model name\s*:\s*(.+)", cpuinfo)
        info["cpu_model"] = model.group(1).strip() if model else "Unknown"
        info["cpu_cores"] = cpuinfo.count("processor")
    except Exception as e:
        info["cpu_model"] = f"Error: {e}"
        info["cpu_cores"] = os.cpu_count() or 0
    
    # Memory info
    try:
        with open("/proc/meminfo", "r") as f:
            meminfo = f.read()
        total = re.search(r"MemTotal:\s*(\d+)\s*kB", meminfo)
        info["ram_gb"] = int(total.group(1)) / (1024 * 1024) if total else 0
    except Exception as e:
        info["ram_gb"] = f"Error: {e}"
    
    # OS info
    try:
        import platform
        info["os"] = platform.system()
        info["os_release"] = platform.release()
        info["python_version"] = platform.python_version()
    except Exception as e:
        info["os"] = f"Error: {e}"
    
    # Rust/SARdine version
    try:
        import sardine
        info["sardine_version"] = getattr(sardine, "__version__", "unknown")
    except Exception as e:
        info["sardine_version"] = f"Error: {e}"
    
    return info


def get_peak_memory_mb() -> float:
    """Get peak memory usage in MB (Linux only)."""
    try:
        usage = resource.getrusage(resource.RUSAGE_SELF)
        return usage.ru_maxrss / 1024  # Convert KB to MB
    except Exception:
        return 0.0


def parse_step_timings(log_lines: list) -> dict:
    """Parse pipeline step timings from log output.
    
    Looks for patterns like:
    - "✅ STEP 1: Read Metadata & Files (0.1s)"
    - "✅ STEP 4: Deburst & Radiometric Calibration (122.9s)"
    """
    step_times = {}
    
    # Pattern: "✅ STEP N: Name (X.Xs)"
    step_pattern = re.compile(r"✅ STEP (\d+): (.+?) \(([\d.]+)s\)")
    
    for line in log_lines:
        match = step_pattern.search(line)
        if match:
            step_num = int(match.group(1))
            step_name = match.group(2).strip()
            duration = float(match.group(3))
            step_times[step_num] = {
                "name": step_name,
                "duration_s": duration,
            }
    
    return step_times


def run_benchmark(output_dir: Path, config: dict, dry_run: bool = False) -> dict:
    """Run the benchmark and collect metrics."""
    
    scene_path = Path(config["scene"])
    if not scene_path.exists():
        raise FileNotFoundError(f"Benchmark scene not found: {scene_path}")
    
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Build command
    cmd = [
        sys.executable, "-m", "sardine.cli", "backscatter",
        str(scene_path),
        str(output_dir / "bench_output"),
        "--polarization", config["polarization"],
        "--multilook", str(config["multilook"][0]), str(config["multilook"][1]),
    ]
    
    # Add threading control
    if config.get("threads"):
        cmd.extend(["--num-threads", str(config["threads"])])
    
    # Processing toggles
    if not config.get("geocode", True):
        cmd.append("--no-geocode")
    if not config.get("terrain_flatten", True):
        cmd.append("--no-terrain-flatten")
    
    # Resolution
    if config.get("resolution"):
        cmd.extend(["--resolution", str(config["resolution"])])
    
    # Disable parallel if sequential requested
    if config.get("sequential", False):
        cmd.append("--sequential")
    
    # Environment setup for reproducibility
    env = os.environ.copy()
    env["RUST_LOG"] = "info"
    env["RAYON_NUM_THREADS"] = str(config.get("threads", 8))
    env["OMP_NUM_THREADS"] = str(config.get("threads", 8))
    env["SARDINE_ORBIT_CACHE"] = str(output_dir / "orbit_cache")
    Path(env["SARDINE_ORBIT_CACHE"]).mkdir(parents=True, exist_ok=True)
    
    if dry_run:
        print(f"DRY RUN - would execute: {' '.join(cmd)}")
        return {"dry_run": True, "command": cmd}
    
    print(f"🚀 Running benchmark...")
    print(f"   Scene: {scene_path.name}")
    print(f"   Threads: {config.get('threads', 'auto')}")
    print(f"   Geocode: {config.get('geocode', True)}")
    print(f"   Terrain flatten: {config.get('terrain_flatten', True)}")
    print(f"   Command: {' '.join(cmd)}")
    print()
    
    # Force garbage collection before benchmark
    gc.collect()
    
    # Track memory before
    mem_before = get_peak_memory_mb()
    
    # Run with timing
    start_time = time.time()
    start_datetime = datetime.now().isoformat()
    
    log_lines = []
    try:
        proc = subprocess.Popen(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            env=env,
            text=True,
            bufsize=1,
        )
        
        for line in proc.stdout:
            log_lines.append(line.rstrip())
            # Show progress
            if "Step" in line or "📌" in line or "✅" in line or "❌" in line:
                print(f"  {line.rstrip()}")
        
        return_code = proc.wait()
    except Exception as e:
        return {
            "error": str(e),
            "success": False,
        }
    
    end_time = time.time()
    end_datetime = datetime.now().isoformat()
    total_duration = end_time - start_time
    
    # Track memory after
    mem_after = get_peak_memory_mb()
    peak_memory_mb = max(mem_before, mem_after)
    
    # Parse step timings from log
    step_timings = parse_step_timings(log_lines)
    
    # Calculate I/O vs compute breakdown (heuristic)
    io_steps = {1, 2, 3, 13}  # Metadata, Orbit, IW Split, Export
    io_time = sum(
        st["duration_s"] for num, st in step_timings.items() if num in io_steps
    )
    compute_time = sum(
        st["duration_s"] for num, st in step_timings.items() if num not in io_steps
    )
    
    # Build result
    result = {
        "success": return_code == 0,
        "return_code": return_code,
        "config": config,
        "timing": {
            "total_wall_time_s": total_duration,
            "start_time": start_datetime,
            "end_time": end_datetime,
            "io_time_s": io_time,
            "compute_time_s": compute_time,
        },
        "memory": {
            "peak_rss_mb": peak_memory_mb,
        },
        "step_timings": step_timings,
        "machine_info": get_machine_info(),
        "log_line_count": len(log_lines),
    }
    
    # Save log
    log_file = output_dir / "benchmark_log.txt"
    with open(log_file, "w") as f:
        f.write("\n".join(log_lines))
    
    return result


def format_timing_table(step_timings: dict) -> str:
    """Format step timings as a Markdown table."""
    if not step_timings:
        return "No step timings captured."
    
    lines = [
        "| Step | Name | Duration (s) | % of Total |",
        "|------|------|-------------|------------|",
    ]
    
    total = sum(st["duration_s"] for st in step_timings.values())
    
    for step_num in sorted(step_timings.keys()):
        st = step_timings[step_num]
        pct = (st["duration_s"] / total * 100) if total > 0 else 0
        lines.append(f"| {step_num} | {st['name']} | {st['duration_s']:.2f} | {pct:.1f}% |")
    
    lines.append(f"| **Total** | - | **{total:.2f}** | **100%** |")
    return "\n".join(lines)


def generate_baseline_report(result: dict, output_path: Path) -> str:
    """Generate Markdown baseline report."""
    
    machine = result.get("machine_info", {})
    timing = result.get("timing", {})
    config = result.get("config", {})
    
    report = f"""# SARdine Baseline Benchmark Report

**Generated:** {datetime.now().isoformat()}

## Machine Configuration

| Property | Value |
|----------|-------|
| CPU | {machine.get('cpu_model', 'N/A')} |
| CPU Cores | {machine.get('cpu_cores', 'N/A')} |
| RAM | {machine.get('ram_gb', 0):.1f} GB |
| OS | {machine.get('os', 'N/A')} {machine.get('os_release', '')} |
| Python | {machine.get('python_version', 'N/A')} |
| SARdine | {machine.get('sardine_version', 'N/A')} |

## Benchmark Configuration

| Parameter | Value |
|-----------|-------|
| Scene | `{Path(config.get('scene', 'N/A')).name}` |
| Polarization | {config.get('polarization', 'N/A')} |
| Multilook | {config.get('multilook', 'N/A')} |
| Threads | {config.get('threads', 'auto')} |
| Geocode | {config.get('geocode', True)} |
| Terrain Flatten | {config.get('terrain_flatten', True)} |
| Resolution | {config.get('resolution', 10.0)} m |

## Summary Results

| Metric | Value |
|--------|-------|
| **Total Wall Time** | **{timing.get('total_wall_time_s', 0):.2f}s** |
| I/O Time | {timing.get('io_time_s', 0):.2f}s |
| Compute Time | {timing.get('compute_time_s', 0):.2f}s |
| Peak Memory | {result.get('memory', {}).get('peak_rss_mb', 0):.1f} MB |
| Success | {'✅ Yes' if result.get('success') else '❌ No'} |

## Per-Stage Breakdown

{format_timing_table(result.get('step_timings', {}))}

## Hotspot Analysis

Top 5 slowest stages:

"""
    
    # Sort by duration
    sorted_steps = sorted(
        result.get("step_timings", {}).items(),
        key=lambda x: x[1]["duration_s"],
        reverse=True
    )[:5]
    
    for i, (step_num, st) in enumerate(sorted_steps, 1):
        total = timing.get("total_wall_time_s", 1)
        pct = st["duration_s"] / total * 100 if total > 0 else 0
        report += f"{i}. **Step {step_num}: {st['name']}** - {st['duration_s']:.2f}s ({pct:.1f}%)\n"
    
    report += f"""

## Reproducibility Notes

To reproduce this benchmark:

```bash
cd /home/datacube/apps/SARdine/SARdine
export RAYON_NUM_THREADS={config.get('threads', 8)}
export RUST_LOG=info
python -m sardine.cli backscatter \\
    "{config.get('scene')}" \\
    ./benchmark_output \\
    --polarization {config.get('polarization', 'VV')} \\
    --multilook {config.get('multilook', (2,2))[0]} {config.get('multilook', (2,2))[1]} \\
    --num-threads {config.get('threads', 8)} \\
    {'--no-geocode' if not config.get('geocode', True) else ''} \\
    {'--no-terrain-flatten' if not config.get('terrain_flatten', True) else ''}
```

## Raw Data

The complete benchmark data is saved to `benchmark_results.json`.
"""
    
    with open(output_path, "w") as f:
        f.write(report)
    
    return report


def main():
    import argparse
    parser = argparse.ArgumentParser(description="SARdine Baseline Benchmark")
    parser.add_argument("--output-dir", "-o", default="/home/datacube/apps/SARdine/SARdine/benchmark/results",
                        help="Directory for benchmark outputs")
    parser.add_argument("--dry-run", action="store_true", help="Show command without running")
    parser.add_argument("--threads", type=int, default=8, help="Thread count for benchmark")
    parser.add_argument("--geocode", action="store_true", help="Enable geocoding (slower)")
    parser.add_argument("--terrain", action="store_true", help="Enable terrain flattening (slower)")
    parser.add_argument("--full", action="store_true", help="Run full pipeline (geocode + terrain)")
    args = parser.parse_args()
    
    output_dir = Path(args.output_dir)
    
    config = BENCHMARK_CONFIG.copy()
    config["threads"] = args.threads
    if args.geocode or args.full:
        config["geocode"] = True
    if args.terrain or args.full:
        config["terrain_flatten"] = True
    
    print("=" * 70)
    print("🔬 SARdine Baseline Benchmark")
    print("=" * 70)
    print()
    
    # Run benchmark
    result = run_benchmark(output_dir, config, dry_run=args.dry_run)
    
    if args.dry_run:
        return 0
    
    # Save results
    results_file = output_dir / "benchmark_results.json"
    with open(results_file, "w") as f:
        json.dump(result, f, indent=2, default=str)
    print(f"\n📊 Results saved to: {results_file}")
    
    # Generate report
    report_file = output_dir / "bench_baseline.md"
    report = generate_baseline_report(result, report_file)
    print(f"📝 Report saved to: {report_file}")
    
    # Print summary
    print("\n" + "=" * 70)
    print("📈 BENCHMARK SUMMARY")
    print("=" * 70)
    timing = result.get("timing", {})
    print(f"Total wall time: {timing.get('total_wall_time_s', 0):.2f}s")
    print(f"I/O time: {timing.get('io_time_s', 0):.2f}s")
    print(f"Compute time: {timing.get('compute_time_s', 0):.2f}s")
    print(f"Peak memory: {result.get('memory', {}).get('peak_rss_mb', 0):.1f} MB")
    
    if result.get("step_timings"):
        print("\nTop 3 hotspots:")
        sorted_steps = sorted(
            result["step_timings"].items(),
            key=lambda x: x[1]["duration_s"],
            reverse=True
        )[:3]
        for step_num, st in sorted_steps:
            print(f"  Step {step_num}: {st['name']} - {st['duration_s']:.2f}s")
    
    return 0 if result.get("success") else 1


if __name__ == "__main__":
    sys.exit(main())
