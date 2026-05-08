#!/usr/bin/env python3
import os
import sys
import time
import subprocess
from pathlib import Path

SAFE_DEFAULT = "/home/datacube/apps/SARdine/data/S1A_IW_SLC__1SDV_20201230T165244_20201230T165311_035918_0434F0_6788.SAFE"


def main():
    safe = sys.argv[1] if len(sys.argv) > 1 else SAFE_DEFAULT
    outdir = Path(sys.argv[2]) if len(sys.argv) > 2 else Path("/home/datacube/apps/SARdine/SARdine/rtc_geotiff_output/parallel_timed")
    outdir.mkdir(parents=True, exist_ok=True)

    env = os.environ.copy()
    env.setdefault("RUST_LOG", "info")
    env["SARDINE_ORBIT_CACHE"] = str(outdir / "orbit_cache")
    Path(env["SARDINE_ORBIT_CACHE"]).mkdir(parents=True, exist_ok=True)

    log_path = outdir / "cli_run.log"
    print(f"Output: {outdir}")
    print(f"RUST_LOG={env['RUST_LOG']}")

    cmd = [
        sys.executable, "-m", "sardine.cli", "backscatter",
        safe, str(outdir),
        "--polarization", "VV",
        "--parallel",
        "--no-geocode",
        "--no-terrain-flatten",
    ]

    # Timing markers from Rust logs
    t = {
        "noise_lut_start": None,
        "noise_lut_end": None,
        "denoise_start": None,
        "denoise_end": None,
    }

    start = time.time()
    with open(log_path, "w") as lf:
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, env=env, text=True, bufsize=1)
        for line in proc.stdout:
            ts = time.time()
            line_stripped = line.rstrip()
            lf.write(f"{ts:.3f} {line_stripped}\n")
            lf.flush()
            # Detect markers
            if "Pre-computing noise LUT for" in line_stripped and t["noise_lut_start"] is None:
                t["noise_lut_start"] = ts
            elif "Noise LUT pre-computation completed" in line_stripped and t["noise_lut_end"] is None:
                t["noise_lut_end"] = ts
            elif "Applying thermal noise removal to" in line_stripped and t["denoise_start"] is None:
                t["denoise_start"] = ts
            elif "Thermal noise removal complete" in line_stripped and t["denoise_end"] is None:
                t["denoise_end"] = ts
        ret = proc.wait()

    end = time.time()

    # Summarize timings
    def dur(a, b):
        return None if a is None or b is None else b - a

    summary = {
        "total_seconds": end - start,
        "noise_lut_seconds": dur(t["noise_lut_start"], t["noise_lut_end"]),
        "denoise_seconds": dur(t["denoise_start"], t["denoise_end"]),
        "markers": t,
        "output_dir": str(outdir),
        "return_code": ret,
    }

    import json
    with open(outdir / "thermal_noise_timing.json", "w") as f:
        json.dump(summary, f, indent=2)

    print("Timing summary:")
    print(json.dumps(summary, indent=2))
    return ret


if __name__ == "__main__":
    raise SystemExit(main())
