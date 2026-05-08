import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

import pytest


def _repo_root() -> Path:
    # This file lives at: SARdine/SARdine/tests/...
    return Path(__file__).resolve().parents[2]


@pytest.mark.integration
@pytest.mark.slow
def test_snap_sigma0_deburst_against_sardine(tmp_path: Path) -> None:
    safe_dir = os.environ.get("SARDINE_TEST_SAFE_DIR")
    if not safe_dir:
        pytest.skip("Set SARDINE_TEST_SAFE_DIR to run this integration test")

    snap_gpt = os.environ.get("SARDINE_TEST_SNAP_GPT", "gpt")
    if snap_gpt == "gpt":
        if shutil.which("gpt") is None:
            pytest.skip("SNAP GPT not found in PATH; set SARDINE_TEST_SNAP_GPT")
    else:
        if not Path(snap_gpt).exists():
            pytest.skip("SARDINE_TEST_SNAP_GPT points to missing executable")

    orbit_cache = os.environ.get("SARDINE_TEST_ORBIT_CACHE")

    script = _repo_root() / "validation" / "run_snap_validation.py"
    if not script.exists():
        pytest.fail(f"Missing validation script: {script}")

    out_dir = tmp_path / "snap_validation"

    cmd = [
        sys.executable,
        str(script),
        str(Path(safe_dir).resolve()),
        "--mode",
        "sigma0_deburst",
        "--subswath",
        os.environ.get("SARDINE_TEST_SUBSWATH", "IW1"),
        "--polarization",
        os.environ.get("SARDINE_TEST_POL", "VH"),
        "--snap-gpt",
        snap_gpt,
        "--output-dir",
        str(out_dir),
    ]
    if orbit_cache:
        cmd += ["--orbit-cache", str(Path(orbit_cache).resolve())]

    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        raise AssertionError(
            "SNAP validation script failed\n"
            f"cmd: {' '.join(cmd)}\n"
            f"stdout:\n{proc.stdout}\n"
            f"stderr:\n{proc.stderr}\n"
        )

    report_path = out_dir / "validation_report.json"
    assert report_path.exists(), "Expected validation_report.json"

    report = json.loads(report_path.read_text())
    cmp = report["comparison"]

    # Basic sanity checks
    assert cmp["matching_pixels"] >= 1000

    # Prefer dB-space metrics (more meaningful for backscatter)
    # These thresholds mirror the defaults in sardine.validation.reference_comparator
    assert cmp["mae_db"] is not None
    assert cmp["rmse_db"] is not None
    assert cmp["correlation_db"] is not None

    assert cmp["mae_db"] <= float(os.environ.get("SARDINE_TEST_MAX_MAE_DB", "1.0"))
    assert cmp["rmse_db"] <= float(os.environ.get("SARDINE_TEST_MAX_RMSE_DB", "1.5"))
    assert cmp["correlation_db"] >= float(os.environ.get("SARDINE_TEST_MIN_CORR", "0.95"))
