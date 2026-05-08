"""Integration test for the full SARdine backscatter pipeline.

This test drives the 14-step pipeline against real Sentinel-1 metadata to
validate that the processor completes without resorting to synthetic fallbacks.

To run, provide the environment variable ``SARDINE_TEST_SAFE_PATH`` pointing to
an unpacked ``.SAFE`` directory (or zip) and optionally ``SARDINE_TEST_ORBIT_CACHE``
with pre-downloaded precise orbit files to avoid network traffic during CI.
"""

from __future__ import annotations

import os
from pathlib import Path

import pytest

pytest.importorskip("sardine")

from sardine.processors import BackscatterProcessor


@pytest.mark.slow
@pytest.mark.integration
def test_backscatter_pipeline_real_metadata(tmp_path: Path) -> None:
    safe_input = os.environ.get("SARDINE_TEST_SAFE_PATH")
    if safe_input is None:
        pytest.skip(
            "SARDINE_TEST_SAFE_PATH not set; provide real SAFE data to run integration test"
        )

    orbit_cache_override = os.environ.get("SARDINE_TEST_ORBIT_CACHE")
    if orbit_cache_override:
        os.environ["SARDINE_ORBIT_CACHE"] = orbit_cache_override

    output_dir = tmp_path / "backscatter"
    options = {
        "verbose": False,
        "terrain_flatten": True,
        "geocode": True,
        "resolution": 30.0,
        "optimization_mode": "complete",
        "allow_synthetic": False,
        "use_real_orbit": True,
    }

    processor = BackscatterProcessor(safe_input, output_dir, options)
    processor.process_backscatter()

    statuses = {entry["step"]: entry["status"] for entry in processor.processing_log}
    assert statuses.get(2) == "success", "Precise orbit application must succeed"
    assert statuses.get(8) == "success", "Terrain flattening must succeed with real metadata"
    assert statuses.get(10) == "success", "Terrain correction must succeed with valid orbit data"

    assert processor.geo_transform is not None, "Geocoded output should provide geo-transform metadata"

    # Regression guard: IW products must perform true multilooking (range looks > 1)
    assert (
        processor.actual_range_looks is not None and processor.actual_range_looks > 1.0
    ), "Step 8 multilooking regressed to single-look behaviour"
