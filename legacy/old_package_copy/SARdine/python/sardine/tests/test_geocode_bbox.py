"""Unit tests for geocoding bounding box refinement.

These tests live inside the `sardine` package, so we make `sardine.tests` a
package to avoid pytest import-name collisions with other `test_geocode_bbox.py`
modules elsewhere in the repo.
"""

from sardine.processors.backscatter.geometry import refine_geocoding_bbox


def test_refine_geocoding_bbox_shrinks_large_bbox_with_cap():
    # A metadata bbox that is far larger than the implied footprint from
    # (rows, cols, spacings) should be shrunk, but never by more than 30%.
    original_bbox = [10.0, 50.0, 16.0, 53.0]  # 6 deg lon x 3 deg lat
    rows = 500
    cols = 700
    range_spacing_m = 80.0
    azimuth_spacing_m = 100.0

    result = refine_geocoding_bbox(
        original_bbox,
        rows,
        cols,
        range_spacing_m,
        azimuth_spacing_m,
    )

    assert result is not None
    refined_bbox, metrics = result

    # Bbox stays well-formed
    assert len(refined_bbox) == 4
    assert refined_bbox[0] < refined_bbox[2]
    assert refined_bbox[1] < refined_bbox[3]

    # For this case, we expect shrinking to occur (but be capped).
    assert 0.0 < metrics["shrink_lat_pct"] <= 0.30
    assert 0.0 < metrics["shrink_lon_pct"] <= 0.30

    assert metrics["new_lat_extent_deg"] <= metrics["original_lat_extent_deg"]
    assert metrics["new_lon_extent_deg"] <= metrics["original_lon_extent_deg"]


def test_refine_geocoding_bbox_invalid_inputs_return_none():
    original_bbox = [10.0, 50.0, 16.0, 53.0]

    assert refine_geocoding_bbox(original_bbox, 0, 100, 90.0, 110.0) is None
    assert refine_geocoding_bbox(original_bbox, 100, 0, 90.0, 110.0) is None
    assert refine_geocoding_bbox(original_bbox, 100, 100, -1.0, 110.0) is None
    assert refine_geocoding_bbox(original_bbox, 100, 100, 90.0, float("nan")) is None
