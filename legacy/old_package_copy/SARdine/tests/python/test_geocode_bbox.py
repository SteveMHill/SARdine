"""Unit tests for geocoding bounding box refinement."""

from sardine.processors.backscatter.geometry import refine_geocoding_bbox


def test_refine_geocoding_bbox_returns_valid_result():
    """Test that refine_geocoding_bbox returns valid bbox and metrics."""
    # Original bbox: 6 deg lon x 3 deg lat
    original_bbox = [10.0, 50.0, 16.0, 53.0]
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

    # Verify bbox format [min_lon, min_lat, max_lon, max_lat]
    assert len(refined_bbox) == 4
    assert refined_bbox[0] < refined_bbox[2]  # min_lon < max_lon
    assert refined_bbox[1] < refined_bbox[3]  # min_lat < max_lat

    # Verify metrics contain expected keys
    assert "original_lat_extent_deg" in metrics
    assert "original_lon_extent_deg" in metrics
    assert "expected_lat_extent_deg" in metrics
    assert "expected_lon_extent_deg" in metrics
    assert "new_lat_extent_deg" in metrics
    assert "new_lon_extent_deg" in metrics
    assert "shrink_lat_pct" in metrics
    assert "shrink_lon_pct" in metrics

    # Shrink percentages should be non-negative
    assert metrics["shrink_lat_pct"] >= 0.0
    assert metrics["shrink_lon_pct"] >= 0.0


def test_refine_geocoding_bbox_no_shrink_when_within_threshold():
    """Test that bbox is not shrunk when expected extent is close to original."""
    # Create a bbox that matches the expected image dimensions
    # At lat 51.5, 1 deg lat ≈ 111320m, 1 deg lon ≈ 69000m
    rows = 2000
    cols = 3000
    range_spacing_m = 20.0   # 3000 * 20 = 60000m ≈ 0.87 deg lon
    azimuth_spacing_m = 14.0  # 2000 * 14 = 28000m ≈ 0.25 deg lat

    # Set bbox to match expected dimensions (within threshold)
    # Expected with 5% margin: 0.26 deg lat, 0.91 deg lon
    # Threshold for shrinking: expected * 1.2
    original_bbox = [10.0, 51.0, 11.0, 51.3]  # 1 deg lon x 0.3 deg lat

    result = refine_geocoding_bbox(
        original_bbox,
        rows,
        cols,
        range_spacing_m,
        azimuth_spacing_m,
    )

    assert result is not None
    refined_bbox, metrics = result

    # Should not shrink much since bbox is close to expected
    assert metrics["shrink_lat_pct"] < 0.5  # Less than 50% shrink
    assert metrics["shrink_lon_pct"] < 0.5


def test_refine_geocoding_bbox_invalid_inputs_return_none():
    """Test that invalid inputs return None."""
    original_bbox = [10.0, 50.0, 16.0, 53.0]

    # Zero rows
    assert (
        refine_geocoding_bbox(original_bbox, 0, 100, 90.0, 110.0)
        is None
    )
    # Zero cols
    assert (
        refine_geocoding_bbox(original_bbox, 100, 0, 90.0, 110.0)
        is None
    )
    # Negative range spacing
    assert (
        refine_geocoding_bbox(original_bbox, 100, 100, -1.0, 110.0)
        is None
    )
    # NaN azimuth spacing
    assert (
        refine_geocoding_bbox(original_bbox, 100, 100, 90.0, float("nan"))
        is None
    )
