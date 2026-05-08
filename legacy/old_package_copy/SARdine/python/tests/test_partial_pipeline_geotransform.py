"""Unit tests for geotransform validation utilities in the partial pipeline."""

from __future__ import annotations

import pytest

from sardine.partial_pipeline import (
    _compute_slice_geotransform,
    _validate_geotransform_uniqueness,
)


def test_validate_geotransform_uniqueness_allows_distinct_transforms() -> None:
    transforms = {
        "IW1": (0.0, 10.0, 0.0, 100.0, 0.0, -10.0),
        "IW2": (10.0, 10.0, 0.0, 90.0, 0.0, -10.0),
        "IW3": (20.0, 10.0, 0.0, 80.0, 0.0, -10.0),
    }

    # Should not raise when all transforms are distinct within tolerance
    _validate_geotransform_uniqueness(transforms)


def test_validate_geotransform_uniqueness_detects_duplicates() -> None:
    transforms = {
        "IW1": (0.0, 10.0, 0.0, 100.0, 0.0, -10.0),
        "IW2": (0.0, 10.0, 0.0, 100.0, 0.0, -10.0),
    }

    with pytest.raises(RuntimeError) as excinfo:
        _validate_geotransform_uniqueness(transforms)

    assert "Duplicate geotransform" in str(excinfo.value)
    assert "IW1" in str(excinfo.value)
    assert "IW2" in str(excinfo.value)


def test_compute_slice_geotransform_reflects_slice_offsets() -> None:
    base_transform = (0.0, 10.0, 0.0, 100.0, 0.0, -10.0)
    first_slice = slice(0, 256)
    second_slice = slice(256, 512)
    first_sample_slice = slice(0, 512)
    second_sample_slice = slice(512, 1024)

    first_geo = _compute_slice_geotransform(base_transform, first_slice, first_sample_slice)
    second_geo = _compute_slice_geotransform(base_transform, second_slice, second_sample_slice)

    assert first_geo != second_geo
    # Sanity-check that the validation logic still accepts the distinct transforms
    _validate_geotransform_uniqueness({"IW1": first_geo, "IW2": second_geo})
