"""Tests for error propagation at PyO3 binding boundaries.

These tests ensure that Rust errors are surfaced to Python as exceptions
rather than causing process aborts (e.g., via `unwrap()` panics).
"""

import pytest

import sardine


@pytest.mark.integration
def test_create_masking_workflow_invalid_gamma_range_raises_valueerror() -> None:
    """create_masking_workflow should raise ValueError on invalid gamma range.

    This exercises a Rust -> Py error path at the binding boundary and
    verifies that the error is mapped to a Python ValueError instead of
    causing a panic.
    """

    # gamma0_min >= gamma0_max is invalid and should trigger a ValueError
    with pytest.raises(ValueError) as excinfo:
        sardine._core.create_masking_workflow(
            lia_threshold=45.0,
            dem_threshold=0.0,
            gamma0_min=2.0,
            gamma0_max=1.0,
        )

    msg = str(excinfo.value)
    assert "gamma0_min must be less than gamma0_max" in msg
