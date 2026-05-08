"""
Quality-check helpers for the backscatter pipeline.

These utilities are intentionally free of side effects; the caller is
responsible for logging or emitting validation messages.
"""

from typing import List, Tuple
import numpy as np


def validate_backscatter_range(data: np.ndarray) -> Tuple[bool, List[str]]:
    """
    Validate that backscatter values are within expected physical ranges.

    Returns a tuple of (is_valid, warnings).
    """
    if not isinstance(data, np.ndarray) or data.size == 0:
        return False, ["Invalid input: data is not a numpy array or is empty"]

    finite_data = data[np.isfinite(data)]
    if finite_data.size == 0:
        return False, ["All values are non-finite"]

    data_min = float(np.min(finite_data))
    data_max = float(np.max(finite_data))

    warnings: List[str] = []

    # Heuristic: if min is negative and max is modest, treat as dB
    if data_min < 0 and data_max < 50:
        if data_min < -60:
            warnings.append(f"Very low values ({data_min:.1f} dB) may indicate noise floor issues")
        if data_max > 20:
            warnings.append(f"Very high values ({data_max:.1f} dB) may indicate calibration issues")
    else:
        # Linear power validation
        if data_min > 0 and data_max > 100:
            warnings.append(f"Linear power values ({data_max:.2f}) seem unusually high")
        if data_min < 0:
            warnings.append(f"Negative linear power values detected (min={data_min:.2e})")

    return len(warnings) == 0, warnings


def validate_power_preservation(
    pre_data: np.ndarray, post_data: np.ndarray, tolerance: float = 0.05
) -> Tuple[bool, float, float, float]:
    """
    Validate that mean power is preserved within tolerance.

    Uses >= 0 filter (not > 0) to match Rust multilooking behavior. This is
    important because multilooking correctly averages all pixels in a block
    including zeros at edges. Using > 0 filter creates false positives when
    edge blocks mix zeros with valid data.

    Returns (is_valid, relative_change, pre_mean, post_mean). If insufficient
    data or near-zero power, relative_change will be 0.0 and validation passes.
    """
    if not isinstance(pre_data, np.ndarray) or not isinstance(post_data, np.ndarray):
        return False, 0.0, 0.0, 0.0  # Invalid input types

    # Use >= 0 to match Rust multilooking behavior (includes zeros in averaging)
    # This avoids false power preservation failures when edge blocks contain
    # scattered zeros mixed with valid data
    pre_finite = pre_data[np.isfinite(pre_data) & (pre_data >= 0)]
    post_finite = post_data[np.isfinite(post_data) & (post_data >= 0)]

    if pre_finite.size < 100 or post_finite.size < 100:
        # Insufficient data to validate - pass with warning
        # Note: This is intentionally True because we can't validate, not because validation failed
        return True, 0.0, 0.0, 0.0

    pre_mean = float(np.mean(pre_finite))
    post_mean = float(np.mean(post_finite))

    if pre_mean < 1e-12:
        # Near-zero power - can't compute relative change meaningfully
        return True, 0.0, pre_mean, post_mean

    relative_change = abs(post_mean - pre_mean) / pre_mean
    is_valid = relative_change <= tolerance
    return is_valid, relative_change, pre_mean, post_mean
