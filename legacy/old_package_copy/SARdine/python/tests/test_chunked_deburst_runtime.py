"""Runtime tests for chunked deburst processing with actual Sentinel-1 data.

NOTE (Jan 2026): Chunked deburst + separable LUT APIs were removed because IW TOPS
mode requires a dense 2D LUT. This file is kept for historical reference and will
skip unless the deprecated APIs are present.
"""

from __future__ import annotations

import os
import time
import logging
from pathlib import Path

import numpy as np
import pytest

pytest.importorskip("sardine")

import sardine


if not hasattr(sardine, "deburst_topsar_chunked") or not hasattr(sardine, "prepare_calibration_luts"):
    pytest.skip(
        "Chunked deburst APIs (deburst_topsar_chunked / prepare_calibration_luts) were removed (Jan 2026).",
        allow_module_level=True,
    )

logger = logging.getLogger(__name__)


def get_test_product_path():
    """Get test product path from environment variable."""
    safe_path = os.environ.get("SARDINE_TEST_SAFE_PATH")
    if safe_path is None:
        pytest.skip(
            "SARDINE_TEST_SAFE_PATH not set. "
            "Set it to a Sentinel-1 .SAFE directory or zip file to run runtime tests."
        )
    return safe_path


def get_orbit_cache():
    """Get orbit cache path from environment variable."""
    orbit_cache = os.environ.get("SARDINE_TEST_ORBIT_CACHE")
    if orbit_cache:
        os.environ["SARDINE_ORBIT_CACHE"] = orbit_cache
    return orbit_cache


@pytest.mark.slow
@pytest.mark.integration
def test_chunked_vs_standard_numerical_equivalence():
    """
    Test that chunked deburst produces identical results to standard deburst.
    
    NOTE: This test requires:
    - Valid Sentinel-1 product (set SARDINE_TEST_SAFE_PATH)
    - Orbit cache with precise orbit files (set SARDINE_ORBIT_CACHE)
    - May need SARDINE_CAL_ALLOW_SEPARABLE_ERR=1 if separable LUT precomputation fails
    """
    """
    Test that chunked deburst produces identical results to standard deburst.
    
    This is the core correctness test - chunked and standard deburst should
    produce identical calibrated output (within floating point precision).
    """
    product_path = get_test_product_path()
    get_orbit_cache()
    
    logger.info(f"Testing numerical equivalence with product: {product_path}")
    
    # Create reader
    reader = sardine.create_cached_slc_reader(str(product_path))
    
    # Test with first available subswath and polarization
    # In practice, you might want to test multiple subswaths
    subswath = "IW1"
    polarization = "VV"
    
    # Get calibration data (needed for chunked deburst)
    # Note: prepare_calibration_luts may fail if separable LUT precomputation fails
    # In that case, we'll skip the test with a clear message
    logger.info("Preparing calibration LUTs...")
    cal_result = sardine.prepare_calibration_luts(
        reader, subswath, polarization, "sigma0", apply_noise_removal=False
    )
    
    if not isinstance(cal_result, dict) or cal_result.get("status") != "success":
        error_msg = cal_result.get("message", "Unknown error") if isinstance(cal_result, dict) else str(cal_result)
        pytest.skip(
            f"Could not prepare calibration LUTs: {error_msg}. "
            f"This may indicate an issue with separable LUT precomputation. "
            f"Try setting SARDINE_CAL_ALLOW_SEPARABLE_ERR=1 or check logs for details."
        )
    
    cal_azimuth = np.asarray(cal_result["calibration_azimuth"], dtype=np.float32)
    cal_range = np.asarray(cal_result["calibration_range"], dtype=np.float32)
    noise_lut = cal_result.get("noise_lut")
    if noise_lut is not None:
        noise_lut = np.asarray(noise_lut, dtype=np.float32)
    
    # Run standard deburst + calibration
    logger.info("Running standard deburst...")
    standard_start = time.time()
    standard_deburst_result = sardine.deburst_topsar_cached(reader, subswath, polarization)
    
    if not isinstance(standard_deburst_result, dict) or standard_deburst_result.get("status") != "success":
        error_msg = standard_deburst_result.get("message", "Unknown error") if isinstance(standard_deburst_result, dict) else str(standard_deburst_result)
        if "orbit" in error_msg.lower() or "Orbit" in error_msg:
            pytest.skip(
                f"Orbit data required but not available: {error_msg}. "
                f"Set SARDINE_ORBIT_CACHE to a directory containing precise orbit files, "
                f"or use DownloadManager to download orbit files first."
            )
        pytest.fail(f"Standard deburst failed: {error_msg}")
    
    standard_deburst_data = np.asarray(standard_deburst_result["data"], dtype=np.complex64)
    standard_deburst_time = time.time() - standard_start
    
    # Calibrate standard deburst result
    logger.info("Calibrating standard deburst result...")
    cal_start = time.time()
    # Check function signature - may need apply_noise_removal instead
    standard_cal_result = sardine.radiometric_calibration_with_denoising_cached(
        reader,
        subswath,
        polarization,
        "sigma0",  # Use sigma0 for comparison
        standard_deburst_data,
        apply_noise_removal=(noise_lut is not None),
    )
    calibration_time = time.time() - cal_start
    
    if not isinstance(standard_cal_result, dict) or standard_cal_result.get("status") != "success":
        pytest.fail(f"Standard calibration failed: {standard_cal_result}")
    
    standard_calibrated = np.asarray(standard_cal_result["calibrated_data"], dtype=np.float32)
    standard_total_time = standard_deburst_time + calibration_time
    
    # Run chunked deburst (includes calibration)
    logger.info("Running chunked deburst...")
    chunked_start = time.time()
    chunked_result = sardine.deburst_topsar_chunked(
        reader,
        subswath,
        polarization,
        cal_azimuth.tolist(),
        cal_range.tolist(),
        noise_lut=noise_lut,
        chunk_lines=512,
    )
    chunked_time = time.time() - chunked_start
    
    if not isinstance(chunked_result, dict) or chunked_result.get("status") != "success":
        pytest.fail(f"Chunked deburst failed: {chunked_result}")
    
    chunked_calibrated = np.asarray(chunked_result["calibrated_data"], dtype=np.float32)
    
    # Compare results
    logger.info("Comparing results...")
    logger.info(f"Standard shape: {standard_calibrated.shape}")
    logger.info(f"Chunked shape: {chunked_calibrated.shape}")
    
    # Shapes must match
    assert standard_calibrated.shape == chunked_calibrated.shape, (
        f"Shape mismatch: standard {standard_calibrated.shape} vs "
        f"chunked {chunked_calibrated.shape}"
    )
    
    # Compare values (accounting for floating point precision)
    # Use relatively tight tolerances since both should use same algorithms
    finite_mask = np.isfinite(standard_calibrated) & np.isfinite(chunked_calibrated)
    finite_count = np.count_nonzero(finite_mask)
    
    if finite_count == 0:
        pytest.fail("No finite values to compare")
    
    standard_finite = standard_calibrated[finite_mask]
    chunked_finite = chunked_calibrated[finite_mask]
    
    # Calculate statistics
    max_diff = np.max(np.abs(standard_finite - chunked_finite))
    mean_diff = np.mean(np.abs(standard_finite - chunked_finite))
    rel_diff = np.abs(standard_finite - chunked_finite) / (np.abs(standard_finite) + 1e-10)
    max_rel_diff = np.max(rel_diff)
    mean_rel_diff = np.mean(rel_diff)
    
    logger.info(f"Comparison statistics:")
    logger.info(f"  Finite pixels: {finite_count}/{standard_calibrated.size}")
    logger.info(f"  Max absolute difference: {max_diff:.6e}")
    logger.info(f"  Mean absolute difference: {mean_diff:.6e}")
    logger.info(f"  Max relative difference: {max_rel_diff:.6e}")
    logger.info(f"  Mean relative difference: {mean_rel_diff:.6e}")
    logger.info(f"  Standard processing time: {standard_total_time:.2f}s")
    logger.info(f"  Chunked processing time: {chunked_time:.2f}s")
    
    # Assert numerical equivalence
    # Use reasonable tolerances for floating point comparison
    rtol = 1e-5  # Relative tolerance
    atol = 1e-6  # Absolute tolerance
    
    assert np.allclose(
        standard_calibrated,
        chunked_calibrated,
        rtol=rtol,
        atol=atol,
        equal_nan=True,
    ), (
        f"Results differ beyond tolerance (rtol={rtol}, atol={atol}). "
        f"Max diff: {max_diff:.6e}, Max rel diff: {max_rel_diff:.6e}"
    )
    
    logger.info("✅ Numerical equivalence test passed!")


@pytest.mark.slow
@pytest.mark.integration
def test_chunked_deburst_performance():
    """
    Test performance of chunked deburst vs standard deburst.
    
    Measures:
    - Processing time
    - Memory usage (if available)
    """
    product_path = get_test_product_path()
    get_orbit_cache()
    
    logger.info(f"Testing performance with product: {product_path}")
    
    reader = sardine.create_cached_slc_reader(str(product_path))
    subswath = "IW1"
    polarization = "VV"
    
    # Get calibration data
    cal_result = sardine.prepare_calibration_luts(
        reader, subswath, polarization, "sigma0", apply_noise_removal=False
    )
    if not isinstance(cal_result, dict) or cal_result.get("status") != "success":
        pytest.skip(f"Could not prepare calibration LUTs: {cal_result}")
    
    cal_azimuth = np.asarray(cal_result["calibration_azimuth"], dtype=np.float32)
    cal_range = np.asarray(cal_result["calibration_range"], dtype=np.float32)
    noise_lut = cal_result.get("noise_lut")
    if noise_lut is not None:
        noise_lut = np.asarray(noise_lut, dtype=np.float32)
    
    # Standard deburst timing
    logger.info("Timing standard deburst...")
    standard_times = []
    for i in range(2):  # Run twice, use second for timing (warm cache)
        start = time.time()
        deburst_result = sardine.deburst_topsar_cached(reader, subswath, polarization)
        deburst_data = np.asarray(deburst_result["data"], dtype=np.complex64)
        cal_result = sardine.radiometric_calibration_with_denoising_cached(
            reader, subswath, polarization, "sigma0", deburst_data,
            apply_noise_removal=(noise_lut is not None),
        )
        standard_times.append(time.time() - start)
    
    standard_time = standard_times[-1]  # Use second run (warm cache)
    
    # Chunked deburst timing
    logger.info("Timing chunked deburst...")
    chunked_times = []
    for i in range(2):  # Run twice, use second for timing
        start = time.time()
        chunked_result = sardine.deburst_topsar_chunked(
            reader, subswath, polarization,
            cal_azimuth.tolist(), cal_range.tolist(),
            noise_lut=noise_lut, chunk_lines=512,
        )
        chunked_times.append(time.time() - start)
    
    chunked_time = chunked_times[-1]  # Use second run
    
    # Report results
    speedup = standard_time / chunked_time if chunked_time > 0 else 0
    logger.info(f"Performance results:")
    logger.info(f"  Standard time: {standard_time:.2f}s")
    logger.info(f"  Chunked time: {chunked_time:.2f}s")
    logger.info(f"  Speedup: {speedup:.2f}x")
    
    # Performance assertion (chunked should be similar or better)
    # We don't fail if it's slower, just log it
    if speedup < 1.0:
        logger.warning(f"Chunked deburst is {1/speedup:.2f}x slower than standard")
    else:
        logger.info(f"Chunked deburst is {speedup:.2f}x faster than standard")
    
    # Test passes as long as both complete successfully
    assert standard_time > 0 and chunked_time > 0, "Both methods must complete"


@pytest.mark.slow
@pytest.mark.integration
def test_chunked_deburst_chunk_size_variation():
    """
    Test that chunked deburst produces consistent results across different chunk sizes.
    
    Different chunk sizes should produce identical results (within floating point precision).
    """
    product_path = get_test_product_path()
    get_orbit_cache()
    
    logger.info(f"Testing chunk size variation with product: {product_path}")
    
    reader = sardine.create_cached_slc_reader(str(product_path))
    subswath = "IW1"
    polarization = "VV"
    
    # Get calibration data
    cal_result = sardine.prepare_calibration_luts(
        reader, subswath, polarization, "sigma0", apply_noise_removal=False
    )
    if not isinstance(cal_result, dict) or cal_result.get("status") != "success":
        pytest.skip(f"Could not prepare calibration LUTs: {cal_result}")
    
    cal_azimuth = np.asarray(cal_result["calibration_azimuth"], dtype=np.float32)
    cal_range = np.asarray(cal_result["calibration_range"], dtype=np.float32)
    noise_lut = cal_result.get("noise_lut")
    if noise_lut is not None:
        noise_lut = np.asarray(noise_lut, dtype=np.float32)
    
    # Test different chunk sizes
    chunk_sizes = [256, 512, 1024, 2048]
    results = {}
    
    for chunk_size in chunk_sizes:
        logger.info(f"Testing chunk size: {chunk_size} lines...")
        result = sardine.deburst_topsar_chunked(
            reader, subswath, polarization,
            cal_azimuth.tolist(), cal_range.tolist(),
            noise_lut=noise_lut, chunk_lines=chunk_size,
        )
        
        if not isinstance(result, dict) or result.get("status") != "success":
            pytest.fail(f"Chunked deburst failed for chunk_size={chunk_size}: {result}")
        
        results[chunk_size] = np.asarray(result["calibrated_data"], dtype=np.float32)
    
    # Compare all results to first one (baseline)
    baseline_size = chunk_sizes[0]
    baseline = results[baseline_size]
    
    logger.info(f"Comparing all chunk sizes to baseline ({baseline_size} lines)...")
    
    for chunk_size in chunk_sizes[1:]:
        data = results[chunk_size]
        
        assert baseline.shape == data.shape, (
            f"Shape mismatch for chunk_size={chunk_size}: "
            f"{baseline.shape} vs {data.shape}"
        )
        
        # Compare values
        finite_mask = np.isfinite(baseline) & np.isfinite(data)
        if np.count_nonzero(finite_mask) == 0:
            pytest.fail(f"No finite values to compare for chunk_size={chunk_size}")
        
        max_diff = np.max(np.abs(baseline[finite_mask] - data[finite_mask]))
        logger.info(f"  Chunk size {chunk_size}: max diff = {max_diff:.6e}")
        
        # Results should be identical (same algorithm, different chunking)
        assert np.allclose(
            baseline, data, rtol=1e-6, atol=1e-7, equal_nan=True
        ), f"Results differ for chunk_size={chunk_size}, max diff={max_diff:.6e}"
    
    logger.info("✅ All chunk sizes produce identical results!")


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    pytest.main([__file__, "-v", "-s"])

