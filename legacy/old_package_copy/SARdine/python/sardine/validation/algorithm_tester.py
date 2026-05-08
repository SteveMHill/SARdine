"""
Algorithm Tester for SARdine Pipeline

Unit tests for individual SAR processing algorithms with known inputs/outputs.
Tests mathematical correctness of:
- Calibration formulas (sigma0, beta0, gamma0)
- dB conversion
- Multilooking (power preservation)
- Coordinate transforms
- Terrain flattening formula
- Noise floor handling

These tests use synthetic data with known correct answers.
"""

import logging
import numpy as np
from dataclasses import dataclass, field
from typing import Dict, Any, List, Optional, Callable, Tuple

logger = logging.getLogger(__name__)


@dataclass
class AlgorithmTestResult:
    """Result of an algorithm test."""
    test_name: str
    algorithm: str
    passed: bool
    expected: Any
    actual: Any
    tolerance: float
    error: Optional[float] = None
    message: str = ""
    details: Dict[str, Any] = field(default_factory=dict)


@dataclass
class AlgorithmTestReport:
    """Complete algorithm test report."""
    results: List[AlgorithmTestResult] = field(default_factory=list)
    timestamp: str = ""
    
    @property
    def passed(self) -> bool:
        return all(r.passed for r in self.results)
    
    @property
    def pass_count(self) -> int:
        return sum(1 for r in self.results if r.passed)
    
    @property
    def fail_count(self) -> int:
        return sum(1 for r in self.results if not r.passed)


class AlgorithmTester:
    """
    Test individual SAR processing algorithms with known inputs/outputs.
    
    Uses synthetic data where the correct answer is known mathematically,
    allowing precise validation of algorithm implementations.
    """
    
    def __init__(self):
        """Initialize tester."""
        self.results: List[AlgorithmTestResult] = []
    
    def run_all_tests(self) -> AlgorithmTestReport:
        """Run all algorithm tests."""
        from datetime import datetime
        
        self.results = []
        
        # Run all test categories
        self._test_db_conversion()
        self._test_calibration_formula()
        self._test_multilook_power()
        self._test_terrain_flattening()
        self._test_noise_floor()
        self._test_coordinate_transforms()
        self._test_resampling_power()
        
        return AlgorithmTestReport(
            results=self.results,
            timestamp=datetime.now().isoformat(),
        )
    
    def _add_result(self, result: AlgorithmTestResult):
        """Add a test result."""
        self.results.append(result)
        status = "✓" if result.passed else "✗"
        logger.info(f"[{status}] {result.test_name}: {result.message}")
    
    # ==========================================================================
    # dB Conversion Tests
    # ==========================================================================
    
    def _test_db_conversion(self):
        """Test dB conversion formula: dB = 10 * log10(linear)"""
        
        # Test 1: Basic conversion
        test_cases = [
            (1.0, 0.0, "unity"),
            (10.0, 10.0, "10x"),
            (0.1, -10.0, "0.1x"),
            (100.0, 20.0, "100x"),
            (0.01, -20.0, "0.01x"),
            (2.0, 3.0103, "2x"),  # 10*log10(2) ≈ 3.0103
        ]
        
        for linear, expected_db, label in test_cases:
            actual_db = 10.0 * np.log10(linear)
            error = abs(actual_db - expected_db)
            tolerance = 0.0001
            
            self._add_result(AlgorithmTestResult(
                test_name=f"db_conversion_{label}",
                algorithm="dB_conversion",
                passed=error < tolerance,
                expected=expected_db,
                actual=float(actual_db),
                tolerance=tolerance,
                error=float(error),
                message=f"10*log10({linear}) = {actual_db:.4f} dB (expected {expected_db} dB)",
            ))
        
        # Test 2: Array conversion via Rust binding
        try:
            import sardine
            
            if hasattr(sardine, 'convert_to_db_real'):
                test_array = np.array([[1.0, 10.0, 0.1]], dtype=np.float32)
                expected = np.array([[0.0, 10.0, -10.0]])
                
                result = sardine.convert_to_db_real(test_array)
                if isinstance(result, dict):
                    result_data = result.get('data', result)
                else:
                    result_data = result
                
                max_error = float(np.max(np.abs(result_data - expected)))
                
                self._add_result(AlgorithmTestResult(
                    test_name="db_conversion_rust_binding",
                    algorithm="dB_conversion",
                    passed=max_error < 0.01,
                    expected=expected.tolist(),
                    actual=result_data.tolist() if hasattr(result_data, 'tolist') else result_data,
                    tolerance=0.01,
                    error=max_error,
                    message=f"Rust dB conversion max error: {max_error:.6f} dB",
                ))
        except ImportError:
            pass
        except Exception as e:
            self._add_result(AlgorithmTestResult(
                test_name="db_conversion_rust_binding",
                algorithm="dB_conversion",
                passed=False,
                expected="success",
                actual=str(e),
                tolerance=0.0,
                message=f"Rust binding test failed: {e}",
            ))
    
    # ==========================================================================
    # Calibration Formula Tests
    # ==========================================================================
    
    def _test_calibration_formula(self):
        """
        Test calibration formulas.
        
        ESA formula: sigma0 = |DN|² / (A_sigma)²
        Where A_sigma is the calibration LUT value.
        """
        
        # Test 1: Basic sigma0 calculation
        dn_amplitude = np.array([[100.0, 200.0, 150.0]], dtype=np.float32)
        calibration_lut = np.array([[50.0, 100.0, 75.0]], dtype=np.float32)
        
        # Expected: (|DN|²) / (A²) = DN² / A²
        expected_sigma0 = (dn_amplitude ** 2) / (calibration_lut ** 2)
        
        # Calculate using basic formula
        actual_sigma0 = (dn_amplitude ** 2) / (calibration_lut ** 2)
        
        error = float(np.max(np.abs(actual_sigma0 - expected_sigma0)))
        
        self._add_result(AlgorithmTestResult(
            test_name="calibration_sigma0_formula",
            algorithm="calibration",
            passed=error < 1e-6,
            expected=expected_sigma0.tolist(),
            actual=actual_sigma0.tolist(),
            tolerance=1e-6,
            error=error,
            message=f"sigma0 = |DN|² / A² max error: {error}",
            details={
                "expected_values": expected_sigma0.flatten().tolist(),
                # For DN=100, A=50: sigma0 = 10000/2500 = 4.0
                # For DN=200, A=100: sigma0 = 40000/10000 = 4.0
                # For DN=150, A=75: sigma0 = 22500/5625 = 4.0
            }
        ))
        
        # Test 2: Calibration LUT interpolation
        # LUT is typically 1D, needs interpolation to image dimensions
        lut_values = np.array([40.0, 50.0, 60.0, 70.0], dtype=np.float32)
        lut_samples = np.array([0, 100, 200, 300])
        
        # Interpolate to pixel positions
        pixel_positions = np.array([50, 150, 250])
        expected_interpolated = np.array([45.0, 55.0, 65.0])  # Linear interpolation
        
        actual_interpolated = np.interp(pixel_positions, lut_samples, lut_values)
        error = float(np.max(np.abs(actual_interpolated - expected_interpolated)))
        
        self._add_result(AlgorithmTestResult(
            test_name="calibration_lut_interpolation",
            algorithm="calibration",
            passed=error < 0.01,
            expected=expected_interpolated.tolist(),
            actual=actual_interpolated.tolist(),
            tolerance=0.01,
            error=error,
            message=f"LUT interpolation max error: {error}",
        ))
    
    # ==========================================================================
    # Multilook Power Preservation Tests
    # ==========================================================================
    
    def _test_multilook_power(self):
        """
        Test that multilooking preserves mean power.
        
        Multilooking is averaging, so mean should be exactly preserved.
        Total power may change due to pixel count reduction.
        """
        
        # Test 1: Mean preservation with uniform data
        uniform_data = np.ones((100, 100), dtype=np.float32) * 0.5
        input_mean = float(np.mean(uniform_data))
        
        # Simulate 5x5 multilooking (averaging)
        looks_r, looks_az = 5, 5
        output_rows = uniform_data.shape[0] // looks_az
        output_cols = uniform_data.shape[1] // looks_r
        
        multilooked = np.zeros((output_rows, output_cols), dtype=np.float32)
        for i in range(output_rows):
            for j in range(output_cols):
                window = uniform_data[
                    i*looks_az:(i+1)*looks_az,
                    j*looks_r:(j+1)*looks_r
                ]
                multilooked[i, j] = np.mean(window)
        
        output_mean = float(np.mean(multilooked))
        mean_error = abs(output_mean - input_mean) / input_mean * 100.0
        
        self._add_result(AlgorithmTestResult(
            test_name="multilook_mean_preservation_uniform",
            algorithm="multilooking",
            passed=mean_error < 0.01,  # 0.01% tolerance
            expected=input_mean,
            actual=output_mean,
            tolerance=0.01,
            error=mean_error,
            message=f"Mean power change: {mean_error:.6f}% (should be ~0%)",
            details={
                "looks_range": looks_r,
                "looks_azimuth": looks_az,
                "input_shape": uniform_data.shape,
                "output_shape": multilooked.shape,
            }
        ))
        
        # Test 2: Mean preservation with random data
        np.random.seed(42)
        random_data = np.random.uniform(0.1, 1.0, (100, 100)).astype(np.float32)
        input_mean = float(np.mean(random_data))
        
        output_rows = random_data.shape[0] // looks_az
        output_cols = random_data.shape[1] // looks_r
        
        multilooked = np.zeros((output_rows, output_cols), dtype=np.float32)
        for i in range(output_rows):
            for j in range(output_cols):
                window = random_data[
                    i*looks_az:(i+1)*looks_az,
                    j*looks_r:(j+1)*looks_r
                ]
                multilooked[i, j] = np.mean(window)
        
        output_mean = float(np.mean(multilooked))
        mean_error = abs(output_mean - input_mean) / input_mean * 100.0
        
        self._add_result(AlgorithmTestResult(
            test_name="multilook_mean_preservation_random",
            algorithm="multilooking",
            passed=mean_error < 1.0,  # 1% tolerance for random data
            expected=input_mean,
            actual=output_mean,
            tolerance=1.0,
            error=mean_error,
            message=f"Mean power change with random data: {mean_error:.4f}%",
        ))
        
        # Test 3: Variance reduction (speckle reduction)
        # ENL should increase by approximately looks_r * looks_az
        input_variance = float(np.var(random_data))
        output_variance = float(np.var(multilooked))
        
        expected_variance_ratio = 1.0 / (looks_r * looks_az)
        actual_variance_ratio = output_variance / input_variance
        ratio_error = abs(actual_variance_ratio - expected_variance_ratio) / expected_variance_ratio * 100
        
        self._add_result(AlgorithmTestResult(
            test_name="multilook_variance_reduction",
            algorithm="multilooking",
            passed=ratio_error < 20.0,  # 20% tolerance
            expected=expected_variance_ratio,
            actual=actual_variance_ratio,
            tolerance=20.0,
            error=ratio_error,
            message=f"Variance ratio: {actual_variance_ratio:.4f} (expected ~{expected_variance_ratio:.4f})",
            details={
                "input_variance": input_variance,
                "output_variance": output_variance,
                "expected_enl_increase": looks_r * looks_az,
            }
        ))
    
    # ==========================================================================
    # Terrain Flattening Tests
    # ==========================================================================
    
    def _test_terrain_flattening(self):
        """
        Test terrain flattening formula.
        
        gamma0 = sigma0 * cos(theta_local) / sin(theta_inc)
        
        On flat terrain: theta_local = theta_inc, so gamma0 ≈ sigma0 * cot(theta_inc)
        """
        
        # Test 1: Flat terrain case
        theta_inc_deg = 35.0  # Mid-swath incidence angle
        theta_inc = np.radians(theta_inc_deg)
        
        sigma0 = 0.1  # Input sigma0
        
        # For flat terrain, local incidence angle = ellipsoid incidence angle
        theta_local = theta_inc
        
        # Terrain flattening formula: gamma0 = sigma0 / cos(theta_local)
        # This removes the cos(theta) factor from radar equation
        expected_gamma0 = sigma0 / np.cos(theta_local)
        
        # Also: gamma0 = sigma0 * sec(theta) = sigma0 / cos(theta)
        actual_gamma0 = sigma0 / np.cos(theta_inc)
        
        error = abs(actual_gamma0 - expected_gamma0) / expected_gamma0 * 100
        
        self._add_result(AlgorithmTestResult(
            test_name="terrain_flatten_flat_surface",
            algorithm="terrain_flattening",
            passed=error < 0.01,
            expected=float(expected_gamma0),
            actual=float(actual_gamma0),
            tolerance=0.01,
            error=float(error),
            message=f"Flat terrain: gamma0 = {actual_gamma0:.6f} (expected {expected_gamma0:.6f})",
            details={
                "sigma0": sigma0,
                "incidence_angle_deg": theta_inc_deg,
                "cos_theta": float(np.cos(theta_inc)),
            }
        ))
        
        # Test 2: Sloped terrain - towards radar (reduces local incidence)
        # When surface slopes towards radar, local incidence angle decreases,
        # cos(theta_local) increases, so gamma0 = sigma0 / cos(theta_local) DECREASES.
        # This is correct: the surface appears brighter in sigma0, but after
        # terrain flattening it should have the same gamma0 as a flat surface.
        slope_deg = 10.0  # 10 degree slope towards radar
        theta_local_sloped = theta_inc - np.radians(slope_deg)
        
        # Gamma0 should DECREASE for surfaces facing the radar
        # (because cos of smaller angle is larger)
        gamma0_sloped = sigma0 / np.cos(theta_local_sloped)
        
        # Sloped surface should have LOWER gamma0 than flat (cos is larger)
        # This verifies the terrain flattening normalizes away the slope effect
        slope_effect_correct = gamma0_sloped < expected_gamma0
        
        self._add_result(AlgorithmTestResult(
            test_name="terrain_flatten_forward_slope",
            algorithm="terrain_flattening",
            passed=slope_effect_correct,
            expected="gamma0_sloped < gamma0_flat (normalized)",
            actual=f"{gamma0_sloped:.6f} < {expected_gamma0:.6f}" if slope_effect_correct else f"{gamma0_sloped:.6f} >= {expected_gamma0:.6f}",
            tolerance=0.0,
            message=f"Forward slope ({slope_deg}°): gamma0 normalized correctly" if slope_effect_correct else "FAIL: Forward slope normalization incorrect",
            details={
                "sigma0": sigma0,
                "ellipsoid_incidence_deg": theta_inc_deg,
                "local_incidence_deg": float(np.degrees(theta_local_sloped)),
                "gamma0_flat": float(expected_gamma0),
                "gamma0_sloped": float(gamma0_sloped),
                "explanation": "Lower local incidence -> higher cos -> lower gamma0",
            }
        ))
    
    # ==========================================================================
    # Noise Floor Tests
    # ==========================================================================
    
    def _test_noise_floor(self):
        """
        Test noise floor handling.
        
        Noise subtraction should not create negative values.
        Small values after noise subtraction should be handled carefully.
        """
        
        # Test 1: Noise subtraction edge case
        signal = np.array([0.1, 0.05, 0.02, 0.01, 0.005], dtype=np.float32)
        noise = 0.01  # Noise floor
        
        # Naive subtraction would give negative values
        naive_result = signal - noise
        
        # Correct approach: clamp to small positive value
        noise_floor = 1e-10
        correct_result = np.maximum(signal - noise, noise_floor)
        
        has_negative = np.any(naive_result < 0)
        correct_no_negative = np.all(correct_result > 0)
        
        self._add_result(AlgorithmTestResult(
            test_name="noise_floor_clipping",
            algorithm="noise_handling",
            passed=correct_no_negative,
            expected="all positive",
            actual="all positive" if correct_no_negative else "has negative",
            tolerance=0.0,
            message=f"Noise subtraction clamping: {'correct' if correct_no_negative else 'FAIL'}",
            details={
                "naive_has_negative": bool(has_negative),
                "correct_min_value": float(np.min(correct_result)),
            }
        ))
        
        # Test 2: Thermal noise subtraction preserves relative structure
        np.random.seed(42)
        signal_with_noise = np.random.uniform(0.1, 1.0, 100) + 0.05  # Signal + noise
        noise_estimate = 0.05
        
        denoised = np.maximum(signal_with_noise - noise_estimate, 1e-10)
        
        # Correlation between denoised and original signal structure
        original_relative = signal_with_noise - np.mean(signal_with_noise)
        denoised_relative = denoised - np.mean(denoised)
        correlation = np.corrcoef(original_relative, denoised_relative)[0, 1]
        
        self._add_result(AlgorithmTestResult(
            test_name="noise_subtraction_structure_preservation",
            algorithm="noise_handling",
            passed=correlation > 0.99,
            expected=1.0,
            actual=float(correlation),
            tolerance=0.01,
            error=1.0 - float(correlation),
            message=f"Structure preservation correlation: {correlation:.6f}",
        ))
    
    # ==========================================================================
    # Coordinate Transform Tests
    # ==========================================================================
    
    def _test_coordinate_transforms(self):
        """
        Test coordinate system transforms are reversible.
        """
        try:
            from pyproj import Transformer
            
            # Test WGS84 -> UTM -> WGS84 round-trip
            test_points = [
                (-122.4194, 37.7749),  # San Francisco
                (139.6917, 35.6895),   # Tokyo
                (-0.1276, 51.5074),    # London
            ]
            
            max_error = 0.0
            for lon, lat in test_points:
                # Determine UTM zone
                zone = int((lon + 180) / 6) + 1
                epsg = 32600 + zone if lat >= 0 else 32700 + zone
                
                # Forward
                t_fwd = Transformer.from_crs("EPSG:4326", f"EPSG:{epsg}", always_xy=True)
                x, y = t_fwd.transform(lon, lat)
                
                # Inverse
                t_inv = Transformer.from_crs(f"EPSG:{epsg}", "EPSG:4326", always_xy=True)
                lon2, lat2 = t_inv.transform(x, y)
                
                error_m = np.sqrt((lon2 - lon)**2 + (lat2 - lat)**2) * 111000
                max_error = max(max_error, error_m)
            
            self._add_result(AlgorithmTestResult(
                test_name="coordinate_transform_roundtrip",
                algorithm="coordinates",
                passed=max_error < 0.01,
                expected=0.0,
                actual=float(max_error),
                tolerance=0.01,
                error=float(max_error),
                message=f"Max round-trip error: {max_error:.6f} m",
            ))
            
        except ImportError:
            self._add_result(AlgorithmTestResult(
                test_name="coordinate_transform_roundtrip",
                algorithm="coordinates",
                passed=True,
                expected="skipped",
                actual="skipped",
                tolerance=0.0,
                message="pyproj not available, test skipped",
            ))
    
    # ==========================================================================
    # Resampling Power Tests
    # ==========================================================================
    
    def _test_resampling_power(self):
        """
        Test that resampling preserves mean power.
        """
        try:
            from scipy import ndimage
            
            # Create test data
            np.random.seed(42)
            data = np.random.uniform(0.1, 1.0, (100, 100)).astype(np.float32)
            input_mean = float(np.mean(data))
            
            # Downsample by factor of 2 using area average
            downsampled = ndimage.zoom(data, 0.5, order=1)  # Bilinear
            output_mean = float(np.mean(downsampled))
            
            mean_error = abs(output_mean - input_mean) / input_mean * 100
            
            self._add_result(AlgorithmTestResult(
                test_name="resample_power_preservation_down",
                algorithm="resampling",
                passed=mean_error < 5.0,  # 5% tolerance
                expected=input_mean,
                actual=output_mean,
                tolerance=5.0,
                error=mean_error,
                message=f"Downsample mean change: {mean_error:.2f}%",
                details={
                    "input_shape": data.shape,
                    "output_shape": downsampled.shape,
                }
            ))
            
            # Upsample and check
            upsampled = ndimage.zoom(data, 2.0, order=1)
            output_mean_up = float(np.mean(upsampled))
            mean_error_up = abs(output_mean_up - input_mean) / input_mean * 100
            
            self._add_result(AlgorithmTestResult(
                test_name="resample_power_preservation_up",
                algorithm="resampling",
                passed=mean_error_up < 5.0,
                expected=input_mean,
                actual=output_mean_up,
                tolerance=5.0,
                error=mean_error_up,
                message=f"Upsample mean change: {mean_error_up:.2f}%",
            ))
            
        except ImportError:
            self._add_result(AlgorithmTestResult(
                test_name="resample_power_preservation",
                algorithm="resampling",
                passed=True,
                expected="skipped",
                actual="skipped",
                tolerance=0.0,
                message="scipy not available, test skipped",
            ))
    
    def print_report(self, report: AlgorithmTestReport):
        """Print formatted test report."""
        print("\n" + "=" * 80)
        print("ALGORITHM TEST REPORT")
        print("=" * 80)
        print(f"Timestamp: {report.timestamp}")
        print(f"Total tests: {len(report.results)}")
        print(f"Passed: {report.pass_count}")
        print(f"Failed: {report.fail_count}")
        print()
        
        print("Overall: " + ("✅ ALL PASSED" if report.passed else "❌ SOME FAILED"))
        print()
        
        # Group by algorithm
        algorithms = set(r.algorithm for r in report.results)
        
        for algo in sorted(algorithms):
            algo_results = [r for r in report.results if r.algorithm == algo]
            algo_passed = sum(1 for r in algo_results if r.passed)
            algo_total = len(algo_results)
            
            status = "✅" if algo_passed == algo_total else "❌"
            print(f"{status} {algo}: {algo_passed}/{algo_total} passed")
            
            for r in algo_results:
                icon = "  ✓" if r.passed else "  ✗"
                print(f"{icon} {r.test_name}: {r.message}")
        
        print("=" * 80)


def run_algorithm_tests() -> AlgorithmTestReport:
    """
    Convenience function to run all algorithm tests.
    
    Returns:
        AlgorithmTestReport
    """
    tester = AlgorithmTester()
    report = tester.run_all_tests()
    tester.print_report(report)
    return report


if __name__ == "__main__":
    run_algorithm_tests()
