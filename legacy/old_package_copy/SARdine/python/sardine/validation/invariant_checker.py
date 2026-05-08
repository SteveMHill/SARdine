"""
Invariant Checker for SARdine Pipeline

Validates that processing invariants are maintained throughout the pipeline:
- Power/energy preservation across processing steps
- Spatial extent preservation (no unexplained data loss)
- Radiometric consistency (calibration preserved)
- Metadata consistency

These invariants should hold regardless of input data.
"""

import logging
import numpy as np
from dataclasses import dataclass, field
from typing import Dict, Any, List, Optional, Tuple, Callable
from pathlib import Path

logger = logging.getLogger(__name__)


@dataclass
class InvariantCheckResult:
    """Result of an invariant check."""
    invariant: str
    step: str
    passed: bool
    expected: float
    actual: float
    tolerance: float
    severity: str  # critical, major, minor
    message: str
    details: Dict[str, Any] = field(default_factory=dict)


@dataclass
class InvariantReport:
    """Complete invariant validation report."""
    results: List[InvariantCheckResult] = field(default_factory=list)
    
    @property
    def passed(self) -> bool:
        """All critical and major invariants passed."""
        for r in self.results:
            if not r.passed and r.severity in ("critical", "major"):
                return False
        return True
    
    @property
    def violations(self) -> List[InvariantCheckResult]:
        """Get all violations."""
        return [r for r in self.results if not r.passed]


class InvariantChecker:
    """
    Check that processing invariants are maintained.
    
    Key invariants for SAR backscatter:
    1. Power preservation: Mean power should be preserved (not lost or amplified)
    2. Valid pixel preservation: Should not lose valid pixels unexpectedly
    3. Radiometric accuracy: Calibration values must be applied correctly
    4. Spatial consistency: Geographic bounds should be preserved
    5. Value range: Data should remain in physically plausible range
    """
    
    DEFAULT_TOLERANCES = {
        # Power can change by at most this % per step
        "power_change_percent": 5.0,
        
        # Valid pixels can decrease by at most this % per step
        # (Some loss is expected from resampling, masking)
        "valid_pixel_loss_percent": 10.0,
        
        # dB values should stay in this range
        "min_db_value": -50.0,
        "max_db_value": 30.0,
        
        # Linear values should stay in this range
        "min_linear_value": 0.0,
        "max_linear_value": 100.0,
        
        # Spatial extent change (0.1% tolerance)
        "extent_change_percent": 0.1,
    }
    
    def __init__(self, tolerances: Optional[Dict[str, float]] = None):
        """
        Initialize checker.
        
        Args:
            tolerances: Custom tolerances (overrides defaults)
        """
        self.tolerances = self.DEFAULT_TOLERANCES.copy()
        if tolerances:
            self.tolerances.update(tolerances)
        
        self.step_data: Dict[str, Dict[str, Any]] = {}
        self.results: List[InvariantCheckResult] = []
    
    def record_step(
        self,
        step_name: str,
        data: np.ndarray,
        metadata: Optional[Dict[str, Any]] = None
    ):
        """
        Record data at a processing step for invariant checking.
        
        Call this before and after each processing step to enable
        cross-step invariant validation.
        
        Args:
            step_name: Name of the processing step (e.g., "after_multilook")
            data: The numpy array at this step
            metadata: Optional metadata dict
        """
        finite_mask = np.isfinite(data)
        finite_data = data[finite_mask]
        
        self.step_data[step_name] = {
            "shape": data.shape,
            "dtype": str(data.dtype),
            "total_pixels": data.size,
            "valid_pixels": int(np.sum(finite_mask)),
            "mean": float(np.mean(finite_data)) if len(finite_data) > 0 else np.nan,
            "std": float(np.std(finite_data)) if len(finite_data) > 0 else np.nan,
            "min": float(np.min(finite_data)) if len(finite_data) > 0 else np.nan,
            "max": float(np.max(finite_data)) if len(finite_data) > 0 else np.nan,
            "sum": float(np.sum(finite_data)) if len(finite_data) > 0 else 0.0,
            "metadata": metadata or {},
        }
    
    def check_power_preservation(
        self,
        step_before: str,
        step_after: str,
        expected_factor: float = 1.0,
    ) -> InvariantCheckResult:
        """
        Check that mean power is preserved between two steps.
        
        For most steps, expected_factor = 1.0 (power preserved).
        For multilooking, power should be exactly preserved (averaging).
        For terrain flattening, power can change based on local incidence angle.
        
        Args:
            step_before: Name of step before processing
            step_after: Name of step after processing
            expected_factor: Expected power ratio (after/before)
            
        Returns:
            InvariantCheckResult
        """
        if step_before not in self.step_data or step_after not in self.step_data:
            return InvariantCheckResult(
                invariant="power_preservation",
                step=f"{step_before} -> {step_after}",
                passed=False,
                expected=expected_factor,
                actual=np.nan,
                tolerance=self.tolerances["power_change_percent"] / 100.0,
                severity="critical",
                message=f"Missing step data for comparison",
            )
        
        before = self.step_data[step_before]
        after = self.step_data[step_after]
        
        mean_before = before["mean"]
        mean_after = after["mean"]
        
        if np.isnan(mean_before) or np.isnan(mean_after):
            return InvariantCheckResult(
                invariant="power_preservation",
                step=f"{step_before} -> {step_after}",
                passed=False,
                expected=expected_factor,
                actual=np.nan,
                tolerance=self.tolerances["power_change_percent"] / 100.0,
                severity="critical",
                message="NaN values in power calculation",
            )
        
        actual_factor = mean_after / mean_before if mean_before != 0 else np.inf
        deviation = abs(actual_factor - expected_factor) / expected_factor
        tolerance = self.tolerances["power_change_percent"] / 100.0
        
        passed = deviation <= tolerance
        
        result = InvariantCheckResult(
            invariant="power_preservation",
            step=f"{step_before} -> {step_after}",
            passed=passed,
            expected=expected_factor,
            actual=actual_factor,
            tolerance=tolerance,
            severity="critical",
            message=f"Power ratio: {actual_factor:.4f} (expected {expected_factor:.4f}, "
                   f"deviation: {deviation*100:.2f}%)",
            details={
                "mean_before": mean_before,
                "mean_after": mean_after,
                "deviation_percent": deviation * 100,
            }
        )
        
        self.results.append(result)
        return result
    
    def check_valid_pixel_preservation(
        self,
        step_before: str,
        step_after: str,
        allow_increase: bool = False,
    ) -> InvariantCheckResult:
        """
        Check that valid pixels are not unexpectedly lost.
        
        Some pixel loss is expected from:
        - Edge effects in filtering
        - Masking operations
        - Resampling to different resolution
        
        But excessive loss indicates a bug.
        
        Args:
            step_before: Name of step before processing
            step_after: Name of step after processing
            allow_increase: If True, also allows valid pixels to increase
            
        Returns:
            InvariantCheckResult
        """
        if step_before not in self.step_data or step_after not in self.step_data:
            return InvariantCheckResult(
                invariant="valid_pixel_preservation",
                step=f"{step_before} -> {step_after}",
                passed=False,
                expected=0.0,
                actual=np.nan,
                tolerance=self.tolerances["valid_pixel_loss_percent"],
                severity="major",
                message="Missing step data for comparison",
            )
        
        before = self.step_data[step_before]
        after = self.step_data[step_after]
        
        valid_before = before["valid_pixels"]
        valid_after = after["valid_pixels"]
        
        if valid_before == 0:
            return InvariantCheckResult(
                invariant="valid_pixel_preservation",
                step=f"{step_before} -> {step_after}",
                passed=False,
                expected=0.0,
                actual=100.0,
                tolerance=self.tolerances["valid_pixel_loss_percent"],
                severity="critical",
                message="Zero valid pixels before step",
            )
        
        loss_percent = (valid_before - valid_after) / valid_before * 100.0
        tolerance = self.tolerances["valid_pixel_loss_percent"]
        
        # Check if loss is within tolerance
        if loss_percent < 0:  # Pixels increased
            passed = allow_increase
            message = f"Valid pixels increased by {-loss_percent:.2f}%"
        else:
            passed = loss_percent <= tolerance
            message = f"Valid pixel loss: {loss_percent:.2f}% (tolerance: {tolerance}%)"
        
        result = InvariantCheckResult(
            invariant="valid_pixel_preservation",
            step=f"{step_before} -> {step_after}",
            passed=passed,
            expected=0.0,
            actual=loss_percent,
            tolerance=tolerance,
            severity="major",
            message=message,
            details={
                "valid_before": valid_before,
                "valid_after": valid_after,
                "loss_percent": loss_percent,
            }
        )
        
        self.results.append(result)
        return result
    
    def check_value_range(
        self,
        step_name: str,
        is_db: bool = False,
    ) -> InvariantCheckResult:
        """
        Check that values are in physically plausible range.
        
        Args:
            step_name: Name of the step to check
            is_db: True if values are in dB, False if linear
            
        Returns:
            InvariantCheckResult
        """
        if step_name not in self.step_data:
            return InvariantCheckResult(
                invariant="value_range",
                step=step_name,
                passed=False,
                expected=0.0,
                actual=np.nan,
                tolerance=0.0,
                severity="critical",
                message="Missing step data",
            )
        
        data = self.step_data[step_name]
        data_min = data["min"]
        data_max = data["max"]
        
        if is_db:
            expected_min = self.tolerances["min_db_value"]
            expected_max = self.tolerances["max_db_value"]
            unit = "dB"
        else:
            expected_min = self.tolerances["min_linear_value"]
            expected_max = self.tolerances["max_linear_value"]
            unit = "linear"
        
        in_range = data_min >= expected_min and data_max <= expected_max
        
        result = InvariantCheckResult(
            invariant="value_range",
            step=step_name,
            passed=in_range,
            expected=expected_max,
            actual=data_max,
            tolerance=0.0,
            severity="major" if not in_range else "info",
            message=f"Value range: [{data_min:.2f}, {data_max:.2f}] {unit} "
                   f"(expected [{expected_min}, {expected_max}])",
            details={
                "min": data_min,
                "max": data_max,
                "expected_min": expected_min,
                "expected_max": expected_max,
                "is_db": is_db,
            }
        )
        
        self.results.append(result)
        return result
    
    def check_multilook_invariant(
        self,
        input_data: np.ndarray,
        output_data: np.ndarray,
        range_looks: int,
        azimuth_looks: int,
    ) -> InvariantCheckResult:
        """
        Check multilooking preserves mean power.
        
        Multilooking is averaging, so mean should be exactly preserved.
        
        Args:
            input_data: Data before multilooking
            output_data: Data after multilooking
            range_looks: Number of range looks
            azimuth_looks: Number of azimuth looks
            
        Returns:
            InvariantCheckResult
        """
        finite_in = input_data[np.isfinite(input_data)]
        finite_out = output_data[np.isfinite(output_data)]
        
        if len(finite_in) == 0 or len(finite_out) == 0:
            return InvariantCheckResult(
                invariant="multilook_power_preservation",
                step="multilooking",
                passed=False,
                expected=1.0,
                actual=np.nan,
                tolerance=0.01,
                severity="critical",
                message="No valid pixels for comparison",
            )
        
        mean_in = np.mean(finite_in)
        mean_out = np.mean(finite_out)
        
        if mean_in == 0:
            return InvariantCheckResult(
                invariant="multilook_power_preservation",
                step="multilooking",
                passed=False,
                expected=1.0,
                actual=np.nan,
                tolerance=0.01,
                severity="critical",
                message="Zero mean input",
            )
        
        ratio = mean_out / mean_in
        tolerance = 0.01  # 1% tolerance
        passed = abs(ratio - 1.0) <= tolerance
        
        result = InvariantCheckResult(
            invariant="multilook_power_preservation",
            step="multilooking",
            passed=passed,
            expected=1.0,
            actual=ratio,
            tolerance=tolerance,
            severity="critical",
            message=f"Mean power ratio: {ratio:.6f} (should be 1.0 ±{tolerance})",
            details={
                "mean_input": float(mean_in),
                "mean_output": float(mean_out),
                "range_looks": range_looks,
                "azimuth_looks": azimuth_looks,
                "expected_pixel_reduction": range_looks * azimuth_looks,
                "actual_pixel_reduction": len(finite_in) / len(finite_out) if len(finite_out) > 0 else np.inf,
            }
        )
        
        self.results.append(result)
        return result
    
    def check_calibration_invariant(
        self,
        raw_amplitude: np.ndarray,
        calibrated_data: np.ndarray,
        calibration_lut: np.ndarray,
    ) -> InvariantCheckResult:
        """
        Check that calibration formula was applied correctly.
        
        For sigma0 calibration: sigma0 = |DN|² / (calibration_lut)²
        
        Args:
            raw_amplitude: Raw amplitude data (|DN|)
            calibrated_data: Calibrated output (should be sigma0)
            calibration_lut: Calibration LUT values
            
        Returns:
            InvariantCheckResult
        """
        # Calculate expected calibrated values
        with np.errstate(divide='ignore', invalid='ignore'):
            expected = (raw_amplitude ** 2) / (calibration_lut ** 2)
        
        # Compare at valid locations
        valid_mask = np.isfinite(expected) & np.isfinite(calibrated_data) & (calibration_lut > 0)
        
        if np.sum(valid_mask) == 0:
            return InvariantCheckResult(
                invariant="calibration_formula",
                step="calibration",
                passed=False,
                expected=0.0,
                actual=np.nan,
                tolerance=0.01,
                severity="critical",
                message="No valid pixels for calibration check",
            )
        
        expected_valid = expected[valid_mask]
        actual_valid = calibrated_data[valid_mask]
        
        # Calculate relative error
        with np.errstate(divide='ignore', invalid='ignore'):
            relative_error = np.abs(actual_valid - expected_valid) / expected_valid
        
        mean_error = float(np.nanmean(relative_error))
        max_error = float(np.nanmax(relative_error))
        
        tolerance = 0.01  # 1% error tolerance
        passed = mean_error <= tolerance
        
        result = InvariantCheckResult(
            invariant="calibration_formula",
            step="calibration",
            passed=passed,
            expected=0.0,
            actual=mean_error,
            tolerance=tolerance,
            severity="critical",
            message=f"Calibration error: mean={mean_error*100:.4f}%, max={max_error*100:.2f}%",
            details={
                "mean_relative_error": mean_error,
                "max_relative_error": max_error,
                "valid_pixels_checked": int(np.sum(valid_mask)),
            }
        )
        
        self.results.append(result)
        return result
    
    def check_db_conversion_invariant(
        self,
        linear_data: np.ndarray,
        db_data: np.ndarray,
    ) -> InvariantCheckResult:
        """
        Check that dB conversion follows: dB = 10 * log10(linear).
        
        Args:
            linear_data: Linear power data
            db_data: dB-converted data
            
        Returns:
            InvariantCheckResult
        """
        # Calculate expected dB values
        with np.errstate(divide='ignore', invalid='ignore'):
            expected_db = 10.0 * np.log10(linear_data)
        
        # Compare at valid locations
        valid_mask = (
            np.isfinite(linear_data) & 
            np.isfinite(db_data) & 
            np.isfinite(expected_db) &
            (linear_data > 0)
        )
        
        if np.sum(valid_mask) == 0:
            return InvariantCheckResult(
                invariant="db_conversion",
                step="db_conversion",
                passed=False,
                expected=0.0,
                actual=np.nan,
                tolerance=0.01,
                severity="critical",
                message="No valid pixels for dB conversion check",
            )
        
        expected_valid = expected_db[valid_mask]
        actual_valid = db_data[valid_mask]
        
        diff = np.abs(actual_valid - expected_valid)
        mean_diff = float(np.mean(diff))
        max_diff = float(np.max(diff))
        
        tolerance = 0.001  # 0.001 dB tolerance
        passed = mean_diff <= tolerance
        
        result = InvariantCheckResult(
            invariant="db_conversion",
            step="db_conversion",
            passed=passed,
            expected=0.0,
            actual=mean_diff,
            tolerance=tolerance,
            severity="critical",
            message=f"dB conversion error: mean={mean_diff:.6f} dB, max={max_diff:.4f} dB",
            details={
                "mean_error_db": mean_diff,
                "max_error_db": max_diff,
                "valid_pixels_checked": int(np.sum(valid_mask)),
            }
        )
        
        self.results.append(result)
        return result
    
    def get_report(self) -> InvariantReport:
        """Get the invariant check report."""
        return InvariantReport(results=self.results)
    
    def print_report(self, report: Optional[InvariantReport] = None):
        """Print formatted invariant report."""
        if report is None:
            report = self.get_report()
        
        print("\n" + "=" * 80)
        print("INVARIANT CHECK REPORT")
        print("=" * 80)
        print(f"Total checks: {len(report.results)}")
        print(f"Passed: {sum(1 for r in report.results if r.passed)}")
        print(f"Failed: {len(report.violations)}")
        print()
        
        print("Overall: " + ("✅ PASSED" if report.passed else "❌ FAILED"))
        print()
        
        if report.violations:
            print("VIOLATIONS:")
            print("-" * 60)
            for r in report.violations:
                severity_icon = "🔴" if r.severity == "critical" else "🟡" if r.severity == "major" else "🔵"
                print(f"{severity_icon} [{r.invariant}] {r.step}")
                print(f"   {r.message}")
                print(f"   Expected: {r.expected:.6f}, Actual: {r.actual:.6f}, Tolerance: {r.tolerance}")
                print()
        
        print("ALL CHECKS:")
        print("-" * 60)
        for r in report.results:
            status = "✅" if r.passed else "❌"
            print(f"{status} [{r.invariant}] {r.step}: {r.message}")
        
        print("=" * 80)
    
    def reset(self):
        """Reset all recorded data and results."""
        self.step_data.clear()
        self.results.clear()
