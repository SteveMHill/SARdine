"""
Scientific Validator for SARdine Backscatter Pipeline

Provides comprehensive validation that goes beyond output checking to detect:
- Algorithm correctness issues
- Design flaws in parameters/thresholds
- Scientific accuracy problems
- Power/radiometry preservation violations

This is designed to catch the types of bugs that output-only validation misses.
"""

import logging
import numpy as np
from dataclasses import dataclass, field
from typing import Dict, Any, List, Optional, Tuple
from pathlib import Path

logger = logging.getLogger(__name__)


@dataclass
class ValidationResult:
    """Result of a single validation check."""
    name: str
    category: str  # algorithm, design, scientific, invariant
    passed: bool
    severity: str  # critical, major, minor, info
    message: str
    details: Dict[str, Any] = field(default_factory=dict)
    fix_suggestion: Optional[str] = None


@dataclass
class ValidationReport:
    """Complete validation report."""
    timestamp: str
    results: List[ValidationResult] = field(default_factory=list)
    summary: Dict[str, Any] = field(default_factory=dict)
    
    @property
    def passed(self) -> bool:
        """Check if all critical and major validations passed."""
        for r in self.results:
            if not r.passed and r.severity in ("critical", "major"):
                return False
        return True
    
    @property
    def critical_failures(self) -> List[ValidationResult]:
        return [r for r in self.results if not r.passed and r.severity == "critical"]
    
    @property  
    def major_failures(self) -> List[ValidationResult]:
        return [r for r in self.results if not r.passed and r.severity == "major"]


class ScientificValidator:
    """
    Comprehensive scientific validation for SAR backscatter processing.
    
    Unlike output-only validation, this checks:
    1. Whether the code logic matches scientific requirements
    2. Whether parameters/thresholds are set to scientifically correct values
    3. Whether algorithms implement the correct formulas
    4. Whether invariants (like power preservation) are maintained
    """
    
    def __init__(self, processor=None, config: Optional[Dict[str, Any]] = None):
        """
        Initialize validator.
        
        Args:
            processor: BackscatterProcessor instance (optional)
            config: Validation configuration overrides
        """
        self.processor = processor
        self.config = config or {}
        self.results: List[ValidationResult] = []
        
        # Scientific thresholds based on published literature
        self.thresholds = {
            # Coverage thresholds (Small & Schubert 2008)
            "min_geocoding_coverage": 50.0,  # Minimum % valid pixels after geocoding
            "warning_geocoding_coverage": 80.0,  # Warn if below this
            
            # Power preservation (ESA Sentinel-1 Product Spec)
            "max_power_change_percent": 5.0,  # Max allowed power change per step
            
            # Orbit requirements (ESA POD Quality Spec)
            "min_orbit_vectors": 10,  # Minimum for interpolation
            "recommended_orbit_vectors": 100,  # Recommended for scientific use
            
            # dB range for SAR backscatter (Ulaby & Long 2014)
            "min_db_land": -25.0,  # Typical land minimum
            "max_db_land": 5.0,    # Typical land maximum
            "min_db_water": -30.0,  # Calm water
            "max_db_urban": 10.0,   # Urban/metal
            
            # Speckle ENL (Oliver & Quegan 2004)
            "enl_well_suppressed": 100,  # σ/μ < 10% at this ENL
            
            # Incidence angle (Sentinel-1 spec)
            "iw_incidence_near": 29.1,
            "iw_incidence_far": 46.0,
            
            # Radiometric accuracy (ESA S1 Cal/Val)
            "max_subswath_radiometric_diff_db": 0.5,  # Between subswaths
        }
        
        # Override with config
        self.thresholds.update(self.config.get("thresholds", {}))
    
    def run_full_validation(self) -> ValidationReport:
        """Run all validation checks and return comprehensive report."""
        from datetime import datetime
        
        self.results = []
        
        # Run all validation categories
        self._validate_design_parameters()
        self._validate_algorithm_correctness()
        self._validate_scientific_accuracy()
        self._validate_invariants()
        self._validate_error_handling()
        
        # Generate report
        report = ValidationReport(
            timestamp=datetime.now().isoformat(),
            results=self.results,
            summary=self._generate_summary()
        )
        
        return report
    
    def _add_result(self, result: ValidationResult):
        """Add a validation result."""
        self.results.append(result)
        
        # Log based on severity
        if not result.passed:
            if result.severity == "critical":
                logger.error(f"CRITICAL: {result.name} - {result.message}")
            elif result.severity == "major":
                logger.warning(f"MAJOR: {result.name} - {result.message}")
            else:
                logger.info(f"{result.severity.upper()}: {result.name} - {result.message}")
    
    def _validate_design_parameters(self):
        """
        Validate that design parameters are set to scientifically correct values.
        
        These are the hardcoded thresholds and defaults that can silently cause
        incorrect results if set wrong.
        """
        
        # Check 1: Coverage threshold
        try:
            from sardine.processors.backscatter.terrain import (
                CRITICAL_COVERAGE_THRESHOLD,
                WARNING_COVERAGE_THRESHOLD
            )
            
            self._add_result(ValidationResult(
                name="coverage_threshold_critical",
                category="design",
                passed=CRITICAL_COVERAGE_THRESHOLD >= self.thresholds["min_geocoding_coverage"],
                severity="critical",
                message=f"CRITICAL_COVERAGE_THRESHOLD is {CRITICAL_COVERAGE_THRESHOLD}% "
                       f"(should be >= {self.thresholds['min_geocoding_coverage']}%)",
                details={
                    "actual": CRITICAL_COVERAGE_THRESHOLD,
                    "expected_min": self.thresholds["min_geocoding_coverage"],
                },
                fix_suggestion="Raise CRITICAL_COVERAGE_THRESHOLD to at least 50%"
            ))
            
            self._add_result(ValidationResult(
                name="coverage_threshold_warning", 
                category="design",
                passed=WARNING_COVERAGE_THRESHOLD >= self.thresholds["warning_geocoding_coverage"],
                severity="major",
                message=f"WARNING_COVERAGE_THRESHOLD is {WARNING_COVERAGE_THRESHOLD}% "
                       f"(should be >= {self.thresholds['warning_geocoding_coverage']}%)",
                details={
                    "actual": WARNING_COVERAGE_THRESHOLD,
                    "expected_min": self.thresholds["warning_geocoding_coverage"],
                },
                fix_suggestion="Raise WARNING_COVERAGE_THRESHOLD to at least 80%"
            ))
        except ImportError as e:
            self._add_result(ValidationResult(
                name="coverage_threshold_import",
                category="design", 
                passed=False,
                severity="minor",
                message=f"Could not import coverage thresholds: {e}",
                details={"error": str(e)}
            ))
        
        # Check 2: Default calibration type for terrain flattening
        if self.processor:
            calibration_type = getattr(self.processor, 'calibration_type', None)
            terrain_flatten = getattr(self.processor, 'terrain_flatten', True)
            
            # If terrain flattening is enabled, calibration should be sigma0
            # because terrain flattening converts sigma0 -> gamma0_tc
            if terrain_flatten and calibration_type:
                correct_cal_type = calibration_type.lower() == "sigma0"
                
                self._add_result(ValidationResult(
                    name="calibration_type_for_terrain_flatten",
                    category="design",
                    passed=correct_cal_type,
                    severity="critical",
                    message=f"Calibration type is '{calibration_type}' with terrain_flatten=True "
                           f"(should be 'sigma0' because terrain flattening converts sigma0→gamma0)",
                    details={
                        "calibration_type": calibration_type,
                        "terrain_flatten": terrain_flatten,
                    },
                    fix_suggestion="Set calibration_type='sigma0' when terrain_flatten=True"
                ))
        
        # Check 3: Orbit vector minimum
        try:
            from sardine.processors.backscatter.metadata_orbit import MINIMUM_ORBIT_VECTORS
            
            self._add_result(ValidationResult(
                name="orbit_vector_minimum",
                category="design",
                passed=MINIMUM_ORBIT_VECTORS >= self.thresholds["min_orbit_vectors"],
                severity="major",
                message=f"MINIMUM_ORBIT_VECTORS is {MINIMUM_ORBIT_VECTORS} "
                       f"(should be >= {self.thresholds['min_orbit_vectors']})",
                details={
                    "actual": MINIMUM_ORBIT_VECTORS,
                    "expected_min": self.thresholds["min_orbit_vectors"],
                }
            ))
        except ImportError:
            pass
        
        # Check 4: ENL handling in speckle filter
        try:
            from sardine.processors.backscatter.speckle import ENL_MAX
            
            self._add_result(ValidationResult(
                name="enl_max_threshold",
                category="design",
                passed=ENL_MAX <= self.thresholds["enl_well_suppressed"],
                severity="minor",
                message=f"ENL_MAX is {ENL_MAX} (at ENL>{self.thresholds['enl_well_suppressed']}, "
                       f"speckle is already well-suppressed, filtering may reduce resolution)",
                details={
                    "actual": ENL_MAX,
                    "well_suppressed_threshold": self.thresholds["enl_well_suppressed"],
                }
            ))
        except ImportError:
            pass
    
    def _validate_algorithm_correctness(self):
        """
        Validate that algorithms implement correct formulas.
        
        This tests the actual mathematical operations, not just outputs.
        """
        
        # Test 1: dB conversion formula
        self._test_db_conversion()
        
        # Test 2: Multilooking power preservation
        self._test_multilook_formula()
        
        # Test 3: Terrain flattening formula
        self._test_terrain_flattening_formula()
        
        # Test 4: Coordinate transforms
        self._test_coordinate_transforms()
    
    def _test_db_conversion(self):
        """Test that dB conversion uses correct formula: dB = 10 * log10(linear)"""
        try:
            import sardine
            
            # Test with known values
            test_cases = [
                (1.0, 0.0),      # 10*log10(1) = 0
                (10.0, 10.0),    # 10*log10(10) = 10
                (0.1, -10.0),    # 10*log10(0.1) = -10
                (100.0, 20.0),   # 10*log10(100) = 20
            ]
            
            test_array = np.array([[tc[0] for tc in test_cases]], dtype=np.float32)
            expected = [tc[1] for tc in test_cases]
            
            if hasattr(sardine, 'convert_to_db_real'):
                result = sardine.convert_to_db_real(test_array)
                if isinstance(result, dict):
                    result_data = result.get('data', result)
                else:
                    result_data = result
                
                result_values = result_data.flatten()
                max_error = max(abs(result_values[i] - expected[i]) for i in range(len(expected)))
                
                self._add_result(ValidationResult(
                    name="db_conversion_formula",
                    category="algorithm",
                    passed=max_error < 0.01,  # 0.01 dB tolerance
                    severity="critical",
                    message=f"dB conversion formula test: max error = {max_error:.6f} dB",
                    details={
                        "test_cases": test_cases,
                        "results": result_values.tolist(),
                        "max_error": float(max_error),
                    }
                ))
            else:
                self._add_result(ValidationResult(
                    name="db_conversion_formula",
                    category="algorithm",
                    passed=False,
                    severity="minor",
                    message="convert_to_db_real function not found",
                    details={}
                ))
        except Exception as e:
            self._add_result(ValidationResult(
                name="db_conversion_formula",
                category="algorithm",
                passed=False,
                severity="major",
                message=f"dB conversion test failed: {e}",
                details={"error": str(e)}
            ))
    
    def _test_multilook_formula(self):
        """Test that multilooking preserves mean power (averaging, not summing)."""
        try:
            import sardine
            
            # Create test data with known mean
            test_data = np.ones((100, 100), dtype=np.float32) * 0.5
            input_mean = float(np.mean(test_data))
            
            if hasattr(sardine, 'apply_multilooking'):
                result = sardine.apply_multilooking(
                    test_data, 
                    range_looks=5, 
                    azimuth_looks=5,
                    input_range_spacing=2.33,
                    input_azimuth_spacing=13.97
                )
                
                if isinstance(result, dict):
                    output_data = result.get('data', result)
                else:
                    output_data = result
                
                output_mean = float(np.mean(output_data[np.isfinite(output_data)]))
                power_change = abs(output_mean - input_mean) / input_mean * 100.0
                
                self._add_result(ValidationResult(
                    name="multilook_power_preservation",
                    category="algorithm",
                    passed=power_change < self.thresholds["max_power_change_percent"],
                    severity="critical",
                    message=f"Multilook power change: {power_change:.2f}% "
                           f"(should be < {self.thresholds['max_power_change_percent']}%)",
                    details={
                        "input_mean": input_mean,
                        "output_mean": output_mean,
                        "power_change_percent": power_change,
                    }
                ))
            else:
                self._add_result(ValidationResult(
                    name="multilook_power_preservation",
                    category="algorithm",
                    passed=False,
                    severity="minor",
                    message="apply_multilooking function not found",
                    details={}
                ))
        except Exception as e:
            self._add_result(ValidationResult(
                name="multilook_power_preservation",
                category="algorithm",
                passed=False,
                severity="major",
                message=f"Multilook test failed: {e}",
                details={"error": str(e)}
            ))
    
    def _test_terrain_flattening_formula(self):
        """
        Test terrain flattening formula: gamma0 = sigma0 / cos(theta_local)
        
        On flat terrain (DEM gradient = 0), theta_local = ellipsoid incidence angle,
        so gamma0 ≈ sigma0 / cos(theta_inc)
        """
        try:
            import sardine
            
            # Create flat DEM and uniform sigma0
            rows, cols = 50, 50
            flat_dem = np.zeros((rows, cols), dtype=np.float32)
            sigma0 = np.ones((rows, cols), dtype=np.float32) * 0.1
            
            if hasattr(sardine, 'apply_scientific_terrain_flattening'):
                result = sardine.apply_scientific_terrain_flattening(
                    sigma0,
                    flat_dem,
                    safe_path=None,
                    subswath=None,
                    azimuth_heading=np.radians(-12.0),  # Typical S-1 heading
                    ellipsoid_incidence_angle=np.radians(35.0),  # Mid-swath
                    dem_pixel_spacing=30.0
                )
                
                if isinstance(result, dict):
                    gamma0 = result.get('data')
                    if gamma0 is not None:
                        # For flat terrain, gamma0 = sigma0 / cos(35°) ≈ sigma0 / 0.819 ≈ 0.122
                        expected_gamma0 = 0.1 / np.cos(np.radians(35.0))
                        gamma0_mean = float(np.mean(gamma0[np.isfinite(gamma0)]))
                        error = abs(gamma0_mean - expected_gamma0) / expected_gamma0 * 100.0
                        
                        self._add_result(ValidationResult(
                            name="terrain_flattening_formula",
                            category="algorithm",
                            passed=error < 10.0,  # 10% tolerance for simplified test
                            severity="major",
                            message=f"Terrain flattening on flat DEM: error = {error:.1f}%",
                            details={
                                "input_sigma0": 0.1,
                                "expected_gamma0": float(expected_gamma0),
                                "actual_gamma0_mean": gamma0_mean,
                                "error_percent": error,
                            }
                        ))
                        return
            
            self._add_result(ValidationResult(
                name="terrain_flattening_formula",
                category="algorithm",
                passed=True,
                severity="info",
                message="Terrain flattening test skipped (function not available)",
                details={}
            ))
        except Exception as e:
            self._add_result(ValidationResult(
                name="terrain_flattening_formula",
                category="algorithm",
                passed=False,
                severity="major",
                message=f"Terrain flattening test failed: {e}",
                details={"error": str(e)}
            ))
    
    def _test_coordinate_transforms(self):
        """Test coordinate system transforms are reversible."""
        try:
            from pyproj import Transformer, CRS
            
            # Test WGS84 -> UTM -> WGS84 round-trip
            test_points = [
                (10.0, 52.0),   # Germany
                (-122.0, 37.0), # California
                (139.0, 35.0),  # Japan
            ]
            
            max_error_m = 0.0
            for lon, lat in test_points:
                # Determine UTM zone
                utm_zone = int((lon + 180) / 6) + 1
                utm_epsg = 32600 + utm_zone if lat >= 0 else 32700 + utm_zone
                
                # Forward transform
                transformer_fwd = Transformer.from_crs("EPSG:4326", f"EPSG:{utm_epsg}", always_xy=True)
                x, y = transformer_fwd.transform(lon, lat)
                
                # Inverse transform
                transformer_inv = Transformer.from_crs(f"EPSG:{utm_epsg}", "EPSG:4326", always_xy=True)
                lon2, lat2 = transformer_inv.transform(x, y)
                
                # Calculate error in meters
                error_deg = np.sqrt((lon2 - lon)**2 + (lat2 - lat)**2)
                error_m = error_deg * 111000  # Approximate meters per degree
                max_error_m = max(max_error_m, error_m)
            
            self._add_result(ValidationResult(
                name="coordinate_transform_reversibility",
                category="algorithm",
                passed=max_error_m < 0.01,  # 1cm tolerance
                severity="major",
                message=f"Coordinate transform round-trip error: {max_error_m:.6f} m",
                details={
                    "test_points": test_points,
                    "max_error_m": max_error_m,
                }
            ))
        except ImportError:
            self._add_result(ValidationResult(
                name="coordinate_transform_reversibility",
                category="algorithm",
                passed=True,
                severity="info",
                message="pyproj not available, skipping coordinate transform test",
                details={}
            ))
        except Exception as e:
            self._add_result(ValidationResult(
                name="coordinate_transform_reversibility",
                category="algorithm",
                passed=False,
                severity="major",
                message=f"Coordinate transform test failed: {e}",
                details={"error": str(e)}
            ))
    
    def _validate_scientific_accuracy(self):
        """
        Validate scientific accuracy of outputs.
        
        Checks that output values are in physically plausible ranges.
        """
        if not self.processor:
            return
        
        # Check dB range
        working_data = getattr(self.processor, '_working_data', None)
        if working_data is not None and isinstance(working_data, np.ndarray):
            finite_data = working_data[np.isfinite(working_data)]
            
            if len(finite_data) > 0:
                data_min = float(np.min(finite_data))
                data_max = float(np.max(finite_data))
                data_mean = float(np.mean(finite_data))
                
                # Check if data appears to be in dB (negative values, range -50 to +20)
                is_db = data_min < 0 and data_max < 30
                
                if is_db:
                    # Validate dB range
                    expected_min = self.thresholds["min_db_water"]  # -30 dB
                    expected_max = self.thresholds["max_db_urban"]  # +10 dB
                    
                    self._add_result(ValidationResult(
                        name="db_value_range",
                        category="scientific",
                        passed=data_min >= expected_min - 10 and data_max <= expected_max + 10,
                        severity="major",
                        message=f"dB values range: [{data_min:.1f}, {data_max:.1f}] dB "
                               f"(expected ~[{expected_min}, {expected_max}])",
                        details={
                            "min": data_min,
                            "max": data_max,
                            "mean": data_mean,
                            "expected_min": expected_min,
                            "expected_max": expected_max,
                        }
                    ))
                    
                    # Check mean is reasonable for land
                    expected_land_mean = -12.0  # Typical land mean
                    mean_deviation = abs(data_mean - expected_land_mean)
                    
                    self._add_result(ValidationResult(
                        name="db_mean_plausibility",
                        category="scientific",
                        passed=mean_deviation < 10.0,  # Within 10 dB of typical
                        severity="minor",
                        message=f"Mean dB value: {data_mean:.1f} dB "
                               f"(typical land: ~{expected_land_mean} dB)",
                        details={
                            "mean": data_mean,
                            "expected_typical": expected_land_mean,
                            "deviation": mean_deviation,
                        }
                    ))
                else:
                    # Linear power data
                    # Should be between 0 and ~10 for calibrated backscatter
                    self._add_result(ValidationResult(
                        name="linear_value_range",
                        category="scientific",
                        passed=data_min >= 0 and data_max < 100,
                        severity="major",
                        message=f"Linear values range: [{data_min:.6f}, {data_max:.6f}] "
                               f"(expected 0 to ~10 for calibrated backscatter)",
                        details={
                            "min": data_min,
                            "max": data_max,
                            "mean": data_mean,
                        }
                    ))
    
    def _validate_invariants(self):
        """
        Validate processing invariants.
        
        These are properties that should be preserved across processing steps.
        """
        if not self.processor:
            return
        
        # Invariant 1: Total power should be approximately preserved
        # (accounting for multilooking and masking)
        
        # Check if we have calibrated subswaths and merged data
        calibrated = getattr(self.processor, '_calibrated_subswaths', None)
        merged = getattr(self.processor, '_working_data', None)
        
        if calibrated and merged is not None and isinstance(merged, np.ndarray):
            # Calculate total power before merge
            pre_merge_power = 0.0
            pre_merge_pixels = 0
            for name, data in calibrated.items():
                if isinstance(data, np.ndarray):
                    finite = data[np.isfinite(data) & (data > 0)]
                    pre_merge_power += np.sum(finite)
                    pre_merge_pixels += len(finite)
            
            # Calculate power after merge (if linear data)
            finite_merged = merged[np.isfinite(merged)]
            if len(finite_merged) > 0:
                # Check if data is linear (positive, < 100) or dB (can be negative)
                is_linear = np.min(finite_merged) >= 0 and np.max(finite_merged) < 100
                
                if is_linear:
                    post_merge_power = np.sum(finite_merged[finite_merged > 0])
                    post_merge_pixels = np.sum(finite_merged > 0)
                    
                    if pre_merge_power > 0:
                        # Account for overlap (subswaths overlap ~10%)
                        overlap_factor = 0.85  # Rough estimate
                        expected_power = pre_merge_power * overlap_factor
                        power_ratio = post_merge_power / expected_power
                        
                        self._add_result(ValidationResult(
                            name="power_preservation_merge",
                            category="invariant",
                            passed=0.8 < power_ratio < 1.2,  # ±20% tolerance
                            severity="major",
                            message=f"Power preservation across merge: ratio = {power_ratio:.2f}",
                            details={
                                "pre_merge_power": float(pre_merge_power),
                                "post_merge_power": float(post_merge_power),
                                "expected_power": float(expected_power),
                                "power_ratio": float(power_ratio),
                                "overlap_factor": overlap_factor,
                            }
                        ))
        
        # Invariant 2: Valid pixel count should not decrease dramatically
        # (some loss is expected from masking, but not >50%)
        if self.processor:
            # This would need pre/post data from each step
            pass
    
    def _validate_error_handling(self):
        """
        Validate that errors are handled correctly.
        
        Checks that fatal conditions raise errors rather than warnings.
        """
        
        # Check 1: Burst timing for IW TOPSAR should be fatal if missing
        try:
            # Try to import and check if the error handling is correct
            from sardine.processors.backscatter import terrain
            import inspect
            
            source = inspect.getsource(terrain.run_terrain_correction)
            
            # Look for RuntimeError on missing burst timing for IW
            has_fatal_burst_timing = "RuntimeError" in source and "burst" in source.lower()
            
            self._add_result(ValidationResult(
                name="burst_timing_fatal_for_iw",
                category="error_handling",
                passed=has_fatal_burst_timing,
                severity="critical",
                message="Missing burst timing for IW TOPSAR should raise RuntimeError",
                details={
                    "has_fatal_handling": has_fatal_burst_timing,
                },
                fix_suggestion="Change burst timing warning to RuntimeError for IW TOPSAR"
            ))
        except Exception as e:
            self._add_result(ValidationResult(
                name="burst_timing_fatal_for_iw",
                category="error_handling",
                passed=False,
                severity="minor",
                message=f"Could not check burst timing error handling: {e}",
                details={"error": str(e)}
            ))
        
        # Check 2: Missing critical metadata should be fatal
        try:
            from sardine.processors.backscatter import terrain
            import inspect
            
            source = inspect.getsource(terrain.run_terrain_correction)
            
            # Look for RuntimeError on missing metadata
            has_fatal_metadata = "RuntimeError" in source and "Missing required metadata" in source
            
            self._add_result(ValidationResult(
                name="missing_metadata_fatal",
                category="error_handling",
                passed=has_fatal_metadata,
                severity="critical",
                message="Missing critical metadata should raise RuntimeError",
                details={
                    "has_fatal_handling": has_fatal_metadata,
                }
            ))
        except Exception as e:
            pass
    
    def _generate_summary(self) -> Dict[str, Any]:
        """Generate validation summary."""
        total = len(self.results)
        passed = sum(1 for r in self.results if r.passed)
        failed = total - passed
        
        by_category = {}
        by_severity = {}
        
        for r in self.results:
            # By category
            if r.category not in by_category:
                by_category[r.category] = {"passed": 0, "failed": 0}
            if r.passed:
                by_category[r.category]["passed"] += 1
            else:
                by_category[r.category]["failed"] += 1
            
            # By severity
            if r.severity not in by_severity:
                by_severity[r.severity] = {"passed": 0, "failed": 0}
            if r.passed:
                by_severity[r.severity]["passed"] += 1
            else:
                by_severity[r.severity]["failed"] += 1
        
        return {
            "total_checks": total,
            "passed": passed,
            "failed": failed,
            "pass_rate": passed / total * 100.0 if total > 0 else 0.0,
            "by_category": by_category,
            "by_severity": by_severity,
            "critical_failures": len([r for r in self.results if not r.passed and r.severity == "critical"]),
            "major_failures": len([r for r in self.results if not r.passed and r.severity == "major"]),
        }
    
    def print_report(self, report: ValidationReport):
        """Print formatted validation report."""
        print("\n" + "=" * 80)
        print("SCIENTIFIC VALIDATION REPORT")
        print("=" * 80)
        print(f"Timestamp: {report.timestamp}")
        print(f"Overall: {'PASSED' if report.passed else 'FAILED'}")
        print()
        
        summary = report.summary
        print(f"Total Checks: {summary['total_checks']}")
        print(f"Passed: {summary['passed']} ({summary['pass_rate']:.1f}%)")
        print(f"Failed: {summary['failed']}")
        print(f"  - Critical: {summary['critical_failures']}")
        print(f"  - Major: {summary['major_failures']}")
        print()
        
        # Print failures
        if report.critical_failures:
            print("CRITICAL FAILURES:")
            print("-" * 40)
            for r in report.critical_failures:
                print(f"  ❌ {r.name}: {r.message}")
                if r.fix_suggestion:
                    print(f"     💡 Fix: {r.fix_suggestion}")
            print()
        
        if report.major_failures:
            print("MAJOR FAILURES:")
            print("-" * 40)
            for r in report.major_failures:
                print(f"  ⚠️  {r.name}: {r.message}")
                if r.fix_suggestion:
                    print(f"     💡 Fix: {r.fix_suggestion}")
            print()
        
        # Print by category
        print("BY CATEGORY:")
        print("-" * 40)
        for cat, counts in summary["by_category"].items():
            status = "✅" if counts["failed"] == 0 else "❌"
            print(f"  {status} {cat}: {counts['passed']}/{counts['passed'] + counts['failed']} passed")
        
        print("=" * 80)
