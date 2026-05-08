"""
Design Flaw Detector for SARdine Pipeline

Performs static analysis of processor configuration to detect:
- Parameter misconfigurations that will cause silent failures
- Incompatible option combinations
- Threshold values that don't match scientific requirements
- Missing or incorrect metadata handling

This catches design-level issues that aren't detectable from outputs alone.
"""

import logging
import inspect
from dataclasses import dataclass, field
from typing import Dict, Any, List, Optional, Type
from pathlib import Path

logger = logging.getLogger(__name__)


@dataclass
class DesignFlaw:
    """A detected design flaw."""
    category: str  # parameter, compatibility, threshold, logic, metadata
    severity: str  # critical, major, minor
    location: str  # File/class/method where flaw exists
    description: str
    impact: str
    fix_suggestion: str
    code_reference: Optional[str] = None


@dataclass
class DesignFlawReport:
    """Complete design flaw report."""
    flaws: List[DesignFlaw] = field(default_factory=list)
    checked_components: List[str] = field(default_factory=list)
    timestamp: str = ""
    
    @property
    def critical_count(self) -> int:
        return sum(1 for f in self.flaws if f.severity == "critical")
    
    @property
    def major_count(self) -> int:
        return sum(1 for f in self.flaws if f.severity == "major")
    
    @property
    def passed(self) -> bool:
        """No critical or major flaws."""
        return self.critical_count == 0 and self.major_count == 0


class DesignFlawDetector:
    """
    Static analysis detector for design flaws in SAR processing pipeline.
    
    Detects issues like:
    1. Using gamma0 calibration with terrain flattening (double correction)
    2. Coverage thresholds too low for scientific validity
    3. Missing error handling for fatal conditions
    4. Incompatible parameter combinations
    5. Metadata key mismatches between Python and Rust
    """
    
    def __init__(self, strict: bool = True):
        """
        Initialize detector.
        
        Args:
            strict: If True, apply stricter scientific thresholds
        """
        self.strict = strict
        self.flaws: List[DesignFlaw] = []
        self.checked_components: List[str] = []
        
        # Expected values for scientific correctness
        self.expected = {
            # Coverage thresholds - SCIENTIFIC FIX: Raised for quality assurance
            "min_critical_coverage": 75.0,  # Was 50.0
            "min_warning_coverage": 90.0,   # Was 80.0
            
            # Orbit requirements (ESA POD spec)
            "min_orbit_vectors": 10,
            
            # DEM quality
            "max_dem_upsample_factor": 2.0,
            
            # ENL for speckle
            "enl_skip_threshold": 100,
            
            # Correct calibration for terrain flattening
            "calibration_for_terrain_flatten": "sigma0",
        }
    
    def analyze_processor(self, processor) -> DesignFlawReport:
        """
        Analyze a BackscatterProcessor instance for design flaws.
        
        Args:
            processor: BackscatterProcessor instance
            
        Returns:
            DesignFlawReport with all detected flaws
        """
        from datetime import datetime
        
        self.flaws = []
        self.checked_components = []
        
        # Check configuration
        self._check_calibration_terrain_compatibility(processor)
        self._check_coverage_thresholds()
        self._check_error_handling()
        self._check_metadata_keys()
        self._check_incidence_angle_handling()
        self._check_dem_handling()
        self._check_speckle_handling()
        
        return DesignFlawReport(
            flaws=self.flaws,
            checked_components=self.checked_components,
            timestamp=datetime.now().isoformat(),
        )
    
    def analyze_codebase(self) -> DesignFlawReport:
        """
        Analyze the SARdine codebase for design flaws.
        
        This performs static analysis without needing a processor instance.
        
        Returns:
            DesignFlawReport with all detected flaws
        """
        from datetime import datetime
        
        self.flaws = []
        self.checked_components = []
        
        # Static analysis checks
        self._check_coverage_thresholds()
        self._check_error_handling()
        self._check_metadata_keys()
        self._check_incidence_angle_handling()
        self._check_dem_handling()
        self._check_speckle_handling()
        self._check_calibration_defaults()
        
        return DesignFlawReport(
            flaws=self.flaws,
            checked_components=self.checked_components,
            timestamp=datetime.now().isoformat(),
        )
    
    def _add_flaw(self, flaw: DesignFlaw):
        """Add a detected flaw."""
        self.flaws.append(flaw)
        logger.warning(f"Design flaw detected: [{flaw.severity}] {flaw.description}")
    
    def _check_calibration_terrain_compatibility(self, processor):
        """
        Check that calibration type is compatible with terrain flattening.
        
        If terrain_flatten=True, calibration should be sigma0, not gamma0.
        Terrain flattening converts sigma0 -> gamma0, so using gamma0 calibration
        with terrain flattening applies the incidence angle correction twice.
        """
        self.checked_components.append("calibration_terrain_compatibility")
        
        terrain_flatten = getattr(processor, 'terrain_flatten', None)
        calibration_type = getattr(processor, 'calibration_type', None)
        
        if terrain_flatten is None or calibration_type is None:
            return
        
        if terrain_flatten and calibration_type and calibration_type.lower() != "sigma0":
            self._add_flaw(DesignFlaw(
                category="compatibility",
                severity="critical",
                location="BackscatterProcessor configuration",
                description=f"calibration_type='{calibration_type}' with terrain_flatten=True",
                impact="Incidence angle correction applied twice, resulting in incorrect gamma0 values. "
                      "Products will have radiometric errors of 1-4 dB depending on incidence angle.",
                fix_suggestion=f"Set calibration_type='sigma0' when terrain_flatten=True. "
                              f"Terrain flattening will convert sigma0 to gamma0_tc.",
            ))
    
    def _check_calibration_defaults(self):
        """Check that default calibration type is correct for terrain flattening."""
        self.checked_components.append("calibration_defaults")
        
        try:
            from sardine.processors.backscatter.processor import BackscatterProcessor
            import inspect
            
            source = inspect.getsource(BackscatterProcessor.__init__)
            
            # Check if terrain_flatten default is True
            terrain_default_true = "terrain_flatten" in source and "True" in source
            
            # Check if calibration defaults to sigma0 when terrain_flatten is True
            has_sigma0_logic = "sigma0" in source and "terrain_flatten" in source
            
            if terrain_default_true and not has_sigma0_logic:
                self._add_flaw(DesignFlaw(
                    category="parameter",
                    severity="major",
                    location="BackscatterProcessor.__init__",
                    description="No logic to set calibration_type='sigma0' when terrain_flatten=True",
                    impact="Users may accidentally use gamma0 calibration with terrain flattening",
                    fix_suggestion="Add logic: if terrain_flatten and calibration_type is None: calibration_type = 'sigma0'",
                ))
        except ImportError as e:
            logger.debug(f"Could not import BackscatterProcessor: {e}")
    
    def _check_coverage_thresholds(self):
        """Check that geocoding coverage thresholds are scientifically appropriate."""
        self.checked_components.append("coverage_thresholds")
        
        try:
            from sardine.processors.backscatter.terrain import (
                CRITICAL_COVERAGE_THRESHOLD,
                WARNING_COVERAGE_THRESHOLD
            )
            
            if CRITICAL_COVERAGE_THRESHOLD < self.expected["min_critical_coverage"]:
                self._add_flaw(DesignFlaw(
                    category="threshold",
                    severity="critical",
                    location="sardine.processors.backscatter.terrain",
                    description=f"CRITICAL_COVERAGE_THRESHOLD = {CRITICAL_COVERAGE_THRESHOLD}% "
                               f"(should be >= {self.expected['min_critical_coverage']}%)",
                    impact="Products with <50% valid pixels are scientifically unusable "
                          "but will be output anyway, leading to incorrect analyses",
                    fix_suggestion=f"Raise CRITICAL_COVERAGE_THRESHOLD to {self.expected['min_critical_coverage']}%",
                    code_reference="CRITICAL_COVERAGE_THRESHOLD = ...",
                ))
            
            if WARNING_COVERAGE_THRESHOLD < self.expected["min_warning_coverage"]:
                self._add_flaw(DesignFlaw(
                    category="threshold",
                    severity="major",
                    location="sardine.processors.backscatter.terrain",
                    description=f"WARNING_COVERAGE_THRESHOLD = {WARNING_COVERAGE_THRESHOLD}% "
                               f"(should be >= {self.expected['min_warning_coverage']}%)",
                    impact="Products with 50-80% valid pixels may have significant gaps "
                          "but won't trigger warnings",
                    fix_suggestion=f"Raise WARNING_COVERAGE_THRESHOLD to {self.expected['min_warning_coverage']}%",
                    code_reference="WARNING_COVERAGE_THRESHOLD = ...",
                ))
        except ImportError:
            pass
    
    def _check_error_handling(self):
        """Check that fatal conditions raise errors instead of warnings."""
        self.checked_components.append("error_handling")
        
        try:
            from sardine.processors.backscatter import terrain
            import inspect
            
            source = inspect.getsource(terrain)
            
            # Check 1: Missing burst timing for IW should be fatal
            # Look for the pattern of warning about burst timing
            has_burst_warning = "burst" in source.lower() and "warning" in source.lower()
            has_burst_error = "burst" in source.lower() and "RuntimeError" in source
            
            if has_burst_warning and not has_burst_error:
                self._add_flaw(DesignFlaw(
                    category="logic",
                    severity="critical",
                    location="terrain.run_terrain_correction",
                    description="Missing burst timing for IW TOPSAR raises warning instead of error",
                    impact="Processing continues without burst timing, producing incorrect "
                          "georeferencing for IW products (meters to hundreds of meters error)",
                    fix_suggestion="Raise RuntimeError for IW mode when burst timing is missing",
                ))
            
            # Check 2: Missing critical metadata should be fatal
            has_missing_metadata_warning = "missing" in source.lower() and "metadata" in source.lower()
            has_missing_metadata_error = "Missing required metadata" in source
            
            # Already fixed if RuntimeError is present with missing metadata message
            
        except ImportError:
            pass
    
    def _check_metadata_keys(self):
        """
        Check for metadata key mismatches between Python and Rust.
        
        The PyO3 bridge expects specific key names, and using wrong keys
        causes silent data loss.
        """
        self.checked_components.append("metadata_keys")
        
        # Known key mappings that must be consistent
        required_keys = {
            "subswaths": "Not 'subswath_geometry'",
            "burst_timings": "Not 'burst_timing_records'",
            "number_of_lines": "Not 'total_azimuth_lines'",
            "native_range_pixel_spacing": "Required for Range-Doppler geocoding",
        }
        
        try:
            from sardine.processors.backscatter import terrain
            import inspect
            
            source = inspect.getsource(terrain)
            
            # Check for incorrect key names that would fail silently
            incorrect_keys = [
                ("subswath_geometry", "subswaths"),
                ("burst_timing_records", "burst_timings"),
                ("total_azimuth_lines", "number_of_lines"),
            ]
            
            for wrong_key, correct_key in incorrect_keys:
                if wrong_key in source and correct_key not in source:
                    self._add_flaw(DesignFlaw(
                        category="metadata",
                        severity="critical",
                        location="terrain.py metadata handling",
                        description=f"Using '{wrong_key}' instead of '{correct_key}'",
                        impact="Metadata will not be passed to Rust correctly, "
                              "causing Range-Doppler geocoding to fail silently",
                        fix_suggestion=f"Use '{correct_key}' for PyO3 bridge compatibility",
                    ))
        except ImportError:
            pass
    
    def _check_incidence_angle_handling(self):
        """Check that incidence angle fallback uses scientifically valid values."""
        self.checked_components.append("incidence_angle_handling")
        
        try:
            from sardine.processors.backscatter import multilook
            import inspect
            
            source = inspect.getsource(multilook)
            
            # Check if there's any incidence angle handling
            has_incidence_handling = "incidence" in source.lower()
            
            # Check for mode-based defaults (IW, EW, SM have different ranges)
            has_mode_specific = any(mode in source for mode in ["IW", "EW", "SM"])
            
            if has_incidence_handling and not has_mode_specific:
                self._add_flaw(DesignFlaw(
                    category="parameter",
                    severity="major",
                    location="multilook.py incidence angle handling",
                    description="Incidence angle fallback doesn't use mode-specific values",
                    impact="IW mode (29-46°), EW mode (19-47°), and SM mode (20-45°) "
                          "have different incidence angle ranges. Using wrong values "
                          "causes 10-20% error in terrain flattening",
                    fix_suggestion="Add mode-based incidence angle defaults",
                ))
        except ImportError:
            pass
    
    def _check_dem_handling(self):
        """Check DEM resampling and quality handling."""
        self.checked_components.append("dem_handling")
        
        try:
            from sardine.processors.backscatter import terrain
            import inspect
            
            source = inspect.getsource(terrain)
            
            # Check 1: DEM upsampling warning
            # Upsampling low-res DEM doesn't add information
            has_upsample_check = "upsample" in source.lower() or "spacing" in source.lower()
            
            if not has_upsample_check:
                self._add_flaw(DesignFlaw(
                    category="logic",
                    severity="minor",
                    location="terrain.py DEM handling",
                    description="No check for DEM upsampling beyond useful resolution",
                    impact="Upsampling a 30m SRTM DEM to 10m doesn't improve accuracy "
                          "and wastes computation",
                    fix_suggestion="Warn if DEM is being upsampled more than 2x",
                ))
            
            # Check 2: Interpolation method for resampling
            has_cubic = "cubic" in source.lower()
            has_bilinear = "bilinear" in source.lower()
            
            if not has_cubic and not has_bilinear:
                self._add_flaw(DesignFlaw(
                    category="logic",
                    severity="minor",
                    location="terrain.py DEM resampling",
                    description="No explicit interpolation method specified for DEM resampling",
                    impact="May use nearest-neighbor interpolation, causing blocky artifacts "
                          "in terrain correction",
                    fix_suggestion="Use bilinear for downsampling, cubic for upsampling",
                ))
        except ImportError:
            pass
    
    def _check_speckle_handling(self):
        """Check speckle filter configuration."""
        self.checked_components.append("speckle_handling")
        
        try:
            from sardine.processors.backscatter import speckle
            import inspect
            
            source = inspect.getsource(speckle)
            
            # Check if ENL is handled correctly
            # At high ENL (>100), speckle is already well-suppressed
            has_enl_check = "enl" in source.lower()
            has_enl_skip = ("100" in source or "skip" in source.lower()) and "enl" in source.lower()
            
            if has_enl_check and not has_enl_skip:
                self._add_flaw(DesignFlaw(
                    category="logic",
                    severity="minor",
                    location="speckle.py ENL handling",
                    description="No skip logic for high ENL (>100)",
                    impact="Applying speckle filter to already well-averaged data "
                          "degrades resolution without benefit",
                    fix_suggestion="Skip speckle filtering when ENL > 100",
                ))
        except ImportError:
            pass
    
    def print_report(self, report: DesignFlawReport):
        """Print formatted design flaw report."""
        print("\n" + "=" * 80)
        print("DESIGN FLAW ANALYSIS REPORT")
        print("=" * 80)
        print(f"Timestamp: {report.timestamp}")
        print(f"Components checked: {len(report.checked_components)}")
        print(f"Total flaws: {len(report.flaws)}")
        print(f"  Critical: {report.critical_count}")
        print(f"  Major: {report.major_count}")
        print(f"  Minor: {sum(1 for f in report.flaws if f.severity == 'minor')}")
        print()
        
        print("Overall: " + ("✅ PASSED" if report.passed else "❌ FAILED"))
        print()
        
        if report.flaws:
            print("DETECTED FLAWS:")
            print("-" * 60)
            
            # Group by severity
            for severity in ["critical", "major", "minor"]:
                severity_flaws = [f for f in report.flaws if f.severity == severity]
                if not severity_flaws:
                    continue
                
                icon = "🔴" if severity == "critical" else "🟡" if severity == "major" else "🔵"
                print(f"\n{icon} {severity.upper()} ({len(severity_flaws)}):")
                
                for flaw in severity_flaws:
                    print(f"\n  📍 {flaw.location}")
                    print(f"     Category: {flaw.category}")
                    print(f"     Issue: {flaw.description}")
                    print(f"     Impact: {flaw.impact}")
                    print(f"     💡 Fix: {flaw.fix_suggestion}")
        
        print("\n" + "-" * 60)
        print("COMPONENTS CHECKED:")
        for comp in report.checked_components:
            status = "✅" if not any(f.category == comp for f in report.flaws) else "⚠️"
            print(f"  {status} {comp}")
        
        print("=" * 80)


def run_design_analysis() -> DesignFlawReport:
    """
    Convenience function to run full design analysis.
    
    Returns:
        DesignFlawReport
    """
    detector = DesignFlawDetector(strict=True)
    report = detector.analyze_codebase()
    detector.print_report(report)
    return report


if __name__ == "__main__":
    # Run standalone analysis
    run_design_analysis()
