"""
Validation utilities for SARdine pipelines.

This module provides comprehensive validation capabilities:

1. ValidationGates - Basic validation gates for pipeline outputs
2. ScientificValidator - Scientific validation of algorithm correctness
3. ReferenceComparator - Compare against ASF/SNAP reference products
4. InvariantChecker - Check processing invariants (power preservation, etc.)
5. DesignFlawDetector - Static analysis to detect design flaws
6. AlgorithmTester - Unit tests for individual algorithms

Usage:
    from sardine.validation import (
        ScientificValidator,
        ReferenceComparator,
        InvariantChecker,
        DesignFlawDetector,
        AlgorithmTester,
    )
    
    # Run scientific validation
    validator = ScientificValidator(processor)
    report = validator.run_full_validation()
    validator.print_report(report)
    
    # Compare with ASF RTC product
    comparator = ReferenceComparator()
    report = comparator.compare_asf_rtc(sardine_output, asf_rtc_path)
    
    # Check invariants during processing
    checker = InvariantChecker()
    checker.record_step("before_multilook", data_before)
    checker.record_step("after_multilook", data_after)
    checker.check_power_preservation("before_multilook", "after_multilook")
    
    # Detect design flaws
    detector = DesignFlawDetector()
    report = detector.analyze_codebase()
    
    # Run algorithm unit tests
    tester = AlgorithmTester()
    report = tester.run_all_tests()
"""

from .gates import ValidationGates, ValidationRecord, ValidationError

# Scientific validation (new)
from .scientific_validator import (
    ScientificValidator,
    ValidationResult,
    ValidationReport,
)

# Reference product comparison (new)
from .reference_comparator import (
    ReferenceComparator,
    ComparisonResult,
    ReferenceComparisonReport,
)

# Invariant checking (new)
from .invariant_checker import (
    InvariantChecker,
    InvariantCheckResult,
    InvariantReport,
)

# Design flaw detection (new)
from .design_flaw_detector import (
    DesignFlawDetector,
    DesignFlaw,
    DesignFlawReport,
    run_design_analysis,
)

# Algorithm testing (new)
from .algorithm_tester import (
    AlgorithmTester,
    AlgorithmTestResult,
    AlgorithmTestReport,
    run_algorithm_tests,
)

__all__ = [
    # Original exports
    "ValidationGates",
    "ValidationRecord", 
    "ValidationError",
    
    # Scientific validation
    "ScientificValidator",
    "ValidationResult",
    "ValidationReport",
    
    # Reference comparison
    "ReferenceComparator",
    "ComparisonResult",
    "ReferenceComparisonReport",
    
    # Invariant checking
    "InvariantChecker",
    "InvariantCheckResult",
    "InvariantReport",
    
    # Design flaw detection
    "DesignFlawDetector",
    "DesignFlaw",
    "DesignFlawReport",
    "run_design_analysis",
    
    # Algorithm testing
    "AlgorithmTester",
    "AlgorithmTestResult",
    "AlgorithmTestReport",
    "run_algorithm_tests",
]
