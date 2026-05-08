"""
Quality Flags for SARdine Processing
====================================

This module implements the quality flag system required by the scientific audit
(NOISE-1, DEBURST-1, RTC-1). Quality flags are collected during processing and
embedded in output metadata.

Usage:
    from sardine.quality_flags import QualityFlagCollector, QualityFlag

    # Create collector at start of processing
    collector = QualityFlagCollector()
    
    # Add flags during processing
    if noise_skipped:
        collector.add_flag(QualityFlag.NOISE_REMOVAL_SKIPPED, 
                          reason="Zero ratio exceeded threshold",
                          zero_ratio=0.96, threshold=0.95)
    
    # Get flags for metadata
    flags_dict = collector.to_dict()
"""

import os
import logging
from dataclasses import dataclass, field, asdict
from datetime import datetime, timezone
from enum import Enum, auto
from typing import Any, Dict, List, Optional, Set
import json

logger = logging.getLogger(__name__)


class QualityFlag(Enum):
    """Quality flag types for scientific processing validation."""
    
    # NOISE-1: Noise removal conditions
    NOISE_REMOVAL_SKIPPED = auto()
    NOISE_THRESHOLDS_MODIFIED = auto()
    
    # DEBURST-1: Deburst fallback conditions  
    DEBURST_T0_FALLBACK = auto()
    DEBURST_TIME_DOMAIN_MISMATCH = auto()
    
    # RTC-1: Terrain correction conditions
    RTC_DEFAULT_INCIDENCE_ANGLE = auto()
    RTC_INCIDENCE_PROFILE_MISSING = auto()
    
    # Calibration conditions
    CALIBRATION_LUT_WARNING = auto()
    
    # General
    WARNING = auto()
    ERROR = auto()


# Severity mapping
FLAG_SEVERITY: Dict[QualityFlag, str] = {
    QualityFlag.NOISE_REMOVAL_SKIPPED: "WARNING",
    QualityFlag.NOISE_THRESHOLDS_MODIFIED: "INFO",
    QualityFlag.DEBURST_T0_FALLBACK: "WARNING",
    QualityFlag.DEBURST_TIME_DOMAIN_MISMATCH: "ERROR",
    QualityFlag.RTC_DEFAULT_INCIDENCE_ANGLE: "WARNING",
    QualityFlag.RTC_INCIDENCE_PROFILE_MISSING: "WARNING",
    QualityFlag.CALIBRATION_LUT_WARNING: "WARNING",
    QualityFlag.WARNING: "WARNING",
    QualityFlag.ERROR: "ERROR",
}


@dataclass
class QualityFlagEntry:
    """Individual quality flag entry with context."""
    flag_type: str
    severity: str
    message: str
    timestamp: str
    details: Dict[str, Any] = field(default_factory=dict)
    
    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)


class QualityFlagCollector:
    """
    Collects quality flags during processing.
    
    Thread-safe collection of quality flags that can be embedded in output metadata.
    """
    
    def __init__(self):
        self._flags: List[QualityFlagEntry] = []
        self._flag_types_seen: Set[str] = set()
    
    def add_flag(
        self,
        flag: QualityFlag,
        message: str = "",
        **details: Any
    ) -> None:
        """
        Add a quality flag.
        
        Args:
            flag: The QualityFlag type
            message: Human-readable description
            **details: Additional context (e.g., threshold values, ratios)
        """
        entry = QualityFlagEntry(
            flag_type=flag.name,
            severity=FLAG_SEVERITY.get(flag, "WARNING"),
            message=message,
            timestamp=datetime.now(timezone.utc).isoformat(),
            details=details,
        )
        
        self._flags.append(entry)
        self._flag_types_seen.add(flag.name)
        
        # Also log the flag
        log_msg = f"🚩 Quality flag: {flag.name}"
        if message:
            log_msg += f" - {message}"
        if details:
            log_msg += f" ({details})"
            
        if entry.severity == "ERROR":
            logger.error(log_msg)
        elif entry.severity == "WARNING":
            logger.warning(log_msg)
        else:
            logger.info(log_msg)
    
    def has_flag(self, flag: QualityFlag) -> bool:
        """Check if a specific flag type is present."""
        return flag.name in self._flag_types_seen
    
    def has_errors(self) -> bool:
        """Check if any error-level flags are present."""
        return any(f.severity == "ERROR" for f in self._flags)
    
    def has_warnings(self) -> bool:
        """Check if any warning-level flags are present."""
        return any(f.severity == "WARNING" for f in self._flags)
    
    def count(self) -> int:
        """Get total number of flags."""
        return len(self._flags)
    
    def get_flags(self) -> List[QualityFlagEntry]:
        """Get all flag entries."""
        return list(self._flags)
    
    def get_flags_by_severity(self, severity: str) -> List[QualityFlagEntry]:
        """Get flags filtered by severity."""
        return [f for f in self._flags if f.severity == severity]
    
    def to_dict(self) -> Dict[str, Any]:
        """
        Convert to dictionary for embedding in metadata.
        
        Returns a dict with:
        - flags: List of all flag entries
        - summary: Dict mapping flag types to counts
        - has_errors: Boolean
        - has_warnings: Boolean
        """
        summary: Dict[str, int] = {}
        for f in self._flags:
            summary[f.flag_type] = summary.get(f.flag_type, 0) + 1
        
        return {
            "flags": [f.to_dict() for f in self._flags],
            "summary": summary,
            "has_errors": self.has_errors(),
            "has_warnings": self.has_warnings(),
            "total_count": len(self._flags),
        }
    
    def to_json(self, indent: int = 2) -> str:
        """Convert to JSON string."""
        return json.dumps(self.to_dict(), indent=indent)
    
    def clear(self) -> None:
        """Clear all flags."""
        self._flags.clear()
        self._flag_types_seen.clear()


def check_noise_thresholds() -> Optional[QualityFlagEntry]:
    """
    Check if noise thresholds are modified by environment variables.
    
    NOISE-1: This surfaces when env vars are modifying noise removal behavior.
    
    Returns:
        QualityFlagEntry if thresholds are modified, None otherwise.
    """
    default_zero = 0.95
    default_ratio = 0.9
    
    zero_threshold = os.environ.get("SARDINE_NOISE_ZERO_THRESHOLD")
    ratio_threshold = os.environ.get("SARDINE_NOISE_RATIO_THRESHOLD")
    
    if zero_threshold is not None or ratio_threshold is not None:
        try:
            actual_zero = float(zero_threshold) if zero_threshold else default_zero
            actual_ratio = float(ratio_threshold) if ratio_threshold else default_ratio
        except ValueError:
            actual_zero = default_zero
            actual_ratio = default_ratio
        
        return QualityFlagEntry(
            flag_type=QualityFlag.NOISE_THRESHOLDS_MODIFIED.name,
            severity="INFO",
            message=f"Noise thresholds modified by environment variables",
            timestamp=datetime.now(timezone.utc).isoformat(),
            details={
                "zero_threshold": actual_zero,
                "ratio_threshold": actual_ratio,
                "default_zero_threshold": default_zero,
                "default_ratio_threshold": default_ratio,
                "SARDINE_NOISE_ZERO_THRESHOLD": zero_threshold,
                "SARDINE_NOISE_RATIO_THRESHOLD": ratio_threshold,
            }
        )
    
    return None


def check_rtc_incidence_available(incidence_near: Optional[float], incidence_far: Optional[float]) -> Optional[QualityFlagEntry]:
    """
    Check if RTC incidence angle is available from annotation.
    
    RTC-1: This surfaces when default incidence angle would be used.
    
    Args:
        incidence_near: Near-range incidence angle from annotation
        incidence_far: Far-range incidence angle from annotation
        
    Returns:
        QualityFlagEntry if incidence angles are missing, None otherwise.
    """
    if incidence_near is None and incidence_far is None:
        return QualityFlagEntry(
            flag_type=QualityFlag.RTC_DEFAULT_INCIDENCE_ANGLE.name,
            severity="WARNING",
            message="Using default incidence angle (35°) - annotation values not available",
            timestamp=datetime.now(timezone.utc).isoformat(),
            details={
                "default_angle_deg": 35.0,
                "reason": "Incidence angle not available from annotation metadata",
            }
        )
    
    return None


# Global collector for the current processing run
_current_collector: Optional[QualityFlagCollector] = None


def get_collector() -> QualityFlagCollector:
    """Get or create the current quality flag collector."""
    global _current_collector
    if _current_collector is None:
        _current_collector = QualityFlagCollector()
    return _current_collector


def reset_collector() -> QualityFlagCollector:
    """Reset and return a new quality flag collector."""
    global _current_collector
    _current_collector = QualityFlagCollector()
    return _current_collector


if __name__ == "__main__":
    # Self-test
    collector = QualityFlagCollector()
    
    # Test adding flags
    collector.add_flag(
        QualityFlag.NOISE_REMOVAL_SKIPPED,
        message="Zero ratio exceeded threshold",
        zero_ratio=0.96,
        threshold=0.95
    )
    
    collector.add_flag(
        QualityFlag.DEBURST_T0_FALLBACK,
        message="Missing polynomial t0 for burst 3",
        burst_id=3
    )
    
    # Test checks
    flag = check_noise_thresholds()
    if flag:
        collector._flags.append(flag)
    
    print("=== Quality Flags ===")
    print(collector.to_json())
    
    print(f"\nHas errors: {collector.has_errors()}")
    print(f"Has warnings: {collector.has_warnings()}")
    print(f"Total count: {collector.count()}")
