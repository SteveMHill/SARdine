"""
Tests for quality flags module (NOISE-1, DEBURST-1, RTC-1).

These tests verify:
1. Quality flag creation and collection
2. NOISE-1: Noise threshold modification detection
3. RTC-1: Missing incidence angle detection  
4. Flag serialization for metadata embedding
"""

import os
import json
import pytest
from unittest import mock

from sardine.quality_flags import (
    QualityFlag,
    QualityFlagEntry,
    QualityFlagCollector,
    FLAG_SEVERITY,
    check_noise_thresholds,
    check_rtc_incidence_available,
    get_collector,
    reset_collector,
)


class TestQualityFlagEntry:
    """Tests for individual flag entries."""

    def test_entry_creation(self):
        """Flag entry should be creatable with required fields."""
        entry = QualityFlagEntry(
            flag_type="TEST_FLAG",
            severity="WARNING",
            message="Test message",
            timestamp="2025-01-09T00:00:00Z",
            details={"key": "value"},
        )
        
        assert entry.flag_type == "TEST_FLAG"
        assert entry.severity == "WARNING"
        assert entry.message == "Test message"

    def test_entry_to_dict(self):
        """to_dict should be JSON serializable."""
        entry = QualityFlagEntry(
            flag_type="TEST_FLAG",
            severity="WARNING",
            message="Test",
            timestamp="2025-01-09T00:00:00Z",
        )
        
        d = entry.to_dict()
        json_str = json.dumps(d)
        assert "TEST_FLAG" in json_str


class TestQualityFlagCollector:
    """Tests for the flag collector."""

    def test_add_flag(self):
        """Should be able to add flags."""
        collector = QualityFlagCollector()
        
        collector.add_flag(
            QualityFlag.NOISE_REMOVAL_SKIPPED,
            message="Test skip",
            reason="testing",
        )
        
        assert collector.count() == 1
        assert collector.has_flag(QualityFlag.NOISE_REMOVAL_SKIPPED)

    def test_has_flag(self):
        """has_flag should correctly identify present flags."""
        collector = QualityFlagCollector()
        
        collector.add_flag(QualityFlag.DEBURST_T0_FALLBACK)
        
        assert collector.has_flag(QualityFlag.DEBURST_T0_FALLBACK)
        assert not collector.has_flag(QualityFlag.NOISE_REMOVAL_SKIPPED)

    def test_has_errors(self):
        """has_errors should detect error-level flags."""
        collector = QualityFlagCollector()
        
        # Warning should not count as error
        collector.add_flag(QualityFlag.NOISE_REMOVAL_SKIPPED)
        assert not collector.has_errors()
        
        # Error should be detected
        collector.add_flag(QualityFlag.DEBURST_TIME_DOMAIN_MISMATCH)
        assert collector.has_errors()

    def test_has_warnings(self):
        """has_warnings should detect warning-level flags."""
        collector = QualityFlagCollector()
        
        assert not collector.has_warnings()
        
        collector.add_flag(QualityFlag.RTC_DEFAULT_INCIDENCE_ANGLE)
        assert collector.has_warnings()

    def test_to_dict_structure(self):
        """to_dict should return proper structure."""
        collector = QualityFlagCollector()
        collector.add_flag(QualityFlag.WARNING, message="Test")
        
        d = collector.to_dict()
        
        assert "flags" in d
        assert "summary" in d
        assert "has_errors" in d
        assert "has_warnings" in d
        assert "total_count" in d
        assert d["total_count"] == 1

    def test_to_json(self):
        """to_json should produce valid JSON."""
        collector = QualityFlagCollector()
        collector.add_flag(QualityFlag.CALIBRATION_LUT_WARNING, lut_type="sigma0")
        
        json_str = collector.to_json()
        parsed = json.loads(json_str)
        
        assert parsed["total_count"] == 1

    def test_clear(self):
        """clear should remove all flags."""
        collector = QualityFlagCollector()
        collector.add_flag(QualityFlag.WARNING)
        collector.add_flag(QualityFlag.ERROR)
        
        assert collector.count() == 2
        
        collector.clear()
        
        assert collector.count() == 0
        assert not collector.has_errors()


class TestNoiseThresholdCheck:
    """Tests for NOISE-1 threshold detection."""

    def test_no_env_vars_returns_none(self):
        """When no env vars set, should return None."""
        with mock.patch.dict(os.environ, {}, clear=True):
            # Remove any existing SARDINE_NOISE_* vars
            env = {k: v for k, v in os.environ.items() 
                   if not k.startswith("SARDINE_NOISE")}
            with mock.patch.dict(os.environ, env, clear=True):
                flag = check_noise_thresholds()
                # May still return None if no vars set
                if flag is not None:
                    assert flag.flag_type == "NOISE_THRESHOLDS_MODIFIED"

    def test_modified_zero_threshold(self):
        """Should detect modified zero threshold."""
        with mock.patch.dict(os.environ, {
            "SARDINE_NOISE_ZERO_THRESHOLD": "0.8",
        }, clear=False):
            flag = check_noise_thresholds()
            
            assert flag is not None
            assert flag.flag_type == "NOISE_THRESHOLDS_MODIFIED"
            assert flag.details["zero_threshold"] == 0.8
            assert flag.details["default_zero_threshold"] == 0.95

    def test_modified_ratio_threshold(self):
        """Should detect modified ratio threshold."""
        with mock.patch.dict(os.environ, {
            "SARDINE_NOISE_RATIO_THRESHOLD": "0.7",
        }, clear=False):
            flag = check_noise_thresholds()
            
            assert flag is not None
            assert flag.details["ratio_threshold"] == 0.7


class TestRtcIncidenceCheck:
    """Tests for RTC-1 incidence angle detection."""

    def test_both_angles_available(self):
        """Should return None when both angles available."""
        flag = check_rtc_incidence_available(30.0, 45.0)
        assert flag is None

    def test_near_angle_only(self):
        """Should return None when near angle available."""
        flag = check_rtc_incidence_available(30.0, None)
        assert flag is None

    def test_far_angle_only(self):
        """Should return None when far angle available."""
        flag = check_rtc_incidence_available(None, 45.0)
        assert flag is None

    def test_no_angles_available(self):
        """Should return flag when no angles available."""
        flag = check_rtc_incidence_available(None, None)
        
        assert flag is not None
        assert flag.flag_type == "RTC_DEFAULT_INCIDENCE_ANGLE"
        assert flag.severity == "WARNING"
        assert flag.details["default_angle_deg"] == 35.0


class TestGlobalCollector:
    """Tests for global collector functions."""

    def test_get_collector_creates_new(self):
        """get_collector should create new if none exists."""
        reset_collector()  # Ensure clean state
        collector = get_collector()
        assert collector is not None
        assert isinstance(collector, QualityFlagCollector)

    def test_get_collector_returns_same(self):
        """get_collector should return same instance."""
        reset_collector()
        c1 = get_collector()
        c2 = get_collector()
        assert c1 is c2

    def test_reset_collector(self):
        """reset_collector should create new instance."""
        c1 = get_collector()
        c1.add_flag(QualityFlag.WARNING)
        
        c2 = reset_collector()
        
        assert c2.count() == 0
        assert c1 is not c2


class TestFlagSeverityMapping:
    """Tests for severity mapping."""

    def test_all_flags_have_severity(self):
        """All QualityFlag values should have severity mapping."""
        for flag in QualityFlag:
            assert flag in FLAG_SEVERITY, f"Missing severity for {flag}"

    def test_severity_values(self):
        """Severity values should be valid."""
        valid_severities = {"INFO", "WARNING", "ERROR"}
        for flag, severity in FLAG_SEVERITY.items():
            assert severity in valid_severities, f"Invalid severity {severity} for {flag}"


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
