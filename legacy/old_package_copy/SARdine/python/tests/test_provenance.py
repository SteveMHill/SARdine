"""
Tests for provenance module (PROV-1, PROV-2).

These tests verify:
1. Package version is captured
2. Git commit hash is captured when available
3. All SARDINE_* env vars are captured
4. Credentials are properly redacted
5. Provenance dict is suitable for JSON serialization
"""

import os
import json
import pytest
from unittest import mock

from sardine.provenance import (
    capture_provenance,
    get_package_version,
    get_git_commit_hash,
    capture_sardine_env_vars,
    get_provenance_dict,
    validate_provenance_completeness,
    ProvenanceInfo,
    KNOWN_SARDINE_ENV_VARS,
    REDACTED_VARS,
)


class TestPackageVersion:
    """Tests for package version capture (PROV-2)."""

    def test_version_is_string(self):
        """Version should be a non-empty string."""
        version = get_package_version()
        assert isinstance(version, str)
        assert len(version) > 0

    def test_version_not_error(self):
        """Version should not be an error message."""
        version = get_package_version()
        # Even if unknown, it should be "unknown" not an exception
        assert "error" not in version.lower() or version == "unknown"


class TestGitCommitHash:
    """Tests for git commit hash capture (PROV-2)."""

    def test_git_hash_is_string_or_none(self):
        """Git hash should be a string or None."""
        git_hash = get_git_commit_hash()
        assert git_hash is None or isinstance(git_hash, str)

    def test_git_hash_format_if_present(self):
        """If git hash is present, it should be 12 hex characters."""
        git_hash = get_git_commit_hash()
        if git_hash is not None:
            assert len(git_hash) == 12
            # Check it's hex
            int(git_hash, 16)


class TestEnvironmentVariables:
    """Tests for SARDINE_* env var capture (PROV-1)."""

    def test_captures_sardine_env_vars(self):
        """Should capture SARDINE_* environment variables."""
        with mock.patch.dict(os.environ, {
            "SARDINE_ORBIT_CACHE": "/tmp/orbit",
            "SARDINE_STRICT": "1",
            "OTHER_VAR": "ignored",
        }, clear=False):
            env_vars = capture_sardine_env_vars()
            
            assert "SARDINE_ORBIT_CACHE" in env_vars
            assert env_vars["SARDINE_ORBIT_CACHE"] == "/tmp/orbit"
            assert "SARDINE_STRICT" in env_vars
            assert env_vars["SARDINE_STRICT"] == "1"
            # Non-SARDINE vars should not be captured
            assert "OTHER_VAR" not in env_vars

    def test_credentials_redacted(self):
        """Credentials should be redacted."""
        with mock.patch.dict(os.environ, {
            "SARDINE_ESA_USERNAME": "secret_user",
            "SARDINE_ESA_PASSWORD": "secret_pass",
            "SARDINE_ASF_USERNAME": "asf_user",
            "SARDINE_ASF_PASSWORD": "asf_pass",
        }, clear=False):
            env_vars = capture_sardine_env_vars()
            
            for var in REDACTED_VARS:
                if var in env_vars:
                    assert env_vars[var] == "[REDACTED]", f"{var} should be redacted"

    def test_unknown_sardine_vars_captured(self):
        """Unknown SARDINE_* vars should still be captured."""
        with mock.patch.dict(os.environ, {
            "SARDINE_NEW_FEATURE_FLAG": "enabled",
        }, clear=False):
            env_vars = capture_sardine_env_vars()
            assert "SARDINE_NEW_FEATURE_FLAG" in env_vars

    def test_password_in_name_redacted(self):
        """Any var with PASSWORD in name should be redacted."""
        with mock.patch.dict(os.environ, {
            "SARDINE_CUSTOM_PASSWORD_VAR": "secret",
        }, clear=False):
            env_vars = capture_sardine_env_vars()
            if "SARDINE_CUSTOM_PASSWORD_VAR" in env_vars:
                assert env_vars["SARDINE_CUSTOM_PASSWORD_VAR"] == "[REDACTED]"


class TestProvenanceInfo:
    """Tests for ProvenanceInfo dataclass."""

    def test_capture_returns_provenance_info(self):
        """capture_provenance should return ProvenanceInfo."""
        prov = capture_provenance()
        assert isinstance(prov, ProvenanceInfo)

    def test_required_fields_present(self):
        """All required fields should be present."""
        prov = capture_provenance()
        
        assert hasattr(prov, "sardine_version")
        assert hasattr(prov, "git_commit")
        assert hasattr(prov, "git_dirty")
        assert hasattr(prov, "sardine_env_vars")
        assert hasattr(prov, "python_version")
        assert hasattr(prov, "platform")
        assert hasattr(prov, "capture_timestamp")
        assert hasattr(prov, "processing_options")
        assert hasattr(prov, "quality_flags")

    def test_to_dict_json_serializable(self):
        """to_dict output should be JSON serializable."""
        prov = capture_provenance()
        prov.set_processing_options({"test": True, "value": 42})
        
        d = prov.to_dict()
        # Should not raise
        json_str = json.dumps(d)
        assert len(json_str) > 0

    def test_to_json_valid(self):
        """to_json should produce valid JSON."""
        prov = capture_provenance()
        json_str = prov.to_json()
        
        # Should parse without error
        parsed = json.loads(json_str)
        assert "sardine_version" in parsed

    def test_add_quality_flag(self):
        """Quality flags should be addable."""
        prov = capture_provenance()
        prov.add_quality_flag("DEBURST_T0_FALLBACK")
        prov.add_quality_flag("RTC_DEFAULT_INCIDENCE")
        prov.add_quality_flag("DEBURST_T0_FALLBACK")  # Duplicate - should not add
        
        assert "DEBURST_T0_FALLBACK" in prov.quality_flags
        assert "RTC_DEFAULT_INCIDENCE" in prov.quality_flags
        assert prov.quality_flags.count("DEBURST_T0_FALLBACK") == 1


class TestProvenanceValidation:
    """Tests for provenance validation."""

    def test_validate_incomplete_provenance(self):
        """Should warn about incomplete provenance."""
        prov = ProvenanceInfo(
            sardine_version="unknown",
            git_commit=None,
            git_dirty=True,
            sardine_env_vars={},
            python_version="3.10",
            platform="Linux",
            capture_timestamp="2025-01-09T00:00:00Z",
            processing_options={},  # Empty
            quality_flags=[],
        )
        
        warnings = validate_provenance_completeness(prov)
        
        # Should have warnings for:
        # - unknown version
        # - no git commit
        # - dirty repo
        # - no processing options
        assert len(warnings) >= 3

    def test_validate_complete_provenance(self):
        """Complete provenance should have minimal warnings."""
        prov = ProvenanceInfo(
            sardine_version="0.3.0",
            git_commit="abc123def456",
            git_dirty=False,
            sardine_env_vars={"SARDINE_STRICT": "1"},
            python_version="3.10",
            platform="Linux",
            capture_timestamp="2025-01-09T00:00:00Z",
            processing_options={"polarization": "VV"},
            quality_flags=[],
        )
        
        warnings = validate_provenance_completeness(prov)
        assert len(warnings) == 0


class TestGetProvenanceDict:
    """Tests for the convenience get_provenance_dict function."""

    def test_returns_dict(self):
        """Should return a dictionary."""
        d = get_provenance_dict()
        assert isinstance(d, dict)

    def test_has_required_keys(self):
        """Should have required keys for metadata embedding."""
        d = get_provenance_dict()
        
        required_keys = [
            "sardine_version",
            "git_commit",
            "sardine_env_vars",
            "capture_timestamp",
        ]
        for key in required_keys:
            assert key in d, f"Missing required key: {key}"


class TestKnownEnvVarsList:
    """Tests for the known env vars list."""

    def test_known_vars_is_list(self):
        """KNOWN_SARDINE_ENV_VARS should be a list."""
        assert isinstance(KNOWN_SARDINE_ENV_VARS, list)

    def test_known_vars_all_start_with_sardine(self):
        """All known vars should start with SARDINE_."""
        for var in KNOWN_SARDINE_ENV_VARS:
            assert var.startswith("SARDINE_"), f"{var} should start with SARDINE_"

    def test_noise_threshold_vars_included(self):
        """NOISE-1 audit: noise threshold vars should be in known list."""
        assert "SARDINE_NOISE_ZERO_THRESHOLD" in KNOWN_SARDINE_ENV_VARS
        assert "SARDINE_NOISE_RATIO_THRESHOLD" in KNOWN_SARDINE_ENV_VARS


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
