"""
Provenance and Reproducibility Module
=====================================

This module captures all configuration sources affecting SARdine outputs:
- Package version and git commit hash
- All SARDINE_* environment variables
- Processing options used
- Timestamp and system info

Implements audit requirements PROV-1 and PROV-2 from SCIENTIFIC_AUDIT_REPORT_V2.md.

Usage:
    from sardine.provenance import capture_provenance, ProvenanceInfo
    
    # Capture all provenance at runtime
    prov = capture_provenance()
    print(prov.to_dict())
"""

import os
import sys
import subprocess
from dataclasses import dataclass, asdict, field
from datetime import datetime, timezone
from typing import Dict, Any, Optional, List
from pathlib import Path
import json

# Known SARDINE_* environment variables that affect behavior
# This list must be kept in sync with the Rust code
KNOWN_SARDINE_ENV_VARS = [
    # Orbit/cache
    "SARDINE_ORBIT_CACHE",
    "SARDINE_DOWNLOAD_CACHE",
    "SARDINE_DOWNLOAD_SOURCES",
    
    # Calibration
    "SARDINE_CAL_TILE_COLS",
    "SARDINE_SKIP_NOISE",
    
    # Parsing
    "SARDINE_SERDE_ONLY",
    "SARDINE_REQUIRE_SUBSWATHS",
    "SARDINE_STRICT",
    "SARDINE_STRICT_BBOX",
    "SARDINE_STRICT_CLAMP",
    
    # Terrain correction
    "SARDINE_ENABLE_FAST_SEEDING",
    "SARDINE_SEED_STRIDE",
    "SARDINE_USE_ORBIT_CACHE",
    "SARDINE_SERIAL_TERRAIN",
    "SARDINE_DEM_NODATA_SEA_LEVEL",
    "SARDINE_LOG_CLAMP_BT",
    
    # I/O
    "SARDINE_TIFF_CHUNK_LINES",
    "SARDINE_CHUNK_SIZE",
    
    # Credentials (values redacted)
    "SARDINE_ESA_USERNAME",
    "SARDINE_ESA_PASSWORD",
    "SARDINE_ASF_USERNAME",
    "SARDINE_ASF_PASSWORD",
    
    # Noise thresholds (NOISE-1 audit item)
    "SARDINE_NOISE_ZERO_THRESHOLD",
    "SARDINE_NOISE_RATIO_THRESHOLD",
]

# Credentials that should be redacted in provenance
REDACTED_VARS = {
    "SARDINE_ESA_USERNAME",
    "SARDINE_ESA_PASSWORD",
    "SARDINE_ASF_USERNAME",
    "SARDINE_ASF_PASSWORD",
}


def get_package_version() -> str:
    """Get the sardine package version."""
    try:
        # Try importlib.metadata first (Python 3.8+)
        from importlib.metadata import version
        return version("sardine")
    except Exception:
        pass
    
    try:
        # Fallback to pkg_resources
        import pkg_resources
        return pkg_resources.get_distribution("sardine").version
    except Exception:
        pass
    
    # Fallback: read from pyproject.toml
    try:
        pyproject_path = Path(__file__).parent.parent.parent / "pyproject.toml"
        if pyproject_path.exists():
            content = pyproject_path.read_text()
            for line in content.split("\n"):
                if line.strip().startswith("version"):
                    # version = "0.2.1"
                    return line.split("=")[1].strip().strip('"').strip("'")
    except Exception:
        pass
    
    return "unknown"


def get_git_commit_hash() -> Optional[str]:
    """
    Get the current git commit hash.
    Returns None if not in a git repository or git not available.
    """
    try:
        # Find the repository root
        repo_root = Path(__file__).parent.parent.parent.parent  # SARdine/SARdine/python/sardine -> SARdine
        
        result = subprocess.run(
            ["git", "rev-parse", "HEAD"],
            cwd=repo_root,
            capture_output=True,
            text=True,
            timeout=5,
        )
        if result.returncode == 0:
            return result.stdout.strip()[:12]  # Short hash (12 chars)
    except (subprocess.SubprocessError, FileNotFoundError, OSError):
        pass
    
    return None


def get_git_dirty_status() -> bool:
    """Check if the git repository has uncommitted changes."""
    try:
        repo_root = Path(__file__).parent.parent.parent.parent
        result = subprocess.run(
            ["git", "status", "--porcelain"],
            cwd=repo_root,
            capture_output=True,
            text=True,
            timeout=5,
        )
        if result.returncode == 0:
            return len(result.stdout.strip()) > 0
    except (subprocess.SubprocessError, FileNotFoundError, OSError):
        pass
    return False


def capture_sardine_env_vars() -> Dict[str, str]:
    """
    Capture all SARDINE_* environment variables.
    Credentials are redacted for security.
    """
    env_vars = {}
    
    # Capture known variables (with explicit list)
    for var in KNOWN_SARDINE_ENV_VARS:
        value = os.environ.get(var)
        if value is not None:
            if var in REDACTED_VARS:
                env_vars[var] = "[REDACTED]"
            else:
                env_vars[var] = value
    
    # Also capture any SARDINE_* that we might have missed
    for key, value in os.environ.items():
        if key.startswith("SARDINE_") and key not in env_vars:
            if key in REDACTED_VARS or "PASSWORD" in key or "SECRET" in key or "TOKEN" in key:
                env_vars[key] = "[REDACTED]"
            else:
                env_vars[key] = value
    
    return env_vars


@dataclass
class ProvenanceInfo:
    """
    Complete provenance information for reproducibility.
    
    This dataclass captures all information needed to reproduce
    a SARdine processing run.
    """
    
    # Package info
    sardine_version: str
    git_commit: Optional[str]
    git_dirty: bool
    
    # Environment
    sardine_env_vars: Dict[str, str]
    python_version: str
    platform: str
    
    # Timestamps
    capture_timestamp: str
    
    # Processing options (populated by caller)
    processing_options: Dict[str, Any] = field(default_factory=dict)
    
    # Quality flags (populated by processing stages)
    quality_flags: List[str] = field(default_factory=list)
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary for JSON serialization."""
        return asdict(self)
    
    def to_json(self, indent: int = 2) -> str:
        """Convert to JSON string."""
        return json.dumps(self.to_dict(), indent=indent, default=str)
    
    def add_quality_flag(self, flag: str) -> None:
        """Add a quality flag (e.g., 'DEBURST_T0_FALLBACK', 'RTC_DEFAULT_INCIDENCE')."""
        if flag not in self.quality_flags:
            self.quality_flags.append(flag)
    
    def set_processing_options(self, options: Dict[str, Any]) -> None:
        """Set the processing options used."""
        self.processing_options = options


def capture_provenance() -> ProvenanceInfo:
    """
    Capture complete provenance information at runtime.
    
    Returns:
        ProvenanceInfo dataclass with all provenance data.
    
    Example:
        prov = capture_provenance()
        prov.set_processing_options({
            "polarization": "VV",
            "multilook_range": 3,
            "terrain_flatten": True,
        })
        print(prov.to_json())
    """
    import platform as plat
    
    return ProvenanceInfo(
        sardine_version=get_package_version(),
        git_commit=get_git_commit_hash(),
        git_dirty=get_git_dirty_status(),
        sardine_env_vars=capture_sardine_env_vars(),
        python_version=sys.version.split()[0],
        platform=plat.platform(),
        capture_timestamp=datetime.now(timezone.utc).isoformat(),
        processing_options={},
        quality_flags=[],
    )


def validate_provenance_completeness(prov: ProvenanceInfo) -> List[str]:
    """
    Validate that provenance is complete for reproducibility.
    
    Returns:
        List of warning messages for incomplete provenance.
    """
    warnings = []
    
    if prov.sardine_version == "unknown":
        warnings.append("PROV-2: Package version unknown - reproducibility impaired")
    
    if prov.git_commit is None:
        warnings.append("PROV-2: Git commit hash unavailable - not in git repository")
    
    if prov.git_dirty:
        warnings.append("PROV-2: Git repository has uncommitted changes - exact reproducibility not possible")
    
    if not prov.processing_options:
        warnings.append("PROV-1: Processing options not set - configuration incomplete")
    
    return warnings


# Convenience function for embedding in output metadata
def get_provenance_dict() -> Dict[str, Any]:
    """
    Get provenance as a dictionary suitable for embedding in output metadata.
    This is the primary API for integrating provenance into outputs.
    """
    prov = capture_provenance()
    return prov.to_dict()


if __name__ == "__main__":
    # Self-test
    prov = capture_provenance()
    prov.set_processing_options({
        "polarization": "VV",
        "test": True,
    })
    print("=== SARdine Provenance ===")
    print(prov.to_json())
    print("\n=== Validation ===")
    for warning in validate_provenance_completeness(prov):
        print(f"  ⚠️  {warning}")
