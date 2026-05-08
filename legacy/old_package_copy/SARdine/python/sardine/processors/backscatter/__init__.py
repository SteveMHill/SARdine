"""Backscatter processing package."""

from .processor import BackscatterProcessor

# Expose helper modules for consumers that want finer-grained control
from . import qc  # noqa: F401
from . import export  # noqa: F401

# Expose metadata normalization utilities for Python→Rust bridge compatibility
from .metadata_normalize import (  # noqa: F401
    normalize_metadata_for_rust,
    validate_metadata_completeness,
    validate_orbit_vector_count,
    prepare_metadata_for_terrain_correction,
    validate_timing_consistency,
    validate_doppler_centroid,
    validate_geometry_metadata,
    log_metadata_summary,
    KEY_MAPPINGS,
    REQUIRED_KEYS_TERRAIN_CORRECTION,
    MINIMUM_ORBIT_VECTORS,
    MissingMetadataError,
)

__all__ = [
    "BackscatterProcessor",
    "normalize_metadata_for_rust",
    "validate_metadata_completeness",
    "validate_orbit_vector_count",
    "prepare_metadata_for_terrain_correction",
    "validate_timing_consistency",
    "validate_doppler_centroid",
    "validate_geometry_metadata",
    "log_metadata_summary",
    "MissingMetadataError",
]

