"""
Memory-adaptive I/O and caching strategies for SARdine.

This module provides adaptive behavior based on available system RAM,
ensuring the pipeline works efficiently on both high-end workstations
and memory-constrained systems.

Memory Tiers:
- LOW (< 4GB): Sequential processing, no prefetch, aggressive cleanup
- MEDIUM (4-8GB): Single subswath prefetch, moderate chunks
- HIGH (8-16GB): Dual prefetch, larger chunks, some caching  
- VERY_HIGH (> 16GB): Full parallel prefetch, aggressive caching
"""

import logging
from dataclasses import dataclass
from enum import Enum
from typing import Optional, Dict, Any, Tuple

logger = logging.getLogger(__name__)


class MemoryTier(Enum):
    """Memory tier classification for adaptive behavior."""
    LOW = "low"           # < 4 GB available
    MEDIUM = "medium"     # 4-8 GB available
    HIGH = "high"         # 8-16 GB available
    VERY_HIGH = "very_high"  # > 16 GB available


@dataclass
class MemoryConfig:
    """Configuration parameters based on available memory."""
    tier: MemoryTier
    available_gb: float
    
    # I/O settings
    enable_prefetch: bool
    max_prefetch_subswaths: int
    chunk_size_lines: int
    
    # Cache settings
    enable_orbit_cache: bool
    enable_antenna_cache: bool
    enable_lut_cache: bool
    max_cache_mb: int
    
    # Processing settings
    max_parallel_bursts: int
    enable_memory_cleanup: bool
    cleanup_frequency: int  # Cleanup every N stages


def get_available_memory_gb() -> float:
    """Get available system memory in GB.
    
    Returns:
        Available memory in gigabytes. Returns 4.0 as safe default if psutil unavailable.
    """
    try:
        import psutil
        return psutil.virtual_memory().available / (1024**3)
    except ImportError:
        logger.warning("psutil not available, assuming 4GB available memory")
        return 4.0
    except Exception as e:
        logger.warning(f"Failed to get memory info: {e}, assuming 4GB")
        return 4.0


def classify_memory_tier(available_gb: Optional[float] = None) -> MemoryTier:
    """Classify system into a memory tier.
    
    Args:
        available_gb: Override available memory (for testing)
        
    Returns:
        MemoryTier classification
    """
    if available_gb is None:
        available_gb = get_available_memory_gb()
    
    if available_gb < 4:
        return MemoryTier.LOW
    elif available_gb < 8:
        return MemoryTier.MEDIUM
    elif available_gb < 16:
        return MemoryTier.HIGH
    else:
        return MemoryTier.VERY_HIGH


def get_memory_config(available_gb: Optional[float] = None) -> MemoryConfig:
    """Get optimal configuration based on available memory.
    
    Args:
        available_gb: Override available memory (for testing)
        
    Returns:
        MemoryConfig with tier-appropriate settings
    """
    if available_gb is None:
        available_gb = get_available_memory_gb()
    
    tier = classify_memory_tier(available_gb)
    
    configs = {
        MemoryTier.LOW: MemoryConfig(
            tier=MemoryTier.LOW,
            available_gb=available_gb,
            # I/O: Very conservative
            enable_prefetch=False,
            max_prefetch_subswaths=0,
            chunk_size_lines=256,
            # Caching: Minimal
            enable_orbit_cache=True,  # Small, always beneficial
            enable_antenna_cache=False,
            enable_lut_cache=False,
            max_cache_mb=256,
            # Processing: Sequential
            max_parallel_bursts=1,
            enable_memory_cleanup=True,
            cleanup_frequency=1,  # Every stage
        ),
        MemoryTier.MEDIUM: MemoryConfig(
            tier=MemoryTier.MEDIUM,
            available_gb=available_gb,
            # I/O: Single prefetch
            enable_prefetch=True,
            max_prefetch_subswaths=1,
            chunk_size_lines=512,
            # Caching: Selective
            enable_orbit_cache=True,
            enable_antenna_cache=False,
            enable_lut_cache=True,
            max_cache_mb=512,
            # Processing: Limited parallelism
            max_parallel_bursts=2,
            enable_memory_cleanup=True,
            cleanup_frequency=2,
        ),
        MemoryTier.HIGH: MemoryConfig(
            tier=MemoryTier.HIGH,
            available_gb=available_gb,
            # I/O: Dual prefetch
            enable_prefetch=True,
            max_prefetch_subswaths=2,
            chunk_size_lines=1024,
            # Caching: Most caches enabled
            enable_orbit_cache=True,
            enable_antenna_cache=True,
            enable_lut_cache=True,
            max_cache_mb=1024,
            # Processing: Good parallelism
            max_parallel_bursts=4,
            enable_memory_cleanup=True,
            cleanup_frequency=3,
        ),
        MemoryTier.VERY_HIGH: MemoryConfig(
            tier=MemoryTier.VERY_HIGH,
            available_gb=available_gb,
            # I/O: Full prefetch
            enable_prefetch=True,
            max_prefetch_subswaths=3,
            chunk_size_lines=2048,
            # Caching: Everything
            enable_orbit_cache=True,
            enable_antenna_cache=True,
            enable_lut_cache=True,
            max_cache_mb=4096,
            # Processing: Maximum parallelism
            max_parallel_bursts=9,  # All bursts in parallel
            enable_memory_cleanup=False,
            cleanup_frequency=0,
        ),
    }
    
    return configs[tier]


def log_memory_config(config: MemoryConfig) -> None:
    """Log the active memory configuration."""
    tier_emoji = {
        MemoryTier.LOW: "⚠️",
        MemoryTier.MEDIUM: "📊",
        MemoryTier.HIGH: "✅",
        MemoryTier.VERY_HIGH: "🚀",
    }
    
    emoji = tier_emoji.get(config.tier, "📊")
    logger.info(f"{emoji} Memory tier: {config.tier.value} ({config.available_gb:.1f} GB available)")
    logger.info(f"   Prefetch: {config.max_prefetch_subswaths} subswaths, chunks: {config.chunk_size_lines} lines")
    logger.info(f"   Caching: orbit={config.enable_orbit_cache}, antenna={config.enable_antenna_cache}, lut={config.enable_lut_cache}")
    logger.info(f"   Parallelism: {config.max_parallel_bursts} bursts, cleanup every {config.cleanup_frequency} stages")


def estimate_subswath_memory_mb(lines: int, samples: int, include_intermediate: bool = True) -> float:
    """Estimate memory needed for a subswath.
    
    Args:
        lines: Number of azimuth lines
        samples: Number of range samples
        include_intermediate: Include intermediate arrays (deramp, accumulator, etc.)
        
    Returns:
        Estimated memory in MB
    """
    # Complex64 SLC data: 8 bytes per pixel
    slc_mb = lines * samples * 8 / (1024 * 1024)
    
    # Float32 output: 4 bytes per pixel
    output_mb = lines * samples * 4 / (1024 * 1024)
    
    if include_intermediate:
        # Deramp ramps (complex64): ~same as SLC for range-dependent
        deramp_mb = slc_mb
        # Weight accumulator (float32)
        weight_mb = output_mb
        # Hit count (uint16)
        hit_mb = lines * samples * 2 / (1024 * 1024)
        
        return slc_mb + output_mb + deramp_mb + weight_mb + hit_mb
    
    return slc_mb + output_mb


def can_prefetch_subswath(config: MemoryConfig, subswath_lines: int, subswath_samples: int) -> bool:
    """Check if we have enough memory to prefetch a subswath.
    
    Args:
        config: Current memory configuration
        subswath_lines: Lines in the subswath
        subswath_samples: Samples in the subswath
        
    Returns:
        True if prefetching is safe
    """
    if not config.enable_prefetch:
        return False
    
    required_mb = estimate_subswath_memory_mb(subswath_lines, subswath_samples)
    available_mb = config.available_gb * 1024 * 0.5  # Use at most 50% of available
    
    return required_mb < available_mb


def get_optimal_chunk_size(config: MemoryConfig, samples_per_line: int) -> int:
    """Get optimal chunk size based on memory and data dimensions.
    
    Args:
        config: Memory configuration
        samples_per_line: Number of samples per line
        
    Returns:
        Optimal number of lines per chunk
    """
    # Target ~100MB per chunk as baseline
    target_chunk_mb = min(100, config.max_cache_mb // 4)
    bytes_per_line = samples_per_line * 8  # complex64
    lines_for_target = int(target_chunk_mb * 1024 * 1024 / bytes_per_line)
    
    # Clamp to config limits
    return max(64, min(lines_for_target, config.chunk_size_lines))


# Global config instance (set once at startup)
_global_config: Optional[MemoryConfig] = None


def initialize_memory_config(override_gb: Optional[float] = None) -> MemoryConfig:
    """Initialize global memory configuration.
    
    Call this at pipeline startup to set adaptive parameters.
    
    Args:
        override_gb: Override detected memory (for testing)
        
    Returns:
        The active MemoryConfig
    """
    global _global_config
    _global_config = get_memory_config(override_gb)
    log_memory_config(_global_config)
    return _global_config


def get_global_config() -> MemoryConfig:
    """Get the global memory configuration.
    
    Returns:
        Current MemoryConfig, initializes with defaults if not set
    """
    global _global_config
    if _global_config is None:
        _global_config = get_memory_config()
    return _global_config
