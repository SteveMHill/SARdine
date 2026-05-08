"""
SARdine Utility Functions

Common utility functions for SAR processing including orbit handling,
coordinate conversions, file operations, and quality assessment tools.
"""

from __future__ import annotations

import os
import json
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Union
import numpy as np


def get_orbit_info(orbit_cache_dir: Union[str, Path]) -> Dict[str, Any]:
    """Get information about available precise orbit files.
    
    Args:
        orbit_cache_dir: Path to the orbit cache directory
        
    Returns:
        Dictionary containing:
            - count: Number of orbit files found
            - files: List of first 5 orbit file names
            - status: Human-readable status message
    """
    if not os.path.exists(orbit_cache_dir):
        return {"count": 0, "files": [], "status": "No cache directory"}
    
    orbit_files = [f for f in os.listdir(orbit_cache_dir) if f.endswith('.EOF')]
    return {
        "count": len(orbit_files),
        "files": orbit_files[:5],  # Show first 5 files
        "status": f"Found {len(orbit_files)} orbit file(s)"
    }


def get_sar_info(zip_path: Union[str, Path], polarization: str = "VV") -> Dict[str, str]:
    """Get basic information about a Sentinel-1 SLC product.
    
    Args:
        zip_path: Path to the Sentinel-1 product (ZIP or SAFE)
        polarization: Polarization to check (default: "VV")
        
    Returns:
        Dictionary containing product_type, mission, mode, polarization, status
    """
    try:
        import sardine
        reader = sardine.SlcReader.new_with_full_cache(str(zip_path))
        metadata = reader.get_cached_metadata()
        
        return {
            "product_type": metadata.get("product_type", "Unknown"),
            "mission": metadata.get("mission", "Unknown"),
            "mode": metadata.get("mode", "Unknown"),
            "polarization": polarization,
            "status": "Valid Sentinel-1 SLC"
        }
    except Exception as e:
        return {
            "status": f"Error reading file: {e}",
            "product_type": "Unknown",
            "mission": "Unknown",
            "mode": "Unknown",
            "polarization": polarization
        }


def validate_processing_options(options: Dict[str, Any]) -> Dict[str, Any]:
    """Validate processing options and provide recommendations.
    
    Args:
        options: Dictionary of processing options to validate
        
    Returns:
        Dictionary containing:
            - valid: Boolean indicating if all options are valid
            - warnings: List of warning messages
            - recommendations: List of suggested improvements
    """
    warnings: List[str] = []
    recommendations: List[str] = []
    
    # Check speckle filter
    valid_filters = ['lee', 'enhanced_lee', 'gamma_map', 'refined_lee', 'frost']
    if options.get('speckle_filter') not in valid_filters:
        warnings.append(f"Unknown speckle filter: {options.get('speckle_filter')}")
        recommendations.append("Use 'enhanced_lee' for best quality")
    
    # Check multilook factors
    if options.get('multilook_range', 1) < 1 or options.get('multilook_azimuth', 1) < 1:
        warnings.append("Multilook factors must be >= 1")
    
    # Check resolution
    resolution = options.get('resolution', 10.0)
    if resolution < 5.0:
        warnings.append("Resolution < 5m may be computationally expensive")
    elif resolution > 100.0:
        warnings.append("Resolution > 100m may lose important details")
    
    return {
        "valid": len(warnings) == 0,
        "warnings": warnings,
        "recommendations": recommendations
    }


def calculate_processing_time_estimate(
    slc_dimensions: Optional[Tuple[int, int]], 
    options: Dict[str, Any]
) -> Dict[str, Any]:
    """Estimate processing time based on data size and options.
    
    Args:
        slc_dimensions: Tuple of (rows, cols) for the SLC data, or None
        options: Processing options dictionary
        
    Returns:
        Dictionary containing estimated time and contributing factors
    """
    rows, cols = slc_dimensions if slc_dimensions else (1000, 1000)
    total_pixels = rows * cols
    
    # Base time per million pixels
    base_time_per_mpx = 2.0  # seconds
    
    # Factors that affect processing time
    speckle_factor = 1.5 if options.get('speckle_filter') == 'enhanced_lee' else 1.0
    terrain_factor = 2.0 if options.get('terrain_flatten') else 1.0
    geocode_factor = 3.0 if options.get('geocode') else 1.0
    
    estimated_seconds = (total_pixels / 1e6) * base_time_per_mpx * speckle_factor * terrain_factor * geocode_factor
    
    return {
        "estimated_time_seconds": estimated_seconds,
        "estimated_time_minutes": estimated_seconds / 60,
        "factors": {
            "data_size_mpx": total_pixels / 1e6,
            "speckle_filter_factor": speckle_factor,
            "terrain_correction_factor": terrain_factor,
            "geocoding_factor": geocode_factor
        }
    }


def create_processing_summary(
    input_path: Union[str, Path, None], 
    output_dir: Union[str, Path], 
    options: Dict[str, Any], 
    success: bool = True
) -> Dict[str, Any]:
    """Create a summary of processing parameters and outputs.
    
    Args:
        input_path: Path to input file (may be None)
        output_dir: Output directory path
        options: Processing options used
        success: Whether processing was successful
        
    Returns:
        Summary dictionary with input, output, and status information
    """
    
    # List output files
    output_files = []
    if os.path.exists(output_dir):
        for file in os.listdir(output_dir):
            if file.endswith(('.tif', '.npy', '.json', '.xml')):
                file_path = Path(output_dir) / file
                output_files.append({
                    "name": file,
                    "size_mb": file_path.stat().st_size / (1024 * 1024) if file_path.exists() else 0,
                    "type": file.split('.')[-1].upper()
                })
    
    summary = {
        "input": {
            "path": str(input_path),
            "exists": os.path.exists(input_path) if input_path else False
        },
        "output": {
            "directory": str(output_dir),
            "file_count": len(output_files),
            "files": output_files
        },
        "processing_options": options,
        "status": "SUCCESS" if success else "FAILED"
    }
    
    return summary


def check_system_resources() -> Dict[str, Any]:
    """Check available system resources for processing.
    
    Returns:
        Dictionary containing memory, disk, and CPU information
        with recommendations for processing feasibility.
        
    Note:
        Requires psutil package to be installed.
    """
    import psutil
    
    # Get memory info
    memory = psutil.virtual_memory()
    
    # Get disk space for common paths
    disk_usage = {}
    for path in ["/tmp", "/home", "."]:
        try:
            if os.path.exists(path):
                usage = psutil.disk_usage(path)
                disk_usage[path] = {
                    "free_gb": usage.free / (1024**3),
                    "total_gb": usage.total / (1024**3),
                    "used_percent": (usage.used / usage.total) * 100
                }
        except (OSError, PermissionError):
            # Skip paths that are inaccessible (network drives, permission issues)
            pass
    
    return {
        "memory": {
            "available_gb": memory.available / (1024**3),
            "total_gb": memory.total / (1024**3),
            "used_percent": memory.percent
        },
        "disk": disk_usage,
        "cpu_count": psutil.cpu_count(),
        "recommendations": {
            "memory_warning": memory.available < 4 * (1024**3),  # Less than 4GB
            "disk_warning": any(d.get("free_gb", 0) < 10 for d in disk_usage.values())  # Less than 10GB
        }
    }


def format_file_size(size_bytes: int) -> str:
    """Format file size in human-readable format.
    
    Args:
        size_bytes: Size in bytes
        
    Returns:
        Formatted string (e.g., "1.5 MB", "2.3 GB")
    """
    if size_bytes < 1024:
        return f"{size_bytes} B"
    elif size_bytes < 1024**2:
        return f"{size_bytes / 1024:.1f} KB"
    elif size_bytes < 1024**3:
        return f"{size_bytes / (1024**2):.1f} MB"
    else:
        return f"{size_bytes / (1024**3):.1f} GB"


def create_cache_directories(base_dir: Union[str, Path]) -> Dict[str, Any]:
    """Create necessary cache directories for processing.
    
    Args:
        base_dir: Base directory for cache structure
        
    Returns:
        Dictionary with base_cache path, created dirs list, and count
    """
    cache_dirs = [
        "orbit",
        "dem", 
        "temp",
        "logs"
    ]
    
    created_dirs = []
    base_path = Path(base_dir) / "cache"
    
    for cache_dir in cache_dirs:
        dir_path = base_path / cache_dir
        dir_path.mkdir(parents=True, exist_ok=True)
        created_dirs.append(str(dir_path))
    
    return {
        "base_cache": str(base_path),
        "created": created_dirs,
        "count": len(created_dirs)
    }


def cleanup_temp_files(
    output_dir: Union[str, Path], 
    keep_logs: bool = True
) -> Dict[str, Any]:
    """Clean up temporary processing files.
    
    Args:
        output_dir: Directory to clean
        keep_logs: If True, preserve log files (default: True)
        
    Returns:
        Dictionary with removed file count and list
    """
    temp_patterns = [
        "*.tmp",
        "*.temp", 
        "temp_*",
        "*.partial"
    ]
    
    if not keep_logs:
        temp_patterns.extend(["*.log", "processing_log*.json"])
    
    removed_files = []
    output_path = Path(output_dir)
    
    for pattern in temp_patterns:
        for file_path in output_path.rglob(pattern):
            try:
                file_path.unlink()
                removed_files.append(str(file_path))
            except (OSError, PermissionError):
                # Skip files that can't be removed (in use, permission denied)
                pass
    
    return {
        "removed_count": len(removed_files),
        "removed_files": removed_files[:10],  # Show first 10
        "patterns_cleaned": temp_patterns
    }
