"""
I/O helpers for backscatter processing (SAFE discovery, caches, etc.).
"""

import os
from pathlib import Path
from typing import Any, Dict, Tuple


def configure_environment(orbit_cache: Path, dem_cache: Path) -> None:
    """
    Set environment variables for cache directories to keep downstream code consistent.
    """
    os.environ.setdefault("SARDINE_ORBIT_CACHE", str(orbit_cache))
    os.environ.setdefault("SARDINE_DEM_CACHE", str(dem_cache))
