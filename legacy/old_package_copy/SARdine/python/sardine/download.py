"""High-level Python interface for SARdine download module"""

from typing import Optional, List, Tuple
from pathlib import Path
import datetime

try:
    from sardine._core import DownloadManager as _DownloadManager
except ImportError:
    _DownloadManager = None


class DownloadManager:
    """Unified download manager for Sentinel-1 products, orbit files, and DEM data.
    
    This class provides a user-friendly interface for downloading all types of
    data required for SAR processing.
    
    Example:
        >>> manager = DownloadManager(
        ...     cache_dir="/path/to/cache",
        ...     sources=["aws", "esa"],
        ...     auto_download=True
        ... )
        >>> product_path = manager.download_product("S1A_IW_SLC__1SDV_20200103T170815_...")
        >>> orbit_path = manager.download_orbit("S1A_IW_SLC__1SDV_...", "2020-01-03T17:08:15Z")
        >>> dem_paths = manager.download_dem(bbox=(lon_min, lat_min, lon_max, lat_max))
    """
    
    def __init__(
        self,
        cache_dir: Optional[str] = None,
        sources: Optional[List[str]] = None,
        auto_download: bool = False,
    ):
        """Initialize download manager.
        
        Args:
            cache_dir: Directory for caching downloaded files. Defaults to ~/.sardine/cache
            sources: List of enabled data sources. Options: ["esa", "asf", "aws", "cds"]
                    Defaults to ["aws", "esa"]
            auto_download: If True, automatically download missing data in pipeline
        """
        if _DownloadManager is None:
            raise ImportError(
                "DownloadManager not available. Make sure SARdine is properly installed."
            )
        
        self._manager = _DownloadManager(
            cache_dir=cache_dir,
            sources=sources,
            auto_download=auto_download,
        )
    
    def download_product(self, product_id: str) -> str:
        """Download a Sentinel-1 product by ID.
        
        Args:
            product_id: Sentinel-1 product identifier (e.g., "S1A_IW_SLC__1SDV_20200103T170815_...")
        
        Returns:
            Path to downloaded product file
        """
        return self._manager.download_product(product_id)
    
    def download_orbit(
        self,
        product_id: str,
        start_time: str,
    ) -> str:
        """Download orbit file for a product.
        
        Args:
            product_id: Sentinel-1 product identifier
            start_time: Product start time in ISO 8601 format (e.g., "2020-01-03T17:08:15Z")
        
        Returns:
            Path to downloaded orbit file
        """
        return self._manager.download_orbit(product_id, start_time)
    
    def download_dem(
        self,
        bbox: Tuple[float, float, float, float],
    ) -> List[str]:
        """Download DEM tiles for a bounding box.
        
        Args:
            bbox: Bounding box as (min_lon, min_lat, max_lon, max_lat)
        
        Returns:
            List of paths to downloaded DEM tile files
        """
        min_lon, min_lat, max_lon, max_lat = bbox
        return self._manager.download_dem(min_lon, min_lat, max_lon, max_lat)
    
    @property
    def auto_download(self) -> bool:
        """Check if auto-download is enabled."""
        return self._manager.auto_download()

