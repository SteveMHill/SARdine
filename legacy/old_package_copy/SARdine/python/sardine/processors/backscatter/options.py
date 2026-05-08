"""
Option parsing and normalization for the backscatter pipeline.
"""

from pathlib import Path
from typing import Any, Dict, Tuple

DEFAULT_TARGET_RESOLUTION_M = 20.0  # Sentinel-1 RTC default output pixel spacing


def normalize_input_path(input_path: Path) -> str:
    """
    Validate and normalize input path (SAFE dir or ZIP). Returns string path.
    Raises ValueError on invalid input.
    """
    input_path = Path(input_path)
    if not input_path.exists():
        raise ValueError(f"❌ SCIENTIFIC MODE: Input file does not exist: {input_path}")

    if input_path.is_dir() and input_path.suffix == ".SAFE":
        return str(input_path)
    if input_path.suffix.lower() == ".zip":
        return str(input_path)
    raise ValueError("❌ SCIENTIFIC MODE: Input must be either a .zip file or .SAFE directory")


def validate_sentinel_product(input_path: str) -> None:
    """
    Ensure product naming matches Sentinel-1 SLC conventions.
    """
    filename = Path(input_path).name
    if not filename.startswith(("S1A_", "S1B_")):
        raise ValueError(
            "❌ SCIENTIFIC MODE: Not a valid Sentinel-1 product. "
            f"Filename must start with S1A_ or S1B_: {filename}"
        )
    if "_SLC_" not in filename:
        raise ValueError(f"❌ SCIENTIFIC MODE: Only SLC products supported. Found: {filename}")


def resolve_caches(input_path: str, output_dir: Path, options: Dict[str, Any]) -> Tuple[Path, Path]:
    """
    Resolve orbit and DEM cache directories with per-product isolation.
    Returns (orbit_cache_dir, dem_cache_dir).
    """
    orbit_dir_opt = options.get("orbit_dir")
    configured_orbit_cache = orbit_dir_opt or None
    orbit_cache_root = Path(configured_orbit_cache) if configured_orbit_cache else output_dir / "orbit_cache"
    orbit_cache_root.mkdir(parents=True, exist_ok=True)
    orbit_cache_path = orbit_cache_root / Path(input_path).stem if orbit_cache_root.name != Path(input_path).stem else orbit_cache_root
    orbit_cache_path.mkdir(parents=True, exist_ok=True)

    dem_dir_opt = options.get("dem_cache")
    dem_cache_root = Path(dem_dir_opt) if dem_dir_opt else output_dir / "dem_cache"
    dem_cache_root.mkdir(parents=True, exist_ok=True)
    dem_cache_path = dem_cache_root / Path(input_path).stem if dem_cache_root.name != Path(input_path).stem else dem_cache_root
    dem_cache_path.mkdir(parents=True, exist_ok=True)

    return orbit_cache_path, dem_cache_path


def normalize_options(options: Dict[str, Any]) -> Dict[str, Any]:
    """
    Normalize common numeric options (resolution) and defaults.
    """
    normalized = dict(options)
    target_resolution = normalized.get("resolution", DEFAULT_TARGET_RESOLUTION_M)
    try:
        normalized["resolution"] = float(target_resolution)
    except (TypeError, ValueError):
        raise ValueError(f"Invalid resolution option: {target_resolution!r}")
    return normalized
