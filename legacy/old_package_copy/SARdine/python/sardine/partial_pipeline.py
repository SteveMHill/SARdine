"""Utilities to run the calibrated+merged portion of the SARdine pipeline.

This module exposes a light-weight driver that executes the IO, deburst,
calibration and TOPSAR merge stages using the cached reader infrastructure.
The calibrated per-subs-wath rasters and the merged stack can optionally be
exported as GeoTIFF files prior to multilooking, terrain flattening or
geocoding.
"""

from __future__ import annotations

import argparse
import json
import math
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Optional, Tuple

import numpy as np

import sardine
from sardine.export import export_to_geotiff

_WGS84_SEMI_MAJOR_AXIS_M = 6_378_137.0
_WGS84_ECCENTRICITY_SQUARED = 6.694_379_990_14e-3


@dataclass
class CalibrationSummary:
    subswath: str
    rows: int
    cols: int
    deburst_seconds: float
    calibration_seconds: float
    output_path: Optional[Path]
    geotransform: Tuple[float, float, float, float, float, float]


def _normalise_product_name(path: Path) -> str:
    name = path.name
    for suffix in (".SAFE", ".safe", ".zip", ".ZIP"):
        if name.endswith(suffix):
            name = name[: -len(suffix)]
            break
    return name


def _parse_subswaths(metadata: Mapping[str, object]) -> List[str]:
    subs = metadata.get("subswaths")
    if not subs or not isinstance(subs, str):
        raise RuntimeError("Metadata is missing the 'subswaths' list")

    raw = [item.strip() for item in subs.split(",") if item.strip()]
    if not raw:
        raise RuntimeError("No subswaths found in metadata")

    def sort_key(sw: str) -> Tuple[int, str]:
        sw_upper = sw.upper()
        if sw_upper.startswith("IW") and sw_upper[2:].isdigit():
            return (0, f"{int(sw_upper[2:]):02d}")
        return (1, sw_upper)

    ordered: List[str] = []
    seen = set()
    for sw in raw:
        sw_upper = sw.upper()
        if sw_upper not in seen:
            seen.add(sw_upper)
            ordered.append(sw_upper)

    return sorted(ordered, key=sort_key)


def _extract_array(result: object, preferred_key: str = "data") -> np.ndarray:
    if isinstance(result, dict):
        status = result.get("status")
        if status == "error":
            message = result.get("message", "operation failed")
            raise RuntimeError(str(message))
        if preferred_key != "data" and preferred_key in result:
            return np.asarray(result[preferred_key])
        if "calibrated_data" in result:
            return np.asarray(result["calibrated_data"])
        if "data" in result:
            return np.asarray(result["data"])
        raise RuntimeError("Operation result does not contain array data")

    return np.asarray(result)


def _compute_geotransform(
    metadata: Mapping[str, object], rows: int, cols: int
) -> Tuple[float, float, float, float, float, float]:
    required_bbox = ("min_longitude", "max_longitude", "min_latitude", "max_latitude")
    missing_bbox = [key for key in required_bbox if key not in metadata]
    if missing_bbox:
        raise RuntimeError(
            "Metadata is missing bounding box fields required for export: "
            + ", ".join(sorted(missing_bbox))
        )

    required_spacing = ("range_pixel_spacing", "azimuth_pixel_spacing")
    missing_spacing = [key for key in required_spacing if key not in metadata]
    if missing_spacing:
        raise RuntimeError(
            "Metadata is missing pixel spacing fields required for export: "
            + ", ".join(sorted(missing_spacing))
        )

    try:
        min_lon = float(metadata["min_longitude"])  # type: ignore[index]
        max_lon = float(metadata["max_longitude"])  # type: ignore[index]
        min_lat = float(metadata["min_latitude"])  # type: ignore[index]
        max_lat = float(metadata["max_latitude"])  # type: ignore[index]
        range_spacing_m = float(metadata["range_pixel_spacing"])  # type: ignore[index]
        azimuth_spacing_m = float(metadata["azimuth_pixel_spacing"])  # type: ignore[index]
    except (TypeError, ValueError) as exc:
        raise RuntimeError("Metadata fields contain non-numeric values") from exc

    if cols <= 0 or rows <= 0:
        raise ValueError("Rows/cols must be positive to compute geotransform")
    if not math.isfinite(range_spacing_m) or range_spacing_m <= 0.0:
        raise RuntimeError(
            f"Invalid range pixel spacing: {range_spacing_m!r} (must be positive)"
        )
    if not math.isfinite(azimuth_spacing_m) or azimuth_spacing_m <= 0.0:
        raise RuntimeError(
            f"Invalid azimuth pixel spacing: {azimuth_spacing_m!r} (must be positive)"
        )

    scene_centre_lat = 0.5 * (min_lat + max_lat)
    lat_rad = math.radians(scene_centre_lat)
    sin_lat = math.sin(lat_rad)
    cos_lat = math.cos(lat_rad)

    denom = math.sqrt(1.0 - _WGS84_ECCENTRICITY_SQUARED * sin_lat * sin_lat)
    if denom == 0.0:
        raise RuntimeError(
            "Cannot compute WGS84 radii at latitude resulting in zero denominator"
        )

    prime_vertical_radius = _WGS84_SEMI_MAJOR_AXIS_M / denom
    meters_per_degree_lon = prime_vertical_radius * cos_lat * math.pi / 180.0

    meridional_radius = (
        _WGS84_SEMI_MAJOR_AXIS_M
        * (1.0 - _WGS84_ECCENTRICITY_SQUARED)
        / (1.0 - _WGS84_ECCENTRICITY_SQUARED * sin_lat * sin_lat) ** 1.5
    )
    meters_per_degree_lat = meridional_radius * math.pi / 180.0

    if meters_per_degree_lon == 0.0 or meters_per_degree_lat == 0.0:
        raise RuntimeError("Derived meters-per-degree values are zero; cannot compute geotransform")

    pixel_width_deg = range_spacing_m / meters_per_degree_lon
    pixel_height_deg = azimuth_spacing_m / meters_per_degree_lat

    if not math.isfinite(pixel_width_deg) or not math.isfinite(pixel_height_deg):
        raise RuntimeError("Computed degree-per-pixel spacing is not finite")

    # GDAL geotransform: (origin_x, pixel_width, x_skew, origin_y, y_skew, pixel_height_negative)
    return (min_lon, pixel_width_deg, 0.0, max_lat, 0.0, -pixel_height_deg)


def _compute_slice_geotransform(
    base_transform: Tuple[float, float, float, float, float, float],
    line_slice: slice,
    sample_slice: slice,
) -> Tuple[float, float, float, float, float, float]:
    """Derive a child geotransform for the specified row/column slice."""

    if line_slice.start is None or sample_slice.start is None:
        raise ValueError("Subswath slices must have explicit start indices for georeferencing")

    origin_x = base_transform[0] + base_transform[1] * float(sample_slice.start)
    origin_y = base_transform[3] + base_transform[5] * float(line_slice.start)

    return (
        origin_x,
        base_transform[1],
        base_transform[2],
        origin_y,
        base_transform[4],
        base_transform[5],
    )


def _validate_geotransform_uniqueness(
    transforms: Mapping[str, Tuple[float, float, float, float, float, float]],
    *,
    decimals: int = 12,
) -> None:
    """Ensure each subswath geotransform is unique within rounding tolerance.

    Parameters
    ----------
    transforms:
        Mapping of subswath identifier to the derived geotransform tuple.
    decimals:
        Number of decimal places to use when quantising values for comparison.

    Raises
    ------
    RuntimeError
        If two subswaths share the same quantised geotransform, indicating the
        derived offsets are not unique.
    """

    seen: Dict[Tuple[float, float, float, float, float, float], str] = {}
    for subswath, transform in transforms.items():
        if len(transform) != 6:
            raise ValueError(
                f"Geotransform for subswath {subswath} must contain 6 values, "
                f"got {len(transform)}"
            )

        quantised = tuple(round(float(value), decimals) for value in transform)
        duplicate = seen.get(quantised)
        if duplicate is not None:
            raise RuntimeError(
                "Duplicate geotransform detected for subswaths "
                f"{duplicate} and {subswath}; derived offsets must be unique"
            )

        seen[quantised] = subswath


def _export_geotiff(
    data: np.ndarray,
    output_path: Path,
    geotransform: Tuple[float, float, float, float, float, float],
    metadata: Mapping[str, object],
) -> Path:
    """Export calibrated data to GeoTIFF with approximate geographic coordinates.
    
    WARNING: This function uses approximate coordinate conversion (meters-to-degrees
    at scene center) and does NOT perform proper SAR geocoding with Range-Doppler
    terrain correction. The output CRS is marked as undefined to prevent misuse
    in spatial analysis.
    
    For scientifically valid geocoding, use the full backscatter processor with
    terrain correction enabled.
    """
    import logging
    logger = logging.getLogger(__name__)
    
    logger.warning(
        "⚠️  APPROXIMATE COORDINATES: Exporting with simple meters-to-degrees conversion.\n"
        "   This is NOT proper SAR geocoding (no Range-Doppler terrain correction).\n"
        "   Output should only be used for quick visual inspection, not scientific analysis.\n"
        "   For accurate geocoding, use the full backscatter pipeline with terrain correction."
    )
    
    data = np.asarray(data, dtype=np.float32)
    export_meta = {
        "POLARIZATION": metadata.get("polarization", "unknown"),
        "PRODUCT_ID": metadata.get("product_id", "unknown"),
        "COORDINATE_TYPE": "APPROXIMATE_SAR_GEOMETRY",
        "WARNING": "Not properly geocoded - meters-to-degrees approximation only",
        "RECOMMENDED_USE": "Visual inspection only - not for spatial analysis",
    }
    
    # Use undefined CRS to prevent misinterpretation as proper geocoding
    # Users who need proper geocoding should use the full terrain correction pipeline
    return export_to_geotiff(
        data=data,
        output_path=output_path,
        geotransform=geotransform,
        crs=None,  # Changed from "EPSG:4326" - undefined CRS for approximate coordinates
        metadata=export_meta,
    )


def _collect_subswath_slices(
    reader: "sardine.SlcReader", polarization: str
) -> Dict[str, Tuple[slice, slice]]:
    geometry = reader.get_all_iw_subswaths()
    by_pol = geometry.get(polarization.upper())
    if not by_pol:
        raise RuntimeError(
            f"No IW subswath geometry available for polarization {polarization}"
        )

    slices: Dict[str, Tuple[slice, slice]] = {}
    for swath_id, info in by_pol.items():
        try:
            first_line = int(info["first_line_global"])
            last_line = int(info["last_line_global"])
            first_sample = int(info["first_sample_global"])
            last_sample = int(info["last_sample_global"])
        except KeyError as exc:  # pragma: no cover - defensive guard
            raise RuntimeError(
                f"Subswath geometry missing required key {exc.args[0]!r}"
            ) from exc

        line_slice = slice(first_line, last_line + 1)
        sample_slice = slice(first_sample, last_sample + 1)
        slices[swath_id.upper()] = (line_slice, sample_slice)

    return slices


def _calibrate_subswath(
    reader: "sardine.SlcReader",
    subswath: str,
    polarization: str,
    calibration_type: str,
    subswath_slices: Mapping[str, Tuple[slice, slice]],
) -> Tuple[np.ndarray, int, int, float, float, slice, slice, Tuple[int, int]]:
    deburst_start = time.perf_counter()
    deburst_result = sardine.deburst_topsar_cached(reader, subswath, polarization)
    deburst_seconds = time.perf_counter() - deburst_start

    deburst_array = _extract_array(deburst_result)
    full_shape = deburst_array.shape

    try:
        line_slice, sample_slice = subswath_slices[subswath]
    except KeyError as exc:
        available = ", ".join(sorted(subswath_slices))
        raise RuntimeError(
            f"No geometry slice for subswath {subswath}. Available: {available}"
        ) from exc

    # NOTE: The deburst output is already local to this subswath - do NOT slice
    # using global coordinates. The slices are kept for positioning during merge.
    # The deburst array represents the complete subswath data.
    complex_data = np.ascontiguousarray(deburst_array, dtype=np.complex64)
    rows, cols = complex_data.shape

    calibration_job = sardine.prepare_calibration_job_cached(
        reader,
        subswath,
        polarization,
        calibration_type,
        int(rows),
        int(cols),
        True,
    )

    calibration_start = time.perf_counter()
    cal_result = sardine.run_calibration_job(calibration_job, complex_data)
    calibration_seconds = time.perf_counter() - calibration_start

    calibrated = _extract_array(cal_result, preferred_key="calibrated_data")
    calibrated = np.asarray(calibrated, dtype=np.float32)

    return (
        calibrated,
        rows,
        cols,
        deburst_seconds,
        calibration_seconds,
        line_slice,
        sample_slice,
        full_shape,
    )


def _merge_subswaths(
    calibrated: Mapping[str, np.ndarray],
    reader: "sardine.SlcReader",
    polarization: str,
) -> Tuple[np.ndarray, Dict[str, object], Optional[np.ndarray]]:
    if not calibrated:
        raise RuntimeError("No calibrated subswaths provided for merge")

    # Normalise arrays to float32 and enforce canonical ordering
    ordered_items = sorted(calibrated.items(), key=lambda item: item[0])
    arrays = {name: np.asarray(data, dtype=np.float32) for name, data in ordered_items}
    upper_names = {name.upper(): arr for name, arr in arrays.items()}

    if {"IW1", "IW2", "IW3"}.issubset(upper_names):
        merged_result = sardine.merge_subswaths_cached(
            upper_names["IW1"],
            upper_names["IW2"],
            upper_names["IW3"],
            reader,
            polarization,
        )
    else:
        merged_result = sardine.topsar_merge_cached(
            arrays,
            polarization,
            reader,
            None,
        )

    merged_array = _extract_array(merged_result)
    info: Dict[str, object] = dict(merged_result) if isinstance(merged_result, dict) else {}
    hit_count = None
    if isinstance(merged_result, dict) and "hit_count" in merged_result:
        hit_count = _extract_array(merged_result, preferred_key="hit_count")

    return np.asarray(merged_array, dtype=np.float32), info, hit_count


def export_calibrated_and_merged_geotiffs(
    input_path: Path,
    output_dir: Path,
    polarization: str = "VV",
    calibration_type: str = "sigma0",
    export_subswaths: bool = True,
    export_merged: bool = True,
) -> Dict[str, object]:
    input_path = input_path.expanduser().resolve()
    output_dir = output_dir.expanduser().resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    sardine.setup_default_orbit_cache(output_dir)

    reader = sardine.create_cached_slc_reader(str(input_path))
    metadata = reader.get_cached_metadata()
    metadata["polarization"] = polarization

    subswaths = _parse_subswaths(metadata)
    product_name = _normalise_product_name(input_path)
    geotiff_records: List[CalibrationSummary] = []
    calibrated_arrays: Dict[str, np.ndarray] = {}
    subswath_slices = _collect_subswath_slices(reader, polarization)

    base_transform: Optional[Tuple[float, float, float, float, float, float]] = None
    subswath_geotransforms: Dict[str, Tuple[float, float, float, float, float, float]] = {}

    for subswath in subswaths:
        (
            calibrated,
            rows,
            cols,
            deburst_s,
            calibration_s,
            line_slice,
            sample_slice,
            full_shape,
        ) = _calibrate_subswath(
            reader, subswath, polarization, calibration_type, subswath_slices
        )
        calibrated_arrays[subswath] = calibrated

        if base_transform is None:
            base_transform = _compute_geotransform(metadata, full_shape[0], full_shape[1])

        geo_transform = _compute_slice_geotransform(base_transform, line_slice, sample_slice)
        output_path = None
        if export_subswaths:
            output_name = f"{product_name}_{subswath.lower()}_calibrated.tif"
            output_path = _export_geotiff(
                calibrated,
                output_dir / output_name,
                geo_transform,
                metadata,
            )

        subswath_geotransforms[subswath] = geo_transform

        geotiff_records.append(
            CalibrationSummary(
                subswath=subswath,
                rows=rows,
                cols=cols,
                deburst_seconds=deburst_s,
                calibration_seconds=calibration_s,
                output_path=output_path,
                geotransform=geo_transform,
            )
        )

    _validate_geotransform_uniqueness(subswath_geotransforms)

    summary: Dict[str, object] = {
        "input": str(input_path),
        "output_dir": str(output_dir),
        "polarization": polarization,
        "calibration_type": calibration_type,
        "subswaths": [record.subswath for record in geotiff_records],
        "subswath_exports": [
            {
                "subswath": record.subswath,
                "rows": record.rows,
                "cols": record.cols,
                "deburst_seconds": record.deburst_seconds,
                "calibration_seconds": record.calibration_seconds,
                "path": str(record.output_path) if record.output_path else None,
                "geotransform": [float(value) for value in record.geotransform],
            }
            for record in geotiff_records
        ],
        "subswath_geotransforms": {
            subswath: [float(value) for value in transform]
            for subswath, transform in subswath_geotransforms.items()
        },
    }

    if export_merged and len(calibrated_arrays) >= 2:
        merged_start = time.perf_counter()
        merged_array, merged_info, merged_hit_count = _merge_subswaths(
            calibrated_arrays, reader, polarization
        )
        summary["merge_seconds"] = time.perf_counter() - merged_start

        geo_transform = _compute_geotransform(
            metadata, merged_array.shape[0], merged_array.shape[1]
        )
        merged_name = f"{product_name}_merged_calibrated.tif"
        merged_path = _export_geotiff(
            merged_array,
            output_dir / merged_name,
            geo_transform,
            metadata,
        )

        summary["merged_raster"] = {
            "rows": int(merged_array.shape[0]),
            "cols": int(merged_array.shape[1]),
            "path": str(merged_path),
            "coverage_percent": merged_info.get("coverage_percent"),
        }

        if merged_hit_count is not None:
            hit_output = output_dir / f"{product_name}_merged_hitcount.tif"
            _export_geotiff(merged_hit_count, hit_output, geo_transform, metadata)
            summary["merged_hitcount"] = str(hit_output)
    else:
        summary["merged_raster"] = None

    return summary


def _dump_summary(summary: Mapping[str, object], output_dir: Path, product_name: str) -> Path:
    summary_path = output_dir / f"{product_name}_merge_summary.json"
    with summary_path.open("w", encoding="utf-8") as fp:
        json.dump(summary, fp, indent=2)
    return summary_path


def _build_argument_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Run the SARdine pipeline through calibration and merge, exporting GeoTIFFs",
    )
    parser.add_argument("input", help="Path to Sentinel-1 SLC ZIP or .SAFE product")
    parser.add_argument("output", help="Directory to write intermediate GeoTIFFs")
    parser.add_argument(
        "--polarization",
        choices=["VV", "VH", "HH", "HV"],
        default="VV",
        help="Polarization to process (default: VV)",
    )
    parser.add_argument(
        "--calibration-type",
        default="sigma0",
        help="Calibration type to request from the Rust backend (default: sigma0)",
    )
    parser.add_argument(
        "--skip-subs",
        action="store_true",
        help="Do not export individual calibrated subswaths",
    )
    parser.add_argument(
        "--skip-merged",
        action="store_true",
        help="Skip the merged GeoTIFF export",
    )
    return parser


def main(argv: Optional[Iterable[str]] = None) -> int:
    parser = _build_argument_parser()
    args = parser.parse_args(argv)

    input_path = Path(args.input)
    output_dir = Path(args.output)
    export_subswaths = not args.skip_subs
    export_merged = not args.skip_merged

    start_time = time.perf_counter()
    summary = export_calibrated_and_merged_geotiffs(
        input_path=input_path,
        output_dir=output_dir,
        polarization=args.polarization,
        calibration_type=args.calibration_type,
        export_subswaths=export_subswaths,
        export_merged=export_merged,
    )
    summary["total_seconds"] = time.perf_counter() - start_time

    product_name = _normalise_product_name(input_path)
    summary_path = _dump_summary(summary, output_dir, product_name)

    print("✅ Calibration + merge exports completed")
    print(f"   • Input: {summary['input']}")
    if summary.get("merged_raster"):
        merged_info = summary["merged_raster"] or {}
        print(
            f"   • Merged GeoTIFF: {merged_info.get('path')} ({merged_info.get('rows')}×{merged_info.get('cols')} pixels)"
        )
    if export_subswaths:
        for record in summary["subswath_exports"]:
            print(
                f"   • {record['subswath']}: {record['path']} ({record['rows']}×{record['cols']} pixels)"
            )
    print(f"   • Summary JSON: {summary_path}")
    print(f"   • Total execution time: {summary['total_seconds']:.1f}s")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
