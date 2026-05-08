"""
Export helpers for backscatter products.

Encapsulates GeoTIFF/NumPy/JSON/text export logic so the processor stays lean.
"""

import json
import time
from datetime import datetime
from pathlib import Path
from typing import List, Optional, Tuple

import numpy as np

from sardine.export import create_cog_with_stac


def export_final_products(processor) -> List[str]:
    """
    Write final products (GeoTIFF/NumPy/JSON) to disk using processor context.

    Returns list of exported filenames.
    """
    # Support both dB and linear output modes
    is_linear_output = getattr(processor, "_linear_output", False)
    db_data = getattr(processor, "_db_data", None)
    
    if not isinstance(db_data, np.ndarray):
        raise RuntimeError("SCIENTIFIC MODE FAILURE: Output data unavailable for export")
    
    # Use appropriate naming based on output type
    output_suffix = "linear" if is_linear_output else "final"
    data_unit = "linear power" if is_linear_output else "dB"

    valid_percentage = getattr(processor, "_valid_percentage", 0.0) or 0.0
    valid_pixel_count = getattr(processor, "_valid_pixel_count", 0) or 0
    total_pixels = getattr(processor, "_total_pixels", int(db_data.size)) or int(db_data.size)

    exported_files: List[str] = []

    try:
        range_looks_used = (
            processor.actual_range_looks
            if getattr(processor, "actual_range_looks", None) is not None
            else processor.multilook_range
        )
        azimuth_looks_used = (
            processor.actual_azimuth_looks
            if getattr(processor, "actual_azimuth_looks", None) is not None
            else processor.multilook_azimuth
        )
        if range_looks_used is None:
            range_looks_used = 1
        if azimuth_looks_used is None:
            azimuth_looks_used = 1
        range_looks_used_f = float(range_looks_used)
        azimuth_looks_used_f = float(azimuth_looks_used)

        if processor.geo_transform is not None:
            output_geotiff = processor.output_dir / f"backscatter_{processor.polarization}_{output_suffix}.tif"

            # Validate CRS before export
            crs_str = processor.coordinate_system
            if not crs_str:
                print("   ⚠️  No coordinate system specified, defaulting to EPSG:4326")
                crs_str = "EPSG:4326"
                processor.coordinate_system = crs_str
            elif not crs_str.upper().startswith("EPSG:"):
                # Try to normalize CRS string
                if crs_str.isdigit():
                    crs_str = f"EPSG:{crs_str}"
                    processor.coordinate_system = crs_str
                    print(f"   ℹ️  Normalized CRS to: {crs_str}")
                else:
                    print(f"   ⚠️  Non-standard CRS format: {crs_str}")

            px_x = abs(processor.geo_transform[1]) if processor.geo_transform else None
            px_y = abs(processor.geo_transform[5]) if processor.geo_transform else None
            range_spacing_m = float(processor.get_current_range_spacing())
            azimuth_spacing_m = float(processor.get_current_azimuth_spacing())
            print(
                "   ✅ Export metadata pixel spacing (m): "
                f"range={range_spacing_m:.3f}, azimuth={azimuth_spacing_m:.3f}"
            )
            if processor.output_epsg == 4326 or px_x is None or px_y is None:
                px_lon_deg = px_x
                px_lat_deg = px_y
                px_x_m = range_spacing_m
                px_y_m = azimuth_spacing_m
            else:
                px_lon_deg = None
                px_lat_deg = None
                px_x_m = px_x
                px_y_m = px_y

            # OPTIMIZATION #97: Single-pass statistics calculation
            # Compute all stats from finite values in one pass (except median)
            try:
                finite_vals = db_data[np.isfinite(db_data)]
                if finite_vals.size > 0:
                    data_min = float(finite_vals.min())
                    data_max = float(finite_vals.max())
                    data_mean = float(finite_vals.mean())
                    data_median = float(np.median(finite_vals))
                else:
                    data_min = data_max = data_mean = data_median = None
            except Exception:
                logging.getLogger(__name__).warning(
                    "Failed to compute dB statistics", exc_info=True
                )
                data_min = data_max = data_mean = data_median = None

            # Extract actual acquisition time from metadata instead of using current time
            acquisition_time = None
            if hasattr(processor, 'metadata') and isinstance(processor.metadata, dict):
                acquisition_time = (
                    processor.metadata.get('start_time') or
                    processor.metadata.get('acquisition_start_time') or
                    processor.metadata.get('startTime')
                )
            if acquisition_time is None:
                acquisition_time = datetime.now().isoformat()
            elif hasattr(acquisition_time, 'isoformat'):
                acquisition_time = acquisition_time.isoformat()
            else:
                acquisition_time = str(acquisition_time)

            # Build metadata keys with appropriate naming for linear vs dB output
            quality_keys = {
                f"quality_{data_unit.replace(' ', '_')}_min": data_min,
                f"quality_{data_unit.replace(' ', '_')}_max": data_max,
                f"quality_{data_unit.replace(' ', '_')}_mean": data_mean,
                f"quality_{data_unit.replace(' ', '_')}_median": data_median,
            }

            sar_metadata = {
                "platform": "sentinel-1",
                "polarization": processor.polarization,
                "processing_level": "COMPLETE_12_STEP_PIPELINE",
                "acquisition_start_time": acquisition_time,
                "orbit_direction": "unknown",
                "range_pixel_spacing": range_spacing_m,
                "azimuth_pixel_spacing": azimuth_spacing_m,
                "pixel_spacing_lon_degrees": px_lon_deg,
                "pixel_spacing_lat_degrees": px_lat_deg,
                "pixel_spacing_x_m": px_x_m,
                "pixel_spacing_y_m": px_y_m,
                "multilook_range": range_looks_used_f,
                "multilook_azimuth": azimuth_looks_used_f,
                "projection_epsg": processor.output_epsg,
                "quality_valid_pixel_percentage": float(valid_percentage),
                "output_unit": data_unit,
                **quality_keys,
            }

            geotiff_path, stac_path = create_cog_with_stac(
                data=db_data,
                output_dir=processor.output_dir,
                filename_base=f"backscatter_{processor.polarization}_{output_suffix}",
                geotransform=tuple(processor.geo_transform),
                sar_metadata=sar_metadata,
                crs=processor.coordinate_system,
                nodata_value=np.nan,  # Use NaN for invalid pixels
                compress="lzw",
            )

            exported_files.extend([geotiff_path.name, stac_path.name])
            print(f"   ✅ PRIMARY: GeoTIFF exported ({data_unit}): {geotiff_path.name}")
        else:
            # No geo_transform (geocoding disabled) - create SAR geometry GeoTIFF
            # Use pixel coordinates with 1:1 pixel spacing in a local Cartesian system
            height, width = db_data.shape
            sar_geotransform = (0.0, 1.0, 0.0, 0.0, 0.0, -1.0)  # Pixel coordinates
            sar_crs = "LOCAL_CS[\"SAR Geometry\"]"  # Local coordinate system
            
            print("   ℹ️  Creating SAR geometry GeoTIFF (geocoding disabled)")

            # Extract acquisition time from metadata
            acquisition_time = None
            if hasattr(processor, 'metadata') and isinstance(processor.metadata, dict):
                acquisition_time = (
                    processor.metadata.get('start_time') or
                    processor.metadata.get('acquisition_start_time') or
                    processor.metadata.get('startTime')
                )
            if acquisition_time is None:
                acquisition_time = datetime.now().isoformat()
            elif hasattr(acquisition_time, 'isoformat'):
                acquisition_time = acquisition_time.isoformat()
            else:
                acquisition_time = str(acquisition_time)

            # Get spacing values
            range_spacing_m = float(processor.get_current_range_spacing()) if processor.get_current_range_spacing() is not None else 1.0
            azimuth_spacing_m = float(processor.get_current_azimuth_spacing()) if processor.get_current_azimuth_spacing() is not None else 1.0

            # Base metadata for SAR geometry output
            sar_metadata = {
                "platform": "sentinel-1",
                "polarization": processor.polarization,
                "processing_level": "COMPLETE_12_STEP_PIPELINE",
                "acquisition_start_time": acquisition_time,
                "orbit_direction": "unknown",
                "range_pixel_spacing": range_spacing_m,
                "azimuth_pixel_spacing": azimuth_spacing_m,
                "multilook_range": float(range_looks_used),
                "multilook_azimuth": float(azimuth_looks_used),
                "quality_valid_pixel_percentage": float(valid_percentage),
                "output_unit": data_unit,
            }
            
            # Update metadata for SAR geometry output
            sar_metadata_local = {
                **sar_metadata,
                "projection_epsg": None,  # No EPSG for local coordinates
                "coordinate_system": "SAR_GEOMETRY",
                "note": "Data in SAR slant-range/azimuth geometry (not geocoded)"
            }
            
            geotiff_path, stac_path = create_cog_with_stac(
                data=db_data,
                output_dir=processor.output_dir,
                filename_base=f"backscatter_{processor.polarization}_{output_suffix}_sar_geometry",
                geotransform=sar_geotransform,
                sar_metadata=sar_metadata_local,
                crs=sar_crs,
                nodata_value=np.nan,
                compress="lzw",
            )
            
            exported_files.extend([geotiff_path.name, stac_path.name])
            print(f"   ✅ SAR GEOMETRY: GeoTIFF exported ({data_unit}): {geotiff_path.name}")

        output_npy = processor.output_dir / f"backscatter_{processor.polarization}_{output_suffix}.npy"
        np.save(output_npy, db_data)
        exported_files.append(output_npy.name)

        try:
            data_valid = np.isfinite(db_data)
            data_vals = db_data[data_valid]
            stats_json = {
                "input": Path(processor.input_path).name,
                "polarization": processor.polarization,
                "calibration_type": processor.calibration_type,
                "output_unit": data_unit,
                "shape": [int(db_data.shape[0]), int(db_data.shape[1])],
                "valid_pixels": int(valid_pixel_count),
                "total_pixels": int(total_pixels),
                "valid_percentage": float(valid_percentage),
                "data_min": float(np.nanmin(data_vals)) if data_vals.size else None,
                "data_max": float(np.nanmax(data_vals)) if data_vals.size else None,
                "data_mean": float(np.nanmean(data_vals)) if data_vals.size else None,
                "data_median": float(np.nanmedian(data_vals)) if data_vals.size else None,
                "range_spacing_m": float(processor.get_current_range_spacing())
                if processor.get_current_range_spacing() is not None
                else None,
                "azimuth_spacing_m": float(processor.get_current_azimuth_spacing())
                if processor.get_current_azimuth_spacing() is not None
                else None,
                "multilook_range": float(range_looks_used_f),
                "multilook_azimuth": float(azimuth_looks_used_f),
                "crs": processor.coordinate_system,
                "epsg": int(processor.output_epsg) if processor.output_epsg is not None else None,
            }
            output_json = processor.output_dir / f"backscatter_{processor.polarization}_summary.json"
            with open(output_json, "w", encoding="utf-8") as jf:
                json.dump(stats_json, jf, indent=2, ensure_ascii=False)
            exported_files.append(output_json.name)
        except Exception as json_error:
            print(f"   ⚠️  Failed to write JSON summary: {json_error}")

        output_txt = processor.output_dir / f"backscatter_{processor.polarization}_summary.txt"
        with open(output_txt, "w", encoding="utf-8") as f:
            f.write("SARdine Backscatter Processing Results\n")
            f.write("=====================================\n")
            f.write(f"Input: {Path(processor.input_path).name}\n")
            f.write(f"Polarization: {processor.polarization}\n")
            f.write(f"Output unit: {data_unit}\n")
            f.write(f"Output shape: {db_data.shape}\n")
            f.write(f"Data range: {np.nanmin(db_data):.6g} to {np.nanmax(db_data):.6g} ({data_unit})\n")
            f.write(
                "Valid pixels: {:,} / {:,} ({:.1f}%)\n".format(
                    valid_pixel_count, total_pixels, valid_percentage
                )
            )
            f.write(f"Processing time: {time.time() - processor.start_time:.1f}s\n")
            f.write(f"Coordinate system: {processor.coordinate_system}\n")
            if processor.geo_transform:
                f.write(f"Geotransform: {processor.geo_transform}\n")
                f.write(
                    f"Pixel size: {abs(processor.geo_transform[1]):.8f}° x {abs(processor.geo_transform[5]):.8f}°\n"
                )

        exported_files.append(output_txt.name)

        print(f"   ✅ Exported files: {', '.join(exported_files)}")
        print(f"   ℹ️  Files saved to: {processor.output_dir}")

        return exported_files

    except Exception:  # pragma: no cover - defensive
        raise
