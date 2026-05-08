#!/usr/bin/env python3
"""
Test early pipeline steps: Read → Orbit → Deburst → Calibrate → Merge → (Multilook)

Exports intermediate GeoTIFFs for visual inspection.
"""

import os
import sys
import time
import logging
from pathlib import Path
from datetime import datetime

# Set up logging BEFORE importing sardine
os.environ['RUST_LOG'] = 'info'
os.environ['PYTHONUNBUFFERED'] = '1'
os.environ['SARDINE_ORBIT_CACHE'] = '/home/datacube/apps/SARdine/orbit_cache'

import numpy as np

# Configure Python logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s [%(levelname)s] %(message)s',
    handlers=[logging.StreamHandler(sys.stdout)]
)
logger = logging.getLogger(__name__)

def export_array_as_geotiff(array: np.ndarray, output_path: str, description: str = ""):
    """Export a numpy array as a GeoTIFF for visual inspection."""
    try:
        import rasterio
        from rasterio.transform import from_bounds
        
        # Handle complex data - convert to intensity
        if np.iscomplexobj(array):
            logger.info(f"Converting complex data to intensity for {description}")
            array = np.abs(array) ** 2
        
        # Ensure 2D
        if array.ndim == 3:
            array = array[0]  # Take first band/polarization
        
        # Basic stats
        valid_mask = np.isfinite(array) & (array > 0)
        if valid_mask.any():
            valid_data = array[valid_mask]
            logger.info(f"  Array stats: min={valid_data.min():.4e}, max={valid_data.max():.4e}, "
                       f"mean={valid_data.mean():.4e}, valid={valid_mask.sum()}/{array.size}")
        
        # Create a simple transform (pixel coordinates)
        height, width = array.shape
        transform = from_bounds(0, 0, width, height, width, height)
        
        # Write GeoTIFF
        with rasterio.open(
            output_path,
            'w',
            driver='GTiff',
            height=height,
            width=width,
            count=1,
            dtype=array.dtype,
            crs=None,  # Pixel coordinates
            transform=transform,
            compress='lzw'
        ) as dst:
            dst.write(array, 1)
            dst.update_tags(description=description)
        
        file_size = Path(output_path).stat().st_size / (1024 * 1024)
        logger.info(f"✅ Exported: {output_path} ({file_size:.1f} MB)")
        return True
        
    except Exception as e:
        logger.error(f"❌ Failed to export {output_path}: {e}")
        import traceback
        traceback.print_exc()
        return False


def main():
    # Configuration
    SAFE_PATH = Path("/home/datacube/apps/SARdine/data/SLC/S1A_IW_SLC__1SDV_20201005T170824_20201005T170851_034664_04098A_1E66.SAFE")
    ORBIT_DIR = Path("/home/datacube/apps/SARdine/orbit_cache")
    OUTPUT_DIR = Path("/home/datacube/apps/SARdine/SARdine/outputs/early_pipeline_test")
    POLARIZATION = "VV"
    
    # Create output directory
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    log_file = OUTPUT_DIR / f"early_pipeline_{timestamp}.log"
    
    # Add file handler for logging
    file_handler = logging.FileHandler(log_file)
    file_handler.setFormatter(logging.Formatter('%(asctime)s [%(levelname)s] %(message)s'))
    logger.addHandler(file_handler)
    
    logger.info("=" * 70)
    logger.info("🧪 EARLY PIPELINE TEST: Read → Deburst → Calibrate → Merge")
    logger.info("=" * 70)
    logger.info(f"Input: {SAFE_PATH}")
    logger.info(f"Output: {OUTPUT_DIR}")
    logger.info(f"Polarization: {POLARIZATION}")
    
    # Import sardine after setting env vars
    import sardine
    
    total_start = time.time()
    
    # =========================================================================
    # STEP 0: Apply Precise Orbit (required before deburst)
    # =========================================================================
    logger.info("\n" + "=" * 50)
    logger.info("🛰️ STEP 0: Apply Precise Orbit")
    logger.info("=" * 50)
    
    step_start = time.time()
    
    # Extract product ID and start time from SAFE path
    product_id = SAFE_PATH.name.replace('.SAFE', '')
    # Parse start time from product ID: S1A_IW_SLC__1SDV_20201005T170824_...
    # Split by underscore: ['S1A', 'IW', 'SLC', '', '1SDV', '20201005T170824', ...]
    parts = product_id.split('_')
    start_time_str = parts[5]  # 20201005T170824 (index 5, not 4)
    start_time_iso = f"{start_time_str[:4]}-{start_time_str[4:6]}-{start_time_str[6:8]}T{start_time_str[9:11]}:{start_time_str[11:13]}:{start_time_str[13:15]}Z"
    logger.info(f"  Product: {product_id}")
    logger.info(f"  Start time: {start_time_iso}")
    
    try:
        orbit_result = sardine.apply_precise_orbit_file(product_id, start_time_iso, str(ORBIT_DIR))
        logger.info(f"  Orbit vectors: {orbit_result.get('orbit_vectors_count', 'N/A')}")
        logger.info(f"  ✅ Orbit applied in {time.time() - step_start:.2f}s")
    except Exception as e:
        logger.warning(f"  ⚠️ Orbit application failed: {e}")
        logger.warning(f"  Continuing without precise orbit...")
    
    # =========================================================================
    # STEP 1: Create cached reader
    # =========================================================================
    logger.info("\n" + "=" * 50)
    logger.info("📖 STEP 1: Create Cached SLC Reader")
    logger.info("=" * 50)
    
    step_start = time.time()
    reader = sardine.create_cached_slc_reader(str(SAFE_PATH))
    metadata = reader.get_cached_metadata()
    
    logger.info(f"  Product ID: {metadata.get('product_id', 'Unknown')}")
    logger.info(f"  Polarizations: {metadata.get('polarizations', [])}")
    logger.info(f"  Subswaths: {metadata.get('subswaths', [])}")
    logger.info(f"  ✅ Reader created in {time.time() - step_start:.2f}s")
    
    # =========================================================================
    # STEP 2: Get IW subswaths
    # =========================================================================
    logger.info("\n" + "=" * 50)
    logger.info("📡 STEP 2: Get IW Subswaths")
    logger.info("=" * 50)
    
    step_start = time.time()
    subswaths_by_pol = reader.get_all_iw_subswaths()
    
    if POLARIZATION not in subswaths_by_pol:
        logger.error(f"Polarization {POLARIZATION} not found!")
        return 1
    
    # Structure is: {pol: {swath_name: swath_info_dict}}
    subswaths_dict = subswaths_by_pol[POLARIZATION]
    logger.info(f"  Found {len(subswaths_dict)} subswaths for {POLARIZATION}:")
    for sw_name, sw_info in subswaths_dict.items():
        logger.info(f"    {sw_name}: {sw_info.get('burst_count', '?')} bursts, "
                   f"{sw_info.get('range_samples', '?')}x{sw_info.get('azimuth_samples', '?')} samples")
    logger.info(f"  ✅ Subswaths retrieved in {time.time() - step_start:.2f}s")
    
    # Get pixel spacings from metadata for multilook later
    range_spacing = float(metadata.get('range_pixel_spacing', 2.33))
    azimuth_spacing = float(metadata.get('azimuth_pixel_spacing', 13.97))
    logger.info(f"  Pixel spacing: {range_spacing:.2f} m (range) x {azimuth_spacing:.2f} m (azimuth)")
    
    # =========================================================================
    # STEP 3: Deburst each subswath
    # =========================================================================
    logger.info("\n" + "=" * 50)
    logger.info("🔄 STEP 3: Deburst Each Subswath")
    logger.info("=" * 50)
    
    deburst_results = {}  # Store raw SLC data (complex) for calibration
    
    for sw_name in sorted(subswaths_dict.keys()):
        sw_info = subswaths_dict[sw_name]
        logger.info(f"\n  Processing {sw_name}...")
        step_start = time.time()
        
        try:
            # Deburst the subswath using cached reader
            deburst_result = sardine.deburst_topsar_cached(
                reader,
                sw_name,
                POLARIZATION
            )
            
            # Check if result is a dict (error or data with metadata)
            if isinstance(deburst_result, dict):
                if deburst_result.get('status') == 'error':
                    logger.error(f"    ❌ Deburst error: {deburst_result.get('message', 'Unknown error')}")
                    continue
                # Extract the actual data array from the dict
                if 'data' in deburst_result:
                    deburst_data = deburst_result['data']
                else:
                    logger.info(f"    Dict keys: {list(deburst_result.keys())}")
                    # Try to find the numpy array in the dict
                    for k, v in deburst_result.items():
                        if isinstance(v, np.ndarray) and v.ndim >= 2:
                            deburst_data = v
                            logger.info(f"    Found data in key '{k}'")
                            break
                    else:
                        logger.error(f"    ❌ No data array found in deburst result")
                        continue
            else:
                deburst_data = deburst_result
            
            logger.info(f"    Deburst shape: {deburst_data.shape}, dtype: {deburst_data.dtype}")
            
            # Keep complex data for calibration
            deburst_results[sw_name] = deburst_data
            
            # Export deburst result (convert to intensity for visualization)
            if np.iscomplexobj(deburst_data):
                intensity = np.abs(deburst_data) ** 2
            else:
                intensity = deburst_data
            
            output_path = OUTPUT_DIR / f"01_deburst_{sw_name}_{POLARIZATION}.tif"
            export_array_as_geotiff(
                intensity.astype(np.float32),
                str(output_path),
                f"Deburst {sw_name} {POLARIZATION} (intensity)"
            )
            
            logger.info(f"    ✅ {sw_name} deburst in {time.time() - step_start:.2f}s")
            
        except Exception as e:
            logger.error(f"    ❌ Failed to deburst {sw_name}: {e}")
            import traceback
            traceback.print_exc()
    
    # =========================================================================
    # STEP 4: Radiometric Calibration (Sigma0)
    # =========================================================================
    logger.info("\n" + "=" * 50)
    logger.info("📊 STEP 4: Radiometric Calibration (Sigma0)")
    logger.info("=" * 50)
    
    calibrated_results = {}
    
    for sw_name, deburst_data in deburst_results.items():
        logger.info(f"\n  Calibrating {sw_name}...")
        step_start = time.time()
        
        try:
            # Radiometric calibration - expects complex64 SLC data
            # Returns dict with 'data' key containing float32 calibrated values
            cal_result = sardine.radiometric_calibration(
                str(SAFE_PATH),
                sw_name,
                POLARIZATION,
                'sigma0',  # calibration type
                deburst_data  # complex64 SLC data
            )
            
            # Handle dict result
            if isinstance(cal_result, dict):
                if cal_result.get('status') == 'error':
                    logger.error(f"    ❌ Calibration error: {cal_result.get('message', 'Unknown')}")
                    continue
                calibrated_data = cal_result.get('data')
                if calibrated_data is None:
                    logger.error(f"    ❌ No data in calibration result")
                    continue
                logger.info(f"    Valid pixels: {cal_result.get('valid_percentage', 0):.1f}%")
            else:
                calibrated_data = cal_result
            
            logger.info(f"    Calibrated shape: {calibrated_data.shape}, dtype: {calibrated_data.dtype}")
            
            # Stats
            valid = calibrated_data[np.isfinite(calibrated_data) & (calibrated_data > 0)]
            if len(valid) > 0:
                # Convert to dB for display
                db_values = 10 * np.log10(valid + 1e-10)
                logger.info(f"    Sigma0 dB range: [{db_values.min():.1f}, {db_values.max():.1f}] dB")
            
            calibrated_results[sw_name] = calibrated_data
            
            # Export calibrated result
            output_path = OUTPUT_DIR / f"02_calibrated_{sw_name}_{POLARIZATION}.tif"
            export_array_as_geotiff(
                calibrated_data.astype(np.float32),
                str(output_path),
                f"Sigma0 {sw_name} {POLARIZATION}"
            )
            
            logger.info(f"    ✅ {sw_name} calibrated in {time.time() - step_start:.2f}s")
            
        except Exception as e:
            logger.error(f"    ❌ Failed to calibrate {sw_name}: {e}")
            import traceback
            traceback.print_exc()
    
    # =========================================================================
    # STEP 5: TOPSAR Merge
    # =========================================================================
    logger.info("\n" + "=" * 50)
    logger.info("🔗 STEP 5: TOPSAR Merge (IW1 + IW2 + IW3)")
    logger.info("=" * 50)
    
    step_start = time.time()
    merged_data = None
    
    try:
        # Prepare data for merge - dict of subswath name -> data
        subswath_data = {}
        for sw_name in ['IW1', 'IW2', 'IW3']:
            if sw_name in calibrated_results:
                subswath_data[sw_name] = calibrated_results[sw_name]
        
        logger.info(f"  Merging {len(subswath_data)} subswaths: {list(subswath_data.keys())}")
        
        # Call merge function - returns dict with 'data' key
        merge_result = sardine.topsar_merge_cached(
            subswath_data,
            POLARIZATION,
            reader
        )
        
        # Handle dict result
        if isinstance(merge_result, dict):
            if merge_result.get('status') == 'error':
                logger.error(f"  ❌ Merge error: {merge_result.get('message', 'Unknown')}")
            else:
                merged_data = merge_result.get('data')
                if merged_data is None:
                    logger.error(f"  ❌ No data in merge result. Keys: {list(merge_result.keys())}")
        else:
            merged_data = merge_result
        
        if merged_data is not None:
            logger.info(f"  Merged shape: {merged_data.shape}, dtype: {merged_data.dtype}")
            
            # Stats
            valid = merged_data[np.isfinite(merged_data) & (merged_data > 0)]
            if len(valid) > 0:
                db_values = 10 * np.log10(valid + 1e-10)
                logger.info(f"  Merged Sigma0 dB range: [{db_values.min():.1f}, {db_values.max():.1f}] dB")
                logger.info(f"  Valid pixels: {len(valid)}/{merged_data.size} ({100*len(valid)/merged_data.size:.1f}%)")
            
            # Export merged result (linear scale)
            output_path = OUTPUT_DIR / f"03_merged_{POLARIZATION}_linear.tif"
            export_array_as_geotiff(
                merged_data.astype(np.float32),
                str(output_path),
                f"Merged Sigma0 {POLARIZATION} (linear)"
            )
            
            # Export merged result (dB scale for visualization)
            merged_db = 10 * np.log10(np.maximum(merged_data, 1e-10))
            merged_db = np.clip(merged_db, -30, 10)  # Clip to reasonable dB range
            output_path_db = OUTPUT_DIR / f"03_merged_{POLARIZATION}_dB.tif"
            export_array_as_geotiff(
                merged_db.astype(np.float32),
                str(output_path_db),
                f"Merged Sigma0 {POLARIZATION} (dB, clipped -30 to 10)"
            )
            
            logger.info(f"  ✅ Merge completed in {time.time() - step_start:.2f}s")
        
    except Exception as e:
        logger.error(f"  ❌ Merge failed: {e}")
        import traceback
        traceback.print_exc()
    
    # =========================================================================
    # STEP 6: Multilook (optional)
    # =========================================================================
    if merged_data is not None:
        logger.info("\n" + "=" * 50)
        logger.info("👁️ STEP 6: Multilook (3x3)")
        logger.info("=" * 50)
        
        step_start = time.time()
        
        try:
            ml_result = sardine.apply_multilooking(
                merged_data,
                3,  # range_looks
                3,  # azimuth_looks
                range_spacing,  # input_range_spacing
                azimuth_spacing  # input_azimuth_spacing
            )
            
            # Handle dict result
            if isinstance(ml_result, dict):
                multilooked = ml_result.get('data')
                if multilooked is None:
                    logger.error(f"  ❌ No data in multilook result. Keys: {list(ml_result.keys())}")
            else:
                multilooked = ml_result
            
            if multilooked is not None:
                logger.info(f"  Input shape: {merged_data.shape}")
                logger.info(f"  Multilooked shape: {multilooked.shape}")
                
                # Export multilooked (linear)
                output_path = OUTPUT_DIR / f"04_multilooked_{POLARIZATION}_linear.tif"
                export_array_as_geotiff(
                    multilooked.astype(np.float32),
                    str(output_path),
                    f"Multilooked Sigma0 {POLARIZATION} (3x3, linear)"
                )
                
                # Export multilooked (dB)
                ml_db = 10 * np.log10(np.maximum(multilooked, 1e-10))
                ml_db = np.clip(ml_db, -30, 10)
                output_path_db = OUTPUT_DIR / f"04_multilooked_{POLARIZATION}_dB.tif"
                export_array_as_geotiff(
                    ml_db.astype(np.float32),
                    str(output_path_db),
                    f"Multilooked Sigma0 {POLARIZATION} (3x3, dB)"
                )
                
                logger.info(f"  ✅ Multilook completed in {time.time() - step_start:.2f}s")
            
        except Exception as e:
            logger.error(f"  ❌ Multilook failed: {e}")
            import traceback
            traceback.print_exc()
    
    # =========================================================================
    # Summary
    # =========================================================================
    total_time = time.time() - total_start
    
    logger.info("\n" + "=" * 70)
    logger.info("📋 SUMMARY")
    logger.info("=" * 70)
    logger.info(f"Total processing time: {total_time:.1f}s ({total_time/60:.1f} min)")
    logger.info(f"Output directory: {OUTPUT_DIR}")
    logger.info(f"Log file: {log_file}")
    
    # List outputs
    logger.info("\n📁 Generated files:")
    for f in sorted(OUTPUT_DIR.glob("*.tif")):
        size_mb = f.stat().st_size / (1024 * 1024)
        logger.info(f"  {f.name} ({size_mb:.1f} MB)")
    
    logger.info("\n✅ Early pipeline test complete!")
    logger.info("   Open the _dB.tif files in QGIS or any image viewer for visual inspection.")
    
    return 0


if __name__ == "__main__":
    sys.exit(main())
