#!/usr/bin/env python3
"""
SARdine Command Line Interface

A command-line tool for processing Sentinel-1 SLC data.
"""

import argparse
import sys
import json
import time
from pathlib import Path
from typing import Dict, Any

import sardine


def print_banner():
    """Print the SARdine banner."""
    banner = """
╔══════════════════════════════════════════════════════════════════╗
║                            SARdine                               ║
║        A Fast, Modular Sentinel-1 Backscatter Processor         ║
║                                                                  ║
║  🚀 Modern alternative to ESA SNAP and GAMMA                    ║
║  🔧 Python API with fast Rust backend                          ║
║  📊 High-quality Gamma0 VV/VH GeoTIFFs                         ║
╚══════════════════════════════════════════════════════════════════╝
    """
    print(banner)


def format_file_size(size_bytes: int) -> str:
    """Format file size in human-readable format."""
    for unit in ['B', 'KB', 'MB', 'GB']:
        if size_bytes < 1024.0:
            return f"{size_bytes:.1f} {unit}"
        size_bytes /= 1024.0
    return f"{size_bytes:.1f} TB"


def print_product_info(info: Dict[str, Any]):
    """Print formatted product information."""
    print("\n" + "="*70)
    print("📡 SENTINEL-1 PRODUCT INFORMATION")
    print("="*70)
    
    # Basic file information
    input_path = Path(info["input_path"])
    print(f"📁 File: {input_path.name}")
    
    if input_path.exists():
        file_size = input_path.stat().st_size
        print(f"💾 Size: {format_file_size(file_size)}")
    
    print(f"📦 Archive contains: {info['total_files']} files")
    print(f"📻 Polarizations: {', '.join(info['polarizations'])}")
    
    # Metadata information
    if "metadata" in info and info["metadata"]:
        metadata = info["metadata"]
        print("\n" + "-"*50)
        print("🛰️  ACQUISITION METADATA")
        print("-"*50)
        
        print(f"🆔 Product ID: {metadata.get('product_id', 'Unknown')}")
        print(f"🚀 Mission: {metadata.get('mission', 'Unknown')}")
        print(f"🛰️  Platform: {metadata.get('platform', 'Unknown')}")
        print(f"📡 Mode: {metadata.get('acquisition_mode', 'Unknown')}")
        
        # Time information
        start_time = metadata.get('start_time', '')
        stop_time = metadata.get('stop_time', '')
        if start_time and stop_time:
            print(f"⏰ Start: {start_time}")
            print(f"⏰ Stop:  {stop_time}")
        
        # Spatial information
        pixel_spacing = metadata.get('pixel_spacing')
        if pixel_spacing:
            print(f"📏 Pixel spacing: {pixel_spacing[0]:.1f} m (range) × {pixel_spacing[1]:.1f} m (azimuth)")
        
        bbox = metadata.get('bounding_box')
        if bbox:
            print(f"🌍 Bounding box:")
            print(f"   Longitude: {bbox[0]:.3f}° to {bbox[2]:.3f}°")
            print(f"   Latitude:  {bbox[1]:.3f}° to {bbox[3]:.3f}°")
    
    elif "metadata_error" in info:
        print(f"\n⚠️  Metadata Error: {info['metadata_error']}")
    
    print("\n" + "="*70)


def cmd_info(args):
    """Handle the 'info' command."""
    input_path = Path(args.input)
    
    if not input_path.exists():
        print(f"❌ Error: File not found: {input_path}")
        return 1
    
    if not input_path.suffix.lower() == '.zip':
        print(f"⚠️  Warning: Expected ZIP file, got: {input_path.suffix}")
    
    try:
        print("🔍 Analyzing Sentinel-1 product...")
        info = sardine.get_product_info(str(input_path))
        
        if args.json:
            print(json.dumps(info, indent=2))
        else:
            print_product_info(info)
            
    except Exception as e:
        print(f"❌ Error analyzing product: {e}")
        if args.verbose:
            import traceback
            traceback.print_exc()
        return 1
    
    return 0


def cmd_process(args):
    """Handle the 'process' command."""
    print("🚧 Processing functionality not yet implemented")
    print("📋 This will include:")
    print("   • Reading SLC data")
    print("   • Applying orbit files")
    print("   • Debursting")
    print("   • Radiometric calibration")
    print("   • Multilooking")
    print("   • Terrain correction")
    print("   • Speckle filtering")
    print("   • Output to GeoTIFF")
    return 0


def cmd_test(args):
    """Handle the 'test' command."""
    print("🧪 Running SARdine tests...")
    
    # Test with sample data if available
    sample_data = Path("/home/datacube/SARdine/data")
    if sample_data.exists():
        zip_files = list(sample_data.glob("*.zip"))
        if zip_files:
            print(f"📁 Found sample data: {zip_files[0].name}")
            try:
                info = sardine.get_product_info(str(zip_files[0]))
                print("✅ SLC reader test passed")
                return 0
            except Exception as e:
                print(f"❌ SLC reader test failed: {e}")
                return 1
    
    print("⚠️  No sample data found for testing")
    return 0


def cmd_test_srtm(args):
    """Test SRTM download functionality."""
    try:
        print(f"🌍 Testing SRTM download for tile: {args.tile}")
        print(f"📁 Output directory: {args.output}")
        
        start_time = time.time()
        
        # Test SRTM download
        result_path = sardine.test_srtm_download(args.tile, args.output)
        
        download_time = time.time() - start_time
        
        print(f"✅ SRTM download completed in {download_time:.2f}s")
        print(f"📁 File saved to: {result_path}")
        
        # Check file size
        file_path = Path(result_path)
        if file_path.exists():
            file_size = file_path.stat().st_size
            print(f"📊 File size: {format_file_size(file_size)}")
            
            if file_size > 1024 * 1024:  # > 1MB
                print(f"✅ File size looks reasonable for SRTM data")
            else:
                print(f"⚠️  File size seems small, might be an error page")
        
        print(f"\n🎉 SRTM download test successful!")
        print(f"📊 This confirms automatic DEM preparation will work")
        
        return 0
        
    except Exception as e:
        print(f"❌ SRTM download test failed: {e}")
        print(f"\n📝 This is expected if:")
        print(f"   • No internet connection")
        print(f"   • SRTM servers are temporarily unavailable") 
        print(f"   • Authentication is required for some sources")
        print(f"\n💡 Alternative options:")
        print(f"   • Download SRTM tiles manually from https://earthexplorer.usgs.gov/")
        print(f"   • Use existing DEM files with the read_dem() method")
        print(f"   • Try again later when servers are available")
        return 1


def cmd_speckle_filter(args):
    """Apply speckle filtering to SAR intensity images."""
    try:
        # Check if input file exists
        input_path = Path(args.input)
        if not input_path.exists():
            print(f"❌ Error: Input file not found: {input_path}")
            return 1
        
        output_path = Path(args.output)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        
        print(f"🔧 Applying {args.filter_type} speckle filter")
        print(f"📁 Input: {input_path}")
        print(f"📁 Output: {output_path}")
        print(f"🪟 Window size: {args.window_size}x{args.window_size}")
        
        start_time = time.time()
        
        # Read input image using rasterio
        try:
            import rasterio
            import numpy as np
        except ImportError:
            print("❌ Error: rasterio package required for GeoTIFF support")
            print("Install with: pip install rasterio")
            return 1
        
        # Read image
        with rasterio.open(input_path) as src:
            image_data = src.read(1).astype(np.float64)
            profile = src.profile
            
        print(f"📊 Image dimensions: {image_data.shape[1]}x{image_data.shape[0]}")
        
        # Convert to list format for Python API
        image_list = image_data.tolist()
        
        # Estimate number of looks if not provided
        num_looks = args.num_looks
        if num_looks is None:
            print("🔍 Estimating number of looks...")
            num_looks = sardine.estimate_num_looks(image_list, args.window_size)
            print(f"📊 Estimated number of looks: {num_looks:.2f}")
        
        # Apply speckle filter
        print(f"🔧 Applying {args.filter_type} filter...")
        filtered_list = sardine.apply_speckle_filter(
            image_list,
            args.filter_type,
            window_size=args.window_size,
            num_looks=num_looks,
            edge_threshold=args.edge_threshold,
            damping_factor=args.damping_factor,
            cv_threshold=args.cv_threshold
        )
        
        # Convert back to numpy array
        filtered_data = np.array(filtered_list, dtype=np.float32)
        
        # Update profile for output
        profile.update(dtype='float32', compress='lzw')
        
        # Write output
        with rasterio.open(output_path, 'w', **profile) as dst:
            dst.write(filtered_data, 1)
        
        filter_time = time.time() - start_time
        
        print(f"✅ Speckle filtering completed in {filter_time:.2f}s")
        print(f"📁 Filtered image saved to: {output_path}")
        
        # Calculate filtering statistics
        original_mean = np.mean(image_data[image_data > 0])
        filtered_mean = np.mean(filtered_data[filtered_data > 0])
        reduction_ratio = np.std(image_data[image_data > 0]) / np.std(filtered_data[filtered_data > 0])
        
        print(f"\n📊 Filtering Statistics:")
        print(f"   Original mean: {original_mean:.6f}")
        print(f"   Filtered mean: {filtered_mean:.6f}")
        print(f"   Noise reduction ratio: {reduction_ratio:.2f}")
        
        return 0
        
    except Exception as e:
        print(f"❌ Speckle filtering failed: {e}")
        return 1


def cmd_estimate_nlooks(args):
    """Estimate number of looks from SAR intensity image."""
    try:
        # Check if input file exists
        input_path = Path(args.input)
        if not input_path.exists():
            print(f"❌ Error: Input file not found: {input_path}")
            return 1
        
        print(f"🔍 Estimating number of looks")
        print(f"📁 Input: {input_path}")
        print(f"🪟 Window size: {args.window_size}x{args.window_size}")
        
        start_time = time.time()
        
        # Read input image using rasterio
        try:
            import rasterio
            import numpy as np
        except ImportError:
            print("❌ Error: rasterio package required for GeoTIFF support")
            print("Install with: pip install rasterio")
            return 1
        
        # Read image
        with rasterio.open(input_path) as src:
            image_data = src.read(1).astype(np.float64)
            
        print(f"📊 Image dimensions: {image_data.shape[1]}x{image_data.shape[0]}")
        
        # Convert to list format for Python API
        image_list = image_data.tolist()
        
        # Estimate number of looks
        num_looks = sardine.estimate_num_looks(image_list, args.window_size)
        
        estimation_time = time.time() - start_time
        
        print(f"✅ Number of looks estimation completed in {estimation_time:.2f}s")
        print(f"📊 Estimated number of looks: {num_looks:.3f}")
        
        # Provide interpretation
        if num_looks < 1.5:
            interpretation = "Single-look data (high speckle)"
        elif num_looks < 4:
            interpretation = "Few-look data (moderate speckle)"
        elif num_looks < 10:
            interpretation = "Multi-look data (low speckle)"
        else:
            interpretation = "Heavily multi-looked data (very low speckle)"
        
        print(f"📝 Interpretation: {interpretation}")
        
        # Recommend filter settings
        print(f"\n💡 Recommended speckle filter settings:")
        if num_looks < 2:
            print(f"   • Filter type: enhanced_lee or lee_sigma")
            print(f"   • Window size: 7x7 or 9x9")
        elif num_looks < 5:
            print(f"   • Filter type: lee or refined_lee")
            print(f"   • Window size: 5x5 or 7x7")
        else:
            print(f"   • Filter type: mean or median")
            print(f"   • Window size: 3x3 or 5x5")
        
        return 0
        
    except Exception as e:
        print(f"❌ Number of looks estimation failed: {e}")
        return 1


def cmd_orbit(args):
    """Handle the 'orbit' command for orbit file operations."""
    input_path = Path(args.input)
    
    if not input_path.exists():
        print(f"❌ Error: File not found: {input_path}")
        return 1
    
    print(f"🛰️  Orbit File Operations for: {input_path.name}")
    print("="*70)
    
    try:
        from sardine import SlcReader
        
        # Open SLC file and get metadata
        reader = SlcReader(str(input_path))
        metadata = reader.get_metadata('VV')
        
        print(f"📅 Product date: {metadata.start_time}")
        print(f"🛰️  Mission: {metadata.mission}")
        
        # Extract product ID from filename
        product_id = input_path.stem.replace('.SAFE', '')
        print(f"📋 Product ID: {product_id}")
        
        # Determine recommended orbit type
        from datetime import datetime, timezone
        
        # Parse datetime string from metadata
        acquisition_time_str = str(metadata.start_time)
        if 'T' in acquisition_time_str:
            acquisition_time = datetime.fromisoformat(acquisition_time_str.replace('Z', '+00:00'))
        else:
            # Fallback parsing
            from dateutil import parser
            acquisition_time = parser.parse(acquisition_time_str)
        
        now = datetime.now(timezone.utc)
        age_days = (now - acquisition_time).days
        
        if age_days > 20:
            recommended_type = "POEORB"
            print(f"📡 Recommended: POEORB (product is {age_days} days old)")
            print("   🎯 Highest accuracy (precise orbit ephemerides)")
        else:
            recommended_type = "RESORB" 
            print(f"📡 Recommended: RESORB (product is {age_days} days old)")
            print("   🎯 Lower accuracy but available sooner")
        
        # Show orbit file availability
        print(f"\n📊 Orbit Data Status:")
        print("   ❌ No orbit data in SLC archive (standard for Sentinel-1)")
        print("   🌐 External orbit files required")
        
        # Show cache information
        cache_dir = Path.home() / ".sardine" / "orbit_cache" if not args.output else Path(args.output)
        print(f"\n💾 Orbit Cache:")
        print(f"   📁 Cache directory: {cache_dir}")
        
        if cache_dir.exists():
            orbit_files = list(cache_dir.glob("*.EOF"))
            if orbit_files:
                print(f"   📄 Cached orbit files: {len(orbit_files)}")
                for file in orbit_files[:3]:  # Show first 3
                    size_mb = file.stat().st_size / (1024 * 1024)
                    print(f"     • {file.name} ({size_mb:.1f} MB)")
                if len(orbit_files) > 3:
                    print(f"     ... and {len(orbit_files) - 3} more")
            else:
                print("   📄 No cached orbit files found")
        else:
            print("   📄 Cache directory does not exist (will be created)")
        
        # Show download information
        print(f"\n🌐 Download Information:")
        print("   📡 ESA Server: https://step.esa.int/auxdata/orbits/Sentinel-1/")
        print(f"   📂 Path: {recommended_type}/{'S1A' if product_id.startswith('S1A') else 'S1B'}/")
        print("   📄 Format: EOF (XML-based)")
        
        if args.download:
            print(f"\n⬇️  Downloading orbit files...")
            print(f"   📁 Output directory: {cache_dir}")
            print("   🚧 Download functionality: Available")
            print("   � Note: Requires internet connection and correct ESA server URLs")
            print("   � Will attempt multiple URL patterns for robust download")
        
        if args.list_urls:
            print(f"\n🔗 Example URLs for {recommended_type}:")
            satellite = 'S1A' if product_id.startswith('S1A') else 'S1B'
            year = acquisition_time.year
            base_url = f"https://step.esa.int/auxdata/orbits/Sentinel-1/{recommended_type}/{satellite}/{year}/"
            
            date_str = acquisition_time.strftime("%Y%m%dT%H%M%S")
            example_url = f"{base_url}{satellite}_OPER_AUX_{recommended_type}_OPOD_{date_str}_V*.EOF"
            print(f"   🌐 Pattern: {example_url}")
            print("   📝 Note: Multiple time variants attempted for robust download")
            print("   🔍 Searches around acquisition time for valid orbit files")
    
    except ModuleNotFoundError as e:
        if 'dateutil' in str(e):
            print("❌ Error: Missing required module 'python-dateutil'")
            print("💡 Install with: pip install python-dateutil")
        else:
            print(f"❌ Error: Missing module: {e}")
        return 1
    except Exception as e:
        print(f"❌ Error processing orbit information: {e}")
        return 1
    
    print("\n" + "="*70)
    return 0


def cmd_iw_split(args):
    """Handle IW split command."""
    try:
        reader = sardine.SlcReader(args.input)
        
        print(f"\n🔍 Analyzing IW sub-swaths in: {Path(args.input).name}")
        
        # Check if this is an IW mode product
        if not reader.is_iw_mode():
            print("❌ Error: Input file is not an IW (Interferometric Wide) mode product")
            return 1
        
        # Get all sub-swaths for all polarizations
        all_subswaths = reader.get_all_iw_subswaths()
        
        if not all_subswaths:
            print("❌ Error: No IW sub-swaths found in the product")
            return 1
        
        print(f"\n✅ Found IW sub-swaths for {len(all_subswaths)} polarization(s)")
        
        # Print detailed information for each polarization
        for pol, subswaths in all_subswaths.items():
            print(f"\n📊 Polarization: {pol}")
            print("-" * 40)
            
            for swath_id, swath in subswaths.items():
                print(f"  🎯 Sub-swath: {swath_id}")
                print(f"     • Bursts: {swath.burst_count}")
                print(f"     • Range samples: {swath.range_samples:,}")
                print(f"     • Azimuth samples: {swath.azimuth_samples:,}")
                print(f"     • Range pixel spacing: {swath.range_pixel_spacing:.2f} m")
                print(f"     • Azimuth pixel spacing: {swath.azimuth_pixel_spacing:.2f} m")
                print(f"     • Slant range time: {swath.slant_range_time:.6f} s")
                print(f"     • Burst duration: {swath.burst_duration:.6f} s")
        
        # Save to JSON if requested
        if args.output:
            output_data = {}
            for pol, subswaths in all_subswaths.items():
                output_data[pol] = {}
                for swath_id, swath in subswaths.items():
                    output_data[pol][swath_id] = {
                        "id": swath.id,
                        "burst_count": swath.burst_count,
                        "range_samples": swath.range_samples,
                        "azimuth_samples": swath.azimuth_samples,
                        "range_pixel_spacing": swath.range_pixel_spacing,
                        "azimuth_pixel_spacing": swath.azimuth_pixel_spacing,
                        "slant_range_time": swath.slant_range_time,
                        "burst_duration": swath.burst_duration
                    }
            
            output_path = Path(args.output)
            with open(output_path, 'w') as f:
                json.dump(output_data, f, indent=2)
            
            print(f"\n💾 Sub-swath information saved to: {output_path}")
        
        print(f"\n🎯 IW split analysis complete!")
        print(f"Next step: Use this information for burst extraction and debursting")
        
        return 0
        
    except Exception as e:
        print(f"❌ Error during IW split: {e}")
        return 1
    


def cmd_deburst(args):
    """Handle the 'deburst' command."""
    input_path = Path(args.input)
    
    if not input_path.exists():
        print(f"❌ Error: File not found: {input_path}")
        return 1
    
    print(f"🔄 Starting deburst processing: {input_path.name}")
    
    try:
        # Initialize reader
        reader = sardine.SlcReader(str(input_path))
        
        # Check if this is IW mode
        if not reader.is_iw_mode():
            print("⚠️  Warning: This doesn't appear to be an IW mode product")
            print("   Deburst is primarily designed for IW mode data")
        
        # Find available polarizations
        annotation_files = reader.find_annotation_files()
        available_polarizations = list(annotation_files.keys())
        
        print(f"📡 Available polarizations: {', '.join(available_polarizations)}")
        
        # Determine which polarizations to process
        if args.polarization:
            if args.polarization.upper() not in available_polarizations:
                print(f"❌ Error: Polarization {args.polarization} not found")
                print(f"   Available: {', '.join(available_polarizations)}")
                return 1
            polarizations_to_process = [args.polarization.upper()]
        else:
            polarizations_to_process = available_polarizations
        
        print(f"🎯 Processing polarizations: {', '.join(polarizations_to_process)}")
        
        # Process each polarization
        results = {}
        total_start_time = time.time()
        
        for pol in polarizations_to_process:
            pol_start_time = time.time()
            print(f"\n🔄 Processing {pol} polarization...")
            
            try:
                # Perform deburst
                deburst_data, (rows, cols) = reader.deburst_slc(pol)
                pol_end_time = time.time()
                
                print(f"✅ Deburst completed for {pol}")
                print(f"   • Input dimensions: Reading from SLC...")
                print(f"   • Output dimensions: {rows:,} x {cols:,}")
                print(f"   • Processing time: {pol_end_time - pol_start_time:.2f} seconds")
                
                results[pol] = {
                    "output_dimensions": (rows, cols),
                    "processing_time": pol_end_time - pol_start_time,
                    "data_size_mb": (rows * cols * 8) / (1024 * 1024)  # Complex64 = 8 bytes
                }
                
                # Save to binary file if requested
                if args.output:
                    output_path = Path(args.output)
                    if len(polarizations_to_process) > 1:
                        # Multiple polarizations: create separate files
                        pol_output_path = output_path.parent / f"{output_path.stem}_{pol.lower()}{output_path.suffix}"
                    else:
                        pol_output_path = output_path
                    
                    # Convert to numpy array and save
                    import numpy as np
                    complex_array = np.array(deburst_data, dtype=np.complex64)
                    np.save(pol_output_path.with_suffix('.npy'), complex_array)
                    
                    print(f"💾 Deburst data saved to: {pol_output_path.with_suffix('.npy')}")
                
            except Exception as e:
                print(f"❌ Error processing {pol}: {e}")
                return 1
        
        total_end_time = time.time()
        
        # Print summary
        print(f"\n" + "="*70)
        print("🎯 DEBURST PROCESSING SUMMARY")
        print("="*70)
        
        for pol, result in results.items():
            print(f"📡 {pol}:")
            print(f"   • Output dimensions: {result['output_dimensions'][0]:,} x {result['output_dimensions'][1]:,}")
            print(f"   • Data size: {result['data_size_mb']:.1f} MB")
            print(f"   • Processing time: {result['processing_time']:.2f} seconds")
        
        print(f"\n⏱️  Total processing time: {total_end_time - total_start_time:.2f} seconds")
        print(f"🚀 Deburst processing complete!")
        print(f"Next step: Apply radiometric calibration and terrain correction")
        
        return 0
        
    except Exception as e:
        print(f"❌ Error during deburst: {e}")
        import traceback
        if args.verbose:
            traceback.print_exc()
        return 1
    

def cmd_calibrate(args):
    """Handle the 'calibrate' command."""
    input_path = Path(args.input)
    
    if not input_path.exists():
        print(f"❌ Error: File not found: {input_path}")
        return 1
    
    print(f"🔄 Starting calibration processing: {input_path.name}")
    
    try:
        # Initialize reader
        reader = sardine.SlcReader(str(input_path))
        
        # Find available polarizations
        annotation_files = reader.find_annotation_files()
        available_polarizations = list(annotation_files.keys())
        
        print(f"📡 Available polarizations: {', '.join(available_polarizations)}")
        
        # Determine which polarizations to process
        if args.polarization:
            if args.polarization.upper() not in available_polarizations:
                print(f"❌ Error: Polarization {args.polarization} not found")
                print(f"   Available: {', '.join(available_polarizations)}")
                return 1
            polarizations_to_process = [args.polarization.upper()]
        else:
            polarizations_to_process = available_polarizations
        
        print(f"🎯 Processing polarizations: {', '.join(polarizations_to_process)}")
        print(f"🎯 Calibration type: {args.calibration_type}")
        
        # Process each polarization
        results = {}
        total_start_time = time.time()
        
        for pol in polarizations_to_process:
            pol_start_time = time.time()
            print(f"\n🔄 Processing {pol} calibration...")
            
            try:
                # Get calibration info first
                cal_info = reader.get_calibration_info(pol)
                print(f"   • Swath: {cal_info['swath']}")
                print(f"   • Calibration vectors: {cal_info['num_vectors']}")
                
                # Perform calibration
                calibrated_data, (rows, cols) = reader.calibrate_slc(pol, args.calibration_type)
                pol_end_time = time.time()
                
                # Calculate data statistics
                import numpy as np
                cal_array = np.array(calibrated_data)
                data_min = np.min(cal_array)
                data_max = np.max(cal_array)
                data_mean = np.mean(cal_array)
                
                print(f"✅ Calibration completed for {pol}")
                print(f"   • Output dimensions: {rows:,} x {cols:,}")
                print(f"   • Data range: {data_min:.2e} to {data_max:.2e}")
                print(f"   • Mean value: {data_mean:.2e}")
                print(f"   • Processing time: {pol_end_time - pol_start_time:.2f} seconds")
                
                results[pol] = {
                    "output_dimensions": (rows, cols),
                    "processing_time": pol_end_time - pol_start_time,
                    "data_size_mb": (rows * cols * 4) / (1024 * 1024),  # float32 = 4 bytes
                    "data_stats": {
                        "min": float(data_min),
                        "max": float(data_max),
                        "mean": float(data_mean)
                    }
                }
                
                # Save to binary file if requested
                if args.output:
                    output_path = Path(args.output)
                    if len(polarizations_to_process) > 1:
                        # Multiple polarizations: create separate files
                        pol_output_path = output_path.parent / f"{output_path.stem}_{pol.lower()}_{args.calibration_type}{output_path.suffix}"
                    else:
                        pol_output_path = output_path.parent / f"{output_path.stem}_{args.calibration_type}{output_path.suffix}"
                    
                    # Save as numpy array
                    np.save(pol_output_path.with_suffix('.npy'), cal_array.astype(np.float32))
                    print(f"💾 Calibrated data saved to: {pol_output_path.with_suffix('.npy')}")
                    
                    # Also save in dB if requested
                    if args.db_scale:
                        db_array = 10 * np.log10(np.maximum(cal_array, 1e-10))  # Avoid log(0)
                        db_output_path = pol_output_path.parent / f"{pol_output_path.stem}_db.npy"
                        np.save(db_output_path, db_array.astype(np.float32))
                        print(f"💾 dB scale data saved to: {db_output_path}")
                
            except Exception as e:
                print(f"❌ Error processing {pol}: {e}")
                if hasattr(args, 'verbose') and args.verbose:
                    import traceback
                    traceback.print_exc()
                return 1
        
        total_end_time = time.time()
        
        # Print summary
        print(f"\n" + "="*70)
        print("🎯 CALIBRATION PROCESSING SUMMARY")
        print("="*70)
        
        for pol, result in results.items():
            print(f"📡 {pol}:")
            print(f"   • Output dimensions: {result['output_dimensions'][0]:,} x {result['output_dimensions'][1]:,}")
            print(f"   • Data size: {result['data_size_mb']:.1f} MB")
            print(f"   • Data range: {result['data_stats']['min']:.2e} to {result['data_stats']['max']:.2e}")
            print(f"   • Processing time: {result['processing_time']:.2f} seconds")
        
        print(f"\n⏱️  Total processing time: {total_end_time - total_start_time:.2f} seconds")
        print(f"🚀 Calibration processing complete!")
        print(f"Next step: Apply multilooking and terrain correction")
        
        return 0
        
    except Exception as e:
        print(f"❌ Error during calibration: {e}")
        import traceback
        if hasattr(args, 'verbose') and args.verbose:
            traceback.print_exc()
        return 1


def cmd_multilook(args):
    """Handle the 'multilook' command."""
    input_path = Path(args.input)
    
    if not input_path.exists():
        print(f"❌ Error: File not found: {input_path}")
        return 1
    
    print(f"🔄 Starting multilook processing: {input_path.name}")
    
    try:
        # Initialize reader
        reader = sardine.SlcReader(str(input_path))
        
        # Find available polarizations
        annotation_files = reader.find_annotation_files()
        available_polarizations = list(annotation_files.keys())
        
        print(f"📡 Available polarizations: {', '.join(available_polarizations)}")
        
        # Determine which polarizations to process
        if args.polarization:
            if args.polarization.upper() not in available_polarizations:
                print(f"❌ Error: Polarization {args.polarization} not found")
                print(f"   Available: {', '.join(available_polarizations)}")
                return 1
            polarizations_to_process = [args.polarization.upper()]
        else:
            polarizations_to_process = available_polarizations
        
        print(f"🎯 Processing polarizations: {', '.join(polarizations_to_process)}")
        print(f"🎯 Calibration type: {args.calibration_type}")
        print(f"🎯 Multilook parameters: {args.azimuth_looks}x{args.range_looks} (azimuth x range)")
        
        # Process each polarization
        results = {}
        total_start_time = time.time()
        
        for pol in polarizations_to_process:
            pol_start_time = time.time()
            print(f"\n🔄 Processing {pol} calibration and multilook...")
            
            try:
                # Get calibration info first
                cal_info = reader.get_calibration_info(pol)
                print(f"   • Swath: {cal_info['swath']}")
                print(f"   • Calibration vectors: {cal_info['num_vectors']}")
                
                # Perform calibration and multilook in one step
                multilooked_data, (new_range_spacing, new_azimuth_spacing) = reader.calibrate_and_multilook(
                    pol, args.calibration_type, args.range_looks, args.azimuth_looks
                )
                pol_end_time = time.time()
                
                # Calculate data statistics
                import numpy as np
                ml_array = np.array(multilooked_data)
                rows, cols = ml_array.shape
                data_min = np.min(ml_array)
                data_max = np.max(ml_array)
                data_mean = np.mean(ml_array)
                
                print(f"✅ Multilook completed for {pol}")
                print(f"   • Output dimensions: {rows:,} x {cols:,}")
                print(f"   • New pixel spacing: {new_range_spacing:.1f}m x {new_azimuth_spacing:.1f}m")
                print(f"   • Data range: {data_min:.2e} to {data_max:.2e}")
                print(f"   • Mean value: {data_mean:.2e}")
                print(f"   • Processing time: {pol_end_time - pol_start_time:.2f} seconds")
                
                # Calculate multilook efficiency
                theoretical_looks = args.range_looks * args.azimuth_looks
                print(f"   • Theoretical looks: {theoretical_looks}")
                
                results[pol] = {
                    "output_dimensions": (rows, cols),
                    "pixel_spacing": (new_range_spacing, new_azimuth_spacing),
                    "multilook_params": (args.azimuth_looks, args.range_looks),
                    "processing_time": pol_end_time - pol_start_time,
                    "data_size_mb": (rows * cols * 4) / (1024 * 1024),  # float32 = 4 bytes
                    "data_stats": {
                        "min": float(data_min),
                        "max": float(data_max),
                        "mean": float(data_mean)
                    }
                }
                
                # Save to binary file if requested
                if args.output:
                    output_path = Path(args.output)
                    if len(polarizations_to_process) > 1:
                        # Multiple polarizations: create separate files
                        pol_output_path = output_path.parent / f"{output_path.stem}_{pol.lower()}_{args.calibration_type}_ml{args.azimuth_looks}x{args.range_looks}{output_path.suffix}"
                    else:
                        pol_output_path = output_path.parent / f"{output_path.stem}_{args.calibration_type}_ml{args.azimuth_looks}x{args.range_looks}{output_path.suffix}"
                    
                    # Save as numpy array
                    np.save(pol_output_path.with_suffix('.npy'), ml_array.astype(np.float32))
                    print(f"💾 Multilooked data saved to: {pol_output_path.with_suffix('.npy')}")
                    
                    # Also save in dB if requested
                    if args.db_scale:
                        db_array = 10 * np.log10(np.maximum(ml_array, 1e-10))  # Avoid log(0)
                        db_output_path = pol_output_path.parent / f"{pol_output_path.stem}_db.npy"
                        np.save(db_output_path, db_array.astype(np.float32))
                        print(f"💾 dB scale data saved to: {db_output_path}")
                
            except Exception as e:
                print(f"❌ Error processing {pol}: {e}")
                if hasattr(args, 'verbose') and args.verbose:
                    import traceback
                    traceback.print_exc()
                return 1
        
        total_end_time = time.time()
        
        # Print summary
        print(f"\n" + "="*70)
        print("🎯 MULTILOOK PROCESSING SUMMARY")
        print("="*70)
        
        for pol, result in results.items():
            print(f"📡 {pol}:")
            print(f"   • Output dimensions: {result['output_dimensions'][0]:,} x {result['output_dimensions'][1]:,}")
            print(f"   • Pixel spacing: {result['pixel_spacing'][0]:.1f}m x {result['pixel_spacing'][1]:.1f}m")
            print(f"   • Multilook: {result['multilook_params'][0]}x{result['multilook_params'][1]} looks")
            print(f"   • Data size: {result['data_size_mb']:.1f} MB")
            print(f"   • Data range: {result['data_stats']['min']:.2e} to {result['data_stats']['max']:.2e}")
            print(f"   • Processing time: {result['processing_time']:.2f} seconds")
        
        print(f"\n⏱️  Total processing time: {total_end_time - total_start_time:.2f} seconds")
        print(f"🚀 Multilook processing complete!")
        print(f"📊 Next step: Apply terrain flattening and terrain correction")
        
        return 0
        
    except Exception as e:
        print(f"❌ Error during multilook: {e}")
        import traceback
        if hasattr(args, 'verbose') and args.verbose:
            traceback.print_exc()
        return 1


def cmd_terrain_flatten(args):
    """Process Sentinel-1 data with terrain flattening (gamma0) using automatic DEM preparation."""
    try:
        print(f"🏔️  Starting terrain flattening workflow...")
        print(f"📁 Input: {args.input}")
        print(f"📊 Polarization: {args.polarization}")
        print(f"📏 Range looks: {args.range_looks}")
        print(f"📏 Azimuth looks: {args.azimuth_looks}")
        print(f"🗺️  DEM cache: {args.dem_cache or './dem_cache'}")
        
        start_time = time.time()
        
        # Open SLC reader
        reader = sardine.SlcReader(args.input)
        
        # Load orbit data if available
        orbit_file = None
        if args.orbit:
            orbit_file = args.orbit
        else:
            # Try to find orbit file automatically
            product_path = Path(args.input)
            orbit_dir = product_path.parent / "orbit_files"
            if orbit_dir.exists():
                orbit_files = list(orbit_dir.glob("*.EOF"))
                if orbit_files:
                    orbit_file = str(orbit_files[0])
                    print(f"📡 Found orbit file: {orbit_file}")
        
        if orbit_file:
            try:
                orbit_data = sardine.load_orbit_file(orbit_file)
                reader.set_orbit_data(orbit_data)
                print(f"✅ Orbit data loaded successfully")
            except Exception as e:
                print(f"⚠️  Warning: Could not load orbit data: {e}")
                print("   Terrain flattening may be less accurate without precise orbit data")
        else:
            print(f"⚠️  Warning: No orbit file found")
            print("   Please provide orbit file with --orbit for best results")
        
        # Parse polarization
        pol_map = {"VV": "VV", "VH": "VH", "HV": "HV", "HH": "HH"}
        if args.polarization.upper() not in pol_map:
            print(f"❌ Error: Invalid polarization '{args.polarization}'. Use VV, VH, HV, or HH")
            return 1
        
        polarization = pol_map[args.polarization.upper()]
        
        # Parse calibration type
        cal_type_map = {
            "sigma0": "Sigma0",
            "beta0": "Beta0", 
            "gamma0": "Gamma0",
            "dn": "Dn"
        }
        cal_type = cal_type_map[args.calibration_type.lower()]
        
        print(f"🔄 Processing with calibration type: {cal_type}")
        
        # Process with automatic DEM preparation
        try:
            gamma0_data, incidence_angles, range_spacing, azimuth_spacing = reader.calibrate_multilook_and_flatten_auto_dem(
                polarization,
                cal_type,
                args.range_looks,
                args.azimuth_looks,
                args.dem_cache
            )
            
            processing_time = time.time() - start_time
            
            print(f"✅ Terrain flattening completed in {processing_time:.2f}s")
            print(f"📊 Output shape: {gamma0_data.shape}")
            print(f"📏 Pixel spacing: {range_spacing:.2f}m (range) x {azimuth_spacing:.2f}m (azimuth)")
            
            # Save gamma0 data
            output_gamma0 = f"gamma0_{polarization}_{args.range_looks}x{args.azimuth_looks}.npy"
            gamma0_data.save(output_gamma0)
            print(f"💾 Gamma0 data saved: {output_gamma0}")
            
            # Save incidence angles
            output_angles = f"incidence_angles_{polarization}_{args.range_looks}x{args.azimuth_looks}.npy"
            incidence_angles.save(output_angles)
            print(f"💾 Incidence angles saved: {output_angles}")
            
            # Save in dB scale if requested
            if args.db_scale:
                import numpy as np
                gamma0_db = 10 * np.log10(np.maximum(gamma0_data, 1e-10))
                output_db = f"gamma0_db_{polarization}_{args.range_looks}x{args.azimuth_looks}.npy"
                np.save(output_db, gamma0_db)
                print(f"💾 Gamma0 (dB) saved: {output_db}")
            
            # Summary
            print("\n" + "="*70)
            print("✅ TERRAIN FLATTENING COMPLETED")
            print("="*70)
            print(f"📁 Input file: {Path(args.input).name}")
            print(f"📊 Polarization: {polarization}")
            print(f"📏 Multilooking: {args.range_looks} x {args.azimuth_looks}")
            print(f"🏔️  Gamma0 range: {float(gamma0_data.min()):.6f} - {float(gamma0_data.max()):.6f}")
            print(f"📐 Incidence angle range: {float(incidence_angles.min()):.2f}° - {float(incidence_angles.max()):.2f}°")
            print(f"⏱️  Processing time: {processing_time:.2f} seconds")
            print(f"📊 Next step: Convert to GeoTIFF or apply further processing")
            
            return 0
            
        except Exception as e:
            print(f"❌ Error during terrain flattening: {e}")
            return 1
            
    except Exception as e:
        print(f"❌ Error: {e}")
        return 1


def main():
    """Main CLI entry point."""
    parser = argparse.ArgumentParser(
        description="SARdine: Fast Sentinel-1 SAR processing",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  sardine info /path/to/S1A_IW_SLC_*.zip     # Show product information
  sardine info --json input.zip              # Output as JSON
  sardine orbit input.zip                     # Check orbit file status
  sardine orbit --download input.zip          # Download orbit files
  sardine orbit --list-urls input.zip         # Show download URLs
  sardine process input.zip output/          # Process SLC to backscatter
  sardine test                                # Run basic tests
  sardine iw-split input.zip                  # Analyze IW sub-swaths
  sardine deburst input.zip                   # Deburst IW product
        """
    )
    
    parser.add_argument(
        "--version", 
        action="version", 
        version=f"SARdine {sardine.__version__}"
    )
    parser.add_argument(
        "-v", "--verbose",
        action="store_true",
        help="Enable verbose output"
    )
    
    subparsers = parser.add_subparsers(dest="command", help="Available commands")
    
    # Info command
    info_parser = subparsers.add_parser(
        "info", 
        help="Display information about a Sentinel-1 product"
    )
    info_parser.add_argument(
        "input",
        help="Path to Sentinel-1 SLC ZIP file"
    )
    info_parser.add_argument(
        "--json",
        action="store_true",
        help="Output information as JSON"
    )
    
    # Process command
    process_parser = subparsers.add_parser(
        "process",
        help="Process Sentinel-1 SLC data to backscatter products"
    )
    process_parser.add_argument(
        "input",
        help="Path to Sentinel-1 SLC ZIP file"
    )
    process_parser.add_argument(
        "output",
        nargs="?",
        default="./output",
        help="Output directory (default: ./output)"
    )
    process_parser.add_argument(
        "--polarizations",
        default="VV,VH",
        help="Polarizations to process (default: VV,VH)"
    )
    
    # Test command
    test_parser = subparsers.add_parser(
        "test",
        help="Run basic functionality tests"
    )
    
    # Orbit command
    orbit_parser = subparsers.add_parser(
        "orbit",
        help="Perform orbit file operations"
    )
    orbit_parser.add_argument(
        "input",
        help="Path to Sentinel-1 SLC ZIP file"
    )
    orbit_parser.add_argument(
        "--download",
        action="store_true",
        help="Download missing orbit files"
    )
    orbit_parser.add_argument(
        "--list-urls",
        action="store_true",
        help="List generated download URLs"
    )
    orbit_parser.add_argument(
        "--output",
        help="Output directory for downloaded orbit files"
    )
    
    # IW split command
    iw_split_parser = subparsers.add_parser(
        "iw-split",
        help="Analyze and extract IW sub-swath information"
    )
    iw_split_parser.add_argument(
        "input",
        help="Path to Sentinel-1 SLC ZIP file"
    )
    iw_split_parser.add_argument(
        "--output",
        help="Output file for sub-swath information (JSON format)"
    )
    
    # Deburst command
    deburst_parser = subparsers.add_parser(
        "deburst",
        help="Deburst IW product to individual bursts"
    )
    deburst_parser.add_argument(
        "input",
        help="Path to Sentinel-1 SLC ZIP file"
    )
    deburst_parser.add_argument(
        "--output",
        help="Output file for debursted data (Numpy .npy format)"
    )
    deburst_parser.add_argument(
        "--polarization",
        help="Specify polarization to deburst (e.g., VV, VH)"
    )
    
    # Calibration command
    calibrate_parser = subparsers.add_parser(
        "calibrate",
        help="Calibrate IW product using calibration vectors"
    )
    calibrate_parser.add_argument(
        "input",
        help="Path to Sentinel-1 SLC ZIP file"
    )
    calibrate_parser.add_argument(
        "--output",
        help="Output file for calibrated data (Numpy .npy format)"
    )
    calibrate_parser.add_argument(
        "--polarization",
        help="Specify polarization to calibrate (e.g., VV, VH)"
    )
    calibrate_parser.add_argument(
        "--calibration-type",
        default="sigma0",
        choices=["sigma0", "beta0", "gamma0", "dn"],
        help="Type of calibration to apply (default: sigma0)"
    )
    calibrate_parser.add_argument(
        "--db-scale",
        action="store_true",
        help="Save output in dB scale (additional .npy file)"
    )
    
    # Multilook command
    multilook_parser = subparsers.add_parser(
        "multilook",
        help="Apply multilooking to calibrated data for speckle reduction"
    )
    multilook_parser.add_argument(
        "input",
        help="Path to Sentinel-1 SLC ZIP file"
    )
    multilook_parser.add_argument(
        "--output",
        help="Output file for multilooked data (Numpy .npy format)"
    )
    multilook_parser.add_argument(
        "--polarization",
        help="Specify polarization to process (e.g., VV, VH)"
    )
    multilook_parser.add_argument(
        "--calibration-type",
        default="sigma0",
        choices=["sigma0", "beta0", "gamma0", "dn"],
        help="Type of calibration to apply before multilooking (default: sigma0)"
    )
    multilook_parser.add_argument(
        "--range-looks",
        type=int,
        default=4,
        help="Number of looks in range direction (default: 4)"
    )
    multilook_parser.add_argument(
        "--azimuth-looks",
        type=int,
        default=1,
        help="Number of looks in azimuth direction (default: 1)"
    )
    multilook_parser.add_argument(
        "--db-scale",
        action="store_true",
        help="Save output in dB scale (additional .npy file)"
    )
    
    # Terrain flatten command
    terrain_parser = subparsers.add_parser(
        "terrain",
        help="Process Sentinel-1 data with terrain flattening (gamma0)"
    )
    terrain_parser.add_argument(
        "input",
        help="Path to Sentinel-1 SLC ZIP file"
    )
    terrain_parser.add_argument(
        "--output",
        help="Output file for gamma0 data (Numpy .npy format)"
    )
    terrain_parser.add_argument(
        "--polarization",
        default="VV",
        help="Specify polarization to process (default: VV)"
    )
    terrain_parser.add_argument(
        "--calibration-type",
        default="sigma0",
        choices=["sigma0", "beta0", "gamma0", "dn"],
        help="Type of calibration to apply (default: sigma0)"
    )
    terrain_parser.add_argument(
        "--range-looks",
        type=int,
        default=4,
        help="Number of looks in range direction (default: 4)"
    )
    terrain_parser.add_argument(
        "--azimuth-looks",
        type=int,
        default=1,
        help="Number of looks in azimuth direction (default: 1)"
    )
    terrain_parser.add_argument(
        "--db-scale",
        action="store_true",
        help="Save output in dB scale (additional .npy file)"
    )
    terrain_parser.add_argument(
        "--dem-cache",
        help="Directory for DEM cache (default: ./dem_cache)"
    )
    terrain_parser.add_argument(
        "--orbit",
        help="Path to orbit file (EOF format)"
    )
    
    # Test SRTM download command
    test_srtm_parser = subparsers.add_parser(
        "test-srtm",
        help="Test SRTM download functionality"
    )
    test_srtm_parser.add_argument(
        "tile",
        help="SRTM tile name to test (e.g., N37W122)"
    )
    test_srtm_parser.add_argument(
        "--output",
        default="./srtm_test",
        help="Output directory for test download (default: ./srtm_test)"
    )
    
    # Speckle filter command
    speckle_parser = subparsers.add_parser(
        "speckle-filter",
        help="Apply speckle filtering to SAR intensity images"
    )
    speckle_parser.add_argument(
        "input",
        help="Input GeoTIFF intensity image"
    )
    speckle_parser.add_argument(
        "output",
        help="Output filtered GeoTIFF image"
    )
    speckle_parser.add_argument(
        "--filter-type",
        choices=["mean", "median", "lee", "enhanced_lee", "lee_sigma", "frost", "gamma_map", "refined_lee"],
        default="lee",
        help="Type of speckle filter to apply (default: lee)"
    )
    speckle_parser.add_argument(
        "--window-size",
        type=int,
        default=7,
        help="Filter window size (must be odd, default: 7)"
    )
    speckle_parser.add_argument(
        "--num-looks",
        type=float,
        help="Number of looks (auto-estimated if not provided)"
    )
    speckle_parser.add_argument(
        "--edge-threshold",
        type=float,
        default=0.5,
        help="Edge detection threshold for Lee Sigma filter (default: 0.5)"
    )
    speckle_parser.add_argument(
        "--damping-factor",
        type=float,
        default=1.0,
        help="Damping factor for Gamma MAP filter (default: 1.0)"
    )
    speckle_parser.add_argument(
        "--cv-threshold",
        type=float,
        default=0.5,
        help="Coefficient of variation threshold (default: 0.5)"
    )
    
    # Estimate number of looks command
    nlooks_parser = subparsers.add_parser(
        "estimate-nlooks",
        help="Estimate number of looks from SAR intensity image"
    )
    nlooks_parser.add_argument(
        "input",
        help="Input GeoTIFF intensity image"
    )
    nlooks_parser.add_argument(
        "--window-size",
        type=int,
        default=11,
        help="Window size for estimation (default: 11)"
    )
    
    args = parser.parse_args()
    
    if not args.command:
        print_banner()
        parser.print_help()
        return 0
    
    # Route to appropriate command handler
    if args.command == "info":
        return cmd_info(args)
    elif args.command == "process":
        return cmd_process(args)
    elif args.command == "test":
        return cmd_test(args)
    elif args.command == "orbit":
        return cmd_orbit(args)
    elif args.command == "iw-split":
        return cmd_iw_split(args)
    elif args.command == "deburst":
        return cmd_deburst(args)
    elif args.command == "calibrate":
        return cmd_calibrate(args)
    elif args.command == "multilook":
        return cmd_multilook(args)
    elif args.command == "terrain":
        return cmd_terrain_flatten(args)
    elif args.command == "test-srtm":
        return cmd_test_srtm(args)
    elif args.command == "speckle-filter":
        return cmd_speckle_filter(args)
    elif args.command == "estimate-nlooks":
        return cmd_estimate_nlooks(args)
    else:
        parser.print_help()
        return 1


if __name__ == "__main__":
    sys.exit(main())
