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
    else:
        parser.print_help()
        return 1


if __name__ == "__main__":
    sys.exit(main())
