#!/usr/bin/env python3
"""
SARdine Command Line Interface

A command-line tool for processing Sentinel-1 SLC data.
"""

import argparse
import sys
import json
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


def main():
    """Main CLI entry point."""
    parser = argparse.ArgumentParser(
        description="SARdine: Fast Sentinel-1 SAR processing",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  sardine info /path/to/S1A_IW_SLC_*.zip     # Show product information
  sardine info --json input.zip              # Output as JSON
  sardine process input.zip output/          # Process SLC to backscatter
  sardine test                                # Run basic tests
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
    else:
        parser.print_help()
        return 1


if __name__ == "__main__":
    sys.exit(main())
