#!/usr/bin/env python3
"""
Debug TOPSAR Burst Information Extraction
==========================================

This script examines the annotation XML structure to understand
why burst information extraction is failing and fix the regex patterns.
"""

import sys
import os
from pathlib import Path

# Add the project directory to the path
project_root = Path(__file__).parent
sys.path.insert(0, str(project_root))

try:
    import sardine
    print("✅ SARdine library imported successfully")
except ImportError as e:
    print(f"❌ Failed to import sardine: {e}")
    sys.exit(1)

def debug_annotation_structure():
    """Debug annotation XML to understand burst structure"""
    
    slc_path = "/home/datacube/apps/SARdine/data/S1A_IW_SLC__1SDV_20200103T170815_20200103T170842_030639_0382D5_DADE.zip"
    
    print(f"🔍 Debugging annotation structure for: {Path(slc_path).name}")
    
    try:
        # Create SLC reader
        reader = sardine.SlcReader(slc_path)
        
        # Try to read metadata which should show annotation structure
        print("\n📊 Extracting metadata to examine annotation structure...")
        metadata = reader.read_metadata("VV")
        
        print("✅ Metadata extraction successful")
        print(f"📈 Available fields: {len(metadata)} total")
        
        # Look for burst-related fields
        burst_fields = [key for key in metadata.keys() if 'burst' in key.lower()]
        if burst_fields:
            print(f"\n🎯 Burst-related fields found:")
            for field in burst_fields:
                print(f"   • {field}: {metadata[field]}")
        else:
            print(f"\n❌ No burst-related fields found in metadata")
            
        # Look for timing fields that might indicate burst structure
        timing_fields = [key for key in metadata.keys() if any(term in key.lower() for term in ['time', 'sensing', 'azimuth'])]
        if timing_fields:
            print(f"\n⏰ Timing-related fields found:")
            for field in timing_fields[:10]:  # Show first 10 to avoid clutter
                print(f"   • {field}: {metadata[field]}")
                
        return True
        
    except Exception as e:
        print(f"❌ Error during debugging: {e}")
        import traceback
        traceback.print_exc()
        return False

if __name__ == "__main__":
    success = debug_annotation_structure()
    sys.exit(0 if success else 1)
