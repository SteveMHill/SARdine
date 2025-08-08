#!/usr/bin/env python3
"""
Debug SlcReader Methods
========================

Check what methods are available on SlcReader to understand how to read annotation data.
"""

import sys
from pathlib import Path

# Add the project directory to the path
project_root = Path(__file__).parent
sys.path.insert(0, str(project_root))

import sardine

def check_slc_reader_methods():
    """Check available methods on SlcReader"""
    
    slc_path = "/home/datacube/apps/SARdine/data/S1A_IW_SLC__1SDV_20200103T170815_20200103T170842_030639_0382D5_DADE.zip"
    
    try:
        reader = sardine.SlcReader(slc_path)
        
        print("🔍 Available methods on SlcReader:")
        methods = [method for method in dir(reader) if not method.startswith('_')]
        
        for i, method in enumerate(methods, 1):
            print(f"   {i:2d}. {method}")
            
        print(f"\n📊 Total methods: {len(methods)}")
        
        # Try to call some methods to see what they return
        print(f"\n🧪 Testing some methods:")
        
        try:
            files = reader.list_files()
            print(f"   • list_files(): {len(files)} files found")
        except Exception as e:
            print(f"   • list_files(): Failed - {e}")
            
        try:
            annotations = reader.find_annotation_files()
            print(f"   • find_annotation_files(): {len(annotations)} annotations found")
            for pol, file in annotations.items():
                print(f"     - {pol}: {file}")
        except Exception as e:
            print(f"   • find_annotation_files(): Failed - {e}")
            
    except Exception as e:
        print(f"❌ Error: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    check_slc_reader_methods()
