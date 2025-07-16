#!/usr/bin/env python3
"""
Quick inspection of calibration files in the Sentinel-1 dataset.
"""

import zipfile
import xml.etree.ElementTree as ET
from pathlib import Path

def inspect_calibration_files():
    """Inspect calibration files in the SLC dataset."""
    
    slc_path = "/home/datacube/SARdine/data/S1A_IW_SLC__1SDV_20200103T170815_20200103T170842_030639_0382D5_DADE.zip"
    
    if not Path(slc_path).exists():
        print(f"❌ Dataset not found: {slc_path}")
        return False
    
    print("🔍 Inspecting Calibration Files")
    print("=" * 50)
    
    try:
        with zipfile.ZipFile(slc_path, 'r') as archive:
            # List all files in the archive
            all_files = archive.namelist()
            
            # Find calibration files
            calibration_files = []
            for filename in all_files:
                if 'annotation/calibration/' in filename and filename.endswith('.xml'):
                    calibration_files.append(filename)
            
            print(f"📋 Total files in archive: {len(all_files)}")
            print(f"📊 Calibration files found: {len(calibration_files)}")
            
            if calibration_files:
                print("\n📄 Calibration files:")
                for i, cal_file in enumerate(calibration_files, 1):
                    print(f"   {i}. {cal_file}")
                
                # Examine the first calibration file
                first_cal_file = calibration_files[0]
                print(f"\n🔍 Examining: {first_cal_file}")
                
                with archive.open(first_cal_file, 'r') as f:
                    content = f.read().decode('utf-8')
                    
                    # Parse XML
                    root = ET.fromstring(content)
                    print(f"   • Root element: {root.tag}")
                    
                    # Look for key elements
                    swath = root.find('.//swath')
                    pol = root.find('.//polarisation')
                    
                    if swath is not None:
                        print(f"   • Swath: {swath.text}")
                    if pol is not None:
                        print(f"   • Polarization: {pol.text}")
                    
                    # Count calibration vectors
                    cal_vectors = root.findall('.//calibrationVector')
                    print(f"   • Calibration vectors: {len(cal_vectors)}")
                    
                    if cal_vectors:
                        # Examine first vector
                        first_vector = cal_vectors[0]
                        print(f"   • First vector elements:")
                        
                        for child in first_vector:
                            if child.tag == 'azimuthTime':
                                print(f"     - {child.tag}: {child.text}")
                            elif child.tag == 'line':
                                print(f"     - {child.tag}: {child.text}")
                            elif child.tag == 'pixel':
                                values = child.text.strip() if child.text else ""
                                pixel_count = len(values.split()) if values else 0
                                print(f"     - {child.tag}: {pixel_count} pixel values")
                            elif child.tag in ['sigmaNought', 'betaNought', 'gamma', 'dn']:
                                values = child.text.strip() if child.text else ""
                                value_count = len(values.split()) if values else 0
                                print(f"     - {child.tag}: {value_count} calibration values")
                                if value_count > 0:
                                    first_few = " ".join(values.split()[:3])
                                    print(f"       First 3: {first_few}")
                
                print(f"\n✅ Calibration files are present and properly structured!")
                return True
            else:
                print("❌ No calibration files found in the dataset!")
                
                # Let's see what annotation files are available
                annotation_files = [f for f in all_files if 'annotation/' in f and f.endswith('.xml')]
                print(f"\n📄 Available annotation files ({len(annotation_files)}):")
                for ann_file in annotation_files[:10]:  # Show first 10
                    print(f"   • {ann_file}")
                if len(annotation_files) > 10:
                    print(f"   ... and {len(annotation_files) - 10} more")
                
                return False
            
    except Exception as e:
        print(f"❌ Error inspecting calibration files: {e}")
        import traceback
        traceback.print_exc()
        return False

if __name__ == "__main__":
    inspect_calibration_files()
