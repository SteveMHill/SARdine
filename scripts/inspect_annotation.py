#!/usr/bin/env python3
"""
Script to inspect Sentinel-1 annotation XML files for burst information
"""
import sys
import os
import zipfile
import xml.etree.ElementTree as ET

def count_burst_elements(xml_content):
    """Count burst elements in XML content"""
    root = ET.fromstring(xml_content)
    
    # Look for burst elements
    burst_elements = root.findall(".//burst")
    if burst_elements:
        return len(burst_elements)
    
    # Try with namespace
    namespaces = {'s1': 'http://www.esa.int/safe/sentinel-1/1.1'}
    burst_elements = root.findall(".//s1:burst", namespaces)
    
    return len(burst_elements)

def extract_burst_info(xml_content):
    """Extract basic burst information"""
    root = ET.fromstring(xml_content)
    
    # Try to find burstList element
    burst_list = root.find(".//burstList")
    if burst_list is None:
        # Try with namespace
        namespaces = {'s1': 'http://www.esa.int/safe/sentinel-1/1.1'}
        burst_list = root.find(".//s1:burstList", namespaces)
    
    if burst_list is None:
        print("No burstList element found")
        return []
    
    # Get burst elements
    bursts = burst_list.findall("./burst") or burst_list.findall("./s1:burst", {'s1': 'http://www.esa.int/safe/sentinel-1/1.1'})
    
    burst_info = []
    for i, burst in enumerate(bursts):
        info = {'index': i}
        
        # Extract basic properties
        for prop in ['azimuthTime', 'sensingTime', 'byteOffset']:
            elem = burst.find(f"./{prop}") or burst.find(f"./s1:{prop}", {'s1': 'http://www.esa.int/safe/sentinel-1/1.1'})
            if elem is not None and elem.text:
                info[prop] = elem.text.strip()
            else:
                info[prop] = None
        
        # Extract first/last valid sample
        for prop in ['firstValidSample', 'lastValidSample']:
            elem = burst.find(f"./{prop}") or burst.find(f"./s1:{prop}", {'s1': 'http://www.esa.int/safe/sentinel-1/1.1'})
            if elem is not None and elem.text:
                # Just count the values
                values = elem.text.strip().split()
                info[f"{prop}_count"] = len(values)
                # Get first few samples for reference
                info[f"{prop}_examples"] = values[:5] if len(values) > 0 else []
            else:
                info[f"{prop}_count"] = 0
                info[f"{prop}_examples"] = []
        
        burst_info.append(info)
    
    return burst_info

def main():
    if len(sys.argv) < 2:
        print(f"Usage: {sys.argv[0]} <path_to_S1_zip>")
        sys.exit(1)
    
    slc_path = sys.argv[1]
    if not os.path.exists(slc_path):
        print(f"File not found: {slc_path}")
        sys.exit(1)
    
    print(f"Inspecting: {os.path.basename(slc_path)}")
    
    try:
        with zipfile.ZipFile(slc_path, 'r') as z:
            # Find annotation files
            annotation_files = [f for f in z.namelist() 
                               if 'annotation/' in f and f.endswith('.xml') 
                               and 'calibration' not in f and 'noise' not in f]
            
            print(f"Found {len(annotation_files)} annotation files")
            
            for ann_file in annotation_files:
                print("\n" + "="*60)
                print(f"Inspecting: {ann_file}")
                print("="*60)
                
                try:
                    with z.open(ann_file) as xml_file:
                        content = xml_file.read().decode('utf-8')
                        
                        # Count burst elements
                        burst_count = count_burst_elements(content)
                        print(f"Found {burst_count} burst elements")
                        
                        # Extract burst information
                        burst_info = extract_burst_info(content)
                        
                        print(f"Extracted information for {len(burst_info)} bursts:")
                        
                        for info in burst_info:
                            print(f"\nBurst #{info['index']}:")
                            print(f"  Azimuth time: {info.get('azimuthTime')}")
                            print(f"  Sensing time: {info.get('sensingTime')}")
                            print(f"  Byte offset: {info.get('byteOffset')}")
                            print(f"  First valid sample count: {info.get('firstValidSample_count')}")
                            print(f"  First valid sample examples: {info.get('firstValidSample_examples')}")
                            print(f"  Last valid sample count: {info.get('lastValidSample_count')}")
                            print(f"  Last valid sample examples: {info.get('lastValidSample_examples')}")
                            
                except Exception as e:
                    print(f"Error processing {ann_file}: {e}")
                    
    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()
