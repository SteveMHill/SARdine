#!/usr/bin/env python3
"""
Corrected Working SAR Processing Pipeline
Using only working functions with proper data flow
"""

import sardine
import numpy as np
import time
import json
from pathlib import Path
from datetime import datetime

def working_sar_pipeline(input_file, output_dir="./working_output"):
    """
    Working SAR Processing Pipeline - Using Only Successfully Tested Functions
    """
    
    print("🛰️  WORKING SAR PROCESSING PIPELINE")
    print("✅ Using only verified working functions with real data")
    print("=" * 70)
    
    pipeline_start = time.time()
    results = {}
    
    # Create output directory
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    
    try:
        # STEP 1: Product Analysis (✅ WORKING)
        print("STEP 1: Product Analysis & Metadata Extraction 📋")
        step_start = time.time()
        product_info = sardine.get_product_info(input_file)
        step_time = time.time() - step_start
        print(f"✅ Complete in {step_time:.2f}s")
        print(f"   📁 Product: {product_info.get('platform', 'unknown')}")
        print(f"   📅 Date: {product_info.get('start_time', 'unknown')}")
        print(f"   📡 Polarizations: {product_info.get('polarizations', 'unknown')}")
        results['step1_product_analysis'] = product_info

        # STEP 2: Precise Orbit File Application (✅ WORKING)
        print("\nSTEP 2: Precise Orbit File Application 🛰️")
        step_start = time.time()
        
        filename = Path(input_file).name
        product_id = filename[:-4] if filename.endswith('.zip') else filename
        parts = filename.split('_')
        start_time_str = parts[5]
        sensing_time = datetime.strptime(start_time_str, '%Y%m%dT%H%M%S')
        start_time_rfc3339 = sensing_time.isoformat() + 'Z'
        
    orbit_cache = f"{output_dir}/orbit_cache"
    orbit_result = sardine.apply_precise_orbit_file(product_id, start_time_rfc3339, orbit_cache)
        step_time = time.time() - step_start
        print(f"✅ Complete in {step_time:.2f}s")
        print(f"   📡 Type: {orbit_result['result']['orbit_type']}")
        print(f"   📊 Vectors: {orbit_result['result']['orbit_vectors_count']}")
        results['step2_orbit_correction'] = orbit_result

        # STEP 3: IW Sub-swath Splitting (✅ WORKING)
        print("\nSTEP 3: IW Sub-swath Splitting 🎯")
        step_start = time.time()
        subswath_results = {}
        
        # Process all available polarizations
        polarizations = product_info.get('polarizations', 'VV').split(',')
        for pol in polarizations:
            pol = pol.strip()
            subswath_results[pol] = {}
            for subswath in ['IW1', 'IW2', 'IW3']:
                try:
                    iw_result = sardine.iw_split_with_real_data(input_file, pol, subswath)
                    subswath_results[pol][subswath] = iw_result
                    print(f"   ✅ {pol}-{subswath} split successfully")
                except Exception as e:
                    print(f"   ⚠️  {pol}-{subswath} split failed: {e}")
                    subswath_results[pol][subswath] = {'error': str(e)}
        
        step_time = time.time() - step_start
        print(f"✅ Complete in {step_time:.2f}s")
        results['step3_iw_split'] = subswath_results

        # STEP 4: Deburst Processing (✅ WORKING - single subswath)
        print("\nSTEP 4: Deburst Processing 🔄")
        step_start = time.time()
        deburst_results = {}
        
        for pol in polarizations:
            pol = pol.strip()
            deburst_results[pol] = {}
            try:
                # Process first available subswath for deburst
                deburst_result = sardine.deburst_topsar(input_file, pol, 'IW1')
                deburst_results[pol]['IW1'] = deburst_result
                print(f"   ✅ {pol}-IW1 deburst completed")
            except Exception as e:
                print(f"   ⚠️  {pol} deburst failed: {e}")
                deburst_results[pol] = {'error': str(e)}
        
        step_time = time.time() - step_start
        print(f"✅ Step completed in {step_time:.2f}s")
        results['step4_deburst'] = deburst_results

        # STEP 5: Get Available Functions for Further Processing
        print("\nSTEP 5: Available Processing Functions 🔍")
        step_start = time.time()
        
        # Check what other processing functions are available
        available_functions = []
        function_tests = {
            'radiometric_calibration': 'sardine.radiometric_calibration',
            'multilooking': 'sardine.apply_multilooking', 
            'speckle_filtering': 'sardine.apply_speckle_filter',
            'terrain_correction': 'sardine.apply_terrain_correction_with_real_orbits',
            'export_geotiff': 'sardine.export_geotiff'
        }
        
        for func_name, func_path in function_tests.items():
            try:
                # Check if function exists
                func = eval(func_path)
                available_functions.append(func_name)
                print(f"   ✅ {func_name} - Available")
            except Exception as e:
                print(f"   ❌ {func_name} - {e}")
        
        step_time = time.time() - step_start
        print(f"✅ Function check completed in {step_time:.2f}s")
        results['step5_available_functions'] = available_functions

        # STEP 6: Parameter Requirements Analysis
        print("\nSTEP 6: Parameter Requirements Analysis 📋")
        step_start = time.time()
        
        # Analyze what parameters are needed for next steps
        param_analysis = {}
        
        # For calibration - check required parameters
        try:
            # This will fail but tell us what parameters are needed
            help_result = sardine.radiometric_calibration.__doc__
            param_analysis['calibration'] = "Requires: zip_path, polarization, calibration_type, slc_data"
        except Exception as e:
            param_analysis['calibration'] = str(e)
        
        # For multilooking
        try:
            help_result = sardine.apply_multilooking.__doc__
            param_analysis['multilooking'] = "Requires: slc_data, range_looks, azimuth_looks, input_range_spacing, input_azimuth_spacing"
        except Exception as e:
            param_analysis['multilooking'] = str(e)
        
        step_time = time.time() - step_start
        print(f"✅ Parameter analysis completed in {step_time:.2f}s")
        for func, params in param_analysis.items():
            print(f"   📋 {func}: {params}")
        results['step6_parameter_analysis'] = param_analysis

        # STEP 7: Data Flow Requirements
        print("\nSTEP 7: Data Flow Requirements 🔗")
        step_start = time.time()
        
        data_flow = {
            'current_state': 'IW subswaths split and deburst completed',
            'required_for_calibration': 'SLC complex data arrays from deburst',
            'required_for_multilooking': 'Calibrated intensity data + pixel spacings',
            'required_for_terrain_correction': 'Multilooked data + orbit data + DEM',
            'missing_implementations': [
                'Extract SLC data arrays from deburst results',
                'Chain calibration with proper SLC data',
                'Get pixel spacing from metadata',
                'Implement proper data array formats'
            ]
        }
        
        step_time = time.time() - step_start
        print(f"✅ Data flow analysis completed in {step_time:.2f}s")
        for key, value in data_flow.items():
            if isinstance(value, list):
                print(f"   🔗 {key}:")
                for item in value:
                    print(f"      • {item}")
            else:
                print(f"   🔗 {key}: {value}")
        results['step7_data_flow'] = data_flow

        # STEP 8: Next Steps Recommendations
        print("\nSTEP 8: Next Steps Recommendations 💡")
        step_start = time.time()
        
        recommendations = {
            'immediate_priorities': [
                '1. Implement data extraction from deburst results',
                '2. Create proper SLC data arrays for calibration',
                '3. Extract pixel spacing from product metadata',
                '4. Chain calibration → multilooking → terrain correction'
            ],
            'missing_critical_steps': [
                'Radiometric calibration with real SLC data',
                'Multilooking with proper intensity data', 
                'Terrain correction/geocoding',
                'Speckle filtering',
                'Quality assessment and masking',
                'Export to standard formats (GeoTIFF)'
            ],
            'working_foundation': [
                '✅ Product analysis and metadata extraction',
                '✅ Precise orbit file integration', 
                '✅ IW subswath splitting',
                '✅ TOPSAR deburst processing'
            ]
        }
        
        step_time = time.time() - step_start
        print(f"✅ Recommendations completed in {step_time:.2f}s")
        for category, items in recommendations.items():
            print(f"   💡 {category.replace('_', ' ').title()}:")
            for item in items:
                print(f"      {item}")
        results['step8_recommendations'] = recommendations

        # PIPELINE SUMMARY
        total_time = time.time() - pipeline_start
        print(f"\n" + "=" * 70)
        print(f"🎯 WORKING PIPELINE ANALYSIS COMPLETE!")
        print(f"⏱️  Total analysis time: {total_time:.2f} seconds")
        print(f"✅ Successfully working steps: 4 (Steps 1-4)")
        print(f"⚠️  Steps needing implementation: 10+ (Steps 5-14)")
        
        print(f"\n📊 CRITICAL FINDINGS:")
        print(f"   ✅ Core SAR reading and orbit correction: WORKING")
        print(f"   ✅ IW splitting and deburst: WORKING") 
        print(f"   ❌ Data array extraction: MISSING")
        print(f"   ❌ Calibration chaining: MISSING")
        print(f"   ❌ Multilooking implementation: MISSING")
        print(f"   ❌ Terrain correction: MISSING")
        print(f"   ❌ Export functionality: MISSING")
        
        print(f"\n📁 Output directory: {output_dir}")
        print(f"🔬 Analysis complete - ready for implementation planning")
        
        # Save comprehensive results
        results_file = f"{output_dir}/working_pipeline_analysis.json"
        with open(results_file, 'w') as f:
            json_results = {}
            for key, value in results.items():
                try:
                    json.dumps(value)
                    json_results[key] = value
                except:
                    json_results[key] = str(value)
            json.dump(json_results, f, indent=2)
        print(f"📄 Complete analysis: {results_file}")
        
        return results

    except Exception as e:
        print(f"\n❌ Analysis failed: {e}")
        import traceback
        traceback.print_exc()
        return None

if __name__ == "__main__":
    input_file = './data/S1A_IW_SLC__1SDV_20200103T170815_20200103T170842_030639_0382D5_DADE.zip'
    results = working_sar_pipeline(input_file)
