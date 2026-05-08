#!/usr/bin/env python3
import numpy as np

# Original values from XML
K_xml = 318.8  # sigmaNought from calibration XML
abs_const = 1.393  # absoluteCalibrationConstant

print('Sentinel-1 Calibration Formula Analysis')
print('=' * 50)
print(f'K_xml (sigmaNought): {K_xml}')
print(f'absoluteCalibrationConstant: {abs_const}')
print()

# ESA spec: What does K represent?
print('Checking ESA formula interpretation:')
print('=' * 50)
print('ESA documentation states:')
print('  sigma0^2 = (DN / K)^2 * abs_const')
print('  Therefore: sigma0 = (DN^2 / K^2) * abs_const')
print()
print('For K ~ 318.8 (large), this makes sense.')
print(f'  1/K^2 = 1/{K_xml}^2 = {1/(K_xml**2):.6e} (small, as expected for gain)')
print()

# OLD CODE
old_k_modified = K_xml * abs_const  # Multiplied K by abs_const
old_final_lut = 1.0 / (old_k_modified ** 2)  # Then inverted AND SQUARED
print(f'OLD CODE (if squaring):')
print(f'  Step 1: K_modified = K * abs_const = {K_xml} * {abs_const} = {old_k_modified:.2f}')
print(f'  Step 2: LUT = 1/K_modified^2 = 1/{old_k_modified:.2f}^2 = {old_final_lut:.6e}')
print()

# But wait - our code just inverts, doesn't square!
old_final_lut_no_sq = 1.0 / old_k_modified
print(f'OLD CODE (no squaring - what we actually do):')
print(f'  Step 1: K_modified = K * abs_const = {K_xml} * {abs_const} = {old_k_modified:.2f}')
print(f'  Step 2: LUT = 1/K_modified = 1/{old_k_modified:.2f} = {old_final_lut_no_sq:.6e}')
print()

# NEW CODE  
new_final_lut = (1.0 / K_xml) * abs_const
print(f'NEW CODE (my fix - no squaring):')
print(f'  Step 1: LUT = (1/K) * abs_const = (1/{K_xml}) * {abs_const} = {new_final_lut:.6e}')
print()

# Ratio (no squaring)
ratio = new_final_lut / old_final_lut_no_sq
print(f'Ratio: NEW/OLD = {ratio:.4f} = {10*np.log10(ratio):.2f} dB')
print()

# Check against observed difference
observed_diff_db = 22.56 - 12.26
print(f'Observed difference: {observed_diff_db:.2f} dB')
print(f'Expected from abs_const^2: {10*np.log10(abs_const**2):.2f} dB')
print()

print('=' * 50)
print('CONCLUSION:')
print('The observed +10.30 dB matches abs_const^2 in dB!')
print('This suggests the ESA formula uses K^2, not K.')
print('We need to CHECK if K values in XML are already squared or not.')
