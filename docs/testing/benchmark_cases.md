# Benchmark Cases

## Purpose
This file tracks Sentinel-1 scenes that may later be used to diagnose, compare, and verify the pipeline.

## Current status
A complete benchmark suite does not yet exist.

## Rule
Do not invent trusted benchmark cases.
If a case is incomplete, mark it as incomplete.

---

## Candidate Case 001

### Case ID
- ID: case_001
- Short name: duplicated bursts and missing area in VV backscatter output

### Current availability
- Available locally: yes
- Raw input still accessible: yes
- Orbit file available: unknown
- DEM choice known: unknown

### What is known
- Sentinel-1 product: `S1A_IW_SLC__1SDV_20201005T170824_20201005T170851_034664_04098A_1E66.SAFE`
- Product type: SLC
- Acquisition mode: IW
- Polarization: VV
- Calibration type: sigma0
- Output unit: dB
- AOI / crop: unknown
- Output CRS: EPSG:32632
- Notes: processed VV backscatter output exists and shows duplicated bursts plus missing scene area

### Output metadata
- Shape: 10078 x 12263
- Valid pixels: 91,113,091
- Total pixels: 123,586,514
- Valid percentage: 73.72413708505444
- Data min: -15.885869026184082
- Data max: -0.8239160180091858
- Data mean: -9.08894157409668
- Data median: -9.139272689819336
- Range spacing (m): 11.64781
- Azimuth spacing (m): 13.89642
- Multilook range: 5.0
- Multilook azimuth: 1.0

### Why this case may matter
- Suspected failure symptom: duplicated bursts and missing part of the final scene
- Where the symptom was observed: processed VV backscatter product
- Confidence in recollection: high
- Diagnostic value: high, because this case directly exposes one of the main known failure classes

### What should later be inspected
- Whether duplication appears before or after geocoding
- Whether the missing area originates during deburst / merge, masking, or terrain correction
- Whether burst overlap handling is causing repeated spatial coverage
- Whether invalid mask propagation removes valid data incorrectly
- Whether target grid definition or resampling causes dropped coverage

### Missing information
- Missing files: exact orbit file unknown, DEM source unknown
- Missing metadata: exact processing configuration, intermediate step outputs
- Missing reproduction notes: exact stage where the defect first becomes visible

### Status
- Reproducible now: partly
- Suitable for verification now: yes, as a forensic baseline
- Notes: this case should be used early for repository forensics and later for differential debugging

### Visual / symptom notes
- Duplicated bursts appear: repeated narrow upper strips / duplicated-looking swath fragments near the upper part of the output
- Missing area appears as: several long band-like blank gaps near the upper portion of the scene
- Earliest stage where artifact is known: final geocoded backscatter output
- Screenshot available: yes

### Screenshot interpretation
- The output footprint is rotated in projected space
- Artifact pattern is structured, not random noise
- Blank regions look like coverage loss or masking gaps, not simple display problems
- Upper scene area suggests possible burst overlap / merge / masking failure