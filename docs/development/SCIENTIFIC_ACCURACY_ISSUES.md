# Scientific Accuracy Issues

## Critical Scientific Problems Found During Audit

### 1. Hardcoded Values Instead of Real Data Extraction

#### Issue: Hardcoded Pixel Spacing  
**Status**: ❌ SCIENTIFIC ERROR
**Location**: Multiple locations using hardcoded 10.0m spacing
**Problem**: Real Sentinel-1 pixel spacing varies by mode and subswath:
- **IW1**: Range ~2.7m, Azimuth ~22m  
- **IW2**: Range ~3.1m, Azimuth ~22m
- **IW3**: Range ~3.5m, Azimuth ~22m

**Impact**: Wrong pixel spacing leads to:
- Incorrect multilooking ratios
- Wrong output resolution calculations  
- Invalid geospatial coordinates
- Scientifically meaningless results

#### Issue: Hardcoded Scene Bounding Box
**Status**: ❌ SCIENTIFIC ERROR  
**Location**: Line 288 in backscatter.py
```python
sar_bbox = [10.0, 45.0, 11.0, 46.0]  # Northern Italy region
```

**Problem**: This hardcoded bbox only works for data over Northern Italy
**Impact**: 
- Cannot process data from other global regions
- Wrong DEM data downloaded
- Invalid terrain correction for non-Italy scenes
- Geographic coordinates completely wrong

#### Issue: Hardcoded GeoTransform
**Status**: ❌ SCIENTIFIC ERROR
**Location**: Line 371 in backscatter.py  
```python
geo_transform = [10.0, 0.0001, 0.0, 46.0, 0.0, -0.0001]
```

**Problem**: Fixed coordinates for Northern Italy (46°N, 10°E)
**Impact**:
- GeoTIFF files have wrong geographic coordinates
- Cannot be used in GIS systems correctly
- Invalid for any scene outside Northern Italy

### 2. Missing Scientific Methodologies

#### Issue: No Doppler Centroid Extraction
**Status**: ❌ MISSING CRITICAL PARAMETER
**Location**: Deburst processing
**Problem**: Uses hardcoded doppler centroid (150.0 Hz)
**Scientific Requirement**: Must extract from annotation XML per burst
**Impact**: 
- Incorrect burst timing
- Invalid azimuth focusing  
- Wrong geometric accuracy

#### Issue: No Real FM Rate Extraction  
**Status**: ❌ MISSING CRITICAL PARAMETER
**Location**: Deburst processing  
**Problem**: Uses hardcoded FM rate (-2300.0 Hz/s)
**Scientific Requirement**: Must extract from annotation XML per burst
**Impact**:
- Incorrect range cell migration correction
- Invalid focusing quality
- Wrong geometric fidelity

#### Issue: No Incidence Angle Calculation
**Status**: ❌ MISSING SCIENTIFIC PARAMETER
**Location**: Terrain flattening
**Problem**: No local incidence angle computation from DEM
**Scientific Requirement**: Essential for terrain flattening (γ⁰ calculation)
**Impact**:
- Cannot perform proper terrain flattening
- Invalid γ⁰ (gamma nought) values
- No topographic normalization

### 3. Invalid Calibration Methodology

#### Issue: Hardcoded Calibration Scaling Factors
**Status**: ❌ SCIENTIFIC ERROR  
**Location**: Line 1e-4 scaling in calibrate_slc
**Problem**: Uses fixed scaling factors instead of calibration vectors
**Scientific Requirement**: Must use real calibration LUT from annotation XML

**Calibration Formula (CORRECT)**:
```
σ⁰ = |DN|² / (A_s · |calibrationVector|²)
```
Where:
- DN = Digital Number (SLC pixel value)  
- A_s = Calibration constant from annotation
- calibrationVector = Range-varying calibration from XML

**Current Implementation (WRONG)**:
```rust
intensity * 1e-4  // Fixed scaling factor - scientifically invalid
```

**Impact**: 
- Backscatter values completely unrealistic
- Cannot compare with other SAR data
- Invalid for scientific analysis

### 4. Missing Geometric Corrections

#### Issue: No Range Cell Migration Correction
**Status**: ❌ MISSING CRITICAL PROCESSING
**Location**: SLC reading and focusing
**Problem**: No correction for satellite motion during pulse transmission
**Impact**: Geometric distortions, reduced focusing quality

#### Issue: No Atmospheric Path Delay Correction  
**Status**: ❌ MISSING PROCESSING
**Location**: Range-Doppler terrain correction
**Problem**: No correction for atmospheric propagation delays
**Impact**: Geometric errors up to several meters

#### Issue: No Earth Rotation Correction
**Status**: ❌ MISSING PROCESSING  
**Location**: Zero-Doppler calculation
**Problem**: No correction for Earth rotation during data acquisition
**Impact**: Along-track positioning errors

### 5. Incorrect Reference Systems

#### Issue: Wrong Coordinate Reference System Usage
**Status**: ⚠️ PARTIALLY INCORRECT
**Location**: Terrain correction and export
**Problem**: Assumes WGS84 geographic (EPSG:4326) for all products
**Scientific Requirement**: Should support appropriate UTM zones for metric calculations

#### Issue: No Ellipsoid Height vs Orthometric Height Distinction
**Status**: ❌ MISSING GEODETIC ACCURACY
**Location**: DEM processing and terrain correction  
**Problem**: No distinction between ellipsoid heights (GPS) and orthometric heights (mean sea level)
**Impact**: Elevation errors up to 100+ meters in some regions

### 6. Missing Quality Assessment

#### Issue: No Radiometric Quality Validation
**Status**: ❌ NO QUALITY CONTROL
**Scientific Requirements Missing**:
- Cross-calibration with known targets
- Noise equivalent sigma zero (NESZ) assessment
- Antenna pattern correction validation
- Polarimetric leakage assessment (for dual-pol)

#### Issue: No Geometric Quality Validation  
**Status**: ❌ NO QUALITY CONTROL
**Scientific Requirements Missing**:
- Ground control point accuracy assessment
- Digital elevation model accuracy validation
- Orthorectification accuracy quantification
- Multi-temporal co-registration validation

### 7. Incorrect Processing Sequence

#### Issue: Wrong Multilooking Implementation
**Status**: ❌ SCIENTIFICALLY INCORRECT
**Problem**: Simple averaging instead of proper intensity multilooking

**Correct Multilooking (Intensity Domain)**:
```
I_ML = (1/N) × Σ|SLC_i|²
```

**Current Implementation**: May be doing amplitude averaging instead of intensity

#### Issue: Missing Speckle Statistics Validation
**Status**: ❌ NO VALIDATION
**Problem**: No validation that multilooking produces expected speckle statistics
**Required**: Verify coefficient of variation = 1/√N_looks

### 8. Missing Scientific References and Validation

#### Issue: No Algorithm Validation Against Literature
**Status**: ❌ NOT VALIDATED
**Missing Validations**:
- Comparison with ESA official processors (SNAP)
- Validation against published test datasets  
- Cross-comparison with other SAR processors
- Reproducibility testing with identical inputs

#### Issue: Insufficient Scientific Documentation
**Status**: ❌ INADEQUATE
**Missing Documentation**:
- Mathematical formulation of each processing step
- Error propagation analysis
- Accuracy specifications and limitations
- Processing parameter sensitivity analysis

## Impact on Research Validity

### Current Status: ❌ NOT SUITABLE FOR RESEARCH
1. **Geographic Limitation**: Only works for Northern Italy
2. **Invalid Calibration**: Backscatter values meaningless  
3. **Wrong Geometry**: Pixel spacing and coordinates incorrect
4. **Missing Critical Steps**: No proper terrain correction or flattening
5. **No Quality Control**: No validation of processing accuracy

### Required for Scientific Validity: ✅ RESEARCH GRADE
1. **Global Processing**: Works for any Sentinel-1 scene worldwide
2. **Proper Calibration**: Uses real calibration vectors from ESA
3. **Accurate Geometry**: Real pixel spacing and coordinates from annotation
4. **Complete Processing**: All 14 steps implemented correctly
5. **Quality Validation**: Comprehensive accuracy assessment

## Fixing Priority (Critical for Research)

### Phase 1: Fix Geographic Limitations (URGENT)
1. Remove all hardcoded coordinates  
2. Implement real scene geometry extraction
3. Calculate proper pixel spacing from metadata
4. Fix bounding box calculation for global coverage

### Phase 2: Fix Calibration (CRITICAL)  
1. Implement proper calibration vector extraction
2. Apply correct radiometric calibration formula
3. Validate against known targets
4. Add calibration quality metrics

### Phase 3: Complete Missing Processing Steps
1. Add proper terrain flattening with incidence angles
2. Implement complete range-Doppler terrain correction  
3. Add atmospheric and geometric corrections
4. Validate geometric accuracy

### Phase 4: Add Scientific Quality Control
1. Implement comprehensive validation framework
2. Add comparison with reference processors
3. Create accuracy assessment tools
4. Document all scientific methodologies

**Without these fixes, the package produces scientifically invalid results unsuitable for research publications.**
