# Step-by-Step Implementation Plan

## Implementation Strategy: From Non-Functional to Research-Grade

### Phase 1: Critical Foundation (Week 1) 🔥 URGENT

#### Step 1.1: Fix Import Issues and Basic Structure
**Priority**: CRITICAL - Package currently broken
**Time**: 1 day

1. **Fix BackscatterProcessor Import**
   ```python
   # In cli.py line 1483 - CURRENTLY BROKEN:
   from backscatter_cli import BackscatterProcessor  # FILE DOES NOT EXIST
   
   # FIX TO:
   from sardine.processors.backscatter import BackscatterProcessor
   ```

2. **Create Missing Python Module Structure**
   ```
   sardine/
   ├── processors/
   │   ├── __init__.py
   │   └── backscatter.py  ✅ EXISTS
   ├── core/
   │   ├── __init__.py
   │   ├── orbit.py        ❌ CREATE
   │   ├── annotation.py   ❌ CREATE  
   │   └── validation.py   ❌ CREATE
   └── utils/
       ├── __init__.py
       ├── xml_parser.py   ❌ CREATE
       └── orbit_urls.py   ❌ CREATE
   ```

#### Step 1.2: Implement Core Data Structures in Rust
**Priority**: CRITICAL - Blocks all processing
**Time**: 2 days

1. **Add AnnotationData Structure** (`src/types.rs`)
   ```rust
   #[derive(Debug, Clone)]
   pub struct AnnotationData {
       pub bursts: Vec<BurstData>,
       pub calibration_vectors: Vec<CalibrationVector>,
       pub product_info: ProductInfo,
       // ... complete structure
   }
   ```

2. **Implement Annotation XML Parser** (`src/io/annotation.rs`)
   ```rust
   pub fn parse_annotation_xml(xml_path: &str) -> Result<AnnotationData, Box<dyn std::error::Error>> {
       // Parse Sentinel-1 annotation XML
       // Extract burst parameters
       // Extract calibration vectors
       // Return structured data
   }
   ```

#### Step 1.3: Fix Critical Missing SlcReader Methods  
**Priority**: CRITICAL - Blocks scientific processing
**Time**: 2 days

1. **Implement `get_pixel_spacing()`**
   ```rust
   fn get_pixel_spacing(&mut self, polarization: String) -> PyResult<std::collections::HashMap<String, f64>> {
       // Extract from annotation XML
       // Return real range/azimuth pixel spacing
   }
   ```

2. **Implement `get_scene_bbox()`**  
   ```rust
   fn get_scene_bbox(&mut self) -> PyResult<Vec<f64>> {
       // Calculate from annotation coordinates
       // Return [min_lon, min_lat, max_lon, max_lat]
   }
   ```

3. **Implement `read_annotation_data()`**
   ```rust
   fn read_annotation_data(&mut self, pol: Polarization) -> Result<AnnotationData, Box<dyn std::error::Error>> {
       // Parse annotation XML for polarization
       // Return complete annotation data
   }
   ```

### Phase 2: Scientific Processing Core (Week 2) 🔬 SCIENCE

#### Step 2.1: Fix Calibration to Use Real Vectors
**Priority**: CRITICAL - Invalid backscatter values
**Time**: 2 days

1. **Implement Real Calibration Vector Extraction**
   ```rust
   // In calibrate_slc() - REMOVE hardcoded scaling:
   intensity * 1e-4  // ❌ WRONG
   
   // ADD proper calibration:
   let calibration_lut = extract_calibration_vectors(&annotation_data)?;
   let sigma0 = |dn|² / (calibration_constant * |calibration_vector|²);
   ```

2. **Validate Calibration Against ESA Standards**
   - Compare with SNAP processor output
   - Validate with known calibration targets
   - Check backscatter value ranges

#### Step 2.2: Fix Deburst to Use Real Parameters
**Priority**: HIGH - Wrong geometric accuracy  
**Time**: 2 days

1. **Remove ALL Hardcoded Burst Parameters**
   ```rust
   // REMOVE these hardcoded values:
   azimuth_fm_rate: -2300.0,     // ❌ HARDCODED
   doppler_centroid: 150.0,      // ❌ HARDCODED  
   burst_duration: 2.758,        // ❌ HARDCODED
   
   // REPLACE with real extraction:
   let burst_params = extract_burst_parameters_from_annotation(&annotation_data)?;
   ```

2. **Implement Complete `extract_burst_parameters_from_annotation()`**
   - Parse all burst timing parameters
   - Extract FM rates per burst
   - Get real doppler centroids
   - Calculate proper burst boundaries

#### Step 2.3: Fix Multilooking to Use Real Pixel Spacing
**Priority**: HIGH - Wrong output resolution
**Time**: 1 day

1. **Update Function Signature** (ALREADY DONE ✅)
   ```rust
   fn apply_multilooking(
       data: PyReadonlyArray2<f32>,
       range_looks: usize,
       azimuth_looks: usize,
       input_range_spacing: f64,    // ✅ ADDED
       input_azimuth_spacing: f64,  // ✅ ADDED  
   )
   ```

2. **Fix Python Calls** (ALREADY DONE ✅)
   ```python
   multilooked_data = sardine.apply_multilooking(
       cal_array, 
       range_looks, 
       azimuth_looks,
       range_pixel_spacing,  # ✅ REAL FROM SLC
       azimuth_pixel_spacing  # ✅ REAL FROM SLC
   )
   ```

### Phase 3: Terrain Processing Implementation (Week 3) 🏔️ GEOGRAPHY

#### Step 3.1: Complete Terrain Correction Implementation
**Priority**: HIGH - Step 9 partially broken
**Time**: 3 days

1. **Implement Missing DEM Reader**
   ```rust
   // Create src/io/dem.rs
   pub struct DemReader;
   impl DemReader {
       pub fn prepare_dem_for_scene(bbox: &BoundingBox, resolution: f64, cache_dir: &str) 
           -> Result<(Array2<f32>, GeoTransform), Box<dyn std::error::Error>> {
           // Download SRTM/ASTER DEM data
           // Crop to scene bbox  
           // Resample to target resolution
       }
   }
   ```

2. **Fix Terrain Correction to Use Real Metadata** (ALREADY DONE ✅)
   - Function signature updated to require real metadata
   - Uses real pixel spacing, PRF, wavelength from SLC
   - No more hardcoded Sentinel-1 "typical" values

3. **Implement `apply_terrain_correction_with_real_orbits()`**
   ```rust
   #[pyfunction]  
   fn apply_terrain_correction_with_real_orbits(
       sar_image: PyReadonlyArray2<f32>,
       sar_bbox: Vec<f64>,
       orbit_data: std::collections::HashMap<String, Vec<f64>>,
       cache_dir: String,
       dem_resolution: f64,
   ) -> PyResult<PyObject> {
       // Wrapper around apply_terrain_correction
       // Handles orbit data conversion
       // Downloads DEM automatically
   }
   ```

#### Step 3.2: Complete Terrain Flattening (ALREADY IMPLEMENTED ✅)
**Status**: ✅ DONE - Step 8 fixed with Small & Schubert (2008) methodology

#### Step 3.3: Fix Scene Geometry (PARTIALLY DONE ✅)
**Status**: ✅ Code updated to extract real scene bbox, needs implementation

### Phase 4: Export and Validation (Week 4) 📊 QUALITY

#### Step 4.1: Fix GeoTIFF Export (PARTIALLY DONE ✅)  
**Status**: ✅ Code updated to use real geospatial transform, needs implementation

1. **Implement `get_output_geotransform()`**
   ```rust
   fn get_output_geotransform(&mut self, data_shape: (usize, usize), bbox: Vec<f64>, resolution: f64) 
       -> PyResult<Vec<f64>> {
       // Calculate real geotransform from processed data
       // Handle different coordinate systems (UTM, Geographic)  
       // Return proper GDAL geotransform
   }
   ```

2. **Implement `export_geotiff()`**
   ```rust
   #[pyfunction]
   fn export_geotiff(
       data: PyReadonlyArray2<f32>,
       output_path: String,
       geotransform: Vec<f64>,
       crs: i32,
       nodata: Option<f32>,
   ) -> PyResult<PyObject> {
       // Create properly georeferenced GeoTIFF
       // Set coordinate system
       // Handle nodata values  
   }
   ```

#### Step 4.2: Implement Orbit File Processing
**Priority**: HIGH - Required for scientific processing
**Time**: 2 days

1. **Implement `load_orbit_file()`**
   ```rust
   #[pyfunction]
   fn load_orbit_file(orbit_file_path: String) -> PyResult<PyObject> {
       // Parse ESA .EOF orbit files
       // Extract state vectors (time, position, velocity)
       // Validate orbit data quality
       // Return in proper Python format
   }
   ```

2. **Implement Orbit File Download in Python**
   ```python
   def download_orbit_files(self, cache_dir):
       """Download real .EOF orbit files from ESA servers"""
       # Calculate orbit file URLs from product metadata  
       # Download from ESA servers
       # Validate downloaded files
       return downloaded_file_paths
   ```

#### Step 4.3: Add Scientific Validation Framework
**Priority**: MEDIUM - Quality assurance
**Time**: 2 days

1. **Create Validation Functions**
   ```python
   def validate_orbit_data(orbit_data):
       """Validate orbit data accuracy and completeness"""
   
   def validate_calibration_quality(cal_data):
       """Validate calibration against expected ranges"""
   
   def validate_geometric_accuracy(processed_data):  
       """Validate geometric accuracy of processed data"""
   ```

### Phase 5: Complete Missing CLI Commands (Week 5) 🖥️ INTERFACE  

#### Step 5.1: Implement Empty CLI Commands
**Priority**: MEDIUM - User functionality  
**Time**: 3 days

1. **Implement `cmd_speckle_filter()`**
2. **Implement `cmd_estimate_nlooks()`** 
3. **Implement `cmd_orbit()`**
4. **Implement `cmd_iw_split()`**
5. **Complete `cmd_deburst()`**

#### Step 5.2: Add Comprehensive Error Handling
**Priority**: MEDIUM - User experience
**Time**: 1 day

1. **Create Scientific Error Classes**
   ```python
   class RealDataRequiredError(Exception): pass
   class OrbitDataError(Exception): pass
   class AnnotationParsingError(Exception): pass
   ```

### Phase 6: Testing and Validation (Week 6) ✅ TESTING

#### Step 6.1: Create Test Framework  
**Priority**: HIGH - Scientific validation
**Time**: 3 days

1. **Unit Tests for Core Functions**
2. **Integration Tests for 14-Step Pipeline**
3. **Scientific Validation Tests**
4. **Comparison Tests with SNAP Processor**

#### Step 6.2: Performance Optimization
**Priority**: LOW - Functionality first
**Time**: 2 days

1. **Optimize Rust Core Functions**
2. **Add Parallel Processing**  
3. **Memory Usage Optimization**

## Implementation Milestones

### Week 1 Milestone: ✅ FUNCTIONAL PACKAGE
- All imports working
- Basic commands functional  
- Core data structures implemented
- Can process at least one complete scene

### Week 2 Milestone: 🔬 SCIENTIFIC ACCURACY
- Real calibration vectors used
- Real burst parameters extracted
- Real pixel spacing calculations
- No more hardcoded values

### Week 3 Milestone: 🏔️ COMPLETE PROCESSING  
- All 14 steps implemented
- Terrain correction working globally
- DEM processing functional
- Export to GeoTIFF working

### Week 4 Milestone: 📊 RESEARCH GRADE
- Orbit file processing complete
- Scientific validation implemented
- Quality assessment working
- Global scene coverage

### Week 5 Milestone: 🖥️ COMPLETE INTERFACE
- All CLI commands working
- Comprehensive error handling
- User-friendly feedback
- Documentation complete

### Week 6 Milestone: ✅ PRODUCTION READY
- Full test coverage
- Performance optimized
- Scientific validation complete
- Ready for research use

## Resource Requirements

### Development Time: 6 weeks (1 developer)
### Critical Dependencies:
- Rust XML parsing crate (quick-xml)
- DEM download capability (SRTM/ASTER APIs)
- GDAL for GeoTIFF export
- ESA orbit file server access

### Success Criteria:
1. ✅ All 14 processing steps functional  
2. ✅ Works for any global Sentinel-1 scene
3. ✅ Uses only real Sentinel-1 data (no synthetic)
4. ✅ Produces research-grade results
5. ✅ Comprehensive test coverage
6. ✅ Scientific validation complete

**Current Status**: Week 0 - Package partially functional with critical scientific issues
**Target**: Week 6 - Production-ready research-grade SAR processing package
