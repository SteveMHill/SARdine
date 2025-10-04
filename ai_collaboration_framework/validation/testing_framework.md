# Testing Framework and Guidelines

## Overview
This document establishes comprehensive testing strategies and guidelines for the SARdine project, ensuring that all code meets scientific accuracy requirements while maintaining software quality standards.

## Testing Philosophy

### Core Principles
1. **Scientific Rigor**: All scientific algorithms must be validated against known standards
2. **Comprehensive Coverage**: Testing covers functionality, performance, and edge cases
3. **Automated Validation**: Tests run automatically to catch regressions early
4. **Real-World Validation**: Testing includes realistic SAR data scenarios
5. **Continuous Integration**: Tests are integrated into the development workflow

### Testing Pyramid Strategy
```
    /\     Unit Tests (70%)
   /  \    - Fast execution
  /____\   - High coverage
 /      \  - Isolated testing
/__________\ 
Integration Tests (20%)  System Tests (10%)
- Module interaction     - End-to-end validation
- API testing           - Performance testing
- Data flow validation  - Scientific accuracy
```

## Unit Testing Framework

### Test Structure and Organization
```python
"""
Standard test file structure for SARdine
"""
import pytest
import numpy as np
from numpy.testing import assert_allclose, assert_array_equal
from sardine.processing import backscatter_calibration
from sardine.testing import TestDataGenerator, ScientificValidators

class TestBackscatterCalibration:
    """
    Comprehensive test suite for backscatter calibration
    """
    
    def setup_method(self):
        """Setup test data and validators for each test"""
        self.test_data = TestDataGenerator()
        self.validators = ScientificValidators()
        
    def test_calibration_accuracy_known_target(self):
        """Test calibration accuracy with known point target"""
        # Generate synthetic point target with known RCS
        target_data = self.test_data.generate_point_target(
            rcs_dbm2=10.0,  # Known radar cross section
            position=(1000, 500),
            noise_level=-30  # dB
        )
        
        # Apply calibration
        calibrated = backscatter_calibration(
            raw_data=target_data.raw,
            calibration_lut=target_data.cal_lut,
            noise_equivalent_sigma0=target_data.nesz
        )
        
        # Validate against known theoretical value
        measured_rcs = self.validators.extract_target_rcs(
            calibrated, target_data.position
        )
        
        # Scientific accuracy requirement: ±0.5 dB
        expected_rcs = 10.0  # dB
        assert abs(measured_rcs - expected_rcs) < 0.5, \
            f"RCS error {measured_rcs - expected_rcs:.2f} dB exceeds ±0.5 dB tolerance"
    
    def test_radiometric_stability(self):
        """Test radiometric stability across scene"""
        # Generate homogeneous distributed target
        scene_data = self.test_data.generate_distributed_target(
            sigma0_db=-15.0,  # Constant backscatter
            scene_size=(2000, 1000),
            correlation_length=50.0
        )
        
        # Apply calibration
        calibrated = backscatter_calibration(
            raw_data=scene_data.raw,
            calibration_lut=scene_data.cal_lut,
            noise_equivalent_sigma0=scene_data.nesz
        )
        
        # Validate radiometric stability
        stability = self.validators.assess_radiometric_stability(calibrated)
        
        # Requirement: <0.3 dB standard deviation
        assert stability.std_db < 0.3, \
            f"Radiometric stability {stability.std_db:.2f} dB exceeds 0.3 dB limit"
    
    def test_edge_cases_and_boundaries(self):
        """Test handling of edge cases and boundary conditions"""
        test_cases = [
            # Very low signal (noise limited)
            {'signal_level': -40, 'noise_level': -35, 'expected_behavior': 'noise_flag'},
            # Very high signal (saturation)
            {'signal_level': 40, 'noise_level': -30, 'expected_behavior': 'saturation_check'},
            # Zero signal
            {'signal_level': -np.inf, 'noise_level': -30, 'expected_behavior': 'handle_gracefully'},
            # Missing calibration data
            {'signal_level': -10, 'cal_data': None, 'expected_behavior': 'raise_exception'}
        ]
        
        for case in test_cases:
            with self.subTest(case=case):
                self._test_edge_case(case)
```

### Scientific Validation Functions
```python
class ScientificValidators:
    """
    Collection of scientific validation methods
    """
    
    def __init__(self):
        self.tolerance_db = 0.5  # Default tolerance for radiometric accuracy
        
    def validate_speckle_statistics(self, sigma0_image, theoretical_looks):
        """
        Validate that speckle statistics match theoretical expectations
        """
        # Calculate empirical statistics
        mean_linear = np.mean(10**(sigma0_image/10))
        var_linear = np.var(10**(sigma0_image/10))
        
        # Theoretical statistics for gamma distribution
        theoretical_cov = 1.0 / np.sqrt(theoretical_looks)
        empirical_cov = np.sqrt(var_linear) / mean_linear
        
        # Validate coefficient of variation
        assert abs(empirical_cov - theoretical_cov) < 0.1, \
            f"Speckle statistics mismatch: {empirical_cov:.3f} vs {theoretical_cov:.3f}"
    
    def validate_geometric_accuracy(self, processed_image, ground_control_points):
        """
        Validate geometric accuracy using ground control points
        """
        geometric_errors = []
        
        for gcp in ground_control_points:
            predicted_pixel = self.geo_to_pixel(
                processed_image.geotransform,
                gcp.latitude, gcp.longitude
            )
            actual_pixel = gcp.pixel_coordinates
            
            error_pixels = np.sqrt(
                (predicted_pixel[0] - actual_pixel[0])**2 +
                (predicted_pixel[1] - actual_pixel[1])**2
            )
            geometric_errors.append(error_pixels)
        
        # Requirement: <1 pixel RMS error
        rms_error = np.sqrt(np.mean(np.array(geometric_errors)**2))
        assert rms_error < 1.0, \
            f"Geometric accuracy {rms_error:.2f} pixels exceeds 1.0 pixel limit"
    
    def validate_phase_preservation(self, complex_image_before, complex_image_after):
        """
        Validate that processing preserves phase information appropriately
        """
        # Calculate coherence between before and after
        coherence = self.calculate_coherence(complex_image_before, complex_image_after)
        
        # For stable targets, coherence should be high
        stable_threshold = 0.9
        stable_fraction = np.sum(coherence > stable_threshold) / coherence.size
        
        # Requirement: >80% of stable targets maintain coherence >0.9
        assert stable_fraction > 0.8, \
            f"Phase preservation insufficient: {stable_fraction:.2f} < 0.8"
```

## Integration Testing Framework

### End-to-End Pipeline Testing
```python
class TestSARProcessingPipeline:
    """
    Integration tests for complete SAR processing pipeline
    """
    
    @pytest.fixture
    def real_sentinel1_data(self):
        """Fixture providing real Sentinel-1 test data"""
        return {
            'slc_path': 'tests/data/S1A_IW_SLC_test.SAFE',
            'orbit_path': 'tests/data/S1A_OPER_AUX_POEORB_test.EOF',
            'dem_path': 'tests/data/srtm_test_area.tif',
            'expected_results': {
                'radiometric_accuracy': 0.3,  # dB
                'geometric_accuracy': 0.5,    # pixels
                'processing_time': 120,       # seconds for test scene
            }
        }
    
    def test_full_rtc_pipeline(self, real_sentinel1_data):
        """Test complete RTC processing pipeline with real data"""
        from sardine.cli import process_rtc
        
        # Configure processing parameters
        config = {
            'input_slc': real_sentinel1_data['slc_path'],
            'orbit_file': real_sentinel1_data['orbit_path'],
            'dem_file': real_sentinel1_data['dem_path'],
            'output_dir': 'tests/output/rtc_integration',
            'polarizations': ['VV', 'VH'],
            'dem_buffer': 0.1,  # degrees
            'multilook_factors': {'range': 4, 'azimuth': 1},
            'speckle_filter': {'type': 'lee', 'window_size': 7}
        }
        
        # Execute processing
        start_time = time.time()
        results = process_rtc(config)
        processing_time = time.time() - start_time
        
        # Validate processing results
        self._validate_rtc_outputs(results, real_sentinel1_data['expected_results'])
        
        # Check processing time
        assert processing_time < real_sentinel1_data['expected_results']['processing_time'], \
            f"Processing too slow: {processing_time:.1f}s > {real_sentinel1_data['expected_results']['processing_time']}s"
    
    def _validate_rtc_outputs(self, results, expected):
        """Validate RTC processing outputs"""
        # Check that all expected outputs exist
        required_outputs = [
            'backscatter_VV.tif',
            'backscatter_VH.tif',
            'metadata.json',
            'processing_log.json',
            'quality_report.json'
        ]
        
        for output_file in required_outputs:
            output_path = os.path.join(results['output_dir'], output_file)
            assert os.path.exists(output_path), f"Missing output: {output_file}"
        
        # Validate geospatial properties
        for pol in ['VV', 'VH']:
            geotiff_path = os.path.join(results['output_dir'], f'backscatter_{pol}.tif')
            self._validate_geotiff_properties(geotiff_path)
            
        # Validate scientific accuracy
        self._validate_scientific_accuracy(results, expected)
```

### Performance Testing Framework
```python
class TestPerformanceRequirements:
    """
    Performance and scalability testing
    """
    
    def test_processing_speed_requirements(self):
        """Test that processing meets speed requirements"""
        test_scenes = [
            {'size': (1000, 1000), 'max_time': 30},    # Small scene
            {'size': (5000, 5000), 'max_time': 300},   # Medium scene  
            {'size': (10000, 10000), 'max_time': 900}, # Large scene
        ]
        
        for scene in test_scenes:
            synthetic_data = self.generate_synthetic_scene(scene['size'])
            
            start_time = time.time()
            processed = process_backscatter(synthetic_data)
            processing_time = time.time() - start_time
            
            pixels_per_second = (scene['size'][0] * scene['size'][1]) / processing_time
            
            # Log performance metrics
            print(f"Scene {scene['size']}: {pixels_per_second:.0f} pixels/second")
            
            assert processing_time < scene['max_time'], \
                f"Processing too slow for {scene['size']}: {processing_time:.1f}s > {scene['max_time']}s"
    
    def test_memory_usage_limits(self):
        """Test memory usage stays within acceptable limits"""
        import psutil
        import gc
        
        # Monitor memory during processing
        process = psutil.Process()
        initial_memory = process.memory_info().rss / 1024 / 1024  # MB
        
        # Process large synthetic dataset
        large_scene = self.generate_synthetic_scene((8000, 8000))
        processed = process_backscatter(large_scene)
        
        peak_memory = process.memory_info().rss / 1024 / 1024  # MB
        memory_usage = peak_memory - initial_memory
        
        # Clean up
        del large_scene, processed
        gc.collect()
        
        # Memory requirement: <4GB for large scenes
        max_memory_mb = 4000
        assert memory_usage < max_memory_mb, \
            f"Memory usage {memory_usage:.0f}MB exceeds {max_memory_mb}MB limit"
```

## Test Data Management

### Synthetic Data Generation
```python
class TestDataGenerator:
    """
    Generate synthetic SAR data for testing
    """
    
    def __init__(self):
        self.random_seed = 42  # Reproducible random data
        np.random.seed(self.random_seed)
    
    def generate_point_target(self, rcs_dbm2, position, noise_level=-30):
        """Generate synthetic point target with known properties"""
        # Create scene with speckle background
        scene_size = (2000, 1000)
        background = self._generate_speckle_background(scene_size, sigma0_db=-20)
        
        # Add point target
        target_amplitude = self._rcs_to_amplitude(rcs_dbm2, position)
        background[position[1], position[0]] += target_amplitude
        
        # Add thermal noise
        noise = self._generate_thermal_noise(scene_size, noise_level)
        raw_data = background + noise
        
        # Generate calibration LUT
        cal_lut = self._generate_calibration_lut(scene_size)
        
        return SyntheticTargetData(
            raw=raw_data,
            cal_lut=cal_lut,
            nesz=noise_level,
            position=position,
            true_rcs=rcs_dbm2
        )
    
    def generate_distributed_target(self, sigma0_db, scene_size, correlation_length):
        """Generate synthetic distributed target scene"""
        # Generate correlated speckle
        speckle = self._generate_correlated_speckle(
            scene_size, correlation_length, sigma0_db
        )
        
        # Add thermal noise
        noise = self._generate_thermal_noise(scene_size, -35)
        raw_data = speckle + noise
        
        # Generate calibration data
        cal_lut = self._generate_calibration_lut(scene_size)
        
        return SyntheticSceneData(
            raw=raw_data,
            cal_lut=cal_lut,
            nesz=-35,
            true_sigma0=sigma0_db
        )
```

### Reference Data Validation
```python
class ReferenceDataValidator:
    """
    Validate processing results against reference data
    """
    
    def __init__(self):
        self.reference_database = self._load_reference_database()
    
    def validate_against_reference(self, processed_result, reference_id):
        """Compare processed result with reference implementation"""
        reference = self.reference_database[reference_id]
        
        # Radiometric comparison
        radiometric_diff = self._compare_radiometric(
            processed_result.sigma0, reference.sigma0
        )
        
        # Geometric comparison
        geometric_diff = self._compare_geometric(
            processed_result.geotransform, reference.geotransform
        )
        
        # Statistical comparison
        statistical_diff = self._compare_statistics(
            processed_result.sigma0, reference.sigma0
        )
        
        return ValidationReport(
            radiometric=radiometric_diff,
            geometric=geometric_diff,
            statistical=statistical_diff,
            overall_quality=self._assess_overall_quality(
                radiometric_diff, geometric_diff, statistical_diff
            )
        )
```

## Continuous Integration and Testing

### Automated Test Execution
```yaml
# Example CI configuration for automated testing
name: SARdine Quality Assurance

on: [push, pull_request]

jobs:
  unit-tests:
    runs-on: ubuntu-latest
    strategy:
      matrix:
        python-version: [3.8, 3.9, '3.10']
    
    steps:
    - uses: actions/checkout@v3
    - name: Set up Python ${{ matrix.python-version }}
      uses: actions/setup-python@v3
      with:
        python-version: ${{ matrix.python-version }}
    
    - name: Install dependencies
      run: |
        pip install -e .[dev]
        pip install pytest pytest-cov pytest-xdist
    
    - name: Run unit tests with coverage
      run: |
        pytest tests/unit/ -v --cov=sardine --cov-report=xml --cov-report=html
        
    - name: Upload coverage reports
      uses: codecov/codecov-action@v3
      with:
        file: ./coverage.xml

  integration-tests:
    runs-on: ubuntu-latest
    needs: unit-tests
    
    steps:
    - uses: actions/checkout@v3
    - name: Set up Python
      uses: actions/setup-python@v3
      with:
        python-version: '3.10'
        
    - name: Download test data
      run: |
        mkdir -p tests/data
        wget -O tests/data/test_scene.zip "https://test-data-repo/sentinel1_test.zip"
        unzip tests/data/test_scene.zip -d tests/data/
        
    - name: Run integration tests
      run: |
        pytest tests/integration/ -v --timeout=600
        
  scientific-validation:
    runs-on: ubuntu-latest
    needs: integration-tests
    
    steps:
    - uses: actions/checkout@v3
    - name: Set up Python
      uses: actions/setup-python@v3
      with:
        python-version: '3.10'
        
    - name: Run scientific validation suite
      run: |
        pytest tests/scientific/ -v --scientific-accuracy
        
    - name: Generate validation report
      run: |
        python scripts/generate_validation_report.py
        
    - name: Upload validation artifacts
      uses: actions/upload-artifact@v3
      with:
        name: validation-report
        path: validation_report.html
```

### Test Reporting and Metrics
```python
class TestMetricsCollector:
    """
    Collect and analyze test execution metrics
    """
    
    def __init__(self):
        self.metrics = {
            'coverage': {},
            'performance': {},
            'scientific_accuracy': {},
            'reliability': {}
        }
    
    def collect_coverage_metrics(self, coverage_report):
        """Extract coverage metrics from test report"""
        self.metrics['coverage'] = {
            'line_coverage': coverage_report.line_coverage_percent,
            'branch_coverage': coverage_report.branch_coverage_percent,
            'function_coverage': coverage_report.function_coverage_percent,
            'uncovered_lines': coverage_report.uncovered_lines,
            'missing_tests': self._identify_missing_tests(coverage_report)
        }
    
    def collect_performance_metrics(self, performance_results):
        """Extract performance metrics from test results"""
        self.metrics['performance'] = {
            'processing_speed': performance_results.pixels_per_second,
            'memory_usage': performance_results.peak_memory_mb,
            'scalability': performance_results.scaling_factor,
            'bottlenecks': performance_results.identified_bottlenecks
        }
    
    def generate_quality_dashboard(self):
        """Generate comprehensive quality dashboard"""
        dashboard = QualityDashboard()
        dashboard.add_coverage_section(self.metrics['coverage'])
        dashboard.add_performance_section(self.metrics['performance'])
        dashboard.add_scientific_section(self.metrics['scientific_accuracy'])
        dashboard.add_reliability_section(self.metrics['reliability'])
        
        return dashboard.generate_html_report()
```

---

**Framework Version**: 1.0
**Test Coverage Target**: >95%
**Scientific Accuracy Tolerance**: ±0.5 dB
**Performance Requirement**: >1M pixels/second