# Contributing to SARdine 🤝

Thank you for your interest in contributing to SARdine! This guide will help you get started with contributing to our high-performance SAR processing library.

## 🎯 Ways to Contribute

### 1. **Code Contributions**
- Core SAR processing algorithms
- Performance optimizations
- New output formats
- Bug fixes and improvements

### 2. **Documentation**
- Example scripts and tutorials
- Algorithm explanations
- User guides and FAQ
- API documentation

### 3. **Testing and Validation**
- Test cases with real SAR data
- Performance benchmarks
- Scientific validation
- Cross-platform testing

### 4. **Community Support**
- Answering questions in issues
- Reviewing pull requests
- Sharing use cases and applications
- Bug reports and feature requests

## 🚀 Getting Started

### Prerequisites

```bash
# Required tools
- Rust 1.70+ (https://rustup.rs/)
- Python 3.8+
- Git

# Install Rust
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh

# Verify installation
rustc --version
python --version
```

### Development Setup

```bash
# 1. Fork the repository on GitHub
# 2. Clone your fork
git clone https://github.com/YOUR_USERNAME/SARdine.git
cd SARdine

# 3. Set up development environment
cd SARdine
cargo build --release

# 4. Install in development mode
pip install -e .

# 5. Run tests to verify setup
cargo test
python -m pytest tests/
```

## 📝 Development Workflow

### 1. **Create a Feature Branch**

```bash
# Create and switch to a new branch
git checkout -b feature/your-feature-name

# Or for bug fixes
git checkout -b fix/bug-description
```

### 2. **Make Your Changes**

#### For Rust Code (`SARdine/src/`)
- Follow Rust naming conventions
- Add comprehensive documentation
- Include unit tests
- Use `cargo fmt` for formatting
- Run `cargo clippy` for linting

#### For Python Code (`examples/`, `tests/`)
- Follow PEP 8 style guidelines
- Add docstrings for functions
- Include type hints where helpful
- Add comprehensive examples

### 3. **Testing Your Changes**

```bash
# Run Rust tests
cd SARdine
cargo test

# Run Python tests
cd ..
python -m pytest tests/ -v

# Test examples with real data
python examples/enhanced_geotiff_pipeline_v2.py
```

### 4. **Commit Your Changes**

```bash
# Stage your changes
git add .

# Commit with descriptive message
git commit -m "feat: add new SAR processing function

- Implement advanced speckle filtering
- Add comprehensive documentation
- Include unit tests and validation
- Performance: 20% faster processing"
```

#### Commit Message Format

```
<type>: <description>

[optional body]

[optional footer]
```

**Types**:
- `feat`: New features
- `fix`: Bug fixes  
- `docs`: Documentation changes
- `test`: Test additions/modifications
- `perf`: Performance improvements
- `refactor`: Code refactoring
- `style`: Formatting changes

### 5. **Submit a Pull Request**

```bash
# Push to your fork
git push origin feature/your-feature-name

# Create pull request on GitHub
# Include:
# - Clear description of changes
# - Test results
# - Performance impact
# - Breaking changes (if any)
```

## 🧪 Testing Guidelines

### Unit Tests

```rust
// Rust unit tests
#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_sar_processing_function() {
        // Test with realistic SAR data
        let test_data = create_test_slc_data();
        let result = process_sar_data(test_data);
        
        assert!(result.is_ok());
        assert_eq!(result.unwrap().len(), expected_length);
    }
}
```

```python
# Python tests
def test_python_wrapper():
    """Test Python bindings work correctly."""
    import sardine
    
    # Test with sample data
    result = sardine.some_function(test_params)
    assert result['status'] == 'success'
    assert len(result['data']) > 0
```

### Integration Tests

```python
# Test complete workflows
def test_complete_pipeline():
    """Test end-to-end processing pipeline."""
    from examples.enhanced_geotiff_pipeline_v2 import enhanced_complete_geotiff_pipeline_v2
    
    # Use test data
    results = enhanced_complete_geotiff_pipeline_v2("test_data.zip")
    
    # Validate outputs
    assert results is not None
    assert Path("output/backscatter_final.tif").exists()
```

### Performance Tests

```python
import time

def test_processing_performance():
    """Ensure processing meets performance benchmarks."""
    start_time = time.time()
    
    # Run processing
    process_large_dataset()
    
    processing_time = time.time() - start_time
    assert processing_time < 120  # Must complete within 2 minutes
```

## 📚 Documentation Standards

### Code Documentation

```rust
/// Calculate SAR backscatter coefficients from SLC data.
///
/// This function implements the ESA Sentinel-1 calibration procedure
/// as specified in the Product Definition document (2016).
///
/// # Arguments
/// * `slc_data` - Complex SLC data array
/// * `calibration_lut` - Calibration lookup table
/// * `noise_lut` - Noise lookup table
///
/// # Returns
/// * `Result<Array2<f32>, ProcessingError>` - Calibrated backscatter
///
/// # Example
/// ```rust
/// let backscatter = calculate_backscatter(slc, cal_lut, noise_lut)?;
/// ```
///
/// # References
/// - ESA Sentinel-1 Product Definition (2016)
/// - Torres et al. (2012): "GMES Sentinel-1 mission"
pub fn calculate_backscatter(
    slc_data: &Array2<Complex<f32>>,
    calibration_lut: &Array1<f32>,
    noise_lut: &Array1<f32>,
) -> Result<Array2<f32>, ProcessingError> {
    // Implementation...
}
```

### Example Documentation

```python
"""
Enhanced SAR Processing Pipeline with GeoTIFF Export

This example demonstrates complete processing of Sentinel-1 SLC data
to terrain-corrected GeoTIFF products suitable for GIS analysis.

Processing Steps:
1. Product analysis and metadata extraction
2. Precise orbit file integration
3. TOPSAR deburst with orbit corrections
4. Multilooking for speckle reduction
5. Advanced speckle filtering
6. Terrain flattening (γ⁰ correction)
7. Geocoding and terrain correction
8. Professional GeoTIFF export

Performance:
- Input: 4.5GB Sentinel-1 SLC
- Output: 208MB GeoTIFF
- Processing time: ~80 seconds
- Memory usage: <3GB

Usage:
    python enhanced_geotiff_pipeline_v2.py

Output:
    complete_output/backscatter_VH_final_wgs84.tif
    complete_output/backscatter_VH_final_utm.tif
"""
```

## 🔬 Scientific Contributions

### Algorithm Implementation

When implementing new SAR processing algorithms:

1. **Literature Review**
   - Cite relevant papers
   - Follow established procedures
   - Validate against known results

2. **Scientific Rigor**
   - Use proper physical constants
   - Implement peer-reviewed methods
   - Document assumptions and limitations

3. **Validation**
   - Test with real SAR data
   - Compare with reference implementations
   - Document accuracy and performance

### Example Algorithm Contribution

```rust
/// Implement the Lee speckle filter for SAR data.
///
/// Based on Lee (1980): "Digital image enhancement and noise filtering
/// by use of local statistics", IEEE Trans. PAMI.
///
/// # Algorithm
/// The Lee filter estimates the local mean and variance to adaptively
/// filter speckle noise while preserving image structure.
///
/// # Parameters
/// * `window_size` - Filter window size (typically 3, 5, or 7)
/// * `num_looks` - Number of looks in the data
///
/// # Performance
/// - Time complexity: O(n²) for n×n image
/// - Memory: O(window_size²) additional storage
///
/// # References
/// - Lee, J.S. (1980). IEEE Trans. PAMI, 2(2), 165-168
/// - Frost et al. (1982). IEEE Trans. PAMI, 4(2), 157-166
pub fn lee_filter(
    data: &Array2<f32>,
    window_size: usize,
    num_looks: f32,
) -> Array2<f32> {
    // Scientifically rigorous implementation...
}
```

## 🐛 Bug Reports

### Good Bug Report Format

```markdown
**Bug Description**
Clear description of the issue

**To Reproduce**
1. Install SARdine version X.Y.Z
2. Run command: `python example.py`
3. Use input file: `test_data.zip`
4. Observe error at step 3

**Expected Behavior**
Processing should complete successfully

**Actual Behavior**
```
Error: orbit integration failed
Traceback: ...
```

**Environment**
- OS: Ubuntu 22.04
- Python: 3.9.7
- Rust: 1.70.0
- SARdine: 0.2.0

**Additional Context**
- Input file size: 4.5GB
- Available memory: 8GB
- This worked in version 0.1.9
```

## 🚀 Feature Requests

### Feature Request Template

```markdown
**Feature Description**
Add support for Radarsat-2 data processing

**Use Case**
Enable researchers to process Canadian SAR data using the same
scientifically rigorous pipeline developed for Sentinel-1.

**Proposed Implementation**
1. Add Radarsat-2 metadata parser
2. Implement RS2-specific calibration
3. Add format conversion utilities
4. Update examples and documentation

**Additional Information**
- Similar to existing Sentinel-1 support
- Would benefit agricultural monitoring in Canada
- Compatible with existing GeoTIFF output
```

## 📊 Performance Guidelines

### Performance Expectations

- **Memory usage**: <4GB for standard processing
- **Processing time**: <2 minutes for 4.5GB SLC
- **Output quality**: Maintain scientific accuracy
- **Scalability**: Support files up to 10GB

### Optimization Tips

```rust
// Use efficient data structures
use ndarray::prelude::*;
use rayon::prelude::*;  // For parallelization

// Parallel processing example
fn parallel_processing(data: &Array2<f32>) -> Array2<f32> {
    data.axis_iter(Axis(0))
        .into_par_iter()
        .map(|row| process_row(&row))
        .collect()
}
```

## 🎨 Code Style

### Rust Style

```rust
// Use descriptive names
fn calculate_doppler_centroid_frequency(
    azimuth_time: f64,
    range_time: f64,
    orbit_state: &OrbitState,
) -> f64 {
    // Implementation follows ESA specifications
    let velocity = orbit_state.velocity_magnitude();
    let wavelength = SPEED_OF_LIGHT / CARRIER_FREQUENCY;
    
    // Calculate Doppler frequency
    2.0 * velocity / wavelength * azimuth_time.sin()
}
```

### Python Style

```python
def enhanced_complete_geotiff_pipeline_v2(
    input_file: str, 
    output_dir: str = "./complete_output"
) -> Optional[Dict[str, Any]]:
    """
    Process Sentinel-1 SLC data to terrain-corrected GeoTIFF.
    
    Args:
        input_file: Path to Sentinel-1 SLC ZIP file
        output_dir: Output directory for processed data
        
    Returns:
        Processing results dictionary or None if failed
        
    Raises:
        ProcessingError: If critical processing step fails
    """
    # Clear, documented implementation
    results = {}
    
    try:
        # Step 1: Product analysis
        product_info = sardine.get_product_info(input_file)
        results['product_analysis'] = product_info
        
        # Continue with remaining steps...
        
    except Exception as e:
        logger.error(f"Processing failed: {e}")
        return None
        
    return results
```

## 🏆 Recognition

Contributors will be recognized in:

- **CHANGELOG.md**: All contributions documented
- **README.md**: Major contributors acknowledged  
- **Academic papers**: Significant algorithm contributions cited
- **Community**: Recognition in discussions and presentations

## 📞 Getting Help

### Communication Channels

- **GitHub Issues**: Bug reports, feature requests
- **GitHub Discussions**: General questions, ideas
- **Pull Request Reviews**: Code-specific discussions
- **Email**: Direct contact for sensitive issues

### Before Asking for Help

1. **Search existing issues** for similar problems
2. **Check documentation** and examples
3. **Try the latest version** of SARdine
4. **Provide complete information** about your environment

---

**Thank you for contributing to SARdine!** 🛰️ 

Your contributions help make SAR processing more accessible and scientifically rigorous for researchers worldwide.
