# SARdine Validation Gates Framework

## Overview

SARdine now includes a comprehensive validation gates framework that enforces critical quality control measures throughout the SAR processing pipeline. This framework ensures robust, production-grade processing with automatic quality validation at each stage.

## Validation Gates

The framework implements 6 critical validation checks as specified for mission-critical SAR processing:

### 1. Power Preservation Test
- **Purpose**: Validates energy conservation during deburst/merge operations
- **Threshold**: ±1% tolerance (configurable)
- **Validation**: `mean(|S|²)` in overlap windows pre/post should match within tolerance
- **Benefits**: Prevents radiometric artifacts from power leakage

### 2. Uncovered Pixel Rate Validation
- **Purpose**: Ensures complete coverage after deburst operations
- **Threshold**: <0.2% uncovered pixels (configurable)
- **Validation**: Counts gaps and holes in debursted imagery
- **Benefits**: Early detection of timing/overlap issues

### 3. NESZ Compliance Check
- **Purpose**: Validates noise removal and calibration accuracy
- **Threshold**: ±1-2 dB tolerance from mission NESZ (configurable)
- **Validation**: Far-range ocean σ⁰ compared to expected noise floor
- **Benefits**: Confirms proper noise removal and calibration

### 4. Point Target Analysis (PTA)
- **Purpose**: Validates focusing quality and spatial resolution
- **Thresholds**: Sentinel-1 mission specifications
  - PSLR: <-13 dB
  - ISLR: <-10 dB  
  - IRW: <20 m
  - Radiometric accuracy: ±0.5 dB
- **Benefits**: Ensures focusing quality and spatial resolution compliance

### 5. Geometry Accuracy Validation
- **Purpose**: Validates orbit, timing, and DEM accuracy
- **Threshold**: <10m CE90 error (configurable)
- **Validation**: Known corner reflectors or coastline edges position accuracy
- **Benefits**: Validates orbit, timing, and DEM accuracy

### 6. LUT Domain Safety Checks
- **Purpose**: Prevents silent extrapolation errors
- **Validation**: Asserts all lookups are inside min/max line/pixel bounds
- **Benefits**: Stops processing on invalid lookups - prevents silent errors

## Usage

### Basic Validation Gates

```rust
use SARdine::core::{ValidationGates, ValidationResult};

// Create validation gates with default thresholds
let validation_gates = ValidationGates::default();

// Or create with custom thresholds
let validation_gates = ValidationGates::with_thresholds(
    0.005,  // ±0.5% power preservation tolerance
    0.001,  // 0.1% uncovered pixel threshold  
    1.5,    // ±1.5 dB NESZ tolerance
    8.0,    // 8m CE90 geometry threshold
);

// Validate power preservation
let result = validation_gates.validate_power_preservation(
    &pre_deburst_power, 
    &post_deburst_power, 
    &overlap_regions
)?;

if result.passed {
    println!("✅ Power preservation: {:.4}% error", result.measured_value * 100.0);
} else {
    println!("❌ {}", result.diagnostic_message);
    for rec in &result.recommendations {
        println!("💡 {}", rec);
    }
}
```

### Integrated Processing Pipeline

```rust
use SARdine::core::{ValidatedProcessingPipeline, ProcessingParameters};

// Create validated processing pipeline
let mut pipeline = ValidatedProcessingPipeline::with_validation_gates(validation_gates);

// Configure for fail-fast mode (stop on first validation failure)
pipeline.fail_fast = true;
pipeline.log_detailed_results = true;

// Process with comprehensive validation
match pipeline.process_sar_with_validation(
    &slc_data,
    &calibration_lut,
    &dem_data,
    &processing_params
) {
    Ok(processed_data) => {
        println!("🎉 Processing completed successfully!");
        
        let summary = &processed_data.validation_summary;
        if summary.overall_success {
            println!("✅ All validation gates PASSED");
        } else {
            println!("⚠️ Some validation gates failed");
        }
    },
    Err(e) => {
        println!("❌ Processing failed: {}", e);
        // Critical validation failures stop processing
    }
}
```

## Key Features

- **Fail-Fast Mode**: Stops processing on critical validation failures
- **Detailed Diagnostics**: Provides specific error messages and actionable recommendations
- **Configurable Thresholds**: Adapt validation criteria for different quality requirements
- **Stage-by-Stage Validation**: Validates at each processing step for early error detection
- **LUT Domain Safety**: Prevents silent extrapolation errors that could corrupt results
- **Production Ready**: Comprehensive error handling with clear failure modes

## Critical Error Handling

The validation framework provides robust error handling with specific guidance:

- **LUT Domain Violations**: Immediate stop with clear error message indicating data bounds issues
- **Power Preservation Failures**: Diagnostic recommendations for algorithm parameter tuning
- **High Uncovered Pixel Rates**: Targeted troubleshooting for timing and overlap parameters
- **Geometry Errors**: Analysis separating azimuth vs range error sources

## Validation Results

Each validation gate returns a `ValidationResult` with:

```rust
pub struct ValidationResult {
    pub passed: bool,
    pub test_name: String,
    pub measured_value: f64,
    pub threshold_value: f64,
    pub diagnostic_message: String,
    pub recommendations: Vec<String>,
}
```

The complete pipeline provides a `ValidationSummary`:

```rust
pub struct ValidationSummary {
    pub total_tests: usize,
    pub passed_tests: usize,
    pub failed_tests: usize,
    pub critical_failures: Vec<String>,
    pub warnings: Vec<String>,
    pub overall_success: bool,
}
```

## Integration with SARdine Core

The validation gates are fully integrated into SARdine's core processing pipeline:

```rust
// Available exports
use SARdine::core::{
    ValidationGates,
    ValidationResult, 
    ValidationSummary,
    ValidatedProcessingPipeline,
    ProcessingParameters,
    ProcessedSarData,
    DeburstMetrics,
};
```

## Performance

The validation framework is designed for production use with minimal performance impact:

- Efficient validation algorithms with O(n) complexity
- Optional validation levels (strict vs permissive)
- Configurable diagnostic detail levels
- Parallel validation when possible

## Testing

Comprehensive test suite validates all validation gates:

```bash
# Test validation gates
cargo test validation_gates --lib -- --nocapture

# Test validated processing pipeline  
cargo test validated_processing_pipeline --lib -- --nocapture

# Run comprehensive validation test
python3 test_validation_gates.py
```

## Benefits for Production Use

1. **Quality Assurance**: Automatic detection of processing artifacts and errors
2. **Reproducibility**: Consistent quality validation across different datasets
3. **Debugging Support**: Clear diagnostic messages accelerate troubleshooting
4. **Mission Compliance**: Validation against official Sentinel-1 specifications
5. **Silent Error Prevention**: LUT domain checks prevent corrupted outputs
6. **Early Error Detection**: Stage-by-stage validation stops bad processing early

The validation gates framework transforms SARdine into a production-grade SAR processing system with comprehensive quality control, ensuring reliable and accurate results for mission-critical applications.