# SARdine Project: Scientific Software Development Standards

## Project Mission Statement

**Creating reliable, research-grade SAR processing software that researchers, scientists, and students can trust for real-world remote sensing applications.**

SARdine is committed to providing scientifically accurate, peer-reviewable Synthetic Aperture Radar processing tools that meet the rigorous standards required for climate research, environmental monitoring, and geoscience applications.

---

## Project Overview

### **Target Domain**: Synthetic Aperture Radar (SAR) Processing
- **Primary Focus**: Sentinel-1 backscatter processing pipeline
- **Target Users**: Graduate students, research scientists, operational analysts
- **Use Cases**: Land use monitoring, deforestation tracking, disaster response, climate studies
- **Quality Standard**: Publication-ready, peer-reviewable scientific accuracy

### **Current Project Status**
- **Working Features**: 14-step processing pipeline, speckle filtering, terrain correction
- **Critical Issues**: Georeferencing pixel size calculation bug (73,000x error)
- **Validation Status**: Partial - coordinates extraction validated, geometric accuracy compromised

---

## Scientific Accuracy Requirements

### **Geometric Accuracy Standards**
- **Target Resolution**: 20m pixel spacing for operational products
- **Coordinate Accuracy**: <1 pixel geometric error (better than 20m)
- **Validation Method**: Comparison against ESA SNAP software and ground control points
- **Acceptance Criteria**: RMS error <10m for well-defined targets

### **Radiometric Accuracy Standards**
- **Calibration Accuracy**: <0.5dB absolute radiometric accuracy
- **Speckle Preservation**: Maintain natural speckle statistics
- **Dynamic Range**: Support full Sentinel-1 dynamic range (16-bit input)
- **Validation Method**: Corner reflector analysis and cross-calibration

### **Scientific Compliance**
- **ESA Standards**: Full compliance with Sentinel-1 User Handbook specifications
- **CEOS Standards**: Follow Committee on Earth Observation Satellites guidelines
- **Reproducibility**: Bit-exact results for identical processing parameters
- **Traceability**: Full processing chain documentation and version control

---

## AI Collaboration Standards for SARdine

### **Mandatory Session Initialization Protocol**

Every AI collaboration session MUST begin with:

```
SARDINE SCIENTIFIC SOFTWARE DEVELOPMENT SESSION

Project: SARdine SAR Processing Pipeline
Domain: Remote Sensing / Synthetic Aperture Radar
Users: Research scientists, graduate students, operational analysts
Quality Standard: Research-grade, publication-ready

AI Assistant Requirements:
1. Read and acknowledge Scientific Software Development Protocol
2. Follow systematic analysis with evidence-based claims
3. Declare analysis scope before any technical work
4. Provide specific file paths and line numbers for all claims
5. State limitations and confidence levels
6. No code modifications without explicit approval
7. Prioritize scientific accuracy over development speed

Current Critical Issues:
- Georeferencing pixel size calculation bug (factor 73,000x error)
- Mathematical validation required for all algorithms
- ESA Sentinel-1 standard compliance verification needed

Session Objective: [Specify specific goal]
Success Criteria: [Define measurable outcomes]

Please confirm understanding of these scientific software requirements.
```

### **Analysis Scope Declaration Requirements**

Before ANY technical analysis, provide:

```
SARDINE ANALYSIS SCOPE DECLARATION:
- Target: [Specific component/function being analyzed]
- Scientific Context: [SAR processing step and scientific purpose]
- Methodology: [Exact tools and validation approach]
- Coverage: [What percentage will be examined]
- Limitations: [What will NOT be verified]
- ESA Standard Reference: [Relevant section of Sentinel-1 documentation]
- Validation Criteria: [How accuracy will be assessed]
- Confidence Level: [High/Medium/Low with justification]
```

### **Evidence Requirements for SAR Processing Claims**

All technical statements must include:
- **File paths and line numbers**: Exact code locations
- **Mathematical basis**: Equations and literature references
- **ESA compliance**: Reference to official Sentinel-1 standards
- **Validation status**: What has/hasn't been tested
- **Scientific limitations**: Assumptions and applicable ranges

---

## SAR-Specific Validation Standards

### **Algorithm Validation Requirements**

#### **TOPSAR Deburst Processing**
- **Literature Reference**: ESA Sentinel-1 Level 1 Detailed Algorithm Definition
- **Validation Method**: Phase continuity across burst boundaries
- **Test Cases**: Multi-burst IW scenes with overlap analysis
- **Acceptance Criteria**: <π/4 phase discontinuity at burst edges

#### **Radiometric Calibration**
- **Mathematical Basis**: σ⁰ = |DN|² / (LUT_σ⁰)²
- **Literature Reference**: ESA Sentinel-1 User Handbook Section 2.3.4
- **Validation Method**: Corner reflector analysis
- **Test Cases**: Known RCS targets and Amazon rainforest
- **Acceptance Criteria**: <0.5dB deviation from reference

#### **Terrain Correction (CRITICAL - CURRENTLY BUGGY)**
- **Mathematical Basis**: Range-Doppler geocoding with DEM
- **Literature Reference**: Small & Schubert, "Guide to SAR Geocoding"
- **Current Issue**: Pixel size calculation wrong by factor ~73,000x
- **Validation Method**: Ground control point analysis
- **Test Cases**: Mountainous terrain with known coordinates
- **Acceptance Criteria**: <20m RMS geometric error

#### **Speckle Filtering**
- **Mathematical Basis**: Gamma MAP filter implementation
- **Literature Reference**: Lopes et al., "Adaptive Speckle Filters"
- **Validation Method**: Speckle statistics preservation
- **Test Cases**: Homogeneous areas with known speckle parameters
- **Acceptance Criteria**: ENL within 10% of theoretical value

### **Data Quality Validation**

#### **Input Data Validation**
- **Sentinel-1 Format**: SLC product compliance checking
- **Orbit Data**: Precise orbit validation against ESA files
- **DEM Quality**: SRTM consistency and void handling
- **Metadata Integrity**: Annotation file parsing accuracy

#### **Output Data Validation**
- **Coordinate System**: WGS84 UTM zone verification
- **Pixel Size**: Target resolution achievement (20m)
- **Data Range**: Physical value ranges for backscatter
- **Metadata**: Complete processing provenance

---

## Critical Issues and Resolution Plan

### **Priority 1: Georeferencing Bug (CRITICAL)**

**Issue**: Pixel size calculation produces 0.0000000024° instead of ~0.00018° for 20m resolution

**Impact**: All geocoded outputs have wrong geographic scaling (factor 73,000x error)

**Root Cause**: Terrain correction function using image dimensions instead of target resolution

**Resolution Plan**:
1. Locate pixel size calculation in Rust terrain_correction function
2. Verify mathematical derivation against geocoding literature
3. Implement resolution-based calculation instead of dimension-based
4. Validate against known ground control points
5. Test with multiple scenes and terrain types

**Validation Criteria**: 
- Pixel size matches target resolution ±5%
- Geographic coordinates accurate to <20m RMS
- Comparison with ESA SNAP processing results

### **Priority 2: Mathematical Function Validation**

**Scope**: Systematic review of all scientific algorithms

**Components**:
- Calibration equations against ESA standards
- Terrain correction mathematics
- Speckle filter implementations
- Coordinate transformations

**Validation Method**:
- Literature review and citation verification
- Mathematical derivation checking
- Comparison with reference implementations
- Real data validation testing

### **Priority 3: Fallback Mechanism Removal**

**Issue**: Fallback orbit mechanisms compromise scientific accuracy

**Resolution**: Disable all fallback modes in scientific processing

**Validation**: Ensure processing fails gracefully when precise data unavailable

---

## Quality Assurance Framework

### **Automated Testing Requirements**

#### **Unit Tests (Required for all scientific functions)**
```rust
#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_radiometric_calibration_known_values() {
        // Test against ESA-provided reference data
        let input_dn = 1000.0;
        let lut_value = 12345.67;
        let expected_sigma0 = 0.00812; // From ESA reference
        
        let result = apply_radiometric_calibration(input_dn, lut_value)?;
        assert!((result - expected_sigma0).abs() < 1e-5);
    }
    
    #[test]
    fn test_geocoding_pixel_size() {
        // Test pixel size calculation for 20m target resolution
        let target_resolution = 20.0; // meters
        let expected_pixel_size = 0.00018; // degrees at mid-latitudes
        
        let result = calculate_pixel_size(target_resolution, 49.5)?; // Test latitude
        assert!((result - expected_pixel_size).abs() < 1e-6);
    }
}
```

#### **Integration Tests**
- End-to-end processing with known input/output pairs
- Multi-scene processing consistency
- Memory usage and performance monitoring
- Cross-platform compatibility verification

#### **Validation Tests**
- Comparison with ESA SNAP processing results
- Ground truth validation with field measurements
- Inter-comparison with other SAR processors
- Scientific literature benchmark reproduction

### **Documentation Standards**

#### **Mathematical Documentation Required**
Every scientific function must include:
```rust
/// Applies terrain flattening using local incidence angle normalization
/// 
/// # Mathematical Basis
/// σ⁰_flat = σ⁰_obs * sin(θ_local) / sin(θ_ref)
/// where θ_local is local incidence angle, θ_ref is reference angle
/// 
/// # Literature References
/// - Ulaby et al., "Microwave Remote Sensing", Vol 2, Chapter 12
/// - Small, "Flattening Gamma: Radiometric Terrain Correction", IEEE TGRS 2011
/// 
/// # ESA Compliance
/// - Follows Sentinel-1 Level 2 specification Section 4.3.2
/// 
/// # Limitations
/// - Assumes single-scattering approximation
/// - Valid for incidence angles 20° to 45°
/// - DEM resolution should be comparable to SAR resolution
/// 
/// # Error Propagation
/// δσ⁰_flat/σ⁰_flat ≈ δσ⁰_obs/σ⁰_obs + cot(θ_local)*δθ_local
```

#### **User Documentation**
- Installation and setup guides
- Processing workflow examples
- Validation result documentation
- Troubleshooting and FAQ
- Citation and attribution guidelines

### **Release and Validation Protocol**

#### **Pre-Release Checklist**
- [ ] All unit tests passing
- [ ] Integration tests with multiple scenes
- [ ] Validation against reference datasets
- [ ] Documentation completeness review
- [ ] Performance benchmarking
- [ ] Expert scientific review
- [ ] User acceptance testing

#### **Release Validation**
- [ ] ESA SNAP comparison report
- [ ] Ground truth validation results
- [ ] Scientific accuracy assessment
- [ ] Performance and memory analysis
- [ ] Cross-platform testing results
- [ ] User feedback integration

---

## Long-term Sustainability Plan

### **Scientific Community Engagement**
- Regular presentations at SAR/remote sensing conferences
- Collaboration with ESA validation teams
- Integration with educational curricula
- Open source community building

### **Maintenance and Evolution**
- Quarterly validation updates
- Annual algorithm review cycle
- Version control and change documentation
- Long-term support planning

### **Crisis Management Protocol**
If critical bugs are discovered in released versions:
1. **Impact Assessment**: Evaluate effect on published results
2. **User Notification**: Contact all known users immediately
3. **Correction Distribution**: Provide fixes and validation tools
4. **Documentation**: Update all documentation and examples
5. **Prevention Analysis**: Implement additional safeguards

---

## Success Metrics

### **Scientific Quality Indicators**
- **Geometric Accuracy**: <20m RMS error vs. ground control points
- **Radiometric Accuracy**: <0.5dB vs. ESA SNAP reference
- **Processing Consistency**: <1% variance in repeated processing
- **Literature Compliance**: All algorithms match published standards

### **Community Impact Indicators**
- **Scientific Publications**: Citations in peer-reviewed papers
- **User Adoption**: Download and usage statistics
- **Educational Use**: Integration in university courses
- **Expert Validation**: Endorsements from SAR community leaders

### **Technical Quality Indicators**
- **Code Coverage**: >90% for all scientific functions
- **Documentation Completeness**: All functions fully documented
- **Performance**: Processing time competitive with commercial tools
- **Reliability**: <1% failure rate on valid input data

---

## Commitment Statement

**SARdine prioritizes scientific integrity and accuracy over development speed and convenience.**

We commit to:
- **Never compromising scientific accuracy for performance or ease of use**
- **Maintaining transparent documentation of all limitations and assumptions**
- **Providing complete validation evidence for all scientific claims**
- **Supporting reproducible research through reliable, well-documented software**
- **Engaging with the scientific community for continuous validation and improvement**

**Every line of code in SARdine carries the responsibility of supporting real-world research that may influence environmental policy, climate understanding, and scientific discovery.**

---

## Project Contact and Governance

**Primary Maintainer**: [Your contact information]
**Scientific Advisory**: [Expert reviewers and advisors]
**User Community**: [Support channels and forums]
**Version Control**: [Repository and release management]

**This document serves as the definitive guide for all SARdine development activities and must be referenced in every AI collaboration session to ensure consistent scientific-grade development standards.**

---

*Last Updated: September 17, 2025*
*Document Version: 1.0*
*Next Review Date: December 17, 2025*