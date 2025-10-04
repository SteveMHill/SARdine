# Scientific Algorithms in SARdine

## Mathematical Foundations and Literature References

This document provides comprehensive mathematical documentation for all scientific algorithms implemented in SARdine, with literature references and validation requirements.

---

## 1. Range-Doppler Geocoding (Terrain Correction)

### Mathematical Basis

The Range-Doppler geocoding algorithm transforms SAR slant-range coordinates to geographic coordinates using Digital Elevation Models (DEM).

#### Core Equations

**Zero-Doppler Condition:**
```
f_d(t) = -2 * (R⃗(t) · V⃗(t)) / (λ |R⃗(t)|) = 0
```

Where:
- `f_d(t)`: Doppler frequency at time t
- `R⃗(t) = X_target - X_sat(t)`: Range vector from satellite to target
- `V⃗(t)`: Satellite velocity vector at time t  
- `λ`: Radar wavelength
- `|R⃗(t)|`: Magnitude of range vector (slant range)

**Range Equation:**
```
R = |X_target - X_sat(t_azimuth)|
```

**Pixel Coordinate Mapping:**
```
range_pixel = (τ_computed - τ_0) / Δτ
azimuth_pixel = (t_azimuth - t_start) * PRF
```

Where:
- `τ_computed = 2R/c`: Two-way travel time
- `τ_0`: Slant range time to first pixel (from annotation)
- `Δτ = 2 * range_pixel_spacing / c`: Range pixel spacing in time
- `PRF`: Pulse Repetition Frequency

### Algorithm Steps

1. **Coordinate Transformation**: Convert (lat, lon, height) → ECEF coordinates
2. **Zero-Doppler Solution**: Solve f_d(t) = 0 using Newton-Raphson iteration
3. **Range Calculation**: Compute slant range from satellite to target
4. **Pixel Mapping**: Transform to SAR image coordinates

### Literature References

- **Bamler & Hartl (1998)**: "Synthetic aperture radar interferometry", *Inverse Problems* 14, R1-R54
- **Small & Schubert (2008)**: "Guide to SAR Geocoding", *Remote Sensing Tutorial*
- **ESA (2019)**: "Sentinel-1 Level 1 Detailed Algorithm Definition", GMES-S1OP-EOPG-TN-13-0007
- **Curlander & McDonough (1991)**: "Synthetic Aperture Radar Systems and Signal Processing"

### Implementation Files
- `src/core/terrain_correction.rs`: Main implementation
- `src/lib.rs`: Python interface (lines 2065-2315)

### Validation Requirements
- Geometric accuracy: <1 pixel RMS error vs known targets
- Comparison against ESA SNAP geocoding results
- Cross-validation with published corner reflector positions

---

## 2. Newton-Raphson Zero-Doppler Solver

### Mathematical Basis

Iterative solution of the zero-Doppler condition using Newton-Raphson method:

```
t_{n+1} = t_n - f(t_n) / f'(t_n)
```

Where:
```
f(t) = -2 * (R⃗(t) · V⃗(t)) / (λ |R⃗(t)|)
f'(t) = df/dt (computed via satellite acceleration)
```

### Convergence Properties
- **Initial Guess**: Closest approach from orbit state vectors
- **Tolerance**: |f(t)| < 1e-9 Hz
- **Typical Iterations**: 3-5 for most cases
- **Convergence Rate**: Quadratic (when well-conditioned)

### Literature References
- **Press et al. (2007)**: "Numerical Recipes: The Art of Scientific Computing"
- **Franceschetti & Lanari (1999)**: "Synthetic Aperture Radar Processing", CRC Press

### Implementation
- File: `src/core/terrain_correction.rs`, function `newton_raphson_zero_doppler`
- Lines: 1943-2050

---

## 3. Radiometric Calibration

### Mathematical Basis

Conversion of Digital Numbers (DN) to calibrated backscatter coefficients:

**Sigma-nought (σ⁰):**
```
σ⁰ = |DN|² / |LUT_σ⁰|²
```

**Beta-nought (β⁰):**
```
β⁰ = |DN|² / |LUT_β⁰|²
```

**Gamma-nought (γ⁰):**
```
γ⁰ = |DN|² / |LUT_γ⁰|²
```

### Unit Conversion Logic

The calibration system automatically detects and converts between dB and linear units:

**dB to Linear Conversion:**
```
linear_value = 10^(dB_value / 10)
```

**Detection Heuristics:**
- **Gamma LUT**: If median ∈ [10, 50] dB → convert to linear
- **Sigma/Beta LUT**: If median ∈ [5, 50] dB and non-integer → convert to linear
- **Explicit Units**: If XML units contain "db" → convert to linear

### Literature References
- **ESA (2019)**: "Sentinel-1 User Handbook", Section 2.3.4
- **Small (2011)**: "Flattening Gamma: Radiometric Terrain Correction for SAR Imagery"
- **Ulaby & Long (2014)**: "Microwave Radar and Radiometric Remote Sensing"

### Implementation Files
- `src/core/calibrate.rs`: Main calibration logic
- `src/io/slc_reader.rs`: XML parsing and LUT extraction

### Validation Requirements
- Absolute radiometric accuracy: <0.5 dB
- Cross-calibration with corner reflectors
- Validation against Amazon rainforest reference values

---

## 4. TOPSAR Deburst Processing

### Mathematical Basis

TOPSAR (Terrain Observation with Progressive Scans SAR) requires seamless merging of burst data with phase continuity.

**Phase Continuity Condition:**
```
Δφ = φ(burst_n+1, overlap_start) - φ(burst_n, overlap_end) < π/4
```

**Burst Timing:**
```
t_burst = t_start + (line_number / PRF)
azimuth_antenna_pattern = sinc²(β * (t - t_burst_center))
```

### Algorithm Steps
1. **Burst Identification**: Extract burst boundaries from annotation
2. **Overlap Analysis**: Identify overlap regions between bursts  
3. **Phase Alignment**: Ensure phase continuity across burst boundaries
4. **Seamless Merging**: Weighted combination in overlap regions

### Literature References
- **ESA (2013)**: "TOPS Sentinel-1 Level-1 Detailed Algorithm Definition"
- **Prats-Iraola et al. (2012)**: "TOPS Interferometry with TerraSAR-X"

### Implementation Files
- `src/core/deburst.rs`: Main deburst implementation
- `src/core/deburst_optimized.rs`: Optimized version

### Validation Requirements
- Phase discontinuity: <π/4 radians at burst boundaries
- Amplitude continuity: <10% variation in overlap regions

---

## 5. Speckle Filtering

### Mathematical Basis

Speckle filtering reduces multiplicative speckle noise while preserving structural information.

**Enhanced Lee Filter:**
```
I_filtered = I_mean + k * (I_center - I_mean)
k = (1 - C_u/C_i) * (C_i/C_u) / (1 + C_i/C_u)
```

Where:
- `C_i`: Coefficient of variation in local window
- `C_u`: Coefficient of variation for pure speckle (≈0.523)
- `I_mean`: Local window mean intensity

**Gamma MAP Filter:**
```
I_filtered = α * I_mean
α = (1 + 1/L) / (1 + (C_i²/C_u²))
```

Where `L` is the number of looks.

### Literature References
- **Lopes et al. (1993)**: "Adaptive Speckle Filters and Scene Heterogeneity"
- **Lee (1980)**: "Digital Image Enhancement and Noise Filtering by Use of Local Statistics"
- **Frost et al. (1982)**: "A Model for Radar Images and Its Application to Adaptive Digital Filtering"

### Implementation Files
- `src/core/speckle_filter.rs`: All speckle filtering algorithms

### Validation Requirements
- Speckle reduction: >50% noise variance reduction
- Edge preservation: <10% degradation of sharp edges
- Statistical validation: Maintain natural speckle statistics in homogeneous areas

---

## 6. Coordinate Transformations

### Mathematical Basis

**WGS84 Ellipsoid to ECEF:**
```
X = (N + h) * cos(φ) * cos(λ)
Y = (N + h) * cos(φ) * sin(λ)  
Z = (N * (1 - e²) + h) * sin(φ)
```

Where:
- `N = a / √(1 - e² sin²φ)`: Prime vertical radius
- `a = 6378137.0 m`: WGS84 semi-major axis
- `e² = 0.00669438002290`: First eccentricity squared
- `φ, λ`: Latitude, longitude (radians)
- `h`: Height above ellipsoid (meters)

**Pixel Size Calculation:**
```
pixel_size_degrees = resolution_meters / (N * cos(φ) * π / 180°)
```

### Literature References
- **Hofmann-Wellenhof et al. (2007)**: "GNSS - Global Navigation Satellite Systems"
- **NIMA (2000)**: "Department of Defense World Geodetic System 1984"

### Implementation Files
- `src/core/terrain_correction.rs`: Coordinate transformation functions

### Validation Requirements
- ECEF accuracy: <1 mm vs reference implementations
- Pixel size accuracy: <0.0001% relative error

---

## Validation Framework

### Scientific Testing Requirements

All algorithms must pass comprehensive validation tests:

1. **Mathematical Consistency**: Compare against analytical solutions
2. **Literature Benchmarks**: Cross-validate with published results  
3. **Reference Implementation**: Compare with ESA SNAP when possible
4. **Ground Truth**: Validate with known corner reflector positions
5. **Error Propagation**: Quantify uncertainty through processing chain

### Validation Implementation
- File: `src/core/terrain_correction.rs`
- Method: `validate_scientific_implementation()`
- Tests: 5 comprehensive validation categories

### Success Criteria
- Geometric accuracy: <1 pixel RMS error
- Radiometric accuracy: <0.5 dB absolute error  
- Phase continuity: <π/4 radian discontinuity
- Convergence rate: >95% for Newton-Raphson
- Coordinate accuracy: <1 mm transformation error

---

## References

### Primary Literature
1. Bamler, R., & Hartl, P. (1998). Synthetic aperture radar interferometry. *Inverse problems*, 14(4), R1.
2. Curlander, J. C., & McDonough, R. N. (1991). *Synthetic aperture radar: systems and signal processing*. Wiley.
3. Small, D. (2011). Flattening gamma: Radiometric terrain correction for SAR imagery. *IEEE TGRS*, 49(8), 3081-3093.
4. Lopes, A., Touzi, R., & Nezry, E. (1993). Adaptive speckle filters and scene heterogeneity. *IEEE TGRS*, 31(6), 1392-1404.

### ESA Documentation
- ESA (2019). Sentinel-1 User Handbook. ESA-GMES-S1OP-EOPG-TN-13-0001
- ESA (2019). Sentinel-1 Level 1 Detailed Algorithm Definition. GMES-S1OP-EOPG-TN-13-0007
- ESA (2013). TOPS Sentinel-1 Level-1 Detailed Algorithm Definition. ESA-EOPG-CSCOP-TN-0014

### Validation Standards
- Committee on Earth Observation Satellites (CEOS) SAR Calibration Standards
- IEEE Standard for SAR Terminology and Definitions (IEEE 686-2017)
- ISO 19115 Geographic Information Metadata Standards

---

*This document should be updated whenever new algorithms are implemented or existing ones are modified. All implementations must include literature references and validation requirements.*