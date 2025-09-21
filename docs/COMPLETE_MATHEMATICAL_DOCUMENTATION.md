# Complete Mathematical Documentation for SARdine
## Scientific SAR Processing Pipeline - Mathematical Foundations

---

## 1. Introduction and Scope

This document provides comprehensive mathematical documentation for all algorithms implemented in the SARdine SAR processing pipeline. SARdine processes Sentinel-1 Single Look Complex (SLC) data through a scientifically rigorous pipeline that maintains mathematical accuracy throughout all processing steps.

### 1.1 Mathematical Precision Requirements
- All calculations use double-precision floating point (f64) for maximum accuracy
- Physical constants are defined to full CODATA 2018 precision
- Coordinate transformations maintain sub-meter accuracy
- Temporal calculations preserve nanosecond precision

### 1.2 Coordinate Systems
SARdine operates with three primary coordinate systems:
1. **Radar Coordinates** (range, azimuth) - Native SAR acquisition geometry
2. **Zero-Doppler Coordinates** - Range-Doppler geocoding reference
3. **Geographic Coordinates** (WGS84) - Global reference system

---

## 2. Fundamental Constants and Parameters

### 2.1 Physical Constants (CODATA 2018)
```rust
/// Speed of light in vacuum (exact definition)
const SPEED_OF_LIGHT_M_S: f64 = 299_792_458.0; // m/s

/// WGS84 Earth ellipsoid parameters
const WGS84_SEMI_MAJOR_AXIS_M: f64 = 6_378_137.0;           // m
const WGS84_SEMI_MINOR_AXIS_M: f64 = 6_356_752.314245;      // m
const WGS84_ECCENTRICITY_SQUARED: f64 = 6.69437999014e-3;   // dimensionless
const WGS84_FLATTENING: f64 = 1.0 / 298.257223563;         // dimensionless
```

### 2.2 Sentinel-1 System Parameters
```rust
/// Sentinel-1 C-band center frequency (extracted from annotation XML)
/// Typical value: 5.405000454334350e9 Hz
/// Wavelength λ = c/f ≈ 0.055465763 m

/// Sentinel-1 IW mode parameters (extracted from annotation XML):
/// - Range pixel spacing: 2.329561... m (varies by subswath)
/// - Azimuth pixel spacing: 13.86708... m (varies by processing)
/// - Pulse repetition frequency: ~486 Hz (varies by mode)
```

---

## 3. Range-Doppler Coordinate Transformations

### 3.1 Zero-Doppler Time Calculation

The zero-Doppler time t₀ for a ground point P is the time when the satellite-to-point vector is perpendicular to the satellite velocity vector:

```
(R⃗ₛₐₜ(t₀) - P⃗) · V⃗ₛₐₜ(t₀) = 0
```

Where:
- `R⃗ₛₐₜ(t₀)` = satellite position at zero-Doppler time
- `P⃗` = target point position
- `V⃗ₛₐₜ(t₀)` = satellite velocity at zero-Doppler time

### 3.2 Slant Range Calculation

The slant range from satellite to ground point:

```
ρ = |R⃗ₛₐₜ(t₀) - P⃗|
```

### 3.3 Range-Doppler Geocoding Equations

Forward transformation (radar → geographic):

```
r = (2 × ρ / c - t_near) / Δt_range
a = (t₀ - t_azimuth_start) / Δt_azimuth
```

Where:
- `r` = range pixel coordinate
- `a` = azimuth pixel coordinate  
- `ρ` = slant range (m)
- `c` = speed of light (m/s)
- `t_near` = near range time (s)
- `Δt_range` = range sampling interval (s)
- `t₀` = zero-Doppler time (s)
- `t_azimuth_start` = azimuth start time (s)
- `Δt_azimuth` = azimuth time interval (s)

Inverse transformation (geographic → radar):

```
ρ = c × (t_near + r × Δt_range) / 2
t₀ = t_azimuth_start + a × Δt_azimuth
```

---

## 4. Orbit Interpolation Mathematics

### 4.1 Cubic Spline Interpolation

SARdine uses cubic spline interpolation for satellite position and velocity:

```
P(t) = a₀ + a₁t + a₂t² + a₃t³
V(t) = a₁ + 2a₂t + 3a₃t²
```

Coefficients are calculated using the constraint:
```
P(tᵢ) = Pᵢ  and  V(tᵢ) = Vᵢ
```

### 4.2 Orbit State Vector Interpolation

For high-precision geocoding, satellite position accuracy must be:
- Position: < 5 cm RMS
- Velocity: < 0.5 mm/s RMS

---

## 5. Terrain Correction Mathematics

### 5.1 Local Incidence Angle Calculation

The local incidence angle θᵢ between radar look vector and surface normal:

```
cos(θᵢ) = L⃗ · N⃗ / (|L⃗| × |N⃗|)
```

Where:
- `L⃗` = radar look vector (satellite to ground point)
- `N⃗` = surface normal vector from DEM

### 5.2 Surface Normal from DEM

Surface normal calculation using finite differences:

```
∂z/∂x = (z[i+1,j] - z[i-1,j]) / (2 × Δx)
∂z/∂y = (z[i,j+1] - z[i,j-1]) / (2 × Δy)

N⃗ = (-∂z/∂x, -∂z/∂y, 1) / √(1 + (∂z/∂x)² + (∂z/∂y)²)
```

### 5.3 Radiometric Terrain Correction

Corrected backscatter coefficient:

```
σ₀_corrected = σ₀_observed × cos(θᵢ) / cos(θᵣ)
```

Where:
- `θᵢ` = local incidence angle
- `θᵣ` = reference incidence angle (typically sensor boresight)

---

## 6. Calibration Mathematics

### 6.1 Radiometric Calibration

Convert digital numbers to calibrated backscatter:

```
σ₀ = (DN² - noise) / LUT²
```

Where:
- `DN` = digital number (complex amplitude)
- `noise` = thermal noise power
- `LUT` = calibration lookup table value

### 6.2 Complex Data Calibration

For complex SLC data:

```
γ₀ = |SLC|² / (A² × LUT²)
```

Where:
- `SLC` = complex pixel value
- `A` = antenna pattern compensation
- `LUT` = calibration coefficient

---

## 7. Multilooking Mathematics

### 7.1 Spatial Averaging

Multilooking reduces speckle through spatial averaging:

```
I_ML = (1/(Nr × Na)) × Σᵢ₌₀^(Nr-1) Σⱼ₌₀^(Na-1) |SLC[i,j]|²
```

Where:
- `Nr` = number of range looks
- `Na` = number of azimuth looks
- `I_ML` = multilooked intensity

### 7.2 Equivalent Number of Looks

Effective number of independent samples:

```
ENL = μ² / σ²
```

Where:
- `μ` = mean intensity
- `σ²` = intensity variance

---

## 8. Speckle Filtering Mathematics

### 8.1 Lee Filter

Adaptive speckle filter preserving edges:

```
I_filtered = μ + k × (I - μ)
```

Where:
```
k = (1 - Cᵤ²) / (1 + Cᵤ²)
Cᵤ² = σ²/μ² (coefficient of variation)
```

### 8.2 Gamma MAP Filter

Maximum a posteriori estimator:

```
I_filtered = μ × (α / (α + 1))
```

Where `α` is estimated from local statistics.

---

## 9. TOPSAR Burst Mode Processing

### 9.1 Azimuth Antenna Pattern Correction

TOPSAR azimuth antenna pattern:

```
A(η) = sinc²(β × η) × exp(-π × (η/θ₃dB)²)
```

Where:
- `η` = azimuth angle
- `β` = steering rate
- `θ₃dB` = 3dB beamwidth

### 9.2 Burst Synchronization

Precise timing for burst alignment:

```
t_burst = t_start + n × T_burst + δt_fine
```

Where:
- `T_burst` = burst repetition interval
- `δt_fine` = fine timing correction

---

## 10. Error Analysis and Accuracy

### 10.1 Geocoding Accuracy

Target geocoding accuracy:
- Horizontal: < 7m (90% confidence)  
- Vertical: < 10m (90% confidence)

### 10.2 Radiometric Accuracy

Target radiometric accuracy:
- Absolute: < 1 dB (σ₀)
- Relative: < 0.5 dB (σ₀)

### 10.3 Phase Accuracy

Complex data phase preservation:
- Phase stability: < 0.1 radians
- Coherence preservation: > 0.95

---

## 11. Quality Assessment Mathematics

### 11.1 Signal-to-Noise Ratio

```
SNR_dB = 10 × log₁₀(P_signal / P_noise)
```

### 11.2 Effective Scattering Area

```
ESA = σ₀ × A_pixel × sin(θ)
```

Where:
- `A_pixel` = pixel area in ground range
- `θ` = incidence angle

---

## 12. Implementation Notes

### 12.1 Numerical Stability
- All divisions check for near-zero denominators
- Trigonometric functions use stable implementations
- Matrix operations use condition number checking

### 12.2 Performance Considerations
- Vectorized operations where possible
- Memory-efficient chunked processing
- Parallel processing for independent operations

### 12.3 Validation Requirements
- All mathematical functions include unit tests
- Cross-validation against reference implementations
- Continuous integration testing

---

## References

This mathematical framework implements algorithms from:
1. ESA Sentinel-1 Product Definition Documents
2. "Digital Processing of Synthetic Aperture Radar Data" - Ian Cumming
3. "Synthetic Aperture Radar: Principles and Applications" - Soumekh
4. ESA SNAP Scientific Exploitation Tools
5. ASF DAAC SAR Processing Documentation

---

*Document Version: 1.0*  
*Last Updated: September 17, 2025*  
*Mathematical Review: Pending*