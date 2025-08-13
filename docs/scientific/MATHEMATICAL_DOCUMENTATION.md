# SARdine Mathematical Documentation
# Complete Equations and Derivations for SAR Processing

## Table of Contents
1. [SAR Fundamentals](#sar-fundamentals)
2. [Radiometric Calibration](#radiometric-calibration)
3. [Range-Doppler Geometry](#range-doppler-geometry)
4. [Terrain Flattening](#terrain-flattening)
5. [TOPSAR Processing](#topsar-processing)
6. [Coordinate Transformations](#coordinate-transformations)
7. [Statistical Processing](#statistical-processing)

---

## SAR Fundamentals

### Basic SAR Equation
The fundamental SAR equation relates the received signal to the target reflectivity:

```
s(τ,η) = ∫∫ σ(x,y) · w_r(τ - 2R(x,y)/c) · w_a(η - η_c(x,y)) dx dy
```

where:
- `s(τ,η)` = received SAR signal (range time τ, azimuth time η)
- `σ(x,y)` = radar cross section at ground position (x,y)
- `w_r(·)` = range window function (chirp pulse)
- `w_a(·)` = azimuth window function (synthetic aperture)
- `R(x,y)` = slant range to point (x,y)
- `c` = speed of light (299,792,458 m/s)
- `η_c(x,y)` = zero-Doppler time for point (x,y)

### Radar Cross Section
The radar cross section per unit area (σ⁰) is defined as:

```
σ⁰ = (4π · R²) / A · |s|² / |s_ref|²
```

where:
- `R` = slant range to target
- `A` = illuminated area on ground
- `|s|²` = measured signal power
- `|s_ref|²` = reference signal power

---

## Radiometric Calibration

### ESA Sentinel-1 Calibration Formula

The ESA-approved calibration formula implemented in SARdine:

```
σ⁰ = |DN|² / (LUT_σ⁰)²
```

**Detailed Components:**

1. **Digital Number (DN)**: Complex SLC value
   ```
   DN = I + jQ = |DN| · e^(jφ)
   |DN|² = I² + Q²
   ```

2. **Lookup Table Interpolation**: 
   ```
   LUT_σ⁰(r,a) = Σᵢ Σⱼ wᵢⱼ · LUT[i,j]
   ```
   where `wᵢⱼ` are bilinear interpolation weights:
   ```
   w₀₀ = (1-α)(1-β)
   w₀₁ = (1-α)β  
   w₁₀ = α(1-β)
   w₁₁ = αβ
   ```
   with `α = (r-r₀)/(r₁-r₀)` and `β = (a-a₀)/(a₁-a₀)`

3. **Calibration in dB**:
   ```
   σ⁰_dB = 10 · log₁₀(σ⁰) = 20 · log₁₀(|DN|) - 20 · log₁₀(LUT_σ⁰)
   ```

### Derivation from First Principles

Starting from the radar equation and antenna gain patterns:

```
P_r = (P_t · G_t · G_r · λ² · σ) / ((4π)³ · R⁴)
```

Converting to digital numbers through ADC quantization:
```
DN = √(P_r · G_system · K_cal)
```

Therefore:
```
σ = |DN|² / (G_t · G_r · λ² · G_system · K_cal) · ((4π)³ · R⁴) / P_t
```

The ESA calibration vectors contain the combined term:
```
LUT_σ⁰ = √(G_t · G_r · λ² · G_system · K_cal · P_t) / (4π)^(3/2) · R²
```

---

## Range-Doppler Geometry

### Fundamental Equations

1. **Slant Range Equation**:
   ```
   R(t) = c · τ / 2
   ```
   where `τ` is the two-way travel time.

2. **Doppler Frequency**:
   ```
   f_d(t) = (2/λ) · d/dt[R_vec(t) · v_rel(t) / |R_vec(t)|]
   ```
   
   Simplified for small angles:
   ```
   f_d(t) ≈ (2/λ) · (v_r(t))
   ```
   where `v_r(t)` is the radial velocity component.

3. **Zero-Doppler Condition**:
   ```
   f_d(t₀) = 0  ⟹  v_rel(t₀) ⊥ R_vec(t₀)
   ```

### Newton-Raphson Solution

To find zero-Doppler time `t₀`, we solve:
```
f_d(t) = 0
```

Newton-Raphson iteration:
```
t_{n+1} = t_n - f_d(t_n) / f'_d(t_n)
```

where the derivative is:
```
f'_d(t) = (2/λ) · d²/dt²[R_vec(t) · v_rel(t) / |R_vec(t)|]
```

**Convergence criteria**: |f_d(t_n)| < 1e-6 Hz

### Coordinate System Transformations

1. **Geographic to ECEF**:
   ```
   X = (N + h) · cos(φ) · cos(λ)
   Y = (N + h) · cos(φ) · sin(λ)
   Z = (N(1-e²) + h) · sin(φ)
   ```
   
   where:
   ```
   N = a / √(1 - e² · sin²(φ))  (radius of curvature)
   e² = 2f - f²                  (first eccentricity squared)
   ```
   
   WGS84 parameters:
   - `a = 6,378,137.0 m` (semi-major axis)
   - `f = 1/298.257223563` (flattening)

2. **ECEF to Range/Azimuth**:
   ```
   R = |P_target - P_satellite|
   
   θ_look = arccos((P_target - P_satellite) · v_satellite / 
                   (|P_target - P_satellite| · |v_satellite|))
   ```

---

## Terrain Flattening

### Basic Terrain Flattening Equation

```
γ⁰ = σ⁰ / cos(θ_lia)
```

where `θ_lia` is the local incidence angle.

### Local Incidence Angle Calculation

The local incidence angle is the angle between:
1. **Radar look vector**: `l̂ = (P_target - P_satellite) / |P_target - P_satellite|`
2. **Surface normal vector**: `n̂` computed from DEM

```
cos(θ_lia) = n̂ · l̂
```

### Surface Normal from DEM

Given elevation `z(x,y)`, the surface normal is:

```
n̂ = [-∂z/∂x, -∂z/∂y, 1] / √(1 + (∂z/∂x)² + (∂z/∂y)²)
```

**Gradient computation using central differences**:
```
∂z/∂x = (z[i,j+1] - z[i,j-1]) / (2·Δx)
∂z/∂y = (z[i+1,j] - z[i-1,j]) / (2·Δy)
```

**Edge handling**: Use forward/backward differences at boundaries:
```
∂z/∂x|_{x=x_min} = (z[i,j+1] - z[i,j]) / Δx     (forward)
∂z/∂x|_{x=x_max} = (z[i,j] - z[i,j-1]) / Δx     (backward)
```

### Physical Interpretation

- **σ⁰ (sigma-naught)**: Backscatter per unit area in slant range geometry
- **γ⁰ (gamma-naught)**: Backscatter per unit area in ground range geometry
- **β⁰ (beta-naught)**: Backscatter coefficient (not area-normalized)

The relationship between these:
```
β⁰ = σ⁰ · sin(θ_i) = γ⁰ · sin(θ_i) · cos(θ_lia)
```

where `θ_i` is the radar incidence angle.

### Terrain Correction Validation

**Valid angle ranges**:
- Local incidence angle: 10° < θ_lia < 80°
- Radar incidence angle: 20° < θ_i < 70° (Sentinel-1)

**Quality checks**:
```
if θ_lia < 10°:  # Too steep (forward slope)
    flag as unreliable
if θ_lia > 80°:  # Grazing angle (shadow)
    flag as unreliable
```

---

## TOPSAR Processing

### TOPSAR Burst Structure

Sentinel-1 TOPSAR mode acquires data in bursts with steering:

```
burst_duration = 2.758277 seconds  (typical)
burst_cycle = 3.2 seconds
```

Each burst has time-varying Doppler centroid:
```
f_dc(t) = f_dc0 + K_a · (t - t_ref)
```

### Azimuth Deramp Formula

```
s_deramp(η) = s_burst(η) · exp(-j · π · K_a · η²)
```

where:
- `s_burst(η)` = complex burst signal
- `K_a` = azimuth FM rate (from annotation)
- `η = t - t_burst_center` = azimuth time relative to burst center

**Detailed Implementation**:
```
K_a = 2 · PRF² / (V_s · L_antenna)
```
where:
- `PRF` = Pulse Repetition Frequency
- `V_s` = satellite velocity
- `L_antenna` = antenna length

### Burst Merging

To create seamless wide-swath images, overlapping burst regions must be merged:

```
s_merged(η) = w₁(η) · s_burst1(η) + w₂(η) · s_burst2(η)
```

**Weight functions** (linear taper):
```
w₁(η) = (η_end - η) / (η_end - η_start)  for η ∈ [η_start, η_end]
w₂(η) = (η - η_start) / (η_end - η_start)  for η ∈ [η_start, η_end]
```

### Phase Continuity

Ensure phase continuity between bursts:
```
Δφ = arg(s_burst2(η_overlap)) - arg(s_burst1(η_overlap))
```

If |Δφ| > π/4, apply phase correction:
```
s_burst2_corrected = s_burst2 · exp(-j · Δφ)
```

---

## Coordinate Transformations

### WGS84 Ellipsoid Parameters

```
a = 6,378,137.0 m           (semi-major axis)
b = 6,356,752.314245 m      (semi-minor axis)  
f = 1/298.257223563         (flattening)
e² = 0.00669437999014       (first eccentricity squared)
e'² = 0.00673949674228      (second eccentricity squared)
```

### Geographic to Cartesian (ECEF)

```
X = (N + h) · cos(φ) · cos(λ)
Y = (N + h) · cos(φ) · sin(λ)
Z = (N · (1 - e²) + h) · sin(φ)
```

where:
```
N = a / √(1 - e² · sin²(φ))  (prime vertical radius of curvature)
```

### Cartesian to Geographic (ECEF to WGS84)

Iterative solution for latitude:

```
p = √(X² + Y²)
λ = atan2(Y, X)

# Iterative latitude calculation
φ₀ = atan2(Z, p · (1 - e²))
for i in range(max_iterations):
    N = a / √(1 - e² · sin²(φᵢ))
    h = p / cos(φᵢ) - N
    φᵢ₊₁ = atan2(Z, p · (1 - e² · N/(N + h)))
    if |φᵢ₊₁ - φᵢ| < tolerance:
        break

h = p / cos(φ) - N
```

### UTM Projection

**Central meridian**: `λ₀ = -183° + 6° × zone_number`

**Scale factor**: `k₀ = 0.9996`

**False easting**: `E₀ = 500,000 m`

**False northing**: `N₀ = 0` (Northern hemisphere), `N₀ = 10,000,000 m` (Southern)

**Forward transformation**:
```
ν = a / √(1 - e² · sin²(φ))
ρ = a · (1 - e²) / (1 - e² · sin²(φ))^(3/2)
η² = e'² · cos²(φ)
t = tan(φ)
c = e'² · cos²(φ)
A = cos(φ) · (λ - λ₀)

T₁ = k₀ · ν · [A + (1-t²+c) · A³/6 + (5-18t²+t⁴+72c-58e'²) · A⁵/120]
T₂ = k₀ · [M + ν·t·(A²/2 + (5-t²+9c+4c²)·A⁴/24 + (61-58t²+t⁴+600c-330e'²)·A⁶/720)]

Easting = E₀ + T₁
Northing = N₀ + T₂
```

---

## Statistical Processing

### Speckle Statistics

SAR images exhibit multiplicative speckle noise. For `L` independent looks:

**Probability density function**:
```
p(I) = (L^L / Γ(L)) · (I/⟨I⟩)^(L-1) · (1/⟨I⟩) · exp(-L·I/⟨I⟩)
```

**Moments**:
```
⟨I⟩ = mean intensity
var(I) = ⟨I⟩² / L
std(I) = ⟨I⟩ / √L
CV = std(I) / ⟨I⟩ = 1 / √L  (coefficient of variation)
```

### Lee Filter Mathematics

**Enhanced Lee Filter**:
```
Î = ⟨I⟩ + k · (I - ⟨I⟩)
```

where the Wiener-like coefficient is:
```
k = (1 - Cu²/Ci²) / (1 + Cu²)
```

**Coefficient of variation**:
```
Cu² = 1/L               (theoretical)
Ci² = var(I)/⟨I⟩²       (observed in local window)
```

**Refined Lee Filter**:
```
k = max(0, (Ci² - Cu²)/(Ci² + Cu²))
```

### Multilooking

**Incoherent averaging**:
```
I_ML = (1/L) · Σᵢ₌₁ᴸ |sᵢ|²
```

**Complex multilooking** (preserves phase):
```
s_ML = (1/L) · Σᵢ₌₁ᴸ sᵢ
```

**Effective number of looks**:
```
L_eff = ⟨I⟩² / var(I)
```

### Noise Equivalent Sigma Zero (NESZ)

```
NESZ = k_B · T_sys · F_noise / (P_tx · G · λ² · τ · Δf / (4π)³)
```

where:
- `k_B = 1.38 × 10⁻²³ J/K` (Boltzmann constant)
- `T_sys` = system noise temperature
- `F_noise` = noise figure
- `P_tx` = transmitted power
- `G` = antenna gain
- `τ` = pulse duration
- `Δf` = processed bandwidth

---

## Error Analysis and Validation

### Geometric Accuracy Assessment

**Root Mean Square Error**:
```
RMSE = √[(1/n) · Σᵢ₌₁ⁿ (x_measured,i - x_true,i)²]
```

**Circular Error Probable (CEP)**:
```
CEP₉₀ = 2.146 · √(σ_E² + σ_N²)
```
where `σ_E` and `σ_N` are standard deviations in easting and northing.

### Radiometric Accuracy

**Relative calibration accuracy**:
```
ε_rel = |σ⁰_measured - σ⁰_reference| / σ⁰_reference
```

**Absolute calibration accuracy** (dB):
```
ε_abs = 10 · log₁₀(σ⁰_measured / σ⁰_reference)
```

**Target specifications**:
- Relative accuracy: < 0.5 dB (1σ)
- Absolute accuracy: < 1.0 dB (1σ)
- Stability: < 0.3 dB/year

### Quality Flags

**Geometric quality indicators**:
```
QI_geo = {
    EXCELLENT:  RMSE < 0.5 pixels
    GOOD:       0.5 ≤ RMSE < 1.0 pixels  
    ACCEPTABLE: 1.0 ≤ RMSE < 2.0 pixels
    POOR:       RMSE ≥ 2.0 pixels
}
```

**Radiometric quality indicators**:
```
QI_rad = {
    EXCELLENT:  |ε_abs| < 0.5 dB
    GOOD:       0.5 ≤ |ε_abs| < 1.0 dB
    ACCEPTABLE: 1.0 ≤ |ε_abs| < 2.0 dB  
    POOR:       |ε_abs| ≥ 2.0 dB
}
```

---

## Implementation Validation

### Unit Test Equations

Each mathematical formula should be validated against:

1. **Synthetic test cases** with known analytical solutions
2. **Reference implementations** (ESA SNAP, ASF MapReady)
3. **Published validation datasets** with ground truth

### Numerical Precision Requirements

- **Coordinate transformations**: ε < 1e-6 degrees
- **Range calculations**: ε < 1e-3 meters  
- **Doppler calculations**: ε < 1e-6 Hz
- **Calibration**: ε < 1e-3 dB
- **Angle calculations**: ε < 1e-4 radians

### Convergence Criteria

- **Newton-Raphson iterations**: max 50 iterations, tolerance 1e-6
- **DEM interpolation**: bilinear with bounds checking
- **Orbit interpolation**: 8th-order polynomial minimum
- **Phase unwrapping**: branch-cut algorithm with quality maps

---

This mathematical documentation provides the complete theoretical foundation for all SAR processing algorithms implemented in SARdine, ensuring research-grade accuracy and reproducibility.
