# 📐 MATHEMATICAL METHODS DOCUMENTATION
## SARdine SAR Processing Algorithms - Scientific Reference

**Document Version:** 1.0  
**Date:** August 13, 2025  
**Compliance:** Research-grade SAR processing standards

---

## 🔬 MATHEMATICAL FOUNDATIONS

### Physical Constants (CODATA 2018)
```rust
const SPEED_OF_LIGHT: f64 = 299_792_458.0; // m/s (exact)
const WGS84_SEMI_MAJOR_AXIS: f64 = 6_378_137.0; // m
const WGS84_FLATTENING: f64 = 1.0 / 298.257223563; // dimensionless
```

### Sentinel-1 Specifications
```rust
const S1_RADAR_FREQUENCY: f64 = 5.405e9; // Hz (C-band)
const S1_WAVELENGTH: f64 = 0.055465763; // m (calculated: c/f)
```

---

## 1. ORBITAL MECHANICS

### 1.1 State Vector Interpolation
**Reference:** Lagrange, J.L. (1795) + ESA S1-TN-MDA-52-7445

**Mathematical Foundation:**
```
Given: N state vectors (tᵢ, rᵢ, vᵢ) where i = 0...N-1
Target: State at time t

Lagrange Basis Functions:
Lⱼ(t) = ∏ₖ≠ⱼ [(t - tₖ) / (tⱼ - tₖ)]

Interpolated Position:
r(t) = Σⱼ rⱼ × Lⱼ(t)

Interpolated Velocity:
v(t) = Σⱼ vⱼ × Lⱼ(t)
```

**Implementation Requirements:**
- Use 8th-order interpolation (9 state vectors)
- Center interpolation points around target time
- Validate temporal bounds before interpolation

### 1.2 Satellite Velocity Calculation
**Reference:** Orbital Mechanics (Curtis, 2020)

**Mathematical Formula:**
```
For consecutive state vectors at times t₁, t₂:
Δr = r₂ - r₁
Δt = t₂ - t₁
v_avg = |Δr| / Δt

For higher accuracy, use centered differences:
v(tᵢ) = |r(tᵢ₊₁) - r(tᵢ₋₁)| / [2 × (tᵢ₊₁ - tᵢ₋₁)]
```

---

## 2. RADIOMETRIC CALIBRATION

### 2.1 Basic Calibration Equation
**Reference:** ESA S1-TN-MDA-52-7448

**Mathematical Foundation:**
```
Normalized Radar Cross Section:
σ⁰(i,j) = |DN(i,j)|² / [A² × LUT(i,j)]

Where:
- σ⁰ = normalized radar cross section (linear units)
- DN = digital number (complex SLC value)
- A = absolute calibration constant (dimensionless)
- LUT = calibration look-up table value
```

### 2.2 Bilinear Interpolation of Calibration Vectors
**Reference:** Numerical Recipes (Press et al., 2007)

**Mathematical Implementation:**
```
For calibration at fractional pixel (x, y):

1. Find bracketing calibration vectors in azimuth:
   line_lower = floor(y)
   line_upper = ceil(y)
   weight_y = y - line_lower

2. Interpolate in range for each vector:
   cal_lower = interpolate_range(x, vector[line_lower])
   cal_upper = interpolate_range(x, vector[line_upper])

3. Final bilinear interpolation:
   cal_factor = cal_lower × (1 - weight_y) + cal_upper × weight_y
```

---

## 3. MULTILOOKING

### 3.1 Spatial Averaging
**Reference:** Lee & Pottier (2009) "Polarimetric Radar Imaging"

**Mathematical Foundation:**
```
Input: Complex SLC data s(i,j)
Output: Multi-looked intensity I_ml(p,q)

Intensity Calculation:
I_ml(p,q) = (1/N) × Σₖ₌₀^(nr-1) Σₗ₌₀^(na-1) |s(p×nr+k, q×na+l)|²

Where:
- nr = number of range looks
- na = number of azimuth looks
- N = nr × na (total number of looks)
```

### 3.2 Output Pixel Spacing
**Mathematical Formula:**
```
Output Range Spacing = Input Range Spacing × nr
Output Azimuth Spacing = Input Azimuth Spacing × na
```

**Scientific Validation:**
- Preserves radiometric accuracy
- Reduces speckle by factor √N
- Maintains geometric relationships

---

## 4. SPECKLE FILTERING

### 4.1 Lee Filter
**Reference:** Lee, J.S. (1980) IEEE Transactions on Geoscience and Remote Sensing

**Mathematical Foundation:**
```
Local Statistics:
μ = (1/N) × Σ I(i,j)                    (local mean)
σ² = (1/N) × Σ [I(i,j) - μ]²          (local variance)
Cv = σ/μ                               (coefficient of variation)

Filter Weight:
k = max(0, (1 - Cu²/Cv²))

Where:
Cu = 1/√L                              (theoretical CV for L looks)

Filtered Value:
I_filtered = μ + k × (I_center - μ)
```

### 4.2 Enhanced Lee Filter
**Reference:** Lopes et al. (1990)

**Mathematical Enhancement:**
```
Classification-based filtering with edge preservation:

1. Classify pixel environment:
   - Homogeneous: Cv < Cu
   - Heterogeneous: Cu ≤ Cv < Cmax
   - Point target: Cv ≥ Cmax

2. Apply appropriate filter:
   - Homogeneous: Standard Lee filter
   - Heterogeneous: Reduced smoothing
   - Point target: No filtering
```

### 4.3 Frost Filter
**Reference:** Frost et al. (1982)

**Mathematical Implementation:**
```
Exponential Weighting Function:
w(i,j) = exp(-A × d(i,j))

Where:
A = damping factor based on local statistics
d(i,j) = distance from center pixel

Filtered Value:
I_filtered = Σ w(i,j) × I(i,j) / Σ w(i,j)
```

---

## 5. TERRAIN CORRECTION

### 5.1 Local Incidence Angle Calculation
**Reference:** Ulaby & Long (2014) "Microwave Radar and Radiometric Remote Sensing"

**Mathematical Foundation:**
```
Surface Normal Calculation (Horn's Method):
∂z/∂x = [z(i,j+1) - z(i,j-1)] / (2×Δx)
∂z/∂y = [z(i+1,j) - z(i-1,j)] / (2×Δy)

Surface Normal Vector:
n̂ = [-∂z/∂x, -∂z/∂y, 1] / |[-∂z/∂x, -∂z/∂y, 1]|

Radar Look Vector:
l̂ = (P_ground - P_satellite) / |P_ground - P_satellite|

Local Incidence Angle:
θ_local = arccos(l̂ · n̂)
```

### 5.2 Terrain Flattening
**Reference:** Small & Schubert (2008)

**Mathematical Formula:**
```
Terrain Flattened Backscatter:
γ⁰ = σ⁰ × cos(θ_local) / cos(θ_ref)

Where:
- γ⁰ = terrain flattened backscatter
- σ⁰ = original backscatter
- θ_local = local incidence angle
- θ_ref = reference incidence angle (typically 30-40°)
```

### 5.3 Range-Doppler Geocoding
**Reference:** Franceschetti & Lanari (1999)

**Fundamental Equations:**
```
Range Equation:
R = (c/2) × (τ - τ₀)

Where:
- R = slant range distance (m)
- c = speed of light (m/s)
- τ = two-way travel time (s)
- τ₀ = time to first range sample (s)

Doppler Equation:
f_d = (2/λ) × (v⃗ · û_los)

Where:
- f_d = Doppler frequency (Hz)
- λ = radar wavelength (m)
- v⃗ = satellite velocity vector (m/s)
- û_los = unit line-of-sight vector
```

---

## 6. TOPSAR DEBURSTING

### 6.1 Azimuth Deramping
**Reference:** Meta et al. (2010), De Zan & Guarnieri (2006)

**Mathematical Foundation:**
```
TOPSAR Deramping Phase:
φ_deramp(t) = π × k_a × (t - t_ref)²

Where:
- k_a = azimuth FM rate (Hz/s)
- t = azimuth time (s)
- t_ref = burst center time (s)

Complex Deramping Factor:
D(t) = exp(-j × φ_deramp(t))

Deramped Signal:
s_deramped(t) = s_original(t) × D(t)
```

### 6.2 Burst Boundary Handling
**Mathematical Approach:**
```
Burst Overlap Processing:
1. Identify overlap regions between adjacent bursts
2. Apply smooth transition weighting:
   w(x) = 0.5 × [1 + cos(π × (x - x_center) / overlap_width)]
3. Combine overlapping areas:
   s_merged = w(x) × s_burst1 + (1 - w(x)) × s_burst2
```

---

## 7. COORDINATE TRANSFORMATIONS

### 7.1 WGS84 Geodetic to ECEF
**Reference:** NIMA TR8350.2 (2000)

**Mathematical Implementation:**
```
Earth Parameters:
a = 6,378,137.0 m                      (semi-major axis)
f = 1/298.257223563                    (flattening)
e² = 2f - f²                           (first eccentricity squared)

Prime Vertical Radius:
N = a / √(1 - e² × sin²φ)

ECEF Coordinates:
X = (N + h) × cos φ × cos λ
Y = (N + h) × cos φ × sin λ
Z = (N × (1 - e²) + h) × sin φ

Where:
- φ = latitude (radians)
- λ = longitude (radians)
- h = ellipsoidal height (m)
```

### 7.2 Map Projection Calculations
**Reference:** Snyder (1987) "Map Projections: A Working Manual"

**UTM Projection:**
```
For UTM coordinates from geographic:
k₀ = 0.9996                            (scale factor)
λ₀ = central meridian of zone
E₀ = 500,000 m                         (false easting)
N₀ = 0 m (NH) or 10,000,000 m (SH)     (false northing)

Easting: E = k₀ × N × [A + (1-T+C) × A³/6 + ...]
Northing: N = k₀ × [M + N × tan φ × A²/2 + ...]
```

---

## 8. QUALITY ASSESSMENT

### 8.1 Signal-to-Noise Ratio
**Reference:** IEEE Standard 686-2017

**Mathematical Definition:**
```
SNR_dB = 10 × log₁₀(P_signal / P_noise)

Equivalent Number of Looks:
ENL = μ² / σ²

Where μ and σ are mean and standard deviation of homogeneous area
```

### 8.2 Geometric Accuracy Assessment
**Mathematical Metrics:**
```
Root Mean Square Error:
RMSE = √[(1/N) × Σ(observed - expected)²]

Circular Error Probable (50% confidence):
CEP = 0.5887 × (σ_x + σ_y)

Where σ_x, σ_y are standard deviations in x and y directions
```

---

## 9. ERROR PROPAGATION

### 9.1 Uncertainty Propagation
**Reference:** Taylor (1997) "An Introduction to Error Analysis"

**Linear Error Propagation:**
```
For function f(x₁, x₂, ..., xₙ):
σ_f² = Σᵢ (∂f/∂xᵢ)² × σᵢ²

For correlated variables:
σ_f² = Σᵢ (∂f/∂xᵢ)² × σᵢ² + 2 × ΣᵢΣⱼ (∂f/∂xᵢ)(∂f/∂xⱼ) × σᵢⱼ
```

### 9.2 SAR-Specific Error Sources
```
Orbit Error: ±0.03 m (precise orbits)
Timing Error: ±0.1 ms
Calibration Error: ±0.3 dB
DEM Error: ±16 m (SRTM)
Atmospheric Delay: ±0.1 m
```

---

## 10. VALIDATION METHODS

### 10.1 Corner Reflector Analysis
**Reference:** IEEE Standard 686-2017

**Theoretical Response:**
```
RCS_theoretical = 4π × a⁴ / λ²

Where:
- a = corner reflector leg length
- λ = radar wavelength

Point Target Response:
PSLR = Peak Side Lobe Ratio
ISLR = Integrated Side Lobe Ratio
```

### 10.2 Cross-Platform Validation
**Comparison Metrics:**
```
Relative Difference:
Δ_rel = |value_test - value_ref| / value_ref × 100%

Acceptable Thresholds:
- Radiometric: <0.3 dB
- Geometric: <0.5 pixels
- Phase: <0.1 radians
```

---

## 📚 COMPLETE REFERENCE BIBLIOGRAPHY

### Essential SAR Processing Literature
1. **Franceschetti, G., & Lanari, R. (1999).** Synthetic aperture radar processing. CRC press.
2. **Curlander, J. C., & McDonough, R. N. (1991).** Synthetic aperture radar systems and signal processing.
3. **Lee, J. S., & Pottier, E. (2009).** Polarimetric radar imaging: from basics to applications.
4. **Ulaby, F. T., & Long, D. G. (2014).** Microwave radar and radiometric remote sensing.
5. **Small, D., & Schubert, A. (2008).** Guide to ASAR geocoding.

### TOPSAR and Sentinel-1 Specific
6. **Meta, A., et al. (2010).** TOPSAR: Terrain Observation by Progressive Scans.
7. **De Zan, F., & Guarnieri, A. M. (2006).** TOPSAR: Terrain observation by progressive scans.
8. **Torres, R., et al. (2012).** GMES Sentinel-1 mission.

### Mathematical and Geodetic References
9. **Horn, B. K. P. (1981).** Hill shading and the reflectance map.
10. **Snyder, J. P. (1987).** Map projections: A working manual.
11. **Press, W. H., et al. (2007).** Numerical recipes: The art of scientific computing.
12. **Taylor, J. R. (1997).** An Introduction to error analysis.

### Standards and Technical Documents
13. **ESA S1-TN-MDA-52-7448:** Sentinel-1 Radiometric Calibration
14. **ESA S1-TN-MDA-52-7445:** Sentinel-1 Orbit State Vectors
15. **ESA S1-RS-MDA-52-7440:** Sentinel-1 Product Specification
16. **IEEE Standard 686-2017:** Radar Definitions and Terminology

---

*This mathematical documentation provides the scientific foundation for all SAR processing algorithms implemented in SARdine, ensuring research-grade accuracy and reproducibility.*
