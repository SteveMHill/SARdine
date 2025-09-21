# SARdine Scientific Documentation
## Mathematical Reference for SAR Processing Algorithms

### Table of Contents
1. [Range-Doppler Geocoding](#range-doppler-geocoding)
2. [Terrain Correction](#terrain-correction)
3. [Coordinate System Transformations](#coordinate-transformations)
4. [Orbit State Vector Interpolation](#orbit-interpolation)
5. [Scientific References](#references)

---

## Range-Doppler Geocoding

### Mathematical Foundation

The Range-Doppler geocoding method is the standard approach for SAR image geolocation. It solves the intersection of two surfaces:
- **Range surface**: constant distance from satellite to ground point
- **Doppler surface**: constant Doppler frequency

### Range Equation

The range equation relates the two-way travel time to slant range:

```
R = c·τ/2
```

Where:
- `R` = slant range (meters)
- `c` = speed of light (299,792,458 m/s)
- `τ` = two-way travel time (seconds)

The range pixel index is calculated as:

```
n_r = (τ - τ₀) / Δτ
```

Where:
- `n_r` = range pixel index
- `τ₀` = two-way travel time to first range sample
- `Δτ` = range sampling interval = 2·Δr/c
- `Δr` = range pixel spacing

### Doppler Equation

The Doppler frequency equation for SAR is:

```
f_d = -2·(v⃗ · r̂) / λ
```

Where:
- `f_d` = Doppler frequency (Hz)
- `v⃗` = satellite velocity vector (m/s)
- `r̂` = unit range vector from satellite to target
- `λ` = radar wavelength (meters)

### Zero-Doppler Condition

For geocoding, we solve for the zero-Doppler condition:

```
f_d = 0  ⟹  v⃗ · r̂ = 0
```

This is solved using Newton-Raphson iteration:

```
t_{k+1} = t_k - f_d(t_k) / f'_d(t_k)
```

Where the derivative is:

```
f'_d = -2/λ · d(v⃗ · r̂)/dt
```

**Reference**: Cumming & Wong (2005), Chapter 4, "SAR Signal Processing"

---

## Terrain Correction

### Local Incidence Angle Correction

Terrain flattening corrects for the local topographic effects on backscatter:

```
σ⁰ = σ_orig · cos(θ_local) / cos(θ_ellipsoid)
```

Where:
- `σ⁰` = terrain-flattened backscatter coefficient
- `σ_orig` = original backscatter coefficient
- `θ_local` = local incidence angle
- `θ_ellipsoid` = ellipsoidal incidence angle

### Local Incidence Angle Calculation

The local incidence angle is computed using the surface normal:

```
θ_local = arccos(n⃗ · ŝ)
```

Where:
- `n⃗` = local surface normal vector
- `ŝ` = unit vector from target to satellite

### Surface Normal from DEM

The surface normal is calculated from DEM gradients:

```
n⃗ = normalize([-∂h/∂x, -∂h/∂y, 1])
```

Where:
- `h` = elevation from DEM
- `∂h/∂x`, `∂h/∂y` = elevation gradients

**Reference**: Small & Schubert (2019), "Guide to ALOS PALSAR Interferometry", ESA TM-19

---

## Coordinate System Transformations

### WGS84 to ECEF Conversion

The transformation from geodetic coordinates to Earth-Centered Earth-Fixed (ECEF):

```
X = (N + h) · cos(φ) · cos(λ)
Y = (N + h) · cos(φ) · sin(λ)
Z = (N(1-e²) + h) · sin(φ)
```

Where:
- `φ` = latitude (radians)
- `λ` = longitude (radians)
- `h` = ellipsoidal height (meters)
- `N` = prime vertical radius of curvature
- `e²` = first eccentricity squared

### Prime Vertical Radius

```
N = a / √(1 - e² · sin²(φ))
```

Where:
- `a` = semi-major axis (6,378,137 m for WGS84)
- `e²` = 0.00669437999014 for WGS84

**Reference**: Hofmann-Wellenhof et al. (2008), "GNSS - Global Navigation Satellite Systems"

---

## Orbit State Vector Interpolation

### Cubic Spline Interpolation

For precise orbit interpolation, we use Catmull-Rom cubic splines:

```
f(t) = 0.5 · [(2P₁) + (-P₀ + P₂)t + (2P₀ - 5P₁ + 4P₂ - P₃)t² + (-P₀ + 3P₁ - 3P₂ + P₃)t³]
```

Where:
- `P₀, P₁, P₂, P₃` = four consecutive orbit state vectors
- `t` = interpolation parameter [0,1]

This provides smooth first and second derivatives for accurate Doppler calculations.

**Reference**: Schaub & Junkins (2003), "Analytical Mechanics of Space Systems"

---

## Scientific References

### Primary SAR Processing References

1. **Cumming, I. G., & Wong, F. H. (2005)**  
   *Digital Processing of Synthetic Aperture Radar Data*  
   Artech House, Boston  
   ISBN: 978-1-58053-058-3

2. **Small, D., & Schubert, A. (2019)**  
   *Guide to ALOS PALSAR Interferometry*  
   ESA Technical Memorandum TM-19  
   Remote Sensing Laboratories, University of Zurich

3. **Hanssen, R. F. (2001)**  
   *Radar Interferometry: Data Interpretation and Error Analysis*  
   Kluwer Academic Publishers  
   ISBN: 978-0-7923-6945-5

### Coordinate Systems and Geodesy

4. **Hofmann-Wellenhof, B., Lichtenegger, H., & Wasle, E. (2008)**  
   *GNSS - Global Navigation Satellite Systems*  
   Springer-Verlag  
   ISBN: 978-3-211-73012-6

5. **Torge, W., & Müller, J. (2012)**  
   *Geodesy*  
   4th Edition, De Gruyter  
   ISBN: 978-3-11-020718-7

### Orbital Mechanics

6. **Schaub, H., & Junkins, J. L. (2003)**  
   *Analytical Mechanics of Space Systems*  
   AIAA Education Series  
   ISBN: 978-1-56347-563-1

### SAR Mission Documentation

7. **ESA Sentinel-1 User Handbook (2013)**  
   European Space Agency  
   ESA-GMES-S1OP-EOPG-TN-13-0001

8. **Torres, R., et al. (2012)**  
   *GMES Sentinel-1 mission*  
   Remote Sensing of Environment, 120, 9-24  
   DOI: 10.1016/j.rse.2011.05.028

---

## Implementation Standards

### Numerical Precision Requirements

- **Doppler frequency convergence**: 1 × 10⁻⁶ Hz
- **Orbit interpolation**: Sub-meter accuracy
- **Coordinate transformations**: Millimeter precision
- **Range timing**: Nanosecond accuracy

### Validation Requirements

All implementations must be validated against:
1. **Known test cases** with analytical solutions
2. **Established SAR processors** (SNAP, GAMMA, etc.)
3. **Peer-reviewed scientific literature**
4. **Sensor specification documents**

### Error Handling

- **No fallback approximations** - fail explicitly when conditions are not met
- **No synthetic data generation** - use only real sensor data
- **Comprehensive validation** - check all physical constraints
- **Traceable errors** - clear error messages with scientific context

---

*This documentation provides the mathematical foundation for all SAR processing algorithms implemented in SARdine. Each equation and method is referenced to peer-reviewed scientific literature to ensure scientific integrity and reproducibility.*
