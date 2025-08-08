# SARdine Mathematical Methods Documentation

## Scientific Implementation References

### Range Doppler Geocoding

**Description:** Range-Doppler geocoding following ESA methodology

**Scientific Reference:** ESA Sentinel-1 Level 1 Detailed Algorithm Definition

**Mathematical Equation:** `R = c * τ/2, where τ is two-way travel time`

**Implementation:** Iterative solution of range-Doppler equations

**Required Parameters:** slant_range, doppler_frequency, orbit_state_vectors

---

### Topsar Debursting

**Description:** TOPSAR debursting with azimuth antenna pattern correction

**Scientific Reference:** Meta, Mittermayer, Steinbrecher (2010) - TOPSAR Signal Processing

**Mathematical Equation:** `Debursted = SLC * exp(-j*φ_antenna)`

**Implementation:** Burst boundary detection and azimuth pattern correction

**Required Parameters:** burst_start_time, burst_end_time, azimuth_fm_rate

---

### Radiometric Calibration

**Description:** Conversion to backscatter coefficient (σ⁰, γ⁰, β⁰)

**Scientific Reference:** ESA Radiometric Calibration of SAR Data (Rosich & Meadows)

**Mathematical Equation:** `σ⁰ = |DN|² / (A² * sin(θ_inc))`

**Implementation:** LUT-based calibration with incidence angle correction

**Required Parameters:** calibration_lut, incidence_angle, range_spreading_loss

---

### Terrain Correction

**Description:** Terrain correction using DEM and precise orbits

**Scientific Reference:** Small, Schubert (2008) - Guide to ASAR Geocoding

**Mathematical Equation:** `σ⁰_flat = σ⁰ * cos(θ_local) / cos(θ_ref)`

**Implementation:** DEM-based local incidence angle calculation

**Required Parameters:** dem_elevation, local_incidence_angle, reference_angle

---

### Orbit Interpolation

**Description:** Precise orbit interpolation using state vectors

**Scientific Reference:** Curlander & McDonough (1991) - SAR Data Processing

**Mathematical Equation:** `P(t) = Σ L_i(t) * P_i (Lagrange interpolation)`

**Implementation:** 8th-order Lagrange polynomial interpolation

**Required Parameters:** orbit_state_vectors, interpolation_time, polynomial_order

---

