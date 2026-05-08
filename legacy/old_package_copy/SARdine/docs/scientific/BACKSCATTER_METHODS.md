# Backscatter Processing: Methods and Equations

This document summarizes the scientific methods implemented in SARdine for Sentinel-1 IW SLC backscatter processing. It provides the core equations, assumptions, and references to the SAR literature. Variables follow ESA Sentinel-1 conventions where possible.

## Radiometric Calibration

Goal: convert raw complex SLC samples to calibrated backscatter in sigma0, beta0, or gamma0.

- Complex intensity: I = |S|^2, where S is the focused complex SLC sample.
- Calibration gain G(r, a) from annotation (per-pixel or per-line LUTs) is applied as a multiplicative factor.

- Beta nought (β0):
  β0 = I · G(r, a)

- Sigma nought (σ0):
  σ0 = β0 · sin(θ_i)
  
  where θ_i is the local incidence angle at the ground point.

- Gamma nought (γ0):
  γ0 = β0 · cos(θ_i)

- Decibel conversion:
  x_dB = 10 · log10(x)

Notes:
- Terrain flattening (RTC) reshapes β0 to γ0 using local incidence angle from DEM and orbit geometry.
- LUTs are interpolated in range/azimuth. We preserve radiometric power during resampling (see Multilooking below).

References:
- Cumming & Wong (2005), Digital Processing of Synthetic Aperture Radar Data, ch. 4–6
- ESA Sentinel-1 Product Specification (S1-RS-MDA-52-7441)
- Small (2011) Flattening Gamma: Radiometric Terrain Correction for SAR Imagery

## TOPS Deburst and IW Merge

Sentinel‑1 IW uses TOPS bursts. Deburst stitches bursts within each subswath; merge aligns IW1–IW3 using slant range timing offsets Δτ.

- Range alignment between subswaths uses near-range two-way travel time τ = 2R/c.
- Merged grid defined by reference subswath; other subswaths are shifted by Δsamples = (Δτ)/Δτ_sample, with Δτ_sample derived from native range pixel spacing.

We use per-burst azimuth timing from the annotation plus Doppler centroids to preserve phase continuity and geometry during merge.

References:
- De Zan & Guarnieri (2006) TOPSAR: Terrain Observation by Progressive Scans
- ESA Sentinel‑1 Level‑1 detailed algorithm definition

## Multilooking (Power-Preserving)

We multilook in linear power domain to preserve radiometry:

- Given a window Lr × La, multilooked power:
  P_mlk(r, a) = (1/(Lr·La)) · Σ_{i=0}^{Lr-1} Σ_{j=0}^{La-1} P(r+i, a+j)

- Ground-range spacing for Sentinel-1 IW SLC is derived from annotation slant spacing Δr_slant via
  Δr_ground = Δr_slant / sin(θ_incidence). We select range looks as
  L_r = clamp(round(r_target / Δr_ground), 1, 8) and azimuth looks as
  L_a = clamp(round(r_target / Δa_native), 1, 4), where r_target is the requested RTC ground
  sampling distance and Δa_native is the native azimuth spacing (~13.9 m). This follows the
  Sentinel-1 Product Specification S1-RS-MDA-52-7441 (IW SLC spacing) and the ASF RTC
  Algorithm Theoretical Basis Document (§3.3).

- Number of looks ENL is estimated to validate speckle reduction.

We avoid dB-domain averaging. Interpolation kernels for resampling are power-preserving (area-weighted where applicable).

References:
- Oliver & Quegan (2004) Understanding Synthetic Aperture Radar Images
- Bamler & Hartl (1998) Synthetic Aperture Radar Interferometry

## Range-Doppler Geocoding (Zero-Doppler Projection)

We map SAR image coordinates (range r, azimuth a) to geographic coordinates (ϕ, λ, h) by solving:

1. Zero-Doppler condition:
   f_D(t) = (λ^−1) · (v(t) · ρ(t)) = 0
   
   where v(t) is the satellite velocity, ρ(t) the look vector from satellite to ground point at time t, and λ the wavelength.

2. Range equation (two-way travel time):
   τ = 2 · |ρ(t)| / c

3. Terrain constraint: h from DEM, Earth model WGS84.

We solve for t (azimuth time) with a robust secant/Newton method with damping and bisection safeguards. The native line time per pixel is the annotation’s azimuthTimeInterval, not 1/PRF for TOPS merged data.

- Convert from relative product time to absolute UTC using orbit_ref_epoch_utc and product_start_rel_s.
- Convert slant range to native range pixel using τ and native range spacing.

References:
- Curlander & McDonough (1991) Synthetic Aperture Radar: Systems and Signal Processing
- Cumming & Wong (2005), ch. 7 (Geocoding)
- Bamler & Hartl (1998)

## Terrain Flattening (RTC)

We convert β0 to γ0 using local incidence angle θ_i from DEM and orbit geometry (surface normal n̂ and look vector l̂):

- cos θ_i = n̂ · l̂
- γ0 = β0 · cos θ_i

For masking, we use thresholds on cos θ_i and DEM validity.

References:
- Small (2011) Flattening Gamma
- Hoekman & Reiche (2015) Multi-Temporal Forest Monitoring Using Coherence and Backscatter

## Coordinate Systems and Units

- Earth model: WGS84 ellipsoid with latitude-aware meters-per-degree conversions
- Times: UTC seconds; orbit-relative times referenced via orbit_ref_epoch_utc + product_start_rel_s
- Pixel spacings: native range and azimuth from annotation; no hardcoded values
- Output grid: default projection resolves to the scene's UTM zone (EPSG:326xx/327xx) with polar
  stereographic fallbacks above 80° latitude; EPSG:4326 is provided as an optional quick-look export.

## Implementation Checks (Scientific Integrity)

- No synthetic orbit data; require .EOF precise or restituted files
- No hardcoded radar frequency; derive wavelength from annotation
- No 1/PRF substitution for azimuthTimeInterval in TOPS; use annotation value
- Strict mode (SARDINE_STRICT=1) forbids heuristic fallbacks

## Suggested Citations

- Curlander, J. C., & McDonough, R. N. (1991). Synthetic Aperture Radar: Systems and Signal Processing. Wiley.
- Cumming, I. G., & Wong, F. H. (2005). Digital Processing of Synthetic Aperture Radar Data. Artech House.
- Bamler, R., & Hartl, P. (1998). Synthetic Aperture Radar Interferometry. Inverse Problems, 14(4), R1–R54.
- Small, D. (2011). Flattening Gamma: Radiometric Terrain Correction for SAR Imagery. IEEE TGRS, 49(8), 3081–3093.
