
# Sentinel‑1 Metadata & Annotation Map (for SARdine pipeline)

This document summarizes **what to read** from Sentinel‑1 products, **which function/module uses it**, and **what derived values** each stage produces that should be passed downstream. It is tailored to the files/modules you shared (`slc_reader.rs`, `sentinel1.rs`, `annotation.rs`, `orbit.rs`, `calibrate.rs`, `deburst.rs`, `topsar_merge.rs`, `terrain_correction.rs`, `scientific_terrain_flatten.rs`, `dem.rs`, `iw_split_optimized.rs`, `multilook.rs`, `validation_gates.rs`, etc.).

> Scope: IW SLC/TOPS (VV/VH, dual or single-pol), but most fields generalize to SM/EW.  
> Time base: Treat `azimuth_time` from solvers as **relative to `orbit_ref_epoch`**, not product start; convert before indexing lines/bursts.

---

## 1) Product & Platform Core

**Read from**: `manifest.safe`, `annotation/s1?-iw?-slc-*.xml`, `measurement/*.tiff`, `calibration/calibration*.xml`, `noise/*.xml`

| Field | Why / Used by | Module/Function | Derived/Pass downstream |
|---|---|---|---|
| `missionId` / `platform` (S1A/S1B) | Consistency, antenna patterns | `sentinel1.rs` | — |
| `acquisitionMode` (IW/SM/EW) | Processing branch | `slc_reader.rs` / `iw_split_optimized.rs` | swath layout |
| `swath` (IW1/IW2/IW3) | Range grid specifics | `iw_split_optimized.rs` | swath → burst list |
| `polarisationChannels` (VV/VH) | Band selection, noise LUT | `slc_reader.rs` | selected pol |
| `lookDirection` (Right) & `passDirection` (ASC/DESC) | Geometry signs, heading estimate | `sentinel1.rs` | azimuth heading sign hint |
| `processingLevel` (SLC) | Choose SLC path | global | — |
| `productType` (SLC) | — | global | — |
| `imageDimensions` (lines, pixels) | Buffer sizes, bounds | `slc_reader.rs` | — |
| `sampleType` (COMPLEX) | Calibration path (SLC vs intensity) | `calibrate.rs` | — |

---

## 2) Timing Models (Line/Pixel Grids)

**Read from**: `annotation/*.xml` (`generalAnnotation`)

| Field | Why / Used by | Module/Function | Derived/Pass downstream |
|---|---|---|---|
| `startTime` (product) | Convert absolute times to line indices | `slc_reader.rs` | `t_product_start` |
| `azimuthTimeInterval` (PRI equiv) | line ⇄ time mapping | deburst, TC | `dt_az` |
| `rangeSamplingRate` | pixel ⇄ fast-time | multiple | `fs`, `dr = c/(2 fs)` |
| `lineTimePolynomial` / `azimuthTime` per line (if present) | more accurate mapping | TC solver | az time prediction |
| `slantRangeTime` of first pixel | base range | geocoding, TC | `τ0` |
| `pixelSpacingRange` (ELLIPSOID) | range spacing on ground proxies | multilook | spacing hints |

**Important conversions**:
- `t_abs = t_orbit_ref_epoch + azimuth_time` (solver)  
- `line = round((t_abs - t_product_start)/azimuthTimeInterval)`  
- `τ(pix) = τ0 + pix * (1/fs) * 2` ⇒ `R = c τ / 2`

---

## 3) Bursts (TOPS)

**Read from**: `annotation/*.xml` (`swathTiming`)

| Field | Why / Used by | Module/Function | Derived/Pass downstream |
|---|---|---|---|
| Per‑burst `azimuthTime` (`sensingStart`, `sensingStop`) | Line window per burst | `deburst.rs` | `burst_line_start/stop` |
| `firstValidSample`/`lastValidSample` per line | Remove invalid edges | `deburst.rs` | valid range windows |
| `linesPerBurst`, `samplesPerBurst` | Tiling & memory | `iw_split_optimized.rs` | burst grid |
| `burstIndex` → `IW1..IW3` ordering | Merge order | `topsar_merge.rs` | burst stacking |
| `subswathGap` (if present) | Seam handling | `topsar_merge.rs` | stitch offsets |

**Derived**: `burst_to_global_line` mapping; overlap windows for power preservation check (→ `validation_gates.rs`).

---

## 4) Doppler (Centroid & Rate)

**Read from**: `annotation/*.xml` (`dcEstimate`, polynomials over range & azimuth)

| Field | Why / Used by | Module/Function | Derived/Pass downstream |
|---|---|---|---|
| Doppler centroid polynomial(s) | Zero‑Doppler time solve; deskew | `terrain_correction.rs` forward model; deburst (optional deskew) | `f_dc(line,pix)` |
| Doppler rate / FM rate | Focusing consistency; optional check | TC, validation | `k_az` |
| Reference epoch for Doppler model | time base | TC | use **same** time base as orbit |

**Note**: In your logs the **derivative** was near‑constant ~1.05e3 Hz/s; be consistent with time bases to keep bounds sane.

---

## 5) Orbit State Vectors (OSV)

**Read from**: `annotation/orbit.xml` or external Precise/Restituted orbits

| Field | Why / Used by | Module/Function | Derived/Pass downstream |
|---|---|---|---|
| `t_i, pos_i, vel_i` time‑tagged | Interpolate sat position/velocity | `orbit.rs` | `sat_pos(t), sat_vel(t)` |
| `ref_epoch` (your `orbit_ref_time`) | Time base for solvers | `terrain_correction.rs` | requires consistent conversion |
| Interp method (Lagrange/Cubic Hermite) | Smooth derivatives | `orbit.rs` | robust `sat_pos/vel` |

**Derived**: `wavelength = c / f0`; `heading` estimate from velocity projected to ENU (useful fallback).

---

## 6) Radar & Instrument

**Read from**: `annotation/*.xml` (`radarParameters`)

| Field | Why / Used by | Module/Function | Derived/Pass downstream |
|---|---|---|---|
| `carrierFrequency` `f0` | Doppler→radar geometry; wavelength | `calibrate.rs`, TC | `λ = c / f0` |
| `chirpBandwidth`, `rangeBandwidth` | Resolution, IRW checks | focusing/validation | expected IRW |
| `PRF` | time grids, ambiguity checks | timing, validation | `PRI = 1/PRF` |
| `azimuthFMRate` (if provided) | optional | validation | — |

---

## 7) Radiometric Calibration (β⁰/σ⁰/γ⁰)

**Read from**: `calibration/calibration-*.xml`

| Field | Why / Used by | Module/Function | Derived/Pass downstream |
|---|---|---|---|
| `betaNought` LUT(s) | SLC intensity to β⁰ | `calibrate.rs` | β⁰ map |
| `sigmaNought` LUT(s) | SLC intensity to σ⁰ | `calibrate.rs` | σ⁰ map |
| `gammaNought` LUT(s) | sometimes provided | `calibrate.rs` | γ⁰ map |
| `calibrationConstant` | absolute gain | `calibrate.rs` | linear gain |
| LUT axes (line/pixel domains) | interpolation bounds | `calibrate.rs` + `validation_gates.rs::validate_lut_domain_bounds` | — |

**Derived**: per‑pixel β⁰/σ⁰; radiometric QA stats; carry **which scale** your pipeline standardizes to (recommend: produce σ⁰ linear, and compute γ⁰ only in terrain‑flattening).

---

## 8) Thermal Noise

**Read from**: `noise/noise-*.xml`

| Field | Why / Used by | Module/Function | Derived/Pass downstream |
|---|---|---|---|
| Noise range profiles / 2‑D LUT | Noise removal pre‑calibration | `calibrate.rs` | noise‑subtracted I,Q or power |
| Axes & domains | Bounds checking | calibration + validation | — |

**Derived**: NESZ estimate for validation (→ `validation_gates.rs::validate_nesz_compliance`).

---

## 9) Geolocation/Tie‑Point Grid

**Read from**: `annotation/*.xml` (`geolocationGrid`)

| Field | Why / Used by | Module/Function | Derived/Pass downstream |
|---|---|---|---|
| `lat, lon, h` on coarse grid | Seed for forward model / grid warping | `terrain_correction.rs` | tie‑points |
| `incidenceAngle`, `elevationAngle` (ellipsoid) | Terrain flatten seed | `scientific_terrain_flatten.rs` | use as θᵢ input |
| `demElevation` (if present) | coarse elevation sanity | TC/RTC | — |

**Derived**: dense lat/lon if you interpolate pre‑TC for quick geocoding previews; **but** final should use DEM‑aware TC. Keep tie‑point residuals for QA.

---

## 10) DEM & Geodesy

**Read from**: external DEM; geoid/ellipsoid models

| Field | Why / Used by | Module/Function | Derived/Pass downstream |
|---|---|---|---|
| DEM grid, pixel spacing, CRS | Gradients & normals | `dem.rs`, `scientific_terrain_flatten.rs` | `p=∂z/∂x`, `q=∂z/∂y`, normals |
| Vertical datum (EGM96) awareness | Elevation bias | TC/RTC | offset handling |
| WGS84 constants (a, e²) | LLH↔ECEF | `type_safe_units.rs::angle_ops::latlon_to_ecef` | ECEF coords |

**Derived**: shadow/layover masks from local incidence; keep `cos(θ_local)` for RTC.

---

## 11) Processing Choices & Switches

| Switch / Param | Why / Used by | Module |
|---|---|---|
| Deskew on/off | TOPS burst merging | `topsar_merge.rs` |
| Interp kernels (sinc, cubic) | Radiometry preservation | `deburst.rs`, merge |
| Parallel chunk sizes | Performance | multiple |
| Clamp/strict flags | Prevent silent inversions | `type_safe_units.rs` + `global_clamp_debug` |

---

## 12) Module‑by‑Module Dependency & Outputs

### `slc_reader.rs`
- **Needs**: image dims, data type, polarization, timing basics (`startTime`, `azimuthTimeInterval`), `slantRangeTime` first pixel, `rangeSamplingRate`.
- **Produces**: complex array tiles; per‑tile stats; pass: `dt_az`, `τ0`, `fs`.

### `iw_split_optimized.rs` / `deburst.rs`
- **Needs**: burst list (`sensingStart/Stop`, `linesPerBurst`, `samplesPerBurst`), `first/lastValidSample` per line.
- **Produces**: debursted complex grid, **overlap windows** metadata, mapping burst→global line. Pass overlap to `validation_gates::validate_power_preservation`.

### `topsar_merge.rs`
- **Needs**: per‑burst tiles, overlaps, deskew option, Doppler centroid (optional for deskew).
- **Produces**: seamless azimuth mosaic; pass merged Doppler summary and uncovered‑pixel mask (holes) to `validation_gates`.

### `calibrate.rs`
- **Needs**: calibration LUTs (β⁰/σ⁰), calibration constant, noise LUT; wavelength from radar params.
- **Produces**: β⁰ or σ⁰ in **linear** units + NESZ estimate; pass σ⁰ to RTC; pass QA (mean ocean σ⁰).

### `scientific_terrain_flatten.rs` (RTC)
- **Needs**: σ⁰ (linear), DEM grid & spacing, **azimuth heading ψ**, **ellipsoid incidence θᵢ** (from annotation), optional orbit data.
- **Produces**: γ⁰ (linear), `cos(θ_local)`, `shadow_mask`, quality mask. Pass γ⁰→multilook; keep `cos(θ_local)` for diagnostics.

### `terrain_correction.rs` (Geocoding / forward model)
- **Needs**: orbit (pos/vel vs time, **ref epoch**), Doppler centroid model & time base, timing grids, DEM.  
- **Produces**: map from output grid (LLH) → input (line, pixel), validity flags; derivative diagnostics. Pass geocode transform + residuals to metadata provenance.

### `multilook.rs`
- **Needs**: selected looks (Lr, La), pixel spacing hints.
- **Produces**: averaged power in radiometrically correct manner; pass updated spacing & ENL to provenance.

### `validation_gates.rs`
- **Needs**: overlap power stats, sigma0/γ0 ocean subset, NESZ LUT, geometry checkpoints, LUT domains.
- **Produces**: structured results + summary; fail‑fast on critical domain errors.

---

## 13) Derived Quantities (compute once, reuse)

- **Wavelength**: `λ = c / f0` (store in provenance).
- **Range pixel spacing**: `ΔR = c/(2 fs)`; **slant range** for pixel p: `R(p) = R0 + p·ΔR`.
- **Azimuth time per line**: `dt_az` (from annotation); absolute time: `t_abs = t_product_start + line·dt_az`.
- **Heading**: from orbit velocity projected to ENU near scene center; fallback if annotation heading absent.
- **Ellipsoid incidence θᵢ**: from geolocation grid at scene/burst; average if per‑pixel not needed.
- **Look vector ŝ** (geometric): used in RTC; store once per tile if ψ, θᵢ constant.
- **Local incidence cosine** `cos(θ_local)` and **shadow mask**: keep for QA & mask propagation.
- **Deskew phase** (if applied): store deskew parameters used during merge.

---

## 14) Provenance & “Pass‑Along” Checklist

When a function finishes, pass these forward (serialize into processing metadata):

- **Time bases & epochs**: `orbit_ref_epoch`, `t_product_start`.
- **Grids**: `dt_az`, `τ0`, `fs`, `ΔR`.
- **Per‑burst windows**: start/stop lines, valid sample masks, overlap regions.
- **Calibration choice**: which LUT used (β⁰ or σ⁰), calibration constant value, noise LUT ID.
- **Wavelength**.
- **RTC inputs**: ψ, θᵢ, DEM spacing, DEM CRS, DEM vertical datum.
- **RTC outputs**: `cos(θ_local)` stats, shadow fraction, γ⁰ min/max (linear & dB) for QA.
- **Geocoding**: solver type (secant/bisection), iterations, residuals, failure flags.
- **Validation**: pass/fail per gate, key metrics (PSLR/ISLR/IRW, NESZ diff, CE90).

---

## 15) Common Pitfalls & Guards

- **Time base mismatch**: Always convert solver `azimuth_time` (relative to **orbit epoch**) to product time before computing line indices.  
- **LUT domain overrun**: check with `validate_lut_domain_bounds` before sampling.
- **SLC vs Intensity**: calibration normalization differs; for **SLC** preserve complex; for **intensity** (power) apply LUTs to power, not amplitude.
- **TOPS deskew**: if applied, do it **before** merging; keep radiometric normalization across overlaps.
- **DEM spacing assumption**: if not explicit, infer carefully (see `infer_dem_spacing_if_missing`); wrong spacing corrupts normals and RTC.
- **Angle units**: use `type_safe_units` to avoid deg↔rad mixups; clamp inversions guarded by `SARDINE_STRICT_CLAMP` and `dbg_clamp`.

---

## 16) Minimal “Required Reads” by Stage (IW SLC)

1. **Open product**: dimensions, pols, start time, `dt_az`, `fs`, `τ0`, swath.  
2. **Bursts**: per‑burst times, valid sample windows.  
3. **Orbit**: OSV + ref epoch.  
4. **Doppler**: centroid polynomial(s) + reference time base.  
5. **Calibration**: β⁰/σ⁰ LUT, calibration constant, noise LUT.  
6. **Geolocation grid**: θᵢ and optional heading if available.  
7. **DEM**: raster + spacing + CRS.

---

## 17) Quick Dependency Graph (text)

```
SLC/Annotation ─┬─> Timing (dt_az, τ0, fs) ───────────────────────────┐
                ├─> Bursts (windows) ──> Deburst ──> Merge ───────────┤
                ├─> Doppler (centroid) ──> TC forward model           │
                ├─> Orbit (pos/vel, epoch) ─┐                         │
                │                            └─> TC (LLH→SAR mapping) ├─> Geocoded coordinates
                ├─> Calibration LUTs + Noise ──> Calibrate (σ⁰/β⁰) ───┤
                ├─> GeolocationGrid (θᵢ, heading) ──> RTC (γ⁰) ───────┤
                └─> DEM (spacing, CRS) ───────┬───────────────────────┘
                                               └─> Normals, cosθ, shadow
```

---

## 18) Recommended “carry” struct

Define a `ProcessingContext` that bundles and persists across stages:

```rust
struct ProcessingContext {
  // time bases
  orbit_ref_epoch: f64,
  t_product_start: f64,

  // timing & range
  dt_az: f64, fs: f64, tau0: f64, dr: f64,

  // geometry
  wavelength: f64, heading_rad: f64, theta_i_rad: f64,

  // bursts
  bursts: Vec<BurstMeta>, overlap_windows: Vec<Overlap>,

  // calibration
  cal_kind: CalKind, cal_const: f64, lut_ids: LutIds, nesz_profile: Vec<f32>,

  // DEM
  dem_spacing_m: f32, dem_crs: String, geoid: String,

  // RTC
  cos_theta_stats: (f32,f32,f32), shadow_frac: f32,

  // validation
  gate_results: Vec<ValidationResult>,
}
```

This prevents re‑reading annotations and enforces a single source of truth.
