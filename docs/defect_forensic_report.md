# SARdine Legacy Pipeline — Defect-Oriented Forensic Report

**Date:** 2025-01-XX  
**Scope:** Broad algorithmic, architectural, and numerical defect identification across all processing stages  
**Method:** Code inspection and pattern analysis (no runtime testing)

---

## 1. Defect Hypothesis Table (Stage-by-Stage)

### Stage 1–2: Read Metadata / Apply Orbit

| ID | Hypothesis | Severity | Evidence |
|----|-----------|----------|----------|
| D01 | **Orbit time reference epoch mismatch** — `orbit_ref_epoch_utc` vs actual orbit state vector reference time may diverge if computed differently in Python vs Rust | HIGH | `types.rs:182` marks orbit ref epoch as "CRITICAL for orbit interpolation"; `types.rs:405` warns `azimuth_time_rel_orbit` defaults to `0.0` which causes "WRONG geocoding" |
| D02 | **RESORB fallback silently degrades accuracy** — `orbit.rs:78-85` falls back to RESORB when POEORB is missing, with only a log message | MEDIUM | Orbit accuracy degrades from ~5cm to ~10cm; no quality flag propagated to output metadata |

### Stage 3: IW Split (Subswath Geometry Extraction)

| ID | Hypothesis | Severity | Evidence |
|----|-----------|----------|----------|
| D03 | **Subswath geometry validation is permissive** — `processor.py:870` logs missing fields as "optional" via `debug()` and returns early without failing | MEDIUM | Missing `slant_range_time`, `burst_count`, or `range_pixel_spacing` would cause downstream silent failures in terrain correction |
| D04 | **Polarization key matching is fragile** — `terrain.py:1459-1470` tries 3 different key formats (exact, `Polarization::` prefix, substring match) to find subswath data | LOW | Could match wrong key if polarization names overlap; the `Polarization::` prefix suggests Rust Debug format leaking into Python |

### Stage 4: TOPSAR Deburst

| ID | Hypothesis | Severity | Evidence |
|----|-----------|----------|----------|
| D05 | **Polynomial time offset domain detection is heuristic** — `iw_deburst.rs:107-140` uses threshold `86_400` to decide if t0 is "relative to day start" vs absolute epoch | HIGH | If a burst has absolute t0 in UTC seconds-of-day (0–86400) but the burst reference time is also UTC seconds-of-day (not UNIX epoch), the heuristic produces wrong result. The threshold was previously 1e6, suggesting this logic has been patched multiple times |
| D06 | **Missing polynomial t0 falls back to 0.0** — `iw_deburst.rs:155-161` uses 0.0 when both DC and FM t0 are missing, with only a warning | HIGH | Zero offset means deramp phase is computed relative to burst center rather than actual polynomial reference time; introduces phase error proportional to burst duration (~2.7s for S1 IW) |
| D07 | **Chunked deburst accumulation changed from power to complex** — `iw_deburst.rs:76-81` has `CRITICAL FIX` comment showing the accumulator was changed from `power_acc: Array2<f32>` to `complex_acc: Array2<SarComplex>` | MEDIUM | If the standard (non-chunked) path still uses a different accumulation order, the two paths produce different results |

### Stage 5: Radiometric Calibration

| ID | Hypothesis | Severity | Evidence |
|----|-----------|----------|----------|
| D08 | **Calibration LUT is azimuth-invariant assumption** — `antenna.rs:210` states "S1 SLC calibration LUTs are azimuth-invariant (sigma0 varies < 0.001 dB across bursts)" | LOW | This is correct for sigma0 but NOT for the noise LUT, which varies per burst. If the same azimuth-invariant assumption is applied to noise subtraction, thermal noise correction will be wrong at burst boundaries |
| D09 | **Radiometric domain transitions are implicit** — No explicit domain tag tracks whether data is in complex SLC, intensity, sigma0, beta0, or gamma0 domain at any given pipeline stage | MEDIUM | The pipeline relies on stage ordering to know the domain. If stages are reordered or skipped, data can be in the wrong radiometric domain without any check |

### Stage 6: Subswath Merge

| ID | Hypothesis | Severity | Evidence |
|----|-----------|----------|----------|
| D10 | **Azimuth overlap uses UNION instead of INTERSECTION** — `detect.rs:251-252` explicitly changed from `max/min` (intersection) to `min/max` (union) to fix "bright stripes" | HIGH | The original intersection logic caused doubled intensity when both swaths contributed without weights. But UNION means the overlap region extends beyond where BOTH swaths have valid data — rows where only ONE swath has data still get a blending weight applied to that one swath, potentially dimming it |
| D11 | **feather_width parameter is silently ignored** — `weights.rs:42-44` states "feather_width parameter is kept for backwards compatibility but ignored as we always use the full overlap width" | LOW | API contract violation; callers passing custom feather widths get full-width cosine regardless |
| D12 | **Hardcoded gap row diagnostics** — `intensity.rs:92,152` hardcode `gap_row_start=12426, gap_row_end=12470` for diagnostic tracking | LOW | Scene-specific debugging code left in production; these exact rows are meaningless for other scenes |
| D13 | **Overlap weight computation assumes range-only blending** — Weights are 1D in range direction only, broadcast to all azimuth rows uniformly | MEDIUM | At burst boundaries within the overlap region, different bursts from adjacent subswaths may not be temporally aligned. Uniform azimuth weighting doesn't account for this |

### Stage 7: Multilook

| ID | Hypothesis | Severity | Evidence |
|----|-----------|----------|----------|
| D14 | **Incidence angle fallback introduces 5–15% error** — `multilook.py:84-102` falls back to mode-based incidence angles (e.g., 29.1°–46.0° for IW) when metadata parsing fails | HIGH | Ground-range spacing = slant_range_spacing / sin(θ). If θ is wrong by ±5°, the multilook factors change, affecting the output resolution and downstream terrain correction coordinate scaling |
| D15 | **Multilook preserve_power bug was fixed but regression risk remains** — `apply/tests.rs:269-293` contains regression tests for a bug where `preserve_power=true` multiplied by fill_factor, causing partial-block dimming | MEDIUM | The fix is tested, but `preserve_power` is disabled by default. If enabled, the interaction with dB conversion needs validation |

### Stage 8–9: Terrain Flatten / Speckle Filter

| ID | Hypothesis | Severity | Evidence |
|----|-----------|----------|----------|
| D16 | **No terrain flattening implementation found** — Stage 8 "Terrain Flatten" is listed in the pipeline but the actual implementation delegates to terrain correction (Stage 10) | HIGH | If terrain flattening is supposed to apply an RTC factor before speckle filtering, its absence means the speckle filter operates on uncorrected σ⁰ instead of γ⁰, which changes the speckle statistics |

### Stage 10: Terrain Correction (Range-Doppler Geocoding)

| ID | Hypothesis | Severity | Evidence |
|----|-----------|----------|----------|
| D17 | **diag_clamp returns NaN on inverted bounds in production** — `mod.rs:30-51` returns `f64::NAN` when min > max in non-strict mode, silently poisoning downstream calculations | CRITICAL | Any pixel where the clamp bounds are inverted (e.g., due to an upstream geometry error) produces NaN that propagates through all arithmetic. The env-var gated strict mode means this is disabled by default |
| D18 | **Negative range pixel clamping to 0.0** — `mod.rs:932-940` clamps negative range pixels to 0.0 within a 50-pixel tolerance | HIGH | Clamping to 0.0 means pixels that should map near-range get assigned to sample 0, which is IW1's first sample. This is geometrically wrong — these should be rejected, not reassigned to the edge |
| D19 | **slant_range_to_native_pixel fallback tolerance is 1ms (~150km)** — `range_doppler.rs:650-656` uses `MAX_FALLBACK_TOLERANCE = 0.001s` | HIGH | As the comment itself acknowledges, "150km tolerance is larger than subswath spacing (~46-58km)". This means a pixel that belongs to IW1 could be falsely matched to IW3, causing a ~26k pixel range shift |
| D20 | **Burst segment miss falls back to closest segment or linear mapping** — `mod.rs:1277-1322` has a cascade of fallbacks when no burst segment matches the azimuth time | HIGH | The "closest segment" fallback clamps to segment boundaries (start_line or end_line), which is spatially incorrect. The "simple linear mapping" fallback `time / azimuth_time_interval` assumes uniform line timing, which is wrong for TOPSAR |
| D21 | **azimuth_time_interval mandatory check only for merged TOPSAR** — `mod.rs:1072-1080` returns None for merged TOPSAR but allows PRF fallback for single-swath | MEDIUM | Even for single-swath processing, using PRF instead of annotation azimuth_time_interval introduces timing errors |
| D22 | **Dimension override in Python uses current working_data shape** — `terrain.py:1415-1426` overrides `total_azimuth_lines` with `working_data.shape[0]`, and computes `native_azimuth_lines = actual * azimuth_ml` | HIGH | If the multilook factor is wrong (see D14), the native dimension computation is wrong, which propagates to `max_valid_range_pixel` in the Rust Range-Doppler solver, causing valid pixels to be rejected |
| D23 | **Geoid conversion is always applied without validation** — `mod.rs:706-710` calls `orthometric_to_ellipsoidal()` on all elevation values | MEDIUM | If the DEM is already in ellipsoidal heights (e.g., some Copernicus DEM products), the conversion is applied twice, shifting all elevations by the geoid undulation (~20-45m in Europe) |
| D24 | **range_doppler_with_geometry clamps elevation to max(0.0)** — `mod.rs:1533` does `elevation_clamped = elevation.max(0.0)` | MEDIUM | This incorrectly zeroes out negative bathymetry for coastal scenes and Dead Sea-type locations. The main `scientific_range_doppler_transformation` does NOT have this clamp (D25), so the RTC path gets different geometry than the geocoding path |
| D25 | **Two Range-Doppler paths with subtle differences** — `scientific_range_doppler_transformation` and `range_doppler_with_geometry` are separate ~200-line functions | HIGH | D24 shows one difference (elevation clamp). Other potential differences: the RTC path (D24) doesn't have the per-subswath burst segment lookup shown in the main path. Any fix applied to one path but not the other creates inconsistency |

### Stage 11: Apply Mask

| ID | Hypothesis | Severity | Evidence |
|----|-----------|----------|----------|
| D26 | **DEM void detection uses limited sentinels** — `terrain.py:194-197` checks `{-32768, -9999, <-9000}` | LOW | Some DEMs use 0 or -99999 or 32767 as nodata. Missing sentinels would pass voids as valid elevation |

### Stage 12: dB Conversion

| ID | Hypothesis | Severity | Evidence |
|----|-----------|----------|----------|
| D27 | **"Emergency bottleneck fix" label on unit conversion** — `unit_conversion.rs:3,82,95,113` has multiple "FIXED: In-place parallel" and "emergency bottleneck fix" comments | MEDIUM | The "broken SIMD versions" comment and "fixes broken export functions that ignored return values" suggest prior implementations silently failed to apply the conversion |

### Stage 13–15: Export / Metadata / QA

| ID | Hypothesis | Severity | Evidence |
|----|-----------|----------|----------|
| D28 | **North-up flip occurs after terrain correction** — `terrain.py:141-168` flips data and adjusts geotransform if pixel_height > 0 | MEDIUM | The flip is applied to the geocoded output, but the geotransform adjustment `transform[3] + transform[5] * (height - 1)` should use `height` not `height - 1` for correct pixel registration |

---

## 2. Top 10 Most Suspicious Code Regions

| Rank | Location | Issue | Risk |
|------|----------|-------|------|
| 1 | `mod.rs:30-51` (diag_clamp) | Returns NaN on inverted bounds in production | Silently poisons entire output tiles |
| 2 | `range_doppler.rs:650` (MAX_FALLBACK_TOLERANCE) | 1ms tolerance larger than subswath spacing | Cross-subswath pixel misassignment |
| 3 | `iw_deburst.rs:107-140` (polynomial_time_offset) | Heuristic domain detection for DC/FM t0 | Wrong deramp phase → phase noise in overlap weighting |
| 4 | `mod.rs:1415-1426` + `multilook.py:84-102` | Dimension override chain depends on incidence angle fallback | Cascading error: wrong θ → wrong ML factor → wrong native dims → wrong TC bounds |
| 5 | `mod.rs:932-940` (negative range pixel clamping) | Clamps to 0.0 instead of rejecting | Geometrically wrong pixel placement at near-range edge |
| 6 | `detect.rs:251-252` (azimuth UNION overlap) | UNION extent may extend beyond valid data for one swath | Potential dimming at subswath azimuth boundaries |
| 7 | `mod.rs:1277-1322` (burst segment miss fallbacks) | Cascade of increasingly wrong fallbacks for azimuth mapping | Wrong azimuth pixel → wrong geocoding position |
| 8 | `iw_deburst.rs:155-161` (missing polynomial t0) | Falls back to 0.0 with only a warning | Phase error across burst |
| 9 | `mod.rs:1533` vs `mod.rs:706-710` | Two Range-Doppler paths with different elevation handling | Inconsistent geometry between geocoding and RTC |
| 10 | `types.rs:405-408` (azimuth_time_rel_orbit default 0.0) | Deserialization defaults to 0.0 for missing timing | Silent geocoding failure: all bursts map to orbit epoch |

---

## 3. False Confidence Zones

These are areas where the code *appears* correct but masks underlying fragility:

### 3A. Validation Theatre
The codebase has extensive validation infrastructure (`ValidationGateway`, `ScientificValidator`, `DesignFlawDetector`) that checks *parameter ranges* but not *algorithmic correctness*. For example:
- Coverage threshold validation ensures 45% fill → but doesn't check if the filled pixels are in the right geographic location
- Wavelength hardcoded-value detection → catches 0.055465... but can't detect if the wavelength is being applied in the wrong equation
- Metadata key normalization → ensures keys exist but not that the values are semantically correct (e.g., `azimuth_time_rel_orbit` present but zero)

### 3B. Instrumentation Masking
The code has ~50 `AtomicBool` "log once" guards (`LOGGED: AtomicBool::new(false)`) that suppress repeated diagnostics. This means:
- Only the FIRST pixel's computation path is logged
- If the first pixel succeeds but subsequent pixels fail differently, the failure mode is invisible
- The first-pixel diagnostic can be misleading if processing order varies across Rayon threads

### 3C. Tolerance Stacking
Multiple layers of tolerance are applied sequentially:
1. Subswath matching: 5% of range_pixel_spacing_time (~0.15ms)
2. Fallback matching: 1ms (MAX_FALLBACK_TOLERANCE)
3. Negative pixel tolerance: -50 to 0 clamped to 0
4. Time window margin: ±10 seconds
5. Azimuth bounds: derived from duration/interval + 50 lines

Each tolerance is individually reasonable, but stacked together they can mask large errors. A pixel that is 0.15ms off in subswath matching, -40 pixels in range (clamped to 0), and 8 seconds outside product window could still be accepted.

### 3D. Test Surface Area
Tests visible in the repo test:
- Weight monotonicity and sum-to-one properties
- Multilook partial-block power preservation
- Burst timing record serialization
- Hardcoded value detection

Tests NOT visible for:
- End-to-end Range-Doppler coordinate accuracy
- Merged subswath temporal alignment
- Geoid conversion correctness
- Cross-subswath pixel boundary continuity
- Burst segment lookup coverage completeness
- Python→Rust metadata round-trip fidelity

---

## 4. Missing Invariants and Diagnostics

| Category | Missing Invariant | Impact |
|----------|-------------------|--------|
| **Coordinate domain** | No type-level tracking of radiometric domain (SLC→intensity→σ⁰→γ⁰→dB) | Cannot statically prevent dB-domain data from being squared or intensity from being log'd twice |
| **Timing domain** | No assertion that all time values in a single processing pass use the same epoch | Mix of UNiX epoch, orbit-relative, product-relative, and UTC seconds-of-day times |
| **Pixel coordinate system** | No tagged type distinguishing native vs multilooked vs geocoded pixels | Wrong scaling factor applied silently |
| **Subswath coverage** | No check that burst segments fully tile the product duration without gaps | Pixels falling in burst gaps get fallback mapping |
| **Range pixel continuity** | No check that overlapping subswath range pixels are contiguous in global coordinates | Gap or double-counted pixels at subswath boundaries |
| **Elevation datum** | No tag or flag indicating whether DEM is orthometric or ellipsoidal | Geoid correction applied unconditionally (D23) |
| **Merge completeness** | No row-level assertion that every output pixel has weight_sum > 0 | Zero-weight pixels produce 0/0 = NaN silently |
| **Export sanity** | No post-geocoding check that output footprint overlaps input bounding box | Complete spatial mismatch would produce all-NaN output without error |

---

## 5. Best Next Diagnostic Steps

In priority order:

1. **Controlled end-to-end coordinate test** — Process a synthetic scene with known ground control points. Compare geocoded positions against truth. This validates the entire Range-Doppler chain including geoid, burst segments, multilook scaling, and epoch conventions. **Estimated value: disproves or confirms D17-D25 in one experiment.**

2. **Burst segment coverage audit** — For a real product, enumerate all `burst_segments` and verify they tile `[product_start_rel_s, product_start_rel_s + product_duration]` without gaps. Print any azimuth time ranges that fall in gaps. **Estimated value: directly tests D20.**

3. **Cross-subswath seam diagnostic** — At the IW1-IW2 and IW2-IW3 overlap boundaries, extract a 100-pixel-wide strip from each subswath BEFORE merge. Cross-correlate in range to verify alignment. Plot intensity profiles across the seam AFTER merge. **Estimated value: tests D10, D13.**

4. **Metadata round-trip test** — From Python, construct a known metadata dict, pass it through `normalize_metadata_for_rust()`, serialize to JSON, deserialize in Rust, and verify all fields match. Check `azimuth_time_rel_orbit`, subswath `slant_range_time`, and `azimuth_time_interval` specifically. **Estimated value: tests D01, D03, D04.**

5. **Dual Range-Doppler path comparison** — For a subset of pixels, call both `scientific_range_doppler_transformation` and `range_doppler_with_geometry` and diff the returned (range, azimuth) coordinates. Any difference reveals D25 inconsistency. **Estimated value: tests D24, D25.**

---

## 6. Remaining Uncertainties

- **Deburst phase correctness** — Cannot assess from code alone whether the deramp/reramp phase is physically correct. Requires test data with known phase response.
- **Calibration LUT range interpolation** — The 1D interpolation in `lut.rs` appears correct for per-row application, but verification against SNAP's calibration requires comparison testing.
- **RTC (Small 2011) implementation** — The RTC module imports in `rtc.rs` exist but the actual area-projection formula was not fully inspected. The `reference_incidence_angle_deg` defaulting to 35° (typical IW mid-swath) is coarse but unlikely to cause large errors.
- **Speckle filter internals** — Not inspected in this review.
- **GDAL write ordering** — Whether the final GeoTIFF export preserves band/pixel order correctly was not inspected.
