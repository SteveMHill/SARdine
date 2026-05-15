# Changelog

All notable changes to SARdine are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html)
once it reaches 1.0.0.

## [Unreleased]

### Fixed

#### Audit round 2 — robustness hardening (2026-05-15)

- **`terrain_correction`**: Added `is_finite()` early-return guard to
  `bilinear_sample_slice`, `bicubic_sample_slice`, and `lanczos_sample_slice`.
  Previously, NaN input coordinates silently produced wrong samples because
  `NaN as isize` evaluates to 0 in Rust.
- **`geoid`**: Added `is_finite()` guard to `Egm96Grid::undulation_m` and
  `Egm2008Grid::undulation_m`. NaN lat/lon previously produced row_f = NaN
  which clamped to row 0 (North Pole cell) and returned a wrong geoid height.
- **`ground_range`**: Three fixes:
  - `slant_col_f` bounds check now includes `!is_finite()` so NaN values are
    correctly rejected (NaN comparisons always return false, bypassing the
    existing `< 0.0` guard).
  - `solve_range_doppler_h0` now returns `None` instead of `Some(NaN, NaN)`
    when the Newton solver produces non-finite coordinates.
  - `RangeInterpolator::interpolate` guards against duplicate slant-range times
    (`t0 == t1`) that previously caused a division-by-zero producing `Inf`.
- **`parse`**: `sort_by(|a, b| a.partial_cmp(b).unwrap())` in burst-delta
  sorting replaced with `.unwrap_or(Equal)` to avoid a panic when a
  malformed annotation contains a NaN burst time.
- **`scene_prep`**: Multilook output buffer allocation uses `checked_mul` for
  `out_lines × out_samples`, returning a descriptive error instead of silently
  overflowing on pathological metadata.
- **`export`**: Two fixes:
  - `write_geotiff_raw_inner` and `write_bigtiff_raw_inner` now clean up the
    `.tmp` file when `write_all` or `flush` fails, not only on rename failure.
  - COG overview pyramid loop termination changed from `||` to `&&`: overviews
    are now generated until *both* dimensions are ≤ tile size, not when either
    is (which stopped too early for non-square images).
- **`dem`**: Longitude normalisation in `elevation_at` and `elevation_at_linear`
  now uses `rem_euclid` to correctly handle values below −180° (e.g. dateline
  wrap-around from the east), replacing the one-sided `>= 180` check.
- **`slice_assembly`**: `burst_index × lines_per_burst` uses `checked_mul` and
  propagates a `SliceAssemblyError::ArithmeticOverflow` instead of overflowing
  silently on corrupted metadata. `merge_bursts` return type updated to
  `Result` accordingly.

#### Audit round 1 — correctness & validation (2026-05-06)

- **`apply_calibration`**: Added warning when K² is suspiciously small
  (`< 1e-6`) to surface miscalibrated products; fixed a truncating cast for
  `col_offset` that could produce wrong LUT indices on large swaths.
- **`calibration`**: Added pixel-index and value validation for range/azimuth
  noise and calibration LUT entries at parse time.
- **`deburst`**: Guard `ati_s <= 0` in burst midpoint selection returns an
  explicit error instead of producing division-by-zero geometry.
- **`output_crs`**: Longitude normalised with `rem_euclid` to handle scenes
  crossing the antimeridian.
- **`orbit`**: Skip (with a warning) orbit state vectors that carry a zero
  position vector rather than propagating NaN into downstream interpolation.
- **`stac`**: STAC item `href` fields written as relative paths so sidecar
  files remain portable when the output directory is moved.
- **`validate`**: Added PRF range check, range pixel spacing sanity check, and
  burst index monotonicity validation.
- **`insar/interferogram`**: `flat_earth_fallback_count` counter added;
  logged at the end of interferogram formation so silent fallbacks are visible.
- **`scene_prep`**: `first_line.checked_sub` avoids underflow on the first
  burst; polarisation channel presence validated before use.
- **`merge_subswaths`**: `NegativeOutputColumns` error returned instead of
  panicking when subswath offsets produce a negative output width.
- **`ground_range`**: `EmptyRangeProfile` error for an empty geolocation grid;
  incidence angle validated to be in (0°, 90°) for every grid point.
- **`parse`**: `first().ok_or(…)` replaces `first().unwrap()` on potentially
  empty burst lists.

### Added
- Apache-2.0 `LICENSE` and `NOTICE` files.
- `CHANGELOG.md` and `CONTRIBUTING.md`.
- GitHub Actions CI workflow running `cargo fmt`, `cargo clippy`, `cargo test`,
  and the `check_no_silent_fallbacks.sh` guard on every push and pull request.
- Workspace-level package metadata (license, repository, edition) inherited
  by member crates.
- `sardine-validation` workspace crate exposing the `sardine-validate`
  binary: a TOML-driven end-to-end radiometric regression harness that
  drives `sardine process` plus `scripts/multilook_compare.py` and
  asserts hard thresholds on median bias, multi-look std deviation and
  joint-valid sample count. Ships with a `munich_s1b_utm32n` scene
  pinned to the post-fix baseline (median bias `+0.0283 dB`, std
  `2.9572 dB`, n = 4 098 547). Six unit tests cover output parsing
  and threshold logic.
- `terrain_correction::cardinal_neighbour_step_deg` helper that converts a
  CRS-units pixel spacing into the equivalent degree step for cardinal-
  neighbour DEM lookups during terrain flattening / LIA computation.
  Three new unit tests pin the WGS84 identity and the UTM
  metres → degrees conversion.

### Changed
- N/A

### Fixed
- **Terrain flattening in metric (UTM) output CRSs.** The cardinal-
  neighbour DEM lookups and `compute_terrain_geometry` previously
  consumed `pixel_spacing_deg` directly, which carries CRS units. In
  UTM mode at 10 m spacing this asked the DEM for elevations 10°
  away from each output pixel, returning out-of-bounds for every
  neighbour and silently masking every pixel as `flat_masked`. The
  S1B Munich UTM 32 N 10 m baseline scene now reports
  `valid = 426 625 590, flat_masked = 2017` (was `valid = 0,
  flat_masked = 426 627 607`). The WGS84 path is bit-identical to
  the previous behaviour.

## [0.1.0] - 2026-04-26

### Added
- Initial pre-release of the `sardine-scene` Rust crate and `sardine-py`
  Python extension. Sentinel-1 IW SLC → terrain-corrected σ⁰/γ⁰ GeoTIFF
  pipeline. Validated to ±0.016 dB median linear bias vs ASF RTC10
  (S1B Munich, 2019-01-23, after 10×10 multilook).
- Phase 0–4 performance optimisations: per-stage timing, DEM HashMap
  spatial index, anchor-interpolated row projection, chunked row
  scheduling, and orbit walk-from-hint cache (terrain correction
  ~2.2× faster, end-to-end ~1.7× faster).

[Unreleased]: https://example.invalid/compare/v0.1.0...HEAD
[0.1.0]: https://example.invalid/releases/tag/v0.1.0
