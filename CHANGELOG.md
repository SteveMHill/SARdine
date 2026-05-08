# Changelog

All notable changes to SARdine are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html)
once it reaches 1.0.0.

## [Unreleased]

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
