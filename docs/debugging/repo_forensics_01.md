# Repo Forensics 01

## Status
Initial repository forensics completed.

## Key provisional conclusions
- The pipeline is substantial and non-trivial.
- The canonical path appears to be the 15-stage `BackscatterProcessor`.
- Python likely orchestrates; Rust likely contains most numerical heavy lifting.
- The merge -> terrain correction boundary is currently the highest-risk interface.
- Burst timing domain transitions are likely a major defect source.
- Existing comments labeled as critical fixes suggest unresolved structural fragility rather than isolated bugs.
- Case 001 is currently most consistent with merge / overlap / merge-to-geocoding coupling defects.

## Important note
These are provisional conclusions from static inspection, not verified truth.

What I Inspected
All files under docs (architecture, known_failures, benchmark_cases, rebuild_plan, verification_strategy)
AGENTS.md, LEGACY_STATUS.md, README.md, README.md
pyproject.toml, Cargo.toml
Full Python package sardine/ (all 30+ .py files including processors, validation, metadata)
Full Rust src/ tree (182 files): lib.rs, types.rs, core/mod.rs, io/mod.rs, bindings/mod.rs, constants/mod.rs
Key Rust implementations: terrain_correction, topsar_merge, deburst, calibration, geometry, range_doppler
All scripts (26), all tests (23 Python+Rust), benchmark suite (3 files)
What I Think Is True
The pipeline is a real, substantial SAR processing system — not a toy or prototype.
The 15-stage BackscatterProcessor is the canonical pipeline path.
All heavy computation happens in Rust; Python orchestrates stages and manages state.
The merge→terrain_correction coupling is the single highest-risk interface.
Burst timing domain transitions (orbit-relative vs absolute vs product-start) are a documented failure class.
The "CRITICAL FIX" comments in the code suggest these bugs were hit in practice.
Case 001 symptoms are consistent with merge/overlap or merge→TC dimension coupling failures.
Constants, types, and I/O parsing are likely correct and reusable.
The validation framework (gates, flags, scientific validator) is well-designed infrastructure.
Strict science mode is the only path with adequate error checking.
What I Am Unsure About
Whether the UNION overlap fix in detect.rs fully resolves the bright-stripe / burst-duplication issue, or only partially addresses it.
Whether the "CRITICAL FIX" dimension override in terrain.py:2485 is always reached in every pipeline path.
Whether the burst timing JSON deserialization actually preserves azimuth_time_rel_orbit correctly through the Python→JSON→Rust round-trip.
Whether the seed-grid Newton-Raphson solver in terrain correction converges reliably for all geometries, or has edge cases.
Whether the geoid conversion (EGM96↔WGS84) implementation matches the expected precision.
Whether the calibration LUT bilinear interpolation handles burst boundaries correctly.
Whether the partial_pipeline.py secondary entry point exercises the same coupling paths as the main BackscatterProcessor.
The exact state of the Rust code — whether it compiles against current dependencies.
Smallest Safe Next Step
Reproduce Case 001 diagnostics without running the full pipeline.

Specifically: read the existing output files for Case 001 (the processed VV backscatter output at outputs), and inspect:

Output array shape vs expected shape from annotation metadata
NaN/nodata pattern (spatial structure of missing regions)
Whether duplicated strips have identical or near-identical pixel values
The geotransform and CRS of the existing output
This requires no code changes, only reading existing artifact files with standard tools (rasterio/GDAL). It would confirm whether the failure pattern matches the merge or TC coupling hypothesis, and narrow the next investigation step