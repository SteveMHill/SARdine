//! InSAR processing pipeline for Sentinel-1 IW SLC data.
//!
//! # Processing overview
//!
//! The InSAR pipeline produces coherence and optionally wrapped interferometric
//! phase from two co-temporal Sentinel-1 IW SLC scenes.  All computations are
//! performed in the reference scene's slant-range geometry; the final product
//! is geocoded via the existing range-Doppler terrain correction.
//!
//! ## Pipeline stages
//!
//! 1. **Parse** — read both SAFE directories; verify DC estimates and FM rates
//!    are present.
//! 2. **Deramp** — apply TOPS azimuth deramping to remove steering phase ramp
//!    (module [`deramp`]).
//! 3. **Co-registration** — geometric sparse-grid co-registration using the
//!    reference and secondary orbits (module [`coreg`]).
//! 4. **Interferogram + coherence** — flat-earth phase removal, complex
//!    cross-product, windowed coherence estimation (module [`interferogram`]).
//! 5. **Geocoding** — pass the coherence raster through the existing terrain
//!    correction pipeline.
//!
//! # Domain notes
//!
//! - All azimuth times are absolute UTC to prevent epoch-coupling bugs.
//! - Slant-range times are two-way (seconds).
//! - DC polynomial `t0` is per-estimate (not global); always evaluate via
//!   [`DcEstimate::evaluate`] rather than raw coefficient indexing.
//! - Deramping does NOT change the power of any pixel (it is a unitary rotation).
//! - Coherence is the amplitude of the normalised complex coherence `γ`; it is
//!   in [0, 1] and is phase-invariant.

pub mod coreg;
pub mod deramp;
pub mod error;
pub mod interferogram;
