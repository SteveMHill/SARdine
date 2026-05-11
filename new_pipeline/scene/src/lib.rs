//! Sentinel-1 scene metadata model for the SARdine rebuild.
//!
//! This crate defines the validated metadata representation for a Sentinel-1 SLC scene.
//! All scientifically required fields are non-optional and validated on construction.
//!
//! # Design Principles
//!
//! - **No optional fields** for scientifically required values.
//! - **Explicit time domains**: burst times stored as absolute UTC; convert to
//!   [`OrbitRelativeSeconds`] via [`BurstEntry::azimuth_time_rel`] at point of use.
//! - **Units in field names**: `_m`, `_s`, `_hz`, `_deg` suffixes prevent unit confusion.
//! - **Exclusive-end bounds**: all `[first, last)` ranges match Rust range semantics.
//! - **Validation collects all errors**: [`validate::check`] reports every violation, not just the first.

pub mod apply_calibration;
pub mod calibration;
pub mod cog;
pub mod deburst;
pub mod dem;
pub mod dem_fetch;
pub mod export;
pub mod geodesy;
pub mod geoid;
pub mod geoid_fetch;
pub mod ground_range;
pub mod insar;
pub mod merge_subswaths;
pub mod multi_slc_reader;
pub mod orbit;
pub mod multi_pol;
pub mod orbit_fetch;
pub mod output_crs;
pub mod parse;
pub mod pipeline_options;
pub mod provenance;
pub mod run;
pub mod run_provenance;
pub mod scene_prep;
pub mod stac;
pub mod slc_fetch;
pub mod slc_reader;
pub mod slice_assembly;
pub mod speckle;
pub mod terrain_correction;
pub mod types;
pub mod validate;

pub use types::*;
pub use validate::{ValidationError, ValidationErrors};
