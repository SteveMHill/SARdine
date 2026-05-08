//! SLC Reader module with correct time domain handling
//!
//! This module provides a scientifically correct SLC reader that:
//! - Uses orbit-reference-relative seconds for all timing computations
//! - Correctly handles TOPS burst range geometry (no walking range windows)
//! - Uses exclusive upper bounds for valid sample ranges
//! - Properly aligns Doppler centroid evaluation with orbit timing

pub mod burst_records;
pub mod time_utils;

pub use burst_records::{build_burst_records_for_subswath_fixed, BurstTimingDomain};
pub use time_utils::{parse_iso8601_to_datetime_utc, seconds_since_epoch};

#[cfg(test)]
mod tests;
