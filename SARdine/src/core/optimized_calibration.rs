//! Deprecated optimized calibration module.
//! The high-performance calibration path has been merged into `core::calibrate`.
//!
//! This stub re-exports the canonical calibration APIs to preserve backwards compatibility
//! while avoiding duplicate implementations.

pub use crate::core::calibrate::*;