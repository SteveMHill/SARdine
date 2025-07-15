//! Core SAR processing modules

pub mod deburst;
pub mod calibrate;

// Re-export main types
pub use deburst::{DeburstProcessor, BurstInfo};
pub use calibrate::{CalibrationProcessor, CalibrationCoefficients, CalibrationType};
