//! Core SAR processing modules

pub mod deburst;
pub mod calibrate;
pub mod multilook;
pub mod terrain_flatten;
pub mod speckle_filter;

// Re-export main types
pub use deburst::{DeburstProcessor, BurstInfo};
pub use calibrate::{CalibrationProcessor, CalibrationCoefficients, CalibrationType};
pub use multilook::{MultilookProcessor, MultilookParams};
pub use terrain_flatten::{TerrainFlattener, TerrainFlatteningParams};
pub use speckle_filter::{SpeckleFilter, SpeckleFilterParams, SpeckleFilterType};
