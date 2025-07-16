//! Core SAR processing modules

pub mod deburst;
pub mod calibrate;
pub mod multilook;
pub mod terrain_flatten;
pub mod speckle_filter;
pub mod terrain_correction;
pub mod topsar_merge;

// Re-export main types
pub use deburst::{DeburstProcessor, BurstInfo};
pub use calibrate::{CalibrationProcessor, CalibrationCoefficients, CalibrationType};
pub use multilook::{MultilookProcessor, MultilookParams};
pub use terrain_flatten::{TerrainFlattener, TerrainFlatteningParams};
pub use speckle_filter::{SpeckleFilter, SpeckleFilterParams, SpeckleFilterType};
pub use terrain_correction::{TerrainCorrector, RangeDopplerParams, GroundPoint, GeocodedPixel};
pub use topsar_merge::{TopsarMerge, MergedSwathData, OverlapRegion, OutputGrid, merge_iw_subswaths};
