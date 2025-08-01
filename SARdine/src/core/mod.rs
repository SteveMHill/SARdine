//! Core SAR processing modules

pub mod deburst;
pub mod calibrate;
pub mod multilook;
pub mod terrain_flatten;
pub mod speckle_filter;
pub mod terrain_correction;
pub mod topsar_merge;
pub mod iw_merge;  // New IW merge module for TOPSAR processing
pub mod advanced_masking;
pub mod quality_assessment;  // Comprehensive quality assessment module
pub mod metadata_provenance;  // Metadata and provenance tracking module

// Re-export main types
pub use deburst::{DeburstProcessor, BurstInfo, TopSarDeburstProcessor, DeburstConfig};
pub use calibrate::{CalibrationProcessor, CalibrationCoefficients, CalibrationType};
pub use multilook::{MultilookProcessor, MultilookParams};
pub use terrain_flatten::{TerrainFlattener, TerrainFlatteningParams};
pub use speckle_filter::{SpeckleFilter, SpeckleFilterParams, SpeckleFilterType};
pub use terrain_correction::{TerrainCorrector, RangeDopplerParams, GroundPoint, GeocodedPixel};
pub use topsar_merge::{TopsarMerge, MergedSwathData, OverlapRegion, OutputGrid, merge_iw_subswaths};
pub use iw_merge::{IwMergeProcessor, SubSwathInfo, IwMergeConfig};  // New IW merge exports
pub use advanced_masking::{
    AdvancedMaskingProcessor, AdvancedMaskingConfig, AdvancedMaskResult, 
    MaskingMethod, apply_advanced_masking, apply_water_masking, apply_terrain_masking
};
pub use quality_assessment::{QualityAssessor, QualityConfig, QualityAssessment, QualityStatistics};
pub use metadata_provenance::{MetadataManager, ProcessingMetadata, InputProvenance, QualityMetrics, UncertaintyEstimates};
