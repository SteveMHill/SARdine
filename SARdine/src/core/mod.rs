//! Core SAR processing modules

pub mod deburst;
pub mod calibrate;
pub mod multilook;
#[cfg(test)]
pub mod multilook_validation_tests;
pub mod terrain_flatten;  // Legacy complex orbital approach
pub mod scientific_terrain_flatten;  // Industry-standard geometric approach
pub mod speckle_filter;
#[cfg(test)]
pub mod speckle_filter_validation_tests;
pub mod terrain_correction;
pub mod topsar_merge;
pub mod iw_merge;  // New IW merge module for TOPSAR processing
pub mod iw_split_optimized;  // Optimized IW split module with parallel processing
pub mod advanced_masking;
pub mod quality_assessment;  // Comprehensive quality assessment module
pub mod metadata_provenance;  // Metadata and provenance tracking module
pub mod simd_optimizations;  // SIMD-accelerated mathematical operations
pub mod memory_optimizations;  // Memory allocation and cache optimization utilities

// Re-export main types
pub use deburst::{DeburstProcessor, BurstInfo, TopSarDeburstProcessor, DeburstConfig};
pub use calibrate::{CalibrationProcessor, CalibrationCoefficients, CalibrationType};
pub use multilook::{MultilookProcessor, MultilookParams};
pub use terrain_flatten::{TerrainFlattener as LegacyTerrainFlattener, TerrainFlatteningParams as LegacyTerrainFlatteningParams};
pub use scientific_terrain_flatten::{
    TerrainFlattener, TerrainFlatteningParams, ProcessingMode,
    terrain_flatten, terrain_flatten_orbital
};
pub use speckle_filter::{SpeckleFilter, SpeckleFilterParams, SpeckleFilterType};
pub use terrain_correction::{TerrainCorrector, RangeDopplerParams, GroundPoint, GeocodedPixel};
pub use topsar_merge::{TopsarMerge, MergedSwathData, OverlapRegion, OutputGrid, merge_iw_subswaths};
pub use iw_merge::{IwMergeProcessor, SubSwathInfo, IwMergeConfig};  // New IW merge exports
pub use iw_split_optimized::{IwSplitOptimized, SplitResult, SubswathGeometry};  // Optimized IW split exports
pub use advanced_masking::{
    AdvancedMaskingProcessor, AdvancedMaskingConfig, AdvancedMaskResult, 
    MaskingMethod, apply_advanced_masking, apply_water_masking, apply_terrain_masking
};
pub use quality_assessment::{QualityAssessor, QualityConfig, QualityAssessment, QualityStatistics};
pub use metadata_provenance::{MetadataManager, ProcessingMetadata, InputProvenance, QualityMetrics, UncertaintyEstimates};
