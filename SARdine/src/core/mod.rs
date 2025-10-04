//! Core SAR processing modules

pub mod advanced_masking;
pub mod calibrate;
pub mod context_extraction; // Extract ProcessingContext from annotation data
pub mod coordinate_frames; // Type-safe coordinate frame handling
pub mod dc_fm_provider; // DC/FM-rate provider trait with mocking support
pub mod deburst;
// pub mod deburst_algorithmic_fixes; // Algorithmic fixes for TOPS deburst - MODULE NOT FOUND
pub mod deburst_optimized; // Ultra-high performance TOPS deburst with targeted optimizations
pub mod f64_numerics; // High-precision f64 numerics for geometry/Doppler/NR
pub mod fast_unit_conversion;
pub mod global_clamp_debug; // Global diagnostic clamp instrumentation
pub mod mask_propagation; // Scientifically-valid mask propagation through processing chain
pub mod memory_optimizations; // Memory allocation and cache optimization utilities
pub mod metadata_parser; // Single-pass metadata parser with SoA arrays
pub mod metadata_provenance; // Metadata and provenance tracking module
pub mod metadata_strictness; // Strict metadata validation with no fallbacks
pub mod multilook;
pub mod processing_context; // Unified processing context and structured provenance
// Note: multilook_validation_tests temporarily removed during refactoring
pub mod power_preserving_resample; // Power-preserving azimuth resampling (deramp-FFT-reramp)
pub mod precision_standards; // Numeric precision and reproducibility standards
pub mod quality_assessment; // Comprehensive quality assessment module
pub mod robust_doppler_solver; // Robust zero-Doppler solver (bracket + secant method)
pub mod scientific_terrain_flatten; // Industry-standard geometric approach
pub mod simd_optimizations; // SIMD-accelerated mathematical operations
pub mod speckle_filter;
#[cfg(test)]
pub mod speckle_filter_validation_tests;
pub mod terrain_correction;
pub mod terrain_flatten; // Legacy complex orbital approach
pub mod topsar_merge;
pub mod topsar_merge_optimized; // High-performance TOPSAR merge with bottleneck fixes
pub mod type_safe_units; // Type-safe units to prevent unit confusion
pub mod validation_gates; // Comprehensive validation framework for processing quality gates
pub mod validated_processing_pipeline; // Complete processing pipeline with integrated validation

// Re-export main types
pub use advanced_masking::{
    apply_advanced_masking, apply_terrain_masking, apply_water_masking, AdvancedMaskResult,
    AdvancedMaskingConfig, AdvancedMaskingProcessor, MaskingMethod,
};
pub use calibrate::{CalibrationCoefficients, CalibrationProcessor, CalibrationType};
pub use context_extraction::{
    extract_from_annotation, update_with_calibration, update_with_dem,
};
pub use coordinate_frames::{
    BurstIdx, SubswathIdx, StitchedIdx, BurstToSubswathConverter, SubswathToStitchedConverter,
    CoordinatePosition, BurstFrame, SubswathFrame, StitchedFrame,
};
pub use dc_fm_provider::{
    DcFmRateProvider, PolynomialDcFmProvider, MockDcFmProvider, CachedDcFmProvider,
    DcFmProviderFactory, DcPolynomial, FmPolynomial,
};
pub use deburst::{BurstInfo, DeburstConfig, DeburstProcessor, TopSarDeburstProcessor};
// pub use deburst_algorithmic_fixes::{
//     AlgorithmicDeburstProcessor, BurstOverlap, ComplexSeamWeights, FractionalFIR,
// };
pub use deburst_optimized::{
    deburst_optimized, deburst_optimized_enhanced, OptimizedDeburstResult,
    OptimizedTopsDeburstProcessor, PerformanceMetrics,
};
pub use f64_numerics::{
    Point3D, Point2D, GeographicCoord, RadarCoord, ImageCoord, GeometryF64, DopplerF64,
    NewtonRaphsonF64, RangeDopplerToGeo, PrecisionConverter,
};
// Use canonical unit conversion from constants module
pub use crate::constants::unit_conversion::{
    db_to_linear_inplace, linear_to_db_inplace, 
    db_to_linear_parallel_chunked, linear_to_db_parallel_chunked,
    db_to_linear_inplace_par, linear_to_db_inplace_par,
    export_db_parallel,
};
// Keep IncidenceTerms and Units from fast_unit_conversion (not duplicated)
pub use fast_unit_conversion::{IncidenceTerms, Units};
pub use metadata_parser::{
    SinglePassXmlParser, CompactCalibrationData, CompactNoiseData, CompactGeometricData,
    CompactTimingData, CalibrationUnits, UnitType,
};
pub use metadata_provenance::{
    InputProvenance, MetadataManager, ProcessingMetadata, QualityMetrics, UncertaintyEstimates,
};
pub use metadata_strictness::{
    BetaVariabilityResult, CoverageResult, LineIndexingResult, PowerSeamResult, 
    StrictMetadataValidator, StrictValidationResult, UnitsValidationResult,
};
pub use multilook::{MultilookParams, MultilookProcessor};
pub use power_preserving_resample::{
    PowerPreservingResampler, calculate_power_ratio, power_preserving_resample_2d,
};
pub use processing_context::{
    ProcessingContext, ProcessingContextBuilder, ProductInfo, TimingContext,
    RadarParameters, GeometryContext, BurstContext, CalibrationContext,
    DemContext, ProcessingStage, ProvenanceValidationResult, StageStatus,
};
pub use precision_standards::{
    DeterministicRng, PrecisionStandards, initialize_precision_standards,
    geometry_precision, deterministic_parallel, validation,
};
pub use quality_assessment::{
    QualityAssessment, QualityAssessor, QualityConfig, QualityStatistics,
};
pub use scientific_terrain_flatten::{
    terrain_flatten, terrain_flatten_orbital, ProcessingMode, TerrainFlattener,
    TerrainFlatteningParams,
};
pub use speckle_filter::{SpeckleFilter, SpeckleFilterParams, SpeckleFilterType};
pub use terrain_correction::{
    DopplerCentroidModel, GeocodedPixel, GroundPoint, RangeDopplerParams, TerrainCorrector,
};
pub use terrain_flatten::{
    TerrainFlattener as LegacyTerrainFlattener,
    TerrainFlatteningParams as LegacyTerrainFlatteningParams,
};
pub use topsar_merge::{
    merge_iw_subswaths, MergedSwathData, OutputGrid, OverlapRegion, TopsarMerge,
};
pub use type_safe_units::{
    Radians, Degrees, Meters, Seconds, Hertz, MetersPerSecond, angle_ops, distance_ops, time_ops,
};
pub use validation_gates::{
    ValidationGates, ValidationResult, ValidationSummary,
};
pub use validated_processing_pipeline::{
    ValidatedProcessingPipeline, ProcessingParameters, ProcessedSarData, 
    DeburstMetrics,
};
