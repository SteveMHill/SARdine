//! Core SAR processing modules

pub mod calibration;
pub mod context; // Grouped context utilities (extraction + processing context)
pub mod deburst;
pub mod geometry; // Type-safe coordinate frame handling, geodesy helpers, DC/FM provider
pub mod global_clamp_debug; // Global diagnostic clamp instrumentation
pub mod metadata; // Grouped metadata parsing/validation/provenance helpers
pub mod multilook;
pub mod perf;
// Validation tests now live under multilook/tests
pub mod masking;
pub mod quality; // Grouped quality helpers
pub mod quality_flags; // Scientific audit quality flags (NOISE-1, DEBURST-1, RTC-1)
pub mod speckle_filter;
pub mod terrain_correction;
pub mod topsar_merge;

// Re-export main types
pub use calibration::{CalibrationCoefficients, CalibrationType};
pub use context::context_extraction::{
    extract_from_annotation, update_with_calibration, update_with_dem,
};
pub use deburst::{BurstInfo, DeburstConfig, DeburstProcessor, TopSarDeburstProcessor};
pub use geometry::coordinate_frames::{
    BurstFrame, BurstIdx, BurstToSubswathConverter, CoordinatePosition, StitchedFrame, StitchedIdx,
    SubswathFrame, SubswathIdx, SubswathToStitchedConverter,
};
pub use geometry::dc_fm_provider::{
    CachedDcFmProvider, DcFmProviderFactory, DcFmRateProvider, DcPolynomial, FmPolynomial,
    MockDcFmProvider, PolynomialDcFmProvider,
};
pub use masking::advanced_masking::{
    apply_advanced_masking, apply_terrain_masking, apply_water_masking, AdvancedMaskResult,
    AdvancedMaskingConfig, AdvancedMaskingProcessor, MaskingMethod,
};
// pub use deburst_algorithmic_fixes::{
//     AlgorithmicDeburstProcessor, BurstOverlap, ComplexSeamWeights, FractionalFIR,
// };
pub use geometry::f64_numerics::{
    DopplerF64, GeographicCoord, GeometryF64, ImageCoord, NewtonRaphsonF64, Point2D, Point3D,
    PrecisionConverter, RadarCoord, RangeDopplerToGeo,
};
// Use canonical unit conversion from constants module
pub use crate::constants::unit_conversion::{
    db_to_linear_inplace, db_to_linear_inplace_par, db_to_linear_parallel_chunked,
    export_db_parallel, linear_to_db_inplace, linear_to_db_inplace_par,
    linear_to_db_parallel_chunked,
};
// Keep IncidenceTerms and Units from fast_unit_conversion (not duplicated)
pub use context::processing_context::{
    BurstContext, CalibrationContext, DemContext, GeometryContext, ProcessingContext,
    ProcessingContextBuilder, ProcessingStage, ProductInfo, ProvenanceValidationResult,
    RadarParameters, StageStatus, TimingContext,
};
pub use deburst::mask_propagation::{
    fill_small_gaps, flags as mask_flags, MaskStats, PixelValidityMask,
    ProcessingStage as MaskProcessingStage,
};
pub use geometry::type_safe_units::{
    angle_ops, distance_ops, time_ops, Degrees, Hertz, Meters, MetersPerSecond, Radians, Seconds,
};
pub use metadata::metadata_parser::{
    CalibrationUnits, CompactCalibrationData, CompactGeometricData, CompactNoiseData,
    CompactTimingData, SinglePassXmlParser, UnitType,
};
pub use metadata::metadata_provenance::{
    InputProvenance, MetadataManager, ProcessingMetadata, QualityMetrics, UncertaintyEstimates,
};
pub use metadata::metadata_strictness::{
    BetaVariabilityResult, CoverageResult, LineIndexingResult, PowerSeamResult,
    StrictMetadataValidator, StrictValidationResult, UnitsValidationResult,
};
pub use metadata::validation_gates::{ValidationGates, ValidationResult, ValidationSummary};
pub use multilook::{MultilookParams, MultilookProcessor};
pub use perf::fast_unit_conversion::{IncidenceTerms, Units};
pub use perf::precision_standards::{
    deterministic_parallel, geometry_precision, initialize_precision_standards, validation,
    DeterministicRng, PrecisionStandards,
};
pub use perf::{
    fast_unit_conversion, memory_optimizations, memory_optimized, precision_standards,
    simd_optimizations,
};
pub use quality::validated_processing_pipeline::{
    DeburstMetrics, ProcessedSarData, ProcessingParameters, ValidatedProcessingPipeline,
};
pub use quality::{QualityAssessment, QualityAssessor, QualityConfig, QualityStatistics};
pub use quality_flags::{
    global_quality_flags, reset_global_quality_flags, GlobalQualityFlags, QualityFlag, QualityFlags,
};
pub use speckle_filter::{SpeckleFilter, SpeckleFilterParams, SpeckleFilterType};
pub use terrain_correction::{
    compute_gradients,
    get_terrain_correction_stats, // Counter reset functions
    gradient_to_normal,
    reset_terrain_correction_counters,
    slope_magnitude,
    DopplerCentroidModel,
    GeocodedPixel,
    GradientOperator, // DUP-4: Centralized gradient ops
    GroundPoint,
    RangeDopplerParams,
    RtcMode,
    RtcQualityFlag,
    TerrainCorrector,
};
pub use topsar_merge::power_preserving_resample::{
    calculate_power_ratio, power_preserving_resample_2d, PowerPreservingResampler,
};
pub use topsar_merge::{
    BlendingMethod, MergeParameters, MergedSwathData, OutputGrid, OverlapQuality, OverlapRegion,
    PerformanceMetrics as MergePerformanceMetrics, QualityControl, QualityResults, TopsarMerge,
};

// Coordinate convention consistency tests (Critical Issue C-3 from forensic audit)
#[cfg(test)]
mod coordinate_convention_test;
