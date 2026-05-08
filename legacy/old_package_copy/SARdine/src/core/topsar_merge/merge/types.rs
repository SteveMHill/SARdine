use crate::core::deburst::iw_deburst::{BurstTimingInfo, RowRangeProvenance};
use ndarray::Array2;
use std::fmt;

/// Enhanced merge parameters for scientific processing
#[derive(Debug, Clone)]
pub struct MergeParameters {
    /// Blending method for overlap regions
    pub blending_method: BlendingMethod,
    /// Phase preservation for complex data
    pub preserve_phase: bool,
    /// Overlap region optimization
    pub optimize_overlaps: bool,
    /// Enable parallel processing
    pub enable_parallel: bool,
    /// Chunk size for memory-efficient processing
    pub chunk_size: usize,
    /// Edge feathering width (pixels)
    pub feather_width: usize,
    /// Skip destriping filter (saves ~40s on large scenes)
    /// Destriping is cosmetic - only needed for visual products
    /// Set SARDINE_SKIP_DESTRIPING=1 to skip at runtime
    pub skip_destriping: bool,
}

/// Blending methods for overlap regions
#[derive(Debug, Clone, Default)]
pub enum BlendingMethod {
    /// Linear blending with distance weighting
    #[default]
    Linear,
    /// Gaussian weighted blending
    Gaussian { sigma: f32 },
    /// Advanced multi-scale blending
    MultiScale { scales: Vec<usize> },
    /// Phase-coherent blending for complex data
    PhaseCoherent,
    /// ESA SNAP compatible blending
    SnapCompatible,
}

/// Quality control parameters
#[derive(Debug, Clone)]
pub struct QualityControl {
    /// Enable quality validation
    pub enable_validation: bool,
    /// Maximum phase discontinuity (radians)
    pub max_phase_discontinuity: f32,
    /// Minimum valid pixel ratio in overlaps
    pub min_valid_pixel_ratio: f32,
    /// Enable seamline optimization
    pub optimize_seamlines: bool,
    /// Radiometric consistency tolerance
    pub radiometric_tolerance: f32,
}

/// Weighting strategy for a merge segment
#[derive(Clone, Debug)]
pub enum MergeWeight {
    /// Constant weight applied to every sample in the segment
    Constant(f32),
    /// Weight derived from an overlap region weight matrix
    Overlap {
        overlap_index: usize,
        row: usize,
        col_offset: usize,
        inverse: bool,
    },
}

/// Row segment describing how to copy/blend a portion of a subswath into the output grid
#[derive(Clone, Debug)]
pub struct MergeRowSegment {
    pub swath_idx: usize,
    pub src_row: usize,
    pub src_col_start: usize,
    pub dst_col_start: usize,
    pub len: usize,
    pub weight: MergeWeight,
}

/// Precomputed merge plan decomposing the merge into copy/blend row segments
#[derive(Clone, Debug)]
pub struct MergePlan {
    pub rows: usize,
    pub cols: usize,
    pub rows_plan: Vec<Vec<MergeRowSegment>>,
    /// MANDATORY: Tracks which segments are gap-filled (fabricated) for scientific integrity
    /// Tuple: (dst_row, dst_col_start, len) - identifies gap-filled pixel ranges
    pub gap_filled_segments: Vec<(usize, usize, usize)>,
}

/// Overlap region between two adjacent sub-swaths
#[derive(Debug, Clone)]
pub struct OverlapRegion {
    /// First sub-swath ID (e.g., "IW1")
    pub swath1_id: String,
    /// Second sub-swath ID (e.g., "IW2")
    pub swath2_id: String,
    /// Range extent of overlap in swath1 coordinates
    pub swath1_range_start: usize,
    pub swath1_range_end: usize,
    /// Range extent of overlap in swath2 coordinates
    pub swath2_range_start: usize,
    pub swath2_range_end: usize,
    /// Azimuth extent (common for both swaths)
    pub azimuth_start: usize,
    pub azimuth_end: usize,
    /// Optimized blending weights for smooth transition
    pub weights: Array2<f32>,
    /// Quality metrics for this overlap
    pub quality_metrics: OverlapQuality,
}

/// Quality metrics for overlap regions
#[derive(Debug, Clone, Default)]
pub struct OverlapQuality {
    /// Phase coherence in overlap region
    pub phase_coherence: f32,
    /// Radiometric consistency
    pub radiometric_consistency: f32,
    /// Valid pixel percentage
    pub valid_pixel_percentage: f32,
    /// Seamline quality score
    pub seamline_quality: f32,
}

/// Enhanced azimuth time modeling for scientific-grade TOPSAR processing
/// Implements precise per-line timing and Doppler polynomial evaluation
#[derive(Clone)]
pub struct AzimuthTimingModel {
    /// Pulse Repetition Frequency (Hz)
    pub prf: f64,
    /// Azimuth time interval between lines (seconds)
    pub azimuth_time_interval: f64,
    /// Per-burst timing information
    pub burst_timing: Vec<BurstTiming>,
    /// Reference azimuth time for polynomial evaluation
    pub reference_azimuth_time: f64,
}

/// Optional deburst-derived timing/provenance to bypass annotation reconstruction
#[derive(Clone, Debug)]
pub struct DeburstTimingOverride {
    pub subswath_id: String,
    pub timing_reference: Option<f64>,
    pub burst_timing: Vec<BurstTimingInfo>,
    pub row_provenance: Vec<RowRangeProvenance>,
    pub azimuth_index_origin: usize,
    pub range_sample_origin: usize,
}

impl fmt::Debug for AzimuthTimingModel {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("AzimuthTimingModel")
            .field("prf", &self.prf)
            .field("azimuth_time_interval", &self.azimuth_time_interval)
            .field("burst_timing", &self.burst_timing)
            .field("reference_azimuth_time", &self.reference_azimuth_time)
            .finish()
    }
}

/// Burst-specific timing parameters for enhanced azimuth modeling
#[derive(Debug, Clone)]
pub struct BurstTiming {
    pub subswath_id: String,
    pub burst_id: usize,
    pub burst_index: usize,
    pub azimuth_time_start: f64,
    pub azimuth_time_end: f64,
    pub prf_hz: f64,
    pub azimuth_time_interval: f64,
    pub first_line_merged: usize,
    pub last_line_merged: usize,
    pub sensing_time_center: f64,
    pub azimuth_steering_rate: f64,
}

/// Doppler centroid polynomial for precise frequency modeling
/// Implements: f_dc(η) = c0 + c1*η + c2*η² + c3*η³
/// where η is azimuth time relative to reference
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct DopplerPolynomial {
    pub coefficients: Vec<f64>,
    pub reference_time: f64,
    pub validity_start: f64,
    pub validity_end: f64,
}

/// Inter-subswath alignment information (DC-aware)
/// CRITICAL: Must be calculated AFTER DC-aware deburst to ensure proper alignment
#[derive(Debug, Clone)]
pub struct SubswathAlignment {
    pub azimuth_offset: f64,
    pub range_offset: f64,
    pub confidence: f32,
    pub method: AlignmentMethod,
    pub swath1_id: String,
    pub swath2_id: String,
}

/// Method used for inter-subswath alignment
#[derive(Debug, Clone, PartialEq)]
pub enum AlignmentMethod {
    DopplerBased,
    CrossCorrelation,
    Combined,
    None,
}

/// Output grid definition with enhanced azimuth time modeling
#[derive(Debug, Clone)]
pub struct OutputGrid {
    pub range_samples: usize,
    pub azimuth_samples: usize,
    pub range_pixel_spacing: f64,
    pub azimuth_pixel_spacing: f64,
    pub range_time_start: f64,
    pub azimuth_time_start: f64,
    pub azimuth_timing: AzimuthTimingModel,
}

/// Results from the TOPSAR merge operation
#[derive(Debug)]
pub struct MergedSwathData {
    pub merged_intensity: Array2<f32>,
    pub merged_complex: Option<Array2<num_complex::Complex32>>,
    pub output_grid: OutputGrid,
    pub overlap_regions: Vec<OverlapRegion>,
    pub quality_results: QualityResults,
    pub processing_metadata: ProcessingMetadata,
    pub merged_hitcount: Array2<f32>,
    pub uncovered_mask: Array2<u8>,
    /// MANDATORY: Tracks gap-filled (fabricated) pixels for scientific integrity
    /// Value 1 = pixel created via gap-filling (row duplication/interpolation)
    /// Value 0 = pixel from actual SAR data
    /// CRITICAL: Must be exported to Python for analysis exclusion
    pub gap_filled_mask: Array2<u8>,
}

/// Quality assessment results
#[derive(Debug, Clone)]
pub struct QualityResults {
    pub overall_quality: f32,
    pub phase_preservation: Option<f32>,
    pub radiometric_consistency: f32,
    pub overlap_qualities: Vec<OverlapQuality>,
    pub validation_passed: bool,
    pub warnings: Vec<String>,
}

/// Processing metadata
#[derive(Debug, Clone)]
pub struct ProcessingMetadata {
    pub processing_time_ms: f64,
    pub subswaths_merged: usize,
    pub overlap_regions_processed: usize,
    pub blending_method: BlendingMethod,
    pub total_azimuth_lines: usize,
    pub azimuth_index_origin: usize,
    pub performance_metrics: PerformanceMetrics,
}

/// Performance metrics
#[derive(Debug, Clone)]
pub struct PerformanceMetrics {
    pub total_time_seconds: f64,
    pub peak_memory_mb: f64,
    pub pixels_per_second: f64,
    pub overlap_efficiency: f32,
}

impl Default for MergeParameters {
    fn default() -> Self {
        // Check environment variable for destriping skip
        let skip_destriping = std::env::var("SARDINE_SKIP_DESTRIPING")
            .map(|v| v == "1" || v.to_lowercase() == "true")
            .unwrap_or(false);

        // Check environment variable for merge blending mode
        // SARDINE_MERGE_MODE=sharp|snap → SnapCompatible (SNAP-style midpoint selection)
        // SARDINE_MERGE_MODE=cosine|feather → Linear (cosine feathering - default)
        let blending_method = std::env::var("SARDINE_MERGE_MODE")
            .map(|v| {
                let v_lower = v.to_lowercase();
                if v_lower == "sharp" || v_lower == "snap" || v_lower == "cutoff" {
                    log::info!(
                        "🔧 SARDINE_MERGE_MODE={}: Using SNAP-compatible sharp cutoff for IW merge",
                        v
                    );
                    BlendingMethod::SnapCompatible
                } else if v_lower == "cosine" || v_lower == "feather" || v_lower == "linear" {
                    log::info!(
                        "🔧 SARDINE_MERGE_MODE={}: Using cosine feathering for IW merge",
                        v
                    );
                    BlendingMethod::Linear
                } else {
                    log::warn!(
                        "🔧 SARDINE_MERGE_MODE={}: Unknown mode, using default (cosine feathering)",
                        v
                    );
                    BlendingMethod::Linear
                }
            })
            .unwrap_or(BlendingMethod::Linear);

        Self {
            blending_method,
            preserve_phase: true,
            optimize_overlaps: true,
            enable_parallel: true,
            chunk_size: 4096,
            feather_width: 16,
            skip_destriping,
        }
    }
}

impl Default for QualityControl {
    fn default() -> Self {
        Self {
            enable_validation: true,
            max_phase_discontinuity: 0.5,
            min_valid_pixel_ratio: 0.8,
            optimize_seamlines: true,
            radiometric_tolerance: 0.3, // Relaxed from 0.1: typical S1 IW subswath differences are 0.2-0.4 dB
        }
    }
}

impl Default for QualityResults {
    fn default() -> Self {
        Self {
            overall_quality: 0.0,
            phase_preservation: None,
            radiometric_consistency: 0.0,
            overlap_qualities: Vec::new(),
            validation_passed: false,
            warnings: Vec::new(),
        }
    }
}
