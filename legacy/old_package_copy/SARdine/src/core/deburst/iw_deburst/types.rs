#![allow(dead_code, unused_variables)]
//! Type definitions for IW TOPSAR deburst processing
//!
//! This module contains the core data structures used throughout the deburst pipeline:
//! - `BurstInfo`: Per-burst metadata including timing, geometry, and polynomial coefficients
//! - `DeburstConfig`: Configuration options for the deburst processor
//! - `DeburstResult`: Output container with debursted image and quality metrics
//! - Supporting types for provenance tracking and internal processing

use crate::types::SarComplex;
use ndarray::Array2;
use serde::Serialize;

use super::RangePolynomial;

/// Region-of-interest descriptor for power diagnostics (fractions of burst extents)
#[derive(Debug, Clone)]
pub struct PowerRoi {
    pub name: String,
    pub line_start_frac: f32,
    pub line_end_frac: f32,
    pub sample_start_frac: f32,
    pub sample_end_frac: f32,
}

impl PowerRoi {
    /// Convert fractional bounds into exclusive index ranges within a burst
    pub fn resolved_ranges(
        &self,
        lines: usize,
        samples: usize,
    ) -> Option<(usize, usize, usize, usize)> {
        if lines == 0 || samples == 0 {
            return None;
        }

        let clamp = |v: f32| v.max(0.0).min(1.0);

        let mut ls = (clamp(self.line_start_frac) * lines as f32).floor() as usize;
        if ls >= lines {
            ls = lines.saturating_sub(1);
        }
        let mut le = (clamp(self.line_end_frac) * lines as f32).ceil() as usize;
        if le <= ls {
            le = (ls + 1).min(lines);
        } else if le > lines {
            le = lines;
        }

        let mut ss = (clamp(self.sample_start_frac) * samples as f32).floor() as usize;
        if ss >= samples {
            ss = samples.saturating_sub(1);
        }
        let mut se = (clamp(self.sample_end_frac) * samples as f32).ceil() as usize;
        if se <= ss {
            se = (ss + 1).min(samples);
        } else if se > samples {
            se = samples;
        }

        Some((ls, le, ss, se))
    }
}

// ============================================================================
// BurstInfo - Core burst metadata structure
// ============================================================================

/// TOPSAR Burst information for proper debursting
/// Based on ESA Sentinel-1 Level 1 Detailed Algorithm Definition
#[derive(Debug, Clone)]
pub struct BurstInfo {
    pub burst_id: usize,
    pub start_line: usize,
    pub end_line: usize,
    pub start_sample: usize,
    pub end_sample: usize,
    pub azimuth_time: String,
    pub sensing_time: String,
    pub first_valid_sample: Vec<i32>,
    pub last_valid_sample: Vec<i32>,
    pub byte_offset: u64,

    // TOPSAR-specific parameters
    pub azimuth_fm_rate: f64, // Azimuth FM rate (Hz/s)

    /// Azimuth steering rate (rad/s) - CRITICAL: Units clarification
    ///
    /// In Sentinel-1 annotation XML, <azimuthSteeringRate> represents the BEAM STEERING ANGLE RATE,
    /// not the phase rate. This is the rate at which the antenna beam sweeps across the swath.
    ///
    /// **Physical Interpretation:**
    /// - Beam angle: θ(t) = θ₀ + (dθ/dt) × t
    /// - Doppler shift: f_doppler = (2v/λ) × sin(θ) ≈ (2v/λ) × θ  (small angle)
    /// - Phase: φ(t) = 2π × ∫f_doppler dt = 2π × (2v/λ) × ∫θ dt
    ///
    /// **Usage in Phase Correction:**
    /// The steering contribution to phase is:
    ///   φ_steering(t) = 2π × (2v/λ) × (dθ/dt) × t²/2
    ///
    /// However, in practice, the Doppler centroid polynomial already accounts for most
    /// of the steering effect. This term captures residual steering not in DC polynomial.
    ///
    /// **Current Implementation:**
    /// We use `steering_rate * t` directly in phase, which assumes the annotation
    /// value has been pre-scaled or that the linear approximation is sufficient for
    /// the residual correction.
    ///
    /// **Typical Values:**
    /// - Sentinel-1 IW: ~0.0015-0.0020 rad/s (beam angle rate)
    /// - If converted to phase rate: ~2.6-3.4 MHz (at v=7600 m/s, λ=0.0555m)
    ///
    /// **References:**
    /// - S1-TN-MDA-52-7445 "TOPSAR Debursting Algorithm" Section 2.3
    /// - De Zan & Guarnieri (2006) "TOPSAR Processing"
    pub azimuth_steering_rate: f64,

    pub slant_range_time: f64,      // Slant range time (s)
    pub doppler_centroid: f64,      // Doppler centroid frequency (Hz)
    pub azimuth_bandwidth: f64,     // Processed azimuth bandwidth (Hz)
    pub range_sampling_rate: f64,   // Range sampling rate (Hz)
    pub range_pixel_spacing: f64,   // Range pixel spacing (m)
    pub azimuth_pixel_spacing: f64, // Azimuth pixel spacing (m)

    // NEW: Enhanced timing parameters for scientifically correct processing
    pub azimuth_time_interval: f64, // PRF interval (seconds) from annotation
    pub dc_polynomial: Vec<f64>,    // Doppler centroid polynomial coefficients
    pub fm_polynomial: Vec<f64>,    // FM rate polynomial coefficients

    // Optional range-dependent polynomial tables (slant-range grid)
    pub(crate) dc_range_poly: Option<RangePolynomial>,
    pub(crate) fm_range_poly: Option<RangePolynomial>,

    // CRITICAL FIX (Issue #1, #8): Polynomial reference time for proper phase calculation
    /// Reference time for DC/FM polynomial evaluation (seconds since epoch)
    /// This is the t0 from annotation <dcEstimateList><dcEstimate><t0>
    /// MUST be used to compute: time_offset = burst_sensing_time - polynomial_ref_time
    pub dc_polynomial_t0: Option<f64>,

    /// Reference time for FM polynomial evaluation (seconds since epoch)
    /// Parsed from <fmRate><t0> when available; usually matches DC t0
    pub fm_polynomial_t0: Option<f64>,

    /// Burst reference time for sensing (seconds since epoch)
    /// Used to compute polynomial time offset: sensing_time - polynomial_t0
    pub burst_reference_time_seconds: Option<f64>,

    /// Optional azimuth timestamp parsed directly from <azimuthTime>
    /// Logged separately so we can prove ProductStart vs OrbitEpoch handling
    pub burst_azimuth_time_seconds: Option<f64>,

    /// Burst start time from <azimuthTime> in UTC seconds (for burst span calculation)
    /// Used to compute accurate burst azimuth duration for deramp timing
    pub burst_start_time_utc: Option<f64>,

    /// Next burst start time from <azimuthTime> in UTC seconds
    /// Used with burst_start_time_utc to compute actual burst azimuth span
    pub next_burst_start_time_utc: Option<f64>,

    /// Diagnostics: dcEstimate indices used for this burst after azimuth interpolation
    pub dc_selection_lower_idx: Option<usize>,
    pub dc_selection_upper_idx: Option<usize>,
    pub dc_selection_weight: Option<f64>,

    /// Diagnostics: fmRate indices used for this burst after azimuth interpolation
    pub fm_selection_lower_idx: Option<usize>,
    pub fm_selection_upper_idx: Option<usize>,
    pub fm_selection_weight: Option<f64>,
}

impl BurstInfo {
    /// Calculate the number of lines in this burst
    pub fn lines(&self) -> usize {
        self.end_line.saturating_sub(self.start_line) + 1
    }

    // ========================================================================
    // REMOVED: calculate_dc_aware_overlap (physically incorrect - overlap comes from timing)
    // REMOVED: calculate_deramp_phase (use precompute_deramp_2d instead)
    // REMOVED: calculate_overlap_weight (use overlap_weight function instead)
    // ========================================================================

    /// Evaluate DC polynomial at azimuth time
    pub(crate) fn eval_dc_poly(poly: &[f64], t_az: f64) -> f64 {
        poly.iter()
            .enumerate()
            .fold(0.0, |acc, (i, &coeff)| acc + coeff * t_az.powi(i as i32))
    }

    /// Calculate the number of valid samples for a given line (OPTIMIZATION 2: Enhanced)
    pub fn valid_samples_for_line(&self, line: usize) -> (usize, usize) {
        if line >= self.first_valid_sample.len() || line >= self.last_valid_sample.len() {
            return (0, 0);
        }

        let first = self.first_valid_sample[line].max(0) as usize;
        let last = self.last_valid_sample[line].max(0) as usize;

        if first <= last {
            (first, last)
        } else {
            (0, 0)
        }
    }

    /// Create a BurstInfo with enhanced timing parameters
    pub fn with_enhanced_timing(
        mut self,
        azimuth_time_interval: f64,
        dc_polynomial: Vec<f64>,
        fm_polynomial: Vec<f64>,
    ) -> Self {
        self.azimuth_time_interval = azimuth_time_interval;
        self.dc_polynomial = dc_polynomial;
        self.fm_polynomial = fm_polynomial;
        self
    }
}

// ============================================================================
// DeburstConfig - Processing configuration
// ============================================================================

/// Configuration for TOPSAR deburst processing
/// Based on ESA Sentinel-1 processing specifications with scientific enhancements
#[derive(Debug, Clone)]
pub struct DeburstConfig {
    pub blend_overlap: bool,       // Enable overlap blending between bursts
    pub blend_lines: usize,        // Number of lines to blend (typically 100-200)
    pub remove_invalid_data: bool, // Remove invalid data regions
    pub seamless_stitching: bool,  // Enable seamless stitching with phase continuity
    pub apply_deramp: bool,        // Apply azimuth deramp for TOPSAR
    pub preserve_phase: bool,      // Preserve interferometric phase information

    /// Use SNAP-compatible midpoint selection for burst overlap zones.
    ///
    /// When enabled, each row in the overlap zone is assigned to the burst
    /// whose center is closer (hard 0/1 selection, no blending). This avoids
    /// TOPSAR scalloping artifacts because data is always taken from nearer
    /// the burst center where antenna gain is highest.
    ///
    /// The S1 calibration LUT does NOT include azimuth scalloping correction
    /// (sigma0 values are constant to within 0.001 dB across azimuth), so
    /// cos² blending of burst edges produces visible dark bands.
    ///
    /// **Default:** true (recommended for intensity/backscatter applications)
    /// **Set to false:** for InSAR where phase continuity is needed
    pub use_midpoint_selection: bool,

    /// Use range-dependent deramp for TOPS processing (experimental)
    ///
    /// When enabled, computes per-pixel deramp phase accounting for range-dependent
    /// Doppler centroid and FM rate polynomials: f_DC(t,r) and K_az(t,r).
    ///
    /// **Benefits:**
    /// - Eliminates faint residual stripes in overlap regions
    /// - Required for precise interferometry with range-varying DC
    /// - Improves IW sub-swath alignment by 0.1-0.3 dB RMS
    ///
    /// **Cost:**
    /// - +10-20% processing time (per-pixel phase computation)
    /// - +memory overhead for 2D deramp LUT
    ///
    /// **When to Enable:**
    /// - DC polynomial has range terms (degree > 2 with range coefficients)
    /// - Visible striping artifacts in overlap regions
    /// - Precision applications (InSAR, change detection)
    ///
    /// **When to Disable:**
    /// - Time-only DC polynomials (for quick-look only)
    /// - Quick-look processing
    /// - Testing/debugging (not recommended for production)
    ///
    /// **Default:** true (enabled for IW merging - essential for stripe elimination)
    /// **Changed:** 2025-10-04 - Now default ON for production-grade IW merging
    ///
    /// # References
    /// - De Zan & Guarnieri (2006): "TOPSAR Processing" Section 4.2
    /// - S1-TN-MDA-52-7445: "TOPSAR Debursting" Section 3.3
    pub use_range_dependent_deramp: bool,

    // NEW: Scientific enhancements
    pub use_annotation_timing: bool, // Use annotation PRF instead of velocity estimates
    pub enable_bilinear_interp: bool, // Enable sub-pixel bilinear interpolation
    pub enable_hit_count_mask: bool, // Generate hit-count quality mask
    pub power_preservation_check: bool, // Verify power conservation in non-overlap regions

    // Radiometric continuity
    pub enable_row_equalization: bool, // Equalize per-row (burst) mean power to suppress striping
    pub fill_small_gaps: bool,         // Fill thin zero-coverage seams (burst/subswath gaps)

    // PERFORMANCE: Phase 1 optimizations
    pub enable_loop_unrolling: bool, // Unroll hot loops for better vectorization (10% faster)

    // PERFORMANCE: Phase 2+3 optimizations
    pub enable_simd: bool, // Use SIMD for blend/normalize (6-8× faster on AVX2)
    pub enable_tiled_processing: bool, // Process in cache-friendly tiles (10-15% faster)
    pub tile_lines: usize, // Tile height (default: 2048 for L3 cache)
    pub tile_samples: usize, // Tile width (default: 4096)

    /// Optional ROIs for focused per-burst power diagnostics (fractions 0-1)
    pub power_rois: Vec<PowerRoi>,
}

impl Default for DeburstConfig {
    fn default() -> Self {
        Self {
            blend_overlap: true,
            blend_lines: 150, // Typical TOPSAR overlap
            remove_invalid_data: true,
            seamless_stitching: true,
            apply_deramp: true,               // Essential for TOPSAR
            preserve_phase: true,             // Important for interferometry
            use_midpoint_selection: true,     // SNAP-compatible: select from closer burst center (avoids scalloping)
            use_range_dependent_deramp: true, // Default ON - essential for IW stripe elimination (changed 2025-10-04)

            // NEW: Scientific defaults
            use_annotation_timing: true, // Prefer annotation PRF over velocity estimates
            enable_bilinear_interp: false, // Conservative default for stability
            enable_hit_count_mask: true, // Always generate quality masks
            power_preservation_check: true, // Always verify scientific correctness
            enable_row_equalization: false, // OFF: midpoint selection handles scalloping; equalization may distort radiometry
            fill_small_gaps: false,         // OFF: gap filling can introduce artifacts

            // PERFORMANCE: Phase 1 optimizations (disabled by default for backward compat)
            enable_loop_unrolling: true, // Safe to enable - compiler optimization hint

            // PERFORMANCE: Phase 2+3 optimizations (enabled by default - safe SIMD with fallback)
            enable_simd: true, // Uses AVX2 if available, scalar fallback otherwise
            enable_tiled_processing: false, // Disabled by default, enable for large scenes
            tile_lines: 2048,  // L3 cache optimized (2048 × 4096 × 4 bytes = 32MB)
            tile_samples: 4096,
            power_rois: Vec::new(),
        }
    }
}

// ============================================================================
// DeburstResult - Output container
// ============================================================================

/// Enhanced deburst result with quality metrics and coverage information
#[derive(Debug)]
pub struct DeburstResult {
    pub image: Array2<SarComplex>,
    pub hit_count: Array2<u16>, // Coverage mask: 0 = uncovered, >0 = hit count
    pub power_ratio: f64,       // Power conservation ratio (should be ~1.0)
    pub uncovered_pixels: usize, // Number of pixels with no coverage
    pub blend_quality_score: f64, // Quality score for seamless blending (0-1)
    pub total_azimuth_lines: usize, // Normalized total azimuth lines in output grid
    pub total_range_samples: usize, // Total range samples in output grid
    pub azimuth_index_origin: usize, // Original azimuth index offset removed during deburst
    pub range_sample_origin: usize, // Range origin (first valid sample) removed during deburst
    pub timing_reference: Option<f64>, // Reference time (seconds since epoch) used for relative timing
    pub burst_timing: Vec<BurstTimingInfo>, // Per-burst timing metadata (relative to timing_reference)
    pub row_provenance: Vec<RowRangeProvenance>, // Piecewise mapping from output rows to burst/line
    pub step2_diagnostics: Option<Step2Diagnostics>, // Optional STEP-2 diagnostics summary
}

// ============================================================================
// Supporting types
// ============================================================================

/// Compact mapping from contiguous output rows to a source burst and starting line
#[derive(Clone, Debug)]
pub struct RowRangeProvenance {
    pub out_row_start: usize,
    pub out_row_end: usize, // exclusive
    pub burst_id: usize,
    pub burst_line_start: usize,
}

// ============================================================================
// STEP-2 diagnostics payload (JSON serializable)
// ============================================================================

#[derive(Clone, Debug, Serialize)]
pub struct DerampPhaseStats {
    pub pre_slope: f64,
    pub post_slope: f64,
    pub reduction_factor: f64,
    pub samples_used: usize,
    pub lines_used: usize,
}

#[derive(Clone, Debug, Serialize)]
pub struct Step2Diagnostics {
    pub coverage_fraction: f64,
    pub covered_pixels: usize,
    pub uncovered_pixels: usize,
    pub contrib_histogram: [usize; 3],
    pub overlap_metrics: Vec<crate::core::deburst::diagnostics::OverlapMetrics>,
    pub burst_means: Vec<f64>,
    pub cv: f64,
    pub deramp_stats: Option<DerampPhaseStats>,
}

/// Timing metadata carried forward to merge to avoid reconstructing time axes
#[derive(Clone, Debug)]
pub struct BurstTimingInfo {
    pub burst_id: usize,
    pub prf_hz: f64,
    pub dt: f64,
    pub t_start_rel: f64,
    pub t_end_rel: f64,
    pub line_count_emitted: u32,
}

// ============================================================================
// Internal types (pub(crate) visibility)
// ============================================================================

/// Row-level copy/blend operation for deburst execution
#[derive(Clone, Debug)]
pub(crate) struct DeburstRowSegment {
    pub(crate) burst_idx: usize,
    pub(crate) line_in_burst: usize,
    pub(crate) src_line: usize,
    pub(crate) src_col_start: usize,
    pub(crate) dst_col_start: usize,
    pub(crate) len: usize,
    pub(crate) weight: f32,
}

/// Precomputed deburst plan describing how each output row is populated
#[derive(Clone, Debug)]
pub(crate) struct DeburstPlan {
    pub(crate) rows: usize,
    pub(crate) cols: usize,
    pub(crate) rows_plan: Vec<Vec<DeburstRowSegment>>,
}
