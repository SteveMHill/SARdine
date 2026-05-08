#![allow(dead_code, unused_variables, deprecated)]
//! Multilook types and configuration
//!
//! Defines the enums, structs, and parameters used by multilook processing.

/// Multilook processing mode
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MultilookMode {
    /// Average intensity: <|I+jQ|²> (for amplitude products σ⁰, β⁰, γ⁰)
    Intensity,
    /// Average complex values: <I+jQ> (for interferometry, coherence)
    Complex,
}

/// Border handling strategy for incomplete blocks at image edges
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BorderMode {
    /// Drop incomplete blocks (loses edge data)
    Drop,
    /// Include partial blocks with normalization (default, most data-preserving)
    Partial,
    /// Zero-pad to complete blocks
    ZeroPad,
    /// Replicate edge values to complete blocks
    EdgeReplicate,
}

/// Comprehensive metadata for multilook processing results
#[derive(Debug, Clone)]
pub struct MultilookMetadata {
    /// Original image dimensions (azimuth, range)
    pub input_dims: (usize, usize),
    /// Output image dimensions (azimuth, range)
    pub output_dims: (usize, usize),
    /// Number of looks applied (azimuth, range)
    pub looks: (usize, usize),
    /// Original pixel spacing in meters (azimuth, range)
    pub input_spacing: (f64, f64),
    /// Output pixel spacing in meters (azimuth, range)
    pub output_spacing: (f64, f64),
    /// Processing mode used
    pub mode: MultilookMode,
    /// Border handling mode
    pub border_mode: BorderMode,
    /// Fraction of valid (non-NaN) pixels in output [0.0-1.0]
    pub valid_pixel_fraction: f64,
    /// Theoretical number of looks
    pub theoretical_enl: usize,
    /// Estimated ENL from output statistics
    pub estimated_enl: f32,
    /// Power preservation ratio (output_power / input_power, should be ~1.0)
    pub power_ratio: f64,
    /// Number of output blocks with all-NaN (invalid) samples
    pub nan_blocks: usize,
    /// Average number of valid samples per block
    pub avg_samples_per_block: f64,
    /// Processing time in seconds
    pub processing_time_s: f64,
}

/// Multilooking parameters for speckle reduction
#[derive(Debug, Clone)]
pub struct MultilookParams {
    /// Number of looks in range direction
    pub range_looks: usize,
    /// Number of looks in azimuth direction
    pub azimuth_looks: usize,
    /// Processing mode (intensity or complex)
    pub mode: MultilookMode,
    /// DEPRECATED for intensity multilook (no effect since v0.8.0).
    /// For intensity data, standard multilook always uses simple mean of valid
    /// samples, which is the scientifically correct behavior. The previous
    /// fill_factor scaling was incorrect (reduced power for partial blocks).
    ///
    /// For COMPLEX multilook (InSAR/PolSAR), this flag still applies:
    /// If true, scales by sqrt(N) to preserve expected power for random-phase
    /// (incoherent) targets after coherent averaging.
    ///
    /// See audit issue ML-1 and CHANGELOG for details.
    #[deprecated(
        since = "0.8.0",
        note = "No effect for intensity multilook. Use default."
    )]
    pub preserve_power: bool,
    /// Border handling strategy
    pub border_mode: BorderMode,
    /// Include partial windows at edges (keeps more data)
    pub include_partial: bool,
}

impl Default for MultilookParams {
    fn default() -> Self {
        Self {
            range_looks: 4,
            azimuth_looks: 1,
            mode: MultilookMode::Intensity,
            preserve_power: true, // ESA/SNAP standard for radiometric accuracy
            border_mode: BorderMode::Partial,
            include_partial: true, // Default to keeping edge data
        }
    }
}
