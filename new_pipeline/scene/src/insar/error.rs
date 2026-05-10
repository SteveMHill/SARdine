//! Typed error enum for the InSAR pipeline.

use thiserror::Error;

/// Errors produced by InSAR operations.
#[derive(Debug, Error)]
pub enum InsarError {
    /// The subswath has no Doppler centroid estimates, which are required for
    /// TOPS deramping.  This means the scene was parsed from a product type
    /// that does not include a `<dopplerCentroid>` block (e.g. GRD).
    #[error(
        "subswath {swath} has no Doppler centroid estimates; \
         InSAR requires an IW/EW SLC product with a <dopplerCentroid> block"
    )]
    NoDcEstimates { swath: String },

    /// The subswath has no azimuth FM rate estimates.
    #[error(
        "subswath {swath} has no azimuth FM rate estimates; \
         InSAR requires an IW/EW SLC product with an <azimuthFmRateList> block"
    )]
    NoFmRates { swath: String },

    /// The azimuth time of a pixel is before all DC estimates.  This indicates
    /// either a mis-parsed azimuth time or a burst outside the expected span.
    #[error(
        "azimuth time {az_time_s:.6} s (relative) is before the first DC estimate \
         at {first_est_time_s:.6} s"
    )]
    AzimuthTimeBeforeDcEstimates {
        az_time_s: f64,
        first_est_time_s: f64,
    },

    /// The deramp input and output buffers have inconsistent dimensions.
    #[error(
        "deramp: input has {input_lines} lines × {input_samples} samples \
         but output has {output_len} elements (expected {expected})"
    )]
    DimensionMismatch {
        input_lines: usize,
        input_samples: usize,
        output_len: usize,
        expected: usize,
    },
}
