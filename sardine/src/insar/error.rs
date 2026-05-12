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

    /// Forward geocoding (slant-range + zero-Doppler → ECEF on ellipsoid)
    /// failed to converge for a sparse grid point.
    #[error(
        "forward geocoding failed to converge at grid point \
         (line {line}, sample {sample})"
    )]
    ForwardGeocodingFailed { line: usize, sample: usize },

    /// Too few valid sparse grid points to fit the co-registration polynomial.
    /// Need at least 6 for degree-2 2D fit; got `n_valid`.
    #[error(
        "co-registration: need at least 6 valid grid points for polynomial fit, \
         got {n_valid}"
    )]
    CoregGridTooSmall { n_valid: usize },

    /// Zero-Doppler solver failed to converge on the secondary orbit for a
    /// sparse grid point.
    #[error(
        "secondary zero-Doppler solver did not converge at grid point \
         (line {line}, sample {sample})"
    )]
    SecondaryZeroDopplerFailed { line: usize, sample: usize },

    /// Orbit interpolation error propagated from the orbit module.
    #[error("orbit interpolation error: {0}")]
    Orbit(#[from] crate::orbit::OrbitError),

    /// Coherence estimation window size is zero.
    #[error(
        "coherence estimation window size must be > 0, got az={az} rg={rg}"
    )]
    InvalidWindowSize { az: usize, rg: usize },
}
