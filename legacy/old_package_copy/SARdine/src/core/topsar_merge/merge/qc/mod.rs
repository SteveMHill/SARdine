//! Quality control modules for TOPSAR merge.

pub mod metrics;
pub mod phase;
pub mod radiometric;

#[allow(unused_imports)]
pub use metrics::calculate_quality_metrics;
#[allow(unused_imports)]
pub use phase::{
    compute_overlap_coherence, compute_overlap_phase_offsets, log_overlap_phase_diagnostics,
    phase_rotation, precompute_overlap_phase_ramps,
};
#[allow(unused_imports)]
pub use radiometric::{
    compute_overlap_gains, isotonic_non_decreasing, smooth_gain_curve,
    validate_radiometric_consistency,
};
