pub mod burst_ops;
pub mod diagnostics;
pub mod diagnostics_integration;
pub mod ew_deburst;
pub mod geometry;
pub mod iw_deburst;
pub mod mask_propagation;
pub mod quality;

// Public API surface - types from the new modular structure
pub use iw_deburst::{
    extract_subswath_complex_data, BurstInfo, DeburstConfig, DeburstProcessor, DeburstResult,
    TopSarDeburstProcessor,
};

// STEP-2 diagnostics
pub use diagnostics::{
    DiagnosticsConfig, Step2Diagnostics, SubswathDiagnostics, Step2DiagnosticsSummary,
    compute_coverage_metrics, compute_contribution_histogram, compute_overlap_metrics,
    compute_burst_power_stats, compute_burst_power_cv, validate_dc_fm_bracket,
    measure_deramp_effectiveness, validate_calibration_lut,
};

pub use diagnostics_integration::DeburstDiagnosticsCollector;

#[cfg(test)]
mod tests;
