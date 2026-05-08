//! DC/FM rate providers for subswath processing.

pub mod dc_fm;

#[allow(unused_imports)]
pub use dc_fm::{
    build_dc_fm_provider_for_swath, build_dc_fm_provider_map, calculate_azimuth_fm_rate,
    validate_dc_fm_providers, validate_dc_prerequisites,
};
