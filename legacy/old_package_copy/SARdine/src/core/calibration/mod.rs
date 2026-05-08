//! Calibration module broken into focused components: model, units, XML IO, apply, noise, and antenna.

pub mod antenna;
pub mod apply;
pub mod audit;
pub mod io_xml;
pub mod lut;
pub mod model;
pub mod noise;
pub mod parsing;
pub mod units;

pub use model::{
    sane_gain, AntennaPatternLUT, AntennaPatternVector, CalibrationCoefficients,
    CalibrationCoordinateMapper, CalibrationLUT, CalibrationScale, CalibrationUnits,
    CalibrationVector, CompactCalibrationData, CompactGeometricData, CompactNoiseData,
    CompactTimingData, EllipsoidIncidenceModel, IncidenceAngleModel, NoiseAzimuthVector,
    NoiseCoefficients, NoiseCoordinateMapper, NoiseLUT, NoiseLutMode, NoiseVector,
    ValidSampleRanges, MIN_VALID_POWER, POWER_ALERT_THRESHOLD,
};

pub use crate::core::metadata::UnitType;

pub use apply::{
    apply_beta0_calibration, apply_calibration_to_denoised, apply_dn_calibration,
    apply_gamma0_calibration, apply_sigma0_calibration, CalibrationType,
};

pub use audit::sigma_audit_and_maybe_write_json;

pub use noise::{
    apply_fused_noise_calibration, apply_fused_slc_calibration, apply_thermal_noise_removal,
    interpolate_noise_row_for_vector, precompute_noise_lut, precompute_noise_lut_with_strategy,
    NoiseStrategy,
};

pub use antenna::{parse_antenna_pattern_from_xml, precompute_antenna_pattern_lut};

pub use io_xml::{
    parse_antenna_pattern_from_xml as parse_antenna_pattern_from_xml_via_io,
    parse_calibration_from_xml, parse_noise_from_xml,
};

pub use lut::{get_calibration_value, precompute_calibration_lut};

#[cfg(test)]
mod tests;
