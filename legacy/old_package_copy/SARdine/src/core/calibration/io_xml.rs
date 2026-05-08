use crate::core::calibration::parsing::{
    parse_antenna_pattern_from_xml as parse_antenna_pattern,
    parse_calibration_from_xml as parse_calibration, parse_noise_from_xml as parse_noise,
};
use crate::types::SarResult;

use super::model::{AntennaPatternVector, CalibrationCoefficients, NoiseCoefficients};

/// Parse calibration XML into `CalibrationCoefficients` using the legacy robust parser.
pub fn parse_calibration_from_xml(xml_content: &str) -> SarResult<CalibrationCoefficients> {
    parse_calibration(xml_content)
}

/// Parse noise XML into `NoiseCoefficients` using the legacy robust parser.
pub fn parse_noise_from_xml(xml_content: &str) -> SarResult<NoiseCoefficients> {
    parse_noise(xml_content)
}

/// Parse antenna pattern vectors from annotation XML using the legacy robust parser.
pub fn parse_antenna_pattern_from_xml(
    xml_content: &str,
    range_sampling_rate_hz: Option<f64>,
) -> SarResult<Vec<AntennaPatternVector>> {
    parse_antenna_pattern(xml_content, range_sampling_rate_hz)
}
