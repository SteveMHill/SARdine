use crate::core::metadata::UnitType;

/// Detect unit type from a string, normalized to lowercase.
pub fn detect_unit_type(units_str: &str) -> UnitType {
    let units_lower = units_str.trim().to_ascii_lowercase();
    match units_lower.as_str() {
        "db" | "decibel" | "decibels" => UnitType::Decibel,
        "db_tenth" | "decibel_tenth" | "decibels_tenth" => UnitType::DecibelTenths,
        "" => UnitType::Linear,
        _ => {
            if units_lower.contains("db") || units_lower.contains("decibel") {
                UnitType::DecibelTenths
            } else {
                UnitType::Linear
            }
        }
    }
}

/// Convert a value to linear units using the detected unit type and an optional scale factor.
pub fn convert_to_linear(value: f64, units: UnitType, scale: f64) -> f64 {
    let scaled = value * scale;
    match units {
        UnitType::Linear => scaled,
        UnitType::Decibel => 10.0f64.powf(scaled / 10.0),
        UnitType::DecibelTenths => 10.0f64.powf((scaled * 0.1) / 10.0),
        UnitType::Unknown => scaled,
    }
}
