use crate::types::{SarError, SarResult};
use chrono::{DateTime, NaiveDateTime, Utc};
use quick_xml::events::Event;
use quick_xml::Reader;
use std::sync::atomic::{AtomicUsize, Ordering as AtomicOrdering};

use super::{
    AntennaPatternVector, CalibrationCoefficients, CalibrationVector, NoiseAzimuthVector,
    NoiseCoefficients, NoiseVector,
};

static ALL_BETA_WARN_COUNT: AtomicUsize = AtomicUsize::new(0);
const ALL_BETA_WARN_LIMIT: usize = 25;

/// Numerical stability floor - DEPRECATED, kept for reference
/// Denoising now preserves negative values (noise > signal) so that
/// calibration can properly mark them as NaN (invalid data).
/// The MIN_VALID_POWER threshold in calibration handles invalid pixels.
const _NOISE_FLOOR_DEPRECATED: f32 = 1e-4;

/// Swath-specific incidence angle ranges for Sentinel-1 IW mode
///
/// These are typical values derived from ESA documentation and real product analysis:
/// - IW1 (near range): 29.1° - 36.8°
/// - IW2 (mid range):  36.8° - 42.5°
/// - IW3 (far range):  42.5° - 46.0°
///
/// Reference: ESA Sentinel-1 Product Specification S1-RS-MDA-52-7441
pub struct SwathIncidenceAngles {
    pub min_deg: f32,
    pub max_deg: f32,
}

impl SwathIncidenceAngles {
    /// Get incidence angle range for a given swath ID
    ///
    /// Returns swath-specific angles if recognized, otherwise returns IW-wide defaults
    pub fn for_swath(swath_id: &str) -> Self {
        match swath_id.to_uppercase().as_str() {
            "IW1" => Self {
                min_deg: 29.1,
                max_deg: 36.8,
            },
            "IW2" => Self {
                min_deg: 36.8,
                max_deg: 42.5,
            },
            "IW3" => Self {
                min_deg: 42.5,
                max_deg: 46.0,
            },
            // EW mode swaths
            "EW1" => Self {
                min_deg: 19.0,
                max_deg: 27.0,
            },
            "EW2" => Self {
                min_deg: 27.0,
                max_deg: 32.0,
            },
            "EW3" => Self {
                min_deg: 32.0,
                max_deg: 37.0,
            },
            "EW4" => Self {
                min_deg: 37.0,
                max_deg: 41.0,
            },
            "EW5" => Self {
                min_deg: 41.0,
                max_deg: 47.0,
            },
            // Stripmap mode (single swath with full range)
            s if s.starts_with("S") => Self {
                min_deg: 20.0,
                max_deg: 45.0,
            },
            // Default: full IW range for unknown swaths
            _ => {
                log::debug!(
                    "Unknown swath '{}', using full IW range for incidence angles",
                    swath_id
                );
                Self {
                    min_deg: 29.1,
                    max_deg: 46.0,
                }
            }
        }
    }

    /// Get incidence angle range from calibration file path (extracts swath from filename)
    ///
    /// Expected format: *-iw1-* or *-iw2-* or *-iw3-* (case insensitive)
    pub fn from_path(path: &str) -> Self {
        let lower = path.to_lowercase();
        if lower.contains("-iw1-") || lower.contains("/iw1/") {
            Self::for_swath("IW1")
        } else if lower.contains("-iw2-") || lower.contains("/iw2/") {
            Self::for_swath("IW2")
        } else if lower.contains("-iw3-") || lower.contains("/iw3/") {
            Self::for_swath("IW3")
        } else if lower.contains("-ew") {
            // Try to extract EW swath number
            for i in 1..=5 {
                if lower.contains(&format!("-ew{}-", i)) {
                    return Self::for_swath(&format!("EW{}", i));
                }
            }
            Self::for_swath("")
        } else {
            Self::for_swath("")
        }
    }
}

/// Robustly parse ISO 8601 timestamp to seconds since Unix epoch
/// Supports variants with/without fractional seconds and timezone suffix
fn parse_azimuth_time_to_seconds(time_str: &str) -> SarResult<f64> {
    let s = time_str.trim();

    // Detect if timezone info is present after the 'T' separator
    let post_t = s.split('T').nth(1).unwrap_or("");
    let has_tz = post_t.contains('Z') || post_t.contains('+') || post_t.contains('-');

    // Prefer handling naive timestamps first by assuming UTC (append Z)
    if !has_tz {
        // Try with fractional seconds (or none) and Z suffix
        if let Ok(dt) = DateTime::parse_from_str(&format!("{}Z", s), "%Y-%m-%dT%H:%M:%S%.fZ") {
            let utc = dt.with_timezone(&Utc);
            return Ok(utc.timestamp() as f64 + utc.timestamp_subsec_nanos() as f64 / 1e9);
        }

        // Fallback: naive datetime with fractional seconds (assume UTC)
        if let Ok(naive) = NaiveDateTime::parse_from_str(s, "%Y-%m-%dT%H:%M:%S%.f") {
            let utc: DateTime<Utc> = naive.and_utc();
            return Ok(utc.timestamp() as f64 + utc.timestamp_subsec_nanos() as f64 / 1e9);
        }

        // Fallback: naive datetime without fractional seconds (assume UTC)
        if let Ok(naive) = NaiveDateTime::parse_from_str(s, "%Y-%m-%dT%H:%M:%S") {
            let utc: DateTime<Utc> = naive.and_utc();
            return Ok(utc.timestamp() as f64 + utc.timestamp_subsec_nanos() as f64 / 1e9);
        }
    }

    // Try RFC3339 (with timezone like Z or offsets)
    if let Ok(dt) = DateTime::parse_from_rfc3339(s) {
        let utc = dt.with_timezone(&Utc);
        return Ok(utc.timestamp() as f64 + utc.timestamp_subsec_nanos() as f64 / 1e9);
    }

    // Explicit formats commonly seen in Sentinel-1 XMLs
    if let Ok(dt) = DateTime::parse_from_str(s, "%Y-%m-%dT%H:%M:%S%.f+00:00") {
        let utc = dt.with_timezone(&Utc);
        return Ok(utc.timestamp() as f64 + utc.timestamp_subsec_nanos() as f64 / 1e9);
    }
    if let Ok(dt) = DateTime::parse_from_str(&format!("{}Z", s), "%Y-%m-%dT%H:%M:%SZ") {
        let utc = dt.with_timezone(&Utc);
        return Ok(utc.timestamp() as f64 + utc.timestamp_subsec_nanos() as f64 / 1e9);
    }

    Err(SarError::Processing(format!(
        "Failed to parse azimuth time '{}' with supported ISO-8601 variants",
        time_str
    )))
}

/// Helper function to parse comma/whitespace-separated float list with robust formatting
fn parse_float_list(raw: &str) -> Vec<f64> {
    // tolerate commas and mixed whitespace, scientific notation, etc.
    raw.replace(',', " ")
        .split_ascii_whitespace()
        .filter_map(|s| s.parse::<f64>().ok())
        .collect()
}

/// Helper function to parse comma/whitespace-separated usize list
fn parse_usize_list(raw: &str) -> Vec<usize> {
    raw.replace(',', " ")
        .split_ascii_whitespace()
        .filter_map(|s| s.parse::<usize>().ok())
        .collect()
}

/// Fast min/max/median (median via nth selection) without full sort.
fn stats_min_max_median(vals: &[f64]) -> (f64, f64, f64) {
    if vals.is_empty() {
        return (0.0, 0.0, 0.0);
    }
    let mut min_val = f64::INFINITY;
    let mut max_val = f64::NEG_INFINITY;
    for &v in vals.iter() {
        if v < min_val {
            min_val = v;
        }
        if v > max_val {
            max_val = v;
        }
    }
    let mut buf = vals.to_vec();
    let mid = buf.len() / 2;
    buf.select_nth_unstable_by(mid, |a, b| a.total_cmp(b));
    let median = buf[mid];
    (median, min_val, max_val)
}

/// Robust unit conversion with automatic detection of tenths-dB format and negative conversion
fn convert_units_auto(tag: &str, vals: &mut Vec<f64>, units: Option<&str>) {
    let u = units.map(|s| s.to_ascii_lowercase());
    let verbose = std::env::var("SARDINE_CAL_VERBOSE")
        .ok()
        .as_deref()
        .map(|v| v == "1" || v.eq_ignore_ascii_case("true"))
        .unwrap_or(false);

    if vals.is_empty() {
        return;
    }

    // Calculate value statistics for debugging
    let (median, min_val, max_val) = stats_min_max_median(vals);

    if verbose {
        log::info!(
            "📊 {} units={:?}: median={:.3e}, range=[{:.3e}, {:.3e}], count={}",
            tag,
            units,
            median,
            min_val,
            max_val,
            vals.len()
        );
    }

    // FIXED: Only convert when unit attribute explicitly says dB
    // Remove dangerous heuristic that misclassifies S1 LUTs
    if matches!(
        u.as_deref(),
        Some("db") | Some("decibel") | Some("decibels")
    ) {
        // BOTTLENECK FIX: Use fast exp() with standard 10^(x/10) formula
        const K: f64 = std::f64::consts::LN_10 / 10.0;
        for v in vals.iter_mut() {
            *v = (K * *v).exp();
        }
        if verbose {
            log::info!(
                "{}: converted from dB to linear ({} values)",
                tag,
                vals.len()
            );
        }
        return;
    }

    // Special case for gamma: Sentinel-1 gamma LUT values are also LINEAR calibration constants
    // Similar to sigmaNought/betaNought - typical range 200-1000
    if tag == "gamma" {
        if median >= 1e6 {
            log::error!("🚨 {} has astronomical median={:.3e} - likely double conversion error, NOT converting", tag, median);
            return;
        } else if median >= 50.0 && median <= 5000.0 {
            // Normal range for Sentinel-1 calibration constants - keep as linear
            if verbose {
                log::info!(
                    "✅ {} values are linear calibration constants (median={:.1}), keeping as-is",
                    tag,
                    median
                );
            }
            return;
        } else if median > 5.0 && median < 50.0 && median != median.round() {
            // Possible dB values
            log::warn!(
                "⚠️ {} has suspicious median {:.2} - possible dB encoding, converting to linear",
                tag,
                median
            );
            const K: f64 = std::f64::consts::LN_10 / 10.0;
            for v in vals.iter_mut() {
                *v = (K * *v).exp();
            }
            return;
        } else {
            // Keep as-is
            log::warn!(
                "⚠️ {} has unusual median {:.3e} - keeping as linear, verify calibration",
                tag,
                median
            );
            return;
        }
    }

    // Special case for sigma0/beta0: Sentinel-1 XML sigma/beta values are calibration constants (K values)
    //
    // IMPORTANT: ESA Sentinel-1 Level-1 Product Specification states:
    // - sigmaNought, betaNought values are LINEAR calibration constants (not dB!)
    // - Typical range: 200-2000 for Sentinel-1 IW mode
    // - Formula: σ⁰ = |DN|² / A_σ² where A_σ is the LUT value
    //
    // Previously this code incorrectly detected values in range 100-500 as "tenths-dB"
    // and converted them, causing ~12 dB calibration error.
    if tag == "sigmaNought" || tag == "betaNought" {
        if median >= 1e6 {
            log::error!("🚨 {} has astronomical median={:.3e} - likely double conversion error, NOT converting", tag, median);
            return;
        } else if median >= 50.0 && median <= 5000.0 {
            // Normal range for Sentinel-1 calibration constants - keep as linear
            if verbose {
                log::info!(
                    "✅ {} values are linear calibration constants (median={:.1}), keeping as-is",
                    tag,
                    median
                );
            }
            return;
        } else if median > 5.0 && median < 50.0 && median != median.round() {
            // Possible dB values (unlikely for Sentinel-1 but handle edge case)
            log::warn!(
                "⚠️ {} has suspicious median {:.2} - possible dB encoding, converting to linear",
                tag,
                median
            );
            const K: f64 = std::f64::consts::LN_10 / 10.0;
            for v in vals.iter_mut() {
                *v = (K * *v).exp();
            }
            return;
        } else {
            // Keep other values as-is with a warning
            log::warn!(
                "⚠️ {} has unusual median {:.3e} - keeping as linear, verify calibration",
                tag,
                median
            );
            return;
        }
    }

    // REMOVED: Dangerous heuristic that converts linear LUTs incorrectly
    // For Sentinel-1, treat missing unit attribute as linear (default)
    if verbose {
        log::info!("{}: keeping as linear domain ({} values)", tag, vals.len());
    }
}

/// Reject impossible LUT scales early to prevent astronomical outputs
fn sanity_check_lut(name: &str, vals: &[f64]) -> Result<(), SarError> {
    if vals.is_empty() {
        return Ok(());
    }
    let (med_raw, min_val, max_val) = stats_min_max_median(vals);
    let med = med_raw.max(1e-30);

    // Different thresholds for different calibration types
    // Note: These are calibration constants (K values), not backscatter values
    let threshold = match name {
        "gamma" => 1e6, // Gamma coefficients can be very large after conversion
        "sigmaNought" | "betaNought" => 1e4, // Sigma/Beta calibration constants up to ~10000 are normal
        _ => 1e4,                            // Allow larger calibration constants by default
    };

    log::info!(
        "🔍 Sanity check for {}: median={:.3e}, min={:.3e}, max={:.3e}, threshold={:.0e}",
        name,
        med,
        min_val,
        max_val,
        threshold
    );

    if med > threshold {
        return Err(SarError::Processing(format!(
            "Unrealistic {} LUT scale (median ~ {:.3e}) — expected < {:.0e} (units likely dB/tenths-dB). Refusing to proceed.",
            name, med, threshold
        )));
    }
    Ok(())
}

/// Detect if a calibration LUT is "flat" (constant values indicating parsing failure)
fn is_flat(vals: &[f32]) -> bool {
    if vals.len() < 16 {
        return false;
    }
    
    // Check if all values are zero or near-zero (corrupted data)
    // DO NOT flag constant non-zero values (e.g., beta=237 in IW mode is CORRECT!)
    let mean = vals.iter().copied().sum::<f32>() / vals.len() as f32;
    
    // If mean is essentially zero, it's flat/corrupted
    if mean.abs() < 1e-6 {
        return true;
    }
    
    // If values are constant but non-zero, they're VALID (not flat)
    // Example: beta=237.0 everywhere in Sentinel-1 IW SLC is correct
    false
}

/// Validate calibration vector sanity (lengths, monotonicity)
fn vector_sanity(vector: &CalibrationVector) -> Result<(), SarError> {
    // Pixel alignment
    if vector.pixels.len() < 2 {
        return Err(SarError::Processing("Too few pixel knots".to_string()));
    }
    if vector.pixels.windows(2).any(|w| w[0] >= w[1]) {
        return Err(SarError::Processing("Non-monotonic pixels".to_string()));
    }

    // Length equality (only check those that exist)
    for (name, len) in [
        ("sigma", vector.sigma_nought.len()),
        ("beta", vector.beta_nought.len()),
        ("gamma", vector.gamma.len()),
    ] {
        if len != 0 && len != vector.pixels.len() {
            return Err(SarError::Processing(format!("{} length mismatch", name)));
        }
    }
    Ok(())
}

/// Quick extraction of swath ID from calibration XML header
/// Returns the swath (IW1, IW2, IW3, etc.) without full parsing
fn extract_swath_from_xml(xml_content: &str) -> Option<String> {
    // Quick regex-like search for <swath>XXX</swath> in adsHeader
    if let Some(start) = xml_content.find("<swath>") {
        let content_start = start + 7; // len("<swath>")
        if let Some(end) = xml_content[content_start..].find("</swath>") {
            let swath = xml_content[content_start..content_start + end].trim();
            if !swath.is_empty() {
                return Some(swath.to_string());
            }
        }
    }
    None
}

/// Parse calibration data from Sentinel-1 XML file with proper unit handling
pub fn parse_calibration_from_xml(xml_content: &str) -> SarResult<CalibrationCoefficients> {
    log::info!("🔬 Parsing calibration XML ({} bytes)", xml_content.len());

    // Extract swath ID first for accurate incidence angle recovery
    let swath_id = extract_swath_from_xml(xml_content);
    if let Some(ref swath) = swath_id {
        log::debug!(
            "📡 Detected swath: {} (will use swath-specific incidence angles)",
            swath
        );
    }

    // Use the new robust XML parser with swath info
    match parse_calibration_vectors_robust(xml_content, true, swath_id.as_deref()) {
        Ok(robust_vectors) => {
            log::info!(
                "✅ NEW ROBUST PARSER SUCCESS: {} vectors",
                robust_vectors.len()
            );

            // Convert to CalibrationCoefficients structure
            let mut calibration = CalibrationCoefficients::new();
            calibration.vectors = robust_vectors;

            // Set swath from extracted value
            if let Some(swath) = swath_id {
                calibration.swath = swath;
            }

            // Parse header metadata using the old parser (header parsing is less complex)
            let mut reader = Reader::from_str(xml_content);
            reader.trim_text(true);
            let mut buf = Vec::new();
            let mut in_ads_header = false;
            #[allow(unused_assignments)]
            let mut current_tag = String::new();
            let mut text_content = String::new();

            loop {
                match reader.read_event_into(&mut buf) {
                    Ok(Event::Start(ref e)) => {
                        current_tag = String::from_utf8_lossy(e.name().as_ref()).to_string();
                        text_content.clear();

                        if current_tag == "adsHeader" {
                            in_ads_header = true;
                        }
                    }
                    Ok(Event::End(ref e)) => {
                        let tag_name = String::from_utf8_lossy(e.name().as_ref()).to_string();

                        if in_ads_header {
                            match tag_name.as_str() {
                                "polarisation" => {
                                    calibration.polarization = text_content.trim().to_string()
                                }
                                "swath" => calibration.swath = text_content.trim().to_string(),
                                "startTime" => {
                                    calibration.product_first_line_utc_time =
                                        text_content.trim().to_string()
                                }
                                "stopTime" => {
                                    calibration.product_last_line_utc_time =
                                        text_content.trim().to_string()
                                }
                                "absoluteCalibrationConstant" => {
                                    if let Ok(abs_const) = text_content.trim().parse::<f64>() {
                                        calibration.abs_const = abs_const;
                                        log::info!(
                                            "📊 Parsed absoluteCalibrationConstant: {}",
                                            abs_const
                                        );
                                    }
                                }
                                "adsHeader" => in_ads_header = false,
                                _ => {}
                            }
                        }
                    }
                    Ok(Event::Text(ref e)) => {
                        let new_text = e
                            .unescape()
                            .map_err(|e| {
                                SarError::Processing(format!("Failed to unescape XML text: {}", e))
                            })?
                            .to_string();
                        text_content.push_str(&new_text);
                    }
                    Ok(Event::Eof) => break,
                    Err(e) => {
                        log::warn!("XML parsing warning during header extraction: {}", e);
                        break;
                    }
                    _ => {}
                }
                buf.clear();
            }

            if calibration.abs_const != 0.0 && calibration.abs_const != 1.0 {
                log::info!(
                    "📊 Loaded absoluteCalibrationConstant = {} (applied during LUT inversion)",
                    calibration.abs_const
                );
            }

            log::info!(
                "🎯 ROBUST PARSER: Parsed {} vectors for swath {} polarization {}",
                calibration.vectors.len(),
                calibration.swath,
                calibration.polarization
            );

            return Ok(calibration);
        }
        Err(e) => {
            log::error!("❌ NEW ROBUST PARSER FAILED: {}", e);

            // Check if the failure is due to identical values detected by robust parser
            if e.to_string().contains("identical beta coefficients")
                || e.to_string().contains("XML parsing regression")
            {
                // This is actually the robust parser working correctly - it's detecting the issue!
                log::error!("🚨 ROBUST PARSER CORRECTLY DETECTED XML PARSING BUG!");
                log::error!("🚨 Root cause: XML parsing producing identical coefficients (β⁰ = 237 problem)");
                log::error!(
                    "🚨 This indicates either corrupted XML or a fundamental parsing issue"
                );

                // Don't fall back to legacy parser - it has the same bug
                return Err(e);
            }

            // Enhanced error handling: don't fall back to legacy parser
            // The comprehensive fixes make the robust parser reliable
            log::error!(
                "❌ Robust XML parser failed, and legacy parser has been removed for safety"
            );
            return Err(e);
        }
    }
}

/// NEW ROBUST XML PARSER - addresses root cause of coefficient collapse (237 everywhere)
///
/// # Arguments
/// * `xml_content` - The calibration XML content
/// * `debug` - Enable debug logging
/// * `swath_id` - Optional swath ID (IW1, IW2, IW3) for swath-specific incidence angle recovery
fn parse_calibration_vectors_robust(
    xml_content: &str,
    debug: bool,
    swath_id: Option<&str>,
) -> Result<Vec<CalibrationVector>, SarError> {
    let mut reader = Reader::from_str(xml_content);
    reader.trim_text(true);

    let mut vectors = Vec::new();
    let mut buf = Vec::new();

    // Per-element state tracking
    let mut current_vector: Option<CalibrationVector> = None;
    let mut in_calibration_vector = false;

    // Per-element tag tracking with independent unit/scale scoping
    let mut in_pixel = false;
    let mut in_sigma = false;
    let mut in_beta = false;
    let mut in_gamma = false;
    let mut in_line = false;
    let mut in_azimuth_time = false;

    // Per-element units and scale (reset for each element)
    let mut current_units: Option<String> = None;
    let mut current_scale: Option<f64> = None;

    // Get swath-specific incidence angles for coefficient recovery
    let inc_angles = swath_id
        .map(SwathIncidenceAngles::for_swath)
        .unwrap_or_else(|| SwathIncidenceAngles::for_swath(""));

    // Text accumulation for robust parsing
    let mut accumulated_text = String::new();

    loop {
        match reader.read_event_into(&mut buf) {
            Ok(Event::Start(ref e)) => {
                let name_bytes = e.name().as_ref().to_vec();
                let name = String::from_utf8_lossy(&name_bytes);
                match name.as_ref() {
                    "calibrationVector" => {
                        in_calibration_vector = true;
                        current_vector = Some(CalibrationVector {
                            azimuth_time: String::new(),
                            line: 0,
                            pixels: Vec::new(),
                            sigma_nought: Vec::new(),
                            beta_nought: Vec::new(),
                            gamma: Vec::new(),
                            dn: Vec::new(),
                            beta_flat: false,
                            sigma_flat: false,
                            gamma_flat: false,
                        });

                        // Note: In actual S1 calibration XML, azimuthTime and line are child elements, not attributes
                        // They will be parsed in the Text event handler below
                    }
                    // Per-element parsing with independent unit/scale scoping
                    "pixel" => {
                        in_pixel = true;
                        current_units = None; // Reset for this element
                        current_scale = None; // Reset for this element
                        accumulated_text.clear();

                        // Parse per-element attributes
                        for attr in e.attributes() {
                            if let Ok(attr) = attr {
                                let key = String::from_utf8_lossy(attr.key.as_ref());
                                let value = String::from_utf8_lossy(&attr.value);
                                match key.as_ref() {
                                    "units" => current_units = Some(value.to_string()),
                                    "scale" => {
                                        current_scale = value.parse().ok();
                                    }
                                    _ => {}
                                }
                            }
                        }
                    }
                    "line" => {
                        in_line = true;
                        accumulated_text.clear();
                    }
                    "azimuthTime" => {
                        in_azimuth_time = true;
                        accumulated_text.clear();
                    }
                    "sigmaNought" => {
                        in_sigma = true;
                        current_units = None; // Reset for this element
                        current_scale = None; // Reset for this element
                        accumulated_text.clear();

                        // Parse per-element attributes
                        for attr in e.attributes() {
                            if let Ok(attr) = attr {
                                let key = String::from_utf8_lossy(attr.key.as_ref());
                                let value = String::from_utf8_lossy(&attr.value);
                                match key.as_ref() {
                                    "units" => current_units = Some(value.to_string()),
                                    "scale" => {
                                        current_scale = value.parse().ok();
                                    }
                                    _ => {}
                                }
                            }
                        }
                    }
                    "betaNought" => {
                        in_beta = true;
                        current_units = None; // Reset for this element
                        current_scale = None; // Reset for this element
                        accumulated_text.clear();

                        // Parse per-element attributes
                        for attr in e.attributes() {
                            if let Ok(attr) = attr {
                                let key = String::from_utf8_lossy(attr.key.as_ref());
                                let value = String::from_utf8_lossy(&attr.value);
                                match key.as_ref() {
                                    "units" => current_units = Some(value.to_string()),
                                    "scale" => {
                                        current_scale = value.parse().ok();
                                    }
                                    _ => {}
                                }
                            }
                        }
                    }
                    "gamma" => {
                        in_gamma = true;
                        current_units = None; // Reset for this element
                        current_scale = None; // Reset for this element
                        accumulated_text.clear();

                        // Parse per-element attributes
                        for attr in e.attributes() {
                            if let Ok(attr) = attr {
                                let key = String::from_utf8_lossy(attr.key.as_ref());
                                let value = String::from_utf8_lossy(&attr.value);
                                match key.as_ref() {
                                    "units" => current_units = Some(value.to_string()),
                                    "scale" => {
                                        current_scale = value.parse().ok();
                                    }
                                    _ => {}
                                }
                            }
                        }
                    }
                    _ => {}
                }
            }
            Ok(Event::Text(e)) => {
                // Accumulate text content across all text events
                if let Ok(text) = e.unescape() {
                    accumulated_text.push_str(&text);
                }
            }
            Ok(Event::End(ref e)) => {
                // Only process accumulated text at End tag
                let name_bytes = e.name().as_ref().to_vec();
                let name = String::from_utf8_lossy(&name_bytes);
                match name.as_ref() {
                    "pixel" if in_pixel => {
                        if let Some(ref mut vector) = current_vector {
                            vector.pixels = parse_usize_list(&accumulated_text);
                        }
                        in_pixel = false;
                    }
                    "line" if in_line => {
                        if let Some(ref mut vector) = current_vector {
                            vector.line = accumulated_text.trim().parse::<i32>().map_err(|e| {
                                SarError::Processing(format!("Failed to parse line number: {}", e))
                            })?;
                        }
                        in_line = false;
                    }
                    "azimuthTime" if in_azimuth_time => {
                        if let Some(ref mut vector) = current_vector {
                            vector.azimuth_time = accumulated_text.trim().to_string();
                        }
                        in_azimuth_time = false;
                    }
                    "sigmaNought" if in_sigma => {
                        if let Some(ref mut vector) = current_vector {
                            let mut values = parse_float_list(&accumulated_text);

                            // Debug: Show original values before any processing
                            if !values.is_empty() {
                                let preview: Vec<_> = values.iter().take(3).collect();
                                log::warn!(
                                    "🔍 SIGMA PROCESSING START | Raw values: {:?} | units: {:?} | scale: {:?}",
                                    preview, current_units, current_scale
                                );
                            }

                            // Apply per-element scale if present
                            if let Some(scale) = current_scale {
                                for v in &mut values {
                                    *v *= scale;
                                }
                            }

                            // Handle robust unit conversion with tenths-dB detection
                            convert_units_auto(
                                "sigmaNought",
                                &mut values,
                                current_units.as_ref().map(|s| s.as_str()),
                            );

                            // Debug: Show values and statistics after unit conversion
                            if !values.is_empty() {
                                let preview: Vec<_> = values.iter().take(3).collect();
                                let mut sorted = values.clone();
                                sorted.sort_by(|a, b| a.total_cmp(b));
                                let median = sorted[sorted.len() / 2];
                                let vmax = *sorted.last().unwrap_or(&median);
                                let vmin = *sorted.first().unwrap_or(&median);
                                log::warn!(
                                    "🔍 SIGMA AFTER UNIT CONVERSION: {:?} | median={:.3e} | range=[{:.3e},{:.3e}] | n={}",
                                    preview, median, vmin, vmax, values.len()
                                );
                            }

                            vector.sigma_nought = values.into_iter().map(|v| v as f32).collect();

                            // Parsed sigma nought values
                        }
                        in_sigma = false;
                    }
                    "betaNought" if in_beta => {
                        if let Some(ref mut vector) = current_vector {
                            let mut values = parse_float_list(&accumulated_text);

                            // Debug: Show original values before any processing
                            if !values.is_empty() {
                                let preview: Vec<_> = values.iter().take(3).collect();
                                log::warn!(
                                    "🔍 BETA PROCESSING START | Raw values: {:?} | units: {:?} | scale: {:?}",
                                    preview, current_units, current_scale
                                );
                            }

                            // Apply per-element scale if present
                            if let Some(scale) = current_scale {
                                for v in &mut values {
                                    *v *= scale;
                                }
                            }

                            // Handle robust unit conversion with tenths-dB detection
                            convert_units_auto(
                                "betaNought",
                                &mut values,
                                current_units.as_ref().map(|s| s.as_str()),
                            );

                            // Basic regression detection for critical issues
                            if values.len() > 1 {
                                let min_val = values.iter().fold(f64::INFINITY, |a, &b| a.min(b));
                                let max_val =
                                    values.iter().fold(f64::NEG_INFINITY, |a, &b| a.max(b));
                                let ratio = if min_val > 0.0 {
                                    max_val / min_val
                                } else {
                                    0.0
                                };

                                // Critical regression detection: warn if all values identical (β⁰ = 237 problem)
                                if ratio < 1.001 && values.len() > 50 {
                                    let warn_index =
                                        ALL_BETA_WARN_COUNT.fetch_add(1, AtomicOrdering::Relaxed);
                                    if warn_index < ALL_BETA_WARN_LIMIT {
                                        log::warn!(
                                            "All beta values identical - potential XML parsing issue"
                                        );
                                    } else if warn_index == ALL_BETA_WARN_LIMIT {
                                        log::warn!(
                                            "All beta values identical - additional occurrences suppressed (>{})",
                                            ALL_BETA_WARN_LIMIT
                                        );
                                    }
                                    // Flatness flag will be set later via is_flat; no double signalling here.
                                }
                            }

                            // Debug: Show values and statistics after unit conversion
                            if !values.is_empty() {
                                let preview: Vec<_> = values.iter().take(3).collect();
                                let mut sorted = values.clone();
                                sorted.sort_by(|a, b| a.total_cmp(b));
                                let median = sorted[sorted.len() / 2];
                                let vmax = *sorted.last().unwrap_or(&median);
                                let vmin = *sorted.first().unwrap_or(&median);
                                log::warn!(
                                    "🔍 BETA AFTER UNIT CONVERSION: {:?} | median={:.3e} | range=[{:.3e},{:.3e}] | n={}",
                                    preview, median, vmin, vmax, values.len()
                                );
                            }

                            vector.beta_nought = values.into_iter().map(|v| v as f32).collect();

                            // Beta nought values parsed
                        }
                        in_beta = false;
                    }
                    "gamma" if in_gamma => {
                        if let Some(ref mut vector) = current_vector {
                            let mut values = parse_float_list(&accumulated_text);

                            // Debug: Show original values before any processing
                            if !values.is_empty() {
                                let first_few: Vec<_> = values.iter().take(3).collect();
                                log::warn!("🔍 GAMMA PROCESSING START | Raw values: {:?} | units: {:?} | scale: {:?}", first_few, current_units, current_scale);
                            }

                            // Apply per-element scale if present BEFORE unit conversion (consistent with sigma/beta)
                            if let Some(scale) = current_scale {
                                log::warn!(
                                    "🔍 Gamma SCALE APPLICATION: applying {} to {} values",
                                    scale,
                                    values.len()
                                );
                                for v in &mut values {
                                    *v *= scale;
                                }
                                // Debug: Show values after scaling
                                if !values.is_empty() {
                                    let first_few: Vec<_> = values.iter().take(3).collect();
                                    log::warn!("🔍 Gamma after scaling: {:?}", first_few);
                                }
                            } else {
                                log::info!("🔍 Gamma: NO scale factor to apply");
                            }

                            // Convert units with robust auto-detection AFTER scaling
                            log::info!(
                                "🔍 Gamma: calling convert_units_auto with units: {:?}",
                                current_units
                            );
                            convert_units_auto(
                                "gamma",
                                &mut values,
                                current_units.as_ref().map(|s| s.as_str()),
                            );

                            // Show final values after unit conversion
                            if !values.is_empty() {
                                let first_few: Vec<_> = values.iter().take(3).collect();
                                log::warn!("🔍 Gamma AFTER UNIT CONVERSION: {:?}", first_few);

                                // Calculate statistics
                                let median = if values.len() > 0 {
                                    let mut sorted = values.clone();
                                    sorted.sort_by(|a, b| a.total_cmp(b));
                                    sorted[values.len() / 2]
                                } else {
                                    0.0
                                };

                                let max_value = values
                                    .iter()
                                    .fold(f64::NEG_INFINITY, |a, &b| a.max(b as f64))
                                    as f64;
                                log::warn!(
                                    "🔍 Gamma STATISTICS: count={}, median={:.3e}, max={:.3e}",
                                    values.len(),
                                    median,
                                    max_value
                                );
                            }

                            // Early sanity check to prevent astronomical outputs
                            sanity_check_lut("gamma", &values)?;

                            vector.gamma = values.into_iter().map(|v| v as f32).collect();

                            // Parsed gamma values
                        }
                        in_gamma = false;
                    }
                    "calibrationVector" => {
                        if in_calibration_vector {
                            if let Some(mut vector) = current_vector.take() {
                                // Validate vector has required data
                                let has_any = !vector.sigma_nought.is_empty()
                                    || !vector.beta_nought.is_empty()
                                    || !vector.gamma.is_empty();

                                if !vector.pixels.is_empty() && has_any {
                                    // Run vector sanity checks
                                    vector_sanity(&vector)?;

                                    // Set flatness flags
                                    vector.beta_flat = is_flat(&vector.beta_nought);
                                    vector.sigma_flat = is_flat(&vector.sigma_nought);
                                    vector.gamma_flat = is_flat(&vector.gamma);

                                    if vector.beta_flat {
                                        log::warn!(
                                            "β⁰ flat at line {} (attempting recovery from σ⁰)",
                                            vector.line
                                        );
                                    }
                                    if vector.sigma_flat {
                                        log::warn!(
                                            "σ⁰ flat at line {} (attempting recovery from β⁰/γ⁰)",
                                            vector.line
                                        );
                                    }
                                    if vector.gamma_flat {
                                        log::warn!(
                                            "γ⁰ flat at line {} (attempting recovery from σ⁰)",
                                            vector.line
                                        );
                                    }

                                    // CRITICAL FIX: Recover corrupted coefficients using swath-specific incidence angles
                                    if vector.beta_flat || vector.sigma_flat || vector.gamma_flat {
                                        // Use swath-specific incidence angles for accurate recovery
                                        let min_inc = inc_angles.min_deg;
                                        let max_inc = inc_angles.max_deg;

                                        // Estimate swath width from pixel array
                                        let swath_width =
                                            vector.pixels.last().copied().unwrap_or(25000);

                                        log::info!(
                                            "🔧 Recovering coefficients at line {} using swath-specific angles: {:.1}°-{:.1}°",
                                            vector.line, min_inc, max_inc
                                        );

                                        let recovered = vector.recover_corrupted_coefficients(
                                            min_inc,
                                            max_inc,
                                            swath_width,
                                        );

                                        if !recovered {
                                            log::error!(
                                                "❌ Failed to recover any calibration coefficients at line {}",
                                                vector.line
                                            );
                                        }
                                    }

                                    // Completed vector (possibly recovered)
                                    vectors.push(vector);
                                } else {
                                    // Skipping incomplete vector
                                }
                            }
                            in_calibration_vector = false;
                        }
                    }
                    _ => {}
                }
            }
            Ok(Event::Eof) => break,
            Err(e) => {
                return Err(SarError::Processing(format!("XML parsing error: {}", e)));
            }
            _ => {}
        }
        buf.clear();
    }

    if debug {
        // Parsed calibration vectors
    }

    if vectors.is_empty() {
        return Err(SarError::Processing(
            "No valid calibration vectors found in XML".to_string(),
        ));
    }

    Ok(vectors)
}

/// Derive Beta0 from Sigma0 using radar equation: β⁰ = σ⁰ / sin(θ)
///
/// # Physical Basis
/// The radar equation relates different calibration types through incidence angle:
/// - Beta0 (β⁰): Radar brightness coefficient (perpendicular to sensor LOS)
/// - Sigma0 (σ⁰): Radar backscatter cross-section (normalized to ground area)
/// - Relationship: β⁰ = σ⁰ / sin(θ_inc)
pub fn derive_beta_from_sigma(
    sigma_nought: &[f32],
    pixels: &[usize],
    min_incidence_deg: f32,
    max_incidence_deg: f32,
    swath_width: usize,
) -> Vec<f32> {
    let min_inc_rad = min_incidence_deg.to_radians();
    let max_inc_rad = max_incidence_deg.to_radians();

    sigma_nought
        .iter()
        .zip(pixels.iter())
        .map(|(&sigma, &pixel)| {
            let range_fraction = (pixel as f32) / (swath_width as f32).max(1.0);
            let range_fraction = range_fraction.clamp(0.0, 1.0);
            let incidence_rad = min_inc_rad + range_fraction * (max_inc_rad - min_inc_rad);
            let sin_inc = incidence_rad.sin();

            if sin_inc > 0.01 {
                sigma / sin_inc
            } else {
                sigma
            }
        })
        .collect()
}

/// Derive Gamma0 from Sigma0 using radar equation: γ⁰ = σ⁰ / cos(θ)
pub fn derive_gamma_from_sigma(
    sigma_nought: &[f32],
    pixels: &[usize],
    min_incidence_deg: f32,
    max_incidence_deg: f32,
    swath_width: usize,
) -> Vec<f32> {
    let min_inc_rad = min_incidence_deg.to_radians();
    let max_inc_rad = max_incidence_deg.to_radians();

    sigma_nought
        .iter()
        .zip(pixels.iter())
        .map(|(&sigma, &pixel)| {
            let range_fraction = (pixel as f32) / (swath_width as f32).max(1.0);
            let range_fraction = range_fraction.clamp(0.0, 1.0);
            let incidence_rad = min_inc_rad + range_fraction * (max_inc_rad - min_inc_rad);
            let cos_inc = incidence_rad.cos();

            // Physical relation: γ⁰ = σ⁰ / cos(θ).
            // For very shallow or grazing angles (cosθ → 0), this can explode numerically.
            // Instead of applying an arbitrary gain (e.g. ×10), fall back to σ⁰ when the
            // geometry is unreliable, and log this at the call site if needed.
            if cos_inc.abs() > 0.01 {
                sigma / cos_inc
            } else {
                sigma
            }
        })
        .collect()
}

/// Derive Sigma0 from Beta0 (reverse transformation): σ⁰ = β⁰ × sin(θ)
pub fn derive_sigma_from_beta(
    beta_nought: &[f32],
    pixels: &[usize],
    min_incidence_deg: f32,
    max_incidence_deg: f32,
    swath_width: usize,
) -> Vec<f32> {
    let min_inc_rad = min_incidence_deg.to_radians();
    let max_inc_rad = max_incidence_deg.to_radians();

    beta_nought
        .iter()
        .zip(pixels.iter())
        .map(|(&beta, &pixel)| {
            let range_fraction = (pixel as f32) / (swath_width as f32).max(1.0);
            let range_fraction = range_fraction.clamp(0.0, 1.0);
            let incidence_rad = min_inc_rad + range_fraction * (max_inc_rad - min_inc_rad);
            let sin_inc = incidence_rad.sin();
            beta * sin_inc
        })
        .collect()
}

/// Parse thermal noise data from Sentinel-1 noise XML file
/// References: ESA Sentinel-1 Level 1 Detailed Algorithm Definition, Section 2.2.3
/// Supports both noiseRangeVectorList (all IPF) and noiseAzimuthVectorList (IPF 3.x+)
pub fn parse_noise_from_xml(xml_content: &str) -> SarResult<NoiseCoefficients> {
    let mut reader = Reader::from_str(xml_content);
    reader.trim_text(true);

    let mut noise = NoiseCoefficients::new();
    let mut buf = Vec::new();

    // State tracking for range noise parsing
    let mut in_ads_header = false;
    let mut in_noise_vector = false;
    let mut _in_noise_vector_list = false;
    let mut current_vector: Option<NoiseVector> = None;

    // State tracking for azimuth noise parsing (IPF 3.x+)
    let mut in_azimuth_noise_vector = false;
    let mut _in_azimuth_noise_list = false;
    let mut current_azimuth_vector: Option<NoiseAzimuthVector> = None;

    #[allow(unused_assignments)]
    let mut current_tag = String::new();
    #[allow(unused_assignments)]
    let mut text_content = String::new();

    loop {
        match reader.read_event_into(&mut buf) {
            Ok(Event::Start(ref e)) => {
                current_tag = String::from_utf8_lossy(e.name().as_ref()).to_string();
                text_content.clear();

                match current_tag.as_str() {
                    "adsHeader" => in_ads_header = true,
                    "noiseRangeVectorList" => _in_noise_vector_list = true,
                    "noiseRangeVector" => {
                        in_noise_vector = true;
                        current_vector = Some(NoiseVector {
                            azimuth_time: String::new(),
                            azimuth_time_seconds: 0.0,
                            line: 0.0,
                            range_pixels: Vec::new(),
                            noise_range_lut: Vec::new(),
                        });
                    }
                    // IPF 3.x+ azimuth noise vectors
                    "noiseAzimuthVectorList" => _in_azimuth_noise_list = true,
                    "noiseAzimuthVector" => {
                        in_azimuth_noise_vector = true;
                        current_azimuth_vector = Some(NoiseAzimuthVector {
                            swath: String::new(),
                            first_azimuth_line: 0,
                            first_range_sample: 0,
                            last_azimuth_line: 0,
                            last_range_sample: 0,
                            lines: Vec::new(),
                            noise_azimuth_lut: Vec::new(),
                        });
                    }
                    _ => {}
                }
            }
            Ok(Event::Text(e)) => {
                text_content.push_str(
                    &e.unescape()
                        .map_err(|e| {
                            SarError::Processing(format!(
                                "Failed to unescape XML text in noise parsing: {}",
                                e
                            ))
                        })?
                        .to_string(),
                );
            }
            Ok(Event::End(ref e)) => {
                let tag_name = String::from_utf8_lossy(e.name().as_ref()).to_string();

                if in_ads_header {
                    match tag_name.as_str() {
                        "polarisation" => noise.polarization = text_content.clone(),
                        "swath" => noise.swath = text_content.clone(),
                        "adsHeader" => in_ads_header = false,
                        _ => {}
                    }
                } else if in_noise_vector && current_vector.is_some() {
                    let vector = current_vector
                        .as_mut()
                        .expect("current_vector confirmed as Some above");

                    match tag_name.as_str() {
                        "azimuthTime" => {
                            vector.azimuth_time = text_content.clone();
                            // Parse azimuth time to seconds for numerical processing
                            match parse_azimuth_time_to_seconds(&text_content) {
                                Ok(seconds) => vector.azimuth_time_seconds = seconds,
                                Err(e) => {
                                    log::warn!(
                                        "Failed to parse azimuth time '{}': {}",
                                        text_content,
                                        e
                                    );
                                    vector.azimuth_time_seconds = 0.0;
                                }
                            }
                        }
                        "line" => {
                            // Parse line number as f64 to preserve negative values and fractional precision
                            if let Ok(parsed_line) = text_content.trim().parse::<f64>() {
                                vector.line = parsed_line;
                            } else {
                                log::warn!("Failed to parse line value as f64: '{}'", text_content);
                            }
                        }
                        "pixel" => {
                            // Parse space-separated pixel indices as f64 to preserve sub-pixel precision
                            vector.range_pixels = text_content
                                .split_whitespace()
                                .filter_map(|s| s.parse::<f64>().ok())
                                .collect();
                        }
                        "noiseRangeLut" => {
                            // Parse space-separated noise values (scientific notation)
                            vector.noise_range_lut = text_content
                                .split_whitespace()
                                .filter_map(|s| s.parse::<f32>().ok())
                                .collect();
                        }
                        "noiseRangeVector" => {
                            // Vector complete, add to coefficients
                            if let Some(vector) = current_vector.take() {
                                if !vector.range_pixels.is_empty()
                                    && !vector.noise_range_lut.is_empty()
                                {
                                    if vector.range_pixels.len() == vector.noise_range_lut.len() {
                                        // Noise vector parsed

                                        // Validate monotonicity of range axis
                                        let is_monotonic =
                                            vector.range_pixels.windows(2).all(|w| w[0] <= w[1]);
                                        if !is_monotonic {
                                            log::warn!("Non-monotonic range pixels in noise vector at line {}", vector.line);
                                        }

                                        noise.vectors.push(vector);
                                    } else {
                                        log::warn!(
                                            "Noise vector pixel/LUT length mismatch: {} vs {}",
                                            vector.range_pixels.len(),
                                            vector.noise_range_lut.len()
                                        );
                                    }
                                } else {
                                    log::warn!("Empty noise vector data for line {}", vector.line);
                                }
                            }
                            in_noise_vector = false;
                        }
                        _ => {}
                    }
                } else if in_azimuth_noise_vector && current_azimuth_vector.is_some() {
                    // Handle azimuth noise vector elements (IPF 3.x+)
                    let vector = current_azimuth_vector
                        .as_mut()
                        .expect("current_azimuth_vector confirmed as Some above");

                    match tag_name.as_str() {
                        "swath" => vector.swath = text_content.clone(),
                        "firstAzimuthLine" => {
                            if let Ok(v) = text_content.trim().parse::<usize>() {
                                vector.first_azimuth_line = v;
                            }
                        }
                        "firstRangeSample" => {
                            if let Ok(v) = text_content.trim().parse::<usize>() {
                                vector.first_range_sample = v;
                            }
                        }
                        "lastAzimuthLine" => {
                            if let Ok(v) = text_content.trim().parse::<usize>() {
                                vector.last_azimuth_line = v;
                            }
                        }
                        "lastRangeSample" => {
                            if let Ok(v) = text_content.trim().parse::<usize>() {
                                vector.last_range_sample = v;
                            }
                        }
                        "line" => {
                            // Parse space-separated line indices
                            vector.lines = text_content
                                .split_whitespace()
                                .filter_map(|s| s.parse::<usize>().ok())
                                .collect();
                        }
                        "noiseAzimuthLut" => {
                            // Parse space-separated noise values
                            vector.noise_azimuth_lut = text_content
                                .split_whitespace()
                                .filter_map(|s| s.parse::<f32>().ok())
                                .collect();
                        }
                        "noiseAzimuthVector" => {
                            // Vector complete, add to coefficients
                            if let Some(vector) = current_azimuth_vector.take() {
                                if !vector.lines.is_empty() && !vector.noise_azimuth_lut.is_empty()
                                {
                                    if vector.lines.len() == vector.noise_azimuth_lut.len() {
                                        log::debug!(
                                            "Parsed azimuth noise vector for {}: {} values, lines {}-{}",
                                            vector.swath,
                                            vector.noise_azimuth_lut.len(),
                                            vector.first_azimuth_line,
                                            vector.last_azimuth_line
                                        );
                                        noise.azimuth_vectors.push(vector);
                                    } else {
                                        log::warn!(
                                            "Azimuth noise vector line/LUT length mismatch: {} vs {}",
                                            vector.lines.len(),
                                            vector.noise_azimuth_lut.len()
                                        );
                                    }
                                }
                            }
                            in_azimuth_noise_vector = false;
                        }
                        _ => {}
                    }
                }
            }
            Ok(Event::Eof) => break,
            Err(e) => {
                log::error!("XML parsing error: {}", e);
                return Err(SarError::Processing(format!("XML parsing error: {}", e)));
            }
            _ => {}
        }

        buf.clear();
    }

    if noise.vectors.is_empty() {
        return Err(SarError::Processing(
            "No noise vectors found in XML".to_string(),
        ));
    }

    // Log parsed noise data summary
    if noise.azimuth_vectors.is_empty() {
        log::info!(
            "Successfully parsed {} range noise vectors from XML (no azimuth noise - pre-IPF 3.x product)",
            noise.vectors.len()
        );
    } else {
        log::info!(
            "Successfully parsed {} range noise vectors and {} azimuth noise vectors from XML (IPF 3.x+)",
            noise.vectors.len(),
            noise.azimuth_vectors.len()
        );
    }

    // Validate monotonicity and consistency after parsing is complete
    noise.validate_vectors()?;

    Ok(noise)
}

/// Parse antenna pattern from Sentinel-1 annotation XML
///
/// The antenna pattern correction removes azimuth scalloping caused by the rotating
/// antenna beam in TOPS mode (typically 0.3-0.5 dB amplitude).
///
/// Reference: ESA-EOPG-CSCOP-TN-0010 "Sentinel-1 TOPS Radiometric Calibration"
pub fn parse_antenna_pattern_from_xml(
    xml_content: &str,
    range_sampling_rate_hz: Option<f64>,
) -> SarResult<Vec<AntennaPatternVector>> {
    log::info!(
        "📡 Parsing antenna pattern from annotation XML ({} bytes)",
        xml_content.len()
    );

    let mut reader = Reader::from_str(xml_content);
    reader.trim_text(true);
    let mut buf = Vec::new();

    let mut vectors = Vec::new();
    let mut current_vector: Option<AntennaPatternVector> = None;
    let mut current_tag = String::new();
    let mut text_content = String::new();
    let mut in_antenna_pattern_list = false;

    // Track expected sample counts from count="" attributes so we can enforce alignment
    let mut expected_pixel_count: Option<usize> = None;
    let mut expected_slant_count: Option<usize> = None;
    let mut expected_pattern_count: Option<usize> = None;
    let mut pattern_is_complex: Option<bool> = None;

    fn resample_values(values: &[f64], target_len: usize) -> Vec<f64> {
        if values.is_empty() || target_len == 0 {
            return Vec::new();
        }
        if values.len() == target_len {
            return values.to_vec();
        }

        if target_len == 1 {
            return vec![values[0]];
        }

        let scale = (values.len() - 1) as f64 / (target_len - 1) as f64;
        let mut out = Vec::with_capacity(target_len);
        for i in 0..target_len {
            let pos = i as f64 * scale;
            let idx0 = pos.floor() as usize;
            let idx1 = pos.ceil() as usize;
            if idx0 == idx1 {
                out.push(values[idx0]);
            } else {
                let w = pos - idx0 as f64;
                let v0 = values[idx0];
                let v1 = values[idx1.min(values.len() - 1)];
                out.push(v0 + (v1 - v0) * w);
            }
        }
        out
    }

    fn normalize_pattern_samples(
        mut raw_values: Vec<f64>,
        target_len: Option<usize>,
        pattern_is_complex: Option<bool>,
        line: i32,
    ) -> Vec<f64> {
        if raw_values.is_empty() {
            return Vec::new();
        }

        let mut treat_as_complex = pattern_is_complex.unwrap_or(false);
        let mut target = target_len.unwrap_or(0);

        log::debug!(
            "Normalizing antenna pattern at line {}: raw_len={}, target_hint={:?}, treat_complex={}",
            line,
            raw_values.len(),
            target_len,
            treat_as_complex
        );

        if !treat_as_complex {
            if let Some(t) = target_len {
                if raw_values.len() == t * 2 {
                    treat_as_complex = true;
                    target = t;
                }
            }
        }

        if treat_as_complex {
            if target == 0 {
                target = raw_values.len() / 2;
            }
            let available_pairs = raw_values.len() / 2;
            if available_pairs == 0 {
                log::warn!(
                    "⚠️  Antenna pattern at line {} has insufficient I/Q samples ({} floats)",
                    line,
                    raw_values.len()
                );
                return Vec::new();
            }

            if available_pairs < target {
                log::warn!(
                    "⚠️  Antenna pattern at line {} expected {} I/Q pairs but only found {}; using available pairs",
                    line,
                    target,
                    available_pairs
                );
                target = available_pairs;
            }

            if available_pairs > target {
                log::debug!(
                    "Truncating antenna pattern I/Q pairs at line {} from {} to {} pairs",
                    line,
                    available_pairs,
                    target
                );
                raw_values.truncate(target * 2);
            }

            let mut magnitudes = Vec::with_capacity(target);
            for chunk in raw_values.chunks_exact(2) {
                let i = chunk[0];
                let q = chunk[1];
                magnitudes.push((i * i + q * q).sqrt());
            }

            log::debug!(
                "Converted antenna pattern I/Q pairs at line {} into {} magnitudes (target_hint={:?})",
                line,
                magnitudes.len(),
                target_len
            );

            if let Some(t) = target_len {
                if t > 0 && magnitudes.len() != t {
                    log::debug!(
                        "Resampling antenna pattern magnitudes at line {} from {} to {} samples",
                        line,
                        magnitudes.len(),
                        t
                    );
                    return resample_values(&magnitudes, t);
                }
            }

            magnitudes
        } else {
            if let Some(t) = target_len {
                if t > 0 && raw_values.len() != t {
                    log::debug!(
                        "Resampling antenna pattern amplitudes at line {} from {} to {} samples",
                        line,
                        raw_values.len(),
                        t
                    );
                    return resample_values(&raw_values, t);
                }
            }
            log::debug!(
                "Keeping antenna pattern amplitudes at line {} with {} samples (target_hint={:?})",
                line,
                raw_values.len(),
                target_len
            );
            raw_values
        }
    }

    loop {
        match reader.read_event_into(&mut buf) {
            Ok(Event::Start(ref e)) => {
                current_tag = String::from_utf8_lossy(e.name().as_ref()).to_string();
                text_content.clear();

                if current_tag == "antennaPatternList" {
                    in_antenna_pattern_list = true;
                    log::debug!("Found antennaPatternList");
                } else if current_tag == "antennaPattern" && in_antenna_pattern_list {
                    current_vector = Some(AntennaPatternVector {
                        azimuth_time: String::new(),
                        line: 0,
                        pixels: Vec::new(),
                        values: Vec::new(),
                    });
                }

                match current_tag.as_str() {
                    "pixel" => {
                        expected_pixel_count = None;
                        for attr in e.attributes().flatten() {
                            if attr.key.as_ref() == b"count" {
                                if let Ok(val) = std::str::from_utf8(&attr.value) {
                                    expected_pixel_count = val.trim().parse::<usize>().ok();
                                }
                            }
                        }
                    }
                    "slantRangeTime" => {
                        expected_slant_count = None;
                        for attr in e.attributes().flatten() {
                            if attr.key.as_ref() == b"count" {
                                if let Ok(val) = std::str::from_utf8(&attr.value) {
                                    expected_slant_count = val.trim().parse::<usize>().ok();
                                }
                            }
                        }
                    }
                    "elevationPattern" => {
                        pattern_is_complex = Some(true);
                        expected_pattern_count = None;
                        for attr in e.attributes().flatten() {
                            if attr.key.as_ref() == b"count" {
                                if let Ok(val) = std::str::from_utf8(&attr.value) {
                                    expected_pattern_count = val.trim().parse::<usize>().ok();
                                }
                            }
                        }
                    }
                    "pattern" => {
                        pattern_is_complex = Some(false);
                        expected_pattern_count = None;
                        for attr in e.attributes().flatten() {
                            if attr.key.as_ref() == b"count" {
                                if let Ok(val) = std::str::from_utf8(&attr.value) {
                                    expected_pattern_count = val.trim().parse::<usize>().ok();
                                }
                            }
                        }
                    }
                    _ => {}
                }
            }
            Ok(Event::Text(e)) => {
                if !current_tag.is_empty() {
                    text_content.push_str(&e.unescape().unwrap_or_default());
                }
            }
            Ok(Event::End(ref e)) => {
                let tag_name = String::from_utf8_lossy(e.name().as_ref()).to_string();

                if tag_name == "antennaPatternList" {
                    in_antenna_pattern_list = false;
                } else if let Some(ref mut vector) = current_vector {
                    match tag_name.as_str() {
                        "azimuthTime" => {
                            vector.azimuth_time = text_content.trim().to_string();
                        }
                        "line" => {
                            if let Ok(line_val) = text_content.trim().parse::<i32>() {
                                vector.line = line_val;
                            }
                        }
                        "pixel" => {
                            let mut pixels = parse_usize_list(&text_content);
                            if let Some(expected) = expected_pixel_count {
                                if pixels.len() > expected {
                                    log::debug!(
                                        "Truncating antenna pattern pixel list at line {} from {} to {} entries",
                                        vector.line,
                                        pixels.len(),
                                        expected
                                    );
                                    pixels.truncate(expected);
                                } else if pixels.len() < expected {
                                    log::warn!(
                                        "Warning: antenna pattern pixel list at line {} has {} entries but expected {}",
                                        vector.line,
                                        pixels.len(),
                                        expected
                                    );
                                }
                            }
                            vector.pixels = pixels;
                            expected_pixel_count = None;
                        }
                        "slantRangeTime" => {
                            let mut times: Vec<f64> = text_content
                                .split_whitespace()
                                .filter_map(|s| s.parse::<f64>().ok())
                                .collect();
                            if let Some(expected) = expected_slant_count {
                                if times.len() > expected {
                                    log::debug!(
                                        "Truncating slant range times at line {} from {} to {} entries",
                                        vector.line,
                                        times.len(),
                                        expected
                                    );
                                    times.truncate(expected);
                                } else if times.len() < expected {
                                    log::warn!(
                                        "Warning: slant range time list at line {} has {} entries but expected {}",
                                        vector.line,
                                        times.len(),
                                        expected
                                    );
                                }
                            }

                            if !times.is_empty() {
                                let pixels: Vec<usize> = if let Some(rate) = range_sampling_rate_hz
                                {
                                    if rate > 0.0 {
                                        let dt = 1.0 / rate;
                                        let t0 = times[0];
                                        times
                                            .iter()
                                            .map(|t| {
                                                let idx = if dt > 0.0 {
                                                    ((t - t0) / dt).round()
                                                } else {
                                                    0.0
                                                };
                                                idx.max(0.0) as usize
                                            })
                                            .collect()
                                    } else {
                                        (0..times.len()).collect()
                                    }
                                } else {
                                    (0..times.len()).collect()
                                };
                                if vector.pixels.is_empty() || vector.pixels.len() != pixels.len() {
                                    if !vector.pixels.is_empty()
                                        && vector.pixels.len() != pixels.len()
                                    {
                                        log::debug!(
                                            "Overriding existing pixel list at line {} ({} vs {} samples)",
                                            vector.line,
                                            vector.pixels.len(),
                                            pixels.len()
                                        );
                                    }
                                    vector.pixels = pixels;
                                }
                            }
                            expected_slant_count = None;
                        }
                        "elevationPattern" => {
                            let raw_values = parse_float_list(&text_content);
                            vector.values = normalize_pattern_samples(
                                raw_values,
                                expected_pattern_count.or_else(|| {
                                    if !vector.pixels.is_empty() {
                                        Some(vector.pixels.len())
                                    } else {
                                        None
                                    }
                                }),
                                pattern_is_complex,
                                vector.line,
                            );
                            expected_pattern_count = None;
                            pattern_is_complex = None;
                        }
                        "pattern" => {
                            let raw_values = parse_float_list(&text_content);
                            vector.values = normalize_pattern_samples(
                                raw_values,
                                expected_pattern_count.or_else(|| {
                                    if !vector.pixels.is_empty() {
                                        Some(vector.pixels.len())
                                    } else {
                                        None
                                    }
                                }),
                                pattern_is_complex,
                                vector.line,
                            );
                            expected_pattern_count = None;
                            pattern_is_complex = None;
                        }
                        "antennaPattern" => {
                            // Vector complete, validate and add
                            if let Some(vector) = current_vector.take() {
                                if !vector.pixels.is_empty() && !vector.values.is_empty() {
                                    let mut normalized = vector;

                                    // Handle common mismatch scenarios before giving up
                                    if normalized.pixels.len() != normalized.values.len() {
                                        if normalized.values.len() == normalized.pixels.len() * 2 {
                                            // Safety net: pattern values still in I/Q pairs (should have been handled earlier)
                                            let mut magnitudes =
                                                Vec::with_capacity(normalized.values.len() / 2);
                                            for chunk in normalized.values.chunks_exact(2) {
                                                let i = chunk[0];
                                                let q = chunk[1];
                                                magnitudes.push((i * i + q * q).sqrt());
                                            }
                                            log::warn!(
                                                "Antenna pattern vector at line {} had I/Q pairs; converted {} pairs to magnitudes",
                                                normalized.line,
                                                magnitudes.len()
                                            );
                                            normalized.values = magnitudes;
                                        }

                                        if normalized.pixels.len() != normalized.values.len() {
                                            let new_len = normalized
                                                .pixels
                                                .len()
                                                .min(normalized.values.len());
                                            if new_len == 0 {
                                                log::warn!(
                                                    "Discarding antenna pattern vector at line {} due to empty pixel/value arrays",
                                                    normalized.line
                                                );
                                                continue;
                                            }
                                            log::warn!(
                                                "Antenna pattern vector length mismatch: {} pixels vs {} values; truncating to {}",
                                                normalized.pixels.len(),
                                                normalized.values.len(),
                                                new_len
                                            );
                                            normalized.pixels.truncate(new_len);
                                            normalized.values.truncate(new_len);
                                        }
                                    }

                                    if normalized.pixels.len() == normalized.values.len() {
                                        let mut mags: Vec<f64> =
                                            normalized.values.iter().map(|v| v.abs()).collect();
                                        mags.retain(|m| m.is_finite() && *m > 0.0);
                                        let scale = if mags.is_empty() {
                                            1.0
                                        } else {
                                            let mid = mags.len() / 2;
                                            mags.sort_by(|a, b| {
                                                a.partial_cmp(b)
                                                    .unwrap_or(std::cmp::Ordering::Equal)
                                            });
                                            mags[mid]
                                        };
                                        let scale = if scale.is_finite() && scale > 0.0 {
                                            scale
                                        } else {
                                            1.0
                                        };
                                        normalized
                                            .values
                                            .iter_mut()
                                            .for_each(|v| *v = (v.abs() / scale).max(1e-6));
                                        log::debug!(
                                            "Parsed antenna pattern vector: line={}, {} points",
                                            normalized.line,
                                            normalized.pixels.len()
                                        );
                                        vectors.push(normalized);
                                    } else {
                                        log::warn!(
                                            "Skipping antenna pattern vector at line {} due to irreconcilable length mismatch ({} pixels vs {} values)",
                                            normalized.line,
                                            normalized.pixels.len(),
                                            normalized.values.len()
                                        );
                                    }
                                }
                            }
                        }
                        _ => {}
                    }
                }

                current_tag.clear();
                text_content.clear();
            }
            Ok(Event::Eof) => break,
            Err(e) => {
                return Err(SarError::Processing(format!(
                    "Error parsing antenna pattern XML at position {}: {:?}",
                    reader.buffer_position(),
                    e
                )));
            }
            _ => {}
        }

        buf.clear();
    }

    if vectors.is_empty() {
        log::warn!("⚠️  No antenna pattern vectors found in XML - scalloping correction will not be applied");
        return Ok(Vec::new()); // Return empty vector, not an error (some products may not have this)
    }

    log::info!(
        "✅ Successfully parsed {} antenna pattern vectors",
        vectors.len()
    );

    // Validate pattern values
    let mut all_values = Vec::new();
    for vector in &vectors {
        all_values.extend(vector.values.iter().copied());
    }

    if !all_values.is_empty() {
        let min_val = all_values.iter().fold(f64::INFINITY, |a, &b| a.min(b));
        let max_val = all_values.iter().fold(f64::NEG_INFINITY, |a, &b| a.max(b));
        let mean_val = all_values.iter().sum::<f64>() / all_values.len() as f64;

        log::info!(
            "Antenna pattern statistics: min={:.3}, max={:.3}, mean={:.3}",
            min_val,
            max_val,
            mean_val
        );

        if max_val > 4.0 || min_val < 0.2 {
            return Err(SarError::Processing(format!(
                "Antenna pattern amplitudes out of bounds [{:.3}, {:.3}] (expected ~[0.5, 2.0]); aborting to prevent 8–20× gain spikes",
                min_val, max_val
            )));
        }

        if min_val < 0.5 || max_val > 2.0 {
            // Downgraded to debug: S1 IW calibration LUTs already include antenna pattern effects,
            // so this parsed pattern won't be used (disabled in precompute_antenna_pattern_lut_separable).
            // Showing a warning for a feature we don't use is confusing.
            log::debug!(
                "Antenna pattern range [{:.3}, {:.3}] outside safe band [0.5, 2.0]; would clamp (but S1 IW uses LUT-included patterns)",
                min_val,
                max_val
            );
            for vec in vectors.iter_mut() {
                for v in vec.values.iter_mut() {
                    if *v < 0.5 {
                        *v = 0.5;
                    } else if *v > 2.0 {
                        *v = 2.0;
                    }
                }
            }
        }
    }

    Ok(vectors)
}
