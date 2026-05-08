#![allow(dead_code)]
use chrono::{NaiveDateTime, TimeZone, Utc};
use regex::Regex;
use std::sync::OnceLock;

use crate::types::{OrbitData, SarError, SarResult, StateVector};

// OPTIMIZATION: Lazy static regex to avoid recompilation per call
static OSV_BLOCK_RE: OnceLock<Regex> = OnceLock::new();
static TAI_UTCID_RE: OnceLock<Regex> = OnceLock::new();
static UTC_REGEX: OnceLock<Regex> = OnceLock::new();

fn get_osv_block_regex() -> &'static Regex {
    OSV_BLOCK_RE.get_or_init(|| Regex::new(r"(?s)<OSV>(.*?)</OSV>").unwrap())
}

/// Parse EOF content, preferring block-based OSV extraction with fallback to legacy UTC lines
pub(crate) fn parse_eof_content_enhanced(content: &str) -> SarResult<OrbitData> {
    let mut state_vectors = Vec::new();
    let mut reference_time = None;

    log::info!("Parsing EOF orbit file ({} bytes)", content.len());

    // Use cached regex instead of compiling each time
    let osv_block_re = get_osv_block_regex();

    let mut found_blocks = 0;
    for captures in osv_block_re.captures_iter(content) {
        if let Some(block_match) = captures.get(1) {
            found_blocks += 1;
            let block_content = block_match.as_str();

            match parse_osv_block(block_content) {
                Ok(Some(sv)) => {
                    if reference_time.is_none() {
                        reference_time = Some(sv.time);
                    }
                    log::debug!(
                        "Parsed OSV block {}: position/velocity at {}",
                        found_blocks,
                        sv.time.format("%Y-%m-%d %H:%M:%S")
                    );
                    state_vectors.push(sv);
                }
                Ok(None) => {
                    log::warn!("Skipped incomplete OSV block {}", found_blocks);
                }
                Err(e) => {
                    log::warn!("Failed to parse OSV block {}: {}", found_blocks, e);
                }
            }
        }
    }

    if found_blocks > 0 {
        log::info!(
            "Block-based parsing found {} OSV blocks, {} valid state vectors",
            found_blocks,
            state_vectors.len()
        );
    }

    // Fallback: legacy text-based parsing only when no OSV blocks exist
    if state_vectors.is_empty() {
        log::info!("No OSV blocks found; falling back to UTC= line parsing only");

        let mut line_count = 0;

        for line in content.lines() {
            line_count += 1;
            let line = line.trim();

            if line.starts_with('#') || line.is_empty() {
                continue;
            }

            if line.starts_with("UTC=") {
                match parse_state_vector_line_enhanced(line) {
                    Ok(Some(sv)) => {
                        if reference_time.is_none() {
                            reference_time = Some(sv.time);
                        }
                        state_vectors.push(sv);
                    }
                    Ok(None) => continue,
                    Err(e) => {
                        log::warn!("Failed to parse state vector at line {}: {}", line_count, e);
                        continue;
                    }
                }
            }
        }
    }

    if state_vectors.is_empty() {
        return Err(SarError::Processing(
            "No valid state vectors found in orbit file. Expected XML format with <OSV> blocks or UTC= lines.".to_string(),
        ));
    }

    state_vectors.sort_by_key(|sv| sv.time);

    let time_span = if state_vectors.len() > 1 {
        (state_vectors
            .last()
            .expect("Vector has at least 2 elements")
            .time
            - state_vectors
                .first()
                .expect("Vector has at least 2 elements")
                .time)
            .num_seconds()
    } else {
        0
    };

    log::info!(
        "Successfully parsed {} state vectors spanning {} seconds",
        state_vectors.len(),
        time_span
    );

    if !state_vectors.is_empty() {
        log::info!(
            "Time range: {} to {}",
            state_vectors
                .first()
                .expect("Vector checked for non-empty")
                .time
                .format("%Y-%m-%d %H:%M:%S"),
            state_vectors
                .last()
                .expect("Vector checked for non-empty")
                .time
                .format("%Y-%m-%d %H:%M:%S")
        );
    }

    validate_orbit_data(&state_vectors)?;

    let reference_time = reference_time.ok_or_else(|| {
        SarError::Processing(
            "Orbit file missing reference time - cannot use current time as fallback!".to_string(),
        )
    })?;

    let orbit_data = OrbitData {
        state_vectors,
        reference_time,
    };

    if let Err(validation_error) = orbit_data.validate_time_precision() {
        log::error!("❌ Orbit validation failed: {}", validation_error);
        return Err(SarError::DataProcessingError(format!(
            "Orbit time validation failed: {}",
            validation_error
        )));
    }

    log::info!(
        "✅ Orbit file loaded: reference_time={}, {} state vectors",
        reference_time,
        orbit_data.state_vectors.len()
    );

    Ok(orbit_data)
}

/// Extract XML value from block-structured content (simple string-based approach)
pub(crate) fn extract_xml_value_block(content: &str, tag: &str) -> Option<String> {
    let start_tag = format!("<{}", tag);
    let end_tag = format!("</{}>", tag);

    let start_pos = content.find(&start_tag)?;
    let tag_end = content[start_pos..].find('>')? + start_pos + 1;
    let end_pos = content[tag_end..].find(&end_tag)? + tag_end;

    let value = content[tag_end..end_pos].trim();

    if let Some(equals_pos) = value.find('=') {
        let after_equals = value[equals_pos + 1..].trim();
        if after_equals.starts_with('"') && after_equals.ends_with('"') {
            Some(after_equals[1..after_equals.len() - 1].to_string())
        } else {
            Some(after_equals.to_string())
        }
    } else {
        Some(value.to_string())
    }
}

/// Parse a complete OSV block (handles multi-line XML structure)
pub(crate) fn parse_osv_block(block: &str) -> SarResult<Option<StateVector>> {
    let utc_str = extract_xml_value_block(block, "UTC");
    let x_str = extract_xml_value_block(block, "X");
    let y_str = extract_xml_value_block(block, "Y");
    let z_str = extract_xml_value_block(block, "Z");
    let vx_str = extract_xml_value_block(block, "VX");
    let vy_str = extract_xml_value_block(block, "VY");
    let vz_str = extract_xml_value_block(block, "VZ");

    if let (Some(utc), Some(x), Some(y), Some(z), Some(vx), Some(vy), Some(vz)) =
        (utc_str, x_str, y_str, z_str, vx_str, vy_str, vz_str)
    {
        let time_str = if utc.starts_with("UTC=") {
            utc.strip_prefix("UTC=").unwrap_or(&utc)
        } else {
            &utc
        };

        let time = NaiveDateTime::parse_from_str(time_str, "%Y-%m-%dT%H:%M:%S%.f")
            .map_err(|e| SarError::Processing(format!("Invalid UTC time '{}': {}", time_str, e)))?;
        let time = Utc.from_utc_datetime(&time);

        let parse_coord = |s: &str, prefix: &str| -> SarResult<f64> {
            let value_str = if s.starts_with(prefix) {
                s.strip_prefix(prefix).unwrap_or(s)
            } else {
                s
            };
            value_str.parse::<f64>().map_err(|e| {
                SarError::Processing(format!("Invalid coordinate '{}': {}", value_str, e))
            })
        };

        let position = [
            parse_coord(&x, "X=")?,
            parse_coord(&y, "Y=")?,
            parse_coord(&z, "Z=")?,
        ];

        let velocity = [
            parse_coord(&vx, "VX=")?,
            parse_coord(&vy, "VY=")?,
            parse_coord(&vz, "VZ=")?,
        ];

        if position == [0.0, 0.0, 0.0] || velocity == [0.0, 0.0, 0.0] {
            log::warn!("Rejecting OSV with zero position or velocity vectors");
            return Ok(None);
        }

        Ok(Some(StateVector {
            time,
            position,
            velocity,
        }))
    } else {
        log::warn!("Incomplete OSV block - missing required fields");
        Ok(None)
    }
}

/// Validate orbit data quality
pub(crate) fn validate_orbit_data(state_vectors: &[StateVector]) -> SarResult<()> {
    if state_vectors.is_empty() {
        return Err(SarError::Processing(
            "No state vectors to validate".to_string(),
        ));
    }

    for sv in state_vectors {
        let velocity_magnitude =
            (sv.velocity[0].powi(2) + sv.velocity[1].powi(2) + sv.velocity[2].powi(2)).sqrt();

        if velocity_magnitude < 6000.0 || velocity_magnitude > 9000.0 {
            log::warn!(
                "Unusual orbital velocity: {:.1} m/s at {}",
                velocity_magnitude,
                sv.time.format("%Y-%m-%d %H:%M:%S")
            );
        }
    }

    for sv in state_vectors {
        let position_magnitude =
            (sv.position[0].powi(2) + sv.position[1].powi(2) + sv.position[2].powi(2)).sqrt();

        if position_magnitude < 6_500_000.0 || position_magnitude > 7_500_000.0 {
            log::warn!(
                "Unusual orbital radius: {:.1} km at {}",
                position_magnitude / 1000.0,
                sv.time.format("%Y-%m-%d %H:%M:%S")
            );
        }
    }

    log::debug!("Orbit data validation completed successfully");
    Ok(())
}

/// Enhanced state vector parsing with better error handling
pub(crate) fn parse_state_vector_line_enhanced(line: &str) -> SarResult<Option<StateVector>> {
    let mut time_opt = None;
    let mut position = [0.0; 3];
    let mut velocity = [0.0; 3];
    let mut found_fields = 0;

    for part in line.split_whitespace() {
        if let Some((key, value)) = part.split_once('=') {
            match key {
                "UTC" => {
                    if let Ok(naive_dt) =
                        NaiveDateTime::parse_from_str(value, "%Y-%m-%dT%H:%M:%S%.f")
                    {
                        time_opt = Some(Utc.from_utc_datetime(&naive_dt));
                        found_fields += 1;
                    } else {
                        log::warn!("Could not parse UTC time: {}", value);
                    }
                }
                "X" => {
                    position[0] = value.parse().map_err(|e| {
                        SarError::Processing(format!("Invalid X coordinate: {} ({})", value, e))
                    })?;
                    found_fields += 1;
                }
                "Y" => {
                    position[1] = value.parse().map_err(|e| {
                        SarError::Processing(format!("Invalid Y coordinate: {} ({})", value, e))
                    })?;
                    found_fields += 1;
                }
                "Z" => {
                    position[2] = value.parse().map_err(|e| {
                        SarError::Processing(format!("Invalid Z coordinate: {} ({})", value, e))
                    })?;
                    found_fields += 1;
                }
                "VX" => {
                    velocity[0] = value.parse().map_err(|e| {
                        SarError::Processing(format!("Invalid VX velocity: {} ({})", value, e))
                    })?;
                    found_fields += 1;
                }
                "VY" => {
                    velocity[1] = value.parse().map_err(|e| {
                        SarError::Processing(format!("Invalid VY velocity: {} ({})", value, e))
                    })?;
                    found_fields += 1;
                }
                "VZ" => {
                    velocity[2] = value.parse().map_err(|e| {
                        SarError::Processing(format!("Invalid VZ velocity: {} ({})", value, e))
                    })?;
                    found_fields += 1;
                }
                _ => {}
            }
        }
    }

    if found_fields < 7 {
        return Ok(None);
    }

    let time = time_opt
        .ok_or_else(|| SarError::Processing("Missing UTC time in state vector".to_string()))?;

    Ok(Some(StateVector {
        time,
        position,
        velocity,
    }))
}
