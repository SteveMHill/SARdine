//! Time conversion utilities for SAR processing
//!
//! Provides robust time parsing and conversion functions that maintain
//! consistency between different time domains:
//! - Orbit-reference-relative seconds (primary computational domain)
//! - Unix timestamps (legacy/interchange format)
//! - ISO8601/RFC3339 strings (annotation format)

use chrono::{DateTime, NaiveDateTime, Utc};

/// Convert a time difference to seconds with nanosecond precision
///
/// This is the primary time conversion helper for all solver-relevant timing.
/// Returns seconds relative to the given epoch, preserving sign and sub-second precision.
///
/// # Arguments
/// * `epoch` - The reference epoch (e.g., orbit reference time)
/// * `t` - The time to convert
///
/// # Returns
/// Seconds since epoch (can be negative if t < epoch)
///
/// # Example
/// ```
/// use chrono::{DateTime, Utc, TimeZone};
/// use sardine::io::sentinel::slc_reader::seconds_since_epoch;
///
/// let epoch = Utc.with_ymd_and_hms(2020, 1, 1, 0, 0, 0).unwrap();
/// let t = Utc.with_ymd_and_hms(2020, 1, 1, 0, 0, 10).unwrap();
/// let delta = seconds_since_epoch(epoch, t);
/// assert!((delta - 10.0).abs() < 1e-9);
/// ```
pub fn seconds_since_epoch(epoch: DateTime<Utc>, t: DateTime<Utc>) -> f64 {
    let delta = t - epoch;
    // Use nanoseconds for maximum precision when available
    if let Some(ns) = delta.num_nanoseconds() {
        ns as f64 * 1e-9
    } else {
        // Fallback for very large time differences (>292 years)
        delta.num_seconds() as f64
    }
}

/// Parse ISO8601/RFC3339 timestamp string to DateTime<Utc>
///
/// Accepts multiple formats commonly found in Sentinel-1 annotations:
/// - RFC3339 with timezone: "2020-01-01T12:00:00.123456Z"
/// - RFC3339 with offset: "2020-01-01T12:00:00.123456+00:00"
/// - Naive (treated as UTC): "2020-01-01T12:00:00.123456"
/// - Naive without fractional: "2020-01-01T12:00:00"
///
/// # Arguments
/// * `value` - The timestamp string to parse
///
/// # Returns
/// `Some(DateTime<Utc>)` if parsing succeeds, `None` otherwise
pub fn parse_iso8601_to_datetime_utc(value: &str) -> Option<DateTime<Utc>> {
    let trimmed = value.trim();
    if trimmed.is_empty() {
        return None;
    }

    // Try RFC3339 first (handles Z and offsets)
    if let Ok(dt) = DateTime::parse_from_rfc3339(trimmed) {
        return Some(dt.with_timezone(&Utc));
    }

    // Try naive formats (treat as UTC per Sentinel-1 convention)
    // With fractional seconds
    if let Ok(naive) = NaiveDateTime::parse_from_str(trimmed, "%Y-%m-%dT%H:%M:%S%.f") {
        return Some(DateTime::<Utc>::from_naive_utc_and_offset(naive, Utc));
    }

    // Without fractional seconds
    if let Ok(naive) = NaiveDateTime::parse_from_str(trimmed, "%Y-%m-%dT%H:%M:%S") {
        return Some(DateTime::<Utc>::from_naive_utc_and_offset(naive, Utc));
    }

    None
}

/// Convert DateTime<Utc> to Unix timestamp with nanosecond precision
///
/// This is used for legacy interchange formats that expect Unix seconds.
/// New code should prefer orbit-reference-relative times.
pub fn datetime_to_unix_seconds(dt: DateTime<Utc>) -> f64 {
    dt.timestamp() as f64 + (dt.timestamp_subsec_nanos() as f64) * 1e-9
}

#[cfg(test)]
mod tests {
    use super::*;
    use chrono::TimeZone;

    #[test]
    fn test_seconds_since_epoch_positive() {
        let epoch = Utc.with_ymd_and_hms(2020, 1, 1, 0, 0, 0).unwrap();
        let t = Utc.with_ymd_and_hms(2020, 1, 1, 0, 0, 10).unwrap();
        let delta = seconds_since_epoch(epoch, t);
        assert!((delta - 10.0).abs() < 1e-9, "Expected 10.0, got {}", delta);
    }

    #[test]
    fn test_seconds_since_epoch_negative() {
        let epoch = Utc.with_ymd_and_hms(2020, 1, 1, 0, 0, 10).unwrap();
        let t = Utc.with_ymd_and_hms(2020, 1, 1, 0, 0, 0).unwrap();
        let delta = seconds_since_epoch(epoch, t);
        assert!(
            (delta - (-10.0)).abs() < 1e-9,
            "Expected -10.0, got {}",
            delta
        );
    }

    #[test]
    fn test_seconds_since_epoch_fractional() {
        let epoch = Utc.with_ymd_and_hms(2020, 1, 1, 0, 0, 0).unwrap();
        // Add 500 milliseconds
        let t = epoch + chrono::Duration::milliseconds(500);
        let delta = seconds_since_epoch(epoch, t);
        assert!((delta - 0.5).abs() < 1e-9, "Expected 0.5, got {}", delta);
    }

    #[test]
    fn test_parse_iso8601_rfc3339_z() {
        let dt = parse_iso8601_to_datetime_utc("2020-01-15T12:30:45.123456Z").unwrap();
        assert_eq!(dt.year(), 2020);
        assert_eq!(dt.month(), 1);
        assert_eq!(dt.day(), 15);
        assert_eq!(dt.hour(), 12);
        assert_eq!(dt.minute(), 30);
        assert_eq!(dt.second(), 45);
    }

    #[test]
    fn test_parse_iso8601_naive() {
        let dt = parse_iso8601_to_datetime_utc("2020-01-15T12:30:45.123456").unwrap();
        assert_eq!(dt.hour(), 12);
    }

    #[test]
    fn test_parse_iso8601_no_fractional() {
        let dt = parse_iso8601_to_datetime_utc("2020-01-15T12:30:45").unwrap();
        assert_eq!(dt.second(), 45);
    }

    #[test]
    fn test_parse_iso8601_empty() {
        assert!(parse_iso8601_to_datetime_utc("").is_none());
        assert!(parse_iso8601_to_datetime_utc("   ").is_none());
    }

    #[test]
    fn test_parse_iso8601_invalid() {
        assert!(parse_iso8601_to_datetime_utc("not-a-date").is_none());
    }

    use chrono::Datelike;
    use chrono::Timelike;
}
