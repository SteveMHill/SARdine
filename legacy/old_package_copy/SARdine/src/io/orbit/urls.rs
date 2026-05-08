#![allow(dead_code, unused_variables)]
use chrono::{DateTime, NaiveDateTime, TimeZone, Timelike, Utc};
use regex::Regex;

use crate::types::{SarError, SarResult};

use super::OrbitType;

pub fn generate_orbit_urls(
    product_id: &str,
    start_time: DateTime<Utc>,
    orbit_type: OrbitType,
) -> SarResult<Vec<String>> {
    let satellite = if product_id.starts_with("S1A") {
        "S1A"
    } else if product_id.starts_with("S1B") {
        "S1B"
    } else {
        return Err(SarError::Processing(format!(
            "Cannot determine satellite from product ID: {}",
            product_id
        )));
    };

    let mut urls = Vec::new();
    let base_url = "https://step.esa.int/auxdata/orbits/Sentinel-1";

    let search_dates = vec![
        start_time - chrono::Duration::days(1),
        start_time,
        start_time + chrono::Duration::days(1),
    ];

    for search_date in search_dates {
        let year = search_date.format("%Y");
        let month = search_date.format("%m");

        let dir_url = format!(
            "{}/{}/{}/{}/{}/",
            base_url, orbit_type, satellite, year, month
        );

        log::debug!("Checking directory: {}", dir_url);

        if let Ok(file_urls) = get_orbit_files_from_directory(&dir_url, start_time, orbit_type) {
            urls.extend(file_urls);
        }
    }

    if urls.is_empty() {
        log::warn!("No orbit files found via directory listing, trying fallback patterns");
        urls.extend(generate_fallback_orbit_urls(
            product_id, start_time, orbit_type,
        )?);
    }

    if urls.is_empty() {
        return Err(SarError::Processing(format!(
            "No orbit file URLs generated for {} {} at {}",
            satellite,
            orbit_type,
            start_time.format("%Y-%m-%d")
        )));
    }

    log::info!("Generated {} potential orbit file URLs", urls.len());
    for (i, url) in urls.iter().enumerate() {
        log::debug!("URL {}: {}", i + 1, url);
    }

    Ok(urls)
}

pub(crate) fn get_orbit_files_from_directory(
    dir_url: &str,
    target_time: DateTime<Utc>,
    orbit_type: OrbitType,
) -> SarResult<Vec<String>> {
    let client = reqwest::blocking::Client::builder()
        .timeout(std::time::Duration::from_secs(30))
        .build()
        .map_err(|e| SarError::Processing(format!("Failed to create HTTP client: {}", e)))?;

    let response = client
        .get(dir_url)
        .send()
        .map_err(|e| SarError::Processing(format!("Failed to fetch directory listing: {}", e)))?;

    if !response.status().is_success() {
        return Err(SarError::Processing(format!(
            "Directory listing failed: {}",
            response.status()
        )));
    }

    let html = response
        .text()
        .map_err(|e| SarError::Processing(format!("Failed to read response: {}", e)))?;

    let mut matching_files = Vec::new();

    for line in html.lines() {
        if line.contains(".EOF.zip\"") {
            if let Some(filename) = extract_filename_from_html(line) {
                if is_orbit_file_relevant(&filename, target_time, orbit_type) {
                    let file_url = format!("{}{}", dir_url, filename);
                    matching_files.push(file_url);
                }
            }
        }
    }

    matching_files.sort_by(|a, b| {
        let score_a = score_orbit_file_relevance(a, target_time);
        let score_b = score_orbit_file_relevance(b, target_time);
        score_a.total_cmp(&score_b)
    });

    Ok(matching_files)
}

pub(crate) fn extract_filename_from_html(html_line: &str) -> Option<String> {
    if let Some(start) = html_line.find("href=\"") {
        let start = start + 6;
        if let Some(end) = html_line[start..].find('"') {
            let filename = &html_line[start..start + end];
            if filename.ends_with(".EOF.zip") && !filename.contains('/') {
                return Some(filename.to_string());
            }
        }
    }
    None
}

pub fn is_orbit_file_relevant(
    filename: &str,
    target_time: DateTime<Utc>,
    _orbit_type: OrbitType,
) -> bool {
    if let Some((start_time, end_time)) = extract_validity_from_name(filename) {
        return target_time >= start_time && target_time <= end_time;
    }

    let target_date = target_time.format("%Y%m%d").to_string();
    filename.contains(&target_date)
}

pub(crate) fn parse_orbit_filename_time(
    time_str: &str,
) -> Result<DateTime<Utc>, chrono::ParseError> {
    NaiveDateTime::parse_from_str(time_str, "%Y%m%dT%H%M%S").map(|dt| Utc.from_utc_datetime(&dt))
}

pub(crate) fn score_orbit_file_relevance(file_url: &str, target_time: DateTime<Utc>) -> f64 {
    let filename = file_url.split('/').next_back().unwrap_or("");

    if let Some((start_time, end_time)) = extract_validity_from_name(filename) {
        let mid_time = start_time + (end_time - start_time) / 2;
        return (target_time - mid_time).num_seconds().abs() as f64;
    }

    let target_date = target_time.format("%Y%m%d").to_string();
    if filename.contains(&target_date) {
        0.0
    } else {
        1_000_000.0
    }
}

pub fn generate_fallback_orbit_urls(
    product_id: &str,
    start_time: DateTime<Utc>,
    orbit_type: OrbitType,
) -> SarResult<Vec<String>> {
    let satellite = if product_id.starts_with("S1A") {
        "S1A"
    } else {
        "S1B"
    };
    let base_url = "https://step.esa.int/auxdata/orbits/Sentinel-1";

    let mut urls = Vec::new();
    let dates = vec![
        start_time - chrono::Duration::days(1),
        start_time,
        start_time + chrono::Duration::days(1),
    ];

    for date in dates {
        let year = date.format("%Y");
        let month = date.format("%m");

        for hour in [0, 6, 12, 18] {
            let validity_start = date.with_hour(hour).unwrap_or(date);
            let validity_end = validity_start + chrono::Duration::days(1);
            let production_time = validity_start + chrono::Duration::hours(3);

            let filename = format!(
                "S1{}_OPER_AUX_{}_OPOD_{}_V{}_{}.EOF.zip",
                satellite.chars().last().unwrap_or('A'),
                orbit_type,
                production_time.format("%Y%m%dT%H%M%S"),
                validity_start.format("%Y%m%dT%H%M%S"),
                validity_end.format("%Y%m%dT%H%M%S")
            );

            let url = format!(
                "{}/{}/{}/{}/{}/{}",
                base_url, orbit_type, satellite, year, month, filename
            );
            urls.push(url);
        }
    }

    Ok(urls)
}

/// Extract validity times from orbit filename using regex
pub fn extract_validity_from_name(filename: &str) -> Option<(DateTime<Utc>, DateTime<Utc>)> {
    let re = Regex::new(r"_V(\d{8}T\d{6})_(\d{8}T\d{6})\.EOF(?:\.zip)?$").ok()?;

    if let Some(captures) = re.captures(filename) {
        let start_str = captures.get(1)?.as_str();
        let end_str = captures.get(2)?.as_str();

        let parse_time = |s: &str| -> Option<DateTime<Utc>> {
            NaiveDateTime::parse_from_str(s, "%Y%m%dT%H%M%S")
                .map(|dt| Utc.from_utc_datetime(&dt))
                .ok()
        };

        if let (Some(start_time), Some(end_time)) = (parse_time(start_str), parse_time(end_str)) {
            return Some((start_time, end_time));
        }
    }
    None
}
