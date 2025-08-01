use crate::types::{OrbitData, SarError, SarResult, StateVector, BurstOrbitData, EofHeader};
use chrono::{DateTime, Utc, NaiveDateTime, Timelike};
use std::path::{Path, PathBuf};
use std::fs;

/// Orbit file types available from ESA
#[derive(Debug, Clone, Copy)]
pub enum OrbitType {
    /// Precise Orbit Ephemerides (best accuracy, ~20 days delay)
    POEORB,
    /// Restituted Orbit Ephemerides (lower accuracy, ~3 hours delay)  
    RESORB,
}

impl std::fmt::Display for OrbitType {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            OrbitType::POEORB => write!(f, "POEORB"),
            OrbitType::RESORB => write!(f, "RESORB"),
        }
    }
}

/// Precise orbit file reader for Sentinel-1
pub struct OrbitReader;

impl OrbitReader {
    /// Read precise orbit file (EOF format)
    pub fn read_orbit_file<P: AsRef<Path>>(path: P) -> SarResult<OrbitData> {
        log::info!("Reading orbit file: {}", path.as_ref().display());
        
        let content = fs::read_to_string(&path)
            .map_err(|e| SarError::Io(e))?;
            
        Self::parse_eof_content_enhanced(&content)
    }

    /// Automatically retrieve the best available orbit file for a product
    pub fn get_orbit_for_product(
        product_id: &str,
        start_time: DateTime<Utc>,
        output_dir: Option<&Path>,
    ) -> SarResult<OrbitData> {
        log::info!("Retrieving orbit file for product: {}", product_id);
        
        // Determine the best orbit type based on product age
        let orbit_type = Self::determine_orbit_type(start_time);
        log::info!("Selected orbit type: {}", orbit_type);
        
        // Try to download the selected orbit type
        let output_path = output_dir.map(|dir| {
            let filename = format!("{}_{}.EOF", product_id, orbit_type);
            dir.join(filename)
        });
        
        match Self::download_orbit_file(product_id, start_time, orbit_type, output_path.as_deref()) {
            Ok(orbit_data) => Ok(orbit_data),
            Err(e) => {
                // If POEORB fails, fallback to RESORB
                if matches!(orbit_type, OrbitType::POEORB) {
                    log::warn!("POEORB download failed: {}. Falling back to RESORB", e);
                    let fallback_path = output_dir.map(|dir| {
                        let filename = format!("{}_{}.EOF", product_id, OrbitType::RESORB);
                        dir.join(filename)
                    });
                    Self::download_orbit_file(product_id, start_time, OrbitType::RESORB, fallback_path.as_deref())
                } else {
                    Err(e)
                }
            }
        }
    }

    /// Generate orbit file URL for ESA download (legacy single URL version)
    pub fn generate_orbit_url(
        product_id: &str,
        orbit_type: OrbitType,
        start_time: DateTime<Utc>,
    ) -> String {
        // This is a simplified version for backwards compatibility
        let satellite = if product_id.starts_with("S1A") { "S1A" } else { "S1B" };
        let year = start_time.format("%Y");
        let base_url = "https://step.esa.int/auxdata/orbits/Sentinel-1";
        
        match orbit_type {
            OrbitType::POEORB => {
                format!("{}/{}/{}/{}/S1{}_OPER_AUX_POEORB_OPOD_{}_V{}.EOF",
                       base_url, orbit_type, satellite, year,
                       satellite.chars().last().unwrap(),
                       start_time.format("%Y%m%dT%H%M%S"),
                       (start_time + chrono::Duration::hours(24)).format("%Y%m%dT%H%M%S"))
            },
            OrbitType::RESORB => {
                format!("{}/{}/{}/{}/S1{}_OPER_AUX_RESORB_OPOD_{}_V{}.EOF",
                       base_url, orbit_type, satellite, year,
                       satellite.chars().last().unwrap(),
                       start_time.format("%Y%m%dT%H%M%S"),
                       (start_time + chrono::Duration::hours(6)).format("%Y%m%dT%H%M%S"))
            }
        }
    }

    /// Download orbit file from ESA servers
    pub fn download_orbit_file(
        product_id: &str,
        start_time: DateTime<Utc>,
        orbit_type: OrbitType,
        output_path: Option<&Path>,
    ) -> SarResult<OrbitData> {
        log::info!("Downloading {} orbit file for product: {}", orbit_type, product_id);
        
        // Try multiple URL patterns since ESA orbit file naming can vary
        let urls = Self::generate_orbit_urls(product_id, start_time, orbit_type)?;
        
        for (i, url) in urls.iter().enumerate() {
            log::info!("Attempting download from URL {}/{}: {}", i + 1, urls.len(), url);
            
            match Self::download_from_url(url, output_path) {
                Ok(content) => {
                    log::info!("Successfully downloaded orbit file from: {}", url);
                    return Self::parse_eof_content_enhanced(&content);
                },
                Err(e) => {
                    log::warn!("Failed to download from {}: {}", url, e);
                    continue;
                }
            }
        }
        
        Err(SarError::Processing(
            format!("Failed to download {} orbit file for {} from all attempted URLs", 
                   orbit_type, product_id)
        ))
    }
    
    /// Generate multiple possible orbit file URLs based on actual ESA structure
    pub fn generate_orbit_urls(
        product_id: &str,
        start_time: DateTime<Utc>,
        orbit_type: OrbitType,
    ) -> SarResult<Vec<String>> {
        // Extract satellite identifier from product ID
        let satellite = if product_id.starts_with("S1A") { 
            "S1A" 
        } else if product_id.starts_with("S1B") { 
            "S1B" 
        } else {
            return Err(SarError::Processing(
                format!("Cannot determine satellite from product ID: {}", product_id)
            ));
        };
        
        let mut urls = Vec::new();
        let base_url = "https://step.esa.int/auxdata/orbits/Sentinel-1";
        
        // Generate date range to search for orbit files
        // Orbit files typically cover 1-3 day periods
        let search_dates = vec![
            start_time - chrono::Duration::days(1),
            start_time,
            start_time + chrono::Duration::days(1),
        ];
        
        for search_date in search_dates {
            let year = search_date.format("%Y");
            let month = search_date.format("%m");
            
            // Create directory listing URL to find actual files
            let dir_url = format!("{}/{}/{}/{}/{}/", 
                                 base_url, orbit_type, satellite, year, month);
            
            log::debug!("Checking directory: {}", dir_url);
            
            // Try to get directory listing and find matching files
            if let Ok(file_urls) = Self::get_orbit_files_from_directory(&dir_url, start_time, orbit_type) {
                urls.extend(file_urls);
            }
        }
        
        // If no specific files found, try some common patterns
        if urls.is_empty() {
            log::warn!("No orbit files found via directory listing, trying fallback patterns");
            urls.extend(Self::generate_fallback_orbit_urls(product_id, start_time, orbit_type)?);
        }
        
        if urls.is_empty() {
            return Err(SarError::Processing(
                format!("No orbit file URLs generated for {} {} at {}", 
                       satellite, orbit_type, start_time.format("%Y-%m-%d"))
            ));
        }
        
        log::info!("Generated {} potential orbit file URLs", urls.len());
        for (i, url) in urls.iter().enumerate() {
            log::debug!("URL {}: {}", i + 1, url);
        }
        
        Ok(urls)
    }
    
    /// Get orbit files from ESA directory listing
    fn get_orbit_files_from_directory(
        dir_url: &str,
        target_time: DateTime<Utc>,
        orbit_type: OrbitType,
    ) -> SarResult<Vec<String>> {
        let client = reqwest::blocking::Client::builder()
            .timeout(std::time::Duration::from_secs(30))
            .build()
            .map_err(|e| SarError::Processing(format!("Failed to create HTTP client: {}", e)))?;
        
        let response = client.get(dir_url)
            .send()
            .map_err(|e| SarError::Processing(format!("Failed to fetch directory listing: {}", e)))?;
        
        if !response.status().is_success() {
            return Err(SarError::Processing(
                format!("Directory listing failed: {}", response.status())
            ));
        }
        
        let html = response.text()
            .map_err(|e| SarError::Processing(format!("Failed to read response: {}", e)))?;
        
        let mut matching_files = Vec::new();
        
        // Parse HTML to find .EOF.zip files
        for line in html.lines() {
            if line.contains(".EOF.zip\"") {
                if let Some(filename) = Self::extract_filename_from_html(line) {
                    if Self::is_orbit_file_relevant(&filename, target_time, orbit_type) {
                        let file_url = format!("{}{}", dir_url, filename);
                        matching_files.push(file_url);
                    }
                }
            }
        }
        
        // Sort by relevance (closest to target time)
        matching_files.sort_by(|a, b| {
            let score_a = Self::score_orbit_file_relevance(a, target_time);
            let score_b = Self::score_orbit_file_relevance(b, target_time);
            score_a.partial_cmp(&score_b).unwrap()
        });
        
        Ok(matching_files)
    }
    
    /// Extract filename from HTML link
    fn extract_filename_from_html(html_line: &str) -> Option<String> {
        // Look for pattern: href="filename.EOF.zip">
        if let Some(start) = html_line.find("href=\"") {
            let start = start + 6; // Skip 'href="'
            if let Some(end) = html_line[start..].find('\"') {
                let filename = &html_line[start..start + end];
                if filename.ends_with(".EOF.zip") && !filename.contains('/') {
                    return Some(filename.to_string());
                }
            }
        }
        None
    }
    
    /// Check if orbit file is relevant for the target time
    fn is_orbit_file_relevant(
        filename: &str,
        target_time: DateTime<Utc>,
        _orbit_type: OrbitType,
    ) -> bool {
        // Parse orbit file validity period from filename
        // Format: S1A_OPER_AUX_POEORB_OPOD_[prod]_V[start]_[end].EOF.zip
        let parts: Vec<&str> = filename.split('_').collect();
        if parts.len() >= 7 {
            // Extract validity start and end times
            if let (Ok(start_time), Ok(end_time)) = (
                Self::parse_orbit_filename_time(&parts[5][1..]), // Skip 'V' prefix
                Self::parse_orbit_filename_time(&parts[6].replace(".EOF.zip", ""))
            ) {
                // Check if target time falls within validity period
                return target_time >= start_time && target_time <= end_time;
            }
        }
        
        // Fallback: check if filename contains target date
        let target_date = target_time.format("%Y%m%d").to_string();
        filename.contains(&target_date)
    }
    
    /// Parse time from orbit filename format (YYYYMMDDTHHMMSS)
    fn parse_orbit_filename_time(time_str: &str) -> Result<DateTime<Utc>, chrono::ParseError> {
        NaiveDateTime::parse_from_str(time_str, "%Y%m%dT%H%M%S")
            .map(|dt| DateTime::from_naive_utc_and_offset(dt, Utc))
    }
    
    /// Score orbit file relevance (lower is better)
    fn score_orbit_file_relevance(file_url: &str, target_time: DateTime<Utc>) -> f64 {
        // Extract filename from URL
        let filename = file_url.split('/').last().unwrap_or("");
        
        // Try to parse validity period
        let parts: Vec<&str> = filename.split('_').collect();
        if parts.len() >= 7 {
            if let (Ok(start_time), Ok(end_time)) = (
                Self::parse_orbit_filename_time(&parts[5][1..]),
                Self::parse_orbit_filename_time(&parts[6].replace(".EOF.zip", ""))
            ) {
                let mid_time = start_time + (end_time - start_time) / 2;
                return (target_time - mid_time).num_seconds().abs() as f64;
            }
        }
        
        // Fallback scoring based on filename similarity
        let target_date = target_time.format("%Y%m%d").to_string();
        if filename.contains(&target_date) {
            0.0
        } else {
            1000000.0
        }
    }
    
    /// Generate fallback URLs when directory listing fails
    fn generate_fallback_orbit_urls(
        product_id: &str,
        start_time: DateTime<Utc>,
        orbit_type: OrbitType,
    ) -> SarResult<Vec<String>> {
        let satellite = if product_id.starts_with("S1A") { "S1A" } else { "S1B" };
        let base_url = "https://step.esa.int/auxdata/orbits/Sentinel-1";
        
        let mut urls = Vec::new();
        let dates = vec![
            start_time - chrono::Duration::days(1), 
            start_time, 
            start_time + chrono::Duration::days(1)
        ];
        
        for date in dates {
            let year = date.format("%Y");
            let month = date.format("%m");
            let _day = date.format("%d");
            
            // Generate potential orbit file names based on ESA naming conventions
            // Format: S1A_OPER_AUX_POEORB_OPOD_[prod]_V[start]_[end].EOF.zip
            for hour in [0, 6, 12, 18] {
                let start_time = date.with_hour(hour).unwrap_or(date);
                let end_time = start_time + chrono::Duration::days(1);
                
                let filename = format!(
                    "S1{}_OPER_AUX_{}_OPOD_{}_V{}__{}.EOF.zip",
                    satellite.chars().last().unwrap(),
                    orbit_type,
                    start_time.format("%Y%m%dT%H%M%S"),
                    start_time.format("%Y%m%dT%H%M%S"),
                    end_time.format("%Y%m%dT%H%M%S")
                );
                
                let url = format!("{}/{}/{}/{}/{}/{}", 
                                 base_url, orbit_type, satellite, year, month, filename);
                urls.push(url);
            }
        }
        
        Ok(urls)
    }
    
    /// Download content from a URL
    fn download_from_url(url: &str, output_path: Option<&Path>) -> SarResult<String> {
        log::debug!("Downloading from: {}", url);
        
        let response = reqwest::blocking::get(url)
            .map_err(|e| SarError::Processing(format!("HTTP request failed: {}", e)))?;
            
        if !response.status().is_success() {
            return Err(SarError::Processing(
                format!("HTTP request failed with status: {}", response.status())
            ));
        }
        
        let bytes = response.bytes()
            .map_err(|e| SarError::Processing(format!("Failed to read response bytes: {}", e)))?;
        
        // Check if this is a ZIP file based on URL or content
        let content = if url.ends_with(".EOF.zip") || Self::is_zip_content(&bytes) {
            log::debug!("Processing ZIP file from: {}", url);
            Self::extract_eof_from_zip(&bytes)?
        } else {
            // Plain text EOF file
            String::from_utf8(bytes.to_vec())
                .map_err(|e| SarError::Processing(format!("Invalid UTF-8 content: {}", e)))?
        };
        
        // Save to file if output path is provided
        if let Some(path) = output_path {
            if let Some(parent) = path.parent() {
                std::fs::create_dir_all(parent)
                    .map_err(|e| SarError::Io(e))?;
            }
            
            std::fs::write(path, &content)
                .map_err(|e| SarError::Io(e))?;
                
            log::info!("Orbit file saved to: {}", path.display());
        }
        
        Ok(content)
    }
    
    /// Check if content is a ZIP file by examining magic bytes
    fn is_zip_content(bytes: &[u8]) -> bool {
        bytes.len() >= 4 && bytes[0..4] == [0x50, 0x4B, 0x03, 0x04] // ZIP magic signature
    }
    
    /// Extract EOF content from ZIP file
    fn extract_eof_from_zip(zip_bytes: &[u8]) -> SarResult<String> {
        use std::io::Cursor;
        use zip::ZipArchive;
        
        let cursor = Cursor::new(zip_bytes);
        let mut archive = ZipArchive::new(cursor)
            .map_err(|e| SarError::Processing(format!("Failed to read ZIP archive: {}", e)))?;
        
        // Look for .EOF file in the archive
        for i in 0..archive.len() {
            let mut file = archive.by_index(i)
                .map_err(|e| SarError::Processing(format!("Failed to read ZIP entry {}: {}", i, e)))?;
            
            if file.name().ends_with(".EOF") {
                log::debug!("Found EOF file in ZIP: {}", file.name());
                
                use std::io::Read;
                let mut contents = String::new();
                file.read_to_string(&mut contents)
                    .map_err(|e| SarError::Processing(format!("Failed to read EOF file: {}", e)))?;
                
                return Ok(contents);
            }
        }
        
        Err(SarError::Processing(
            "No .EOF file found in ZIP archive".to_string()
        ))
    }

    /// Interpolate orbit position at specific time using Lagrange interpolation
    pub fn interpolate_position(
        orbit: &OrbitData,
        target_time: DateTime<Utc>,
    ) -> SarResult<[f64; 3]> {
        if orbit.state_vectors.is_empty() {
            return Err(SarError::Processing(
                "No state vectors in orbit data".to_string()
            ));
        }
        
        let target_timestamp = target_time.timestamp_millis() as f64 / 1000.0;
        let selected_svs = Self::find_interpolation_vectors(&orbit.state_vectors, target_timestamp)?;
        Self::lagrange_interpolate(&selected_svs, target_timestamp)
    }

    /// Fast binary search to find the best state vectors for interpolation
    fn find_interpolation_vectors(
        state_vectors: &[StateVector],
        target_timestamp: f64,
    ) -> SarResult<Vec<&StateVector>> {
        if state_vectors.is_empty() {
            return Err(SarError::Processing("No state vectors available".to_string()));
        }

        // Binary search to find the closest state vector
        let closest_idx = Self::binary_search_closest_time(state_vectors, target_timestamp);
        
        // Select 4 vectors around the target time for interpolation
        let num_points = std::cmp::min(4, state_vectors.len());
        if num_points < 2 {
            log::warn!("Not enough state vectors for interpolation, using nearest");
            return Ok(vec![&state_vectors[closest_idx]]);
        }

        // Get indices for interpolation points (try to center around target)
        let start_idx = if closest_idx >= num_points / 2 {
            std::cmp::min(closest_idx - num_points / 2, state_vectors.len() - num_points)
        } else {
            0
        };
        
        let selected_svs: Vec<&StateVector> = (start_idx..start_idx + num_points)
            .map(|i| &state_vectors[i])
            .collect();

        Ok(selected_svs)
    }

    /// Binary search to find the state vector closest in time to target
    fn binary_search_closest_time(
        state_vectors: &[StateVector],
        target_timestamp: f64,
    ) -> usize {
        if state_vectors.is_empty() {
            return 0;
        }

        let mut left = 0;
        let mut right = state_vectors.len() - 1;
        
        while left < right {
            let mid = left + (right - left) / 2;
            let mid_timestamp = state_vectors[mid].time.timestamp_millis() as f64 / 1000.0;
            
            if mid_timestamp < target_timestamp {
                left = mid + 1;
            } else {
                right = mid;
            }
        }

        // Check if left-1 is closer than left
        if left > 0 {
            let left_dist = (state_vectors[left].time.timestamp_millis() as f64 / 1000.0 - target_timestamp).abs();
            let left_minus_1_dist = (state_vectors[left - 1].time.timestamp_millis() as f64 / 1000.0 - target_timestamp).abs();
            
            if left_minus_1_dist < left_dist {
                return left - 1;
            }
        }

        left
    }
    
    /// Lagrange polynomial interpolation for orbit position
    fn lagrange_interpolate(
        state_vectors: &[&StateVector], 
        target_time: f64
    ) -> SarResult<[f64; 3]> {
        let n = state_vectors.len();
        let mut result = [0.0; 3];
        
        for i in 0..n {
            let sv_i = state_vectors[i];
            let ti = sv_i.time.timestamp_millis() as f64 / 1000.0;
            let mut li = 1.0;
            
            // Calculate Lagrange basis polynomial
            for j in 0..n {
                if i != j {
                    let tj = state_vectors[j].time.timestamp_millis() as f64 / 1000.0;
                    li *= (target_time - tj) / (ti - tj);
                }
            }
            
            // Add contribution to each coordinate
            for coord in 0..3 {
                result[coord] += li * sv_i.position[coord];
            }
        }
        
        log::debug!("Interpolated position at {}: {:?}", target_time, result);
        Ok(result)
    }
    
    /// Interpolate orbit velocity at specific time using Lagrange interpolation
    pub fn interpolate_velocity(
        orbit: &OrbitData,
        target_time: DateTime<Utc>,
    ) -> SarResult<[f64; 3]> {
        if orbit.state_vectors.is_empty() {
            return Err(SarError::Processing(
                "No state vectors in orbit data".to_string()
            ));
        }
        
        let target_timestamp = target_time.timestamp_millis() as f64 / 1000.0;
        let selected_svs = Self::find_interpolation_vectors(&orbit.state_vectors, target_timestamp)?;
        Self::lagrange_interpolate_velocity(&selected_svs, target_timestamp)
    }
    
    /// Lagrange interpolation for velocity
    fn lagrange_interpolate_velocity(
        state_vectors: &[&StateVector],
        target_time: f64
    ) -> SarResult<[f64; 3]> {
        let n = state_vectors.len();
        let mut result = [0.0; 3];
        
        for i in 0..n {
            let sv_i = state_vectors[i];
            let ti = sv_i.time.timestamp_millis() as f64 / 1000.0;
            let mut li = 1.0;
            
            for j in 0..n {
                if i != j {
                    let tj = state_vectors[j].time.timestamp_millis() as f64 / 1000.0;
                    li *= (target_time - tj) / (ti - tj);
                }
            }
            
            for coord in 0..3 {
                result[coord] += li * sv_i.velocity[coord];
            }
        }
        
        Ok(result)
    }
    
    /// Determine which orbit type to use based on product age
    pub fn determine_orbit_type(acquisition_time: DateTime<Utc>) -> OrbitType {
        let now = Utc::now();
        let age_days = (now - acquisition_time).num_days();
        
        if age_days > 20 {
            OrbitType::POEORB // Precise orbits should be available
        } else {
            OrbitType::RESORB // Use restituted orbits for recent data
        }
    }
    
    /// Get satellite position and velocity for a burst at specific azimuth times
    /// This is essential for SAR processing - each pixel row has a different azimuth time
    /// OPTIMIZED VERSION: Uses batch processing and pre-sorted state vectors
    pub fn interpolate_burst_orbit(
        orbit: &OrbitData,
        burst_start_time: DateTime<Utc>,
        azimuth_time_interval: f64, // seconds between azimuth samples
        num_azimuth_lines: usize,
    ) -> SarResult<BurstOrbitData> {
        log::info!("Interpolating orbit for burst: {} lines, {} s interval", 
                  num_azimuth_lines, azimuth_time_interval);
        
        if orbit.state_vectors.is_empty() {
            return Err(SarError::Processing("No state vectors in orbit data".to_string()));
        }

        // Pre-sort state vectors by time for efficient binary search
        let mut sorted_state_vectors = orbit.state_vectors.clone();
        sorted_state_vectors.sort_by(|a, b| a.time.cmp(&b.time));

        let mut positions = Vec::with_capacity(num_azimuth_lines);
        let mut velocities = Vec::with_capacity(num_azimuth_lines);
        let mut azimuth_times = Vec::with_capacity(num_azimuth_lines);
        
        // Batch processing: find interpolation range once for the entire burst
        let burst_start_timestamp = burst_start_time.timestamp_millis() as f64 / 1000.0;
        let burst_end_timestamp = burst_start_timestamp + (num_azimuth_lines as f64 * azimuth_time_interval);
        
        // Find the range of state vectors we'll need for this entire burst
        let start_search_idx = Self::binary_search_closest_time(&sorted_state_vectors, burst_start_timestamp);
        let end_search_idx = Self::binary_search_closest_time(&sorted_state_vectors, burst_end_timestamp);
        
        // Expand search range to ensure we have enough vectors for interpolation
        let search_start = if start_search_idx >= 2 { start_search_idx - 2 } else { 0 };
        let search_end = std::cmp::min(end_search_idx + 3, sorted_state_vectors.len());
        
        let relevant_vectors = &sorted_state_vectors[search_start..search_end];
        
        log::debug!("Using {} state vectors for {} azimuth lines (optimization: {:.1}x reduction)",
                   relevant_vectors.len(), num_azimuth_lines, 
                   sorted_state_vectors.len() as f64 / relevant_vectors.len() as f64);

        for line_idx in 0..num_azimuth_lines {
            // Calculate azimuth time for this pixel row
            let azimuth_time = burst_start_time + 
                chrono::Duration::milliseconds((line_idx as f64 * azimuth_time_interval * 1000.0) as i64);
            
            let target_timestamp = azimuth_time.timestamp_millis() as f64 / 1000.0;
            
            // Find interpolation vectors from the reduced set
            let selected_svs = Self::find_interpolation_vectors_from_subset(relevant_vectors, target_timestamp)?;
            
            // Interpolate satellite position and velocity at this time
            let position = Self::lagrange_interpolate(&selected_svs, target_timestamp)?;
            let velocity = Self::lagrange_interpolate_velocity(&selected_svs, target_timestamp)?;
            
            positions.push(position);
            velocities.push(velocity);
            azimuth_times.push(azimuth_time);
            
            if line_idx % 1000 == 0 {
                log::debug!("Line {}: pos=[{:.1}, {:.1}, {:.1}] km", 
                           line_idx, position[0]/1000.0, position[1]/1000.0, position[2]/1000.0);
            }
        }
        
        Ok(BurstOrbitData {
            positions,
            velocities,
            azimuth_times,
            burst_start_time,
            azimuth_time_interval,
        })
    }

    /// Find interpolation vectors from a pre-filtered subset (for batch processing)
    fn find_interpolation_vectors_from_subset(
        state_vectors: &[StateVector],
        target_timestamp: f64,
    ) -> SarResult<Vec<&StateVector>> {
        if state_vectors.is_empty() {
            return Err(SarError::Processing("No state vectors available".to_string()));
        }

        // Binary search within the subset
        let closest_idx = Self::binary_search_closest_time(state_vectors, target_timestamp);
        
        // Select 4 vectors around the target time for interpolation
        let num_points = std::cmp::min(4, state_vectors.len());
        if num_points < 2 {
            log::warn!("Not enough state vectors for interpolation, using nearest");
            return Ok(vec![&state_vectors[closest_idx]]);
        }

        // Get indices for interpolation points (try to center around target)
        let start_idx = if closest_idx >= num_points / 2 {
            std::cmp::min(closest_idx - num_points / 2, state_vectors.len() - num_points)
        } else {
            0
        };
        
        let selected_svs: Vec<&StateVector> = (start_idx..start_idx + num_points)
            .map(|i| &state_vectors[i])
            .collect();

        Ok(selected_svs)
    }
    
    /// Calculate Doppler centroid frequency for a given satellite velocity and look direction
    /// This is useful for accurate SAR focusing and geolocation
    pub fn calculate_doppler_centroid(
        satellite_velocity: [f64; 3],
        look_direction: [f64; 3], // unit vector pointing towards ground target
        wavelength: f64, // radar wavelength in meters (C-band ≈ 0.055 m)
    ) -> f64 {
        // Doppler frequency = 2 * (v_sat · look_dir) / λ
        let velocity_dot_look = satellite_velocity[0] * look_direction[0] +
                               satellite_velocity[1] * look_direction[1] +
                               satellite_velocity[2] * look_direction[2];
        
        2.0 * velocity_dot_look / wavelength
    }
    
    /// Create a mock orbit file for testing purposes
    /// This generates synthetic but realistic state vectors for a Sentinel-1-like orbit
    pub fn create_mock_orbit_file(
        start_time: DateTime<Utc>,
        duration_hours: f64,
        output_path: &Path,
    ) -> SarResult<()> {
        log::info!("Creating mock orbit file: {}", output_path.display());
        
        // Create directory if it doesn't exist
        if let Some(parent) = output_path.parent() {
            std::fs::create_dir_all(parent)?;
        }
        
        // Sentinel-1 orbital parameters (approximate)
        let orbit_period = 5940.0; // seconds (about 99 minutes)
        let orbit_radius = 7_070_000.0; // meters (700 km altitude)
        let _earth_radius = 6_371_000.0; // meters
        let inclination = 98.18_f64.to_radians(); // degrees to radians
        
        let mut content = String::new();
        content.push_str("<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n");
        content.push_str("<Earth_Explorer_File>\n");
        content.push_str("  <Earth_Fixed_Coordinate_System>ITRF</Earth_Fixed_Coordinate_System>\n");
        content.push_str("  <Time_Reference>UTC</Time_Reference>\n");
        content.push_str("  <Data_Block type=\"xml\">\n");
        content.push_str("    <List_of_OSVs count=\"72\">\n");
        
        // Generate state vectors every 5 minutes for the specified duration
        let vector_interval = 300.0; // seconds
        let num_vectors = ((duration_hours * 3600.0) / vector_interval) as usize;
        
        for i in 0..num_vectors {
            let time_offset = i as f64 * vector_interval;
            let current_time = start_time + chrono::Duration::seconds(time_offset as i64);
            
            // Simple circular orbit model for testing
            let angle = (time_offset / orbit_period) * 2.0 * std::f64::consts::PI;
            
            // Position in orbital plane
            let x_orbit = orbit_radius * angle.cos();
            let y_orbit = orbit_radius * angle.sin();
            let z_orbit = 0.0;
            
            // Apply inclination (rotate around x-axis)
            let x = x_orbit;
            let y = y_orbit * inclination.cos() - z_orbit * inclination.sin();
            let z = y_orbit * inclination.sin() + z_orbit * inclination.cos();
            
            // Velocity (derivative of position)
            let orbital_velocity = 2.0 * std::f64::consts::PI * orbit_radius / orbit_period;
            let vx = -orbital_velocity * angle.sin();
            let vy = orbital_velocity * angle.cos() * inclination.cos();
            let vz = orbital_velocity * angle.cos() * inclination.sin();
            
            content.push_str(&format!(
                "      <OSV>\n"));
            content.push_str(&format!(
                "        <TAI>TAI={}</TAI>\n", 
                current_time.format("%Y-%m-%dT%H:%M:%S.%6f")));
            content.push_str(&format!(
                "        <UTC>UTC={}</UTC>\n", 
                current_time.format("%Y-%m-%dT%H:%M:%S.%6f")));
            content.push_str(&format!(
                "        <UT1>UT1={}</UT1>\n", 
                current_time.format("%Y-%m-%dT%H:%M:%S.%6f")));
            content.push_str(&format!(
                "        <Absolute_Orbit>+{}</Absolute_Orbit>\n", 30000 + i));
            content.push_str(&format!(
                "        <X unit=\"m\">{:.6}</X>\n", x));
            content.push_str(&format!(
                "        <Y unit=\"m\">{:.6}</Y>\n", y));
            content.push_str(&format!(
                "        <Z unit=\"m\">{:.6}</Z>\n", z));
            content.push_str(&format!(
                "        <VX unit=\"m/s\">{:.6}</VX>\n", vx));
            content.push_str(&format!(
                "        <VY unit=\"m/s\">{:.6}</VY>\n", vy));
            content.push_str(&format!(
                "        <VZ unit=\"m/s\">{:.6}</VZ>\n", vz));
            content.push_str("      </OSV>\n");
        }
        
        content.push_str("    </List_of_OSVs>\n");
        content.push_str("  </Data_Block>\n");
        content.push_str("</Earth_Explorer_File>\n");
        
        std::fs::write(output_path, content)?;
        log::info!("Mock orbit file created with {} state vectors", num_vectors);
        Ok(())
    }
    
    /// Enhanced EOF parsing that handles XML format and simple format
    fn parse_eof_content_enhanced(content: &str) -> SarResult<OrbitData> {
        let mut state_vectors = Vec::new();
        let mut reference_time = None;
        let mut header_info = EofHeader::default();
        
        log::info!("Parsing EOF orbit file ({} bytes)", content.len());
        
        let mut in_osv_block = false;
        let mut current_osv_data: Option<(DateTime<Utc>, [f64; 3], [f64; 3])> = None;
        let mut line_count = 0;
        
        for line in content.lines() {
            line_count += 1;
            let line = line.trim();
            
            // Skip comments and empty lines
            if line.starts_with('#') || line.is_empty() {
                continue;
            }
            
            // Parse header information
            if line.starts_with("<Earth_Fixed_Coordinate_System>") {
                header_info.coordinate_system = Some("Earth_Fixed".to_string());
            } else if line.starts_with("<Time_Reference>") {
                if let Some(time_ref) = Self::extract_xml_value(line, "Time_Reference") {
                    header_info.time_reference = Some(time_ref);
                }
            } else if line.starts_with("<Orbit_File_Name>") {
                if let Some(file_name) = Self::extract_xml_value(line, "Orbit_File_Name") {
                    header_info.file_name = Some(file_name);
                }
            }
            
            // Look for OSV (Orbit State Vector) blocks
            if line.contains("<OSV>") {
                in_osv_block = true;
                current_osv_data = None;
                log::debug!("Found OSV block at line {}", line_count);
                continue;
            }
            
            if line.contains("</OSV>") {
                // Finalize current OSV if we have complete data
                if let Some((time, position, velocity)) = current_osv_data {
                    if reference_time.is_none() {
                        reference_time = Some(time);
                    }
                    state_vectors.push(StateVector {
                        time,
                        position,
                        velocity,
                    });
                }
                in_osv_block = false;
                current_osv_data = None;
                continue;
            }
            
            // Parse OSV data elements
            if in_osv_block {
                // Handle XML format: <UTC>UTC=2020-01-03T17:00:00.000000</UTC>
                if line.starts_with("<UTC>") && line.contains("UTC=") {
                    if let Some(utc_value) = Self::extract_xml_value(line, "UTC") {
                        // Extract the actual time from "UTC=2020-01-03T17:00:00.000000"
                        if let Some(time_str) = utc_value.strip_prefix("UTC=") {
                            if let Ok(naive_dt) = NaiveDateTime::parse_from_str(time_str, "%Y-%m-%dT%H:%M:%S%.f") {
                                let time = DateTime::from_naive_utc_and_offset(naive_dt, Utc);
                                current_osv_data = Some((time, [0.0; 3], [0.0; 3]));
                            }
                        }
                    }
                } else if line.starts_with("<X ") || line.starts_with("<X>") {
                    if let Some(x_value) = Self::extract_xml_value(line, "X") {
                        // Handle both "X=value" and just "value" formats
                        let x_str = if x_value.starts_with("X=") {
                            x_value.strip_prefix("X=").unwrap_or(&x_value)
                        } else {
                            &x_value
                        };
                        if let Ok(x) = x_str.parse::<f64>() {
                            if let Some((time, ref mut position, velocity)) = current_osv_data {
                                position[0] = x;
                                current_osv_data = Some((time, *position, velocity));
                            }
                        }
                    }
                } else if line.starts_with("<Y ") || line.starts_with("<Y>") {
                    if let Some(y_value) = Self::extract_xml_value(line, "Y") {
                        let y_str = if y_value.starts_with("Y=") {
                            y_value.strip_prefix("Y=").unwrap_or(&y_value)
                        } else {
                            &y_value
                        };
                        if let Ok(y) = y_str.parse::<f64>() {
                            if let Some((time, ref mut position, velocity)) = current_osv_data {
                                position[1] = y;
                                current_osv_data = Some((time, *position, velocity));
                            }
                        }
                    }
                } else if line.starts_with("<Z ") || line.starts_with("<Z>") {
                    if let Some(z_value) = Self::extract_xml_value(line, "Z") {
                        let z_str = if z_value.starts_with("Z=") {
                            z_value.strip_prefix("Z=").unwrap_or(&z_value)
                        } else {
                            &z_value
                        };
                        if let Ok(z) = z_str.parse::<f64>() {
                            if let Some((time, ref mut position, velocity)) = current_osv_data {
                                position[2] = z;
                                current_osv_data = Some((time, *position, velocity));
                            }
                        }
                    }
                } else if line.starts_with("<VX ") || line.starts_with("<VX>") {
                    if let Some(vx_value) = Self::extract_xml_value(line, "VX") {
                        let vx_str = if vx_value.starts_with("VX=") {
                            vx_value.strip_prefix("VX=").unwrap_or(&vx_value)
                        } else {
                            &vx_value
                        };
                        if let Ok(vx) = vx_str.parse::<f64>() {
                            if let Some((time, position, ref mut velocity)) = current_osv_data {
                                velocity[0] = vx;
                                current_osv_data = Some((time, position, *velocity));
                            }
                        }
                    }
                } else if line.starts_with("<VY ") || line.starts_with("<VY>") {
                    if let Some(vy_value) = Self::extract_xml_value(line, "VY") {
                        let vy_str = if vy_value.starts_with("VY=") {
                            vy_value.strip_prefix("VY=").unwrap_or(&vy_value)
                        } else {
                            &vy_value
                        };
                        if let Ok(vy) = vy_str.parse::<f64>() {
                            if let Some((time, position, ref mut velocity)) = current_osv_data {
                                velocity[1] = vy;
                                current_osv_data = Some((time, position, *velocity));
                            }
                        }
                    }
                } else if line.starts_with("<VZ ") || line.starts_with("<VZ>") {
                    if let Some(vz_value) = Self::extract_xml_value(line, "VZ") {
                        let vz_str = if vz_value.starts_with("VZ=") {
                            vz_value.strip_prefix("VZ=").unwrap_or(&vz_value)
                        } else {
                            &vz_value
                        };
                        if let Ok(vz) = vz_str.parse::<f64>() {
                            if let Some((time, position, ref mut velocity)) = current_osv_data {
                                velocity[2] = vz;
                                current_osv_data = Some((time, position, *velocity));
                            }
                        }
                    }
                }
            }
            
            // Also try old format for backward compatibility
            if line.starts_with("UTC=") {
                match Self::parse_state_vector_line_enhanced(line) {
                    Ok(Some(sv)) => {
                        if reference_time.is_none() {
                            reference_time = Some(sv.time);
                        }
                        state_vectors.push(sv);
                    },
                    Ok(None) => continue,
                    Err(e) => {
                        log::warn!("Failed to parse state vector at line {}: {}", line_count, e);
                        continue;
                    }
                }
            }
        }
        
        if state_vectors.is_empty() {
            return Err(SarError::Processing(
                format!("No valid state vectors found in {} lines of orbit file. Expected XML format with <OSV> blocks or UTC= lines.", line_count)
            ));
        }
        
        // Sort state vectors by time to ensure proper ordering
        state_vectors.sort_by_key(|sv| sv.time);
        
        let time_span = if state_vectors.len() > 1 {
            (state_vectors.last().unwrap().time - state_vectors.first().unwrap().time).num_seconds()
        } else {
            0
        };
        
        log::info!("Successfully parsed {} state vectors spanning {} seconds", 
                  state_vectors.len(), time_span);
        log::info!("Time range: {} to {}", 
                  state_vectors.first().unwrap().time.format("%Y-%m-%d %H:%M:%S"),
                  state_vectors.last().unwrap().time.format("%Y-%m-%d %H:%M:%S"));
        
        // Validate orbit data quality
        Self::validate_orbit_data(&state_vectors)?;
        
        Ok(OrbitData {
            state_vectors,
            reference_time: reference_time.unwrap_or_else(|| Utc::now()),
        })
    }
    
    /// Extract value from XML element (handles tags with attributes)
    fn extract_xml_value(line: &str, tag: &str) -> Option<String> {
        // Look for opening tag (with or without attributes)
        let start_pattern = format!("<{}", tag);
        let end_tag = format!("</{}>", tag);
        
        if let Some(start) = line.find(&start_pattern) {
            // Find the end of the opening tag
            if let Some(tag_end) = line[start..].find('>') {
                let content_start = start + tag_end + 1;
                if let Some(end) = line[content_start..].find(&end_tag) {
                    return Some(line[content_start..content_start + end].to_string());
                }
            }
        }
        None
    }
    
    /// Validate orbit data quality
    fn validate_orbit_data(state_vectors: &[StateVector]) -> SarResult<()> {
        if state_vectors.is_empty() {
            return Err(SarError::Processing("No state vectors to validate".to_string()));
        }
        
        // Check for reasonable orbital velocities (should be ~7.5 km/s for LEO)
        for sv in state_vectors {
            let velocity_magnitude = (sv.velocity[0].powi(2) + 
                                    sv.velocity[1].powi(2) + 
                                    sv.velocity[2].powi(2)).sqrt();
            
            if velocity_magnitude < 6000.0 || velocity_magnitude > 9000.0 {
                log::warn!("Unusual orbital velocity: {:.1} m/s at {}", 
                          velocity_magnitude, sv.time.format("%Y-%m-%d %H:%M:%S"));
            }
        }
        
        // Check for reasonable orbital radius (should be ~7000 km for Sentinel-1)
        for sv in state_vectors {
            let position_magnitude = (sv.position[0].powi(2) + 
                                    sv.position[1].powi(2) + 
                                    sv.position[2].powi(2)).sqrt();
            
            if position_magnitude < 6_500_000.0 || position_magnitude > 7_500_000.0 {
                log::warn!("Unusual orbital radius: {:.1} km at {}", 
                          position_magnitude / 1000.0, sv.time.format("%Y-%m-%d %H:%M:%S"));
            }
        }
        
        log::debug!("Orbit data validation completed successfully");
        Ok(())
    }
    
    /// Enhanced state vector parsing with better error handling
    fn parse_state_vector_line_enhanced(line: &str) -> SarResult<Option<StateVector>> {
        let mut time_opt = None;
        let mut position = [0.0; 3];
        let mut velocity = [0.0; 3];
        let mut found_fields = 0;
        
        // Parse all key=value pairs
        for part in line.split_whitespace() {
            if let Some((key, value)) = part.split_once('=') {
                match key {
                    "UTC" => {
                        // Try multiple time formats
                        if let Ok(naive_dt) = NaiveDateTime::parse_from_str(value, "%Y-%m-%dT%H:%M:%S%.f") {
                            time_opt = Some(DateTime::from_naive_utc_and_offset(naive_dt, Utc));
                            found_fields += 1;
                        } else {
                            log::warn!("Could not parse UTC time: {}", value);
                        }
                    },
                    "X" => {
                        position[0] = value.parse()
                            .map_err(|e| SarError::Processing(format!("Invalid X coordinate: {} ({})", value, e)))?;
                        found_fields += 1;
                    },
                    "Y" => {
                        position[1] = value.parse()
                            .map_err(|e| SarError::Processing(format!("Invalid Y coordinate: {} ({})", value, e)))?;
                        found_fields += 1;
                    },
                    "Z" => {
                        position[2] = value.parse()
                            .map_err(|e| SarError::Processing(format!("Invalid Z coordinate: {} ({})", value, e)))?;
                        found_fields += 1;
                    },
                    "VX" => {
                        velocity[0] = value.parse()
                            .map_err(|e| SarError::Processing(format!("Invalid VX velocity: {} ({})", value, e)))?;
                        found_fields += 1;
                    },
                    "VY" => {
                        velocity[1] = value.parse()
                            .map_err(|e| SarError::Processing(format!("Invalid VY velocity: {} ({})", value, e)))?;
                        found_fields += 1;
                    },
                    "VZ" => {
                        velocity[2] = value.parse()
                            .map_err(|e| SarError::Processing(format!("Invalid VZ velocity: {} ({})", value, e)))?;
                        found_fields += 1;
                    },
                    _ => {} // Ignore other fields like TAI, UT1, GPS
                }
            }
        }
        
        // Require minimum fields: time + position + velocity
        if found_fields < 7 {
            return Ok(None);
        }
        
        let time = time_opt.ok_or_else(|| {
            SarError::Processing("Missing UTC time in state vector".to_string())
        })?;
        
        Ok(Some(StateVector {
            time,
            position,
            velocity,
        }))
    }
}

/// OrbitCache for managing locally stored orbit files
pub struct OrbitCache {
    cache_dir: PathBuf,
}

impl OrbitCache {
    pub fn new(cache_dir: PathBuf) -> Self {
        Self { cache_dir }
    }
    
    pub fn get_cached_orbit(&self, product_id: &str, orbit_type: OrbitType) -> Option<PathBuf> {
        let filename = format!("{}_{}.EOF", product_id, orbit_type);
        let path = self.cache_dir.join(filename);
        if path.exists() { Some(path) } else { None }
    }
    
    pub fn cache_orbit(&self, product_id: &str, orbit_type: OrbitType, content: &str) -> SarResult<PathBuf> {
        std::fs::create_dir_all(&self.cache_dir)?;
        let filename = format!("{}_{}.EOF", product_id, orbit_type);
        let path = self.cache_dir.join(filename);
        std::fs::write(&path, content)?;
        Ok(path)
    }
}

/// OrbitManager provides high-level orbit file management
pub struct OrbitManager {
    cache: OrbitCache,
}

impl OrbitManager {
    pub fn new(cache_dir: PathBuf) -> Self {
        Self {
            cache: OrbitCache::new(cache_dir),
        }
    }

    pub fn get_orbit_data(&self, product_id: &str, start_time: DateTime<Utc>) -> SarResult<OrbitData> {
        // Try cache first, then download
        let orbit_type = OrbitReader::determine_orbit_type(start_time);
        
        if let Some(cached_path) = self.cache.get_cached_orbit(product_id, orbit_type) {
            log::info!("Using cached orbit file: {}", cached_path.display());
            return OrbitReader::read_orbit_file(cached_path);
        }
        
        // Download and cache
        log::info!("Downloading orbit file for {}", product_id);
        let orbit_data = OrbitReader::download_orbit_file(product_id, start_time, orbit_type, None)?;
        
        // Note: For now we skip caching and download each time
        // TODO: Implement proper text-based caching
        
        Ok(orbit_data)
    }
    
    pub fn has_orbit_cached(&self, product_id: &str, orbit_type: OrbitType) -> bool {
        self.cache.get_cached_orbit(product_id, orbit_type).is_some()
    }
    
    pub fn get_cache_dir(&self) -> &Path {
        &self.cache.cache_dir
    }
    
    pub fn download_and_cache_orbit_public(
        &self, 
        product_id: &str, 
        start_time: DateTime<Utc>, 
        orbit_type: OrbitType
    ) -> SarResult<PathBuf> {
        let orbit_data = OrbitReader::download_orbit_file(product_id, start_time, orbit_type, None)?;
        
        // Create a simple text representation for caching
        let mut content = String::new();
        for sv in &orbit_data.state_vectors {
            content.push_str(&format!(
                "UTC={} X={} Y={} Z={} VX={} VY={} VZ={}\n",
                sv.time.format("%Y-%m-%dT%H:%M:%S%.6f"),
                sv.position[0], sv.position[1], sv.position[2],
                sv.velocity[0], sv.velocity[1], sv.velocity[2]
            ));
        }
        
        self.cache.cache_orbit(product_id, orbit_type, &content)
    }
}
