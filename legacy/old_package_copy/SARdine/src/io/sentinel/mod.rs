#![allow(dead_code, unused_variables)]
use crate::constants::physical::SPEED_OF_LIGHT_M_S;
use crate::core::calibration::CalibrationCoefficients;
use crate::io::orbit::OrbitManager;
use crate::types::{
    AcquisitionMode, BurstOrbitData, BurstRecord, CoordinateSystem, DopplerCentroidPolynomial,
    GeoTransform, OrbitData, OrbitStatus, Polarization, SarComplex, SarError, SarImage,
    SarMetadata, SarResult, SubSwath,
};
mod bbox;
mod calibration;
mod deburst;
mod noise;
pub mod slc_reader;
use chrono::{DateTime, NaiveDateTime, Utc};
use ndarray::Array2;
use quick_xml::de::from_str;
use serde::Deserialize;
use std::collections::HashMap;
use std::fs::File;
use std::io::Read;
use std::path::{Path, PathBuf};
use std::sync::Arc;
use zip::ZipArchive;

/// Simplified XML structures for Sentinel-1 annotation parsing
#[derive(Debug, Deserialize)]
struct Product {
    #[serde(rename = "adsHeader")]
    #[allow(dead_code)]
    ads_header: AdsHeader,
    #[serde(rename = "generalAnnotation")]
    #[allow(dead_code)]
    general_annotation: GeneralAnnotation,
    #[serde(rename = "imageAnnotation")]
    #[allow(dead_code)]
    image_annotation: ImageAnnotation,
}

#[derive(Debug, Deserialize)]
struct AdsHeader {
    #[serde(rename = "missionId")]
    #[allow(dead_code)]
    mission_id: String,
    #[serde(rename = "productType")]
    #[allow(dead_code)]
    product_type: String,
    #[serde(rename = "startTime")]
    #[allow(dead_code)]
    start_time: String,
    #[serde(rename = "stopTime")]
    #[allow(dead_code)]
    stop_time: String,
}

#[derive(Debug, Deserialize)]
struct GeneralAnnotation {
    #[serde(rename = "productInformation")]
    #[allow(dead_code)]
    product_information: ProductInformation,
}

#[derive(Debug, Deserialize)]
struct ProductInformation {
    #[serde(rename = "missionDataTakeId")]
    #[allow(dead_code)]
    mission_data_take_id: String,
    #[serde(rename = "productComposition")]
    #[allow(dead_code)]
    product_composition: String,
}

#[derive(Debug, Deserialize)]
struct ImageAnnotation {
    #[serde(rename = "imageInformation")]
    #[allow(dead_code)]
    image_information: ImageInformation,
}

#[derive(Debug, Deserialize)]
struct ImageInformation {
    #[serde(rename = "numberOfSamples")]
    #[allow(dead_code)]
    number_of_samples: usize,
    #[serde(rename = "numberOfLines")]
    #[allow(dead_code)]
    number_of_lines: usize,
    #[serde(rename = "slantRangeTime")]
    #[allow(dead_code)]
    slant_range_time: f64,
    #[serde(rename = "rangePixelSpacing")]
    #[allow(dead_code)]
    range_pixel_spacing: f64,
    #[serde(rename = "azimuthPixelSpacing")]
    #[allow(dead_code)]
    azimuth_pixel_spacing: f64,
}

/// Supported Sentinel-1 product formats
#[derive(Debug, Clone)]
pub enum ProductFormat {
    Zip,
    Safe,
}

/// Cache performance statistics for monitoring caching effectiveness
#[derive(Debug, Clone)]
pub struct CacheStats {
    pub annotations_cached: usize,
    pub xml_files_cached: usize,
    pub calibration_files_cached: usize,
    pub total_memory_usage_bytes: usize,
    pub cache_initialized: bool,
}

/// Sentinel-1 SLC reader that supports both ZIP and SAFE formats
/// Enhanced with comprehensive caching for optimal performance
pub struct SlcReader {
    product_path: PathBuf,
    format: ProductFormat,
    archive: Option<ZipArchive<File>>,
    orbit_data: Option<OrbitData>,

    /// Legacy calibration cache
    calibration_cache: std::collections::HashMap<String, CalibrationCoefficients>,

    /// NEW: Comprehensive caching system for optimal performance
    /// Cached metadata extracted once during initialization
    cached_metadata: Option<SarMetadata>,

    /// Cached annotations per polarization to eliminate XML re-parsing
    cached_annotations: HashMap<Polarization, Vec<Arc<crate::io::annotation::ProductRoot>>>,

    /// Cached calibration coefficients per polarization  
    cached_calibration: HashMap<Polarization, CalibrationCoefficients>,

    /// Cached raw XML content to eliminate file I/O
    cached_xml_content: HashMap<String, String>,

    /// Cache initialization status
    cache_initialized: bool,
}

impl SlcReader {
    /// Return the path to the underlying Sentinel-1 product
    pub fn product_path(&self) -> PathBuf {
        self.product_path.clone()
    }

    /// Internal constructor that performs format detection and wiring but does NOT warm caches
    fn new_internal<P: AsRef<Path>>(product_path: P) -> SarResult<Self> {
        let product_path = product_path.as_ref().to_path_buf();

        if !product_path.exists() {
            return Err(SarError::Io(std::io::Error::new(
                std::io::ErrorKind::NotFound,
                format!("Product not found: {}", product_path.display()),
            )));
        }

        // Determine format based on file extension and structure
        let format = if product_path.is_file()
            && product_path.extension() == Some(std::ffi::OsStr::new("zip"))
        {
            ProductFormat::Zip
        } else if product_path.is_dir()
            && (product_path
                .file_name()
                .map(|name| name.to_string_lossy().contains(".SAFE"))
                .unwrap_or_else(|| {
                    log::debug!(
                        "Could not determine filename for SAFE detection - checking manifest.safe"
                    );
                    false
                })
                || product_path.join("manifest.safe").exists())
        {
            ProductFormat::Safe
        } else {
            return Err(SarError::InvalidFormat(format!(
                "Unsupported format: {}. Must be .zip file or .SAFE directory",
                product_path.display()
            )));
        };

        log::info!(
            "Detected Sentinel-1 product format: {:?} at {}",
            format,
            product_path.display()
        );

        Ok(Self {
            product_path,
            format,
            archive: None,
            orbit_data: None,
            calibration_cache: std::collections::HashMap::new(),

            // Initialize new caching system
            cached_metadata: None,
            cached_annotations: HashMap::new(),
            cached_calibration: HashMap::new(),
            cached_xml_content: HashMap::new(),
            cache_initialized: false,
        })
    }

    /// Create a new SLC reader for a Sentinel-1 product (ZIP or SAFE) with full cache warm-up
    pub fn new<P: AsRef<Path>>(product_path: P) -> SarResult<Self> {
        Self::new_with_full_cache(product_path)
    }

    // ============================================================================
    // PATH CHECKING HELPERS - Cross-platform and case-insensitive
    // ============================================================================

    /// Check if a file path is an annotation XML (robust)
    fn is_annotation_xml(file: &str) -> bool {
        use std::path::Path;
        let p = Path::new(file);
        let lf = file.to_ascii_lowercase();

        p.components().any(|c| c.as_os_str() == "annotation")
            && p.extension()
                .and_then(|e| e.to_str())
                .map(|s| s.eq_ignore_ascii_case("xml"))
                .unwrap_or(false)
            && !lf.contains("calibration")
            && !lf.contains("noise")
    }

    /// Check if a file path is a measurement TIFF (robust)
    fn is_measurement_tiff(file: &str) -> bool {
        use std::path::Path;
        let p = Path::new(file);

        p.components().any(|c| c.as_os_str() == "measurement")
            && p.extension()
                .and_then(|e| e.to_str())
                .map(|s| s.eq_ignore_ascii_case("tiff"))
                .unwrap_or(false)
    }

    /// Check if a file path is a calibration XML (robust)
    fn is_calibration_xml(file: &str) -> bool {
        use std::path::Path;
        let p = Path::new(file);
        let lf = file.to_ascii_lowercase();

        // Check for annotation/calibration path using components
        let has_calibration_path = {
            let components: Vec<_> = p
                .components()
                .map(|c| c.as_os_str().to_string_lossy().to_lowercase())
                .collect();
            components
                .windows(2)
                .any(|w| w[0] == "annotation" && w[1] == "calibration")
        };

        has_calibration_path
            && p.extension()
                .and_then(|e| e.to_str())
                .map(|s| s.eq_ignore_ascii_case("xml"))
                .unwrap_or(false)
            && lf.contains("calibration-")
            && !lf.contains("noise-")
    }

    /// Check if a file path is a noise XML (robust)
    fn is_noise_xml(file: &str) -> bool {
        use std::path::Path;
        let p = Path::new(file);
        let lf = file.to_ascii_lowercase();

        // Check for annotation/calibration path using components
        let has_calibration_path = {
            let components: Vec<_> = p
                .components()
                .map(|c| c.as_os_str().to_string_lossy().to_lowercase())
                .collect();
            components
                .windows(2)
                .any(|w| w[0] == "annotation" && w[1] == "calibration")
        };

        has_calibration_path
            && p.extension()
                .and_then(|e| e.to_str())
                .map(|s| s.eq_ignore_ascii_case("xml"))
                .unwrap_or(false)
            && lf.contains("noise-")
    }

    /// Extract polarization from filename (robust)
    fn extract_polarization(file: &str) -> Option<Polarization> {
        let lf = file.to_ascii_lowercase();
        if lf.contains("-vv-") || lf.contains("_vv_") || lf.ends_with("vv.tiff") {
            Some(Polarization::VV)
        } else if lf.contains("-vh-") || lf.contains("_vh_") || lf.ends_with("vh.tiff") {
            Some(Polarization::VH)
        } else if lf.contains("-hv-") || lf.contains("_hv_") || lf.ends_with("hv.tiff") {
            Some(Polarization::HV)
        } else if lf.contains("-hh-") || lf.contains("_hh_") || lf.ends_with("hh.tiff") {
            Some(Polarization::HH)
        } else {
            None
        }
    }

    /// Extract subswath from filename (e.g., "iw1", "iw2", "iw3" -> "IW1", "IW2", "IW3")
    /// Handles Sentinel-1 IW TOPSAR naming convention: s1b-iw1-slc-vv-...tiff
    fn extract_subswath(file: &str) -> Option<String> {
        let lf = file.to_ascii_lowercase();
        // Pattern: -iw1-, -iw2-, -iw3- (Sentinel-1 IW mode)
        if lf.contains("-iw1-") || lf.contains("_iw1_") {
            Some("IW1".to_string())
        } else if lf.contains("-iw2-") || lf.contains("_iw2_") {
            Some("IW2".to_string())
        } else if lf.contains("-iw3-") || lf.contains("_iw3_") {
            Some("IW3".to_string())
        // Also support EW mode subswaths
        } else if lf.contains("-ew1-") || lf.contains("_ew1_") {
            Some("EW1".to_string())
        } else if lf.contains("-ew2-") || lf.contains("_ew2_") {
            Some("EW2".to_string())
        } else if lf.contains("-ew3-") || lf.contains("_ew3_") {
            Some("EW3".to_string())
        } else if lf.contains("-ew4-") || lf.contains("_ew4_") {
            Some("EW4".to_string())
        } else if lf.contains("-ew5-") || lf.contains("_ew5_") {
            Some("EW5".to_string())
        } else {
            None
        }
    }

    // ============================================================================
    // ARCHIVE AND FILE ACCESS
    // ============================================================================

    /// Open the ZIP archive (only for ZIP format)
    fn open_archive(&mut self) -> SarResult<&mut ZipArchive<File>> {
        match self.format {
            ProductFormat::Zip => {
                if self.archive.is_none() {
                    let file = File::open(&self.product_path)?;
                    let archive = ZipArchive::new(file).map_err(|e| {
                        SarError::InvalidFormat(format!("Failed to open ZIP: {}", e))
                    })?;
                    self.archive = Some(archive);
                }
                Ok(self
                    .archive
                    .as_mut()
                    .expect("Archive should exist after successful creation"))
            }
            ProductFormat::Safe => Err(SarError::InvalidFormat(
                "Cannot open ZIP archive on SAFE directory".to_string(),
            )),
        }
    }

    /// List all files in the product (works for both ZIP and SAFE)
    pub fn list_files(&mut self) -> SarResult<Vec<String>> {
        match self.format {
            ProductFormat::Zip => {
                let archive = self.open_archive()?;
                let mut files = Vec::new();

                for i in 0..archive.len() {
                    let file = archive.by_index(i).map_err(|e| {
                        SarError::Io(std::io::Error::new(
                            std::io::ErrorKind::Other,
                            format!("Failed to access file {}: {}", i, e),
                        ))
                    })?;
                    files.push(file.name().to_string());
                }

                Ok(files)
            }
            ProductFormat::Safe => {
                let mut files = Vec::new();
                Self::list_safe_files_recursive(&self.product_path, "", &mut files)?;
                Ok(files)
            }
        }
    }

    /// Recursively list files in SAFE directory structure
    fn list_safe_files_recursive(
        dir: &Path,
        prefix: &str,
        files: &mut Vec<String>,
    ) -> SarResult<()> {
        let entries = std::fs::read_dir(dir).map_err(|e| SarError::Io(e))?;

        for entry in entries {
            let entry = entry.map_err(|e| SarError::Io(e))?;
            let path = entry.path();
            let file_name = entry.file_name();
            let file_name_str = file_name.to_string_lossy();

            let relative_path = if prefix.is_empty() {
                file_name_str.to_string()
            } else {
                format!("{}/{}", prefix, file_name_str)
            };

            if path.is_dir() {
                Self::list_safe_files_recursive(&path, &relative_path, files)?;
            } else {
                files.push(relative_path);
            }
        }

        Ok(())
    }

    /// Read file content from product (works for both ZIP and SAFE)
    pub fn read_file(&mut self, file_path: &str) -> SarResult<Vec<u8>> {
        match self.format {
            ProductFormat::Zip => {
                let archive = self.open_archive()?;
                let mut file = archive.by_name(file_path).map_err(|e| {
                    SarError::Io(std::io::Error::new(
                        std::io::ErrorKind::Other,
                        format!("Failed to read {}: {}", file_path, e),
                    ))
                })?;

                let mut content = Vec::with_capacity(file.size() as usize);
                file.read_to_end(&mut content)?;
                Ok(content)
            }
            ProductFormat::Safe => {
                let full_path = self.product_path.join(file_path);
                let meta = std::fs::metadata(&full_path).ok();
                let mut buf = Vec::with_capacity(meta.map(|m| m.len() as usize).unwrap_or(0));
                let mut f = std::fs::File::open(&full_path)?;
                f.read_to_end(&mut buf)?;
                Ok(buf)
            }
        }
    }

    /// Read file content as string from product
    pub fn read_file_as_string(&mut self, file_path: &str) -> SarResult<String> {
        let content = self.read_file(file_path)?;
        String::from_utf8(content)
            .map_err(|e| SarError::InvalidFormat(format!("Invalid UTF-8 in {}: {}", file_path, e)))
    }

    /// Set orbit data for the reader
    pub fn set_orbit_data(&mut self, orbit_data: OrbitData) {
        self.orbit_data = Some(orbit_data);
    }

    /// Find annotation files for each polarization - FIXED to properly collect all files  
    pub fn find_annotation_files(&mut self) -> SarResult<HashMap<Polarization, String>> {
        let files = self.list_files()?;
        let mut all_annotations: HashMap<Polarization, Vec<String>> = HashMap::new();

        // First, collect ALL annotation files per polarization (don't overwrite)
        for file in files {
            if Self::is_annotation_xml(&file) {
                // Parse polarization from filename - Sentinel-1 standard naming
                if let Some(pol) = Self::extract_polarization(&file) {
                    // Add to the list for this polarization (no overwriting!)
                    all_annotations
                        .entry(pol)
                        .or_insert_with(Vec::new)
                        .push(file);
                }
            }
        }

        if all_annotations.is_empty() {
            return Err(SarError::InvalidFormat(
                "No annotation files found. Expected files like 's1a-iw-slc-vv-*.xml' in annotation/ directory".to_string(),
            ));
        }

        // For backward compatibility, return just the first file per polarization
        // Logging is toned down because full subswath handling happens elsewhere
        let mut result = HashMap::new();
        for (pol, files) in &all_annotations {
            if let Some(first_file) = files.first() {
                result.insert(*pol, first_file.clone());
                if files.len() > 1 {
                    log::debug!(
                        "Found {} annotation files for {}. Using first for compatibility APIs: {}",
                        files.len(),
                        pol,
                        first_file
                    );
                    for (i, file) in files.iter().enumerate() {
                        log::debug!("  Candidate {}: {}", i + 1, file);
                    }
                } else {
                    log::info!("Found 1 annotation file for {}: {}", pol, first_file);
                }
            }
        }

        // Note: Full subswath discovery is handled by find_all_annotation_files()/cache init.
        // Avoid alarming warnings here to prevent log spam in strict workflows.
        let total_files: usize = all_annotations.values().map(|v| v.len()).sum();
        log::debug!(
            "Annotation overview: {} total files across {} polarizations; representative selection returned {}",
            total_files,
            all_annotations.len(),
            result.len()
        );

        Ok(result)
    }

    /// Find ALL annotation files for each polarization (all subswaths)
    pub fn find_all_annotation_files(&mut self) -> SarResult<HashMap<Polarization, Vec<String>>> {
        let files = self.list_files()?;
        let mut annotations: HashMap<Polarization, Vec<String>> = HashMap::new();

        for file in files {
            // Only include main annotation files, exclude calibration, noise, and RFI subdirectories
            if Self::is_annotation_xml(&file) && !file.contains("rfi") {
                // Additional check: only include files directly in annotation/ directory, not subdirectories
                use std::path::Path;
                let path = Path::new(&file);
                // Check if there are subdirectories after 'annotation' component
                let mut found_annotation = false;
                let mut has_subdir = false;
                for comp in path.components() {
                    if found_annotation && comp.as_os_str() != path.file_name().unwrap_or_default()
                    {
                        has_subdir = true;
                        break;
                    }
                    if comp
                        .as_os_str()
                        .to_str()
                        .map(|s| s.to_lowercase() == "annotation")
                        .unwrap_or(false)
                    {
                        found_annotation = true;
                    }
                }
                if has_subdir {
                    // Skip files in subdirectories (like calibration/, noise/, rfi/)
                    continue;
                }

                // Parse polarization from filename - Sentinel-1 standard naming
                if let Some(pol) = Self::extract_polarization(&file) {
                    // Add to the list for this polarization
                    annotations.entry(pol).or_insert_with(Vec::new).push(file);
                }
            }
        }

        if annotations.is_empty() {
            return Err(SarError::InvalidFormat(
                "No annotation files found. Expected files like 's1a-iw-slc-vv-*.xml' in annotation/ directory".to_string(),
            ));
        }

        // Validate per-polarization IW coverage and ordering early to fail fast on malformed SAFE structures
        Self::validate_annotation_coverage(&annotations)?;

        let total_files: usize = annotations.values().map(|v| v.len()).sum();
        log::info!(
            "Found {} total annotation files across {} polarizations",
            total_files,
            annotations.len()
        );
        for (pol, files) in &annotations {
            log::info!("  {}: {} files", pol, files.len());
        }

        Ok(annotations)
    }

    fn validate_annotation_coverage(
        annotations: &HashMap<Polarization, Vec<String>>,
    ) -> SarResult<()> {
        const EXPECTED_IW: [&str; 3] = ["IW1", "IW2", "IW3"];

        for (pol, files) in annotations {
            let mut swaths: Vec<String> = files
                .iter()
                .filter_map(|f| Self::extract_subswath(f))
                .collect();
            swaths.sort();
            swaths.dedup();

            let iw_swaths: Vec<&String> = swaths.iter().filter(|s| s.starts_with("IW")).collect();

            if !iw_swaths.is_empty() {
                for expected in EXPECTED_IW {
                    if !iw_swaths.iter().any(|s| s.as_str() == expected) {
                        return Err(SarError::Processing(format!(
                            "Missing {} annotation for {:?} (found {:?})",
                            expected, pol, swaths
                        )));
                    }
                }

                let mut last_idx = 0;
                for s in iw_swaths {
                    let idx = match s.as_str() {
                        "IW1" => 1,
                        "IW2" => 2,
                        "IW3" => 3,
                        _ => continue,
                    };
                    if idx < last_idx {
                        return Err(SarError::Processing(format!(
                            "IW ordering invalid for {:?}: {:?}",
                            pol, swaths
                        )));
                    }
                    last_idx = idx;
                }
            }
        }

        Ok(())
    }

    /// Read and parse annotation XML for a polarization, merging across all subswaths (IW1/2/3)
    pub fn read_annotation(&mut self, pol: Polarization) -> SarResult<SarMetadata> {
        let product_id = Self::product_id_from_path(&self.product_path);

        if !self.cache_initialized {
            self.initialize_all_caches()?;
        }

        if let Some(cached_roots) = self.cached_annotations.get(&pol) {
            let metas = cached_roots
                .iter()
                .map(|root| Self::metadata_from_annotation_root(root.as_ref(), pol, &product_id))
                .collect::<SarResult<Vec<_>>>()?;
            return Self::merge_metadata(metas);
        }

        let all = self.find_all_annotation_files()?;
        let files = all.get(&pol).ok_or_else(|| {
            SarError::InvalidFormat(format!(
                "No annotation files found for polarization {}",
                pol
            ))
        })?;

        if files.is_empty() {
            return Err(SarError::InvalidFormat(format!(
                "Empty annotation file list for polarization {}",
                pol
            )));
        }

        log::info!(
            "Reading {} annotation file(s) for {} and merging subswaths",
            files.len(),
            pol
        );

        let mut metas: Vec<SarMetadata> = Vec::with_capacity(files.len());
        for f in files {
            let xml = self.read_file_as_string(f)?;
            let annotation_root = self.parse_annotation_root_xml(&xml)?;
            let meta = Self::metadata_from_annotation_root(&annotation_root, pol, &product_id)?;
            metas.push(meta);
        }

        Self::merge_metadata(metas)
    }

    fn extract_real_value_f64(xml: &str, tag: &str) -> Option<f64> {
        Self::extract_value(xml, tag).and_then(|s| {
            let cleaned = s
                .trim()
                .replace(" m", "")
                .replace(" Hz", "")
                .replace(" s", "")
                .replace(" deg", "");

            match cleaned.parse::<f64>() {
                Ok(val) => {
                    log::debug!("Extracted {} = {} from XML", tag, val);
                    Some(val)
                }
                Err(e) => {
                    log::warn!("Failed to parse {} value '{}': {}", tag, cleaned, e);
                    None
                }
            }
        })
    }

    /// Flexible time parsing that handles multiple formats
    fn parse_time_flexible(time_str: &str) -> SarResult<DateTime<Utc>> {
        let s = time_str.trim();
        let post_t = s.split('T').nth(1).unwrap_or("");
        let has_tz = post_t.contains('Z') || post_t.contains('+') || post_t.contains('-');

        // Prefer handling naive timestamps first by assuming UTC (append Z)
        if !has_tz {
            if let Ok(dt) = DateTime::parse_from_str(&format!("{}Z", s), "%Y-%m-%dT%H:%M:%S%.fZ") {
                return Ok(dt.with_timezone(&Utc));
            }
            if let Ok(naive_dt) = chrono::NaiveDateTime::parse_from_str(s, "%Y-%m-%dT%H:%M:%S%.f") {
                return Ok(DateTime::<Utc>::from_naive_utc_and_offset(naive_dt, Utc));
            }
            if let Ok(naive_dt) = chrono::NaiveDateTime::parse_from_str(s, "%Y-%m-%dT%H:%M:%S") {
                return Ok(DateTime::<Utc>::from_naive_utc_and_offset(naive_dt, Utc));
            }
        }

        // RFC3339 with timezone (Z or offsets)
        if let Ok(time) = DateTime::parse_from_rfc3339(s) {
            return Ok(time.with_timezone(&Utc));
        }
        // Explicit UTC offset with fractional seconds
        if let Ok(time) = DateTime::parse_from_str(s, "%Y-%m-%dT%H:%M:%S%.f+00:00") {
            return Ok(time.with_timezone(&Utc));
        }
        // With Z but no fractional seconds
        if let Ok(time) = DateTime::parse_from_str(&format!("{}Z", s), "%Y-%m-%dT%H:%M:%SZ") {
            return Ok(time.with_timezone(&Utc));
        }

        Err(SarError::InvalidFormat(
            format!("Could not parse time '{}' - no valid time format recognized. Real Sentinel-1 annotation required.", time_str)
        ))
    }

    /// Simple XML value extraction helper
    fn extract_value(xml: &str, tag: &str) -> Option<String> {
        let start_tag = format!("<{}>", tag);
        let end_tag = format!("</{}>", tag);

        if let Some(start) = xml.find(&start_tag) {
            let content_start = start + start_tag.len();
            if let Some(end) = xml[content_start..].find(&end_tag) {
                return Some(xml[content_start..content_start + end].to_string());
            }
        }
        None
    }

    /// Extract specific value with unit handling
    #[allow(dead_code)]
    fn extract_value_with_unit(xml: &str, tag: &str) -> Option<String> {
        if let Some(value) = Self::extract_value(xml, tag) {
            // Remove common units and clean the value
            let cleaned = value
                .replace(" Hz", "")
                .replace(" m", "")
                .replace(" deg", "")
                .replace(" s", "")
                .trim()
                .to_string();
            Some(cleaned)
        } else {
            None
        }
    }

    /// Get annotation data for a specific polarization
    pub fn get_annotation_for_polarization(
        &mut self,
        pol: Polarization,
    ) -> SarResult<Arc<crate::io::annotation::ProductRoot>> {
        if !self.cache_initialized {
            self.initialize_all_caches()?;
        }

        self.cached_annotations
            .get(&pol)
            .and_then(|annotations| annotations.first().cloned())
            .ok_or_else(|| {
                SarError::Processing(format!(
                    "Annotation not cached for {:?}. Available: {:?}",
                    pol,
                    self.cached_annotations.keys().collect::<Vec<_>>()
                ))
            })
    }

    /// Retrieve all cached subswath annotations for a polarization without re-parsing XML
    pub fn get_all_subswath_annotations(
        &mut self,
        pol: Polarization,
    ) -> SarResult<Vec<Arc<crate::io::annotation::ProductRoot>>> {
        if !self.cache_initialized {
            self.initialize_all_caches()?;
        }

        self.get_all_cached_annotations(pol)
    }

    /// **MIGRATION HELPER**: Get primary annotation with subswath coverage awareness
    ///
    /// This method provides a migration path from single-annotation processing to
    /// multi-subswath awareness. It returns the first annotation but logs warnings
    /// about any additional subswaths that would be missed.
    pub fn get_primary_annotation_with_coverage_check(
        &mut self,
        pol: Polarization,
    ) -> SarResult<Arc<crate::io::annotation::ProductRoot>> {
        let all_annotations = self.get_all_subswath_annotations(pol)?;

        if all_annotations.is_empty() {
            return Err(SarError::Processing(format!(
                "No annotations found for polarization {:?}",
                pol
            )));
        }

        if all_annotations.len() > 1 {
            log::warn!(
                "🚨 SCIENTIFIC DATA LOSS WARNING: Using only 1 of {} available subswaths for {}",
                all_annotations.len(),
                pol
            );
            log::warn!("🚨 This may result in incomplete scientific analysis!");
            log::warn!("🚨 Consider migrating to get_all_subswath_annotations() for full coverage");
        }

        Ok(all_annotations
            .into_iter()
            .next()
            .expect("Vector guaranteed non-empty after validation"))
    }

    /// **PHASE 2A: UNIFIED PARSING INFRASTRUCTURE**
    ///
    /// This method demonstrates the unified parsing approach that consolidates
    /// regex-based and serde-based XML parsing for scientific consistency.
    ///
    /// Features:
    /// - Validation between parsing approaches
    /// - Graceful fallback with scientific warnings
    /// - Migration path for existing code
    ///
    /// # Arguments
    /// * `pol` - Polarization to parse
    /// * `validate_consistency` - Whether to validate both parsing methods
    ///
    /// # Returns
    /// * `Ok(ProductRoot)` - Parsed annotation with validation
    /// * `Err(SarError)` - Parsing failure or validation issues
    pub fn get_annotation_unified_validated(
        &mut self,
        pol: Polarization,
        validate_consistency: bool,
    ) -> SarResult<crate::io::annotation::ProductRoot> {
        // Get annotation files (all subswaths) and select the first for validation
        let annotations = self.find_all_annotation_files()?;
        let files = annotations.get(&pol).ok_or_else(|| {
            SarError::Processing(format!(
                "No annotation files found for polarization {:?}",
                pol
            ))
        })?;
        let first = files.first().ok_or_else(|| {
            SarError::Processing(format!("Empty annotation list for polarization {:?}", pol))
        })?;

        // Read XML content
        let xml_content = self.read_file_as_string(first)?;

        if validate_consistency {
            log::info!(
                "🔬 SCIENTIFIC VALIDATION: Comparing parsing methods for {}",
                pol
            );

            // Validate parsing consistency
            let validator = crate::io::parsing_validation::ParsingValidator::new();
            match validator.validate_parsing_equivalence(&xml_content) {
                Ok(validation_result) => {
                    if !validation_result.equivalent {
                        log::warn!("🚨 PARSING INCONSISTENCY DETECTED for {}", pol);
                        log::warn!("🚨 Details: {}", validation_result.details);
                    } else {
                        log::info!(
                            "✅ Parsing validation passed for {}: {}",
                            pol,
                            validation_result.details
                        );
                    }
                }
                Err(e) => {
                    log::warn!("⚠️  Parsing validation failed for {}: {}", pol, e);
                }
            }
        }

        // Primary parsing attempt - use serde (recommended)
        match crate::io::annotation::parse_annotation_xml(&xml_content) {
            Ok(product_root) => {
                log::info!("✅ Serde parsing successful for {}", pol);
                Ok(product_root)
            }
            Err(serde_err) => Err(SarError::Processing(format!(
                "Unified annotation parser failed for {}: {}",
                pol, serde_err
            ))),
        }
    }

    /// Find measurement data files (TIFF files containing SLC data)
    pub fn find_measurement_files(&mut self) -> SarResult<HashMap<Polarization, String>> {
        let files = self.list_files()?;
        let mut measurements = HashMap::new();

        for file in files {
            if Self::is_measurement_tiff(&file) {
                // Parse polarization from filename - handle both ZIP and SAFE naming
                if let Some(pol) = Self::extract_polarization(&file) {
                    measurements.insert(pol, file);
                }
            }
        }

        if measurements.is_empty() {
            return Err(SarError::InvalidFormat(
                "No measurement files found. Expected TIFF files in measurement/ directory"
                    .to_string(),
            ));
        }

        log::info!(
            "Found {} measurement files: {:?}",
            measurements.len(),
            measurements.keys().collect::<Vec<_>>()
        );

        Ok(measurements)
    }

    /// Find measurement file for a specific subswath and polarization
    ///
    /// For IW TOPSAR products, there are separate TIFF files per subswath:
    /// - s1b-iw1-slc-vv-...tiff (IW1, VV)
    /// - s1b-iw2-slc-vv-...tiff (IW2, VV)  
    /// - s1b-iw3-slc-vv-...tiff (IW3, VV)
    ///
    /// This method finds the correct file for a given subswath/polarization combination.
    pub fn find_measurement_file_for_subswath(
        &mut self,
        subswath: &str,
        pol: Polarization,
    ) -> SarResult<String> {
        let files = self.list_files()?;
        let subswath_upper = subswath.to_uppercase();

        for file in files {
            if Self::is_measurement_tiff(&file) {
                let file_pol = Self::extract_polarization(&file);
                let file_sw = Self::extract_subswath(&file);

                if file_pol == Some(pol) {
                    // Check if subswath matches (or if file has no subswath marker, accept it)
                    if let Some(sw) = file_sw {
                        if sw == subswath_upper {
                            log::info!("Found measurement file for {} {}: {}", subswath, pol, file);
                            return Ok(file);
                        }
                    }
                }
            }
        }

        // Fallback: if no subswath-specific file found, try generic (for non-TOPSAR products)
        let measurements = self.find_measurement_files()?;
        if let Some(generic_file) = measurements.get(&pol) {
            log::warn!(
                "No subswath-specific file found for {} {}, using generic: {}",
                subswath,
                pol,
                generic_file
            );
            return Ok(generic_file.clone());
        }

        Err(SarError::InvalidFormat(format!(
            "No measurement file found for subswath {} polarization {}",
            subswath, pol
        )))
    }

    fn measurement_raster_size_from_file(
        &mut self,
        measurement_file: &str,
    ) -> SarResult<(usize, usize)> {
        match self.format {
            ProductFormat::Safe => {
                let full_path = self.product_path.join(measurement_file);
                if !full_path.exists() {
                    return Err(SarError::Io(std::io::Error::new(
                        std::io::ErrorKind::NotFound,
                        format!("Measurement file not found: {:?}", full_path),
                    )));
                }

                let dataset = gdal::Dataset::open(&full_path).map_err(|e| {
                    SarError::Io(std::io::Error::new(
                        std::io::ErrorKind::Other,
                        format!("Failed to open TIFF with GDAL: {}", e),
                    ))
                })?;

                Ok(dataset.raster_size())
            }
            ProductFormat::Zip => {
                let vsizip_path = format!(
                    "/vsizip/{}/{}",
                    self.product_path.display(),
                    measurement_file
                );

                match gdal::Dataset::open(&vsizip_path) {
                    Ok(dataset) => Ok(dataset.raster_size()),
                    Err(err) => {
                        log::warn!(
                            "⚠️  VSIZIP open failed for {} ({}); falling back to temp extraction",
                            measurement_file,
                            err
                        );

                        use tempfile::NamedTempFile;

                        let (width, height) = {
                            let archive = self.open_archive()?;
                            let mut zip_file = archive.by_name(measurement_file).map_err(|e| {
                                SarError::Io(std::io::Error::new(
                                    std::io::ErrorKind::Other,
                                    format!("Failed to access {}: {}", measurement_file, e),
                                ))
                            })?;

                            let mut temp_file = NamedTempFile::new().map_err(|e| SarError::Io(e))?;
                            std::io::copy(&mut zip_file, &mut temp_file)
                                .map_err(|e| SarError::Io(e))?;

                            let dataset = gdal::Dataset::open(temp_file.path()).map_err(|e| {
                                SarError::Io(std::io::Error::new(
                                    std::io::ErrorKind::Other,
                                    format!("Failed to open TIFF with GDAL: {}", e),
                                ))
                            })?;

                            dataset.raster_size()
                        };

                        Ok((width, height))
                    }
                }
            }
        }
    }

    fn verify_measurement_geometry(
        &mut self,
        pol: Polarization,
        swath_id: &str,
        swath: &SubSwath,
    ) -> SarResult<()> {
        let measurement_file = self.find_measurement_file_for_subswath(swath_id, pol)?;
        let (width, height) = self.measurement_raster_size_from_file(&measurement_file)?;

        if swath.range_samples != width {
            return Err(SarError::Processing(format!(
                "Measurement raster width mismatch for {} {}: TIFF width {} vs annotation range_samples {}",
                swath_id, pol, width, swath.range_samples
            )));
        }

        if swath.azimuth_samples != height {
            return Err(SarError::Processing(format!(
                "Measurement raster height mismatch for {} {}: TIFF height {} vs annotation azimuth_samples {}",
                swath_id, pol, height, swath.azimuth_samples
            )));
        }

        let expected_height = swath
            .lines_per_burst
            .saturating_mul(swath.burst_count.max(1));
        if expected_height != height {
            return Err(SarError::Processing(format!(
                "Burst geometry invariant violated for {} {}: {} bursts × {} lines = {}, but raster height is {}",
                swath_id,
                pol,
                swath.burst_count,
                swath.lines_per_burst,
                expected_height,
                height
            )));
        }

        log::info!(
            "🔐 Verified measurement geometry for {} {}: {} x {}",
            swath_id, pol, width, height
        );

        Ok(())
    }

    /// Read SLC data for a specific subswath and polarization
    ///
    /// This is the correct method to use for IW TOPSAR products where each
    /// subswath has its own TIFF file.
    pub fn read_slc_data_for_subswath(
        &mut self,
        subswath: &str,
        pol: Polarization,
    ) -> SarResult<SarImage> {
        use tempfile::NamedTempFile;

        let measurement_file = self.find_measurement_file_for_subswath(subswath, pol)?;

        log::info!(
            "📂 Reading SLC data for subswath {} {} from {}",
            subswath,
            pol,
            measurement_file
        );
        let start_time = std::time::Instant::now();

        // Handle different formats (same logic as read_slc_data)
        let dataset = match self.format {
            ProductFormat::Zip => {
                let archive = self.open_archive()?;
                let mut zip_file = archive.by_name(&measurement_file).map_err(|e| {
                    SarError::Io(std::io::Error::new(
                        std::io::ErrorKind::Other,
                        format!("Failed to access {}: {}", measurement_file, e),
                    ))
                })?;

                let mut temp_file = NamedTempFile::new().map_err(|e| SarError::Io(e))?;
                std::io::copy(&mut zip_file, &mut temp_file).map_err(|e| SarError::Io(e))?;

                let extract_time = start_time.elapsed();
                log::debug!("TIFF extraction took: {:?}", extract_time);

                gdal::Dataset::open(temp_file.path()).map_err(|e| {
                    SarError::Io(std::io::Error::new(
                        std::io::ErrorKind::Other,
                        format!("Failed to open TIFF with GDAL: {}", e),
                    ))
                })?
            }
            ProductFormat::Safe => {
                let full_path = self.product_path.join(&measurement_file);
                if !full_path.exists() {
                    return Err(SarError::Io(std::io::Error::new(
                        std::io::ErrorKind::NotFound,
                        format!("Measurement file not found: {:?}", full_path),
                    )));
                }

                gdal::Dataset::open(&full_path).map_err(|e| {
                    SarError::Io(std::io::Error::new(
                        std::io::ErrorKind::Other,
                        format!("Failed to open TIFF with GDAL: {}", e),
                    ))
                })?
            }
        };

        // Get raster dimensions
        let (width, height) = dataset.raster_size();
        log::info!(
            "📐 SLC raster size for {} {}: {} x {} (width x height)",
            subswath,
            pol,
            width,
            height
        );

        // Read complex data using the same approach as read_slc_data
        let band = dataset.rasterband(1).map_err(|e| {
            SarError::Io(std::io::Error::new(
                std::io::ErrorKind::Other,
                format!("Failed to get rasterband: {}", e),
            ))
        })?;

        // Read CInt16 complex data in chunks to avoid large single allocations
        let array = Self::read_complex_band_chunked(&band, width, height)?;

        let elapsed = start_time.elapsed();
        log::info!(
            "✅ Read SLC for {} {} in {:?}: {} lines x {} samples",
            subswath,
            pol,
            elapsed,
            array.nrows(),
            array.ncols()
        );

        Ok(array)
    }

    /// Read SLC data for a specific polarization (supports both ZIP and SAFE)
    pub fn read_slc_data(&mut self, pol: Polarization) -> SarResult<SarImage> {
        use tempfile::NamedTempFile;

        let measurements = self.find_measurement_files()?;
        let measurement_file = measurements.get(&pol).ok_or_else(|| {
            SarError::InvalidFormat(format!("No measurement found for polarization {}", pol))
        })?;

        log::info!("Reading SLC data for {} from {}", pol, measurement_file);
        let start_time = std::time::Instant::now();

        // Handle different formats
        let dataset = match self.format {
            ProductFormat::Zip => {
                // Extract the TIFF file to a temporary location for ZIP format
                let archive = self.open_archive()?;
                let mut zip_file = archive.by_name(measurement_file).map_err(|e| {
                    SarError::Io(std::io::Error::new(
                        std::io::ErrorKind::Other,
                        format!("Failed to access {}: {}", measurement_file, e),
                    ))
                })?;

                // Create a temporary file to extract the TIFF
                let mut temp_file = NamedTempFile::new().map_err(|e| SarError::Io(e))?;

                // Copy data from ZIP to temporary file
                std::io::copy(&mut zip_file, &mut temp_file).map_err(|e| SarError::Io(e))?;

                let extract_time = start_time.elapsed();
                log::debug!("TIFF extraction took: {:?}", extract_time);

                // Read with GDAL from temporary file
                gdal::Dataset::open(temp_file.path()).map_err(|e| {
                    SarError::Io(std::io::Error::new(
                        std::io::ErrorKind::Other,
                        format!("Failed to open TIFF with GDAL: {}", e),
                    ))
                })?
            }
            ProductFormat::Safe => {
                // Read directly from SAFE directory structure
                let full_path = self.product_path.join(measurement_file);
                gdal::Dataset::open(full_path).map_err(|e| {
                    SarError::Io(std::io::Error::new(
                        std::io::ErrorKind::Other,
                        format!("Failed to open TIFF with GDAL: {}", e),
                    ))
                })?
            }
        };

        let gdal_start = std::time::Instant::now();

        // Get raster dimensions and band count
        let raster_size = dataset.raster_size();
        let (width, height) = (raster_size.0, raster_size.1);
        let band_count = dataset.raster_count();
        log::debug!(
            "TIFF dimensions: {} x {}, bands: {}",
            width,
            height,
            band_count
        );

        // Handle different band configurations
        let slc_data = if band_count == 1 {
            // Single band containing complex data (Sentinel-1 SLC format: CInt16)
            let band = dataset.rasterband(1).map_err(|e| {
                SarError::Io(std::io::Error::new(
                    std::io::ErrorKind::Other,
                    format!("Failed to get band 1: {}", e),
                ))
            })?;

            // **CRITICAL GDAL SAFETY CHECK**: Verify band data type
            let band_type = band.band_type();
            log::info!("🔬 GDAL band type detected: {:?}", band_type);

            // Validate GDAL band type for complex SAR data
            // Sentinel-1 SLC uses CInt16 (complex 16-bit integers)
            use gdal::raster::GdalDataType;
            match band_type {
                GdalDataType::Int16 | GdalDataType::Int32 => {
                    // CInt16/CInt32 - Complex integer types (real/imag interleaved)
                    log::info!(
                        "✅ GDAL band type {:?} - valid for complex SAR data",
                        band_type
                    );
                }
                GdalDataType::Float32 | GdalDataType::Float64 => {
                    // CFloat32/CFloat64 - Complex float types
                    log::info!(
                        "✅ GDAL band type {:?} - valid for complex SAR data",
                        band_type
                    );
                }
                GdalDataType::UInt16 | GdalDataType::UInt32 => {
                    // Unsigned types - may be magnitude-only data
                    log::warn!(
                        "⚠️ GDAL band type {:?} - may not contain phase information",
                        band_type
                    );
                }
                _ => {
                    log::warn!(
                        "⚠️ Unexpected GDAL band type {:?} - proceeding with complex read attempt",
                        band_type
                    );
                }
            }

            // Streamed complex read
            Self::read_complex_band_chunked(&band, width, height)?
        } else if band_count >= 2 {
            // Separate I and Q bands (less common for Sentinel-1)
            let band1 = dataset.rasterband(1).map_err(|e| {
                SarError::Io(std::io::Error::new(
                    std::io::ErrorKind::Other,
                    format!("Failed to get band 1: {}", e),
                ))
            })?;

            let band2 = dataset.rasterband(2).map_err(|e| {
                SarError::Io(std::io::Error::new(
                    std::io::ErrorKind::Other,
                    format!("Failed to get band 2: {}", e),
                ))
            })?;

            let window = (0, 0);
            let window_size = (width, height);
            let buffer_size = (width, height);

            let i_data = band1
                .read_as::<f32>(window, window_size, buffer_size, None)
                .map_err(|e| {
                    SarError::Io(std::io::Error::new(
                        std::io::ErrorKind::Other,
                        format!("Failed to read I data: {}", e),
                    ))
                })?;

            let q_data = band2
                .read_as::<f32>(window, window_size, buffer_size, None)
                .map_err(|e| {
                    SarError::Io(std::io::Error::new(
                        std::io::ErrorKind::Other,
                        format!("Failed to read Q data: {}", e),
                    ))
                })?;

            // Convert to complex SLC data
            Self::convert_to_complex_parallel(i_data, q_data, width, height)?
        } else {
            return Err(SarError::InvalidFormat(format!(
                "Unexpected number of bands in TIFF: {}",
                band_count
            )));
        };

        let gdal_time = gdal_start.elapsed();
        log::debug!("GDAL read took: {:?}", gdal_time);

        let total_time = start_time.elapsed();

        log::info!(
            "SLC data read complete for {}: {} x {} pixels",
            pol,
            width,
            height
        );
        log::info!("Performance timing:");
        log::info!("  - Total: {:?}", total_time);
        log::info!("  - GDAL read: {:?}", gdal_time);

        let mb_size = (width * height * 8) as f64 / (1024.0 * 1024.0); // 8 bytes per complex float
        let mb_per_sec = mb_size / total_time.as_secs_f64();
        log::info!("  - Data size: {:.1} MB", mb_size);
        log::info!("  - Throughput: {:.1} MB/s", mb_per_sec);

        Ok(slc_data)
    }

    /// Read SLC data with optimized parallel processing
    pub fn read_slc_data_parallel(&mut self, pol: Polarization) -> SarResult<SarImage> {
        use tempfile::NamedTempFile;

        let measurements = self.find_measurement_files()?;
        let measurement_file = measurements.get(&pol).ok_or_else(|| {
            SarError::InvalidFormat(format!("No measurement found for polarization {}", pol))
        })?;

        log::info!(
            "Reading SLC data (PARALLEL) for {} from {}",
            pol,
            measurement_file
        );
        let start_time = std::time::Instant::now();

        // Handle TIFF access based on format
        let (dataset, extract_time) = match self.format {
            ProductFormat::Zip => {
                // Extract the TIFF file to a temporary location
                let archive = self.open_archive()?;
                let mut zip_file = archive.by_name(measurement_file).map_err(|e| {
                    SarError::Io(std::io::Error::new(
                        std::io::ErrorKind::Other,
                        format!("Failed to access {}: {}", measurement_file, e),
                    ))
                })?;

                // Create a temporary file to extract the TIFF
                let mut temp_file = NamedTempFile::new().map_err(|e| SarError::Io(e))?;

                // Copy data from ZIP to temporary file
                std::io::copy(&mut zip_file, &mut temp_file).map_err(|e| SarError::Io(e))?;

                let extract_time = start_time.elapsed();
                log::debug!("TIFF extraction took: {:?}", extract_time);

                // Open with GDAL
                let dataset = gdal::Dataset::open(temp_file.path()).map_err(|e| {
                    SarError::Io(std::io::Error::new(
                        std::io::ErrorKind::Other,
                        format!("Failed to open TIFF with GDAL: {}", e),
                    ))
                })?;

                (dataset, extract_time)
            }
            ProductFormat::Safe => {
                // Open TIFF directly from SAFE directory
                let extract_time = start_time.elapsed(); // No extraction needed for SAFE
                log::debug!("SAFE direct access took: {:?}", extract_time);

                let full_path = self.product_path.join(measurement_file);
                let dataset = gdal::Dataset::open(full_path).map_err(|e| {
                    SarError::Io(std::io::Error::new(
                        std::io::ErrorKind::Other,
                        format!("Failed to open TIFF with GDAL: {}", e),
                    ))
                })?;

                (dataset, extract_time)
            }
        };

        // Read with GDAL
        let gdal_start = std::time::Instant::now();

        // Get raster dimensions and band count
        let raster_size = dataset.raster_size();
        let (width, height) = (raster_size.0, raster_size.1);
        let band_count = dataset.raster_count();
        log::debug!(
            "TIFF dimensions: {} x {}, bands: {}",
            width,
            height,
            band_count
        );

        // Handle different band configurations
        let slc_data = if band_count == 1 {
            // Single band containing complex data (Sentinel-1 SLC format: CInt16)
            let band = dataset.rasterband(1).map_err(|e| {
                SarError::Io(std::io::Error::new(
                    std::io::ErrorKind::Other,
                    format!("Failed to get band 1: {}", e),
                ))
            })?;

            // **CRITICAL GDAL SAFETY CHECK**: Verify band data type
            let band_type = band.band_type();
            log::debug!("🔬 GDAL band type detected: {:?}", band_type);

            // Note: Complex data validation - proceeding with complex read
            log::debug!(
                "🔬 Proceeding with complex CInt16 read for type: {:?}",
                band_type
            );

            let window = (0, 0);
            let window_size = (width, height);

            // **PROPER COMPLEX READ**: Read complex data directly without width*2 trick
            let complex_data = band
                .read_as::<i16>(window, window_size, (width, height), None)
                .map_err(|e| {
                    SarError::Io(std::io::Error::new(
                        std::io::ErrorKind::Other,
                        format!("Failed to read complex CInt16 data from GDAL: {}", e),
                    ))
                })?;

            // **SAFETY VALIDATION**: Verify data size for complex format
            let expected_samples = width * height * 2;
            if complex_data.data.len() != expected_samples {
                return Err(SarError::InvalidFormat(format!(
                    "Complex data size validation failed: expected {} samples, got {}",
                    expected_samples,
                    complex_data.data.len()
                )));
            }

            log::debug!(
                "Read {} complex values as interleaved i16 (parallel)",
                complex_data.data.len() / 2
            );

            // Convert interleaved i16 data to complex f32 array with parallel processing
            Self::convert_cint16_to_complex_parallel(complex_data, width, height)?
        } else if band_count >= 2 {
            // Separate I and Q bands (less common for Sentinel-1)
            let band1 = dataset.rasterband(1).map_err(|e| {
                SarError::Io(std::io::Error::new(
                    std::io::ErrorKind::Other,
                    format!("Failed to get band 1: {}", e),
                ))
            })?;

            let band2 = dataset.rasterband(2).map_err(|e| {
                SarError::Io(std::io::Error::new(
                    std::io::ErrorKind::Other,
                    format!("Failed to get band 2: {}", e),
                ))
            })?;

            let window = (0, 0);
            let window_size = (width, height);
            let buffer_size = (width, height);

            let i_data = band1
                .read_as::<f32>(window, window_size, buffer_size, None)
                .map_err(|e| {
                    SarError::Io(std::io::Error::new(
                        std::io::ErrorKind::Other,
                        format!("Failed to read I data: {}", e),
                    ))
                })?;

            let q_data = band2
                .read_as::<f32>(window, window_size, buffer_size, None)
                .map_err(|e| {
                    SarError::Io(std::io::Error::new(
                        std::io::ErrorKind::Other,
                        format!("Failed to read Q data: {}", e),
                    ))
                })?;

            // Convert to complex SLC data with parallel processing
            Self::convert_to_complex_parallel(i_data, q_data, width, height)?
        } else {
            return Err(SarError::InvalidFormat(format!(
                "Unexpected number of bands in TIFF: {}",
                band_count
            )));
        };
        let gdal_time = gdal_start.elapsed();
        log::debug!("GDAL read took: {:?}", gdal_time);

        let total_time = start_time.elapsed();

        log::info!(
            "SLC data read complete (PARALLEL) for {}: {} x {} pixels",
            pol,
            width,
            height
        );
        log::info!("Performance timing (PARALLEL):");
        log::info!("  - Extraction: {:?}", extract_time);
        log::info!("  - GDAL read: {:?}", gdal_time);
        log::info!("  - Total: {:?}", total_time);

        let mb_size = (width * height * 8) as f64 / (1024.0 * 1024.0);
        let mb_per_sec = mb_size / total_time.as_secs_f64();
        log::info!("  - Data size: {:.1} MB", mb_size);
        log::info!("  - Throughput: {:.1} MB/s", mb_per_sec);

        Ok(slc_data)
    }

    /// Read SLC data with streaming (no temp files) - OPTIMIZED
    pub fn parse_annotation_xml_comprehensive(
        xml_content: &str,
        pol: Polarization,
        product_path: &Path,
    ) -> SarResult<SarMetadata> {
        log::info!(
            "Starting comprehensive XML parser for polarization {:?}",
            pol
        );

        let product_id = Self::product_id_from_path(product_path);
        log::debug!("Product ID: {}", product_id);

        let annotation = crate::io::annotation::parse_annotation_xml(xml_content).map_err(|e| {
            log::error!("Failed to parse annotation XML: {}", e);
            SarError::Processing(
                "Failed to parse annotation XML - cannot extract required parameters".to_string(),
            )
        })?;

        log::info!("Annotation parsing completed successfully");

        Self::metadata_from_annotation_root(&annotation, pol, &product_id)
    }

    fn metadata_from_annotation_root(
        annotation: &crate::io::annotation::ProductRoot,
        pol: Polarization,
        product_id: &str,
    ) -> SarResult<SarMetadata> {
        let (incidence_angle_near, incidence_angle_far) =
            Self::derive_incidence_angle_range(annotation)?;

        // Extract mid-swath incidence angle from annotation for RTC normalization
        let incidence_angle_mid_swath = annotation
            .image_annotation
            .as_ref()
            .and_then(|img| img.image_information.as_ref())
            .and_then(|info| info.incidence_angle_mid_swath);

        let range_pixel_spacing = annotation
            .image_annotation
            .as_ref()
            .and_then(|img| img.image_information.as_ref())
            .and_then(|info| info.range_pixel_spacing)
            .ok_or_else(|| SarError::Processing("Range pixel spacing not found in annotation XML - required for scientific accuracy".to_string()))?;

        let azimuth_pixel_spacing = annotation
            .image_annotation
            .as_ref()
            .and_then(|img| img.image_information.as_ref())
            .and_then(|info| info.azimuth_pixel_spacing)
            .ok_or_else(|| SarError::Processing("Azimuth pixel spacing not found in annotation XML - required for scientific accuracy".to_string()))?;

        let slant_range_time = annotation
            .image_annotation
            .as_ref()
            .and_then(|img| img.image_information.as_ref())
            .and_then(|info| info.slant_range_time)
            .ok_or_else(|| {
                SarError::Processing(
                    "Slant range time not found in annotation XML - required".to_string(),
                )
            })?;

        let range_sampling_rate = annotation
            .general_annotation
            .as_ref()
            .and_then(|gen| gen.product_information.as_ref())
            .map(|prod| prod.range_sampling_rate)
            .ok_or_else(|| {
                SarError::Processing(
                    "Range sampling rate not found in annotation XML - required".to_string(),
                )
            })?;

        let prf = annotation.get_pulse_repetition_frequency().map_err(|e| {
            log::error!("PRF extraction failed: {}", e);
            SarError::Processing(format!(
                "Failed to extract consistent PRF from annotation XML: {}",
                e
            ))
        })?;

        log::info!("Final PRF value extracted (strict): {} Hz", prf);

        let radar_frequency = annotation
            .general_annotation
            .as_ref()
            .and_then(|gen| gen.product_information.as_ref())
            .map(|prod| prod.radar_frequency)
            .ok_or_else(|| {
                SarError::Processing(
                    "Radar frequency not found in annotation XML - required for wavelength"
                        .to_string(),
                )
            })?;

        let wavelength = SPEED_OF_LIGHT_M_S / radar_frequency;

        log::info!("Successfully extracted parameters with comprehensive parser:");
        log::info!("  - Range pixel spacing: {:.4} m", range_pixel_spacing);
        log::info!("  - Azimuth pixel spacing: {:.4} m", azimuth_pixel_spacing);
        log::info!("  - Slant range time: {:.6} s", slant_range_time);
        log::info!("  - Range sampling rate: {:.0} Hz", range_sampling_rate);
        log::info!("  - PRF: {:.1} Hz", prf);
        log::info!("  - Radar frequency: {:.3} GHz", radar_frequency / 1e9);
        log::info!("  - Wavelength: {:.4} m", wavelength);
        log::info!(
            "  - Incidence angles: near {:.2}°, far {:.2}°",
            incidence_angle_near,
            incidence_angle_far
        );

        let mission = if product_id.starts_with("S1A") {
            "Sentinel-1A".to_string()
        } else if product_id.starts_with("S1B") {
            "Sentinel-1B".to_string()
        } else {
            "Sentinel-1".to_string()
        };

        let platform = mission.clone();
        let acquisition_mode = AcquisitionMode::IW;

        let (start_time, stop_time) = if let Some(ads_header) = &annotation.ads_header {
            let start_time = ads_header.start_time.as_deref().ok_or_else(|| {
                SarError::Processing("Start time missing in ADS header".to_string())
            })?;
            let stop_time = ads_header.stop_time.as_deref().ok_or_else(|| {
                SarError::Processing("Stop time missing in ADS header".to_string())
            })?;
            (
                Self::parse_time_flexible(start_time)?,
                Self::parse_time_flexible(stop_time)?,
            )
        } else {
            return Err(SarError::Processing(
                "ADS header with start/stop time is required".to_string(),
            ));
        };

        let mut sub_swaths = crate::io::annotation::ProductRoot::extract_subswaths(annotation).map_err(|e| {
            log::error!(
                "❌ SCIENTIFIC MODE FAILURE: Failed to extract subswaths from annotation: {}",
                e
            );
            SarError::XmlParsing(format!(
                "Failed to extract subswaths from annotation XML. No fallback parsers allowed for scientific accuracy. Error: {}",
                e
            ))
        })?;

        let bounding_box = match crate::io::annotation::ProductRoot::extract_bounding_box(
            annotation,
        ) {
            Ok(real_bbox) => {
                log::info!(
                    "Extracted real bounding box: [{:.6}, {:.6}, {:.6}, {:.6}]",
                    real_bbox.min_lon,
                    real_bbox.min_lat,
                    real_bbox.max_lon,
                    real_bbox.max_lat
                );
                real_bbox
            }
            Err(e) => {
                log::warn!("Failed to extract bounding box from annotation: {}", e);
                return Err(SarError::XmlParsing(
                    "Failed to extract valid bounding box from annotation XML with structured parser. No fallback parsers allowed as requested.".to_string(),
                ));
            }
        };

        let orbit_data = if let Some(orbit_list) = &annotation.orbit_list {
            if !orbit_list.is_empty() {
                log::info!(
                    "Extracted {} orbit state vectors for velocity calculations",
                    orbit_list.len()
                );

                // Use earliest orbit time as reference rather than product start
                let reference_time = orbit_list
                    .iter()
                    .min_by_key(|sv| sv.time)
                    .map(|sv| sv.time)
                    .unwrap_or(start_time);
                if reference_time == start_time {
                    log::warn!(
                        "Orbit reference time matches product start; verify timing offsets upstream"
                    );
                } else {
                    let delta = start_time.signed_duration_since(reference_time);
                    let delta_secs = delta.num_nanoseconds().unwrap_or(0) as f64 * 1e-9;
                    log::info!(
                        "Orbit reference epoch {:?} precedes product start by {:.3} s",
                        reference_time,
                        delta_secs
                    );
                }
                Some(OrbitData {
                    state_vectors: orbit_list.clone(),
                    reference_time,
                })
            } else {
                log::warn!("No orbit state vectors found in annotation");
                None
            }
        } else {
            log::warn!("No orbit data found in annotation");
            None
        };

        // FIX A1: Pass DateTime<Utc> directly instead of converting to f64
        let orbit_ref_epoch: Option<DateTime<Utc>> =
            orbit_data.as_ref().map(|data| data.reference_time);

        if orbit_ref_epoch.is_none() {
            log::warn!(
                "⚠️  Unable to derive orbit reference epoch; azimuthFmRate entries remain without orbit-relative timestamps"
            );
        }

        Self::annotate_fm_rate_relative_times(&mut sub_swaths, orbit_ref_epoch)?;

        log::info!("Creating SarMetadata with PRF = {} Hz", prf);

        let burst_records = Self::collect_burst_records(annotation, &sub_swaths, orbit_ref_epoch);

        let doppler_centroid = {
            let est = annotation
                .general_annotation
                .as_ref()
                .and_then(|ga| ga.dc_estimate_list.as_ref())
                .and_then(|dl| dl.dc_estimates.as_ref())
                .and_then(|v| v.first().cloned())
                .or_else(|| {
                    annotation
                        .doppler_centroid
                        .as_ref()
                        .and_then(|dc| dc.dc_estimate_list.as_ref())
                        .and_then(|dl| dl.dc_estimates.as_ref())
                        .and_then(|v| v.first().cloned())
                });

            if let Some(ref estimate) = est {
                crate::io::parsing_validation::validate_dc_polynomial(
                    &estimate.data_dc_polynomial,
                )?;
            }

            est.map(|estimate| DopplerCentroidPolynomial {
                t0: estimate.t0,
                coefficients: estimate.data_dc_polynomial,
            })
        };

        Ok(SarMetadata {
            product_id: product_id.to_string(),
            mission,
            platform,
            instrument: "C-SAR".to_string(),
            acquisition_mode,
            polarizations: vec![pol],
            start_time,
            stop_time,
            radar_frequency: Some(radar_frequency),
            wavelength: Some(wavelength),
            slant_range_time: Some(slant_range_time),
            prf: Some(prf),
            range_sampling_rate: Some(range_sampling_rate), // ✅ Fixed: Now extracted from annotation XML
            incidence_angle_near: Some(incidence_angle_near),
            incidence_angle_far: Some(incidence_angle_far),
            incidence_angle_mid_swath,
            radar_frequency_extracted: true,
            bounding_box,
            coordinate_system: CoordinateSystem::Radar,
            sub_swaths,
            orbit_data,
            range_looks: 1,
            azimuth_looks: 1,
            pixel_spacing: (range_pixel_spacing, azimuth_pixel_spacing),
            burst_records,
            doppler_centroid,
        })
    }

    fn derive_incidence_angle_range(
        annotation: &crate::io::annotation::ProductRoot,
    ) -> SarResult<(f64, f64)> {
        let mut min_angle = f64::INFINITY;
        let mut max_angle = f64::NEG_INFINITY;

        if let Some(geo) = &annotation.geolocation_grid {
            if let Some(list) = &geo.geolocation_grid_point_list {
                if let Some(points) = &list.geolocation_grid_points {
                    for point in points {
                        let angle = point.incidence_angle;
                        if angle.is_finite() {
                            min_angle = min_angle.min(angle);
                            max_angle = max_angle.max(angle);
                        }
                    }
                }
            }
        }

        if !(min_angle.is_finite() && max_angle.is_finite() && max_angle > min_angle) {
            min_angle = f64::INFINITY;
            max_angle = f64::NEG_INFINITY;

            match annotation.get_antenna_patterns() {
                Ok(patterns) => {
                    for (_, incidence, _) in patterns.iter() {
                        for &angle in incidence {
                            if angle.is_finite() {
                                min_angle = min_angle.min(angle);
                                max_angle = max_angle.max(angle);
                            }
                        }
                    }
                }
                Err(err) => {
                    log::warn!(
                        "⚠️  Failed to extract antenna-based incidence angles: {}",
                        err
                    );
                }
            }
        }

        if !(min_angle.is_finite() && max_angle.is_finite() && max_angle > min_angle) {
            if let Some(mid) = annotation
                .image_annotation
                .as_ref()
                .and_then(|img| img.image_information.as_ref())
                .and_then(|info| info.incidence_angle_mid_swath)
            {
                if mid.is_finite() {
                    min_angle = min_angle.min(mid);
                    max_angle = max_angle.max(mid);
                }
            }
        }

        if min_angle.is_finite() && max_angle.is_finite() && max_angle > min_angle {
            log::info!(
                "Derived incidence angle range from annotation: near {:.2}°, far {:.2}°",
                min_angle,
                max_angle
            );
            Ok((min_angle, max_angle))
        } else {
            Err(SarError::Metadata(
                "Unable to derive incidence angle range from annotation".to_string(),
            ))
        }
    }

    /// Collect burst records from annotation with correct TOPS geometry
    ///
    /// **Scientific Fix (2025-01):** Uses orbit-relative timing and constant range windows
    fn collect_burst_records(
        annotation: &crate::io::annotation::ProductRoot,
        sub_swaths: &HashMap<String, SubSwath>,
        orbit_ref_epoch: Option<DateTime<Utc>>,
    ) -> Vec<BurstRecord> {
        let mut records = Vec::new();

        if orbit_ref_epoch.is_none() {
            log::warn!(
                "⚠️  Orbit reference epoch unavailable; burst timings will use annotation azimuth_anx_time values"
            );
        } else {
            log::info!(
                "🕐 Using orbit reference epoch: {}",
                orbit_ref_epoch.unwrap().to_rfc3339()
            );
        }

        // In strict scientific mode, delegate burst record construction to the dedicated
        // `slc_reader::burst_records` implementation, which enforces constant range
        // windows, orbit-relative timing, and additional tiling invariants. This keeps
        // all burst geometry logic in one place and avoids silent fallbacks.
        let strict = std::env::var("SARDINE_STRICT")
            .map(|v| v == "1" || v.to_lowercase() == "true")
            .unwrap_or(false);

        for (swath_id, subswath) in sub_swaths {
            let mut per_swath_records = if strict {
                crate::io::sentinel::slc_reader::burst_records::build_burst_records_for_subswath_fixed(
                    annotation,
                    swath_id,
                    subswath,
                    orbit_ref_epoch,
                )
            } else {
                Self::build_burst_records_for_subswath(
                    annotation,
                    swath_id,
                    subswath,
                    orbit_ref_epoch,
                )
            };

            if per_swath_records.is_empty() {
                log::warn!(
                    "⚠️  No burst records built for subswath {} (strict_mode={})",
                    swath_id,
                    strict
                );
            }

            records.append(&mut per_swath_records);
        }

        records
    }

    /// Build burst records for a subswath with CORRECT TOPS geometry
    ///
    /// **Critical fixes applied (2025-01):**
    /// 1. **A2: Range window is CONSTANT** - bursts tile in azimuth only, not range
    /// 2. **A3: Valid sample bounds are EXCLUSIVE** - last_valid_sample = lastValidSample + 1
    /// 3. **A1: Timing uses orbit-reference-relative seconds** when epoch available
    fn build_burst_records_for_subswath(
        annotation: &crate::io::annotation::ProductRoot,
        subswath_id: &str,
        subswath: &SubSwath,
        orbit_ref_epoch: Option<DateTime<Utc>>,
    ) -> Vec<BurstRecord> {
        use slc_reader::time_utils::{parse_iso8601_to_datetime_utc, seconds_since_epoch};

        let swath_timing = match &annotation.swath_timing {
            Some(timing) => timing,
            None => return Vec::new(),
        };
        let bursts = match swath_timing
            .burst_list
            .as_ref()
            .and_then(|list| list.bursts.as_ref())
        {
            Some(list) if !list.is_empty() => list,
            _ => return Vec::new(),
        };

        let total_lines = subswath
            .last_line_global
            .saturating_sub(subswath.first_line_global)
            .saturating_add(1);
        if total_lines == 0 {
            return Vec::new();
        }

        let total_bursts = bursts.len().max(1);
        let lines_per_burst = swath_timing
            .lines_per_burst
            .map(|v| v as usize)
            .filter(|&v| v > 0)
            .unwrap_or_else(|| std::cmp::max(1, total_lines / total_bursts));

        // =========================================================================
        // FIX A2: RANGE WINDOW IS CONSTANT FOR ALL BURSTS IN A SUBSWATH
        // In TOPS/IW SLC, bursts tile in AZIMUTH, not range.
        // All bursts share the same range extent.
        // =========================================================================
        let constant_first_sample = subswath.first_sample_global;
        let constant_last_sample = subswath.last_sample_global;

        log::debug!(
            "📐 {} burst geometry: {} bursts, {} lines/burst, range=[{}, {}] (constant)",
            subswath_id,
            total_bursts,
            lines_per_burst,
            constant_first_sample,
            constant_last_sample
        );

        let mut current_line = subswath.first_line_global;
        let mut records = Vec::with_capacity(bursts.len());

        for (idx, burst) in bursts.iter().enumerate() {
            // Azimuth (line) window advances per burst
            let first_line = current_line;
            let mut last_line = first_line.saturating_add(lines_per_burst).saturating_sub(1);
            if idx == total_bursts - 1 || last_line > subswath.last_line_global {
                last_line = subswath.last_line_global;
            }
            current_line = last_line.saturating_add(1);

            // =====================================================================
            // FIX A3: VALID SAMPLE BOUNDS USE EXCLUSIVE UPPER BOUND
            // first_valid_sample: inclusive (first valid sample index)
            // last_valid_sample: exclusive (one past last valid sample)
            // =====================================================================
            let first_valid_sample = burst
                .first_valid_sample
                .iter()
                .filter(|&&v| v >= 0)
                .map(|&v| constant_first_sample.saturating_add(v as usize))
                .min()
                .unwrap_or(constant_first_sample);

            // Convert to EXCLUSIVE upper bound: lastValidSample + 1
            let last_valid_sample_exclusive = burst
                .last_valid_sample
                .iter()
                .filter(|&&v| v >= 0)
                .map(|&v| {
                    constant_first_sample
                        .saturating_add(v as usize)
                        .saturating_add(1)
                })
                .max()
                .unwrap_or(constant_last_sample.saturating_add(1));

            // Clamp to subswath range window
            let first_valid_clamped = first_valid_sample.max(constant_first_sample);
            let last_valid_exclusive_clamped =
                last_valid_sample_exclusive.min(constant_last_sample.saturating_add(1));

            // =====================================================================
            // FIX A1: TIMING USES ORBIT-REFERENCE-RELATIVE SECONDS
            // =====================================================================
            let burst_datetime = burst
                .azimuth_time
                .as_ref()
                .or(burst.sensing_time.as_ref())
                .and_then(|t| parse_iso8601_to_datetime_utc(t));

            let azimuth_time_absolute = burst_datetime
                .map(|dt| dt.timestamp() as f64 + (dt.timestamp_subsec_nanos() as f64) * 1e-9);

            let azimuth_time_rel_orbit = match (orbit_ref_epoch, burst_datetime) {
                (Some(epoch), Some(burst_dt)) => {
                    let rel = seconds_since_epoch(epoch, burst_dt);
                    if rel.is_finite() {
                        Some(rel)
                    } else {
                        log::warn!(
                            "⚠️  Non-finite orbit-relative time; falling back to azimuth_anx_time"
                        );
                        if burst.azimuth_anx_time.is_finite() {
                            Some(burst.azimuth_anx_time)
                        } else {
                            None
                        }
                    }
                }
                _ => {
                    // Fallback to azimuth_anx_time
                    if burst.azimuth_anx_time.is_finite() && burst.azimuth_anx_time > 0.0 {
                        if burst_datetime.is_none() && idx == 0 {
                            log::warn!(
                                "⚠️  Burst datetime unavailable; using annotation azimuth_anx_time"
                            );
                        }
                        Some(burst.azimuth_anx_time)
                    } else {
                        None
                    }
                }
            };

            records.push(BurstRecord {
                subswath_id: subswath_id.to_string(),
                burst_index: idx,
                // Azimuth bounds (advance per burst)
                first_line_global: first_line,
                last_line_global: last_line,
                // FIX A2: Range bounds are CONSTANT for all bursts
                start_sample_global: constant_first_sample,
                end_sample_global: constant_last_sample,
                // FIX A3: Valid bounds with EXCLUSIVE upper bound
                first_valid_sample: Some(first_valid_clamped),
                last_valid_sample: Some(last_valid_exclusive_clamped),
                first_valid_line: Some(first_line),
                last_valid_line: Some(last_line),
                // FIX A1: Orbit-relative timing
                azimuth_time_rel_orbit,
                azimuth_time_absolute,
            });
        }

        log::info!(
            "✅ {} built {} burst records: lines={}..{}, samples={}..{} (constant)",
            subswath_id,
            records.len(),
            subswath.first_line_global,
            subswath.last_line_global,
            constant_first_sample,
            constant_last_sample
        );

        records
    }

    fn parse_iso8601_to_unix_seconds(value: &str) -> Option<f64> {
        let trimmed = value.trim();
        if trimmed.is_empty() {
            return None;
        }

        if let Ok(dt) = DateTime::parse_from_rfc3339(trimmed) {
            let utc = dt.with_timezone(&Utc);
            return Some(utc.timestamp() as f64 + (utc.timestamp_subsec_nanos() as f64) * 1e-9);
        }

        // ASF SAFE products often omit an explicit timezone suffix. Treat those as UTC.
        if let Ok(naive) = NaiveDateTime::parse_from_str(trimmed, "%Y-%m-%dT%H:%M:%S%.f") {
            let utc = DateTime::<Utc>::from_naive_utc_and_offset(naive, Utc);
            return Some(utc.timestamp() as f64 + (utc.timestamp_subsec_nanos() as f64) * 1e-9);
        }

        if let Ok(naive) = NaiveDateTime::parse_from_str(trimmed, "%Y-%m-%dT%H:%M:%S") {
            let utc = DateTime::<Utc>::from_naive_utc_and_offset(naive, Utc);
            return Some(utc.timestamp() as f64 + (utc.timestamp_subsec_nanos() as f64) * 1e-9);
        }

        None
    }

    fn validate_burst_timing_continuity(
        pol: Polarization,
        annotation: &crate::io::annotation::ProductRoot,
    ) -> SarResult<()> {
        let swath_id = Self::guess_swath_id(annotation);
        let swath_timing = match &annotation.swath_timing {
            Some(timing) => timing,
            None => {
                log::warn!(
                    "⚠️  No swathTiming block for {} {}; skipping burst continuity check",
                    swath_id,
                    pol
                );
                if crate::types::strict_mode() {
                    return Err(SarError::Metadata(format!(
                        "Missing swathTiming block for {} {}",
                        swath_id, pol
                    )));
                }
                return Ok(());
            }
        };
        let bursts = match swath_timing
            .burst_list
            .as_ref()
            .and_then(|bl| bl.bursts.as_ref())
        {
            Some(b) if !b.is_empty() => b,
            _ => {
                log::warn!(
                    "⚠️  No bursts found in swathTiming for {} {}; cannot validate timing continuity",
                    swath_id,
                    pol
                );
                if crate::types::strict_mode() {
                    return Err(SarError::Metadata(format!(
                        "Missing burst records for {} {}",
                        swath_id, pol
                    )));
                }
                return Ok(());
            }
        };

        if let Some(expected) = swath_timing.burst_list.as_ref().and_then(|bl| bl.count) {
            if expected as usize != bursts.len() {
                log::warn!(
                    "⚠️  Burst count mismatch for {} {}: declared {} vs parsed {}",
                    swath_id,
                    pol,
                    expected,
                    bursts.len()
                );
            }
        }

        let mut prev_time = None;
        let mut max_deviation = 0.0f64;
        let mut discontinuities = 0usize;
        let mut deltas: Vec<f64> = Vec::new();

        for burst in bursts {
            let current = burst
                .azimuth_time
                .as_ref()
                .and_then(|t| Self::parse_iso8601_to_unix_seconds(t))
                .or_else(|| {
                    if burst.azimuth_anx_time.is_finite() {
                        Some(burst.azimuth_anx_time)
                    } else {
                        None
                    }
                });

            if let (Some(prev), Some(cur)) = (prev_time, current) {
                let dt: f64 = cur - prev;
                if dt <= 0.0 {
                    log::warn!(
                        "⚠️  Non-increasing burst time detected in {} {}: dt={:.6}s",
                        swath_id,
                        pol,
                        dt
                    );
                    discontinuities += 1;
                } else {
                    deltas.push(dt);
                }
            }
            prev_time = current;
        }

        let expected_dt = if deltas.is_empty() {
            None
        } else {
            let mut sorted = deltas.clone();
            sorted.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
            let median = if sorted.len() % 2 == 1 {
                sorted[sorted.len() / 2]
            } else {
                let mid = sorted.len() / 2;
                0.5 * (sorted[mid - 1] + sorted[mid])
            };

            for dt in &deltas {
                let dev = (dt - median).abs();
                max_deviation = max_deviation.max(dev);
                if dev > median * 0.02 {
                    discontinuities += 1;
                    log::warn!(
                        "⚠️  Burst spacing drift in {} {}: dt={:.6}s median≈{:.6}s (dev {:.6}s)",
                        swath_id,
                        pol,
                        dt,
                        median,
                        dev
                    );
                }
            }

            Some(median)
        };

        if discontinuities == 0 {
            if let Some(exp) = expected_dt {
                log::info!(
                    "✅ Burst timing continuity OK for {} {} (median spacing {:.6}s, max dev {:.6}s)",
                    swath_id,
                    pol,
                    exp,
                    max_deviation
                );
            } else {
                log::info!(
                    "✅ Burst timing continuity OK for {} {} (no expected spacing available)",
                    swath_id,
                    pol
                );
            }
            return Ok(());
        }

        log::error!(
            "❌ Burst timing discontinuities detected for {} {} ({} anomalies)",
            swath_id,
            pol,
            discontinuities
        );
        if crate::types::strict_mode() {
            return Err(SarError::Metadata(format!(
                "Burst timing continuity failed for {} {} ({} anomalies)",
                swath_id, pol, discontinuities
            )));
        }

        Ok(())
    }

    fn guess_swath_id(annotation: &crate::io::annotation::ProductRoot) -> String {
        if let Ok(subs) = crate::io::annotation::ProductRoot::extract_subswaths(annotation) {
            if let Some(first) = subs.keys().next() {
                return first.clone();
            }
        }
        annotation
            .image_annotation
            .as_ref()
            .and_then(|ia| ia.processing_information.as_ref())
            .and_then(|pi| {
                pi.swath_proc_params_list
                    .as_ref()
                    .and_then(|l| l.swath_proc_params.as_ref())
                    .and_then(|v| v.first())
                    .and_then(|p| p.swath.clone())
            })
            .unwrap_or_else(|| "UNKNOWN".to_string())
    }

    fn summarize_lut_presence(
        &self,
        pol: Polarization,
        calibration: Option<&CalibrationCoefficients>,
        noise_present: bool,
    ) -> SarResult<()> {
        let (sigma_vecs, beta_vecs, gamma_vecs) =
            calibration.map_or((0usize, 0usize, 0usize), |c| {
                c.vectors
                    .iter()
                    .fold((0usize, 0usize, 0usize), |mut acc, v| {
                        if !v.sigma_nought.is_empty() {
                            acc.0 += 1;
                        }
                        if !v.beta_nought.is_empty() {
                            acc.1 += 1;
                        }
                        if !v.gamma.is_empty() {
                            acc.2 += 1;
                        }
                        acc
                    })
            });

        let mut antenna_present = calibration
            .map(|c| c.antenna_pattern_lut.is_some() || !c.antenna_pattern_vectors.is_empty())
            .unwrap_or(false);

        if !antenna_present {
            match self.cached_annotations.get(&pol) {
                Some(annotation_roots) => {
                    for root in annotation_roots {
                        match root.get_antenna_patterns() {
                            Ok(patterns) => {
                                if !patterns.is_empty() {
                                    antenna_present = true;
                                    log::info!(
                                        "📡 Antenna pattern entries available in cached annotations for {:?} ({} entries)",
                                        pol,
                                        patterns.len()
                                    );
                                    break;
                                }
                            }
                            Err(e) => {
                                let msg = format!(
                                    "Failed to extract antenna patterns from annotation for {:?}: {}",
                                    pol,
                                    e
                                );
                                if crate::types::strict_mode() {
                                    return Err(SarError::MissingCalibrationData(msg));
                                } else {
                                    log::warn!("⚠️  {}", msg);
                                }
                            }
                        }
                    }

                    if !antenna_present {
                        log::warn!(
                            "⚠️  Cached annotations for {:?} do not include antennaPattern entries",
                            pol
                        );
                    }
                }
                None => {
                    log::warn!(
                        "⚠️  No cached annotations available for {:?}; antenna pattern availability unknown",
                        pol
                    );
                }
            }
        }

        log::info!(
            "Calibration/LUT presence for {}: sigma={} beta={} gamma={} antenna={} noise_file={}",
            pol,
            sigma_vecs,
            beta_vecs,
            gamma_vecs,
            antenna_present,
            noise_present
        );

        if sigma_vecs == 0 {
            log::warn!("⚠️  Missing sigma⁰ calibration vectors for {}", pol);
        }
        if beta_vecs == 0 {
            log::warn!("⚠️  Missing beta⁰ calibration vectors for {}", pol);
        }
        if gamma_vecs == 0 {
            log::warn!("⚠️  Missing gamma⁰ calibration vectors for {}", pol);
        }
        if !antenna_present {
            log::warn!("⚠️  Missing antenna pattern vectors/LUT for {}", pol);
            if crate::types::strict_mode() {
                return Err(SarError::MissingCalibrationData(format!(
                    "Antenna pattern LUT missing for {}",
                    pol
                )));
            }
        }
        if !noise_present {
            log::warn!("⚠️  Missing thermal noise XML for {}", pol);
            if crate::types::strict_mode() {
                return Err(SarError::MissingCalibrationData(format!(
                    "Thermal noise removal metadata missing for {}",
                    pol
                )));
            }
        }

        Ok(())
    }

    fn merge_metadata(mut metas: Vec<SarMetadata>) -> SarResult<SarMetadata> {
        if metas.is_empty() {
            return Err(SarError::InvalidFormat(
                "No annotation metadata available to merge".to_string(),
            ));
        }

        let mut base = metas.remove(0);
        let mut global_prf_conflict = false;
        let strict = std::env::var("SARDINE_REQUIRE_SUBSWATHS").ok().as_deref() == Some("1");
        for mut m in metas {
            base.bounding_box.min_lat = base.bounding_box.min_lat.min(m.bounding_box.min_lat);
            base.bounding_box.max_lat = base.bounding_box.max_lat.max(m.bounding_box.max_lat);
            base.bounding_box.min_lon = base.bounding_box.min_lon.min(m.bounding_box.min_lon);
            base.bounding_box.max_lon = base.bounding_box.max_lon.max(m.bounding_box.max_lon);
            // Normalize after each merge step to ensure invariants
            base.bounding_box.normalize();

            let other_prf = m.prf;
            let other_radar_frequency = m.radar_frequency;
            let other_incidence_near = m.incidence_angle_near;
            let other_incidence_far = m.incidence_angle_far;

            log::debug!(
                "Merging metadata: base_prf={:?}, other_prf={:?}, base_freq={:?}, other_freq={:?}",
                base.prf,
                other_prf,
                base.radar_frequency,
                other_radar_frequency
            );

            // Fill in missing core parameters from other subswaths if available
            if base.prf.is_none() && !global_prf_conflict {
                base.prf = other_prf;
            }
            if base.radar_frequency.is_none() {
                base.radar_frequency = other_radar_frequency;
            }
            if base.incidence_angle_near.is_none() {
                base.incidence_angle_near = other_incidence_near;
            }
            if base.incidence_angle_far.is_none() {
                base.incidence_angle_far = other_incidence_far;
            }

            for (k, v) in m.sub_swaths.drain() {
                base.sub_swaths.insert(k, v);
            }

            base.burst_records.extend(m.burst_records);

            for pol in m.polarizations {
                if !base.polarizations.contains(&pol) {
                    base.polarizations.push(pol);
                }
            }

            if strict {
                let approx_equal = |a: f64, b: f64, rel_tol: f64, abs_tol: f64| {
                    let diff = (a - b).abs();
                    if diff <= abs_tol {
                        true
                    } else {
                        let max_mag = a.abs().max(b.abs());
                        diff <= rel_tol * max_mag
                    }
                };

                let prf_mismatch = match (base.prf, other_prf) {
                    (Some(a), Some(b)) => !approx_equal(a, b, 1e-6, 1e-3),
                    (None, None) => false,
                    _ => true,
                };

                let freq_mismatch = match (base.radar_frequency, other_radar_frequency) {
                    (Some(a), Some(b)) => !approx_equal(a, b, 1e-9, 1e-3),
                    (None, None) => false,
                    _ => true,
                };

                if prf_mismatch || freq_mismatch {
                    log::warn!(
                        "Subswath metadata mismatch: base_prf={:?}, other_prf={:?}, base_freq={:?}, other_freq={:?}",
                        base.prf, other_prf, base.radar_frequency, other_radar_frequency
                    );
                }

                if prf_mismatch {
                    base.prf = None;
                    global_prf_conflict = true;
                }
            }

            base.incidence_angle_near = match (base.incidence_angle_near, other_incidence_near) {
                (Some(a), Some(b)) => Some(a.min(b)),
                (Some(a), None) => Some(a),
                (None, Some(b)) => Some(b),
                (None, None) => None,
            };

            base.incidence_angle_far = match (base.incidence_angle_far, other_incidence_far) {
                (Some(a), Some(b)) => Some(a.max(b)),
                (Some(a), None) => Some(a),
                (None, Some(b)) => Some(b),
                (None, None) => None,
            };
        }

        // Do NOT derive a single representative PRF across subswaths; keep per-subswath PRF.
        if base.sub_swaths.len() > 1 {
            let mut prfs: Vec<f64> = base
                .sub_swaths
                .values()
                .filter_map(|sw| sw.prf_hz)
                .collect();

            prfs.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
            prfs.dedup_by(|a, b| {
                let diff = (*a - *b).abs();
                let max_mag = a.abs().max(b.abs()).max(1.0);
                diff <= 1e-6 * max_mag || diff <= 1e-3
            });

            base.prf = if prfs.len() == 1 {
                prfs.first().copied()
            } else {
                log::warn!(
                    "⚠️  PRF differs across subswaths; skipping global PRF selection (use per-subswath PRF instead)"
                );
                None
            };
        }

        base.polarizations.sort_by_key(|p| match p {
            Polarization::VV => 0,
            Polarization::VH => 1,
            Polarization::HV => 2,
            Polarization::HH => 3,
        });
        base.polarizations.dedup();

        // Final safety normalization (idempotent) before returning
        base.bounding_box.normalize();

        Ok(base)
    }

    fn annotate_fm_rate_relative_times(
        sub_swaths: &mut HashMap<String, SubSwath>,
        orbit_ref_epoch: Option<DateTime<Utc>>,
    ) -> SarResult<()> {
        let strict = crate::types::strict_mode();

        for (swath_id, swath) in sub_swaths.iter_mut() {
            if let Some(models) = swath.fm_rate_estimates.as_mut() {
                let mut missing_time = false;

                for model in models.iter_mut() {
                    match (model.azimuth_time_utc.as_ref(), orbit_ref_epoch.as_ref()) {
                        (Some(az_utc), Some(ref_epoch)) => {
                            let delta = *az_utc - *ref_epoch;
                            let rel = delta.num_seconds() as f64
                                + (delta.subsec_nanos() as f64) * 1e-9;
                            model.azimuth_time_rel_orbit = Some(rel);
                        }
                        (Some(_), None) => {
                            // No orbit reference yet; leave relative time unset.
                            model.azimuth_time_rel_orbit = None;
                        }
                        (None, Some(_)) => {
                            missing_time = true;
                        }
                        (None, None) => {
                            // No orbit reference and no azimuth time; nothing to report yet.
                        }
                    }
                }

                if missing_time && orbit_ref_epoch.is_some() {
                    if strict {
                        return Err(SarError::Processing(format!(
                            "STRICT MODE: azimuthFmRate entries for {} lack azimuthTime even though precise orbit data is available. Cannot form orbit-relative FM polynomials.",
                            swath_id
                        )));
                    } else {
                        log::warn!(
                            "⚠️  azimuthFmRate entries for {} are missing azimuthTime with a known orbit reference; orbit-relative timestamps unavailable",
                            swath_id
                        );
                    }
                }
            }
        }

        Ok(())
    }

    fn product_id_from_path(product_path: &Path) -> String {
        product_path
            .file_stem()
            .and_then(|s| s.to_str())
            .unwrap_or("S1A_UNKNOWN")
            .replace(".SAFE", "")
            .replace(".zip", "")
    }

    /// Convert I/Q data to complex SLC data in parallel
    fn convert_to_complex_parallel(
        i_data: gdal::raster::Buffer<f32>,
        q_data: gdal::raster::Buffer<f32>,
        width: usize,
        height: usize,
    ) -> SarResult<SarImage> {
        use num_complex::Complex;

        let mut slc_array = Array2::zeros((height, width));

        // Convert buffer data to complex numbers
        for i in 0..height {
            for j in 0..width {
                let idx = i * width + j;
                if idx < i_data.data.len() && idx < q_data.data.len() {
                    slc_array[[i, j]] = Complex::new(i_data.data[idx], q_data.data[idx]);
                }
            }
        }

        Ok(slc_array)
    }

    /// Convert CInt16 interleaved data to complex SLC data with parallel processing
    fn convert_cint16_to_complex_parallel(
        complex_data: gdal::raster::Buffer<i16>,
        width: usize,
        height: usize,
    ) -> SarResult<SarImage> {
        use num_complex::Complex;
        #[cfg(feature = "parallel")]
        use rayon::prelude::*;

        let total_pixels = width * height;

        // Data is interleaved as [real, imag, real, imag, ...]
        if complex_data.data.len() < total_pixels * 2 {
            return Err(SarError::InvalidFormat(format!(
                "Insufficient data: expected {} i16 values, got {}",
                total_pixels * 2,
                complex_data.data.len()
            )));
        }

        let mut slc_array = Array2::zeros((height, width));

        #[cfg(feature = "parallel")]
        {
            // Process in parallel chunks by dividing into row chunks
            let chunk_size = std::cmp::max(1, height / rayon::current_num_threads());

            slc_array
                .axis_chunks_iter_mut(ndarray::Axis(0), chunk_size)
                .into_par_iter()
                .enumerate()
                .for_each(|(chunk_idx, mut chunk)| {
                    let start_row = chunk_idx * chunk_size;
                    for (local_row, mut row) in chunk.axis_iter_mut(ndarray::Axis(0)).enumerate() {
                        let global_row = start_row + local_row;
                        if global_row < height {
                            for col in 0..width {
                                let pixel_idx = global_row * width + col;
                                let data_idx = pixel_idx * 2;

                                if data_idx + 1 < complex_data.data.len() {
                                    let real_i16 = complex_data.data[data_idx];
                                    let imag_i16 = complex_data.data[data_idx + 1];

                                    let real_f32 = real_i16 as f32;
                                    let imag_f32 = imag_i16 as f32;

                                    row[col] = Complex::new(real_f32, imag_f32);
                                }
                            }
                        }
                    }
                });
        }

        #[cfg(not(feature = "parallel"))]
        {
            // Sequential processing fallback
            for row in 0..height {
                for col in 0..width {
                    let pixel_idx = row * width + col;
                    let data_idx = pixel_idx * 2;

                    if data_idx + 1 < complex_data.data.len() {
                        let real_i16 = complex_data.data[data_idx];
                        let imag_i16 = complex_data.data[data_idx + 1];

                        let real_f32 = real_i16 as f32;
                        let imag_f32 = imag_i16 as f32;

                        slc_array[[row, col]] = Complex::new(real_f32, imag_f32);
                    }
                }
            }
        }

        // Log some statistics
        let sample_pixel = slc_array[[0, 0]];
        log::debug!(
            "First pixel (parallel): real={}, imag={}, magnitude={:.2}",
            sample_pixel.re,
            sample_pixel.im,
            sample_pixel.norm()
        );

        Ok(slc_array)
    }

    /// Read SLC data for a subswath with streaming callback (for I/O optimization)
    ///
    /// This function reads SLC data in chunks and calls a callback for each chunk,
    /// enabling processing to happen while reading. This provides finer-grained
    /// I/O overlap than pre-reading entire subswaths.
    ///
    /// # Arguments
    /// * `subswath` - Subswath identifier (e.g., "IW1")
    /// * `pol` - Polarization
    /// * `chunk_callback` - Callback function called for each chunk of data
    /// * `chunk_lines` - Number of lines per chunk (default: 512)
    ///
    /// # Returns
    /// Ok(()) if all chunks processed successfully
    pub fn read_slc_data_streaming<F>(
        &mut self,
        subswath: &str,
        pol: Polarization,
        mut chunk_callback: F,
        chunk_lines: Option<usize>,
    ) -> SarResult<()>
    where
        F: FnMut(Array2<SarComplex>, usize, usize) -> SarResult<()>, // (chunk_data, start_line, end_line)
    {
        use tempfile::NamedTempFile;

        let measurement_file = self.find_measurement_file_for_subswath(subswath, pol)?;

        log::info!(
            "📖 Streaming SLC data for subswath {} {} from {}",
            subswath,
            pol,
            measurement_file
        );
        let start_time = std::time::Instant::now();

        // Handle different formats (same logic as read_slc_data_for_subswath)
        let dataset = match self.format {
            ProductFormat::Zip => {
                let archive = self.open_archive()?;
                let mut zip_file = archive.by_name(&measurement_file).map_err(|e| {
                    SarError::Io(std::io::Error::new(
                        std::io::ErrorKind::Other,
                        format!("Failed to access {}: {}", measurement_file, e),
                    ))
                })?;

                let mut temp_file = NamedTempFile::new().map_err(|e| SarError::Io(e))?;
                std::io::copy(&mut zip_file, &mut temp_file).map_err(|e| SarError::Io(e))?;

                let extract_time = start_time.elapsed();
                log::debug!("TIFF extraction took: {:?}", extract_time);

                gdal::Dataset::open(temp_file.path()).map_err(|e| {
                    SarError::Io(std::io::Error::new(
                        std::io::ErrorKind::Other,
                        format!("Failed to open TIFF with GDAL: {}", e),
                    ))
                })?
            }
            ProductFormat::Safe => {
                let full_path = self.product_path.join(&measurement_file);
                if !full_path.exists() {
                    return Err(SarError::Io(std::io::Error::new(
                        std::io::ErrorKind::NotFound,
                        format!("Measurement file not found: {:?}", full_path),
                    )));
                }

                gdal::Dataset::open(&full_path).map_err(|e| {
                    SarError::Io(std::io::Error::new(
                        std::io::ErrorKind::Other,
                        format!("Failed to open TIFF with GDAL: {}", e),
                    ))
                })?
            }
        };

        // Get raster dimensions
        let (width, height) = dataset.raster_size();
        log::info!(
            "📐 SLC raster size for {} {}: {} x {} (width x height)",
            subswath,
            pol,
            width,
            height
        );

        // Read complex data using chunked approach with callbacks
        let band = dataset.rasterband(1).map_err(|e| {
            SarError::Io(std::io::Error::new(
                std::io::ErrorKind::Other,
                format!("Failed to get rasterband: {}", e),
            ))
        })?;

        let chunk_size = chunk_lines.unwrap_or(512);
        let mut row = 0usize;

        while row < height {
            let h = std::cmp::min(chunk_size, height - row);
            let window = (0, row as isize);
            let window_size = (width, h);
            let buffer_size = (width * 2, h);

            // Read chunk
            let buf = band
                .read_as::<i16>(window, window_size, buffer_size, None)
                .map_err(|e| {
                    SarError::Io(std::io::Error::new(
                        std::io::ErrorKind::Other,
                        format!(
                            "Failed to read CInt16 chunk at row {} ({} lines): {}",
                            row, h, e
                        ),
                    ))
                })?;

            if buf.data.len() != width * h * 2 {
                return Err(SarError::InvalidFormat(format!(
                    "Chunk size mismatch at row {}: expected {} samples, got {}",
                    row,
                    width * h * 2,
                    buf.data.len()
                )));
            }

            // Convert chunk to complex array
            let mut chunk = Array2::<SarComplex>::zeros((h, width));
            let mut idx = 0usize;
            for r in 0..h {
                for c in 0..width {
                    let re = buf.data[idx] as f32;
                    let im = buf.data[idx + 1] as f32;
                    chunk[[r, c]] = SarComplex::new(re, im);
                    idx += 2;
                }
            }

            // Call callback with chunk
            log::info!("📤 Calling chunk callback for lines [{}..{})", row, row + h);
            chunk_callback(chunk, row, row + h)?;
            log::info!(
                "📥 Chunk callback completed for lines [{}..{})",
                row,
                row + h
            );

            row += h;
        }

        let elapsed = start_time.elapsed();
        log::info!(
            "✅ Streamed SLC for {} {} in {:?}: {} lines x {} samples",
            subswath,
            pol,
            elapsed,
            height,
            width
        );

        Ok(())
    }

    /// Streamed read of complex band in chunks to reduce peak memory.
    /// OPTIMIZATION: Parallel conversion from i16 to Complex<f32> within each chunk
    fn read_complex_band_chunked(
        band: &gdal::raster::RasterBand,
        width: usize,
        height: usize,
    ) -> SarResult<SarImage> {
        use num_complex::Complex;
        use rayon::prelude::*;

        let chunk_lines: usize = std::env::var("SARDINE_TIFF_CHUNK_LINES")
            .ok()
            .and_then(|v| v.parse().ok())
            .filter(|&v| v > 0)
            .unwrap_or(512);

        let mut out = Array2::<Complex<f32>>::zeros((height, width));
        let mut row = 0usize;
        while row < height {
            let h = std::cmp::min(chunk_lines, height - row);
            let window = (0, row as isize);
            let window_size = (width, h);
            let buffer_size = (width * 2, h);
            let buf = band
                .read_as::<i16>(window, window_size, buffer_size, None)
                .map_err(|e| {
                    SarError::Io(std::io::Error::new(
                        std::io::ErrorKind::Other,
                        format!(
                            "Failed to read CInt16 chunk at row {} ({} lines): {}",
                            row, h, e
                        ),
                    ))
                })?;
            if buf.data.len() != width * h * 2 {
                return Err(SarError::InvalidFormat(format!(
                    "Chunk size mismatch at row {}: expected {} samples, got {}",
                    row,
                    width * h * 2,
                    buf.data.len()
                )));
            }

            // OPTIMIZATION: Parallel conversion from i16 pairs to Complex<f32>
            // This is safe because we're writing to non-overlapping rows
            let mut chunk_slice = out.slice_mut(ndarray::s![row..row + h, ..]);
            let raw_data = &buf.data;

            // Process rows in parallel - each row is independent
            chunk_slice
                .axis_iter_mut(ndarray::Axis(0))
                .into_par_iter()
                .enumerate()
                .for_each(|(r, mut row_slice)| {
                    let row_start = r * width * 2;
                    for c in 0..width {
                        let idx = row_start + c * 2;
                        let re = raw_data[idx] as f32;
                        let im = raw_data[idx + 1] as f32;
                        row_slice[c] = Complex::new(re, im);
                    }
                });

            row += h;
        }
        Ok(out)
    }

    /// Get comprehensive orbit data for the product
    /// Checks: 1) SLC embedded data, 2) Local cache, 3) Download from ESA
    pub fn rget_orbit_data(&mut self, orbit_cache_dir: Option<&Path>) -> SarResult<OrbitData> {
        log::info!("Getting orbit data for SLC product");

        // Get product metadata first
        let annotations = self.find_all_annotation_files()?;
        let first_pol = annotations
            .keys()
            .next()
            .cloned()
            .ok_or_else(|| SarError::Processing("No polarizations found".to_string()))?;
        let metadata = self.read_annotation(first_pol)?;

        // Extract product ID from filename
        let product_id = self
            .product_path
            .file_stem()
            .and_then(|s| s.to_str())
            .ok_or_else(|| SarError::Processing("Invalid filename".to_string()))?
            .replace(".SAFE", "")
            .replace(".zip", "");

        // Check if SLC contains orbit data (it won't for Sentinel-1, but check anyway)
        let slc_orbit_data = metadata.orbit_data.as_ref();
        if slc_orbit_data.is_some() {
            log::info!("Orbit data found embedded in SLC");
        } else {
            log::info!("No orbit data in SLC - will use external sources");
        }

        // Use orbit manager to get orbit data
        // SCIENTIFIC REQUIREMENT: Orbit cache directory must be explicitly specified.
        // Treat SARDINE_ORBIT_CACHE env var as explicit configuration if argument is None.
        let cache_dir = match orbit_cache_dir {
            Some(dir) => dir.to_path_buf(),
            None => {
                if let Ok(env_path) = std::env::var("SARDINE_ORBIT_CACHE") {
                    let p = std::path::PathBuf::from(env_path);
                    if p.as_path().exists() {
                        p
                    } else {
                        return Err(SarError::Processing(format!(
                            "❌ SCIENTIFIC ERROR: SARDINE_ORBIT_CACHE points to a non-existent path: {}. Set SARDINE_ORBIT_CACHE to a valid directory or pass orbit_cache_dir.",
                            p.display()
                        )));
                    }
                } else {
                    return Err(SarError::Processing(
                        "❌ SCIENTIFIC ERROR: Orbit cache directory must be explicitly specified. Set SARDINE_ORBIT_CACHE or pass orbit_cache_dir."
                            .to_string(),
                    ));
                }
            }
        };
        let orbit_manager = OrbitManager::new(cache_dir);
        orbit_manager.get_orbit_data(&product_id, metadata.start_time)
    }

    /// Check orbit data availability status
    pub fn check_orbit_status(&mut self, orbit_cache_dir: Option<&Path>) -> SarResult<OrbitStatus> {
        let annotations = self.find_all_annotation_files()?;
        let first_pol = annotations
            .keys()
            .next()
            .cloned()
            .ok_or_else(|| SarError::Processing("No polarizations found".to_string()))?;
        let metadata = self.read_annotation(first_pol)?;

        // Extract product ID
        let product_id = self
            .product_path
            .file_stem()
            .and_then(|s| s.to_str())
            .ok_or_else(|| SarError::Processing("Invalid filename".to_string()))?
            .replace(".SAFE", "")
            .replace(".zip", "");

        // Check SLC embedded data
        let has_embedded = metadata.orbit_data.is_some();

        // SCIENTIFIC REQUIREMENT: Orbit cache directory must be explicitly specified.
        // Treat SARDINE_ORBIT_CACHE env var as explicit configuration if argument is None.
        let cache_dir = match orbit_cache_dir {
            Some(dir) => dir.to_path_buf(),
            None => {
                if let Ok(env_path) = std::env::var("SARDINE_ORBIT_CACHE") {
                    let p = std::path::PathBuf::from(env_path);
                    if p.as_path().exists() {
                        p
                    } else {
                        return Err(SarError::Processing(format!(
                            "❌ SCIENTIFIC ERROR: SARDINE_ORBIT_CACHE points to a non-existent path: {}. Set SARDINE_ORBIT_CACHE to a valid directory or pass orbit_cache_dir.",
                            p.display()
                        )));
                    }
                } else {
                    return Err(SarError::Processing(
                        "❌ SCIENTIFIC ERROR: Orbit cache directory must be explicitly specified. Set SARDINE_ORBIT_CACHE or pass orbit_cache_dir."
                            .to_string(),
                    ));
                }
            }
        };
        let orbit_manager = OrbitManager::new(cache_dir);
        let primary_orbit_type =
            crate::io::orbit::OrbitReader::determine_orbit_type(metadata.start_time);
        let fallback_orbit_type = match primary_orbit_type {
            crate::io::orbit::OrbitType::POEORB => crate::io::orbit::OrbitType::RESORB,
            crate::io::orbit::OrbitType::RESORB => crate::io::orbit::OrbitType::POEORB,
        };

        let has_primary_cached = orbit_manager.has_orbit_cached(&product_id, primary_orbit_type);
        let has_fallback_cached = orbit_manager.has_orbit_cached(&product_id, fallback_orbit_type);

        Ok(OrbitStatus {
            product_id,
            start_time: metadata.start_time,
            has_embedded,
            primary_orbit_type,
            has_primary_cached,
            fallback_orbit_type,
            has_fallback_cached,
            cache_dir: orbit_manager.get_cache_dir().to_path_buf(),
        })
    }

    /// Download orbit files to cache without loading
    pub fn download_orbit_files(
        &mut self,
        orbit_cache_dir: Option<&Path>,
    ) -> SarResult<Vec<PathBuf>> {
        let annotations = self.find_all_annotation_files()?;
        let first_pol = annotations
            .keys()
            .next()
            .cloned()
            .ok_or_else(|| SarError::Processing("No polarizations found".to_string()))?;
        let metadata = self.read_annotation(first_pol)?;

        let product_id = self
            .product_path
            .file_stem()
            .and_then(|s| s.to_str())
            .ok_or_else(|| SarError::Processing("Invalid filename".to_string()))?
            .replace(".SAFE", "")
            .replace(".zip", "");

        // SCIENTIFIC REQUIREMENT: Orbit cache directory must be explicitly specified.
        // Treat SARDINE_ORBIT_CACHE env var as explicit configuration if argument is None.
        let cache_dir = match orbit_cache_dir {
            Some(dir) => dir.to_path_buf(),
            None => {
                if let Ok(env_path) = std::env::var("SARDINE_ORBIT_CACHE") {
                    let p = std::path::PathBuf::from(env_path);
                    if p.as_path().exists() {
                        p
                    } else {
                        return Err(SarError::Processing(format!(
                            "❌ SCIENTIFIC ERROR: SARDINE_ORBIT_CACHE points to a non-existent path: {}. Set SARDINE_ORBIT_CACHE to a valid directory or pass orbit_cache_dir.",
                            p.display()
                        )));
                    }
                } else {
                    return Err(SarError::Processing(
                        "❌ SCIENTIFIC ERROR: Orbit cache directory must be explicitly specified. Set SARDINE_ORBIT_CACHE or pass orbit_cache_dir."
                            .to_string(),
                    ));
                }
            }
        };
        let orbit_manager = OrbitManager::new(cache_dir);
        let primary_orbit_type =
            crate::io::orbit::OrbitReader::determine_orbit_type(metadata.start_time);

        let mut downloaded_files = Vec::new();

        // Try to download primary orbit type only (no silent fallbacks in scientific mode)
        match orbit_manager.download_and_cache_orbit_public(
            &product_id,
            metadata.start_time,
            primary_orbit_type,
        ) {
            Ok(_) => {
                let path = orbit_manager
                    .get_cache_dir()
                    .join(format!("{}_{}.EOF", product_id, primary_orbit_type));
                downloaded_files.push(path);
                log::info!("Downloaded {} orbit file", primary_orbit_type);
            }
            Err(e) => {
                return Err(SarError::Processing(format!(
                    "Failed to download {} orbit file: {}. Scientific mode forbids fallback orbits.",
                    primary_orbit_type, e
                )));
            }
        }

        Ok(downloaded_files)
    }

    /// Get satellite position and velocity for each pixel row in a burst
    /// This is essential for accurate geolocation and Doppler processing
    pub fn get_burst_orbit_data(
        &mut self,
        pol: Polarization,
        orbit_cache_dir: Option<&Path>,
    ) -> SarResult<BurstOrbitData> {
        log::info!("Computing burst orbit data for polarization {}", pol);

        // Get orbit data
        let orbit_data = self.rget_orbit_data(orbit_cache_dir)?;

        // Get annotation data to extract timing information
        let annotation = self.read_annotation(pol)?;

        // Extract burst timing parameters from annotation
        // These would typically come from the annotation XML, but for now use defaults
        let burst_start_time = annotation.start_time;

        // Get image dimensions to determine number of azimuth lines
        let measurements = self.find_measurement_files()?;
        let measurement_file = measurements.get(&pol).ok_or_else(|| {
            SarError::InvalidFormat(format!("No measurement found for polarization {}", pol))
        })?;

        // Get image dimensions from TIFF
        let (width, height) = match self.format {
            ProductFormat::Zip => {
                use tempfile::NamedTempFile;
                let archive = self.open_archive()?;
                let mut zip_file = archive.by_name(measurement_file).map_err(|e| {
                    SarError::Io(std::io::Error::new(
                        std::io::ErrorKind::Other,
                        format!("Failed to access {}: {}", measurement_file, e),
                    ))
                })?;
                let mut temp_file = NamedTempFile::new().map_err(SarError::Io)?;
                std::io::copy(&mut zip_file, &mut temp_file).map_err(SarError::Io)?;
                let dataset = gdal::Dataset::open(temp_file.path()).map_err(|e| {
                    SarError::Io(std::io::Error::new(
                        std::io::ErrorKind::Other,
                        format!("Failed to open TIFF with GDAL: {}", e),
                    ))
                })?;
                let rs = dataset.raster_size();
                (rs.0, rs.1)
            }
            ProductFormat::Safe => {
                let full_path = self.product_path.join(measurement_file);
                let dataset = gdal::Dataset::open(full_path).map_err(|e| {
                    SarError::Io(std::io::Error::new(
                        std::io::ErrorKind::Other,
                        format!("Failed to open TIFF with GDAL: {}", e),
                    ))
                })?;
                let rs = dataset.raster_size();
                (rs.0, rs.1)
            }
        };

        log::info!("Image dimensions: {} x {} pixels", width, height);

        // Calculate azimuth time interval (prefer annotation swath value, then PRF, then duration/height)
        let azimuth_time_interval = {
            // Prefer per-swath azimuthTimeInterval from cached metadata
            let swath_interval = annotation
                .sub_swaths
                .values()
                .find_map(|sw| sw.azimuth_time_interval);

            if let Some(dt) = swath_interval {
                log::info!(
                    "Azimuth timing: using swath azimuthTimeInterval from annotation ({:.9}s)",
                    dt
                );
                dt
            } else if let Some(prf_dt) = annotation.sub_swaths.values().find_map(|sw| {
                sw.prf_hz
                    .and_then(|p| if p > 0.0 { Some(1.0 / p) } else { None })
            }) {
                log::info!(
                    "Azimuth timing: using PRF-derived interval from annotation (1/PRF = {:.9}s)",
                    prf_dt
                );
                prf_dt
            } else {
                // Fallback: derive from product duration and number of azimuth lines
                let acquisition_duration = (annotation.stop_time - annotation.start_time)
                    .num_milliseconds() as f64
                    / 1000.0;
                let dt = acquisition_duration / height as f64;
                log::warn!(
                    "Azimuth timing: falling back to duration/height derivation ({:.9}s); annotation interval/PRF missing",
                    dt
                );
                dt
            }
        };

        // Report which timing source was used for orbit interpolation
        {
            let used = if annotation
                .sub_swaths
                .values()
                .any(|sw| sw.azimuth_time_interval.is_some())
            {
                "annotation_azimuthTimeInterval"
            } else if annotation.sub_swaths.values().any(|sw| sw.prf_hz.is_some()) {
                "annotation_prf"
            } else {
                "fallback_duration"
            };
            log::info!(
                "Azimuth timing source: {} (interval={:.9}s)",
                used,
                azimuth_time_interval
            );
        };

        log::info!(
            "Azimuth timing: duration={:.3}s, interval={:.6}s, lines={}",
            (annotation.stop_time - annotation.start_time).num_milliseconds() as f64 / 1000.0,
            azimuth_time_interval,
            height
        );

        // Interpolate orbit for burst
        crate::io::orbit::OrbitReader::interpolate_burst_orbit(
            &orbit_data,
            burst_start_time,
            azimuth_time_interval,
            height,
        )
    }

    /// Calculate satellite position at a specific pixel (azimuth line, range sample)
    pub fn get_satellite_position_at_pixel(
        &mut self,
        pol: Polarization,
        azimuth_line: usize,
        range_sample: usize,
        orbit_cache_dir: Option<&Path>,
    ) -> SarResult<([f64; 3], [f64; 3])> {
        // Returns (position, velocity)
        let burst_orbit = self.get_burst_orbit_data(pol, orbit_cache_dir)?;

        let position = burst_orbit
            .get_position_at_line(azimuth_line)
            .ok_or_else(|| {
                SarError::Processing(format!(
                    "Azimuth line {} out of range (max: {})",
                    azimuth_line,
                    burst_orbit.num_lines()
                ))
            })?;

        let velocity = burst_orbit
            .get_velocity_at_line(azimuth_line)
            .ok_or_else(|| {
                SarError::Processing(format!(
                    "Azimuth line {} out of range for velocity",
                    azimuth_line
                ))
            })?;

        log::debug!("Satellite at pixel [{}, {}]: pos=[{:.1}, {:.1}, {:.1}] km, vel=[{:.1}, {:.1}, {:.1}] m/s",
                   azimuth_line, range_sample,
                   position[0]/1000.0, position[1]/1000.0, position[2]/1000.0,
                   velocity[0], velocity[1], velocity[2]);

        Ok((position, velocity))
    }

    /// Calculate Doppler centroid frequency for SAR focusing
    ///
    /// **FIX A6:** Uses orbit-reference-relative time when available for consistency
    /// with orbit interpolation. Falls back to product-start-relative time with warning.
    pub fn calculate_doppler_centroid(
        &mut self,
        pol: Polarization,
        azimuth_line: usize,
        _range_sample: usize,
        _orbit_cache_dir: Option<&Path>,
    ) -> SarResult<f64> {
        use slc_reader::time_utils::seconds_since_epoch;

        let annotation_root = self.get_annotation_for_polarization(pol)?;

        // Determine azimuth time interval (s/line)
        let az_time_interval = if let Some(image) = &annotation_root.image_annotation {
            if let Some(info) = &image.image_information {
                if let Some(dt) = info.azimuth_time_interval {
                    dt
                } else {
                    // SCIENTIFIC FIX: Always require azimuthTimeInterval from annotation
                    // 1/PRF fallback is only correct for SM mode, wrong for IW/EW merged products
                    return Err(SarError::Metadata(
                        "azimuthTimeInterval missing in annotation. \
                        This field is required for accurate timing calculations. \
                        Check product integrity or use single-burst SM products."
                            .to_string(),
                    ));
                }
            } else {
                return Err(SarError::Metadata(
                    "Missing imageInformation in annotation".to_string(),
                ));
            }
        } else {
            return Err(SarError::Metadata(
                "Missing imageAnnotation in annotation".to_string(),
            ));
        };

        // FIX A6: Use orbit-relative time when available
        // This ensures consistency with orbit state vector interpolation
        let (eval_time, time_domain) = {
            // Try to get orbit epoch from cached metadata
            let orbit_epoch = self
                .cached_metadata
                .as_ref()
                .and_then(|m| m.orbit_data.as_ref())
                .map(|o| o.reference_time);

            let product_start = self.cached_metadata.as_ref().map(|m| m.start_time);

            match (orbit_epoch, product_start) {
                (Some(epoch), Some(start)) => {
                    // Compute product_start_rel_s (seconds from orbit epoch to product start)
                    let product_start_rel_s = seconds_since_epoch(epoch, start);
                    // t_rel_orbit = product_start_rel_s + line * azimuth_time_interval
                    let t_rel_orbit = product_start_rel_s + azimuth_line as f64 * az_time_interval;
                    (t_rel_orbit, "orbit-relative")
                }
                _ => {
                    // Fallback to product-start-relative (reduced rigor)
                    log::warn!(
                        "⚠️  Doppler evaluation using product-start-relative time (orbit epoch unavailable)"
                    );
                    let t_since_start = azimuth_line as f64 * az_time_interval;
                    (t_since_start, "product-start-relative")
                }
            }
        };

        let doppler_hz = annotation_root.evaluate_doppler_centroid(eval_time)?;
        log::debug!(
            "Doppler at line {} (t={:.6}s {}): {:.2} Hz",
            azimuth_line,
            eval_time,
            time_domain,
            doppler_hz
        );
        Ok(doppler_hz)
    }

    /// Extract IW sub-swaths information from annotation files
    /// This implements the "IW split" step in the SAR processing pipeline
    pub fn extract_iw_subswaths(
        &mut self,
        pol: Polarization,
    ) -> SarResult<HashMap<String, SubSwath>> {
        // Use all-subs discovery and merge subswaths from all files for this polarization
        let all = self.find_all_annotation_files()?;
        let files = all.get(&pol).ok_or_else(|| {
            SarError::InvalidFormat(format!(
                "No annotation files found for polarization {}",
                pol
            ))
        })?;

        let mut result: HashMap<String, SubSwath> = HashMap::new();
        for f in files {
            let xml = self.read_file_as_string(f)?;
            let annotation = crate::io::annotation::parse_annotation_xml(&xml).map_err(|e| {
                SarError::XmlParsing(format!("Failed to parse annotation XML {}: {}", f, e))
            })?;
            let subs = crate::io::annotation::ProductRoot::extract_subswaths(&annotation).map_err(
                |e| {
                    SarError::Processing(format!(
                        "Failed to extract IW sub-swaths from {}: {}",
                        f, e
                    ))
                },
            )?;
            for (k, v) in subs {
                result.insert(k, v);
            }
        }
        Ok(result)
    }

    // Removed fallback string-based IW subswath extraction to enforce scientific parser only.

    // Deprecated: string-based numeric extraction removed for scientific rigor

    /// Find ALL annotation files for all IW subswaths (not just one per polarization)
    pub fn find_all_iw_annotation_files(
        &mut self,
    ) -> SarResult<HashMap<Polarization, Vec<String>>> {
        let files = self.list_files()?;
        let mut annotations: HashMap<Polarization, Vec<String>> = HashMap::new();

        for file in files {
            // Check for annotation directory using path components
            let has_annotation = std::path::Path::new(&file).components().any(|c| {
                c.as_os_str()
                    .to_str()
                    .map(|s| s.to_lowercase() == "annotation")
                    .unwrap_or(false)
            });

            if has_annotation && file.ends_with(".xml") && !file.contains("calibration") {
                // Parse polarization from filename
                let pol = if file.contains("-vv-") {
                    Some(Polarization::VV)
                } else if file.contains("-vh-") {
                    Some(Polarization::VH)
                } else if file.contains("-hv-") {
                    Some(Polarization::HV)
                } else if file.contains("-hh-") {
                    Some(Polarization::HH)
                } else {
                    None
                };

                if let Some(polarization) = pol {
                    annotations.entry(polarization).or_default().push(file);
                }
            }
        }

        // Sort the files for each polarization to ensure consistent order (IW1, IW2, IW3)
        for files in annotations.values_mut() {
            files.sort();
        }

        if annotations.is_empty() {
            return Err(SarError::InvalidFormat(
                "No IW annotation files found".to_string(),
            ));
        }

        Ok(annotations)
    }

    /// Extract ALL IW sub-swaths for a specific polarization (IW1, IW2, IW3)
    pub fn extract_all_iw_subswaths(
        &mut self,
        pol: Polarization,
    ) -> SarResult<HashMap<String, SubSwath>> {
        log::info!(
            "🔍 DEBUGGING: Starting extract_all_iw_subswaths for {:?}",
            pol
        );
        let all_annotations = self.find_all_iw_annotation_files()?;
        log::info!(
            "🔍 DEBUGGING: Found annotation files: {:?}",
            all_annotations
        );
        let annotation_files = all_annotations.get(&pol).ok_or_else(|| {
            SarError::InvalidFormat(format!(
                "No annotation files found for polarization {}",
                pol
            ))
        })?;
        log::info!(
            "🔍 DEBUGGING: Annotation files for {:?}: {:?}",
            pol,
            annotation_files
        );

        let mut all_subswaths = HashMap::new();

        for annotation_file in annotation_files {
            log::info!(
                "🔍 DEBUGGING: Processing annotation file: {}",
                annotation_file
            );

            let xml_content = self.read_file_as_string(annotation_file)?;
            log::info!(
                "🔍 DEBUGGING: Read XML content, length: {}",
                xml_content.len()
            );

            // Extract the IW subswath from this annotation file
            let annotation =
                crate::io::annotation::parse_annotation_xml(&xml_content).map_err(|e| {
                    SarError::XmlParsing(format!("Failed to parse {}: {}", annotation_file, e))
                })?;
            log::info!("🔍 DEBUGGING: Parsed annotation XML successfully");

            log::info!("🔍 DEBUGGING: Extracting subswaths from parsed annotation");
            let subswaths = crate::io::annotation::ProductRoot::extract_subswaths(&annotation)
                .map_err(|e| {
                    SarError::Processing(format!(
                        "Failed to extract IW sub-swaths from {}: {}",
                        annotation_file, e
                    ))
                })?;
            log::info!(
                "🔍 DEBUGGING: extract_subswaths returned {} subswaths",
                subswaths.len()
            );

            // Add all extracted subswaths to our collection
            for (swath_id, subswath) in subswaths {
                log::info!("🔍 DEBUGGING: Adding subswath: {}", swath_id);
                all_subswaths.insert(swath_id, subswath);
            }
        }

        log::info!(
            "🔍 DEBUGGING: extract_all_iw_subswaths returning {} total subswaths",
            all_subswaths.len()
        );
        Ok(all_subswaths)
    }

    /// Get available IW sub-swaths for all polarizations (fixed to extract ALL subswaths)
    ///
    /// **CRITICAL FIX**: This function now returns subswaths with corrected global range positions
    /// derived from slant range time differences. The cached metadata (which contains the
    /// corrected positions) is used when available.
    pub fn get_all_iw_subswaths(
        &mut self,
    ) -> SarResult<HashMap<Polarization, HashMap<String, SubSwath>>> {
        log::info!("Starting get_all_iw_subswaths (with global range position fix)");

        // CRITICAL FIX: Use cached metadata which has the corrected global range positions.
        // The cached metadata calculates global positions from slant range times, ensuring
        // that IW1, IW2, IW3 are properly offset in range direction for merge.
        let cached_meta = self.get_cached_metadata()?;

        // If we have corrected subswaths in cache, use them
        if !cached_meta.sub_swaths.is_empty() {
            log::info!(
                "Using cached subswaths with corrected global range positions: {:?}",
                cached_meta.sub_swaths.keys().collect::<Vec<_>>()
            );

            // Log the corrected positions for verification
            for (name, sw) in &cached_meta.sub_swaths {
                log::info!(
                    "  {}: samples {}..{} (corrected from slant_range_time={:.6}s)",
                    name,
                    sw.first_sample_global,
                    sw.last_sample_global,
                    sw.slant_range_time
                );
            }

            // Build result map: duplicate subswaths for each available polarization
            let mut all_subswaths = HashMap::new();
            for pol in &cached_meta.polarizations {
                let mut pol_subswaths = HashMap::new();
                for (swath_id, subswath) in &cached_meta.sub_swaths {
                    pol_subswaths.insert(swath_id.clone(), subswath.clone());
                }
                all_subswaths.insert(*pol, pol_subswaths);
            }

            log::info!(
                "Returning {} polarizations with corrected subswath positions",
                all_subswaths.len()
            );
            return Ok(all_subswaths);
        }

        // Fallback to original extraction if cache is empty (shouldn't happen with new_with_full_cache)
        log::warn!("Cached subswaths empty, falling back to direct annotation extraction");
        let all_annotations = self.find_all_iw_annotation_files()?;
        log::info!(
            "Found annotation files for {} polarizations",
            all_annotations.len()
        );
        let mut all_subswaths = HashMap::new();
        let strict_mode = std::env::var("SARDINE_STRICT").is_ok();

        for pol in all_annotations.keys().cloned() {
            log::info!("Processing polarization: {:?}", pol);
            match self.extract_all_iw_subswaths(pol) {
                Ok(subswaths) => {
                    if subswaths.is_empty() {
                        log::warn!(
                            "No subswaths extracted for {:?} - annotation may be empty",
                            pol
                        );
                        // In strict mode, fail on empty results
                        if strict_mode {
                            return Err(SarError::Processing(format!(
                                "No subswaths extracted for polarization {:?}",
                                pol
                            )));
                        }
                    }
                    log::info!(
                        "Extracted {} subswaths for polarization {:?}",
                        subswaths.len(),
                        pol
                    );
                    all_subswaths.insert(pol, subswaths);
                }
                Err(e) => {
                    // In strict mode, fail immediately on any extraction error
                    if strict_mode {
                        return Err(SarError::Processing(format!(
                            "Failed to extract sub-swaths for {:?}: {}",
                            pol, e
                        )));
                    }
                    log::error!("Failed to extract sub-swaths for {:?}: {}", pol, e);
                    // Continue with other polarizations in non-strict mode
                }
            }
        }

        // Validate we extracted at least one polarization
        if all_subswaths.is_empty() {
            return Err(SarError::Processing(
                "No subswaths extracted for any polarization. \
                 This may indicate missing or corrupted annotation files."
                    .to_string(),
            ));
        }

        log::info!(
            "Total polarizations with subswaths: {}",
            all_subswaths.len()
        );
        Ok(all_subswaths)
    }

    /// Get subswath names as a vector of strings
    pub fn subswath_names(metadata: &SarMetadata) -> Vec<String> {
        metadata.sub_swaths.keys().cloned().collect()
    }

    /// Check if this is an IW mode SLC product
    pub fn is_iw_mode(&mut self) -> SarResult<bool> {
        // Check if annotation files contain IW mode indicators
        let files = self.list_files()?;

        for file in files {
            // Check for annotation directory using path components
            let has_annotation = std::path::Path::new(&file).components().any(|c| {
                c.as_os_str()
                    .to_str()
                    .map(|s| s.to_lowercase() == "annotation")
                    .unwrap_or(false)
            });

            if has_annotation && file.ends_with(".xml") && file.contains("-iw") {
                return Ok(true);
            }
        }

        Ok(false)
    }

    /// Extract polarization from filename
    fn extract_polarization_from_filename(&self, filename: &str) -> Option<Polarization> {
        if filename.contains("-vv-") {
            Some(Polarization::VV)
        } else if filename.contains("-vh-") {
            Some(Polarization::VH)
        } else if filename.contains("-hv-") {
            Some(Polarization::HV)
        } else if filename.contains("-hh-") {
            Some(Polarization::HH)
        } else {
            None
        }
    }

    /// Apply multilooking to calibrated intensity data
    ///
    /// This method takes calibrated intensity data and applies multilooking
    /// to reduce speckle noise.
    pub fn multilook_intensity(
        &mut self,
        intensity_data: &Array2<f32>,
        pol: Polarization,
        range_looks: usize,
        azimuth_looks: usize,
    ) -> SarResult<(Array2<f32>, f64, f64)> {
        log::info!(
            "Applying multilook: {}x{} looks to intensity data",
            azimuth_looks,
            range_looks
        );

        // Get metadata for pixel spacing
        let metadata = self.read_annotation(pol)?;

        // Extract pixel spacings from metadata
        let (range_spacing, azimuth_spacing) = metadata.pixel_spacing;

        log::info!(
            "Input pixel spacing: range={}m, azimuth={}m",
            range_spacing,
            azimuth_spacing
        );

        // Create multilook processor
        #[allow(deprecated)]
        let params = crate::core::multilook::MultilookParams {
            range_looks,
            azimuth_looks,
            mode: crate::core::multilook::MultilookMode::Intensity,
            preserve_power: true, // ESA/SNAP standard
            border_mode: crate::core::multilook::BorderMode::Partial,
            include_partial: true, // Include partial windows for edge preservation
        };

        let processor = crate::core::multilook::MultilookProcessor::new(params);

        // Apply multilooking
        let (multilooked_data, new_range_spacing, new_azimuth_spacing) =
            processor.apply_multilook_filtered(intensity_data, range_spacing, azimuth_spacing)?;

        log::info!(
            "Multilooking complete: {}x{} -> {}x{}",
            intensity_data.nrows(),
            intensity_data.ncols(),
            multilooked_data.nrows(),
            multilooked_data.ncols()
        );
        log::info!(
            "Output pixel spacing: range={}m, azimuth={}m",
            new_range_spacing,
            new_azimuth_spacing
        );

        Ok((multilooked_data, new_range_spacing, new_azimuth_spacing))
    }

    /// Create a GeoTransform for SAR data based on metadata
    fn create_sar_geotransform(
        &self,
        metadata: &SarMetadata,
        range_spacing: f64,
        azimuth_spacing: f64,
    ) -> SarResult<GeoTransform> {
        // Use latitude-aware meters-per-degree based on WGS84 at scene center
        let bbox = &metadata.bounding_box;
        let mid_lat = (bbox.min_lat + bbox.max_lat) / 2.0;

        // WGS84 ellipsoid parameters
        let a = crate::constants::geodetic::WGS84_SEMI_MAJOR_AXIS_M;
        let e2 = crate::constants::geodetic::WGS84_ECCENTRICITY_SQUARED;

        let lat_rad = mid_lat.to_radians();
        let sin_lat = lat_rad.sin();
        let denom = (1.0 - e2 * sin_lat * sin_lat).sqrt();
        let prime_vertical_radius = a / denom; // N(φ)
        let meridional_radius = a * (1.0 - e2) / (1.0 - e2 * sin_lat * sin_lat).powf(1.5); // M(φ)

        let meters_per_degree_lon =
            prime_vertical_radius * lat_rad.cos() * std::f64::consts::PI / 180.0;
        let meters_per_degree_lat = meridional_radius * std::f64::consts::PI / 180.0;

        // Convert desired spacings (meters) to degrees at mid-latitude
        let pixel_width_deg = if meters_per_degree_lon.is_finite() && meters_per_degree_lon > 0.0 {
            range_spacing / meters_per_degree_lon
        } else {
            // SCIENTIFIC FIX: Always fail on undefined coordinate conversion
            return Err(SarError::Processing(
                "Cannot compute pixel width: meters-per-degree (lon) is undefined at scene center. \
                Check that the scene bounding box has valid latitude coordinates.".to_string(),
            ));
        };

        let pixel_height_deg = if meters_per_degree_lat.is_finite() && meters_per_degree_lat > 0.0 {
            -(azimuth_spacing / meters_per_degree_lat) // negative for north-up
        } else {
            // SCIENTIFIC FIX: Always fail on undefined coordinate conversion
            return Err(SarError::Processing(
                "Cannot compute pixel height: meters-per-degree (lat) is undefined at scene center. \
                Check that the scene bounding box has valid latitude coordinates.".to_string(),
            ));
        };

        Ok(GeoTransform {
            top_left_x: bbox.min_lon,
            pixel_width: pixel_width_deg,
            rotation_x: 0.0,
            top_left_y: bbox.max_lat,
            rotation_y: 0.0,
            pixel_height: pixel_height_deg,
        })
    }

    /// Extract product type from XML or infer from product ID
    fn extract_product_type(xml_content: &str, product_id: &str) -> Option<String> {
        // Try to extract from XML first using serde (handles multi-line tags)
        #[derive(Deserialize)]
        struct ProductTag {
            #[serde(rename = "productType")]
            product_type: Option<String>,
        }

        #[derive(Deserialize)]
        struct ProductWrapper {
            product: Option<ProductTag>,
        }

        if let Ok(wrapper) = from_str::<ProductWrapper>(xml_content) {
            if let Some(product_type) = wrapper.product.and_then(|p| p.product_type) {
                return Some(product_type);
            }
        }

        // Infer from product ID pattern: S1A_IW_SLC__1SDV_...
        if let Some(start) = product_id.find('_') {
            if let Some(second_underscore) = product_id[start + 1..].find('_') {
                let start_pos = start + 1;
                let end_pos = start_pos + second_underscore;
                let mode_and_type = &product_id[start_pos..end_pos];

                // Extract the type part (SLC, GRD, etc.)
                if mode_and_type.len() >= 6 {
                    let product_type = &mode_and_type[3..]; // Skip "IW_" or similar
                    return Some(product_type.to_string());
                }
            }
        }

        None
    }

    /// Parse orbit information from Sentinel-1 product ID
    /// Format: S1A_IW_SLC__1SDV_20200103T170815_20200103T170842_030639_0382D5_DADE
    /// Returns: (orbit_number, relative_orbit, orbit_direction)
    pub(crate) fn parse_orbit_from_product_id(product_id: &str) -> Option<(u32, u32, String)> {
        let parts: Vec<&str> = product_id.split('_').collect();
        // Debug: log the product ID parsing (info -> debug to reduce noise)
        log::debug!(
            "Parsing orbit from product ID: {} ({} segments)",
            product_id,
            parts.len()
        );

        // Typical segmentation (note the empty field after the double underscore before the level/class field):
        // S1A_IW_SLC__1SDV_20200103T170815_20200103T170842_030639_0382D5_DADE
        //  0  1  2   3    4              5               6      7       8
        // 0:S1A 1:IW 2:SLC 3:"" 4:1SDV 5:START_TIME 6:STOP_TIME 7:ABS_ORB 8:MISSION_DATA 9:UNIQUE_ID
        // Some variants (e.g. GRD) may not contain the double underscore; we therefore *search* for the
        // first purely numeric 6-digit field that plausibly represents the absolute orbit number.

        // Find first part that looks like a 6-digit absolute orbit number
        let mut abs_orbit_idx = None;
        for (i, p) in parts.iter().enumerate() {
            if p.len() == 6 && p.chars().all(|c| c.is_ascii_digit()) {
                abs_orbit_idx = Some(i);
                break;
            }
        }

        let idx = match abs_orbit_idx {
            Some(i) => i,
            None => {
                log::warn!(
                    "No 6-digit numeric segment found in product ID '{}'",
                    product_id
                );
                return None;
            }
        };

        let orbit_number: u32 = match parts[idx].parse() {
            Ok(v) => v,
            Err(e) => {
                log::warn!(
                    "Failed to parse numeric orbit field '{}' at idx {}: {}",
                    parts[idx],
                    idx,
                    e
                );
                return None;
            }
        };

        // Sentinel-1 relative orbit formula (ESA doc): (ABS_ORBIT - 73) mod 175 + 1
        // (The previous implementation used (orbit_number - 1) which shifts all relative orbits by 72.)
        let relative_orbit = ((orbit_number as i64 - 73).rem_euclid(175) as u32) + 1;

        // Orbit direction: Sentinel-1 uses sun-synchronous orbit with 175-orbit repeat cycle.
        // Ascending passes occur when orbit_number mod 175 falls in certain ranges.
        // More accurate: check z-velocity from orbit state vectors (positive = ascending).
        // Fallback heuristic: relative orbits 1-87 are typically ascending, 88-175 descending.
        // This is approximate but more accurate than simple parity.
        let orbit_direction = if relative_orbit <= 87 {
            "ASCENDING"
        } else {
            "DESCENDING"
        };

        log::debug!(
            "Orbit parse result: abs_orbit={} (idx {}), rel_orbit={}, dir={}",
            orbit_number,
            idx,
            relative_orbit,
            orbit_direction
        );
        Some((orbit_number, relative_orbit, orbit_direction.to_string()))
    }

    /// Extract incidence angle range from the array of values
    fn extract_incidence_angle_range(xml_content: &str) -> Option<(f64, f64)> {
        let start_tag = "<incidenceAngle count=";

        if let Some(start_pos) = xml_content.find(start_tag) {
            if let Some(end_bracket) = xml_content[start_pos..].find('>') {
                let content_start = start_pos + end_bracket + 1;
                if let Some(end_tag_pos) = xml_content[content_start..].find("</incidenceAngle>") {
                    let angles_str = &xml_content[content_start..content_start + end_tag_pos];
                    let angles: Vec<f64> = angles_str
                        .split_whitespace()
                        .filter_map(|s| s.parse().ok())
                        .collect();

                    if !angles.is_empty() {
                        let near = angles[0];
                        let far = angles[angles.len() - 1];
                        return Some((near, far));
                    }
                }
            }
        }
        None
    }

    /// Extract calibration constant from annotation XML or calibration files
    fn extract_calibration_constant(xml_content: &str) -> Option<f64> {
        // First try to find global calibration constants in annotation XML
        if let Some(value) =
            Self::extract_real_value_f64(xml_content, "absoluteCalibrationConstant")
                .or_else(|| Self::extract_real_value_f64(xml_content, "calibrationConstant"))
                .or_else(|| Self::extract_real_value_f64(xml_content, "rescalingFactor"))
        {
            return Some(value);
        }

        // If not found in annotation, return a default reference constant for SLC data
        // Sentinel-1 SLC data uses digital numbers with no specific calibration constant
        // Return 1.0 as SLC data is already calibrated relative to amplitude
        Some(1.0)
    }

    /// Read annotation file directly without full metadata parsing
    fn read_annotation_file_raw(&mut self, annotation_file: &str) -> SarResult<SarMetadata> {
        let xml_content = self.read_file_as_string(annotation_file)?;

        // Extract actual polarization from filename instead of hardcoding VV
        let polarization = Self::extract_polarization(annotation_file).unwrap_or_else(|| {
            log::warn!(
                "Could not extract polarization from filename '{}', defaulting to VV",
                annotation_file
            );
            Polarization::VV
        });

        Self::parse_annotation_xml_comprehensive(&xml_content, polarization, &self.product_path)
    }

    // ============================================================================
    // COMPREHENSIVE CACHING SYSTEM - Phase 1 Implementation
    // ============================================================================

    /// Create a new SlcReader with immediate comprehensive metadata caching
    /// This eliminates all redundant metadata extraction during processing
    pub fn new_with_full_cache<P: AsRef<Path>>(product_path: P) -> SarResult<Self> {
        let mut reader = Self::new_internal(product_path)?;

        log::info!("🚀 Initializing comprehensive metadata cache for optimal performance");
        reader.initialize_all_caches()?;

        Ok(reader)
    }

    /// Initialize ALL metadata caches for maximum performance
    /// This method extracts and caches ALL metadata once during initialization
    pub fn initialize_all_caches(&mut self) -> SarResult<()> {
        if self.cache_initialized {
            log::debug!("Cache already initialized, skipping");
            return Ok(());
        }

        log::info!("🔄 Extracting and caching ALL metadata for efficient processing");

        // Step 1: Cache all annotation files and their content
        // Use ALL annotation files per polarization (IW1/2/3), not just first
        let all_annotation_files = self.find_all_annotation_files()?;
        // Cross-check measurement coverage vs annotations to fail fast on incomplete pols
        if let Ok(measurements) = self.find_measurement_files() {
            let missing: Vec<_> = measurements
                .keys()
                .filter(|pol| !all_annotation_files.contains_key(pol))
                .copied()
                .collect();
            if !missing.is_empty() {
                return Err(SarError::Processing(format!(
                    "Measurement files present for {:?} but no matching annotations found",
                    missing
                )));
            }
        }
        // Also keep first-per-pol for legacy expectations
        let mut first_annotation_per_pol: HashMap<Polarization, String> = HashMap::new();
        for (pol, files) in &all_annotation_files {
            if let Some(first) = files.first() {
                first_annotation_per_pol.insert(*pol, first.clone());
            }
        }

        log::info!(
            "📁 Found {} polarizations to cache",
            all_annotation_files.len()
        );

        // I/O Optimization: Parallel annotation parsing for better performance
        // Read all files first (I/O), then parse XML in parallel (CPU-bound)
        let mut annotation_data: Vec<(Polarization, String, String)> = Vec::new();
        for (pol, files) in &all_annotation_files {
            for annotation_path in files {
                match self.read_file_as_string(annotation_path) {
                    Ok(xml_content) => {
                        annotation_data.push((*pol, annotation_path.clone(), xml_content));
                    }
                    Err(e) => {
                        log::warn!("⚠️  Failed to read annotation {}: {}", annotation_path, e);
                    }
                }
            }
        }

        log::info!(
            "🔄 Parsing {} annotation files in parallel",
            annotation_data.len()
        );

        // Parse XML in parallel (CPU-bound and thread-safe)
        use rayon::prelude::*;
        let parsed_results: Vec<(
            Polarization,
            String,
            SarResult<Arc<crate::io::annotation::ProductRoot>>,
            String,
        )> = annotation_data
            .into_par_iter()
            .map(|(pol, annotation_path, xml_content)| {
                // Parse XML (CPU-bound - safe to parallelize)
                let parse_result = self
                    .parse_annotation_root_xml(&xml_content)
                    .map(|root| Arc::new(root));
                (pol, annotation_path, parse_result, xml_content)
            })
            .collect();

        // Store parsed results sequentially (to avoid race conditions)
        for (pol, annotation_path, parse_result, xml_content) in parsed_results {
            match parse_result {
                Ok(annotation_root) => {
                    // Cache XML content
                    self.cached_xml_content
                        .insert(annotation_path.clone(), xml_content);
                    // Store parsed annotation
                    self.cached_annotations
                        .entry(pol)
                        .or_insert_with(Vec::new)
                        .push(annotation_root);
                    log::debug!("✅ Cached annotation for {:?}: {}", pol, annotation_path);
                }
                Err(e) => {
                    log::warn!("⚠️  Failed to parse annotation {}: {}", annotation_path, e);
                }
            }
        }

        // Cache calibration data (sequential - depends on parsed annotations)
        for pol in all_annotation_files.keys() {
            if let Ok(calibration) = self.read_calibration_data(*pol) {
                self.cached_calibration.insert(*pol, calibration);
                log::debug!("✅ Cached calibration for {:?}", pol);
            } else {
                log::warn!("⚠️  Could not cache calibration for {:?}", pol);
            }
        }

        // Early validation: burst timing continuity and DC model sanity, plus LUT presence
        let noise_files = self.find_noise_files().unwrap_or_default();
        for (pol, roots) in &self.cached_annotations {
            for root in roots {
                Self::validate_burst_timing_continuity(*pol, root.as_ref())?;
                root.validate_doppler_model_against_bursts(*pol);
            }
            self.summarize_lut_presence(
                *pol,
                self.cached_calibration.get(pol),
                noise_files.get(pol).is_some(),
            )?;
        }

        // Get first polarization for metadata initialization
        if first_annotation_per_pol.is_empty() {
            return Err(SarError::Processing(
                "No annotation files found for metadata cache initialization".to_string(),
            ));
        }

        // Step 2: Build comprehensive metadata directly from cached annotation roots
        // Gather metadata for every cached annotation (all polarizations and subswaths)
        let product_id = Self::product_id_from_path(&self.product_path);
        let mut metas = Vec::new();
        for (pol, roots) in &self.cached_annotations {
            for root in roots {
                metas.push(Self::metadata_from_annotation_root(
                    root.as_ref(),
                    *pol,
                    &product_id,
                )?);
            }
        }

        let mut comprehensive_metadata = Self::merge_metadata(metas)?;

        let measurement_geometry_checks_enabled = match self.list_files() {
            Ok(files) => {
                let has_measurements = files.iter().any(|f| Self::is_measurement_tiff(f));
                if !has_measurements {
                    log::warn!(
                        "Skipping measurement geometry verification: no measurement TIFFs discovered"
                    );
                }
                has_measurements
            }
            Err(err) => {
                log::warn!(
                    "Skipping measurement geometry verification: unable to enumerate files ({})",
                    err
                );
                false
            }
        };

        // SCIENTIFIC FIX: Extract subswaths from ALL cached annotations, not just first polarization
        let mut all_subswaths = HashMap::new();
        log::info!("🔍 Extracting subswaths from all cached annotation files");
        let mut extracted_subswaths: Vec<(Polarization, String, SubSwath)> = Vec::new();
        for (pol, roots) in &self.cached_annotations {
            for root in roots {
                match crate::io::annotation::ProductRoot::extract_subswaths(root.as_ref()) {
                    Ok(pol_subswaths) => {
                        for (swath_id, swath_data) in pol_subswaths {
                            extracted_subswaths.push((*pol, swath_id, swath_data));
                        }
                    }
                    Err(e) => {
                        log::warn!(
                            "Failed to extract subswaths from cached annotation for {:?}: {}",
                            pol,
                            e
                        );
                    }
                }
            }
        }

        for (pol, swath_id, swath_data) in extracted_subswaths.into_iter() {
            if measurement_geometry_checks_enabled {
                self.verify_measurement_geometry(pol, &swath_id, &swath_data)?;
            }
            log::info!("✅ Found subswath {} from {:?}", swath_id, pol);

            // Preserve existing metadata when merging across polarizations.
            // In particular, keep any dc_polynomial_t0 we already captured
            // so a later polarization without t0 cannot overwrite it.
            use std::collections::hash_map::Entry;
            match all_subswaths.entry(swath_id) {
                Entry::Vacant(v) => {
                    v.insert(swath_data);
                }
                Entry::Occupied(mut entry) => {
                    let existing = entry.get_mut();

                    // Only adopt new values when the existing metadata is missing.
                    if existing.dc_polynomial.is_none() && swath_data.dc_polynomial.is_some() {
                        existing.dc_polynomial = swath_data.dc_polynomial;
                    }

                    if existing.dc_polynomial_t0.is_none()
                        && swath_data.dc_polynomial_t0.is_some()
                    {
                        existing.dc_polynomial_t0 = swath_data.dc_polynomial_t0;
                    }

                    if existing.azimuth_time_interval.is_none()
                        && swath_data.azimuth_time_interval.is_some()
                    {
                        existing.azimuth_time_interval = swath_data.azimuth_time_interval;
                    }

                    if existing.fm_rate_estimates.is_none()
                        && swath_data.fm_rate_estimates.is_some()
                    {
                        existing.fm_rate_estimates = swath_data.fm_rate_estimates;
                    }
                }
            }
        }

        log::info!(
            "✅ Combined subswaths from all annotations: {:?}",
            all_subswaths.keys().collect::<Vec<_>>()
        );

        comprehensive_metadata.sub_swaths = all_subswaths;

        // CRITICAL FIX: Calculate correct global range positions from slant range times
        // This must be done AFTER combining all subswaths from multiple annotation files
        if comprehensive_metadata.sub_swaths.len() > 1 {
            // Find reference subswath (IW1) with earliest slant range time
            let mut sorted_swaths: Vec<_> = comprehensive_metadata.sub_swaths.iter().collect();
            sorted_swaths.sort_by(|a, b| {
                a.1.slant_range_time
                    .partial_cmp(&b.1.slant_range_time)
                    .unwrap()
            });

            if let Some((ref_name, ref_swath)) = sorted_swaths.first() {
                let ref_tau = ref_swath.slant_range_time;
                let ref_spacing = ref_swath.range_pixel_spacing;

                use crate::constants::physical::SPEED_OF_LIGHT_M_S;

                log::info!("🔧 TOPSAR MERGE FIX: Calculating global range positions from slant range times");
                log::info!(
                    "   Reference: {} @ τ₀={:.9}s, spacing={:.3}m",
                    ref_name,
                    ref_tau,
                    ref_spacing
                );

                // Update all subswaths with corrected global positions
                for (swath_name, swath) in comprehensive_metadata.sub_swaths.iter_mut() {
                    let delta_tau = swath.slant_range_time - ref_tau;
                    let delta_range_m = delta_tau * (SPEED_OF_LIGHT_M_S / 2.0);
                    let offset_samples = (delta_range_m / ref_spacing).round() as usize;

                    // Preserve azimuth positions (they're correct from burst geometry)
                    let azimuth_start = swath.first_line_global;
                    let azimuth_end = swath.last_line_global;

                    // FIX: Use EXCLUSIVE last_sample_global = first + width (one past last valid index)
                    // This matches the documented convention in types.rs
                    // Only update if range_samples is valid (> 0)
                    let prev_first = swath.first_sample_global;
                    let prev_last = swath.last_sample_global;

                    if swath.range_samples > 0 {
                        swath.first_sample_global = offset_samples;
                        // EXCLUSIVE upper bound: last = first + count (one past the last valid sample)
                        swath.last_sample_global = offset_samples + swath.range_samples;

                        // CRITICAL FIX: Also update valid_first_sample and valid_last_sample
                        // to match the new global coordinate system. These were initially set
                        // to local coordinates (0..range_samples) but must be updated when
                        // global positions are corrected for TOPSAR merge.
                        swath.valid_first_sample = Some(swath.first_sample_global);
                        swath.valid_last_sample = Some(swath.last_sample_global); // both now exclusive

                        log::info!(
                            "   {} @ τ={:.9}s: Δτ={:.6}ms → {:.1}m → offset={} samples",
                            swath_name,
                            swath.slant_range_time,
                            delta_tau * 1000.0,
                            delta_range_m,
                            offset_samples
                        );
                        log::info!(
                            "      Global coords: lines {}..{}, samples {}..{} (was {}..{}, width={})",
                            azimuth_start,
                            azimuth_end,
                            swath.first_sample_global,
                            swath.last_sample_global,
                            prev_first,
                            prev_last,
                            swath.range_samples
                        );
                    } else {
                        log::warn!(
                            "⚠️  {} has range_samples=0; preserving existing range window {}..{}",
                            swath_name,
                            prev_first,
                            prev_last
                        );
                    }
                }

                log::info!("✅ Global range positions corrected for TOPSAR merge");
                
                // SCIENTIFIC VALIDATION (Jan 2026): Check for overlaps/gaps between PHYSICALLY ADJACENT swaths
                // CRITICAL: Must sort by range position (first_sample_global) to check actual neighbors
                let mut swaths_sorted: Vec<(&String, &crate::types::SubSwath)> = 
                    comprehensive_metadata.sub_swaths.iter().collect();
                swaths_sorted.sort_by_key(|(_, sw)| sw.first_sample_global);
                
                for i in 0..swaths_sorted.len().saturating_sub(1) {
                    let (curr_name, curr) = swaths_sorted[i];
                    let (next_name, next) = swaths_sorted[i + 1];
                    
                    {
                        let curr_end = curr.last_sample_global;
                        let next_start = next.first_sample_global;
                        
                        if next_start < curr_end {
                            // Overlap detected
                            let overlap_samples = curr_end - next_start;
                            let overlap_m = overlap_samples as f64 * ref_spacing;
                            let overlap_fraction = overlap_samples as f64 / curr.range_samples as f64;
                            
                            log::info!(
                                "🔍 TOPSAR range overlap: {}↔{} = {} samples ({:.1}m, {:.1}% of {})",
                                curr_name, next_name, overlap_samples, overlap_m,
                                overlap_fraction * 100.0, curr_name
                            );
                            
                            // TOPSAR IW typically has 10-20% range overlap for interferometry
                            if overlap_fraction < 0.05 || overlap_fraction > 0.30 {
                                log::warn!(
                                    "⚠️  SCIENTIFIC WARNING: Overlap fraction {:.1}% is outside typical TOPSAR IW range (10-20%). \
                                     Verify slant_range_time values and merge calculation.",
                                    overlap_fraction * 100.0
                                );
                            } else {
                                log::info!("✅ Overlap fraction {:.1}% within expected TOPSAR IW range", overlap_fraction * 100.0);
                            }
                        } else if next_start > curr_end {
                            // Gap detected - this is physically wrong for TOPSAR
                            let gap_samples = next_start - curr_end;
                            let gap_m = gap_samples as f64 * ref_spacing;
                            log::error!(
                                "❌ CRITICAL ERROR: Gap detected between {} (ends at {}) and {} (starts at {}): {} samples ({:.1}m). \
                                 TOPSAR swaths must overlap or be contiguous. Check slant_range_time extraction.",
                                curr_name, curr_end, next_name, next_start, gap_samples, gap_m
                            );
                        } else {
                            log::info!("✅ Swaths {} and {} are perfectly contiguous (no overlap, no gap)", curr_name, next_name);
                        }
                    }
                }

                // CRITICAL: Also update burst_records to match the corrected subswath positions
                // Each burst in a subswath shares the same range extent as the subswath
                log::info!("🔧 Updating burst_records with corrected global range positions");
                for burst in comprehensive_metadata.burst_records.iter_mut() {
                    if let Some(subswath) =
                        comprehensive_metadata.sub_swaths.get(&burst.subswath_id)
                    {
                        let old_start = burst.start_sample_global;
                        let old_end = burst.end_sample_global;
                        burst.start_sample_global = subswath.first_sample_global;
                        burst.end_sample_global = subswath.last_sample_global;
                        log::debug!(
                            "   {} burst {}: samples {}..{} → {}..{}",
                            burst.subswath_id,
                            burst.burst_index,
                            old_start,
                            old_end,
                            burst.start_sample_global,
                            burst.end_sample_global
                        );
                    }
                }
                log::info!("✅ Burst records updated with corrected range positions");
            }
        }

        Self::annotate_fm_rate_relative_times(
            &mut comprehensive_metadata.sub_swaths,
            comprehensive_metadata
                .orbit_data
                .as_ref()
                .map(|orbit| orbit.reference_time),
        )?;

        // Aggregate polarization list across all cached annotations
        let mut all_polarizations: Vec<Polarization> =
            all_annotation_files.keys().copied().collect();
        all_polarizations.sort_by_key(|p| p.to_string());
        all_polarizations.dedup();
        comprehensive_metadata.polarizations = all_polarizations;

        self.cached_metadata = Some(comprehensive_metadata);

        // Mark cache as initialized
        self.cache_initialized = true;

        log::info!("✅ Comprehensive metadata cache initialized successfully");
        log::info!(
            "   📊 Cached {} polarizations",
            first_annotation_per_pol.len()
        );
        log::info!("   📊 Cached {} XML files", self.cached_xml_content.len());
        log::info!(
            "   📊 Cached {} calibration datasets",
            self.cached_calibration.len()
        );

        Ok(())
    }

    /// Get cached metadata (replaces all get_metadata() calls)
    /// This eliminates redundant metadata extraction throughout processing
    pub fn get_cached_metadata(&self) -> SarResult<&SarMetadata> {
        self.ensure_cache_initialized()?;

        self.cached_metadata.as_ref().ok_or_else(|| {
            SarError::Processing(
                "Metadata not cached. Call initialize_all_caches() (construct with SlcReader::new())"
                    .to_string(),
            )
        })
    }

    /// Get cached annotation for specific polarization (no XML re-parsing)
    /// Returns all subswaths already held in memory. Requires cache to be initialized.
    pub fn get_all_cached_annotations(
        &self,
        pol: Polarization,
    ) -> SarResult<Vec<Arc<crate::io::annotation::ProductRoot>>> {
        self.ensure_cache_initialized()?;

        self.cached_annotations
            .get(&pol)
            .map(|annotations| annotations.iter().cloned().collect())
            .ok_or_else(|| {
                SarError::Processing(format!(
                    "Annotations not cached for {:?}. Available: {:?}",
                    pol,
                    self.cached_annotations.keys().collect::<Vec<_>>()
                ))
            })
    }

    /// Get cached calibration coefficients for a polarization
    pub fn get_cached_calibration(&self, pol: Polarization) -> SarResult<&CalibrationCoefficients> {
        self.ensure_cache_initialized()?;

        self.cached_calibration.get(&pol).ok_or_else(|| {
            SarError::Processing(format!(
                "Calibration not cached for {:?}. Available: {:?}",
                pol,
                self.cached_calibration.keys().collect::<Vec<_>>()
            ))
        })
    }

    /// Get cached XML content for specific file (eliminates file re-reading)
    pub fn get_cached_xml_content(&self, file_path: &str) -> SarResult<&String> {
        self.ensure_cache_initialized()?;

        self.cached_xml_content.get(file_path).ok_or_else(|| {
            SarError::Processing(format!("XML content not cached for {}", file_path))
        })
    }

    /// Get all available polarizations from cache
    pub fn get_available_polarizations(&self) -> SarResult<Vec<Polarization>> {
        self.ensure_cache_initialized()?;

        Ok(self.cached_annotations.keys().copied().collect())
    }

    /// Convert cached metadata to legacy HashMap format for compatibility
    pub fn get_cached_metadata_as_map(&self) -> SarResult<HashMap<String, String>> {
        let metadata = self.get_cached_metadata()?;

        let mut map = HashMap::new();

        // Convert SarMetadata to legacy string map format
        map.insert("product_id".to_string(), metadata.product_id.clone());
        map.insert("mission".to_string(), metadata.mission.clone());
        map.insert("platform".to_string(), metadata.platform.clone());
        map.insert("start_time".to_string(), metadata.start_time.to_rfc3339());
        map.insert("stop_time".to_string(), metadata.stop_time.to_rfc3339());

        // Add bounding box coordinates
        map.insert(
            "min_latitude".to_string(),
            metadata.bounding_box.min_lat.to_string(),
        );
        map.insert(
            "max_latitude".to_string(),
            metadata.bounding_box.max_lat.to_string(),
        );
        map.insert(
            "min_longitude".to_string(),
            metadata.bounding_box.min_lon.to_string(),
        );
        map.insert(
            "max_longitude".to_string(),
            metadata.bounding_box.max_lon.to_string(),
        );

        // Add radar parameters (extracted from annotation, not hardcoded)
        if let Some(radar_freq) = metadata.radar_frequency {
            map.insert("radar_frequency".to_string(), radar_freq.to_string());
        }

        if let Some(wavelength) = metadata.wavelength {
            map.insert("wavelength".to_string(), wavelength.to_string());
        }

        // CRITICAL: Add slant_range_time for terrain correction
        if let Some(slant_range_time) = metadata.slant_range_time {
            map.insert("slant_range_time".to_string(), slant_range_time.to_string());
        }

        // CRITICAL: Add range_sampling_rate for processing
        if let Some(range_sampling_rate) = metadata.range_sampling_rate {
            map.insert(
                "range_sampling_rate".to_string(),
                range_sampling_rate.to_string(),
            );
        }

        // CRITICAL: Add PRF for terrain correction (per subswath, then global if consistent)
        log::info!("Checking PRF in SarMetadata struct: {:?}", metadata.prf);

        let per_subswath_prf: HashMap<String, f64> = metadata
            .sub_swaths
            .iter()
            .filter_map(|(name, sw)| sw.prf_hz.map(|prf| (name.clone(), prf)))
            .collect();

        if !per_subswath_prf.is_empty() {
            if let Ok(json) = serde_json::to_string(&per_subswath_prf) {
                map.insert("prf_by_subswath".to_string(), json);
            }
        }

        let mut unique_prfs: Vec<f64> = per_subswath_prf.values().copied().collect();
        unique_prfs.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
        unique_prfs.dedup_by(|a, b| {
            let diff = (*a - *b).abs();
            let max_mag = a.abs().max(b.abs()).max(1.0);
            diff <= 1e-6 * max_mag || diff <= 1e-3
        });

        let global_prf = metadata.prf.or_else(|| {
            if unique_prfs.len() == 1 {
                unique_prfs.first().copied()
            } else if unique_prfs.len() > 1 {
                // IW mode: subswaths have different PRFs. Use median as representative
                // value for terrain correction timing. The exact per-subswath PRFs are
                // still available via prf_by_subswath for code that needs them.
                let mid = unique_prfs.len() / 2;
                let median_prf = unique_prfs[mid];
                log::info!(
                    "PRF differs across subswaths ({:?}); using median {:.3} Hz as global PRF",
                    unique_prfs, median_prf
                );
                Some(median_prf)
            } else {
                None
            }
        });

        if let Some(prf) = global_prf {
            log::info!("Adding PRF to metadata cache: {} Hz", prf);
            map.insert("prf".to_string(), prf.to_string());

            // CRITICAL FIX: Export azimuth_time_interval for terrain correction
            // azimuth_time_interval = 1/PRF (seconds per azimuth line)
            // Also try to get from subswath annotation first (more accurate)
            let azimuth_time_interval = metadata
                .sub_swaths
                .values()
                .find_map(|sw| sw.azimuth_time_interval)
                .unwrap_or(1.0 / prf);
            map.insert(
                "azimuth_time_interval".to_string(),
                azimuth_time_interval.to_string(),
            );
            log::info!(
                "Adding azimuth_time_interval to metadata cache: {:.9}s",
                azimuth_time_interval
            );
        } else {
            log::warn!("PRF is None in cached SarMetadata - extraction failed!");
        }

        let per_subswath_ati: HashMap<String, f64> = metadata
            .sub_swaths
            .iter()
            .filter_map(|(name, sw)| sw.azimuth_time_interval.map(|dt| (name.clone(), dt)))
            .collect();
        if !per_subswath_ati.is_empty() {
            if let Ok(json) = serde_json::to_string(&per_subswath_ati) {
                map.insert("azimuth_time_interval_by_subswath".to_string(), json);
            }
        }

        // CRITICAL FIX: Add absolute product timing for correct terrain correction
        // Convert start/stop times (DateTime<Utc>) to Unix timestamp in seconds
        let product_start_time_abs = metadata.start_time.timestamp() as f64
            + (metadata.start_time.timestamp_subsec_nanos() as f64) * 1e-9;
        map.insert(
            "product_start_time_abs".to_string(),
            product_start_time_abs.to_string(),
        );

        let product_stop_time_abs = metadata.stop_time.timestamp() as f64
            + (metadata.stop_time.timestamp_subsec_nanos() as f64) * 1e-9;
        map.insert(
            "product_stop_time_abs".to_string(),
            product_stop_time_abs.to_string(),
        );

        let product_duration = (product_stop_time_abs - product_start_time_abs).max(0.0);
        map.insert("product_duration".to_string(), product_duration.to_string());

        log::info!(
            "🔬 SCIENTIFIC FIX: Added product timing: start={:.6}s stop={:.6}s duration={:.3}s",
            product_start_time_abs,
            product_stop_time_abs,
            product_duration
        );

        if !metadata.burst_records.is_empty() {
            match serde_json::to_string(&metadata.burst_records) {
                Ok(json) => {
                    map.insert("burst_timing_json".to_string(), json);
                }
                Err(err) => {
                    log::warn!(
                        "⚠️  Unable to serialize burst timing metadata to JSON: {}",
                        err
                    );
                }
            }
        }

        // Orbit diagnostics for downstream QA and STEP 1 validation (always export counts)
        if let Some(orbit) = &metadata.orbit_data {
            let vector_count = orbit.state_vectors.len();
            map.insert("orbit_vectors_count".to_string(), vector_count.to_string());

            if vector_count < 10 {
                log::warn!(
                    "⚠️  Only {} orbit vectors available; terrain correction may be unstable",
                    vector_count
                );
                if crate::types::strict_mode() {
                    return Err(SarError::OrbitError(format!(
                        "Insufficient orbit vectors ({}) for scientific processing",
                        vector_count
                    )));
                }
            }

            if let (Some(first), Some(last)) =
                (orbit.state_vectors.first(), orbit.state_vectors.last())
            {
                let to_secs = |dt: chrono::DateTime<chrono::Utc>| {
                    dt.timestamp() as f64 + (dt.timestamp_subsec_nanos() as f64) * 1e-9
                };
                let first_sec = to_secs(first.time);
                let last_sec = to_secs(last.time);
                map.insert("orbit_first_vector_time".to_string(), first_sec.to_string());
                map.insert("orbit_last_vector_time".to_string(), last_sec.to_string());
                map.insert(
                    "orbit_span_seconds".to_string(),
                    (last_sec - first_sec).max(0.0).to_string(),
                );
            }
        } else {
            map.insert("orbit_vectors_count".to_string(), "0".to_string());
            map.insert("orbit_first_vector_time".to_string(), "0".to_string());
            map.insert("orbit_last_vector_time".to_string(), "0".to_string());
            map.insert("orbit_span_seconds".to_string(), "0".to_string());
        }

        // -----------------------------------------------------------------
        // EXPERT FIX: Derive realistic azimuth bounds and product duration
        // from finalized burst records (annotation-derived). Many SAFE
        // files contain more conservative start/stop times (~25s window)
        // which are larger than the actual acquired burst span. Using the
        // per-subswhath valid first/last lines (when available) yields a
        // tighter native azimuth span. Compute total_azimuth_lines from
        // those bounds and override product_duration = span / PRF + guard.
        // -----------------------------------------------------------------
        if !metadata.sub_swaths.is_empty() {
            let mut min_first_line: Option<usize> = None;
            let mut max_last_line: Option<usize> = None;

            for (_name, sw) in metadata.sub_swaths.iter() {
                let first = sw.valid_first_line.unwrap_or(sw.first_line_global);
                let last = sw.valid_last_line.unwrap_or(sw.last_line_global);

                min_first_line = Some(match min_first_line {
                    Some(prev) => std::cmp::min(prev, first),
                    None => first,
                });
                max_last_line = Some(match max_last_line {
                    Some(prev) => std::cmp::max(prev, last),
                    None => last,
                });
            }

            if let (Some(min_l), Some(max_l)) = (min_first_line, max_last_line) {
                if max_l >= min_l {
                    let total_azimuth_lines = max_l.saturating_sub(min_l).saturating_add(1);
                    map.insert(
                        "total_azimuth_lines".to_string(),
                        total_azimuth_lines.to_string(),
                    );

                    // NOTE: Do NOT override product_duration from lines here!
                    // The product_duration must be the FULL time window (start_time to stop_time)
                    // not the per-subswath line count / PRF, which gives incorrect (short) values.
                    // The correct duration was already set above from product_stop_time_abs - product_start_time_abs.
                    log::info!(
                        "📐 total_azimuth_lines={} (subswath bounds); product_duration={:.3}s (from manifest times)",
                        total_azimuth_lines,
                        product_duration
                    );
                }
            }
        }

        if let Some(dc) = &metadata.doppler_centroid {
            map.insert("doppler_centroid_t0".to_string(), dc.t0.to_string());
            match serde_json::to_string(&dc.coefficients) {
                Ok(json) => {
                    map.insert("doppler_centroid_coeffs".to_string(), json);
                }
                Err(err) => {
                    log::warn!(
                        "⚠️  Unable to serialize doppler centroid coefficients: {}",
                        err
                    );
                }
            }
        }

        if let Some(orbit_data) = &metadata.orbit_data {
            let orbit_ref_epoch_utc = orbit_data.reference_time.timestamp() as f64
                + (orbit_data.reference_time.timestamp_subsec_nanos() as f64) * 1e-9;
            map.insert(
                "orbit_ref_epoch_utc".to_string(),
                orbit_ref_epoch_utc.to_string(),
            );

            let product_start_rel_s = product_start_time_abs - orbit_ref_epoch_utc;
            map.insert(
                "product_start_rel_s".to_string(),
                product_start_rel_s.to_string(),
            );

            log::info!(
                "🛰️  Preserving orbit timing: orbit_ref_epoch_utc={:.6}s, product_start_rel_s={:.6}s",
                orbit_ref_epoch_utc,
                product_start_rel_s
            );
        } else {
            log::warn!(
                "⚠️  Orbit data unavailable in cached metadata; product_start_rel_s will require fallback"
            );
        }

        // Add pixel spacing - the SLC annotation values are the NATIVE spacings
        // before any multilooking. Export both names for clarity.
        map.insert(
            "range_pixel_spacing".to_string(),
            metadata.pixel_spacing.0.to_string(),
        );
        map.insert(
            "native_range_pixel_spacing".to_string(),
            metadata.pixel_spacing.0.to_string(),
        );
        map.insert(
            "azimuth_pixel_spacing".to_string(),
            metadata.pixel_spacing.1.to_string(),
        );
        map.insert(
            "native_azimuth_pixel_spacing".to_string(),
            metadata.pixel_spacing.1.to_string(),
        );

        if let Some(near) = metadata.incidence_angle_near {
            map.insert("incidence_angle_near".to_string(), near.to_string());
        } else {
            log::warn!("Incidence angle (near) missing from cached metadata");
        }

        if let Some(far) = metadata.incidence_angle_far {
            map.insert("incidence_angle_far".to_string(), far.to_string());
        } else {
            log::warn!("Incidence angle (far) missing from cached metadata");
        }

        // Mid-swath incidence angle - critical for RTC normalization (Small 2011)
        if let Some(mid) = metadata.incidence_angle_mid_swath {
            map.insert("incidence_angle_mid_swath".to_string(), mid.to_string());
            log::info!("RTC reference incidence angle: {:.2}°", mid);
        } else {
            log::warn!(
                "Incidence angle (mid-swath) missing from cached metadata - RTC will use fallback"
            );
        }

        // Add acquisition mode
        map.insert(
            "acquisition_mode".to_string(),
            metadata.acquisition_mode.to_string(),
        );

        // CRITICAL: Add subswaths information for scientific processing
        let subswaths_list: Vec<String> = metadata.sub_swaths.keys().cloned().collect();
        if !subswaths_list.is_empty() {
            map.insert("subswaths".to_string(), subswaths_list.join(","));
            log::info!(
                "Added subswaths to metadata cache: {}",
                subswaths_list.join(",")
            );
        } else {
            log::warn!(
                "No subswaths found in metadata - this will cause scientific processing failure"
            );
        }

        // Add polarizations as comma-separated string
        let polarizations_str: Vec<String> = metadata
            .polarizations
            .iter()
            .map(|p| p.to_string())
            .collect();
        map.insert("polarizations".to_string(), polarizations_str.join(","));

        // CRITICAL FIX: Add orbit state vectors for terrain correction
        // The Python code expects keys like: orbit_state_vector_{i}_time, orbit_state_vector_{i}_x_position, etc.
        // FIX A5: Export times both as Unix seconds AND orbit-relative seconds
        if let Some(orbit_data) = &metadata.orbit_data {
            use slc_reader::time_utils::seconds_since_epoch;
            let orbit_ref = orbit_data.reference_time;

            // Export orbit reference epoch
            let orbit_ref_unix =
                orbit_ref.timestamp() as f64 + (orbit_ref.timestamp_subsec_nanos() as f64) * 1e-9;
            map.insert(
                "orbit_ref_epoch_utc_unix_s".to_string(),
                orbit_ref_unix.to_string(),
            );

            for (i, state_vector) in orbit_data.state_vectors.iter().enumerate() {
                // FIX A5: Export time as BOTH Unix seconds and orbit-relative seconds
                let sv_unix_s = state_vector.time.timestamp() as f64
                    + (state_vector.time.timestamp_subsec_nanos() as f64) * 1e-9;
                let sv_rel_orbit_s = seconds_since_epoch(orbit_ref, state_vector.time);

                // Legacy key (Unix seconds) for backward compatibility
                map.insert(
                    format!("orbit_state_vector_{}_time", i),
                    sv_unix_s.to_string(),
                );
                // New explicit keys for clarity
                map.insert(
                    format!("orbit_state_vector_{}_time_utc_unix_s", i),
                    sv_unix_s.to_string(),
                );
                map.insert(
                    format!("orbit_state_vector_{}_time_rel_orbit_s", i),
                    sv_rel_orbit_s.to_string(),
                );

                // Position components (in meters)
                map.insert(
                    format!("orbit_state_vector_{}_x_position", i),
                    state_vector.position[0].to_string(),
                );
                map.insert(
                    format!("orbit_state_vector_{}_y_position", i),
                    state_vector.position[1].to_string(),
                );
                map.insert(
                    format!("orbit_state_vector_{}_z_position", i),
                    state_vector.position[2].to_string(),
                );

                // Velocity components (in m/s)
                map.insert(
                    format!("orbit_state_vector_{}_x_velocity", i),
                    state_vector.velocity[0].to_string(),
                );
                map.insert(
                    format!("orbit_state_vector_{}_y_velocity", i),
                    state_vector.velocity[1].to_string(),
                );
                map.insert(
                    format!("orbit_state_vector_{}_z_velocity", i),
                    state_vector.velocity[2].to_string(),
                );
            }

            log::info!(
                "Added {} orbit state vectors to metadata cache (with both UTC and orbit-relative times)",
                orbit_data.state_vectors.len()
            );
            
            // SCIENTIFIC VALIDATION (Jan 2026): Ensure orbit vectors span product duration with margin
            if let (Some(first_osv), Some(last_osv)) = (
                orbit_data.state_vectors.first(),
                orbit_data.state_vectors.last()
            ) {
                use slc_reader::time_utils::seconds_since_epoch;
                let orbit_span_s = seconds_since_epoch(first_osv.time, last_osv.time);
                let product_start_s = metadata.start_time.timestamp() as f64 + (metadata.start_time.timestamp_subsec_nanos() as f64) * 1e-9;
                let product_stop_s = metadata.stop_time.timestamp() as f64 + (metadata.stop_time.timestamp_subsec_nanos() as f64) * 1e-9;
                let product_duration_s = product_stop_s - product_start_s;
                let margin_ratio = orbit_span_s / product_duration_s;
                
                log::info!(
                    "🛰️  Orbit temporal coverage: {:.1}s (product duration: {:.1}s, margin ratio: {:.2}x)",
                    orbit_span_s, product_duration_s, margin_ratio
                );
                
                if margin_ratio < 1.0 {
                    log::error!(
                        "❌ CRITICAL ERROR: Orbit vectors ({:.1}s) do not span product duration ({:.1}s). \
                         Interpolation will fail at product boundaries.",
                        orbit_span_s, product_duration_s
                    );
                } else if margin_ratio < 1.1 {
                    log::warn!(
                        "⚠️  Low orbit margin: {:.1}% beyond product duration. Recommend ≥10% for safe interpolation.",
                        (margin_ratio - 1.0) * 100.0
                    );
                } else {
                    log::info!("✅ Orbit vectors provide adequate temporal coverage ({:.1}% margin)", (margin_ratio - 1.0) * 100.0);
                }
            }
        } else {
            log::warn!("No orbit data available in comprehensive metadata for terrain correction");
        }

        Ok(map)
    }

    /// Ensure cache is initialized before accessing cached data
    fn ensure_cache_initialized(&self) -> SarResult<()> {
        if !self.cache_initialized {
            return Err(SarError::Processing(
                "Cache not initialized. Call initialize_all_caches() (construct with SlcReader::new())"
                    .to_string(),
            ));
        }
        Ok(())
    }

    /// Get cache performance statistics
    pub fn get_cache_stats(&self) -> CacheStats {
        let total_annotations = self
            .cached_annotations
            .values()
            .map(|v| v.len())
            .sum::<usize>();
        CacheStats {
            annotations_cached: total_annotations,
            xml_files_cached: self.cached_xml_content.len(),
            calibration_files_cached: self.cached_calibration.len(),
            total_memory_usage_bytes: self.estimate_cache_memory_usage(),
            cache_initialized: self.cache_initialized,
        }
    }

    /// Estimate memory usage of all caches
    fn estimate_cache_memory_usage(&self) -> usize {
        let xml_size: usize = self
            .cached_xml_content
            .values()
            .map(|content| content.len())
            .sum();

        let total_annotations = self
            .cached_annotations
            .values()
            .map(|v| v.len())
            .sum::<usize>();
        let annotation_size =
            total_annotations * std::mem::size_of::<crate::io::annotation::ProductRoot>();

        let calibration_size =
            self.cached_calibration.len() * std::mem::size_of::<CalibrationCoefficients>();

        xml_size + annotation_size + calibration_size + 1024 // Approximate overhead
    }

    /// Helper method to parse annotation XML (used by caching system)
    fn parse_annotation_root_xml(
        &self,
        xml_content: &str,
    ) -> SarResult<crate::io::annotation::ProductRoot> {
        crate::io::annotation::parse_annotation_xml(xml_content)
            .map_err(|e| SarError::Processing(format!("Failed to parse annotation XML: {}", e)))
    }
}

// ============================================================================
// TESTS (orbit parsing)
// ============================================================================
#[cfg(test)]
mod tests_orbit_parse {
    use super::SlcReader;

    #[test]
    fn test_parse_standard_slc_id() {
        let id = "S1A_IW_SLC__1SDV_20200103T170815_20200103T170842_030639_0382D5_DADE";
        let (abs, rel, dir) =
            SlcReader::parse_orbit_from_product_id(id).expect("Should parse standard SLC ID");
        assert_eq!(abs, 30639, "Absolute orbit mismatch");
        assert_eq!(rel, 117, "Relative orbit formula incorrect (expected 117)");
        assert!(dir == "ASCENDING" || dir == "DESCENDING");
    }

    #[test]
    fn test_parse_variant_without_double_underscore() {
        let id = "S1A_IW_GRDH_1SDV_20200103T170815_20200103T170842_030639_0382D5";
        let (abs, rel, _) =
            SlcReader::parse_orbit_from_product_id(id).expect("Should parse variant ID");
        assert_eq!(abs, 30639);
        assert_eq!(rel, 117);
    }

    #[test]
    fn test_parse_failure() {
        let id = "S1A_IW_INVALID_PRODUCT";
        assert!(SlcReader::parse_orbit_from_product_id(id).is_none());
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::types::BoundingBox;
    use std::fs;
    use std::io::Write;
    use std::path::PathBuf;
    use tempfile::tempdir;

    #[test]
    fn test_slc_reader_creation() {
        // This test would need actual test data
        let dummy_path = PathBuf::from("nonexistent.zip");
        let result = SlcReader::new(&dummy_path);
        assert!(result.is_err());
    }

    // --- Small helpers -------------------------------------------------------

    fn write(path: &Path, content: &str) {
        fs::create_dir_all(path.parent().unwrap()).unwrap();
        let mut f = fs::File::create(path).unwrap();
        f.write_all(content.as_bytes()).unwrap();
    }

    /// Build a minimal yet representative Sentinel-1 annotation XML.
    /// - `swath_id`: "IW1" | "IW2" | "IW3"
    /// - `include_downlink_prf`: when false, the PRF must be taken from imageInformation.azimuthFrequency
    /// - geogrid corner box controls (two points) so we can test bbox union
    fn make_annotation_xml(
        swath_id: &str,
        include_downlink_prf: bool,
        prf_hz: f64,
        az_freq_hz: f64,
        srt_s: f64,
        range_ps_m: f64,
        az_ps_m: f64,
        line0_lat: f64,
        line0_lon: f64,
        line1_lat: f64,
        line1_lon: f64,
    ) -> String {
        let downlink = if include_downlink_prf {
            format!(
                r#"
      <downlinkInformationList count="1">
        <downlinkInformation><prf>{}</prf></downlinkInformation>
      </downlinkInformationList>
            "#,
                prf_hz
            )
        } else {
            String::new()
        };

        format!(
            r#"
<product>
  <adsHeader>
    <missionId>S1A</missionId>
    <productType>SLC</productType>
    <startTime>2020-12-28T21:59:42.123456Z</startTime>
    <stopTime>2020-12-28T21:59:52Z</stopTime>
    <absoluteOrbitNumber>30639</absoluteOrbitNumber>
    <missionDataTakeId>123456</missionDataTakeId>
    <imageNumber>1</imageNumber>
  </adsHeader>

  <generalAnnotation>
    <productInformation>
      <platformHeading>12.34</platformHeading>
      <rangeSamplingRate>4.991e+07</rangeSamplingRate>
      <radarFrequency>5.405e+09</radarFrequency>
      <azimuthSteeringRate>0.001</azimuthSteeringRate>
      <rangePixelSpacing>{range_ps_m}</rangePixelSpacing>
      <azimuthPixelSpacing>{az_ps_m}</azimuthPixelSpacing>
    </productInformation>
    {downlink}
    <dcEstimateList count="1">
      <dcEstimate>
        <t0>0.0</t0>
        <dataDcPolynomial>100 2 -0.5</dataDcPolynomial>
      </dcEstimate>
    </dcEstimateList>
        <azimuthFmRateList count="2">
            <azimuthFmRate>
                <azimuthTime>2020-12-28T21:59:42Z</azimuthTime>
                <t0>0.0</t0>
                <azimuthFmRatePolynomial>-3000.0 12.0 -0.002</azimuthFmRatePolynomial>
            </azimuthFmRate>
            <azimuthFmRate>
                <azimuthTime>2020-12-28T21:59:46Z</azimuthTime>
                <t0>0.0</t0>
                <azimuthFmRatePolynomial>-3001.0 11.5 -0.0018</azimuthFmRatePolynomial>
            </azimuthFmRate>
        </azimuthFmRateList>
  </generalAnnotation>

  <imageAnnotation>
    <imageInformation>
      <slantRangeTime>{srt_s}</slantRangeTime>
      <rangePixelSpacing>{range_ps_m}</rangePixelSpacing>
      <azimuthPixelSpacing>{az_ps_m}</azimuthPixelSpacing>
      <numberOfSamples>100</numberOfSamples>
      <numberOfLines>200</numberOfLines>
      <azimuthTimeInterval>0.00028</azimuthTimeInterval>
      <azimuthFrequency>{az_freq_hz}</azimuthFrequency>
      <productFirstLineUtcTime>2020-12-28T21:59:42.123456Z</productFirstLineUtcTime>
      <productLastLineUtcTime>2020-12-28T21:59:52Z</productLastLineUtcTime>
    </imageInformation>
    <processingInformation>
      <swathProcParamsList>
        <swathProcParams><swath>{swath_id}</swath></swathProcParams>
      </swathProcParamsList>
    </processingInformation>
  </imageAnnotation>

  <swathTiming>
    <linesPerBurst>1507</linesPerBurst>
    <samplesPerBurst>21632</samplesPerBurst>
    <burstList count="1">
      <burst>
        <azimuthTime>2020-12-28T21:59:42Z</azimuthTime>
        <azimuthAnxTime>1234.5</azimuthAnxTime>
        <sensingTime>2020-12-28T21:59:42Z</sensingTime>
        <byteOffset>5000000000</byteOffset>
        <firstValidSample>10 10 10</firstValidSample>
        <lastValidSample>90 90 90</lastValidSample>
      </burst>
    </burstList>
  </swathTiming>

  <geolocationGrid>
    <geolocationGridPointList count="2">
      <geolocationGridPoint>
        <slantRangeTime>{srt_s}</slantRangeTime>
        <line>0</line><pixel>0</pixel>
        <latitude>{line0_lat}</latitude><longitude>{line0_lon}</longitude><height>0</height>
        <incidenceAngle>30</incidenceAngle><elevationAngle>5</elevationAngle>
      </geolocationGridPoint>
      <geolocationGridPoint>
        <slantRangeTime>{srt_s}</slantRangeTime>
        <line>1</line><pixel>1</pixel>
        <latitude>{line1_lat}</latitude><longitude>{line1_lon}</longitude><height>0</height>
        <incidenceAngle>40</incidenceAngle><elevationAngle>6</elevationAngle>
      </geolocationGridPoint>
    </geolocationGridPointList>
  </geolocationGrid>

  <antennaPattern>
    <antennaPattern>
      <elevationAngle>1 2 3</elevationAngle>
      <incidenceAngle>30 35 40</incidenceAngle>
      <slantRangeTime>{srt_s} {srt_s} {srt_s}</slantRangeTime>
      <elevationPattern>1 2 3</elevationPattern>
      <terrainHeight>0 0 0</terrainHeight>
      <roll>0</roll>
    </antennaPattern>
  </antennaPattern>

  <!-- SAFE-style orbit blocks -->
  <orbit>
    <time>2020-12-28T21:59:42.123456Z</time>
    <position><x>1</x><y>2</y><z>3</z></position>
    <velocity><x>7000</x><y>0</y><z>0</z></velocity>
  </orbit>
  <orbit>
    <time>2020-12-28T21:59:43Z</time>
    <position><x>2</x><y>3</y><z>4</z></position>
    <velocity><x>7000</x><y>10</y><z>0</z></velocity>
  </orbit>
</product>
"#
        )
    }

    /// Create a temp SAFE directory with three IW annotations for VV.
    /// Returns the temp dir (to keep it alive), product dir path and the expected merged bbox.
    fn build_fake_safe() -> (tempfile::TempDir, PathBuf, BoundingBox) {
        let td = tempdir().unwrap();
        // Name must include .SAFE so SlcReader detects ProductFormat::Safe
        let product_dir = td
            .path()
            .join("S1A_IW_SLC__1SDV_20201228T215942_20201228T215952_030639_0382D5_DADE.SAFE");
        fs::create_dir_all(&product_dir).unwrap();

        let ann_dir = product_dir.join("annotation");
        fs::create_dir_all(&ann_dir).unwrap();
        let calib_dir = ann_dir.join("calibration");
        fs::create_dir_all(&calib_dir).unwrap();

        // IW1: has explicit <prf> (1710 Hz)
        let iw1_xml = make_annotation_xml(
            "IW1", true, 1710.0, 1710.0, 0.0038, 2.329560, 13.943035, 10.0, 20.0, 12.0, 22.0,
        );
        write(&ann_dir.join("s1a-iw1-slc-vv-001.xml"), &iw1_xml);

        // IW2: omit downlink PRF -> must fallback to azimuthFrequency (also 1710 Hz)
        let iw2_xml = make_annotation_xml(
            "IW2", false, 0.0, 1710.0, 0.0039, 2.329560, 13.943035, 11.0, 21.0, 13.0, 23.0,
        );
        write(&ann_dir.join("s1a-iw2-slc-vv-001.xml"), &iw2_xml);

        // IW3: explicit <prf>
        let iw3_xml = make_annotation_xml(
            "IW3", true, 1710.0, 1710.0, 0.0040, 2.329560, 13.943035, 9.0, 19.0, 11.0, 21.0,
        );
        write(&ann_dir.join("s1a-iw3-slc-vv-001.xml"), &iw3_xml);

        // Strict cache initialization requires at least one thermal-noise XML per polarization.
        // A tiny placeholder file is sufficient because tests never parse it.
        write(
            &calib_dir.join("s1a-calibration-noise-vv-001.xml"),
            "<noise xmlns=\"sardine:test\"></noise>",
        );

        // Expected merged bbox across the three swaths
        let merged_bbox = BoundingBox {
            min_lat: 9.0,
            max_lat: 13.0,
            min_lon: 19.0,
            max_lon: 23.0,
        };

        (td, product_dir, merged_bbox)
    }

    fn approx(a: f64, b: f64, eps: f64) {
        assert!((a - b).abs() <= eps, "approx failed: |{a} - {b}| > {eps}");
    }

    // --- Tests ---------------------------------------------------------------

    #[test]
    fn metadata_discovery_and_parsing_end_to_end() {
        let (_temp_dir, product_dir, merged_bbox) = build_fake_safe();

        let mut reader = SlcReader::new(&product_dir).expect("create reader");
        // Find all annotations and ensure 3 IW files for VV
        let all = reader
            .find_all_annotation_files()
            .expect("list annotations");
        let vv_files = all.get(&Polarization::VV).expect("VV present");
        assert_eq!(vv_files.len(), 3);

        // Read + merge metadata for VV across IW1/2/3
        let meta = reader
            .read_annotation(Polarization::VV)
            .expect("read annotation");

        // Product id inferred from SAFE dir name (no .SAFE)
        assert!(meta
            .product_id
            .starts_with("S1A_IW_SLC__1SDV_20201228T215942"));

        // Core radar params pulled from XML (not hardcoded)
        approx(meta.radar_frequency.unwrap(), 5.405e9, 1e3);
        approx(
            meta.wavelength.unwrap(),
            (SPEED_OF_LIGHT_M_S as f64) / 5.405e9,
            5e-6,
        );
        approx(meta.prf.unwrap(), 1710.0, 1e-6);

        // Pixel spacing from imageInformation
        approx(meta.pixel_spacing.0, 2.329560, 1e-9);
        approx(meta.pixel_spacing.1, 13.943035, 1e-6);

        // Subswaths merged
        let names = SlcReader::subswath_names(&meta);
        assert!(names.contains(&"IW1".to_string()));
        assert!(names.contains(&"IW2".to_string()));
        assert!(names.contains(&"IW3".to_string()));

        // Bounding box is the union of the 3 files
        approx(meta.bounding_box.min_lat, merged_bbox.min_lat, 1e-9);
        approx(meta.bounding_box.max_lat, merged_bbox.max_lat, 1e-9);
        approx(meta.bounding_box.min_lon, merged_bbox.min_lon, 1e-9);
        approx(meta.bounding_box.max_lon, merged_bbox.max_lon, 1e-9);

        // Orbit data enriched from <orbit> blocks
        let od = meta.orbit_data.as_ref().expect("orbit data present");
        assert!(od.state_vectors.len() >= 2);
        let vmag = {
            let v = od.state_vectors[0].velocity;
            (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt()
        };
        assert!(vmag > 6900.0 && vmag < 7050.0, "vmag={vmag}");
    }

    #[test]
    fn doppler_centroid_polynomial_via_reader_api() {
        let (_temp_dir, product_dir, _) = build_fake_safe();
        let mut reader = SlcReader::new(&product_dir).expect("create reader");

        // compute doppler at azimuth line 100: t = 100 * 0.00028 = 0.028
        // f(t) = 100 + 2*t - 0.5*t^2 ≈ 100.055608
        let doppler = reader
            .calculate_doppler_centroid(Polarization::VV, 100, 0, None)
            .expect("doppler");
        approx(doppler, 100.055608, 1e-6);
    }

    #[test]
    fn all_iw_subswaths_helpers_work() {
        let (_temp_dir, product_dir, _) = build_fake_safe();
        let mut reader = SlcReader::new(&product_dir).expect("create reader");

        let subs_by_pol = reader.get_all_iw_subswaths().expect("all subs");
        let vv_map = subs_by_pol.get(&Polarization::VV).expect("vv subs");
        assert!(vv_map.contains_key("IW1"));
        assert!(vv_map.contains_key("IW2"));
        assert!(vv_map.contains_key("IW3"));

        let vv_only = reader
            .extract_all_iw_subswaths(Polarization::VV)
            .expect("vv subs only");
        assert_eq!(vv_only.len(), 3);
    }

    #[test]
    fn merged_bbox_utility_matches_union() {
        let (_temp_dir, product_dir, merged_bbox) = build_fake_safe();
        let mut reader = SlcReader::new(&product_dir).expect("create reader");

        let merged = reader
            .extract_merged_bounding_box_all_subswaths(Polarization::VV)
            .expect("merged bbox");

        approx(merged.min_lat, merged_bbox.min_lat, 1e-9);
        approx(merged.max_lat, merged_bbox.max_lat, 1e-9);
        approx(merged.min_lon, merged_bbox.min_lon, 1e-9);
        approx(merged.max_lon, merged_bbox.max_lon, 1e-9);
    }

    #[test]
    fn cache_warmup_and_getters() {
        let (_temp_dir, product_dir, merged_bbox) = build_fake_safe();
        let reader = SlcReader::new_with_full_cache(&product_dir).expect("new_with_full_cache");

        // Cached metadata present
        let meta = reader.get_cached_metadata().expect("cached meta");
        approx(meta.bounding_box.min_lat, merged_bbox.min_lat, 1e-9);
        assert!(
            reader
                .get_all_cached_annotations(Polarization::VV)
                .unwrap()
                .len()
                == 3
        );
        assert!(reader
            .get_available_polarizations()
            .unwrap()
            .contains(&Polarization::VV));

        // Stats are sane
        let stats = reader.get_cache_stats();
        assert!(stats.cache_initialized);
        assert!(stats.annotations_cached >= 3);
        assert!(stats.xml_files_cached >= 3);
    }

    #[test]
    fn error_on_nonexistent_product_path() {
        let dummy = PathBuf::from("really_does_not_exist.zip");
        let res = SlcReader::new(&dummy);
        assert!(res.is_err());
    }

    #[test]
    fn meters_per_degree_wgs84_sanity() {
        // At equator (lat=0):
        let a = crate::constants::geodetic::WGS84_SEMI_MAJOR_AXIS_M;
        let e2 = crate::constants::geodetic::WGS84_ECCENTRICITY_SQUARED;
        let lat_rad = 0.0f64;
        let sin_lat = lat_rad.sin();
        let denom = (1.0 - e2 * sin_lat * sin_lat).sqrt();
        let n = a / denom;
        let m = a * (1.0 - e2) / (1.0 - e2 * sin_lat * sin_lat).powf(1.5);
        let mpd_lon = n * lat_rad.cos() * std::f64::consts::PI / 180.0;
        let mpd_lat = m * std::f64::consts::PI / 180.0;

        // Known reference magnitudes near equator
        approx(mpd_lon, 111_319.5, 300.0); // ~111.32 km (allow generous tolerance)
        approx(mpd_lat, 110_574.0, 300.0); // ~110.57 km
    }
}
