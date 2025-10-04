use crate::constants::physical::SPEED_OF_LIGHT_M_S;
use crate::core::calibrate::CalibrationCoefficients;
use crate::core::deburst::DeburstConfig;
use crate::io::orbit::OrbitManager;
use crate::types::{
    AcquisitionMode, BoundingBox, BurstOrbitData, CoordinateSystem, GeoTransform, OrbitData,
    OrbitStatus, Polarization, SarError, SarImage, SarMetadata, SarResult, SubSwath,
};
use chrono::{DateTime, Utc};
use ndarray::Array2;
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
    /// Create a new SLC reader for a Sentinel-1 product (ZIP or SAFE)
    pub fn new<P: AsRef<Path>>(product_path: P) -> SarResult<Self> {
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

        lf.contains("annotation/calibration/")
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

        lf.contains("annotation/calibration/")
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
                        SarError::Io(std::io::Error::new(std::io::ErrorKind::Other, format!(
                            "Failed to access file {}: {}",
                            i, e
                        )))
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
                    SarError::Io(std::io::Error::new(std::io::ErrorKind::Other, format!(
                        "Failed to read {}: {}",
                        file_path, e
                    )))
                })?;

                let mut content = Vec::new();
                file.read_to_end(&mut content)?;
                Ok(content)
            }
            ProductFormat::Safe => {
                let full_path = self.product_path.join(file_path);
                std::fs::read(full_path).map_err(|e| SarError::Io(e))
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
                let annotation_relative = file.split("annotation/").nth(1).unwrap_or("");
                if annotation_relative.contains('/') {
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
                    SarError::Io(std::io::Error::new(std::io::ErrorKind::Other, format!(
                        "Failed to access {}: {}",
                        measurement_file, e
                    )))
                })?;

                // Create a temporary file to extract the TIFF
                let mut temp_file = NamedTempFile::new().map_err(|e| SarError::Io(e))?;

                // Copy data from ZIP to temporary file
                std::io::copy(&mut zip_file, &mut temp_file).map_err(|e| SarError::Io(e))?;

                let extract_time = start_time.elapsed();
                log::debug!("TIFF extraction took: {:?}", extract_time);

                // Read with GDAL from temporary file
                gdal::Dataset::open(temp_file.path()).map_err(|e| {
                    SarError::Io(std::io::Error::new(std::io::ErrorKind::Other, format!(
                        "Failed to open TIFF with GDAL: {}",
                        e
                    )))
                })?
            }
            ProductFormat::Safe => {
                // Read directly from SAFE directory structure
                let full_path = self.product_path.join(measurement_file);
                gdal::Dataset::open(full_path).map_err(|e| {
                    SarError::Io(std::io::Error::new(std::io::ErrorKind::Other, format!(
                        "Failed to open TIFF with GDAL: {}",
                        e
                    )))
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
                SarError::Io(std::io::Error::new(std::io::ErrorKind::Other, format!(
                    "Failed to get band 1: {}",
                    e
                )))
            })?;

            // **CRITICAL GDAL SAFETY CHECK**: Verify band data type
            let band_type = band.band_type();
            log::info!("🔬 GDAL band type detected: {:?}", band_type);

            // Verify this is actually complex data before attempting complex read
            // Note: GDAL CInt16 corresponds to complex 16-bit integers (real/imag pairs)
            log::info!(
                "🔬 GDAL band type detected: {:?} - proceeding with complex read",
                band_type
            );
            // TODO: Add proper type validation once GDAL enum variants are confirmed

            let window = (0, 0);
            let window_size = (width, height);

            // **IMPROVED COMPLEX READ**: Read with explicit type validation
            // Sentinel-1 SLC uses complex 16-bit integers (CInt16) stored as interleaved i16 pairs
            let complex_data = band
                .read_as::<i16>(window, window_size, (width * 2, height), None)
                .map_err(|e| {
                    SarError::Io(std::io::Error::new(std::io::ErrorKind::Other, format!(
                        "Failed to read complex CInt16 data from GDAL band: {}",
                        e
                    )))
                })?;

            // **SAFETY VALIDATION**: Verify data size matches expected complex format
            let expected_samples = (width * height * 2) as usize; // 2 samples per complex pixel
            if complex_data.data.len() != expected_samples {
                return Err(SarError::InvalidFormat(format!(
                    "GDAL complex data size mismatch: expected {} samples ({}x{}x2), got {} samples. This indicates incorrect complex interpretation.",
                    expected_samples, width, height, complex_data.data.len()
                )));
            }

            log::debug!(
                "✅ Read {} complex values as interleaved i16 ({}x{} pixels)",
                complex_data.data.len() / 2,
                width,
                height
            );

            // Convert interleaved i16 data to complex f32 array
            Self::convert_cint16_to_complex_parallel(complex_data, width, height)?
        } else if band_count >= 2 {
            // Separate I and Q bands (less common for Sentinel-1)
            let band1 = dataset.rasterband(1).map_err(|e| {
                SarError::Io(std::io::Error::new(std::io::ErrorKind::Other, format!(
                    "Failed to get band 1: {}",
                    e
                )))
            })?;

            let band2 = dataset.rasterband(2).map_err(|e| {
                SarError::Io(std::io::Error::new(std::io::ErrorKind::Other, format!(
                    "Failed to get band 2: {}",
                    e
                )))
            })?;

            let window = (0, 0);
            let window_size = (width, height);
            let buffer_size = (width, height);

            let i_data = band1
                .read_as::<f32>(window, window_size, buffer_size, None)
                .map_err(|e| {
                    SarError::Io(std::io::Error::new(std::io::ErrorKind::Other, format!(
                        "Failed to read I data: {}",
                        e
                    )))
                })?;

            let q_data = band2
                .read_as::<f32>(window, window_size, buffer_size, None)
                .map_err(|e| {
                    SarError::Io(std::io::Error::new(std::io::ErrorKind::Other, format!(
                        "Failed to read Q data: {}",
                        e
                    )))
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
                    SarError::Io(std::io::Error::new(std::io::ErrorKind::Other, format!(
                        "Failed to access {}: {}",
                        measurement_file, e
                    )))
                })?;

                // Create a temporary file to extract the TIFF
                let mut temp_file = NamedTempFile::new().map_err(|e| SarError::Io(e))?;

                // Copy data from ZIP to temporary file
                std::io::copy(&mut zip_file, &mut temp_file).map_err(|e| SarError::Io(e))?;

                let extract_time = start_time.elapsed();
                log::debug!("TIFF extraction took: {:?}", extract_time);

                // Open with GDAL
                let dataset = gdal::Dataset::open(temp_file.path()).map_err(|e| {
                    SarError::Io(std::io::Error::new(std::io::ErrorKind::Other, format!(
                        "Failed to open TIFF with GDAL: {}",
                        e
                    )))
                })?;

                (dataset, extract_time)
            }
            ProductFormat::Safe => {
                // Open TIFF directly from SAFE directory
                let extract_time = start_time.elapsed(); // No extraction needed for SAFE
                log::debug!("SAFE direct access took: {:?}", extract_time);

                let full_path = self.product_path.join(measurement_file);
                let dataset = gdal::Dataset::open(full_path).map_err(|e| {
                    SarError::Io(std::io::Error::new(std::io::ErrorKind::Other, format!(
                        "Failed to open TIFF with GDAL: {}",
                        e
                    )))
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
                SarError::Io(std::io::Error::new(std::io::ErrorKind::Other, format!(
                    "Failed to get band 1: {}",
                    e
                )))
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

            // **IMPROVED COMPLEX READ**: Sentinel-1 SLC uses complex 16-bit integers (CInt16)
            let complex_data = band
                .read_as::<i16>(window, window_size, (width * 2, height), None)
                .map_err(|e| {
                    SarError::Io(std::io::Error::new(std::io::ErrorKind::Other, format!(
                        "Failed to read complex CInt16 data from GDAL: {}",
                        e
                    )))
                })?;

            // **SAFETY VALIDATION**: Verify data size for complex format
            let expected_samples = (width * height * 2) as usize;
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
                SarError::Io(std::io::Error::new(std::io::ErrorKind::Other, format!(
                    "Failed to get band 1: {}",
                    e
                )))
            })?;

            let band2 = dataset.rasterband(2).map_err(|e| {
                SarError::Io(std::io::Error::new(std::io::ErrorKind::Other, format!(
                    "Failed to get band 2: {}",
                    e
                )))
            })?;

            let window = (0, 0);
            let window_size = (width, height);
            let buffer_size = (width, height);

            let i_data = band1
                .read_as::<f32>(window, window_size, buffer_size, None)
                .map_err(|e| {
                    SarError::Io(std::io::Error::new(std::io::ErrorKind::Other, format!(
                        "Failed to read I data: {}",
                        e
                    )))
                })?;

            let q_data = band2
                .read_as::<f32>(window, window_size, buffer_size, None)
                .map_err(|e| {
                    SarError::Io(std::io::Error::new(std::io::ErrorKind::Other, format!(
                        "Failed to read Q data: {}",
                        e
                    )))
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

        let prf_from_downlink = annotation.get_pulse_repetition_frequency();
        log::info!(
            "PRF extraction from <prf> field (DownlinkInformation): {:?}",
            prf_from_downlink
        );

        let prf = prf_from_downlink
            .or_else(|| {
                log::warn!(
                    "<prf> field not found in DownlinkInformation, trying <azimuthFrequency> from ImageInformation"
                );
                annotation
                    .image_annotation
                    .as_ref()
                    .and_then(|img| img.image_information.as_ref())
                    .and_then(|info| info.azimuth_frequency)
            })
            .ok_or_else(|| {
                SarError::Processing(
                    "Neither PRF (DownlinkInformation) nor azimuthFrequency (ImageInformation) found in annotation XML - required".to_string(),
                )
            })?;

        log::info!("Final PRF value extracted: {} Hz", prf);

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

        let sub_swaths = crate::io::annotation::ProductRoot::extract_subswaths(annotation).map_err(|e| {
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
                Some(OrbitData {
                    state_vectors: orbit_list.clone(),
                    reference_time: start_time,
                })
            } else {
                log::warn!("No orbit state vectors found in annotation");
                None
            }
        } else {
            log::warn!("No orbit data found in annotation");
            None
        };

        log::info!("Creating SarMetadata with PRF = {} Hz", prf);

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
            radar_frequency_extracted: true,
            bounding_box,
            coordinate_system: CoordinateSystem::Radar,
            sub_swaths,
            orbit_data,
            range_looks: 1,
            azimuth_looks: 1,
            pixel_spacing: (range_pixel_spacing, azimuth_pixel_spacing),
        })
    }

    fn merge_metadata(mut metas: Vec<SarMetadata>) -> SarResult<SarMetadata> {
        if metas.is_empty() {
            return Err(SarError::InvalidFormat(
                "No annotation metadata available to merge".to_string(),
            ));
        }

        let mut base = metas.remove(0);
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

            for (k, v) in m.sub_swaths.drain() {
                base.sub_swaths.insert(k, v);
            }

            for pol in m.polarizations {
                if !base.polarizations.contains(&pol) {
                    base.polarizations.push(pol);
                }
            }

            if strict {
                if base.prf != other_prf || base.radar_frequency != other_radar_frequency {
                    return Err(SarError::Processing(
                        "Inconsistent PRF or radar frequency across subswaths".to_string(),
                    ));
                }
            }
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

    /// Get comprehensive orbit data for the product
    /// Checks: 1) SLC embedded data, 2) Local cache, 3) Download from ESA
    pub fn get_orbit_data(&mut self, orbit_cache_dir: Option<&Path>) -> SarResult<OrbitData> {
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
        // SCIENTIFIC REQUIREMENT: Orbit cache directory must be explicitly specified
        let cache_dir = match orbit_cache_dir {
            Some(dir) => dir.to_path_buf(),
            None => {
                return Err(SarError::Processing(
                    "❌ SCIENTIFIC ERROR: Orbit cache directory must be explicitly specified! 
                    No fallback to temp directory permitted for scientific accuracy."
                        .to_string(),
                ));
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

        // SCIENTIFIC REQUIREMENT: Orbit cache directory must be explicitly specified
        let cache_dir = match orbit_cache_dir {
            Some(dir) => dir.to_path_buf(),
            None => {
                return Err(SarError::Processing(
                    "❌ SCIENTIFIC ERROR: Orbit cache directory must be explicitly specified! 
                    No fallback to temp directory permitted for scientific accuracy."
                        .to_string(),
                ));
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

        // SCIENTIFIC REQUIREMENT: Orbit cache directory must be explicitly specified
        let cache_dir = match orbit_cache_dir {
            Some(dir) => dir.to_path_buf(),
            None => {
                return Err(SarError::Processing(
                    "❌ SCIENTIFIC ERROR: Orbit cache directory must be explicitly specified! 
                    No fallback to temp directory permitted for scientific accuracy."
                        .to_string(),
                ));
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
        let orbit_data = self.get_orbit_data(orbit_cache_dir)?;

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
                    SarError::Io(std::io::Error::new(std::io::ErrorKind::Other, format!(
                        "Failed to access {}: {}",
                        measurement_file, e
                    )))
                })?;
                let mut temp_file = NamedTempFile::new().map_err(SarError::Io)?;
                std::io::copy(&mut zip_file, &mut temp_file).map_err(SarError::Io)?;
                let dataset = gdal::Dataset::open(temp_file.path()).map_err(|e| {
                    SarError::Io(std::io::Error::new(std::io::ErrorKind::Other, format!(
                        "Failed to open TIFF with GDAL: {}",
                        e
                    )))
                })?;
                let rs = dataset.raster_size();
                (rs.0, rs.1)
            }
            ProductFormat::Safe => {
                let full_path = self.product_path.join(measurement_file);
                let dataset = gdal::Dataset::open(full_path).map_err(|e| {
                    SarError::Io(std::io::Error::new(std::io::ErrorKind::Other, format!(
                        "Failed to open TIFF with GDAL: {}",
                        e
                    )))
                })?;
                let rs = dataset.raster_size();
                (rs.0, rs.1)
            }
        };

        log::info!("Image dimensions: {} x {} pixels", width, height);

        // Calculate azimuth time interval
        // For Sentinel-1 IW mode, typical azimuth time interval is ~2.8e-4 seconds
        let acquisition_duration =
            (annotation.stop_time - annotation.start_time).num_milliseconds() as f64 / 1000.0;
        let azimuth_time_interval = acquisition_duration / height as f64;

        log::info!(
            "Azimuth timing: duration={:.3}s, interval={:.6}s, lines={}",
            acquisition_duration,
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
    /// This uses satellite velocity and geometry to estimate Doppler shift
    pub fn calculate_doppler_centroid(
        &mut self,
        pol: Polarization,
        azimuth_line: usize,
        _range_sample: usize,
        _orbit_cache_dir: Option<&Path>,
    ) -> SarResult<f64> {
        // Use real Doppler centroid model from annotation (dcEstimate.dataDcPolynomial)
        // Map azimuth_line to time since product start using azimuth_time_interval derived from PRF or annotation
        let annotation_root = self.get_annotation_for_polarization(pol)?;
        let _sar_meta = self.read_annotation(pol)?;

        // Determine azimuth time interval (s/line)
        let az_time_interval = if let Some(image) = &annotation_root.image_annotation {
            if let Some(info) = &image.image_information {
                if let Some(dt) = info.azimuth_time_interval {
                    dt
                } else {
                    // Compute from PRF if available: time per line ~ 1/PRF
                    if let Some(prf) = info.azimuth_frequency {
                        1.0 / prf
                    } else {
                        return Err(SarError::Metadata(
                            "Missing azimuthTimeInterval/PRF in annotation".to_string(),
                        ));
                    }
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

        // Time since start in seconds
        let az_time_since_start = azimuth_line as f64 * az_time_interval;

        let doppler_hz = annotation_root.evaluate_doppler_centroid(az_time_since_start)?;
        log::debug!(
            "Doppler (dcPolynomial) at line {} (t={:.6}s): {:.2} Hz",
            azimuth_line,
            az_time_since_start,
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
            if file.contains("annotation/")
                && file.ends_with(".xml")
                && !file.contains("calibration")
            {
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
    pub fn get_all_iw_subswaths(
        &mut self,
    ) -> SarResult<HashMap<Polarization, HashMap<String, SubSwath>>> {
        log::info!("Starting get_all_iw_subswaths");
        let all_annotations = self.find_all_iw_annotation_files()?;
        log::info!(
            "Found annotation files for {} polarizations",
            all_annotations.len()
        );
        let mut all_subswaths = HashMap::new();

        for pol in all_annotations.keys().cloned() {
            log::info!("Processing polarization: {:?}", pol);
            match self.extract_all_iw_subswaths(pol) {
                Ok(subswaths) => {
                    log::info!(
                        "Extracted {} subswaths for polarization {:?}",
                        subswaths.len(),
                        pol
                    );
                    all_subswaths.insert(pol, subswaths);
                }
                Err(e) => {
                    log::error!("Failed to extract sub-swaths for {:?}: {}", pol, e);
                }
            }
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
            if file.contains("annotation/") && file.ends_with(".xml") && file.contains("-iw") {
                return Ok(true);
            }
        }

        Ok(false)
    }

    /// Deburst SLC data for a specific polarization
    pub fn deburst_slc(&mut self, pol: Polarization) -> SarResult<SarImage> {
        use crate::core::deburst::{DeburstConfig, DeburstProcessor};

        log::info!("Starting deburst processing for polarization {:?}", pol);

        // Read SLC data
        let slc_data = self.read_slc_data(pol)?;
        log::debug!("Read SLC data with dimensions: {:?}", slc_data.dim());

        // Get annotation for burst information
        let annotation_files = self.find_all_annotation_files()?;
        let files = annotation_files.get(&pol).ok_or_else(|| {
            SarError::Processing(format!("No annotation files found for {:?}", pol))
        })?;
        let first = files
            .first()
            .ok_or_else(|| SarError::Processing(format!("Empty annotation list for {:?}", pol)))?;

        // Read annotation content
        let annotation_content = self.read_file_as_string(first)?;

        // Extract burst information
        let (total_lines, total_samples) = slc_data.dim();
        let burst_info = DeburstProcessor::extract_burst_info_from_annotation(
            &annotation_content,
            total_lines,
            total_samples,
        )?;

        log::info!(
            "Extracted {} bursts for deburst processing",
            burst_info.len()
        );

        // Create deburst processor
        // Scientific requirement: derive satellite velocity from real orbit data (no hardcoded values)
        let cache_dir = std::env::var("SARDINE_ORBIT_CACHE")
            .map(std::path::PathBuf::from)
            .map_err(|_| {
                SarError::Processing(
                    "SARDINE_ORBIT_CACHE not set. Required for orbit-based velocity in deburst."
                        .to_string(),
                )
            })?;
        let orbit_data = self.get_orbit_data(Some(cache_dir.as_path()))?;
        if orbit_data.state_vectors.is_empty() {
            return Err(SarError::Processing(
                "No orbit state vectors available for velocity calculation".to_string(),
            ));
        }
        let v = &orbit_data.state_vectors[0].velocity;
        let satellite_velocity = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
        if satellite_velocity < 7000.0 || satellite_velocity > 8000.0 {
            return Err(SarError::Processing(format!(
                "Invalid satellite velocity: {:.1} m/s (expected 7000-8000 m/s)",
                satellite_velocity
            )));
        }
        let processor = DeburstProcessor::new(burst_info, satellite_velocity);

        // Configure deburst processing
        let config = DeburstConfig {
            blend_overlap: true,
            blend_lines: 10,
            remove_invalid_data: true,
            seamless_stitching: true,
            apply_deramp: true,
            preserve_phase: true,
            antenna_pattern_correction: false,
            use_range_dependent_deramp: false, // Use fast time-only deramp

            // NEW: Scientific enhancements (use conservative defaults)
            use_annotation_timing: true,
            enable_bilinear_interp: false,
            enable_hit_count_mask: true,
            power_preservation_check: true,
        };

        // Perform deburst
        let deburst_data = processor.deburst(&slc_data, &config)?;

        log::info!(
            "Deburst completed. Output dimensions: {:?}",
            deburst_data.dim()
        );
        Ok(deburst_data)
    }

    /// Deburst SLC data with custom configuration for TOPSAR processing
    pub fn deburst_slc_with_config(
        &mut self,
        pol: Polarization,
        config: DeburstConfig,
    ) -> SarResult<SarImage> {
        use crate::core::deburst::DeburstProcessor;

        log::info!(
            "Starting TOPSAR deburst processing for polarization {:?}",
            pol
        );

        // Read SLC data
        let slc_data = self.read_slc_data(pol)?;
        log::debug!("Read SLC data with dimensions: {:?}", slc_data.dim());

        // Get annotation for burst information
        let annotation_files = self.find_all_annotation_files()?;
        let files = annotation_files.get(&pol).ok_or_else(|| {
            SarError::Processing(format!("No annotation files found for {:?}", pol))
        })?;
        let first = files
            .first()
            .ok_or_else(|| SarError::Processing(format!("Empty annotation list for {:?}", pol)))?;

        // Read annotation content
        let annotation_content = self.read_file_as_string(first)?;

        // Extract burst information with TOPSAR parameters
        let (total_lines, total_samples) = slc_data.dim();
        let burst_info = DeburstProcessor::extract_burst_info_from_annotation(
            &annotation_content,
            total_lines,
            total_samples,
        )?;

        log::info!(
            "Extracted {} bursts with TOPSAR parameters for deburst processing",
            burst_info.len()
        );

        // Create deburst processor
        // Scientific requirement: derive satellite velocity from real orbit data (no hardcoded values)
        let cache_dir = std::env::var("SARDINE_ORBIT_CACHE")
            .map(std::path::PathBuf::from)
            .map_err(|_| {
                SarError::Processing(
                    "SARDINE_ORBIT_CACHE not set. Required for orbit-based velocity in deburst."
                        .to_string(),
                )
            })?;
        let orbit_data = self.get_orbit_data(Some(cache_dir.as_path()))?;
        if orbit_data.state_vectors.is_empty() {
            return Err(SarError::Processing(
                "No orbit state vectors available for velocity calculation".to_string(),
            ));
        }
        let v = &orbit_data.state_vectors[0].velocity;
        let satellite_velocity = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
        if satellite_velocity < 7000.0 || satellite_velocity > 8000.0 {
            return Err(SarError::Processing(format!(
                "Invalid satellite velocity: {:.1} m/s (expected 7000-8000 m/s)",
                satellite_velocity
            )));
        }
        let processor = DeburstProcessor::new(burst_info, satellite_velocity);

        // Perform deburst with custom configuration
        let deburst_data = processor.deburst(&slc_data, &config)?;

        log::info!(
            "TOPSAR deburst completed. Output dimensions: {:?}",
            deburst_data.dim()
        );
        Ok(deburst_data)
    }

    /// Find ALL calibration files for each polarization in the ZIP
    pub fn find_all_calibration_files(&mut self) -> SarResult<HashMap<Polarization, Vec<String>>> {
        log::debug!("Finding all calibration files (multi-subswath support)");

        let files = self.list_files()?;
        let mut calibration_files: HashMap<Polarization, Vec<String>> = HashMap::new();

        for file in files {
            // Only include actual calibration files, not noise files
            if Self::is_calibration_xml(&file) {
                if let Some(pol) = Self::extract_polarization(&file) {
                    // Collect ALL calibration files for each polarization (IW1, IW2, IW3)
                    calibration_files
                        .entry(pol)
                        .or_insert_with(Vec::new)
                        .push(file);
                }
            }
        }

        log::info!(
            "Found calibration files: {:?}",
            calibration_files
                .iter()
                .map(|(pol, files)| (pol, files.len()))
                .collect::<Vec<_>>()
        );
        Ok(calibration_files)
    }

    /// Find calibration files for each polarization in the ZIP (legacy - single file per pol)
    pub fn find_calibration_files(&mut self) -> SarResult<HashMap<Polarization, String>> {
        log::debug!("Finding calibration files (legacy single-file mode)");

        let files = self.list_files()?;
        let mut calibration_files = HashMap::new();

        for file in files {
            // Only include actual calibration files, not noise files
            if Self::is_calibration_xml(&file) {
                if let Some(pol) = Self::extract_polarization(&file) {
                    // For now, just take the first calibration file found for each polarization
                    // Future enhancement: handle multiple subswaths per polarization
                    if !calibration_files.contains_key(&pol) {
                        calibration_files.insert(pol, file);
                    }
                }
            }
        }

        log::info!("Found {} calibration files", calibration_files.len());
        Ok(calibration_files)
    }

    /// Read calibration data for a specific polarization (handles multiple subswaths)
    pub fn read_calibration_data(
        &mut self,
        pol: Polarization,
    ) -> SarResult<crate::core::calibrate::CalibrationCoefficients> {
        log::debug!("Reading calibration data for polarization {:?}", pol);

        let all_calibration_files = self.find_all_calibration_files()?;
        let file_paths = all_calibration_files.get(&pol).ok_or_else(|| {
            SarError::Processing(format!(
                "No calibration files found for polarization {:?}",
                pol
            ))
        })?;

        log::info!(
            "Found {} calibration files for {:?}: {:?}",
            file_paths.len(),
            pol,
            file_paths
        );

        // Parse and merge calibration data from all subswaths
        let mut merged_calibration = None;

        for file_path in file_paths {
            let xml_content = self.read_file_as_string(file_path)?;
            let calibration_data =
                crate::core::calibrate::parse_calibration_from_xml(&xml_content)?;

            log::info!(
                "Parsed {} calibration vectors from {} (swath: {})",
                calibration_data.vectors.len(),
                file_path,
                calibration_data.swath
            );

            match merged_calibration {
                None => {
                    // First calibration file - use as base and set swath to "IW" for IW mode
                    let mut base_calibration = calibration_data;
                    base_calibration.swath = "IW".to_string();
                    merged_calibration = Some(base_calibration);
                }
                Some(ref mut merged) => {
                    // Merge vectors from additional subswaths
                    merged.vectors.extend(calibration_data.vectors);
                    log::info!(
                        "Merged calibration vectors - total now: {}",
                        merged.vectors.len()
                    );

                    // Swath field already set to "IW" from first file
                }
            }
        }

        let final_calibration = merged_calibration.ok_or_else(|| {
            SarError::Processing(format!(
                "No valid calibration data found for polarization {:?}",
                pol
            ))
        })?;

        log::info!(
            "Successfully merged calibration data for {:?} - total {} vectors from {} files",
            pol,
            final_calibration.vectors.len(),
            file_paths.len()
        );
        Ok(final_calibration)
    }

    /// Find noise files for all available polarizations
    pub fn find_noise_files(&mut self) -> SarResult<HashMap<Polarization, String>> {
        log::debug!("Finding noise files");

        let files = self.list_files()?;
        let mut noise_files = HashMap::new();

        for file in files {
            // Look for noise files in annotation/calibration/ directory
            if file.contains("annotation/calibration/")
                && file.ends_with(".xml")
                && file.contains("noise-")
            {
                if let Some(pol) = Self::extract_polarization(&file) {
                    // For now, just take the first noise file found for each polarization
                    // Future enhancement: handle multiple subswaths per polarization
                    if !noise_files.contains_key(&pol) {
                        noise_files.insert(pol, file);
                    }
                }
            }
        }

        log::info!("Found {} noise files", noise_files.len());
        Ok(noise_files)
    }

    /// Read thermal noise data for a specific polarization
    pub fn read_noise_data(
        &mut self,
        pol: Polarization,
    ) -> SarResult<crate::core::calibrate::NoiseCoefficients> {
        log::debug!("Reading noise data for polarization {:?}", pol);

        let noise_files = self.find_noise_files()?;
        let file_path = noise_files.get(&pol).ok_or_else(|| {
            SarError::Processing(format!("No noise file found for polarization {:?}", pol))
        })?;

        let xml_content = self.read_file_as_string(file_path)?;

        // Parse thermal noise data from XML
        let noise_data = crate::core::calibrate::parse_noise_from_xml(&xml_content)?;

        log::info!("Successfully parsed noise data for polarization {:?}", pol);
        Ok(noise_data)
    }

    /// Read calibration data with caching
    pub fn read_calibration_data_cached(
        &mut self,
        pol: Polarization,
    ) -> SarResult<CalibrationCoefficients> {
        let cache_key = format!("{}", pol);

        // Check if already cached
        if let Some(cached) = self.calibration_cache.get(&cache_key) {
            log::debug!("Using cached calibration data for {}", pol);
            return Ok(cached.clone());
        }

        // Read and cache calibration data
        log::debug!("Reading calibration data for {} (not cached)", pol);
        let cal_data = self.read_calibration_data(pol)?;
        self.calibration_cache.insert(cache_key, cal_data.clone());

        Ok(cal_data)
    }

    /// Read calibration data for all available polarizations
    pub fn read_all_calibration_data(
        &mut self,
    ) -> SarResult<HashMap<Polarization, crate::core::calibrate::CalibrationCoefficients>> {
        let mut all_cal_data = HashMap::new();

        let available_pols = self
            .find_all_annotation_files()?
            .keys()
            .cloned()
            .collect::<Vec<_>>();
        for pol in available_pols {
            match self.read_calibration_data(pol) {
                Ok(cal_data) => {
                    all_cal_data.insert(pol, cal_data);
                }
                Err(e) => {
                    log::warn!("Failed to read calibration data for {:?}: {}", pol, e);
                }
            }
        }

        Ok(all_cal_data)
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

    /// Complete workflow: calibrate and multilook SLC data
    ///
    /// This is a convenience method that combines calibration and multilooking
    pub fn calibrate_and_multilook(
        &mut self,
        pol: Polarization,
        cal_type: crate::core::calibrate::CalibrationType,
        range_looks: usize,
        azimuth_looks: usize,
    ) -> SarResult<(Array2<f32>, f64, f64)> {
        log::info!("Starting calibrate and multilook workflow for {:?}", pol);

        // First, get deburst data
        let deburst_data = self.deburst_slc(pol)?;
        log::info!(
            "Deburst data: {} x {}",
            deburst_data.nrows(),
            deburst_data.ncols()
        );

        // Get calibration coefficients
        let cal_data = self.read_calibration_data(pol)?;

        // Create calibration processor
        let processor = crate::core::calibrate::CalibrationProcessor::new(cal_data, cal_type);

        // Apply calibration to get intensity data
        let intensity_data = processor.calibrate(&deburst_data)?;
        log::info!(
            "Calibrated data: {} x {}",
            intensity_data.nrows(),
            intensity_data.ncols()
        );

        // Apply multilooking
        let (multilooked_data, new_range_spacing, new_azimuth_spacing) =
            self.multilook_intensity(&intensity_data, pol, range_looks, azimuth_looks)?;

        log::info!(
            "Complete workflow finished: final dimensions {}x{}",
            multilooked_data.nrows(),
            multilooked_data.ncols()
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
        // Simplified geotransform creation
        // In practice, this would use precise geo-location from metadata
        let bbox = &metadata.bounding_box;

        Ok(GeoTransform {
            top_left_x: bbox.min_lon,
            pixel_width: range_spacing / 111320.0, // Approximate meters to degrees
            rotation_x: 0.0,
            top_left_y: bbox.max_lat,
            rotation_y: 0.0,
            pixel_height: -azimuth_spacing / 111320.0, // Negative for north-up
        })
    }

    /// Extract product type from XML or infer from product ID
    fn extract_product_type(xml_content: &str, product_id: &str) -> Option<String> {
        // Try to extract from XML first - use the correct Sentinel-1 tag
        if let Some(product_type) = Self::extract_xml_value(xml_content, "productType") {
            return Some(product_type);
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
                log::warn!("No 6-digit numeric segment found in product ID '{}'", product_id);
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

        // Orbit direction heuristic (placeholder): true determination needs precise state vectors.
        // Keep simple parity-based heuristic but document limitation.
        let orbit_direction = if orbit_number % 2 == 0 {
            "DESCENDING"
        } else {
            "ASCENDING"
        };

        log::debug!(
            "Orbit parse result: abs_orbit={} (idx {}), rel_orbit={}, dir={}",
            orbit_number, idx, relative_orbit, orbit_direction
        );
        Some((orbit_number, relative_orbit, orbit_direction.to_string()))
    }


// (tests moved to end of file after impl block)
    /// Extract XML value using simple string parsing (fallback method)
    fn extract_xml_value(xml_content: &str, tag_name: &str) -> Option<String> {
        let start_tag = format!("<{}>", tag_name);
        let end_tag = format!("</{}>", tag_name);

        if let Some(start_pos) = xml_content.find(&start_tag) {
            let start_content = start_pos + start_tag.len();
            if let Some(end_pos) = xml_content[start_content..].find(&end_tag) {
                let value = &xml_content[start_content..start_content + end_pos];
                return Some(value.trim().to_string());
            }
        }
        None
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

    /// Extract merged bounding box from ALL annotation files for complete scene coverage
    fn extract_merged_bounding_box_all_subswaths(
        &mut self,
        pol: Polarization,
    ) -> SarResult<crate::types::BoundingBox> {
        println!(
            "🌐 DEBUG: Starting merged bounding box extraction for polarization {}",
            pol
        );
        log::info!("🌐 Extracting merged bounding box from ALL subswaths for complete coverage");

        // Get ALL annotation files for this polarization
        println!("📁 DEBUG: Getting all annotation files...");
        let all_annotation_files = match self.find_all_annotation_files() {
            Ok(files) => files,
            Err(e) => {
                println!("❌ DEBUG: Failed to find annotation files: {}", e);
                return Err(e);
            }
        };

        println!(
            "📁 DEBUG: Found annotation files for {} polarizations",
            all_annotation_files.len()
        );

        let pol_files = all_annotation_files.get(&pol).ok_or_else(|| {
            println!(
                "❌ DEBUG: No annotation files found for polarization {}",
                pol
            );
            SarError::Processing(format!(
                "No annotation files found for polarization {}",
                pol
            ))
        })?;

        println!(
            "📁 DEBUG: Found {} annotation files for {}",
            pol_files.len(),
            pol
        );
        log::info!(
            "📁 Found {} annotation files for {}: merging bounding boxes",
            pol_files.len(),
            pol
        );

        let mut merged_min_lat = f64::INFINITY;
        let mut merged_max_lat = f64::NEG_INFINITY;
        let mut merged_min_lon = f64::INFINITY;
        let mut merged_max_lon = f64::NEG_INFINITY;

        let mut processed_count = 0;

        // Process each annotation file and merge bounding boxes
        for annotation_file in pol_files {
            log::debug!("📄 Processing annotation file: {}", annotation_file);

            match self.read_annotation_file_raw(annotation_file) {
                Ok(sar_metadata) => {
                    let bbox = &sar_metadata.bounding_box;

                    // Validate this subswath's bounding box
                    if bbox.min_lat < bbox.max_lat
                        && bbox.min_lon < bbox.max_lon
                        && bbox.min_lat != 0.0
                        && bbox.max_lat != 0.0
                        && bbox.min_lon != 0.0
                        && bbox.max_lon != 0.0
                    {
                        // Merge with overall bounds
                        merged_min_lat = merged_min_lat.min(bbox.min_lat);
                        merged_max_lat = merged_max_lat.max(bbox.max_lat);
                        merged_min_lon = merged_min_lon.min(bbox.min_lon);
                        merged_max_lon = merged_max_lon.max(bbox.max_lon);

                        processed_count += 1;

                        log::info!(
                            "✅ Subswath {}: [{:.6}, {:.6}, {:.6}, {:.6}]",
                            processed_count,
                            bbox.min_lon,
                            bbox.min_lat,
                            bbox.max_lon,
                            bbox.max_lat
                        );
                    } else {
                        log::warn!("⚠️  Invalid bounding box in {}, skipping", annotation_file);
                    }
                }
                Err(e) => {
                    log::warn!(
                        "⚠️  Failed to read annotation file {}: {}",
                        annotation_file,
                        e
                    );
                }
            }
        }

        if processed_count == 0 {
            return Err(SarError::Processing(
                "No valid bounding boxes found in any annotation files".to_string(),
            ));
        }

        let merged_bbox = crate::types::BoundingBox {
            min_lat: merged_min_lat,
            max_lat: merged_max_lat,
            min_lon: merged_min_lon,
            max_lon: merged_max_lon,
        };

        log::info!("🎯 MERGED BOUNDING BOX from {} subswaths:", processed_count);
        log::info!(
            "   📍 [{:.6}, {:.6}, {:.6}, {:.6}]",
            merged_bbox.min_lon,
            merged_bbox.min_lat,
            merged_bbox.max_lon,
            merged_bbox.max_lat
        );
        log::info!(
            "   📏 Coverage: {:.3}deg x {:.3}deg ({:.1}km x {:.1}km)",
            merged_bbox.max_lon - merged_bbox.min_lon,
            merged_bbox.max_lat - merged_bbox.min_lat,
            (merged_bbox.max_lon - merged_bbox.min_lon) * 111.0,
            (merged_bbox.max_lat - merged_bbox.min_lat) * 111.0
        );

        Ok(merged_bbox)
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
        let mut reader = Self::new(product_path)?;

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

        for (pol, files) in &all_annotation_files {
            for annotation_path in files {
                log::debug!("📄 Caching annotation for {:?}: {}", pol, annotation_path);
                let xml_content = self.read_file_as_string(annotation_path)?;
                self.cached_xml_content
                    .insert(annotation_path.clone(), xml_content.clone());
                // Parse and collect all subswaths per polarization (IW1, IW2, IW3)
                let annotation_root = Arc::new(self.parse_annotation_root_xml(&xml_content)?);

                // Ensure we have a vector for this polarization and append the annotation
                self.cached_annotations
                    .entry(*pol)
                    .or_insert_with(Vec::new)
                    .push(annotation_root);
            }
            // Cache calibration data once per polarization
            if let Ok(calibration) = self.read_calibration_data(*pol) {
                self.cached_calibration.insert(*pol, calibration);
                log::debug!("✅ Cached calibration for {:?}", pol);
            } else {
                log::warn!("⚠️  Could not cache calibration for {:?}", pol);
            }
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

        // SCIENTIFIC FIX: Extract subswaths from ALL cached annotations, not just first polarization
        let mut all_subswaths = HashMap::new();
        log::info!("🔍 Extracting subswaths from all cached annotation files");
        for (pol, roots) in &self.cached_annotations {
            for root in roots {
                match crate::io::annotation::ProductRoot::extract_subswaths(root.as_ref()) {
                    Ok(pol_subswaths) => {
                        for (swath_id, swath_data) in pol_subswaths {
                            log::info!("✅ Found subswath {} from {:?}", swath_id, pol);
                            all_subswaths.insert(swath_id, swath_data);
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

        log::info!(
            "✅ Combined subswaths from all annotations: {:?}",
            all_subswaths.keys().collect::<Vec<_>>()
        );

        comprehensive_metadata.sub_swaths = all_subswaths;

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
                "Metadata not cached. Use new_with_full_cache() instead of new()".to_string(),
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

        // CRITICAL: Add PRF for terrain correction
        log::info!("Checking PRF in SarMetadata struct: {:?}", metadata.prf);
        if let Some(prf) = metadata.prf {
            log::info!("Adding PRF to metadata cache: {} Hz", prf);
            map.insert("prf".to_string(), prf.to_string());
        } else {
            log::warn!("PRF is None in cached SarMetadata - extraction failed!");
        }

        // CRITICAL FIX: Add product_start_time_abs for correct terrain correction
        // Convert start_time (DateTime<Utc>) to Unix timestamp in seconds
        let product_start_time_abs = metadata.start_time.timestamp() as f64 
            + (metadata.start_time.timestamp_subsec_nanos() as f64) * 1e-9;
        map.insert("product_start_time_abs".to_string(), product_start_time_abs.to_string());
        log::info!(
            "🔬 SCIENTIFIC FIX: Added product_start_time_abs={:.6}s to metadata for coordinate calculations", 
            product_start_time_abs
        );

        // Add pixel spacing
        map.insert(
            "range_pixel_spacing".to_string(),
            metadata.pixel_spacing.0.to_string(),
        );
        map.insert(
            "azimuth_pixel_spacing".to_string(),
            metadata.pixel_spacing.1.to_string(),
        );

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
        if let Some(orbit_data) = &metadata.orbit_data {
            for (i, state_vector) in orbit_data.state_vectors.iter().enumerate() {
                // Convert time to timestamp (seconds since Unix epoch)
                let timestamp = state_vector.time.timestamp() as f64;
                map.insert(
                    format!("orbit_state_vector_{}_time", i),
                    timestamp.to_string(),
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
                "Added {} orbit state vectors to metadata cache",
                orbit_data.state_vectors.len()
            );
        } else {
            log::warn!("No orbit data available in comprehensive metadata for terrain correction");
        }

        Ok(map)
    }

    /// Ensure cache is initialized before accessing cached data
    fn ensure_cache_initialized(&self) -> SarResult<()> {
        if !self.cache_initialized {
            return Err(SarError::Processing(
                "Cache not initialized. Use SlcReader::new_with_full_cache() instead of new()"
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
        let (abs, rel, dir) = SlcReader::parse_orbit_from_product_id(id)
            .expect("Should parse standard SLC ID");
        assert_eq!(abs, 30639, "Absolute orbit mismatch");
        assert_eq!(rel, 117, "Relative orbit formula incorrect (expected 117)");
        assert!(dir == "ASCENDING" || dir == "DESCENDING");
    }

    #[test]
    fn test_parse_variant_without_double_underscore() {
        let id = "S1A_IW_GRDH_1SDV_20200103T170815_20200103T170842_030639_0382D5";
        let (abs, rel, _) = SlcReader::parse_orbit_from_product_id(id)
            .expect("Should parse variant ID");
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
        let mut reader = SlcReader::new_with_full_cache(&product_dir).expect("new_with_full_cache");

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
}
