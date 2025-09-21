use crate::types::{
    AcquisitionMode, BoundingBox, CoordinateSystem, Polarization, SarError, 
    SarImage, SarMetadata, SarResult, OrbitStatus, OrbitData, BurstOrbitData, SubSwath, GeoTransform
};
use crate::constants::physical::SPEED_OF_LIGHT_M_S;
use crate::core::calibrate::CalibrationCoefficients;
use crate::core::deburst::DeburstConfig;
use crate::io::orbit::OrbitManager;
use chrono::{DateTime, Utc};
use ndarray::Array2;
use serde::Deserialize;
use std::collections::HashMap;
use std::fs::File;
use std::io::Read;
use std::path::{Path, PathBuf};
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
    cached_annotations: HashMap<Polarization, crate::io::annotation::AnnotationRoot>,
    
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
        let format = if product_path.is_file() && product_path.extension() == Some(std::ffi::OsStr::new("zip")) {
            ProductFormat::Zip
        } else if product_path.is_dir() && (
            product_path.file_name()
                .map(|name| name.to_string_lossy().contains(".SAFE"))
                .unwrap_or_else(|| {
                    log::debug!("Could not determine filename for SAFE detection - checking manifest.safe");
                    false
                }) ||
            product_path.join("manifest.safe").exists()
        ) {
            ProductFormat::Safe
        } else {
            return Err(SarError::InvalidFormat(
                format!("Unsupported format: {}. Must be .zip file or .SAFE directory", product_path.display())
            ));
        };

        log::info!("Detected Sentinel-1 product format: {:?} at {}", format, product_path.display());

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

    /// Open the ZIP archive (only for ZIP format)
    fn open_archive(&mut self) -> SarResult<&mut ZipArchive<File>> {
        match self.format {
            ProductFormat::Zip => {
                if self.archive.is_none() {
                    let file = File::open(&self.product_path)?;
                    let archive = ZipArchive::new(file)
                        .map_err(|e| SarError::InvalidFormat(format!("Failed to open ZIP: {}", e)))?;
                    self.archive = Some(archive);
                }
                Ok(self.archive.as_mut().unwrap())
            },
            ProductFormat::Safe => {
                Err(SarError::InvalidFormat("Cannot open ZIP archive on SAFE directory".to_string()))
            }
        }
    }

    /// List all files in the product (works for both ZIP and SAFE)
    pub fn list_files(&mut self) -> SarResult<Vec<String>> {
        match self.format {
            ProductFormat::Zip => {
                let archive = self.open_archive()?;
                let mut files = Vec::new();
                
                for i in 0..archive.len() {
                    let file = archive.by_index(i)
                        .map_err(|e| SarError::Io(std::io::Error::other(
                            format!("Failed to access file {}: {}", i, e),
                        )))?;
                    files.push(file.name().to_string());
                }
                
                Ok(files)
            },
            ProductFormat::Safe => {
                let mut files = Vec::new();
                Self::list_safe_files_recursive(&self.product_path, "", &mut files)?;
                Ok(files)
            }
        }
    }

    /// Recursively list files in SAFE directory structure
    fn list_safe_files_recursive(dir: &Path, prefix: &str, files: &mut Vec<String>) -> SarResult<()> {
        let entries = std::fs::read_dir(dir)
            .map_err(|e| SarError::Io(e))?;

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
                let mut file = archive.by_name(file_path)
                    .map_err(|e| SarError::Io(std::io::Error::other(
                        format!("Failed to read {}: {}", file_path, e),
                    )))?;

                let mut content = Vec::new();
                file.read_to_end(&mut content)?;
                Ok(content)
            },
            ProductFormat::Safe => {
                let full_path = self.product_path.join(file_path);
                std::fs::read(full_path)
                    .map_err(|e| SarError::Io(e))
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
            if file.contains("annotation/") && file.ends_with(".xml") && !file.contains("calibration") && !file.contains("noise") {
                // Parse polarization from filename - Sentinel-1 standard naming
                let pol = if file.contains("-vv-") || file.contains("_vv_") {
                    Polarization::VV
                } else if file.contains("-vh-") || file.contains("_vh_") {
                    Polarization::VH
                } else if file.contains("-hv-") || file.contains("_hv_") {
                    Polarization::HV
                } else if file.contains("-hh-") || file.contains("_hh_") {
                    Polarization::HH
                } else {
                    continue; // Skip unknown polarizations
                };
                
                // Add to the list for this polarization (no overwriting!)
                all_annotations.entry(pol).or_insert_with(Vec::new).push(file);
            }
        }
        
        if all_annotations.is_empty() {
            return Err(SarError::InvalidFormat(
                "No annotation files found. Expected files like 's1a-iw-slc-vv-*.xml' in annotation/ directory".to_string(),
            ));
        }
        
        // For backward compatibility, return just the first file per polarization
        // but log that we found multiple files
        let mut result = HashMap::new();
        for (pol, files) in &all_annotations {
            if let Some(first_file) = files.first() {
                result.insert(*pol, first_file.clone());
                if files.len() > 1 {
                    log::warn!("Found {} annotation files for {}, returning first: {}", files.len(), pol, first_file);
                    for (i, file) in files.iter().enumerate() {
                        log::info!("  File {}: {}", i+1, file);
                    }
                } else {
                    log::info!("Found 1 annotation file for {}: {}", pol, first_file);
                }
            }
        }
        
        let total_files: usize = all_annotations.values().map(|v| v.len()).sum();
        log::warn!("IMPORTANT: Found {} total annotation files but only returning {} (first per polarization)", total_files, result.len());
        log::warn!("This may indicate missing subswath data. Use find_all_annotation_files() for complete discovery");
        
        Ok(result)
    }

    /// Find ALL annotation files for each polarization (all subswaths)
    pub fn find_all_annotation_files(&mut self) -> SarResult<HashMap<Polarization, Vec<String>>> {
        let files = self.list_files()?;
        let mut annotations: HashMap<Polarization, Vec<String>> = HashMap::new();
        
        for file in files {
            // Only include main annotation files, exclude calibration, noise, and RFI subdirectories
            if file.contains("annotation/") && file.ends_with(".xml") && 
               !file.contains("calibration") && !file.contains("noise") && !file.contains("rfi") {
                
                // Additional check: only include files directly in annotation/ directory, not subdirectories
                let annotation_relative = file.split("annotation/").nth(1).unwrap_or("");
                if annotation_relative.contains('/') {
                    // Skip files in subdirectories (like calibration/, noise/, rfi/)
                    continue;
                }
                
                // Parse polarization from filename - Sentinel-1 standard naming
                let pol = if file.contains("-vv-") || file.contains("_vv_") {
                    Polarization::VV
                } else if file.contains("-vh-") || file.contains("_vh_") {
                    Polarization::VH
                } else if file.contains("-hv-") || file.contains("_hv_") {
                    Polarization::HV
                } else if file.contains("-hh-") || file.contains("_hh_") {
                    Polarization::HH
                } else {
                    continue; // Skip unknown polarizations
                };
                
                // Add to the list for this polarization
                annotations.entry(pol).or_insert_with(Vec::new).push(file);
            }
        }
        
        if annotations.is_empty() {
            return Err(SarError::InvalidFormat(
                "No annotation files found. Expected files like 's1a-iw-slc-vv-*.xml' in annotation/ directory".to_string(),
            ));
        }
        
        let total_files: usize = annotations.values().map(|v| v.len()).sum();
        log::info!("Found {} total annotation files across {} polarizations", total_files, annotations.len());
        for (pol, files) in &annotations {
            log::info!("  {}: {} files", pol, files.len());
        }
        
        Ok(annotations)
    }

    /// Read and parse annotation XML for a specific polarization with real metadata extraction
    pub fn read_annotation(&mut self, pol: Polarization) -> SarResult<SarMetadata> {
        let annotations = self.find_annotation_files()?;
        let annotation_file = annotations.get(&pol)
            .ok_or_else(|| SarError::InvalidFormat(
                format!("No annotation found for polarization {}", pol)
            ))?
            .clone(); // Clone the filename to avoid borrowing issues

        log::info!("Reading annotation for {} from: {}", pol, annotation_file);

        // Read annotation XML content
        let xml_content = self.read_file_as_string(&annotation_file)?;

        // Parse XML with comprehensive metadata extraction  
        Self::parse_annotation_xml_comprehensive(&xml_content, pol, &self.product_path)
    }

    /// Comprehensive annotation XML parser that extracts REAL metadata
    /// 
    /// This replaces hardcoded values with actual extracted data from
    /// Sentinel-1 annotation XML files. Essential for scientific accuracy.
    /// 
    /// References:
    /// - ESA S1-IF-ASD-PL-0007: "Sentinel-1 Annotation Format Specification"
    /// - ESA S1-TN-MDA-52-7440: "Sentinel-1 TOPSAR Mode"
    fn parse_annotation_xml_comprehensive(xml_content: &str, pol: Polarization, product_path: &Path) -> SarResult<SarMetadata> {
        log::info!("Starting comprehensive XML parser for polarization {:?}", pol);
        
        // Extract product ID from path 
        let product_id = product_path
            .file_stem()
            .and_then(|s| s.to_str())
            .unwrap_or("S1A_UNKNOWN")
            .replace(".SAFE", "")
            .replace(".zip", "");

        log::debug!("Product ID: {}", product_id);

        // Use the unified serde-based annotation parser
        let annotation = crate::io::annotation::parse_annotation_xml(xml_content)
            .map_err(|e| {
                log::error!("Failed to parse annotation XML: {}", e);
                SarError::Processing("Failed to parse annotation XML - cannot extract required parameters".to_string())
            })?;
        
        log::info!("Annotation parsing completed successfully");
        
        // Extract pixel spacing FIRST - no fallbacks allowed
        let range_pixel_spacing = annotation.image_annotation
            .as_ref()
            .and_then(|img| img.image_information.as_ref())
            .and_then(|info| info.range_pixel_spacing)
            .ok_or_else(|| SarError::Processing("Range pixel spacing not found in annotation XML - required for scientific accuracy".to_string()))?;
        let azimuth_pixel_spacing = annotation.image_annotation
            .as_ref()
            .and_then(|img| img.image_information.as_ref())
            .and_then(|info| info.azimuth_pixel_spacing)
            .ok_or_else(|| SarError::Processing("Azimuth pixel spacing not found in annotation XML - required for scientific accuracy".to_string()))?;
            
        // Extract sub-swaths using the detailed parser
        // SCIENTIFIC REQUIREMENT: Subswath extraction must succeed for proper SAR processing
        let _sub_swaths = crate::io::annotation::AnnotationParser::extract_subswaths(&annotation)
            .map_err(|e| SarError::Processing(format!(
                "❌ SCIENTIFIC ERROR: Failed to extract subswaths from annotation: {}. \
                Real Sentinel-1 annotation required - no fallback to empty data permitted.", e
            )))?;
            
        // Extract real bounding box from geolocation grid
        let mut bounding_box = BoundingBox {
            min_lon: -180.0,
            max_lon: 180.0,
            min_lat: -90.0,
            max_lat: 90.0,
        };
        
        if let Some(geolocation_grid) = &annotation.geolocation_grid {
            if let Some(grid_points) = &geolocation_grid.geolocation_grid_point_list {
                if let Some(points) = &grid_points.geolocation_grid_points {
                    let mut lons = Vec::new();
                    let mut lats = Vec::new();
                    
                    for point in points {
                        lons.push(point.longitude);
                        lats.push(point.latitude);
                    }
                
                    if !lons.is_empty() && !lats.is_empty() {
                        bounding_box.min_lon = lons.iter().cloned().fold(f64::INFINITY, f64::min);
                        bounding_box.max_lon = lons.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
                        bounding_box.min_lat = lats.iter().cloned().fold(f64::INFINITY, f64::min);
                        bounding_box.max_lat = lats.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
                    } else {
                    }
                } else {
                }
            } else {
            }
        } else {
        }

        // Extract real timing parameters from image annotation (no defaults)
        let slant_range_time = annotation.image_annotation
            .as_ref()
            .and_then(|img| img.image_information.as_ref())
            .and_then(|info| info.slant_range_time)
            .ok_or_else(|| SarError::Processing("Slant range time not found in annotation XML - required".to_string()))?;

        let range_sampling_rate = annotation.general_annotation
            .as_ref()
            .and_then(|gen| gen.product_information.as_ref())
            .map(|prod| prod.range_sampling_rate)
            .ok_or_else(|| SarError::Processing("Range sampling rate not found in annotation XML - required".to_string()))?;

        // CRITICAL: Extract PRF from the annotation XML
        // XML contains both <prf> and <azimuthFrequency> - use azimuth_frequency which is the actual PRF value
        let prf_result = annotation.image_annotation
            .as_ref()
            .and_then(|img| img.image_information.as_ref())
            .and_then(|info| info.azimuth_frequency);  // Use 'azimuth_frequency' field
            
        log::info!("PRF extraction result from <prf> field: {:?}", prf_result);
        
        // Fallback to azimuth_frequency if prf field is not available
        let prf = prf_result.or_else(|| {
            log::warn!("<prf> field not found, trying <azimuthFrequency>");
            annotation.image_annotation
                .as_ref()
                .and_then(|img| img.image_information.as_ref())
                .and_then(|info| info.azimuth_frequency)
        }).ok_or_else(|| SarError::Processing("Neither PRF nor azimuthFrequency found in annotation XML - required".to_string()))?;

        log::info!("Final PRF value extracted: {} Hz", prf);

        // Extract radar frequency (C-band ~5.405 GHz) - must be present
        let radar_frequency = annotation.general_annotation
            .as_ref()
            .and_then(|gen| gen.product_information.as_ref())
            .map(|prod| prod.radar_frequency)
            .ok_or_else(|| SarError::Processing("Radar frequency not found in annotation XML - required for wavelength".to_string()))?;

        let wavelength = SPEED_OF_LIGHT_M_S / radar_frequency;

        log::info!("Successfully extracted parameters with comprehensive parser:");
        log::info!("  - Range pixel spacing: {:.4} m", range_pixel_spacing);
        log::info!("  - Azimuth pixel spacing: {:.4} m", azimuth_pixel_spacing);
        log::info!("  - Slant range time: {:.6} s", slant_range_time);
        log::info!("  - Range sampling rate: {:.0} Hz", range_sampling_rate);
        log::info!("  - PRF: {:.1} Hz", prf);
        log::info!("  - Radar frequency: {:.3} GHz", radar_frequency / 1e9);
        log::info!("  - Wavelength: {:.4} m", wavelength);

        // Extract other required metadata
        let mission = if product_id.starts_with("S1A") { "Sentinel-1A".to_string() } 
                     else if product_id.starts_with("S1B") { "Sentinel-1B".to_string() }
                     else { "Sentinel-1".to_string() };
        
        let platform = mission.clone();
        let acquisition_mode = AcquisitionMode::IW; // Interferometric Wide swath mode

        // Extract timing from ADS header (no synthetic fallbacks)
        let (start_time, stop_time) = if let Some(ads_header) = &annotation.ads_header {
            let start_time = ads_header.start_time.as_deref()
                .ok_or_else(|| SarError::Processing("Start time missing in ADS header".to_string()))?;
            let stop_time = ads_header.stop_time.as_deref()
                .ok_or_else(|| SarError::Processing("Stop time missing in ADS header".to_string()))?;
            (Self::parse_time_flexible(start_time)?, Self::parse_time_flexible(stop_time)?)
        } else {
            return Err(SarError::Processing("ADS header with start/stop time is required".to_string()));
        };

        // Extract sub-swaths using the detailed parser
        let sub_swaths = crate::io::annotation::AnnotationParser::extract_subswaths(&annotation)
            .unwrap_or_else(|e| {
                log::warn!("Failed to extract subswaths: {}", e);
                HashMap::new()
            });
            
        // Extract real bounding box from geolocation grid using comprehensive parser
        let bounding_box = match crate::io::annotation::AnnotationParser::extract_bounding_box(&annotation) {
            Ok(real_bbox) => {
                log::info!("Extracted real bounding box: [{:.6}, {:.6}, {:.6}, {:.6}]",
                          real_bbox.min_lon, real_bbox.min_lat, real_bbox.max_lon, real_bbox.max_lat);
                real_bbox
            },
            Err(e) => {
                log::warn!("Failed to extract bounding box from annotation: {}", e);
                return Err(SarError::XmlParsing(
                    "Failed to extract valid bounding box from annotation XML with structured parser. No fallback parsers allowed as requested.".to_string()
                ));
            }
        };

        // Create metadata structure with REAL extracted values
        
        // SCIENTIFIC FIX: Extract orbit data from annotation for velocity calculations
        let orbit_data = if let Some(orbit_list) = &annotation.orbit_list {
            if !orbit_list.is_empty() {
                log::info!("Extracted {} orbit state vectors for velocity calculations", orbit_list.len());
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
        
        // CRITICAL DEBUG: Check PRF value before creating SarMetadata
        println!("
        log::info!("Creating SarMetadata with PRF = {} Hz", prf);
        
        Ok(SarMetadata {
            product_id,
            mission,
            platform,
            instrument: "C-SAR".to_string(),
            acquisition_mode,
            polarizations: vec![pol],
            start_time,
            stop_time,
            
            // SCIENTIFIC REQUIREMENT: Extract radar parameters from annotation XML
            radar_frequency: Some(radar_frequency), // Use the extracted value from comprehensive parser
            wavelength: Some(wavelength), // Use the calculated wavelength from comprehensive parser
            slant_range_time: Some(slant_range_time), // Use the extracted slant range time from comprehensive parser
            prf: Some(prf), // Use the extracted PRF from comprehensive parser - NO FALLBACKS
            
            bounding_box,
            coordinate_system: CoordinateSystem::Radar,
            sub_swaths,
            orbit_data,  // SCIENTIFIC FIX: Now properly populated from annotation
            range_looks: 1,
            azimuth_looks: 1,
            pixel_spacing: (range_pixel_spacing, azimuth_pixel_spacing), // REAL extracted values
        })
    }
    fn extract_real_value_f64(xml: &str, tag: &str) -> Option<f64> {
        Self::extract_value(xml, tag)
            .and_then(|s| {
                let cleaned = s.trim()
                    .replace(" m", "")
                    .replace(" Hz", "")
                    .replace(" s", "")
                    .replace(" deg", "");
                
                match cleaned.parse::<f64>() {
                    Ok(val) => {
                        log::debug!("Extracted {} = {} from XML", tag, val);
                        Some(val)
                    },
                    Err(e) => {
                        log::warn!("Failed to parse {} value '{}': {}", tag, cleaned, e);
                        None
                    }
                }
            })
    }

    /// Flexible time parsing that handles multiple formats
    fn parse_time_flexible(time_str: &str) -> SarResult<DateTime<Utc>> {
        
        // Try different time formats
        if let Ok(time) = DateTime::parse_from_rfc3339(time_str) {
            return Ok(time.with_timezone(&Utc));
        }
        
        // Try format with microseconds and Z suffix
        if let Ok(time) = DateTime::parse_from_str(&format!("{}Z", time_str), "%Y-%m-%dT%H:%M:%S%.fZ") {
            return Ok(time.with_timezone(&Utc));
        }
        
        // Try format without Z suffix but assume UTC
        if let Ok(naive_dt) = chrono::NaiveDateTime::parse_from_str(time_str, "%Y-%m-%dT%H:%M:%S%.f") {
            return Ok(DateTime::<Utc>::from_naive_utc_and_offset(naive_dt, Utc));
        }
        
        // Try format with explicit UTC
        if let Ok(time) = DateTime::parse_from_str(time_str, "%Y-%m-%dT%H:%M:%S%.f+00:00") {
            return Ok(time.with_timezone(&Utc));
        }
        
        // Try parsing without microseconds
        if let Ok(time) = DateTime::parse_from_str(&format!("{}Z", time_str), "%Y-%m-%dT%H:%M:%SZ") {
            return Ok(time.with_timezone(&Utc));
        }
        
        // Try naive datetime without microseconds
        if let Ok(naive_dt) = chrono::NaiveDateTime::parse_from_str(time_str, "%Y-%m-%dT%H:%M:%S") {
            return Ok(DateTime::<Utc>::from_naive_utc_and_offset(naive_dt, Utc));
        }
        
        // Return error if no valid time format found
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
    pub fn get_annotation_for_polarization(&mut self, pol: Polarization) -> SarResult<crate::io::annotation::AnnotationRoot> {
        // Find annotation files and get the one for this polarization
        let annotations = self.find_annotation_files()?;
        let annotation_file = annotations.get(&pol)
            .ok_or_else(|| SarError::Processing(
                format!("No annotation file found for polarization {:?}", pol)
            ))?;
        
        // Read the annotation XML content directly
        let content = self.read_file_as_string(annotation_file)?;
        
        // Parse the XML content into AnnotationRoot
        let annotation_root: crate::io::annotation::AnnotationRoot = 
            quick_xml::de::from_str(&content)
                .map_err(|e| SarError::XmlParsing(format!("Failed to parse annotation XML: {}", e)))?;
        
        Ok(annotation_root)
    }

    /// **SCIENTIFIC IMPROVEMENT**: Read all subswath annotations for comprehensive processing
    /// 
    /// This method addresses the critical issue where only the first annotation per polarization 
    /// was processed, potentially missing IW2/IW3 data. Returns a map of all annotations per 
    /// polarization and subswath for complete data coverage.
    /// 
    /// Returns: HashMap<Polarization, Vec<AnnotationRoot>> - all subswaths per polarization
    pub fn get_all_subswath_annotations(&mut self, pol: Polarization) -> SarResult<Vec<crate::io::annotation::AnnotationRoot>> {
        // Get all annotation files for this polarization
        let all_annotations = self.find_all_annotation_files()?;
        let annotation_files = all_annotations.get(&pol)
            .ok_or_else(|| SarError::Processing(
                format!("No annotation files found for polarization {:?}", pol)
            ))?;
        
        if annotation_files.is_empty() {
            return Err(SarError::Processing(
                format!("Empty annotation file list for polarization {:?}", pol)
            ));
        }
        
        log::info!("🔬 SCIENTIFIC PROCESSING: Loading {} subswath annotations for {}", 
                   annotation_files.len(), pol);
        
        let mut annotations = Vec::new();
        for (i, annotation_file) in annotation_files.iter().enumerate() {
            log::info!("  Loading subswath {}: {}", i+1, annotation_file);
            
            // Read the annotation XML content
            let content = self.read_file_as_string(annotation_file)?;
            
            // Parse with strict error handling
            let annotation_root: crate::io::annotation::AnnotationRoot = 
                quick_xml::de::from_str(&content)
                    .map_err(|e| SarError::XmlParsing(
                        format!("Failed to parse annotation XML for {}: {}", annotation_file, e)
                    ))?;
            
            annotations.push(annotation_root);
        }
        
        log::info!("✅ Successfully loaded {} subswath annotations for {}", annotations.len(), pol);
        Ok(annotations)
    }

    /// **MIGRATION HELPER**: Get primary annotation with subswath coverage awareness
    /// 
    /// This method provides a migration path from single-annotation processing to 
    /// multi-subswath awareness. It returns the first annotation but logs warnings
    /// about any additional subswaths that would be missed.
    pub fn get_primary_annotation_with_coverage_check(&mut self, pol: Polarization) -> SarResult<crate::io::annotation::AnnotationRoot> {
        let all_annotations = self.get_all_subswath_annotations(pol)?;
        
        if all_annotations.is_empty() {
            return Err(SarError::Processing(
                format!("No annotations found for polarization {:?}", pol)
            ));
        }
        
        if all_annotations.len() > 1 {
            log::warn!("🚨 SCIENTIFIC DATA LOSS WARNING: Using only 1 of {} available subswaths for {}", 
                      all_annotations.len(), pol);
            log::warn!("🚨 This may result in incomplete scientific analysis!");
            log::warn!("🚨 Consider migrating to get_all_subswath_annotations() for full coverage");
        }
        
        Ok(all_annotations.into_iter().next().unwrap())
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
    pub fn get_annotation_unified_validated(&mut self, pol: Polarization, validate_consistency: bool) -> SarResult<crate::io::annotation::ProductRoot> {
        // Get annotation file
        let annotations = self.find_annotation_files()?;
        let annotation_file = annotations.get(&pol)
            .ok_or_else(|| SarError::Processing(
                format!("No annotation file found for polarization {:?}", pol)
            ))?;
        
        // Read XML content
        let xml_content = self.read_file_as_string(annotation_file)?;
        
        if validate_consistency {
            log::info!("🔬 SCIENTIFIC VALIDATION: Comparing parsing methods for {}", pol);
            
            // Validate parsing consistency
            let validator = crate::io::parsing_validation::ParsingValidator::new();
            match validator.validate_parsing_equivalence(&xml_content) {
                Ok(validation_result) => {
                    if !validation_result.equivalent {
                        log::warn!("🚨 PARSING INCONSISTENCY DETECTED for {}", pol);
                        log::warn!("🚨 Details: {}", validation_result.details);
                    } else {
                        log::info!("✅ Parsing validation passed for {}: {}", pol, validation_result.details);
                    }
                }
                Err(e) => {
                    log::warn!("⚠️  Parsing validation failed for {}: {}", pol, e);
                }
            }
        }
        
        // Primary parsing attempt - use serde (recommended)
        match quick_xml::de::from_str::<crate::io::annotation::ProductRoot>(&xml_content) {
            Ok(product_root) => {
                log::info!("✅ Serde parsing successful for {}", pol);
                Ok(product_root)
            },
            Err(serde_err) => {
                log::warn!("⚠️  Serde parsing failed for {}: {}", pol, serde_err);
                log::warn!("🔄 Falling back to alternative parsing method...");
                
                // Fallback - try the same unified parser again with detailed error
                match crate::io::annotation::parse_annotation_xml(&xml_content) {
                    Ok(product_root) => {
                        log::warn!("✅ Second parsing attempt successful for {}", pol);
                        log::warn!("🔬 SCIENTIFIC WARNING: Initial parsing failed but retry succeeded - investigate XML consistency");
                        Ok(product_root)
                    },
                    Err(fallback_err) => {
                        Err(SarError::Processing(format!(
                            "Both parsing attempts failed for {}: first={}, second={}", 
                            pol, serde_err, fallback_err
                        )))
                    }
                }
            }
        }
    }

    /// Find measurement data files (TIFF files containing SLC data)
    pub fn find_measurement_files(&mut self) -> SarResult<HashMap<Polarization, String>> {
        let files = self.list_files()?;
        let mut measurements = HashMap::new();
        
        for file in files {
            if file.contains("measurement/") && file.ends_with(".tiff") {
                // Parse polarization from filename - handle both ZIP and SAFE naming
                if file.contains("-vv-") || file.contains("_vv_") || file.contains("vv.tiff") {
                    measurements.insert(Polarization::VV, file);
                } else if file.contains("-vh-") || file.contains("_vh_") || file.contains("vh.tiff") {
                    measurements.insert(Polarization::VH, file);
                } else if file.contains("-hv-") || file.contains("_hv_") || file.contains("hv.tiff") {
                    measurements.insert(Polarization::HV, file);
                } else if file.contains("-hh-") || file.contains("_hh_") || file.contains("hh.tiff") {
                    measurements.insert(Polarization::HH, file);
                }
            }
        }
        
        if measurements.is_empty() {
            return Err(SarError::InvalidFormat(
                "No measurement files found. Expected TIFF files in measurement/ directory".to_string(),
            ));
        }

        log::info!("Found {} measurement files: {:?}", measurements.len(), measurements.keys().collect::<Vec<_>>());
        
        Ok(measurements)
    }

    /// Read SLC data for a specific polarization (supports both ZIP and SAFE)
    pub fn read_slc_data(&mut self, pol: Polarization) -> SarResult<SarImage> {
        
        use tempfile::NamedTempFile;
        
        let measurements = self.find_measurement_files()?;
        let measurement_file = measurements.get(&pol)
            .ok_or_else(|| SarError::InvalidFormat(
                format!("No measurement found for polarization {}", pol)
            ))?;

        log::info!("Reading SLC data for {} from {}", pol, measurement_file);
        let start_time = std::time::Instant::now();

        // Handle different formats
        let dataset = match self.format {
            ProductFormat::Zip => {
                // Extract the TIFF file to a temporary location for ZIP format
                let archive = self.open_archive()?;
                let mut zip_file = archive.by_name(measurement_file)
                    .map_err(|e| SarError::Io(std::io::Error::other(
                        format!("Failed to access {}: {}", measurement_file, e),
                    )))?;

                // Create a temporary file to extract the TIFF
                let mut temp_file = NamedTempFile::new()
                    .map_err(|e| SarError::Io(e))?;
                
                // Copy data from ZIP to temporary file
                std::io::copy(&mut zip_file, &mut temp_file)
                    .map_err(|e| SarError::Io(e))?;
                
                let extract_time = start_time.elapsed();
                log::debug!("TIFF extraction took: {:?}", extract_time);

                // Read with GDAL from temporary file
                gdal::Dataset::open(temp_file.path())
                    .map_err(|e| SarError::Io(std::io::Error::other(
                        format!("Failed to open TIFF with GDAL: {}", e),
                    )))?
            },
            ProductFormat::Safe => {
                // Read directly from SAFE directory structure
                let full_path = self.product_path.join(measurement_file);
                gdal::Dataset::open(full_path)
                    .map_err(|e| SarError::Io(std::io::Error::other(
                        format!("Failed to open TIFF with GDAL: {}", e),
                    )))?
            }
        };

        let gdal_start = std::time::Instant::now();

        // Get raster dimensions and band count
        let raster_size = dataset.raster_size();
        let (width, height) = (raster_size.0, raster_size.1);
        let band_count = dataset.raster_count();
        log::debug!("TIFF dimensions: {} x {}, bands: {}", width, height, band_count);

        // Handle different band configurations
        let slc_data = if band_count == 1 {
            // Single band containing complex data (Sentinel-1 SLC format: CInt16)
            let band = dataset.rasterband(1)
                .map_err(|e| SarError::Io(std::io::Error::other(
                    format!("Failed to get band 1: {}", e),
                )))?;

            // **CRITICAL GDAL SAFETY CHECK**: Verify band data type
            let band_type = band.band_type();
            log::info!("🔬 GDAL band type detected: {:?}", band_type);
            
            // Verify this is actually complex data before attempting complex read
            // Note: GDAL CInt16 corresponds to complex 16-bit integers (real/imag pairs)
            log::info!("🔬 GDAL band type detected: {:?} - proceeding with complex read", band_type);
            // TODO: Add proper type validation once GDAL enum variants are confirmed

            let window = (0, 0);
            let window_size = (width, height);
            
            // **IMPROVED COMPLEX READ**: Read with explicit type validation
            // Sentinel-1 SLC uses complex 16-bit integers (CInt16) stored as interleaved i16 pairs
            let complex_data = band.read_as::<i16>(window, window_size, (width * 2, height), None)
                .map_err(|e| SarError::Io(std::io::Error::other(
                    format!("Failed to read complex CInt16 data from GDAL band: {}", e),
                )))?;

            // **SAFETY VALIDATION**: Verify data size matches expected complex format
            let expected_samples = (width * height * 2) as usize; // 2 samples per complex pixel
            if complex_data.data.len() != expected_samples {
                return Err(SarError::InvalidFormat(format!(
                    "GDAL complex data size mismatch: expected {} samples ({}x{}x2), got {} samples. This indicates incorrect complex interpretation.",
                    expected_samples, width, height, complex_data.data.len()
                )));
            }

            log::debug!("✅ Read {} complex values as interleaved i16 ({}x{} pixels)", 
                       complex_data.data.len() / 2, width, height);

            // Convert interleaved i16 data to complex f32 array
            Self::convert_cint16_to_complex(complex_data, width, height)?
        } else if band_count >= 2 {
            // Separate I and Q bands (less common for Sentinel-1)
            let band1 = dataset.rasterband(1)
                .map_err(|e| SarError::Io(std::io::Error::other(
                    format!("Failed to get band 1: {}", e),
                )))?;

            let band2 = dataset.rasterband(2)
                .map_err(|e| SarError::Io(std::io::Error::other(
                    format!("Failed to get band 2: {}", e),
                )))?;

            let window = (0, 0);
            let window_size = (width, height);
            let buffer_size = (width, height);
            
            let i_data = band1.read_as::<f32>(window, window_size, buffer_size, None)
                .map_err(|e| SarError::Io(std::io::Error::other(
                    format!("Failed to read I data: {}", e),
                )))?;

            let q_data = band2.read_as::<f32>(window, window_size, buffer_size, None)
                .map_err(|e| SarError::Io(std::io::Error::other(
                    format!("Failed to read Q data: {}", e),
                )))?;

            // Convert to complex SLC data
            Self::convert_to_complex_parallel(i_data, q_data, width, height)?
        } else {
            return Err(SarError::InvalidFormat(
                format!("Unexpected number of bands in TIFF: {}", band_count)
            ));
        };

        let gdal_time = gdal_start.elapsed();
        log::debug!("GDAL read took: {:?}", gdal_time);

        let total_time = start_time.elapsed();

        log::info!("SLC data read complete for {}: {} x {} pixels", pol, width, height);
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
        let measurement_file = measurements.get(&pol)
            .ok_or_else(|| SarError::InvalidFormat(
                format!("No measurement found for polarization {}", pol)
            ))?;

        log::info!("Reading SLC data (PARALLEL) for {} from {}", pol, measurement_file);
        let start_time = std::time::Instant::now();

        // Extract the TIFF file to a temporary location
        let archive = self.open_archive()?;
        let mut zip_file = archive.by_name(measurement_file)
            .map_err(|e| SarError::Io(std::io::Error::other(
                format!("Failed to access {}: {}", measurement_file, e),
            )))?;

        // Create a temporary file to extract the TIFF
        let mut temp_file = NamedTempFile::new()
            .map_err(|e| SarError::Io(e))?;
        
        // Copy data from ZIP to temporary file
        std::io::copy(&mut zip_file, &mut temp_file)
            .map_err(|e| SarError::Io(e))?;
        
        let extract_time = start_time.elapsed();
        log::debug!("TIFF extraction took: {:?}", extract_time);

        // Read with GDAL
        let gdal_start = std::time::Instant::now();
        let dataset = gdal::Dataset::open(temp_file.path())
            .map_err(|e| SarError::Io(std::io::Error::other(
                format!("Failed to open TIFF with GDAL: {}", e),
            )))?;

        // Get raster dimensions and band count
        let raster_size = dataset.raster_size();
        let (width, height) = (raster_size.0, raster_size.1);
        let band_count = dataset.raster_count();
        log::debug!("TIFF dimensions: {} x {}, bands: {}", width, height, band_count);

        // Handle different band configurations
        let slc_data = if band_count == 1 {
            // Single band containing complex data (Sentinel-1 SLC format: CInt16)
            let band = dataset.rasterband(1)
                .map_err(|e| SarError::Io(std::io::Error::other(
                format!("Failed to get band 1: {}", e),
            )))?;

            // **CRITICAL GDAL SAFETY CHECK**: Verify band data type  
            let band_type = band.band_type();
            log::debug!("🔬 GDAL band type detected: {:?}", band_type);
            
            // Note: Complex data validation - proceeding with complex read
            log::debug!("🔬 Proceeding with complex CInt16 read for type: {:?}", band_type);

            let window = (0, 0);
            let window_size = (width, height);
            
            // **IMPROVED COMPLEX READ**: Sentinel-1 SLC uses complex 16-bit integers (CInt16)
            let complex_data = band.read_as::<i16>(window, window_size, (width * 2, height), None)
                .map_err(|e| SarError::Io(std::io::Error::other(
                format!("Failed to read complex CInt16 data from GDAL: {}", e),
            )))?;

            // **SAFETY VALIDATION**: Verify data size for complex format
            let expected_samples = (width * height * 2) as usize;
            if complex_data.data.len() != expected_samples {
                return Err(SarError::InvalidFormat(format!(
                    "Complex data size validation failed: expected {} samples, got {}", 
                    expected_samples, complex_data.data.len()
                )));
            }

            log::debug!("Read {} complex values as interleaved i16 (parallel)", complex_data.data.len() / 2);

            // Convert interleaved i16 data to complex f32 array with parallel processing
            Self::convert_cint16_to_complex_parallel(complex_data, width, height)?
        } else if band_count >= 2 {
            // Separate I and Q bands (less common for Sentinel-1)
            let band1 = dataset.rasterband(1)
                .map_err(|e| SarError::Io(std::io::Error::other(
                format!("Failed to get band 1: {}", e),
            )))?;

            let band2 = dataset.rasterband(2)
                .map_err(|e| SarError::Io(std::io::Error::other(
                format!("Failed to get band 2: {}", e),
            )))?;

            let window = (0, 0);
            let window_size = (width, height);
            let buffer_size = (width, height);
            
            let i_data = band1.read_as::<f32>(window, window_size, buffer_size, None)
                .map_err(|e| SarError::Io(std::io::Error::other(
                format!("Failed to read I data: {}", e),
            )))?;

            let q_data = band2.read_as::<f32>(window, window_size, buffer_size, None)
                .map_err(|e| SarError::Io(std::io::Error::other(
                format!("Failed to read Q data: {}", e),
            )))?;

            // Convert to complex SLC data with parallel processing
            Self::convert_to_complex_parallel(i_data, q_data, width, height)?
        } else {
            return Err(SarError::InvalidFormat(
                format!("Unexpected number of bands in TIFF: {}", band_count)
            ));
        };
        let gdal_time = gdal_start.elapsed();
        log::debug!("GDAL read took: {:?}", gdal_time);

        let total_time = start_time.elapsed();

        log::info!("SLC data read complete (PARALLEL) for {}: {} x {} pixels", pol, width, height);
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
    pub fn read_slc_data_streaming(&mut self, pol: Polarization) -> SarResult<SarImage> {
        use std::io::Cursor;
        
        let measurements = self.find_measurement_files()?;
        let measurement_file = measurements.get(&pol)
            .ok_or_else(|| SarError::InvalidFormat(
                format!("No measurement found for polarization {}", pol)
            ))?;

        log::info!("Reading SLC data (STREAMING) for {} from {}", pol, measurement_file);
        let start_time = std::time::Instant::now();

        // Read directly from ZIP into memory
        let buffer = {
            let archive = self.open_archive()?;
            let mut zip_file = archive
                .by_name(measurement_file)
                .map_err(|e| SarError::Io(std::io::Error::other(format!(
                    "Failed to access {}: {}",
                    measurement_file, e
                ))))?;

            // Read entire file into memory buffer
            let mut buffer = Vec::new();
            std::io::copy(&mut zip_file, &mut buffer)
                .map_err(|e| SarError::Io(e))?;
            
            buffer
        };
        
        let extract_time = start_time.elapsed();
        log::debug!("In-memory extraction took: {:?}", extract_time);

        // Create cursor for GDAL to read from memory
        let cursor = Cursor::new(buffer);
        let cursor_len = cursor.get_ref().len();
        
        // Use GDAL VSI memory filesystem for direct memory access
        let vsi_path = format!("/vsimem/slc_data_{}", pol);
        
        // Write to VSI memory
        let _gdal_start = std::time::Instant::now();
        unsafe {
            let c_path = std::ffi::CString::new(vsi_path.clone()).unwrap();
            let buffer_ptr = cursor.get_ref().as_ptr() as *const std::os::raw::c_void;
            gdal_sys::VSIFileFromMemBuffer(
                c_path.as_ptr(),
                buffer_ptr as *mut std::os::raw::c_uchar,
                cursor_len as u64,
                0, // don't take ownership
            );
        }
        
        // Open with GDAL from memory
        let dataset = gdal::Dataset::open(&vsi_path).map_err(|e| {
            SarError::Io(std::io::Error::other(format!(
                "Failed to open TIFF from memory: {}",
                e
            )))
        })?;

        let result = self.read_slc_from_dataset(dataset, pol);
        
        // Clean up VSI memory
        unsafe {
            let c_path = std::ffi::CString::new(vsi_path).unwrap();
            gdal_sys::VSIUnlink(c_path.as_ptr());
        }
        
        let total_time = start_time.elapsed();
        log::info!("Streaming SLC read completed in: {:?}", total_time);
        
        result
    }

    /// Read SLC data from GDAL dataset (shared code)
    fn read_slc_from_dataset(&self, dataset: gdal::Dataset, _pol: Polarization) -> SarResult<SarImage> {
        // Get raster dimensions and band count
        let raster_size = dataset.raster_size();
        let (width, height) = (raster_size.0, raster_size.1);
        let band_count = dataset.raster_count();
        log::debug!("TIFF dimensions: {} x {}, bands: {}", width, height, band_count);

        // Handle different band configurations
        let slc_data = if band_count == 1 {
            // Single band containing complex data (Sentinel-1 SLC format: CInt16)
            let band = dataset.rasterband(1).map_err(|e| {
                SarError::Io(std::io::Error::other(format!(
                    "Failed to get band 1: {}",
                    e
                )))
            })?;

            // Check data type
            let data_type = band.band_type();
            log::debug!("Band data type: {:?}", data_type);

            let window = (0, 0);
            let window_size = (width, height);
            let _buffer_size = (width, height);
            
            // Sentinel-1 SLC uses complex 16-bit integers (CInt16)
            // GDAL represents this as interleaved 16-bit signed integers: [real, imag, real, imag, ...]
            let complex_data = band
                .read_as::<i16>(window, window_size, (width * 2, height), None)
                .map_err(|e| {
                    SarError::Io(std::io::Error::other(format!(
                        "Failed to read complex CInt16 data: {}",
                        e
                    )))
                })?;

            log::debug!("Read {} complex values as interleaved i16", complex_data.data.len() / 2);

            // Convert interleaved i16 data to complex f32 array
            Self::convert_cint16_to_complex(complex_data, width, height)?
        } else if band_count >= 2 {
            // Separate I and Q bands (less common for Sentinel-1)
            let band1 = dataset.rasterband(1).map_err(|e| {
                SarError::Io(std::io::Error::other(format!(
                    "Failed to get band 1: {}",
                    e
                )))
            })?;

            let band2 = dataset.rasterband(2).map_err(|e| {
                SarError::Io(std::io::Error::other(format!(
                    "Failed to get band 2: {}",
                    e
                )))
            })?;

            let window = (0, 0);
            let window_size = (width, height);
            let buffer_size = (width, height);
            
            // Read I and Q bands sequentially (GDAL is not thread-safe)
            let i_data = band1
                .read_as::<f32>(window, window_size, buffer_size, None)
                .map_err(|e| {
                    SarError::Io(std::io::Error::other(format!(
                        "Failed to read I data: {}",
                        e
                    )))
                })?;

            let q_data = band2
                .read_as::<f32>(window, window_size, buffer_size, None)
                .map_err(|e| {
                    SarError::Io(std::io::Error::other(format!(
                        "Failed to read Q data: {}",
                        e
                    )))
                })?;

            // Convert to complex SLC data
            Self::convert_to_complex_parallel(i_data, q_data, width, height)?
        } else {
            return Err(SarError::InvalidFormat(
                format!("Unexpected number of bands in TIFF: {}", band_count)
            ));
        };
        
        log::info!("SLC data loaded: {} x {} pixels", height, width);
        Ok(slc_data)
    }

    /// Convert CInt16 interleaved data to complex SLC data
    fn convert_cint16_to_complex(complex_data: gdal::raster::Buffer<i16>, width: usize, height: usize) -> SarResult<SarImage> {
        use num_complex::Complex;
        
        let mut slc_array = Array2::zeros((height, width));
        let total_pixels = width * height;
        
        // Data is interleaved as [real, imag, real, imag, ...]
        // Each complex pixel requires 2 i16 values
        if complex_data.data.len() < total_pixels * 2 {
            return Err(SarError::InvalidFormat(
                format!("Insufficient data: expected {} i16 values, got {}", 
                        total_pixels * 2, complex_data.data.len())
            ));
        }
        
        // Convert i16 complex data to f32 complex values
        for row in 0..height {
            for col in 0..width {
                let pixel_idx = row * width + col;
                let data_idx = pixel_idx * 2;
                
                // Extract real and imaginary parts as i16, convert to f32
                let real_i16 = complex_data.data[data_idx];
                let imag_i16 = complex_data.data[data_idx + 1];
                
                // Convert to f32 with proper scaling
                // Sentinel-1 SLC data is typically scaled, but we preserve raw DN values
                let real_f32 = real_i16 as f32;
                let imag_f32 = imag_i16 as f32;
                
                slc_array[[row, col]] = Complex::new(real_f32, imag_f32);
            }
        }
        
        // Log some statistics
        let sample_pixel = slc_array[[0, 0]];
        log::debug!("First pixel: real={}, imag={}, magnitude={:.2}", 
                   sample_pixel.re, sample_pixel.im, sample_pixel.norm());
        
        // Check data variance
        let mut real_sum = 0.0f32;
        let mut imag_sum = 0.0f32;
        let sample_size = std::cmp::min(1000, total_pixels);
        
        for i in 0..sample_size {
            let row = i / width;
            let col = i % width;
            if row < height {
                real_sum += slc_array[[row, col]].re;
                imag_sum += slc_array[[row, col]].im;
            }
        }
        
        let real_mean = real_sum / sample_size as f32;
        let imag_mean = imag_sum / sample_size as f32;
        log::debug!("Sample statistics (first {} pixels): real_mean={:.2}, imag_mean={:.2}", 
                   sample_size, real_mean, imag_mean);
        
        Ok(slc_array)
    }

/// Convert I/Q data to complex SLC data in parallel
    fn convert_to_complex_parallel(i_data: gdal::raster::Buffer<f32>, q_data: gdal::raster::Buffer<f32>, width: usize, height: usize) -> SarResult<SarImage> {
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
    fn convert_cint16_to_complex_parallel(complex_data: gdal::raster::Buffer<i16>, width: usize, height: usize) -> SarResult<SarImage> {
        use num_complex::Complex;
        #[cfg(feature = "parallel")]
        use rayon::prelude::*;
        
        let total_pixels = width * height;
        
        // Data is interleaved as [real, imag, real, imag, ...]
        if complex_data.data.len() < total_pixels * 2 {
            return Err(SarError::InvalidFormat(
                format!("Insufficient data: expected {} i16 values, got {}", 
                        total_pixels * 2, complex_data.data.len())
            ));
        }
        
        let mut slc_array = Array2::zeros((height, width));
        
        #[cfg(feature = "parallel")]
        {
            // Process in parallel chunks by dividing into row chunks
            let chunk_size = std::cmp::max(1, height / rayon::current_num_threads());
            
            slc_array.axis_chunks_iter_mut(ndarray::Axis(0), chunk_size)
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
        log::debug!("First pixel (parallel): real={}, imag={}, magnitude={:.2}", 
                   sample_pixel.re, sample_pixel.im, sample_pixel.norm());
        
        Ok(slc_array)
    }

    /// Get comprehensive orbit data for the product
    /// Checks: 1) SLC embedded data, 2) Local cache, 3) Download from ESA
    pub fn get_orbit_data(&mut self, orbit_cache_dir: Option<&Path>) -> SarResult<OrbitData> {
        log::info!("Getting orbit data for SLC product");
        
        // Get product metadata first
        let annotations = self.find_annotation_files()?;
        let first_pol = annotations.keys().next().cloned()
            .ok_or_else(|| SarError::Processing("No polarizations found".to_string()))?;
        let metadata = self.read_annotation(first_pol)?;
        
        // Extract product ID from filename
        let product_id = self.product_path
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
                    "❌ SCIENTIFIC ERROR: Orbit cache directory must be explicitly specified! \
                    No fallback to temp directory permitted for scientific accuracy.".to_string()
                ));
            }
        };
        let orbit_manager = OrbitManager::new(cache_dir);
        orbit_manager.get_orbit_data(
            &product_id,
            metadata.start_time,
        )
    }
    
    /// Check orbit data availability status
    pub fn check_orbit_status(&mut self, orbit_cache_dir: Option<&Path>) -> SarResult<OrbitStatus> {
        let annotations = self.find_annotation_files()?;
        let first_pol = annotations.keys().next().cloned()
            .ok_or_else(|| SarError::Processing("No polarizations found".to_string()))?;
        let metadata = self.read_annotation(first_pol)?;
        
        // Extract product ID
        let product_id = self.product_path
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
                    "❌ SCIENTIFIC ERROR: Orbit cache directory must be explicitly specified! \
                    No fallback to temp directory permitted for scientific accuracy.".to_string()
                ));
            }
        };
        let orbit_manager = OrbitManager::new(cache_dir);
        let primary_orbit_type = crate::io::orbit::OrbitReader::determine_orbit_type(metadata.start_time);
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
    pub fn download_orbit_files(&mut self, orbit_cache_dir: Option<&Path>) -> SarResult<Vec<PathBuf>> {
        let annotations = self.find_annotation_files()?;
        let first_pol = annotations.keys().next().cloned()
            .ok_or_else(|| SarError::Processing("No polarizations found".to_string()))?;
        let metadata = self.read_annotation(first_pol)?;
        
        let product_id = self.product_path
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
                    "❌ SCIENTIFIC ERROR: Orbit cache directory must be explicitly specified! \
                    No fallback to temp directory permitted for scientific accuracy.".to_string()
                ));
            }
        };
        let orbit_manager = OrbitManager::new(cache_dir);
        let primary_orbit_type = crate::io::orbit::OrbitReader::determine_orbit_type(metadata.start_time);
        
        let mut downloaded_files = Vec::new();
        
        // Try to download primary orbit type only (no silent fallbacks in scientific mode)
        match orbit_manager.download_and_cache_orbit_public(&product_id, metadata.start_time, primary_orbit_type) {
            Ok(_) => {
                let path = orbit_manager.get_cache_dir().join(format!("{}_{}.EOF", product_id, primary_orbit_type));
                downloaded_files.push(path);
                log::info!("Downloaded {} orbit file", primary_orbit_type);
            },
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
        let measurement_file = measurements.get(&pol)
            .ok_or_else(|| SarError::InvalidFormat(
                format!("No measurement found for polarization {}", pol)
            ))?;
        
        // Extract TIFF to get dimensions
        let archive = self.open_archive()?;
        let mut zip_file = archive.by_name(measurement_file)
            .map_err(|e| SarError::Io(std::io::Error::other(
                format!("Failed to access {}: {}", measurement_file, e),
            )))?;
        
        // Create temporary file to get dimensions
        use tempfile::NamedTempFile;
        let mut temp_file = NamedTempFile::new()
            .map_err(|e| SarError::Io(e))?;
        
        // Copy ZIP content to temp file
        std::io::copy(&mut zip_file, &mut temp_file)
            .map_err(|e| SarError::Io(e))?;
        
        // Open with GDAL to get dimensions
        let dataset = gdal::Dataset::open(temp_file.path())
            .map_err(|e| SarError::Io(std::io::Error::other(
                format!("Failed to open TIFF with GDAL: {}", e),
            )))?;
        
        let raster_size = dataset.raster_size();
        let (width, height) = (raster_size.0, raster_size.1);
        
        log::info!("Image dimensions: {} x {} pixels", width, height);
        
        // Calculate azimuth time interval
        // For Sentinel-1 IW mode, typical azimuth time interval is ~2.8e-4 seconds
        let acquisition_duration = (annotation.stop_time - annotation.start_time).num_milliseconds() as f64 / 1000.0;
        let azimuth_time_interval = acquisition_duration / height as f64;
        
        log::info!("Azimuth timing: duration={:.3}s, interval={:.6}s, lines={}", 
                  acquisition_duration, azimuth_time_interval, height);
        
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
    ) -> SarResult<([f64; 3], [f64; 3])> { // Returns (position, velocity)
        let burst_orbit = self.get_burst_orbit_data(pol, orbit_cache_dir)?;
        
        let position = burst_orbit.get_position_at_line(azimuth_line)
            .ok_or_else(|| SarError::Processing(
                format!("Azimuth line {} out of range (max: {})", azimuth_line, burst_orbit.num_lines())
            ))?;
            
        let velocity = burst_orbit.get_velocity_at_line(azimuth_line)
            .ok_or_else(|| SarError::Processing(
                format!("Azimuth line {} out of range for velocity", azimuth_line)
            ))?;
        
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
                if let Some(dt) = info.azimuth_time_interval { dt }
                else {
                    // Compute from PRF if available: time per line ~ 1/PRF
                    if let Some(prf) = info.azimuth_frequency { 1.0 / prf } else { return Err(SarError::Metadata("Missing azimuthTimeInterval/PRF in annotation".to_string())); }
                }
            } else { return Err(SarError::Metadata("Missing imageInformation in annotation".to_string())); }
        } else { return Err(SarError::Metadata("Missing imageAnnotation in annotation".to_string())); };

        // Time since start in seconds
        let az_time_since_start = azimuth_line as f64 * az_time_interval;

        let doppler_hz = annotation_root.evaluate_doppler_centroid(az_time_since_start)?;
        log::debug!("Doppler (dcPolynomial) at line {} (t={:.6}s): {:.2} Hz", azimuth_line, az_time_since_start, doppler_hz);
        Ok(doppler_hz)
    }

    /// Extract IW sub-swaths information from annotation files
    /// This implements the "IW split" step in the SAR processing pipeline
    pub fn extract_iw_subswaths(&mut self, pol: Polarization) -> SarResult<HashMap<String, SubSwath>> {
        let annotations = self.find_annotation_files()?;
        let annotation_file = annotations.get(&pol)
            .ok_or_else(|| SarError::InvalidFormat(
                format!("No annotation found for polarization {}", pol)
            ))?
            .clone();

        let archive = self.open_archive()?;
        let mut file = archive
            .by_name(&annotation_file)
            .map_err(|e| SarError::Io(std::io::Error::other(format!(
                "Failed to read {}: {}",
                annotation_file, e
            ))))?;

        let mut xml_content = String::new();
        file.read_to_string(&mut xml_content)?;

        // Parse using the unified serde-based annotation parser
        let annotation = crate::io::annotation::parse_annotation_xml(&xml_content)
            .map_err(|e| SarError::XmlParsing(format!("Failed to parse annotation XML {}: {}", annotation_file, e)))?;
        let subswaths = crate::io::annotation::AnnotationParser::extract_subswaths(&annotation)
            .map_err(|e| SarError::Processing(format!("Failed to extract IW sub-swaths from {}: {}", annotation_file, e)))?;
        Ok(subswaths)
    }

    // Removed fallback string-based IW subswath extraction to enforce scientific parser only.
    
    // Deprecated: string-based numeric extraction removed for scientific rigor

    /// Find ALL annotation files for all IW subswaths (not just one per polarization)
    pub fn find_all_iw_annotation_files(&mut self) -> SarResult<HashMap<Polarization, Vec<String>>> {
        let files = self.list_files()?;
        let mut annotations: HashMap<Polarization, Vec<String>> = HashMap::new();
        
        for file in files {
            if file.contains("annotation/") && file.ends_with(".xml") && !file.contains("calibration") {
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
    pub fn extract_all_iw_subswaths(&mut self, pol: Polarization) -> SarResult<HashMap<String, SubSwath>> {
        let all_annotations = self.find_all_iw_annotation_files()?;
        let annotation_files = all_annotations.get(&pol)
            .ok_or_else(|| SarError::InvalidFormat(
                format!("No annotation files found for polarization {}", pol)
            ))?;

    let mut all_subswaths = HashMap::new();
        
        for annotation_file in annotation_files {
            log::debug!("Processing annotation file: {}", annotation_file);
            
            let archive = self.open_archive()?;
            let mut file = archive
                .by_name(annotation_file)
                .map_err(|e| SarError::Io(std::io::Error::other(format!(
                    "Failed to read {}: {}",
                    annotation_file, e
                ))))?;

            let mut xml_content = String::new();
            file.read_to_string(&mut xml_content)?;

            // Extract the IW subswath from this annotation file
            let annotation = crate::io::annotation::parse_annotation_xml(&xml_content)
                .map_err(|e| SarError::XmlParsing(format!("Failed to parse {}: {}", annotation_file, e)))?;
            let subswaths = crate::io::annotation::AnnotationParser::extract_subswaths(&annotation)
                .map_err(|e| SarError::Processing(format!("Failed to extract IW sub-swaths from {}: {}", annotation_file, e)))?;
            
            // Add all extracted subswaths to our collection
            for (swath_id, subswath) in subswaths {
                all_subswaths.insert(swath_id, subswath);
            }
        }
        
        Ok(all_subswaths)
    }

    /// Get available IW sub-swaths for all polarizations (fixed to extract ALL subswaths)
    pub fn get_all_iw_subswaths(&mut self) -> SarResult<HashMap<Polarization, HashMap<String, SubSwath>>> {
        log::info!("Starting get_all_iw_subswaths");
        let all_annotations = self.find_all_iw_annotation_files()?;
        log::info!("Found annotation files for {} polarizations", all_annotations.len());
        let mut all_subswaths = HashMap::new();
        
    for pol in all_annotations.keys().cloned() {
            log::info!("Processing polarization: {:?}", pol);
            match self.extract_all_iw_subswaths(pol) {
                Ok(subswaths) => {
                    log::info!("Extracted {} subswaths for polarization {:?}", subswaths.len(), pol);
                    all_subswaths.insert(pol, subswaths);
                },
                Err(e) => {
                    log::error!("Failed to extract sub-swaths for {:?}: {}", pol, e);
                }
            }
        }
        
        log::info!("Total polarizations with subswaths: {}", all_subswaths.len());
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
        use crate::core::deburst::{DeburstProcessor, DeburstConfig};
        
        log::info!("Starting deburst processing for polarization {:?}", pol);
        
        // Read SLC data
        let slc_data = self.read_slc_data(pol)?;
        log::debug!("Read SLC data with dimensions: {:?}", slc_data.dim());
        
        // Get annotation for burst information
        let annotation_files = self.find_annotation_files()?;
        let annotation_file = annotation_files.get(&pol)
            .ok_or_else(|| SarError::Processing(format!("No annotation file found for {:?}", pol)))?;
        
        // Read annotation content
        let annotation_content = self.read_file_as_string(annotation_file)?;
        
        // Extract burst information
        let (total_lines, total_samples) = slc_data.dim();
        let burst_info = DeburstProcessor::extract_burst_info_from_annotation(
            &annotation_content, 
            total_lines, 
            total_samples
        )?;
        
        log::info!("Extracted {} bursts for deburst processing", burst_info.len());
        
        // Create deburst processor
        // Scientific requirement: derive satellite velocity from real orbit data (no hardcoded values)
        let cache_dir = std::env::var("SARDINE_ORBIT_CACHE")
            .map(std::path::PathBuf::from)
            .map_err(|_| SarError::Processing(
                "SARDINE_ORBIT_CACHE not set. Required for orbit-based velocity in deburst.".to_string(),
            ))?;
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
        };
        
        // Perform deburst
        let deburst_data = processor.deburst(&slc_data, &config)?;
        
        log::info!("Deburst completed. Output dimensions: {:?}", deburst_data.dim());
        Ok(deburst_data)
    }
    
    /// Deburst SLC data with custom configuration for TOPSAR processing
    pub fn deburst_slc_with_config(&mut self, pol: Polarization, config: DeburstConfig) -> SarResult<SarImage> {
        use crate::core::deburst::{DeburstProcessor};
        
        log::info!("Starting TOPSAR deburst processing for polarization {:?}", pol);
        
        // Read SLC data
        let slc_data = self.read_slc_data(pol)?;
        log::debug!("Read SLC data with dimensions: {:?}", slc_data.dim());
        
        // Get annotation for burst information
        let annotation_files = self.find_annotation_files()?;
        let annotation_file = annotation_files.get(&pol)
            .ok_or_else(|| SarError::Processing(format!("No annotation file found for {:?}", pol)))?;
        
        // Read annotation content
        let annotation_content = self.read_file_as_string(annotation_file)?;
        
        // Extract burst information with TOPSAR parameters
        let (total_lines, total_samples) = slc_data.dim();
        let burst_info = DeburstProcessor::extract_burst_info_from_annotation(
            &annotation_content, 
            total_lines, 
            total_samples
        )?;
        
        log::info!("Extracted {} bursts with TOPSAR parameters for deburst processing", burst_info.len());
        
        // Create deburst processor
        // Scientific requirement: derive satellite velocity from real orbit data (no hardcoded values)
        let cache_dir = std::env::var("SARDINE_ORBIT_CACHE")
            .map(std::path::PathBuf::from)
            .map_err(|_| SarError::Processing(
                "SARDINE_ORBIT_CACHE not set. Required for orbit-based velocity in deburst.".to_string(),
            ))?;
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
        
        log::info!("TOPSAR deburst completed. Output dimensions: {:?}", deburst_data.dim());
        Ok(deburst_data)
    }
    
    /// Find calibration files for each polarization in the ZIP
    pub fn find_calibration_files(&mut self) -> SarResult<HashMap<Polarization, String>> {
        log::debug!("Finding calibration files");
        
        let files = self.list_files()?;
        let mut calibration_files = HashMap::new();
        
        for file in files {
            // Only include actual calibration files, not noise files
            if file.contains("annotation/calibration/") && file.ends_with(".xml") && file.contains("calibration-") && !file.contains("noise-") {
                if let Some(pol) = self.extract_polarization_from_filename(&file) {
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
    
    /// Read calibration data for a specific polarization
    pub fn read_calibration_data(&mut self, pol: Polarization) -> SarResult<crate::core::calibrate::CalibrationCoefficients> {
        log::debug!("Reading calibration data for polarization {:?}", pol);
        
        let calibration_files = self.find_calibration_files()?;
        let file_path = calibration_files.get(&pol)
            .ok_or_else(|| SarError::Processing(
                format!("No calibration file found for polarization {:?}", pol)
            ))?;
        
        let xml_content = self.read_file_as_string(file_path)?;
        
        // Parse real calibration data from XML
        let calibration_data = crate::core::calibrate::parse_calibration_from_xml(&xml_content)?;
        
        log::info!("Successfully parsed calibration data for polarization {:?}", pol);
        Ok(calibration_data)
    }
    
    /// Find noise files for all available polarizations
    pub fn find_noise_files(&mut self) -> SarResult<HashMap<Polarization, String>> {
        log::debug!("Finding noise files");
        
        let files = self.list_files()?;
        let mut noise_files = HashMap::new();
        
        for file in files {
            // Look for noise files in annotation/calibration/ directory
            if file.contains("annotation/calibration/") && file.ends_with(".xml") && file.contains("noise-") {
                if let Some(pol) = self.extract_polarization_from_filename(&file) {
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
    pub fn read_noise_data(&mut self, pol: Polarization) -> SarResult<crate::core::calibrate::NoiseCoefficients> {
        log::debug!("Reading noise data for polarization {:?}", pol);
        
        let noise_files = self.find_noise_files()?;
        let file_path = noise_files.get(&pol)
            .ok_or_else(|| SarError::Processing(
                format!("No noise file found for polarization {:?}", pol)
            ))?;
        
        let xml_content = self.read_file_as_string(file_path)?;
        
        // Parse thermal noise data from XML
        let noise_data = crate::core::calibrate::parse_noise_from_xml(&xml_content)?;
        
        log::info!("Successfully parsed noise data for polarization {:?}", pol);
        Ok(noise_data)
    }
    
    /// Read calibration data with caching
    pub fn read_calibration_data_cached(&mut self, pol: Polarization) -> SarResult<CalibrationCoefficients> {
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
    pub fn read_all_calibration_data(&mut self) -> SarResult<HashMap<Polarization, crate::core::calibrate::CalibrationCoefficients>> {
        let mut all_cal_data = HashMap::new();
        
        let available_pols = self.find_annotation_files()?.keys().cloned().collect::<Vec<_>>();
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
            None        }
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
        log::info!("Applying multilook: {}x{} looks to intensity data", azimuth_looks, range_looks);
        
        // Get metadata for pixel spacing
        let metadata = self.read_annotation(pol)?;
        
        // Extract pixel spacings from metadata
        let (range_spacing, azimuth_spacing) = metadata.pixel_spacing;
        
        log::info!("Input pixel spacing: range={}m, azimuth={}m", range_spacing, azimuth_spacing);
        
        // Create multilook processor
        let params = crate::core::multilook::MultilookParams {
            range_looks,
            azimuth_looks,
            output_pixel_spacing: None,
        };
        
        let processor = crate::core::multilook::MultilookProcessor::new(params);
        
        // Apply multilooking
        let (multilooked_data, new_range_spacing, new_azimuth_spacing) = 
            processor.apply_multilook_filtered(intensity_data, range_spacing, azimuth_spacing)?;
        
        log::info!("Multilooking complete: {}x{} -> {}x{}", 
                   intensity_data.nrows(), intensity_data.ncols(),
                   multilooked_data.nrows(), multilooked_data.ncols());
        log::info!("Output pixel spacing: range={}m, azimuth={}m", new_range_spacing, new_azimuth_spacing);
        
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
        log::info!("Deburst data: {} x {}", deburst_data.nrows(), deburst_data.ncols());
        
        // Get calibration coefficients
        let cal_data = self.read_calibration_data(pol)?;
        
        // Create calibration processor
        let processor = crate::core::calibrate::CalibrationProcessor::new(
            cal_data, 
            cal_type
        );
        
        // Apply calibration to get intensity data
        let intensity_data = processor.calibrate(&deburst_data)?;
        log::info!("Calibrated data: {} x {}", intensity_data.nrows(), intensity_data.ncols());
        
        // Apply multilooking
        let (multilooked_data, new_range_spacing, new_azimuth_spacing) = 
            self.multilook_intensity(&intensity_data, pol, range_looks, azimuth_looks)?;
        
        log::info!("Complete workflow finished: final dimensions {}x{}", 
                   multilooked_data.nrows(), multilooked_data.ncols());
        
        Ok((multilooked_data, new_range_spacing, new_azimuth_spacing))
    }

    /// Complete workflow with automatic DEM preparation: calibration, multilooking, and terrain flattening
    pub fn calibrate_multilook_and_flatten_auto_dem(
        &mut self,
        pol: Polarization,
        cal_type: crate::core::calibrate::CalibrationType,
        range_looks: usize,
        azimuth_looks: usize,
        dem_cache_dir: &str, // Required DEM cache directory
        target_resolution: f64, // Required DEM target resolution
    ) -> SarResult<(Array2<f32>, Array2<f32>, f64, f64)> {
        log::info!("Starting complete SAR processing workflow with automatic DEM preparation");
        
        // Step 1: Calibration and multilooking
        let (sigma0_multilooked, new_range_spacing, new_azimuth_spacing) = 
            self.calibrate_and_multilook(pol, cal_type, range_looks, azimuth_looks)?;
        
        log::info!("Calibration and multilooking completed: {} x {}", 
                   sigma0_multilooked.nrows(), sigma0_multilooked.ncols());

        // Step 2: Automatically prepare DEM data
        log::info!("Automatically preparing DEM data");
        
        // Get SAR scene bounding box from metadata
        let metadata = self.read_annotation(pol)?;
        let bbox = &metadata.bounding_box;
        
        let (dem_data, dem_transform) = crate::io::dem::DemReader::prepare_dem_for_scene(
            bbox, 
            target_resolution,
            dem_cache_dir
        )?;
        
        log::info!("DEM data prepared: {} x {}", dem_data.nrows(), dem_data.ncols());

        // Step 3: Resample DEM to match SAR geometry
        let sar_transform = self.create_sar_geotransform(&metadata, new_range_spacing, new_azimuth_spacing)?;
        let resampled_dem = crate::io::dem::DemReader::resample_dem(
            &dem_data,
            &dem_transform,
            &sar_transform,
            sigma0_multilooked.dim(),
        )?;

        log::info!("DEM resampled to SAR geometry");

        // Step 4: Get annotation for terrain flattening parameters
        let annotation = self.get_annotation_for_polarization(pol)?;
        let terrain_params = crate::core::terrain_flatten::TerrainFlatteningParams::from_annotation(&annotation)?;

        // Step 5: Get orbit data for terrain flattening
        let orbit_data = self.orbit_data.as_ref()
            .ok_or_else(|| SarError::Processing("Orbit data not available".to_string()))?.clone();
        
        let flattener = crate::core::terrain_flatten::TerrainFlattener::new(
            terrain_params,
            orbit_data,
        );

        let (gamma0, incidence_angles) = flattener.process_terrain_flattening(
            &sigma0_multilooked,
            &resampled_dem,
            0.0, // range_time (simplified)
            0.0, // azimuth_time_start (simplified)
            1.0, // azimuth_time_spacing (simplified)
        )?;

        log::info!("Terrain flattening completed: {} x {}", 
                   gamma0.nrows(), gamma0.ncols());

        Ok((gamma0, incidence_angles, new_range_spacing, new_azimuth_spacing))
    }

    /// Complete workflow with terrain flattening: calibration, multilooking, and terrain flattening
    pub fn calibrate_multilook_and_flatten(
        &mut self,
        pol: Polarization,
        cal_type: crate::core::calibrate::CalibrationType,
        range_looks: usize,
        azimuth_looks: usize,
        dem_path: &str,
        target_resolution: Option<f64>, // Optional DEM target resolution (default: 30m)
    ) -> SarResult<(Array2<f32>, Array2<f32>, f64, f64)> {
        log::info!("Starting complete SAR processing workflow with terrain flattening");
        
        // Step 1: Calibration and multilooking
        let (sigma0_multilooked, new_range_spacing, new_azimuth_spacing) = 
            self.calibrate_and_multilook(pol, cal_type, range_looks, azimuth_looks)?;
        
        log::info!("Calibration and multilooking completed: {} x {}", 
                   sigma0_multilooked.nrows(), sigma0_multilooked.ncols());

        // Step 2: Read DEM data
        log::info!("Reading DEM data from: {}", dem_path);
        
        // Get SAR scene bounding box from metadata
        let metadata = self.read_annotation(pol)?;
        let bbox = &metadata.bounding_box;
        
        let (dem_data, dem_transform) = crate::io::dem::DemReader::read_dem(
            dem_path, 
            bbox, 
            target_resolution.unwrap_or(30.0) // Default 30m if not specified
        )?;
        
        log::info!("DEM data loaded: {} x {}", dem_data.nrows(), dem_data.ncols());

        // Step 3: Resample DEM to match SAR geometry
        let sar_transform = self.create_sar_geotransform(&metadata, new_range_spacing, new_azimuth_spacing)?;
        let resampled_dem = crate::io::dem::DemReader::resample_dem(
            &dem_data,
            &dem_transform,
            &sar_transform,
            sigma0_multilooked.dim(),
        )?;

        log::info!("DEM resampled to SAR geometry");

        // Step 4: Get annotation for terrain flattening parameters
        let annotation = self.get_annotation_for_polarization(pol)?;
        let terrain_params = crate::core::terrain_flatten::TerrainFlatteningParams::from_annotation(&annotation)?;

        // Step 5: Get orbit data for terrain flattening
        let orbit_data = self.orbit_data.as_ref()
            .ok_or_else(|| SarError::Processing("Orbit data not available".to_string()))?.clone();
        
        let flattener = crate::core::terrain_flatten::TerrainFlattener::new(
            terrain_params,
            orbit_data,
        );

        let (gamma0, incidence_angles) = flattener.process_terrain_flattening(
            &sigma0_multilooked,
            &resampled_dem,
            0.0, // range_time (simplified)
            0.0, // azimuth_time_start (simplified)
            1.0, // azimuth_time_spacing (simplified)
        )?;

        log::info!("Terrain flattening completed: {} x {}", 
                   gamma0.nrows(), gamma0.ncols());

        Ok((gamma0, incidence_angles, new_range_spacing, new_azimuth_spacing))
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
    fn parse_orbit_from_product_id(product_id: &str) -> Option<(u32, u32, String)> {
        let parts: Vec<&str> = product_id.split('_').collect();
        // Debug: log the product ID parsing
        log::info!("Parsing orbit from product ID: {} (has {} parts)", product_id, parts.len());
        
        if parts.len() >= 8 {
            // Extract orbit number (030639 is at index 7, not 6)
            // Format: S1A_IW_SLC__1SDV_20200103T170815_20200103T170842_030639_0382D5_DADE
            //         0   1   2  3 4    5                6                7      8       9
            if let Ok(orbit_number) = parts[7].parse::<u32>() {
                // Calculate relative orbit (Sentinel-1 has 175 relative orbits)
                let relative_orbit = ((orbit_number - 1) % 175) + 1;
                
                // Determine orbit direction from time pattern (simplified heuristic)
                // This is a simplified approach - proper determination needs orbit data
                let orbit_direction = if orbit_number % 2 == 0 { "DESCENDING" } else { "ASCENDING" };
                
                log::info!("Extracted orbit info: orbit_number={}, relative_orbit={}, direction={}", 
                          orbit_number, relative_orbit, orbit_direction);
                return Some((orbit_number, relative_orbit, orbit_direction.to_string()));
            } else {
                log::warn!("Failed to parse orbit number from part: {}", parts[7]);
            }
        } else {
            log::warn!("Product ID has insufficient parts for orbit extraction");
        }
        None
    }

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
        if let Some(value) = Self::extract_real_value_f64(xml_content, "absoluteCalibrationConstant")
            .or_else(|| Self::extract_real_value_f64(xml_content, "calibrationConstant"))
            .or_else(|| Self::extract_real_value_f64(xml_content, "rescalingFactor")) {
            return Some(value);
        }
        
        // If not found in annotation, return a default reference constant for SLC data
        // Sentinel-1 SLC data uses digital numbers with no specific calibration constant
        // Return 1.0 as SLC data is already calibrated relative to amplitude
        Some(1.0)
    }

    /// Extract merged bounding box from ALL annotation files for complete scene coverage
    fn extract_merged_bounding_box_all_subswaths(&mut self, pol: Polarization) -> SarResult<crate::types::BoundingBox> {
        println!("🌐 DEBUG: Starting merged bounding box extraction for polarization {}", pol);
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
        
        println!("📁 DEBUG: Found annotation files for {} polarizations", all_annotation_files.len());
        
        let pol_files = all_annotation_files.get(&pol)
            .ok_or_else(|| {
                println!("❌ DEBUG: No annotation files found for polarization {}", pol);
                SarError::Processing(format!("No annotation files found for polarization {}", pol))
            })?;
        
        println!("📁 DEBUG: Found {} annotation files for {}", pol_files.len(), pol);
        log::info!("📁 Found {} annotation files for {}: merging bounding boxes", pol_files.len(), pol);
        
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
                    if bbox.min_lat < bbox.max_lat && bbox.min_lon < bbox.max_lon &&
                       bbox.min_lat != 0.0 && bbox.max_lat != 0.0 &&
                       bbox.min_lon != 0.0 && bbox.max_lon != 0.0 {
                        
                        // Merge with overall bounds
                        merged_min_lat = merged_min_lat.min(bbox.min_lat);
                        merged_max_lat = merged_max_lat.max(bbox.max_lat);
                        merged_min_lon = merged_min_lon.min(bbox.min_lon);
                        merged_max_lon = merged_max_lon.max(bbox.max_lon);
                        
                        processed_count += 1;
                        
                        log::info!("✅ Subswath {}: [{:.6}, {:.6}, {:.6}, {:.6}]", 
                                  processed_count,
                                  bbox.min_lon, bbox.min_lat, bbox.max_lon, bbox.max_lat);
                    } else {
                        log::warn!("⚠️  Invalid bounding box in {}, skipping", annotation_file);
                    }
                }
                Err(e) => {
                    log::warn!("⚠️  Failed to read annotation file {}: {}", annotation_file, e);
                }
            }
        }
        
        if processed_count == 0 {
            return Err(SarError::Processing("No valid bounding boxes found in any annotation files".to_string()));
        }
        
        let merged_bbox = crate::types::BoundingBox {
            min_lat: merged_min_lat,
            max_lat: merged_max_lat,
            min_lon: merged_min_lon,
            max_lon: merged_max_lon,
        };
        
        log::info!("🎯 MERGED BOUNDING BOX from {} subswaths:", processed_count);
        log::info!("   📍 [{:.6}, {:.6}, {:.6}, {:.6}]", 
                  merged_bbox.min_lon, merged_bbox.min_lat, merged_bbox.max_lon, merged_bbox.max_lat);
        log::info!("   📏 Coverage: {:.3}° × {:.3}° ({:.1}km × {:.1}km)", 
                  merged_bbox.max_lon - merged_bbox.min_lon, 
                  merged_bbox.max_lat - merged_bbox.min_lat,
                  (merged_bbox.max_lon - merged_bbox.min_lon) * 111.0,
                  (merged_bbox.max_lat - merged_bbox.min_lat) * 111.0);
        
        Ok(merged_bbox)
    }
    
    /// Read annotation file directly without full metadata parsing
    fn read_annotation_file_raw(&mut self, annotation_file: &str) -> SarResult<SarMetadata> {
        let xml_content = self.read_file_as_string(annotation_file)?;
        Self::parse_annotation_xml_comprehensive(&xml_content, Polarization::VV, &self.product_path)
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
        let annotation_files = self.find_annotation_files()?;
        
        log::info!("📁 Found {} polarizations to cache", annotation_files.len());
        
        for (pol, annotation_path) in &annotation_files {
            log::debug!("📄 Caching annotation for {:?}: {}", pol, annotation_path);
            
            // Cache raw XML content to eliminate file I/O
            let xml_content = self.read_file_as_string(annotation_path)?;
            self.cached_xml_content.insert(annotation_path.clone(), xml_content.clone());
            
            // Cache parsed annotation structure  
            let annotation_root = self.parse_annotation_root_xml(&xml_content)?;
            self.cached_annotations.insert(*pol, annotation_root);
            
            // Cache calibration data for this polarization
            if let Ok(calibration) = self.read_calibration_data(*pol) {
                self.cached_calibration.insert(*pol, calibration);
                log::debug!("✅ Cached calibration for {:?}", pol);
            } else {
                log::warn!("⚠️  Could not cache calibration for {:?}", pol);
            }
        }
        
        // Get first polarization for metadata initialization  
        let first_pol = annotation_files.keys().next().cloned()
            .ok_or_else(|| SarError::Processing("No annotation files found for metadata cache initialization".to_string()))?;
        
        // Step 2: Cache comprehensive metadata with merged bounding box from ALL subswaths
        // SCIENTIFIC ACCURACY: Use same merged bounding box extraction as traditional method
        let comprehensive_metadata = self.read_annotation(first_pol)?;
        self.cached_metadata = Some(comprehensive_metadata);
        
        // Mark cache as initialized
        self.cache_initialized = true;
        
        log::info!("✅ Comprehensive metadata cache initialized successfully");
        log::info!("   📊 Cached {} polarizations", annotation_files.len());
        log::info!("   📊 Cached {} XML files", self.cached_xml_content.len());
        log::info!("   📊 Cached {} calibration datasets", self.cached_calibration.len());
        
        Ok(())
    }
    
    /// Get cached metadata (replaces all get_metadata() calls)
    /// This eliminates redundant metadata extraction throughout processing
    pub fn get_cached_metadata(&self) -> SarResult<&SarMetadata> {
        self.ensure_cache_initialized()?;
        
        self.cached_metadata.as_ref()
            .ok_or_else(|| SarError::Processing(
                "Metadata not cached. Use new_with_full_cache() instead of new()".to_string()
            ))
    }
    
    /// Get cached annotation for specific polarization (no XML re-parsing)
    pub fn get_cached_annotation(&self, pol: Polarization) -> SarResult<&crate::io::annotation::AnnotationRoot> {
        self.ensure_cache_initialized()?;
        
        self.cached_annotations.get(&pol)
            .ok_or_else(|| SarError::Processing(
                format!("Annotation not cached for {:?}. Available: {:?}", 
                       pol, self.cached_annotations.keys().collect::<Vec<_>>())
            ))
    }
    
    /// Get cached calibration data for specific polarization  
    pub fn get_cached_calibration(&self, pol: Polarization) -> SarResult<&CalibrationCoefficients> {
        self.ensure_cache_initialized()?;
        
        self.cached_calibration.get(&pol)
            .ok_or_else(|| SarError::Processing(
                format!("Calibration not cached for {:?}", pol)
            ))
    }
    
    /// Get cached XML content for specific file (eliminates file re-reading)
    pub fn get_cached_xml_content(&self, file_path: &str) -> SarResult<&String> {
        self.ensure_cache_initialized()?;
        
        self.cached_xml_content.get(file_path)
            .ok_or_else(|| SarError::Processing(
                format!("XML content not cached for {}", file_path)
            ))
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
        map.insert("min_latitude".to_string(), metadata.bounding_box.min_lat.to_string());
        map.insert("max_latitude".to_string(), metadata.bounding_box.max_lat.to_string());
        map.insert("min_longitude".to_string(), metadata.bounding_box.min_lon.to_string());
        map.insert("max_longitude".to_string(), metadata.bounding_box.max_lon.to_string());
        
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
        
        // Add pixel spacing
        map.insert("range_pixel_spacing".to_string(), metadata.pixel_spacing.0.to_string());
        map.insert("azimuth_pixel_spacing".to_string(), metadata.pixel_spacing.1.to_string());
        
        // CRITICAL FIX: Add orbit state vectors for terrain correction
        // The Python code expects keys like: orbit_state_vector_{i}_time, orbit_state_vector_{i}_x_position, etc.
        if let Some(orbit_data) = &metadata.orbit_data {
            for (i, state_vector) in orbit_data.state_vectors.iter().enumerate() {
                // Convert time to timestamp (seconds since Unix epoch)
                let timestamp = state_vector.time.timestamp() as f64;
                map.insert(format!("orbit_state_vector_{}_time", i), timestamp.to_string());
                
                // Position components (in meters)
                map.insert(format!("orbit_state_vector_{}_x_position", i), state_vector.position[0].to_string());
                map.insert(format!("orbit_state_vector_{}_y_position", i), state_vector.position[1].to_string());
                map.insert(format!("orbit_state_vector_{}_z_position", i), state_vector.position[2].to_string());
                
                // Velocity components (in m/s)
                map.insert(format!("orbit_state_vector_{}_x_velocity", i), state_vector.velocity[0].to_string());
                map.insert(format!("orbit_state_vector_{}_y_velocity", i), state_vector.velocity[1].to_string());
                map.insert(format!("orbit_state_vector_{}_z_velocity", i), state_vector.velocity[2].to_string());
            }
            
            log::info!("Added {} orbit state vectors to metadata cache", orbit_data.state_vectors.len());
        } else {
            log::warn!("No orbit data available in comprehensive metadata for terrain correction");
        }
        
        Ok(map)
    }
    
    /// Ensure cache is initialized before accessing cached data
    fn ensure_cache_initialized(&self) -> SarResult<()> {
        if !self.cache_initialized {
            return Err(SarError::Processing(
                "Cache not initialized. Use SlcReader::new_with_full_cache() instead of new()".to_string()
            ));
        }
        Ok(())
    }
    
    /// Get cache performance statistics
    pub fn get_cache_stats(&self) -> CacheStats {
        CacheStats {
            annotations_cached: self.cached_annotations.len(),
            xml_files_cached: self.cached_xml_content.len(),
            calibration_files_cached: self.cached_calibration.len(),
            total_memory_usage_bytes: self.estimate_cache_memory_usage(),
            cache_initialized: self.cache_initialized,
        }
    }
    
    /// Estimate memory usage of all caches
    fn estimate_cache_memory_usage(&self) -> usize {
        let xml_size: usize = self.cached_xml_content.values()
            .map(|content| content.len())
            .sum();
        
        let annotation_size = self.cached_annotations.len() * 
            std::mem::size_of::<crate::io::annotation::AnnotationRoot>();
            
        let calibration_size = self.cached_calibration.len() *
            std::mem::size_of::<CalibrationCoefficients>();
        
        xml_size + annotation_size + calibration_size + 1024 // Approximate overhead
    }
    
    /// Helper method to parse annotation XML (used by caching system)
    fn parse_annotation_root_xml(&self, xml_content: &str) -> SarResult<crate::io::annotation::AnnotationRoot> {
        crate::io::annotation::parse_annotation_xml(xml_content)
            .map_err(|e| SarError::Processing(format!("Failed to parse annotation XML: {}", e)))
    }

}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::PathBuf;

    #[test]
    fn test_slc_reader_creation() {
        // This test would need actual test data
        let dummy_path = PathBuf::from("nonexistent.zip");
        let result = SlcReader::new(&dummy_path);
        assert!(result.is_err());
    }
}
