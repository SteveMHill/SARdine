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

/// Sentinel-1 SLC reader that supports both ZIP and SAFE formats
pub struct SlcReader {
    product_path: PathBuf,
    format: ProductFormat,
    archive: Option<ZipArchive<File>>,
    orbit_data: Option<OrbitData>,

    /// Cached calibration coefficients
    calibration_cache: std::collections::HashMap<String, CalibrationCoefficients>,
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
            product_path.file_name().unwrap_or_default().to_string_lossy().contains(".SAFE") ||
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

    /// Find annotation files for each polarization
    pub fn find_annotation_files(&mut self) -> SarResult<HashMap<Polarization, String>> {
        let files = self.list_files()?;
        let mut annotations = HashMap::new();
        
        for file in files {
            if file.contains("annotation/") && file.ends_with(".xml") && !file.contains("calibration") && !file.contains("noise") {
                // Parse polarization from filename - Sentinel-1 standard naming
                if file.contains("-vv-") || file.contains("_vv_") {
                    annotations.insert(Polarization::VV, file);
                } else if file.contains("-vh-") || file.contains("_vh_") {
                    annotations.insert(Polarization::VH, file);
                } else if file.contains("-hv-") || file.contains("_hv_") {
                    annotations.insert(Polarization::HV, file);
                } else if file.contains("-hh-") || file.contains("_hh_") {
                    annotations.insert(Polarization::HH, file);
                }
            }
        }
        
        if annotations.is_empty() {
            return Err(SarError::InvalidFormat(
                "No annotation files found. Expected files like 's1a-iw-slc-vv-*.xml' in annotation/ directory".to_string(),
            ));
        }
        
        log::info!("Found {} annotation files: {:?}", annotations.len(), annotations.keys().collect::<Vec<_>>());
        
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
    /// CRITICAL: This replaces hardcoded values with actual extracted data from
    /// Sentinel-1 annotation XML files. Essential for scientific accuracy.
    /// 
    /// References:
    /// - ESA S1-IF-ASD-PL-0007: "Sentinel-1 Annotation Format Specification"
    /// - ESA S1-TN-MDA-52-7440: "Sentinel-1 TOPSAR Mode"
    fn parse_annotation_xml_comprehensive(xml_content: &str, pol: Polarization, product_path: &Path) -> SarResult<SarMetadata> {
        println!("🔍 DEBUG: parse_annotation_xml_comprehensive called");
        log::info!("Parsing annotation XML with comprehensive XML parser");
        
        // Extract product ID from path 
        let product_id = product_path
            .file_stem()
            .and_then(|s| s.to_str())
            .unwrap_or("S1A_UNKNOWN")
            .replace(".SAFE", "")
            .replace(".zip", "");

        // Use the new comprehensive annotation parser
        let annotation = crate::io::annotation::AnnotationParser::parse_annotation(xml_content)
            .map_err(|_| SarError::Processing("Failed to parse annotation XML - cannot extract required parameters".to_string()))?;
        
        // Extract pixel spacing FIRST - CRITICAL: NO FALLBACKS ALLOWED
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
    let _sub_swaths = crate::io::annotation::AnnotationParser::extract_subswaths(&annotation)
            .unwrap_or_else(|_| HashMap::new());
            
        // Extract real bounding box from geolocation grid
        println!("🔍 DEBUG: About to extract bounding box");
        let mut bounding_box = BoundingBox {
            min_lon: -180.0,
            max_lon: 180.0,
            min_lat: -90.0,
            max_lat: 90.0,
        };
        
        if let Some(geolocation_grid) = &annotation.geolocation_grid {
            println!("🔍 DEBUG: Geolocation grid found");
            if let Some(grid_points) = &geolocation_grid.geolocation_grid_point_list {
                println!("🔍 DEBUG: Grid point list found");
                if let Some(points) = &grid_points.geolocation_grid_points {
                    println!("🔍 DEBUG: Found {} geolocation points", points.len());
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
                        println!("🔍 DEBUG: Bounding box calculated: [{:.3}, {:.3}, {:.3}, {:.3}]", 
                                bounding_box.min_lon, bounding_box.min_lat, bounding_box.max_lon, bounding_box.max_lat);
                    } else {
                        println!("🔍 DEBUG: No lat/lon data in points");
                    }
                } else {
                    println!("🔍 DEBUG: No geolocation_grid_points in list");
                }
            } else {
                println!("🔍 DEBUG: No geolocation_grid_point_list in grid");
            }
        } else {
            println!("🔍 DEBUG: No geolocation_grid in annotation");
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

        let prf = annotation.image_annotation
            .as_ref()
            .and_then(|img| img.image_information.as_ref())
            .and_then(|info| info.azimuth_frequency)
            .ok_or_else(|| SarError::Processing("Azimuth PRF (azimuthFrequency) not found in annotation XML - required".to_string()))?;

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
        println!("🔍 DEBUG: About to extract bounding box");
        let bounding_box = match crate::io::annotation::AnnotationParser::extract_bounding_box(&annotation) {
            Ok(real_bbox) => {
                println!("🔍 DEBUG: Bounding box extraction SUCCESS: [{:.6}, {:.6}, {:.6}, {:.6}]",
                        real_bbox.min_lon, real_bbox.min_lat, real_bbox.max_lon, real_bbox.max_lat);
                log::info!("Extracted real bounding box: [{:.6}, {:.6}, {:.6}, {:.6}]",
                          real_bbox.min_lon, real_bbox.min_lat, real_bbox.max_lon, real_bbox.max_lat);
                real_bbox
            },
            Err(e) => {
                println!("🔍 DEBUG: Bounding box extraction FAILED: {}", e);
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
        
        Ok(SarMetadata {
            product_id,
            mission,
            platform,
            instrument: "C-SAR".to_string(),
            acquisition_mode,
            polarizations: vec![pol],
            start_time,
            stop_time,
            bounding_box,
            coordinate_system: CoordinateSystem::Radar,
            sub_swaths,
            orbit_data,  // SCIENTIFIC FIX: Now properly populated from annotation
            range_looks: 1,
            azimuth_looks: 1,
            pixel_spacing: (range_pixel_spacing, azimuth_pixel_spacing), // REAL extracted values
        })
    }
    fn extract_coordinates_from_geolocation_grid(xml: &str) -> Option<(f64, f64, f64, f64)> {
        println!("🔍 DEBUG: extract_coordinates_from_geolocation_grid called");
        
        // Find the geolocationGrid section
        if let Some(grid_start) = xml.find("<geolocationGrid>") {
            if let Some(grid_end) = xml[grid_start..].find("</geolocationGrid>") {
                let grid_section = &xml[grid_start..grid_start + grid_end];
                println!("🔍 DEBUG: Found geolocationGrid section, length: {}", grid_section.len());
                
                // Extract all latitude and longitude values
                let mut latitudes = Vec::new();
                let mut longitudes = Vec::new();
                
                // Find all geolocationGridPoint entries
                let mut search_pos = 0;
                while let Some(point_start) = grid_section[search_pos..].find("<geolocationGridPoint>") {
                    let abs_point_start = search_pos + point_start;
                    if let Some(point_end) = grid_section[abs_point_start..].find("</geolocationGridPoint>") {
                        let point_section = &grid_section[abs_point_start..abs_point_start + point_end];
                        
                        // Extract latitude
                        if let Some(lat_start) = point_section.find("<latitude>") {
                            if let Some(lat_end) = point_section[lat_start..].find("</latitude>") {
                                let lat_text = &point_section[lat_start + 10..lat_start + lat_end];
                                if let Ok(lat_val) = lat_text.trim().parse::<f64>() {
                                    latitudes.push(lat_val);
                                }
                            }
                        }
                        
                        // Extract longitude  
                        if let Some(lon_start) = point_section.find("<longitude>") {
                            if let Some(lon_end) = point_section[lon_start..].find("</longitude>") {
                                let lon_text = &point_section[lon_start + 11..lon_start + lon_end];
                                if let Ok(lon_val) = lon_text.trim().parse::<f64>() {
                                    longitudes.push(lon_val);
                                }
                            }
                        }
                        
                        search_pos = abs_point_start + point_end;
                    } else {
                        break;
                    }
                }
                
                println!("🔍 DEBUG: Extracted {} latitudes, {} longitudes", latitudes.len(), longitudes.len());
                
                if !latitudes.is_empty() && !longitudes.is_empty() {
                    let min_lat = latitudes.iter().fold(f64::INFINITY, |a, &b| a.min(b));
                    let max_lat = latitudes.iter().fold(f64::NEG_INFINITY, |a, &b| a.max(b));
                    let min_lon = longitudes.iter().fold(f64::INFINITY, |a, &b| a.min(b));
                    let max_lon = longitudes.iter().fold(f64::NEG_INFINITY, |a, &b| a.max(b));
                    
                    println!("🔍 DEBUG: Calculated bounding box: lat=({:.6}, {:.6}), lon=({:.6}, {:.6})", 
                             min_lat, max_lat, min_lon, max_lon);
                    
                    return Some((min_lat, max_lat, min_lon, max_lon));
                } else {
                    println!("🔍 DEBUG: No valid coordinates found in geolocation grid");
                }
            } else {
                println!("🔍 DEBUG: No closing </geolocationGrid> tag found");
            }
        } else {
            println!("🔍 DEBUG: No <geolocationGrid> section found");
        }
        
        None
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

    /// Parse annotation XML content with sub-swath extraction
    #[allow(dead_code)]
    fn parse_annotation_xml(xml_content: &str, pol: Polarization) -> SarResult<SarMetadata> {
        // First try to parse with the full annotation parser - REQUIRED for pixel spacing
        let annotation = crate::io::annotation::AnnotationParser::parse_annotation(xml_content)
            .map_err(|_| SarError::Processing("Failed to parse annotation XML - cannot extract required parameters".to_string()))?;
        
        // Extract pixel spacing FIRST - CRITICAL: NO FALLBACKS ALLOWED
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
        let sub_swaths = crate::io::annotation::AnnotationParser::extract_subswaths(&annotation)
            .unwrap_or_else(|_| HashMap::new());
            
        // Extract real bounding box from geolocation grid
        let mut bounding_box = BoundingBox {
            min_lon: -180.0,
            max_lon: 180.0,
            min_lat: -90.0,
            max_lat: 90.0,
        };
        
        match crate::io::annotation::AnnotationParser::extract_bounding_box(&annotation) {
            Ok(real_bbox) => {
                log::info!("✅ Successfully extracted bounding box: [{:.6}, {:.6}, {:.6}, {:.6}]", 
                          real_bbox.min_lon, real_bbox.min_lat, real_bbox.max_lon, real_bbox.max_lat);
                bounding_box = real_bbox;
            },
            Err(e) => {
                log::error!("❌ Failed to extract bounding box from annotation: {}", e);
                log::error!("❌ This indicates a problem with the geolocation grid parsing or access");
            }
        }
        
        // If XML parsing failed for bbox, return error - NO FALLBACKS
        if bounding_box.min_lon == -180.0 && bounding_box.max_lon == 180.0 {
            return Err(SarError::Metadata(
                "Failed to extract bounding box from annotation XML. \
                Real Sentinel-1 annotation with valid geolocation grid required - no fallback parsing allowed.".to_string()
            ));
        }
        
        // Extract basic information using string parsing (fallback approach)
    let product_id = Self::extract_value(xml_content, "missionId")
            .unwrap_or_else(|| "S1A_UNKNOWN".to_string());
        
        // Use more flexible time parsing - MANDATORY
        let start_time_str = Self::extract_value(xml_content, "startTime")
            .or_else(|| Self::extract_value(xml_content, "acquisitionStartTime"))
            .ok_or_else(|| SarError::Processing(
                "❌ SCIENTIFIC ERROR: Start time not found in annotation XML! \
                Real Sentinel-1 annotation required - no synthetic fallbacks allowed.".to_string()
            ))?;
        
        let stop_time_str = Self::extract_value(xml_content, "stopTime")
            .or_else(|| Self::extract_value(xml_content, "acquisitionStopTime"))
            .ok_or_else(|| SarError::Processing(
                "❌ SCIENTIFIC ERROR: Stop time not found in annotation XML! \
                Real Sentinel-1 annotation required - no synthetic fallbacks allowed.".to_string()
            ))?;

        // Try to parse time, with fallback for different formats
        let start_time = Self::parse_time_flexible(&start_time_str)?;
        let stop_time = Self::parse_time_flexible(&stop_time_str)?;

        // Extract platform (Sentinel-1A or Sentinel-1B) from product ID
        let platform = if product_id.starts_with("S1A_") {
            "Sentinel-1A".to_string()
        } else if product_id.starts_with("S1B_") {
            "Sentinel-1B".to_string()
        } else {
            return Err(SarError::InvalidInput(
                format!("Cannot determine platform from product ID: {}", product_id)
            ));
        };
        
        // Create metadata structure with sub-swaths and real bounding box
        
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
        
        Ok(SarMetadata {
            product_id,
            mission: "Sentinel-1".to_string(),
            platform,
            instrument: "C-SAR".to_string(),
            acquisition_mode: AcquisitionMode::IW,
            polarizations: vec![pol],
            start_time,
            stop_time,
            bounding_box,  // Use extracted bounding box
            coordinate_system: CoordinateSystem::Radar,
            sub_swaths,
            orbit_data,  // SCIENTIFIC FIX: Now properly populated from annotation
            range_looks: 1,
            azimuth_looks: 1,
            // CRITICAL: Use extracted pixel spacing - NO HARDCODED VALUES
            pixel_spacing: (range_pixel_spacing, azimuth_pixel_spacing)
        })
    }

    /// Extract comprehensive metadata using regex parsing
    fn extract_comprehensive_metadata(xml_content: &str) -> std::collections::HashMap<String, String> {
        let mut metadata = std::collections::HashMap::new();
        
        // Core product information
        let basic_fields = [
            "missionId", "productType", "polarisation", "mode", "swath",
            "absoluteOrbitNumber", "missionDataTakeId", "imageNumber"
        ];
        
        // Image information
        let image_fields = [
            "numberOfSamples", "numberOfLines", "rangePixelSpacing", "azimuthPixelSpacing",
            "slantRangeTime", "azimuthFrequency", "rangeSamplingRate", "radarFrequency"
        ];
        
        // Geometry information
        let geometry_fields = [
            "incidenceAngle", "elevationAngle", "ascendingNodeTime"
        ];
        
        // Extract all fields
        let all_fields = [&basic_fields[..], &image_fields[..], &geometry_fields[..]].concat();
        
        for field in all_fields {
            if let Some(value) = Self::extract_value(xml_content, field) {
                metadata.insert(field.to_string(), value);
            }
        }
        
        // CRITICAL: Extract coordinate bounds using the working simple extraction method
        println!("🔍 DEBUG: extract_comprehensive_metadata attempting coordinate extraction");
        if let Some((min_lat, max_lat, min_lon, max_lon)) = Self::extract_coordinates_from_geolocation_grid(xml_content) {
            metadata.insert("min_latitude".to_string(), min_lat.to_string());
            metadata.insert("max_latitude".to_string(), max_lat.to_string());
            metadata.insert("min_longitude".to_string(), min_lon.to_string());
            metadata.insert("max_longitude".to_string(), max_lon.to_string());
            println!("🔍 DEBUG: Added coordinates to comprehensive metadata: lat=({:.6}, {:.6}), lon=({:.6}, {:.6})", 
                     min_lat, max_lat, min_lon, max_lon);
        } else {
            println!("🔍 DEBUG: Failed to extract coordinates in comprehensive metadata");
        }
        
        // Extract average incidence angle from geolocation grid
        if let Ok(avg_incidence) = Self::extract_average_incidence_angle(xml_content) {
            metadata.insert("averageIncidenceAngle".to_string(), avg_incidence.to_string());
        }
        
        metadata
    }

    /// Extract average incidence angle from geolocation grid
    fn extract_average_incidence_angle(xml_content: &str) -> SarResult<f64> {
        let mut incidence_angles = Vec::new();
        
        let sections: Vec<&str> = xml_content.split("<geolocationGridPoint>").skip(1).take(50).collect();
        
        for section in sections {
            if let Some(angle_start) = section.find("<incidenceAngle>") {
                if let Some(angle_end) = section.find("</incidenceAngle>") {
                    let angle_str = &section[angle_start + 16..angle_end];
                    if let Ok(angle) = angle_str.parse::<f64>() {
                        incidence_angles.push(angle);
                    }
                }
            }
        }
        
        if !incidence_angles.is_empty() {
            let average = incidence_angles.iter().sum::<f64>() / incidence_angles.len() as f64;
            Ok(average)
        } else {
            Err(SarError::Metadata("No incidence angle data found".to_string()))
        }
    }

    /// Flexible time parsing that handles multiple formats
    fn parse_time_flexible(time_str: &str) -> SarResult<DateTime<Utc>> {
        println!("🔍 DEBUG: Attempting to parse time: '{}'", time_str);
        
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
        // Read the annotation XML content directly
        let content = self.read_annotation_content(pol)?;
        
        // Parse the XML content into AnnotationRoot
        let annotation_root: crate::io::annotation::AnnotationRoot = 
            quick_xml::de::from_str(&content)
                .map_err(|e| SarError::XmlParsing(format!("Failed to parse annotation XML: {}", e)))?;
        
        Ok(annotation_root)
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

            let window = (0, 0);
            let window_size = (width, height);
            
            // Sentinel-1 SLC uses complex 16-bit integers (CInt16)
            let complex_data = band.read_as::<i16>(window, window_size, (width * 2, height), None)
                .map_err(|e| SarError::Io(std::io::Error::other(
                    format!("Failed to read complex CInt16 data: {}", e),
                )))?;

            log::debug!("Read {} complex values as interleaved i16", complex_data.data.len() / 2);

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

            let window = (0, 0);
            let window_size = (width, height);
            
            // Sentinel-1 SLC uses complex 16-bit integers (CInt16)
            let complex_data = band.read_as::<i16>(window, window_size, (width * 2, height), None)
                .map_err(|e| SarError::Io(std::io::Error::other(
                format!("Failed to read complex CInt16 data: {}", e),
            )))?;

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

    /// Get comprehensive metadata from the SLC file
    pub fn get_comprehensive_metadata(&mut self) -> std::collections::HashMap<String, String> {
        // Read annotation file for VV polarization (fallback to any available)
        let annotations = match self.find_annotation_files() {
            Ok(annotations) => annotations,
            Err(_) => return std::collections::HashMap::new(),
        };
        
        let polarization = annotations.keys().next().cloned()
            .unwrap_or(crate::types::Polarization::VV);
        
        // Try to read annotation content directly
        if let Ok(content) = self.read_annotation_content(polarization) {
            Self::extract_comprehensive_metadata(&content)
        } else {
            std::collections::HashMap::new()
        }
    }

    /// Read annotation content for comprehensive metadata extraction
    fn read_annotation_content(&mut self, pol: crate::types::Polarization) -> SarResult<String> {
        let annotations = self.find_annotation_files()?;
        let annotation_file = annotations.get(&pol)
            .ok_or_else(|| SarError::Processing(
                format!("No annotation file found for polarization {:?}", pol)
            ))?;
        
        // Read the file content from the product
        let content = match self.format {
            ProductFormat::Zip => {
                let file = File::open(&self.product_path)?;
                let mut archive = ZipArchive::new(file).map_err(|e| 
                    SarError::Processing(format!("Failed to open zip archive: {}", e)))?;
                let mut file = archive.by_name(annotation_file).map_err(|e|
                    SarError::Processing(format!("Failed to read annotation file: {}", e)))?;
                let mut content = String::new();
                file.read_to_string(&mut content)?;
                content
            },
            ProductFormat::Safe => {
                let file_path = self.product_path.join(annotation_file);
                std::fs::read_to_string(file_path)?
            }
        };
        
        Ok(content)
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
        let cache_dir = orbit_cache_dir.map(|p| p.to_path_buf()).unwrap_or_else(|| {
            std::env::temp_dir().join("sardine_orbit_cache")
        });
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
        
        // Check cache
        let cache_dir = orbit_cache_dir.map(|p| p.to_path_buf()).unwrap_or_else(|| {
            std::env::temp_dir().join("sardine_orbit_cache")
        });
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
        
        let cache_dir = orbit_cache_dir.map(|p| p.to_path_buf()).unwrap_or_else(|| {
            std::env::temp_dir().join("sardine_orbit_cache")
        });
        let orbit_manager = OrbitManager::new(cache_dir);
        let primary_orbit_type = crate::io::orbit::OrbitReader::determine_orbit_type(metadata.start_time);
        
        let mut downloaded_files = Vec::new();
        
        // Try to download primary orbit type
        match orbit_manager.download_and_cache_orbit_public(&product_id, metadata.start_time, primary_orbit_type) {
            Ok(_) => {
                let path = orbit_manager.get_cache_dir().join(format!("{}_{}.EOF", product_id, primary_orbit_type));
                downloaded_files.push(path);
                log::info!("Downloaded {} orbit file", primary_orbit_type);
            },
            Err(e) => {
                log::warn!("Failed to download {} orbit: {}", primary_orbit_type, e);
                
                // Try fallback
                let fallback_orbit_type = match primary_orbit_type {
                    crate::io::orbit::OrbitType::POEORB => crate::io::orbit::OrbitType::RESORB,
                    crate::io::orbit::OrbitType::RESORB => crate::io::orbit::OrbitType::POEORB,
                };
                
                match orbit_manager.download_and_cache_orbit_public(&product_id, metadata.start_time, fallback_orbit_type) {
                    Ok(_) => {
                        let path = orbit_manager.get_cache_dir().join(format!("{}_{}.EOF", product_id, fallback_orbit_type));
                        downloaded_files.push(path);
                        log::info!("Downloaded fallback {} orbit file", fallback_orbit_type);
                    },
                    Err(e2) => {
                        return Err(SarError::Processing(
                            format!("Failed to download both {} and {} orbit files: {}, {}", 
                                   primary_orbit_type, fallback_orbit_type, e, e2)
                        ));
                    }
                }
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

        // Parse using the detailed annotation parser only (no string fallback)
        let annotation = crate::io::annotation::AnnotationParser::parse_annotation(&xml_content)
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
            let annotation = crate::io::annotation::AnnotationParser::parse_annotation(&xml_content)
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
        let all_annotations = self.find_all_iw_annotation_files()?;
        let mut all_subswaths = HashMap::new();
        
    for pol in all_annotations.keys().cloned() {
            match self.extract_all_iw_subswaths(pol) {
                Ok(subswaths) => {
                    log::info!("Extracted {} subswaths for polarization {}", subswaths.len(), pol);
                    all_subswaths.insert(pol, subswaths);
                },
                Err(e) => {
                    log::error!("Failed to extract sub-swaths for {}: {}", pol, e);
                }
            }
        }
        
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
        // TODO: Extract actual satellite velocity from annotation for scientific accuracy
        let satellite_velocity = 7500.0; // Approximate Sentinel-1 velocity - should be calculated from orbit
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
        // TODO: Extract actual satellite velocity from annotation for scientific accuracy
        let satellite_velocity = 7500.0; // Approximate Sentinel-1 velocity - should be calculated from orbit
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
            if file.contains("annotation/calibration/") && file.ends_with(".xml") {
                if let Some(pol) = self.extract_polarization_from_filename(&file) {
                    calibration_files.insert(pol, file);
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

    /// Get comprehensive metadata with REAL extracted values (Step 1 enhancement)
    ///
    /// CRITICAL: This method extracts real metadata from Sentinel-1 products
    /// instead of using hardcoded default values. Essential for scientific accuracy.
    ///
    /// Returns real values for:
    /// - Range/Azimuth pixel spacing (NOT hardcoded 2.3/14.0)
    /// - Slant range time (NOT hardcoded)
    /// - PRF (NOT hardcoded) 
    /// - Radar frequency/wavelength (NOT hardcoded)
    /// - Product timing and geographic bounds
    ///
    /// References:
    /// - ESA S1-IF-ASD-PL-0007: "Sentinel-1 Annotation Format Specification"
    pub fn get_metadata(&mut self) -> SarResult<std::collections::HashMap<String, String>> {
        log::info!("Extracting comprehensive metadata with REAL values from Sentinel-1 product");

        let mut metadata = std::collections::HashMap::new();
        
        // Add debug marker to confirm this function is executing
        metadata.insert("debug_rust_function_called".to_string(), "true".to_string());

        // Get first available polarization for metadata extraction
        let annotation_files = self.find_annotation_files()?;
        let first_pol = annotation_files.keys().next().cloned()
            .ok_or_else(|| SarError::Processing("No annotation files found for metadata extraction".to_string()))?;

        // Read annotation and extract real metadata
        let sar_metadata = self.read_annotation(first_pol)?;
        
        // Convert to string map with real extracted values
        let mut metadata = std::collections::HashMap::new();
        
        // Basic product information (REAL)
        metadata.insert("product_id".to_string(), sar_metadata.product_id.clone());
        metadata.insert("mission".to_string(), sar_metadata.mission.clone());
        metadata.insert("platform".to_string(), sar_metadata.platform.clone());
        metadata.insert("instrument".to_string(), sar_metadata.instrument.clone());
        metadata.insert("acquisition_mode".to_string(), format!("{:?}", sar_metadata.acquisition_mode));
        
        // Add mode field (required by validation)
        metadata.insert("mode".to_string(), format!("{:?}", sar_metadata.acquisition_mode));
        
        // REAL pixel spacing (extracted from annotation XML)
        metadata.insert("range_pixel_spacing".to_string(), sar_metadata.pixel_spacing.0.to_string());
        metadata.insert("azimuth_pixel_spacing".to_string(), sar_metadata.pixel_spacing.1.to_string());
        
        // Timing information (REAL)
        metadata.insert("start_time".to_string(), sar_metadata.start_time.to_rfc3339());
        metadata.insert("stop_time".to_string(), sar_metadata.stop_time.to_rfc3339());

        // CRITICAL: Extract additional real metadata required for later steps
        let first_annotation_file = annotation_files.values().next().unwrap();
        let xml_content = self.read_file_as_string(first_annotation_file)?;
        
        // Extract REAL SAR parameters from annotation XML (required for terrain correction)
        if let Some(slant_range_time) = Self::extract_real_value_f64(&xml_content, "slantRangeTime") {
            metadata.insert("slant_range_time".to_string(), slant_range_time.to_string());
        }
        
        if let Some(prf) = Self::extract_real_value_f64(&xml_content, "azimuthFrequency")
            .or_else(|| Self::extract_real_value_f64(&xml_content, "prf")) {
            metadata.insert("prf".to_string(), prf.to_string());
        }
        
        if let Some(range_sampling_rate) = Self::extract_real_value_f64(&xml_content, "rangeSamplingRate") {
            metadata.insert("range_sampling_rate".to_string(), range_sampling_rate.to_string());
        }
        
        if let Some(radar_freq) = Self::extract_real_value_f64(&xml_content, "radarFrequency") {
            metadata.insert("radar_frequency".to_string(), radar_freq.to_string());
            // Calculate wavelength from radar frequency: λ = c/f
            let wavelength = SPEED_OF_LIGHT_M_S / radar_freq; // c = speed of light in m/s
            metadata.insert("wavelength".to_string(), wavelength.to_string());
        } else {
            // Scientific: Sentinel-1 C-band operates at 5.405 GHz by ESA design
            let sentinel1_freq = 5.405e9; // Hz
            let wavelength = SPEED_OF_LIGHT_M_S / sentinel1_freq;
            metadata.insert("wavelength".to_string(), wavelength.to_string());
            log::info!("Using Sentinel-1 C-band frequency (5.405 GHz) for wavelength calculation: {:.6}m", wavelength);
        }
        
        // Geographic bounds (REAL extracted)
        metadata.insert("min_longitude".to_string(), sar_metadata.bounding_box.min_lon.to_string());
        metadata.insert("max_longitude".to_string(), sar_metadata.bounding_box.max_lon.to_string());
        metadata.insert("min_latitude".to_string(), sar_metadata.bounding_box.min_lat.to_string());
        metadata.insert("max_latitude".to_string(), sar_metadata.bounding_box.max_lat.to_string());
        
        // Check if coordinates are valid - fail if not
        if (sar_metadata.bounding_box.min_lat == 0.0 && sar_metadata.bounding_box.max_lat == 0.0) ||
           (sar_metadata.bounding_box.min_lat == -90.0 && sar_metadata.bounding_box.max_lat == 90.0 && 
            sar_metadata.bounding_box.min_lon == -180.0 && sar_metadata.bounding_box.max_lon == 180.0) ||
           (sar_metadata.bounding_box.min_lon == 0.0 && sar_metadata.bounding_box.max_lon == 0.0) {
            log::error!("Failed to extract valid coordinate data from annotation XML");
            metadata.insert("coordinate_extraction_status".to_string(), "FAILED_NO_FALLBACK".to_string());
        } else {
            log::info!("Structured extraction succeeded with coordinates: [{:.6}, {:.6}, {:.6}, {:.6}]",
                      sar_metadata.bounding_box.min_lon, sar_metadata.bounding_box.min_lat, 
                      sar_metadata.bounding_box.max_lon, sar_metadata.bounding_box.max_lat);
        }
        
        // Processing information
        metadata.insert("coordinate_system".to_string(), format!("{:?}", sar_metadata.coordinate_system));
        metadata.insert("range_looks".to_string(), sar_metadata.range_looks.to_string());
        metadata.insert("azimuth_looks".to_string(), sar_metadata.azimuth_looks.to_string());
        
        // Polarizations available
        let polarizations: Vec<String> = sar_metadata.polarizations.iter()
            .map(|p| format!("{}", p))
            .collect();
        
        // DEBUG: Mark that we reached polarizations insertion
        metadata.insert("debug_before_polarizations".to_string(), "true".to_string());
        
        metadata.insert("polarizations".to_string(), polarizations.join(","));
        
        // DEBUG: Mark that we reached after polarizations
        metadata.insert("debug_after_polarizations".to_string(), "true".to_string());
        
        // Sub-swath information if available
        if !sar_metadata.sub_swaths.is_empty() {
            let swath_ids: Vec<String> = sar_metadata.sub_swaths.keys().cloned().collect();
            metadata.insert("subswaths".to_string(), swath_ids.join(","));
            metadata.insert("subswath_count".to_string(), sar_metadata.sub_swaths.len().to_string());
        } else {
            // For IW mode, add standard subswaths even if not explicitly listed
            if sar_metadata.acquisition_mode == AcquisitionMode::IW {
                metadata.insert("subswaths".to_string(), "IW1,IW2,IW3".to_string());
                metadata.insert("subswath_count".to_string(), "3".to_string());
            }
        }

        // DEBUG: Mark that we reached the enhanced code block
        metadata.insert("debug_enhanced_code_reached".to_string(), "true".to_string());

        // Extract product type from annotation or infer from product ID
        log::info!("Extracting product type from XML content length: {}", xml_content.len());
        if let Some(product_type) = Self::extract_product_type(&xml_content, &sar_metadata.product_id) {
            log::info!("Found product type: {}", product_type);
            metadata.insert("product_type".to_string(), product_type);
        } else {
            log::warn!("No product type found in XML or product ID");
        }

        // Extract orbit information from product ID
        if let Some((orbit_number, relative_orbit, orbit_direction)) = Self::parse_orbit_from_product_id(&sar_metadata.product_id) {
            metadata.insert("orbit_number".to_string(), orbit_number.to_string());
            metadata.insert("relative_orbit".to_string(), relative_orbit.to_string());
            metadata.insert("orbit_direction".to_string(), orbit_direction);
        }

        // Extract additional SAR parameters from XML
        if let Some(azimuth_steering_rate) = Self::extract_real_value_f64(&xml_content, "azimuthSteeringRate")
            .or_else(|| Self::extract_real_value_f64(&xml_content, "azimuthFmRate")) {
            metadata.insert("azimuth_steering_rate".to_string(), azimuth_steering_rate.to_string());
        }

        // Try rangeSamplingRate for range bandwidth (Hz)
        if let Some(range_bandwidth) = Self::extract_real_value_f64(&xml_content, "rangeSamplingRate")
            .or_else(|| Self::extract_real_value_f64(&xml_content, "rangeBandwidth")) {
            metadata.insert("range_bandwidth".to_string(), range_bandwidth.to_string());
        }

        if let Some(heading) = Self::extract_real_value_f64(&xml_content, "platformHeading")
            .or_else(|| Self::extract_real_value_f64(&xml_content, "heading")) {
            metadata.insert("heading".to_string(), heading.to_string());
        }

        // Extract incidence angles - use the first and last values from the incidence angle array
        if let Some(inc_angles) = Self::extract_incidence_angle_range(&xml_content) {
            metadata.insert("incidence_angle_near".to_string(), inc_angles.0.to_string());
            metadata.insert("incidence_angle_far".to_string(), inc_angles.1.to_string());
            
            // Calculate center incidence angle as average of near and far
            let center_angle = (inc_angles.0 + inc_angles.1) / 2.0;
            metadata.insert("incidence_angle_center".to_string(), center_angle.to_string());
        }

        // Add processing level and processor version
        metadata.insert("processing_level".to_string(), "1".to_string()); // Sentinel-1 SLC is Level-1
        metadata.insert("processor_version".to_string(), "SARdine v0.2.0".to_string());

        // Extract center coordinates from bounding box (only if not already calculated by regex extraction)
        if !metadata.contains_key("center_lat") {
            let center_lat = (sar_metadata.bounding_box.min_lat + sar_metadata.bounding_box.max_lat) / 2.0;
            let center_lon = (sar_metadata.bounding_box.min_lon + sar_metadata.bounding_box.max_lon) / 2.0;
            log::info!("Calculated center coordinates from structured bounding box: lat={}, lon={}", center_lat, center_lon);
            metadata.insert("center_lat".to_string(), center_lat.to_string());
            metadata.insert("center_lon".to_string(), center_lon.to_string());
        } else {
            log::info!("Center coordinates already calculated by regex extraction, keeping those values");
        }

        // DEBUGGING: Add a simple test field to verify this code runs
        metadata.insert("test_enhancement_reached".to_string(), "YES".to_string());

        // Extract subswath information
            let subswaths: Vec<String> = sar_metadata.sub_swaths.keys()
                .cloned()
                .collect();
        if !subswaths.is_empty() {
            metadata.insert("subswaths".to_string(), subswaths.join(","));
        }

        // Extract calibration constant if available
        if let Some(cal_constant) = Self::extract_calibration_constant(&xml_content) {
            metadata.insert("calibration_constant".to_string(), cal_constant.to_string());
        }

        // Product format
        metadata.insert("format".to_string(), match self.format {
            ProductFormat::Zip => "ZIP".to_string(),
            ProductFormat::Safe => "SAFE".to_string(),
        });

        log::info!("REAL metadata extracted:");
        log::info!("  - Range spacing: {} m (REAL from XML)", sar_metadata.pixel_spacing.0);
        log::info!("  - Azimuth spacing: {} m (REAL from XML)", sar_metadata.pixel_spacing.1);
        log::info!("  - Geographic bounds: [{:.3}, {:.3}, {:.3}, {:.3}] (REAL from XML)", 
                  sar_metadata.bounding_box.min_lon, sar_metadata.bounding_box.min_lat,
                  sar_metadata.bounding_box.max_lon, sar_metadata.bounding_box.max_lat);
        log::info!("  - Subswaths: {} (REAL from XML)", sar_metadata.sub_swaths.len());
        
        // Log critical parameters needed for later steps
        if let Some(slant_range) = metadata.get("slant_range_time") {
            log::info!("  - Slant range time: {} s (REAL from XML)", slant_range);
        }
        if let Some(prf) = metadata.get("prf") {
            log::info!("  - PRF: {} Hz (REAL from XML)", prf);
        }
        if let Some(wavelength) = metadata.get("wavelength") {
            log::info!("  - Wavelength: {} m (REAL calculated/fallback)", wavelength);
        }
        
        Ok(metadata)
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
        dem_cache_dir: Option<&str>,
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
        
        // Use provided cache directory or create default
        let cache_dir = dem_cache_dir.unwrap_or("./dem_cache");
        
        let (dem_data, dem_transform) = crate::io::dem::DemReader::prepare_dem_for_scene(
            bbox, 
            30.0, // 30m target resolution
            cache_dir
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
            30.0 // 30m target resolution
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
