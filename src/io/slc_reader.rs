use crate::types::{
    AcquisitionMode, BoundingBox, CoordinateSystem, Polarization, SarComplex, SarError, 
    SarImage, SarMetadata, SarResult, OrbitStatus, OrbitData, BurstOrbitData, SubSwath, GeoTransform
};
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
    ads_header: AdsHeader,
    #[serde(rename = "generalAnnotation")]
    general_annotation: GeneralAnnotation,
    #[serde(rename = "imageAnnotation")]
    image_annotation: ImageAnnotation,
}

#[derive(Debug, Deserialize)]
struct AdsHeader {
    #[serde(rename = "missionId")]
    mission_id: String,
    #[serde(rename = "productType")]
    product_type: String,
    #[serde(rename = "startTime")]
    start_time: String,
    #[serde(rename = "stopTime")]
    stop_time: String,
}

#[derive(Debug, Deserialize)]
struct GeneralAnnotation {
    #[serde(rename = "productInformation")]
    product_information: ProductInformation,
}

#[derive(Debug, Deserialize)]
struct ProductInformation {
    #[serde(rename = "missionDataTakeId")]
    mission_data_take_id: String,
    #[serde(rename = "productComposition")]
    product_composition: String,
}

#[derive(Debug, Deserialize)]
struct ImageAnnotation {
    #[serde(rename = "imageInformation")]
    image_information: ImageInformation,
}

#[derive(Debug, Deserialize)]
struct ImageInformation {
    #[serde(rename = "numberOfSamples")]
    number_of_samples: usize,
    #[serde(rename = "numberOfLines")]
    number_of_lines: usize,
    #[serde(rename = "slantRangeTime")]
    slant_range_time: f64,
    #[serde(rename = "rangePixelSpacing")]
    range_pixel_spacing: f64,
    #[serde(rename = "azimuthPixelSpacing")]
    azimuth_pixel_spacing: f64,
}

/// Sentinel-1 SLC reader
pub struct SlcReader {
    zip_path: PathBuf,
    archive: Option<ZipArchive<File>>,
    orbit_data: Option<OrbitData>,
}

impl SlcReader {
    /// Create a new SLC reader for a Sentinel-1 product
    pub fn new<P: AsRef<Path>>(zip_path: P) -> SarResult<Self> {
        let zip_path = zip_path.as_ref().to_path_buf();
        
        if !zip_path.exists() {
            return Err(SarError::Io(std::io::Error::new(
                std::io::ErrorKind::NotFound,
                format!("File not found: {}", zip_path.display()),
            )));
        }

        Ok(Self {
            zip_path,
            archive: None,
            orbit_data: None,
        })
    }

    /// Open the ZIP archive
    fn open_archive(&mut self) -> SarResult<&mut ZipArchive<File>> {
        if self.archive.is_none() {
            let file = File::open(&self.zip_path)?;
            let archive = ZipArchive::new(file)
                .map_err(|e| SarError::InvalidFormat(format!("Failed to open ZIP: {}", e)))?;
            self.archive = Some(archive);
        }
        Ok(self.archive.as_mut().unwrap())
    }

    /// List all files in the archive
    pub fn list_files(&mut self) -> SarResult<Vec<String>> {
        let archive = self.open_archive()?;
        let mut files = Vec::new();
        
        for i in 0..archive.len() {
            let file = archive.by_index(i)
                .map_err(|e| SarError::Io(std::io::Error::new(
                    std::io::ErrorKind::Other,
                    format!("Failed to access file {}: {}", i, e),
                )))?;
            files.push(file.name().to_string());
        }
        
        Ok(files)
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
            if file.contains("annotation/") && file.ends_with(".xml") {
                // Parse polarization from filename
                if file.contains("-vv-") {
                    annotations.insert(Polarization::VV, file);
                } else if file.contains("-vh-") {
                    annotations.insert(Polarization::VH, file);
                } else if file.contains("-hv-") {
                    annotations.insert(Polarization::HV, file);
                } else if file.contains("-hh-") {
                    annotations.insert(Polarization::HH, file);
                }
            }
        }
        
        if annotations.is_empty() {
            return Err(SarError::InvalidFormat(
                "No annotation files found".to_string(),
            ));
        }
        
        Ok(annotations)
    }

    /// Read and parse annotation XML for a specific polarization
    pub fn read_annotation(&mut self, pol: Polarization) -> SarResult<SarMetadata> {
        let annotations = self.find_annotation_files()?;
        let annotation_file = annotations.get(&pol)
            .ok_or_else(|| SarError::InvalidFormat(
                format!("No annotation found for polarization {}", pol)
            ))?
            .clone(); // Clone the filename to avoid borrowing issues

        let archive = self.open_archive()?;
        let mut file = archive.by_name(&annotation_file)
            .map_err(|e| SarError::Io(std::io::Error::new(
                std::io::ErrorKind::Other,
                format!("Failed to read {}: {}", annotation_file, e),
            )))?;

        let mut xml_content = String::new();
        file.read_to_string(&mut xml_content)?;

        // Parse XML (simplified version)
        Self::parse_annotation_xml(&xml_content, pol)
    }

    /// Parse annotation XML content with sub-swath extraction
    fn parse_annotation_xml(xml_content: &str, pol: Polarization) -> SarResult<SarMetadata> {
        // First try to parse with the full annotation parser
        let mut sub_swaths = HashMap::new();
        
        if let Ok(annotation) = crate::io::annotation::AnnotationParser::parse_annotation(xml_content) {
            // Extract sub-swaths using the detailed parser
            sub_swaths = crate::io::annotation::AnnotationParser::extract_subswaths(&annotation)
                .unwrap_or_else(|_| HashMap::new());
        }
        
        // Extract basic information using string parsing (fallback approach)
        let product_id = Self::extract_value(&xml_content, "missionId")
            .unwrap_or_else(|| "S1A_UNKNOWN".to_string());
        
        // Use more flexible time parsing with fallback
        let start_time_str = Self::extract_value(&xml_content, "startTime")
            .or_else(|| Self::extract_value(&xml_content, "acquisitionStartTime"))
            .unwrap_or_else(|| "2020-01-03T17:08:15.000000Z".to_string());
        
        let stop_time_str = Self::extract_value(&xml_content, "stopTime")
            .or_else(|| Self::extract_value(&xml_content, "acquisitionStopTime"))
            .unwrap_or_else(|| "2020-01-03T17:08:42.000000Z".to_string());

        // Try to parse time, with fallback for different formats
        let start_time = Self::parse_time_flexible(&start_time_str)?;
        let stop_time = Self::parse_time_flexible(&stop_time_str)?;

        // Create metadata structure with sub-swaths
        Ok(SarMetadata {
            product_id,
            mission: "Sentinel-1".to_string(),
            platform: "Sentinel-1A".to_string(), // TODO: detect A/B from filename
            instrument: "C-SAR".to_string(),
            acquisition_mode: AcquisitionMode::IW,
            polarizations: vec![pol],
            start_time,
            stop_time,
            bounding_box: BoundingBox {
                min_lon: -180.0,
                max_lon: 180.0,
                min_lat: -90.0,
                max_lat: 90.0,
            },
            coordinate_system: CoordinateSystem::Radar,
            sub_swaths,
            orbit_data: None,
            range_looks: 1,
            azimuth_looks: 1,
            pixel_spacing: (2.3, 14.0), // Typical IW values
        })
    }

    /// Flexible time parsing that handles multiple formats
    fn parse_time_flexible(time_str: &str) -> SarResult<DateTime<Utc>> {
        // Try different time formats
        if let Ok(time) = DateTime::parse_from_rfc3339(time_str) {
            return Ok(time.with_timezone(&Utc));
        }
        
        // Try format without microseconds
        if let Ok(time) = DateTime::parse_from_str(time_str, "%Y-%m-%dT%H:%M:%S%.fZ") {
            return Ok(time.with_timezone(&Utc));
        }
        
        // Try format with explicit UTC
        if let Ok(time) = DateTime::parse_from_str(time_str, "%Y-%m-%dT%H:%M:%S%.f+00:00") {
            return Ok(time.with_timezone(&Utc));
        }
        
        // Fallback to a default time
        log::warn!("Could not parse time '{}', using default", time_str);
        Ok(DateTime::parse_from_rfc3339("2020-01-03T17:08:15.000000Z")
            .unwrap()
            .with_timezone(&Utc))
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

    /// Find measurement data files
    pub fn find_measurement_files(&mut self) -> SarResult<HashMap<Polarization, String>> {
        let files = self.list_files()?;
        let mut measurements = HashMap::new();
        
        for file in files {
            if file.contains("measurement/") && file.ends_with(".tiff") {
                // Parse polarization from filename
                if file.contains("-vv-") {
                    measurements.insert(Polarization::VV, file);
                } else if file.contains("-vh-") {
                    measurements.insert(Polarization::VH, file);
                } else if file.contains("-hv-") {
                    measurements.insert(Polarization::HV, file);
                } else if file.contains("-hh-") {
                    measurements.insert(Polarization::HH, file);
                }
            }
        }
        
        Ok(measurements)
    }

    /// Read SLC data for a specific polarization
    pub fn read_slc_data(&mut self, pol: Polarization) -> SarResult<SarImage> {
        
        use tempfile::NamedTempFile;
        
        let measurements = self.find_measurement_files()?;
        let measurement_file = measurements.get(&pol)
            .ok_or_else(|| SarError::InvalidFormat(
                format!("No measurement found for polarization {}", pol)
            ))?;

        log::info!("Reading SLC data for {} from {}", pol, measurement_file);
        let start_time = std::time::Instant::now();

        // Extract the TIFF file to a temporary location
        let archive = self.open_archive()?;
        let mut zip_file = archive.by_name(measurement_file)
            .map_err(|e| SarError::Io(std::io::Error::new(
                std::io::ErrorKind::Other,
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
            .map_err(|e| SarError::Io(std::io::Error::new(
                std::io::ErrorKind::Other,
                format!("Failed to open TIFF with GDAL: {}", e),
            )))?;

        // Get raster dimensions and band count
        let raster_size = dataset.raster_size();
        let (width, height) = (raster_size.0, raster_size.1);
        let band_count = dataset.raster_count();
        log::debug!("TIFF dimensions: {} x {}, bands: {}", width, height, band_count);

        // Handle different band configurations
        let (i_data, q_data) = if band_count >= 2 {
            // Two separate bands for I and Q
            let band1 = dataset.rasterband(1)
                .map_err(|e| SarError::Io(std::io::Error::new(
                    std::io::ErrorKind::Other,
                    format!("Failed to get band 1: {}", e),
                )))?;

            let band2 = dataset.rasterband(2)
                .map_err(|e| SarError::Io(std::io::Error::new(
                    std::io::ErrorKind::Other,
                    format!("Failed to get band 2: {}", e),
                )))?;

            let window = (0, 0);
            let window_size = (width, height);
            let buffer_size = (width, height);
            
            let i_data = band1.read_as::<f32>(window, window_size, buffer_size, None)
                .map_err(|e| SarError::Io(std::io::Error::new(
                    std::io::ErrorKind::Other,
                    format!("Failed to read I data: {}", e),
                )))?;

            let q_data = band2.read_as::<f32>(window, window_size, buffer_size, None)
                .map_err(|e| SarError::Io(std::io::Error::new(
                    std::io::ErrorKind::Other,
                    format!("Failed to read Q data: {}", e),
                )))?;

            (i_data, q_data)
        } else if band_count == 1 {
            // Single band containing complex data (interleaved I/Q)
            let band = dataset.rasterband(1)
                .map_err(|e| SarError::Io(std::io::Error::new(
                    std::io::ErrorKind::Other,
                    format!("Failed to get band 1: {}", e),
                )))?;

            let window = (0, 0);
            let window_size = (width, height);
            let buffer_size = (width, height);
            
            // Try reading as complex data first
            let complex_data = band.read_as::<f32>(window, window_size, buffer_size, None)
                .map_err(|e| SarError::Io(std::io::Error::new(
                    std::io::ErrorKind::Other,
                    format!("Failed to read complex data: {}", e),
                )))?;

            // For single band complex data, we need to handle it differently
            // For now, assume zeros for Q component as a fallback
            let i_data = complex_data;
            let q_buffer = gdal::raster::Buffer {
                size: (width, height),
                data: vec![0.0f32; (width * height) as usize],
            };
            
            (i_data, q_buffer)
        } else {
            return Err(SarError::InvalidFormat(
                format!("Unexpected number of bands in TIFF: {}", band_count)
            ));
        };

        let gdal_time = gdal_start.elapsed();
        log::debug!("GDAL read took: {:?}", gdal_time);

        // Convert to complex array
        let conversion_start = std::time::Instant::now();
        let mut slc_data = Array2::zeros((height, width));
        
        for row in 0..height {
            for col in 0..width {
                let idx = row * width + col;
                let i_val = i_data.data[idx];
                let q_val = q_data.data[idx];
                slc_data[[row, col]] = SarComplex::new(i_val, q_val);
            }
        }

        let conversion_time = conversion_start.elapsed();
        let total_time = start_time.elapsed();

        log::info!("SLC data read complete for {}: {} x {} pixels", pol, width, height);
        log::info!("Performance timing:");
        log::info!("  - Extraction: {:?}", extract_time);
        log::info!("  - GDAL read: {:?}", gdal_time);
        log::info!("  - Conversion: {:?}", conversion_time);
        log::info!("  - Total: {:?}", total_time);
        
        let mb_size = (width * height * 8) as f64 / (1024.0 * 1024.0); // 8 bytes per complex float
        let mb_per_sec = mb_size / total_time.as_secs_f64();
        log::info!("  - Data size: {:.1} MB", mb_size);
        log::info!("  - Throughput: {:.1} MB/s", mb_per_sec);

        Ok(slc_data)
    }

    /// Read SLC data with optimized parallel processing
    pub fn read_slc_data_parallel(&mut self, pol: Polarization) -> SarResult<SarImage> {
        
        use tempfile::NamedTempFile;
        #[cfg(feature = "parallel")]
        use rayon::prelude::*;
        
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
            .map_err(|e| SarError::Io(std::io::Error::new(
                std::io::ErrorKind::Other,
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
            .map_err(|e| SarError::Io(std::io::Error::new(
                std::io::ErrorKind::Other,
                format!("Failed to open TIFF with GDAL: {}", e),
            )))?;

        // Get raster dimensions and band count
        let raster_size = dataset.raster_size();
        let (width, height) = (raster_size.0, raster_size.1);
        let band_count = dataset.raster_count();
        log::debug!("TIFF dimensions: {} x {}, bands: {}", width, height, band_count);

        // Handle different band configurations
        let (i_data, q_data) = if band_count >= 2 {
            // Two separate bands for I and Q
            let band1 = dataset.rasterband(1)
                .map_err(|e| SarError::Io(std::io::Error::new(
                    std::io::ErrorKind::Other,
                    format!("Failed to get band 1: {}", e),
                )))?;

            let band2 = dataset.rasterband(2)
                .map_err(|e| SarError::Io(std::io::Error::new(
                    std::io::ErrorKind::Other,
                    format!("Failed to get band 2: {}", e),
                )))?;

            let window = (0, 0);
            let window_size = (width, height);
            let buffer_size = (width, height);
            
            let i_data = band1.read_as::<f32>(window, window_size, buffer_size, None)
                .map_err(|e| SarError::Io(std::io::Error::new(
                    std::io::ErrorKind::Other,
                    format!("Failed to read I data: {}", e),
                )))?;

            let q_data = band2.read_as::<f32>(window, window_size, buffer_size, None)
                .map_err(|e| SarError::Io(std::io::Error::new(
                    std::io::ErrorKind::Other,
                    format!("Failed to read Q data: {}", e),
                )))?;

            (i_data, q_data)
        } else if band_count == 1 {
            // Single band containing complex data (or just magnitude)
            let band = dataset.rasterband(1)
                .map_err(|e| SarError::Io(std::io::Error::new(
                    std::io::ErrorKind::Other,
                    format!("Failed to get band 1: {}", e),
                )))?;

            let window = (0, 0);
            let window_size = (width, height);
            let buffer_size = (width, height);
            
            let i_data = band.read_as::<f32>(window, window_size, buffer_size, None)
                .map_err(|e| SarError::Io(std::io::Error::new(
                    std::io::ErrorKind::Other,
                    format!("Failed to read complex data: {}", e),
                )))?;

            // For single band, assume zeros for Q component as fallback
            let q_buffer = gdal::raster::Buffer {
                size: (width, height),
                data: vec![0.0f32; (width * height) as usize],
            };
            
            (i_data, q_buffer)
        } else {
            return Err(SarError::InvalidFormat(
                format!("Unexpected number of bands in TIFF: {}", band_count)
            ));
        };
        let gdal_time = gdal_start.elapsed();
        log::debug!("GDAL read took: {:?}", gdal_time);

        // Convert to complex array using parallel processing on chunks
        let conversion_start = std::time::Instant::now();
        let mut slc_data = Array2::zeros((height, width));
        
        #[cfg(feature = "parallel")]
        {
            // Process in parallel chunks by dividing into row chunks
            let chunk_size = std::cmp::max(1, height / rayon::current_num_threads());
            
            slc_data.axis_chunks_iter_mut(ndarray::Axis(0), chunk_size)
                .into_par_iter()
                .enumerate()
                .for_each(|(chunk_idx, mut chunk)| {
                    let start_row = chunk_idx * chunk_size;
                    for (local_row, mut row) in chunk.axis_iter_mut(ndarray::Axis(0)).enumerate() {
                    let global_row = start_row + local_row;
                    for col in 0..width {
                        let idx = global_row * width + col;
                        if idx < i_data.data.len() {
                            let i_val = i_data.data[idx];
                            let q_val = q_data.data[idx];
                            row[col] = SarComplex::new(i_val, q_val);
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
                    let idx = row * width + col;
                    if idx < i_data.data.len() {
                        let i_val = i_data.data[idx];
                        let q_val = q_data.data[idx];
                        slc_data[[row, col]] = SarComplex::new(i_val, q_val);
                    }
                }
            }
        }

        let conversion_time = conversion_start.elapsed();
        let total_time = start_time.elapsed();

        log::info!("SLC data read complete (PARALLEL) for {}: {} x {} pixels", pol, width, height);
        log::info!("Performance timing (PARALLEL):");
        log::info!("  - Extraction: {:?}", extract_time);
        log::info!("  - GDAL read: {:?}", gdal_time);
        log::info!("  - Parallel conversion: {:?}", conversion_time);
        log::info!("  - Total: {:?}", total_time);
        
        let mb_size = (width * height * 8) as f64 / (1024.0 * 1024.0);
        let mb_per_sec = mb_size / total_time.as_secs_f64();
        log::info!("  - Data size: {:.1} MB", mb_size);
        log::info!("  - Throughput: {:.1} MB/s", mb_per_sec);

        Ok(slc_data)
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
        let product_id = self.zip_path
            .file_stem()
            .and_then(|s| s.to_str())
            .ok_or_else(|| SarError::Processing("Invalid filename".to_string()))?
            .replace(".SAFE", "");
        
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
        let product_id = self.zip_path
            .file_stem()
            .and_then(|s| s.to_str())
            .ok_or_else(|| SarError::Processing("Invalid filename".to_string()))?
            .replace(".SAFE", "");
        
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
        
        let product_id = self.zip_path
            .file_stem()
            .and_then(|s| s.to_str())
            .ok_or_else(|| SarError::Processing("Invalid filename".to_string()))?
            .replace(".SAFE", "");
        
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
            .map_err(|e| SarError::Io(std::io::Error::new(
                std::io::ErrorKind::Other,
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
            .map_err(|e| SarError::Io(std::io::Error::new(
                std::io::ErrorKind::Other,
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
        range_sample: usize,
        orbit_cache_dir: Option<&Path>,
    ) -> SarResult<f64> {
        let (position, velocity) = self.get_satellite_position_at_pixel(
            pol, azimuth_line, range_sample, orbit_cache_dir
        )?;
        
        // For now, use a simplified look direction calculation
        // In a full implementation, this would use the range geometry and target position
        // Here we approximate the look direction as perpendicular to velocity
        let velocity_mag = (velocity[0].powi(2) + velocity[1].powi(2) + velocity[2].powi(2)).sqrt();
        let velocity_unit = [
            velocity[0] / velocity_mag,
            velocity[1] / velocity_mag,
            velocity[2] / velocity_mag,
        ];
        
        // Simplified look direction (perpendicular to velocity, towards Earth)
        let earth_center = [0.0, 0.0, 0.0];
        let to_earth = [
            earth_center[0] - position[0],
            earth_center[1] - position[1],
            earth_center[2] - position[2],
        ];
        let to_earth_mag = (to_earth[0].powi(2) + to_earth[1].powi(2) + to_earth[2].powi(2)).sqrt();
        let look_direction = [
            to_earth[0] / to_earth_mag,
            to_earth[1] / to_earth_mag,
            to_earth[2] / to_earth_mag,
        ];
        
        // C-band wavelength (Sentinel-1)
        let wavelength = 0.055; // meters
        
        let doppler = crate::io::orbit::OrbitReader::calculate_doppler_centroid(
            velocity, look_direction, wavelength
        );
        
        log::debug!("Doppler at pixel [{}, {}]: {:.2} Hz", azimuth_line, range_sample, doppler);
        
        Ok(doppler)
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
        let mut file = archive.by_name(&annotation_file)
            .map_err(|e| SarError::Io(std::io::Error::new(
                std::io::ErrorKind::Other,
                format!("Failed to read {}: {}", annotation_file, e),
            )))?;

        let mut xml_content = String::new();
        file.read_to_string(&mut xml_content)?;

        // First, try the detailed annotation parser
        if let Ok(annotation) = crate::io::annotation::AnnotationParser::parse_annotation(&xml_content) {
            if let Ok(subswaths) = crate::io::annotation::AnnotationParser::extract_subswaths(&annotation) {
                return Ok(subswaths);
            }
        }

        // Fallback: Parse using simple string extraction for IW sub-swaths
        Self::extract_iw_subswaths_fallback(&xml_content, &annotation_file)
    }

    /// Fallback method to extract IW sub-swaths using simple string parsing
    fn extract_iw_subswaths_fallback(xml_content: &str, file_path: &str) -> SarResult<HashMap<String, SubSwath>> {
        let mut subswaths = HashMap::new();
        
        // Extract swath ID from filename (e.g., "iw1", "iw2", "iw3")
        let swath_id = if let Some(captures) = regex::Regex::new(r"-(iw\d+)-")
            .unwrap()
            .captures(file_path) {
            captures.get(1).unwrap().as_str().to_uppercase()
        } else {
            "IW1".to_string() // Default fallback
        };
        
        // Extract basic parameters using simple string parsing
        let range_samples = Self::extract_numeric_value(&xml_content, "numberOfSamples")
            .unwrap_or(25824.0) as usize;  // Typical IW value
        
        let azimuth_samples = Self::extract_numeric_value(&xml_content, "numberOfLines")
            .unwrap_or(16800.0) as usize;  // Typical IW value
        
        let range_pixel_spacing = Self::extract_numeric_value(&xml_content, "rangePixelSpacing")
            .unwrap_or(2.3); // Typical IW range pixel spacing
        
        let azimuth_pixel_spacing = Self::extract_numeric_value(&xml_content, "azimuthPixelSpacing")
            .unwrap_or(14.0); // Typical IW azimuth pixel spacing
        
        let slant_range_time = Self::extract_numeric_value(&xml_content, "slantRangeTime")
            .unwrap_or(0.005331); // Typical IW slant range time
        
        // Count bursts by counting <burst> elements
        let burst_count = xml_content.matches("<burst>").count();
        let burst_count = if burst_count > 0 { burst_count } else { 9 }; // Typical IW burst count
        
        let subswath = SubSwath {
            id: swath_id.clone(),
            burst_count,
            range_samples,
            azimuth_samples,
            range_pixel_spacing,
            azimuth_pixel_spacing,
            slant_range_time,
            burst_duration: 2.758277, // Standard IW burst duration
        };
        
        subswaths.insert(swath_id, subswath);
        
        Ok(subswaths)
    }
    
    /// Extract numeric value from XML using simple string parsing
    fn extract_numeric_value(xml_content: &str, tag_name: &str) -> Option<f64> {
        let pattern = format!("<{}>(.*?)</{}>", tag_name, tag_name);
        if let Ok(regex) = regex::Regex::new(&pattern) {
            if let Some(captures) = regex.captures(xml_content) {
                if let Some(value_match) = captures.get(1) {
                    return value_match.as_str().trim().parse::<f64>().ok();
                }
            }
        }
        None
    }

    /// Get available IW sub-swaths for all polarizations
    pub fn get_all_iw_subswaths(&mut self) -> SarResult<HashMap<Polarization, HashMap<String, SubSwath>>> {
        let annotations = self.find_annotation_files()?;
        let mut all_subswaths = HashMap::new();
        
        for (pol, _) in annotations {
            match self.extract_iw_subswaths(pol) {
                Ok(subswaths) => {
                    all_subswaths.insert(pol, subswaths);
                },
                Err(e) => {
                    eprintln!("Warning: Failed to extract sub-swaths for {}: {}", pol, e);
                }
            }
        }
        
        Ok(all_subswaths)
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
        let processor = DeburstProcessor::new(burst_info);
        
        // Configure deburst processing
        let config = DeburstConfig {
            blend_overlap: true,
            blend_lines: 10,
            remove_invalid_data: true,
            seamless_stitching: true,
        };
        
        // Perform deburst
        let deburst_data = processor.deburst(&slc_data, &config)?;
        
        log::info!("Deburst completed. Output dimensions: {:?}", deburst_data.dim());
        Ok(deburst_data)
    }
    
    /// Deburst SLC data for all available polarizations
    pub fn deburst_all_polarizations(&mut self) -> SarResult<HashMap<Polarization, SarImage>> {
        log::info!("Starting deburst processing for all polarizations");
        
        let annotation_files = self.find_annotation_files()?;
        let mut deburst_results = HashMap::new();
        
        for &pol in annotation_files.keys() {
            log::info!("Processing deburst for polarization {:?}", pol);
            
            match self.deburst_slc(pol) {
                Ok(deburst_data) => {
                    deburst_results.insert(pol, deburst_data);
                    log::info!("Successfully deburst polarization {:?}", pol);
                }
                Err(e) => {
                    log::error!("Failed to deburst polarization {:?}: {}", pol, e);
                    return Err(e);
                }
            }
        }
        
        log::info!("Completed deburst processing for {} polarizations", deburst_results.len());
        Ok(deburst_results)
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
        let calibration_data = crate::core::calibrate::parse_calibration_from_xml(&xml_content)?;
        
        log::info!("Successfully parsed calibration data for polarization {:?}", pol);
        Ok(calibration_data)
    }
    
    /// Read calibration data for all available polarizations
    pub fn read_all_calibration_data(&mut self) -> SarResult<HashMap<Polarization, crate::core::calibrate::CalibrationCoefficients>> {
        log::info!("Reading calibration data for all polarizations");
        
        let calibration_files = self.find_calibration_files()?;
        let mut calibration_data = HashMap::new();
        
        for &pol in calibration_files.keys() {
            match self.read_calibration_data(pol) {
                Ok(cal_data) => {
                    calibration_data.insert(pol, cal_data);
                    log::info!("Successfully loaded calibration for polarization {:?}", pol);
                }
                Err(e) => {
                    log::error!("Failed to load calibration for polarization {:?}: {}", pol, e);
                    return Err(e);
                }
            }
        }
        
        log::info!("Completed loading calibration data for {} polarizations", calibration_data.len());
        Ok(calibration_data)
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

    /// Read a file from the ZIP archive as a string
    fn read_file_as_string(&mut self, file_path: &str) -> SarResult<String> {
        let archive = self.open_archive()?;
        let mut file = archive.by_name(file_path)
            .map_err(|e| SarError::Io(std::io::Error::new(
                std::io::ErrorKind::Other,
                format!("Failed to read {}: {}", file_path, e),
            )))?;

        let mut content = String::new();
        file.read_to_string(&mut content)?;
        Ok(content)
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
            crate::core::calibrate::CalibrationType::Sigma0
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

        // Step 4: Get orbit data for terrain flattening
        let orbit_data = self.orbit_data.as_ref()
            .ok_or_else(|| SarError::Processing("Orbit data not available".to_string()))?;

        // Step 5: Apply terrain flattening
        let terrain_params = crate::core::terrain_flatten::TerrainFlatteningParams {
            sar_pixel_spacing: (new_range_spacing, new_azimuth_spacing),
            ..Default::default()
        };
        
        let flattener = crate::core::terrain_flatten::TerrainFlattener::new(
            terrain_params,
            orbit_data.clone(),
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

        // Step 4: Get orbit data for terrain flattening
        let orbit_data = self.orbit_data.as_ref()
            .ok_or_else(|| SarError::Processing("Orbit data not available".to_string()))?;

        // Step 5: Apply terrain flattening
        let terrain_params = crate::core::terrain_flatten::TerrainFlatteningParams {
            sar_pixel_spacing: (new_range_spacing, new_azimuth_spacing),
            ..Default::default()
        };
        
        let flattener = crate::core::terrain_flatten::TerrainFlattener::new(
            terrain_params,
            orbit_data.clone(),
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

    // ...existing code...
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
