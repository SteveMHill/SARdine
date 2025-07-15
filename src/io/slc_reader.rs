use crate::types::{
    AcquisitionMode, BoundingBox, CoordinateSystem, Polarization, SarComplex, SarError, 
    SarImage, SarMetadata, SarResult, SubSwath
};
use chrono::{DateTime, Utc};
use ndarray::Array2;
use quick_xml::de::from_str;
use serde::Deserialize;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufReader, Read};
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

    /// Parse annotation XML content
    fn parse_annotation_xml(xml_content: &str, pol: Polarization) -> SarResult<SarMetadata> {
        // This is a simplified parser. In practice, you'd need more robust XML parsing
        // for the complex Sentinel-1 annotation format
        
        // Extract basic information using string parsing (temporary approach)
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

        // Create minimal metadata structure
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
            sub_swaths: HashMap::new(),
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
        use std::io::Write;
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
        use std::io::Write;
        use tempfile::NamedTempFile;
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
