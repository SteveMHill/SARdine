//! Data extraction functions for IW TOPSAR deburst processing
//!
//! This module provides functions for extracting complex SLC data from Sentinel-1
//! ZIP archives and SAFE directories. These functions handle the CInt16 format
//! used in Sentinel-1 SLC products and convert to Complex<f32> for processing.

use crate::io::SlcReader;
use crate::types::{Polarization, SarComplex, SarError, SarResult};
use ndarray::Array2;

/// Extract complex SLC data from ZIP file or SAFE directory
///
/// This function handles both ZIP archives and extracted SAFE directories,
/// routing to the appropriate extraction method based on the input path.
pub fn extract_subswath_complex_from_zip(
    zip_path: &str,
    subswath: &str,
    polarization: &str,
) -> SarResult<Array2<num_complex::Complex<f32>>> {
    log::info!(
        "📡 Extracting COMPLEX subswath {} {} from {}",
        subswath,
        polarization,
        zip_path
    );

    // Check if input is a ZIP file or SAFE directory
    let input_path = std::path::Path::new(zip_path);
    let is_zip =
        input_path.is_file() && input_path.extension() == Some(std::ffi::OsStr::new("zip"));
    let is_safe = input_path.is_dir() && (
        input_path.file_name()
            .map(|name| name.to_string_lossy().contains(".SAFE"))
            .unwrap_or_else(|| {
                log::warn!("⚠️  Could not determine filename for path validation - checking manifest.safe instead");
                false
            }) ||
        input_path.join("manifest.safe").exists()
    );

    if !is_zip && !is_safe {
        return Err(SarError::Processing(format!(
            "Input must be either a ZIP file or SAFE directory: {}",
            zip_path
        )));
    }

    if is_safe {
        // For SAFE directories, delegate to SlcReader which has full SAFE support
        log::info!("Detected SAFE directory, using SlcReader for complex data extraction");
        return extract_subswath_complex_from_safe(zip_path, subswath, polarization);
    }

    // Original ZIP file processing continues below
    // Open the ZIP file
    let file = std::fs::File::open(zip_path)
        .map_err(|e| SarError::Processing(format!("Failed to open ZIP file: {}", e)))?;

    let mut archive = zip::ZipArchive::new(file)
        .map_err(|e| SarError::Processing(format!("Failed to read ZIP archive: {}", e)))?;

    // Find the appropriate measurement file for the subswath and polarization
    // Support both Sentinel-1A (s1a) and Sentinel-1B (s1b) products
    let subswath_num = subswath.strip_prefix("IW").unwrap_or(subswath);
    let measurement_pattern_a = format!(
        "s1a-iw{}-slc-{}-",
        subswath_num.to_lowercase(),
        polarization.to_lowercase()
    );
    let measurement_pattern_b = format!(
        "s1b-iw{}-slc-{}-",
        subswath_num.to_lowercase(),
        polarization.to_lowercase()
    );

    let mut measurement_file = None;
    for i in 0..archive.len() {
        let file = archive
            .by_index(i)
            .map_err(|e| SarError::Processing(format!("Failed to read ZIP entry: {}", e)))?;

        let name = file.name();
        if (name.contains(&measurement_pattern_a) || name.contains(&measurement_pattern_b))
            && name.ends_with(".tiff")
            && name.contains("measurement/")
        {
            measurement_file = Some(i);
            log::info!("Found measurement file: {}", name);
            break;
        }
    }

    if measurement_file.is_none() {
        return Err(SarError::Processing(format!(
            "No measurement file found for subswath {} and polarization {}",
            subswath, polarization
        )));
    }

    // Extract TIFF to temporary file for GDAL reading
    let measurement_index =
        measurement_file.expect("measurement_file checked above for Some value");

    let mut zip_file = archive
        .by_index(measurement_index)
        .map_err(|e| SarError::Processing(format!("Failed to access measurement file: {}", e)))?;

    let mut temp_file = tempfile::NamedTempFile::new()
        .map_err(|e| SarError::Processing(format!("Failed to create temp file: {}", e)))?;

    std::io::copy(&mut zip_file, &mut temp_file)
        .map_err(|e| SarError::Processing(format!("Failed to extract TIFF to temp file: {}", e)))?;

    // Use GDAL to read the TIFF file
    let dataset = gdal::Dataset::open(temp_file.path())
        .map_err(|e| SarError::Processing(format!("Failed to open TIFF with GDAL: {}", e)))?;

    let raster_size = dataset.raster_size();
    let (width, height) = (raster_size.0, raster_size.1);
    let band_count = dataset.raster_count();

    log::info!(
        "TIFF dimensions: {} x {}, bands: {}",
        width,
        height,
        band_count
    );

    // Handle complex SLC data
    if band_count == 1 {
        let band = dataset
            .rasterband(1)
            .map_err(|e| SarError::Processing(format!("Failed to get band 1: {}", e)))?;

        let window = (0, 0);
        let window_size = (width, height);

        // Read complex 16-bit integers (CInt16 format)
        let complex_data = band
            .read_as::<i16>(window, window_size, (width * 2, height), None)
            .map_err(|e| {
                SarError::Processing(format!("Failed to read complex CInt16 data: {}", e))
            })?;

        // Convert to Complex<f32> array (preserve I/Q values)
        convert_cint16_to_complex(complex_data, width, height)
    } else {
        Err(SarError::Processing(format!(
            "Complex extraction not supported for {} bands",
            band_count
        )))
    }
}

/// Convert CInt16 interleaved data to Complex<f32> (preserving I/Q values)
fn convert_cint16_to_complex(
    complex_data: gdal::raster::Buffer<i16>,
    width: usize,
    height: usize,
) -> SarResult<Array2<num_complex::Complex<f32>>> {
    let mut complex_array = Array2::zeros((height, width));
    let total_pixels = width * height;

    // Data is interleaved as [real, imag, real, imag, ...]
    if complex_data.data.len() < total_pixels * 2 {
        return Err(SarError::Processing(format!(
            "Insufficient data: expected {} i16 values, got {}",
            total_pixels * 2,
            complex_data.data.len()
        )));
    }

    // Convert i16 complex data to Complex<f32> values
    for row in 0..height {
        for col in 0..width {
            let pixel_idx = row * width + col;
            let data_idx = pixel_idx * 2;

            // Extract real and imaginary parts as i16, convert to f32
            let real_i16 = complex_data.data[data_idx];
            let imag_i16 = complex_data.data[data_idx + 1];

            // Convert to f32 and create complex number
            let real_f32 = real_i16 as f32;
            let imag_f32 = imag_i16 as f32;

            complex_array[[row, col]] = SarComplex::new(real_f32, imag_f32);
        }
    }

    log::info!(
        "✅ Successfully extracted COMPLEX SLC data: {}x{} pixels",
        width,
        height
    );

    Ok(complex_array)
}

/// Extract complex SLC data (preserving I/Q values) from ZIP or SAFE directory for proper radiometric calibration
pub fn extract_subswath_complex_data(
    product_path: &str,
    subswath: &str,
    polarization: &str,
) -> SarResult<Array2<num_complex::Complex<f32>>> {
    // Delegate to the ZIP/SAFE-aware function
    extract_subswath_complex_from_zip(product_path, subswath, polarization)
}

/// Helper function to extract complex SLC data from SAFE directory
fn extract_subswath_complex_from_safe(
    safe_path: &str,
    _subswath: &str,
    polarization: &str,
) -> SarResult<Array2<num_complex::Complex<f32>>> {
    log::info!(
        "Extracting complex SLC data from SAFE directory: {}",
        safe_path
    );

    // Parse polarization
    let pol = match polarization.to_uppercase().as_str() {
        "VV" => Polarization::VV,
        "VH" => Polarization::VH,
        "HV" => Polarization::HV,
        "HH" => Polarization::HH,
        _ => {
            return Err(SarError::Processing(format!(
                "Invalid polarization: {}",
                polarization
            )))
        }
    };

    // Create SlcReader for SAFE directory
    let mut reader = SlcReader::new(safe_path)
        .map_err(|e| SarError::Processing(format!("Failed to create SLC reader: {}", e)))?;

    // Read SLC data for the specified polarization
    let sar_image = reader
        .read_slc_data(pol)
        .map_err(|e| SarError::Processing(format!("Failed to read SLC data: {}", e)))?;

    log::info!(
        "Successfully read SLC data with dimensions: {}x{}",
        sar_image.nrows(),
        sar_image.ncols()
    );

    // Convert from Complex64 to Complex<f32> and ensure it's complex data
    let complex_array = sar_image.mapv(|val| SarComplex::new(val.re as f32, val.im as f32));

    // Verify this is actually complex data
    let has_imaginary = complex_array.iter().any(|&val| val.im.abs() > 1e-10);
    if !has_imaginary {
        log::warn!("SLC data appears to have no imaginary component - may be intensity data");
    }

    log::info!(
        "✅ Successfully extracted COMPLEX SLC data from SAFE: {}x{} pixels",
        complex_array.nrows(),
        complex_array.ncols()
    );

    Ok(complex_array)
}
