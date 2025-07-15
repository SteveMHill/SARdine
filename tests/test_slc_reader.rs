use sardine::io::SlcReader;
use sardine::types::Polarization;
use std::path::PathBuf;

#[test]
fn test_slc_reader_with_real_data() {
    // Use the actual Sentinel-1 data file
    let test_data_path = PathBuf::from("/home/datacube/SARdine/data/S1A_IW_SLC__1SDV_20200103T170815_20200103T170842_030639_0382D5_DADE.zip");
    
    // Skip test if file doesn't exist (for CI/CD environments)
    if !test_data_path.exists() {
        println!("Test data not found, skipping test");
        return;
    }

    let mut reader = SlcReader::new(&test_data_path).expect("Failed to create SLC reader");
    
    // Test listing files
    let files = reader.list_files().expect("Failed to list files");
    assert!(!files.is_empty(), "No files found in archive");
    
    println!("Found {} files in archive", files.len());
    for file in &files[..5] { // Print first 5 files
        println!("  - {}", file);
    }
    
    // Test finding annotation files
    let annotations = reader.find_annotation_files().expect("Failed to find annotations");
    assert!(!annotations.is_empty(), "No annotation files found");
    
    println!("Found annotation files:");
    for (pol, file) in &annotations {
        println!("  {} -> {}", pol, file);
    }
    
    // Test finding measurement files
    let measurements = reader.find_measurement_files().expect("Failed to find measurements");
    assert!(!measurements.is_empty(), "No measurement files found");
    
    println!("Found measurement files:");
    for (pol, file) in &measurements {
        println!("  {} -> {}", pol, file);
    }
    
    // Test reading annotation for VV polarization (most common)
    if annotations.contains_key(&Polarization::VV) {
        let metadata = reader.read_annotation(Polarization::VV)
            .expect("Failed to read annotation");
        
        println!("Metadata for VV polarization:");
        println!("  Product ID: {}", metadata.product_id);
        println!("  Mission: {}", metadata.mission);
        println!("  Start time: {}", metadata.start_time);
        println!("  Stop time: {}", metadata.stop_time);
        println!("  Acquisition mode: {:?}", metadata.acquisition_mode);
        println!("  Pixel spacing: {:?}", metadata.pixel_spacing);
    }
}

#[test]
fn test_slc_reader_error_handling() {
    // Test with non-existent file
    let result = SlcReader::new("nonexistent.zip");
    assert!(result.is_err());
    
    // Test with invalid ZIP file (skip on CI if /dev/null doesn't exist)
    let invalid_path = PathBuf::from("/dev/null");
    if invalid_path.exists() {
        let mut reader = SlcReader::new(&invalid_path).unwrap(); // This may succeed
        let result = reader.list_files(); // But this should fail
        assert!(result.is_err());
    }
}

#[test]
fn test_polarization_detection() {
    let mut reader = SlcReader::new("/home/datacube/SARdine/data/S1A_IW_SLC__1SDV_20200103T170815_20200103T170842_030639_0382D5_DADE.zip");
    
    if reader.is_err() {
        println!("Test data not available, skipping polarization test");
        return;
    }
    
    let mut reader = reader.unwrap();
    let files = reader.list_files().unwrap_or_default();
    
    // Check that we can identify different file types
    let mut has_annotation = false;
    let mut has_measurement = false;
    let mut has_calibration = false;
    
    for file in files {
        if file.contains("annotation/") && file.ends_with(".xml") {
            has_annotation = true;
        }
        if file.contains("measurement/") && file.ends_with(".tiff") {
            has_measurement = true;
        }
        if file.contains("annotation/calibration/") {
            has_calibration = true;
        }
    }
    
    println!("File type detection:");
    println!("  Annotation files: {}", has_annotation);
    println!("  Measurement files: {}", has_measurement);
    println!("  Calibration files: {}", has_calibration);
}
