use sardine::io::SlcReader;
use sardine::types::Polarization;
use std::io::Read;

#[test]
fn test_orbit_file_inspection() {
    let test_data_path = "/home/datacube/SARdine/data/S1A_IW_SLC__1SDV_20200103T170815_20200103T170842_030639_0382D5_DADE.zip";
    
    if !std::path::Path::new(test_data_path).exists() {
        println!("Test data not found, skipping orbit inspection test");
        return;
    }

    let mut reader = SlcReader::new(test_data_path).expect("Failed to create SLC reader");
    let files = reader.list_files().expect("Failed to list files");
    
    println!("=== Sentinel-1 SLC Archive File Inspection ===");
    println!("Total files: {}", files.len());
    
    // Check for orbit-related files
    println!("\n🛰️ Orbit-related files:");
    let mut orbit_files = Vec::new();
    for file in &files {
        if file.to_lowercase().contains("orbit") 
            || file.to_lowercase().contains("aux") 
            || file.ends_with(".EOF") 
            || file.ends_with(".eof") {
            orbit_files.push(file);
            println!("  ✓ {}", file);
        }
    }
    
    if orbit_files.is_empty() {
        println!("  ❌ No dedicated orbit files found (.EOF format)");
    }
    
    // Check annotation files (may contain orbit info)
    println!("\n📋 Annotation files:");
    let mut annotation_files = Vec::new();
    for file in &files {
        if file.starts_with("annotation/") && file.ends_with(".xml") && !file.contains("calibration") {
            annotation_files.push(file);
            println!("  ✓ {}", file);
        }
    }
    
    // Check if manifest contains orbit info
    println!("\n📄 Manifest and metadata files:");
    for file in &files {
        if file.contains("manifest") || file.contains("metadata") || file.ends_with("safe-structure") {
            println!("  ✓ {}", file);
        }
    }
    
    // Try to read annotation and look for orbit information
    if !annotation_files.is_empty() {
        println!("\n🔍 Checking annotation files for orbit information...");
        
        let annotations = reader.find_annotation_files().expect("Failed to find annotations");
        
        if let Some(vv_file) = annotations.get(&Polarization::VV) {
            println!("Reading annotation: {}", vv_file);
            
            // Get metadata which internally reads the XML
            let metadata = reader.read_annotation(Polarization::VV).expect("Failed to read annotation");
            println!("  ✓ Successfully parsed annotation metadata");
            println!("  📅 Start time: {}", metadata.start_time);
            println!("  📅 Stop time: {}", metadata.stop_time);
            println!("  🛰️ Mission: {}", metadata.mission);
            
            // For now, we'll check if orbit_data field exists in metadata
            if let Some(_orbit_data) = &metadata.orbit_data {
                println!("  ✅ Orbit data found in metadata!");
            } else {
                println!("  ⚠️  No orbit data in parsed metadata");
            }
        }
    }
    
    // Determine orbit file strategy
    println!("\n📊 Orbit File Strategy Assessment:");
    if !orbit_files.is_empty() {
        println!("  ✅ Strategy: Use included orbit files (.EOF format)");
        println!("  📁 Found {} orbit file(s)", orbit_files.len());
    } else {
        println!("  ⚠️  Strategy: External orbit files needed");
        println!("  🌐 Need to download POEORB/RESORB from ESA");
        println!("  📅 Product date determines availability:");
        println!("     - POEORB: Available ~20 days after acquisition (most accurate)");
        println!("     - RESORB: Available ~3 hours after acquisition (backup)");
    }
}
