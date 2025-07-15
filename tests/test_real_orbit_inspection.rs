use sardine::io::orbit::{OrbitReader, OrbitType};
use sardine::io::SlcReader;
use sardine::types::Polarization;
use chrono::{DateTime, Utc};
use std::fs;

#[test]
fn test_real_orbit_file_inspection() {
    // Initialize logging to see download progress
    env_logger::init();
    
    let test_data_path = "/home/datacube/SARdine/data/S1A_IW_SLC__1SDV_20200103T170815_20200103T170842_030639_0382D5_DADE.zip";
    
    if !std::path::Path::new(test_data_path).exists() {
        println!("Test data not found, skipping orbit file inspection test");
        return;
    }

    // Extract product metadata
    let mut reader = SlcReader::new(test_data_path).expect("Failed to create SLC reader");
    let metadata = reader.read_annotation(Polarization::VV).expect("Failed to read metadata");
    
    println!("=== Real Orbit File Inspection ===");
    println!("Product: {}", test_data_path.split('/').last().unwrap());
    println!("Start time: {}", metadata.start_time);
    
    // Create persistent directory for downloaded files
    let orbit_dir = std::path::Path::new("/home/datacube/SARdine/orbit_test");
    fs::create_dir_all(&orbit_dir).expect("Failed to create orbit directory");
    
    // Extract product ID from filename
    let product_id = test_data_path.split('/').last().unwrap()
        .replace(".zip", "")
        .replace(".SAFE", "");
    
    println!("Product ID: {}", product_id);
    
    // Try to download RESORB orbit file first (smaller and more likely to be available)
    println!("\n=== Downloading RESORB Orbit File ===");
    match OrbitReader::download_orbit_file(&product_id, metadata.start_time, OrbitType::RESORB, Some(&orbit_dir.join(format!("{}_RESORB.EOF", product_id)))) {
        Ok(orbit_data) => {
            println!("✅ Successfully downloaded and parsed RESORB orbit file!");
            println!("  - State vectors: {}", orbit_data.state_vectors.len());
            
            if !orbit_data.state_vectors.is_empty() {
                let first_sv = &orbit_data.state_vectors[0];
                println!("  - First state vector time: {}", first_sv.time);
                println!("  - Position: [{:.3}, {:.3}, {:.3}] m", 
                        first_sv.position[0], first_sv.position[1], first_sv.position[2]);
                println!("  - Velocity: [{:.6}, {:.6}, {:.6}] m/s", 
                        first_sv.velocity[0], first_sv.velocity[1], first_sv.velocity[2]);
            }
        },
        Err(e) => {
            println!("❌ RESORB download failed: {}", e);
            
            // Check if file was downloaded but parsing failed
            let orbit_file_path = orbit_dir.join(format!("{}_RESORB.EOF", product_id));
            if orbit_file_path.exists() {
                println!("File was downloaded, inspecting format...");
                inspect_orbit_file_format(&orbit_file_path);
            }
        }
    }
    
    // List all files in the orbit directory
    println!("\n=== Files in orbit directory ===");
    if let Ok(entries) = fs::read_dir(&orbit_dir) {
        for entry in entries {
            if let Ok(entry) = entry {
                let path = entry.path();
                if let Ok(metadata) = fs::metadata(&path) {
                    println!("  ✓ {} ({} bytes)", 
                            path.file_name().unwrap().to_string_lossy(),
                            metadata.len());
                    
                    // If it's an EOF file, inspect its format
                    if path.extension().and_then(|s| s.to_str()) == Some("EOF") {
                        inspect_orbit_file_format(&path);
                    }
                }
            }
        }
    }
}

fn inspect_orbit_file_format(file_path: &std::path::Path) {
    println!("\n=== Inspecting {} ===", file_path.file_name().unwrap().to_string_lossy());
    
    match fs::read_to_string(file_path) {
        Ok(content) => {
            let lines: Vec<&str> = content.lines().collect();
            println!("File contains {} lines ({} bytes)", lines.len(), content.len());
            
            // Show first 20 lines
            println!("\nFirst 20 lines:");
            for (i, line) in lines.iter().take(20).enumerate() {
                println!("{:3}: {}", i + 1, line);
            }
            
            // Look for specific patterns
            let utc_lines: Vec<(usize, &str)> = lines.iter().enumerate()
                .filter(|(_, line)| line.contains("UTC=") || line.contains("<UTC>"))
                .take(3)
                .map(|(i, &line)| (i, line))
                .collect();
            
            if !utc_lines.is_empty() {
                println!("\nLines containing UTC:");
                for (line_num, line) in utc_lines {
                    println!("{:3}: {}", line_num + 1, line);
                }
            }
            
            let osv_lines: Vec<(usize, &str)> = lines.iter().enumerate()
                .filter(|(_, line)| line.contains("<OSV>") || line.contains("</OSV>"))
                .take(5)
                .map(|(i, &line)| (i, line))
                .collect();
            
            if !osv_lines.is_empty() {
                println!("\nOSV block markers:");
                for (line_num, line) in osv_lines {
                    println!("{:3}: {}", line_num + 1, line);
                }
            }
            
            // Try parsing manually
            println!("\nTrying to parse with current parser...");
            match OrbitReader::read_orbit_file(file_path) {
                Ok(orbit_data) => {
                    println!("✅ Parsing successful: {} state vectors", orbit_data.state_vectors.len());
                },
                Err(e) => {
                    println!("❌ Parsing failed: {}", e);
                }
            }
        },
        Err(e) => {
            println!("❌ Failed to read file: {}", e);
        }
    }
}
