use sardine::io::orbit::{OrbitReader, OrbitType};
use sardine::io::SlcReader;
use sardine::types::Polarization;
use chrono::{DateTime, Utc};
use std::path::Path;
use tempfile::TempDir;

#[test]
fn test_orbit_download_functionality() {
    // Initialize logging to see download progress
    env_logger::init();
    
    let test_data_path = "/home/datacube/SARdine/data/S1A_IW_SLC__1SDV_20200103T170815_20200103T170842_030639_0382D5_DADE.zip";
    
    if !std::path::Path::new(test_data_path).exists() {
        println!("Test data not found, skipping orbit download test");
        return;
    }

    // Extract product metadata
    let mut reader = SlcReader::new(test_data_path).expect("Failed to create SLC reader");
    let metadata = reader.read_annotation(Polarization::VV).expect("Failed to read metadata");
    
    println!("=== Orbit Download Test ===");
    println!("Product: {}", test_data_path.split('/').last().unwrap());
    println!("Start time: {}", metadata.start_time);
    println!("Mission: {}", metadata.mission);
    
    // Create temporary directory for downloads
    let temp_dir = TempDir::new().expect("Failed to create temp directory");
    let orbit_dir = temp_dir.path().join("orbits");
    std::fs::create_dir_all(&orbit_dir).expect("Failed to create orbit directory");
    
    // Extract product ID from filename
    let product_id = test_data_path.split('/').last().unwrap()
        .replace(".zip", "")
        .replace(".SAFE", "");
    
    println!("Product ID: {}", product_id);
    
    // Test orbit type determination
    let recommended_type = OrbitReader::determine_orbit_type(metadata.start_time);
    println!("Recommended orbit type: {}", recommended_type);
    
    // Test URL generation
    println!("\n=== Testing URL Generation ===");
    let poeorb_url = OrbitReader::generate_orbit_url(&product_id, OrbitType::POEORB, metadata.start_time);
    let resorb_url = OrbitReader::generate_orbit_url(&product_id, OrbitType::RESORB, metadata.start_time);
    
    println!("POEORB URL: {}", poeorb_url);
    println!("RESORB URL: {}", resorb_url);
    
    // Test automatic orbit retrieval
    println!("\n=== Testing Automatic Orbit Retrieval ===");
    match OrbitReader::get_orbit_for_product(&product_id, metadata.start_time, Some(&orbit_dir)) {
        Ok(orbit_data) => {
            println!("✅ Successfully retrieved orbit data!");
            println!("  - State vectors: {}", orbit_data.state_vectors.len());
            
            if !orbit_data.state_vectors.is_empty() {
                let first_sv = &orbit_data.state_vectors[0];
                println!("  - First state vector time: {}", first_sv.time);
                println!("  - Position: [{:.3}, {:.3}, {:.3}] m", 
                        first_sv.position[0], first_sv.position[1], first_sv.position[2]);
                println!("  - Velocity: [{:.6}, {:.6}, {:.6}] m/s", 
                        first_sv.velocity[0], first_sv.velocity[1], first_sv.velocity[2]);
                
                // Test interpolation
                println!("\n=== Testing Orbit Interpolation ===");
                let test_time = metadata.start_time + chrono::Duration::minutes(15);
                match OrbitReader::interpolate_position(&orbit_data, test_time) {
                    Ok(pos) => {
                        println!("✅ Interpolated position at {}: [{:.3}, {:.3}, {:.3}] m", 
                                test_time, pos[0], pos[1], pos[2]);
                    },
                    Err(e) => println!("❌ Position interpolation failed: {}", e),
                }
                
                match OrbitReader::interpolate_velocity(&orbit_data, test_time) {
                    Ok(vel) => {
                        println!("✅ Interpolated velocity at {}: [{:.6}, {:.6}, {:.6}] m/s", 
                                test_time, vel[0], vel[1], vel[2]);
                    },
                    Err(e) => println!("❌ Velocity interpolation failed: {}", e),
                }
            }
        },
        Err(e) => {
            println!("⚠️  Orbit retrieval failed: {}", e);
            println!("This is expected if:");
            println!("  - No internet connection");
            println!("  - ESA servers are down");
            println!("  - Orbit files are not yet available for this date");
            println!("  - URL patterns have changed");
        }
    }
    
    // Check what files were downloaded
    println!("\n=== Downloaded Files ===");
    if let Ok(entries) = std::fs::read_dir(&orbit_dir) {
        for entry in entries {
            if let Ok(entry) = entry {
                let path = entry.path();
                if let Ok(metadata) = std::fs::metadata(&path) {
                    println!("  ✓ {} ({} bytes)", 
                            path.file_name().unwrap().to_string_lossy(),
                            metadata.len());
                }
            }
        }
    } else {
        println!("  (No files downloaded)");
    }
}

#[test]
fn test_orbit_url_generation_patterns() {
    println!("=== Testing Orbit URL Generation Patterns ===");
    
    let product_id = "S1A_IW_SLC__1SDV_20200103T170815_20200103T170842_030639_0382D5_DADE";
    let test_time = DateTime::parse_from_rfc3339("2020-01-03T17:08:15Z")
        .unwrap()
        .with_timezone(&Utc);
    
    // Test both orbit types
    for orbit_type in [OrbitType::POEORB, OrbitType::RESORB] {
        println!("\n{} URLs:", orbit_type);
        
        // Test single URL generation
        let single_url = OrbitReader::generate_orbit_url(product_id, orbit_type, test_time);
        println!("  Single URL: {}", single_url);
        
        // Test multiple URL generation (private function, so we'll test via download attempt)
        println!("  Multiple URLs would be generated for robust download attempts");
    }
    
    // Test with S1B product
    let s1b_product = "S1B_IW_SLC__1SDV_20200103T170815_20200103T170842_030639_0382D5_DADE";
    let s1b_url = OrbitReader::generate_orbit_url(s1b_product, OrbitType::POEORB, test_time);
    println!("\nS1B URL: {}", s1b_url);
}
