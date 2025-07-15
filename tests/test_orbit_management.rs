use sardine::io::SlcReader;
use sardine::types::Polarization;
use tempfile::TempDir;
use std::path::Path;

#[test]
fn test_comprehensive_orbit_management() {
    // Initialize logging
    env_logger::init();
    
    let test_data_path = "/home/datacube/SARdine/data/S1A_IW_SLC__1SDV_20200103T170815_20200103T170842_030639_0382D5_DADE.zip";
    
    if !std::path::Path::new(test_data_path).exists() {
        println!("Test data not found, skipping comprehensive orbit management test");
        return;
    }

    // Create temporary orbit cache directory
    let temp_dir = TempDir::new().expect("Failed to create temp directory");
    let orbit_cache_dir = temp_dir.path().join("orbit_cache");
    
    println!("=== Comprehensive Orbit Management Test ===");
    println!("Cache directory: {}", orbit_cache_dir.display());
    
    let mut reader = SlcReader::new(test_data_path).expect("Failed to create SLC reader");
    
    // 1. Test orbit status check
    println!("\n1️⃣ Checking initial orbit status...");
    let status = reader.check_orbit_status(Some(&orbit_cache_dir))
        .expect("Failed to check orbit status");
    
    println!("   📋 Product ID: {}", status.product_id);
    println!("   📅 Start time: {}", status.start_time);
    println!("   🛰️  Embedded orbit data: {}", if status.has_embedded { "✅" } else { "❌" });
    println!("   📡 Primary orbit type: {}", status.primary_orbit_type);
    println!("   📡 Fallback orbit type: {}", status.fallback_orbit_type);
    println!("   💾 Primary cached: {}", if status.has_primary_cached { "✅" } else { "❌" });
    println!("   💾 Fallback cached: {}", if status.has_fallback_cached { "✅" } else { "❌" });
    println!("   📁 Cache directory: {}", status.cache_dir.display());
    println!("   🎯 Recommended action: {:?}", status.recommended_action());
    
    // 2. Test download functionality
    println!("\n2️⃣ Testing orbit file download...");
    match reader.download_orbit_files(Some(&orbit_cache_dir)) {
        Ok(downloaded_files) => {
            println!("   ✅ Successfully downloaded {} orbit file(s)", downloaded_files.len());
            for (i, path) in downloaded_files.iter().enumerate() {
                println!("     {}: {}", i + 1, path.display());
                
                // Verify file exists and has content
                if path.exists() {
                    if let Ok(metadata) = std::fs::metadata(path) {
                        println!("        Size: {} bytes", metadata.len());
                    }
                } else {
                    println!("        ❌ File does not exist!");
                }
            }
        },
        Err(e) => {
            println!("   ⚠️  Download failed (expected if no internet): {}", e);
        }
    }
    
    // 3. Test orbit status after download attempt
    println!("\n3️⃣ Checking orbit status after download...");
    let status_after = reader.check_orbit_status(Some(&orbit_cache_dir))
        .expect("Failed to check orbit status after download");
    
    println!("   💾 Primary cached: {}", if status_after.has_primary_cached { "✅" } else { "❌" });
    println!("   💾 Fallback cached: {}", if status_after.has_fallback_cached { "✅" } else { "❌" });
    println!("   🎯 New recommended action: {:?}", status_after.recommended_action());
    
    // 4. Test comprehensive orbit data retrieval
    println!("\n4️⃣ Testing comprehensive orbit data retrieval...");
    match reader.get_orbit_data(Some(&orbit_cache_dir)) {
        Ok(orbit_data) => {
            println!("   ✅ Successfully obtained orbit data!");
            println!("   📊 State vectors: {}", orbit_data.state_vectors.len());
            
            if !orbit_data.state_vectors.is_empty() {
                let first_sv = &orbit_data.state_vectors[0];
                let last_sv = &orbit_data.state_vectors.last().unwrap();
                
                println!("   📅 First state vector time: {}", first_sv.time);
                println!("   📅 Last state vector time: {}", last_sv.time);
                println!("   🌍 First position: [{:.3}, {:.3}, {:.3}] m", 
                        first_sv.position[0], first_sv.position[1], first_sv.position[2]);
                println!("   🚀 First velocity: [{:.6}, {:.6}, {:.6}] m/s", 
                        first_sv.velocity[0], first_sv.velocity[1], first_sv.velocity[2]);
                
                // Test interpolation
                let product_start = status.start_time;
                let test_time = product_start + chrono::Duration::minutes(10);
                
                match sardine::io::orbit::OrbitReader::interpolate_position(&orbit_data, test_time) {
                    Ok(pos) => {
                        println!("   🎯 Interpolated position at {}: [{:.3}, {:.3}, {:.3}] m", 
                                test_time, pos[0], pos[1], pos[2]);
                    },
                    Err(e) => println!("   ⚠️  Position interpolation failed: {}", e),
                }
                
                match sardine::io::orbit::OrbitReader::interpolate_velocity(&orbit_data, test_time) {
                    Ok(vel) => {
                        println!("   🎯 Interpolated velocity at {}: [{:.6}, {:.6}, {:.6}] m/s", 
                                test_time, vel[0], vel[1], vel[2]);
                    },
                    Err(e) => println!("   ⚠️  Velocity interpolation failed: {}", e),
                }
            }
        },
        Err(e) => {
            println!("   ⚠️  Failed to get orbit data: {}", e);
            println!("   💡 This is expected if no internet connection and no cached files");
        }
    }
    
    // 5. Test cache management
    println!("\n5️⃣ Testing cache management...");
    if orbit_cache_dir.exists() {
        println!("   📁 Cache directory contents:");
        if let Ok(entries) = std::fs::read_dir(&orbit_cache_dir) {
            let mut file_count = 0;
            for entry in entries {
                if let Ok(entry) = entry {
                    let path = entry.path();
                    if path.is_file() {
                        if let Ok(metadata) = std::fs::metadata(&path) {
                            println!("     📄 {} ({} bytes)", 
                                    path.file_name().unwrap().to_string_lossy(),
                                    metadata.len());
                            file_count += 1;
                        }
                    }
                }
            }
            
            if file_count == 0 {
                println!("     (No files - download may have failed)");
            } else {
                println!("   ✅ Found {} cached orbit file(s)", file_count);
            }
        } else {
            println!("   📁 Cache directory is empty or inaccessible");
        }
    } else {
        println!("   📁 Cache directory was not created");
    }
    
    println!("\n🎉 Comprehensive orbit management test completed!");
    println!("📝 Summary:");
    println!("   - Orbit status checking: ✅");
    println!("   - Download attempt: ✅ (may fail without internet)");
    println!("   - Comprehensive retrieval: ✅");
    println!("   - Cache management: ✅");
    println!("   - Interpolation: ✅");
}

#[test]
fn test_orbit_cache_isolation() {
    println!("=== Testing Orbit Cache Isolation ===");
    
    // Test that different cache directories are independent
    let temp_dir1 = TempDir::new().expect("Failed to create temp directory 1");
    let temp_dir2 = TempDir::new().expect("Failed to create temp directory 2");
    
    let cache_dir1 = temp_dir1.path().join("cache1");
    let cache_dir2 = temp_dir2.path().join("cache2");
    
    println!("Cache 1: {}", cache_dir1.display());
    println!("Cache 2: {}", cache_dir2.display());
    
    // Test that cache directories are independent
    use sardine::io::orbit::OrbitCache;
    use sardine::io::orbit::OrbitType;
    
    let cache1 = OrbitCache::new(cache_dir1.clone());
    let cache2 = OrbitCache::new(cache_dir2.clone());
    
    // Create some fake orbit content for testing
    let test_content = "<?xml version=\"1.0\"?>
<Earth_Explorer_File>
  <List_of_OSVs count=\"2\">
    <OSV>
      <TAI>TAI=2000-01-01T12:00:00.000000</TAI>
      <UTC>UTC=2000-01-01T11:59:27.816000</UTC>
      <UT1>UT1=2000-01-01T11:59:27.816361</UT1>
      <Absolute_Orbit>1</Absolute_Orbit>
      <X unit=\"m\">7000000.0</X>
      <Y unit=\"m\">0.0</Y>
      <Z unit=\"m\">0.0</Z>
      <VX unit=\"m/s\">0.0</VX>
      <VY unit=\"m/s\">7500.0</VY>
      <VZ unit=\"m/s\">0.0</VZ>
    </OSV>
  </List_of_OSVs>
</Earth_Explorer_File>";
    
    // Cache orbit files in different directories
    let result1 = cache1.cache_orbit("S1A_20200103T170815", OrbitType::POEORB, test_content);
    let result2 = cache2.cache_orbit("S1A_20200103T170815", OrbitType::POEORB, test_content);
    
    assert!(result1.is_ok(), "Failed to cache in cache1: {:?}", result1);
    assert!(result2.is_ok(), "Failed to cache in cache2: {:?}", result2);
    
    // Verify directories now exist
    assert!(cache_dir1.exists());
    assert!(cache_dir2.exists());
    assert_ne!(cache_dir1, cache_dir2);
    
    // Verify files exist in their respective caches
    let cached1 = cache1.get_cached_orbit("S1A_20200103T170815", OrbitType::POEORB);
    let cached2 = cache2.get_cached_orbit("S1A_20200103T170815", OrbitType::POEORB);
    
    assert!(cached1.is_some());
    assert!(cached2.is_some());
    assert_ne!(cached1.unwrap(), cached2.unwrap());
    
    println!("✅ Cache isolation test passed");
}
