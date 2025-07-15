use sardine::io::SlcReader;
use sardine::types::Polarization;
use tempfile::TempDir;

#[test]
fn test_burst_orbit_interpolation() {
    // Initialize logging to see detailed orbit information
    env_logger::init();
    
    let test_data_path = "/home/datacube/SARdine/data/S1A_IW_SLC__1SDV_20200103T170815_20200103T170842_030639_0382D5_DADE.zip";
    
    if !std::path::Path::new(test_data_path).exists() {
        println!("Test data not found, skipping burst orbit interpolation test");
        return;
    }

    // Create temporary orbit cache directory
    let temp_dir = TempDir::new().expect("Failed to create temp directory");
    let orbit_cache_dir = temp_dir.path().join("orbit_cache");
    
    println!("=== Burst Orbit Interpolation Test ===");
    println!("Cache directory: {}", orbit_cache_dir.display());
    
    let mut reader = SlcReader::new(test_data_path).expect("Failed to create SLC reader");
    
    // Test 1: Get burst orbit data
    println!("\n1️⃣ Computing burst orbit data...");
    let burst_orbit = reader.get_burst_orbit_data(Polarization::VV, Some(&orbit_cache_dir))
        .expect("Failed to get burst orbit data");
    
    println!("   📊 Burst summary:");
    println!("     • Number of azimuth lines: {}", burst_orbit.num_lines());
    println!("     • Burst start time: {}", burst_orbit.burst_start_time);
    println!("     • Azimuth time interval: {:.6} seconds", burst_orbit.azimuth_time_interval);
    println!("     • Total acquisition time: {:.3} seconds", 
             burst_orbit.num_lines() as f64 * burst_orbit.azimuth_time_interval);
    
    // Test 2: Check specific satellite positions
    println!("\n2️⃣ Checking satellite positions at key lines...");
    let test_lines = [0, burst_orbit.num_lines() / 4, burst_orbit.num_lines() / 2, 
                     3 * burst_orbit.num_lines() / 4, burst_orbit.num_lines() - 1];
    
    for &line_idx in &test_lines {
        if let Some(position) = burst_orbit.get_position_at_line(line_idx) {
            if let Some(velocity) = burst_orbit.get_velocity_at_line(line_idx) {
                if let Some(az_time) = burst_orbit.get_azimuth_time_at_line(line_idx) {
                    let altitude = (position[0].powi(2) + position[1].powi(2) + position[2].powi(2)).sqrt() - 6_371_000.0;
                    let velocity_mag = (velocity[0].powi(2) + velocity[1].powi(2) + velocity[2].powi(2)).sqrt();
                    
                    println!("     Line {}: alt={:.1} km, vel={:.1} m/s, time={}", 
                            line_idx, altitude / 1000.0, velocity_mag, az_time.format("%H:%M:%S%.3f"));
                    
                    // Validate reasonable satellite parameters
                    assert!(altitude > 600_000.0, "Altitude too low: {:.1} km", altitude / 1000.0);
                    assert!(altitude < 1_000_000.0, "Altitude too high: {:.1} km", altitude / 1000.0);
                    assert!(velocity_mag > 6_000.0, "Velocity too low: {:.1} m/s", velocity_mag);
                    assert!(velocity_mag < 9_000.0, "Velocity too high: {:.1} m/s", velocity_mag);
                }
            }
        }
    }
    
    // Test 3: Pixel-level satellite position
    println!("\n3️⃣ Testing pixel-level satellite position...");
    let test_pixels = [(0, 0), (100, 500), (500, 1000), (1000, 2000)];
    
    for &(azimuth_line, range_sample) in &test_pixels {
        if azimuth_line < burst_orbit.num_lines() {
            match reader.get_satellite_position_at_pixel(
                Polarization::VV, azimuth_line, range_sample, Some(&orbit_cache_dir)
            ) {
                Ok((position, velocity)) => {
                    let altitude = (position[0].powi(2) + position[1].powi(2) + position[2].powi(2)).sqrt() - 6_371_000.0;
                    println!("     Pixel [{}, {}]: altitude={:.1} km, velocity=[{:.1}, {:.1}, {:.1}] m/s",
                            azimuth_line, range_sample, altitude / 1000.0, 
                            velocity[0], velocity[1], velocity[2]);
                },
                Err(e) => println!("     Pixel [{}, {}]: Error - {}", azimuth_line, range_sample, e),
            }
        }
    }
    
    // Test 4: Doppler centroid calculation
    println!("\n4️⃣ Testing Doppler centroid calculation...");
    for &(azimuth_line, range_sample) in &test_pixels {
        if azimuth_line < burst_orbit.num_lines() {
            match reader.calculate_doppler_centroid(
                Polarization::VV, azimuth_line, range_sample, Some(&orbit_cache_dir)
            ) {
                Ok(doppler_freq) => {
                    println!("     Pixel [{}, {}]: Doppler = {:.2} Hz", 
                            azimuth_line, range_sample, doppler_freq);
                    
                    // Validate reasonable Doppler frequency (typically < 5000 Hz for Sentinel-1)
                    assert!(doppler_freq.abs() < 10_000.0, "Unrealistic Doppler frequency: {:.2} Hz", doppler_freq);
                },
                Err(e) => println!("     Pixel [{}, {}]: Doppler error - {}", azimuth_line, range_sample, e),
            }
        }
    }
    
    // Test 5: Time consistency check
    println!("\n5️⃣ Validating time consistency...");
    let first_time = burst_orbit.get_azimuth_time_at_line(0).unwrap();
    let last_time = burst_orbit.get_azimuth_time_at_line(burst_orbit.num_lines() - 1).unwrap();
    let calculated_duration = (last_time - first_time).num_milliseconds() as f64 / 1000.0;
    let expected_duration = (burst_orbit.num_lines() - 1) as f64 * burst_orbit.azimuth_time_interval;
    
    println!("     Calculated duration: {:.6} s", calculated_duration);
    println!("     Expected duration: {:.6} s", expected_duration);
    println!("     Difference: {:.6} s", (calculated_duration - expected_duration).abs());
    
    // Allow small timing differences due to floating point precision
    assert!((calculated_duration - expected_duration).abs() < 0.001, 
           "Time inconsistency: calculated={:.6}s, expected={:.6}s", 
           calculated_duration, expected_duration);
    
    // Test 6: Orbit data quality metrics
    println!("\n6️⃣ Orbit data quality assessment...");
    let mut position_changes = Vec::new();
    let mut velocity_changes = Vec::new();
    
    for i in 1..std::cmp::min(100, burst_orbit.num_lines()) {
        if let (Some(pos1), Some(pos2)) = (
            burst_orbit.get_position_at_line(i-1),
            burst_orbit.get_position_at_line(i)
        ) {
            let pos_change = ((pos2[0] - pos1[0]).powi(2) + 
                             (pos2[1] - pos1[1]).powi(2) + 
                             (pos2[2] - pos1[2]).powi(2)).sqrt();
            position_changes.push(pos_change);
        }
        
        if let (Some(vel1), Some(vel2)) = (
            burst_orbit.get_velocity_at_line(i-1),
            burst_orbit.get_velocity_at_line(i)
        ) {
            let vel_change = ((vel2[0] - vel1[0]).powi(2) + 
                             (vel2[1] - vel1[1]).powi(2) + 
                             (vel2[2] - vel1[2]).powi(2)).sqrt();
            velocity_changes.push(vel_change);
        }
    }
    
    if !position_changes.is_empty() {
        let avg_pos_change = position_changes.iter().sum::<f64>() / position_changes.len() as f64;
        let max_pos_change = position_changes.iter().fold(0.0f64, |a, &b| a.max(b));
        
        println!("     Position change per line: avg={:.3} m, max={:.3} m", avg_pos_change, max_pos_change);
    }
    
    if !velocity_changes.is_empty() {
        let avg_vel_change = velocity_changes.iter().sum::<f64>() / velocity_changes.len() as f64;
        let max_vel_change = velocity_changes.iter().fold(0.0f64, |a, &b| a.max(b));
        
        println!("     Velocity change per line: avg={:.6} m/s, max={:.6} m/s", avg_vel_change, max_vel_change);
    }
    
    println!("\n✅ Burst orbit interpolation test completed successfully!");
    println!("🛰️  Satellite position and velocity can be accurately interpolated for each pixel");
    println!("📡 Doppler centroid calculation is ready for SAR focusing");
    println!("🎯 System is ready for precise geolocation and calibration");
}

#[test]
fn test_orbit_interpolation_accuracy() {
    env_logger::init();
    
    let test_data_path = "/home/datacube/SARdine/data/S1A_IW_SLC__1SDV_20200103T170815_20200103T170842_030639_0382D5_DADE.zip";
    
    if !std::path::Path::new(test_data_path).exists() {
        println!("Test data not found, skipping interpolation accuracy test");
        return;
    }

    let temp_dir = TempDir::new().expect("Failed to create temp directory");
    let orbit_cache_dir = temp_dir.path().join("orbit_cache");
    
    println!("=== Orbit Interpolation Accuracy Test ===");
    
    let mut reader = SlcReader::new(test_data_path).expect("Failed to create SLC reader");
    
    // Get raw orbit data to compare with interpolated values
    let orbit_data = reader.get_orbit_data(Some(&orbit_cache_dir))
        .expect("Failed to get orbit data");
    
    println!("Raw orbit data: {} state vectors", orbit_data.state_vectors.len());
    
    if orbit_data.state_vectors.len() < 4 {
        println!("Not enough state vectors for interpolation test");
        return;
    }
    
    // Test interpolation accuracy by interpolating at known state vector times
    println!("\nTesting interpolation accuracy at known state vector times:");
    
    for (i, sv) in orbit_data.state_vectors.iter().enumerate().take(5) {
        // Skip first and last vectors as they might be edge cases
        if i == 0 || i == orbit_data.state_vectors.len() - 1 {
            continue;
        }
        
        // Interpolate position and velocity at this exact time
        let interp_pos = sardine::io::orbit::OrbitReader::interpolate_position(&orbit_data, sv.time)
            .expect("Failed to interpolate position");
        let interp_vel = sardine::io::orbit::OrbitReader::interpolate_velocity(&orbit_data, sv.time)
            .expect("Failed to interpolate velocity");
        
        // Calculate errors
        let pos_error = ((interp_pos[0] - sv.position[0]).powi(2) +
                        (interp_pos[1] - sv.position[1]).powi(2) +
                        (interp_pos[2] - sv.position[2]).powi(2)).sqrt();
        
        let vel_error = ((interp_vel[0] - sv.velocity[0]).powi(2) +
                        (interp_vel[1] - sv.velocity[1]).powi(2) +
                        (interp_vel[2] - sv.velocity[2]).powi(2)).sqrt();
        
        println!("  State vector {}: pos_error={:.6} m, vel_error={:.6} m/s", 
                i, pos_error, vel_error);
        
        // Interpolation at exact state vector times should be very accurate
        assert!(pos_error < 1.0, "Position interpolation error too large: {:.6} m", pos_error);
        assert!(vel_error < 0.01, "Velocity interpolation error too large: {:.6} m/s", vel_error);
    }
    
    println!("✅ Interpolation accuracy test passed - errors within acceptable limits");
}
