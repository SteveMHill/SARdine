use sardine::io::SlcReader;
use sardine::types::Polarization;
use std::time::Instant;

#[test]
fn test_slc_read_performance() {
    // Initialize logging to see performance metrics
    env_logger::init();
    
    let test_data_path = "/home/datacube/SARdine/data/S1A_IW_SLC__1SDV_20200103T170815_20200103T170842_030639_0382D5_DADE.zip";
    
    // Skip test if file doesn't exist
    if !std::path::Path::new(test_data_path).exists() {
        println!("Test data not found, skipping performance test");
        return;
    }

    let mut reader = SlcReader::new(test_data_path).expect("Failed to create SLC reader");
    
    // Find available polarizations
    let annotations = reader.find_annotation_files().expect("Failed to find annotations");
    println!("Available polarizations: {:?}", annotations.keys().collect::<Vec<_>>());
    
    // Test reading performance for VV polarization if available
    if annotations.contains_key(&Polarization::VV) {
        println!("\n=== Testing VV Polarization Read Performance ===");
        
        let start_time = Instant::now();
        let slc_data = reader.read_slc_data(Polarization::VV).expect("Failed to read SLC data");
        let total_time = start_time.elapsed();
        
        let (height, width) = slc_data.dim();
        let total_pixels = height * width;
        let data_size_mb = (total_pixels * 8) as f64 / (1024.0 * 1024.0); // 8 bytes per complex pixel
        let throughput_mb_s = data_size_mb / total_time.as_secs_f64();
        let pixels_per_sec = total_pixels as f64 / total_time.as_secs_f64();
        
        println!("Performance Results:");
        println!("  - Image size: {} x {} = {} pixels", width, height, total_pixels);
        println!("  - Data size: {:.2} MB", data_size_mb);
        println!("  - Total time: {:.3} seconds", total_time.as_secs_f64());
        println!("  - Throughput: {:.1} MB/s", throughput_mb_s);
        println!("  - Pixel rate: {:.0} pixels/s", pixels_per_sec);
        
        // Check some basic properties of the data
        println!("\nData validation:");
        println!("  - Non-zero pixels: {}", slc_data.iter().filter(|&&x| x.norm() > 0.0).count());
        println!("  - Max magnitude: {:.2}", slc_data.iter().map(|x| x.norm()).fold(0.0f32, f32::max));
        println!("  - Average magnitude: {:.2}", slc_data.iter().map(|x| x.norm()).sum::<f32>() / total_pixels as f32);
        
        // Performance assessment
        if throughput_mb_s > 100.0 {
            println!("✅ EXCELLENT performance (>100 MB/s)");
        } else if throughput_mb_s > 50.0 {
            println!("✅ GOOD performance (>50 MB/s)");
        } else if throughput_mb_s > 20.0 {
            println!("⚠️  ACCEPTABLE performance (>20 MB/s)");
        } else {
            println!("❌ SLOW performance (<20 MB/s) - needs optimization");
        }
    }

    // Test reading multiple polarizations if available
    if annotations.contains_key(&Polarization::VH) {
        println!("\n=== Testing VH Polarization Read Performance ===");
        
        let start_time = Instant::now();
        let _slc_data = reader.read_slc_data(Polarization::VH).expect("Failed to read VH SLC data");
        let total_time = start_time.elapsed();
        
        println!("VH read time: {:.3} seconds", total_time.as_secs_f64());
    }
}

#[test]
fn test_metadata_read_performance() {
    let test_data_path = "/home/datacube/SARdine/data/S1A_IW_SLC__1SDV_20200103T170815_20200103T170842_030639_0382D5_DADE.zip";
    
    if !std::path::Path::new(test_data_path).exists() {
        println!("Test data not found, skipping metadata performance test");
        return;
    }

    let mut reader = SlcReader::new(test_data_path).expect("Failed to create SLC reader");
    
    println!("\n=== Testing Metadata Read Performance ===");
    
    // Time file listing
    let start_time = Instant::now();
    let files = reader.list_files().expect("Failed to list files");
    let list_time = start_time.elapsed();
    
    // Time annotation finding
    let start_time = Instant::now();
    let annotations = reader.find_annotation_files().expect("Failed to find annotations");
    let annotation_time = start_time.elapsed();
    
    // Time metadata reading
    let start_time = Instant::now();
    if annotations.contains_key(&Polarization::VV) {
        let _metadata = reader.read_annotation(Polarization::VV).expect("Failed to read metadata");
    }
    let metadata_time = start_time.elapsed();
    
    println!("Metadata Performance:");
    println!("  - File listing ({} files): {:.3} ms", files.len(), list_time.as_millis());
    println!("  - Annotation finding: {:.3} ms", annotation_time.as_millis());
    println!("  - Metadata parsing: {:.3} ms", metadata_time.as_millis());
    
    let total_metadata_time = list_time + annotation_time + metadata_time;
    println!("  - Total metadata time: {:.3} ms", total_metadata_time.as_millis());
    
    if total_metadata_time.as_millis() < 100 {
        println!("✅ EXCELLENT metadata performance (<100ms)");
    } else if total_metadata_time.as_millis() < 500 {
        println!("✅ GOOD metadata performance (<500ms)");
    } else {
        println!("⚠️  SLOW metadata performance (>500ms)");
    }
}
