use sardine::io::SlcReader;
use sardine::types::Polarization;
use std::time::Instant;

#[test]
fn test_parallel_vs_sequential_performance() {
    // Initialize logging to see performance metrics
    env_logger::init();
    
    let test_data_path = "/home/datacube/SARdine/data/S1A_IW_SLC__1SDV_20200103T170815_20200103T170842_030639_0382D5_DADE.zip";
    
    // Skip test if file doesn't exist
    if !std::path::Path::new(test_data_path).exists() {
        println!("Test data not found, skipping parallel performance test");
        return;
    }

    let mut reader = SlcReader::new(test_data_path).expect("Failed to create SLC reader");
    
    // Find available polarizations
    let annotations = reader.find_annotation_files().expect("Failed to find annotations");
    
    if !annotations.contains_key(&Polarization::VV) {
        println!("VV polarization not found, skipping test");
        return;
    }

    println!("\n=== Performance Comparison: Sequential vs Parallel ===");
    
    // Test sequential reading
    println!("\n--- Sequential Reading ---");
    let sequential_start = Instant::now();
    let sequential_data = reader.read_slc_data(Polarization::VV).expect("Failed to read SLC data sequentially");
    let sequential_time = sequential_start.elapsed();
    
    let (height, width) = sequential_data.dim();
    let total_pixels = height * width;
    let data_size_mb = (total_pixels * 8) as f64 / (1024.0 * 1024.0);
    let sequential_throughput = data_size_mb / sequential_time.as_secs_f64();
    
    println!("Sequential performance:");
    println!("  - Time: {:.3} seconds", sequential_time.as_secs_f64());
    println!("  - Throughput: {:.1} MB/s", sequential_throughput);
    
    // Test parallel reading
    println!("\n--- Parallel Reading ---");
    let parallel_start = Instant::now();
    let parallel_data = reader.read_slc_data_parallel(Polarization::VV).expect("Failed to read SLC data in parallel");
    let parallel_time = parallel_start.elapsed();
    
    let parallel_throughput = data_size_mb / parallel_time.as_secs_f64();
    
    println!("Parallel performance:");
    println!("  - Time: {:.3} seconds", parallel_time.as_secs_f64());
    println!("  - Throughput: {:.1} MB/s", parallel_throughput);
    
    // Compare results
    println!("\n--- Performance Comparison ---");
    let speedup = sequential_time.as_secs_f64() / parallel_time.as_secs_f64();
    let throughput_improvement = (parallel_throughput / sequential_throughput - 1.0) * 100.0;
    
    println!("  - Speedup: {:.2}x", speedup);
    println!("  - Throughput improvement: {:.1}%", throughput_improvement);
    
    // Verify data consistency (basic check)
    let mut pixel_diffs = 0;
    let mut max_diff = 0.0f32;
    
    for row in 0..height.min(100) { // Check first 100 rows for speed
        for col in 0..width.min(100) {
            let seq_val = sequential_data[[row, col]];
            let par_val = parallel_data[[row, col]];
            let diff = (seq_val - par_val).norm();
            if diff > 1e-6 {
                pixel_diffs += 1;
                max_diff = max_diff.max(diff);
            }
        }
    }
    
    println!("\n--- Data Consistency Check ---");
    println!("  - Pixel differences in sample (100x100): {}", pixel_diffs);
    println!("  - Maximum difference: {:.2e}", max_diff);
    
    if pixel_diffs == 0 {
        println!("✅ Data consistency verified");
    } else if max_diff < 1e-3 {
        println!("⚠️  Minor differences (likely floating point precision)");
    } else {
        println!("❌ Significant differences detected");
    }
    
    // Performance assessment
    if speedup > 1.5 {
        println!("✅ Significant speedup achieved with parallel processing");
    } else if speedup > 1.1 {
        println!("✅ Moderate speedup with parallel processing");
    } else if speedup < 0.9 {
        println!("❌ Parallel processing is slower (overhead dominates)");
    } else {
        println!("⚠️  Minimal difference between sequential and parallel");
    }
}
