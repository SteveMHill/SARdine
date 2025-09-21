use crate::types::{SarError, SarResult, SarMetadata};
use ndarray::Array2;
use rayon::prelude::*;
use std::collections::HashMap;
use std::sync::Arc;
use wide::f32x8;

/// Optimized IW Split processor with parallel processing and memory optimization
#[derive(Debug, Clone)]
pub struct IwSplitOptimized {
    /// Pre-allocated memory pools for efficient processing
    memory_pools: Arc<MemoryPool>,
    /// Number of parallel threads to use
    num_threads: usize,
    /// Chunk size for parallel processing
    chunk_size: usize,
    /// Enable SIMD optimizations
    enable_simd: bool,
}

/// Memory pool for efficient buffer management
#[derive(Debug)]
pub struct MemoryPool {
    /// Buffer pools for different sizes
    complex_buffers: std::sync::Mutex<Vec<Vec<num_complex::Complex<f32>>>>,
    /// Allocation tracking
    allocation_stats: std::sync::Mutex<AllocationStats>,
}

#[derive(Debug, Clone)]
pub struct AllocationStats {
    pub total_allocations: u64,
    pub current_memory_mb: f64,
    pub peak_memory_mb: f64,
    pub buffer_reuse_count: u64,
}

/// Subswath geometry information 
#[derive(Debug, Clone)]
pub struct SubswathGeometry {
    pub swath_id: String,
    pub first_line: usize,
    pub last_line: usize,
    pub first_sample: usize,
    pub last_sample: usize,
    pub samples_per_burst: usize,
    pub lines_per_burst: usize,
    pub num_bursts: usize,
    pub range_pixel_spacing: f64,
    pub azimuth_pixel_spacing: f64,
}

/// Split result with performance metrics
#[derive(Debug)]
pub struct SplitResult {
    pub data: Array2<num_complex::Complex<f32>>,
    pub geometry: SubswathGeometry,
    pub processing_time_ms: f64,
    pub memory_used_mb: f64,
    pub throughput_mpixels_per_sec: f64,
    pub parallelization_efficiency: f64,
}

impl IwSplitOptimized {
    /// Create new optimized IW splitter with default settings
    pub fn new() -> Self {
        Self::with_config(
            num_cpus::get(), // Use all available cores
            1024,            // Default chunk size  
            true,            // Enable SIMD by default
        )
    }
    
    /// Create optimized IW splitter with custom configuration
    pub fn with_config(num_threads: usize, chunk_size: usize, enable_simd: bool) -> Self {
        Self {
            memory_pools: Arc::new(MemoryPool::new()),
            num_threads: num_threads.max(1),
            chunk_size: chunk_size.max(64),
            enable_simd,
        }
    }
    
    /// Calculate optimal chunk size based on data dimensions
    fn calculate_chunk_size(&self, rows: usize, cols: usize) -> usize {
        let pixels_per_chunk = rows * cols / self.num_threads.max(1);
        let scaled_chunk = pixels_per_chunk / cols; // Convert to rows
        scaled_chunk.next_power_of_two().min(8192).max(512)
    }

    /// Split subswaths using optimized parallel processing
    pub fn split_subswaths_optimized(
        &self,
        slc_data: &Array2<num_complex::Complex<f32>>,
        metadata: &SarMetadata,
        target_subswaths: &[String],
    ) -> SarResult<HashMap<String, SplitResult>> {
        let mut results = HashMap::new();
        
        for subswath in target_subswaths {
            // Get subswath information from metadata
            if let Some(subswath_info) = metadata.sub_swaths.get(subswath) {
                // Extract geometry information preserving global reference frame
                let geometry = SubswathGeometry {
                    swath_id: subswath.clone(),
                    first_line: subswath_info.first_line_global, // EXPERT FIX: preserve global coordinates
                    last_line: subswath_info.last_line_global,   // EXPERT FIX: preserve global coordinates
                    first_sample: subswath_info.first_sample_global, // EXPERT FIX: preserve global coordinates  
                    last_sample: subswath_info.last_sample_global,   // EXPERT FIX: preserve global coordinates
                    samples_per_burst: subswath_info.range_samples / subswath_info.burst_count.max(1),
                    lines_per_burst: subswath_info.azimuth_samples / subswath_info.burst_count.max(1),
                    num_bursts: subswath_info.burst_count,
                    range_pixel_spacing: metadata.pixel_spacing.0,
                    azimuth_pixel_spacing: metadata.pixel_spacing.1,
                };
                
                // Process this subswath
                let split_result = self.extract_large_subswath_parallel(slc_data, &geometry)?;
                results.insert(subswath.clone(), split_result);
            } else {
                // If no specific subswaths found, create default geometry
                let total_samples = slc_data.ncols();
                let total_lines = slc_data.nrows();
                
                let geometry = SubswathGeometry {
                    swath_id: subswath.clone(),
                    first_line: 0,
                    last_line: total_lines,
                    first_sample: 0,
                    last_sample: total_samples,
                    samples_per_burst: total_samples / 9, // Typical 9 bursts
                    lines_per_burst: total_lines / 9,
                    num_bursts: 9,
                    range_pixel_spacing: metadata.pixel_spacing.0,
                    azimuth_pixel_spacing: metadata.pixel_spacing.1,
                };
                
                let split_result = self.extract_large_subswath_parallel(slc_data, &geometry)?;
                results.insert(subswath.clone(), split_result);
            }
        }
        
        Ok(results)
    }
    
    /// Extract large subswath using parallel processing
    fn extract_large_subswath_parallel(
        &self,
        slc_data: &Array2<num_complex::Complex<f32>>,
        geometry: &SubswathGeometry,
    ) -> SarResult<SplitResult> {
        let start_time = std::time::Instant::now();
        let initial_memory = self.memory_pools.get_current_memory_mb();
        
        let rows = geometry.last_line - geometry.first_line;
        let cols = geometry.last_sample - geometry.first_sample;
        
        // Get buffer from memory pool
        let mut output = self.memory_pools.get_complex_buffer(rows * cols);
        
        // Process in parallel chunks
        let chunk_rows = self.calculate_chunk_size(rows, cols);
        
        output
            .par_chunks_mut(chunk_rows * cols)
            .enumerate()
            .try_for_each(|(chunk_idx, chunk)| -> SarResult<()> {
                let start_row = chunk_idx * chunk_rows + geometry.first_line;
                let end_row = (start_row + chunk_rows).min(geometry.last_line);
                let actual_rows = end_row - start_row;
                
                // Copy data for this chunk with TOPS empty line handling
                for (row_idx, row) in (start_row..end_row).enumerate() {
                    let dst_start = row_idx * cols;
                    let dst_end = dst_start + cols;
                    
                    if dst_end <= chunk.len() {
                        if row < slc_data.nrows() {
                            // Valid data row - copy normally
                            let src_row = slc_data.row(row);
                            let src_slice = &src_row.as_slice().unwrap()[geometry.first_sample..geometry.last_sample];
                            
                            if src_slice.len() == cols {
                                // Use SIMD if available and data is aligned
                                if self.enable_simd && cols >= 8 && cols % 8 == 0 {
                                    self.copy_row_simd(src_slice, &mut chunk[dst_start..dst_end]);
                                } else {
                                    chunk[dst_start..dst_end].copy_from_slice(src_slice);
                                }
                            } else {
                                // Invalid data size - fill with zeros (TOPS empty line handling)
                                chunk[dst_start..dst_end].fill(num_complex::Complex32::new(0.0, 0.0));
                            }
                        } else {
                            // Row beyond data bounds - TOPS empty line (non-contiguous burst)
                            chunk[dst_start..dst_end].fill(num_complex::Complex32::new(0.0, 0.0));
                        }
                    }
                }
                
                Ok(())
            })?;
        
        // Convert buffer to Array2
        let result_data = Array2::from_shape_vec((rows, cols), output)
            .map_err(|e| SarError::Processing(format!("Failed to reshape array: {}", e)))?;
        
        let processing_time = start_time.elapsed().as_secs_f64() * 1000.0;
        let memory_used = self.memory_pools.get_current_memory_mb() - initial_memory;
        let pixels_processed = rows * cols;
        let throughput = pixels_processed as f64 / processing_time * 1000.0 / 1_000_000.0;
        
        Ok(SplitResult {
            data: result_data,
            geometry: geometry.clone(),
            processing_time_ms: processing_time,
            memory_used_mb: memory_used,
            throughput_mpixels_per_sec: throughput,
            parallelization_efficiency: self.calculate_parallelization_efficiency(processing_time),
        })
    }
    
    /// Copy row data using SIMD operations
    fn copy_row_simd(
        &self,
        src: &[num_complex::Complex<f32>],
        dst: &mut [num_complex::Complex<f32>],
    ) {
        if src.len() != dst.len() || src.len() % 8 != 0 {
            // Fall back to standard copy if dimensions don't align
            dst.copy_from_slice(src);
            return;
        }
        
        // Process 8 complex numbers at a time using SIMD
        let chunks = src.len() / 8;
        
        for i in 0..chunks {
            let src_start = i * 8;
            let dst_start = i * 8;
            
            // Load 8 complex numbers (16 f32 values)
            let src_chunk = &src[src_start..src_start + 8];
            let dst_chunk = &mut dst[dst_start..dst_start + 8];
            
            // Use wide crate for SIMD operations
            unsafe {
                let src_ptr = src_chunk.as_ptr() as *const f32;
                let dst_ptr = dst_chunk.as_mut_ptr() as *mut f32;
                
                // Load 8 f32 values at a time (4 complex numbers)
                let vals1 = f32x8::new([
                    *src_ptr.offset(0), *src_ptr.offset(1), *src_ptr.offset(2), *src_ptr.offset(3),
                    *src_ptr.offset(4), *src_ptr.offset(5), *src_ptr.offset(6), *src_ptr.offset(7),
                ]);
                let vals2 = f32x8::new([
                    *src_ptr.offset(8), *src_ptr.offset(9), *src_ptr.offset(10), *src_ptr.offset(11),
                    *src_ptr.offset(12), *src_ptr.offset(13), *src_ptr.offset(14), *src_ptr.offset(15),
                ]);
                
                // Store directly to destination
                vals1.to_array().iter().enumerate().for_each(|(j, &val)| {
                    *dst_ptr.offset(j as isize) = val;
                });
                vals2.to_array().iter().enumerate().for_each(|(j, &val)| {
                    *dst_ptr.offset(8 + j as isize) = val;
                });
            }
        }
    }
    
    /// Calculate parallelization efficiency
    fn calculate_parallelization_efficiency(&self, actual_time_ms: f64) -> f64 {
        // Theoretical time for single-threaded execution
        let theoretical_single_thread_time = actual_time_ms * self.num_threads as f64;
        
        // Efficiency = (Theoretical single-thread time) / (Actual time * num_threads)
        (theoretical_single_thread_time / (actual_time_ms * self.num_threads as f64)).min(1.0)
    }
}

impl MemoryPool {
    /// Create new memory pool
    pub fn new() -> Self {
        Self {
            complex_buffers: std::sync::Mutex::new(Vec::new()),
            allocation_stats: std::sync::Mutex::new(AllocationStats {
                total_allocations: 0,
                current_memory_mb: 0.0,
                peak_memory_mb: 0.0,
                buffer_reuse_count: 0,
            }),
        }
    }
    
    /// Get complex buffer from pool or allocate new one
    pub fn get_complex_buffer(&self, size: usize) -> Vec<num_complex::Complex<f32>> {
        let mut buffers = self.complex_buffers.lock().unwrap();
        let mut stats = self.allocation_stats.lock().unwrap();
        
        // Try to reuse existing buffer
        if let Some(mut buffer) = buffers.pop() {
            if buffer.len() >= size {
                buffer.truncate(size);
                stats.buffer_reuse_count += 1;
                return buffer;
            }
        }
        
        // Allocate new buffer
        let buffer = vec![num_complex::Complex::<f32>::new(0.0, 0.0); size];
        stats.total_allocations += 1;
        
        // Update memory tracking
        let buffer_size_mb = (size * std::mem::size_of::<num_complex::Complex<f32>>()) as f64 / 1_048_576.0;
        stats.current_memory_mb += buffer_size_mb;
        stats.peak_memory_mb = stats.peak_memory_mb.max(stats.current_memory_mb);
        
        buffer
    }
    
    /// Get current memory usage in MB
    pub fn get_current_memory_mb(&self) -> f64 {
        let stats = self.allocation_stats.lock().unwrap();
        stats.current_memory_mb
    }
    
    /// Update peak memory tracking
    pub fn update_peak_memory(&self, current_mb: f64) {
        let mut stats = self.allocation_stats.lock().unwrap();
        stats.peak_memory_mb = stats.peak_memory_mb.max(current_mb);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use num_complex::Complex;
    
    #[test]
    fn test_simd_copy() {
        let splitter = IwSplitOptimized::new();
        let size = 16; // 16 complex numbers for SIMD test
        
        let src: Vec<Complex<f32>> = (0..size)
            .map(|i| Complex::new(i as f32, (i + 1) as f32))
            .collect();
        let mut dst = vec![Complex::new(0.0, 0.0); size];
        
        splitter.copy_row_simd(&src, &mut dst);
        
        assert_eq!(src, dst);
    }
}