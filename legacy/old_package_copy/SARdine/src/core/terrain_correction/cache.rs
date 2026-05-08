use std::collections::HashMap;
use std::sync::{Arc, Mutex};

use ndarray::Array2;

/// Memory pool for efficient allocation of working arrays
/// Reduces garbage collection overhead and improves cache locality
#[derive(Debug)]
pub struct MemoryPool {
    /// Pool of reusable 2D float arrays for intermediate computations
    float_arrays_2d: Arc<Mutex<Vec<Array2<f32>>>>,
    /// Pool of reusable coordinate buffers for SIMD operations
    coord_buffers: Arc<Mutex<Vec<Vec<(f32, f32)>>>>,
    /// Statistics for pool usage optimization
    allocation_stats: Arc<Mutex<MemoryPoolStats>>,
}

#[derive(Debug, Default)]
struct MemoryPoolStats {
    total_requests: u64,
    cache_hits: u64,
    cache_misses: u64,
    peak_pool_size: usize,
}

impl MemoryPool {
    pub fn new() -> Self {
        Self {
            float_arrays_2d: Arc::new(Mutex::new(Vec::with_capacity(16))),
            coord_buffers: Arc::new(Mutex::new(Vec::with_capacity(32))),
            allocation_stats: Arc::new(Mutex::new(MemoryPoolStats::default())),
        }
    }

    /// Get a reusable 2D float array, creating new if none available
    pub fn get_float_array_2d(&self, height: usize, width: usize) -> Array2<f32> {
        // SAFETY: Use unwrap_or_else to recover from poisoned mutex - caching can proceed
        // even if a previous thread panicked while holding the lock
        let mut stats = self
            .allocation_stats
            .lock()
            .unwrap_or_else(|poisoned| poisoned.into_inner());
        stats.total_requests += 1;

        let mut pool = self
            .float_arrays_2d
            .lock()
            .unwrap_or_else(|poisoned| poisoned.into_inner());

        // Try to find a suitable array from the pool
        for i in 0..pool.len() {
            if pool[i].nrows() == height && pool[i].ncols() == width {
                stats.cache_hits += 1;
                let mut array = pool.swap_remove(i);
                array.fill(f32::NAN); // Reset for reuse
                return array;
            }
        }

        // Create new array if no suitable one found
        stats.cache_misses += 1;
        Array2::<f32>::from_elem((height, width), f32::NAN)
    }

    /// Return a 2D float array to the pool for reuse
    pub fn return_float_array_2d(&self, array: Array2<f32>) {
        let mut pool = self
            .float_arrays_2d
            .lock()
            .unwrap_or_else(|poisoned| poisoned.into_inner());

        // Limit pool size to prevent unbounded growth
        if pool.len() < 32 {
            pool.push(array);

            let mut stats = self
                .allocation_stats
                .lock()
                .unwrap_or_else(|poisoned| poisoned.into_inner());
            stats.peak_pool_size = stats.peak_pool_size.max(pool.len());
        }
    }

    /// Get a reusable coordinate buffer for SIMD operations
    pub fn get_coord_buffer(&self, capacity: usize) -> Vec<(f32, f32)> {
        let mut stats = self
            .allocation_stats
            .lock()
            .unwrap_or_else(|poisoned| poisoned.into_inner());
        stats.total_requests += 1;

        let mut pool = self
            .coord_buffers
            .lock()
            .unwrap_or_else(|poisoned| poisoned.into_inner());

        for i in 0..pool.len() {
            if pool[i].capacity() >= capacity {
                stats.cache_hits += 1;
                let mut buffer = pool.swap_remove(i);
                buffer.clear();
                return buffer;
            }
        }

        stats.cache_misses += 1;
        Vec::with_capacity(capacity)
    }

    /// Return a coordinate buffer to the pool
    pub fn return_coord_buffer(&self, buffer: Vec<(f32, f32)>) {
        let mut pool = self
            .coord_buffers
            .lock()
            .unwrap_or_else(|poisoned| poisoned.into_inner());
        if pool.len() < 64 {
            pool.push(buffer);
        }
    }

    /// Print memory pool statistics for optimization
    pub fn print_stats(&self) {
        let stats = self
            .allocation_stats
            .lock()
            .unwrap_or_else(|poisoned| poisoned.into_inner());
        let hit_rate = if stats.total_requests > 0 {
            100.0 * stats.cache_hits as f64 / stats.total_requests as f64
        } else {
            0.0
        };

        log::info!("🧠 Memory Pool Stats:");
        log::info!("   📊 Total requests: {}", stats.total_requests);
        log::info!("   ✅ Cache hit rate: {:.1}%", hit_rate);
        log::info!("   📈 Peak pool size: {}", stats.peak_pool_size);
    }
}

/// Cache-friendly data layout for optimized memory access
/// Groups related data together to improve CPU cache performance
#[derive(Debug)]
pub struct CacheFriendlyLUT {
    /// Interleaved range and azimuth values for better cache locality
    /// Layout: [range0, azimuth0, range1, azimuth1, ...]
    pub interleaved_data: Vec<f32>,
    /// Validity flags packed for efficient access
    pub validity_mask: Vec<bool>,
    /// Grid dimensions
    pub height: usize,
    pub width: usize,
    /// Base grid spacing for coordinate mapping
    pub grid_spacing: f32,
}

impl CacheFriendlyLUT {
    pub fn new(
        range_lut: &Array2<f32>,
        azimuth_lut: &Array2<f32>,
        valid_lut: &Array2<bool>,
        grid_spacing: f32,
    ) -> Self {
        let height = range_lut.nrows();
        let width = range_lut.ncols();
        let total_pixels = height * width;

        let mut interleaved_data = Vec::with_capacity(total_pixels * 2);
        let mut validity_mask = Vec::with_capacity(total_pixels);

        // Interleave range and azimuth data for better cache locality
        for i in 0..height {
            for j in 0..width {
                interleaved_data.push(range_lut[[i, j]]);
                interleaved_data.push(azimuth_lut[[i, j]]);
                validity_mask.push(valid_lut[[i, j]]);
            }
        }

        Self {
            interleaved_data,
            validity_mask,
            height,
            width,
            grid_spacing,
        }
    }

    /// Get range and azimuth values at grid position (i, j)
    #[inline(always)]
    pub fn get_values(&self, i: usize, j: usize) -> (f32, f32) {
        let index = (i * self.width + j) * 2;
        (
            self.interleaved_data[index],
            self.interleaved_data[index + 1],
        )
    }

    /// Check if position (i, j) has valid data
    #[inline(always)]
    pub fn is_valid(&self, i: usize, j: usize) -> bool {
        let index = i * self.width + j;
        self.validity_mask[index]
    }
}

/// GPU computation context for OpenCL/CUDA acceleration
/// Prepared for future GPU acceleration implementation
#[derive(Debug)]
pub struct GPUContext {
    /// Device information for optimal kernel selection
    pub device_type: String,
    /// Available memory for buffer allocation
    pub available_memory: usize,
    /// Optimal work group size for kernels
    pub work_group_size: usize,
    /// Whether double precision is supported
    pub supports_double_precision: bool,
    /// Enable GPU acceleration if available
    pub enabled: bool,
}

impl Default for GPUContext {
    fn default() -> Self {
        Self {
            device_type: "CPU".to_string(),
            available_memory: 0,
            work_group_size: 1,
            supports_double_precision: true,
            enabled: false, // Disabled by default until implementation complete
        }
    }
}

impl GPUContext {
    /// Initialize GPU context if CUDA/OpenCL is available
    pub fn try_initialize() -> Self {
        // Note: Implement GPU detection and initialization
        // For now, return CPU fallback
        log::info!("🔧 GPU acceleration not yet implemented, using CPU");
        Self::default()
    }

    /// Check if problem size is suitable for GPU acceleration
    pub fn should_use_gpu(&self, total_pixels: usize) -> bool {
        self.enabled && total_pixels > 1_000_000 // Only for large problems
    }
}

/// DEM cache for fast elevation lookups
#[derive(Debug, Clone)]
pub struct DemCache {
    /// Cached DEM values for frequently accessed coordinates
    cache: HashMap<(i32, i32), f32>,
    /// Maximum cache size
    max_size: usize,
}

impl DemCache {
    pub fn new(max_size: usize) -> Self {
        Self {
            cache: HashMap::new(),
            max_size,
        }
    }

    pub fn get(&self, x: i32, y: i32) -> Option<f32> {
        self.cache.get(&(x, y)).copied()
    }

    pub fn insert(&mut self, x: i32, y: i32, value: f32) {
        if self.cache.len() >= self.max_size {
            // Simple eviction: clear half the cache
            let keys_to_remove: Vec<_> = self
                .cache
                .keys()
                .take(self.cache.len() / 2)
                .cloned()
                .collect();
            for key in keys_to_remove {
                self.cache.remove(&key);
            }
        }
        self.cache.insert((x, y), value);
    }
}

// NOTE: DEMTileCache was removed in Dec 2025 cleanup.
// The disabled implementation caused confusion and was never used.
// If DEM tile caching is needed in the future, implement it from scratch
// with proper GDAL thread-safety handling (serial loading or pre-loaded arrays).
