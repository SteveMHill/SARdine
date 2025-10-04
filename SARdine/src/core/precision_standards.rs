//! # Numeric Precision & Reproducibility Standards
//!
//! Enforces production-grade numeric precision and deterministic behavior:
//! - f64 for geometry, Doppler, ECEF, Newton-Raphson calculations
//! - f32 for complex samples and image data
//! - Deterministic RNG seeds for reproducible results
//! - Stable parallel processing with deterministic ordering

use std::sync::atomic::{AtomicU64, Ordering};
use rand::{Rng, SeedableRng};
use rand_chacha::ChaCha8Rng;

/// Global seed counter for deterministic RNG initialization
static GLOBAL_SEED_COUNTER: AtomicU64 = AtomicU64::new(0x12345678_9ABCDEF0);

/// Precision standards for different calculation domains
pub struct PrecisionStandards;

impl PrecisionStandards {
    /// Convergence tolerance for Newton-Raphson iterations (f64)
    pub const NEWTON_RAPHSON_TOLERANCE: f64 = 1e-12;
    
    /// Convergence tolerance for ECEF coordinate calculations (f64)
    pub const ECEF_TOLERANCE: f64 = 1e-9;
    
    /// Convergence tolerance for Doppler calculations (f64)
    pub const DOPPLER_TOLERANCE: f64 = 1e-8;
    
    /// Relative tolerance for geometry calculations (f64)
    pub const GEOMETRY_TOLERANCE: f64 = 1e-10;
    
    /// Absolute tolerance for wavelength validation (f64)
    pub const WAVELENGTH_TOLERANCE: f64 = 1e-6;
    
    /// Tolerance for complex sample comparisons (f32)
    pub const COMPLEX_SAMPLE_TOLERANCE: f32 = 1e-6;
    
    /// Tolerance for radiometric calibration (f32)
    pub const RADIOMETRIC_TOLERANCE: f32 = 1e-5;
    
    /// Maximum relative error for energy conservation checks
    pub const ENERGY_CONSERVATION_TOLERANCE: f64 = 1e-8;
}

/// Deterministic random number generator for reproducible results
pub struct DeterministicRng {
    rng: ChaCha8Rng,
    seed: u64,
}

impl DeterministicRng {
    /// Create a new deterministic RNG with a unique seed
    pub fn new() -> Self {
        let seed = GLOBAL_SEED_COUNTER.fetch_add(1, Ordering::SeqCst);
        let rng = ChaCha8Rng::seed_from_u64(seed);
        
        Self { rng, seed }
    }
    
    /// Create a deterministic RNG with a specific seed
    pub fn with_seed(seed: u64) -> Self {
        let rng = ChaCha8Rng::seed_from_u64(seed);
        Self { rng, seed }
    }
    
    /// Get the current seed (for reproducibility logging)
    pub fn seed(&self) -> u64 {
        self.seed
    }
    
    /// Generate a random f32 in [0, 1)
    pub fn gen_f32(&mut self) -> f32 {
        self.rng.gen()
    }
    
    /// Generate a random f64 in [0, 1)
    pub fn gen_f64(&mut self) -> f64 {
        self.rng.gen()
    }
    
    /// Generate a random value in a range
    pub fn gen_range<T>(&mut self, range: std::ops::Range<T>) -> T
    where
        T: rand::distributions::uniform::SampleUniform + std::cmp::PartialOrd,
    {
        self.rng.gen_range(range)
    }
}

impl Default for DeterministicRng {
    fn default() -> Self {
        Self::new()
    }
}

/// High-precision floating-point utilities for geometry calculations
pub mod geometry_precision {
    use super::PrecisionStandards;
    
    /// High-precision dot product for f64 vectors
    pub fn dot_product_f64(a: &[f64], b: &[f64]) -> f64 {
        assert_eq!(a.len(), b.len(), "Vector lengths must match");
        
        // Use Kahan summation for maximum precision
        let mut sum = 0.0;
        let mut c = 0.0; // Compensation for lost low-order bits
        
        for (ai, bi) in a.iter().zip(b.iter()) {
            let product = ai * bi;
            let y = product - c;
            let t = sum + y;
            c = (t - sum) - y;
            sum = t;
        }
        
        sum
    }
    
    /// High-precision vector norm for f64 vectors
    pub fn vector_norm_f64(v: &[f64]) -> f64 {
        dot_product_f64(v, v).sqrt()
    }
    
    /// Check if two f64 values are approximately equal
    pub fn approx_eq_f64(a: f64, b: f64, tolerance: f64) -> bool {
        let diff = (a - b).abs();
        let max_val = a.abs().max(b.abs());
        
        if max_val < f64::EPSILON {
            diff < tolerance
        } else {
            diff / max_val < tolerance
        }
    }
    
    /// Check if a geometry calculation converged
    pub fn check_convergence(current: f64, previous: f64) -> bool {
        approx_eq_f64(current, previous, PrecisionStandards::GEOMETRY_TOLERANCE)
    }
}

/// Deterministic parallel processing utilities
pub mod deterministic_parallel {
    use rayon::prelude::*;
    use std::sync::Mutex;
    
    /// Process data in deterministic chunks with stable ordering
    pub fn process_chunks_deterministic<T, F, R>(
        data: &[T],
        chunk_size: usize,
        processor: F,
    ) -> Vec<R>
    where
        T: Sync,
        F: Fn(&[T]) -> R + Sync,
        R: Send,
    {
        let chunks: Vec<_> = data.chunks(chunk_size).collect();
        let results = Mutex::new(Vec::with_capacity(chunks.len()));
        
        // Process chunks in parallel but maintain deterministic ordering
        chunks.par_iter().enumerate().for_each(|(index, chunk)| {
            let result = processor(chunk);
            
            // Store result with index to maintain order
            let mut results_guard = results.lock().unwrap();
            if results_guard.len() <= index {
                results_guard.resize_with(index + 1, || panic!("Should not happen"));
            }
            results_guard[index] = result;
        });
        
        results.into_inner().unwrap()
    }
    
    /// Deterministic parallel reduction with Kahan summation
    pub fn parallel_sum_f64(data: &[f64]) -> f64 {
        const CHUNK_SIZE: usize = 1024;
        
        let chunk_sums: Vec<f64> = data
            .par_chunks(CHUNK_SIZE)
            .map(|chunk| {
                // Use Kahan summation within each chunk
                let mut sum = 0.0;
                let mut c = 0.0;
                
                for &value in chunk {
                    let y = value - c;
                    let t = sum + y;
                    c = (t - sum) - y;
                    sum = t;
                }
                
                sum
            })
            .collect();
        
        // Final Kahan summation of chunk sums
        let mut sum = 0.0;
        let mut c = 0.0;
        
        for chunk_sum in chunk_sums {
            let y = chunk_sum - c;
            let t = sum + y;
            c = (t - sum) - y;
            sum = t;
        }
        
        sum
    }
}

/// Precision validation utilities
pub mod validation {
    use super::PrecisionStandards;
    use crate::types::SarError;
    
    /// Validate Newton-Raphson convergence
    pub fn validate_newton_raphson_convergence(
        current: f64,
        previous: f64,
        iteration: usize,
        max_iterations: usize,
    ) -> Result<(), SarError> {
        let diff = (current - previous).abs();
        
        if diff < PrecisionStandards::NEWTON_RAPHSON_TOLERANCE {
            log::debug!(
                "Newton-Raphson converged in {} iterations (diff: {:.2e})",
                iteration,
                diff
            );
            return Ok(());
        }
        
        if iteration >= max_iterations {
            return Err(SarError::NumericalError(format!(
                "Newton-Raphson failed to converge after {} iterations (diff: {:.2e}, tolerance: {:.2e})",
                max_iterations,
                diff,
                PrecisionStandards::NEWTON_RAPHSON_TOLERANCE
            )));
        }
        
        Ok(())
    }
    
    /// Validate ECEF coordinate precision
    pub fn validate_ecef_coordinates(coords: &[f64; 3]) -> Result<(), SarError> {
        let magnitude = (coords[0].powi(2) + coords[1].powi(2) + coords[2].powi(2)).sqrt();
        
        // ECEF coordinates should be Earth-scale (6-7 million meters)
        if magnitude < 6.0e6 || magnitude > 7.0e6 {
            return Err(SarError::InvalidParameter(format!(
                "ECEF coordinates appear invalid: magnitude {:.1}m (expected 6-7 million meters)",
                magnitude
            )));
        }
        
        // Check for precision loss (values too close to zero)
        for (i, &coord) in coords.iter().enumerate() {
            if coord.abs() < f64::EPSILON * 1e6 {
                log::warn!(
                    "ECEF coordinate {} = {:.2e} may have precision loss",
                    i,
                    coord
                );
            }
        }
        
        Ok(())
    }
    
    /// Validate energy conservation in signal processing
    pub fn validate_energy_conservation(
        input_energy: f64,
        output_energy: f64,
        operation: &str,
    ) -> Result<(), SarError> {
        let energy_ratio = output_energy / input_energy;
        let energy_error = (energy_ratio - 1.0).abs();
        
        if energy_error > PrecisionStandards::ENERGY_CONSERVATION_TOLERANCE {
            return Err(SarError::NumericalError(format!(
                "Energy not conserved in {}: input={:.6e}, output={:.6e}, ratio={:.6}, error={:.2e}",
                operation,
                input_energy,
                output_energy,
                energy_ratio,
                energy_error
            )));
        }
        
        log::debug!(
            "Energy conservation verified for {}: ratio={:.6}, error={:.2e}",
            operation,
            energy_ratio,
            energy_error
        );
        
        Ok(())
    }
}

/// Initialize global precision standards
pub fn initialize_precision_standards() {
    log::info!("🎯 Initializing numeric precision standards:");
    log::info!("  Newton-Raphson tolerance: {:.2e}", PrecisionStandards::NEWTON_RAPHSON_TOLERANCE);
    log::info!("  ECEF tolerance: {:.2e}", PrecisionStandards::ECEF_TOLERANCE);
    log::info!("  Doppler tolerance: {:.2e}", PrecisionStandards::DOPPLER_TOLERANCE);
    log::info!("  Geometry tolerance: {:.2e}", PrecisionStandards::GEOMETRY_TOLERANCE);
    log::info!("  Energy conservation tolerance: {:.2e}", PrecisionStandards::ENERGY_CONSERVATION_TOLERANCE);
    
    // Initialize global seed for reproducibility
    let initial_seed = GLOBAL_SEED_COUNTER.load(Ordering::SeqCst);
    log::info!("🎲 Deterministic RNG initialized with seed: 0x{:016X}", initial_seed);
}

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_deterministic_rng_reproducibility() {
        let mut rng1 = DeterministicRng::with_seed(12345);
        let mut rng2 = DeterministicRng::with_seed(12345);
        
        // Same seed should produce same sequence
        for _ in 0..100 {
            assert_eq!(rng1.gen_f64(), rng2.gen_f64());
        }
    }
    
    #[test]
    fn test_high_precision_dot_product() {
        use geometry_precision::*;
        
        let a = vec![1.0, 2.0, 3.0];
        let b = vec![4.0, 5.0, 6.0];
        
        let result = dot_product_f64(&a, &b);
        let expected = 1.0 * 4.0 + 2.0 * 5.0 + 3.0 * 6.0;
        
        assert!(approx_eq_f64(result, expected, 1e-15));
    }
    
    #[test]
    fn test_deterministic_parallel_sum() {
        use deterministic_parallel::*;
        
        let data: Vec<f64> = (0..10000).map(|i| i as f64 * 0.1).collect();
        
        // Multiple runs should give identical results
        let sum1 = parallel_sum_f64(&data);
        let sum2 = parallel_sum_f64(&data);
        let sum3 = parallel_sum_f64(&data);
        
        assert_eq!(sum1, sum2);
        assert_eq!(sum2, sum3);
        
        // Should match sequential sum (within numerical precision)
        let seq_sum: f64 = data.iter().sum();
        let diff = (sum1 - seq_sum).abs();
        assert!(diff < 1e-10, "Parallel sum differs from sequential: {:.2e}", diff);
    }
    
    #[test]
    fn test_newton_raphson_convergence_validation() {
        use validation::*;
        
        // Should pass when converged
        let result = validate_newton_raphson_convergence(
            1.0000000001,
            1.0,
            5,
            100,
        );
        assert!(result.is_ok());
        
        // Should fail when not converged after max iterations
        let result = validate_newton_raphson_convergence(
            1.1,
            1.0,
            100,
            100,
        );
        assert!(result.is_err());
    }
    
    #[test]
    fn test_energy_conservation_validation() {
        use validation::*;
        
        // Should pass when energy is conserved
        let result = validate_energy_conservation(
            1000.0,
            1000.0000001,
            "test_operation",
        );
        assert!(result.is_ok());
        
        // Should fail when energy is not conserved
        let result = validate_energy_conservation(
            1000.0,
            1100.0,
            "test_operation",
        );
        assert!(result.is_err());
    }
}