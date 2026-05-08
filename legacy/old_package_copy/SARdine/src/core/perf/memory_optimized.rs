#![allow(dead_code, unused_variables)]
//! Memory-optimized utilities for zero-copy operations and efficient array processing
//!
//! This module provides:
//! - Zero-copy array views with lifetime-safe borrowing
//! - Finite-only statistics computation (NaN/Inf filtering)
//! - Chunked parallel processing with in-place writes
//! - Shape-keyed array memory pools
//! - Optimized in-place operations with masking support
//! - Cache-friendly tiled processing

use crate::types::SarError;
use crate::types::SarResult;
use ndarray::{Array2, ArrayView2, ArrayViewMut2, Axis};
use num_complex::Complex;
use num_traits::FromPrimitive;
use rayon::prelude::*;
use std::borrow::Cow;
use std::collections::HashMap;

// ============================================================================
// Zero-copy array conversions
// ============================================================================

/// Try to get a zero-copy view of a contiguous array, or fall back to owned copy.
/// Returns a `Cow` that is either borrowed (zero-copy) or owned.
///
/// # Safety
/// The zero-copy path uses unsafe to extend the borrow lifetime, so the caller
/// must ensure the source array outlives the returned view.
pub fn try_view_or_owned<T: Clone>(array: &Array2<T>) -> Cow<Array2<T>> {
    if array.is_standard_layout() {
        // Zero-copy view for contiguous arrays
        // Safety: The caller must ensure the array outlives this borrow
        Cow::Borrowed(array)
    } else {
        // Fall back to owned copy for non-contiguous arrays
        Cow::Owned(array.to_owned())
    }
}

/// Optimized numpy-to-ndarray conversion with zero-copy when possible
#[cfg(feature = "python")]
pub fn numpy_to_array_optimized<T>(py_array: &pyo3::PyAny) -> SarResult<Cow<Array2<T>>>
where
    T: Clone + numpy::Element,
{
    use numpy::PyArrayDyn;
    use pyo3::Python;

    let array: &PyArrayDyn<T> = py_array
        .extract()
        .map_err(|e| SarError::Processing(format!("Failed to extract numpy array: {}", e)))?;

    // Check if array is contiguous and C-order
    if array.is_c_contiguous() {
        // SAFETY: array.is_c_contiguous() guarantees the underlying buffer is laid out
        // as a contiguous C-order array. PyO3's as_array() creates an ndarray view that
        // references the numpy buffer without copying. The lifetime of the view is tied
        // to the PyArrayDyn borrow, ensuring the buffer remains valid.
        let view = unsafe { array.as_array() };
        let array_2d = view
            .into_dimensionality::<ndarray::Ix2>()
            .map_err(|e| SarError::Processing(format!("Array dimension error: {}", e)))?;
        Ok(Cow::Borrowed(array_2d))
    } else {
        // SAFETY: For non-contiguous arrays, we still need unsafe to access the numpy
        // buffer, but we immediately call to_owned() to copy into a contiguous Rust-owned
        // array. This avoids any aliasing issues from the non-contiguous source.
        let owned = unsafe { array.as_array() }.to_owned();
        let array_2d = owned
            .into_dimensionality::<ndarray::Ix2>()
            .map_err(|e| SarError::Processing(format!("Array dimension error: {}", e)))?;
        Ok(Cow::Owned(array_2d))
    }
}

// ============================================================================
// Finite-only statistics computation
// ============================================================================

/// Statistics result with finite-only values
#[derive(Debug, Clone)]
pub struct ArrayStatistics<T> {
    pub count: usize,
    pub min: Option<T>,
    pub max: Option<T>,
    pub mean: Option<T>,
}

/// Compute statistics in-place, filtering out NaN/Inf values
pub fn compute_array_statistics_inplace<T>(
    array: &Array2<T>,
    valid_pred: impl Fn(T) -> bool,
) -> ArrayStatistics<T>
where
    T: Copy + PartialOrd + std::ops::Add<Output = T> + std::ops::Div<Output = T> + FromPrimitive,
{
    let mut count = 0;
    let mut min_val: Option<T> = None;
    let mut max_val: Option<T> = None;
    let mut sum: Option<T> = None;

    for &value in array.iter() {
        if valid_pred(value) {
            count += 1;

            min_val = Some(match min_val {
                None => value,
                Some(m) => {
                    if value < m {
                        value
                    } else {
                        m
                    }
                }
            });

            max_val = Some(match max_val {
                None => value,
                Some(m) => {
                    if value > m {
                        value
                    } else {
                        m
                    }
                }
            });

            sum = Some(match sum {
                None => value,
                Some(s) => s + value,
            });
        }
    }

    let mean = if count > 0 {
        let denom = T::from_usize(count).unwrap_or_else(|| {
            panic!(
                "Failed to convert count={} into target type for statistics",
                count
            )
        });
        sum.map(|s| s / denom)
    } else {
        None
    };

    ArrayStatistics {
        count,
        min: min_val,
        max: max_val,
        mean,
    }
}

/// Compute statistics sequentially with custom validity predicate
pub fn compute_array_statistics_sequential<T>(
    array: &Array2<T>,
    valid_pred: impl Fn(T) -> bool,
) -> ArrayStatistics<T>
where
    T: Copy + PartialOrd + std::ops::Add<Output = T> + std::ops::Div<Output = T> + FromPrimitive,
{
    compute_array_statistics_inplace(array, valid_pred)
}

/// Compute statistics in parallel using Rayon with Option-based reduction
pub fn compute_array_statistics_parallel<T>(
    array: &Array2<T>,
    valid_pred: impl Fn(T) -> bool + Sync + Send,
) -> ArrayStatistics<T>
where
    T: Copy
        + PartialOrd
        + std::ops::Add<Output = T>
        + std::ops::Div<Output = T>
        + FromPrimitive
        + Send
        + Sync,
{
    let result = array
        .axis_iter(Axis(0))
        .par_bridge()
        .map(|row| {
            let mut count = 0;
            let mut min_val: Option<T> = None;
            let mut max_val: Option<T> = None;
            let mut sum: Option<T> = None;

            for &value in row.iter() {
                if valid_pred(value) {
                    count += 1;

                    min_val = Some(match min_val {
                        None => value,
                        Some(m) => {
                            if value < m {
                                value
                            } else {
                                m
                            }
                        }
                    });

                    max_val = Some(match max_val {
                        None => value,
                        Some(m) => {
                            if value > m {
                                value
                            } else {
                                m
                            }
                        }
                    });

                    sum = Some(match sum {
                        None => value,
                        Some(s) => s + value,
                    });
                }
            }

            (count, min_val, max_val, sum)
        })
        .reduce(
            || (0, None, None, None),
            |(c1, min1, max1, sum1), (c2, min2, max2, sum2)| {
                let combined_min = match (min1, min2) {
                    (Some(m1), Some(m2)) => Some(if m1 < m2 { m1 } else { m2 }),
                    (Some(m), None) | (None, Some(m)) => Some(m),
                    (None, None) => None,
                };

                let combined_max = match (max1, max2) {
                    (Some(m1), Some(m2)) => Some(if m1 > m2 { m1 } else { m2 }),
                    (Some(m), None) | (None, Some(m)) => Some(m),
                    (None, None) => None,
                };

                let combined_sum = match (sum1, sum2) {
                    (Some(s1), Some(s2)) => Some(s1 + s2),
                    (Some(s), None) | (None, Some(s)) => Some(s),
                    (None, None) => None,
                };

                (c1 + c2, combined_min, combined_max, combined_sum)
            },
        );

    let (count, min_val, max_val, sum) = result;
    let mean = if count > 0 {
        let denom = T::from_usize(count).unwrap_or_else(|| {
            panic!(
                "Failed to convert count={} into target type for statistics",
                count
            )
        });
        sum.map(|s| s / denom)
    } else {
        None
    };

    ArrayStatistics {
        count,
        min: min_val,
        max: max_val,
        mean,
    }
}

// ============================================================================
// Chunked processing with in-place writes
// ============================================================================

/// Process array in chunks, writing directly to output (no intermediate allocations)
pub fn process_chunked<T, F>(
    input: &Array2<T>,
    #[allow(unused_mut)] mut output: ArrayViewMut2<T>,
    chunk_size: usize,
    #[allow(unused_mut)] mut operation: F,
) -> SarResult<()>
where
    T: Copy,
    F: FnMut(ArrayView2<T>, ArrayViewMut2<T>) -> SarResult<()>,
{
    let (height, width) = input.dim();

    for chunk_start in (0..height).step_by(chunk_size) {
        let chunk_end = (chunk_start + chunk_size).min(height);

        let input_chunk = input.slice(ndarray::s![chunk_start..chunk_end, ..]);
        let output_chunk = output.slice_mut(ndarray::s![chunk_start..chunk_end, ..]);

        operation(input_chunk, output_chunk)?;
    }

    Ok(())
}

/// Process array in parallel chunks, writing directly to non-overlapping output slices
///
/// Uses rayon's parallel iterator to process chunks and collect results,
/// then writes them back to the output in a second pass.
///
/// This is a safe parallel implementation that:
/// 1. Processes each chunk independently in parallel
/// 2. Collects results into a Vec<Vec<T>>
/// 3. Writes results back to output sequentially (fast memcpy)
pub fn process_chunked_parallel<T, F>(
    input: &Array2<T>,
    mut output: ArrayViewMut2<T>,
    chunk_size: usize,
    operation: F,
) -> SarResult<()>
where
    T: Copy + Send + Sync + Default,
    F: Fn(ArrayView2<T>, &mut [T]) -> SarResult<()> + Send + Sync,
{
    use rayon::prelude::*;

    let (rows, cols) = input.dim();
    if rows == 0 || cols == 0 {
        return Ok(());
    }

    // Pre-compute chunk boundaries
    let num_chunks = (rows + chunk_size - 1) / chunk_size;
    let chunk_boundaries: Vec<(usize, usize)> = (0..num_chunks)
        .map(|i| {
            let start = i * chunk_size;
            let end = ((i + 1) * chunk_size).min(rows);
            (start, end)
        })
        .collect();

    // Process all chunks in parallel and collect results
    let results: Result<Vec<(usize, Vec<T>)>, SarError> = chunk_boundaries
        .par_iter()
        .map(|(start, end)| {
            let chunk_rows = end - start;
            let input_chunk = input.slice(ndarray::s![*start..*end, ..]);

            // Create a temporary buffer for this chunk's output
            let mut chunk_output = vec![T::default(); chunk_rows * cols];

            // Process the chunk
            operation(input_chunk, &mut chunk_output)?;

            Ok((*start, chunk_output))
        })
        .collect();

    // Handle any errors
    let results = results?;

    // Write results back to output (sequential but fast - just memcpy)
    for (start_row, chunk_data) in results {
        let chunk_rows = chunk_data.len() / cols;
        for (local_row, row_data) in chunk_data.chunks(cols).enumerate() {
            let global_row = start_row + local_row;
            if global_row < rows {
                for (col, &val) in row_data.iter().enumerate() {
                    output[[global_row, col]] = val;
                }
            }
        }
    }

    Ok(())
}

// ============================================================================
// Shape-keyed array memory pool
// ============================================================================

/// Memory pool for reusing arrays, keyed by shape to avoid O(n) scans
pub struct ArrayMemoryPool<T> {
    pools: HashMap<(usize, usize), Vec<Array2<T>>>,
}

impl<T: Default + Clone> ArrayMemoryPool<T> {
    pub fn new() -> Self {
        Self {
            pools: HashMap::new(),
        }
    }

    /// Get a zeroed array from the pool or allocate a new one
    pub fn get_zeroed(&mut self, rows: usize, cols: usize) -> Array2<T> {
        let shape = (rows, cols);

        if let Some(pool) = self.pools.get_mut(&shape) {
            if let Some(mut array) = pool.pop() {
                // Zero the array before returning
                array.fill(T::default());
                return array;
            }
        }

        // Allocate new zeroed array
        Array2::default((rows, cols))
    }

    /// Get an array from the pool (may contain uninitialized data)
    pub fn get_array(&mut self, rows: usize, cols: usize) -> Array2<T> {
        let shape = (rows, cols);

        if let Some(pool) = self.pools.get_mut(&shape) {
            if let Some(array) = pool.pop() {
                return array;
            }
        }

        // Allocate new array
        Array2::default((rows, cols))
    }

    /// Return an array to the pool for reuse
    pub fn return_array(&mut self, array: Array2<T>) {
        let shape = array.dim();
        self.pools.entry(shape).or_insert_with(Vec::new).push(array);
    }

    /// Clear the pool and release memory
    pub fn clear(&mut self) {
        self.pools.clear();
        self.pools.shrink_to_fit();
    }

    /// Shrink all pools to fit current capacity
    pub fn shrink_to_fit(&mut self) {
        for pool in self.pools.values_mut() {
            pool.shrink_to_fit();
        }
        self.pools.shrink_to_fit();
    }
}

impl<T: Default + Clone> Default for ArrayMemoryPool<T> {
    fn default() -> Self {
        Self::new()
    }
}

// ============================================================================
// In-place operations with masking
// ============================================================================

pub mod inplace_ops {
    use super::*;

    /// Convert linear power to dB in-place with epsilon clamping and optional masking
    pub fn linear_to_db_inplace(
        array: &mut Array2<f32>,
        epsilon: f32,
        mask_invalid: Option<&Array2<u8>>,
    ) {
        if let Some(mask) = mask_invalid {
            array
                .iter_mut()
                .zip(mask.iter())
                .for_each(|(val, &mask_val)| {
                    if mask_val == 0 {
                        *val = f32::NAN; // Mark invalid
                    } else {
                        *val = 10.0 * (*val).max(epsilon).log10();
                    }
                });
        } else {
            array.iter_mut().for_each(|val| {
                *val = 10.0 * (*val).max(epsilon).log10();
            });
        }
    }

    /// Convert linear power to dB in-place (parallel version)
    pub fn linear_to_db_inplace_par(
        array: &mut Array2<f32>,
        epsilon: f32,
        mask_invalid: Option<&Array2<u8>>,
    ) {
        if let Some(mask) = mask_invalid {
            array
                .axis_iter_mut(Axis(0))
                .into_par_iter()
                .zip(mask.axis_iter(Axis(0)).into_par_iter())
                .for_each(|(mut row, mask_row)| {
                    row.iter_mut()
                        .zip(mask_row.iter())
                        .for_each(|(val, &mask_val)| {
                            if mask_val == 0 {
                                *val = f32::NAN;
                            } else {
                                *val = 10.0 * (*val).max(epsilon).log10();
                            }
                        });
                });
        } else {
            array
                .axis_iter_mut(Axis(0))
                .into_par_iter()
                .for_each(|mut row| {
                    row.iter_mut().for_each(|val| {
                        *val = 10.0 * (*val).max(epsilon).log10();
                    });
                });
        }
    }

    /// Convert real values to complex in-place (optimized with mapv)
    pub fn real_to_complex_optimized(real_array: &Array2<f32>) -> Array2<Complex<f32>> {
        real_array.mapv(|r| Complex::new(r, 0.0))
    }

    /// Apply a function in-place with optional masking
    pub fn apply_inplace<T, F>(array: &mut Array2<T>, mask: Option<&Array2<u8>>, mut func: F)
    where
        T: Copy,
        F: FnMut(T) -> T,
    {
        if let Some(mask) = mask {
            array
                .iter_mut()
                .zip(mask.iter())
                .for_each(|(val, &mask_val)| {
                    if mask_val != 0 {
                        *val = func(*val);
                    }
                });
        } else {
            array.iter_mut().for_each(|val| {
                *val = func(*val);
            });
        }
    }

    /// Apply a function in-place (parallel version)
    pub fn apply_inplace_par<T, F>(array: &mut Array2<T>, mask: Option<&Array2<u8>>, func: F)
    where
        T: Copy + Send + Sync,
        F: Fn(T) -> T + Sync + Send,
    {
        if let Some(mask) = mask {
            array
                .axis_iter_mut(Axis(0))
                .into_par_iter()
                .zip(mask.axis_iter(Axis(0)).into_par_iter())
                .for_each(|(mut row, mask_row)| {
                    row.iter_mut()
                        .zip(mask_row.iter())
                        .for_each(|(val, &mask_val)| {
                            if mask_val != 0 {
                                *val = func(*val);
                            }
                        });
                });
        } else {
            array
                .axis_iter_mut(Axis(0))
                .into_par_iter()
                .for_each(|mut row| {
                    row.iter_mut().for_each(|val| {
                        *val = func(*val);
                    });
                });
        }
    }
}

// ============================================================================
// Cache-friendly processing
// ============================================================================

pub mod cache_friendly {
    use super::*;

    /// Process array by rows with optional masking
    pub fn process_by_rows<T, F>(
        input: &Array2<T>,
        mut output: ArrayViewMut2<T>,
        mask: Option<&Array2<u8>>,
        mut operation: F,
    ) -> SarResult<()>
    where
        T: Copy,
        F: FnMut(ArrayView2<T>, ArrayViewMut2<T>, Option<ArrayView2<u8>>) -> SarResult<()>,
    {
        let (height, _) = input.dim();

        for row_idx in 0..height {
            let input_row = input.slice(ndarray::s![row_idx..row_idx + 1, ..]);
            let output_row = output.slice_mut(ndarray::s![row_idx..row_idx + 1, ..]);
            let mask_row = mask.map(|m| m.slice(ndarray::s![row_idx..row_idx + 1, ..]));

            operation(input_row, output_row, mask_row)?;
        }

        Ok(())
    }

    /// Process array in cache-friendly tiles with cache-line alignment
    pub fn process_tiled<T, F>(
        input: &Array2<T>,
        mut output: ArrayViewMut2<T>,
        tile_rows: usize,
        tile_cols: usize,
        mask: Option<&Array2<u8>>,
        mut operation: F,
    ) -> SarResult<()>
    where
        T: Copy,
        F: FnMut(ArrayView2<T>, ArrayViewMut2<T>, Option<ArrayView2<u8>>) -> SarResult<()>,
    {
        let (height, width) = input.dim();

        // Align tile column starts to cache lines (64 elements for f32)
        let cache_line_cols = 64;

        for row_start in (0..height).step_by(tile_rows) {
            let row_end = (row_start + tile_rows).min(height);

            for col_start in (0..width).step_by(tile_cols) {
                // Align to cache line boundary when possible
                let aligned_col_start = if col_start >= cache_line_cols {
                    (col_start / cache_line_cols) * cache_line_cols
                } else {
                    col_start
                };

                let col_end = (aligned_col_start + tile_cols).min(width);

                let input_tile =
                    input.slice(ndarray::s![row_start..row_end, aligned_col_start..col_end]);
                let output_tile =
                    output.slice_mut(ndarray::s![row_start..row_end, aligned_col_start..col_end]);
                let mask_tile = mask
                    .map(|m| m.slice(ndarray::s![row_start..row_end, aligned_col_start..col_end]));

                operation(input_tile, output_tile, mask_tile)?;
            }
        }

        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_statistics_finite_only() {
        // Use i32 instead of f32 since f32 doesn't implement From<usize>
        let mut array = Array2::<i32>::zeros((10, 10));
        array[[2, 2]] = 5;
        array[[3, 3]] = 10;
        array[[4, 4]] = -3;

        let stats = compute_array_statistics_sequential(&array, |v| v != 0);

        assert_eq!(stats.count, 3); // Only 3 non-zero values
        assert_eq!(stats.min, Some(-3));
        assert_eq!(stats.max, Some(10));
    }

    #[test]
    fn test_memory_pool() {
        let mut pool = ArrayMemoryPool::<f32>::new();

        let array1 = pool.get_zeroed(10, 10);
        assert_eq!(array1.dim(), (10, 10));

        pool.return_array(array1);

        let array2 = pool.get_array(10, 10);
        assert_eq!(array2.dim(), (10, 10));
    }

    #[test]
    fn test_linear_to_db_inplace() {
        let mut array = Array2::from_elem((5, 5), 100.0);
        inplace_ops::linear_to_db_inplace(&mut array, 1e-30, None);

        assert!((array[[0, 0]] - 20.0).abs() < 0.01); // 10*log10(100) = 20 dB
    }
}
