//! Parallel and NUMA-optimized filter processing
//!
//! Provides high-performance parallel processing for large images.

use crate::core::speckle_filter::config::TileSize;
use crate::core::speckle_filter::stats::IntegralImage;
use crate::types::{SarError, SarResult};
use ndarray::{s, Array2};

/// Parallel processing implementations for SpeckleFilter
impl super::SpeckleFilter {
    /// NUMA-aware parallel processing for very large images
    #[cfg(feature = "parallel")]
    pub fn apply_filter_numa_optimized(
        &self,
        image: &Array2<f32>,
        filter_type: super::SpeckleFilterType,
    ) -> SarResult<Array2<f32>> {
        use rayon::prelude::*;

        log::info!(
            "Applying {:?} speckle filter with NUMA-optimized parallel processing",
            filter_type
        );

        let (height, width) = image.dim();
        let total_pixels = height * width;

        // Use NUMA optimization only for very large images
        if total_pixels < 50_000_000 {
            return self.apply_filter_tiled(image, filter_type, None);
        }

        let num_threads = rayon::current_num_threads();
        let numa_regions = num_threads.min(8); // Limit NUMA regions
        let rows_per_region = height / numa_regions;

        log::debug!(
            "NUMA processing: {} regions, {} rows per region",
            numa_regions,
            rows_per_region
        );

        // Process in NUMA-aware regions
        let region_results: Vec<_> = (0..numa_regions)
            .into_par_iter()
            .map(|region_id| {
                let start_row = region_id * rows_per_region;
                let end_row = if region_id == numa_regions - 1 {
                    height
                } else {
                    (region_id + 1) * rows_per_region
                };

                // Extract region with padding for window operations
                let half_window = self.params.window_size / 2;
                let padded_start = start_row.saturating_sub(half_window);
                let padded_end = (end_row + half_window).min(height);

                let region_slice = image.slice(s![padded_start..padded_end, ..]).to_owned();

                // Apply filter to region
                let filtered_region =
                    self.apply_filter_tiled(&region_slice, filter_type, Some(256))?;

                // Calculate copy parameters
                let copy_start = half_window.min(start_row - padded_start);
                let copy_height = end_row - start_row;

                Ok::<_, SarError>((start_row, copy_start, copy_height, filtered_region))
            })
            .collect::<Result<Vec<_>, _>>()?;

        // Combine results
        let mut result = Array2::zeros((height, width));
        for (start_row, copy_start, copy_height, region_data) in region_results {
            result
                .slice_mut(s![start_row..start_row + copy_height, ..])
                .assign(&region_data.slice(s![copy_start..copy_start + copy_height, ..]));
        }

        log::info!("NUMA-optimized parallel processing completed successfully");
        Ok(result)
    }

    #[cfg(not(feature = "parallel"))]
    pub fn apply_filter_numa_optimized(
        &self,
        image: &Array2<f32>,
        filter_type: super::SpeckleFilterType,
    ) -> SarResult<Array2<f32>> {
        // Fallback to tiled processing if parallel feature is not available
        self.apply_filter_tiled(image, filter_type, None)
    }

    /// Apply speckle filtering with tiled processing for optimal cache performance
    pub fn apply_filter_tiled(
        &self,
        image: &Array2<f32>,
        filter_type: super::SpeckleFilterType,
        tile_size: Option<usize>,
    ) -> SarResult<Array2<f32>> {
        self.params.validate_window()?;

        log::info!(
            "Applying {:?} speckle filter with tiled processing",
            filter_type
        );

        // Gamma-MAP needs consistent border handling; reuse the non-tiled path for bit-identical output
        if matches!(filter_type, super::SpeckleFilterType::GammaMAP) {
            return self.apply_filter(image, filter_type);
        }

        let (height, width) = image.dim();
        let effective_tile_size = self.determine_optimal_tile_size(height, width, tile_size);

        log::debug!(
            "Using tile size: {}x{}",
            effective_tile_size,
            effective_tile_size
        );

        // For small images, use standard processing
        if height <= effective_tile_size && width <= effective_tile_size {
            return self.apply_filter(image, filter_type);
        }

        let half_window = self.params.window_size / 2;

        // OPTIMIZATION #46: Collect tile coordinates for parallel processing
        let mut tile_coords: Vec<(usize, usize, usize, usize)> = Vec::new();
        for tile_row in (0..height).step_by(effective_tile_size) {
            for tile_col in (0..width).step_by(effective_tile_size) {
                let tile_end_row = (tile_row + effective_tile_size).min(height);
                let tile_end_col = (tile_col + effective_tile_size).min(width);
                tile_coords.push((tile_row, tile_col, tile_end_row, tile_end_col));
            }
        }

        // Process tiles in parallel
        use rayon::prelude::*;
        let tile_results: Vec<_> = tile_coords
            .par_iter()
            .map(|&(tile_row, tile_col, tile_end_row, tile_end_col)| {
                // Extract tile with padding for window operations
                let padded_start_row = tile_row.saturating_sub(half_window);
                let padded_end_row = (tile_end_row + half_window).min(height);
                let padded_start_col = tile_col.saturating_sub(half_window);
                let padded_end_col = (tile_end_col + half_window).min(width);

                let tile = image
                    .slice(s![
                        padded_start_row..padded_end_row,
                        padded_start_col..padded_end_col
                    ])
                    .to_owned();

                // Build per-tile integral image when window >= 9
                let tile_integral = if self.params.window_size >= 9 {
                    Some(IntegralImage::new(&tile))
                } else {
                    None
                };

                // Process tile with its own integral image
                let filtered_tile = self
                    .apply_filter_to_tile(&tile, filter_type, tile_integral.as_ref())
                    .unwrap_or_else(|_| tile.clone());

                // Calculate offsets for copying results back
                let copy_start_row = half_window.min(tile_row - padded_start_row);
                let copy_start_col = half_window.min(tile_col - padded_start_col);
                let copy_height = tile_end_row - tile_row;
                let copy_width = tile_end_col - tile_col;

                (
                    tile_row,
                    tile_col,
                    tile_end_row,
                    tile_end_col,
                    copy_start_row,
                    copy_start_col,
                    copy_height,
                    copy_width,
                    filtered_tile,
                )
            })
            .collect();

        // Merge results (sequential - writes to shared array)
        let mut result = Array2::zeros((height, width));
        for (
            tile_row,
            tile_col,
            tile_end_row,
            tile_end_col,
            copy_start_row,
            copy_start_col,
            copy_height,
            copy_width,
            filtered_tile,
        ) in tile_results
        {
            result
                .slice_mut(s![tile_row..tile_end_row, tile_col..tile_end_col])
                .assign(&filtered_tile.slice(s![
                    copy_start_row..copy_start_row + copy_height,
                    copy_start_col..copy_start_col + copy_width
                ]));
        }

        log::info!("Tiled speckle filtering completed successfully");
        Ok(result)
    }

    /// Determine optimal tile size based on image dimensions and cache considerations
    pub(super) fn determine_optimal_tile_size(
        &self,
        height: usize,
        width: usize,
        user_tile_size: Option<usize>,
    ) -> usize {
        if let Some(size) = user_tile_size {
            return size;
        }

        if self.params.tile_size > 0 {
            return self.params.tile_size;
        }

        // Auto-select tile size based on image size and memory considerations
        let total_pixels = height * width;

        if total_pixels < 1_000_000 {
            // Small images: use smaller tiles for better cache locality
            TileSize::Small as usize
        } else if total_pixels < 10_000_000 {
            // Medium images: balanced approach
            TileSize::Medium as usize
        } else {
            // Large images: use larger tiles to reduce overhead
            TileSize::Large as usize
        }
    }

    /// Apply filter to a single tile with optimized processing
    pub(super) fn apply_filter_to_tile(
        &self,
        tile: &Array2<f32>,
        filter_type: super::SpeckleFilterType,
        integral: Option<&IntegralImage>,
    ) -> SarResult<Array2<f32>> {
        let (tile_height, tile_width) = tile.dim();
        let mut filtered = Array2::zeros((tile_height, tile_width));

        match filter_type {
            super::SpeckleFilterType::Lee => {
                self.apply_lee_filter_to_tile(tile, &mut filtered, integral)
            }
            super::SpeckleFilterType::EnhancedLee => {
                self.apply_enhanced_lee_filter_to_tile(tile, &mut filtered, integral)
            }
            super::SpeckleFilterType::GammaMAP => {
                self.apply_gamma_map_filter_to_tile(tile, &mut filtered, integral)
            }
            super::SpeckleFilterType::Mean => {
                self.apply_mean_filter_to_tile(tile, &mut filtered, integral)
            }
            super::SpeckleFilterType::Median => {
                self.apply_median_filter_to_tile(tile, &mut filtered)
            }
            _ => {
                // For other filters, use standard implementation
                return self.apply_filter(tile, filter_type);
            }
        }?;

        Ok(filtered)
    }
}
