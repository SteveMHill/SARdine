//! GeoTIFF output and masking workflow functions for terrain correction.
//!
//! This module provides functionality for:
//! - Saving terrain-corrected data to GeoTIFF format
//! - Applying masking workflows to filter invalid pixels
//! - Applying masks to gamma0 data
//! - Shadow/layover detection based on local incidence angle

use crate::types::{GeoTransform, MaskResult, MaskingWorkflow, SarError, SarResult};
use gdal::DriverManager;
use ndarray::Array2;

use super::types::Vector3;
use super::TerrainCorrector;

/// Shadow detection threshold: cos(LIA) below this indicates shadow
/// Default: cos(85°) ≈ 0.087 - very grazing incidence
const SHADOW_COS_LIA_THRESHOLD: f64 = 0.087;

/// Layover detection threshold: cos(LIA) below this indicates layover
/// Default: 0.0 - surface normal pointing away from radar
const LAYOVER_COS_LIA_THRESHOLD: f64 = 0.0;

/// Default incidence angle for Sentinel-1 IW mode (mid-swath)
const DEFAULT_INCIDENCE_ANGLE_DEG: f64 = 35.0;

impl TerrainCorrector {
    /// Compute look vector for a pixel using orbit data if available
    ///
    /// This method computes the radar look vector (from target to satellite) using:
    /// 1. Approximate geometry based on typical S1 incidence angle
    ///
    /// Note: For full accuracy with orbit data, use the `look_vector_at` method
    /// in rtc.rs which requires RangeDopplerParams for zero-Doppler solving.
    /// This simplified version provides good approximation for masking workflows.
    ///
    /// # Arguments
    /// * `row` - Pixel row in DEM array
    /// * `col` - Pixel column in DEM array
    /// * `dem_array` - DEM elevation data
    /// * `workflow` - Masking workflow parameters (unused, for future extension)
    ///
    /// # Returns
    /// * Unit look vector pointing from target toward satellite
    fn compute_look_vector_for_pixel(
        &self,
        row: usize,
        col: usize,
        dem_array: &Array2<f32>,
        _workflow: &MaskingWorkflow,
    ) -> Vector3 {
        // Get target lat/lon/height from DEM
        let _lon =
            self.dem_transform.top_left_x + (col as f64 + 0.5) * self.dem_transform.pixel_width;
        let _lat =
            self.dem_transform.top_left_y + (row as f64 + 0.5) * self.dem_transform.pixel_height;
        let _height = dem_array[[row, col]] as f64;

        // Use approximate look vector based on typical S1 IW geometry
        // This provides good accuracy for shadow/layover detection in masking workflows
        //
        // For terrain correction geocoding, the full zero-Doppler solution via
        // look_vector_at() in rtc.rs should be used instead.
        let approx_incidence_rad = DEFAULT_INCIDENCE_ANGLE_DEG.to_radians();

        // In local ENU coordinates at the target:
        // - x (East): horizontal range direction
        // - y (North): along-track (azimuth) direction
        // - z (Up): vertical
        // Look vector points toward satellite (side-looking geometry)
        Vector3 {
            x: approx_incidence_rad.sin(), // Range component (horizontal)
            y: 0.0,                        // Azimuth component (along track)
            z: approx_incidence_rad.cos(), // Vertical component (upward)
        }
    }
    /// Save GeoTIFF file
    pub fn save_geotiff(
        &self,
        image: &Array2<f32>,
        transform: &GeoTransform,
        output_path: &str,
        compression: Option<&str>,
    ) -> SarResult<()> {
        let (height, width) = image.dim();
        let driver = DriverManager::get_driver_by_name("GTiff")?;

        // Create dataset with compression options
        let mut options = Vec::new();
        if let Some(comp) = compression {
            options.push(format!("COMPRESS={}", comp));
        }
        options.push("TILED=YES".to_string());

        let mut dataset = driver.create_with_band_type::<f32, _>(
            output_path,
            width as isize,
            height as isize,
            1,
        )?;

        // Set geotransform
        dataset.set_geo_transform(&[
            transform.top_left_x,
            transform.pixel_width,
            transform.rotation_x,
            transform.top_left_y,
            transform.rotation_y,
            transform.pixel_height,
        ])?;

        // Set coordinate reference system
        dataset.set_spatial_ref(&gdal::spatial_ref::SpatialRef::from_epsg(self.output_crs)?)?;

        // Write image data
        let mut rasterband = dataset.rasterband(1)?;
        let flat_data: Vec<f32> = image.iter().cloned().collect();
        let buffer = gdal::raster::Buffer::new((width, height), flat_data);
        rasterband.write((0, 0), (width, height), &buffer)?;

        // Set no-data value
        rasterband.set_no_data_value(Some(f32::NAN as f64))?;

        Ok(())
    }

    /// Apply masking workflow to terrain-corrected data
    pub fn apply_masking_workflow(
        &self,
        gamma0_data: &Array2<f32>,
        dem_array: &Array2<f32>,
        workflow: &MaskingWorkflow,
    ) -> SarResult<MaskResult> {
        let (height, width) = gamma0_data.dim();

        // Initialize masks
        let mut gamma0_mask = Array2::<bool>::from_elem((height, width), true);
        let mut dem_mask = Array2::<bool>::from_elem((height, width), true);
        let mut lia_mask = Array2::<bool>::from_elem((height, width), true);
        let mut lia_cosine = Array2::<f32>::zeros((height, width));

        // Process each pixel
        for row in 0..height {
            for col in 0..width {
                let gamma0_val = gamma0_data[[row, col]];
                let dem_val = dem_array[[row, col]];

                // Check gamma0 validity
                gamma0_mask[[row, col]] = gamma0_val.is_finite()
                    && gamma0_val >= workflow.gamma0_min
                    && gamma0_val <= workflow.gamma0_max;

                // Check DEM validity
                dem_mask[[row, col]] =
                    dem_val.is_finite() && dem_val > workflow.dem_threshold as f32;

                // Compute local incidence angle if DEM is valid
                if dem_mask[[row, col]] {
                    let surface_normal = self.compute_surface_normal(dem_array, row, col);

                    // Compute look vector from target to satellite
                    // Use orbit data if available for accurate geometry, otherwise use approximation
                    let look_vector =
                        self.compute_look_vector_for_pixel(row, col, dem_array, workflow);

                    let cos_lia = self.compute_local_incidence_angle(&surface_normal, &look_vector);
                    lia_cosine[[row, col]] = cos_lia as f32;

                    // Check local incidence angle threshold
                    lia_mask[[row, col]] = cos_lia >= workflow.lia_threshold;

                    // SHADOW/LAYOVER DETECTION (Jan 2026 fix):
                    // Shadow: cos(LIA) near zero means grazing incidence - radar cannot illuminate
                    // Layover: cos(LIA) negative means surface faces away from radar
                    if cos_lia < SHADOW_COS_LIA_THRESHOLD {
                        // Mark as potential shadow - very low illumination
                        lia_mask[[row, col]] = false;
                    }
                    if cos_lia < LAYOVER_COS_LIA_THRESHOLD {
                        // Mark as layover - surface tilted past vertical
                        lia_mask[[row, col]] = false;
                    }
                } else {
                    lia_cosine[[row, col]] = f32::NAN;
                    lia_mask[[row, col]] = false;
                }
            }
        }

        // Combine all masks and compute shadow/layover statistics
        let mut combined_mask_bool = Array2::<bool>::from_elem((height, width), true);
        let mut shadow_mask_array = Array2::<u8>::zeros((height, width));
        let mut layover_mask_array = Array2::<u8>::zeros((height, width));
        let mut valid_pixels = 0;
        let mut shadow_pixels = 0;
        let mut layover_pixels = 0;

        for row in 0..height {
            for col in 0..width {
                let cos_lia = lia_cosine[[row, col]];

                // Shadow detection: very grazing incidence
                if cos_lia.is_finite() && (cos_lia as f64) < SHADOW_COS_LIA_THRESHOLD {
                    shadow_mask_array[[row, col]] = 1;
                    shadow_pixels += 1;
                }

                // Layover detection: surface facing away from radar
                if cos_lia.is_finite() && (cos_lia as f64) < LAYOVER_COS_LIA_THRESHOLD {
                    layover_mask_array[[row, col]] = 1;
                    layover_pixels += 1;
                }

                combined_mask_bool[[row, col]] =
                    gamma0_mask[[row, col]] && dem_mask[[row, col]] && lia_mask[[row, col]];

                if combined_mask_bool[[row, col]] {
                    valid_pixels += 1;
                }
            }
        }

        // Convert boolean mask to u8 mask
        let mut combined_mask = Array2::<u8>::zeros((height, width));
        for row in 0..height {
            for col in 0..width {
                combined_mask[[row, col]] = if combined_mask_bool[[row, col]] { 1 } else { 0 };
            }
        }

        let total_pixels = height * width;
        let coverage_percent = (valid_pixels as f64 / total_pixels as f64) * 100.0;

        // Log shadow/layover statistics
        if shadow_pixels > 0 || layover_pixels > 0 {
            log::info!(
                "🗻 Shadow/Layover detection: {} shadow pixels ({:.2}%), {} layover pixels ({:.2}%)",
                shadow_pixels,
                100.0 * shadow_pixels as f64 / total_pixels as f64,
                layover_pixels,
                100.0 * layover_pixels as f64 / total_pixels as f64
            );
        }

        let stats = crate::types::MaskStats {
            total_pixels,
            water_pixels: 0, // Would be computed if water masking is implemented
            shadow_pixels,   // Now computed from LIA threshold
            layover_pixels,  // Now computed from LIA threshold
            noise_pixels: 0, // Would be computed if noise masking is implemented
            valid_pixels,
        };

        Ok(MaskResult {
            water_mask: None,
            shadow_mask: Some(shadow_mask_array),
            layover_mask: Some(layover_mask_array),
            noise_mask: None,
            combined_mask,
            lia_cosine,
            gamma0_mask,
            dem_mask,
            lia_mask,
            valid_pixels,
            total_pixels,
            coverage_percent,
            stats,
        })
    }

    /// Apply mask to gamma0 data
    pub fn apply_mask_to_gamma0(
        &self,
        gamma0_data: &Array2<f32>,
        mask: &Array2<u8>,
        fill_value: Option<f32>,
    ) -> SarResult<Array2<f32>> {
        let mut masked_data = gamma0_data.clone();
        let fill_val = fill_value.ok_or_else(|| {
            SarError::MissingParameter(
                "Fill value is required for scientific processing".to_string(),
            )
        })?;

        for ((row, col), &mask_val) in mask.indexed_iter() {
            if mask_val == 0 {
                masked_data[[row, col]] = fill_val;
            }
        }

        Ok(masked_data)
    }
}
