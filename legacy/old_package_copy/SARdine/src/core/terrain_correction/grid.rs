//! Grid generation functions for terrain correction output
//!
//! This module provides methods for calculating output bounds and creating
//! output grids in both geographic (WGS84) and projected (UTM, etc.) coordinate systems.

use crate::constants::geodetic::WGS84_SEMI_MAJOR_AXIS_M;
use crate::types::{BoundingBox, GeoTransform, SarError, SarResult};

use super::types::{AlgorithmStatus, ExecutionMode};
use super::TerrainCorrector;

impl TerrainCorrector {
    /// Calculate output bounds from SAR bounding box with validation
    pub(crate) fn calculate_output_bounds(&self, sar_bbox: &BoundingBox) -> SarResult<BoundingBox> {
        // Validate input bounding box
        if sar_bbox.max_lat <= sar_bbox.min_lat || sar_bbox.max_lon <= sar_bbox.min_lon {
            return Err(SarError::InvalidInput(
                "Invalid bounding box: max values must be greater than min values".to_string(),
            ));
        }

        // Validate bounding box using configuration-based thresholds
        let lat_diff = sar_bbox.max_lat - sar_bbox.min_lat;
        let lon_diff = sar_bbox.max_lon - sar_bbox.min_lon;

        // Record validation attempt
        let start_time = std::time::Instant::now();
        let mut validation_status = AlgorithmStatus {
            algorithm_name: "bounding_box_validation".to_string(),
            execution_mode: ExecutionMode::Primary,
            iterations_used: None,
            convergence_achieved: Some(true),
            fallback_reason: None,
            processing_time_ms: 0.0,
        };

        // Use scientific configuration instead of hardcoded values
        if lat_diff > self.config.max_bounding_box_degrees
            || lon_diff > self.config.max_bounding_box_degrees
        {
            validation_status.execution_mode = ExecutionMode::Failed(format!(
                "Bounding box exceeds maximum: {:.2}° x {:.2}° (max {:.2}° x {:.2}°)",
                lat_diff,
                lon_diff,
                self.config.max_bounding_box_degrees,
                self.config.max_bounding_box_degrees
            ));
            validation_status.processing_time_ms = start_time.elapsed().as_secs_f64() * 1000.0;
            return Err(SarError::InvalidInput(
                format!("Bounding box too large: {:.2}° x {:.2}° (max {:.2}° x {:.2}°) - this may indicate an error in bounding box calculation or multi-scene processing",
                       lat_diff, lon_diff, self.config.max_bounding_box_degrees, self.config.max_bounding_box_degrees)
            ));
        }

        // Warn about large bounding boxes using configuration
        if lat_diff > self.config.warning_bounding_box_degrees
            || lon_diff > self.config.warning_bounding_box_degrees
        {
            let warning_msg = format!("Large bounding box detected: {:.2}° x {:.2}° (warning threshold: {:.2}°) - processing may be slow and require significant memory",
                                    lat_diff, lon_diff, self.config.warning_bounding_box_degrees);
            log::warn!("{}", warning_msg);
            log::warn!("Consider subdividing the processing area or checking if the bounding box is correct");
            validation_status.fallback_reason = Some(warning_msg);
        }

        // Log normal processing info (degrees only; avoid rough km approximations)
        if lat_diff <= self.config.warning_bounding_box_degrees
            && lon_diff <= self.config.warning_bounding_box_degrees
        {
            log::info!("Processing area: {:.2}° x {:.2}°", lat_diff, lon_diff);
        }

        log::debug!(
            "Input bounding box: ({:.6}, {:.6}) to ({:.6}, {:.6})",
            sar_bbox.min_lon,
            sar_bbox.min_lat,
            sar_bbox.max_lon,
            sar_bbox.max_lat
        );
        log::debug!("Bounding box size: {:.6}° x {:.6}°", lon_diff, lat_diff);
        log::debug!(
            "Using configuration limits: max={:.1}°, warning={:.1}°",
            self.config.max_bounding_box_degrees,
            self.config.warning_bounding_box_degrees
        );

        validation_status.processing_time_ms = start_time.elapsed().as_secs_f64() * 1000.0;

        // Scientific mode: Accept bounding box for Range-Doppler terrain correction
        // The actual Range-Doppler geocoding occurs in the pixel-by-pixel projection step
        log::info!("Scientific mode: Using Range-Doppler terrain correction with validated bounds");

        // Return the validated bounding box - the actual Range-Doppler geocoding
        // happens pixel-by-pixel in the main terrain correction loop
        Ok(sar_bbox.clone())
    }

    /// Create output grid from bounds with proper CRS handling
    pub(crate) fn create_output_grid(
        &self,
        bounds: &BoundingBox,
    ) -> SarResult<(usize, usize, GeoTransform)> {
        log::debug!("Creating output grid for CRS EPSG:{}", self.output_crs);
        log::debug!(
            "Input bounds: min_lat={:.6}, max_lat={:.6}, min_lon={:.6}, max_lon={:.6}",
            bounds.min_lat,
            bounds.max_lat,
            bounds.min_lon,
            bounds.max_lon
        );
        log::debug!("Output spacing: {:.2}m", self.output_spacing);

        let (width, height, transform) = if self.output_crs == 4326 {
            // Geographic coordinate system (degrees)
            self.create_geographic_grid(bounds)?
        } else {
            // Projected coordinate system (meters) - UTM, etc.
            self.create_projected_grid(bounds)?
        };

        log::info!("📐 Output grid: {}x{} pixels", width, height);
        log::debug!(
            "GeoTransform: [{:.8}, {:.8}, {:.2}, {:.8}, {:.2}, {:.8}]",
            transform.top_left_x,
            transform.pixel_width,
            transform.rotation_x,
            transform.top_left_y,
            transform.rotation_y,
            transform.pixel_height
        );

        Ok((width, height, transform))
    }

    /// Create grid for geographic coordinate system (WGS84)
    pub(crate) fn create_geographic_grid(
        &self,
        bounds: &BoundingBox,
    ) -> SarResult<(usize, usize, GeoTransform)> {
        // For geographic coordinates, convert output spacing from meters to degrees
        let center_lat = (bounds.min_lat + bounds.max_lat) / 2.0;
        let center_lat_rad = center_lat.to_radians();

        // Accurate degree-to-meter conversion using WGS84 ellipsoid
        let a = WGS84_SEMI_MAJOR_AXIS_M;
        let e2 = crate::constants::geodetic::WGS84_ECCENTRICITY_SQUARED;

        let sin_lat = center_lat_rad.sin();
        let meridional_radius = a * (1.0 - e2) / (1.0 - e2 * sin_lat * sin_lat).powf(1.5);
        let prime_vertical_radius = a / (1.0 - e2 * sin_lat * sin_lat).sqrt();

        let meters_per_degree_lat = meridional_radius * std::f64::consts::PI / 180.0;
        let meters_per_degree_lon =
            prime_vertical_radius * center_lat_rad.cos() * std::f64::consts::PI / 180.0;

        // Convert output spacing from meters to degrees
        let pixel_size_lat = self.output_spacing / meters_per_degree_lat;
        let pixel_size_lon = self.output_spacing / meters_per_degree_lon;

        log::debug!(
            "Geographic pixel sizes: lat={:.8}°, lon={:.8}°",
            pixel_size_lat,
            pixel_size_lon
        );

        // Calculate grid dimensions
        let lat_extent = bounds.max_lat - bounds.min_lat;
        let lon_extent = bounds.max_lon - bounds.min_lon;

        let original_height = (lat_extent / pixel_size_lat).ceil() as usize;
        let original_width = (lon_extent / pixel_size_lon).ceil() as usize;

        // Apply configuration limits
        let max_dimension = self.config.max_output_dimension;
        let height = original_height.min(max_dimension);
        let width = original_width.min(max_dimension);

        // CRITICAL DEBUG: Log grid calculation to verify no row/col swap
        log::debug!("🔍 RUST GRID CALCULATION:");
        log::debug!(
            "  Bbox: lon=[{:.6}, {:.6}], lat=[{:.6}, {:.6}]",
            bounds.min_lon,
            bounds.max_lon,
            bounds.min_lat,
            bounds.max_lat
        );
        log::debug!("  Extents: lon={:.6}°, lat={:.6}°", lon_extent, lat_extent);
        log::debug!(
            "  Pixel sizes: lon={:.8}°/px, lat={:.8}°/px",
            pixel_size_lon,
            pixel_size_lat
        );
        log::debug!(
            "  Grid calc: width = lon_extent/pixel_lon = {:.2}/{:.8} = {}",
            lon_extent,
            pixel_size_lon,
            original_width
        );
        log::debug!(
            "  Grid calc: height = lat_extent/pixel_lat = {:.2}/{:.8} = {}",
            lat_extent,
            pixel_size_lat,
            original_height
        );
        log::debug!("  Output dimensions: width={}, height={}", width, height);
        log::debug!(
            "  → Array will be allocated as: {}×{} (rows × cols = height × width)",
            height,
            width
        );

        if width != original_width || height != original_height {
            log::warn!(
                "Output dimensions clamped: {}x{} -> {}x{}",
                original_width,
                original_height,
                width,
                height
            );
        }

        // Create geotransform for geographic coordinates
        let transform = GeoTransform {
            top_left_x: bounds.min_lon,  // Western longitude
            pixel_width: pixel_size_lon, // Degrees per pixel (east)
            rotation_x: 0.0,
            top_left_y: bounds.max_lat, // Northern latitude
            rotation_y: 0.0,
            pixel_height: -pixel_size_lat, // Negative for north-up orientation
        };

        Ok((width, height, transform))
    }

    /// Create grid for projected coordinate system (UTM, etc.)
    fn create_projected_grid(
        &self,
        bounds: &BoundingBox,
    ) -> SarResult<(usize, usize, GeoTransform)> {
        // CRITICAL FIX: For UTM projections, we must compute ALL FOUR corners
        // and find the actual min/max in projected space.
        //
        // The previous code incorrectly assumed min_x comes from (min_lon, min_lat),
        // but in UTM, easting varies with latitude for the same longitude.
        // For example: lon=8.0 at lat=49.0 has different easting than lon=8.0 at lat=50.4

        // Convert all four corners to projected coordinates
        let (x_ll, y_ll) = self.geographic_to_projected(bounds.min_lon, bounds.min_lat)?; // lower-left
        let (x_lr, y_lr) = self.geographic_to_projected(bounds.max_lon, bounds.min_lat)?; // lower-right
        let (x_ul, y_ul) = self.geographic_to_projected(bounds.min_lon, bounds.max_lat)?; // upper-left
        let (x_ur, y_ur) = self.geographic_to_projected(bounds.max_lon, bounds.max_lat)?; // upper-right

        // Find actual min/max in projected coordinates
        let proj_min_x = x_ll.min(x_lr).min(x_ul).min(x_ur);
        let proj_max_x = x_ll.max(x_lr).max(x_ul).max(x_ur);
        let proj_min_y = y_ll.min(y_lr).min(y_ul).min(y_ur);
        let proj_max_y = y_ll.max(y_lr).max(y_ul).max(y_ur);

        log::debug!(
            "Projected corners: LL=({:.1},{:.1}), LR=({:.1},{:.1}), UL=({:.1},{:.1}), UR=({:.1},{:.1})",
            x_ll, y_ll, x_lr, y_lr, x_ul, y_ul, x_ur, y_ur
        );
        log::debug!(
            "Projected bounds: min_x={:.1}, max_x={:.1}, min_y={:.1}, max_y={:.1}",
            proj_min_x,
            proj_max_x,
            proj_min_y,
            proj_max_y
        );

        // Calculate extent in projected coordinates
        let x_extent = proj_max_x - proj_min_x;
        let y_extent = proj_max_y - proj_min_y;

        // Calculate grid dimensions based on output spacing in meters
        let original_width = (x_extent / self.output_spacing).ceil() as usize;
        let original_height = (y_extent / self.output_spacing).ceil() as usize;

        // Apply configuration limits
        let max_dimension = self.config.max_output_dimension;
        let width = original_width.min(max_dimension);
        let height = original_height.min(max_dimension);

        if width != original_width || height != original_height {
            log::warn!(
                "Output dimensions clamped: {}x{} -> {}x{}",
                original_width,
                original_height,
                width,
                height
            );
        }

        // For projected coordinates (UTM, etc.): use meters directly
        // Note: Geographic coordinates (EPSG:4326) are handled by create_geographic_grid()
        log::debug!(
            "Using meter spacing for projected CRS (EPSG:{}):",
            self.output_crs
        );
        log::debug!("   Pixel spacing: {:.1}m", self.output_spacing);

        let transform = GeoTransform {
            top_left_x: proj_min_x,
            pixel_width: self.output_spacing,
            rotation_x: 0.0,
            top_left_y: proj_max_y, // Start from northern edge
            rotation_y: 0.0,
            pixel_height: -self.output_spacing, // Negative for north-up orientation
        };

        Ok((width, height, transform))
    }
}
