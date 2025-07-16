use crate::types::{SarError, SarResult, BoundingBox, GeoTransform, OrbitData, StateVector};
use gdal::{Dataset, DriverManager, Metadata};
use gdal::raster::Buffer;
use ndarray::Array2;
use std::path::Path;

/// Simple 3D vector for geometric calculations
#[derive(Debug, Clone)]
pub struct Vector3 {
    pub x: f64,
    pub y: f64,
    pub z: f64,
}

/// Terrain correction processor for SAR geocoding
pub struct TerrainCorrector {
    /// Digital Elevation Model data
    pub dem: Array2<f32>,  // Made public for Python API access
    /// DEM geotransform
    dem_transform: GeoTransform,
    /// DEM no-data value
    dem_nodata: f32,
    /// Output coordinate reference system (EPSG code)
    output_crs: u32,
    /// Output pixel spacing in meters
    output_spacing: f64,
}

/// Range-Doppler terrain correction parameters
#[derive(Debug, Clone)]
pub struct RangeDopplerParams {
    /// Range pixel spacing in meters
    pub range_pixel_spacing: f64,
    /// Azimuth pixel spacing in meters  
    pub azimuth_pixel_spacing: f64,
    /// Slant range time to first pixel (seconds)
    pub slant_range_time: f64,
    /// Pulse repetition frequency (Hz)
    pub prf: f64,
    /// Radar wavelength (meters)
    pub wavelength: f64,
    /// Speed of light (m/s)
    pub speed_of_light: f64,
}

impl Default for RangeDopplerParams {
    fn default() -> Self {
        Self {
            range_pixel_spacing: 2.33,  // Sentinel-1 IW typical
            azimuth_pixel_spacing: 14.1, // Sentinel-1 IW typical
            slant_range_time: 5.44e-3,   // Sentinel-1 IW typical
            prf: 486.5,                  // Sentinel-1 IW typical
            wavelength: 0.0555,          // C-band
            speed_of_light: 299_792_458.0,
        }
    }
}

/// Ground point with coordinates and elevation
#[derive(Debug, Clone)]
pub struct GroundPoint {
    pub latitude: f64,
    pub longitude: f64,
    pub elevation: f64,
    pub x: f64,  // Projected X coordinate
    pub y: f64,  // Projected Y coordinate
}

/// Geocoding result for a single pixel
#[derive(Debug, Clone)]
pub struct GeocodedPixel {
    pub sar_range: usize,
    pub sar_azimuth: usize,
    pub ground_point: GroundPoint,
    pub value: f32,
    pub valid: bool,
}

impl TerrainCorrector {
    /// Create new terrain correction processor
    pub fn new(
        dem: Array2<f32>,
        dem_transform: GeoTransform,
        dem_nodata: f32,
        output_crs: u32,
        output_spacing: f64,
    ) -> Self {
        Self {
            dem,
            dem_transform,
            dem_nodata,
            output_crs,
            output_spacing,
        }
    }

    /// Load DEM from file for terrain correction
    pub fn from_dem_file<P: AsRef<Path>>(
        dem_path: P,
        output_crs: u32,
        output_spacing: f64,
    ) -> SarResult<Self> {
        log::info!("Loading DEM for terrain correction: {}", dem_path.as_ref().display());
        
        let dataset = Dataset::open(dem_path.as_ref())?;
        let geo_transform = dataset.geo_transform()?;
        let (width, height) = dataset.raster_size();
        
        // Read DEM data
        let rasterband = dataset.rasterband(1)?;
        let nodata_value = rasterband.no_data_value().unwrap_or(-32768.0) as f32;
        let band_data = rasterband.read_as::<f32>((0, 0), (width, height), (width, height), None)?;
        
        let dem_array = Array2::from_shape_vec((height, width), band_data.data)
            .map_err(|e| SarError::Processing(format!("Failed to reshape DEM data: {}", e)))?;

        let dem_transform_struct = GeoTransform {
            top_left_x: geo_transform[0],
            pixel_width: geo_transform[1],
            rotation_x: geo_transform[2],
            top_left_y: geo_transform[3],
            rotation_y: geo_transform[4],
            pixel_height: geo_transform[5],
        };

        Ok(Self::new(
            dem_array,
            dem_transform_struct,
            nodata_value,
            output_crs,
            output_spacing,
        ))
    }

    /// Perform Range-Doppler terrain correction
    pub fn range_doppler_terrain_correction(
        &self,
        sar_image: &Array2<f32>,
        orbit_data: &OrbitData,
        params: &RangeDopplerParams,
        sar_bbox: &BoundingBox,
    ) -> SarResult<(Array2<f32>, GeoTransform)> {
        log::info!("🗺️  Starting Range-Doppler terrain correction");
        log::debug!("SAR image shape: {:?}", sar_image.dim());
        log::debug!("Output CRS: EPSG:{}", self.output_crs);
        log::debug!("Output spacing: {:.2}m", self.output_spacing);

        // Step 1: Calculate output grid bounds
        let output_bounds = self.calculate_output_bounds(sar_bbox)?;
        log::debug!("Output bounds: {:?}", output_bounds);

        // Step 2: Create output grid
        let (output_width, output_height, output_transform) = 
            self.create_output_grid(&output_bounds)?;
        log::info!("Output grid: {}x{} pixels", output_width, output_height);

        // Step 3: Initialize output image
        let mut output_image = Array2::zeros((output_height, output_width));
        let mut valid_count = 0;

        // Step 4: Backward geocoding - for each output pixel, find corresponding SAR pixel
        for i in 0..output_height {
            for j in 0..output_width {
                // Convert output pixel to geographic coordinates
                let map_x = output_transform.top_left_x + (j as f64) * output_transform.pixel_width;
                let map_y = output_transform.top_left_y + (i as f64) * output_transform.pixel_height;

                // Convert map coordinates to lat/lon
                let (lat, lon) = self.map_to_geographic(map_x, map_y)?;

                // Get elevation from DEM
                if let Some(elevation) = self.get_elevation_at_latlon(lat, lon) {
                    // Find corresponding SAR pixel using range-doppler equations
                    if let Some((sar_range, sar_azimuth)) = self.latlon_to_sar_pixel(
                        lat, lon, elevation, orbit_data, params
                    )? {
                        // Check if SAR pixel is within image bounds
                        if sar_range < sar_image.dim().1 && sar_azimuth < sar_image.dim().0 {
                            // Bilinear interpolation from SAR image
                            let value = self.bilinear_interpolate(
                                sar_image, 
                                sar_range as f64, 
                                sar_azimuth as f64
                            );
                            output_image[[i, j]] = value;
                            valid_count += 1;
                        } else {
                            output_image[[i, j]] = f32::NAN;
                        }
                    } else {
                        output_image[[i, j]] = f32::NAN;
                    }
                } else {
                    output_image[[i, j]] = f32::NAN;
                }
            }

            // Progress reporting
            if i % (output_height / 10).max(1) == 0 {
                let progress = (i as f64 / output_height as f64) * 100.0;
                log::info!("Terrain correction progress: {:.1}%", progress);
            }
        }

        let coverage = (valid_count as f64 / (output_width * output_height) as f64) * 100.0;
        log::info!("✅ Terrain correction completed: {:.1}% coverage", coverage);

        Ok((output_image, output_transform))
    }

    /// Calculate output bounds in the target projection
    fn calculate_output_bounds(&self, sar_bbox: &BoundingBox) -> SarResult<BoundingBox> {
        // For simplicity, assume output CRS is WGS84 (EPSG:4326)
        // In a full implementation, you'd reproject the SAR bbox to the target CRS
        let buffer = 0.01; // degrees
        Ok(BoundingBox {
            min_lat: sar_bbox.min_lat - buffer,
            max_lat: sar_bbox.max_lat + buffer,
            min_lon: sar_bbox.min_lon - buffer,
            max_lon: sar_bbox.max_lon + buffer,
        })
    }

    /// Create output grid with specified spacing
    fn create_output_grid(&self, bounds: &BoundingBox) -> SarResult<(usize, usize, GeoTransform)> {
        // Convert spacing from meters to degrees (approximate)
        let lat_center = (bounds.min_lat + bounds.max_lat) / 2.0;
        let meters_per_degree_lat = 111_320.0; // meters per degree latitude
        let meters_per_degree_lon = 111_320.0 * lat_center.to_radians().cos(); // longitude varies with latitude

        let pixel_size_lat = self.output_spacing / meters_per_degree_lat;
        let pixel_size_lon = self.output_spacing / meters_per_degree_lon;

        let width = ((bounds.max_lon - bounds.min_lon) / pixel_size_lon).ceil() as usize;
        let height = ((bounds.max_lat - bounds.min_lat) / pixel_size_lat).ceil() as usize;

        let transform = GeoTransform {
            top_left_x: bounds.min_lon,
            pixel_width: pixel_size_lon,
            rotation_x: 0.0,
            top_left_y: bounds.max_lat,
            rotation_y: 0.0,
            pixel_height: -pixel_size_lat, // Negative for north-up images
        };

        Ok((width, height, transform))
    }

    /// Convert map coordinates to geographic coordinates
    fn map_to_geographic(&self, map_x: f64, map_y: f64) -> SarResult<(f64, f64)> {
        // For EPSG:4326, map coordinates are already lat/lon
        if self.output_crs == 4326 {
            Ok((map_y, map_x)) // (lat, lon)
        } else {
            // For other CRS, would need proper reprojection using PROJ library
            // For now, simplified implementation
            Err(SarError::Processing(format!(
                "Reprojection from EPSG:{} not yet implemented", self.output_crs
            )))
        }
    }

    /// Get elevation at lat/lon from DEM using bilinear interpolation
    fn get_elevation_at_latlon(&self, lat: f64, lon: f64) -> Option<f64> {
        // Convert lat/lon to DEM pixel coordinates
        let col = (lon - self.dem_transform.top_left_x) / self.dem_transform.pixel_width;
        let row = (lat - self.dem_transform.top_left_y) / self.dem_transform.pixel_height;

        // Check bounds
        if col < 0.0 || row < 0.0 || 
           col >= (self.dem.dim().1 - 1) as f64 || 
           row >= (self.dem.dim().0 - 1) as f64 {
            return None;
        }

        // Bilinear interpolation
        let x1 = col.floor() as usize;
        let y1 = row.floor() as usize;
        let x2 = (x1 + 1).min(self.dem.dim().1 - 1);
        let y2 = (y1 + 1).min(self.dem.dim().0 - 1);

        let dx = col - x1 as f64;
        let dy = row - y1 as f64;

        let v11 = self.dem[[y1, x1]] as f64;
        let v12 = self.dem[[y2, x1]] as f64;
        let v21 = self.dem[[y1, x2]] as f64;
        let v22 = self.dem[[y2, x2]] as f64;

        // Check for no-data values
        if v11 == self.dem_nodata as f64 || v12 == self.dem_nodata as f64 ||
           v21 == self.dem_nodata as f64 || v22 == self.dem_nodata as f64 {
            return None;
        }

        let interpolated = v11 * (1.0 - dx) * (1.0 - dy) +
                          v21 * dx * (1.0 - dy) +
                          v12 * (1.0 - dx) * dy +
                          v22 * dx * dy;

        Some(interpolated)
    }

    /// Convert lat/lon/elevation to SAR pixel coordinates using range-doppler equations
    fn latlon_to_sar_pixel(
        &self,
        lat: f64,
        lon: f64,
        elevation: f64,
        orbit_data: &OrbitData,
        params: &RangeDopplerParams,
    ) -> SarResult<Option<(usize, usize)>> {
        // Convert lat/lon/elevation to ECEF coordinates
        let ground_point = self.latlon_to_ecef(lat, lon, elevation);

        // Find the azimuth time when the satellite was closest to the ground point
        let azimuth_time = self.find_azimuth_time(&ground_point, orbit_data)?;

        // Get satellite position and velocity at azimuth time
        let sat_state = self.interpolate_satellite_state(orbit_data, azimuth_time)?;

        // Calculate slant range
        let range_vector = [
            ground_point[0] - sat_state.position[0],
            ground_point[1] - sat_state.position[1],
            ground_point[2] - sat_state.position[2],
        ];
        let slant_range = (range_vector[0].powi(2) + range_vector[1].powi(2) + range_vector[2].powi(2)).sqrt();

        // Convert slant range to range pixel
        let range_time = slant_range / params.speed_of_light * 2.0; // Two-way travel time
        let range_pixel = (range_time - params.slant_range_time) / 
                         (params.range_pixel_spacing / params.speed_of_light * 2.0);

        // Convert azimuth time to azimuth pixel
        let azimuth_pixel = azimuth_time * params.prf;

        if range_pixel >= 0.0 && azimuth_pixel >= 0.0 {
            Ok(Some((range_pixel as usize, azimuth_pixel as usize)))
        } else {
            Ok(None)
        }
    }

    /// Convert lat/lon/elevation to ECEF coordinates
    fn latlon_to_ecef(&self, lat: f64, lon: f64, elevation: f64) -> [f64; 3] {
        let a = 6_378_137.0; // WGS84 semi-major axis
        let e2 = 0.00669437999014; // WGS84 first eccentricity squared

        let lat_rad = lat.to_radians();
        let lon_rad = lon.to_radians();

        let n = a / (1.0 - e2 * lat_rad.sin().powi(2)).sqrt();

        let x = (n + elevation) * lat_rad.cos() * lon_rad.cos();
        let y = (n + elevation) * lat_rad.cos() * lon_rad.sin();
        let z = (n * (1.0 - e2) + elevation) * lat_rad.sin();

        [x, y, z]
    }

    /// Find azimuth time when satellite was closest to ground point
    fn find_azimuth_time(&self, ground_point: &[f64; 3], orbit_data: &OrbitData) -> SarResult<f64> {
        let mut min_distance = f64::INFINITY;
        let mut best_time = 0.0;

        for (i, state_vector) in orbit_data.state_vectors.iter().enumerate() {
            let distance = self.distance_to_point(&state_vector.position, ground_point);
            
            if distance < min_distance {
                min_distance = distance;
                // Convert to relative time from start of orbit data
                best_time = i as f64; // Simplified - should use actual time differences
            }
        }

        Ok(best_time)
    }

    /// Calculate distance between two 3D points
    fn distance_to_point(&self, point1: &[f64; 3], point2: &[f64; 3]) -> f64 {
        ((point1[0] - point2[0]).powi(2) + 
         (point1[1] - point2[1]).powi(2) + 
         (point1[2] - point2[2]).powi(2)).sqrt()
    }

    /// Interpolate satellite state at given time
    fn interpolate_satellite_state(&self, orbit_data: &OrbitData, time: f64) -> SarResult<StateVector> {
        // Simplified linear interpolation between orbit state vectors
        let index = time as usize;
        
        if index >= orbit_data.state_vectors.len() - 1 {
            return Ok(orbit_data.state_vectors.last().unwrap().clone());
        }

        let sv1 = &orbit_data.state_vectors[index];
        let sv2 = &orbit_data.state_vectors[index + 1];
        
        let fraction = time - index as f64;

        let interpolated = StateVector {
            time: sv1.time, // Simplified
            position: [
                sv1.position[0] + fraction * (sv2.position[0] - sv1.position[0]),
                sv1.position[1] + fraction * (sv2.position[1] - sv1.position[1]),
                sv1.position[2] + fraction * (sv2.position[2] - sv1.position[2]),
            ],
            velocity: [
                sv1.velocity[0] + fraction * (sv2.velocity[0] - sv1.velocity[0]),
                sv1.velocity[1] + fraction * (sv2.velocity[1] - sv1.velocity[1]),
                sv1.velocity[2] + fraction * (sv2.velocity[2] - sv1.velocity[2]),
            ],
        };

        Ok(interpolated)
    }

    /// Bilinear interpolation from SAR image
    fn bilinear_interpolate(&self, image: &Array2<f32>, x: f64, y: f64) -> f32 {
        let x1 = x.floor() as usize;
        let y1 = y.floor() as usize;
        let x2 = (x1 + 1).min(image.dim().1 - 1);
        let y2 = (y1 + 1).min(image.dim().0 - 1);

        if x1 >= image.dim().1 || y1 >= image.dim().0 {
            return f32::NAN;
        }

        let dx = x - x1 as f64;
        let dy = y - y1 as f64;

        let v11 = image[[y1, x1]] as f64;
        let v12 = image[[y2, x1]] as f64;
        let v21 = image[[y1, x2]] as f64;
        let v22 = image[[y2, x2]] as f64;

        let interpolated = v11 * (1.0 - dx) * (1.0 - dy) +
                          v21 * dx * (1.0 - dy) +
                          v12 * (1.0 - dx) * dy +
                          v22 * dx * dy;

        interpolated as f32
    }

    /// Save geocoded result as GeoTIFF
    pub fn save_geotiff<P: AsRef<Path>>(
        &self,
        image: &Array2<f32>,
        transform: &GeoTransform,
        output_path: P,
        compression: Option<&str>,
    ) -> SarResult<()> {
        log::info!("Saving geocoded result as GeoTIFF: {}", output_path.as_ref().display());

        let driver = DriverManager::get_driver_by_name("GTiff")?;
        let (height, width) = image.dim();

        let mut dataset = driver.create_with_band_type::<f32, _>(
            output_path.as_ref(),
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
        let crs = format!("EPSG:{}", self.output_crs);
        dataset.set_spatial_ref(&gdal::spatial_ref::SpatialRef::from_epsg(self.output_crs)?)?;

        // Write image data
        let mut rasterband = dataset.rasterband(1)?;
        let flat_data: Vec<f32> = image.iter().cloned().collect();
        let buffer = gdal::raster::Buffer::new((width, height), flat_data);
        rasterband.write((0, 0), (width, height), &buffer)?;

        // Set no-data value
        rasterband.set_no_data_value(Some(f32::NAN as f64))?;

        // Set compression if specified
        if let Some(compression_type) = compression {
            dataset.set_metadata_item("COMPRESS", compression_type, "")?;
        }

        log::info!("✅ GeoTIFF saved successfully");
        Ok(())
    }

    /// Complete terrain correction pipeline
    pub fn complete_terrain_correction_pipeline(
        sar_image: &Array2<f32>,
        dem_path: &str,
        orbit_data: &OrbitData,
        sar_bbox: &BoundingBox,
        output_path: &str,
        output_crs: u32,
        output_spacing: f64,
    ) -> SarResult<()> {
        log::info!("🌍 Starting complete terrain correction pipeline");
        
        // Step 1: Load terrain corrector with DEM
        let corrector = Self::from_dem_file(dem_path, output_crs, output_spacing)?;
        
        // Step 2: Set up range-doppler parameters
        let params = RangeDopplerParams::default();
        
        // Step 3: Perform terrain correction
        let (geocoded_image, geo_transform) = corrector.range_doppler_terrain_correction(
            sar_image,
            orbit_data,
            &params,
            sar_bbox,
        )?;
        
        // Step 4: Save result as GeoTIFF
        corrector.save_geotiff(
            &geocoded_image,
            &geo_transform,
            output_path,
            Some("LZW"), // Use LZW compression
        )?;
        
        log::info!("🎉 Terrain correction pipeline completed successfully!");
        log::info!("Output saved to: {}", output_path);
        
        Ok(())
    }

    /// Enhanced terrain correction pipeline with integrated masking
    pub fn enhanced_terrain_correction_pipeline(
        sar_image: &Array2<f32>,
        dem_path: &str,
        orbit_data: &OrbitData,
        sar_bbox: &BoundingBox,
        output_path: &str,
        output_crs: u32,
        output_spacing: f64,
        masking_config: Option<&MaskingWorkflow>,
        save_intermediate: bool,
    ) -> SarResult<()> {
        log::info!("🌍 Starting enhanced terrain correction pipeline with masking");
        
        // Step 1: Load terrain corrector with DEM
        let corrector = Self::from_dem_file(dem_path, output_crs, output_spacing)?;
        
        // Step 2: Set up range-doppler parameters
        let params = RangeDopplerParams::default();
        
        // Step 3: Perform terrain correction
        let (mut geocoded_image, geo_transform) = corrector.range_doppler_terrain_correction(
            sar_image,
            orbit_data,
            &params,
            sar_bbox,
        )?;
        
        // Step 4: Apply masking workflow if provided
        if let Some(workflow) = masking_config {
            log::info!("🎭 Applying masking workflow");
            
            // Apply masking to the geocoded image using the corrector's DEM
            let mask_result = corrector.apply_masking_workflow(
                &geocoded_image,
                &corrector.dem,
                workflow,
            )?;
            
            log::info!("📊 Masking statistics:");
            log::info!("  - Valid pixels: {}/{} ({:.1}%)", 
                mask_result.valid_pixels, 
                mask_result.total_pixels, 
                mask_result.coverage_percent
            );
            
            // Apply the combined mask to the geocoded image
            geocoded_image = corrector.apply_mask_to_gamma0(
                &geocoded_image,
                &mask_result.combined_mask,
                Some(f32::NAN),
            )?;
            
            // Save intermediate masking products if requested
            if save_intermediate {
                let base_path = Path::new(output_path).parent().unwrap_or(Path::new("."));
                let stem = Path::new(output_path).file_stem().unwrap().to_str().unwrap();
                
                // Save LIA cosine
                let lia_path = base_path.join(format!("{}_lia_cosine.tif", stem));
                corrector.save_geotiff(
                    &mask_result.lia_cosine,
                    &geo_transform,
                    lia_path.to_str().unwrap(),
                    Some("LZW"),
                )?;
                
                // Save combined mask as byte data
                let mask_u8: Array2<u8> = mask_result.combined_mask.mapv(|v| if v { 1u8 } else { 0u8 });
                corrector.save_mask_geotiff(
                    &mask_u8,
                    &geo_transform,
                    &base_path.join(format!("{}_combined_mask.tif", stem)),
                )?;
                
                log::info!("💾 Saved intermediate masking products");
            }
        }
        
        // Step 5: Save final result as GeoTIFF
        corrector.save_geotiff(
            &geocoded_image,
            &geo_transform,
            output_path,
            Some("LZW"), // Use LZW compression
        )?;
        
        log::info!("🎉 Enhanced terrain correction pipeline completed successfully!");
        log::info!("Output saved to: {}", output_path);
        
        Ok(())
    }

    /// Adaptive terrain correction with quality assessment
    pub fn adaptive_terrain_correction(
        sar_image: &Array2<f32>,
        dem_path: &str,
        orbit_data: &OrbitData,
        sar_bbox: &BoundingBox,
        output_path: &str,
        output_crs: u32,
        output_spacing: f64,
        adaptive_thresholds: bool,
    ) -> SarResult<(Array2<f32>, MaskResult, GeoTransform)> {
        log::info!("🧠 Starting adaptive terrain correction with quality assessment");
        
        // Load terrain corrector
        let corrector = Self::from_dem_file(dem_path, output_crs, output_spacing)?;
        let params = RangeDopplerParams::default();
        
        // Perform initial terrain correction
        let (geocoded_image, geo_transform) = corrector.range_doppler_terrain_correction(
            sar_image,
            orbit_data,
            &params,
            sar_bbox,
        )?;
        
        // Create adaptive masking workflow
        let workflow = if adaptive_thresholds {
            corrector.create_adaptive_masking_workflow(&geocoded_image)?
        } else {
            MaskingWorkflow::default()
        };
        
        // Apply masking
        let mask_result = corrector.apply_masking_workflow(
            &geocoded_image,
            &corrector.dem,
            &workflow,
        )?;
        
        // Quality assessment
        corrector.assess_terrain_correction_quality(&mask_result)?;
        
        Ok((geocoded_image, mask_result, geo_transform))
    }

    /// Create adaptive masking workflow based on data characteristics
    pub fn create_adaptive_masking_workflow(&self, gamma0_data: &Array2<f32>) -> SarResult<MaskingWorkflow> {
        log::info!("🔧 Creating adaptive masking thresholds");
        
        // Compute gamma0 statistics for adaptive thresholding
        let valid_gamma0: Vec<f32> = gamma0_data.iter()
            .filter(|&&x| x.is_finite() && x > 0.0)
            .cloned()
            .collect();
        
        if valid_gamma0.is_empty() {
            return Ok(MaskingWorkflow::default());
        }
        
        // Calculate percentiles for adaptive thresholding
        let mut sorted_gamma0 = valid_gamma0.clone();
        sorted_gamma0.sort_by(|a, b| a.partial_cmp(b).unwrap());
        
        let p1_idx = (sorted_gamma0.len() as f32 * 0.01) as usize;
        let p99_idx = (sorted_gamma0.len() as f32 * 0.99) as usize;
        
        let gamma0_min = sorted_gamma0[p1_idx.min(sorted_gamma0.len() - 1)];
        let gamma0_max = sorted_gamma0[p99_idx.min(sorted_gamma0.len() - 1)];
        
        // Adaptive DEM threshold based on elevation range
        let valid_dem: Vec<f32> = self.dem.iter()
            .filter(|&&x| x.is_finite() && x != self.dem_nodata)
            .cloned()
            .collect();
        
        let dem_threshold = if !valid_dem.is_empty() {
            let min_elevation = valid_dem.iter().fold(f32::INFINITY, |a, &b| a.min(b));
            (min_elevation - 50.0).max(-500.0) // Allow 50m below minimum elevation
        } else {
            -100.0 // Default threshold
        };
        
        log::info!("📏 Adaptive thresholds: γ₀ [{:.3}, {:.3}], DEM > {:.1}m", 
            gamma0_min, gamma0_max, dem_threshold);
        
        Ok(MaskingWorkflow {
            lia_threshold: 0.1, // Keep conservative LIA threshold
            dem_threshold: dem_threshold as f64,
            gamma0_min,
            gamma0_max,
        })
    }

    /// Assess terrain correction quality using masking results
    pub fn assess_terrain_correction_quality(&self, mask_result: &MaskResult) -> SarResult<()> {
        log::info!("📈 Assessing terrain correction quality");
        
        let coverage = mask_result.coverage_percent;
        
        // Quality assessment based on coverage
        let quality_level = if coverage >= 90.0 {
            "Excellent"
        } else if coverage >= 75.0 {
            "Good"
        } else if coverage >= 50.0 {
            "Fair"
        } else if coverage >= 25.0 {
            "Poor"
        } else {
            "Very Poor"
        };
        
        log::info!("🏆 Quality Assessment:");
        log::info!("  - Coverage: {:.1}% ({} quality)", coverage, quality_level);
        log::info!("  - Valid pixels: {}/{}", mask_result.valid_pixels, mask_result.total_pixels);
        
        // Check specific mask components
        let gamma0_coverage = mask_result.gamma0_mask.iter().filter(|&&x| x).count() as f64 
            / mask_result.total_pixels as f64 * 100.0;
        let dem_coverage = mask_result.dem_mask.iter().filter(|&&x| x).count() as f64 
            / mask_result.total_pixels as f64 * 100.0;
        let lia_coverage = mask_result.lia_mask.iter().filter(|&&x| x).count() as f64 
            / mask_result.total_pixels as f64 * 100.0;
        
        log::info!("📊 Component Coverage:");
        log::info!("  - Gamma0 valid: {:.1}%", gamma0_coverage);
        log::info!("  - DEM valid: {:.1}%", dem_coverage);
        log::info!("  - LIA valid: {:.1}%", lia_coverage);
        
        if coverage < 50.0 {
            log::warn!("⚠️  Low coverage detected. Consider adjusting masking thresholds.");
        }
        
        Ok(())
    }

    /// Save mask data as GeoTIFF (byte format)
    pub fn save_mask_geotiff<P: AsRef<Path>>(
        &self,
        mask_data: &Array2<u8>,
        transform: &GeoTransform,
        output_path: P,
    ) -> SarResult<()> {
        let (height, width) = mask_data.dim();
        
        // Create output dataset
        let driver = DriverManager::get_driver_by_name("GTiff")?;
        let mut dataset = driver.create_with_band_type::<u8, _>(
            output_path.as_ref(),
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
        
        // Write mask data
        let mut rasterband = dataset.rasterband(1)?;
        let flat_data: Vec<u8> = mask_data.iter().cloned().collect();
        let buffer = gdal::raster::Buffer::new((width, height), flat_data);
        rasterband.write((0, 0), (width, height), &buffer)?;
        
        // Set no-data value for mask (0 = invalid, 1 = valid)
        rasterband.set_no_data_value(Some(255.0))?;
        
        Ok(())
    }

    /// Compute surface normal from DEM using central differences
    pub fn compute_surface_normal(&self, dem_array: &Array2<f32>, row: usize, col: usize) -> SurfaceNormal {
        let (height, width) = dem_array.dim();
        let pixel_spacing = self.output_spacing;
        
        // Central differences with boundary handling
        let dz_dx = if col == 0 {
            dem_array[[row, col + 1]] - dem_array[[row, col]]
        } else if col == width - 1 {
            dem_array[[row, col]] - dem_array[[row, col - 1]]
        } else {
            (dem_array[[row, col + 1]] - dem_array[[row, col - 1]]) / 2.0
        };
        
        let dz_dy = if row == 0 {
            dem_array[[row + 1, col]] - dem_array[[row, col]]
        } else if row == height - 1 {
            dem_array[[row, col]] - dem_array[[row - 1, col]]
        } else {
            (dem_array[[row + 1, col]] - dem_array[[row - 1, col]]) / 2.0
        };
        
        // Compute normal vector: cross product of tangent vectors
        // Tangent vector in x: (pixel_spacing, 0, dz_dx)
        // Tangent vector in y: (0, pixel_spacing, dz_dy)
        // Normal = tx × ty = (-dz_dx, -dz_dy, pixel_spacing^2)
        let mut normal = SurfaceNormal {
            x: -dz_dx as f64 / pixel_spacing,
            y: -dz_dy as f64 / pixel_spacing,
            z: pixel_spacing,
        };
        normal.normalize();
        normal
    }

    /// Compute local incidence angle cosine
    pub fn compute_local_incidence_angle(&self, surface_normal: &SurfaceNormal, radar_look_vector: &Vector3) -> f64 {
        // Local incidence angle is angle between radar look vector and surface normal
        // cos(θ_lia) = |look_vector · surface_normal|
        let dot_product = radar_look_vector.x * surface_normal.x +
                         radar_look_vector.y * surface_normal.y +
                         radar_look_vector.z * surface_normal.z;
        dot_product.abs()
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
                gamma0_mask[[row, col]] = gamma0_val.is_finite() &&
                    gamma0_val >= workflow.gamma0_min &&
                    gamma0_val <= workflow.gamma0_max;
                
                // Check DEM validity
                dem_mask[[row, col]] = dem_val.is_finite() && dem_val > workflow.dem_threshold as f32;
                
                // Compute local incidence angle if DEM is valid
                if dem_mask[[row, col]] {
                    let surface_normal = self.compute_surface_normal(dem_array, row, col);
                    
                    // Approximate radar look vector from target to sensor
                    // For now, use a simplified vertical look assumption
                    // In practice, this should use actual SAR geometry
                    let look_vector = Vector3 {
                        x: 0.0,
                        y: 0.0,
                        z: 1.0, // Vertical look (nadir)
                    };
                    
                    let cos_lia = self.compute_local_incidence_angle(&surface_normal, &look_vector);
                    lia_cosine[[row, col]] = cos_lia as f32;
                    
                    // Check local incidence angle threshold
                    lia_mask[[row, col]] = cos_lia >= workflow.lia_threshold;
                } else {
                    lia_cosine[[row, col]] = f32::NAN;
                    lia_mask[[row, col]] = false;
                }
            }
        }
        
        // Combine all masks
        let mut combined_mask = Array2::<bool>::from_elem((height, width), true);
        let mut valid_pixels = 0;
        
        for row in 0..height {
            for col in 0..width {
                combined_mask[[row, col]] = gamma0_mask[[row, col]] &&
                    dem_mask[[row, col]] &&
                    lia_mask[[row, col]];
                
                if combined_mask[[row, col]] {
                    valid_pixels += 1;
                }
            }
        }
        
        let total_pixels = height * width;
        let coverage_percent = (valid_pixels as f64 / total_pixels as f64) * 100.0;
        
        Ok(MaskResult {
            combined_mask,
            lia_cosine,
            gamma0_mask,
            dem_mask,
            lia_mask,
            valid_pixels,
            total_pixels,
            coverage_percent,
        })
    }

    /// Apply mask to gamma0 data
    pub fn apply_mask_to_gamma0(
        &self,
        gamma0_data: &Array2<f32>,
        mask: &Array2<bool>,
        fill_value: Option<f32>,
    ) -> SarResult<Array2<f32>> {
        let mut masked_data = gamma0_data.clone();
        let fill_val = fill_value.unwrap_or(f32::NAN);
        
        for ((row, col), &is_valid) in mask.indexed_iter() {
            if !is_valid {
                masked_data[[row, col]] = fill_val;
            }
        }
        
        Ok(masked_data)
    }
}

/// Masking workflow for post-terrain correction processing
#[derive(Debug, Clone)]
pub struct MaskingWorkflow {
    /// Local incidence angle threshold (cosine)
    pub lia_threshold: f64,
    /// DEM validity threshold
    pub dem_threshold: f64,
    /// Gamma0 validity thresholds
    pub gamma0_min: f32,
    pub gamma0_max: f32,
}

impl Default for MaskingWorkflow {
    fn default() -> Self {
        Self {
            lia_threshold: 0.1,     // cos(84°) - steep slopes
            dem_threshold: -100.0,   // Below sea level threshold
            gamma0_min: -50.0,      // Minimum gamma0 in dB
            gamma0_max: 10.0,       // Maximum gamma0 in dB
        }
    }
}

/// Combined mask result
#[derive(Debug, Clone)]
pub struct MaskResult {
    /// Combined validity mask (true = valid pixel)
    pub combined_mask: Array2<bool>,
    /// Local incidence angle cosine values
    pub lia_cosine: Array2<f32>,
    /// Individual mask components
    pub gamma0_mask: Array2<bool>,
    pub dem_mask: Array2<bool>,
    pub lia_mask: Array2<bool>,
    /// Statistics
    pub valid_pixels: usize,
    pub total_pixels: usize,
    pub coverage_percent: f64,
}

/// Surface normal vector for local incidence angle calculation
#[derive(Debug, Clone)]
pub struct SurfaceNormal {
    pub x: f64,
    pub y: f64,
    pub z: f64,
}

impl SurfaceNormal {
    /// Normalize the surface normal vector
    pub fn normalize(&mut self) {
        let magnitude = (self.x * self.x + self.y * self.y + self.z * self.z).sqrt();
        if magnitude > 0.0 {
            self.x /= magnitude;
            self.y /= magnitude;
            self.z /= magnitude;
        }
    }
}
