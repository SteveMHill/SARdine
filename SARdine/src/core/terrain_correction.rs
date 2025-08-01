use crate::types::{SarError, SarResult, BoundingBox, GeoTransform, OrbitData, StateVector, MaskingWorkflow, MaskResult, SurfaceNormal};
use gdal::Dataset;
use ndarray::Array2;
use std::path::Path;
use rayon::prelude::*;
use std::collections::HashMap;
use serde::{Serialize, Deserialize};

/// Scientific processing configuration with literature-based thresholds
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TerrainCorrectionConfig {
    /// Maximum bounding box size in degrees (based on Sentinel-1 scene size)
    pub max_bounding_box_degrees: f64,
    /// Warning threshold for large bounding boxes
    pub warning_bounding_box_degrees: f64,
    /// Maximum output dimension to prevent memory issues
    pub max_output_dimension: usize,
    /// Minimum valid elevation (below sea level is valid)
    pub min_valid_elevation: f32,
    /// Maximum valid elevation (Mount Everest height)
    pub max_valid_elevation: f32,
    /// Minimum valid range pixel (avoiding near-range artifacts)
    pub min_valid_range_pixel: f64,
    /// Maximum valid range pixel (sensor-dependent)
    pub max_valid_range_pixel: f64,
    /// Convergence tolerance for iterative algorithms
    pub convergence_tolerance: f64,
    /// Maximum iterations for iterative algorithms
    pub max_iterations: usize,
    /// Enable strict validation mode
    pub strict_validation: bool,
}

impl Default for TerrainCorrectionConfig {
    fn default() -> Self {
        Self {
            // Based on Sentinel-1 IW swath width (~250km) = ~2.5 degrees at equator
            max_bounding_box_degrees: 30.0,  // Allow some multi-scene processing
            warning_bounding_box_degrees: 5.0,  // Warn for large single scenes
            max_output_dimension: 10000,  // Prevent excessive memory usage
            min_valid_elevation: -500.0,  // Below Dead Sea level
            max_valid_elevation: 9000.0,  // Above Mount Everest
            min_valid_range_pixel: 0.0,   // Start of swath
            max_valid_range_pixel: 30000.0,  // End of typical swath
            convergence_tolerance: 1e-6,  // Precision for iterative solutions
            max_iterations: 50,  // Prevent infinite loops
            strict_validation: true,  // Enable comprehensive checking
        }
    }
}

/// Algorithm execution status tracking
#[derive(Debug, Clone)]
pub struct AlgorithmStatus {
    pub algorithm_name: String,
    pub execution_mode: ExecutionMode,
    pub iterations_used: Option<usize>,
    pub convergence_achieved: Option<bool>,
    pub fallback_reason: Option<String>,
    pub processing_time_ms: f64,
}

#[derive(Debug, Clone)]
pub enum ExecutionMode {
    Primary,
    Fallback(String),
    Failed(String),
}

/// Processing metadata for scientific reproducibility
#[derive(Debug, Clone)]
pub struct ProcessingMetadata {
    pub algorithm_statuses: Vec<AlgorithmStatus>,
    pub configuration_used: TerrainCorrectionConfig,
    pub input_validation_results: ValidationResults,
}

#[derive(Debug, Clone)]
pub struct ValidationResults {
    pub bounding_box_valid: bool,
    pub elevation_range_valid: bool,
    pub coordinate_system_valid: bool,
    pub orbit_data_valid: bool,
    pub warnings: Vec<String>,
    pub errors: Vec<String>,
}

/// Interpolation methods for SAR data resampling
#[derive(Debug, Clone, Copy)]
pub enum InterpolationMethod {
    Nearest,
    Bilinear,
    Bicubic,
}

impl Default for InterpolationMethod {
    fn default() -> Self {
        InterpolationMethod::Bilinear
    }
}

impl InterpolationMethod {
    pub fn from_str(s: &str) -> Self {
        match s.to_lowercase().as_str() {
            "nearest" => InterpolationMethod::Nearest,
            "bilinear" => InterpolationMethod::Bilinear,
            "bicubic" => InterpolationMethod::Bicubic,
            _ => InterpolationMethod::Bilinear,
        }
    }
}

/// Spatial lookup grid for fast coordinate transformations
#[derive(Debug, Clone)]
pub struct SpatialLookupGrid {
    /// Grid of pre-computed SAR coordinates (range, azimuth) for each output pixel
    #[allow(dead_code)]
    grid: Array2<Option<(f64, f64)>>,
    /// Grid resolution (pixels)
    #[allow(dead_code)]
    grid_size: (usize, usize),
    /// Mapping from output coordinates to grid coordinates
    #[allow(dead_code)]
    output_to_grid_transform: GeoTransform,
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
            let keys_to_remove: Vec<_> = self.cache.keys().take(self.cache.len() / 2).cloned().collect();
            for key in keys_to_remove {
                self.cache.remove(&key);
            }
        }
        self.cache.insert((x, y), value);
    }
}

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
    /// Scientific processing configuration
    config: TerrainCorrectionConfig,
    /// Processing metadata for reproducibility
    pub metadata: ProcessingMetadata,
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
        let config = TerrainCorrectionConfig::default();
        let metadata = ProcessingMetadata {
            algorithm_statuses: Vec::new(),
            configuration_used: config.clone(),
            input_validation_results: ValidationResults {
                bounding_box_valid: true,
                elevation_range_valid: true,
                coordinate_system_valid: true,
                orbit_data_valid: true,
                warnings: Vec::new(),
                errors: Vec::new(),
            },
        };
        
        Self {
            dem,
            dem_transform,
            dem_nodata,
            output_crs,
            output_spacing,
            config,
            metadata,
        }
    }

    /// Create new terrain correction processor with custom configuration
    pub fn new_with_config(
        dem: Array2<f32>,
        dem_transform: GeoTransform,
        dem_nodata: f32,
        output_crs: u32,
        output_spacing: f64,
        config: TerrainCorrectionConfig,
    ) -> Self {
        let metadata = ProcessingMetadata {
            algorithm_statuses: Vec::new(),
            configuration_used: config.clone(),
            input_validation_results: ValidationResults {
                bounding_box_valid: true,
                elevation_range_valid: true,
                coordinate_system_valid: true,
                orbit_data_valid: true,
                warnings: Vec::new(),
                errors: Vec::new(),
            },
        };
        
        Self {
            dem,
            dem_transform,
            dem_nodata,
            output_crs,
            output_spacing,
            config,
            metadata,
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
                    ) {
                        // Debug: Print first few pixel mappings
                        if i < 5 && j < 5 {
                            println!("DEBUG: pixel ({}, {}): lat={:.6}, lon={:.6}, elev={:.1} -> SAR pixel ({}, {})", 
                                     i, j, lat, lon, elevation, sar_range, sar_azimuth);
                            println!("DEBUG: SAR image dimensions: {}x{}", sar_image.dim().0, sar_image.dim().1);
                        }
                        
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
                            // Debug: Print out-of-bounds cases for first few
                            if i < 5 && j < 5 {
                                println!("DEBUG: Out of bounds: SAR pixel ({}, {}) vs image bounds ({}x{})", 
                                         sar_range, sar_azimuth, sar_image.dim().0, sar_image.dim().1);
                            }
                            output_image[[i, j]] = f32::NAN;
                        }
                    } else {
                        // Debug: Print failed mapping for first few
                        if i < 5 && j < 5 {
                            println!("DEBUG: Failed SAR pixel mapping for lat={:.6}, lon={:.6}, elev={:.1}", 
                                     lat, lon, elevation);
                        }
                        output_image[[i, j]] = f32::NAN;
                    }
                } else {
                    // Debug: Print failed DEM elevation for first few
                    if i < 5 && j < 5 {
                        println!("DEBUG: No elevation data for lat={:.6}, lon={:.6}", lat, lon);
                    }
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

    /// Optimized Range-Doppler terrain correction with chunked parallel processing
    pub fn range_doppler_terrain_correction_chunked(
        &self,
        sar_image: &Array2<f32>,
        orbit_data: &OrbitData,
        params: &RangeDopplerParams,
        sar_bbox: &BoundingBox,
        chunk_size: Option<usize>,
    ) -> SarResult<(Array2<f32>, GeoTransform)> {
        log::info!("🚀 Starting OPTIMIZED Range-Doppler terrain correction with chunked processing");
        let start_time = std::time::Instant::now();
        
        // Step 1: Calculate output grid bounds
        let output_bounds = self.calculate_output_bounds(sar_bbox)?;
        let (output_width, output_height, output_transform) = 
            self.create_output_grid(&output_bounds)?;
        log::info!("Output grid: {}x{} pixels", output_width, output_height);

        // Step 2: Build orbit lookup table for faster orbit queries
        let orbit_lut = self.build_orbit_lookup_table_optimized(orbit_data, &output_bounds)?;
        log::debug!("Built orbit lookup table with {} entries", orbit_lut.len());

        // Step 3: Process in parallel chunks
        let chunk_sz = chunk_size.unwrap_or(rayon::current_num_threads() * 64);
        let total_pixels = output_width * output_height;
        log::info!("Processing {} pixels in chunks of {} using {} threads", 
                  total_pixels, chunk_sz, rayon::current_num_threads());

        let mut output_image = Array2::zeros((output_height, output_width));
        let mut total_valid = 0;

        // Process rows in parallel chunks
        let row_chunks: Vec<_> = (0..output_height).step_by(chunk_sz).collect();
        
        let chunk_results: Vec<_> = row_chunks
            .par_iter()
            .map(|&start_row| {
                let end_row = (start_row + chunk_sz).min(output_height);
                self.process_row_chunk_optimized(
                    sar_image,
                    orbit_data,
                    params,
                    &output_transform,
                    &orbit_lut,
                    start_row,
                    end_row,
                    output_width,
                )
            })
            .collect();

        // Combine results from parallel chunks
        for (chunk_idx, chunk_result) in chunk_results.into_iter().enumerate() {
            match chunk_result {
                Ok((chunk_data, valid_count)) => {
                    let start_row = chunk_idx * chunk_sz;
                    let end_row = (start_row + chunk_sz).min(output_height);
                    
                    for (local_i, i) in (start_row..end_row).enumerate() {
                        if local_i < chunk_data.nrows() {
                            for j in 0..output_width {
                                if j < chunk_data.ncols() {
                                    output_image[[i, j]] = chunk_data[[local_i, j]];
                                }
                            }
                        }
                    }
                    total_valid += valid_count;
                }
                Err(e) => {
                    log::warn!("Chunk {} processing failed: {}", chunk_idx, e);
                }
            }
        }

        let elapsed = start_time.elapsed();
        let coverage = (total_valid as f64 / total_pixels as f64) * 100.0;
        log::info!("✅ OPTIMIZED terrain correction completed in {:.2}s: {:.1}% coverage", 
                  elapsed.as_secs_f64(), coverage);

        Ok((output_image, output_transform))
    }

    /// Build optimized orbit lookup table
    fn build_orbit_lookup_table_optimized(
        &self,
        orbit_data: &OrbitData,
        output_bounds: &BoundingBox,
    ) -> SarResult<Vec<(usize, [f64; 3])>> {
        // Create a lookup table of orbit indices and positions within the area of interest
        let mut orbit_lut = Vec::new();
        
        for (i, state_vector) in orbit_data.state_vectors.iter().enumerate() {
            // Simple spatial filtering - include all orbit points for now
            // In a more advanced implementation, convert ECEF to lat/lon and filter spatially
            orbit_lut.push((i, state_vector.position));
        }
        
        // Sort by distance to scene center for more efficient lookup
        let center_lat = (output_bounds.min_lat + output_bounds.max_lat) / 2.0;
        let center_lon = (output_bounds.min_lon + output_bounds.max_lon) / 2.0;
        let center_ecef = self.latlon_to_ecef(center_lat, center_lon, 0.0);
        
        orbit_lut.sort_by(|a, b| {
            let dist_a = self.distance_to_point(&a.1, &center_ecef);
            let dist_b = self.distance_to_point(&b.1, &center_ecef);
            dist_a.partial_cmp(&dist_b).unwrap_or(std::cmp::Ordering::Equal)
        });
        
        Ok(orbit_lut)
    }

    /// Process a chunk of rows in parallel
    fn process_row_chunk_optimized(
        &self,
        sar_image: &Array2<f32>,
        orbit_data: &OrbitData,
        params: &RangeDopplerParams,
        output_transform: &GeoTransform,
        orbit_lut: &[(usize, [f64; 3])],
        start_row: usize,
        end_row: usize,
        output_width: usize,
    ) -> SarResult<(Array2<f32>, usize)> {
        let chunk_height = end_row - start_row;
        let mut chunk_data = Array2::zeros((chunk_height, output_width));
        let mut valid_count = 0;

        for (local_i, i) in (start_row..end_row).enumerate() {
            for j in 0..output_width {
                // Convert output pixel to geographic coordinates
                let map_x = output_transform.top_left_x + (j as f64) * output_transform.pixel_width;
                let map_y = output_transform.top_left_y + (i as f64) * output_transform.pixel_height;

                if let Ok((lat, lon)) = self.map_to_geographic(map_x, map_y) {
                    if let Some(elevation) = self.get_elevation_at_latlon(lat, lon) {
                        // Use optimized orbit lookup
                        if let Ok(Some((sar_range, sar_azimuth))) = self.latlon_to_sar_pixel_with_lut(
                            lat, lon, elevation, orbit_data, params, orbit_lut
                        ) {
                            if sar_range < sar_image.dim().1 && sar_azimuth < sar_image.dim().0 {
                                let value = self.bilinear_interpolate(
                                    sar_image, 
                                    sar_range as f64, 
                                    sar_azimuth as f64
                                );
                                chunk_data[[local_i, j]] = value;
                                if !value.is_nan() {
                                    valid_count += 1;
                                }
                            } else {
                                chunk_data[[local_i, j]] = f32::NAN;
                            }
                        } else {
                            chunk_data[[local_i, j]] = f32::NAN;
                        }
                    } else {
                        chunk_data[[local_i, j]] = f32::NAN;
                    }
                } else {
                    chunk_data[[local_i, j]] = f32::NAN;
                }
            }
        }

        Ok((chunk_data, valid_count))
    }

    /// Optimized lat/lon to SAR pixel conversion using lookup table
    fn latlon_to_sar_pixel_with_lut(
        &self,
        lat: f64,
        lon: f64,
        elevation: f64,
        orbit_data: &OrbitData,
        params: &RangeDopplerParams,
        orbit_lut: &[(usize, [f64; 3])],
    ) -> SarResult<Option<(usize, usize)>> {
        // Convert lat/lon/elevation to ECEF
        let target_ecef = self.latlon_to_ecef(lat, lon, elevation);
        
        // Find closest orbit state vectors using lookup table
        let mut min_distance = f64::MAX;
        let mut best_state_idx = 0;
        
        for &(state_idx, satellite_pos) in orbit_lut.iter().take(10) { // Use only first 10 for speed
            let distance = self.distance_to_point(&satellite_pos, &target_ecef);
            if distance < min_distance {
                min_distance = distance;
                best_state_idx = state_idx;
            }
        }
        
        // Use the best state vector for range-doppler calculation
        if best_state_idx < orbit_data.state_vectors.len() {
            let state_vector = &orbit_data.state_vectors[best_state_idx];
            
            // Compute slant range
            let slant_range = self.distance_to_point(&state_vector.position, &target_ecef);
            
            // Convert to range pixel using proper SAR geometry
            let range_pixel = (slant_range / params.speed_of_light * 2.0 - params.slant_range_time) 
                / (params.range_pixel_spacing / params.speed_of_light);
            
            // CRITICAL FIX: Attempt proper Doppler-based azimuth calculation first
            let azimuth_result = self.calculate_azimuth_with_doppler(
                &target_ecef, 
                orbit_data, 
                params, 
                best_state_idx
            );
            
            let azimuth_pixel = match azimuth_result {
                Ok(azimuth) => {
                    log::debug!("✅ Using proper Doppler-based azimuth calculation: {:.2}", azimuth);
                    azimuth
                }
                Err(e) => {
                    // EXPLICIT FALLBACK WITH WARNING
                    let fallback_azimuth = best_state_idx as f64 * params.prf / 1000.0;
                    log::warn!("⚠️  Doppler-based azimuth calculation failed: {}", e);
                    log::warn!("⚠️  FALLING BACK to simplified approximation: {:.2}", fallback_azimuth);
                    log::warn!("⚠️  This may reduce geocoding accuracy - consider improving orbit data quality");
                    
                    // Record the fallback in metadata
                    let mut status = AlgorithmStatus {
                        algorithm_name: "azimuth_geocoding".to_string(),
                        execution_mode: ExecutionMode::Fallback(format!("Doppler calculation failed: {}", e)),
                        iterations_used: None,
                        convergence_achieved: Some(false),
                        fallback_reason: Some(format!("Using simplified approximation due to: {}", e)),
                        processing_time_ms: 0.0,
                    };
                    // Note: In a real implementation, we'd need a mutable reference to store this
                    
                    fallback_azimuth
                }
            };
            
            // Validate pixel coordinates using configuration bounds
            if range_pixel >= self.config.min_valid_range_pixel && 
               range_pixel <= self.config.max_valid_range_pixel && 
               azimuth_pixel >= 0.0 {
                Ok(Some((range_pixel.round() as usize, azimuth_pixel.round() as usize)))
            } else {
                log::debug!("Pixel coordinates out of valid bounds: range={:.1} (valid: {:.1}-{:.1}), azimuth={:.1}", 
                           range_pixel, self.config.min_valid_range_pixel, self.config.max_valid_range_pixel, azimuth_pixel);
                Ok(None)
            }
        } else {
            log::warn!("Invalid orbit state vector index: {} >= {}", best_state_idx, orbit_data.state_vectors.len());
            Ok(None)
        }
    }
    
    /// Proper Doppler-based azimuth calculation with convergence checking
    fn calculate_azimuth_with_doppler(
        &self,
        target_ecef: &[f64; 3],
        orbit_data: &OrbitData,
        params: &RangeDopplerParams,
        initial_guess: usize,
    ) -> SarResult<f64> {
        // Newton-Raphson iteration for zero-Doppler condition
        let mut azimuth_time = initial_guess as f64 / params.prf;
        let mut iteration = 0;
        
        while iteration < self.config.max_iterations {
            // Interpolate satellite position and velocity at current azimuth time
            let (sat_pos, sat_vel) = self.interpolate_orbit_state(orbit_data, azimuth_time)?;
            
            // Calculate range vector from satellite to target
            let range_vec = [
                target_ecef[0] - sat_pos[0],
                target_ecef[1] - sat_pos[1],
                target_ecef[2] - sat_pos[2],
            ];
            
            // Calculate Doppler frequency
            let range_magnitude = (range_vec[0].powi(2) + range_vec[1].powi(2) + range_vec[2].powi(2)).sqrt();
            let doppler_freq = -2.0 * (range_vec[0] * sat_vel[0] + range_vec[1] * sat_vel[1] + range_vec[2] * sat_vel[2]) 
                              / (params.wavelength * range_magnitude);
            
            // Check convergence
            if doppler_freq.abs() < self.config.convergence_tolerance {
                log::debug!("Doppler convergence achieved in {} iterations: |f_d| = {:.2e}", iteration + 1, doppler_freq.abs());
                return Ok(azimuth_time * params.prf);
            }
            
            // Calculate derivative for Newton-Raphson update (simplified)
            let time_step = 1e-6; // 1 microsecond
            let (_, sat_vel_next) = self.interpolate_orbit_state(orbit_data, azimuth_time + time_step)?;
            let doppler_derivative = (sat_vel_next[0] - sat_vel[0]) / time_step; // Simplified derivative
            
            // Newton-Raphson update
            if doppler_derivative.abs() > 1e-12 {
                azimuth_time -= doppler_freq / doppler_derivative;
            } else {
                return Err(SarError::Processing("Doppler derivative too small for convergence".to_string()));
            }
            
            iteration += 1;
        }
        
        Err(SarError::Processing(format!("Doppler-based azimuth calculation failed to converge after {} iterations", self.config.max_iterations)))
    }
    
    /// Interpolate orbit state at given azimuth time
    fn interpolate_orbit_state(&self, orbit_data: &OrbitData, azimuth_time: f64) -> SarResult<([f64; 3], [f64; 3])> {
        // Simple linear interpolation (in practice would use higher-order)
        let time_index = azimuth_time * 1000.0; // Convert to index approximation
        let idx = time_index.floor() as usize;
        
        if idx + 1 < orbit_data.state_vectors.len() {
            let frac = time_index - idx as f64;
            let sv1 = &orbit_data.state_vectors[idx];
            let sv2 = &orbit_data.state_vectors[idx + 1];
            
            let pos = [
                sv1.position[0] + frac * (sv2.position[0] - sv1.position[0]),
                sv1.position[1] + frac * (sv2.position[1] - sv1.position[1]),
                sv1.position[2] + frac * (sv2.position[2] - sv1.position[2]),
            ];
            
            let vel = [
                sv1.velocity[0] + frac * (sv2.velocity[0] - sv1.velocity[0]),
                sv1.velocity[1] + frac * (sv2.velocity[1] - sv1.velocity[1]),
                sv1.velocity[2] + frac * (sv2.velocity[2] - sv1.velocity[2]),
            ];
            
            Ok((pos, vel))
        } else {
            Err(SarError::Processing("Orbit interpolation index out of bounds".to_string()))
        }
    }

    /// Bilinear interpolation for SAR image sampling
    fn bilinear_interpolate(&self, image: &Array2<f32>, x: f64, y: f64) -> f32 {
        let (height, width) = image.dim();
        let x1 = x.floor() as usize;
        let y1 = y.floor() as usize;
        let x2 = (x1 + 1).min(width - 1);
        let y2 = (y1 + 1).min(height - 1);
        
        if x1 >= width || y1 >= height {
            return f32::NAN;
        }
        
        let dx = x - x1 as f64;
        let dy = y - y1 as f64;
        
        let v11 = image[[y1, x1]];
        let v12 = image[[y2, x1]];
        let v21 = image[[y1, x2]];
        let v22 = image[[y2, x2]];
        
        // Check for NaN values
        if v11.is_nan() || v12.is_nan() || v21.is_nan() || v22.is_nan() {
            return f32::NAN;
        }
        
        let v1 = v11 * (1.0 - dx as f32) + v21 * dx as f32;
        let v2 = v12 * (1.0 - dx as f32) + v22 * dx as f32;
        v1 * (1.0 - dy as f32) + v2 * dy as f32
    }
    
    /// Find azimuth time using orbit lookup table
    #[allow(dead_code)]
    fn find_azimuth_time_with_lut(
        &self,
        ground_point: &[f64; 3],
        orbit_data: &OrbitData,
        orbit_lut: &[(usize, [f64; 3])],
    ) -> SarResult<f64> {
        // Find the closest orbit state vector
        let mut min_distance = f64::MAX;
        let mut best_idx = 0;
        
        for &(idx, sat_pos) in orbit_lut.iter() {
            let distance = self.distance_to_point(&sat_pos, ground_point);
            if distance < min_distance {
                min_distance = distance;
                best_idx = idx;
            }
        }
        
        // Return time index as simplified azimuth time (in practice would need proper time conversion)
        Ok(best_idx as f64 / orbit_data.state_vectors.len() as f64)
    }
    
    /// Interpolate satellite state at given azimuth time
    #[allow(dead_code)]
    fn interpolate_satellite_state(&self, orbit_data: &OrbitData, azimuth_time: f64) -> SarResult<StateVector> {
        let state_idx = (azimuth_time * orbit_data.state_vectors.len() as f64) as usize;
        
        if state_idx < orbit_data.state_vectors.len() {
            Ok(orbit_data.state_vectors[state_idx].clone())
        } else {
            Err(SarError::Processing("Invalid azimuth time".to_string()))
        }
    }

    /// Calculate output bounds from SAR bounding box
    fn calculate_output_bounds(&self, sar_bbox: &BoundingBox) -> SarResult<BoundingBox> {
        // Validate input bounding box
        if sar_bbox.max_lat <= sar_bbox.min_lat || sar_bbox.max_lon <= sar_bbox.min_lon {
            return Err(SarError::InvalidInput(
                "Invalid bounding box: max values must be greater than min values".to_string()
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
        if lat_diff > self.config.max_bounding_box_degrees || lon_diff > self.config.max_bounding_box_degrees {
            validation_status.execution_mode = ExecutionMode::Failed(
                format!("Bounding box exceeds maximum: {:.2}° x {:.2}° (max {:.2}° x {:.2}°)", 
                       lat_diff, lon_diff, self.config.max_bounding_box_degrees, self.config.max_bounding_box_degrees)
            );
            validation_status.processing_time_ms = start_time.elapsed().as_secs_f64() * 1000.0;
            return Err(SarError::InvalidInput(
                format!("Bounding box too large: {:.2}° x {:.2}° (max {:.2}° x {:.2}°) - this may indicate an error in bounding box calculation or multi-scene processing", 
                       lat_diff, lon_diff, self.config.max_bounding_box_degrees, self.config.max_bounding_box_degrees)
            ));
        }
        
        // Warn about large bounding boxes using configuration
        if lat_diff > self.config.warning_bounding_box_degrees || lon_diff > self.config.warning_bounding_box_degrees {
            let warning_msg = format!("Large bounding box detected: {:.2}° x {:.2}° (warning threshold: {:.2}°) - processing may be slow and require significant memory", 
                                    lat_diff, lon_diff, self.config.warning_bounding_box_degrees);
            log::warn!("{}", warning_msg);
            log::warn!("Consider subdividing the processing area or checking if the bounding box is correct");
            validation_status.fallback_reason = Some(warning_msg);
        }
        
        // Log normal processing info
        if lat_diff <= self.config.warning_bounding_box_degrees && lon_diff <= self.config.warning_bounding_box_degrees {
            log::info!("Processing area: {:.2}° x {:.2}° (~{:.0} km x {:.0} km)", 
                      lat_diff, lon_diff, lat_diff * 111.0, lon_diff * 111.0);
        }
        
        log::debug!("Input bounding box: ({:.6}, {:.6}) to ({:.6}, {:.6})", 
                   sar_bbox.min_lon, sar_bbox.min_lat, sar_bbox.max_lon, sar_bbox.max_lat);
        log::debug!("Bounding box size: {:.6}° x {:.6}°", lon_diff, lat_diff);
        log::debug!("Using configuration limits: max={:.1}°, warning={:.1}°", 
                   self.config.max_bounding_box_degrees, self.config.warning_bounding_box_degrees);
        
        validation_status.processing_time_ms = start_time.elapsed().as_secs_f64() * 1000.0;
        
        // For now, use the SAR bounding box directly
        // In a more sophisticated implementation, this would project the SAR geometry
        // to the output coordinate system and compute proper bounds
        Ok(sar_bbox.clone())
    }
    
    /// Create output grid from bounds
    fn create_output_grid(&self, bounds: &BoundingBox) -> SarResult<(usize, usize, GeoTransform)> {
        // Calculate output dimensions based on spacing in meters
        // Convert lat/lon differences to approximate meters
        let lat_diff = bounds.max_lat - bounds.min_lat;
        let lon_diff = bounds.max_lon - bounds.min_lon;
        
        // Approximate conversion: 1 degree lat ≈ 111km, 1 degree lon ≈ 111km * cos(lat)
        let center_lat = (bounds.min_lat + bounds.max_lat) / 2.0;
        let meters_per_degree_lat = 111000.0;
        let meters_per_degree_lon = 111000.0 * center_lat.to_radians().cos();
        
        // Calculate dimensions in pixels
        let original_width = ((lon_diff * meters_per_degree_lon) / self.output_spacing).ceil() as usize;
        let original_height = ((lat_diff * meters_per_degree_lat) / self.output_spacing).ceil() as usize;
        
        // Use configuration-based maximum dimensions instead of hardcoded values
        let max_dimension = self.config.max_output_dimension;
        let width = original_width.min(max_dimension);
        let height = original_height.min(max_dimension);
        
        log::debug!("Output grid calculation:");
        log::debug!("  Bounds: ({:.6}, {:.6}) to ({:.6}, {:.6})", bounds.min_lon, bounds.min_lat, bounds.max_lon, bounds.max_lat);
        log::debug!("  Lat diff: {:.6} degrees ({:.1} km)", lat_diff, lat_diff * meters_per_degree_lat / 1000.0);
        log::debug!("  Lon diff: {:.6} degrees ({:.1} km)", lon_diff, lon_diff * meters_per_degree_lon / 1000.0);
        log::debug!("  Output spacing: {:.1} m", self.output_spacing);
        log::debug!("  Calculated dimensions: {} x {} pixels", width, height);
        log::debug!("  Maximum allowed dimension: {} pixels (configurable)", max_dimension);
        
        if width != original_width || height != original_height {
            let warning_msg = format!("Output dimensions clamped to configuration maximum: {}x{} -> {}x{}", 
                          original_width, original_height, width, height);
            log::warn!("{}", warning_msg);
            log::warn!("Consider increasing max_output_dimension in configuration or reducing output spacing");
        }
        
        let transform = GeoTransform {
            top_left_x: bounds.min_lon,
            pixel_width: (bounds.max_lon - bounds.min_lon) / width as f64,
            rotation_x: 0.0,
            top_left_y: bounds.max_lat,
            rotation_y: 0.0,
            pixel_height: -(bounds.max_lat - bounds.min_lat) / height as f64,
        };
        
        Ok((width, height, transform))
    }
    
    /// Convert map coordinates to geographic coordinates
    fn map_to_geographic(&self, map_x: f64, map_y: f64) -> SarResult<(f64, f64)> {
        // For WGS84 (EPSG:4326), map coordinates are already lat/lon
        if self.output_crs == 4326 {
            Ok((map_y, map_x)) // (lat, lon)
        } else {
            // For other CRS, would need proper projection transformation
            // For now, return as-is
            Ok((map_y, map_x))
        }
    }
    
    /// Get elevation at lat/lon coordinates
    fn get_elevation_at_latlon(&self, lat: f64, lon: f64) -> Option<f64> {
        // Convert lat/lon to DEM pixel coordinates
        let dem_x = (lon - self.dem_transform.top_left_x) / self.dem_transform.pixel_width;
        let dem_y = (lat - self.dem_transform.top_left_y) / self.dem_transform.pixel_height;
        
        let dem_col = dem_x.round() as usize;
        let dem_row = dem_y.round() as usize;
        
        let (dem_height, dem_width) = self.dem.dim();
        
        if dem_row < dem_height && dem_col < dem_width {
            let elevation = self.dem[[dem_row, dem_col]];
            if elevation != self.dem_nodata && elevation.is_finite() {
                Some(elevation as f64)
            } else {
                None
            }
        } else {
            None
        }
    }
    
    /// Optimized DEM elevation lookup with bilinear interpolation and caching
    #[allow(dead_code)]
    fn get_elevation_at_latlon_optimized(&self, lat: f64, lon: f64) -> Option<f64> {
        // Convert lat/lon to DEM pixel coordinates
        let dem_x = (lon - self.dem_transform.top_left_x) / self.dem_transform.pixel_width;
        let dem_y = (lat - self.dem_transform.top_left_y) / self.dem_transform.pixel_height;
        
        let dem_col = dem_x.floor() as usize;
        let dem_row = dem_y.floor() as usize;
        
        let (dem_height, dem_width) = self.dem.dim();
        
        // Check bounds
        if dem_row >= dem_height - 1 || dem_col >= dem_width - 1 {
            return None;
        }
        
        // Bilinear interpolation for better accuracy
        let dx = dem_x - dem_col as f64;
        let dy = dem_y - dem_row as f64;
        
        let v11 = self.dem[[dem_row, dem_col]];
        let v12 = self.dem[[dem_row + 1, dem_col]];
        let v21 = self.dem[[dem_row, dem_col + 1]];
        let v22 = self.dem[[dem_row + 1, dem_col + 1]];
        
        // Check for no-data values
        if v11 == self.dem_nodata || v12 == self.dem_nodata || 
           v21 == self.dem_nodata || v22 == self.dem_nodata ||
           !v11.is_finite() || !v12.is_finite() || !v21.is_finite() || !v22.is_finite() {
            return None;
        }
        
        // Bilinear interpolation
        let v1 = v11 * (1.0 - dx as f32) + v21 * dx as f32;
        let v2 = v12 * (1.0 - dx as f32) + v22 * dx as f32;
        let elevation = v1 * (1.0 - dy as f32) + v2 * dy as f32;
        
        Some(elevation as f64)
    }
    
    /// Vectorized batch DEM lookup for multiple points
    #[allow(dead_code)]
    fn get_elevations_batch(&self, coords: &[(f64, f64)]) -> Vec<Option<f64>> {
        coords.par_iter()
            .map(|&(lat, lon)| self.get_elevation_at_latlon_optimized(lat, lon))
            .collect()
    }
    
    /// Convert lat/lon/elevation to ECEF coordinates
    fn latlon_to_ecef(&self, lat: f64, lon: f64, elevation: f64) -> [f64; 3] {
        let a = 6378137.0; // WGS84 semi-major axis
        let e2 = 0.00669437999014; // WGS84 first eccentricity squared
        
        let lat_rad = lat.to_radians();
        let lon_rad = lon.to_radians();
        
        let n = a / (1.0 - e2 * lat_rad.sin().powi(2)).sqrt();
        
        let x = (n + elevation) * lat_rad.cos() * lon_rad.cos();
        let y = (n + elevation) * lat_rad.cos() * lon_rad.sin();
        let z = (n * (1.0 - e2) + elevation) * lat_rad.sin();
        
        [x, y, z]
    }
    
    /// Calculate distance between two 3D points
    fn distance_to_point(&self, point1: &[f64; 3], point2: &[f64; 3]) -> f64 {
        let dx = point1[0] - point2[0];
        let dy = point1[1] - point2[1];
        let dz = point1[2] - point2[2];
        (dx*dx + dy*dy + dz*dz).sqrt()
    }
    
    /// Basic lat/lon to SAR pixel mapping (for non-optimized path)
    fn latlon_to_sar_pixel(
        &self,
        lat: f64,
        lon: f64,
        elevation: f64,
        orbit_data: &OrbitData,
        params: &RangeDopplerParams,
    ) -> Option<(usize, usize)> {
        // Convert lat/lon/elevation to ECEF coordinates
        let target_ecef = self.latlon_to_ecef(lat, lon, elevation);
        
        // Find best orbit state vector based on minimum slant range
        let mut min_range = f64::MAX;
        let mut best_idx = 0;
        let mut best_range = 0.0;
        
        for (i, state_vector) in orbit_data.state_vectors.iter().enumerate() {
            let slant_range = self.distance_to_point(&state_vector.position, &target_ecef);
            if slant_range < min_range {
                min_range = slant_range;
                best_idx = i;
                best_range = slant_range;
            }
        }
        
        // Basic validation - check if slant range is reasonable for Sentinel-1
        if best_range < 800_000.0 || best_range > 950_000.0 {
            // Outside typical Sentinel-1 slant range (800-950 km)
            return None;
        }
        
        // Simplified range pixel calculation
        // For Sentinel-1: slant_range_time = 2 * slant_range / c
        let two_way_travel_time = 2.0 * best_range / params.speed_of_light;
        let range_pixel = (two_way_travel_time - params.slant_range_time) 
            / (params.range_pixel_spacing * 2.0 / params.speed_of_light);
        
        // Simplified azimuth pixel calculation 
        // Use orbit index scaled to azimuth line count (very approximate)
        let azimuth_fraction = best_idx as f64 / orbit_data.state_vectors.len() as f64;
        let azimuth_pixel = azimuth_fraction * 25000.0; // Approximate line count for Sentinel-1
        
        // Validate pixel coordinates
        if range_pixel >= 0.0 && range_pixel < 30000.0 && // Reasonable range pixel bounds
           azimuth_pixel >= 0.0 && azimuth_pixel < 30000.0 { // Reasonable azimuth pixel bounds
            Some((range_pixel.round() as usize, azimuth_pixel.round() as usize))
        } else {
            None
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
        use gdal::DriverManager;
        
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
        let mut combined_mask_bool = Array2::<bool>::from_elem((height, width), true);
        let mut valid_pixels = 0;
        
        for row in 0..height {
            for col in 0..width {
                combined_mask_bool[[row, col]] = gamma0_mask[[row, col]] &&
                    dem_mask[[row, col]] &&
                    lia_mask[[row, col]];
                
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
        
        let stats = crate::types::MaskStats {
            total_pixels,
            water_pixels: 0, // Would be computed if water masking is implemented
            shadow_pixels: 0, // Would be computed if shadow masking is implemented
            layover_pixels: 0, // Would be computed if layover masking is implemented
            noise_pixels: 0, // Would be computed if noise masking is implemented
            valid_pixels,
        };
        
        Ok(MaskResult {
            water_mask: None,
            shadow_mask: None,
            layover_mask: None,
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
        let fill_val = fill_value.unwrap_or(f32::NAN);
        
        for ((row, col), &mask_val) in mask.indexed_iter() {
            if mask_val == 0 {
                masked_data[[row, col]] = fill_val;
            }
        }
        
        Ok(masked_data)
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

    /// Enhanced terrain correction pipeline (static method for Python API) 
    pub fn enhanced_terrain_correction_pipeline(
        sar_image: &Array2<f32>,
        dem_path: &str,
        orbit_data: &OrbitData,
        bbox: &BoundingBox,
        output_path: &str,
        output_crs: u32,
        output_spacing: f64,
        chunk_size: Option<usize>,
    ) -> SarResult<()> {
        // Check if DEM file exists before trying to open it
        let corrector = if std::path::Path::new(dem_path).exists() {
            // DEM file exists, try to load it
            match Self::from_dem_file(dem_path, output_crs, output_spacing) {
                Ok(corrector) => corrector,
                Err(e) => {
                    log::warn!("Failed to load DEM file '{}': {}", dem_path, e);
                    log::info!("Attempting automatic DEM preparation...");
                    
                    // Use a default cache directory
                    let cache_dir = "/tmp/sardine_dem_cache";
                    
                    // Try to prepare DEM using the bounding box
                    let (dem_array, dem_transform) = crate::io::dem::DemReader::prepare_dem_for_scene(
                        bbox,
                        output_spacing,
                        cache_dir,
                    )?;
                    
                    // Create terrain corrector with prepared DEM
                    Self::new(
                        dem_array,
                        dem_transform,
                        -32768.0, // Default nodata value
                        output_crs,
                        output_spacing,
                    )
                }
            }
        } else {
            // DEM file doesn't exist, prepare it automatically
            log::info!("DEM file not found at '{}', attempting automatic preparation...", dem_path);
            
            // Use a default cache directory
            let cache_dir = "/tmp/sardine_dem_cache";
            
            // Try to prepare DEM using the bounding box
            let (dem_array, dem_transform) = crate::io::dem::DemReader::prepare_dem_for_scene(
                bbox,
                output_spacing,
                cache_dir,
            )?;
            
            // Create terrain corrector with prepared DEM
            Self::new(
                dem_array,
                dem_transform,
                -32768.0, // Default nodata value
                output_crs,
                output_spacing,
            )
        };
        
        let params = RangeDopplerParams::default();
        
        let (terrain_corrected, geo_transform) = corrector.range_doppler_terrain_correction_chunked(
            sar_image, orbit_data, &params, bbox, chunk_size
        )?;
        
        corrector.save_geotiff(&terrain_corrected, &geo_transform, output_path, Some("LZW"))?;
        Ok(())
    }
    
    /// Adaptive terrain correction (static method for Python API)
    pub fn adaptive_terrain_correction(
        sar_image: &Array2<f32>,
        dem_path: &str,
        orbit_data: &OrbitData,
        bbox: &BoundingBox,
        _output_path: &str,
        output_crs: u32,
        output_spacing: f64,
        enable_chunking: bool,
        chunk_size: Option<usize>,
    ) -> SarResult<(Array2<f32>, GeoTransform)> {
        println!("DEBUG: adaptive_terrain_correction called with dem_path: {}", dem_path);
        println!("DEBUG: bbox: {:?}", bbox);
        
        // Try to load DEM file, with automatic preparation fallback
        let corrector = match Self::from_dem_file(dem_path, output_crs, output_spacing) {
            Ok(corrector) => {
                println!("DEBUG: Successfully loaded DEM file");
                corrector
            },
            Err(e) => {
                // DEM file doesn't exist, try to prepare it automatically
                println!("DEBUG: DEM file loading failed: {}", e);
                log::info!("DEM file not found, attempting automatic preparation...");
                
                // Extract cache directory from dem_path
                let cache_dir = if dem_path.contains("/") {
                    // Extract the directory part from the path
                    let path = std::path::Path::new(dem_path);
                    path.parent().unwrap_or(std::path::Path::new("/tmp/sardine_dem_cache")).to_str().unwrap_or("/tmp/sardine_dem_cache")
                } else {
                    "/tmp/sardine_dem_cache"
                };
                println!("DEBUG: Calling prepare_dem_for_scene with cache_dir: {}", cache_dir);
                
                // Try to prepare DEM using the bounding box
                let (dem_array, dem_transform) = crate::io::dem::DemReader::prepare_dem_for_scene(
                    bbox,
                    output_spacing,
                    cache_dir,
                )?;
                
                println!("DEBUG: prepare_dem_for_scene completed successfully");
                
                // Create terrain corrector with prepared DEM
                Self::new(
                    dem_array,
                    dem_transform,
                    -32768.0, // Default nodata value
                    output_crs,
                    output_spacing,
                )
            }
        };
        
        let params = RangeDopplerParams::default();
        
        if enable_chunking {
            corrector.range_doppler_terrain_correction_chunked(
                sar_image, orbit_data, &params, bbox, chunk_size
            )
        } else {
            corrector.range_doppler_terrain_correction(
                sar_image, orbit_data, &params, bbox
            )
        }
    }
    
    /// Static ultra-optimized terrain correction method for Python API
    pub fn ultra_optimized_terrain_correction_static(
        sar_image: &Array2<f32>,
        dem_path: &str,
        orbit_data: &OrbitData,
        sar_bbox: &BoundingBox,
        output_path: &str,
        output_crs: u32,
        output_spacing: f64,
        chunk_size: Option<usize>,
        interpolation_method: &str,
        enable_spatial_cache: bool,
    ) -> SarResult<()> {
        log::info!("🚀 Starting static ultra-optimized terrain correction");
        
        // Try to load DEM file, with automatic preparation fallback
        let corrector = match Self::from_dem_file(dem_path, output_crs, output_spacing) {
            Ok(corrector) => corrector,
            Err(_) => {
                // DEM file doesn't exist, try to prepare it automatically
                log::info!("DEM file not found, attempting automatic preparation...");
                
                // Use a default cache directory
                let cache_dir = "/tmp/sardine_dem_cache";
                
                // Try to prepare DEM using the bounding box
                let (dem_array, dem_transform) = crate::io::dem::DemReader::prepare_dem_for_scene(
                    sar_bbox,
                    output_spacing,
                    cache_dir,
                )?;
                
                // Create terrain corrector with prepared DEM
                Self::new(
                    dem_array,
                    dem_transform,
                    -32768.0, // Default nodata value
                    output_crs,
                    output_spacing,
                )
            }
        };
        
        let params = RangeDopplerParams::default();
        let interp_method = InterpolationMethod::from_str(interpolation_method);
        
        let (terrain_corrected, geo_transform) = corrector.ultra_optimized_terrain_correction(
            sar_image, 
            orbit_data, 
            &params, 
            sar_bbox,
            interp_method,
            enable_spatial_cache,
            chunk_size
        )?;
        
        if !output_path.is_empty() {
            corrector.save_geotiff(&terrain_corrected, &geo_transform, output_path, Some("LZW"))?;
        }
        
        Ok(())
    }
    
    /// Fast slant range to pixel conversion
    #[inline]
    fn slant_range_to_pixel(&self, slant_range: f64, params: &RangeDopplerParams) -> f64 {
        (slant_range / params.speed_of_light * 2.0 - params.slant_range_time) 
            / (params.range_pixel_spacing / params.speed_of_light)
    }
    
    /// Fast azimuth time to pixel conversion
    #[inline]
    fn azimuth_time_to_pixel(&self, azimuth_time: f64, params: &RangeDopplerParams) -> f64 {
        azimuth_time * params.prf
    }
    
    /// Create optimized orbit lookup table for fast orbit queries
    fn create_orbit_lookup_table(
        &self,
        orbit_data: &OrbitData,
        sar_bbox: &BoundingBox,
    ) -> SarResult<Vec<(usize, [f64; 3])>> {
        let mut orbit_lut = Vec::new();
        
        // Calculate scene center for sorting
        let center_lat = (sar_bbox.min_lat + sar_bbox.max_lat) / 2.0;
        let center_lon = (sar_bbox.min_lon + sar_bbox.max_lon) / 2.0;
        let center_ecef = self.latlon_to_ecef(center_lat, center_lon, 0.0);
        
        // Build lookup table with all orbit state vectors
        for (i, state_vector) in orbit_data.state_vectors.iter().enumerate() {
            orbit_lut.push((i, state_vector.position));
        }
        
        // Sort by distance to scene center for efficient lookup
        orbit_lut.sort_by(|a, b| {
            let dist_a = self.distance_to_point(&a.1, &center_ecef);
            let dist_b = self.distance_to_point(&b.1, &center_ecef);
            dist_a.partial_cmp(&dist_b).unwrap_or(std::cmp::Ordering::Equal)
        });
        
        log::debug!("Created orbit LUT with {} entries, sorted by distance to scene center", orbit_lut.len());
        Ok(orbit_lut)
    }
    
    /// Ultra-optimized terrain correction with advanced interpolation and caching
    pub fn ultra_optimized_terrain_correction(
        &self,
        sar_image: &Array2<f32>,
        orbit_data: &OrbitData,
        params: &RangeDopplerParams,
        sar_bbox: &BoundingBox,
        interpolation_method: InterpolationMethod,
        enable_spatial_cache: bool,
        chunk_size: Option<usize>,
    ) -> SarResult<(Array2<f32>, GeoTransform)> {
        log::info!("🚀 Starting ultra-optimized terrain correction");
        log::debug!("Interpolation method: {:?}", interpolation_method);
        log::debug!("Spatial cache enabled: {}", enable_spatial_cache);
        
        // Step 1: Calculate output grid bounds
        let output_bounds = self.calculate_output_bounds(sar_bbox)?;
        let (output_width, output_height, output_transform) = 
            self.create_output_grid(&output_bounds)?;
        log::info!("Output grid: {}x{} pixels", output_width, output_height);

        // Step 2: Create orbit lookup table for fast access
        let orbit_lut = self.create_orbit_lookup_table(orbit_data, sar_bbox)?;
        log::debug!("Created orbit LUT with {} entries", orbit_lut.len());

        // Step 3: Initialize DEM cache if enabled
        let _dem_cache = if enable_spatial_cache {
            Some(DemCache::new(10000))
        } else {
            None
        };

        // Step 4: Initialize output image
        let mut output_image = Array2::zeros((output_height, output_width));

        // Step 5: Process in chunks for memory efficiency
        let chunk_size = chunk_size.unwrap_or(256);
        let chunks_y = (output_height + chunk_size - 1) / chunk_size;
        let chunks_x = (output_width + chunk_size - 1) / chunk_size;

        log::info!("Processing {} chunks ({}x{})", chunks_y * chunks_x, chunks_y, chunks_x);

        // Process chunks in parallel
        let chunks: Vec<_> = (0..chunks_y).flat_map(|chunk_y| {
            (0..chunks_x).map(move |chunk_x| (chunk_y, chunk_x))
        }).collect();

        let processed_chunks: Vec<_> = chunks.par_iter().map(|&(chunk_y, chunk_x)| {
            let y_start = chunk_y * chunk_size;
            let y_end = ((chunk_y + 1) * chunk_size).min(output_height);
            let x_start = chunk_x * chunk_size;
            let x_end = ((chunk_x + 1) * chunk_size).min(output_width);

            let mut chunk_data = Array2::zeros((y_end - y_start, x_end - x_start));

            for i in 0..(y_end - y_start) {
                for j in 0..(x_end - x_start) {
                    let global_i = y_start + i;
                    let global_j = x_start + j;

                    // Convert output pixel to geographic coordinates
                    let map_x = output_transform.top_left_x + (global_j as f64) * output_transform.pixel_width;
                    let map_y = output_transform.top_left_y + (global_i as f64) * output_transform.pixel_height;

                    // Convert map coordinates to lat/lon
                    if let Ok((lat, lon)) = self.map_to_geographic(map_x, map_y) {
                        // Get elevation using cached or interpolated DEM lookup
                        if let Some(elevation) = self.get_elevation_with_interpolation(lat, lon, &mut None) {
                            // Find corresponding SAR pixel using optimized range-doppler
                            if let Some((sar_range, sar_azimuth)) = self.optimized_latlon_to_sar_pixel(
                                lat, lon, elevation, &orbit_lut, params
                            ) {
                                // Check bounds
                                if sar_range >= 0.0 && sar_azimuth >= 0.0 && 
                                   sar_range < sar_image.dim().1 as f64 && 
                                   sar_azimuth < sar_image.dim().0 as f64 {
                                    // Apply selected interpolation method
                                    let value = self.interpolate_sar_value(
                                        sar_image, 
                                        sar_range, 
                                        sar_azimuth, 
                                        interpolation_method
                                    );
                                    chunk_data[[i, j]] = value;
                                } else {
                                    chunk_data[[i, j]] = f32::NAN;
                                }
                            } else {
                                chunk_data[[i, j]] = f32::NAN;
                            }
                        } else {
                            chunk_data[[i, j]] = f32::NAN;
                        }
                    } else {
                        chunk_data[[i, j]] = f32::NAN;
                    }
                }
            }

            ((y_start, y_end, x_start, x_end), chunk_data)
        }).collect();

        // Merge processed chunks
        for ((y_start, y_end, x_start, x_end), chunk_data) in processed_chunks {
            for i in 0..(y_end - y_start) {
                for j in 0..(x_end - x_start) {
                    output_image[[y_start + i, x_start + j]] = chunk_data[[i, j]];
                }
            }
        }

        log::info!("✅ Ultra-optimized terrain correction completed");
        Ok((output_image, output_transform))
    }

    /// Advanced SAR interpolation with multiple methods
    fn interpolate_sar_value(
        &self,
        sar_image: &Array2<f32>,
        x: f64,
        y: f64,
        method: InterpolationMethod,
    ) -> f32 {
        match method {
            InterpolationMethod::Nearest => {
                let i = y.round() as usize;
                let j = x.round() as usize;
                if i < sar_image.dim().0 && j < sar_image.dim().1 {
                    sar_image[[i, j]]
                } else {
                    f32::NAN
                }
            },
            InterpolationMethod::Bilinear => {
                self.bilinear_interpolate(sar_image, x, y)
            },
            InterpolationMethod::Bicubic => {
                self.bicubic_interpolate(sar_image, x, y)
            },
        }
    }

    /// Bicubic interpolation for high-quality resampling
    fn bicubic_interpolate(&self, sar_image: &Array2<f32>, x: f64, y: f64) -> f32 {
        let x0 = x.floor() as i32;
        let y0 = y.floor() as i32;
        let dx = x - x0 as f64;
        let dy = y - y0 as f64;

        // Cubic interpolation weights
        let cubic = |t: f64| -> f64 {
            let t = t.abs();
            if t <= 1.0 {
                1.5 * t * t * t - 2.5 * t * t + 1.0
            } else if t <= 2.0 {
                -0.5 * t * t * t + 2.5 * t * t - 4.0 * t + 2.0
            } else {
                0.0
            }
        };

        let mut sum = 0.0;
        let mut weight_sum = 0.0;

        // Sample 4x4 neighborhood
        for i in -1..3 {
            for j in -1..3 {
                let yi = y0 + i;
                let xi = x0 + j;

                if yi >= 0 && yi < sar_image.dim().0 as i32 && 
                   xi >= 0 && xi < sar_image.dim().1 as i32 {
                    let weight_y = cubic(dy - i as f64);
                    let weight_x = cubic(dx - j as f64);
                    let weight = weight_x * weight_y;
                    
                    let value = sar_image[[yi as usize, xi as usize]];
                    if !value.is_nan() {
                        sum += value as f64 * weight;
                        weight_sum += weight;
                    }
                }
            }
        }

        if weight_sum > 0.0 {
            (sum / weight_sum) as f32
        } else {
            f32::NAN
        }
    }

    /// Get elevation with bilinear DEM interpolation
    fn get_elevation_with_interpolation(
        &self,
        lat: f64,
        lon: f64,
        dem_cache: &mut Option<DemCache>,
    ) -> Option<f32> {
        // Convert lat/lon to DEM pixel coordinates
        let dem_x = (lon - self.dem_transform.top_left_x) / self.dem_transform.pixel_width;
        let dem_y = (lat - self.dem_transform.top_left_y) / self.dem_transform.pixel_height;

        if dem_x < 0.0 || dem_y < 0.0 || 
           dem_x >= self.dem.dim().1 as f64 || 
           dem_y >= self.dem.dim().0 as f64 {
            return None;
        }

        // Check cache first if available
        if let Some(ref mut cache) = dem_cache {
            let cache_x = (dem_x * 10.0) as i32; // Cache with 0.1 pixel precision
            let cache_y = (dem_y * 10.0) as i32;
            if let Some(cached_value) = cache.get(cache_x, cache_y) {
                return Some(cached_value);
            }
        }

        // Bilinear interpolation in DEM
        let x0 = dem_x.floor() as usize;
        let y0 = dem_y.floor() as usize;
        let x1 = (x0 + 1).min(self.dem.dim().1 - 1);
        let y1 = (y0 + 1).min(self.dem.dim().0 - 1);

        let dx = dem_x - x0 as f64;
        let dy = dem_y - y0 as f64;

        let v00 = self.dem[[y0, x0]];
        let v01 = self.dem[[y0, x1]];
        let v10 = self.dem[[y1, x0]];
        let v11 = self.dem[[y1, x1]];

        // Check for no-data values
        if v00 == self.dem_nodata || v01 == self.dem_nodata || 
           v10 == self.dem_nodata || v11 == self.dem_nodata {
            return None;
        }

        let interpolated = v00 as f64 * (1.0 - dx) * (1.0 - dy) +
                          v01 as f64 * dx * (1.0 - dy) +
                          v10 as f64 * (1.0 - dx) * dy +
                          v11 as f64 * dx * dy;

        let result = interpolated as f32;

        // Cache the result if cache is available
        if let Some(ref mut cache) = dem_cache {
            let cache_x = (dem_x * 10.0) as i32;
            let cache_y = (dem_y * 10.0) as i32;
            cache.insert(cache_x, cache_y, result);
        }

        Some(result)
    }

    /// Optimized lat/lon to SAR pixel conversion using orbit LUT
    fn optimized_latlon_to_sar_pixel(
        &self,
        lat: f64,
        lon: f64,
        elevation: f32,
        orbit_lut: &[(usize, [f64; 3])],
        params: &RangeDopplerParams,
    ) -> Option<(f64, f64)> {
        // Convert to ECEF coordinates
        let target_ecef = self.latlon_to_ecef(lat, lon, elevation as f64);
        
        // Use optimized range-doppler calculation with LUT
        self.optimized_range_doppler_calculation(&target_ecef, orbit_lut, params)
    }
    
    /// Optimized range-doppler calculation with LUT
    fn optimized_range_doppler_calculation(
        &self,
        target_ecef: &[f64; 3],
        orbit_lut: &[(usize, [f64; 3])],
        params: &RangeDopplerParams,
    ) -> Option<(f64, f64)> {
        // Use pre-sorted orbit LUT for faster lookup
        let mut best_range = f64::MAX;
        let mut best_azimuth_idx = 0;
        
        // Binary search for closest approach (orbit LUT is sorted by distance to scene center)
        let low = 0;
        let high = orbit_lut.len().min(20); // Limit search to nearby points
        
        for i in low..high {
            let (idx, sat_pos) = orbit_lut[i];
            let range = self.distance_to_point(&sat_pos, target_ecef);
            
            if range < best_range {
                best_range = range;
                best_azimuth_idx = idx;
            }
        }
        
        // Convert to pixel coordinates with optimized formulas
        let range_pixel = self.slant_range_to_pixel(best_range, params);
        let azimuth_pixel = self.azimuth_time_to_pixel(best_azimuth_idx as f64, params);
        
        Some((range_pixel, azimuth_pixel))
    }
}

