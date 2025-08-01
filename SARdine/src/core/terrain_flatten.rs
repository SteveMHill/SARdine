use crate::types::{SarError, SarResult, OrbitData};
use ndarray::Array2;
use rayon::prelude::*;
use std::f32::consts::PI;

/// Parameters for terrain flattening computation
#[derive(Debug, Clone)]
pub struct TerrainFlatteningParams {
    /// DEM pixel spacing in meters (range, azimuth)
    pub dem_pixel_spacing: (f64, f64),
    /// SAR pixel spacing in meters (range, azimuth)
    pub sar_pixel_spacing: (f64, f64),
    /// Wavelength in meters (C-band ~= 0.055m)
    pub wavelength: f64,
    /// Whether to apply layover/shadow masking
    pub apply_masking: bool,
    /// Minimum valid local incidence angle (degrees)
    pub min_incidence_angle: f32,
    /// Maximum valid local incidence angle (degrees)
    pub max_incidence_angle: f32,
    /// Chunk size for parallel processing
    pub chunk_size: usize,
    /// Enable parallel processing
    pub enable_parallel: bool,
}

impl Default for TerrainFlatteningParams {
    fn default() -> Self {
        Self {
            dem_pixel_spacing: (30.0, 30.0),  // SRTM 30m
            sar_pixel_spacing: (2.3, 14.0),   // Typical Sentinel-1 IW
            wavelength: 0.055,                 // C-band wavelength
            apply_masking: true,
            min_incidence_angle: 10.0,         // Avoid steep angles
            max_incidence_angle: 80.0,         // Avoid grazing angles
            chunk_size: 64,                   // Default chunk size for parallel processing
            enable_parallel: true,             // Enable parallel processing by default
        }
    }
}

/// Terrain flattening processor
pub struct TerrainFlattener {
    params: TerrainFlatteningParams,
    orbit_data: OrbitData,
}

impl TerrainFlattener {
    /// Create a new terrain flattening processor
    pub fn new(params: TerrainFlatteningParams, orbit_data: OrbitData) -> Self {
        Self {
            params,
            orbit_data,
        }
    }

    /// Create a standard terrain flattening processor
    pub fn standard(orbit_data: OrbitData) -> Self {
        Self::new(TerrainFlatteningParams::default(), orbit_data)
    }

    /// Compute slope and aspect from DEM with parallel processing
    /// 
    /// Returns (slope_radians, aspect_radians)
    pub fn compute_slope_aspect(&self, dem: &Array2<f32>) -> SarResult<(Array2<f32>, Array2<f32>)> {
        let (rows, cols) = dem.dim();
        let mut slope = Array2::<f32>::zeros((rows, cols));
        let mut aspect = Array2::<f32>::zeros((rows, cols));

        let dx_scale = self.params.dem_pixel_spacing.0 as f32;
        let dy_scale = self.params.dem_pixel_spacing.1 as f32;

        if self.params.enable_parallel && rows > self.params.chunk_size {
            // Parallel chunked processing
            self.compute_slope_aspect_parallel(dem, &mut slope, &mut aspect, dx_scale, dy_scale)?;
        } else {
            // Sequential processing for small arrays
            self.compute_slope_aspect_sequential(dem, &mut slope, &mut aspect, dx_scale, dy_scale);
        }

        // Handle edges by copying nearest valid values
        self.fill_edges(&mut slope)?;
        self.fill_edges(&mut aspect)?;

        Ok((slope, aspect))
    }

    /// Sequential slope/aspect computation
    fn compute_slope_aspect_sequential(
        &self,
        dem: &Array2<f32>,
        slope: &mut Array2<f32>,
        aspect: &mut Array2<f32>,
        dx_scale: f32,
        dy_scale: f32,
    ) {
        let (rows, cols) = dem.dim();
        
        for i in 1..rows-1 {
            for j in 1..cols-1 {
                // Calculate gradients using central differences
                let dz_dx = (dem[[i, j+1]] - dem[[i, j-1]]) / (2.0 * dx_scale);
                let dz_dy = (dem[[i+1, j]] - dem[[i-1, j]]) / (2.0 * dy_scale);

                // Compute slope (in radians)
                slope[[i, j]] = (dz_dx * dz_dx + dz_dy * dz_dy).sqrt().atan();

                // Compute aspect (in radians, 0 = North, clockwise positive)
                aspect[[i, j]] = (-dz_dy).atan2(dz_dx);
            }
        }
    }

    /// Parallel slope/aspect computation using rayon
    fn compute_slope_aspect_parallel(
        &self,
        dem: &Array2<f32>,
        slope: &mut Array2<f32>,
        aspect: &mut Array2<f32>,
        dx_scale: f32,
        dy_scale: f32,
    ) -> SarResult<()> {
        let (rows, cols) = dem.dim();
        let chunk_size = self.params.chunk_size;
        
        // Collect all row indices to process
        let row_indices: Vec<usize> = (1..rows-1).collect();
        
        // Process in parallel chunks, collecting results
        let results: Vec<_> = row_indices.par_chunks(chunk_size)
            .map(|row_chunk| {
                let mut local_slope_data = Vec::new();
                let mut local_aspect_data = Vec::new();
                
                for &i in row_chunk {
                    for j in 1..cols-1 {
                        // Calculate gradients using central differences
                        let dz_dx = (dem[[i, j+1]] - dem[[i, j-1]]) / (2.0 * dx_scale);
                        let dz_dy = (dem[[i+1, j]] - dem[[i-1, j]]) / (2.0 * dy_scale);

                        // Compute slope (in radians)
                        let slope_val = (dz_dx * dz_dx + dz_dy * dz_dy).sqrt().atan();
                        
                        // Compute aspect (in radians, 0 = North, clockwise positive)
                        let aspect_val = (-dz_dy).atan2(dz_dx);

                        local_slope_data.push((i, j, slope_val));
                        local_aspect_data.push((i, j, aspect_val));
                    }
                }
                (local_slope_data, local_aspect_data)
            })
            .collect();
        
        // Apply results back to arrays
        for (slope_data, aspect_data) in results {
            for (i, j, val) in slope_data {
                slope[[i, j]] = val;
            }
            for (i, j, val) in aspect_data {
                aspect[[i, j]] = val;
            }
        }

        Ok(())
    }

    /// Compute surface normal vectors from slope and aspect with parallel processing
    pub fn compute_surface_normals(
        &self, 
        slope: &Array2<f32>, 
        aspect: &Array2<f32>
    ) -> Array2<[f32; 3]> {
        let (rows, cols) = slope.dim();
        let mut normals = Array2::<[f32; 3]>::from_elem((rows, cols), [0.0, 0.0, 1.0]);

        if self.params.enable_parallel && rows > self.params.chunk_size {
            // Parallel processing - collect results then apply
            let row_indices: Vec<usize> = (0..rows).collect();
            let results: Vec<_> = row_indices.par_chunks(self.params.chunk_size)
                .map(|row_chunk| {
                    let mut local_normals = Vec::new();
                    for &i in row_chunk {
                        for j in 0..cols {
                            let slope_rad = slope[[i, j]];
                            let aspect_rad = aspect[[i, j]];

                            // Convert slope/aspect to 3D normal vector
                            // Normal vector points upward from surface
                            let nx = -slope_rad.sin() * aspect_rad.sin();
                            let ny = slope_rad.sin() * aspect_rad.cos();
                            let nz = slope_rad.cos();

                            local_normals.push((i, j, [nx, ny, nz]));
                        }
                    }
                    local_normals
                })
                .collect();
            
            // Apply results
            for normal_data in results {
                for (i, j, normal) in normal_data {
                    normals[[i, j]] = normal;
                }
            }
        } else {
            // Sequential processing for small arrays
            for i in 0..rows {
                for j in 0..cols {
                    let slope_rad = slope[[i, j]];
                    let aspect_rad = aspect[[i, j]];

                    // Convert slope/aspect to 3D normal vector
                    // Normal vector points upward from surface
                    let nx = -slope_rad.sin() * aspect_rad.sin();
                    let ny = slope_rad.sin() * aspect_rad.cos();
                    let nz = slope_rad.cos();

                    normals[[i, j]] = [nx, ny, nz];
                }
            }
        }

        normals
    }

    /// Compute radar look vectors for each pixel
    /// 
    /// This computes realistic look vectors using satellite position and ground coordinates
    pub fn compute_radar_look_vectors(
        &self,
        sar_dims: (usize, usize),
        range_time: f64,
        azimuth_time_start: f64,
        azimuth_time_spacing: f64,
    ) -> SarResult<Array2<[f32; 3]>> {
        let (rows, cols) = sar_dims;
        let mut look_vectors = Array2::<[f32; 3]>::from_elem((rows, cols), [0.0, 0.0, 0.0]);

        // Get orbit state vectors
        let state_vectors = &self.orbit_data.state_vectors;
        if state_vectors.is_empty() {
            return Err(SarError::Processing("No orbit state vectors available".to_string()));
        }

        log::debug!("Computing look vectors for {}x{} pixels", rows, cols);

        for i in 0..rows {
            let azimuth_time = azimuth_time_start + (i as f64) * azimuth_time_spacing;
            
            // Interpolate satellite position at this azimuth time
            let sat_pos = self.interpolate_satellite_position(azimuth_time)?;

            for j in 0..cols {
                // Calculate slant range distance for this pixel
                let slant_range = self.calculate_slant_range(range_time, j)?;
                
                // Compute ground position (simplified - assuming ellipsoid earth)
                let ground_pos = self.compute_ground_position(sat_pos, slant_range, j as f64)?;
                
                // Look vector points from satellite to ground
                let dx = ground_pos[0] - sat_pos[0];
                let dy = ground_pos[1] - sat_pos[1]; 
                let dz = ground_pos[2] - sat_pos[2];
                
                // Normalize the look vector
                let length = (dx*dx + dy*dy + dz*dz).sqrt();
                if length > 1e-6 {
                    look_vectors[[i, j]] = [
                        (dx / length) as f32,
                        (dy / length) as f32, 
                        (dz / length) as f32
                    ];
                } else {
                    // Fallback for degenerate cases
                    look_vectors[[i, j]] = [0.0, 0.0, -1.0];
                }
            }
        }

        Ok(look_vectors)
    }

    /// Calculate slant range for a given range sample
    fn calculate_slant_range(&self, range_time: f64, range_pixel: usize) -> SarResult<f64> {
        // Speed of light
        const C: f64 = 299792458.0;
        
        // Range time to first pixel plus pixel offset
        let pixel_range_time = range_time + (range_pixel as f64) * self.params.sar_pixel_spacing.0 / C;
        
        // Slant range = (time * speed_of_light) / 2
        Ok(pixel_range_time * C / 2.0)
    }

    /// Compute approximate ground position for satellite position and slant range
    fn compute_ground_position(&self, sat_pos: [f64; 3], slant_range: f64, _range_pixel: f64) -> SarResult<[f64; 3]> {
        // This is a simplified calculation - for production use, proper
        // SAR geometry and earth ellipsoid models should be used
        
        // Earth radius (simplified spherical earth)
        const EARTH_RADIUS: f64 = 6371000.0;
        
        // Satellite altitude (approximate)
        let sat_altitude = (sat_pos[0]*sat_pos[0] + sat_pos[1]*sat_pos[1] + sat_pos[2]*sat_pos[2]).sqrt() - EARTH_RADIUS;
        
        // For simplification, assume ground point is directly below satellite
        // offset by range pixel spacing
        let _ground_distance = (slant_range*slant_range - sat_altitude*sat_altitude).sqrt();
        
        // Normalize satellite position to get direction to ground
        let sat_length = (sat_pos[0]*sat_pos[0] + sat_pos[1]*sat_pos[1] + sat_pos[2]*sat_pos[2]).sqrt();
        let sat_unit = [sat_pos[0]/sat_length, sat_pos[1]/sat_length, sat_pos[2]/sat_length];
        
        // Ground position (simplified)
        let ground_pos = [
            sat_unit[0] * EARTH_RADIUS,
            sat_unit[1] * EARTH_RADIUS, 
            sat_unit[2] * EARTH_RADIUS
        ];
        
        Ok(ground_pos)
    }

    /// Compute local incidence angle with parallel processing
    /// 
    /// θ_lia = arccos(normal · look_vector)
    pub fn compute_local_incidence_angle(
        &self,
        surface_normals: &Array2<[f32; 3]>,
        look_vectors: &Array2<[f32; 3]>,
    ) -> Array2<f32> {
        let (rows, cols) = surface_normals.dim();
        let mut incidence_angles = Array2::<f32>::zeros((rows, cols));

        if self.params.enable_parallel && rows > self.params.chunk_size {
            // Parallel processing - collect results then apply
            let row_indices: Vec<usize> = (0..rows).collect();
            let results: Vec<_> = row_indices.par_chunks(self.params.chunk_size)
                .map(|row_chunk| {
                    let mut local_angles = Vec::new();
                    for &i in row_chunk {
                        for j in 0..cols {
                            let normal = surface_normals[[i, j]];
                            let look = look_vectors[[i, j]];

                            // Compute dot product
                            let dot_product = normal[0] * look[0] + normal[1] * look[1] + normal[2] * look[2];
                            
                            // Clamp to valid range to avoid numerical issues
                            let dot_clamped = dot_product.clamp(-1.0, 1.0);
                            
                            // Compute local incidence angle in radians
                            let theta_lia = dot_clamped.abs().acos();
                            
                            local_angles.push((i, j, theta_lia));
                        }
                    }
                    local_angles
                })
                .collect();
            
            // Apply results
            for angle_data in results {
                for (i, j, angle) in angle_data {
                    incidence_angles[[i, j]] = angle;
                }
            }
        } else {
            // Sequential processing for small arrays
            for i in 0..rows {
                for j in 0..cols {
                    let normal = surface_normals[[i, j]];
                    let look = look_vectors[[i, j]];

                    // Compute dot product
                    let dot_product = normal[0] * look[0] + normal[1] * look[1] + normal[2] * look[2];
                    
                    // Clamp to valid range to avoid numerical issues
                    let dot_clamped = dot_product.clamp(-1.0, 1.0);
                    
                    // Compute local incidence angle in radians
                    let theta_lia = dot_clamped.abs().acos();
                    
                    incidence_angles[[i, j]] = theta_lia;
                }
            }
        }

        incidence_angles
    }

    /// Apply terrain flattening with parallel processing: gamma0 = sigma0 / cos(theta_lia)
    pub fn apply_terrain_flattening(
        &self,
        sigma0: &Array2<f32>,
        local_incidence_angles: &Array2<f32>,
    ) -> SarResult<Array2<f32>> {
        if sigma0.dim() != local_incidence_angles.dim() {
            return Err(SarError::Processing(
                "Sigma0 and incidence angle arrays must have same dimensions".to_string()
            ));
        }

        let (rows, cols) = sigma0.dim();
        let mut gamma0 = Array2::<f32>::zeros((rows, cols));

        let min_angle_rad = self.params.min_incidence_angle * PI / 180.0;
        let max_angle_rad = self.params.max_incidence_angle * PI / 180.0;

        if self.params.enable_parallel && rows > self.params.chunk_size {
            // Parallel processing - collect results then apply
            let row_indices: Vec<usize> = (0..rows).collect();
            let results: Vec<_> = row_indices.par_chunks(self.params.chunk_size)
                .map(|row_chunk| {
                    let mut local_gamma0 = Vec::new();
                    for &i in row_chunk {
                        for j in 0..cols {
                            let theta_lia = local_incidence_angles[[i, j]];
                            let sigma0_val = sigma0[[i, j]];

                            let gamma0_val = if self.params.apply_masking {
                                if theta_lia < min_angle_rad || theta_lia > max_angle_rad {
                                    f32::NAN // Mark as invalid
                                } else {
                                    // Apply terrain flattening
                                    let cos_theta = theta_lia.cos();
                                    if cos_theta > 1e-6 { // Avoid division by zero
                                        sigma0_val / cos_theta
                                    } else {
                                        f32::NAN // Mark as invalid
                                    }
                                }
                            } else {
                                // Apply terrain flattening without masking
                                let cos_theta = theta_lia.cos();
                                if cos_theta > 1e-6 { // Avoid division by zero
                                    sigma0_val / cos_theta
                                } else {
                                    f32::NAN // Mark as invalid
                                }
                            };

                            local_gamma0.push((i, j, gamma0_val));
                        }
                    }
                    local_gamma0
                })
                .collect();
            
            // Apply results
            for gamma0_data in results {
                for (i, j, val) in gamma0_data {
                    gamma0[[i, j]] = val;
                }
            }
        } else {
            // Sequential processing for small arrays
            for i in 0..rows {
                for j in 0..cols {
                    let theta_lia = local_incidence_angles[[i, j]];
                    let sigma0_val = sigma0[[i, j]];

                    // Apply masking if enabled
                    if self.params.apply_masking {
                        if theta_lia < min_angle_rad || theta_lia > max_angle_rad {
                            gamma0[[i, j]] = f32::NAN; // Mark as invalid
                            continue;
                        }
                    }

                    // Apply terrain flattening
                    let cos_theta = theta_lia.cos();
                    if cos_theta > 1e-6 { // Avoid division by zero
                        gamma0[[i, j]] = sigma0_val / cos_theta;
                    } else {
                        gamma0[[i, j]] = f32::NAN; // Mark as invalid
                    }
                }
            }
        }

        Ok(gamma0)
    }

    /// Complete terrain flattening workflow with improved geometry
    pub fn process_terrain_flattening(
        &self,
        sigma0: &Array2<f32>,
        dem: &Array2<f32>,
        _range_time: f64,
        _azimuth_time_start: f64,
        _azimuth_time_spacing: f64,
    ) -> SarResult<(Array2<f32>, Array2<f32>)> {
        log::info!("Starting terrain flattening process");
        
        let (sigma_rows, sigma_cols) = sigma0.dim();
        let (dem_rows, dem_cols) = dem.dim();
        
        log::debug!("Sigma0 dimensions: {}x{}, DEM dimensions: {}x{}", 
                   sigma_rows, sigma_cols, dem_rows, dem_cols);
        
        // Use simplified geometry for better results
        if sigma_rows * sigma_cols > 10000 {
            log::debug!("Using simplified geometry for large arrays");
            
            // Step 1: Compute slope and aspect from DEM
            let (slope, aspect) = Self::compute_slope_aspect_robust(dem, self.params.dem_pixel_spacing)?;
            
            // Step 2: Compute local incidence angles using realistic Sentinel-1 geometry
            // Use DEM dimensions for slope/aspect, not sigma0 dimensions
            let incidence_angles = Self::compute_incidence_angles_sentinel1(&slope, &aspect, dem_rows, dem_cols)?;
            
            // Step 3: Resample incidence angles to match sigma0 grid if dimensions differ
            let resampled_incidence = if (dem_rows, dem_cols) != (sigma_rows, sigma_cols) {
                log::debug!("Resampling incidence angles from {}x{} to {}x{}", 
                           dem_rows, dem_cols, sigma_rows, sigma_cols);
                Self::resample_to_sigma_grid(&incidence_angles, sigma_rows, sigma_cols)?
            } else {
                incidence_angles
            };
            
            // Step 4: Apply terrain flattening
            let gamma0 = self.apply_terrain_flattening(sigma0, &resampled_incidence)?;
            
            Ok((gamma0, resampled_incidence))
        } else {
            // Original detailed computation for small arrays
            log::debug!("Using detailed geometry computation");
            
            // Step 1: Compute slope and aspect from DEM
            let (slope, aspect) = self.compute_slope_aspect(dem)?;
            
            // Step 2: Compute surface normals
            let surface_normals = self.compute_surface_normals(&slope, &aspect);
            
            // Step 3: Compute radar look vectors (simplified)
            let look_vectors = self.compute_simplified_look_vectors(sigma_rows, sigma_cols)?;
            
            // Step 4: Compute local incidence angles
            let incidence_angles = self.compute_local_incidence_angle(&surface_normals, &look_vectors);
            
            // Step 5: Apply terrain flattening
            let gamma0 = self.apply_terrain_flattening(sigma0, &incidence_angles)?;
            
            Ok((gamma0, incidence_angles))
        }
    }

    /// Compute simplified radar look vectors for Sentinel-1
    fn compute_simplified_look_vectors(
        &self,
        rows: usize,
        cols: usize,
    ) -> SarResult<Array2<[f32; 3]>> {
        let mut look_vectors = Array2::<[f32; 3]>::from_elem((rows, cols), [0.0, 0.0, 0.0]);
        
        // Sentinel-1 look angles (approximate)
        let near_range_angle = 29.1_f32.to_radians();
        let far_range_angle = 46.0_f32.to_radians();
        
        for i in 0..rows {
            for j in 0..cols {
                // Interpolate look angle across range
                let range_fraction = j as f32 / (cols - 1) as f32;
                let look_angle = near_range_angle + range_fraction * (far_range_angle - near_range_angle);
                
                // Right-looking radar geometry
                // Look vector points from satellite to ground
                let look_x = look_angle.sin();  // Range direction (positive = away from radar)
                let look_y = 0.0;               // Azimuth direction
                let look_z = -look_angle.cos(); // Vertical component (negative = downward)
                
                look_vectors[[i, j]] = [look_x, look_y, look_z];
            }
        }
        
        Ok(look_vectors)
    }

    /// Helper function to fill edge pixels
    fn fill_edges(&self, array: &mut Array2<f32>) -> SarResult<()> {
        let (rows, cols) = array.dim();
        
        // Fill top and bottom edges
        for j in 0..cols {
            array[[0, j]] = array[[1, j]];
            array[[rows-1, j]] = array[[rows-2, j]];
        }
        
        // Fill left and right edges
        for i in 0..rows {
            array[[i, 0]] = array[[i, 1]];
            array[[i, cols-1]] = array[[i, cols-2]];
        }
        
        Ok(())
    }

    /// Interpolate satellite position at given time using orbit state vectors
    fn interpolate_satellite_position(&self, time: f64) -> SarResult<[f64; 3]> {
        let state_vectors = &self.orbit_data.state_vectors;
        
        if state_vectors.is_empty() {
            return Err(SarError::Processing("No state vectors available".to_string()));
        }
        
        // If only one state vector, return its position
        if state_vectors.len() == 1 {
            return Ok(state_vectors[0].position);
        }
        
        // Convert target time to seconds since reference
        let reference_time = self.orbit_data.reference_time.timestamp() as f64;
        let target_time_sec = time - reference_time;
        
        // Find surrounding state vectors
        let mut before_idx = 0;
        let mut after_idx = state_vectors.len() - 1;
        
        for (i, sv) in state_vectors.iter().enumerate() {
            let sv_time = sv.time.timestamp() as f64 - reference_time;
            if sv_time <= target_time_sec {
                before_idx = i;
            }
            if sv_time >= target_time_sec && after_idx == state_vectors.len() - 1 {
                after_idx = i;
                break;
            }
        }
        
        // If exact match or single vector
        if before_idx == after_idx {
            return Ok(state_vectors[before_idx].position);
        }
        
        // Linear interpolation between state vectors
        let sv_before = &state_vectors[before_idx];
        let sv_after = &state_vectors[after_idx];
        
        let time_before = sv_before.time.timestamp() as f64 - reference_time;
        let time_after = sv_after.time.timestamp() as f64 - reference_time;
        
        if (time_after - time_before).abs() < 1e-6 {
            return Ok(sv_before.position);
        }
        
        let weight = (target_time_sec - time_before) / (time_after - time_before);
        
        // Interpolate position
        let pos = [
            sv_before.position[0] * (1.0 - weight) + sv_after.position[0] * weight,
            sv_before.position[1] * (1.0 - weight) + sv_after.position[1] * weight,
            sv_before.position[2] * (1.0 - weight) + sv_after.position[2] * weight,
        ];
        
        Ok(pos)
    }

    /// Create a simple terrain flattening processor for testing
    pub fn create_simple(dem_pixel_spacing: (f64, f64)) -> SarResult<Self> {
        // Create minimal orbit data for testing
        let reference_time = chrono::Utc::now();
        let state_vector = crate::types::StateVector {
            time: reference_time,
            position: [7000000.0, 0.0, 0.0], // Typical satellite altitude
            velocity: [0.0, 7500.0, 0.0],    // Typical orbital velocity
        };
        
        let orbit_data = crate::types::OrbitData {
            state_vectors: vec![state_vector],
            reference_time,
        };
        
        let mut params = TerrainFlatteningParams::default();
        params.dem_pixel_spacing = dem_pixel_spacing;
        
        Ok(Self::new(params, orbit_data))
    }

    /// Simplified terrain flattening for production use
    /// This function handles the complete workflow with sensible defaults
    pub fn flatten_terrain_simple(
        sigma0: &Array2<f32>,
        dem: &Array2<f32>,
        _orbit_data: OrbitData,
        dem_pixel_spacing: Option<(f64, f64)>,
    ) -> SarResult<Array2<f32>> {
        log::info!("Starting simplified terrain flattening with robust geometry");
        
        let spacing = dem_pixel_spacing.unwrap_or((30.0, 30.0));
        let (rows, cols) = sigma0.dim();
        
        // Check input data validity
        if rows == 0 || cols == 0 {
            return Err(SarError::Processing("Empty input arrays".to_string()));
        }
        
        if dem.dim() != sigma0.dim() {
            return Err(SarError::Processing("DEM and sigma0 arrays must have same dimensions".to_string()));
        }
        
        // Step 1: Compute slope and aspect from DEM
        log::debug!("Computing terrain slope and aspect");
        let (slope, aspect) = Self::compute_slope_aspect_robust(dem, spacing)?;
        
        // Step 2: Compute local incidence angles using simplified geometry
        log::debug!("Computing local incidence angles");
        let incidence_angles = Self::compute_incidence_angles_sentinel1(&slope, &aspect, rows, cols)?;
        
        // Step 3: Apply terrain flattening correction
        log::debug!("Applying terrain flattening correction");
        let mut gamma0 = Array2::<f32>::zeros((rows, cols));
        
        for i in 0..rows {
            for j in 0..cols {
                let sigma0_val = sigma0[[i, j]];
                let theta_lia = incidence_angles[[i, j]];
                
                // Apply terrain flattening: gamma0 = sigma0 / cos(theta_lia)
                let cos_theta = theta_lia.cos();
                if cos_theta > 1e-3 && theta_lia < 80.0_f32.to_radians() {
                    gamma0[[i, j]] = sigma0_val / cos_theta;
                } else {
                    gamma0[[i, j]] = sigma0_val; // Keep original value for extreme angles
                }
            }
        }
        
        log::info!("Simplified terrain flattening completed successfully");
        Ok(gamma0)
    }

    /// Compute slope and aspect robustly from DEM
    fn compute_slope_aspect_robust(
        dem: &Array2<f32>,
        dem_spacing: (f64, f64),
    ) -> SarResult<(Array2<f32>, Array2<f32>)> {
        let (rows, cols) = dem.dim();
        let mut slope = Array2::<f32>::zeros((rows, cols));
        let mut aspect = Array2::<f32>::zeros((rows, cols));

        let dx_scale = dem_spacing.0 as f32;
        let dy_scale = dem_spacing.1 as f32;

        for i in 1..rows-1 {
            for j in 1..cols-1 {
                // Calculate gradients using central differences
                let dz_dx = (dem[[i, j+1]] - dem[[i, j-1]]) / (2.0 * dx_scale);
                let dz_dy = (dem[[i+1, j]] - dem[[i-1, j]]) / (2.0 * dy_scale);

                // Compute slope (in radians)
                slope[[i, j]] = (dz_dx * dz_dx + dz_dy * dz_dy).sqrt().atan();

                // Compute aspect (in radians, 0 = North, clockwise positive)
                aspect[[i, j]] = (-dz_dy).atan2(dz_dx);
            }
        }

        // Fill edges with nearest neighbor values
        for j in 0..cols {
            slope[[0, j]] = slope[[1, j]];
            slope[[rows-1, j]] = slope[[rows-2, j]];
            aspect[[0, j]] = aspect[[1, j]];
            aspect[[rows-1, j]] = aspect[[rows-2, j]];
        }
        
        for i in 0..rows {
            slope[[i, 0]] = slope[[i, 1]];
            slope[[i, cols-1]] = slope[[i, cols-2]];
            aspect[[i, 0]] = aspect[[i, 1]];
            aspect[[i, cols-1]] = aspect[[i, cols-2]];
        }

        Ok((slope, aspect))
    }

    /// Compute realistic local incidence angles for Sentinel-1 geometry
    fn compute_incidence_angles_sentinel1(
        slope: &Array2<f32>,
        aspect: &Array2<f32>,
        rows: usize,
        cols: usize,
    ) -> SarResult<Array2<f32>> {
        let mut incidence_angles = Array2::<f32>::zeros((rows, cols));
        
        // Sentinel-1 IW mode typical parameters
        let near_range_angle = 29.1_f32.to_radians();  // degrees to radians
        let far_range_angle = 46.0_f32.to_radians();   // degrees to radians
        let look_direction = -90.0_f32.to_radians();   // Right-looking SAR, perpendicular to flight direction
        
        for i in 0..rows {
            for j in 0..cols {
                // Interpolate radar incidence angle across range
                let range_fraction = j as f32 / (cols - 1) as f32;
                let radar_incidence = near_range_angle + range_fraction * (far_range_angle - near_range_angle);
                
                // Get local terrain slope and aspect
                let terrain_slope = slope[[i, j]];
                let terrain_aspect = aspect[[i, j]];
                
                // Compute the angle between radar look direction and terrain aspect
                let aspect_diff = (terrain_aspect - look_direction).abs();
                let aspect_factor = aspect_diff.cos();
                
                // Local incidence angle considers both radar geometry and terrain
                // When terrain faces the radar: smaller incidence angle
                // When terrain faces away: larger incidence angle
                let terrain_correction = terrain_slope * aspect_factor;
                let local_incidence = radar_incidence + terrain_correction;
                
                // Clamp to realistic range
                incidence_angles[[i, j]] = local_incidence.clamp(5.0_f32.to_radians(), 85.0_f32.to_radians());
            }
        }
        
        Ok(incidence_angles)
    }

    /// Enhanced terrain flattening with masking
    /// Returns both flattened data and quality mask
    pub fn flatten_terrain_with_mask(
        sigma0: &Array2<f32>,
        dem: &Array2<f32>,
        orbit_data: OrbitData,
        min_incidence: Option<f32>,
        max_incidence: Option<f32>,
    ) -> SarResult<(Array2<f32>, Array2<bool>)> {
        log::info!("Starting terrain flattening with quality masking");
        
        // Use scientifically-based defaults instead of arbitrary fallbacks
        let default_params = TerrainFlatteningParams::default();
        let min_angle = min_incidence.unwrap_or(default_params.min_incidence_angle);
        let max_angle = max_incidence.unwrap_or(default_params.max_incidence_angle);
        
        log::debug!("Using incidence angle limits: {:.1}° - {:.1}° (configurable via TerrainFlatteningParams)", 
                   min_angle, max_angle);
        
        if min_incidence.is_none() {
            log::info!("ℹ️  Using default minimum incidence angle: {:.1}° (based on SAR literature)", min_angle);
        }
        if max_incidence.is_none() {
            log::info!("ℹ️  Using default maximum incidence angle: {:.1}° (based on SAR literature)", max_angle);
        }
        
        let mut params = TerrainFlatteningParams::default();
        params.min_incidence_angle = min_angle;
        params.max_incidence_angle = max_angle;
        params.apply_masking = true;
        
        let flattener = Self::new(params, orbit_data);
        
        // Use default timing parameters
        let range_time = 0.005;
        let azimuth_time_start = 0.0;
        let azimuth_time_spacing = 1e-4;
        
        let (gamma0, incidence_angles) = flattener.process_terrain_flattening(
            sigma0,
            dem,
            range_time,
            azimuth_time_start,
            azimuth_time_spacing,
        )?;
        
        // Create quality mask
        let (rows, cols) = gamma0.dim();
        let mut quality_mask = Array2::<bool>::from_elem((rows, cols), true);
        
        let min_angle_rad = min_angle * PI / 180.0;
        let max_angle_rad = max_angle * PI / 180.0;
        
        for i in 0..rows {
            for j in 0..cols {
                let is_valid = !gamma0[[i, j]].is_nan() 
                    && incidence_angles[[i, j]] >= min_angle_rad
                    && incidence_angles[[i, j]] <= max_angle_rad;
                quality_mask[[i, j]] = is_valid;
            }
        }
        
        log::info!("Terrain flattening with masking completed");
        Ok((gamma0, quality_mask))
    }

    /// Resample incidence angles from DEM grid to sigma0 grid
    fn resample_to_sigma_grid(
        incidence_angles: &Array2<f32>,
        target_rows: usize,
        target_cols: usize,
    ) -> SarResult<Array2<f32>> {
        let (src_rows, src_cols) = incidence_angles.dim();
        let mut resampled = Array2::<f32>::zeros((target_rows, target_cols));
        
        let row_scale = (src_rows - 1) as f32 / (target_rows - 1) as f32;
        let col_scale = (src_cols - 1) as f32 / (target_cols - 1) as f32;
        
        for i in 0..target_rows {
            for j in 0..target_cols {
                let src_row = (i as f32 * row_scale).min((src_rows - 1) as f32);
                let src_col = (j as f32 * col_scale).min((src_cols - 1) as f32);
                
                let src_row_floor = src_row.floor() as usize;
                let src_col_floor = src_col.floor() as usize;
                let src_row_ceil = (src_row_floor + 1).min(src_rows - 1);
                let src_col_ceil = (src_col_floor + 1).min(src_cols - 1);
                
                // Bilinear interpolation
                let t = src_row - src_row_floor as f32;
                let s = src_col - src_col_floor as f32;
                
                let v00 = incidence_angles[[src_row_floor, src_col_floor]];
                let v01 = incidence_angles[[src_row_floor, src_col_ceil]];
                let v10 = incidence_angles[[src_row_ceil, src_col_floor]];
                let v11 = incidence_angles[[src_row_ceil, src_col_ceil]];
                
                let v0 = v00 * (1.0 - s) + v01 * s;
                let v1 = v10 * (1.0 - s) + v11 * s;
                let value = v0 * (1.0 - t) + v1 * t;
                
                resampled[[i, j]] = value;
            }
        }
        
        Ok(resampled)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use chrono::{DateTime, Utc};
    use crate::types::StateVector;

    fn create_test_orbit_data() -> OrbitData {
        let test_time = DateTime::parse_from_rfc3339("2020-01-03T17:08:15Z")
            .unwrap()
            .with_timezone(&Utc);
        
        OrbitData {
            state_vectors: vec![
                StateVector {
                    time: test_time,
                    position: [7000000.0, 0.0, 0.0],
                    velocity: [0.0, 7500.0, 0.0],
                }
            ],
            reference_time: test_time,
        }
    }

    #[test]
    fn test_slope_aspect_computation() {
        let orbit_data = create_test_orbit_data();
        let flattener = TerrainFlattener::standard(orbit_data);

        // Create a simple tilted plane DEM
        let mut dem = Array2::<f32>::zeros((5, 5));
        for i in 0..5 {
            for j in 0..5 {
                dem[[i, j]] = (i + j) as f32 * 10.0; // Rising elevation
            }
        }

        let result = flattener.compute_slope_aspect(&dem);
        assert!(result.is_ok());
        
        let (slope, aspect) = result.unwrap();
        assert_eq!(slope.dim(), dem.dim());
        assert_eq!(aspect.dim(), dem.dim());
        
        // Check that slopes are positive (tilted plane)
        assert!(slope[[2, 2]] > 0.0);
    }

    #[test]
    fn test_surface_normals() {
        let orbit_data = create_test_orbit_data();
        let flattener = TerrainFlattener::standard(orbit_data);

        let slope = Array2::<f32>::from_elem((3, 3), 0.1); // Small slope
        let aspect = Array2::<f32>::zeros((3, 3)); // Facing north

        let normals = flattener.compute_surface_normals(&slope, &aspect);
        
        assert_eq!(normals.dim(), (3, 3));
        
        // For flat terrain, normal should point mostly upward
        let normal = normals[[1, 1]];
        assert!(normal[2] > 0.9); // Z component should be close to 1
    }

    #[test]
    fn test_terrain_flattening() {
        let orbit_data = create_test_orbit_data();
        let flattener = TerrainFlattener::standard(orbit_data);

        // Create test data
        let sigma0 = Array2::<f32>::from_elem((3, 3), 0.1);
        let incidence_angles = Array2::<f32>::from_elem((3, 3), PI / 4.0); // 45 degrees

        let result = flattener.apply_terrain_flattening(&sigma0, &incidence_angles);
        assert!(result.is_ok());
        
        let gamma0 = result.unwrap();
        
        // Gamma0 should be larger than sigma0 for 45-degree incidence
        assert!(gamma0[[1, 1]] > sigma0[[1, 1]]);
    }
}
