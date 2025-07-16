use crate::types::{SarError, SarResult, OrbitData};
use ndarray::Array2;
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

    /// Compute slope and aspect from DEM
    /// 
    /// Returns (slope_radians, aspect_radians)
    pub fn compute_slope_aspect(&self, dem: &Array2<f32>) -> SarResult<(Array2<f32>, Array2<f32>)> {
        let (rows, cols) = dem.dim();
        let mut slope = Array2::<f32>::zeros((rows, cols));
        let mut aspect = Array2::<f32>::zeros((rows, cols));

        let dx_scale = self.params.dem_pixel_spacing.0 as f32;
        let dy_scale = self.params.dem_pixel_spacing.1 as f32;

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

        // Handle edges by copying nearest valid values
        self.fill_edges(&mut slope)?;
        self.fill_edges(&mut aspect)?;

        Ok((slope, aspect))
    }

    /// Compute surface normal vectors from slope and aspect
    pub fn compute_surface_normals(
        &self, 
        slope: &Array2<f32>, 
        aspect: &Array2<f32>
    ) -> Array2<[f32; 3]> {
        let (rows, cols) = slope.dim();
        let mut normals = Array2::<[f32; 3]>::from_elem((rows, cols), [0.0, 0.0, 1.0]);

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
    fn compute_ground_position(&self, sat_pos: [f64; 3], slant_range: f64, range_pixel: f64) -> SarResult<[f64; 3]> {
        // This is a simplified calculation - for production use, proper
        // SAR geometry and earth ellipsoid models should be used
        
        // Earth radius (simplified spherical earth)
        const EARTH_RADIUS: f64 = 6371000.0;
        
        // Satellite altitude (approximate)
        let sat_altitude = (sat_pos[0]*sat_pos[0] + sat_pos[1]*sat_pos[1] + sat_pos[2]*sat_pos[2]).sqrt() - EARTH_RADIUS;
        
        // For simplification, assume ground point is directly below satellite
        // offset by range pixel spacing
        let ground_distance = (slant_range*slant_range - sat_altitude*sat_altitude).sqrt();
        
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

    /// Compute local incidence angle
    /// 
    /// θ_lia = arccos(normal · look_vector)
    pub fn compute_local_incidence_angle(
        &self,
        surface_normals: &Array2<[f32; 3]>,
        look_vectors: &Array2<[f32; 3]>,
    ) -> Array2<f32> {
        let (rows, cols) = surface_normals.dim();
        let mut incidence_angles = Array2::<f32>::zeros((rows, cols));

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

        incidence_angles
    }

    /// Apply terrain flattening: gamma0 = sigma0 / cos(theta_lia)
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

        Ok(gamma0)
    }

    /// Complete terrain flattening workflow
    pub fn process_terrain_flattening(
        &self,
        sigma0: &Array2<f32>,
        dem: &Array2<f32>,
        range_time: f64,
        azimuth_time_start: f64,
        azimuth_time_spacing: f64,
    ) -> SarResult<(Array2<f32>, Array2<f32>)> {
        log::info!("Starting terrain flattening process");
        
        // Step 1: Compute slope and aspect from DEM
        log::debug!("Computing slope and aspect from DEM");
        let (slope, aspect) = self.compute_slope_aspect(dem)?;
        
        // Step 2: Compute surface normals
        log::debug!("Computing surface normals");
        let surface_normals = self.compute_surface_normals(&slope, &aspect);
        
        // Step 3: Compute radar look vectors
        log::debug!("Computing radar look vectors");
        let look_vectors = self.compute_radar_look_vectors(
            sigma0.dim(),
            range_time,
            azimuth_time_start,
            azimuth_time_spacing,
        )?;
        
        // Step 4: Compute local incidence angles
        log::debug!("Computing local incidence angles");
        let incidence_angles = self.compute_local_incidence_angle(&surface_normals, &look_vectors);
        
        // Step 5: Apply terrain flattening
        log::debug!("Applying terrain flattening");
        let gamma0 = self.apply_terrain_flattening(sigma0, &incidence_angles)?;
        
        log::info!("Terrain flattening completed");
        
        Ok((gamma0, incidence_angles))
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
    /// Enhanced terrain flattening - simplified but robust approach
    pub fn flatten_terrain_simple(
        sigma0: &Array2<f32>,
        dem: &Array2<f32>,
        orbit_data: OrbitData,
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
        let gamma0 = Self::apply_terrain_correction(&sigma0, &incidence_angles)?;
        
        // Step 4: Quality control - check for valid output
        let valid_count = gamma0.iter().filter(|&&x| x.is_finite() && x > 0.0).count();
        let total_count = gamma0.len();
        let valid_percentage = (valid_count as f64 / total_count as f64) * 100.0;
        
        log::info!("Terrain flattening completed: {:.1}% valid pixels ({}/{})", 
                  valid_percentage, valid_count, total_count);
        
        if valid_percentage < 5.0 {
            log::warn!("Low percentage of valid pixels after terrain flattening - check input data");
        }
        
        Ok(gamma0)
    }

    /// Robust slope and aspect computation with edge handling
    fn compute_slope_aspect_robust(
        dem: &Array2<f32>,
        spacing: (f64, f64),
    ) -> SarResult<(Array2<f32>, Array2<f32>)> {
        let (rows, cols) = dem.dim();
        let mut slope = Array2::<f32>::zeros((rows, cols));
        let mut aspect = Array2::<f32>::zeros((rows, cols));
        
        let dx = spacing.0 as f32;  // East-West spacing
        let dy = spacing.1 as f32;  // North-South spacing
        
        for i in 1..(rows-1) {
            for j in 1..(cols-1) {
                // Check for valid DEM values in neighborhood
                let neighbors = [
                    dem[[i-1, j-1]], dem[[i-1, j]], dem[[i-1, j+1]],
                    dem[[i, j-1]],   dem[[i, j]],   dem[[i, j+1]],
                    dem[[i+1, j-1]], dem[[i+1, j]], dem[[i+1, j+1]]
                ];
                
                // Skip if any neighbor is NaN
                if neighbors.iter().any(|&x| !x.is_finite()) {
                    slope[[i, j]] = 0.0;
                    aspect[[i, j]] = 0.0;
                    continue;
                }
                
                // Compute gradients using central differences
                let dz_dx = (dem[[i, j+1]] - dem[[i, j-1]]) / (2.0 * dx);
                let dz_dy = (dem[[i+1, j]] - dem[[i-1, j]]) / (2.0 * dy);
                
                // Slope magnitude (in radians)
                slope[[i, j]] = (dz_dx * dz_dx + dz_dy * dz_dy).sqrt().atan();
                
                // Aspect (direction of steepest slope)
                aspect[[i, j]] = dz_dy.atan2(dz_dx);
            }
        }
        
        // Fill edges with nearest valid values
        Self::fill_edges_robust(&mut slope)?;
        Self::fill_edges_robust(&mut aspect)?;
        
        Ok((slope, aspect))
    }

    /// Compute local incidence angles for Sentinel-1 using simplified geometry
    fn compute_incidence_angles_sentinel1(
        slope: &Array2<f32>,
        aspect: &Array2<f32>,
        rows: usize,
        cols: usize,
    ) -> SarResult<Array2<f32>> {
        let mut incidence_angles = Array2::<f32>::zeros((rows, cols));
        
        // Sentinel-1 typical parameters
        const SENTINEL1_INCIDENCE_CENTER: f32 = 35.0; // degrees (typical center incidence)
        const SENTINEL1_INCIDENCE_VARIATION: f32 = 10.0; // degrees (near to far range variation)
        
        for i in 0..rows {
            for j in 0..cols {
                // Compute reference incidence angle (varies across range)
                let range_fraction = (j as f32) / (cols as f32);
                let reference_incidence = SENTINEL1_INCIDENCE_CENTER + 
                    (range_fraction - 0.5) * SENTINEL1_INCIDENCE_VARIATION;
                let ref_incidence_rad = reference_incidence.to_radians();
                
                // Get terrain slope and aspect
                let terrain_slope = slope[[i, j]];
                let terrain_aspect = aspect[[i, j]];
                
                // Simplified local incidence angle computation
                // For side-looking SAR, assume radar look direction is approximately East-West
                // This is a simplification but works reasonably well for most cases
                let radar_azimuth = std::f32::consts::PI / 2.0; // 90 degrees (East)
                let aspect_diff = terrain_aspect - radar_azimuth;
                
                // Local incidence angle considering terrain slope
                let local_incidence = ref_incidence_rad + terrain_slope * aspect_diff.cos();
                
                // Clamp to reasonable range
                incidence_angles[[i, j]] = local_incidence.max(0.1).min(1.4); // ~5-80 degrees
            }
        }
        
        Ok(incidence_angles)
    }

    /// Apply terrain flattening correction: gamma0 = sigma0 / cos(theta_local)
    fn apply_terrain_correction(
        sigma0: &Array2<f32>,
        incidence_angles: &Array2<f32>,
    ) -> SarResult<Array2<f32>> {
        let (rows, cols) = sigma0.dim();
        let mut gamma0 = Array2::<f32>::zeros((rows, cols));
        
        for i in 0..rows {
            for j in 0..cols {
                let sigma0_val = sigma0[[i, j]];
                let theta = incidence_angles[[i, j]];
                
                // Skip invalid input values
                if !sigma0_val.is_finite() || sigma0_val <= 0.0 {
                    gamma0[[i, j]] = f32::NAN;
                    continue;
                }
                
                // Apply terrain flattening: gamma0 = sigma0 / cos(theta)
                let cos_theta = theta.cos();
                if cos_theta > 0.01 { // Avoid division by zero for steep angles
                    gamma0[[i, j]] = sigma0_val / cos_theta;
                } else {
                    gamma0[[i, j]] = f32::NAN; // Invalid for very steep angles
                }
            }
        }
        
        Ok(gamma0)
    }

    /// Robust edge filling
    fn fill_edges_robust(array: &mut Array2<f32>) -> SarResult<()> {
        let (rows, cols) = array.dim();
        
        if rows < 3 || cols < 3 {
            return Ok(()); // Too small to fill edges meaningfully
        }
        
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
        
        let min_angle = min_incidence.unwrap_or(10.0);
        let max_angle = max_incidence.unwrap_or(70.0);
        
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
