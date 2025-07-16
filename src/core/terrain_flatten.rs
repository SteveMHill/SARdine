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

    /// Compute radar look vector for each pixel
    /// 
    /// This is a simplified version assuming constant look vector per range line
    /// For full accuracy, this should be computed per-pixel using precise orbit
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

        for i in 0..rows {
            let azimuth_time = azimuth_time_start + (i as f64) * azimuth_time_spacing;
            
            // Interpolate satellite position at this azimuth time
            let sat_pos = self.interpolate_satellite_position(azimuth_time)?;

            for j in 0..cols {
                // For simplicity, assume constant look vector per azimuth line
                // In reality, this should account for earth curvature and precise geometry
                
                // Simplified look vector calculation
                // This points from satellite to ground (nadir direction)
                let look_vector = [0.0, 0.0, -1.0]; // Simplified: pointing down
                
                look_vectors[[i, j]] = look_vector;
            }
        }

        Ok(look_vectors)
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

    /// Interpolate satellite position at given time
    fn interpolate_satellite_position(&self, time: f64) -> SarResult<[f64; 3]> {
        let state_vectors = &self.orbit_data.state_vectors;
        
        // For simplicity, use the first state vector
        // In reality, this should interpolate between state vectors
        if let Some(first_sv) = state_vectors.first() {
            Ok(first_sv.position)
        } else {
            Err(SarError::Processing("No state vectors available".to_string()))
        }
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
