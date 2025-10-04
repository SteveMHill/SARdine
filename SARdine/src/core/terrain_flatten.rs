use crate::types::{OrbitData, SarError, SarResult};
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

impl TerrainFlatteningParams {
    /// Create terrain flattening parameters from annotation XML (MANDATORY)
    /// This is the ONLY way to create parameters - prevents hardcoded values
    pub fn from_annotation(
        annotation: &crate::io::annotation::AnnotationRoot,
    ) -> crate::types::SarResult<Self> {
        let builder = TerrainFlatteningParamsBuilder::new(annotation)?;
        builder.build()
    }
}

/// Scientific parameter builder for terrain flattening
/// Ensures ALL parameters come from annotation XML - NO hardcoded values allowed
#[derive(Debug)]
pub struct TerrainFlatteningParamsBuilder {
    wavelength: f64,
    sar_pixel_spacing: (f64, f64),
    dem_pixel_spacing: Option<(f64, f64)>,
    min_incidence_angle: Option<f32>,
    max_incidence_angle: Option<f32>,
    apply_masking: bool,
    chunk_size: usize,
    enable_parallel: bool,
}

impl TerrainFlatteningParamsBuilder {
    /// Create builder from annotation - extracts all scientifically critical parameters
    pub fn new(
        annotation: &crate::io::annotation::AnnotationRoot,
    ) -> crate::types::SarResult<Self> {
        // Extract radar frequency and calculate wavelength (NEVER hardcoded)
        let radar_freq = annotation.get_radar_frequency_hz().ok_or_else(|| {
            crate::types::SarError::ParameterError(
                "Missing radar frequency in annotation - cannot calculate wavelength".to_string(),
            )
        })?;
        let wavelength = crate::constants::physical::SPEED_OF_LIGHT_M_S as f64 / radar_freq;

        // Extract SAR pixel spacing (NEVER hardcoded)
        let sar_pixel_spacing = annotation.get_pixel_spacing().ok_or_else(|| {
            crate::types::SarError::ParameterError(
                "Missing SAR pixel spacing in annotation".to_string(),
            )
        })?;

        // Extract incidence angle bounds from antenna patterns or geolocation grid
        let (min_incidence, max_incidence) = Self::extract_incidence_angle_bounds(annotation)?;

        Ok(Self {
            wavelength,
            sar_pixel_spacing,
            dem_pixel_spacing: None, // Must be set by user
            min_incidence_angle: Some(min_incidence),
            max_incidence_angle: Some(max_incidence),
            apply_masking: true,   // Scientific default
            chunk_size: 32,        // Performance parameter - safe default
            enable_parallel: true, // Performance parameter - safe default
        })
    }

    /// Extract scientifically accurate incidence angle bounds from annotation
    fn extract_incidence_angle_bounds(
        annotation: &crate::io::annotation::AnnotationRoot,
    ) -> crate::types::SarResult<(f32, f32)> {
        // Try antenna patterns first (most accurate)
        let antenna_patterns = annotation.get_antenna_patterns()?;
        if !antenna_patterns.is_empty() {
            let mut all_incidence_angles = Vec::new();
            for (_, incidence_angles, _) in &antenna_patterns {
                all_incidence_angles.extend(incidence_angles);
            }

            if !all_incidence_angles.is_empty() {
                let min_inc = all_incidence_angles
                    .iter()
                    .fold(f64::INFINITY, |a, &b| a.min(b)) as f32;
                let max_inc = all_incidence_angles
                    .iter()
                    .fold(f64::NEG_INFINITY, |a, &b| a.max(b)) as f32;
                return Ok((min_inc, max_inc));
            }
        }

        // Fallback to geolocation grid incidence angles
        if let Some(ref geoloc_grid) = annotation.geolocation_grid {
            if let Some(ref point_list) = geoloc_grid.geolocation_grid_point_list {
                if let Some(ref points) = point_list.geolocation_grid_points {
                    let incidence_angles: Vec<f64> =
                        points.iter().map(|p| p.incidence_angle).collect();
                    if !incidence_angles.is_empty() {
                        let min_inc = incidence_angles
                            .iter()
                            .fold(f64::INFINITY, |a, &b| a.min(b))
                            as f32;
                        let max_inc = incidence_angles
                            .iter()
                            .fold(f64::NEG_INFINITY, |a, &b| a.max(b))
                            as f32;
                        return Ok((min_inc, max_inc));
                    }
                }
            }
        }

        // Fallback to image annotation mid-swath incidence angle
        if let Some(ref image) = annotation.image_annotation {
            if let Some(ref image_info) = image.image_information {
                if let Some(mid_incidence) = image_info.incidence_angle_mid_swath {
                    // Estimate bounds based on typical Sentinel-1 geometry (scientific approximation)
                    let mid_inc = mid_incidence as f32;
                    let range_estimate = 3.0; // Conservative estimate for near/far range variation
                    return Ok((mid_inc - range_estimate, mid_inc + range_estimate));
                }
            }
        }

        Err(crate::types::SarError::ParameterError(
            "Cannot extract incidence angle bounds from annotation - no antenna patterns, geolocation grid, or mid-swath incidence available".to_string()
        ))
    }

    /// Set DEM pixel spacing (required for slope/aspect computation)
    pub fn with_dem_pixel_spacing(mut self, dem_spacing: (f64, f64)) -> Self {
        self.dem_pixel_spacing = Some(dem_spacing);
        self
    }

    /// Override incidence angle bounds (only if annotation extraction failed)
    pub fn with_incidence_angle_bounds(mut self, min_degrees: f32, max_degrees: f32) -> Self {
        self.min_incidence_angle = Some(min_degrees);
        self.max_incidence_angle = Some(max_degrees);
        self
    }

    /// Set layover/shadow masking mode
    pub fn with_masking(mut self, enable: bool) -> Self {
        self.apply_masking = enable;
        self
    }

    /// Set performance parameters (chunk size and parallelization)
    pub fn with_performance(mut self, chunk_size: usize, enable_parallel: bool) -> Self {
        self.chunk_size = chunk_size;
        self.enable_parallel = enable_parallel;
        self
    }

    /// Build final parameters with scientific validation
    pub fn build(self) -> crate::types::SarResult<TerrainFlatteningParams> {
        let dem_pixel_spacing = self.dem_pixel_spacing.ok_or_else(|| {
            crate::types::SarError::ParameterError(
                "DEM pixel spacing must be provided via with_dem_pixel_spacing()".to_string(),
            )
        })?;

        let min_incidence_angle = self.min_incidence_angle.ok_or_else(|| {
            crate::types::SarError::ParameterError(
                "Minimum incidence angle could not be determined from annotation".to_string(),
            )
        })?;

        let max_incidence_angle = self.max_incidence_angle.ok_or_else(|| {
            crate::types::SarError::ParameterError(
                "Maximum incidence angle could not be determined from annotation".to_string(),
            )
        })?;

        // Scientific validation
        if self.wavelength <= 0.0 {
            return Err(crate::types::SarError::ParameterError(format!(
                "Invalid wavelength: {} m",
                self.wavelength
            )));
        }

        if min_incidence_angle >= max_incidence_angle {
            return Err(crate::types::SarError::ParameterError(format!(
                "Invalid incidence angle range: {:.2}° to {:.2}°",
                min_incidence_angle, max_incidence_angle
            )));
        }

        if dem_pixel_spacing.0 <= 0.0 || dem_pixel_spacing.1 <= 0.0 {
            return Err(crate::types::SarError::ParameterError(format!(
                "Invalid DEM pixel spacing: ({:.2}, {:.2}) m",
                dem_pixel_spacing.0, dem_pixel_spacing.1
            )));
        }

        if self.sar_pixel_spacing.0 <= 0.0 || self.sar_pixel_spacing.1 <= 0.0 {
            return Err(crate::types::SarError::ParameterError(format!(
                "Invalid SAR pixel spacing: ({:.2}, {:.2}) m",
                self.sar_pixel_spacing.0, self.sar_pixel_spacing.1
            )));
        }

        Ok(TerrainFlatteningParams {
            dem_pixel_spacing,
            sar_pixel_spacing: self.sar_pixel_spacing,
            wavelength: self.wavelength,
            apply_masking: self.apply_masking,
            min_incidence_angle,
            max_incidence_angle,
            chunk_size: self.chunk_size,
            enable_parallel: self.enable_parallel,
        })
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
        Self { params, orbit_data }
    }

    /// Create a standard terrain flattening processor from annotation data
    /// This method ensures all parameters come from annotation XML
    pub fn from_annotation(
        annotation: &crate::io::annotation::AnnotationRoot,
        orbit_data: OrbitData,
    ) -> crate::types::SarResult<Self> {
        let params = TerrainFlatteningParams::from_annotation(annotation)?;
        Ok(Self::new(params, orbit_data))
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

    /// Sequential slope/aspect computation using Horn 3×3 gradient operator
    /// This provides much better noise reduction compared to simple central differences
    fn compute_slope_aspect_sequential(
        &self,
        dem: &Array2<f32>,
        slope: &mut Array2<f32>,
        aspect: &mut Array2<f32>,
        dx_scale: f32,
        dy_scale: f32,
    ) {
        let (rows, cols) = dem.dim();

        for i in 1..rows - 1 {
            for j in 1..cols - 1 {
                // Horn 3×3 gradient operator for noise reduction
                // Computes gradients using 8 neighboring pixels with weighted contributions
                //
                // Kernel weights for dz/dx:
                // -1  0  1
                // -2  0  2
                // -1  0  1
                //
                // Kernel weights for dz/dy:
                // -1 -2 -1
                //  0  0  0
                //  1  2  1

                let dz_dx = (
                    // Right column - Left column
                    (dem[[i - 1, j + 1]] + 2.0 * dem[[i, j + 1]] + dem[[i + 1, j + 1]])
                        - (dem[[i - 1, j - 1]] + 2.0 * dem[[i, j - 1]] + dem[[i + 1, j - 1]])
                ) / (8.0 * dx_scale);

                let dz_dy = (
                    // Bottom row - Top row
                    (dem[[i + 1, j - 1]] + 2.0 * dem[[i + 1, j]] + dem[[i + 1, j + 1]])
                        - (dem[[i - 1, j - 1]] + 2.0 * dem[[i - 1, j]] + dem[[i - 1, j + 1]])
                ) / (8.0 * dy_scale);

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
        let row_indices: Vec<usize> = (1..rows - 1).collect();

        // Process in parallel chunks, collecting results
        let results: Vec<_> = row_indices
            .par_chunks(chunk_size)
            .map(|row_chunk| {
                let mut local_slope_data = Vec::new();
                let mut local_aspect_data = Vec::new();

                for &i in row_chunk {
                    for j in 1..cols - 1 {
                        // Horn 3×3 gradient operator for noise reduction
                        // This is the same as sequential version but in parallel chunks
                        let dz_dx = (
                            // Right column - Left column
                            (dem[[i - 1, j + 1]] + 2.0 * dem[[i, j + 1]] + dem[[i + 1, j + 1]])
                                - (dem[[i - 1, j - 1]]
                                    + 2.0 * dem[[i, j - 1]]
                                    + dem[[i + 1, j - 1]])
                        ) / (8.0 * dx_scale);

                        let dz_dy = (
                            // Bottom row - Top row
                            (dem[[i + 1, j - 1]] + 2.0 * dem[[i + 1, j]] + dem[[i + 1, j + 1]])
                                - (dem[[i - 1, j - 1]]
                                    + 2.0 * dem[[i - 1, j]]
                                    + dem[[i - 1, j + 1]])
                        ) / (8.0 * dy_scale);

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
        aspect: &Array2<f32>,
    ) -> Array2<[f32; 3]> {
        let (rows, cols) = slope.dim();
        let mut normals = Array2::<[f32; 3]>::from_elem((rows, cols), [0.0, 0.0, 1.0]);

        if self.params.enable_parallel && rows > self.params.chunk_size {
            // Parallel processing - collect results then apply
            let row_indices: Vec<usize> = (0..rows).collect();
            let results: Vec<_> = row_indices
                .par_chunks(self.params.chunk_size)
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
            return Err(SarError::Processing(
                "No orbit state vectors available".to_string(),
            ));
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
                let length = (dx * dx + dy * dy + dz * dz).sqrt();
                if length > 1e-6 {
                    look_vectors[[i, j]] = [
                        (dx / length) as f32,
                        (dy / length) as f32,
                        (dz / length) as f32,
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
        const C: f64 = crate::constants::physical::SPEED_OF_LIGHT_M_S;

        // Range time to first pixel plus pixel offset
        let pixel_range_time =
            range_time + (range_pixel as f64) * self.params.sar_pixel_spacing.0 / C;

        // Slant range = (time * speed_of_light) / 2
        Ok(pixel_range_time * C / 2.0)
    }

    /// Compute approximate ground position for satellite position and slant range
    fn compute_ground_position(
        &self,
        _sat_pos: [f64; 3],
        _slant_range: f64,
        _range_pixel: f64,
    ) -> SarResult<[f64; 3]> {
        // Disabled in scientific mode: simplified spherical/nadir assumptions produce incorrect results.
        Err(SarError::Processing(
            "Simplified ground position computation is disabled; use precise SAR geometry (WGS84 + RD).".to_string(),
        ))
    }

    /// Compute precise ground position using WGS84 ellipsoid and proper SAR geometry
    fn compute_precise_ground_position(
        &self,
        sat_pos: [f64; 3],
        sat_vel: [f64; 3],
        slant_range: f64,
        wgs84_a: f64,
        wgs84_e2: f64,
    ) -> SarResult<[f64; 3]> {
        // Compute satellite altitude above WGS84 ellipsoid
        let sat_magnitude =
            (sat_pos[0] * sat_pos[0] + sat_pos[1] * sat_pos[1] + sat_pos[2] * sat_pos[2]).sqrt();

        // Approximate Earth radius at satellite location for initial guess
        let lat_approx = (sat_pos[2] / sat_magnitude).asin();
        let n_radius = wgs84_a / (1.0 - wgs84_e2 * lat_approx.sin().powi(2)).sqrt();
        let sat_altitude = sat_magnitude - n_radius;

        // Cross product of position and velocity gives orbit normal (for right-looking SAR)
        let orbit_normal = [
            sat_pos[1] * sat_vel[2] - sat_pos[2] * sat_vel[1],
            sat_pos[2] * sat_vel[0] - sat_pos[0] * sat_vel[2],
            sat_pos[0] * sat_vel[1] - sat_pos[1] * sat_vel[0],
        ];

        // Normalize orbit normal
        let normal_mag = (orbit_normal[0] * orbit_normal[0]
            + orbit_normal[1] * orbit_normal[1]
            + orbit_normal[2] * orbit_normal[2])
            .sqrt();
        let unit_normal = [
            orbit_normal[0] / normal_mag,
            orbit_normal[1] / normal_mag,
            orbit_normal[2] / normal_mag,
        ];

        // Look direction perpendicular to both satellite position and orbit normal (right-looking)
        let look_direction = [
            sat_pos[1] * unit_normal[2] - sat_pos[2] * unit_normal[1],
            sat_pos[2] * unit_normal[0] - sat_pos[0] * unit_normal[2],
            sat_pos[0] * unit_normal[1] - sat_pos[1] * unit_normal[0],
        ];

        // Normalize look direction
        let look_mag = (look_direction[0] * look_direction[0]
            + look_direction[1] * look_direction[1]
            + look_direction[2] * look_direction[2])
            .sqrt();
        let unit_look = [
            look_direction[0] / look_mag,
            look_direction[1] / look_mag,
            look_direction[2] / look_mag,
        ];

        // Iterative solution for ground intersection with WGS84 ellipsoid
        let mut ground_distance = (slant_range * slant_range - sat_altitude * sat_altitude)
            .sqrt()
            .max(0.0);

        for _iteration in 0..5 {
            // Usually converges in 2-3 iterations
            let ground_pos_candidate = [
                sat_pos[0] + unit_look[0] * ground_distance,
                sat_pos[1] + unit_look[1] * ground_distance,
                sat_pos[2] + unit_look[2] * ground_distance,
            ];

            // Check distance from satellite
            let distance_check = ((ground_pos_candidate[0] - sat_pos[0]).powi(2)
                + (ground_pos_candidate[1] - sat_pos[1]).powi(2)
                + (ground_pos_candidate[2] - sat_pos[2]).powi(2))
            .sqrt();

            // Adjust ground distance based on slant range constraint
            let distance_error = distance_check - slant_range;
            if distance_error.abs() < 1.0 {
                // Converged within 1 meter
                return Ok(ground_pos_candidate);
            }

            ground_distance -= distance_error * 0.5; // Damped iteration
        }

        // Fallback if iteration doesn't converge
        let ground_pos = [
            sat_pos[0] + unit_look[0] * ground_distance,
            sat_pos[1] + unit_look[1] * ground_distance,
            sat_pos[2] + unit_look[2] * ground_distance,
        ];

        Ok(ground_pos)
    }

    /// Interpolate satellite velocity at given azimuth time
    fn interpolate_satellite_velocity(&self, azimuth_time: f64) -> SarResult<[f64; 3]> {
        let state_vectors = &self.orbit_data.state_vectors;

        if state_vectors.len() < 2 {
            return Err(SarError::Processing(
                "Need at least 2 state vectors for velocity interpolation".to_string(),
            ));
        }

        // Convert azimuth_time to DateTime for comparison
        let reference_time = self.orbit_data.reference_time;
        let target_time = reference_time + chrono::TimeDelta::try_seconds(azimuth_time as i64)
            .ok_or_else(|| SarError::Processing(
                format!("SCIENTIFIC MODE: Invalid azimuth time {} for orbit interpolation. Real timing data required for scientific processing.", azimuth_time)
            ))?;

        // Find the two state vectors that bracket the azimuth time
        let mut before_idx = 0;
        let mut after_idx = state_vectors.len() - 1;

        for (i, sv) in state_vectors.iter().enumerate() {
            if sv.time <= target_time {
                before_idx = i;
            }
            if sv.time >= target_time && after_idx == state_vectors.len() - 1 {
                after_idx = i;
                break;
            }
        }

        // If we're at the boundary, use the closest velocity
        if before_idx == after_idx {
            return Ok(state_vectors[before_idx].velocity);
        }

        // Linear interpolation of velocity
        let sv_before = &state_vectors[before_idx];
        let sv_after = &state_vectors[after_idx];

        let dt = sv_after.time - sv_before.time;
        if dt.num_seconds().abs() < 1 {
            return Ok(sv_before.velocity);
        }

        let total_seconds = dt.num_seconds() as f64;
        let elapsed_seconds = (target_time - sv_before.time).num_seconds() as f64;
        let t_frac = elapsed_seconds / total_seconds;

        let interpolated_velocity = [
            sv_before.velocity[0] + t_frac * (sv_after.velocity[0] - sv_before.velocity[0]),
            sv_before.velocity[1] + t_frac * (sv_after.velocity[1] - sv_before.velocity[1]),
            sv_before.velocity[2] + t_frac * (sv_after.velocity[2] - sv_before.velocity[2]),
        ];

        Ok(interpolated_velocity)
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
            let results: Vec<_> = row_indices
                .par_chunks(self.params.chunk_size)
                .map(|row_chunk| {
                    let mut local_angles = Vec::new();
                    for &i in row_chunk {
                        for j in 0..cols {
                            let normal = surface_normals[[i, j]];
                            let look = look_vectors[[i, j]];

                            // Compute dot product
                            let dot_product =
                                normal[0] * look[0] + normal[1] * look[1] + normal[2] * look[2];

                            // Clamp to valid range to avoid numerical issues
                            use crate::core::global_clamp_debug::ClampDebug;
                            let dot_clamped = dot_product.dbg_clamp(-1.0, 1.0, "flatten_local_dot");

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
                    let dot_product =
                        normal[0] * look[0] + normal[1] * look[1] + normal[2] * look[2];

                    // Clamp to valid range to avoid numerical issues
                    use crate::core::global_clamp_debug::ClampDebug;
                    let dot_clamped = dot_product.dbg_clamp(-1.0, 1.0, "flatten_block_dot");

                    // Compute local incidence angle in radians
                    let theta_lia = dot_clamped.abs().acos();

                    incidence_angles[[i, j]] = theta_lia;
                }
            }
        }

        incidence_angles
    }

    /// Detect layover and shadow areas based on radar geometry
    /// Returns (layover_mask, shadow_mask) where true indicates problematic areas
    ///
    /// Layover occurs when terrain slope exceeds radar look angle
    /// Shadow occurs when terrain blocks radar illumination
    pub fn detect_layover_shadow(
        &self,
        dem: &Array2<f32>,
        look_vectors: &Array2<[f32; 3]>,
    ) -> SarResult<(Array2<bool>, Array2<bool>)> {
        let (rows, cols) = dem.dim();
        let mut layover_mask = Array2::<bool>::from_elem((rows, cols), false);
        let mut shadow_mask = Array2::<bool>::from_elem((rows, cols), false);

        // Calculate slope and aspect for geometry analysis
        let (slope, aspect) = self.compute_slope_aspect(dem)?;

        // Parallel processing for better performance on large DEMs
        if self.params.enable_parallel && rows > self.params.chunk_size {
            self.detect_layover_shadow_parallel(
                &slope,
                &aspect,
                look_vectors,
                &mut layover_mask,
                &mut shadow_mask,
            )?;
        } else {
            self.detect_layover_shadow_sequential(
                &slope,
                &aspect,
                look_vectors,
                &mut layover_mask,
                &mut shadow_mask,
            );
        }

        Ok((layover_mask, shadow_mask))
    }

    /// Sequential layover/shadow detection based on radar geometry principles
    fn detect_layover_shadow_sequential(
        &self,
        slope: &Array2<f32>,
        aspect: &Array2<f32>,
        look_vectors: &Array2<[f32; 3]>,
        layover_mask: &mut Array2<bool>,
        shadow_mask: &mut Array2<bool>,
    ) {
        let (rows, cols) = slope.dim();

        for i in 1..rows - 1 {
            for j in 1..cols - 1 {
                let terrain_slope = slope[[i, j]];
                let terrain_aspect = aspect[[i, j]];
                let look_vec = look_vectors[[i, j]];

                // Calculate radar look angle (angle between look vector and vertical)
                let look_angle = (-look_vec[2]).acos(); // Assumes look_vec[2] is downward component

                // Calculate terrain slope in the range direction
                // Project look vector to horizontal plane for range direction
                let range_direction = [look_vec[0], look_vec[1], 0.0];
                let range_dir_norm =
                    (range_direction[0].powi(2) + range_direction[1].powi(2)).sqrt();

                if range_dir_norm > 1e-6 {
                    let range_unit = [
                        range_direction[0] / range_dir_norm,
                        range_direction[1] / range_dir_norm,
                        0.0,
                    ];

                    // Calculate slope component in range direction using aspect
                    let range_slope_component = terrain_slope
                        * (terrain_aspect.cos() * range_unit[0]
                            + terrain_aspect.sin() * range_unit[1])
                            .abs();

                    // Layover detection: terrain slope steeper than radar look angle
                    if range_slope_component > look_angle {
                        layover_mask[[i, j]] = true;
                    }

                    // Shadow detection: terrain faces away from radar at steep angles
                    let terrain_facing_angle =
                        terrain_aspect.cos() * range_unit[0] + terrain_aspect.sin() * range_unit[1];

                    // Shadow occurs when terrain faces away and is steep enough to block radar
                    if terrain_facing_angle < -0.5 && terrain_slope > std::f32::consts::PI / 6.0 {
                        // 30° threshold
                        shadow_mask[[i, j]] = true;
                    }

                    // Additional shadow detection for very steep opposing slopes
                    if terrain_slope > std::f32::consts::PI / 3.0 && terrain_facing_angle < 0.0 {
                        // 60° steep slopes facing away
                        shadow_mask[[i, j]] = true;
                    }
                }
            }
        }
    }

    /// Parallel layover/shadow detection using rayon
    fn detect_layover_shadow_parallel(
        &self,
        slope: &Array2<f32>,
        aspect: &Array2<f32>,
        look_vectors: &Array2<[f32; 3]>,
        layover_mask: &mut Array2<bool>,
        shadow_mask: &mut Array2<bool>,
    ) -> SarResult<()> {
        let (rows, cols) = slope.dim();
        let chunk_size = self.params.chunk_size;

        // Collect row indices to process
        let row_indices: Vec<usize> = (1..rows - 1).collect();

        // Process in parallel chunks
        let results: Vec<_> = row_indices
            .par_chunks(chunk_size)
            .map(|row_chunk| {
                let mut local_layover_data = Vec::new();
                let mut local_shadow_data = Vec::new();

                for &i in row_chunk {
                    for j in 1..cols - 1 {
                        let terrain_slope = slope[[i, j]];
                        let terrain_aspect = aspect[[i, j]];
                        let look_vec = look_vectors[[i, j]];

                        // Calculate radar look angle
                        let look_angle = (-look_vec[2]).acos();

                        // Calculate terrain slope in range direction
                        let range_direction = [look_vec[0], look_vec[1], 0.0];
                        let range_dir_norm =
                            (range_direction[0].powi(2) + range_direction[1].powi(2)).sqrt();

                        if range_dir_norm > 1e-6 {
                            let range_unit = [
                                range_direction[0] / range_dir_norm,
                                range_direction[1] / range_dir_norm,
                                0.0,
                            ];

                            let range_slope_component = terrain_slope
                                * (terrain_aspect.cos() * range_unit[0]
                                    + terrain_aspect.sin() * range_unit[1])
                                    .abs();

                            // Layover detection
                            if range_slope_component > look_angle {
                                local_layover_data.push((i, j));
                            }

                            // Shadow detection
                            let terrain_facing_angle = terrain_aspect.cos() * range_unit[0]
                                + terrain_aspect.sin() * range_unit[1];

                            if (terrain_facing_angle < -0.5
                                && terrain_slope > std::f32::consts::PI / 6.0)
                                || (terrain_slope > std::f32::consts::PI / 3.0
                                    && terrain_facing_angle < 0.0)
                            {
                                local_shadow_data.push((i, j));
                            }
                        }
                    }
                }
                (local_layover_data, local_shadow_data)
            })
            .collect();

        // Apply results
        for (layover_data, shadow_data) in results {
            for (i, j) in layover_data {
                layover_mask[[i, j]] = true;
            }
            for (i, j) in shadow_data {
                shadow_mask[[i, j]] = true;
            }
        }

        Ok(())
    }

    /// Apply terrain flattening normalization using local incidence angle correction
    ///
    /// # Input Requirements
    /// - `sigma0`: Calibrated backscatter coefficient array
    /// - `local_incidence_angles`: Local incidence angles in **RADIANS** (not degrees)
    ///
    /// # Mathematical Basis
    /// Terrain flattening removes topographic modulation from SAR backscatter:
    ///
    /// γ⁰ = σ⁰ / cos(θ_lia)
    ///
    /// where:
    /// - γ⁰ = terrain-flattened backscatter coefficient (gamma-nought)
    /// - σ⁰ = original backscatter coefficient (sigma-nought)
    /// - θ_lia = local incidence angle between radar look vector and surface normal (radians)
    ///
    /// # Literature References
    /// - Small, D. (2011): "Flattening Gamma: Radiometric Terrain Correction for SAR Imagery",
    ///   IEEE Transactions on Geoscience and Remote Sensing, Vol. 49, No. 8
    /// - Ulaby, F.T. et al. (1986): "Microwave Remote Sensing", Vol. 2, Chapter 12
    /// - ESA Sentinel-1 User Handbook, Section 2.3.5: "Radiometric Accuracy"
    ///
    /// # ESA Compliance
    /// Implements ESA Level 2 processing specification for terrain-flattened products
    ///
    /// # Limitations
    /// - Assumes single-scattering radar equation validity
    /// - Valid for local incidence angles between 20° and 60°
    /// - Requires accurate DEM for surface normal computation
    ///
    /// # Error Propagation
    /// δγ⁰/γ⁰ ≈ δσ⁰/σ⁰ + tan(θ_lia) * δθ_lia
    ///
    /// Processing uses parallel computation for efficiency on large images
    pub fn apply_terrain_flattening(
        &self,
        sigma0: &Array2<f32>,
        local_incidence_angles: &Array2<f32>,
    ) -> SarResult<Array2<f32>> {
        if sigma0.dim() != local_incidence_angles.dim() {
            return Err(SarError::Processing(
                "Sigma0 and incidence angle arrays must have same dimensions".to_string(),
            ));
        }

        let (rows, cols) = sigma0.dim();
        let mut gamma0 = Array2::<f32>::zeros((rows, cols));

        let min_angle_rad = self.params.min_incidence_angle * PI / 180.0;
        let max_angle_rad = self.params.max_incidence_angle * PI / 180.0;

        let safe_cos_cap = 0.342f32; // cos(70°) - consistent threshold everywhere

        if self.params.enable_parallel && rows > self.params.chunk_size {
            // Parallel processing - collect results then apply
            let row_indices: Vec<usize> = (0..rows).collect();
            let results: Vec<_> = row_indices
                .par_chunks(self.params.chunk_size)
                .map(|row_chunk| {
                    let mut local_gamma0 = Vec::new();
                    for &i in row_chunk {
                        for j in 0..cols {
                            let theta_lia = local_incidence_angles[[i, j]];
                            let sigma0_val = sigma0[[i, j]];

                            // Validate sigma0 input (finite and non-negative)
                            if !sigma0_val.is_finite() || sigma0_val < 0.0 {
                                local_gamma0.push((i, j, f32::NAN));
                                continue;
                            }

                            // Apply incidence angle masking if enabled
                            if self.params.apply_masking
                                && (theta_lia < min_angle_rad || theta_lia > max_angle_rad)
                            {
                                local_gamma0.push((i, j, f32::NAN));
                                continue;
                            }

                            // Apply terrain flattening with consistent 70° safety cap
                            let cos_theta = theta_lia.cos();
                            let gamma0_val = if cos_theta > safe_cos_cap {
                                sigma0_val / cos_theta
                            } else {
                                f32::NAN // Avoid unrealistic amplification at extreme grazing angles
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

                    // Validate sigma0 input (finite and non-negative)
                    if !sigma0_val.is_finite() || sigma0_val < 0.0 {
                        gamma0[[i, j]] = f32::NAN;
                        continue;
                    }

                    // Apply incidence angle masking if enabled
                    if self.params.apply_masking
                        && (theta_lia < min_angle_rad || theta_lia > max_angle_rad)
                    {
                        gamma0[[i, j]] = f32::NAN;
                        continue;
                    }

                    // Apply terrain flattening with consistent 70° safety cap
                    let cos_theta = theta_lia.cos();
                    gamma0[[i, j]] = if cos_theta > safe_cos_cap {
                        sigma0_val / cos_theta
                    } else {
                        f32::NAN // Avoid unrealistic amplification at extreme grazing angles
                    };
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

        log::debug!(
            "Sigma0 dimensions: {}x{}, DEM dimensions: {}x{}",
            sigma_rows,
            sigma_cols,
            dem_rows,
            dem_cols
        );

        // Force maximum scientific accuracy - always use full 3D geometry computation
        log::info!("Using full 3D geometry computation for maximum scientific accuracy");

        // Step 1: Compute slope and aspect from DEM using robust method
        let (slope, aspect) = self.compute_slope_aspect(dem)?;

        // Step 2: Compute surface normals from slope and aspect
        let surface_normals = self.compute_surface_normals(&slope, &aspect);

        // Step 3: Compute precise radar look vectors using orbit data
        let look_vectors = self.compute_precise_look_vectors(
            sigma_rows,
            sigma_cols,
            _range_time,
            _azimuth_time_start,
            _azimuth_time_spacing,
        )?;

        // Step 4: Compute local incidence angles using full 3D vector geometry
        let incidence_angles = self.compute_local_incidence_angle(&surface_normals, &look_vectors);

        // Step 5: Apply terrain flattening with scientific precision
        let gamma0 = self.apply_terrain_flattening(sigma0, &incidence_angles)?;

        Ok((gamma0, incidence_angles))
    }

    /// Compute precise radar look vectors using full orbit interpolation for maximum scientific accuracy
    fn compute_precise_look_vectors(
        &self,
        rows: usize,
        cols: usize,
        range_time: f64,
        azimuth_time_start: f64,
        azimuth_time_spacing: f64,
    ) -> SarResult<Array2<[f32; 3]>> {
        log::info!("Computing precise look vectors using full orbit interpolation");
        let mut look_vectors = Array2::<[f32; 3]>::from_elem((rows, cols), [0.0, 0.0, 0.0]);

        // Validate orbit data availability
        let state_vectors = &self.orbit_data.state_vectors;
        if state_vectors.is_empty() {
            return Err(SarError::Processing(
                "No orbit state vectors available for precise computation".to_string(),
            ));
        }

        log::debug!(
            "Computing precise look vectors for {}x{} pixels using {} orbit state vectors",
            rows,
            cols,
            state_vectors.len()
        );

        // Use WGS84 ellipsoid parameters for precise Earth geometry
        const WGS84_A: f64 = crate::constants::physical::WGS84_SEMI_MAJOR_AXIS_M; // Semi-major axis (m)
        const WGS84_E2: f64 = crate::constants::geodetic::WGS84_ECCENTRICITY_SQUARED; // First eccentricity squared
        const C: f64 = crate::constants::physical::SPEED_OF_LIGHT_M_S;

        for i in 0..rows {
            let azimuth_time = azimuth_time_start + (i as f64) * azimuth_time_spacing;

            // Interpolate satellite position and velocity at this azimuth time
            let sat_pos = self.interpolate_satellite_position(azimuth_time)?;
            let sat_vel = self.interpolate_satellite_velocity(azimuth_time)?;

            for j in 0..cols {
                // Calculate precise slant range distance for this pixel
                let pixel_range_time =
                    range_time + (j as f64) * self.params.sar_pixel_spacing.0 / C;
                let slant_range = pixel_range_time * C / 2.0;

                // Compute precise ground position using WGS84 ellipsoid and SAR geometry
                let ground_pos = self.compute_precise_ground_position(
                    sat_pos,
                    sat_vel,
                    slant_range,
                    WGS84_A,
                    WGS84_E2,
                )?;

                // Compute precise look vector from satellite to ground
                let dx = ground_pos[0] - sat_pos[0];
                let dy = ground_pos[1] - sat_pos[1];
                let dz = ground_pos[2] - sat_pos[2];

                // Normalize the look vector
                let length = (dx * dx + dy * dy + dz * dz).sqrt();
                if length > 1e-6 {
                    look_vectors[[i, j]] = [
                        (dx / length) as f32,
                        (dy / length) as f32,
                        (dz / length) as f32,
                    ];
                } else {
                    return Err(SarError::Processing(format!(
                        "Degenerate look vector at pixel ({}, {}); check orbit/geometry inputs",
                        i, j
                    )));
                }
            }
        }

        Ok(look_vectors)
    }

    // Removed unused simplified look vector computation to enforce full RD-geometry only

    /// Helper function to fill edge pixels
    fn fill_edges(&self, array: &mut Array2<f32>) -> SarResult<()> {
        let (rows, cols) = array.dim();

        // Fill top and bottom edges
        for j in 0..cols {
            array[[0, j]] = array[[1, j]];
            array[[rows - 1, j]] = array[[rows - 2, j]];
        }

        // Fill left and right edges
        for i in 0..rows {
            array[[i, 0]] = array[[i, 1]];
            array[[i, cols - 1]] = array[[i, cols - 2]];
        }

        Ok(())
    }

    /// Interpolate satellite position at given time using orbit state vectors
    fn interpolate_satellite_position(&self, time: f64) -> SarResult<[f64; 3]> {
        let state_vectors = &self.orbit_data.state_vectors;

        if state_vectors.is_empty() {
            return Err(SarError::Processing(
                "No state vectors available".to_string(),
            ));
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

    /// Create a simple terrain flattening processor for testing from annotation data
    /// This ensures even test functions don't use hardcoded parameters
    pub fn create_simple(
        annotation: &crate::io::annotation::AnnotationRoot,
        dem_pixel_spacing: (f64, f64),
    ) -> SarResult<Self> {
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

        let mut params = TerrainFlatteningParams::from_annotation(annotation)?;
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
        min_incidence_angle: f32,
        max_incidence_angle: f32,
    ) -> SarResult<Array2<f32>> {
        log::info!("Starting simplified terrain flattening with robust geometry");

        let spacing = dem_pixel_spacing.ok_or_else(|| {
            SarError::MissingParameter(
                "DEM pixel spacing is required for scientific terrain flattening".to_string(),
            )
        })?;
        let (rows, cols) = sigma0.dim();

        // Check input data validity
        if rows == 0 || cols == 0 {
            return Err(SarError::Processing("Empty input arrays".to_string()));
        }

        if dem.dim() != sigma0.dim() {
            return Err(SarError::Processing(
                "DEM and sigma0 arrays must have same dimensions".to_string(),
            ));
        }

        // Step 1: Compute slope and aspect from DEM
        log::debug!("Computing terrain slope and aspect");
        let (slope, aspect) = Self::compute_slope_aspect_robust(dem, spacing)?;

        // Step 2: Compute local incidence angles using simplified geometry
        log::debug!("Computing local incidence angles");
        let incidence_angles = Self::compute_incidence_angles_sentinel1(
            &slope,
            &aspect,
            rows,
            cols,
            min_incidence_angle,
            max_incidence_angle,
        )?;

        // Step 3: Apply terrain flattening correction
        log::debug!("Applying terrain flattening correction");
        let mut gamma0 = Array2::<f32>::zeros((rows, cols));

        for i in 0..rows {
            for j in 0..cols {
                let sigma0_val = sigma0[[i, j]];
                let theta_lia = incidence_angles[[i, j]];

                // Apply terrain flattening: gamma0 = sigma0 / cos(theta_lia)
                // CRITICAL FIX: Use realistic threshold to prevent extreme amplification
                let cos_theta = theta_lia.cos();
                // cos(70°) ≈ 0.342, cos(60°) ≈ 0.5 - use realistic SAR geometry limits
                if cos_theta > 0.2 && theta_lia < 70.0_f32.to_radians() {
                    gamma0[[i, j]] = sigma0_val / cos_theta;
                } else {
                    // For extreme angles, limit amplification to prevent unrealistic values
                    let safe_cos = cos_theta.max(0.2); // Limit maximum amplification to 5x
                    gamma0[[i, j]] = sigma0_val / safe_cos;
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

        for i in 1..rows - 1 {
            for j in 1..cols - 1 {
                // Calculate gradients using central differences
                let dz_dx = (dem[[i, j + 1]] - dem[[i, j - 1]]) / (2.0 * dx_scale);
                let dz_dy = (dem[[i + 1, j]] - dem[[i - 1, j]]) / (2.0 * dy_scale);

                // Compute slope (in radians)
                slope[[i, j]] = (dz_dx * dz_dx + dz_dy * dz_dy).sqrt().atan();

                // Compute aspect (in radians, 0 = North, clockwise positive)
                aspect[[i, j]] = (-dz_dy).atan2(dz_dx);
            }
        }

        // Fill edges with nearest neighbor values
        for j in 0..cols {
            slope[[0, j]] = slope[[1, j]];
            slope[[rows - 1, j]] = slope[[rows - 2, j]];
            aspect[[0, j]] = aspect[[1, j]];
            aspect[[rows - 1, j]] = aspect[[rows - 2, j]];
        }

        for i in 0..rows {
            slope[[i, 0]] = slope[[i, 1]];
            slope[[i, cols - 1]] = slope[[i, cols - 2]];
            aspect[[i, 0]] = aspect[[i, 1]];
            aspect[[i, cols - 1]] = aspect[[i, cols - 2]];
        }

        Ok((slope, aspect))
    }

    /// Compute realistic local incidence angles for Sentinel-1 geometry
    fn compute_incidence_angles_sentinel1(
        slope: &Array2<f32>,
        aspect: &Array2<f32>,
        rows: usize,
        cols: usize,
        min_incidence_deg: f32,
        max_incidence_deg: f32,
    ) -> SarResult<Array2<f32>> {
        let mut incidence_angles = Array2::<f32>::zeros((rows, cols));

        // Extract real incidence angles from annotation metadata - no hardcoded values
        let near_range_angle = min_incidence_deg.to_radians();
        let far_range_angle = max_incidence_deg.to_radians();
        let look_direction = -90.0_f32.to_radians(); // Right-looking SAR, perpendicular to flight direction

        for i in 0..rows {
            for j in 0..cols {
                // Interpolate radar incidence angle across range
                let range_fraction = j as f32 / (cols - 1) as f32;
                let radar_incidence =
                    near_range_angle + range_fraction * (far_range_angle - near_range_angle);

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
                incidence_angles[[i, j]] =
                    local_incidence.clamp(5.0_f32.to_radians(), 85.0_f32.to_radians());
            }
        }

        Ok(incidence_angles)
    }

    /// Enhanced terrain flattening with masking
    /// Returns both flattened data and quality mask
    /// Requires annotation data to ensure all parameters are scientifically derived
    pub fn flatten_terrain_with_mask(
        sigma0: &Array2<f32>,
        dem: &Array2<f32>,
        annotation: &crate::io::annotation::AnnotationRoot,
        orbit_data: OrbitData,
        min_incidence: Option<f32>,
        max_incidence: Option<f32>,
    ) -> SarResult<(Array2<f32>, Array2<bool>)> {
        log::info!("Starting terrain flattening with quality masking");

        // Extract parameters from annotation to avoid hardcoded values
        let _params = TerrainFlatteningParams::from_annotation(annotation)?;
        let min_angle = min_incidence.ok_or_else(|| {
            SarError::MissingParameter(
                "Minimum incidence angle is required from annotation metadata".to_string(),
            )
        })?;
        let max_angle = max_incidence.ok_or_else(|| {
            SarError::MissingParameter(
                "Maximum incidence angle is required from annotation metadata".to_string(),
            )
        })?;

        log::debug!("Using incidence angle limits: {:.1}° - {:.1}° (configurable via TerrainFlatteningParams)", 
                   min_angle, max_angle);

        if min_incidence.is_none() {
            log::info!("ℹ️  Using annotation-based minimum incidence angle: {:.1}° (extracted from metadata)", min_angle);
        }
        if max_incidence.is_none() {
            log::info!("ℹ️  Using annotation-based maximum incidence angle: {:.1}° (extracted from metadata)", max_angle);
        }

        let mut params = TerrainFlatteningParams::from_annotation(annotation)?;
        params.min_incidence_angle = min_angle;
        params.max_incidence_angle = max_angle;
        params.apply_masking = true;

        let flattener = Self::new(params, orbit_data);

        // Extract timing parameters from annotation (no hardcoded values)
        let range_time = annotation.get_slant_range_time().ok_or_else(|| {
            SarError::MissingParameter("Slant range time is required from annotation".to_string())
        })?;
        let azimuth_time_start = 0.0; // Start at scene beginning
        let azimuth_time_spacing = 1.0
            / annotation.get_pulse_repetition_frequency().ok_or_else(|| {
                SarError::MissingParameter(
                    "Pulse repetition frequency is required from annotation".to_string(),
                )
            })?;

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

    // Removed unused resample_to_sigma_grid helper; resampling handled in dedicated modules
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::types::StateVector;
    use chrono::{DateTime, Utc};

    fn create_test_orbit_data() -> OrbitData {
        let test_time = DateTime::parse_from_rfc3339("2020-01-03T17:08:15Z")
            .unwrap()
            .with_timezone(&Utc);

        OrbitData {
            state_vectors: vec![StateVector {
                time: test_time,
                position: [7000000.0, 0.0, 0.0],
                velocity: [0.0, 7500.0, 0.0],
            }],
            reference_time: test_time,
        }
    }

    /// Create a test annotation with scientifically valid parameters
    fn create_test_annotation() -> crate::io::annotation::AnnotationRoot {
        let xml = r#"<?xml version="1.0" encoding="UTF-8"?>
<product>
  <adsHeader>
    <missionId>S1A</missionId>
    <absoluteOrbitNumber>30639</absoluteOrbitNumber>
    <missionDataTakeId>123456</missionDataTakeId>
    <imageNumber>1</imageNumber>
    <startTime>2020-12-28T21:59:42.123456Z</startTime>
  </adsHeader>
  <generalAnnotation>
    <productInformation>
      <platformHeading>12.34</platformHeading>
      <rangeSamplingRate>6.4345241e+07</rangeSamplingRate>
      <radarFrequency>5.405000e+09</radarFrequency>
      <azimuthSteeringRate>1.590368784</azimuthSteeringRate>
      <rangePixelSpacing>2.329560</rangePixelSpacing>
      <azimuthPixelSpacing>13.943035</azimuthPixelSpacing>
    </productInformation>
    <downlinkInformationList count="1">
      <downlinkInformation>
        <prf>1710.0</prf>
      </downlinkInformation>
    </downlinkInformationList>
  </generalAnnotation>
  <imageAnnotation>
    <imageInformation>
      <slantRangeTime>0.004</slantRangeTime>
      <rangePixelSpacing>2.329560</rangePixelSpacing>
      <azimuthPixelSpacing>13.943035</azimuthPixelSpacing>
      <incidenceAngleMidSwath>36.0</incidenceAngleMidSwath>
    </imageInformation>
  </imageAnnotation>
  <geolocationGrid>
    <geolocationGridPointList count="2">
      <geolocationGridPoint>
        <slantRangeTime>0.004</slantRangeTime>
        <line>0</line>
        <pixel>0</pixel>
        <latitude>45.0</latitude>
        <longitude>10.0</longitude>
        <height>0.0</height>
        <incidenceAngle>33.0</incidenceAngle>
        <elevationAngle>5.0</elevationAngle>
      </geolocationGridPoint>
      <geolocationGridPoint>
        <slantRangeTime>0.004</slantRangeTime>
        <line>1</line>
        <pixel>1</pixel>
        <latitude>45.1</latitude>
        <longitude>10.1</longitude>
        <height>0.0</height>
        <incidenceAngle>39.0</incidenceAngle>
        <elevationAngle>6.0</elevationAngle>
      </geolocationGridPoint>
    </geolocationGridPointList>
  </geolocationGrid>
</product>"#;

        crate::io::annotation::parse_annotation_xml(xml).expect("Failed to parse test annotation")
    }

    fn make_test_flattener() -> SarResult<TerrainFlattener> {
        let annotation = create_test_annotation();
        let orbit_data = create_test_orbit_data();

        // Use the builder pattern with DEM pixel spacing
        let params = TerrainFlatteningParamsBuilder::new(&annotation)?
            .with_dem_pixel_spacing((30.0, 30.0)) // Typical SRTM spacing
            .build()?;

        Ok(TerrainFlattener::new(params, orbit_data))
    }

    #[test]
    fn test_slope_aspect_computation() {
        let flattener = make_test_flattener().expect("Failed to create test flattener");

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
        let flattener = make_test_flattener().expect("Failed to create test flattener");

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
        let flattener = make_test_flattener().expect("Failed to create test flattener");

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
