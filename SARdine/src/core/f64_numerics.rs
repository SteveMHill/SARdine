/// F64 numerics module for high-precision geometry, Doppler, and Newton-Raphson calculations
/// 
/// This module implements the "Numerics: All geometry/Doppler/NR in f64" refactor from the practical refactors:
/// All geometry/Doppler/NR in f64; only convert to f32 at the last I/O boundary.

use crate::core::type_safe_units::{Radians, Degrees, Meters, Seconds, Hertz, MetersPerSecond};
use crate::types::{SarResult, SarError};

// ============================================================================
// HIGH-PRECISION GEOMETRIC TYPES
// ============================================================================

/// High-precision 3D point for geometric calculations
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Point3D {
    pub x: f64,
    pub y: f64,
    pub z: f64,
}

/// High-precision 2D point for image coordinates
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Point2D {
    pub x: f64,
    pub y: f64,
}

/// High-precision geographic coordinate
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GeographicCoord {
    pub latitude: Degrees,
    pub longitude: Degrees,
    pub elevation: Meters,
}

/// High-precision radar coordinate
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RadarCoord {
    pub azimuth_time: Seconds,
    pub slant_range: Meters,
    pub incidence_angle: Radians,
}

/// High-precision image coordinate
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ImageCoord {
    pub line: f64,
    pub pixel: f64,
}

// ============================================================================
// GEOMETRIC OPERATIONS - ALL F64
// ============================================================================

/// High-precision geometric operations
pub struct GeometryF64;

impl GeometryF64 {
    /// Calculate distance between two 3D points
    pub fn distance_3d(p1: Point3D, p2: Point3D) -> Meters {
        let dx = p1.x - p2.x;
        let dy = p1.y - p2.y;
        let dz = p1.z - p2.z;
        Meters::new((dx * dx + dy * dy + dz * dz).sqrt())
    }
    
    /// Calculate distance between two 2D points
    pub fn distance_2d(p1: Point2D, p2: Point2D) -> f64 {
        let dx = p1.x - p2.x;
        let dy = p1.y - p2.y;
        (dx * dx + dy * dy).sqrt()
    }
    
    /// Cross product of two 3D vectors
    pub fn cross_product(a: Point3D, b: Point3D) -> Point3D {
        Point3D {
            x: a.y * b.z - a.z * b.y,
            y: a.z * b.x - a.x * b.z,
            z: a.x * b.y - a.y * b.x,
        }
    }
    
    /// Dot product of two 3D vectors
    pub fn dot_product(a: Point3D, b: Point3D) -> f64 {
        a.x * b.x + a.y * b.y + a.z * b.z
    }
    
    /// Normalize a 3D vector
    pub fn normalize(v: Point3D) -> SarResult<Point3D> {
        let length = (v.x * v.x + v.y * v.y + v.z * v.z).sqrt();
        if length < f64::EPSILON {
            return Err(SarError::NumericalError("Cannot normalize zero vector".to_string()));
        }
        
        Ok(Point3D {
            x: v.x / length,
            y: v.y / length,
            z: v.z / length,
        })
    }
    
    /// Calculate angle between two 3D vectors
    pub fn angle_between(a: Point3D, b: Point3D) -> SarResult<Radians> {
        let norm_a = Self::normalize(a)?;
        let norm_b = Self::normalize(b)?;
        let dot = Self::dot_product(norm_a, norm_b);
        
        // Clamp to handle numerical precision issues
    use crate::core::global_clamp_debug::ClampDebug;
    let clamped_dot = dot.dbg_clamp(-1.0, 1.0, "dot_product_normalization");
        Ok(Radians::new(clamped_dot.acos()))
    }
    
    /// Rotate point around Z-axis
    pub fn rotate_z(point: Point3D, angle: Radians) -> Point3D {
        let cos_a = angle.value().cos();
        let sin_a = angle.value().sin();
        
        Point3D {
            x: point.x * cos_a - point.y * sin_a,
            y: point.x * sin_a + point.y * cos_a,
            z: point.z,
        }
    }
    
    /// Convert spherical to Cartesian coordinates
    pub fn spherical_to_cartesian(
        radius: Meters,
        azimuth: Radians,
        elevation: Radians,
    ) -> Point3D {
        let r = radius.value();
        let az = azimuth.value();
        let el = elevation.value();
        
        let cos_el = el.cos();
        
        Point3D {
            x: r * cos_el * az.cos(),
            y: r * cos_el * az.sin(),
            z: r * el.sin(),
        }
    }
    
    /// Convert Cartesian to spherical coordinates
    pub fn cartesian_to_spherical(point: Point3D) -> (Meters, Radians, Radians) {
        let radius = Meters::new((point.x * point.x + point.y * point.y + point.z * point.z).sqrt());
        let azimuth = Radians::new(point.y.atan2(point.x));
        let elevation = Radians::new(point.z.atan2((point.x * point.x + point.y * point.y).sqrt()));
        
        (radius, azimuth, elevation)
    }
}

// ============================================================================
// DOPPLER CALCULATIONS - ALL F64
// ============================================================================

/// High-precision Doppler calculations
pub struct DopplerF64;

impl DopplerF64 {
    /// Calculate Doppler frequency from relative velocity
    /// 
    /// # Arguments
    /// * `relative_velocity` - Relative velocity between radar and target
    /// * `wavelength` - Radar wavelength in meters
    /// 
    /// # Returns
    /// * Doppler frequency in Hz
    pub fn frequency_from_velocity(
        relative_velocity: MetersPerSecond,
        wavelength: Meters,
    ) -> Hertz {
        let doppler_hz = -2.0 * relative_velocity.value() / wavelength.value();
        Hertz::new(doppler_hz)
    }
    
    /// Calculate relative velocity from Doppler frequency
    pub fn velocity_from_frequency(
        doppler_freq: Hertz,
        wavelength: Meters,
    ) -> MetersPerSecond {
        let velocity_ms = -doppler_freq.value() * wavelength.value() / 2.0;
        MetersPerSecond::new(velocity_ms)
    }
    
    /// Calculate range migration for given Doppler parameters
    /// 
    /// This uses the full f64 precision for accurate range migration calculation
    /// which is critical for focusing quality.
    pub fn range_migration(
        azimuth_time: Seconds,
        reference_time: Seconds,
        slant_range: Meters,
        wavelength: Meters,
        doppler_centroid: Hertz,
        fm_rate: f64,
    ) -> Meters {
        let dt = azimuth_time.value() - reference_time.value();
        let lambda = wavelength.value();
        let r0 = slant_range.value();
        let fdc = doppler_centroid.value();
        
        // Range migration formula: ΔR = λ²/(8πR₀) * (fDC + fR*t)² * t²
        let doppler_term = fdc + fm_rate * dt;
        let range_migration_m = (lambda * lambda) / (8.0 * std::f64::consts::PI * r0) * 
                               doppler_term * doppler_term * dt * dt;
        
        Meters::new(range_migration_m)
    }
    
    /// Calculate azimuth frequency modulation rate
    /// 
    /// This is a critical parameter for SAR focusing that requires high precision
    pub fn azimuth_fm_rate(
        velocity: MetersPerSecond,
        slant_range: Meters,
        wavelength: Meters,
    ) -> f64 {
        let v = velocity.value();
        let r = slant_range.value();
        let lambda = wavelength.value();
        
        // FM rate = -2*v²/(λ*R)
        -2.0 * v * v / (lambda * r)
    }
    
    /// Calculate squint angle from Doppler centroid
    pub fn squint_from_doppler(
        doppler_centroid: Hertz,
        velocity: MetersPerSecond,
        wavelength: Meters,
    ) -> Radians {
        let fdc = doppler_centroid.value();
        let v = velocity.value();
        let lambda = wavelength.value();
        
        // Squint angle: θ = asin(λ*fDC/(2*v))
        let sin_squint = (lambda * fdc) / (2.0 * v);
    use crate::core::global_clamp_debug::ClampDebug;
    let clamped_sin = sin_squint.dbg_clamp(-1.0, 1.0, "sin_squint_normalization");
        
        Radians::new(clamped_sin.asin())
    }
    
    /// Doppler bandwidth calculation
    pub fn doppler_bandwidth(
        prf: Hertz,
        processed_bandwidth_fraction: f64,
    ) -> Hertz {
        Hertz::new(prf.value() * processed_bandwidth_fraction)
    }
    
    /// Calculate theoretical resolution in azimuth
    pub fn azimuth_resolution(
        velocity: MetersPerSecond,
        doppler_bandwidth: Hertz,
    ) -> Meters {
        let v = velocity.value();
        let bd = doppler_bandwidth.value();
        
        // Theoretical azimuth resolution: ρa = v/Bd
        Meters::new(v / bd)
    }
}

// ============================================================================
// NEWTON-RAPHSON SOLVER - ALL F64
// ============================================================================

/// High-precision Newton-Raphson solver for SAR coordinate transformations
/// 
/// This solver maintains f64 precision throughout to prevent accumulation
/// of numerical errors in iterative solutions.
pub struct NewtonRaphsonF64 {
    /// Maximum number of iterations
    max_iterations: usize,
    /// Convergence tolerance (absolute)
    tolerance: f64,
    /// Convergence tolerance (relative)
    relative_tolerance: f64,
}

/// Function trait for Newton-Raphson
pub trait NRFunction {
    /// Evaluate function at point x
    fn evaluate(&self, x: f64) -> f64;
    
    /// Evaluate derivative at point x
    fn derivative(&self, x: f64) -> f64;
}

/// Multi-dimensional function trait for Newton-Raphson
pub trait NRFunctionND {
    /// Number of dimensions
    fn dimensions(&self) -> usize;
    
    /// Evaluate function vector at point x
    fn evaluate(&self, x: &[f64]) -> Vec<f64>;
    
    /// Evaluate Jacobian matrix at point x
    fn jacobian(&self, x: &[f64]) -> Vec<Vec<f64>>;
}

impl Default for NewtonRaphsonF64 {
    fn default() -> Self {
        Self {
            max_iterations: 100,
            tolerance: 1e-12,
            relative_tolerance: 1e-10,
        }
    }
}

impl NewtonRaphsonF64 {
    /// Create new Newton-Raphson solver with custom parameters
    pub fn new(max_iterations: usize, tolerance: f64, relative_tolerance: f64) -> Self {
        Self {
            max_iterations,
            tolerance,
            relative_tolerance,
        }
    }
    
    /// Solve single-variable equation f(x) = 0
    pub fn solve_1d<F: NRFunction>(
        &self,
        function: &F,
        initial_guess: f64,
    ) -> SarResult<f64> {
        let mut x = initial_guess;
        let mut iteration = 0;
        
        while iteration < self.max_iterations {
            let f_val = function.evaluate(x);
            let df_val = function.derivative(x);
            
            // Check for zero derivative
            if df_val.abs() < f64::EPSILON {
                return Err(SarError::NumericalError(
                    format!("Zero derivative at x = {:.12}", x)
                ));
            }
            
            // Newton-Raphson step
            let delta = f_val / df_val;
            let x_new = x - delta;
            
            // Check convergence
            let abs_error = delta.abs();
            let rel_error = if x.abs() > f64::EPSILON {
                abs_error / x.abs()
            } else {
                abs_error
            };
            
            if abs_error < self.tolerance || rel_error < self.relative_tolerance {
                return Ok(x_new);
            }
            
            x = x_new;
            iteration += 1;
        }
        
        Err(SarError::NumericalError(
            format!("Newton-Raphson failed to converge after {} iterations", self.max_iterations)
        ))
    }
    
    /// Solve multi-dimensional system F(x) = 0
    pub fn solve_nd<F: NRFunctionND>(
        &self,
        function: &F,
        initial_guess: &[f64],
    ) -> SarResult<Vec<f64>> {
        let n = function.dimensions();
        if initial_guess.len() != n {
            return Err(SarError::InvalidInput(
                "Initial guess dimension mismatch".to_string()
            ));
        }
        
        let mut x = initial_guess.to_vec();
        let mut iteration = 0;
        
        while iteration < self.max_iterations {
            let f_vals = function.evaluate(&x);
            let jacobian = function.jacobian(&x);
            
            // Solve linear system: J * delta = -F
            let delta = self.solve_linear_system(&jacobian, &f_vals)?;
            
            // Update solution
            for i in 0..n {
                x[i] -= delta[i];
            }
            
            // Check convergence
            let norm_delta = Self::vector_norm(&delta);
            let norm_x = Self::vector_norm(&x);
            
            let abs_error = norm_delta;
            let rel_error = if norm_x > f64::EPSILON {
                abs_error / norm_x
            } else {
                abs_error
            };
            
            if abs_error < self.tolerance || rel_error < self.relative_tolerance {
                return Ok(x);
            }
            
            iteration += 1;
        }
        
        Err(SarError::NumericalError(
            format!("Multi-dimensional Newton-Raphson failed to converge after {} iterations", 
                   self.max_iterations)
        ))
    }
    
    /// Solve linear system Ax = b using Gaussian elimination with partial pivoting
    fn solve_linear_system(&self, a: &[Vec<f64>], b: &[f64]) -> SarResult<Vec<f64>> {
        let n = a.len();
        if n == 0 || a[0].len() != n || b.len() != n {
            return Err(SarError::InvalidInput("Invalid matrix dimensions".to_string()));
        }
        
        // Create augmented matrix
        let mut aug = vec![vec![0.0; n + 1]; n];
        for i in 0..n {
            for j in 0..n {
                aug[i][j] = a[i][j];
            }
            aug[i][n] = b[i];
        }
        
        // Forward elimination with partial pivoting
        for i in 0..n {
            // Find pivot
            let mut max_row = i;
            for k in (i + 1)..n {
                if aug[k][i].abs() > aug[max_row][i].abs() {
                    max_row = k;
                }
            }
            
            // Swap rows
            if max_row != i {
                aug.swap(i, max_row);
            }
            
            // Check for zero pivot
            if aug[i][i].abs() < f64::EPSILON {
                return Err(SarError::NumericalError("Singular matrix in linear solve".to_string()));
            }
            
            // Eliminate column
            for k in (i + 1)..n {
                let factor = aug[k][i] / aug[i][i];
                for j in i..=n {
                    aug[k][j] -= factor * aug[i][j];
                }
            }
        }
        
        // Back substitution
        let mut x = vec![0.0; n];
        for i in (0..n).rev() {
            x[i] = aug[i][n];
            for j in (i + 1)..n {
                x[i] -= aug[i][j] * x[j];
            }
            x[i] /= aug[i][i];
        }
        
        Ok(x)
    }
    
    /// Calculate L2 norm of vector
    fn vector_norm(v: &[f64]) -> f64 {
        v.iter().map(|&x| x * x).sum::<f64>().sqrt()
    }
}

// ============================================================================
// SAR-SPECIFIC COORDINATE TRANSFORM FUNCTIONS
// ============================================================================

/// Range-Doppler to geographic coordinate transformation
pub struct RangeDopplerToGeo {
    /// Orbital state vectors
    orbit_position: Point3D,
    orbit_velocity: Point3D,
    /// Radar parameters
    wavelength: Meters,
    /// Earth model parameters
    earth_radius: Meters,
    earth_flattening: f64,
}

impl NRFunctionND for RangeDopplerToGeo {
    fn dimensions(&self) -> usize {
        3 // [latitude, longitude, elevation]
    }
    
    fn evaluate(&self, coords: &[f64]) -> Vec<f64> {
        let lat = Radians::new(coords[0]);
        let lon = Radians::new(coords[1]);
        let elev = Meters::new(coords[2]);
        
        // Convert to ECEF
        let target_ecef = self.geo_to_ecef(lat, lon, elev);
        
        // Range equation: |target - satellite| - measured_range = 0
        let range_residual = GeometryF64::distance_3d(target_ecef, self.orbit_position).value()
                           - self.measured_range().value();
        
        // Doppler equation: dot(relative_velocity, los_unit) * 2/λ - measured_doppler = 0
        let relative_vel = Point3D {
            x: -self.orbit_velocity.x, // Target velocity assumed zero
            y: -self.orbit_velocity.y,
            z: -self.orbit_velocity.z,
        };
        
        let los_vector = Point3D {
            x: target_ecef.x - self.orbit_position.x,
            y: target_ecef.y - self.orbit_position.y,
            z: target_ecef.z - self.orbit_position.z,
        };
        
        let los_unit = GeometryF64::normalize(los_vector).unwrap_or(Point3D { x: 0.0, y: 0.0, z: 1.0 });
        let doppler_velocity = GeometryF64::dot_product(relative_vel, los_unit);
        let computed_doppler = -2.0 * doppler_velocity / self.wavelength.value();
        let doppler_residual = computed_doppler - self.measured_doppler().value();
        
        // Elevation constraint (if using DEM)
        let elevation_residual = elev.value() - self.reference_elevation().value();
        
        vec![range_residual, doppler_residual, elevation_residual]
    }
    
    fn jacobian(&self, coords: &[f64]) -> Vec<Vec<f64>> {
        // Numerical differentiation for Jacobian
        let delta = 1e-8;
        let f0 = self.evaluate(coords);
        let mut jacobian = vec![vec![0.0; 3]; 3];
        
        for j in 0..3 {
            let mut coords_plus = coords.to_vec();
            coords_plus[j] += delta;
            let f_plus = self.evaluate(&coords_plus);
            
            for i in 0..3 {
                jacobian[i][j] = (f_plus[i] - f0[i]) / delta;
            }
        }
        
        jacobian
    }
}

impl RangeDopplerToGeo {
    /// Convert geographic to ECEF coordinates
    fn geo_to_ecef(&self, lat: Radians, lon: Radians, elev: Meters) -> Point3D {
        let a = self.earth_radius.value();
        let f = self.earth_flattening;
        let e2 = 2.0 * f - f * f; // First eccentricity squared
        
        let sin_lat = lat.value().sin();
        let cos_lat = lat.value().cos();
        let sin_lon = lon.value().sin();
        let cos_lon = lon.value().cos();
        
        let n = a / (1.0 - e2 * sin_lat * sin_lat).sqrt(); // Radius of curvature in prime vertical
        
        Point3D {
            x: (n + elev.value()) * cos_lat * cos_lon,
            y: (n + elev.value()) * cos_lat * sin_lon,
            z: (n * (1.0 - e2) + elev.value()) * sin_lat,
        }
    }
    
    // These would be set from actual measurements
    fn measured_range(&self) -> Meters { Meters::new(850000.0) }
    fn measured_doppler(&self) -> Hertz { Hertz::new(-7193.46) }
    fn reference_elevation(&self) -> Meters { Meters::new(0.0) }
}

// ============================================================================
// PRECISION CONVERSION UTILITIES
// ============================================================================

/// Utilities for converting between f64 (internal) and f32 (I/O boundary)
pub struct PrecisionConverter;

impl PrecisionConverter {
    /// Convert f64 array to f32 for I/O (with bounds checking)
    pub fn f64_to_f32_array(input: &[f64]) -> SarResult<Vec<f32>> {
        let mut output = Vec::with_capacity(input.len());
        
        for &value in input {
            if value.is_finite() {
                let f32_val = value as f32;
                // Check for precision loss warning
                if (value - f32_val as f64).abs() / value.abs() > 1e-6 {
                    eprintln!("Warning: Significant precision loss converting {} to f32", value);
                }
                output.push(f32_val);
            } else {
                return Err(SarError::NumericalError(
                    format!("Non-finite value in f64 to f32 conversion: {}", value)
                ));
            }
        }
        
        Ok(output)
    }
    
    /// Convert f32 array to f64 for internal calculations
    pub fn f32_to_f64_array(input: &[f32]) -> Vec<f64> {
        input.iter().map(|&x| x as f64).collect()
    }
    
    /// Convert complex f64 to complex f32 for I/O
    pub fn complex_f64_to_f32(real: f64, imag: f64) -> SarResult<(f32, f32)> {
        if !real.is_finite() || !imag.is_finite() {
            return Err(SarError::NumericalError(
                "Non-finite value in complex conversion".to_string()
            ));
        }
        
        Ok((real as f32, imag as f32))
    }
    
    /// Validate that f64 values are suitable for f32 conversion
    pub fn validate_f32_conversion_range(values: &[f64]) -> SarResult<()> {
        for &value in values {
            if !value.is_finite() {
                return Err(SarError::NumericalError(
                    format!("Non-finite value: {}", value)
                ));
            }
            
            if value.abs() > f32::MAX as f64 {
                return Err(SarError::NumericalError(
                    format!("Value {} exceeds f32 range", value)
                ));
            }
            
            if value != 0.0 && value.abs() < f32::MIN_POSITIVE as f64 {
                return Err(SarError::NumericalError(
                    format!("Value {} below f32 minimum positive", value)
                ));
            }
        }
        
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_point3d_operations() {
        let p1 = Point3D { x: 1.0, y: 2.0, z: 3.0 };
        let p2 = Point3D { x: 4.0, y: 5.0, z: 6.0 };
        
        let distance = GeometryF64::distance_3d(p1, p2);
        assert!((distance.value() - ((3.0_f64.powi(2) * 3.0).sqrt())).abs() < 1e-10);
        
        let dot = GeometryF64::dot_product(p1, p2);
        assert_eq!(dot, 32.0); // 1*4 + 2*5 + 3*6 = 32
        
        let cross = GeometryF64::cross_product(p1, p2);
        assert_eq!(cross.x, -3.0); // 2*6 - 3*5 = -3
        assert_eq!(cross.y, 6.0);  // 3*4 - 1*6 = 6
        assert_eq!(cross.z, -3.0); // 1*5 - 2*4 = -3
    }
    
    #[test]
    fn test_doppler_calculations() {
        let velocity = MetersPerSecond::new(7500.0); // Typical satellite velocity
        let wavelength = Meters::new(0.055); // C-band wavelength
        
        let doppler = DopplerF64::frequency_from_velocity(velocity, wavelength);
        let expected_doppler = -2.0 * 7500.0 / 0.055;
        assert!((doppler.value() - expected_doppler).abs() < 1e-6);
        
        let recovered_velocity = DopplerF64::velocity_from_frequency(doppler, wavelength);
        assert!((recovered_velocity.value() - velocity.value()).abs() < 1e-6);
    }
    
    #[test]
    fn test_newton_raphson_1d() {
        // Solve x^2 - 4 = 0 (should give x = 2)
        struct QuadraticFunction;
        impl NRFunction for QuadraticFunction {
            fn evaluate(&self, x: f64) -> f64 {
                x * x - 4.0
            }
            
            fn derivative(&self, x: f64) -> f64 {
                2.0 * x
            }
        }
        
        let solver = NewtonRaphsonF64::default();
        let function = QuadraticFunction;
        let result = solver.solve_1d(&function, 1.0).unwrap();
        
        assert!((result - 2.0).abs() < 1e-10);
    }
    
    #[test]
    fn test_precision_conversion() {
        let f64_values = vec![1.23456789, -9.87654321, 0.0, 1e-20, 1e20];
        
        // Test conversion with validation
        assert!(PrecisionConverter::validate_f32_conversion_range(&f64_values[0..3]).is_ok());
        
        let f32_converted = PrecisionConverter::f64_to_f32_array(&f64_values[0..3]).unwrap();
        assert_eq!(f32_converted.len(), 3);
        
        let f64_recovered = PrecisionConverter::f32_to_f64_array(&f32_converted);
        assert_eq!(f64_recovered.len(), 3);
        
        // Test out-of-range detection
        let large_values = vec![f64::MAX];
        assert!(PrecisionConverter::validate_f32_conversion_range(&large_values).is_err());
    }
    
    #[test]
    fn test_spherical_conversion() {
        let radius = Meters::new(100.0);
        let azimuth = Radians::new(std::f64::consts::PI / 4.0); // 45 degrees
        let elevation = Radians::new(std::f64::consts::PI / 6.0); // 30 degrees
        
        let cartesian = GeometryF64::spherical_to_cartesian(radius, azimuth, elevation);
        let (r_back, az_back, el_back) = GeometryF64::cartesian_to_spherical(cartesian);
        
        assert!((r_back.value() - radius.value()).abs() < 1e-10);
        assert!((az_back.value() - azimuth.value()).abs() < 1e-10);
        assert!((el_back.value() - elevation.value()).abs() < 1e-10);
    }
}