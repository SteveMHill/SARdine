/// Type-safe units for SAR processing to prevent unit confusion and violations
/// 
/// This module implements the "Type-safe units" refactor from the practical refactors:
/// Newtypes for Radians, Degrees, Meters, Seconds; functions accept only the correct type.

use std::fmt;
use std::ops::{Add, Sub, Mul, Div, AddAssign, SubAssign};

// ============================================================================
// ANGULAR UNITS - Type-safe radians and degrees
// ============================================================================

/// Type-safe radians to prevent angle unit confusion
#[derive(Debug, Clone, Copy, PartialEq, PartialOrd)]
pub struct Radians(pub f64);

/// Type-safe degrees to prevent angle unit confusion
#[derive(Debug, Clone, Copy, PartialEq, PartialOrd)]
pub struct Degrees(pub f64);

impl Radians {
    pub fn new(value: f64) -> Self {
        Self(value)
    }
    
    pub fn value(&self) -> f64 {
        self.0
    }
    
    pub fn to_degrees(&self) -> Degrees {
        Degrees(self.0 * 180.0 / std::f64::consts::PI)
    }
    
    pub fn sin(&self) -> f64 {
        self.0.sin()
    }
    
    pub fn cos(&self) -> f64 {
        self.0.cos()
    }
    
    pub fn tan(&self) -> f64 {
        self.0.tan()
    }
    
    pub fn clamp(self, min: Radians, max: Radians) -> Radians {
        use crate::core::global_clamp_debug::ClampDebug;
        if min.0 > max.0 {
            let strict = std::env::var("SARDINE_STRICT_CLAMP").ok().as_deref() == Some("1");
            let log_bt = std::env::var("SARDINE_LOG_CLAMP_BT").ok().as_deref() == Some("1");
            if log_bt {
                let bt = std::backtrace::Backtrace::force_capture();
                log::error!(
                    "🚨 Radians clamp inversion: value={} min={} max={} strict={}\n{:?}",
                    self.0, min.0, max.0, strict, bt
                );
            } else {
                log::error!(
                    "🚨 Radians clamp inversion: value={} min={} max={} strict={}",
                    self.0, min.0, max.0, strict
                );
            }
            if strict {
                panic!("Radians clamp inversion (strict mode)");
            }
            // Swap using dbg_clamp so global instrumentation still records bounds (label notes swap)
            return Radians(self.0.dbg_clamp(max.0, min.0, "Radians_inversion_swapped"));
        }
        Radians(self.0.dbg_clamp(min.0, max.0, "Radians"))
    }
    
    /// Create radians from degrees (explicit conversion)
    pub fn from_degrees(degrees: Degrees) -> Self {
        Self(degrees.0 * std::f64::consts::PI / 180.0)
    }
}

impl Degrees {
    pub fn new(value: f64) -> Self {
        Self(value)
    }
    
    pub fn value(&self) -> f64 {
        self.0
    }
    
    pub fn to_radians(&self) -> Radians {
        Radians(self.0 * std::f64::consts::PI / 180.0)
    }
    
    pub fn clamp(self, min: Degrees, max: Degrees) -> Degrees {
        use crate::core::global_clamp_debug::ClampDebug;
        if min.0 > max.0 {
            let strict = std::env::var("SARDINE_STRICT_CLAMP").ok().as_deref() == Some("1");
            let log_bt = std::env::var("SARDINE_LOG_CLAMP_BT").ok().as_deref() == Some("1");
            if log_bt {
                let bt = std::backtrace::Backtrace::force_capture();
                log::error!(
                    "🚨 Degrees clamp inversion: value={} min={} max={} strict={}\n{:?}",
                    self.0, min.0, max.0, strict, bt
                );
            } else {
                log::error!(
                    "🚨 Degrees clamp inversion: value={} min={} max={} strict={}",
                    self.0, min.0, max.0, strict
                );
            }
            if strict {
                panic!("Degrees clamp inversion (strict mode)");
            }
            return Degrees(self.0.dbg_clamp(max.0, min.0, "Degrees_inversion_swapped"));
        }
        Degrees(self.0.dbg_clamp(min.0, max.0, "Degrees"))
    }
}

impl fmt::Display for Radians {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{:.6} rad", self.0)
    }
}

impl fmt::Display for Degrees {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{:.3}°", self.0)
    }
}

// Math operations for Radians
impl Add for Radians {
    type Output = Self;
    fn add(self, other: Self) -> Self {
        Self(self.0 + other.0)
    }
}

impl Sub for Radians {
    type Output = Self;
    fn sub(self, other: Self) -> Self {
        Self(self.0 - other.0)
    }
}

impl Mul<f64> for Radians {
    type Output = Self;
    fn mul(self, scalar: f64) -> Self {
        Self(self.0 * scalar)
    }
}

impl Div<f64> for Radians {
    type Output = Self;
    fn div(self, scalar: f64) -> Self {
        Self(self.0 / scalar)
    }
}

// Math operations for Degrees
impl Add for Degrees {
    type Output = Self;
    fn add(self, other: Self) -> Self {
        Self(self.0 + other.0)
    }
}

impl Sub for Degrees {
    type Output = Self;
    fn sub(self, other: Self) -> Self {
        Self(self.0 - other.0)
    }
}

impl Mul<f64> for Degrees {
    type Output = Self;
    fn mul(self, scalar: f64) -> Self {
        Self(self.0 * scalar)
    }
}

impl Div<f64> for Degrees {
    type Output = Self;
    fn div(self, scalar: f64) -> Self {
        Self(self.0 / scalar)
    }
}

// ============================================================================
// DISTANCE UNITS - Type-safe meters
// ============================================================================

/// Type-safe meters to prevent distance unit confusion
#[derive(Debug, Clone, Copy, PartialEq, PartialOrd)]
pub struct Meters(pub f64);

impl Meters {
    pub fn new(value: f64) -> Self {
        Self(value)
    }
    
    pub fn value(&self) -> f64 {
        self.0
    }
    
    pub fn to_kilometers(&self) -> f64 {
        self.0 / 1000.0
    }
    
    pub fn from_kilometers(km: f64) -> Self {
        Self(km * 1000.0)
    }
    
    pub fn abs(&self) -> Self {
        Self(self.0.abs())
    }
    
    pub fn clamp(self, min: Meters, max: Meters) -> Meters {
        use crate::core::global_clamp_debug::ClampDebug;
        if min.0 > max.0 {
            let strict = std::env::var("SARDINE_STRICT_CLAMP").ok().as_deref() == Some("1");
            let log_bt = std::env::var("SARDINE_LOG_CLAMP_BT").ok().as_deref() == Some("1");
            if log_bt {
                let bt = std::backtrace::Backtrace::force_capture();
                log::error!(
                    "🚨 Meters clamp inversion: value={} min={} max={} strict={}\n{:?}",
                    self.0, min.0, max.0, strict, bt
                );
            } else {
                log::error!(
                    "🚨 Meters clamp inversion: value={} min={} max={} strict={}",
                    self.0, min.0, max.0, strict
                );
            }
            if strict {
                panic!("Meters clamp inversion (strict mode)");
            }
            return Meters(self.0.dbg_clamp(max.0, min.0, "Meters_inversion_swapped"));
        }
        Meters(self.0.dbg_clamp(min.0, max.0, "Meters"))
    }
}

impl fmt::Display for Meters {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        if self.0.abs() >= 1000.0 {
            write!(f, "{:.3} km", self.0 / 1000.0)
        } else {
            write!(f, "{:.3} m", self.0)
        }
    }
}

// Math operations for Meters
impl Add for Meters {
    type Output = Self;
    fn add(self, other: Self) -> Self {
        Self(self.0 + other.0)
    }
}

impl Sub for Meters {
    type Output = Self;
    fn sub(self, other: Self) -> Self {
        Self(self.0 - other.0)
    }
}

impl Mul<f64> for Meters {
    type Output = Self;
    fn mul(self, scalar: f64) -> Self {
        Self(self.0 * scalar)
    }
}

impl Div<f64> for Meters {
    type Output = Self;
    fn div(self, scalar: f64) -> Self {
        Self(self.0 / scalar)
    }
}

impl Div for Meters {
    type Output = f64;
    fn div(self, other: Self) -> f64 {
        self.0 / other.0
    }
}

impl AddAssign for Meters {
    fn add_assign(&mut self, other: Self) {
        self.0 += other.0;
    }
}

impl SubAssign for Meters {
    fn sub_assign(&mut self, other: Self) {
        self.0 -= other.0;
    }
}

// ============================================================================
// TIME UNITS - Type-safe seconds
// ============================================================================

/// Type-safe seconds to prevent time unit confusion
#[derive(Debug, Clone, Copy, PartialEq, PartialOrd)]
pub struct Seconds(pub f64);

impl Seconds {
    pub fn new(value: f64) -> Self {
        Self(value)
    }
    
    pub fn value(&self) -> f64 {
        self.0
    }
    
    pub fn to_milliseconds(&self) -> f64 {
        self.0 * 1000.0
    }
    
    pub fn from_milliseconds(ms: f64) -> Self {
        Self(ms / 1000.0)
    }
    
    pub fn abs(&self) -> Self {
        Self(self.0.abs())
    }
    
    pub fn clamp(self, min: Seconds, max: Seconds) -> Seconds {
        use crate::core::global_clamp_debug::ClampDebug;
        if min.0 > max.0 {
            let strict = std::env::var("SARDINE_STRICT_CLAMP").ok().as_deref() == Some("1");
            let log_bt = std::env::var("SARDINE_LOG_CLAMP_BT").ok().as_deref() == Some("1");
            if log_bt {
                let bt = std::backtrace::Backtrace::force_capture();
                log::error!(
                    "🚨 Seconds clamp inversion: value={} min={} max={} strict={}\n{:?}",
                    self.0, min.0, max.0, strict, bt
                );
            } else {
                log::error!(
                    "🚨 Seconds clamp inversion: value={} min={} max={} strict={}",
                    self.0, min.0, max.0, strict
                );
            }
            if strict {
                panic!("Seconds clamp inversion (strict mode)");
            }
            return Seconds(self.0.dbg_clamp(max.0, min.0, "Seconds_inversion_swapped"));
        }
        Seconds(self.0.dbg_clamp(min.0, max.0, "Seconds"))
    }
}

impl fmt::Display for Seconds {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        if self.0.abs() < 1.0 {
            write!(f, "{:.3} ms", self.0 * 1000.0)
        } else {
            write!(f, "{:.6} s", self.0)
        }
    }
}

// Math operations for Seconds
impl Add for Seconds {
    type Output = Self;
    fn add(self, other: Self) -> Self {
        Self(self.0 + other.0)
    }
}

impl Sub for Seconds {
    type Output = Self;
    fn sub(self, other: Self) -> Self {
        Self(self.0 - other.0)
    }
}

impl Mul<f64> for Seconds {
    type Output = Self;
    fn mul(self, scalar: f64) -> Self {
        Self(self.0 * scalar)
    }
}

impl Div<f64> for Seconds {
    type Output = Self;
    fn div(self, scalar: f64) -> Self {
        Self(self.0 / scalar)
    }
}

impl Div for Seconds {
    type Output = f64;
    fn div(self, other: Self) -> f64 {
        self.0 / other.0
    }
}

// ============================================================================
// DERIVED UNITS - Type-safe combinations
// ============================================================================

/// Type-safe meters per second
#[derive(Debug, Clone, Copy, PartialEq, PartialOrd)]
pub struct MetersPerSecond(pub f64);

impl MetersPerSecond {
    pub fn new(value: f64) -> Self {
        Self(value)
    }
    
    pub fn value(&self) -> f64 {
        self.0
    }
}

impl fmt::Display for MetersPerSecond {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{:.3} m/s", self.0)
    }
}

impl Div<Seconds> for Meters {
    type Output = MetersPerSecond;
    fn div(self, time: Seconds) -> MetersPerSecond {
        MetersPerSecond(self.0 / time.0)
    }
}

impl Mul<Seconds> for MetersPerSecond {
    type Output = Meters;
    fn mul(self, time: Seconds) -> Meters {
        Meters(self.0 * time.0)
    }
}

/// Type-safe hertz (frequency)
#[derive(Debug, Clone, Copy, PartialEq, PartialOrd)]
pub struct Hertz(pub f64);

impl Hertz {
    pub fn new(value: f64) -> Self {
        Self(value)
    }
    
    pub fn value(&self) -> f64 {
        self.0
    }
    
    pub fn to_period(&self) -> Seconds {
        Seconds(1.0 / self.0)
    }
    
    pub fn from_period(period: Seconds) -> Self {
        Self(1.0 / period.0)
    }
}

impl fmt::Display for Hertz {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        if self.0.abs() >= 1e6 {
            write!(f, "{:.3} MHz", self.0 / 1e6)
        } else if self.0.abs() >= 1e3 {
            write!(f, "{:.3} kHz", self.0 / 1e3)
        } else {
            write!(f, "{:.3} Hz", self.0)
        }
    }
}

// ============================================================================
// UTILITY FUNCTIONS - Type-safe unit conversions
// ============================================================================

/// Type-safe angle operations
pub mod angle_ops {
    use super::{Radians, Degrees, Meters};
    
    /// Calculate local incidence angle from radar and terrain geometry
    pub fn calculate_local_incidence_angle(
        radar_incidence: Radians,
        terrain_slope: Radians,
        aspect_factor: f64,
    ) -> Radians {
        let terrain_correction = terrain_slope * aspect_factor;
        radar_incidence + Radians::new(terrain_correction.value())
    }
    
    /// Clamp incidence angle to valid SAR range
    pub fn clamp_incidence_angle(angle: Radians) -> Radians {
        let min_incidence = Degrees::new(5.0).to_radians();
        let max_incidence = Degrees::new(85.0).to_radians();
        angle.clamp(min_incidence, max_incidence)
    }
    
    /// Convert WGS84 lat/lon to ECEF coordinates using type-safe angles
    pub fn latlon_to_ecef(lat: Degrees, lon: Degrees, elevation_m: Meters) -> [f64; 3] {
        let lat_rad = lat.to_radians();
        let lon_rad = lon.to_radians();
        let h = elevation_m.value();
        
        // WGS84 parameters
        const A: f64 = 6_378_137.0; // Semi-major axis (m)
        const E2: f64 = 0.006_694_379_990_14; // First eccentricity squared
        
        let sin_lat = lat_rad.sin();
        let cos_lat = lat_rad.cos();
        let sin_lon = lon_rad.sin();
        let cos_lon = lon_rad.cos();
        
        let n = A / (1.0 - E2 * sin_lat * sin_lat).sqrt();
        
        let x = (n + h) * cos_lat * cos_lon;
        let y = (n + h) * cos_lat * sin_lon;
        let z = (n * (1.0 - E2) + h) * sin_lat;
        
        [x, y, z]
    }
}

/// Type-safe distance operations
pub mod distance_ops {
    use super::Meters;
    
    /// Calculate euclidean distance between two 3D points
    pub fn euclidean_distance_3d(p1: &[f64; 3], p2: &[f64; 3]) -> Meters {
        let dx = p1[0] - p2[0];
        let dy = p1[1] - p2[1];
        let dz = p1[2] - p2[2];
        Meters::new((dx * dx + dy * dy + dz * dz).sqrt())
    }
    
    /// Calculate CE90 error from position offsets
    pub fn calculate_ce90_error(position_errors: &[Meters]) -> Meters {
        if position_errors.is_empty() {
            return Meters::new(0.0);
        }
        
        let mut errors: Vec<f64> = position_errors.iter().map(|e| e.value()).collect();
        errors.sort_by(|a, b| a.partial_cmp(b).unwrap());
        
        // 90th percentile
        let index = ((errors.len() as f64) * 0.9).floor() as usize;
        let index = index.min(errors.len() - 1);
        
        Meters::new(errors[index])
    }
}

/// Type-safe time operations
pub mod time_ops {
    use super::{Seconds, Hertz};
    
    /// Convert PRF to azimuth time interval
    pub fn prf_to_time_interval(prf: Hertz) -> Seconds {
        prf.to_period()
    }
    
    /// Calculate relative time between two absolute times
    pub fn time_difference(start: Seconds, end: Seconds) -> Seconds {
        end - start
    }
    
    /// Clamp time to valid burst duration
    pub fn clamp_to_burst_duration(time: Seconds, burst_duration: Seconds) -> Seconds {
        time.clamp(Seconds::new(0.0), burst_duration)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_angle_conversions() {
        let deg = Degrees::new(45.0);
        let rad = deg.to_radians();
        assert!((rad.value() - std::f64::consts::PI / 4.0).abs() < 1e-10);
        
        let back_to_deg = rad.to_degrees();
        assert!((back_to_deg.value() - 45.0).abs() < 1e-10);
    }
    
    #[test]
    fn test_distance_operations() {
        let d1 = Meters::new(1000.0);
        let d2 = Meters::new(500.0);
        
        let sum = d1 + d2;
        assert_eq!(sum.value(), 1500.0);
        
        let diff = d1 - d2;
        assert_eq!(diff.value(), 500.0);
        
        let scaled = d1 * 2.0;
        assert_eq!(scaled.value(), 2000.0);
    }
    
    #[test]
    fn test_time_operations() {
        let t1 = Seconds::new(2.5);
        let t2 = Seconds::new(1.0);
        
        let diff = t1 - t2;
        assert_eq!(diff.value(), 1.5);
        
        let prf = Hertz::new(486.486);
        let period = prf.to_period();
        assert!((period.value() - 0.002055).abs() < 1e-6);
    }
    
    #[test]
    fn test_derived_units() {
        let distance = Meters::new(100.0);
        let time = Seconds::new(10.0);
        
        let velocity = distance / time;
        assert_eq!(velocity.value(), 10.0);
        
        let back_distance = velocity * time;
        assert_eq!(back_distance.value(), 100.0);
    }
    
    #[test]
    fn test_type_safety() {
        // These should compile and work correctly
        let angle_deg = Degrees::new(30.0);
        let angle_rad = angle_deg.to_radians();
        let _incidence = angle_ops::clamp_incidence_angle(angle_rad);
        
        let dist = Meters::new(1000.0);
        let time = Seconds::new(2.0);
        let _speed = dist / time;
        
        // Type safety prevents mixing units incorrectly
        // These would NOT compile:
        // let bad = angle_deg + dist;  // Cannot add degrees and meters
        // let bad = time * angle_rad;  // Cannot multiply time and angle
    }
    
    #[test]
    fn test_coordinate_conversion() {
        let lat = Degrees::new(45.0);
        let lon = Degrees::new(9.0);
        let elevation = Meters::new(1000.0);
        
        let ecef = angle_ops::latlon_to_ecef(lat, lon, elevation);
        
        // Verify reasonable ECEF coordinates for this location
        assert!(ecef[0] > 4e6 && ecef[0] < 5e6); // X coordinate
        assert!(ecef[1] > 0.5e6 && ecef[1] < 1e6); // Y coordinate
        assert!(ecef[2] > 4e6 && ecef[2] < 5e6); // Z coordinate
    }
}