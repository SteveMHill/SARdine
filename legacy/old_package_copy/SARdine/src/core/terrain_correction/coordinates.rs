#![allow(dead_code, unused_variables)]
//! Coordinate transformation utilities for terrain correction
//!
//! Provides coordinate system conversions between:
//! - Geographic (WGS84 lat/lon)
//! - UTM projected coordinates
//! - ECEF (Earth-Centered Earth-Fixed) Cartesian coordinates
//!
//! All transformations use WGS84 ellipsoid parameters.

use wide::f64x4;

use super::types::{LatLon, Vector3};
use super::TerrainCorrector;
use crate::constants::geodetic::WGS84_SEMI_MAJOR_AXIS_M;
use crate::types::{SarError, SarResult};

impl TerrainCorrector {
    /// Convert map coordinates to geographic coordinates with proper CRS handling
    ///
    /// # Coordinate Ordering Convention
    /// - Input: `map_x` = longitude/easting, `map_y` = latitude/northing
    /// - Output: `(latitude, longitude)` - standard geographic ordering
    ///
    /// # Supported CRS
    /// - EPSG:4326 (WGS84 Geographic)
    /// - EPSG:32601-32660 (UTM North zones)
    /// - EPSG:32701-32760 (UTM South zones)
    pub(super) fn map_to_geographic(&self, map_x: f64, map_y: f64) -> SarResult<(f64, f64)> {
        // CRITICAL FIX: Ensure consistent coordinate ordering
        // map_x should be longitude (X/easting)
        // map_y should be latitude (Y/northing)

        // Validate input coordinates are finite
        if !map_x.is_finite() || !map_y.is_finite() {
            return Err(SarError::Processing(format!(
                "Invalid map coordinates: map_x={}, map_y={} (must be finite)",
                map_x, map_y
            )));
        }

        if self.output_crs == 4326 {
            // For Geographic (WGS84): map_x=lon, map_y=lat
            // Validate geographic bounds
            if map_y < -90.0 || map_y > 90.0 {
                return Err(SarError::Processing(format!(
                    "Latitude out of bounds: {} (must be in [-90, 90])",
                    map_y
                )));
            }
            if map_x < -180.0 || map_x > 180.0 {
                return Err(SarError::Processing(format!(
                    "Longitude out of bounds: {} (must be in [-180, 180])",
                    map_x
                )));
            }
            Ok((map_y, map_x)) // Return (latitude, longitude)
        } else if self.output_crs >= 32601 && self.output_crs <= 32760 {
            // UTM zones: need proper projection
            self.utm_to_geographic(map_x, map_y, self.output_crs)
        } else {
            // Unknown CRS - this is a critical error
            Err(SarError::Processing(format!(
                "Unsupported coordinate reference system: EPSG:{}",
                self.output_crs
            )))
        }
    }

    /// Proper UTM to Geographic transformation using WGS84 ellipsoid
    pub(super) fn utm_to_geographic(
        &self,
        easting: f64,
        northing: f64,
        epsg_code: u32,
    ) -> SarResult<(f64, f64)> {
        // Extract UTM zone and hemisphere from EPSG code
        let zone = if epsg_code >= 32601 && epsg_code <= 32660 {
            epsg_code - 32600 // North
        } else if epsg_code >= 32701 && epsg_code <= 32760 {
            epsg_code - 32700 // South
        } else {
            return Err(SarError::Processing("Invalid UTM EPSG code".to_string()));
        };

        let is_north = epsg_code <= 32660;

        // WGS84 ellipsoid parameters
        let a = WGS84_SEMI_MAJOR_AXIS_M;
        let f = crate::constants::geodetic::WGS84_FLATTENING;
        let e2 = 2.0 * f - f * f; // First eccentricity squared

        // UTM projection parameters
        let k0 = 0.9996; // Scale factor
        let false_easting = 500000.0;
        let false_northing = if is_north { 0.0 } else { 10000000.0 };

        // Central meridian for this UTM zone
        let lon0 = ((zone as f64 - 1.0) * 6.0 - 180.0 + 3.0).to_radians();

        // Normalize coordinates
        let x = easting - false_easting;
        let y = northing - false_northing;

        // Iterative solution for latitude using Bowring's method
        let m = y / k0;
        let mu = m / (a * (1.0 - e2 / 4.0 - 3.0 * e2 * e2 / 64.0 - 5.0 * e2 * e2 * e2 / 256.0));

        let e1 = (1.0 - (1.0 - e2).sqrt()) / (1.0 + (1.0 - e2).sqrt());
        let j1 = 3.0 * e1 / 2.0 - 27.0 * e1.powi(3) / 32.0;
        let j2 = 21.0 * e1.powi(2) / 16.0 - 55.0 * e1.powi(4) / 32.0;
        let j3 = 151.0 * e1.powi(3) / 96.0;
        let j4 = 1097.0 * e1.powi(4) / 512.0;

        let fp = mu
            + j1 * (2.0 * mu).sin()
            + j2 * (4.0 * mu).sin()
            + j3 * (6.0 * mu).sin()
            + j4 * (8.0 * mu).sin();

        let e_prime2 = e2 / (1.0 - e2);
        let c1 = e_prime2 * fp.cos().powi(2);
        let t1 = fp.tan().powi(2);
        let r1 = a * (1.0 - e2) / (1.0 - e2 * fp.sin().powi(2)).powf(1.5);
        let n1 = a / (1.0 - e2 * fp.sin().powi(2)).sqrt();
        let d = x / (n1 * k0);

        let lat = fp
            - (n1 * fp.tan() / r1)
                * (d.powi(2) / 2.0
                    - (5.0 + 3.0 * t1 + 10.0 * c1 - 4.0 * c1.powi(2) - 9.0 * e_prime2) * d.powi(4)
                        / 24.0
                    + (61.0 + 90.0 * t1 + 298.0 * c1 + 45.0 * t1.powi(2)
                        - 252.0 * e_prime2
                        - 3.0 * c1.powi(2))
                        * d.powi(6)
                        / 720.0);

        let lon = lon0
            + (d - (1.0 + 2.0 * t1 + c1) * d.powi(3) / 6.0
                + (5.0 - 2.0 * c1 + 28.0 * t1 - 3.0 * c1.powi(2)
                    + 8.0 * e_prime2
                    + 24.0 * t1.powi(2))
                    * d.powi(5)
                    / 120.0)
                / fp.cos();

        Ok((lat.to_degrees(), lon.to_degrees()))
    }

    /// Convert geographic coordinates to projected coordinates
    pub(super) fn geographic_to_projected(&self, lon: f64, lat: f64) -> SarResult<(f64, f64)> {
        if self.output_crs == 4326 {
            // Already geographic
            Ok((lon, lat))
        } else if self.output_crs >= 32601 && self.output_crs <= 32760 {
            // UTM projection - use enhanced transformation with proper coordinate type
            let coords = LatLon::new(lat, lon)?;
            self.enhanced_geographic_to_utm(coords, self.output_crs)
        } else {
            Err(SarError::Processing(format!(
                "Unsupported projection: EPSG:{}",
                self.output_crs
            )))
        }
    }

    /// Convert geographic to UTM coordinates using WGS84 ellipsoid
    ///
    /// Enhanced UTM coordinate transformation with improved accuracy and documentation.
    /// Uses standard UTM projection formulas from NIMA Technical Report TR8350.2
    /// with WGS84 ellipsoid parameters and proper convergence series.
    ///
    /// # Literature References
    /// - NIMA TR8350.2: "Department of Defense World Geodetic System 1984"
    /// - Snyder, J.P. (1987): "Map Projections - A Working Manual", USGS
    pub(super) fn enhanced_geographic_to_utm(
        &self,
        coords: LatLon,
        epsg_code: u32,
    ) -> SarResult<(f64, f64)> {
        // Validate EPSG code and extract zone information
        let (zone, is_north) = if epsg_code >= 32601 && epsg_code <= 32660 {
            (epsg_code - 32600, true) // North
        } else if epsg_code >= 32701 && epsg_code <= 32760 {
            (epsg_code - 32700, false) // South
        } else {
            return Err(SarError::Processing(format!(
                "Unsupported EPSG code: {}. Only UTM zones 1-60 supported",
                epsg_code
            )));
        };

        // Validate coordinates for UTM applicability
        if coords.lat.abs() > 84.0 {
            return Err(SarError::Processing(format!(
                "Latitude {} outside UTM valid range (±84°)",
                coords.lat
            )));
        }

        // Convert to radians
        let lat_rad = coords.lat.to_radians();
        let lon_rad = coords.lon.to_radians();

        // Central meridian for this UTM zone
        let central_meridian = ((zone as f64 - 1.0) * 6.0 - 180.0 + 3.0).to_radians();
        let delta_lon = lon_rad - central_meridian;

        // WGS84 ellipsoid parameters (using constants)
        let a = WGS84_SEMI_MAJOR_AXIS_M;
        let f = crate::constants::geodetic::WGS84_FLATTENING;
        let e2 = 2.0 * f - f * f; // First eccentricity squared
        let ep2 = e2 / (1.0 - e2); // Second eccentricity squared
        let k0 = 0.9996; // UTM scale factor

        let sin_lat = lat_rad.sin();
        let cos_lat = lat_rad.cos();
        let tan_lat = lat_rad.tan();

        // Radius of curvature in the prime vertical
        let nu = a / (1.0 - e2 * sin_lat * sin_lat).sqrt();

        // Terms for UTM projection
        let t = tan_lat * tan_lat;
        let c = ep2 * cos_lat * cos_lat;
        let a_term = cos_lat * delta_lon;

        // Meridional arc length (M)
        let m = a
            * ((1.0 - e2 / 4.0 - 3.0 * e2 * e2 / 64.0 - 5.0 * e2.powi(3) / 256.0) * lat_rad
                - (3.0 * e2 / 8.0 + 3.0 * e2 * e2 / 32.0 + 45.0 * e2.powi(3) / 1024.0)
                    * (2.0 * lat_rad).sin()
                + (15.0 * e2 * e2 / 256.0 + 45.0 * e2.powi(3) / 1024.0) * (4.0 * lat_rad).sin()
                - (35.0 * e2.powi(3) / 3072.0) * (6.0 * lat_rad).sin());

        // UTM Easting with higher-order terms
        let easting = 500000.0
            + k0 * nu
                * (a_term
                    + (1.0 - t + c) * a_term.powi(3) / 6.0
                    + (5.0 - 18.0 * t + t * t + 72.0 * c - 58.0 * ep2) * a_term.powi(5) / 120.0);

        // UTM Northing with higher-order terms
        let northing_base = k0
            * (m + nu
                * tan_lat
                * (a_term.powi(2) / 2.0
                    + (5.0 - t + 9.0 * c + 4.0 * c * c) * a_term.powi(4) / 24.0
                    + (61.0 - 58.0 * t + t * t + 600.0 * c - 330.0 * ep2) * a_term.powi(6)
                        / 720.0));

        // Apply false northing for southern hemisphere
        let northing = if is_north {
            northing_base
        } else {
            northing_base + 10000000.0
        };

        // Validate results are reasonable
        if easting < 166000.0 || easting > 834000.0 {
            log::warn!(
                "UTM easting {} outside typical range [166000, 834000]",
                easting
            );
        }

        Ok((easting, northing))
    }

    /// Transform lat/lon coordinates to DEM coordinate reference system
    pub(super) fn transform_latlon_to_dem_crs(&self, lat: f64, lon: f64) -> SarResult<(f64, f64)> {
        if self.dem_crs == 4326 {
            // DEM is already in WGS84 geographic
            Ok((lon, lat))
        } else if self.dem_crs >= 32601 && self.dem_crs <= 32760 {
            // DEM is in UTM - transform from geographic to UTM with enhanced accuracy
            let coords = LatLon::new(lat, lon)?;
            self.enhanced_geographic_to_utm(coords, self.dem_crs)
        } else {
            Err(SarError::Processing(format!(
                "Unsupported DEM CRS: EPSG:{}",
                self.dem_crs
            )))
        }
    }

    /// Placeholder for future GDAL/PROJ integration
    #[allow(dead_code)]
    pub(super) fn future_gdal_transform(
        &self,
        _coords: LatLon,
        _target_epsg: u32,
    ) -> SarResult<(f64, f64)> {
        Err(SarError::Processing(
            "GDAL integration not yet implemented".to_string(),
        ))
    }

    /// Convert lat/lon/elevation to ECEF coordinates
    pub(super) fn latlon_to_ecef(&self, lat: f64, lon: f64, elevation: f64) -> [f64; 3] {
        let a = crate::constants::geodetic::WGS84_SEMI_MAJOR_AXIS_M;
        let e2 = crate::constants::geodetic::WGS84_ECCENTRICITY_SQUARED;

        let lat_rad = lat.to_radians();
        let lon_rad = lon.to_radians();

        // Pre-compute trigonometric values to avoid redundant calculations
        let (sin_lat, cos_lat) = lat_rad.sin_cos();
        let (sin_lon, cos_lon) = lon_rad.sin_cos();

        // Optimized prime vertical radius calculation
        let sin_lat_sq = sin_lat * sin_lat;
        let n = a / (1.0 - e2 * sin_lat_sq).sqrt();

        // Pre-compute common factor
        let n_plus_h_cos_lat = (n + elevation) * cos_lat;

        let x = n_plus_h_cos_lat * cos_lon;
        let y = n_plus_h_cos_lat * sin_lon;
        let z = (n * (1.0 - e2) + elevation) * sin_lat;

        [x, y, z]
    }

    /// SIMD-optimized batch coordinate transformation for 4 points at once
    pub(super) fn latlon_to_ecef_simd_batch(
        &self,
        lats: &[f64; 4],
        lons: &[f64; 4],
        elevations: &[f64; 4],
    ) -> [[f64; 3]; 4] {
        // Constants for WGS84
        let a = crate::constants::geodetic::WGS84_SEMI_MAJOR_AXIS_M;
        let e2 = crate::constants::geodetic::WGS84_ECCENTRICITY_SQUARED;

        // Load coordinates into SIMD registers
        let lat_simd = f64x4::new([lats[0], lats[1], lats[2], lats[3]]);
        let lon_simd = f64x4::new([lons[0], lons[1], lons[2], lons[3]]);
        let elev_simd = f64x4::new([elevations[0], elevations[1], elevations[2], elevations[3]]);

        // Convert to radians
        let pi_over_180 = f64x4::splat(std::f64::consts::PI / 180.0);
        let lat_rad = lat_simd * pi_over_180;
        let lon_rad = lon_simd * pi_over_180;

        // Calculate trigonometric values using SIMD
        let sin_lat = lat_rad.sin();
        let cos_lat = lat_rad.cos();
        let sin_lon = lon_rad.sin();
        let cos_lon = lon_rad.cos();

        // Calculate prime vertical radius
        let sin_lat_sq = sin_lat * sin_lat;
        let a_simd = f64x4::splat(a);
        let e2_simd = f64x4::splat(e2);
        let one = f64x4::splat(1.0);
        let n = a_simd / (one - e2_simd * sin_lat_sq).sqrt();

        // Calculate ECEF coordinates
        let n_plus_h_cos_lat = (n + elev_simd) * cos_lat;
        let x = n_plus_h_cos_lat * cos_lon;
        let y = n_plus_h_cos_lat * sin_lon;
        let z = (n * (one - e2_simd) + elev_simd) * sin_lat;

        // Extract results from SIMD registers
        let x_array = x.to_array();
        let y_array = y.to_array();
        let z_array = z.to_array();

        [
            [x_array[0], y_array[0], z_array[0]],
            [x_array[1], y_array[1], z_array[1]],
            [x_array[2], y_array[2], z_array[2]],
            [x_array[3], y_array[3], z_array[3]],
        ]
    }

    /// Calculate distance between two 3D points
    pub(super) fn distance_to_point(&self, point1: &[f64; 3], point2: &[f64; 3]) -> f64 {
        let dx = point1[0] - point2[0];
        let dy = point1[1] - point2[1];
        let dz = point1[2] - point2[2];
        (dx * dx + dy * dy + dz * dz).sqrt()
    }

    /// Calculate distance between Vector3 and array point
    pub(super) fn distance_vector3_to_array(&self, point1: &Vector3, point2: &[f64; 3]) -> f64 {
        let dx = point1.x - point2[0];
        let dy = point1.y - point2[1];
        let dz = point1.z - point2[2];
        (dx * dx + dy * dy + dz * dz).sqrt()
    }

    /// Convert Vector3 to array
    pub(super) fn vector3_to_array(&self, vec: &Vector3) -> [f64; 3] {
        [vec.x, vec.y, vec.z]
    }

    /// Convert array to Vector3
    pub(super) fn array_to_vector3(&self, arr: &[f64; 3]) -> Vector3 {
        Vector3 {
            x: arr[0],
            y: arr[1],
            z: arr[2],
        }
    }
}
