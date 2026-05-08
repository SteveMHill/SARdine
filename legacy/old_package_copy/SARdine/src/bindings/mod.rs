//! Python bindings for SARdine - organized by functionality
//!
//! This module contains PyO3 bindings organized into submodules.
//! Currently in incremental migration from lib.rs.
//!
//! ## Migration Status
//! - `utils`: Helper functions (migrated)
//! - `calibration`: Calibration types and helpers (partial)
//! - `reader`: PySlcReader struct definition (partial)
//! - `orbit`: Orbit file loading and application (migrated)
//! - `export`: GeoTIFF export and metadata generation (migrated)
//! - Others: Placeholder modules for future migration

pub mod calibration;
pub mod export;
pub mod orbit;
pub mod reader;
pub mod utils;

// Placeholder modules - will be populated in future refactoring
pub mod deburst;
pub mod dem;
pub mod download;
pub mod merge;
pub mod processing;
pub mod validation;

// Re-export pyclass types
pub use calibration::PyCalibrationJob;
pub use reader::PySlcReader;

// Re-export pyfunctions from utils
pub use utils::{
    convert_to_db_real, db_to_linear_inplace_py, estimate_num_looks, export_db_parallel_py,
    extract_ellipsoid_incidence_angle, extract_platform_heading, get_product_info,
    get_product_info_cached, latlon_to_ecef, linear_to_db_inplace_py,
};

// Re-export pyfunctions from orbit
pub use orbit::{apply_precise_orbit_file, load_orbit_file};

// Re-export pyfunctions from export
pub use export::{export_cog_with_stac, export_geotiff, export_metadata_json, export_metadata_xml};
