//! I/O modules for reading SAR data, orbits, and DEMs

pub mod slc_reader;
pub mod orbit;
pub mod dem;
pub mod annotation;  // Comprehensive XML parser with regex-based approach

pub use slc_reader::SlcReader;
pub use orbit::OrbitReader;
pub use dem::DemReader;
