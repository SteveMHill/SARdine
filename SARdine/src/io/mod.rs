//! I/O modules for reading SAR data, orbits, and DEMs

pub mod annotation; // Comprehensive XML parser with regex-based approach
pub mod dem;
pub mod orbit;
pub mod parsing_validation; // Validation tools for parsing consistency
pub mod sentinel1;
pub mod slc_reader; // Sentinel-1 data download functionality

#[cfg(test)]
pub mod parsing_tests; // Comprehensive test suite for parsing validation

pub use dem::DemReader;
pub use orbit::OrbitReader;
pub use sentinel1::{
    download_by_product_names, DataProvider, SearchParams, Sentinel1Downloader, Sentinel1Product,
};
pub use slc_reader::SlcReader;
