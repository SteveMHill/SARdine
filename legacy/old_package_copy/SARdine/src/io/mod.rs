//! I/O modules for reading SAR data, orbits, and DEMs

pub mod annotation; // Comprehensive XML parser with regex-based approach
pub mod dem;
pub mod download; // Unified download module for products, orbits, and DEMs
pub mod orbit;
pub mod parsing_validation; // Validation tools for parsing consistency
pub mod sentinel;
pub mod sentinel1; // Sentinel-1 SLC reader facade split into submodules

#[cfg(test)]
pub mod parsing_tests; // Comprehensive test suite for parsing validation

pub use dem::DemReader;
pub use download::{DownloadConfig, DownloadManager, DownloadResult};
pub use orbit::OrbitReader;
pub use sentinel::SlcReader;
pub use sentinel1::{
    download_by_product_names, DataProvider, SearchParams, Sentinel1Downloader, Sentinel1Product,
};
