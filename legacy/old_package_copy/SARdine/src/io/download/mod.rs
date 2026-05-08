//! Unified download module for Sentinel-1 products, orbit files, and DEM data
//!
//! This module provides a centralized, user-friendly interface for downloading
//! all types of data required for SAR processing.

pub mod cache;
pub mod dem;
pub mod manager;
pub mod orbit;
pub mod products;
pub mod progress;
pub mod providers;
pub mod queue;
pub mod utils;

#[cfg(test)]
mod tests;

pub use manager::{DownloadConfig, DownloadManager, DownloadResult};
pub use progress::ProgressCallback;
