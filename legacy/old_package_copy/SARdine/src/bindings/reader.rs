//! PySlcReader class - Python wrapper for SlcReader
//!
//! This is a placeholder - the full implementation remains in lib.rs
//! and will be migrated incrementally.

use pyo3::prelude::*;

/// Python wrapper for the SLC reader
#[pyclass(name = "SlcReader", module = "sardine._core")]
pub struct PySlcReader {
    pub inner: crate::io::SlcReader,
}

// NOTE: The #[pymethods] impl block remains in lib.rs for now
// due to its size (~1000 lines) and would be moved in a subsequent refactor
