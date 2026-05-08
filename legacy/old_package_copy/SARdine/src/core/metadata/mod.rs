//! Metadata submodule: grouped parsers, provenance tracking, and validation gates.
pub mod metadata_parser;
pub mod metadata_provenance;
pub mod metadata_strictness;
pub mod validation_gates;

pub use metadata_parser::*;
pub use metadata_provenance::*;
pub use metadata_strictness::*;
pub use validation_gates::*;
