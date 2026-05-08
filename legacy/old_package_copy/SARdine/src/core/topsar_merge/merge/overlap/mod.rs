//! Overlap module: detection, weight generation, and quality metrics.

pub mod detect;
pub mod weights;

pub use super::types::OverlapRegion;
pub use detect::{audit_overlap_regions, detect_subswath_overlaps_with_blending, overlap_range};
pub use weights::create_complementary_cosine_weights;
