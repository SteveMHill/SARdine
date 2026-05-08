//! Merge execution modules for intensity and complex data.

pub mod complex;
pub mod intensity;

#[allow(unused_imports)]
pub use complex::{
    execute_merge_plan_complex, normalize_complex_by_weight_sum, precompute_phase_alignment,
};
#[allow(unused_imports)]
pub use intensity::{
    apply_null_handling, execute_merge_plan, fill_thin_gaps, normalize_by_weight_sum,
};
