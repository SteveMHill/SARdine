//! Extra-Wide (EW) mode deburst placeholder.
//!
//! EW mode deburst is not yet implemented. This module returns a graceful error
//! instead of panicking when EW mode data is encountered.

use crate::types::{SarError, SarResult};

/// Check if the mode is EW and return an appropriate error.
///
/// EW (Extra-Wide) mode uses a different TOPSAR configuration than IW mode:
/// - 5 subswaths instead of 3
/// - Different burst timing and overlap patterns
/// - Wider coverage but lower resolution
///
/// # Returns
/// Always returns `Err` with a descriptive message explaining that EW mode
/// is not yet supported.
pub fn deburst_ew_placeholder() -> SarResult<()> {
    Err(SarError::NotImplemented(
        "EW (Extra-Wide) mode deburst is not yet implemented. \
        Please use IW (Interferometric Wide) mode data, or contribute \
        an EW deburst implementation. See TOPSAR documentation for \
        EW-specific burst patterns and timing."
            .to_string(),
    ))
}

/// Check if a product is EW mode based on its identifier
pub fn is_ew_mode(product_id: &str) -> bool {
    product_id.contains("_EW_") || product_id.contains("_EW1") || product_id.contains("_EW2")
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_ew_mode_returns_error() {
        let result = deburst_ew_placeholder();
        assert!(result.is_err());
        let err_msg = result.unwrap_err().to_string();
        assert!(err_msg.contains("EW"));
        assert!(err_msg.contains("not yet implemented"));
    }

    #[test]
    fn test_is_ew_mode() {
        assert!(is_ew_mode("S1A_EW_GRDH_1SDH_20201005"));
        assert!(is_ew_mode("S1B_EW1_SLC__1SDV_20201005"));
        assert!(!is_ew_mode("S1A_IW_SLC__1SDV_20201005"));
    }
}
