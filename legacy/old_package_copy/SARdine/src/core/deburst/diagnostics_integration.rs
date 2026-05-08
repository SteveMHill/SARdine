#![allow(dead_code, unused_variables)]
/*!
 * STEP-2 Diagnostics Integration Helper
 *
 * This module provides integration points to hook STEP-2 diagnostics
 * into the existing deburst pipeline without changing core algorithms.
 *
 * NOTE: Full integration requires updating DeburstResult and DeburstPlan structures.
 * For now, this provides the basic hooks that can be called manually.
 */

use super::diagnostics::*;
use crate::types::SarResult;

/// Simplified helper for collecting diagnostics
/// Can be used standalone without modifying existing pipeline structures
pub struct DeburstDiagnosticsCollector {
    diagnostics: Step2Diagnostics,
    subswath_name: String,
    polarization: String,
}

impl DeburstDiagnosticsCollector {
    /// Create new collector for a subswath
    pub fn new(
        product_id: String,
        subswath_name: String,
        polarization: String,
        config: DiagnosticsConfig,
    ) -> Self {
        Self {
            diagnostics: Step2Diagnostics::new(product_id, config),
            subswath_name,
            polarization,
        }
    }
    
    /// Check if enabled
    pub fn is_enabled(&self) -> bool {
        self.diagnostics.is_enabled()
    }
    
    /// Get config reference
    pub fn config(&self) -> &DiagnosticsConfig {
        &self.diagnostics.config
    }
    
    /// Access inner diagnostics for manual data collection
    pub fn inner_mut(&mut self) -> &mut Step2Diagnostics {
        &mut self.diagnostics
    }
    
    /// Finalize and write diagnostics
    pub fn finalize(&self, filename: &str) -> SarResult<()> {
        self.diagnostics.write_json(filename)
    }
}

// Helper functions for manual data collection will be added in future iterations
// For now, users can access diagnostics.inner_mut() and call the compute_* functions directly
