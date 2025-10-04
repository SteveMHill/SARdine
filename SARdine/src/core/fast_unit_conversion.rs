/// CONSOLIDATION NOTICE: Unit conversion functions moved to canonical location
///
/// The basic unit conversion functions (`db_to_linear_inplace`, `linear_to_db_inplace`, etc.)
/// have been consolidated into `crate::constants::unit_conversion` to eliminate code duplication.
///
/// This module now contains only the specialized functions:
/// - `IncidenceTerms`: Precomputed incidence angle terms for performance
/// - `Units`: Canonical unit representation enum
///
/// For basic unit conversions, import from `crate::constants::unit_conversion` or `crate::core`
/// which re-exports the canonical implementations.
///
/// # Migration Guide
/// ```ignore
/// // OLD (still works but deprecated):
/// use sardine::core::fast_unit_conversion::{db_to_linear_inplace, linear_to_db_inplace};
///
/// // NEW (recommended):
/// use sardine::constants::unit_conversion::{db_to_linear_inplace, linear_to_db_inplace};
/// // OR:
/// use sardine::core::{db_to_linear_inplace, linear_to_db_inplace};
/// ```
use rayon::prelude::*;

/// Canonical unit representation - convert once, flag forever
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Units {
    Linear, // Already in linear domain
    Db,     // Decibel scale (needs conversion)
}

// ============================================================================
// CONSOLIDATED: Basic unit conversion functions moved to crate::constants::unit_conversion
// ============================================================================
// 
// The following functions are now in `crate::constants::unit_conversion`:
// - db_to_linear_inplace(v: &mut [f32])
// - linear_to_db_inplace(v: &mut [f32])
// - db_to_linear_parallel_chunked(buf: &mut [f32])
// - linear_to_db_parallel_chunked(buf: &mut [f32])
// - db_to_linear_inplace_par(buf: &mut [f32])
// - linear_to_db_inplace_par(buf: &mut [f32])
// - export_db_parallel(image_linear: &mut [f32])
//
// These are re-exported from `crate::core` for convenience.
// ============================================================================

/// Precomputed incidence angle terms to eliminate per-pixel trig
/// Avoids calling sin/cos per pixel during β⁰/γ⁰ derivation
pub struct IncidenceTerms {
    pub inv_sin: Vec<f32>, // 1/sin(alpha) for β⁰ derivation
    pub inv_cos: Vec<f32>, // 1/cos(alpha) for γ⁰ derivation
}

impl IncidenceTerms {
    /// Precompute incidence terms once per column
    /// Eliminates millions of sin/cos calls in hot loops
    pub fn from_incidence_per_column(alpha_rad: &[f32]) -> Self {
        let mut inv_sin = Vec::with_capacity(alpha_rad.len());
        let mut inv_cos = Vec::with_capacity(alpha_rad.len());

        for &a in alpha_rad {
            let (s, c) = a.sin_cos();
            inv_sin.push(1.0 / s.max(1e-6)); // Avoid division by zero
            inv_cos.push(1.0 / c.max(1e-6));
        }

        Self { inv_sin, inv_cos }
    }

    /// Apply β⁰ conversion: β⁰ = σ⁰ / sin(α)
    /// Vectorized elementwise multiplication - no trig in hot loop
    pub fn apply_beta_conversion(&self, sigma0_row: &[f32], beta0_row: &mut [f32]) {
        assert_eq!(sigma0_row.len(), self.inv_sin.len());
        assert_eq!(beta0_row.len(), self.inv_sin.len());

        // Vectorized: β⁰[i] = σ⁰[i] * inv_sin[i]
        for ((sigma, beta), &inv_s) in sigma0_row
            .iter()
            .zip(beta0_row.iter_mut())
            .zip(self.inv_sin.iter())
        {
            *beta = sigma * inv_s;
        }
    }

    /// Apply γ⁰ conversion: γ⁰ = σ⁰ / cos(α)  
    /// Vectorized elementwise multiplication - no trig in hot loop
    pub fn apply_gamma_conversion(&self, sigma0_row: &[f32], gamma0_row: &mut [f32]) {
        assert_eq!(sigma0_row.len(), self.inv_cos.len());
        assert_eq!(gamma0_row.len(), self.inv_cos.len());

        // Vectorized: γ⁰[i] = σ⁰[i] * inv_cos[i]
        for ((sigma, gamma), &inv_c) in sigma0_row
            .iter()
            .zip(gamma0_row.iter_mut())
            .zip(self.inv_cos.iter())
        {
            *gamma = sigma * inv_c;
        }
    }
}

/// Fast unit detection and one-time conversion for calibration vectors
/// Eliminates repeated unit checking in hot loops
pub struct FastUnitProcessor;

impl FastUnitProcessor {
    /// Detect and convert units once at XML load time
    /// Never check units again during LUT fill or pixel lookup
    pub fn normalize_to_linear(values: &mut Vec<f64>, tag: &str, units: Option<&str>) -> Units {
        let detected_units = Self::detect_units(values, units);

        match detected_units {
            Units::Db => {
                Self::convert_f64_db_to_linear(values);
                log::info!(
                    "{}: converted from dB to linear domain ({} values)",
                    tag,
                    values.len()
                );
            }
            Units::Linear => {
                log::info!(
                    "{}: already in linear domain ({} values)",
                    tag,
                    values.len()
                );
            }
        }

        Units::Linear // Always linear after this call
    }

    /// Smart unit detection using heuristics
    fn detect_units(values: &[f64], units_attr: Option<&str>) -> Units {
        // Check explicit unit attribute first
        if let Some(u) = units_attr {
            let u_lower = u.to_ascii_lowercase();
            if matches!(u_lower.as_str(), "db" | "decibel" | "decibels") {
                return Units::Db;
            }
        }

        // Heuristic: deci-dB when mean is in [100..500] range
        // This catches the tenths-dB format from Sentinel-1 calibration vectors
        if !values.is_empty() {
            let mean = values.iter().sum::<f64>() / values.len() as f64;
            if mean > 100.0 && mean < 500.0 {
                return Units::Db; // Likely tenths-dB format
            }
        }

        Units::Linear
    }

    /// Fast f64 dB to linear conversion
    fn convert_f64_db_to_linear(values: &mut [f64]) {
        // For calibration gains: treat as NEGATIVE tenths-dB
        // linear = 10^(-v/100) for proper calibration scaling
        for v in values.iter_mut() {
            *v = 10f64.powf(-*v / 100.0);
        }
    }
}

/// Fast export to dB for final output
/// Keeps pipeline in linear domain until the very end
pub fn export_db_parallel(image_linear: &mut [f32]) {
    crate::constants::unit_conversion::linear_to_db_inplace_par(image_linear);
}

/// Sanity check for flat arrays - avoid converting garbage
pub fn is_flat_array(vals: &[f64], threshold: f64) -> bool {
    if vals.len() < 16 {
        return false;
    }

    let mean = vals.iter().sum::<f64>() / vals.len() as f64;
    let variance: f64 = vals
        .iter()
        .map(|&x| {
            let d = x - mean;
            d * d
        })
        .sum::<f64>()
        / vals.len() as f64;

    let coefficient_of_variation = (variance.sqrt()) / mean.abs().max(1e-12);
    coefficient_of_variation < threshold
}

/// Performance-optimized calibration vector preprocessing
/// Converts all vectors to linear domain once at load time
pub fn preprocess_calibration_vectors(
    sigma_values: &mut Vec<f64>,
    beta_values: &mut Vec<f64>,
    gamma_values: &mut Vec<f64>,
    sigma_units: Option<&str>,
    beta_units: Option<&str>,
    gamma_units: Option<&str>,
) -> (Units, Units, Units) {
    // Check for flat arrays and skip conversion if garbage
    let sigma_flat = is_flat_array(sigma_values, 1e-6);
    let beta_flat = is_flat_array(beta_values, 1e-6);
    let gamma_flat = is_flat_array(gamma_values, 1e-6);

    if sigma_flat {
        log::warn!("Sigma vector is flat - potential XML parsing issue, skipping conversion");
    }
    if beta_flat {
        log::warn!("Beta vector is flat - potential XML parsing issue, forcing σ⁰-derivation path");
    }
    if gamma_flat {
        log::warn!(
            "Gamma vector is flat - potential XML parsing issue, forcing σ⁰-derivation path"
        );
    }

    // Convert once, flag as linear forever
    let sigma_units_final = if !sigma_flat {
        FastUnitProcessor::normalize_to_linear(sigma_values, "sigmaNought", sigma_units)
    } else {
        Units::Linear
    };

    let beta_units_final = if !beta_flat {
        FastUnitProcessor::normalize_to_linear(beta_values, "betaNought", beta_units)
    } else {
        Units::Linear
    };

    let gamma_units_final = if !gamma_flat {
        FastUnitProcessor::normalize_to_linear(gamma_values, "gamma", gamma_units)
    } else {
        Units::Linear
    };

    (sigma_units_final, beta_units_final, gamma_units_final)
}

#[cfg(test)]
mod tests {
    use crate::constants::unit_conversion::{db_to_linear_inplace, linear_to_db_inplace};
    use super::IncidenceTerms;

    #[test]
    fn test_db_to_linear_conversion() {
        let mut values: Vec<f32> = vec![10.0, 20.0, 30.0]; // dB values
        db_to_linear_inplace(&mut values);

        // 10^(10/10) = 10, 10^(20/10) = 100, 10^(30/10) = 1000
        // Using 0.001 tolerance due to f32 numerical precision in power operations
        assert!((values[0] - 10.0_f32).abs() < 0.001);
        assert!((values[1] - 100.0_f32).abs() < 0.001);
        assert!((values[2] - 1000.0_f32).abs() < 0.001);
    }

    #[test]
    fn test_linear_to_db_conversion() {
        let mut values: Vec<f32> = vec![10.0, 100.0, 1000.0]; // linear values
        linear_to_db_inplace(&mut values);

        // 10*log10(10) = 10, 10*log10(100) = 20, 10*log10(1000) = 30
        assert!((values[0] - 10.0_f32).abs() < 1e-5);
        assert!((values[1] - 20.0_f32).abs() < 1e-5);
        assert!((values[2] - 30.0_f32).abs() < 1e-5);
    }

    #[test]
    fn test_incidence_terms() {
        let alpha_rad: Vec<f32> = vec![0.5, 1.0, 1.5]; // radians
        let terms = IncidenceTerms::from_incidence_per_column(&alpha_rad);

        assert_eq!(terms.inv_sin.len(), 3);
        assert_eq!(terms.inv_cos.len(), 3);

        // Check that inv_sin * sin ≈ 1
        for (i, &alpha) in alpha_rad.iter().enumerate() {
            let sin_alpha = alpha.sin();
            assert!((terms.inv_sin[i] * sin_alpha - 1.0_f32).abs() < 1e-5);
        }
    }
}
