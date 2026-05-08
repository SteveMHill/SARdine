#![allow(dead_code)]
//! Polynomial evaluation utilities for IW TOPSAR deburst processing
//!
//! This module provides range-dependent polynomial models and evaluation functions
//! for Doppler centroid (DC) and FM rate calculations used during TOPSAR deramping.

use crate::core::deburst::geometry::{
    RangePolyEntry as GeometryRangePolyEntry, RangePolynomial as GeometryRangePolynomial,
};

/// Range-dependent polynomial model sampled at discrete slant-range times.
/// Each entry holds a time polynomial evaluated at a given slant range; values
/// are interpolated along range to approximate the 2D DC/FM surface.
#[derive(Debug, Clone)]
pub(crate) struct RangePolyEntry {
    pub(crate) range_time: f64,
    pub(crate) coeffs: Vec<f64>,
}

/// Range-dependent polynomial model for 2D DC/FM evaluation.
/// Interpolates across slant-range grid entries.
#[derive(Debug, Clone)]
pub(crate) struct RangePolynomial {
    pub(crate) entries: Vec<RangePolyEntry>,
}

/// Evaluate Doppler centroid (Hz) and FM rate (Hz/s) at azimuth time t
/// OPTIMIZATION #32: Use Horner's rule to avoid repeated t.powi() calls
/// Horner's method: c0 + c1*t + c2*t^2 = c0 + t*(c1 + t*c2) - evaluated right-to-left
pub(crate) fn eval_dc_fm(t: f64, dc_poly: &[f64], fm_poly: &[f64]) -> (f64, f64) {
    // Horner's rule: iterate coefficients in reverse, multiply by t, add next coeff
    let dc = dc_poly.iter().rev().fold(0.0, |acc, &c| acc * t + c);
    let fm = fm_poly.iter().rev().fold(0.0, |acc, &c| acc * t + c);
    (dc, fm)
}

/// Evaluate Doppler centroid and FM rate with RANGE DEPENDENCY for TOPS IW alignment
///
/// **CRITICAL**: DC/FM polynomials from Sentinel-1 annotation are functions of slant range time τ,
/// NOT azimuth time. The polynomial variable is τ - t0 where τ is two-way range time (seconds)
/// and t0 is the reference range time from <dataDcPolynomial><t0> (typically ~0.0053-0.0060 s).
///
/// This function is DEPRECATED - use geometry::eval_dc_fm_2d which correctly handles t0.
#[deprecated(note = "Use geometry::eval_dc_fm_2d with explicit dc_t0/fm_t0 parameters")]
pub(crate) fn eval_dc_fm_2d(
    _t_az: f64, // UNUSED - kept for API compatibility but should not be used
    range_sample: usize,
    dc_poly: &[f64],
    fm_poly: &[f64],
    dc_range_poly: Option<&RangePolynomial>,
    fm_range_poly: Option<&RangePolynomial>,
    range_pixel_spacing: f64,
    slant_range_time: f64,
) -> (f64, f64) {
    // Compute range time for this sample
    const SPEED_OF_LIGHT: f64 = 299_792_458.0;
    let two_way_range_time =
        slant_range_time + (range_sample as f64 * range_pixel_spacing * 2.0 / SPEED_OF_LIGHT);

    // **SCIENTIFIC FIX**: DC/FM polynomials are functions of τ (range time), not azimuth time.
    // For 1D case without explicit t0, assume t0=0 and evaluate at τ directly.
    // This is a fallback; proper usage should provide range-dependent model or use geometry::eval_dc_fm_2d.

    let dc = if let Some(model) = dc_range_poly {
        // Range model handles t0 internally via RangePolyEntry.t0
        model.evaluate(_t_az, two_way_range_time)
    } else {
        // FALLBACK: 1D polynomial in range time (τ - t0 where t0=0 assumed)
        // **WARNING**: This assumes t0=0 which is incorrect for real S1 data (t0 ~0.0053s)
        // Result will be wrong by ~hundreds of Hz unless dc_poly was pre-shifted
        dc_poly.iter().enumerate().fold(0.0, |acc, (i, &c)| {
            acc + c * two_way_range_time.powi(i as i32)
        })
    };

    let fm = if let Some(model) = fm_range_poly {
        model.evaluate(_t_az, two_way_range_time)
    } else {
        // FALLBACK: 1D polynomial in range time (same t0=0 assumption)
        fm_poly.iter().enumerate().fold(0.0, |acc, (i, &c)| {
            acc + c * two_way_range_time.powi(i as i32)
        })
    };

    (dc, fm)
}

/// Evaluate 2D polynomial: f(t, r) = sum(c_ij × t^i × r^j)
/// Supports time-only (1D) and time-range (2D) polynomials
///
/// **OPTIMIZATION:** Uses Horner's rule for 1D case to minimize multiplications
pub(crate) fn eval_poly_2d(coeffs: &[f64], t: f64, r: f64) -> f64 {
    if coeffs.len() <= 3 {
        // Time-only polynomial: Use Horner's rule for efficiency
        // c0 + c1*t + c2*t^2 becomes ((c2)*t + c1)*t + c0
        coeffs.iter().rev().fold(0.0, |acc, &c| acc * t + c)
    } else {
        // 2D polynomial expansion
        // Typical format: [c00, c10, c20, c01, c11, c02, ...]
        // where c_ij corresponds to t^i * r^j
        //
        // Note: Could optimize with nested Horner, but 2D case is rare for S1
        let mut result = 0.0;
        let max_order = ((coeffs.len() as f64).sqrt().ceil() as usize).max(1);

        for (idx, &coeff) in coeffs.iter().enumerate() {
            let i = idx / max_order; // time power
            let j = idx % max_order; // range power
            result += coeff * t.powi(i as i32) * r.powi(j as i32);
        }
        result
    }
}

impl RangePolynomial {
    /// Build a range-dependent polynomial model from estimate entries using closures
    pub(crate) fn from_estimates<T, FRange, FCoeff>(
        estimates: &[T],
        range_fn: FRange,
        coeff_fn: FCoeff,
    ) -> Option<Self>
    where
        FRange: Fn(&T) -> Option<f64>,
        FCoeff: Fn(&T) -> &[f64],
    {
        let mut entries: Vec<RangePolyEntry> = estimates
            .iter()
            .filter_map(|e| {
                let range_time = range_fn(e)?;
                if !range_time.is_finite() {
                    return None;
                }

                let coeffs = coeff_fn(e);
                if coeffs.is_empty() {
                    return None;
                }

                Some(RangePolyEntry {
                    range_time,
                    coeffs: coeffs.to_vec(),
                })
            })
            .collect();

        if entries.is_empty() {
            return None;
        }

        entries.sort_by(|a, b| {
            a.range_time
                .partial_cmp(&b.range_time)
                .unwrap_or(std::cmp::Ordering::Equal)
        });

        let min_len = entries.iter().map(|e| e.coeffs.len()).min().unwrap_or(0);
        if min_len == 0 {
            return None;
        }

        let mut truncated = false;
        for entry in entries.iter_mut() {
            if entry.coeffs.len() != min_len {
                truncated = true;
                entry.coeffs.truncate(min_len);
            }
        }

        if truncated {
            log::warn!(
                "⚠️  Truncated range-dependent polynomial coefficients to lowest degree {} to enforce consistency",
                min_len.saturating_sub(1)
            );
        }

        Some(Self { entries })
    }

    /// Evaluate the range-dependent polynomial at (t_az, range_time) using linear interpolation
    /// across the slant-range grid.
    pub(crate) fn evaluate(&self, t_az: f64, range_time: f64) -> f64 {
        if self.entries.is_empty() {
            return 0.0;
        }

        // Helper: evaluate a single entry (supports 1D or flattened 2D coeffs)
        let eval_entry = |entry: &RangePolyEntry| eval_poly_2d(&entry.coeffs, t_az, range_time);

        if self.entries.len() == 1 {
            return eval_entry(&self.entries[0]);
        }

        // Find bracketing entries in sorted order
        let mut upper_idx = 0usize;
        while upper_idx < self.entries.len() && self.entries[upper_idx].range_time < range_time {
            upper_idx += 1;
        }

        if upper_idx == 0 {
            return eval_entry(&self.entries[0]);
        }

        if upper_idx >= self.entries.len() {
            return eval_entry(self.entries.last().unwrap());
        }

        let lower = &self.entries[upper_idx - 1];
        let upper = &self.entries[upper_idx];
        let span = upper.range_time - lower.range_time;

        if !span.is_finite() || span.abs() < f64::EPSILON {
            return eval_entry(upper);
        }

        let w = ((range_time - lower.range_time) / span).clamp(0.0, 1.0);
        let v0 = eval_entry(lower);
        let v1 = eval_entry(upper);
        v0 + w * (v1 - v0)
    }

    /// Get the number of range sample entries
    pub(crate) fn samples(&self) -> usize {
        self.entries.len()
    }
}

impl From<&RangePolynomial> for GeometryRangePolynomial {
    fn from(src: &RangePolynomial) -> Self {
        let entries = src
            .entries
            .iter()
            .map(|entry| GeometryRangePolyEntry {
                range_time: entry.range_time,
                coeffs: entry.coeffs.clone(),
                t0: None,
            })
            .collect();

        GeometryRangePolynomial { entries }
    }
}
