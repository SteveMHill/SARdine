#![allow(dead_code, unused_variables)]
use crate::core::geometry::type_safe_units::{Hertz, Seconds};
/// DC and FM-rate provider trait for preventing zero polynomial fallbacks
///
/// This module implements the "DC/FM-rate provider" refactor from the practical refactors:
/// DC/FM-rate provider: one trait with 'get_dc(az_time)' & 'get_fm_rate(az_time)', mocked in tests.
use crate::types::{SarError, SarResult};
use std::collections::HashMap;

// ============================================================================
// OUT-OF-RANGE POLICY
// ============================================================================

/// Policy for handling azimuth times outside valid range
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum OutOfRangePolicy {
    /// Return error for out-of-range times (strict, recommended)
    Error,
    /// Clamp to nearest valid time
    Clamp,
    /// Use nearest burst polynomial (extrapolate)
    NearestBurst,
}

impl Default for OutOfRangePolicy {
    fn default() -> Self {
        OutOfRangePolicy::Clamp // Favor robustness by default
    }
}

// ============================================================================
// DC/FM-RATE PROVIDER TRAIT
// ============================================================================

/// Provider trait for Doppler Centroid (DC) and FM rate values
///
/// This trait abstracts the source of DC and FM rate data, allowing for:
/// - Real polynomial-based providers from annotation XML
/// - Mock providers for testing to avoid "zero polynomial" fallbacks
/// - Cached providers for performance optimization
/// - Interpolated providers for smooth temporal variation
pub trait DcFmRateProvider: Send + Sync {
    /// Get Doppler Centroid at the specified azimuth time
    ///
    /// # Arguments
    /// * `az_time` - Azimuth time in type-safe seconds
    ///
    /// # Returns
    /// * `Ok(Hertz)` - DC value in type-safe Hertz units
    /// * `Err(SarError)` - If DC cannot be computed for this time
    fn get_dc(&self, az_time: Seconds) -> SarResult<Hertz>;

    /// Get FM rate at the specified azimuth time
    ///
    /// # Arguments  
    /// * `az_time` - Azimuth time in type-safe seconds
    ///
    /// # Returns
    /// * `Ok(f64)` - FM rate in Hz/s (always f64 for precision)
    /// * `Err(SarError)` - If FM rate cannot be computed for this time
    fn get_fm_rate(&self, az_time: Seconds) -> SarResult<f64>;

    /// Get the valid time range for this provider
    ///
    /// # Returns
    /// * `(start_time, end_time)` - Valid azimuth time range in seconds
    fn get_time_range(&self) -> (Seconds, Seconds);

    /// Check if the provider can handle the given azimuth time
    ///
    /// # Arguments
    /// * `az_time` - Azimuth time to check
    ///
    /// # Returns
    /// * `true` if az_time is within valid range, `false` otherwise
    fn is_time_valid(&self, az_time: Seconds) -> bool {
        let (start, end) = self.get_time_range();
        az_time >= start && az_time <= end
    }

    /// Get provider type for debugging/logging
    fn provider_type(&self) -> &'static str;
}

// ============================================================================
// POLYNOMIAL-BASED PROVIDER (PRODUCTION)
// ============================================================================

/// Real DC/FM rate provider based on polynomial coefficients from annotation XML
#[derive(Debug, Clone)]
pub struct PolynomialDcFmProvider {
    /// DC polynomial coefficients per burst
    dc_polynomials: Vec<DcPolynomial>,
    /// FM rate polynomial coefficients per burst  
    fm_polynomials: Vec<FmPolynomial>,
    /// Burst time ranges for polynomial selection
    burst_time_ranges: Vec<(Seconds, Seconds)>,
    /// Overall valid time range
    time_range: (Seconds, Seconds),
    /// Out-of-range policy
    out_of_range_policy: OutOfRangePolicy,
}

/// DC polynomial representation
#[derive(Debug, Clone)]
pub struct DcPolynomial {
    /// Polynomial coefficients [c0, c1, c2, ...] where DC(t) = c0 + c1*t + c2*t^2 + ...
    coefficients: Vec<f64>,
    /// Reference time for polynomial (typically burst center time)
    reference_time: Seconds,
    /// Valid time range for this polynomial
    time_range: (Seconds, Seconds),
}

/// FM rate polynomial representation
#[derive(Debug, Clone)]
pub struct FmPolynomial {
    /// Polynomial coefficients [c0, c1, c2, ...] where FM(t) = c0 + c1*t + c2*t^2 + ...
    coefficients: Vec<f64>,
    /// Reference time for polynomial (typically burst center time)
    reference_time: Seconds,
    /// Valid time range for this polynomial
    time_range: (Seconds, Seconds),
}

impl DcPolynomial {
    /// Construct a DC polynomial with explicit reference time and validity window
    pub fn new(
        coefficients: Vec<f64>,
        reference_time: Seconds,
        time_range: (Seconds, Seconds),
    ) -> Self {
        Self {
            coefficients,
            reference_time,
            time_range,
        }
    }

    /// Expose coefficients for downstream logging/validation
    pub fn coefficients(&self) -> &[f64] {
        &self.coefficients
    }

    /// Expose reference time for downstream logging/validation
    pub fn reference_time(&self) -> Seconds {
        self.reference_time
    }
}

impl FmPolynomial {
    /// Construct an FM-rate polynomial with explicit reference time and validity window
    pub fn new(
        coefficients: Vec<f64>,
        reference_time: Seconds,
        time_range: (Seconds, Seconds),
    ) -> Self {
        Self {
            coefficients,
            reference_time,
            time_range,
        }
    }
}

impl PolynomialDcFmProvider {
    /// Create new polynomial provider from metadata
    pub fn new(
        dc_polynomials: Vec<DcPolynomial>,
        fm_polynomials: Vec<FmPolynomial>,
        burst_time_ranges: Vec<(Seconds, Seconds)>,
    ) -> SarResult<Self> {
        Self::with_policy(
            dc_polynomials,
            fm_polynomials,
            burst_time_ranges,
            OutOfRangePolicy::default(),
        )
    }

    /// Create new polynomial provider with custom out-of-range policy
    pub fn with_policy(
        dc_polynomials: Vec<DcPolynomial>,
        fm_polynomials: Vec<FmPolynomial>,
        burst_time_ranges: Vec<(Seconds, Seconds)>,
        out_of_range_policy: OutOfRangePolicy,
    ) -> SarResult<Self> {
        if dc_polynomials.is_empty() {
            return Err(SarError::InvalidMetadata(
                "DC polynomials cannot be empty".to_string(),
            ));
        }

        if dc_polynomials.len() != fm_polynomials.len() {
            return Err(SarError::InvalidMetadata(
                "Mismatch between DC and FM polynomial counts".to_string(),
            ));
        }

        if dc_polynomials.len() != burst_time_ranges.len() {
            return Err(SarError::InvalidMetadata(
                "Mismatch between polynomial count and burst time ranges".to_string(),
            ));
        }

        // Validate burst time ranges are ordered and non-overlapping
        for i in 0..burst_time_ranges.len() - 1 {
            let (_, end_i) = burst_time_ranges[i];
            let (start_next, _) = burst_time_ranges[i + 1];
            if end_i > start_next {
                return Err(SarError::InvalidMetadata(format!(
                    "Burst time ranges overlap: burst {} ends at {:.6}s, burst {} starts at {:.6}s",
                    i,
                    end_i.value(),
                    i + 1,
                    start_next.value()
                )));
            }
        }

        // Compute overall time range
        let mut min_time = Seconds::new(f64::INFINITY);
        let mut max_time = Seconds::new(f64::NEG_INFINITY);

        for &(start, end) in &burst_time_ranges {
            min_time = Seconds::new(min_time.value().min(start.value()));
            max_time = Seconds::new(max_time.value().max(end.value()));
        }

        let time_range = (min_time, max_time);

        Ok(Self {
            dc_polynomials,
            fm_polynomials,
            burst_time_ranges,
            time_range,
            out_of_range_policy,
        })
    }

    /// Find the appropriate polynomial index for the given azimuth time
    fn find_polynomial_index(&self, az_time: Seconds) -> SarResult<usize> {
        for (idx, &(start, end)) in self.burst_time_ranges.iter().enumerate() {
            if az_time >= start && az_time <= end {
                return Ok(idx);
            }
        }

        // If not in any burst range, find closest
        let mut closest_idx = 0;
        let mut min_distance = f64::INFINITY;

        for (idx, &(start, end)) in self.burst_time_ranges.iter().enumerate() {
            let burst_center = Seconds::new((start.value() + end.value()) / 2.0);
            let distance = (az_time.value() - burst_center.value()).abs();
            if distance < min_distance {
                min_distance = distance;
                closest_idx = idx;
            }
        }

        Ok(closest_idx)
    }

    /// Evaluate polynomial at given time using Horner's method
    ///
    /// Horner's method is numerically stable and efficient:
    /// p(x) = c0 + c1*x + c2*x² + c3*x³ + ...
    ///      = c0 + x*(c1 + x*(c2 + x*(c3 + ...)))
    ///
    /// Benefits over naive power method:
    /// - O(n) multiplications instead of O(n²)
    /// - Better numerical stability (fewer roundoff errors)
    /// - More cache-friendly
    fn evaluate_polynomial(coefficients: &[f64], time: Seconds, reference_time: Seconds) -> f64 {
        if coefficients.is_empty() {
            return 0.0;
        }

        let dt = time.value() - reference_time.value();

        // Horner's method: evaluate from highest degree down
        let mut result = coefficients[coefficients.len() - 1];
        for &coeff in coefficients[..coefficients.len() - 1].iter().rev() {
            result = result * dt + coeff;
        }

        result
    }
}

impl DcFmRateProvider for PolynomialDcFmProvider {
    fn get_dc(&self, az_time: Seconds) -> SarResult<Hertz> {
        // Check time validity first based on policy
        let time_to_use = match self.out_of_range_policy {
            OutOfRangePolicy::Error => {
                if !self.is_time_valid(az_time) {
                    return Err(SarError::InvalidInput(format!(
                        "Azimuth time {:.6} s is outside valid range [{:.6}, {:.6}] s",
                        az_time.value(),
                        self.time_range.0.value(),
                        self.time_range.1.value()
                    )));
                }
                az_time
            }
            OutOfRangePolicy::Clamp => {
                let clamped = az_time
                    .value()
                    .clamp(self.time_range.0.value(), self.time_range.1.value());
                Seconds::new(clamped)
            }
            OutOfRangePolicy::NearestBurst => az_time, // find_polynomial_index handles this
        };

        let poly_idx = self.find_polynomial_index(time_to_use)?;
        let dc_poly = &self.dc_polynomials[poly_idx];

        let dc_hz =
            Self::evaluate_polynomial(&dc_poly.coefficients, time_to_use, dc_poly.reference_time);

        Ok(Hertz::new(dc_hz))
    }

    fn get_fm_rate(&self, az_time: Seconds) -> SarResult<f64> {
        // Check time validity first based on policy
        let time_to_use = match self.out_of_range_policy {
            OutOfRangePolicy::Error => {
                if !self.is_time_valid(az_time) {
                    return Err(SarError::InvalidInput(format!(
                        "Azimuth time {:.6} s is outside valid range [{:.6}, {:.6}] s",
                        az_time.value(),
                        self.time_range.0.value(),
                        self.time_range.1.value()
                    )));
                }
                az_time
            }
            OutOfRangePolicy::Clamp => {
                let clamped = az_time
                    .value()
                    .clamp(self.time_range.0.value(), self.time_range.1.value());
                Seconds::new(clamped)
            }
            OutOfRangePolicy::NearestBurst => az_time, // find_polynomial_index handles this
        };

        let poly_idx = self.find_polynomial_index(time_to_use)?;
        let fm_poly = &self.fm_polynomials[poly_idx];

        let fm_rate =
            Self::evaluate_polynomial(&fm_poly.coefficients, time_to_use, fm_poly.reference_time);

        Ok(fm_rate)
    }

    fn get_time_range(&self) -> (Seconds, Seconds) {
        self.time_range
    }

    fn provider_type(&self) -> &'static str {
        "PolynomialDcFmProvider"
    }
}

// ============================================================================
// MOCK PROVIDER (TESTING)
// ============================================================================

/// Mock DC/FM rate provider for testing
///
/// This provider eliminates "zero polynomial" fallbacks in tests by providing
/// realistic, controllable DC and FM rate values.
#[derive(Debug, Clone)]
pub struct MockDcFmProvider {
    /// Fixed DC value to return
    dc_value: Hertz,
    /// Fixed FM rate value to return
    fm_rate_value: f64,
    /// Valid time range
    time_range: (Seconds, Seconds),
    /// Optional time-varying behavior
    time_variation: Option<MockTimeVariation>,
}

/// Mock time variation patterns
#[derive(Debug, Clone)]
pub enum MockTimeVariation {
    /// Linear variation: value = base + slope * (time - reference)
    Linear {
        slope_dc: f64,
        slope_fm: f64,
        reference_time: Seconds,
    },
    /// Sinusoidal variation for testing periodic effects
    Sinusoidal {
        amplitude_dc: f64,
        amplitude_fm: f64,
        period: Seconds,
        phase: f64,
    },
    /// Step function for testing discontinuities
    Step { steps: Vec<(Seconds, Hertz, f64)> }, // (time, dc_value, fm_value)
}

impl MockDcFmProvider {
    /// Create constant mock provider
    pub fn constant(dc_hz: f64, fm_rate: f64) -> Self {
        Self {
            dc_value: Hertz::new(dc_hz),
            fm_rate_value: fm_rate,
            time_range: (Seconds::new(0.0), Seconds::new(1000.0)), // Default large range
            time_variation: None,
        }
    }

    /// Create realistic Sentinel-1 mock provider
    pub fn sentinel1_realistic() -> Self {
        Self {
            dc_value: Hertz::new(-7193.46), // Typical S1 DC value
            fm_rate_value: -2313.13,        // Typical S1 FM rate
            time_range: (Seconds::new(0.0), Seconds::new(1000.0)),
            time_variation: None,
        }
    }

    /// Create linear variation mock provider
    pub fn linear_variation(
        base_dc: f64,
        base_fm: f64,
        dc_slope: f64,
        fm_slope: f64,
        reference_time: Seconds,
    ) -> Self {
        Self {
            dc_value: Hertz::new(base_dc),
            fm_rate_value: base_fm,
            time_range: (Seconds::new(0.0), Seconds::new(1000.0)),
            time_variation: Some(MockTimeVariation::Linear {
                slope_dc: dc_slope,
                slope_fm: fm_slope,
                reference_time,
            }),
        }
    }

    /// Create sinusoidal variation mock provider
    pub fn sinusoidal_variation(
        base_dc: f64,
        base_fm: f64,
        amplitude_dc: f64,
        amplitude_fm: f64,
        period: Seconds,
    ) -> Self {
        Self {
            dc_value: Hertz::new(base_dc),
            fm_rate_value: base_fm,
            time_range: (Seconds::new(0.0), Seconds::new(1000.0)),
            time_variation: Some(MockTimeVariation::Sinusoidal {
                amplitude_dc,
                amplitude_fm,
                period,
                phase: 0.0,
            }),
        }
    }

    /// Set custom time range
    pub fn with_time_range(mut self, start: Seconds, end: Seconds) -> Self {
        self.time_range = (start, end);
        self
    }

    /// Get DC value at time (considering variation)
    fn get_dc_at_time(&self, az_time: Seconds) -> Hertz {
        let base_dc = self.dc_value.value();

        match &self.time_variation {
            None => self.dc_value,
            Some(MockTimeVariation::Linear {
                slope_dc,
                reference_time,
                ..
            }) => {
                let dt = az_time.value() - reference_time.value();
                Hertz::new(base_dc + slope_dc * dt)
            }
            Some(MockTimeVariation::Sinusoidal {
                amplitude_dc,
                period,
                phase,
                ..
            }) => {
                let phase_rad =
                    2.0 * std::f64::consts::PI * az_time.value() / period.value() + phase;
                Hertz::new(base_dc + amplitude_dc * phase_rad.sin())
            }
            Some(MockTimeVariation::Step { steps }) => {
                // Find appropriate step (iterate in reverse to get last step <= az_time)
                for &(step_time, dc_val, _) in steps.iter().rev() {
                    if az_time >= step_time {
                        return dc_val;
                    }
                }
                self.dc_value // Default if before all steps
            }
        }
    }

    /// Get FM rate at time (considering variation)
    fn get_fm_rate_at_time(&self, az_time: Seconds) -> f64 {
        let base_fm = self.fm_rate_value;

        match &self.time_variation {
            None => base_fm,
            Some(MockTimeVariation::Linear {
                slope_fm,
                reference_time,
                ..
            }) => {
                let dt = az_time.value() - reference_time.value();
                base_fm + slope_fm * dt
            }
            Some(MockTimeVariation::Sinusoidal {
                amplitude_fm,
                period,
                phase,
                ..
            }) => {
                let phase_rad =
                    2.0 * std::f64::consts::PI * az_time.value() / period.value() + phase;
                base_fm + amplitude_fm * phase_rad.sin()
            }
            Some(MockTimeVariation::Step { steps }) => {
                // Find appropriate step (iterate in reverse to get last step <= az_time)
                for &(step_time, _, fm_val) in steps.iter().rev() {
                    if az_time >= step_time {
                        return fm_val;
                    }
                }
                base_fm // Default if before all steps
            }
        }
    }
}

impl DcFmRateProvider for MockDcFmProvider {
    fn get_dc(&self, az_time: Seconds) -> SarResult<Hertz> {
        if !self.is_time_valid(az_time) {
            return Err(SarError::InvalidInput(format!(
                "Azimuth time {:.6} s is outside valid range [{:.6}, {:.6}] s",
                az_time.value(),
                self.time_range.0.value(),
                self.time_range.1.value()
            )));
        }

        Ok(self.get_dc_at_time(az_time))
    }

    fn get_fm_rate(&self, az_time: Seconds) -> SarResult<f64> {
        if !self.is_time_valid(az_time) {
            return Err(SarError::InvalidInput(format!(
                "Azimuth time {:.6} s is outside valid range [{:.6}, {:.6}] s",
                az_time.value(),
                self.time_range.0.value(),
                self.time_range.1.value()
            )));
        }

        Ok(self.get_fm_rate_at_time(az_time))
    }

    fn get_time_range(&self) -> (Seconds, Seconds) {
        self.time_range
    }

    fn provider_type(&self) -> &'static str {
        "MockDcFmProvider"
    }
}

// ============================================================================
// CACHED PROVIDER (PERFORMANCE OPTIMIZATION)
// ============================================================================

/// Cached DC/FM rate provider for performance optimization
///
/// Wraps another provider and caches results to avoid repeated polynomial evaluations.
pub struct CachedDcFmProvider<T: DcFmRateProvider> {
    /// Underlying provider
    inner: T,
    /// DC cache: time -> DC value (uses signed i64 key for negative times)
    dc_cache: std::sync::RwLock<HashMap<i64, Hertz>>,
    /// FM rate cache: time -> FM rate (uses signed i64 key for negative times)
    fm_cache: std::sync::RwLock<HashMap<i64, f64>>,
    /// Cache time resolution (seconds)
    time_resolution: f64,
}

impl<T: DcFmRateProvider> CachedDcFmProvider<T> {
    /// Create new cached provider
    pub fn new(inner: T, time_resolution: f64) -> Self {
        Self {
            inner,
            dc_cache: std::sync::RwLock::new(HashMap::new()),
            fm_cache: std::sync::RwLock::new(HashMap::new()),
            time_resolution,
        }
    }

    /// Convert time to cache key using signed i64 with floor
    ///
    /// This correctly handles negative times (e.g., orbit-relative times)
    /// Uses floor instead of round for consistent bucketing
    fn time_to_key(&self, time: Seconds) -> i64 {
        (time.value() / self.time_resolution).floor() as i64
    }

    /// Clear caches
    pub fn clear_cache(&self) {
        if let Ok(mut dc_cache) = self.dc_cache.write() {
            dc_cache.clear();
        }
        if let Ok(mut fm_cache) = self.fm_cache.write() {
            fm_cache.clear();
        }
    }

    /// Get cache statistics
    pub fn cache_stats(&self) -> CacheStats {
        let dc_size = self.dc_cache.read().map(|c| c.len()).unwrap_or(0);
        let fm_size = self.fm_cache.read().map(|c| c.len()).unwrap_or(0);

        CacheStats {
            dc_entries: dc_size,
            fm_entries: fm_size,
            time_resolution: self.time_resolution,
        }
    }
}

/// Cache statistics
#[derive(Debug)]
pub struct CacheStats {
    pub dc_entries: usize,
    pub fm_entries: usize,
    pub time_resolution: f64,
}

impl<T: DcFmRateProvider> DcFmRateProvider for CachedDcFmProvider<T> {
    fn get_dc(&self, az_time: Seconds) -> SarResult<Hertz> {
        let key = self.time_to_key(az_time);

        // Try to get from cache first
        if let Ok(cache) = self.dc_cache.read() {
            if let Some(&cached_value) = cache.get(&key) {
                return Ok(cached_value);
            }
        }

        // Not in cache, compute and store
        let dc_value = self.inner.get_dc(az_time)?;

        if let Ok(mut cache) = self.dc_cache.write() {
            cache.insert(key, dc_value);
        }

        Ok(dc_value)
    }

    fn get_fm_rate(&self, az_time: Seconds) -> SarResult<f64> {
        let key = self.time_to_key(az_time);

        // Try to get from cache first
        if let Ok(cache) = self.fm_cache.read() {
            if let Some(&cached_value) = cache.get(&key) {
                return Ok(cached_value);
            }
        }

        // Not in cache, compute and store
        let fm_value = self.inner.get_fm_rate(az_time)?;

        if let Ok(mut cache) = self.fm_cache.write() {
            cache.insert(key, fm_value);
        }

        Ok(fm_value)
    }

    fn get_time_range(&self) -> (Seconds, Seconds) {
        self.inner.get_time_range()
    }

    fn provider_type(&self) -> &'static str {
        "CachedDcFmProvider"
    }
}

// ============================================================================
// PROVIDER FACTORY
// ============================================================================

/// Factory for creating DC/FM rate providers
pub struct DcFmProviderFactory;

impl DcFmProviderFactory {
    /// Create provider from metadata polynomials
    pub fn from_polynomials(
        dc_polynomials: Vec<DcPolynomial>,
        fm_polynomials: Vec<FmPolynomial>,
        burst_time_ranges: Vec<(Seconds, Seconds)>,
    ) -> SarResult<Box<dyn DcFmRateProvider>> {
        let provider =
            PolynomialDcFmProvider::new(dc_polynomials, fm_polynomials, burst_time_ranges)?;
        Ok(Box::new(provider))
    }

    /// Create mock provider for testing
    pub fn mock_constant(dc_hz: f64, fm_rate: f64) -> Box<dyn DcFmRateProvider> {
        Box::new(MockDcFmProvider::constant(dc_hz, fm_rate))
    }

    /// Create realistic Sentinel-1 mock provider
    pub fn mock_sentinel1() -> Box<dyn DcFmRateProvider> {
        Box::new(MockDcFmProvider::sentinel1_realistic())
    }

    /// Create cached provider
    pub fn with_cache<T: DcFmRateProvider + 'static>(
        provider: T,
        time_resolution: f64,
    ) -> Box<dyn DcFmRateProvider> {
        Box::new(CachedDcFmProvider::new(provider, time_resolution))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_mock_constant_provider() {
        let provider = MockDcFmProvider::constant(-7000.0, -2300.0);

        let az_time = Seconds::new(100.0);

        let dc = provider.get_dc(az_time).unwrap();
        assert_eq!(dc.value(), -7000.0);

        let fm_rate = provider.get_fm_rate(az_time).unwrap();
        assert_eq!(fm_rate, -2300.0);

        assert_eq!(provider.provider_type(), "MockDcFmProvider");
    }

    #[test]
    fn test_mock_linear_variation() {
        let provider = MockDcFmProvider::linear_variation(
            -7000.0,             // base DC
            -2300.0,             // base FM rate
            -10.0,               // DC slope (Hz/s)
            -5.0,                // FM slope (Hz/s²)
            Seconds::new(100.0), // reference time
        );

        // At reference time
        let dc_ref = provider.get_dc(Seconds::new(100.0)).unwrap();
        assert_eq!(dc_ref.value(), -7000.0);

        let fm_ref = provider.get_fm_rate(Seconds::new(100.0)).unwrap();
        assert_eq!(fm_ref, -2300.0);

        // At offset time
        let dc_offset = provider.get_dc(Seconds::new(110.0)).unwrap();
        assert_eq!(dc_offset.value(), -7000.0 + (-10.0 * 10.0)); // -7100.0

        let fm_offset = provider.get_fm_rate(Seconds::new(110.0)).unwrap();
        assert_eq!(fm_offset, -2300.0 + (-5.0 * 10.0)); // -2350.0
    }

    #[test]
    fn test_mock_time_validation() {
        let provider = MockDcFmProvider::constant(-7000.0, -2300.0)
            .with_time_range(Seconds::new(50.0), Seconds::new(150.0));

        // Valid time
        assert!(provider.is_time_valid(Seconds::new(100.0)));
        assert!(provider.get_dc(Seconds::new(100.0)).is_ok());

        // Invalid times
        assert!(!provider.is_time_valid(Seconds::new(25.0)));
        assert!(!provider.is_time_valid(Seconds::new(200.0)));
        assert!(provider.get_dc(Seconds::new(25.0)).is_err());
        assert!(provider.get_fm_rate(Seconds::new(200.0)).is_err());
    }

    #[test]
    fn test_polynomial_evaluation() {
        let coeffs = vec![100.0, -50.0, 25.0]; // 100 - 50*t + 25*t^2
        let time = Seconds::new(2.0);
        let ref_time = Seconds::new(0.0);

        let result = PolynomialDcFmProvider::evaluate_polynomial(&coeffs, time, ref_time);

        // Using Horner's method: 100 - 50*2 + 25*4 = 100 - 100 + 100 = 100
        assert_eq!(result, 100.0);
    }

    #[test]
    fn test_horner_vs_naive() {
        // Test that Horner's method gives same result as naive for small values
        let coeffs = vec![1.0, 2.0, 3.0, 4.0]; // 1 + 2*t + 3*t^2 + 4*t^3
        let time = Seconds::new(0.5);
        let ref_time = Seconds::new(0.0);

        let horner_result = PolynomialDcFmProvider::evaluate_polynomial(&coeffs, time, ref_time);

        // Expected: 1 + 2*0.5 + 3*0.25 + 4*0.125 = 1 + 1 + 0.75 + 0.5 = 3.25
        assert!((horner_result - 3.25).abs() < 1e-10);
    }

    #[test]
    fn test_out_of_range_policy_error() {
        let dc_poly = DcPolynomial {
            coefficients: vec![-7000.0],
            reference_time: Seconds::new(100.0),
            time_range: (Seconds::new(50.0), Seconds::new(150.0)),
        };

        let fm_poly = FmPolynomial {
            coefficients: vec![-2300.0],
            reference_time: Seconds::new(100.0),
            time_range: (Seconds::new(50.0), Seconds::new(150.0)),
        };

        let burst_ranges = vec![(Seconds::new(50.0), Seconds::new(150.0))];

        let provider = PolynomialDcFmProvider::with_policy(
            vec![dc_poly],
            vec![fm_poly],
            burst_ranges,
            OutOfRangePolicy::Error,
        )
        .unwrap();

        // Valid time - should work
        assert!(provider.get_dc(Seconds::new(100.0)).is_ok());

        // Out of range - should error
        assert!(provider.get_dc(Seconds::new(25.0)).is_err());
        assert!(provider.get_fm_rate(Seconds::new(200.0)).is_err());
    }

    #[test]
    fn test_out_of_range_policy_clamp() {
        let dc_poly = DcPolynomial {
            coefficients: vec![-7000.0, -10.0], // -7000 - 10*dt
            reference_time: Seconds::new(100.0),
            time_range: (Seconds::new(50.0), Seconds::new(150.0)),
        };

        let fm_poly = FmPolynomial {
            coefficients: vec![-2300.0],
            reference_time: Seconds::new(100.0),
            time_range: (Seconds::new(50.0), Seconds::new(150.0)),
        };

        let burst_ranges = vec![(Seconds::new(50.0), Seconds::new(150.0))];

        let provider = PolynomialDcFmProvider::with_policy(
            vec![dc_poly],
            vec![fm_poly],
            burst_ranges,
            OutOfRangePolicy::Clamp,
        )
        .unwrap();

        // Out of range below - should clamp to 50.0
        let dc_low = provider.get_dc(Seconds::new(25.0)).unwrap();
        let expected_low = -7000.0 - 10.0 * (50.0 - 100.0); // t=50 clamped
        assert_eq!(dc_low.value(), expected_low);

        // Out of range above - should clamp to 150.0
        let dc_high = provider.get_dc(Seconds::new(200.0)).unwrap();
        let expected_high = -7000.0 - 10.0 * (150.0 - 100.0); // t=150 clamped
        assert_eq!(dc_high.value(), expected_high);
    }

    #[test]
    fn test_cache_key_signed() {
        let mock = MockDcFmProvider::constant(-7000.0, -2300.0);
        let cached = CachedDcFmProvider::new(mock, 1.0); // 1s resolution

        // Test positive time
        let key_pos = cached.time_to_key(Seconds::new(100.5));
        assert_eq!(key_pos, 100);

        // Test negative time (orbit-relative)
        let key_neg = cached.time_to_key(Seconds::new(-50.3));
        assert_eq!(key_neg, -51); // floor of -50.3

        // Test zero
        let key_zero = cached.time_to_key(Seconds::new(0.0));
        assert_eq!(key_zero, 0);
    }

    #[test]
    fn test_mock_step_variation_reverse() {
        let steps = vec![
            (Seconds::new(0.0), Hertz::new(-7000.0), -2300.0),
            (Seconds::new(50.0), Hertz::new(-7100.0), -2350.0),
            (Seconds::new(100.0), Hertz::new(-7200.0), -2400.0),
        ];

        let provider = MockDcFmProvider {
            dc_value: Hertz::new(-6900.0),
            fm_rate_value: -2250.0,
            time_range: (Seconds::new(0.0), Seconds::new(200.0)),
            time_variation: Some(MockTimeVariation::Step { steps }),
        };

        // Before all steps
        let dc_early = provider.get_dc_at_time(Seconds::new(-10.0));
        assert_eq!(dc_early.value(), -6900.0); // base value

        // At first step
        let dc_0 = provider.get_dc_at_time(Seconds::new(0.0));
        assert_eq!(dc_0.value(), -7000.0);

        // Between steps
        let dc_75 = provider.get_dc_at_time(Seconds::new(75.0));
        assert_eq!(dc_75.value(), -7100.0); // second step

        // After last step
        let dc_150 = provider.get_dc_at_time(Seconds::new(150.0));
        assert_eq!(dc_150.value(), -7200.0); // third step
    }

    #[test]
    fn test_polynomial_provider_creation() {
        let dc_poly = DcPolynomial {
            coefficients: vec![-7000.0, -10.0],
            reference_time: Seconds::new(100.0),
            time_range: (Seconds::new(50.0), Seconds::new(150.0)),
        };

        let fm_poly = FmPolynomial {
            coefficients: vec![-2300.0, -5.0],
            reference_time: Seconds::new(100.0),
            time_range: (Seconds::new(50.0), Seconds::new(150.0)),
        };

        let burst_ranges = vec![(Seconds::new(50.0), Seconds::new(150.0))];

        let provider =
            PolynomialDcFmProvider::new(vec![dc_poly], vec![fm_poly], burst_ranges).unwrap();

        assert_eq!(provider.provider_type(), "PolynomialDcFmProvider");

        let (start, end) = provider.get_time_range();
        assert_eq!(start.value(), 50.0);
        assert_eq!(end.value(), 150.0);
    }

    #[test]
    fn test_cached_provider() {
        let mock = MockDcFmProvider::constant(-7000.0, -2300.0);
        let cached = CachedDcFmProvider::new(mock, 0.1); // 0.1s resolution

        let az_time = Seconds::new(100.0);

        // First call - not cached
        let dc1 = cached.get_dc(az_time).unwrap();
        assert_eq!(dc1.value(), -7000.0);

        // Second call - should be cached
        let dc2 = cached.get_dc(az_time).unwrap();
        assert_eq!(dc2.value(), -7000.0);

        let stats = cached.cache_stats();
        assert_eq!(stats.dc_entries, 1);
        assert_eq!(stats.time_resolution, 0.1);
    }

    #[test]
    fn test_provider_factory() {
        // Test mock creation
        let mock_provider = DcFmProviderFactory::mock_constant(-7000.0, -2300.0);
        let dc = mock_provider.get_dc(Seconds::new(100.0)).unwrap();
        assert_eq!(dc.value(), -7000.0);

        // Test Sentinel-1 mock
        let s1_provider = DcFmProviderFactory::mock_sentinel1();
        let dc_s1 = s1_provider.get_dc(Seconds::new(100.0)).unwrap();
        assert_eq!(dc_s1.value(), -7193.46);
    }
}
